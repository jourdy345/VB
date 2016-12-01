# Run all code below

# To fit a logistic GLMM in batch mode:
# VAbatchL(y,X,Z,C,vni,centering=2,updateW=0,fit=NULL,pv=FALSE) 

# To fit a logistic GLMM in minibatch mode:
# VAonlineL(y,X,Z,C,vni,centering=2,updateW=0,batchsize=48,pow=0.5,stab=0,fit=NULL,pv=FALSE,MAX) 

# centering = 2 : PNCP 
# centering = 0 : CP 
# centering = 1 : NCP
# updateW = 0 : do not further update tuning parameters
# updateW = 1 : further update tuning parameters 
# pv = FALSE : do not compute confict p-values
# pv = TRUE : compute confict p-values
# batchsize: size of minibatch (minibatches should differ in size by at most 1)
# vni : vector stating size of each cluster


library(mvtnorm)
library(Matrix)
library(fastGHQuad)

b <- function(x) log(1+exp(x))
b1 <- function(x) 1/(1+exp(-x))                             # 1st derivative of b
b2 <- function(x) exp(x)/(1+exp(x))^2                       # 2nd derivative of b
b3 <- function(x) exp(x)*(1-exp(x))/(1+exp(x))^3            # 3rd derivative of b
b4 <- function(x) exp(x)*(1-4*exp(x)+exp(2*x))/(1+exp(x))^4 # 4th derivative of b

# function for calculating sigmahat #
sg <- function(deriv,sigma,mu,x){
ins <- sigma*x+mu
if (deriv==0){SUB <- (b(ins)*b2(ins)-b1(ins)^2)/(b(ins)^2)}
if (deriv==1){SUB <- (b1(ins)*b3(ins)-b2(ins)^2)/(b1(ins)^2)}
if (deriv==2){SUB <- (b2(ins)*b4(ins)-b3(ins)^2)/(b2(ins)^2)}
(1-sigma^2*SUB)^(-0.5)}


 #------------------#
 # Common functions #
 #------------------#

#function to compute t(x)%*%A^(-1)%*%x (x is a vector)
quadinv <- function(x,A) sum(x*solve(A,x))

#function to compute tr(A^(-1)%*%B)
tri <- function(A,B) sum(diag(as.matrix(solve(A,B)))) 

#function to compute trace of A
tr <- function(A) sum(diag(as.matrix(A)))

CT <- function(vni){
ct <- matrix(0,n,2)
for (i in 1:n){ct[i,] <- c( sum(vni[1:i-1])+1,sum(vni[1:i]) )}
return(ct)}

LG <- function(i,vni,t1,t2,muhat,rule){
lg1 <- rep(0,vni[i])
lg2 <- rep(0,vni[i])
for (j in 1:vni[i]){
g1 <- function(x) b1(t1[j]+t2[j]*x)*dnorm(x)
g2 <- function(x) b2(t1[j]+t2[j]*x)*dnorm(x)
lg1[j] <- aghQuad(g1, muhat[[i]][j], sg(1,t2[j],t1[j], muhat[[i]][j]), rule)
lg2[j] <- aghQuad(g2, muhat[[i]][j], sg(2,t2[j],t1[j], muhat[[i]][j]), rule)}
list(lg1=lg1,lg2=lg2)
}

 #--------------------#
 # Lower Bound (PNCP) #
 #--------------------#

# LBfix: only need to compute once at the beginning of the algorithm

LBfix <- function(n,p,r,nuq,nu,S,sig0){
lb <- (sum(lgamma(0.5*(nuq-1:r+1))-lgamma(0.5*(nu-1:r+1)) )
       +0.5*n*r*log(2)+(p+n*r)/2-0.5*determinant(as.matrix(sig0))$modulus[1]
       +0.5*nu*determinant(as.matrix(S))$modulus[1] )
return(lb) }


# LBvar: used in minibatch mode

LBvar <- function(n,r,s,y,Z,V,Wtp,Wtpr,ct,mubq,sigbq,muaq,sigaq,nuq,Sq,S,sig0,rule,muhat){

lb <- 0
Vmb <- as.vector(V%*%mubq)
Vsb <- rowSums((V%*%sigbq)*V)
Wb <- Wtpr%*%mubq[1:s]
ab <- matrix(as.vector(t(muaq))- Wb,n,r,byrow=TRUE)
Ws <- array(as.vector(tcrossprod(sigbq[1:s,1:s],Wtpr)),dim=c(s,r,n))
t3 <- NULL

for (i in 1:n){
Zi <- Z[ct[i,1]:ct[i,2],,drop=FALSE]
t1 <- Vmb[ct[i,1]:ct[i,2]]+Zi%*%muaq[i,]
t2 <- sqrt(Vsb[ct[i,1]:ct[i,2]]+rowSums((Zi%*%sigaq[i,,])*Zi))
t3 <- c(t3,t1)

for (j in 1:vni[i]){
g <- function(x) b(t1[j]+t2[j]*x)*dnorm(x)
lb <- lb-aghQuad(g, muhat[[i]][j], sg(0,t2[j],t1[j],muhat[[i]][j]), rule)}

lb <- lb+0.5*determinant(as.matrix(sigaq[i,,]))$modulus[1] }

lb <- ( lb-0.5*nuq*determinant(as.matrix(Sq))$modulus[1]+0.5*nuq*r
        +0.5*determinant(as.matrix(sigbq))$modulus[1]+crossprod(y,t3)
        -0.5*(quadinv(mubq,sig0)+tri(sig0,sigbq)) 
        -0.5*nuq*tri(Sq,S+crossprod(ab)+colSums(sigaq)
        +matrix(aperm(Wtp,c(2,3,1)),r,n*s)%*%matrix(aperm(Ws,c(1,3,2)),n*s,r))
)
return(lb) }


# LBvar2: used in batch mode, only valid if nuq and Sq are the latest to be updated 
           
LBvar2 <- function(n,r,y,Z,V,ct,mubq,sigbq,muaq,sigaq,nuq,Sq,sig0,rule,muhat){
lb <- 0
Vmb <- as.vector(V%*%mubq)
Vsb <- rowSums((V%*%sigbq)*V)
t3 <- NULL   
for (i in 1:n){
Zi <- Z[ct[i,1]:ct[i,2],,drop=FALSE]
t1 <- Vmb[ct[i,1]:ct[i,2]]+Zi%*%muaq[i,]
t2 <- sqrt(Vsb[ct[i,1]:ct[i,2]]+rowSums((Zi%*%sigaq[i,,])*Zi))
t3 <- c(t3,t1)

for (j in 1:vni[i]){
g <- function(x) b(t1[j]+t2[j]*x)*dnorm(x)
lb <- lb-aghQuad(g, muhat[[i]][j], sg(0,t2[j],t1[j],muhat[[i]][j]), rule) }

lb <- lb+0.5*determinant(as.matrix(sigaq[i,,]))$modulus[1] }
lb <- ( lb-0.5*nuq*determinant(as.matrix(Sq))$modulus[1]
        +0.5*determinant(as.matrix(sigbq))$modulus[1]        
        -0.5*(quadinv(mubq,sig0)+tri(sig0,sigbq))+crossprod(y,t3) 
)
return(lb) }
         


 #---------------------#
 # Diagnostics pvalues #
 #---------------------#

Logisticpvalue <- function(n,r,y,Z,V,ct,mubq,sigbq,Wb,muaq,sigaq,nuq,Sq,rule,muhat){

Vsb <- rowSums((V%*%sigbq)*V)
Vmb <- as.vector(V%*%mubq)
if (r==1){pvalue <- matrix(0,n,3)}
if (r>1){pvalue <- matrix(0,n,r+1)}
for (i in 1:n){
yi <- y[ct[i,1]:ct[i,2]]
Zi <- Z[ct[i,1]:ct[i,2],,drop=FALSE]
t1 <- Vmb[ct[i,1]:ct[i,2]]+Zi%*%muaq[i,]
t2 <- sqrt(Vsb[ct[i,1]:ct[i,2]]+rowSums((Zi%*%sigaq[i,,])*Zi))
lg <- LG(i,vni,t1,t2,muhat,rule)
t21 <- crossprod(sqrt(lg$lg2)*Zi)
sig <- solve(t21)
mu <- muaq[i,]+crossprod(sig,crossprod(Zi,yi-lg$lg1))
mu <- Wb[((i-1)*r+1):(i*r)]-mu
sig <- Sq/nuq+sig
if (r==1){
pvalue[i,1] <- pnorm(0,mean=mu,sd=sqrt(sig),lower.tail=FALSE,log.p=FALSE)
pvalue[i,2] <- pnorm(0,mean=mu,sd=sqrt(sig),lower.tail=TRUE,log.p=FALSE)
pvalue[i,3] <- 2*min(pvalue[i,c(1,2)])
}
if (r>1) {
pvalue[i,1] <- pchisq(sum(mu*solve(sig,mu)),df=r,ncp=0,lower.tail=FALSE,log.p=FALSE)
for (j in 1:r){
pvalue[i,(j+1)] <- 2*min( pnorm(0,mean=mu[j],sd=sqrt(sig[j,j]),lower.tail=TRUE,log.p=FALSE),
                      pnorm(0,mean=mu[j],sd=sqrt(sig[j,j]),lower.tail=FALSE,log.p=FALSE) )
}}} 
return(pvalue) }


#####################       
# Logistic Batch VA #
#####################

VAbatchL <- function(y,X,Z,C,vni,centering=2,updateW=0,fit=NULL,pv=FALSE) {

n <- length(vni)
N <- sum(vni)
r <- dim(Z)[2]
s <- dim(C)[2]
p <- dim(X)[2]
ct <- CT(vni)
rule <- gaussHermiteData(10)
lbc <- c()

# Default Conjugate prior (Kass and natarajan) #
glmfit <- glm(y~X-1, family = binomial(link = "logit"))
linpd <- glmfit$linear.predictors
de <- exp(linpd)/(1+exp(linpd))^2
mn <- 0*diag(r)
for (i in 1:n) {
Zi <- matrix(Z[(sum(vni[1:i-1])+1):(sum(vni[1:i])),],vni[i],r)
dei <- de[(sum(vni[1:i-1])+1):(sum(vni[1:i]))]
mn <- (i-1)*mn/i+crossprod(dei*Zi,Zi)/i }
R <- solve(mn)
nu <- r
S <- r*R
sig0 <- 1000*diag(p)


# Initialization #  
V <- matrix(0,N,s)
W <- array(0,dim=c(n,r,r))
Wtp <- array(0,dim=c(n,r,s))
nuq <- nu+n
muaq <- matrix(0,n,r)
sigaq <- array(0,dim=c(n,r,r))

if (is.null(fit)){

mubq <- as.vector(coef(glmfit))
sigbq <- vcov(glmfit)
ransig <- S/r
for (i in 1:n){
Zi <- Z[ct[i,1]:ct[i,2],,drop=FALSE]
Ci <- C[((i-1)*r+1):(i*r),,drop=FALSE]
if (centering==1){W[i,,] <- diag(r)}
if (centering==2){
eta <- as.vector(glmfit$linear.predictors[ct[i,1]:ct[i,2]])
ff <- exp(eta)/(1+exp(eta))^2
W[i,,] <- solve(crossprod(ransig,crossprod(sqrt(ff)*Zi))+diag(r))}
Wtp[i,,] <- (diag(r)-W[i,,])%*%Ci  
muaq[i,] <- Wtp[i,,]%*%mubq[1:s]
sigaq[i,,] <- ransig
V[ct[i,1]:ct[i,2],] <- Zi%*%W[i,,]%*%Ci}

} else { 
           
mubq <- as.vector(fixef(fit))
sigbq <- as.matrix(vcov(fit))
if (r==1){ransig <- as.numeric(VarCorr(fit)[1:r,1])
} else {ransig <- diag(as.numeric(VarCorr(fit)[1:r,1]))}
ref <- as.matrix(ranef(fit),n,r)
for (i in 1:n){
Zi <- Z[ct[i,1]:ct[i,2],,drop=FALSE]
Xi <- X[ct[i,1]:ct[i,2],,drop=FALSE]
Ci <- C[((i-1)*r+1):(i*r),,drop=FALSE]
if (centering==1){W[i,,] <- diag(r)}
if (centering==2){
eta <- as.vector(Zi%*%ref[i,]+Xi%*%mubq)
ff <- exp(eta)/(1+exp(eta))^2
W[i,,] <- solve(crossprod(ransig,crossprod(sqrt(ff)*Zi))+diag(r))}
Wtp[i,,] <- (diag(r)-W[i,,])%*%Ci 
muaq[i,] <- Wtp[i,,]%*%mubq[1:s]+ref[i,]
sigaq[i,,] <- ransig
V[ct[i,1]:ct[i,2],] <- Zi%*%W[i,,]%*%Ci  }}

Wtpr <- matrix(aperm(Wtp,c(2,1,3)),n*r,s) 
Sq <- (nuq-r-1)*ransig
if (s<p) {
XG <- X[,(s+1):p,drop=FALSE]
V <- cbind(V,XG)}

lb1 <- LBfix(n,p,r,nuq,nu,S,sig0)
count <- 0
dif <- 10
lbold <- -1.0e17
record <- NULL
totaltime <- 0

while (dif > 1.0e-6) {

bef <- proc.time()
count <- count+1

# update W #
if (updateW ==1 & count > 1){
Vmb <- as.vector(V%*%mubq)
for (i in 1:n) {
Zi <- Z[ct[i,1]:ct[i,2],,drop=FALSE]
Ci <- C[((i-1)*r+1):(i*r),,drop=FALSE]
eta <- as.vector(Zi%*%muaq[i,]+Vmb[ct[i,1]:ct[i,2]])
ff <- exp(eta)/(1+exp(eta))^2
W[i,,] <- solve(crossprod(Sq/(nuq-r-1),crossprod(sqrt(ff)*Zi))+diag(r))
Wtp[i,,] <- (diag(r)-W[i,,])%*%Ci
V[ct[i,1]:ct[i,2],1:s] <- Zi%*%W[i,,]%*%Ci }
Wtpr <- matrix(aperm(Wtp,c(2,1,3)),n*r,s)
}

Vsb <- rowSums((V%*%sigbq)*V)
Vmb <- as.vector(V%*%mubq)

muhat <- list(NULL)
for (i in 1:n) {
muhati <- rep(0,vni[i])
Zi <- Z[ct[i,1]:ct[i,2],,drop=FALSE]
t1 <- Vmb[ct[i,1]:ct[i,2]]+Zi%*%muaq[i,]
t2 <- sqrt(Vsb[ct[i,1]:ct[i,2]]+rowSums((Zi%*%sigaq[i,,])*Zi))
for (j in 1:vni[i]){
g1 <- function(x) b1(t1[j]+t2[j]*x)*dnorm(x)
muhati[j] <- optimize(g1,c(-10,10),maximum=TRUE)$maximum}
muhat[[i]] <- muhati }

# update sigaq and muaq #
Wb <- Wtpr%*%mubq[1:s]
ab <- matrix(as.vector(t(muaq))- Wb,n,r,byrow=TRUE)
Sab <- -nuq*solve(Sq,t(ab))
for (i in 1:n){
yi <- y[ct[i,1]:ct[i,2]]
Zi <- Z[ct[i,1]:ct[i,2],,drop=FALSE]
t1 <- Vmb[ct[i,1]:ct[i,2]]+Zi%*%muaq[i,]
t2 <- sqrt(Vsb[ct[i,1]:ct[i,2]]+ rowSums((Zi%*%sigaq[i,,])*Zi))
lg <- LG(i,vni,t1,t2,muhat,rule)
t21 <- crossprod(sqrt(lg$lg2)*Zi)
sigaq[i,,] <- solve(nuq*solve(Sq)+t21)
muaq[i,] <- muaq[i,]+crossprod(sigaq[i,,],crossprod(Zi,yi-lg$lg1)+Sab[,i])}

# update sigbq and mubq #
t11 <- NULL
t12 <- NULL
ab <- matrix(as.vector(t(muaq))- Wb,n,r,byrow=TRUE)
for (i in 1:n){
Zi <- Z[ct[i,1]:ct[i,2],,drop=FALSE]
t1 <- Vmb[ct[i,1]:ct[i,2]]+Zi%*%muaq[i,]
t2 <- sqrt(Vsb[ct[i,1]:ct[i,2]] + rowSums((Zi%*%sigaq[i,,])*Zi))
lg <- LG(i,vni,t1,t2,muhat,rule)
t11 <- c(t11,lg$lg1)
t12 <- c(t12,lg$lg2)}
t15 <- 0*diag(p)
t14 <- array(as.vector(solve(Sq,matrix(aperm(Wtp,c(2,3,1)),r,n*s))),dim=c(r,s,n))
t15[1:s,1:s] <- crossprod(Wtpr,matrix(aperm(t14,c(1,3,2)),n*r,s))
sigbq <- solve(solve(sig0)+nuq*t15+crossprod(sqrt(t12)*V))
mubq <- mubq+crossprod(sigbq,-solve(sig0,mubq)+crossprod(V,y-t11)
         +nuq*c(crossprod(Wtpr,as.vector(solve(Sq,t(ab)))),rep(0,p-s)))

# update Sq #
Wb <- Wtpr%*%mubq[1:s]
Ws <- array(as.vector(tcrossprod(sigbq[1:s,1:s],Wtpr)),dim=c(s,r,n))
ab <- matrix(as.vector(t(muaq))- Wb,n,r,byrow=TRUE)
Sq <- S+matrix(aperm(Wtp,c(2,3,1)),r,n*s)%*%matrix(aperm(Ws,c(1,3,2)),n*s,r)+
      crossprod(ab)+colSums(sigaq)

lb <- lb1+LBvar2(n,r,y,Z,V,ct,mubq,sigbq,muaq,sigaq,nuq,Sq,sig0,rule,muhat)
lbc <- c(lbc,lb)
dif <- abs((lb-lbold)/lb)
lbold <- lb
aft <- proc.time()
time <- (aft-bef)[1]
totaltime <- totaltime+time
cat(count,totaltime,lb,dif,mubq,Sq[1,1],"\n")
record <- rbind(record, c(count,totaltime,lb))  }

if (pv==TRUE) {
bef <- proc.time()
pvalue <- Logisticpvalue(n,r,y,Z,V,ct,mubq,sigbq,Wb,muaq,sigaq,nuq,Sq,rule,muhat)
aft <- proc.time()
time <- (aft-bef)[1]
record <- rbind(record, c(time,0,0))  
} else {pvalue <- NULL}

list(count=count,mubq=mubq,sigbq=sigbq,muaq=muaq,sigaq=sigaq,nuq=nuq,Sq=Sq,
lb=lb,lbc=lbc,W=W,pvalue=pvalue,totaltime=totaltime,record=record)}


######################       
# Logistic Online VA #
######################

VAonlineL <- function(y,X,Z,C,vni,centering=2,updateW=0,batchsize=48,pow=0.5,stab=0,fit=NULL,pv=FALSE,MAX) {

n <- length(vni)
N <- sum(vni)
r <- dim(Z)[2]
s <- dim(C)[2]
p <- dim(X)[2]
ct <- CT(vni)
rule <- gaussHermiteData(10)

# Default Conjugate prior (Kass and natarajan) #
glmfit <- glm(y~X-1, family = binomial(link = "logit"))
linpd <- glmfit$linear.predictors
de <- exp(linpd)/(1+exp(linpd))^2
mn <- 0*diag(r)
for (i in 1:n) {
Zi <- matrix(Z[(sum(vni[1:i-1])+1):(sum(vni[1:i])),],vni[i],r)
dei <- de[(sum(vni[1:i-1])+1):(sum(vni[1:i]))]
mn <- (i-1)*mn/i+crossprod(dei*Zi,Zi)/i }
R <- solve(mn)
nu <- r
S <- r*R
sig0 <- 1000*diag(p)

# Initialization #  
V <- matrix(0,N,s)
W <- array(0,dim=c(n,r,r))
Wtp <- array(0,dim=c(n,r,s))
nuq <- nu+n
muaq <- matrix(0,n,r)
sigaq <- array(0,dim=c(n,r,r))

if (is.null(fit)){

mubq <- as.vector(coef(glmfit))
sigbq <- vcov(glmfit)
ransig <- S/r
for (i in 1:n){
Zi <- Z[ct[i,1]:ct[i,2],,drop=FALSE]
Ci <- C[((i-1)*r+1):(i*r),,drop=FALSE]
if (centering==1){W[i,,] <- diag(r)}
if (centering==2){
eta <- as.vector(glmfit$linear.predictors[ct[i,1]:ct[i,2]])
ff <- exp(eta)/(1+exp(eta))^2
W[i,,] <- solve(crossprod(ransig,crossprod(sqrt(ff)*Zi))+diag(r))}
Wtp[i,,] <- (diag(r)-W[i,,])%*%Ci  
muaq[i,] <- Wtp[i,,]%*%mubq[1:s]
sigaq[i,,] <- ransig
V[ct[i,1]:ct[i,2],] <- Zi%*%W[i,,]%*%Ci}

} else { 
           
mubq <- as.vector(fixef(fit))
sigbq <- as.matrix(vcov(fit))
if (r==1){ransig <- as.numeric(VarCorr(fit)[1:r,1])
} else {ransig <- diag(as.numeric(VarCorr(fit)[1:r,1]))}
ref <- as.matrix(ranef(fit),n,r)
for (i in 1:n){
Zi <- Z[ct[i,1]:ct[i,2],,drop=FALSE]
Xi <- X[ct[i,1]:ct[i,2],,drop=FALSE]
Ci <- C[((i-1)*r+1):(i*r),,drop=FALSE]
if (centering==1){W[i,,] <- diag(r)}
if (centering==2){
eta <- as.vector(Zi%*%ref[i,]+Xi%*%mubq)
ff <- exp(eta)/(1+exp(eta))^2
W[i,,] <- solve(crossprod(ransig,crossprod(sqrt(ff)*Zi))+diag(r))}
Wtp[i,,] <- (diag(r)-W[i,,])%*%Ci 
muaq[i,] <- Wtp[i,,]%*%mubq[1:s]+ref[i,]
sigaq[i,,] <- ransig
V[ct[i,1]:ct[i,2],] <- Zi%*%W[i,,]%*%Ci  }}

Wtpr <- matrix(aperm(Wtp,c(2,1,3)),n*r,s) 
Sq <- (nuq-r-1)*ransig
if (s<p) {
XG <- X[,(s+1):p,drop=FALSE]
V <- cbind(V,XG)}

lb1 <- LBfix(n,p,r,nuq,nu,S,sig0)
swp <- 0
record <- NULL
totaltime <- 0
numold <- 0
perdif <- 1

muhat <- list(NULL)
max <- floor(n/batchsize)
re <- n-batchsize*max
samplesizes <- c(rep(batchsize,max-re),rep(batchsize+1,re))

while (perdif > 1.0e-3) {
ranseq <- sample(seq(1,n,1),n) 
swp <- swp+1

# update W #
if (updateW ==1 & swp > 1){
Vmb <- as.vector(V%*%mubq)
for (i in 1:n) {
Zi <- Z[ct[i,1]:ct[i,2],,drop=FALSE]
Ci <- C[((i-1)*r+1):(i*r),,drop=FALSE]
eta <- as.vector(Zi%*%muaq[i,]+Vmb[ct[i,1]:ct[i,2]])
ff <- exp(eta)/(1+exp(eta))^2
W[i,,] <- solve(crossprod(Sq/(nuq-r-1),crossprod(sqrt(ff)*Zi))+diag(r))
Wtp[i,,] <- (diag(r)-W[i,,])%*%Ci
V[ct[i,1]:ct[i,2],1:s] <- Zi%*%W[i,,]%*%Ci }
Wtpr <- matrix(aperm(Wtp,c(2,1,3)),n*r,s)
}

for (it in 1:max){
bef1 <- proc.time()
step <- (swp+(it-1)/max+stab)^(-pow)

mb <- sort(ranseq[(sum(samplesizes[1:it-1])+1):(sum(samplesizes[1:it]))])
nS <- samplesizes[it]
mba <- which(rep(1:n,vni) %in% mb)
vnimb <- vni[mb]
Vpar <- V[mba,,drop=FALSE]
Wtpmb <- Wtp[mb,,,drop=FALSE]

Vsb <- rowSums((Vpar%*%sigbq)*Vpar)
Vmb <- as.vector(Vpar%*%mubq)
Wtpmbr <- matrix(aperm(Wtpmb,c(2,1,3)),length(mb)*r,s)
Wb <- Wtpmbr%*%mubq[1:s]

ii <- 0
for (i in mb){
ii <- ii+1
muhati <- rep(0,vni[i])
Zi <- Z[ct[i,1]:ct[i,2],,drop=FALSE]
t1 <- Vmb[(sum(vnimb[1:ii-1])+1):(sum(vnimb[1:ii]))]+Zi%*%muaq[i,]
t2 <- sqrt(Vsb[(sum(vnimb[1:ii-1])+1):(sum(vnimb[1:ii]))]+rowSums((Zi%*%sigaq[i,,])*Zi))
for (j in 1:vni[i]){
g1 <- function(x) b1(t1[j]+t2[j]*x)*dnorm(x)
muhati[j] <- optimize(g1,c(-10,10),maximum = TRUE)$maximum}
muhat[[i]] <- muhati }


# update sigaq and muaq #
num <- 0
localdiff <- 1
muaqold <- muaq
while (localdiff>=0.05 & num < MAX) {
num <- num+1
ab <- matrix(as.vector(t(muaq[mb,]))- Wb,length(mb),r,byrow=TRUE)
Sab <- -nuq*solve(Sq,t(ab))
ii <- 0
for (i in mb){
ii <- ii+1
yi <- y[ct[i,1]:ct[i,2]]
Zi <- Z[ct[i,1]:ct[i,2],,drop=FALSE]
t1 <- Vmb[(sum(vnimb[1:ii-1])+1):(sum(vnimb[1:ii]))]+Zi%*%muaq[i,]
t2 <- sqrt(Vsb[(sum(vnimb[1:ii-1])+1):(sum(vnimb[1:ii]))]+rowSums((Zi%*%sigaq[i,,])*Zi))
lg <- LG(i,vni,t1,t2,muhat,rule)
t21 <- crossprod(matrix(rep(sqrt(lg$lg2),r),vni[i],r)*Zi)
sigaq[i,,] <- solve(nuq*solve(Sq)+t21)
muaq[i,] <- muaq[i,]+crossprod(sigaq[i,,],crossprod(Zi,yi-lg$lg1)+Sab[,ii])}
localdiff <-norm(as.matrix((muaq-muaqold)[mb,]),type='F')/ norm(as.matrix(muaq[mb,]),type='F') # Euclidean norm
muaqold <- muaq}


# update sigbq and mubq #
t11 <- NULL
t12 <- NULL
ii <- 0
ab <- matrix(as.vector(t(muaq[mb,]))- Wb,length(mb),r,byrow=TRUE)
for (i in mb){
ii <- ii+1
yi <- y[ct[i,1]:ct[i,2]]
Zi <- Z[ct[i,1]:ct[i,2],,drop=FALSE]
t1 <- Vmb[(sum(vnimb[1:ii-1])+1):(sum(vnimb[1:ii]))]+Zi%*%muaq[i,]
t2 <- sqrt(Vsb[(sum(vnimb[1:ii-1])+1):(sum(vnimb[1:ii]))]+rowSums((Zi%*%sigaq[i,,])*Zi))
lg <- LG(i,vni,t1,t2,muhat,rule)
t11 <- c(t11,lg$lg1)
t12 <- c(t12,lg$lg2) }
t15 <- 0*diag(p)
t14 <- array(as.vector(solve(Sq,matrix(aperm(Wtpmb,c(2,3,1)),r,length(mb)*s))),dim=c(r,s,length(mb)))
t15[1:s,1:s] <- crossprod(Wtpmbr,matrix(aperm(t14,c(1,3,2)),length(mb)*r,s))
sigbq <- solve((1-step)*solve(sigbq)+step*(solve(sig0)
         +n/nS*(nuq*t15+crossprod(sqrt(t12)*Vpar))))
mubq <- mubq+step*crossprod(sigbq,-solve(sig0,mubq)+n/nS*(crossprod(Vpar,y[mba]-t11)
        +nuq*c(crossprod(Wtpmbr, as.vector(solve(Sq,t(ab)))),rep(0,p-s))))

# update Sq #
Wb <- Wtpmbr%*%mubq[1:s]
ab <- matrix(as.vector(t(muaq[mb,]))- Wb,length(mb),r,byrow=TRUE)
Ws <- array(as.vector(tcrossprod(sigbq[1:s,1:s],Wtpmbr)),dim=c(s,r,length(mb)))
if (r==1) {t6 <- sum(sigaq[mb,,])} else {t6 <- colSums(sigaq[mb,,])}
Sq <- (1-step)*Sq+step*( S+ n/nS*(matrix(aperm(Wtpmb,c(2,3,1)),r,length(mb)*s)%*%
       matrix(aperm(Ws,c(1,3,2)),length(mb)*s,r)+t6+crossprod(ab)) )

aft1 <- proc.time()
time <- (aft1-bef1)[1]
totaltime <- totaltime+time
nummax <- max(num,numold)
numold <- num
}

bef <- proc.time()
lb <- lb1+LBvar(n,r,s,y,Z,V,Wtp,Wtpr,ct,mubq,sigbq,muaq,sigaq,nuq,Sq,S,sig0,rule,muhat)
aft <- proc.time()
time <- (aft-bef)[1]
totaltime <- totaltime+time
cat(swp,totaltime,lb,Sq[1,1],nummax,"\n")

record <- rbind(record, c(swp,step,totaltime,lb,nummax))
if (swp>1) {
perdif <- (record[swp,4]-record[(swp-1),4])/abs(record[swp,4])
cat(swp,perdif,"\n")}
}

lbold <- lb
dif <- 1
count <- 0

while (dif > 1.0e-6) {

bef <- proc.time()
swp <- swp+1

# update W #
if (updateW ==1){
Vmb <- as.vector(V%*%mubq)
for (i in 1:n) {
Zi <- Z[ct[i,1]:ct[i,2],,drop=FALSE]
Ci <- C[((i-1)*r+1):(i*r),,drop=FALSE]
eta <- as.vector(Zi%*%muaq[i,]+Vmb[ct[i,1]:ct[i,2]])
ff <- exp(eta)/(1+exp(eta))^2
W[i,,] <- solve(crossprod(Sq/(nuq-r-1),crossprod(sqrt(ff)*Zi))+diag(r))
Wtp[i,,] <- (diag(r)-W[i,,])%*%Ci
V[ct[i,1]:ct[i,2],1:s] <- Zi%*%W[i,,]%*%Ci }
Wtpr <- matrix(aperm(Wtp,c(2,1,3)),n*r,s)
}

Vsb <- rowSums((V%*%sigbq)*V)
Vmb <- as.vector(V%*%mubq)

for (i in 1:n) {
muhati <- rep(0,vni[i])
Zi <- Z[ct[i,1]:ct[i,2],,drop=FALSE]
t1 <- Vmb[ct[i,1]:ct[i,2]]+Zi%*%muaq[i,]
t2 <- sqrt(Vsb[ct[i,1]:ct[i,2]]+rowSums((Zi%*%sigaq[i,,])*Zi))
for (j in 1:vni[i]){
g1 <- function(x) b1(t1[j]+t2[j]*x)*dnorm(x)
muhati[j] <- optimize(g1,c(-10,10),maximum=TRUE)$maximum}
muhat[[i]] <- muhati }

# update sigaq and muaq #
Wb <- Wtpr%*%mubq[1:s]
ab <- matrix(as.vector(t(muaq))- Wb,n,r,byrow=TRUE)
Sab <- -nuq*solve(Sq,t(ab))
for (i in 1:n){
yi <- y[ct[i,1]:ct[i,2]]
Zi <- Z[ct[i,1]:ct[i,2],,drop=FALSE]
t1 <- Vmb[ct[i,1]:ct[i,2]]+Zi%*%muaq[i,]
t2 <- sqrt(Vsb[ct[i,1]:ct[i,2]]+ rowSums((Zi%*%sigaq[i,,])*Zi))
lg <- LG(i,vni,t1,t2,muhat,rule)
t21 <- crossprod(sqrt(lg$lg2)*Zi)
sigaq[i,,] <- solve(nuq*solve(Sq)+t21)
muaq[i,] <- muaq[i,]+crossprod(sigaq[i,,],crossprod(Zi,yi-lg$lg1)+Sab[,i])}

# update sigbq and mubq #
t11 <- NULL
t12 <- NULL
ab <- matrix(as.vector(t(muaq))- Wb,n,r,byrow=TRUE)
for (i in 1:n){
Zi <- Z[ct[i,1]:ct[i,2],,drop=FALSE]
t1 <- Vmb[ct[i,1]:ct[i,2]]+Zi%*%muaq[i,]
t2 <- sqrt(Vsb[ct[i,1]:ct[i,2]] + rowSums((Zi%*%sigaq[i,,])*Zi))
lg <- LG(i,vni,t1,t2,muhat,rule)
t11 <- c(t11,lg$lg1)
t12 <- c(t12,lg$lg2)}
t15 <- 0*diag(p)
t14 <- array(as.vector(solve(Sq,matrix(aperm(Wtp,c(2,3,1)),r,n*s))),dim=c(r,s,n))
t15[1:s,1:s] <- crossprod(Wtpr,matrix(aperm(t14,c(1,3,2)),n*r,s))
sigbq <- solve(solve(sig0)+nuq*t15+crossprod(sqrt(t12)*V))
mubq <- mubq+crossprod(sigbq,-solve(sig0,mubq)+crossprod(V,y-t11)
         +nuq*c(crossprod(Wtpr,as.vector(solve(Sq,t(ab)))),rep(0,p-s)))

# update Sq #
Wb <- Wtpr%*%mubq[1:s]
Ws <- array(as.vector(tcrossprod(sigbq[1:s,1:s],Wtpr)),dim=c(s,r,n))
ab <- matrix(as.vector(t(muaq))- Wb,n,r,byrow=TRUE)
Sq <- S+matrix(aperm(Wtp,c(2,3,1)),r,n*s)%*%matrix(aperm(Ws,c(1,3,2)),n*s,r)+
      crossprod(ab)+colSums(sigaq)


lb <- lb1+LBvar2(n,r,y,Z,V,ct,mubq,sigbq,muaq,sigaq,nuq,Sq,sig0,rule,muhat)
dif <- abs((lb-lbold)/lb)
lbold <- lb
aft <- proc.time()
time <- (aft-bef)[1]
totaltime <- totaltime+time
cat(swp,totaltime,lb,dif,mubq,Sq[1,1],"\n")
record <- rbind(record, c(swp,0,totaltime,lb,0))
}

if (pv==TRUE) {
bef <- proc.time()
pvalue <- Logisticpvalue(n,r,y,Z,V,ct,mubq,sigbq,Wb,muaq,sigaq,nuq,Sq,rule,muhat)
aft <- proc.time()
time <- (aft-bef)[1]
record <- rbind(record, c(time,0,0,0))  
} else {pvalue <- NULL}

list(mubq=mubq,sigbq=sigbq,muaq=muaq,sigaq=sigaq,nuq=nuq,Sq=Sq,
lb=lb,W=W,pvalue=pvalue,totaltime=totaltime,record=record)}



###################
# Bristol example #
###################


set.seed(7)
vni <- c(143, 187, 323, 122, 164, 405, 239, 482, 195, 177, 581, 301)
deaths <- c(41,25,24,23,25,42,24,53,26,25,58,31)
n <- length(vni)
N <- sum(vni)
y <- rep(0,N)
ynew <- NULL
for (i in 1:n) {
yi <- y[(sum(vni[1:i-1])+1):(sum(vni[1:i]))]
sam <- sample(seq(1,vni[i],1),deaths[i])
yi[sam] <- rep(1,deaths[i])
ynew <- c(ynew,yi)
}
y <- ynew
Z <- matrix(1,N,1)
X <- Z
C <- matrix(1,n,1)
id <- rep(seq(1,n,1),vni)
r <- dim(Z)[2]
XR <- Z

fitVA <- VAbatchL(y,X,Z,C,vni,centering=2,updateW=1,fit=NULL,pv=TRUE)

