library(data.table) # replace aggregate function
library(mvtnorm)    # for multivariate normal density
library(MASS)
library(bayesm)     # MCMC for mixed MNL models
library(Ecdat)      # for economic data


# prior for zeta: N(m0,S0)
# variational posterior for Omega: IW(omeg,Ups)


#############################
# Fixed part of lower bound #   # same for all approximations 
#############################

LBFIX <- function(H, K, omeg, b, nu, A, S0){
( (H + 1 + omeg + omeg * log(2))*K/2 + K*(nu+K-1)*log(nu)/2 - 0.5 * determinant(S0)$modulus[1] - K*lgamma(0.5)
  + sum( lgamma((omeg-1:K+1)/2) - lgamma((nu-1:K+K)/2)+ lgamma(b) - log(A) ))}


#########################
# Laplace approximation #
#########################

# Note: OU = omeg*solve(Ups), cOU = chol(OU)

# Function for maximizing beta_h #
# yX:  K by 1, sum over h, j and t
# yXH: H by K, sum over j and t

fbeta <- function(betah, yXh, Xh, muz, cOU, Th, J) {
mp <- matrix(Xh %*% betah, Th, J, byrow=TRUE)
ma <- apply(mp, 1, max)
sum(yXh * betah) - sum(log(rowSums(exp(mp-ma))) + ma) - 0.5*crossprod( cOU %*% (betah-muz) ) }

# gradient for maximizing beta_h #

grbeta <- function(betah, yXh, Xh, muz, cOU, Th, J) {
mp <- matrix(Xh %*% betah, Th, J, byrow=TRUE)
mp <- exp(mp-apply(mp, 1, max))
Xp <- Xh * as.vector(t(mp/rowSums(mp)))
as.vector(yXh - colSums(Xp) - crossprod(cOU, cOU %*% (betah-muz))) }

# Hessian #

Hessbeta <- function(betah, yXh, Xh, cOU, Th, J) {
mp <- matrix(Xh %*% betah, Th, J, byrow=TRUE)
mp <- exp(mp-apply(mp, 1, max))
Xp <- Xh * as.vector(t(mp/rowSums(mp)))
grp <- rep(1:Th,rep(J,Th))
DT <- data.table(Xp,grp)
temp <- as.matrix(DT[,lapply(.SD,sum),by=grp])[,2:(K+1),drop=FALSE] # H by K
H <- - crossprod(Xp,Xh) + crossprod(temp) - crossprod(cOU) 
return(H)}

LBVAR_Laplace <- function(H, K, Sig, muz, Sigz, b, c, omeg, Ups, m0, S0i, cS0i, OU, fvalue){
 L <- sum(fvalue) - H * K / 2 + 0.5 * (- H * sum(Sigz*OU) - crossprod(cS0i %*% (muz-m0)) - sum(S0i * Sigz) 
      - omeg * determinant(Ups)$modulus[1] + determinant(Sigz)$modulus[1]) - sum(b * log(c)) 
 for (h in 1:H){ L <- L + 0.5*determinant(as.matrix(Sig[h,,]))$modulus[1] }
 return(L) 
}


#########
# NCVMP #
#########

LBVAR_NCVMP <- function(y, X, H, J, K, T, sT, Tc, TJc, mu, Sig, muz, Sigz, b, c, omeg, Ups,
                        m0, S0i, cS0i, OU, cOU){

 mm <- mu - matrix(rep(muz,H), H, K, byrow=TRUE)
 L <- 0.5 * ( - sum(colSums(Sig) * OU) - sum((tcrossprod(mm,cOU))^2)  
      - H * sum(Sigz*OU) - crossprod(cS0i %*% (muz-m0)) - sum(S0i * Sigz) 
      - omeg * determinant(Ups)$modulus[1] + determinant(Sigz)$modulus[1]) - sum(b * log(c))
     
 mp <- rowSums(X * mu[rep(1:H,J*T),])  # x_{htj}^T %*% mu_h
 mp1 <- matrix(mp, sT, J, byrow=TRUE)
 ma <- apply(mp1, 1, max)
 mp2 <- exp(mp1 - ma)
 mp3 <- rowSums(mp2)
 Xp <- X * as.vector(t(mp2/mp3))   # each row is p_{htj}*x_{htj}^T
 grp <- rep(1:sT,rep(J,sT))
 DT <- data.table(Xp,grp)
 PTX <- as.matrix(DT[,lapply(.SD,sum),by=grp])[,2:(K+1),drop=FALSE] # sum p_{htj}*x_{htj}^T over j
 for (h in 1:H){ 
 L <- L + (determinant(as.matrix(Sig[h,,]))$modulus[1]/2 
          - 0.5 * sum(Sig[h,,] * ( crossprod(Xp[(TJc[h]+1):TJc[h+1],],X[(TJc[h]+1):TJc[h+1],,drop=FALSE])
                                     -crossprod(PTX[(Tc[h]+1):Tc[h+1],,drop=FALSE] )))) }
 L <- L + sum(y*mp) - sum( log(mp3) + ma )
 return(L) 
}


################################
# Stochastic Linear Regression #
################################

# beta: H by K matrix
# gb: H by K matrix
# Hb: H by K by K array
# compute gradient and Hessian w.r.t. beta

gH <- function(y, X, H, K, J, Tc, TJc, T, sT, beta, muz, OU, yXH) {

Hb <- array(0, dim=c(H,K,K))
mp <- matrix(rowSums(X * beta[rep(1:H,J*T),]), sT, J, byrow=TRUE)
mp <- exp(mp - apply(mp,1,max))
mp <- as.vector(t(mp/rowSums(mp)))
Xp <- X * mp
mm <- t(beta) - muz

grp <- rep(1:sT,rep(J,sT))
DT <- data.table(Xp,grp)
temp1 <- as.matrix(DT[,lapply(.SD,sum),by=grp])[,2:(K+1),drop=FALSE] # sum over j 

grp <- rep(1:H,T)
DT <- data.table(temp1,grp)
temp2 <- as.matrix(DT[,lapply(.SD,sum),by=grp])[,2:(K+1),drop=FALSE] # sum over j and t

gb <- yXH - temp2 - crossprod(mm, OU)
for (h in 1:H) {Hb[h,,] <- -(crossprod(Xp[(TJc[h]+1):(TJc[h+1]),,drop=FALSE],X[(TJc[h]+1):(TJc[h+1]),,drop=FALSE])
                                 -crossprod(temp1[(Tc[h]+1):(Tc[h+1]),,drop=FALSE]))- OU }

list(gb=gb, Hb=Hb) }


LBVAR_SA <- function(y, X, H, K, T, mu, Sig, muz, Sigz, b, omeg, Ups, c, m0, S0i, cS0i, OU, cOU){

 mp <- rowSums(X * mu[rep(1:H,J*T),])   # x_{htj}^T %*% mu_h
 mm <- mu - matrix(rep(muz,H),H,K,byrow=TRUE)
 L <- 0.5 * ( - sum(colSums(Sig) * OU) - sum((tcrossprod(mm,cOU))^2)  
      - H * sum(Sigz*OU) - crossprod(cS0i %*% (muz-m0)) - sum(S0i * Sigz) 
      - omeg * determinant(Ups)$modulus[1] + determinant(Sigz)$modulus[1]) - sum(b * log(c)) + sum(y*mp)
  for (h in 1:H){ L <- L + determinant(as.matrix(Sig[h,,]))$modulus[1]/2 }
 return(L) 
}


  
#########################
# Variational Algorithm #
#########################

# option = "Laplace" 
# option = "NCVMP"  
# option = "SA"      

VA <- function(y, X, T, J, tol=0.005, option="NCVMP", initialfit=NULL, N=40, w=0.25, pre=1){

lb <- 0
K <- dim(X)[2]
H <- length(T)
sT <- sum(T)
Tc <- c(0,cumsum(T))
TJc <- c(0,cumsum(T*J))
Jc <- c(0,cumsum(rep(J,sT)))

yX <- as.vector(crossprod(y,X)) # K by 1, sum over h, j and t
grp <- rep(1:H,J*T)
DT <- data.table(y*X,grp)
yXH <- as.matrix(DT[,lapply(.SD,sum),by=grp])[,2:(K+1),drop=FALSE] # H by K, sum over j and t
if (option == "Laplace") { fvalue <- rep(0,H) }

# Priors #
m0 <- rep(0,K)
S0 <- 1e+06*diag(K)
S0i <- chol2inv(chol(S0))
cS0i <- chol(S0i)
Ob <- as.vector(S0i%*%m0)
nu <- 2
A <- rep(10^3,K)

# Initialization #
omeg <- nu + H + K - 1     # deterministic
b <- rep((nu+K)/2,K)       # deterministic
if ( is.null(initialfit) ){
 Ups <- (omeg-K+1)*diag(K)
 c <- b
 muz <- rep(0,K)
 Sigz <- 0.01*diag(K)
 mu <- matrix(0,H,K)
 Sig <- array(0,dim=c(H,K,K))
 for (h in 1:H) { Sig[h,,] <- 0.01*diag(K) }
 } else {
 Ups <- initialfit$Ups
 c <- initialfit$c
 muz <- initialfit$muz
 Sigz <- initialfit$Sigz
 mu <- initialfit$mu
 Sig <- initialfit$Sig 
}

OU <- omeg * chol2inv(chol(Ups))
cOU <- chol(OU)

if (option == "SA"){
 beta <- matrix(0, H, K)  # for structure
 mb <- mu
 Pb <- array(0, dim=c(H,K,K)); for (h in 1:H){ Pb[h,,] <- chol2inv(chol(Sig[h,,])) }
 gb <- matrix(0,H,K)
}

lbfix <- LBFIX(H, K, omeg, b, nu, A, S0)
count <- 0
dif <- 1
record <- NULL
totaltime <- 0
parold <- c(muz, diag(Ups), c)
if (pre > 1) { parmat <- matrix(0, pre, length(parold)); parmat[pre,] <- parold }

while (dif > tol) {

bef <- proc.time()
count <- count + 1

# Update mu, Sig #

if (option == "Laplace"){
 for (h in 1:H) {
   opt <- optim(mu[h,], fbeta, gr=grbeta, yXh=yXH[h,], Xh=X[(TJc[h]+1):TJc[h+1],,drop=FALSE], muz=muz,  
                cOU=cOU, Th=T[h], J=J, method="BFGS", control = list(fnscale=-1), hessian=TRUE)
   mu[h,]<- opt$par
   Sig[h,,] <- chol2inv(chol(-opt$hessian)) }
}


if (option == "NCVMP"){
mm <- mu - matrix(rep(muz,H), H, K, byrow=TRUE)
mp <- matrix( rowSums(X * mu[rep(1:H,J*T),]), sT, J, byrow=TRUE)
mp <- mp - apply(mp,1,max)
mp <- exp(mp) / rowSums(exp(mp))  # sT by J
mp <- as.vector(t(mp))            # sT * J by 1
Xp <- X * mp  
grp <- rep(1:sT,rep(J,sT))
DT <- data.table(Xp,grp)
PTX <- as.matrix(DT[,lapply(.SD,sum),by=grp])[,2:(K+1),drop=FALSE]  # sum p_{htj}*X_{htj}^T over j
grp <- rep(1:H,T)
DT <- data.table(PTX,grp)
sPTX <- as.matrix(DT[,lapply(.SD,sum),by=grp])[,2:(K+1),drop=FALSE] # sum p_{htj}*X_{htj}^T over j and t
A2 <- yXH - sPTX - mm %*% OU
i <- 0

for (h in 1:H){

SW <- tcrossprod(Sig[h,,],X[(TJc[h]+1):TJc[h+1],,drop=FALSE])

  for (t in 1:T[h]){
  i <- i+1
  WSW <- X[(Jc[i]+1):Jc[i+1],,drop=FALSE] %*% SW[,((t-1)*J+1):(t*J)]
  temp2 <- WSW %*% mp[(Jc[i]+1):Jc[i+1]] - 0.5 * diag(WSW)
  temp3 <- t(Xp[(Jc[i]+1):Jc[i+1],,drop=FALSE]) - as.matrix(PTX[i,]) %*% mp[(Jc[i]+1):Jc[i+1]]
  A2[h,] <- A2[h,] + temp3 %*% temp2 }

Sig[h,,] <- chol2inv(chol(OU + crossprod(Xp[(TJc[h]+1):TJc[h+1],,drop=FALSE], X[(TJc[h]+1):TJc[h+1],,drop=FALSE])
                     - crossprod(PTX[(Tc[h]+1):Tc[h+1],,drop=FALSE]) ))
mu[h,] <- mu[h,] + Sig[h,,] %*% A2[h,]  }
}


if (option == "SA"){
gbbar <- matrix(0, H, K)
mbbar <- matrix(0, H, K)
Pbbar <- array(0, dim=c(H, K, K))
 
 for (num in 1:N){

# Generate random variates:
  Z <- matrix(rnorm(H*K), H, K)
  for (h in 1:H) { 
  U <- backsolve( chol(Pb[h,,]),diag(K) ) 
  beta[h,] <- U %*% (Z[h,] + crossprod(U, gb[h,])) 
  }
  beta <- beta + mb

# Compute gradients and Hessians
  out <- gH(y, X, H, K, J, Tc, TJc, T, sT, beta, muz, OU, yXH)
  mb <- (1-w) * mb + w * beta
  gb <- (1-w) * gb + w * out$gb
  Pb <- (1-w) * Pb - w * out$Hb

  if ( num > (N/2) ){
  mbbar <- mbbar + beta/N*2
  gbbar <- gbbar + out$gb/N*2
  Pbbar <- Pbbar - out$Hb/N*2  
  } 
 }

for (h in 1:H) {
Sig[h,,] <- chol2inv(chol(Pbbar[h,,]))
mu[h,] <- Sig[h,,]%*%gbbar[h,] + mbbar[h,] 
}

mb <- mu
Pb <- Pbbar
gb <- matrix(0,H,K)
}

Sigz <- chol2inv(chol(H * OU + S0i))
muz <- as.vector( Sigz %*% (OU %*% colSums(mu) + Ob) )

mm <- t(mu)- muz                               # K by H matrix
if (K==1) {ct <- as.numeric(2 * nu * b/c + H * Sigz)} else {ct <- 2 * nu * diag(b/c) + H * Sigz}
Ups <- ct + colSums(Sig) + tcrossprod(mm)
OU <- omeg * chol2inv(chol(Ups))       
cOU <- chol(OU)               # Cholesky decomposition of OU
c <- nu * diag(OU) + 1/(A^2)  


par <- c(muz, diag(Ups), c)
if (pre == 1) {dif <- max(abs(par-parold)/(abs(parold) + 1.0e-8))
} else {
 parmat <- rbind(parmat[2:pre,],par)
 parmatnew <- colSums(parmat)/pre
 if (count <= pre) { dif <- max(abs(par-parold)/(abs(parold) + 1.0e-8))
 } else { dif <- max(abs(parmatnew - parmatold)/(abs(parmatold) + 1.0e-8)) }
 parmatold <- parmatnew 
}
parold <- par

aft <- proc.time()
time <- (aft-bef)[1]
totaltime <- totaltime + time

if (option == "NCVMP"){ lb <- lbfix + LBVAR_NCVMP(y, X, H, J, K, T, sT, Tc, TJc, mu, Sig, muz, Sigz, b, c, omeg, Ups, m0, S0i, cS0i, OU, cOU) }
if (option == "Laplace"){ 
 for (h in 1:H) { fvalue[h] <- fbeta(mu[h,], yXH[h,], Xh=X[(TJc[h]+1):TJc[h+1],,drop=FALSE], muz, cOU, Th=T[h], J) }
  lb <- lbfix + LBVAR_Laplace(H, K, Sig, muz, Sigz, b, c, omeg, Ups, m0, S0i, cS0i, OU, fvalue) 
}

record <- rbind(record, c(count, dif, totaltime, lb))
cat(count, dif, lb, totaltime, round(muz,2), round(diag(Ups)), "\n")
}

list(count = count, mu=mu, Sig=Sig, muz = muz, Sigz = Sigz, omeg = omeg, Ups = Ups, 
b = b, c = c, w = w, totaltime = totaltime, lb=lb, record = record)
}





VALARGE <- function(y, X, T, J, tol=0.005, option="NCVMP", initialfit=NULL, N=40, w=0.25, IB=25, fac=20, pre=1){

num=0
lb <- 0
K <- dim(X)[2]
H <- length(T)
sT <- sum(T)
Tc <- c(0,cumsum(T))
TJc <- c(0,cumsum(T*J))
Jc <- c(0,cumsum(rep(J,sT)))

yX <- as.vector(crossprod(y,X)) # K by 1, sum over h, j and t
grp <- rep(1:H,J*T)
DT <- data.table(y*X,grp)
yXH <- as.matrix(DT[,lapply(.SD,sum),by=grp])[,2:(K+1),drop=FALSE] # H by K, sum over j and t
if (option == "Laplace") { fvalue <- rep(0,H) }

STEP <- function(B, H) ( 1 - (1 - 0.4) * (H - B) / (H-25) )
CUT <-  function(B, H) ( 1 - (1 - 0.4) * (H - B) / (H-25) )


# Priors #
m0 <- rep(0,K)
S0 <- 10^6*diag(K)
S0i <- chol2inv(chol(S0))
cS0i <- chol(S0i)
Ob <- as.vector(S0i%*%m0)
nu <- 2
A <- rep(10^3,K)

# Initialization #
omeg <- nu + H + K - 1     # deterministic
b <- rep((nu+K)/2,K)       # deterministic

if ( is.null(initialfit) ){
  Ups <- (omeg-K+1)*diag(K)
  c <- b
  muz <- rep(0,K)
  Sigz <- 0.01*diag(K)
  mu <- matrix(0,H,K)
  Sig <- array(0,dim=c(H,K,K))
  InvSig <- array(0,dim=c(H,K,K))
  for (h in 1:H) { Sig[h,,] <- 0.01*diag(K); InvSig[h,,] <- 100*diag(K) }
} else {
  Ups <- initialfit$Ups
  c <- initialfit$c
  muz <- initialfit$muz
  Sigz <- initialfit$Sigz
  mu <- initialfit$mu
  Sig <- initialfit$Sig
  InvSig <- array(0,dim=c(H,K,K))
  for (h in 1:H) { InvSig[h,,] <- chol2inv(chol(Sig[h,,])) }
}

OU <- omeg * chol2inv(chol(Ups))
cOU <- chol(OU)
lbfix <- LBFIX(H, K, omeg, b, nu, A, S0)

B <- IB
s <- STEP(B, H)
co <- CUT(B, H)
int <- 20
hyp <- FALSE
RPP <- 1
mat <- c(muz, diag(Ups))
record <- NULL
icount <- 0
totaltime <- 0
if (option == "SA") {beta <- matrix(0, B, K)}


while (B <= H | hyp == FALSE) {


if (B == H & hyp == TRUE) break

if (hyp == TRUE) {
B <- min(fac*B,H)
s <- STEP(B, H)
co <- CUT(B,H)
RPP <- 1
hyp <- FALSE
icount <- 0    
mat <- c(muz, diag(Ups))
if (option == "SA") {beta <- matrix(0, B, K)}
}     

bef <- proc.time()
icount <- icount + 1

if (B < H) {

batch <- sort(sample(1:H, B))
rowidX <- which(rep(1:H, T*J) %in% batch)
Xb <- X[rowidX,,drop=FALSE]
yb <- y[rowidX]
Tb <- T[batch]
sTb <- sum(Tb)
Jcb <- c(0, cumsum(rep(J,sTb)))
TJcb <- c(0, cumsum(Tb*J))
Tcb <- c(0, cumsum(Tb))
grp <- rep(1:B, J*Tb)
DT <- data.table(yb * Xb,grp)
yXHb <- as.matrix(DT[,lapply(.SD,sum),by=grp])[,2:(K+1),drop=FALSE] # B by K, sum over j and t


# Update mu, Sig #

if (option == "Laplace"){
 for (h in batch) {
   opt <- optim(mu[h,], fbeta, gr=grbeta, yXh=yXH[h,], Xh=X[(TJc[h]+1):TJc[h+1],,drop=FALSE], muz=muz, 
                cOU=cOU, Th=T[h], J=J, method="BFGS", control = list(fnscale=-1), hessian=TRUE)
   mu[h,]<- opt$par
   Sig[h,,] <- chol2inv(chol(-opt$hessian)) }
}


if (option == "NCVMP"){

num <- 0
localdiff <- 1
muold <- mu
while (localdiff > 0.1 & num < 3) {
num <- num+1

mm <- mu[batch,,drop=FALSE] - matrix(rep(muz,B), B, K, byrow=TRUE)
mp <- mu[batch,,drop=FALSE]
mp <- matrix(rowSums( Xb * mp[rep(1:B,J*Tb),]), sTb, J, byrow=TRUE)
mp <- mp - apply(mp,1,max)
mp <- exp(mp)/rowSums(exp(mp))
mp <- as.vector(t(mp))
Xp <- Xb * mp

grp <- rep(1:sTb,rep(J,sTb))
DT <- data.table(Xp,grp)
PTX <- as.matrix(DT[,lapply(.SD,sum),by=grp])[,2:(K+1),drop=FALSE]  # sum p_{htj}*X_{htj}^T over j

grp <- rep(1:B, Tb)
DT <- data.table(PTX, grp)
sPTX <- as.matrix(DT[,lapply(.SD,sum),by=grp])[,2:(K+1),drop=FALSE]

A2 <- yXHb - sPTX - mm %*% OU  
n1 <- 0
i <- 0

for (h in batch){
i <- i+1

SW <- tcrossprod(Sig[h,,],X[(TJc[h]+1):TJc[h+1],,drop=FALSE])

  for (t in 1:T[h]){
  n1 <- n1+1
  WSW <- Xb[(Jcb[n1]+1):Jcb[n1+1],,drop=FALSE] %*% SW[,((t-1)*J+1):(t*J)]
  temp2 <- WSW %*% mp[(Jcb[n1]+1):Jcb[n1+1]] - 0.5 * diag(WSW)
  temp3 <- t(Xp[(Jcb[n1]+1):Jcb[n1+1],,drop=FALSE]) - as.matrix(PTX[n1,]) %*% mp[(Jcb[n1]+1):Jcb[n1+1]]
  A2[i,] <- A2[i,] + temp3 %*% temp2
  }

Sig[h,,] <- chol2inv(chol(OU + crossprod(Xp[(TJcb[i]+1):TJcb[i+1],,drop=FALSE], X[(TJc[h]+1):TJc[h+1],,drop=FALSE])
                     - crossprod(PTX[(Tcb[i]+1):Tcb[i+1],,drop=FALSE]) ))
mu[h,] <- mu[h,] + Sig[h,,] %*% A2[i,]  }

localdiff <- norm(as.matrix(mu[batch,]-muold[batch,]),type='F')/ norm(as.matrix(mu[batch,]),type='F') # Euclidean norm
muold <- mu}
}


if (option == "SA"){

gbbar <- matrix(0,B,K)
mbbar <- matrix(0,B,K)
Pbbar <- array(0,dim=c(B,K,K))
mb <- mu[batch,,drop=FALSE]
Pb <- InvSig[batch,,,drop=FALSE]
gb <- matrix(0, B, K)

for (num in 1:N){

# Generate random variates:
 Z <- matrix(rnorm(B*K), B, K)
 i <- 0
 for (h in batch) { 
  i <- i+1
  U <- backsolve( chol(Pb[i,,]), diag(K) ) 
  beta[i,] <- U %*% (Z[i,] + crossprod(U, gb[i,])) 
 }
 beta <- beta + mb

# Compute gradients and Hessians
 out <- gH(yb, Xb, B, K, J, Tcb, TJcb, Tb, sTb, beta, muz, OU, yXHb)
 mb <- (1-w) * mb + w * beta
 gb <- (1-w) * gb + w * out$gb
 Pb <- (1-w) * Pb - w * out$Hb

 if ( num > (N/2) ) {
   mbbar <- mbbar + beta/N*2
   gbbar <- gbbar + out$gb/N*2
   Pbbar <- Pbbar - out$Hb/N*2 }
}

 i <- 0
 for (h in batch) {
   i <- i + 1
   InvSig[h,,] <- Pbbar[i,,]
   Sig[h,,] <- chol2inv(chol(Pbbar[i,,]))
   mu[h,] <- Sig[h,,] %*% gbbar[i,] + mbbar[i,] }
}

mubatch <- mu[batch,,drop=FALSE] # B by K matrix
Sigbatch <- Sig[batch,,,drop=FALSE]

Sigz <- chol2inv(chol(H * OU + S0i))
mng <- Sigz %*% (H/B * OU %*% colSums(mubatch) + Ob)
muz <- as.vector((1 - s) * muz + s * mng)

mm <- t(mubatch)- muz        # K by B matrix
if (K==1) { ct <- as.numeric(2 * nu * b/c + H * Sigz) } else { ct <- 2 * nu * diag(b/c) + H * Sigz }
Ups <- (1 - s) * Ups + s * ( H / B * (colSums(Sigbatch) + tcrossprod(mm)) + ct ) 
OU <- omeg * chol2inv(chol(Ups))      
cOU <- chol(OU)               # Cholesky decomposition of OU
c <- nu * diag(OU) + 1/(A^2) 


if (icount <= int) {
   mat <- rbind(mat,  c(muz, diag(Ups)) )
   if (icount > 5) {
   RPP <- abs(mat[1,] - mat[(icount+1),]) / colSums(abs(mat[1 : icount,,drop=FALSE] - mat[2 : (icount+1),,drop=FALSE]))}
} else {
   mat <- rbind(mat[2:int,,drop=FALSE],  c(muz, diag(Ups)) )
   RPP <- abs(mat[1,] - mat[int,]) / colSums(abs(mat[1 : (int - 1),,drop=FALSE] - mat[2 : int,,drop=FALSE]))
}

if (max(RPP) < co) {hyp <- TRUE} else {hyp <- FALSE}


} else {

if (icount==1) { 
 parold <- c(muz, diag(Ups), c)
 if (pre > 1) { parmat <- matrix(0, pre, length(parold)); parmat[pre,] <- parold } 
}


# Update mu, Sig #

if (option == "Laplace"){
 for (h in 1:H) {
   opt <- optim(mu[h,], fbeta, gr=grbeta, yXh=yXH[h,], Xh=X[(TJc[h]+1):TJc[h+1],,drop=FALSE], muz=muz, 
                cOU=cOU, Th=T[h], J=J, method="BFGS", control = list(fnscale=-1), hessian=TRUE)
   mu[h,]<- opt$par
   Sig[h,,] <- chol2inv(chol(-opt$hessian)) }
}


if (option == "NCVMP"){

 if (icount==1){
   num <- 0
   localdiff <- 1
   muold <- mu
   while (localdiff > 0.1 & num < 3) {
   num <- num+1

   mm <- mu - matrix(rep(muz,H), H, K, byrow=TRUE)
   mp <- matrix(rowSums( X * mu[rep(1:H,J*T),]), sT, J, byrow=TRUE)
   mp <- mp - apply(mp,1,max)
   mp <- exp(mp)/rowSums(exp(mp))
   mp <- as.vector(t(mp))
   Xp <- X * mp

   grp <- rep(1:sT,rep(J,sT))
   DT <- data.table(Xp,grp)
   PTX <- as.matrix(DT[,lapply(.SD,sum),by=grp])[,2:(K+1),drop=FALSE]  

   grp <- rep(1:H, T)
   DT <- data.table(PTX, grp)
   sPTX <- as.matrix(DT[,lapply(.SD,sum),by=grp])[,2:(K+1),drop=FALSE]

   A2 <- yXH - sPTX - mm %*% OU  
   n1 <- 0

   for (h in 1:H){

   SW <- tcrossprod(Sig[h,,],X[(TJc[h]+1):TJc[h+1],,drop=FALSE])

   for (t in 1:T[h]){
    n1 <- n1+1
    WSW <- X[(Jc[n1]+1):Jc[n1+1],,drop=FALSE] %*% SW[,((t-1)*J+1):(t*J)]
    temp2 <- WSW %*% mp[(Jc[n1]+1):Jc[n1+1]] - 0.5 * diag(WSW)
    temp3 <- t(Xp[(Jc[n1]+1):Jc[n1+1],,drop=FALSE]) - as.matrix(PTX[n1,]) %*% mp[(Jc[n1]+1):Jc[n1+1]]
    A2[h,] <- A2[h,] + temp3 %*% temp2
   }

  Sig[h,,] <- chol2inv(chol(OU + crossprod(Xp[(TJc[h]+1):TJc[h+1],,drop=FALSE], X[(TJc[h]+1):TJc[h+1],,drop=FALSE])
                     - crossprod(PTX[(Tc[h]+1):Tc[h+1],,drop=FALSE]) ))
  mu[h,] <- mu[h,] + Sig[h,,] %*% A2[h,]  }

  localdiff <- norm(as.matrix(mu-muold),type='F')/ norm(as.matrix(mu),type='F') 
  muold <- mu}

 } else {

  mm <- mu - matrix(rep(muz,H), H, K, byrow=TRUE)
  mp <- matrix( rowSums(X * mu[rep(1:H,J*T),]), sT, J, byrow=TRUE)
  mp <- mp - apply(mp,1,max)
  mp <- exp(mp) / rowSums(exp(mp))  # sT by J
  mp <- as.vector(t(mp))            # sT * J by 1
  Xp <- X * mp  
  grp <- rep(1:sT,rep(J,sT))
  DT <- data.table(Xp,grp)
  PTX <- as.matrix(DT[,lapply(.SD,sum),by=grp])[,2:(K+1),drop=FALSE]  # sum p_{htj}*X_{htj}^T over j

  grp <- rep(1:H,T)
  DT <- data.table(PTX,grp)
  sPTX <- as.matrix(DT[,lapply(.SD,sum),by=grp])[,2:(K+1),drop=FALSE] # sum p_{htj}*X_{htj}^T over j and t

  A2 <- yXH - sPTX - mm %*% OU
  n1 <- 0

  for (h in 1:H){

  SW <- tcrossprod(Sig[h,,],X[(TJc[h]+1):TJc[h+1],,drop=FALSE])

  for (t in 1:T[h]){
    n1 <- n1+1
    WSW <- X[(Jc[n1]+1):Jc[n1+1],,drop=FALSE] %*% SW[,((t-1)*J+1):(t*J)]
    temp2 <- WSW %*% mp[(Jc[n1]+1):Jc[n1+1]] - 0.5 * diag(WSW)
    temp3 <- t(Xp[(Jc[n1]+1):Jc[n1+1],,drop=FALSE]) - as.matrix(PTX[n1,]) %*% mp[(Jc[n1]+1):Jc[n1+1]]
    A2[h,] <- A2[h,] + temp3 %*% temp2 }

  Sig[h,,] <- chol2inv(chol(OU + crossprod(Xp[(TJc[h]+1):TJc[h+1],,drop=FALSE], X[(TJc[h]+1):TJc[h+1],,drop=FALSE])
                     - crossprod(PTX[(Tc[h]+1):Tc[h+1],,drop=FALSE]) ))
  mu[h,] <- mu[h,] + Sig[h,,] %*% A2[h,]  }
}
}


if (option == "SA"){

mb <- mu
Pb <- InvSig
gb <- matrix(0, H, K)
gbbar <- matrix(0,H,K)
mbbar <- matrix(0,H,K)
Pbbar <- array(0,dim=c(H,K,K))
 
for (num in 1:N){

# Generate random variates:
  Z <- matrix(rnorm(H*K), H, K)
  for (h in 1:H) { 
  U <- backsolve( chol(Pb[h,,]),diag(K) ) 
  beta[h,] <- U %*% (Z[h,] + crossprod(U, gb[h,])) }
  beta <- beta + mb

# Compute gradients and Hessians
  out <- gH(y, X, H, K, J, Tc, TJc, T, sT, beta, muz, OU, yXH)
  mb <- (1-w) * mb + w * beta
  gb <- (1-w) * gb + w * out$gb
  Pb <- (1-w) * Pb - w * out$Hb

  if ( num > (N/2) ){
  mbbar <- mbbar + beta/N*2
  gbbar <- gbbar + out$gb/N*2
  Pbbar <- Pbbar - out$Hb/N*2  } 
 }

 for (h in 1:H) {
  InvSig[h,,] <- Pbbar[h,,]
  Sig[h,,] <- chol2inv(chol(Pbbar[h,,]))
  mu[h,] <- Sig[h,,]%*%gbbar[h,] + mbbar[h,]
 }
}

Sigz <- chol2inv(chol(H * OU + S0i))
muz <- as.vector( Sigz %*% (OU %*% colSums(mu) + Ob) )

mm <- t(mu)- muz                               # K by H matrix
if (K==1) {ct <- as.numeric(2 * nu * b/c + H * Sigz)} else {ct <- 2 * nu * diag(b/c) + H * Sigz}
Ups <- ct + colSums(Sig) + tcrossprod(mm)
OU <- omeg * chol2inv(chol(Ups))       
cOU <- chol(OU)             
c <- nu * diag(OU) + 1/(A^2)  

par <- c(muz, diag(Ups), c)
if (pre == 1) {dif <- max(abs(par-parold)/(abs(parold) + 1.0e-8))
} else {
 parmat <- rbind(parmat[2:pre,],par)
 parmatnew <- colSums(parmat)/pre
 if (icount <= pre) { dif <- max(abs(par-parold)/(abs(parold) + 1.0e-8))
 } else { dif <- max(abs(parmatnew - parmatold)/(abs(parmatold) + 1.0e-8)) }
 parmatold <- parmatnew 
}
parold <- par

if ( dif < tol ) { hyp <- TRUE } else { hyp <- FALSE }
}

aft <- proc.time()
time <- (aft-bef)[1]
totaltime <- totaltime + time

if (B < H)  {
cat( icount, B, totaltime, max(RPP), round(muz,2), round(diag(Ups)), co, s, num,"\n")
record <- rbind( record, c(B, icount, totaltime, max(RPP)) )}
if (B == H) {
cat( icount, B, totaltime, dif, round(muz,2), round(diag(Ups)), num, "\n")
record <- rbind( record, c(B, icount, totaltime, dif) )}
}

if (option == "Laplace"){ 
for (h in 1:H) { fvalue[h] <- fbeta(mu[h,], yXH[h,], X[(TJc[h]+1):TJc[h+1],,drop=FALSE], muz, cOU, Th=T[h], J) }
lb <- lbfix + LBVAR_Laplace(H, K, Sig, muz, Sigz, b, c, omeg, Ups, m0, S0i, cS0i, OU, fvalue) 
}

if (option == "NCVMP") { 
lb <- lbfix + LBVAR_NCVMP(y, X, H, J, K, T, sT, Tc, TJc, mu, Sig, muz, Sigz, b, c, omeg, Ups, m0, S0i, cS0i, OU, cOU) }

cat(lb, "\n")
list(mu = mu, Sig = Sig, muz = muz, Sigz = Sigz, omeg = omeg, Ups = Ups, b = b, c = c, totaltime = totaltime, lb=lb, record = record)
}



