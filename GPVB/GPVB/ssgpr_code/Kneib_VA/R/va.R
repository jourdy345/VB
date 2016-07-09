Rdelta <- function(vektor) .C("R_delta",x=as.double(vektor), w=numeric(1), w1=numeric(1),PACKAGE="VA")
### Einführen der Checkfunktion zur Bewertung der Quantile

checkf <- function(u, tau){
check<- u * (tau - (u<0))
return(check)
}

deltadens <- function(x,p,q,r) p*log(x)+q*(x)-r*x^2

###Einführen des MSE zur Bewertung der Mittelwerte
mse <- function(y, ydach){
mean((y-ydach)^2)
}

fu <- function(x,si){
  n<- nrow(x)
  v <- vector(length=n)
  for(i in 1:n){
  v[i] <- x[i,]%*%si[,i]
 }
  return(v)
}
# Rdelta
va <- function(formula,  data=NULL, intercept=TRUE,
              quant = FALSE, ta = 0.5,  k = 100, sto = 1e-04,
              a_0 = 1, b_0 = 1e-05,  
              plotmse = FALSE, norm = TRUE, pred = NULL,...)

{


mu_delt = 1
me <-0
geo <-NULL
gramap <- NULL
gam <- NULL
be <- NULL

deltw <- mu_delt
mse <-rep(0, k)
delli <- rep(0,k)
##Funktion zur schnelleren Matrixhälftenmultiplikation

col <- c("black")

data <- as.data.frame(data)
y <- (model.frame(formula, data)[,1])
n <- length(y)
#### Normierung der Zielgröße
  if(norm){
  yw <- y
  me <- (max(y)+min(y))/2
  se <- max(abs(y-me))
  y <- (y-me)/(max(abs(y-me)))}
print('>')
n <- length(y)
qua <-0
j<-1
a_1_s <- 0
b_1 <- 0
MU0<-0
SIG0<-0
  Z <- list()
  Zt<- list()
  ZtZ <- list()
knots<-list()
mx<-list()
K <- list()
Xlinear <- as.matrix(rep(1,n))
Xnonlinear <- as.matrix(rep(1,n))
colnames(Xlinear)[1] <- "intercept"
type <- vector(length = length(model.frame(formula, data))-1)
for(i in 2:length(model.frame(formula, data=data))){
  b <- model.frame(formula, data)[,i]
   if(is.null(attr(b,"class"))){Xlinear<-cbind(Xlinear,b)
                                colnames(Xlinear)[ncol(Xlinear)]<-paste(get.vars(formula)[i])
                                MU0<-c(MU0,0)
                                SIG0<-c(SIG0,10000)
                                }  
  else{if(attr(b,"class")%in%c("bspline", "random", "geo")){
    Z[[j]]<-as.matrix(b)
    names(Z)[j] <- paste(get.vars(formula)[i])
            K[[j]] <- attr(b, "K")
    
  Zt[[j]]<-t(Z[[j]])                                                    
  ZtZ[[j]] <- t(Z[[j]])%*%Z[[j]]
  rang  <- qr(K[[j]])$rank
  a_1 <- attr(b, "a_1")
  a_1_s<-c(a_1_s, a_1+ rang/2)
  b_1<-c(b_1,attr(b, "b_1"))
    if(attr(b,"class")==c("bspline")){Xnonlinear <- cbind(Xnonlinear, attr(b, "data"))
                                      knots[[j]]<-attr(b,"knots")
                                      mx[[j]]<- attr(b, "mx")}
    if(attr(b,"class")==c("geo")){geo <-attr(b, "data")
                                  names(Z)[j]<- "geo"}
    if(attr(b,"class")==c("random")){names(Z)[j]<-attr(b, "name")}
    j<-j+1}
       if(attr(b,"class")==("linear")){
        Xlinear<-cbind(Xlinear,b)
        colnames(Xlinear)[ncol(Xlinear)]<-names(Z)[j]
        MU0<-c(MU0, attr(b, "mu0"))
        SIG0<-c(SIG0, attr(b, "sig0"))
        
       }
  }
}
print('>>')
Xnonlinear<- Xnonlinear[,-1]
if(length(Z)>0){spli<-TRUE}else{spli<-FALSE}   
if(spli){
   Zg <- Z[[1]]
q <- length(Z)
  z <-matrix(nrow=q, ncol=2)                                                                        
  z[1,1]<-1
  z[1,2] <- ncol(Zg)

   if(q >1){
     for(j in 2:q){
   ### Wert zuweisen, um geschickter auf die Spalten zugreifen zu können
    z[j,1]<-ncol(Zg)+1
    Zg <- cbind(Zg, Z[[j]])
   z[j,2]<-ncol(Zg)
   }
   }
rownames(z)<- names(Z)
   print(z)        
a_1_s <- a_1_s[-1]
b_1 <- b_1[-1]
 Zgt <-(t(Zg))
   ZgtZg <- Zgt%*%Zg
  l <- ncol(Zg)
  mu_gamG <- rep(1, l)
  mu_gamGneu <- rep(1, l)
  sig_gamG <- matrix(nrow=l, ncol=l, data=0)
    sig_gam_invG <-diag(1,l)
  mu_tau_inv <- rep(50, q)
print('>>>')


sigZ<-sig_gamG%*%Zgt                  
   }else{

Zg<- as.matrix(rep(0,n))
mu_gamG<-0
Zgt<-t(Zg)
mu_tau_inv<-0
sig_gamG<-0
  ZgtZg <- 0
sigZ<-0}
if(intercept==FALSE){
  if(ncol(Xlinear)>1){
  Xlinear <- Xlinear[,-1]
  li <- TRUE}else{li<-FALSE}
  }
if(intercept){
  li=TRUE}

### Lineare Effekte
if(intercept){
  SIG0[1] <- 10000
}else{
MU0 <- MU0[-1]
SIG0 <- SIG0[-1]}

if(li){ 
  Xlinear <- as.matrix(Xlinear)  
  mu_bet_prior=MU0
  if(ncol(Xlinear)>1){
  sig_bet_prior =diag(SIG0)}else{sig_bet_prior<-SIG0}
    Xtlin <- t(Xlinear)
  pred_bet_prior <- solve(sig_bet_prior)
    XlintXlin <- Xtlin%*%Xlinear

  pred_bet_prior <- solve(sig_bet_prior)

  mm <- pred_bet_prior%*%mu_bet_prior
  mu_bet <- rep(1,ncol(Xlinear))           
  }else{
     int <- rep(0,n)
  Xlinear <- as.matrix(int)
  XlintXlin <- 0
  Xtlin <- t(Xlinear)
  pred_bet_prior <- 0
  sig_bet_prior<-0
  mm <- 0
  mu_bet<-0
  }
 
print('>>>>')
## Quantilregression  
  if(quant){
  w <- rep(1,n)
  xi <- (1-2*ta)/(ta*(1-ta))
  sig2 <-  2/(ta*(1-ta))
  p <- 2*a_0 +n -2 + nrow(pred_bet_prior)
  print(nrow(pred_bet_prior))
  print(p)
}

res <- y

  quak <- 0.001

  a_0_s <- a_0 + n/2
    sig_bet <- sig_bet_prior

############################################################################################
###  ab hier beginnt das echte Programm   ##################################################
############################################################################################
  for(i in 1:k)
    {
       
    ###########################################
    ## Gewichte für quant #####################
    ###########################################   
    if(quant){        

    fa<-4*ta^2*(1-ta)^2*mu_delt
    resq<- res^2
      if(spli){
        
      zum <- fu(Zg, sigZ)}else{zum <- rep(0,n)}
          if(li){
        xum <- fu(Xlinear, sig_bet%*%Xtlin)}else{xum<-rep(0,n)}
for(j in 1:n){w[j] <- (fa*(resq[j] + zum[j]+ xum[j]))^{-1/2}}
    winv <- 1/w+4*ta*(1-ta)#((1-2*ta)^2+2)
wend<-w
    ###########################################
    ## beta ###################################
    ###########################################   
print('>>>>>')
 if(li){
      sig_bet_inv <- (pred_bet_prior+ mu_delt*t((w)*Xlinear)%*%(Xlinear)/(sig2))
      sig_bet <- solve(sig_bet_inv)
      mu_bet <- (1/sig2)*sig_bet%*%(Xtlin%*%(w*(y*mu_delt-Zg%*%mu_gamG*mu_delt-.5*xi*deltw/w))+mm)
      }
      }

 
  
  ### wieder zurück zu mean
else{
if(li){
      sig_bet_inv <- (pred_bet_prior+ XlintXlin)*mu_delt
      sig_bet     <- solve(sig_bet_inv)
      mu_bet      <- sig_bet%*%(Xtlin%*%(y -Zg%*%mu_gamG)+mm)*mu_delt
      }  
    }

    ###########################################
    ## gamma ##################################
    ########################################### 
    if(spli){
  
        for(j in 1:q){


    ind <- z[j,1]:z[j,2]
    if(quant){
    if(names(Z)[j]=="random"){
      sig_gam_invG[ind, ind]<- diag(mu_tau_inv[j], length(ind)) + mu_delt*(t(w*Z[[j]]))%*%Z[[j]]/sig2
    }else{  
     sig_gam_invG[ind, ind]<- mu_tau_inv[j]*K[[j]] + mu_delt*(t(w*Z[[j]]))%*%Z[[j]]/sig2}
    #if(det(sig_gam_invG[ind, ind])!=0){
    sig_gamG[ind, ind]    <- solve(sig_gam_invG[ind,ind])#}else{print(i)}
    sigZ[ind,]<-sig_gamG[ind, ind]%*%(Zt[[j]])
    mu_gamGneu[ind]       <- mu_delt/sig2*sigZ[ind,]%*%(w*y - w*Zg[,-ind]%*%mu_gamGneu[-ind] - 
                                                          w*Xlinear%*%mu_bet-deltw*.5*xi/(mu_delt))
    mu_gamG[ind]          <- mu_gamGneu[ind]- mean(mu_gamGneu[ind])  
    
    }
    
    else{
    if(names(Z)[j]=="random"){
      sig_gam_invG[ind, ind] <- diag(mu_tau_inv[j], length(ind)) + mu_delt*ZtZ[[j]]
    }else{
      sig_gam_invG[ind, ind] <- mu_tau_inv[j]*K[[j]] + mu_delt*ZtZ[[j]]}
    sig_gamG[ind, ind]     <- solve(sig_gam_invG[ind,ind])
    sigZ[ind,]<-sig_gamG[ind, ind]%*%(Zt[[j]])
    mu_gamGneu[ind]        <- mu_delt*sigZ[ind,] %*%(y - (Zg[,-ind]%*%mu_gamG[-ind] + Xlinear%*%mu_bet))
    mu_gamG[ind]           <- mu_gamGneu[ind]- mean(mu_gamGneu[ind])
    }

print('>>>>>>')
###########################################
## tau^2 ##################################
########################################### 
if(names(Z)[j]=="random"){
  b_1_s <- b_1[j] + 0.5*sum((mu_gamG[ind])^2)+ 0.5*sum((sig_gamG[ind, ind]))
}else{
   b_1_s <- b_1[j] + 0.5*t(mu_gamG[ind])%*%K[[j]]%*%mu_gamG[ind] + 0.5*sum((K[[j]]*sig_gamG[ind, ind]))}

    mu_tau_inv[j] <- (a_1_s[j])/b_1_s   
    
    }  

}

###########################################
## delta ##################################
########################################### 
res <- (y - (Zg%*%mu_gamG+Xlinear%*%mu_bet))
if(quant) {# C3 <- (ta*(1-ta))*(sum(w*(res^2)) + 
  #(sum(diag(t(Xlinear)%*%diag(w)%*%Xlinear)%*%sig_bet))
  #sum(w*(Xlinear%*%sig_bet*Xlinear))+sum(w*(Zgt*sigZ))
#)+b_0 + t(mu_bet)%*%pred_bet_prior%*%(mu_bet)
       #     print(.5*t(mu_bet)%*%pred_bet_prior%*%(mu_bet))
 #C3<-0           
 C3 <- 0.5*(ta*(1-ta))*(sum(w*(res^2)) + 
  sum(w*(Xlinear%*%sig_bet*Xlinear) )
+sum(w*(Zgt*sigZ))
)+b_0

   C2 <- (.5-ta)*sum(res)


print('>>>>>>>')
vektor <- c(p+2,C2,C3)

   
print('>>>>>>>>')
integr <- Rdelta(vektor) 
if(is.na(integr$w)){
  col <- cbind(col, "blue")
lpdf <- function(x) deltadens(x,p,C2,C3)
gen <- ars.new(logpdf=lpdf, lb=0, ub=Inf)
deb <- ur(gen,100)
integr$w <- mean(deb^2)
integr$w1<-mean((deb))
}else{col<-cbind(col, "black")}
print('>>>>>>>>>')

mu_delt <-integr$w
deltw <- integr$w1
delli[i]<-mu_delt
mse[i] <- mean(checkf(res,ta)) 
   if(i>3){

       if(abs(mse[i-1]-mse[i])/mse[i-1] <sto)# & mse[i]<mse[i-1])
	{
			break
			}
			
      }
   if(plotmse){
     par(mfrow=c(1,2))
t <- 1:i
  ind <- which(mse!=0)            
plot(t[ind], mse[ind], col=col,
  ylab="risk", xlab="iterations",
  las=1, pch=20)

     plot(t[ind], delli[ind], col=col,
  ylab="delta", xlab="iterations",
  las=1, pch=20)
 abline(h=mse[i])
}

}
   else{
		b_0_s <- 0.5*sum(res^2) + 0.5*(sum((ZgtZg*sig_gamG))+ sum((XlintXlin*sig_bet))) + b_0
		mu_delt <- as.numeric((a_0_s)/b_0_s)
		mse[i]<-   (sum(res^2)/n)
		   if(i>3){
       if(abs(mse[i-1]-mse[i])/mse[i-1] <sto)# & mse[i]<mse[i-1])
	{
			break
			}
			}
    
    if(plotmse){

t <- 1:i
  ind <- which(mse!=0)
plot(t[ind], mse[ind], col=c(rep("black",length(ind)-1),"red"), las=1, pch=20, xlab="iterations", ylab="mse")
 abline(h=mse[i])
 } 

    }
    }
qua<- qua[-1]
  if(length(qua)==(k)){
  print("There were integrating problems for one parameter, maybe the data is too sparse in this quantile.")
  }
  if(i==k){
  print("The function did not reach the stopping criterion.")
  }

###########################################
## Normierung   ###########################
###########################################	
  
  if(norm){
  mu_bet<-mu_bet*se
  mu_gamG <- mu_gamG*se
  mu_delt <- mu_delt/(se^2)
  sig_bet <- sig_bet*se^2
  sig_gamG <- sig_gamG*se^2
  if(intercept){  mu_bet[1] <- mu_bet[1]+me}
  residuals <- y -(Xlinear%*%mu_bet + Zg%*%mu_gamG)
  }else{residuals <- res}
  
###########################################
## Prediction  ############################
###########################################	

  if(is.null(pred)==FALSE){
datapred <- as.data.frame(pred)  
Zpred <- list()
  j <- 1
n <- nrow(pred)
Xlinearpred <- as.matrix(rep(1,n))
if(length(which(colnames(datapred)==as.character((formula)[2])))==0){
  datapred<- cbind(1, datapred)
  colnames(datapred)[1]<-as.character((formula)[2])
}
if("random"%in%names(Z)){
  datapred<-cbind(datapred, 0)
  colnames(datapred)[ncol(datapred)]<-"ind"
}
type <- vector(length = length(model.frame(formula, datapred))-1)
for(l in 2:length(model.frame(formula, data=datapred))){
  b <- model.frame(formula, datapred)[,l]
   if(is.null(attr(b,"class"))){Xlinearpred<-cbind(Xlinearpred,b)}  
  else{if(attr(b,"class")%in%c("bspline", "geo", "random")){
    if(attr(b,"class")=="bspline"){
   pas<-which(names(Z)== names(get_all_vars(formula)[l]))
    KNOTS<-knots[[pas]]
    MX <- mx[[pas]]
    b <- bsplines(attr(b, "data"), knots=KNOTS, mx=MX)}
      Zpred[[j]]<-as.matrix(b)
  
  j<-j+1}
  }
}

if(spli){
   Zpredg <- Zpred[[1]]
   if(q >1){
 	for(j in c(2:q)){
      	     Zpredg <- cbind(Zpredg, Zpred[[j]])
   }
   
   }
}

if(intercept==FALSE){Xlinearpred <- Xlinearpred[,-1]}
if(spli&li) {
  mu_gamGpred <- mu_gamG
  if("random"%in%rownames(z)){mu_gamGpred[seq(z[rownames(z)=="random",1],z[rownames(z)=="random",2])]<-0}
prediction <- Zpredg%*%mu_gamGpred+ Xlinearpred%*%mu_bet}else{
  if(spli){
  mu_gamGpred <- mu_gamG
  if("random"%in%rownames(z)){mu_gamGpred[seq(z[rownames(z)=="random",1],z[rownames(z)=="random",2])]<-0}
    prediction <- Zpredg%*%mu_gamGpred
  }else{
    if(li){
      prediction <-  Xlinearpred%*%mu_bet
    }
  }
}
}else{
  prediction <-y-res
  
}
if(spli){
  mu_gam <- list()
  sig_gamma <- list()
for(j in 1:nrow(z)){
mu_gam[[j]] <- mu_gamG[seq(z[j,1],z[j,2])]
names(mu_gam)[j] <- names(Z)[j]
sig_gamma[[j]]<- sig_gamG[seq(z[j,1],z[j,2]),seq(z[j,1],z[j,2])]
names(sig_gamma)[j] <- names(Z)[j]
if(names(mu_gam)[j]=="geo"){names(Z)[j]<-"geo"
                            gramap <- K[[j]]}
if(names(mu_gam)[j]=="random"){names(Z)[j]<- "random"}

 }
}else{mu_gam <- NULL
       sig_gamma <-NULL
  }
if(li){
rownames(mu_bet)<-colnames(Xlinear)}
  
###########################################
## Rückgabe der Werte  ####################
###########################################	  
    if(is.null(pred)){prediction <- Zg%*%mu_gamG+Xlinear%*%mu_bet}  
	l<- list(
    "pqr"<-vektor,
    "w"=wend,
    "lambda"=4*ta*(1-ta),
    "delli"=delli,
    "call"=formula,
    "beta"=mu_bet,   
    "Xlinear"=Xlinear,
    "sig_beta"=sig_bet,
	  "Xnonlinear"=Xnonlinear,
    "Z"=Z,
    "gamma"= mu_gam, 
    "sig_gamma" = sig_gamma,
    "knots"=knots,
    "mx"=mx,
		"delta" = mu_delt,
		"tau_inv" = mu_tau_inv,
	  #"K"=K,
    "i"=i,
    "residuals"= residuals,
    #"mse"=mse[1:i],
    "replace"=qua,
    "prediction"=prediction,
    "geo"=geo,
    "gramap"=gramap)
if(quant){
  l$tau <- ta
}
attr(l, "names")<- names(l)
attr(l, "class") <- "va"
		return(l)
	}
