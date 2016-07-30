#--------------------------------------------------------------------------------------#
#--Generalized Inverse Gaussian distribution-------------------------------------------#
#--Usage: dgig(x, lambda = 1, chi = 1, psi = 1, logvalue=FALSE)                        #
#-------- pgig(q, lambda = 1, chi = 1, psi = 1, ...)                                   #
#-------- rgig(n, lambda = 1, chi = 1, psi = 1)                                        #
#-------- Egig(lambda, chi, psi, func = c('x', 'logx', '1/x', 'var'), check.pars=TRUE) #
#--Note:  Egig returns the expected value of each functional form ---------------------#
#---------aq2 = chi, bq2 = psi --------------------------------------------------------#
#--------------------------------------------------------------------------------------#
require('ghyp')
require('quantmod')

#---expectation of folded-normal---#
efnorm <- function(mu,sigma,sigma2) {
  sigma*sqrt(2/pi)*exp(-mu^2/(2*sigma2))+mu*(1-2*pnorm(-mu/sigma))
}

#-------MGF of folded-normal-------#
mfnorm <- function(mu,sigma,sigma2,t) {
  exp(sigma2*t^2/2+mu*t)*(1-pnorm(-mu/sigma-sigma*t))+exp(sigma2*t^2/2-mu*t)*(1-pnorm(mu/sigma-sigma*t))
}

#-------derivative of efnorm wrt mu--------#
dmefnorm <- function(mu,sigma) {
  2*pnorm(mu/sigma)-1
}

#-------derivative of efnorm wrt sigma--------#
dsefnorm <- function(mu,sigma) {
  dnorm(mu/sigma)
}

#-------error function--------#
erf <- function(x) 2*pnorm(x*sqrt(2))-1

#-------derivative of mfnorm wrt mu-------#
dmmfnorm <- function(mu,sigma,sigma2,t) {
  res <- -t*exp(t^2*sigma2/2-t*mu)*(0.5-0.5*erf((mu/sigma-t*sigma)/sqrt(2)))+
         t*exp(t^2*sigma2/2+t*mu)*(0.5+0.5*erf((mu/sigma+t*sigma)/sqrt(2)))-
         1/sqrt(2*pi*sigma2)*(exp(t^2*sigma2/2-0.5*(mu/sigma-t*sigma)^2-t*mu)-exp(t^2*sigma2/2-0.5*(t*sigma+mu/sigma)^2+t*mu))
  
  res
}

#-------derivative of mfnorm wrt sigma-------#
dsmfnorm <- function(mu,sigma,sigma2,t) {
  res <- t^2/2*exp(t^2*sigma2/2-t*mu)*(0.5-0.5*erf((mu/sigma-t*sigma)/sqrt(2)))+
         t^2/2*exp(t^2*sigma2/2+t*mu)*(0.5+0.5*erf((mu/sigma-t*sigma)/sqrt(2)))-
         1/sqrt(2*pi)*exp(t^2*sigma2/2-0.5*(mu/sigma-t*sigma)^2-t*mu)*(-t/(2*sigma)-mu/(2*sigma^3))+
         1/sqrt(2*pi)*exp(t^2*sigma2/2-0.5*(mu/sigma+t*sigma)^2+t*mu)*(t/(2*sigma)-mu/(2*sigma^3))

  res
}

#-----------lower bound------------#
lowerBound <- function(y,W,varphi,quant,J,A,B,Aq,Bq,aq2,bq2,muBeta,SigmaBeta_inv,muBetaq,SigmaBetaq,muThetaq,SigmaThetaq,mug,sigmag,sigmag2) {
  n     <- length(y)
  p     <- dim(W)[2]
  taup2 <- 2/(quant*(1-quant))
  nup   <- (1-2*quant)/(quant*(1-quant))
  t1    <- SigmaBeta_inv%*%SigmaBetaq


  res <- -n/2*log(2*pi*taup2) - 0.5*sum(1/taup2*Egig(rep(0.5,n),aq2,rep(bq2,n),func='1/x')*(y-W%*%muBetaq-varphi%*%muThetaq)^2+diag(W%*%tcrossprod(SigmaBetaq,W))+diag(varphi%*%tcrossprod(SigmaThetaq,varphi)))-
  sum(Egig(rep(0.5,n),aq2,rep(bq2,n),func='x')) - p/2*log(2*pi)+0.5*determinant(t1)$modulus[1]-0.5*sum((muBetaq-muBeta)*(SigmaBeta_inv%*%(muBetaq-muBeta)))-
  J/2*(log(2*pi)+log(Bq)-digamma(Aq))+0.25*J*(J+1)*efnorm(mug,sigmag,sigmag2)-0.5*Aq/Bq*sum(mfnorm(mug,sigmag,sigmag2,1:J)*(muThetaq^2+diag(SigmaThetaq)))+
  log(w0/2)-w0*efnorm(mug,sigmag,sigmag2)+A*log(B)-lgamma(A)-(A+1)*(log(Bq)-digamma(Aq))-B*Aq/Bq+
  0.5*log(sqrt(aq2/bq2))+log(2*besselK(sqrt(aq2*bq2),0.5,expon.scaled=FALSE))+0.5*(sum(aq2*Egig(rep(0.5,n),aq2,rep(bq2,n),func='1/x'))+bq2*sum(Egig(rep(0.5,n),aq2,rep(bq2,n),func='x')))+
  (p+J+1)/2*(1+log(2*pi))+0.5*determinant(SigmaThetaq)$modulus[1]+0.5*log(sigmag2)+Aq+log(Bq)+lgamma(Aq)-(1-Aq)*digamma(Aq)

  res
}

vbQuant <- function(y,x,W,priors,quant,mug.start,J=20,tol=1.0e-05) {
  if (is.matrix(W)==FALSE) W <- as.matrix(W)
  n             <- length(y)
  p             <- ncol(W)
  A             <- priors$A
  B             <- priors$B
  muBeta        <- as.vector(priors$muBeta)
  SigmaBeta     <- as.matrix(priors$SigmaBeta)
  w0            <- priors$w0
  varphi        <- sqrt(2)*cos(outer(x,pi*(1:J)))
  WtW           <- crossprod(W)
  VtV           <- crossprod(varphi)
  SigmaBeta_inv <- solve(SigmaBeta)
  SigmaBeta_inv_muBeta <- as.vector(SigmaBeta_inv%*%muBeta)
  taup2  <- 2/(quant*(1-quant))
  nup    <- (1-2*quant)/(quant*(1-quant))
  const1 <- 0.25*J*(J+1)-w0

  #---initialize variational parameters
  bq2     <- 2+nup^2/taup2
  aq2     <- rep(1,n)
  Aq      <- J/2+A
  Bq      <- B
  mug     <- mug.start
  sigmag2 <- mug^2/100
  sigmag  <- sqrt(sigmag2)
  muBetaq <- muBeta
  SigmaBetaq <- SigmaBeta
  muThetaq   <- rep(1,J)
  SigmaThetaq <- diag(1,J)

  #---create call objects---#
  clb   <- call('lowerBound',y,W,varphi,quant,J,A,B,Aq,quote(Bq),quote(aq2),bq2,muBeta,SigmaBeta_inv,muBetaq,quote(SigmaBetaq),quote(muThetaq),quote(SigmaThetaq),quote(mug),quote(sigmag),quote(sigmag2))
  clb.try <- call('lowerBound',y,W,varphi,quant,J,A,B,Aq,quote(Bq),quote(aq2),bq2,muBeta,SigmaBeta_inv,muBetaq,quote(SigmaBetaq),quote(muThetaq),quote(SigmaThetaq),quote(mug.try),quote(sigmag.try),quote(sigmag2.try))
  dif   <- tol+1
  lb    <- c()
  lbold <- -Inf
  lbnew <- 0
  count <- 0
  while (dif>tol | dif < 0) {
    count       <- count+1
    Estar       <- mfnorm(mug,sigmag,sigmag2,1:J)
    cat('Estar: ', Estar, '\n')
    aq2         <- 1/taup2*sum((y-(as.vector(W%*%muBetaq+varphi%*%muThetaq)))^2+diag(W%*%tcrossprod(SigmaBetaq,W))+diag(varphi%*%tcrossprod(SigmaThetaq,varphi)))
    e1overu     <- Egig(rep(0.5,n),aq2,rep(bq2,n),func='1/x')
    re1overu    <- sqrt(e1overu) # root of e1overu
    SigmaThetaq <- solve(crossprod(varphi*(re1overu%o%rep(1,ncol(varphi))))/taup2+Aq/Bq*diag(Estar))
    muThetaq    <- as.vector(SigmaThetaq%*%colSums(varphi*((e1overu*(y-as.vector(W%*%muBetaq))-nup)%o%rep(1,ncol(varphi)))/taup2))
    # e1overu     <- Egig(rep(0.5,n),aq2,rep(bq2,n),func='1/x')
    # re1overu    <- sqrt(e1overu) # root of e1overu
    SigmaBetaq  <- solve(crossprod(W*(re1overu%o%rep(1,ncol(W))))/taup2+SigmaBeta_inv)
    muBetaq     <- as.vector(SigmaBetaq%*%(colSums(W*((e1overu*as.vector((y-varphi%*%muThetaq)-nup))%o%rep(1,ncol(W))))/taup2+SigmaBeta_inv_muBeta))
    Bq          <- 0.5*(sum(muThetaq^2*Estar)+sum(Estar*diag(SigmaThetaq)))+B
    #---check lower bound---#
    lbtest      <- eval(clb)
    #---NCVMP for gamma---#
    sigmag2.old <- sigmag2
    sigmag2     <- -0.5/(const1*dsefnorm(mug,sigmag)-0.5*Aq/Bq*sum((muThetaq^2+diag(SigmaThetaq))*dsmfnorm(mug,sigmag,sigmag2,1:J)))
    sigmag      <- sqrt(sigmag2)
    cat('sigmag2: ', sigmag2, '\n')
    mug.old     <- mug
    mug         <- mug+sigmag2*(const1*dmefnorm(mug,sigmag)-0.5*Aq/Bq*sum((muThetaq^2+diag(SigmaThetaq))*dmmfnorm(mug,sigmag,sigmag2,1:J)))
    #---check lower bound---#
    lbnew       <- eval(clb)
    cat('lbnew1: ',lbnew,'\n')
    dif         <- lbnew-lbtest

    if (dif<0) {
      step <- 1
      dif.try <- dif
      while (dif.try<0) {
        step <- step*0.5
        sigmag2.try <- 1/(1/sigmag2.old+step*(1/sigmag2-1/sigmag2.old))
        sigmag.try  <- sqrt(sigmag2.try)
        mug.try     <- sigmag2.try*(mug.old/sigmag2.old+step*(mug/sigmag2-mug.old/sigmag2.old))
        lbnew       <- eval(clb.try)
        dif.try     <- lbnew-lbtest
      }
      sigmag  <- sigmag.try
      sigmag2 <- sigmag2.try
      mug     <- mug.try
    }
    dif <- (lbnew-lbold)/abs(lbnew)
    lbold <- lbnew
    lb <- c(lb,lbnew)
    cat('count: ', count, ', lbnew: ', lbnew, ', dif: ', dif, '\n')
  }
  list(lb=lb,muThetaq=muThetaq,SigmaThetaq=SigmaThetaq,Aq=Aq,Bq=Bq,muBetaq=muBetaq,SigmaBetaq=SigmaBetaq,mug=mug,sigmag2=sigmag2,aq2=aq2,bq2=bq2)
}


# x <- runif(1000)
# y <- sin(2*pi*x)-log(x)

london <- new.env()
#-----------from google------------#
status <- getSymbols('LON:HSBA',env=london,src="google",from=as.Date("2005-01-01"))
status <- getSymbols('LON:BARC',env=london,src="google",from=as.Date('2005-01-01'))
HSBC   <- get("LON:HSBA",envir=london)
BARC   <- get("LON:BARC",envir=london)

y      <- as.vector(HSBC[,4])
x2     <- as.vector(BARC[,4])

if (length(y)>length(x2)) {
  y  <- y[1:length(x2)]
} else if (length(x2)>length(y)) {
  x2 <- x2[1:length(y)]
}

xmax <- max(x2)
xmin <- min(x2)
x    <- (x2-xmin)/(xmax-xmin)


A           <- 1
B           <- 1
muBeta      <- 1
SigmaBeta   <- matrix(1)
w0          <- 1
W           <- matrix(1,nr=length(y),nc=1)
priors      <- list(A=A,B=B,muBeta=muBeta,SigmaBeta=SigmaBeta,w0=w0)
res         <- NULL
mug.start   <- 0.1
while(is.null(res)) {
  try(res   <- vbQuant(y,x,W,priors,0.5,mug.start,J=30))
  mug.start <- mug.start + 0.1
}