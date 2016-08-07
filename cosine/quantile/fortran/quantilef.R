dyn.load('quantilef.so')

wrapFortranMCMC <- function(y,x,W,priors,J,quant,nSample,burnIn,thinIn)
{
  n <- length(y)
  p <- ncol(W)
  if (!is.matrix(W)) W <- as.matrix(W)
  A             <- priors$A
  B             <- priors$B
  muBeta        <- priors$muBeta
  SigmaBeta     <- priors$SigmaBeta
  w0            <- priors$w0
  if (any(x<0 | x>1)) {
    xmax   <- max(x)
    xmin   <- min(x)
    x_adj  <- (x-xmin)/(xmax-xmin)
    varphi <- sqrt(2/(xmax-xmin))*cos(outer(x_adj,pi*(1:J)))
  } else {
    varphi <- sqrt(2)*cos(outer(x,pi*(1:J)))
  }
  res <- .Fortran('quantilef',y=as.matrix(y),vphi=as.matrix(varphi),
                  W=as.matrix(W),A=as.double(A),B=as.double(B),muBeta=as.vector(muBeta),
                  SigmaBeta=as.matrix(SigmaBeta),w0=as.double(w0),nbasis=as.integer(J),
                  quant=as.double(quant),nSample=as.integer(nSample),burnIn=as.integer(burnIn),
                  thinIn=as.integer(thinIn),beta=vector(mode='numeric',length=p),
                  theta=vector(mode='numeric',length=J),gamma=numeric(1),tau=numeric(1),n=as.integer(n),p=as.integer(p))
  res
}

N <- 1000
x <- runif(N)
y <- sin(2*pi*x)-log(x)
A         <- 0.01
B         <- 0.01
muBeta    <- 0
SigmaBeta <- matrix(1)
w0        <- 1
W         <- matrix(1,nr=length(y),nc=1)
priors    <- list(A=A,B=B,muBeta=muBeta,SigmaBeta=SigmaBeta,w0=w0)
MCres     <- wrapFortranMCMC(y,x,W,priors,30,0.5,2000,5000,30)