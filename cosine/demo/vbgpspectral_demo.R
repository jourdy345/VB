source('vbgpspectral.R')

# Test example

set.seed(1)
n<-500
x<-0.1+0.8*runif(n)
#Z<-cbind(rep(1,times=n),runif(n))
#y<-sin(x*pi)+Z%*%c(1,1)+0.1*rnorm(n)
Z<-rep(1,times=n)
y<-sin(x*pi)+Z+0.1*rnorm(n)
T<-30
tol<-0.0001

# Set up prior parameters

rsig.0<-0.01
ssig.0<-0.01
rtau.0<-0.01
stau.0<-0.01
w0<-1
#mubeta.0<-c(0,0)
#sigbeta.0<-diag(2)
mubeta.0<-0
sigbeta.0<-matrix(1,nrow=1,ncol=1)
prior.parms<-list(rsig.0=rsig.0,ssig.0=ssig.0,rtau.0=rtau.0,stau.0=stau.0,w0=w0,mubeta.0=mubeta.0,sigbeta.0=sigbeta.0)

# Fit the model

test<-vbgpspectral(y,x,Z,T,tol,prior.parms,1)

vphi<-sqrt(2)*cos(outer(x,pi*(1:T)))
ord<-order(x)
plot(x,y-mean(y),col="red")
fits<-(vphi[,1:length(test$mutheta.q)]%*%test$mutheta.q)
fits<-fits-mean(fits)
lines(x[ord],fits[ord],col="purple")

#rv<-test$rsig.q/2
#sv<-test$ssig.q/2
#cg<-sv/(rv-1)
#lowlim<-max(0.001,cg-2*sqrt(cg^2/(rv-2)))
#uplim<-cg+2*sqrt(cg^2/(rv-2))
#grid<-seq(from=lowlim,to=uplim,length=1000)
#plot(grid,(rv*log(sv)-lgamma(rv)-(rv+1)*log(grid)-sv/grid),type="l")



