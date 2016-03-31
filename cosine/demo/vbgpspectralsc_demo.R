source('vbgpspectralsc.R')

# Test example

set.seed(1)
n<-500
x<-runif(n)
#Z<-cbind(rep(1,times=n),runif(n))
#y<- -sin(x*pi/2)+Z%*%c(1,1)+0.1*rnorm(n)
Z<-matrix(1,nrow=n,ncol=1)
y<-as.vector(-sin(x*pi/2)+Z+0.1*rnorm(n))
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
sig02<-100
prior.parms<-list(rsig.0=rsig.0,ssig.0=ssig.0,rtau.0=rtau.0,stau.0=stau.0,w0=w0,mubeta.0=mubeta.0,sigbeta.0=sigbeta.0,sig02=sig02)
delta<- -1
mupsi.q.start<-1
n.grid<-200

# Fit the model

test<-vbgpspectralsc(y,x,Z,T,tol,prior.parms,delta,mupsi.q.start,n.grid)

# Plot fits

fitfun<-function(psi,mu) {
  return(sum(mu*(psi%*%mu))) 
}
fits<-apply(test$dmats,1,fitfun,test$mutheta.q)
plot(x,y-mean(y),col="red")
ord<-order(test$x.grid)
lines(test$x.grid[ord],delta*fits[ord],col="purple",lwd=2)

# fits<-apply(dmats,1,fitfun,mutheta.q)
# plot(x,y-mean(y),col="red")
# lines(x.grid,fits,col="purple")


