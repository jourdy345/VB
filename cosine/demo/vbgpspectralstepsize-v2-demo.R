source('vbgpspectralstepsize-v2.R')
source('vbgpspectral.R')

# Test example

set.seed(1)
#n<-50
#x<-runif(n)
#Z<-cbind(rep(1,times=n),runif(n))
#y<- -sin(x*pi/2)+Z%*%c(1,1)+0.1*rnorm(n)
#Z<-matrix(1,nrow=n,ncol=1)
#y<-as.vector(-sin(x*pi/2)+Z+0.1*rnorm(n))
#T<-30
#tol<-0.0001

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

#test<-vbgpspectralsc(y,x,Z,T,tol,prior.parms,delta,mupsi.q.start,n.grid)

# Plot fits

#fitfun<-function(psi,mu) {
#  return(sum(mu*(psi%*%mu))) 
#}
#fits<-apply(test$dmats,1,fitfun,test$mutheta.q)
#plot(x,y-mean(y),col="red")
#ord<-order(test$x.grid)
#lines(test$x.grid[ord],delta*fits[ord],col="purple")

#fits<-apply(dmats,1,fitfun,mutheta.q)
#plot(x,y-mean(y),col="red")
#lines(x.grid,fits,col="purple")

nsim=1
n=100      # sample size
nbasis=50  # number of basis

f1=function(x) 200*(x-0.1)*(x-0.6)*(x-0.8)
f2=function(x) 10*(exp(15*(x-0.4))/(exp(15*(x-0.4))+1))+exp(20*(x-0.9))
f3=function(x) 20*exp(5*(x-0.9))
f4=function(x) 20*exp(-15*(x-0.7)^2)
f5=function(x) 10*(exp(10*(x-0.4))/(exp(10*(x-0.4))+1))

gen.data=function(nsim,n,f,sig,seed,equal=T){
	if(!missing(seed)) set.seed(seed)
	dat=list()
	for(isim in 1:nsim){
		xy=list()
		if(equal) x=(2*(1:n)-1)/(2*n) 
		else x=runif(n)
		fx=f(x)
		xy$x=x
		xy$fx=fx
		xy$y=fx+rnorm(n,sd=sig)
		dat[[isim]]=xy
	}
	dat
}
dat=gen.data(nsim,n,f1,sig=1,seed=1,equal=F) # data generation

for(isim in 1:nsim){
	x=dat[[isim]]$x
	y=dat[[isim]]$y
	fx=dat[[isim]]$fx

      sigma2_m0=1
	sigma2_v0=1000
	sigma2_r0=2*(2+sigma2_m0^2/sigma2_v0)
	sigma2_s0=sigma2_m0*(sigma2_r0-2)
	
	tau2_m0=1
	tau2_v0=100
	tau2_r0=2*(2+tau2_m0^2/tau2_v0)
	tau2_s0=tau2_m0*(tau2_r0-2)
	tau2_rn=tau2_r0+nbasis


	prior.parms<-list(rsig.0=sigma2_r0,ssig.0=sigma2_s0,sig02=10000,
					  rtau.0=tau2_m0,stau.0=tau2_v0,w0=2,
					  mubeta.0=0,sigbeta.0=matrix(100,nrow=1,ncol=1))


	x1<-x[order(x)][-1]
      y1<-sqrt(length(y)*diff(y[order(y)]))
      fit4start=vbgpspectral(y=y1,x=x1,Z=rep(1,times=(n-1)),T=nbasis,tol=0.0001,prior.parms=prior.parms,
		 				mupsi.q.start=0.5)
      mutheta.q.start<-rep(0,times=(nbasis+1))
      mutheta.q.start[1:length(fit4start$mutheta.q)]<-fit4start$mutheta.q
			  
      fit4=vbgpspectralsc.step.v2(y=y,x=x,Z=rep(1,times=n),T=nbasis,tol=0.0001,prior.parms=prior.parms,
		 				delta=1,mupsi.q.start=1.0,mutheta.q.start=mutheta.q.start,n.grid=200,stepsize=0.5)

	fitfun<-function(psi,mu,sigma) sum(mu*(psi%*%mu))+sum(diag(psi%*%sigma))
	fit4.fits<-apply(fit4$dmats,1,fitfun,fit4$mutheta.q,fit4$sigtheta.q)+fit4$mubeta.q
	
	RMISEy[isim,'VB-BSARM']=sqrt(mean((y-fit4.fits)^2))
	RMISEf[isim,'VB-BSARM']=sqrt(mean((fx-fit4.fits)^2))
	
	cat('\nData set = ',isim,sep='','\n\n')
}

plot(x,y,col='gray60',pch=20)
lines(x[order(x)],fit4.fits[order(x)],lwd=2,col=5)
#lines(fit4$x.grid,fit4.fits,lwd=2,col=5)

