library(DPpackage)
library(gam)
library(gammSlice)

source('simCosineProbitMonotone.R')
source('../BSAR/GBSAR_v2/Gbsar_v2.R')
source('../BSAR/GBSAR_v2/Gbsar_probit_v2.R')
source('../BSAR/GBSAR_v2/Gbsar_fun_Module_v2.R')



set.seed(11)
n=300
xobs=runif(n)

f1=function(x)10*x-5
f2=function(x)exp(x^2)-2
f3=function(x) 10*(x-0.2)*(x-0.4)*(x-0.7)
f4 = function(x)pi*x-3*sin(pi*x)

plot(xobs,f1(xobs))
plot(xobs,f2(xobs))
plot(xobs,f3(xobs))
plot(xobs,f4(xobs))

y1=rbinom(length(xobs),1,pnorm(f1(xobs)))
y2=rbinom(length(xobs),1,pnorm(f2(xobs)))
y3=rbinom(length(xobs),1,pnorm(f3(xobs)))
y4=rbinom(length(xobs),1,pnorm(f4(xobs)))


n.burn=1000
n.thin=1
n.post=1000
mcmc=list(n.burn=n.burn,n.thin=n.thin,n.post=n.post)

set.seed(123)
fit_VB=simCosineProbitMonotone(f2,30,draw=F)
x=fit_VB$x
y=rbinom(length(x),1,pnorm(f2(x)))
xo=order(x)
fit_Gbsar=Gbsar(y,xobs=x,family='bernoulli',link='probit',method='DataAug',
           fmodel=1,fpm=1,nbasis=30)#mcmc=list(nblow=10000))
fit_Gbsar2=Gbsar(y,xobs=x,family='bernoulli',link='probit',method='DataAug',
                fmodel=2,fpm=1,nbasis=30)
fit_Gbsar3=Gbsar(y,xobs=x,family='bernoulli',link='probit',method='DataAug',
                 fmodel=3,fpm=1,nbasis=30)
fit.hannah=DPhannah.probit(y,x,20)
fhat=colMeans(fit.hannah$cdens)

save.image("C:/Users/Dynamic SOM/Dropbox/MyOwn/GBSAR_family/DPhannah/probittt.RData")


fit.gam=gam(y~s(x),family=binomial(link=probit))
prior <- list(taub1=2.02,taub2=0.02,beta0=rep(0,1),Sbeta0=diag(100,1),tau1=6.01,tau2=2.01)
mcmc <- list(nburn=1000,nsave=1000,nskip=1,ndisplay=1000)
fit.psgam=PSgam(formula=y3~ps(xobs,k=20,degree=3,pord=1),family=binomial(probit),prior=prior,mcmc=mcmc,
                ngrid=20,status=T)
cont=gSlc.control(nBurnin = 100,nIter = 100,nThin = 1,fixedEffPriorVar = 1e10,
                  sdPriorScale = 1e5)
fit.gamsl=gSlc(y~s(x),family='binomial',control=cont)



plot(fit_VB$x[xo],f2(fit_VB$x[xo]),type='l',col=1,lwd=2,xlab='x',ylab='f(x)',main='10(x-0.2)(x-0.4)(x-0.7)')
lines(fit_VB$x[xo],fit_VB$unknown_f[xo]-fit_VB$W%*%fit_VB$fit$mubq,col='blue',lwd=2,lty=2)
lines(x[xo],fit_Gbsar$post.est$muhatm[xo],col=2,lwd=2)
#lines(x[xo],fit_Gbsar2$post.est$muhatm[xo],col=2,lwd=2)
lines(x[xo],fhat[xo],col='orange',lwd=2,lty=5)
lines(x[xo],qnorm(fitted(fit.gam)[xo]),col='darkgreen',lwd=2,lty=6)

RMSE.VB=sqrt(crossprod(f3(fit_VB$x[xo])-fit_VB$unknown_f[xo])/length(x))
RMSE.GBSAR=sqrt(crossprod(f3(fit_VB$x[xo])-fit_Gbsar$post.est$muhatm[xo])/length(x))
#MSE.GBSAR2=crossprod(f3(fit_VB$x[xo])-fit_Gbsar2$post.est$muhatm[xo])/length(x)
RMSE.hannah=sqrt(crossprod(f3(fit_VB$x[xo])-fhat[xo])/length(x))
RMSE.gam=sqrt(crossprod(f3(fit_VB$x[xo])-qnorm(fitted(fit.gam)[xo]))/length(x))
RMSE.psgam=sqrt(crossprod())

legend('topleft',legend=c('True',paste("GBSAR = ",round(RMSE.GBSAR,4),sep=''),
                          paste("Hannah = ",round(RMSE.hannah,4),sep=''),
                          paste('VB = ',round(RMSE.VB,4),sep=''),
                          paste('GAM = ',round(RMSE.gam,4),sep='')),
       col=c(1,2,'orange','blue','darkgreen'),lwd=2,lty=c(1,1,5,2,6),cex=0.8)


plot(fit_VB$x[xo],pnorm(f3(fit_VB$x[xo])),type='l',col=1,lwd=2,xlab='x',ylab='Probability',
     main='10(x-0.2)(x-0.4)(x-0.7)',ylim=c(0,1))
lines(fit_VB$x[xo],pnorm(fit_VB$unknown_f[xo]),col='blue',lwd=2,lty=2)
lines(x[xo],pnorm(fit_Gbsar$post.est$muhatm[xo]),col=2,lwd=2)
lines(x[xo],pnorm(fhat[xo]),col='orange',lwd=2,lty=5)
lines(x[xo],fitted(fit.gam)[xo],col='darkgreen',lwd=2,lty=6)
points(x,y,pch=19,col='grey')

SSE.VB=crossprod(pnorm(f3(fit_VB$x[xo]))-pnorm(fit_VB$unknown_f[xo]))/length(x)
SSE.GBSAR=crossprod(pnorm(f3(fit_VB$x[xo]))-pnorm(fit_Gbsar$post.est$muhatm[xo]))/length(x)
SSE.hannah=crossprod(pnorm(f3(fit_VB$x[xo]))-pnorm(fhat[xo]))/length(x)
SSE.gam=crossprod(pnorm(f3(fit_VB$x[xo]))-fitted(fit.gam)[xo])/length(x)

legend('topleft',legend=c('True',paste("GBSAR = ",round(SSE.GBSAR,4),sep=''),
                          paste("Hannah = ",round(SSE.hannah,4),sep=''),
                          paste('VB = ',round(SSE.VB,4),sep=''),
                          paste('GAM = ',round(SSE.gam,4),sep='')),
       col=c(1,2,'orange','blue','darkgreen'),lwd=2,lty=c(1,1,5,2,6),cex=0.8)
