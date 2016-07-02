source('simCosineProbit.R')
source('../BSAR/Gbsar_v2.R')
source('../BSAR/Gbsar_probit_v2.R')
source('../BSAR/Gbsar_fun_Module_v2.R')
# source('Dropbox/MyOwn/GBSAR_family/GBSAR_VB/VB_probit_sim.R')  #VB
# source('Dropbox/MyOwn/GBSAR_family/GBSAR_v2/Gbsar_v2.R')
# source('Dropbox/MyOwn/GBSAR_family/GBSAR_v2/Gbsar_probit_v2.R')
# source('Dropbox/MyOwn/GBSAR_family/GBSAR_v2/Gbsar_fun_Module_v2.R')



set.seed(11)
n=300
xobs=runif(n)

f1=function(x)10*x-5
# f2=function(x)exp(x^2)-2
f3=function(x) 10*(x-0.2)*(x-0.4)*(x-0.7)
f4 = function(x)pi*x-3*sin(pi*x)
f2 <- function(x) sin(2 * pi * x) * log(x / 5)

plot(xobs,f1(xobs))
plot(xobs,f2(xobs))
plot(xobs,f3(xobs))
plot(xobs,f4(xobs))

y1=rbinom(length(xobs),1,pnorm(f1(xobs)))
y2=rbinom(length(xobs),1,pnorm(f2(xobs)))
y3=rbinom(length(xobs),1,pnorm(f3(xobs)))
y4=rbinom(length(xobs),1,pnorm(f4(xobs)))

nblow0=1000 # B0 iteration adaptive MCMC for shape-restricted model
nblow=1000 # Burn in
smcmc=1000 # Save
nskip=10 # Thin
ndisp=1000 # Display number of iteration
mcmc=list(nblow0=nblow0,nblow=nblow,smcmc=smcmc,nskip=nskip,ndisp=ndisp)

set.seed(123)
fit_VB=sim_cosine_probit(f2,30,draw=F) # fit f2 function in VB
x=fit_VB$x          # draw x from VB
y=rbinom(length(x),1,pnorm(f2(x))) # generate y using x from VB

# fmodel 1: no shape restriction
# fmodel 2: monotone increasing / monotone decreasing
# fmodel 3: increasing convex / decresing concave)
# fmodel 4: increasing concave / decresing convex

# fpm = 1 (increasing), -1 (decreasing)
# fix fpm=1 without shape restriction
fit_Gbsar1=Gbsar(y,xobs=x,family='bernoulli',link='probit',method='DataAug',
           fmodel=1,fpm=1,nbasis=30,mcmc=mcmc) # no shape restriction
# fit_Gbsar2=Gbsar(y,xobs=x,family='bernoulli',link='probit',method='DataAug',
                # fmodel=2,fpm=1,nbasis=30,mcmc=mcmc) # monotone increasing
# fit_Gbsar3=Gbsar(y,xobs=x,family='bernoulli',link='probit',method='DataAug',
                 # fmodel=3,fpm=1,nbasis=30,mcmc=mcmc) # increasing convex

#plot estimated functions
xo=order(x)
plot(fit_VB$x[xo],f2(fit_VB$x[xo]),type='l',col=1,lwd=2,xlab='x',ylab='f(x)',main='sin(2*pi*x)*log(x/5)')
lines(fit_VB$x[xo],fit_VB$unknown_f[xo]-fit_VB$W%*%fit_VB$fit$mubq,col='blue',lwd=2)
lines(x[xo],fit_Gbsar1$post.est$muhatm[xo],col='orange',lwd=5)
# lines(x[xo],fit_Gbsar2$post.est$muhatm[xo],col='blue',lwd=2,lty=2)
# lines(x[xo],fit_Gbsar3$post.est$muhatm[xo],col='darkgreen',lwd=2,lty=6)

RMSE.VB=sqrt(crossprod(f2(fit_VB$x[xo])-fit_VB$unknown_f[xo])/length(x))
RMSE.GBSAR1=sqrt(crossprod(f2(fit_VB$x[xo])-fit_Gbsar1$post.est$muhatm[xo])/length(x))
# RMSE.GBSAR2=sqrt(crossprod(f2(fit_VB$x[xo])-fit_Gbsar2$post.est$muhatm[xo])/length(x))
# RMSE.GBSAR3=sqrt(crossprod(f2(fit_VB$x[xo])-fit_Gbsar3$post.est$muhatm[xo])/length(x))

legend('topleft', legend = c('True', paste('VB = ', round(RMSE.VB,4),sep=''),
                                      paste('GBSAR_free = ', round(RMSE.GBSAR1,4),sep='')),
              col=c(1,2),lwd=2,lty=c(1,1),cex=0.8)

# legend('topleft',legend=c('True',paste('VB = ',round(RMSE.VB,4),sep=''),
#                           paste("GBSAR_free = ",round(RMSE.GBSAR1,4),sep=''),
#                           paste("GBSAR_monotone = ",round(RMSE.GBSAR2,4),sep=''),
#                           paste('GBSAR_convex = ',round(RMSE.GBSAR3,4),sep='')),
#        col=c(1,2,'orange','blue','darkgreen'),lwd=2,lty=c(1,1,5,2,6),cex=0.8)


#plot estimated probabilities
plot(fit_VB$x[xo],pnorm(f2(fit_VB$x[xo])),type='l',col=1,lwd=2,xlab='x',ylab='Probability',
     main='exp(x^2)-2',ylim=c(0,1))
lines(fit_VB$x[xo],pnorm(fit_VB$unknown_f[xo]),col='blue',lwd=2,lty=2)
lines(x[xo],pnorm(fit_Gbsar1$post.est$muhatm[xo]),col=2,lwd=2)
# lines(x[xo],pnorm(fit_Gbsar2$post.est$muhatm[xo]),col='orange',lwd=2,lty=5)
# lines(x[xo],pnorm(fit_Gbsar3$post.est$muhatm[xo]),col='darkgreen',lwd=2,lty=6)
points(x,y,pch=19,col='grey')

SSE.VB=crossprod(pnorm(f2(fit_VB$x[xo]))-pnorm(fit_VB$unknown_f[xo]))/length(x)
SSE.GBSAR=crossprod(pnorm(f2(fit_VB$x[xo]))-pnorm(fit_Gbsar1$post.est$muhatm[xo]))/length(x)
# SSE.hannah=crossprod(pnorm(f2(fit_VB$x[xo]))-pnorm(fit_Gbsar$post.est$muhatm[xo]))/length(x)
# SSE.gam=crossprod(pnorm(f2(fit_VB$x[xo]))-pnorm(fit_Gbsar$post.est$muhatm[xo]))/length(x)

legend('topleft', legend = c('True', paste('VB = ', round(SSE.VB,4),sep=''),
                                      paste('GBSAR_free = ', round(SSE.GBSAR1,4),sep='')),
              col=c(1,2),lwd=2,lty=c(1,1),cex=0.8)

# legend('topleft',legend=c('True',paste('VB = ',round(SSE.VB,4),sep=''),
#                           paste("GBSAR_free = ",round(SSE.GBSAR,4),sep=''),
#                           paste("GBSAr_monotne = ",round(SSE.hannah,4),sep=''),
#                           paste('GBSAR_convex = ',round(SSE.gam,4),sep='')),
#        col=c(1,2,'orange','blue','darkgreen'),lwd=2,lty=c(1,1,5,2,6),cex=0.8)
