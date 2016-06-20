sse=lmarg=numeric(6)
names(sse)=names(lmarg)=c('VBU','VBM','VBMC','BSAR','BSARM','BSARMC')
  
setEPS()
postscript('CI_VB.eps')

mtest <- matrix(c(1,1,1,1,2,2,3,3,0,4,4,0),nrow = 3,ncol = 4,byrow = TRUE)

layout(mat = mtest,heights = c(1,6,6))
  par(mar = c(0,0,0,0))
plot.new()
par(xpd=TRUE)

text = c('Residual','BSAR','VB','95% CI (BSAR)','95% CI (VB)')

legend("center", text.width = max(sapply(text, strwidth))+0.02, legend=text,density=c(0,0,0,25,0),
		col=c('gray60',1,'gray30',1,'gray60'),lty=c(1,2,1,NA,NA),pch=c(20,NA,NA,NA,15),lwd=c(NA,2,2,1.5,1.5),cex=1.2,border=c(NA,NA,NA,1,NA),
, horiz = TRUE,pt.cex=2, x.intersp=0.75,xjust=0.5)

par(xpd=FALSE)

  par(mar =c (5, 4, 3, 2) + 0.1)
  par(mai =c (0.5,0.5,0.5,0.5))
Modeltype = 1
source('VB_Elec.R')
sse[1]=sqrt(mean((ZB+fx-y)^2))
sse[4]=sqrt(mean((semi.fit$yhat$mean-y)^2))
lmarg[1]=semi.vb$lb[length(semi.vb$lb)]
lmarg[4]=semi.mcmc$lmarg.gd

plot(tempm[o],y-semi.fit$wbeta$mean,ylab='Parametric Residual',xlab='Temperature',pch=20,ylim=c(-0.3,0.35),cex=0.6,col='gray60')
polygon(c(tempm[o],rev(tempm[o])),c(CI$CIfx[1,],rev(CI$CIfx[2,])),border=NA,col='gray60')
polygon(c(tempm[o],rev(tempm[o])),c(semi.fit$fxobs$upper,rev(semi.fit$fxobs$lower)),density = c(10),lwd=1.5,col='gray10')

points(tempm[o],y-semi.fit$wbeta$mean,ylab='Parametric Residual',xlab='Temperature',pch=20,cex=0.6,col='gray60')
lines(tempm[o],fx,col='gray30',lwd=3,lty=1)
lines(tempm[o],semi.fit$fxobs$mean,lwd=3,col=1,lty=2)

main=expression(paste("(a) Unrestricted "))
title(main=main, line = 2)

Modeltype = 2
source('VB_Elec.R')
sse[2]=sqrt(mean((ZB+fx-y)^2))
sse[5]=sqrt(mean((semi.fit$yhat$mean-y)^2))
lmarg[2]=semi.vb$lb[length(semi.vb$lb)]
lmarg[5]=semi.mcmc$lmarg.gd

plot(tempm[o],y-semi.fit$wbeta$mean,ylab='Parametric Residual',xlab='Temperature',pch=20,ylim=c(-0.3,0.35),cex=0.6,col='gray60')
polygon(c(tempm[o],rev(tempm[o])),c(CI$CIfx[1,],rev(CI$CIfx[2,])),border=NA,col='gray60')
polygon(c(tempm[o],rev(tempm[o])),c(semi.fit$fxobs$upper,rev(semi.fit$fxobs$lower)),density = c(10),lwd=1.5,col='gray10')

points(tempm[o],y-semi.fit$wbeta$mean,ylab='Parametric Residual',xlab='Temperature',pch=20,cex=0.6,col='gray60')
lines(tempm[o],fx,col='gray30',lwd=3,lty=1)
lines(tempm[o],semi.fit$fxobs$mean,lwd=3,col=1,lty=2)

main=expression(paste("(b) Monotone "))
title(main=main, line = 2)

Modeltype = 3
source('VB_Elec.R')
sse[3]=sqrt(mean((ZB+fx-y)^2))
sse[6]=sqrt(mean((semi.fit$yhat$mean-y)^2))
lmarg[3]=semi.vb$lb[length(semi.vb$lb)]
lmarg[6]=semi.mcmc$lmarg.gd

plot(tempm[o],y-semi.fit$wbeta$mean,ylab='Parametric Residual',xlab='Temperature',pch=20,ylim=c(-0.3,0.35),cex=0.6,col='gray60')
polygon(c(tempm[o],rev(tempm[o])),c(CI$CIfx[1,],rev(CI$CIfx[2,])),border=NA,col='gray60')
polygon(c(tempm[o],rev(tempm[o])),c(semi.fit$fxobs$upper,rev(semi.fit$fxobs$lower)),density = c(10),lwd=1.5,col='gray10')

points(tempm[o],y-semi.fit$wbeta$mean,ylab='Parametric Residual',xlab='Temperature',pch=20,cex=0.6,col='gray60')
lines(tempm[o],fx,col='gray30',lwd=3,lty=1)
lines(tempm[o],semi.fit$fxobs$mean,lwd=3,col=1,lty=2)


main=expression(paste("(c) Monotone Convex/Concave "))
title(main=main, line = 2)

dev.off()
