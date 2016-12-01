Normal_post <- function(obj1,obj2,ylim,title)
{
x = -1000:1000*0.002
obj <- obj1
n <- obj$n
yobs <- obj$y.obs
mean <- n/(n+1)*mean(yobs)
sd <- sqrt(1/(1+n))
hx <- dnorm(x,mean,sd)
if (missing(ylim))
{
ylim = c(min(hx),max(hx))
}

mean <- obj$mu
sd <- sqrt(obj$Sig)
hx <- dnorm(x,mean,sd)
plot(x, hx, type="l", lty=1, xlab=expression(paste(theta)),
  ylab="Density",col="blue",lwd=2,cex=2,ylim=ylim)

obj <- obj2
mean <- obj$mu
sd <- sqrt(obj$Sig)
hx <- dnorm(x,mean,sd)
lines(x, hx, type="l", lty=1, xlab="x value",
  ylab="Density",col="red",lwd=2,cex=2)

mean <- n/(n+1)*mean(yobs)
sd <- sqrt(1/(1+n))
hx <- dnorm(x,mean,sd)
lines(x, hx, type="l", lty=3, xlab="mu",
  ylab="Density", main=title,col="black",lwd=3,cex=2,ylim=ylim)
}

Normal_grad <- function(obj1,obj2,title)
{
ymax <- max(apply(obj1$grad.sd[,1,],2,sd),apply(obj2$grad.sd[,1,],2,sd))
ymin <- min(apply(obj1$grad.sd[,1,],2,sd),apply(obj2$grad.sd[,1,],2,sd))
obj <- obj1
plot(apply(obj$grad.sd[,1,],2,sd), main=title, xlab="No. of iterations",ylab="Standard deviation", type="l",ylim=c(ymin,ymax),col='blue',cex=2)
obj <- obj2
lines(apply(obj$grad.sd[,1,],2,sd), type="l",col='red',cex=2)
}

Normal_LB <- function(obj1,obj2,title,method,ylim)
{
if (missing(ylim))
{
ylim = c(min(obj1$LB)-0.05,max(obj1$LB)+0.05)
}

plot(obj1$LB[1,], main=title, xlab="No. of iterations",ylab="Lower bound",ylim=ylim,lty=1,type="l",col='blue',cex=2)
lines(obj2$LB[1,],lty=1,type="l",col='red',cex=2)
n = obj1$n
if(method ==1)
{
anaLB = -0.5*log(2*pi) - 1/(2*n)* sum(obj1$y.obs^2) - 1/(2*n)*log(n+1) + 1/(2*n*(n+1)) *sum(obj1$y.obs)^2
abline(h = anaLB,lwd=4)
}
if(method ==2)
{
	epsABC = (1 + obj1$eps)
#	anaLB = -0.5*log(2*pi) - (n-1)/(2*n)*log(epsABC^2) - 1/(2*n)*log(n+epsABC^2) - 1/(2*n*epsABC^2) * sum(obj1$y.obs^2) + sum(obj1$y.obs)^2/(2*epsABC^2*(epsABC^2+n)*n)
tau = 0.5
	anaLB = -0.5*log(2*pi) - 0.5*log(epsABC) - 1/(2*n*epsABC) * sum(obj1$y.obs^2) - 1/(2*n)*log(n/epsABC + 1) + mean(obj1$y.obs)^2*((n/epsABC)^2)/(2*n*(1 + n/(epsABC))) - tau/(2*n)
print(anaLB)
abline(h = anaLB,col='red',lwd=4)
tau = 0.1
	anaLB = -0.5*log(2*pi) - 0.5*log(epsABC) - 1/(2*n*epsABC) * sum(obj1$y.obs^2) - 1/(2*n)*log(n/epsABC + 1) + mean(obj1$y.obs)^2*((n/epsABC)^2)/(2*n*(1 + n/(epsABC))) - tau/(2*n)
print(anaLB)
abline(h = anaLB,col='blue',lwd=4)
}
}


Normal_LB_NULL <- function(obj1,obj2,title,method,ylim)
{
  if (missing(ylim))
  {
    ylim = c(min(obj1$LB)-0.05,max(obj1$LB)+0.05)
  }
  
  plot(obj1$LB[1,], main=title, xlab="No. of iterations",ylab="Lower bound",ylim=ylim,lty=1,type="l",col='blue',cex=2)
  lines(obj2$LB[1,],lty=1,type="l",col='red',cex=2)
  n = obj1$n
  if(method ==1)
  {
    anaLB = -0.5*log(2*pi) - 1/(2*n)* sum(obj1$y.obs^2) - 1/(2*n)*log(n+1) + 1/(2*n*(n+1)) *sum(obj1$y.obs)^2
    abline(h = anaLB,lwd=4)
  }
  if(method ==2)
  {
    epsABC = (1 + obj1$eps)
    #	anaLB = -0.5*log(2*pi) - (n-1)/(2*n)*log(epsABC^2) - 1/(2*n)*log(n+epsABC^2) - 1/(2*n*epsABC^2) * sum(obj1$y.obs^2) + sum(obj1$y.obs)^2/(2*epsABC^2*(epsABC^2+n)*n)
    tau = 0.5
    anaLB = -0.5*log(2*pi) - 0.5*log(epsABC) - 1/(2*n*epsABC) * sum(obj1$y.obs^2) - 1/(2*n)*log(n/epsABC + 1) + mean(obj1$y.obs)^2*((n/epsABC)^2)/(2*n*(1 + n/(epsABC))) - tau/(2*n)
    print(anaLB)
    abline(h = anaLB,col='black',lwd=4)
  }
}