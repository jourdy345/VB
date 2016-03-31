source('vbfunctions.R')
library(BSAR)
source('vb_con.R')

is.not.null <- function(x) ! is.null(x)
fitfun<-function(psi,mu) {
  return(sum(mu*(psi%*%mu))) 
}

fitfun2<-function(psi,sigma) {
  return(sum(diag(psi%*%sigma))) 
}

fitfun3<-function(psi,sigma) {
  return(sum(diag(psi%*%sigma%*%psi%*%sigma))) 
}

fitfun4<-function(psi,sigma,mu) {
  return(sum(mu*(psi%*%sigma%*%psi%*%mu))) 
}

for (Modelnumber in 1:4)
{
maxloop = 50
iterall = rep(0,maxloop)
VB_RMSE = rep(0,maxloop)
VB_RMSEy = rep(0,maxloop)
VB_time = rep(0,maxloop)
VBSC_RMSE = rep(0,maxloop)
VBSC_RMSEy= rep(0,maxloop)
VBSC_time = rep(0,maxloop)
VBSCCon_RMSE = rep(0,maxloop)
VBSCCon_RMSEy= rep(0,maxloop)
VBSCCon_time = rep(0,maxloop)
MCMC_RMSE = rep(0,maxloop) 
MCMC_RMSEy= rep(0,maxloop) 
MCMC_time = rep(0,maxloop)
n<-50
T<-30

for (loop in 1:maxloop)
{

set.seed(Modelnumber*2+5*loop+n)

############################
#### Generating dataset ####
############################

if(Modelnumber ==1)
{
	ftn = function(x) exp(6*x-3)
delta_sc <- 1
delta_con <- 1
# fmodel = 3 and fpm =  1:  Increasing, convex f
fmodel = 3
fpm = 1
swap = 0
}
if(Modelnumber ==2)
{
	ftn = function(x) 16*x^2 - 4/(pi^2)*cos(2*pi*x) - 1/(pi^2)*cos(4*pi*x) - 32/(9*pi^2)*cos(3*pi*x) - 32/(pi^2)*cos(pi*x) + 365/(9*pi^2)
delta_sc <- 1
delta_con <- 1
# fmodel = 3 and fpm =  1:  Increasing, convex f
fmodel = 3
fpm = 1
swap = 0
}
if(Modelnumber ==3)
{

	ftn = function(x) log(1 + 10*x)
delta_sc <- 1
delta_con <- -1
# fmodel = 4 and fpm =  1:	Increasing, concave f
fmodel = 4
fpm = 1
swap = 1
}
if(Modelnumber ==4)
{
	ftn = function(x) sqrt(1 - x^2)
delta_sc <- -1
delta_con <- -1
# fmodel = 3 and fpm = -1:	Decreasing, concave f
fmodel = 3
fpm = -1
swap = 0
}

#x<-runif(n)
x <- (2*(1:n)-1)/(2*n) 
Z<-matrix(1,nrow=n,ncol=1)

tol<-0.0001

y = ftn(x) + rnorm(n)

####################################
#### Start of VC Unconstrainted ####
####################################

# Set up prior parameters

sigma2_m0 = 1;				   # @ Prior mean of sigma2 @
sigma2_v0 = 1000;			   # @ Prior var  of sigma2 @
sigma2_r0 = 2*(2 + sigma2_m0^2/sigma2_v0);
sigma2_s0 = sigma2_m0*(sigma2_r0-2);

rsig.0<-sigma2_r0
ssig.0<-sigma2_s0

tau2_m0	= 1;	# Make the prior mean too small to avoid the "dead zone" of IG @
tau2_v0 = 100;
rtau.0<- 2*(2 + tau2_m0^2/tau2_v0)
stau.0<- tau2_m0*(rtau.0 -2)	
w0<-2

mubeta.0<-0
sigbeta.0<-matrix(100,nrow=1,ncol=1)
sig02=100^2

prior.parms<-list(rsig.0=rsig.0,ssig.0=ssig.0,rtau.0=rtau.0,stau.0=stau.0,w0=w0,mubeta.0=mubeta.0,sigbeta.0=sigbeta.0,sig02=sig02)


VB_time[loop] <- system.time(test<-vbgpspectral(y,x,Z,T,tol,prior.parms,1))[3]

# y = Z\beta + \phi \theta + \epsilon 

vphi<-sqrt(2)*cos(outer(x,pi*(1:T)))  #This is eta_j (x) = \sqrt{2} cos(\pi x_j)
ord<-order(x)
plot(x,y,col="red")
fits<-(vphi[,1:length(test$mutheta.q)]%*%test$mutheta.q)
fits <- fits + Z*test$mubeta.q
#fits<-fits-mean(fits)
#lines(x[ord],fits[ord],col="purple")


true = ftn(x[ord])
lines(x[ord],true,col="black")
VB_RMSE[loop] = sqrt(1/n*sum((true - fits[ord])^2))
VB_RMSEy[loop] = sqrt(1/n*sum((y[ord] - fits[ord])^2))

####################################
#### Start of VC Shape Standard ####
####################################

nbasis = T
    mutheta.q.start<-rep(0,times=(nbasis+1))
    mutheta.q.start[1] <- 1
n.grid<-n

# Fit the model

VBSC_time[loop] <- system.time(test<- vbgpspectralsc(y,x,Z,T,tol,prior.parms,delta=delta_sc,mupsi.q.start=0.5,mutheta.q.start=mutheta.q.start,n.grid,stepsize=0.5))[3]

term1 <- apply(test$dmats,1,fitfun2,test$sigtheta.q)
term2 <-apply(test$dmats,1,fitfun,test$mutheta.q) #theta^T phi(x_i) theta (1 to n)
term3 <- Z%*%(test$mubeta.q)
fits <- delta_sc*(term1 + term2) + term3
plot(x,y,col="red")
ord<-order(test$x.grid) 
lines(test$x.grid[ord],fits[ord],col="purple")
test$mubeta.q

true = ftn(test$x.grid[ord])
lines(test$x.grid[ord],true,col="black")
VBSC_RMSE[loop] = sqrt(1/n*sum((true - fits[ord])^2))
VBSC_RMSEy[loop] = sqrt(1/n*sum((y[ord] - fits[ord])^2))

### Calculation of Standard Deviation ###
sdterm1 <- 2*apply(test$dmats,1,fitfun3,test$sigtheta.q)
sdterm2 <- 4*apply(test$dmats,1,fitfun4,test$sigtheta.q,test$mutheta.q)
sdterm3 <- (Z^2)*test$sigbeta.q[1,1]
sd_standard <- sqrt(sdterm1 + sdterm2 + sdterm3)

##################################
#### Start of VC Shape Convex ####
##################################
nbasis = T
    mutheta.q.start<-rep(0,times=(nbasis+1))
    mutheta.q.start[1] <- 1
n.grid<-n

# Fit the model

### Modificaton for convex/concave ###
Z<-matrix(1,nrow=n,ncol=1)
mualpha.0 = 3
sigalpha.0 = 50^2
mutheta.q.start <- c(0.5,mutheta.q.start)
######################################
alpha02 = sig02

prior.parms<-list(rsig.0=rsig.0,ssig.0=ssig.0,rtau.0=rtau.0,stau.0=stau.0,w0=w0,mubeta.0=mubeta.0,sigbeta.0=sigbeta.0,sig02=sig02,alpha02=alpha02)

test = NULL
mupsi.q.start = 0.5
while(is.null(test))
{
VBSCCon_time[loop] <- system.time(try(test<- vbgpspectral_con(y,x,Z,T,tol,prior.parms,delta=delta_con,mupsi.q.start=mupsi.q.start,mutheta.q.start=mutheta.q.start,n.grid,stepsize=0.5,swap=swap)))[3]
mupsi.q.start = mupsi.q.start + 0.5
mutheta.q.start[1] = mutheta.q.start[1] + 0.5
}

term1 <- apply(test$dmats,1,fitfun2,test$sigtheta.q)
term2 <-apply(test$dmats,1,fitfun,test$mutheta.q) #theta^T phi(x_i) theta (1 to n)
term3 <- Z%*%(test$mubeta.q)
fits <- delta_con*(term1 + term2) + term3
ord<-order(test$x.grid) 
lines(test$x.grid[ord],fits[ord],col="blue")

true = ftn(test$x.grid[ord])
lines(test$x.grid[ord],true,col="black")
VBSCCon_RMSE[loop] = sqrt(1/n*sum((true - fits[ord])^2))
VBSCCon_RMSEy[loop] = sqrt(1/n*sum((y[ord] - fits[ord])^2))

### Calculation of Standard Deviation ###
sdterm1 <- 2*apply(test$dmats,1,fitfun3,test$sigtheta.q)
sdterm2 <- 4*apply(test$dmats,1,fitfun4,test$sigtheta.q,test$mutheta.q)
sdterm3 <- (Z^2)*test$sigbeta.q[1,1]
sd_con <- sqrt(sdterm1 + sdterm2 + sdterm3)

#######################
#### Start of MCMC ####
#######################

xmin = 0
xmax = 1
iflagprior = 0		# 1 = Lasso Smoother, 0 = T Smoother

nblow0 = 1000		# Initialization period for adpative metropolis and adjust metm @
nblow = 10000		# Number of MCMC in transition period 	
smcmc = 1000		# Number of MCMC for analysis 			
nskip = 10			# Number of MCMC to skip after nblow    
ndisp = 200			# Giving the number of saved posterior samples to be displayed teon screen 

nint = n			# number of intervals for plotting f 
nbasis = T			# number of cosine basis functions   


fit = NULL
while(is.null(fit))
{
MCMC_time[loop] = system.time(try(fit <- BSAR(y,x,wdata=NULL,fmodel=fmodel,fpm=fpm,xmin=xmin,xmax=xmax,iflagprior=iflagprior,
		 nint=nint,nbasis=nbasis,nblow0=nblow0,nblow=nblow,smcmc=smcmc,nskip=nskip,ndisp=ndisp)))[3]
}
plot(x,y)
 lines(x,fit$yhatm)

MCMC_RMSE[loop] = sqrt(1/n*sum((true - fit$yhatm[ord])^2))
MCMC_RMSEy[loop] = sqrt(1/n*sum((y[ord] - fit$yhatm[ord])^2))

cat("Current loop is",loop,"\n")
}

	if(Modelnumber ==1)
	{
AllRMSE1 <- rbind(VB_RMSE,VBSC_RMSE,VBSCCon_RMSE,MCMC_RMSE)
AllRMSE1y <- rbind(VB_RMSEy,VBSC_RMSEy,VBSCCon_RMSEy,MCMC_RMSEy)
Alltime1 <- rbind(VB_time,VBSC_time,VBSCCon_time,MCMC_time)
apply(AllRMSE1,1,mean)
apply(Alltime1,1,mean)
	}
	if(Modelnumber ==2)
	{
AllRMSE2 <- rbind(VB_RMSE,VBSC_RMSE,VBSCCon_RMSE,MCMC_RMSE)
AllRMSE2y <- rbind(VB_RMSEy,VBSC_RMSEy,VBSCCon_RMSEy,MCMC_RMSEy)
Alltime2 <- rbind(VB_time,VBSC_time,VBSCCon_time,MCMC_time)
apply(AllRMSE2,1,mean)
apply(Alltime2,1,mean)
	}
	if(Modelnumber ==3)
	{
AllRMSE3 <- rbind(VB_RMSE,VBSC_RMSE,VBSCCon_RMSE,MCMC_RMSE)
AllRMSE3y <- rbind(VB_RMSEy,VBSC_RMSEy,VBSCCon_RMSEy,MCMC_RMSEy)
Alltime3 <- rbind(VB_time,VBSC_time,VBSCCon_time,MCMC_time)
apply(AllRMSE3,1,mean)
apply(Alltime3,1,mean)
	}
	if(Modelnumber ==4)
	{
AllRMSE4 <- rbind(VB_RMSE,VBSC_RMSE,VBSCCon_RMSE,MCMC_RMSE)
AllRMSE4y <- rbind(VB_RMSEy,VBSC_RMSEy,VBSCCon_RMSEy,MCMC_RMSEy)
Alltime4 <- rbind(VB_time,VBSC_time,VBSCCon_time,MCMC_time)
apply(AllRMSE4,1,mean)
apply(Alltime4,1,mean)
	}
#save.image("convex_n50.Rdata")
}

apply(AllRMSE1,1,mean)
apply(AllRMSE1,1,sd)/sqrt(50)
apply(Alltime1,1,mean)
apply(AllRMSE2,1,mean)
apply(AllRMSE2,1,sd)/sqrt(50)
apply(Alltime2,1,mean)
apply(AllRMSE3,1,mean)
apply(AllRMSE3,1,sd)/sqrt(50)
apply(Alltime3,1,mean)
apply(AllRMSE4,1,mean)
apply(AllRMSE4,1,sd)/sqrt(50)
apply(Alltime4,1,mean)

