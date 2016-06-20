# .libPaths("C:/Users/Physium/Documents/R/win-library/3.1")
# setwd("C:/Users/Physium/Dropbox/Edge Selection/shape_constrained_smoothing/Shape Restricted VB Jan 2016")

#Don't use Modeltype when using the Figure plot.R file. 
#Modeltype = 1  #Only run Unrestricted model
#Modeltype = 2  #Only run Monotone model
#Modeltype = 3  #Only run Monotone Convex/Concave

###### Functions ########

fitfun<-function(psi,mu) {
  return(sum(mu*(psi%*%mu))) 
}

fitfun2<-function(psi,mu) {
  return(sum(diag(psi%*%mu))) 
}

matfun6<-function(psi,sigma,mu) {
    return(2*sum(diag(psi%*%sigma%*%psi%*%sigma))+4*sum(mu*(psi%*%sigma%*%psi%*%mu)))
}


######################################
## ElecDemand Data in Lenk and Choi.##
##									##
## Ong et al. (2015)                ##
######################################
source('vb_unrestricted.R')
source('vb_monotone.R')
library(BSAM)
source('vb_con.R')
library(mvtnorm)
source('CI.R')
#====================================================================
# Loading Data
#====================================================================
library(xlsx)

Load=read.xlsx("ElecLoad.xls",sheetIndex=1)
nobs=nrow(Load)
lener      <- log(Load$ener)
lenerm     <- log(Load$enerm)
lpelec     <- log(Load$pelec)
lpgas      <- log(Load$pgas)
lgdp       <- log(Load$gdp)
temp       <- Load$cddq  - Load$hddq
tempm      <- Load$cddqm - Load$hddqm
tempmsq    <- tempm^2
lpelecpgas <- log(Load$pelec/Load$pgas)
lenermgdp  <- log(Load$enerm/Load$gdp)

o = order(tempm)
x = tempm[o];
y = lenermgdp[o];
n = length(y)

# @ Get wdata, if any @
w = lpelecpgas[o];
#====================================================================
# Fitting semiparametric model
#====================================================================
# Set up prior parameters
p=2
nbasis=30
T = 30
sigma2_m0=1
sigma2_v0=100
sigma2_r0=2*(2+sigma2_m0^2/sigma2_v0)
sigma2_s0=sigma2_m0*(sigma2_r0-2)

tau2_m0=1
tau2_v0=100
rtau.0 = 2*(2+tau2_m0^2/tau2_v0)
stau.0<- tau2_m0*(rtau.0 -2)	

prior.parms<-list(rsig.0=sigma2_r0,ssig.0=sigma2_s0,
                  rtau.0=tau2_m0,stau.0=tau2_v0,w0=2,
                  mubeta.0=rep(0,p),sigbeta.0=diag(100,p))


Z=cbind(1,w)
xmin=min(x)
xmax=max(x)
x=(x-xmin)/(xmax-xmin)


#####################################################
##### Fit the model using unconstrainted model ######
#####################################################

if(Modeltype == 1)
{

######################
##### MCMC BSAR ######

system.time({
semi.mcmc <- bsar(y=y,w=w,x=x,nbasis=30,shape='Free')
})
semi.mcmc
semi.fit=fitted(semi.mcmc)


############################
##### VB Unrestricted ######

T = nbasis*2

system.time({
semi.vb<-vbgpspectral(y=y,x=x,Z=Z,T=T,tol=0.0001,prior.parms,1)
})

#############################
##### Fitting the model #####

vphi<-sqrt(2)*cos(outer(x,pi*(1:nbasis)))  #This is eta_j (x) = \sqrt{2} cos(\pi x_j)
fx <- (vphi[,1:length(semi.vb$mutheta.q)]%*%semi.vb$mutheta.q)
ZB <-  Z%*%semi.vb$mubeta.q


#############################
##### Credible Interval #####

CI <- CI_noshape(semi.vb,vphi,Z)

}

###############################################
##### Fit the model using monotone model ######
###############################################

if (Modeltype == 2)
{

######################
##### MCMC BSAR ######

system.time({
semi.mcmc <- bsar(y=y,w=w,x=x,nbasis=30,shape='Decreasing')
})
semi.mcmc
semi.fit=fitted(semi.mcmc)

########################
##### VB MONOTONE ######

VBSC_time = 0
tol<-0.0001
delta<- -1 #Decreasing

prior.parms<-list(rsig.0=sigma2_r0,ssig.0=sigma2_s0,
                  rtau.0=tau2_m0,stau.0=tau2_v0,w0=2,
                  mubeta.0=rep(0,p),sigbeta.0=diag(100,p),sig02=10000)



## Setting for larger number of basis ##
mupsi.q.start<-2
n.grid<-n
T = nbasis*2
mutheta.q.start<-rep(1,times=(T+1))

semi.vb = NULL
while(is.null(semi.vb))
{
VBSC_time <- system.time(try(semi.vb<- vbgpspectralsc(y,x,Z,T,tol,prior.parms,delta,mupsi.q.start,mutheta.q.start,n.grid,stepsize=0.5,maxit=100)))[3]
mupsi.q.start<- mupsi.q.start + 0.5
}

#############################
##### Fitting the model #####

fxterm1 <- apply(semi.vb$dmats,1,fitfun2,semi.vb$sigtheta.q) #trace of sigtheta %*% phi(x_i)
fxterm2 <- apply(semi.vb$dmats,1,fitfun,semi.vb$mutheta.q) #theta^T phi(x_i) theta 
ZB <- (Z)%*%semi.vb$mubeta.q
fx <- delta*(fxterm1 + fxterm2) 


#############################
##### Credible Interval #####

CI <- CI_shape(semi.vb,Z,delta=delta)


#######################################
##### VB Decreasing convex model ######
#######################################

if (Modeltype == 3)
{
######################
##### MCMC BSAR ######

system.time({
semi.mcmc <- bsar(y=y,w=w,x=x,nbasis=30,shape='DecreasingConvex')
})
semi.mcmc
semi.fit=fitted(semi.mcmc)

################################
##### VB Decreasing Convex #####

delta_con <- 1
swap <- 1 #For decreasing convex and increasing concave, -1 otherwise

sig02=10000
mualpha.0 = 3
sigalpha.0 = 50^2
alpha02 = sig02
prior.parms<-list(rsig.0=sigma2_r0,ssig.0=sigma2_s0,
                  rtau.0=tau2_m0,stau.0=tau2_v0,w0=2,
                  mubeta.0=rep(0,p),sigbeta.0=diag(100,p),sig02=sig02,alpha02=alpha02)

T = nbasis*2
mutheta.q.start<-rep(5,times=(T+1))
mutheta.q.start <- c(5,mutheta.q.start)

semi.vb = NULL
VBSCCon_time = NULL
mupsi.q.start = 5
tol<-0.0001
n.grid<-n

VBSCCon_time <- system.time(try(semi.vb<- vbgpspectral_con(y,x,Z,T,tol,prior.parms,delta=delta_con,mupsi.q.start=mupsi.q.start,mutheta.q.start=mutheta.q.start,n.grid,stepsize=0.5,swap=swap,maxit=100)))[3]
VBSCCon_time

#############################
##### Fitting the model #####

fxterm1 <- apply(semi.vb$dmats,1,fitfun2,semi.vb$sigtheta.q)
fxterm2 <- apply(semi.vb$dmats,1,fitfun,semi.vb$mutheta.q) #theta^T phi(x_i) theta (1 to n)
ZB <- Z%*%(semi.vb$mubeta.q)
fx <- delta_con*(fxterm1 + fxterm2) 

#############################
##### Credible Interval #####

CI <- CI_shape(semi.vb,Z,delta=delta_con)


}

