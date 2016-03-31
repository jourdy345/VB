#### Loading R package
library(BSAM)

#### Data generation
set.seed(1)
ftn=function(x) 2-5*x+exp(5*(x-0.6))
n=100
xobs=(2*(1:n)-1)/(2*n)   				# nonlinear component
yobs=ftn(xobs)+rnorm(n,sd=0.8)
wdata=NULL								# fixed effects excluding intercept

# Priors  
iflagprior = 1;		# 1 = Lasso Smoother, theta[j] ~ N(0,sigma2*tau2*exp(-j*gamma)), tau2 ~ IGa(r_{0,\tau}/2,s_{0,\tau}/2)
					# 0 = T Smoother, theta[j] ~ N(0,sigma2*tau2*exp(-j*gamma)), tau2 ~ Exp(u0)
					
iflagpsi = 0		# 1 = estimate slope of squish function for S-shaped 
					# 0 = use fixed value of psi
psifixed = 1000		# Starting point for psi:

# Shape-restriction	(fmodel and fpm togles type of shape restriction)
fpm = 1     # If Free f.  fpm is forced to 1
fmodel = 1 	# Free f. 
			# fmodel = 2 and fpm =  1:	Increasing f
			# fmodel = 2 and fpm = -1: 	Decreasing f
			# fmodel = 3 and fpm =  1:  Increasing, convex f
			# fmodel = 3 and fpm = -1:	Decreasing, concave f
			# fmodel = 4 and fpm =  1:	Increasing, concave f
			# fmodel = 4 and fpm = -1:	Decreasing, convex f
			# fmodel = 5 and fpm =  1:	Increasing S or increasing, convex-to-concave f
			# fmodel = 5 and fpm = -1:	Decreasing S or decreasing, concave-to-convex f
			# fmodel = 6 and fpm =  1:	Increasing Rotated S or increasing, concave-to-convex f
			# fmodel = 6 and fpm = -1:	Decreasing Rotated S or decresaing, convex-to-concave f
			# fmodel = 7 and fpm =  1:  Inverted U: increasing to omega and decreasing after omega
			# fmodel = 7 and fpm = -1:	U Shapes: decreasing to omega and increasing after omega

# mcmc parameters
nblow0 = 1000;		# Initialization period for adpative metropolis
nblow = 10000;		# Number of MCMC in transition period 	
smcmc = 1000;		# Number of MCMC for analysis 			
nskip = 50;			# Number of MCMC to skip after nblow    
ndisp = 1000;		# Print on screen
maxmodmet = 5;		# Maximum number of times to modify metropolis for shape-restricted function

nint = 200;			# number of intervals for plotting f 
nbasis = 50;		# number of cosine basis functions   


# fit model
fit=bsar(yobs=yobs,xobs=xobs,wdata=wdata,fmodel=fmodel,fpm=fpm,iflagprior=iflagprior,nint=nint,nbasis=nbasis,
	     nblow0=nblow0,nblow=nblow,smcmc=smcmc,nskip=nskip,ndisp=ndisp,maxmodmet=maxmodmet,
	     psifixed=psifixed,iflagpsi=iflagpsi)
	     
print(fit) # summary
plot(fit)  # plot

