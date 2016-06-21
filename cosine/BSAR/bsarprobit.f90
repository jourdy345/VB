!*****************************************************************************************
!* R-Package for Bayesian Spectral Analysis Regression (BSAR) with Shape Restrictions
!*	   in Probit regression model.
!*----------------------------------------------------------------------------------------
!* Created by Seongil Jo (joseongil@gmail.com)
!* Department of Statistics, Korea University
!* Date. June/6/2015.
!*****************************************************************************************
!*	Model:
!*		Y ~ Ber(mu), E(Y) = mu
!*		mu = g(w'beta + sum_k fk(xk)), g : Normal CDF
!*		xmin < x < xmax
!*		beta has the  intercept
!*		Normalize fk with int_a^b fk(x)dx = 0
!*	Define
!*		Z(x) 	= sum_j=0^J theta_j*phi_j(x)
!*		g^a(x) 	= int_0^x |Z(s)|^2 dx
!*			 	= theta' Phi^a(x)*theta
!*		where Phi^a(x) is matrix with 
!*		phi^a_{j,k}(x) = int_0^x phi_j(s) phi_k(s) ds in (j,k) 
!*		g^b(x)	= int_0^x g^a(s) ds
!*				= theta' Phi^b(x)*theta
!*		where Phi^b(x) is a matrix with
!*		phi^b_{j,k}(x) = int_0^x phi^a_{j,k}(s) ds
!*
!*	f has
!*		1. No shape restriction
!*			f(x) = Z(x)  
!*
!*		2. Increasing f
!*			f(x) = g^a(x) 
!*		2' Decreasing f
!*			f(x) = -g^a(x) 
!*
!*		3. Increasing convex f
!*			f(x) = g^b(x)  + alpha*(x-xmid)
!*			alpha > 0
!*		3' Decreasing, concave f
!*			f(x) = -g^b(x)  + alpha*(x-xmid)
!*			alpha < 0
!*
!*		4. Increasing concave f
!*			f(x) = -g^b(a+b-x)  + alpha*(x-xmid) 
!*			alpha > 0
!*		4' Decreasing convex f
!*			f(x) = g^b(a+b-x)  + alpha*(x-xmid)
!*			alpha < 0
!*
!*		5. Increasing, "S" shaped or increasing convex-to-concave
!*			f(x) = int_a^x int_a^s Z(t)^2 h1(t)dt ds  + alpha*(x-xmid) - xi*(x-xmid)
!*			h1(t) = {1-exp[psi*(x-omega)]}/{1+exp[psi*(x-omega)]} 
!*			psi  > 0, a < omega < b, alpha > 0
!*			xi = min(0,f'(x)) to make sure that f' > 0
!*
!*		5' Decreasing, "S" shaped or decreasing concave-to-convex
!*			f(x) = -int_a^x int_a^s Z(t)^2 h1(t)dt ds  + alpha*(x-xmid) - xi*(x-xmid)
!*			h1(t) = {1-exp[psi*(x-omega)]}/{1+exp[psi*(x-omega)]} 
!*			xi   = max(0,f'(x)) to make sure that f' < 0
!*			psi  > 0, a < omega < b, alpha < 0
!*
!*		6  Increasing concave-to-convex
!*			Concave before omega and convex after omega
!*			f(x) = int_a^x int_a^s Z(t)^2 h2(t)dt ds  + alpha*(x-xmid) - xi*(x-xmin)
!*			h2(t) = {exp[psi*(x-omega)]-1}/{exp[psi*(x-omega)]+1} 
!*			psi  > 0, a < omega < b, alpha > 0
!*			xi = min(0,f'(x)) to make sure that f' > 0
!*
!*		6'  Decreasing convex-to-concave
!*			Concave before omega and convex after omega
!*			f(x) = int_a^x int_a^s Z(t)^2 h2(t)dt ds  + alpha*(x-xmid) - xi*(x-xmin)
!*			h2(t) = {exp[psi*(x-omega)]-1}/{exp[psi*(x-omega)]+1} 
!*			psi  > 0, a < omega < b, alpha > 0
!*			xi = min(0,f'(x)) to make sure that f' > 0
!*
!*
!*		7	Concave increasing to decreasing or inverted U shaped
!*			Increasing before omega and decreasing after omega
!*			f(x) 	= int_0^x Z^2(s)h(s) ds - zeta + alpha_0 
!*			zeta 	= min(0,min int_0^x Z^2(s)h(s)ds)
!*			h(x) 	= {1- exp[psi*(x-omega)]}/{exp[psi*(x-omega)]+1} 
!*
!*		7'	Convex decreasing to increasing or U shaped
!*			f(x) 	= -int_0^x Z^2(s)h(s) ds + zeta + alpha_0 
!*			zeta 	= min(0,min int_0^x Z^2(s)h(s)ds)
!*			h(x) 	= {1 - exp[psi*(x-omega)]}/{exp[psi*(x-omega)]+1} 
!*
!*	Cosine basis 
!*	phi_0(x) = 1/sqrt((xmax-xmin)) for xmin < x < xmax
!*	phi_k(x) = sqrt(2/((xmax-xmin)))*cos(pi*k*(x-xmin)/((xmax-xmin)))
!*
!*	Scale invariant priors:
!*	beta 		~ N(b0,B0)
!*	alpha		~ N(m0,v0)I(delta*alpha>0)  	Truncated normal
!*	theta_0 	~ N(0,v0)					for free f
!*	theta_0		~ N(0,v0)I(theta_0 > 0)		for constrained f
!*	theta_k 	~ N(0,tau^2*exp(-gamma*k) 	for free f
!*	theta_k 	~ N(0,tau^2*exp(-gamma*k)  	for constrained f
!*
!*	Smoothing parameters tau and gamma
!*	Choice of two priors for tau2 
!*		1.) T Smoother: 	tau2 ~ IG(r0/2,s0/2) 
!*		2.) Lasso Smoother:	tau2 ~ Exp(u0)
!*	gamma ~ Exp(w0)
!*   Note: gamma = 0 for Lasso Prior 
!*	Note: posterior of tau and gamma have "banana" contours
!*	zeta = ln(tau2) - kbar*gamma is less dependent with gamma or log(gamma)
!*	I have not been able to fruitfully exploit this fact.
!*
!*	"S" models uses "squish" function (reparameterization of hyperbolic tangent)
!*	that depens on slope psi and location omega (inflection point for S).
!*	psi 	~ N(m0,v0)I(psi > 0)  Truncated Normal
!*       omega 	~ N(m0,v0)I(xmin < omega < xmax) and m0 = (xmin+xmax)/2
!*	
!*****************************************************************************************
!*	fmodel and fpm togles type of shape restriction
!*	  fmodel = 1: 			Free f.  fpm is forced to 1
!*	  fmodel = 2 and fpm =  1:	Increasing f
!*	  fmodel = 2 and fpm = -1: 	Decreasing f
!*	  fmodel = 3 and fpm =  1:  Increasing, convex f
!*	  fmodel = 3 and fpm = -1:	Decreasing, concave f
!*	  fmodel = 4 and fpm =  1:	Increasing, concave f
!*	  fmodel = 4 and fpm = -1:	Decreasing, convex f
!*	  fmodel = 5 and fpm =  1:	Increasing S or increasing, convex-to-concave f
!*	  fmodel = 5 and fpm = -1:	Decreasing S or decreasing, concave-to-convex f
!*	  fmodel = 6 and fpm =  1:	Increasing Rotated S or increasing, concave-to-convex f
!*	  fmodel = 6 and fpm = -1:	Decreasing Rotated S or decresaing, convex-to-concave f
!*	  fmodel = 7 and fpm =  1:  Inverted U: increasing to omega and decreasing after omega
!*	  fmodel = 7 and fpm = -1:	U Shapes: decreasing to omega and increasing after omega
!*****************************************************************************************
subroutine bsarprobit(yobs,wdata,xobs,nobs,nparw,nfun,nbasis,nint,fmodel,fpm,&
                      theta_m0,theta0_m0,theta0_s0,tau2_m0,tau2_v0,w0,beta_m0,beta_v0,&
                      alpha_m0,alpha_s0,psi_m0,psi_s0,psifixed,omega_m0,omega_s0,&
                      xinxgrid,xidelta,iflagprior,iflagpsi,&
                      maxmodmet,nblow0,nblow,smcmc,nskip,ndisp,&
                      zetag,tau2g,gammag,thetag,betag,alphag,psig,omegag,&
                      fxgridg,fxobsg,yestg,phatg,wbg,loglikeg,logpriorg,imodmetg,pmetg)
	use ToolsRfunf
	use bsamTools
	implicit none
	
	!input arguments
	integer,intent(in) :: nobs,nparw,nfun,nbasis,nint,fmodel(nfun),iflagprior,iflagpsi
	integer,intent(in) :: maxmodmet,nblow0,nblow,smcmc,nskip,ndisp,xinxgrid(nobs,nfun)
	real*8, intent(in) :: yobs(nobs),wdata(nobs,nparw),xobs(nobs,nfun),fpm(nfun)
	real*8, intent(in) :: theta_m0(nbasis+1),theta0_m0,theta0_s0,xidelta(nobs,nfun)
	real*8, intent(in) :: tau2_m0,tau2_v0,w0,beta_m0(nparw),beta_v0(nparw,nparw)
	real*8, intent(in) :: alpha_m0,alpha_s0,psi_m0,psi_s0,omega_m0,omega_s0,psifixed
	
	!output arguments
	integer,intent(out) :: imodmetg
	real*8,intent(out) :: zetag(smcmc,nfun),tau2g(smcmc,nfun),gammag(smcmc,nfun)
	real*8,intent(out) :: thetag(nbasis+1,nfun,smcmc),betag(smcmc,nparw)
	real*8,intent(out) :: alphag(smcmc,nfun),psig(smcmc,nfun),omegag(smcmc,nfun)
	real*8,intent(out) :: fxobsg(nobs,nfun,smcmc),fxgridg(nint+1,nfun,smcmc)
	real*8,intent(out) :: wbg(smcmc,nobs),yestg(smcmc,nobs),phatg(smcmc,nobs)
	real*8,intent(out) :: loglikeg(smcmc),logpriorg(smcmc),pmetg(nfun)
	
	!internal arguments
	real*8  :: stime,itime			! cpu time
	integer :: imcmc,isave,nmcmc	! iteration
	
	integer :: iloop,ifun,iobs
	integer :: quadfacts((nbasis+1)*(nbasis+2)/2,3),leftindex((nbasis+1)*(nbasis+2)/2)
	integer :: rightindex((nbasis+1)*(nbasis+2)/2),multfact((nbasis+1)*(nbasis+2)/2)
	integer :: intsimpfacts(nint+1),tmpA(nbasis+1,1),tmpB(1,nbasis+1)
	integer :: tmpI(nbasis+1,nbasis+1),tmpC(nbasis+1,nbasis+1)
	
	real*8 :: xmin(nfun),xmax(nfun),xrange(nfun),xmid(nfun),xdelta(nfun)
	real*8 :: wtw(nparw,nparw),wtwi(nparw,nparw),wdatat(nparw,nobs)
	real*8 :: xgrid(nint+1,nfun),xbar(nfun),xtx(nfun)
	real*8 :: xobs2(nobs),xgrid2(nint+1)
	
	integer :: k,nr
	real*8 :: kall(nbasis),kall0(nbasis+1),kbar,kvec(nbasis)
	
	! trig functions
	integer :: nfree,nfunconstraint,ifree,irest
	real*8,allocatable :: phixobsfree(:,:,:),phixobsfreet(:,:,:)
	real*8,allocatable :: phixgridfree(:,:,:),phi2(:,:,:)
	real*8,allocatable :: phixobs(:,:,:),phixgrid(:,:,:)
	
	! priors
	real*8 :: theta0_v0,theta0_v0i,theta0_v0im0,theta0_lnv0
	real*8 :: tau2_r0,tau2_s0,tau2_rn,tau2_u0,wk
	real*8 :: beta_v0i(nparw,nparw),beta_v0im0(nparw),beta_lnv0
	real*8 :: alpha_v0,alpha_v0i,alpha_v0im0,alpha_lnv0
	real*8 :: psi_v0i,psi_v0,psi_lnp0,psi_lnv0
	real*8 :: omega_v0i,omega_v0,omega_lnp0,omega_lnv0
	
	! Random walk Metropolis-Hastings
	integer :: imodmet,pmet(nfun),iflag_AM,icount_AM,pok
	real*8 :: metw,met_alpha,met_omega_alpha,met_psi_alpha
	real*8 :: theta_mAM(nbasis+1,nfun),theta_vAM(nbasis+1,nfun)
	real*8 :: metm(nfun),mets(nfun),metv(nfun),met_beta(nfun),met_beta_AM(nfun)
	real*8 :: met_mAM(nfun),met_vAM(nfun),met_omega_m(nfun),met_omega_beta(nfun)
	real*8 :: met_psi_m(nfun),met_psi_beta(nfun),met_var_all(nfun),ck
	
	character*65 acceptmessage
	character*6 acceptr
	character*1 func
	
	! Parameters
	real*8 :: tau(nfun),tau2(nfun),tau2i(nfun),gampar(nfun),lngampar(nfun)
	real*8 :: theta(nbasis+1,nfun),theta02(nfun),theta2(nbasis,nfun),beta(nparw),wb(nobs)
	real*8 :: theta0(nfun),gamvec0(nbasis+1,nfun),thv0(nbasis,nfun),thv00(nbasis+1,nfun)
	real*8 :: psi(nfun),omega(nfun),zeta(nfun),alpha(nfun),gamvec(nbasis,nfun)
	real*8 :: fxobs(nobs,nfun),fxgrid(nint+1,nfun),yest(nobs),ehat(nobs),zobs(nobs)
	
	! Rfunctions
	real*8 :: cdfnorm,rndnorm,ltnormrnd,rtnormrnd,binomrnd
	
	call cpu_time(stime) ! start time
	call rndstart()		 ! call R random number generators
	
	! Data
	wdatat=transpose(wdata)
	wtw=matmul(wdatat,wdata)
	call inverse(wtw,nparw,wtwi)
	xbar=sum(xobs,1)/dble(nobs)
	do ifun=1,nfun
	  xmin(ifun)=minval(xobs(:,ifun))
	  xmax(ifun)=maxval(xobs(:,ifun))
	  xmid(ifun)=(xmin(ifun)+xmax(ifun))/2.d0
	  xrange(ifun)=xmax(ifun)-xmin(ifun)
	  xdelta(ifun)=(xrange(ifun))/dble(nint)
	  xtx(ifun)=sum((xobs(:,ifun)-xmid(ifun))**2.d0)
	  xgrid(1,ifun)=xmin(ifun)
	  do iloop=2,nint+1
	    xgrid(iloop,ifun)=xgrid(iloop-1,ifun)+xdelta(ifun)
	  end do
	end do
	
	! for basis
	nr=(nbasis+1)*(nbasis+2)/2
	kall=(/ (dble(k),k=1,nbasis) /)
	kall0=(/ (dble(k),k=0,nbasis) /)
	kbar=sum(kall0)/dble(nbasis+1)
	kvec=w0/(w0+kall)
	
	! ************************************************************************************
	! * Vector used in simpson's integration
	! * 1, 4, 2, 4, 2, ..., 2, 4, 1
	! ************************************************************************************
	intsimpfacts(1)=1	
	intsimpfacts(nint+1)=1
	do iloop=2,nint,2
		intsimpfacts(iloop)=4
	end do
	do iloop=3,nint-1,2
		intsimpfacts(iloop)=2
	end do
	
	! ************************************************************************************
	! * With shape constraint, 
	! *	f(x) includes  theta'Psi^c(xi)*theta for c = a or b
	! *	where Psi^c(xi) is a matrix of "basis" functions from Z^2
	! * QuadMult is a fast way to compute theta'Psi^c(xi)*theta for all x obs or xgrid
	! *	It needs some helper variables 
	! *   Get indicies for theta and multiplcation factor to do vectorized quadratic forms 
	! *   QuadMult = sumc(multfact.*theta[leftindex].*vech(phi).*theta[rightindex])
	! *		multfact   = vector or 1 for theta_j^2 and 2 for theta_j*theta_k terms
	! *		leftindex  = row index for vech, which vectorizes lower diagonal matrix
	! *		rightindex = col index for vech, which vectorizes lower diagonal matrx
	! *  Functional call for Quadratic Form Multiplication: 
	! * 	    Q = QuadMult(x, qvech, quadfacts);	
	! *
	! *  x can be a vector, such as xobs or xgrid.
	! *  QuadMult computes quadartic form for all x without loops!
	! *  The multfact, leftindex, and rightindex work with vech, 
	! *  which vectorizes symmetric matrics
	! ************************************************************************************
	tmpC=2
	call Idiag(1,nbasis+1,tmpI)
	call Ivech(tmpC-tmpI,nbasis+1,nbasis+1,multfact)
	quadfacts(:,1)=multfact
	
	tmpB=1
	tmpA(:,1)=(/ (iloop,iloop=1,nbasis+1) /)
	call Ikron(tmpB,1,nbasis+1,tmpA,nbasis+1,1,tmpC)
	call Ivech(tmpC,nbasis+1,nbasis+1,leftindex)
	quadfacts(:,2)=leftindex
	
	tmpA=1
	tmpB(1,:)=(/ (iloop,iloop=1,nbasis+1) /)
	call Ikron(tmpA,nbasis+1,1,tmpB,1,nbasis+1,tmpC)
	call Ivech(tmpC,nbasis+1,nbasis+1,rightindex)
	quadfacts(:,3)=rightindex
	
	! ************************************************************************************
    ! * Random walk Metropolis-Hastings with t-errors for fmodel > 1
    ! *	pmet = P(accept candidate theta)
    ! *	If pmet is > 0.6, then step size is too small: increase metm and/or mets
    ! *	If pmet is < 0.3, then step size is too large: decrease metm and/or mets
    ! ************************************************************************************
	metw=0.5d0          ! weight*(IG factor) + (1-weight)*(Running Mean)
	metm=0.01d0         ! Mean of inverse gamma
	met_alpha=3.d0
	if (met_alpha.gt.2.d0) then
	  mets=metm/dsqrt(met_alpha-2.d0)
	  metv=mets**2.d0   ! Variance of inverse gamma
	else
	  mets=1.d0
	  metv=1.d0
	end if
	met_beta=(met_alpha-1.d0)*metm
	met_beta_AM=met_beta  ! updated IG parameter for adaptive metropolis
	met_var_all=metm      ! IG values
	met_mAM=metm
	met_vAM=metv
	
	! Metropolis for omega and psi
	met_omega_m=0.0001d0
	met_omega_alpha=3.d0
	met_omega_beta=(met_omega_alpha-1.d0)*met_omega_m
	
	met_psi_m=0.0001d0
	met_psi_alpha=3.d0
	met_psi_beta=(met_psi_alpha-1.d0)*met_psi_m
	
	! ************************************************************************************
	! * Prior distributions
	! ************************************************************************************
	! @ Smoothing Prior: signal to noise ratio @
	! @ Free f: theta ~ N(0,tau2*exp(-k*gamma)) 		for k = 1 ..., nbasis @
	! @ Constrained f: theta ~ N(0,tau^2*exp(-gamma*k)) for k = 1, .., nbasis @
	
	! @ theta0 ~ N(0,v0^2) for free f @
	! @ theta0 ~ N(0,v0^2)I(theta0 > 0)  for shape restricted f @
	theta0_v0=theta0_s0**2.d0
	theta0_v0i=1.d0/theta0_v0
	theta0_v0im0=theta0_v0i*theta0_m0
	theta0_lnv0=dlog(theta0_v0)
	
	! @ IG prior for tau2, the signal to noise ratio: define mean and variance @
	tau2_r0=2.d0*(2.d0+tau2_m0**2.d0/tau2_v0)
	tau2_s0=tau2_m0*(tau2_r0-2.d0)
	tau2_rn=tau2_r0+dble(nbasis)
	
	! @ Exponential prior for tau2 if Lasso prior @
	tau2_u0=1/tau2_m0
	
	! @ gamma ~ Exp(w0) @
	wk=sum(kall)/2.d0-w0
	
	! @ beta ~ N(bm0,bv0) @
	call inverse(beta_v0,nparw,beta_v0i)
	beta_v0im0=matmul(beta_v0i,beta_m0)
	beta_lnv0=dlog(determinant(beta_v0,nparw))
	
	! @ alpha ~ Truncated normal(m0, v0) @
	alpha_v0=alpha_s0**2.d0
	alpha_v0i=1.d0/alpha_v0
	alpha_v0im0=alpha_m0*alpha_v0i
	alpha_lnv0=dlog(alpha_v0)
	
	! @ psi is slope for logisitic function in S shaped f              @
	! @ psi ~ Truncated normal  psi is slope of squish (tanh) function @
	psi_v0=psi_s0**2.d0
	psi_v0i=1.d0/psi_v0
	psi_lnp0=cdfnorm(-psi_m0/psi_s0,0.d0,1.d0,0,1)   ! @ Log Normalizing constant @
	psi_lnv0=dlog(psi_v0)
	
	! @ omega is point of inflection for "S" shaped f @
	! @ omega ~ Truncated normal omega is zero of squish (tanh) function @
	omega_v0=omega_s0**2.d0
	omega_v0i=1.d0/omega_v0
	omega_lnv0=dlog(omega_v0)
	! @ Normal probability xmin < omega < xmax @
	omega_lnp0=dlog(cdfnorm((xmax-omega_m0)/(omega_s0),0.d0,1.d0,1,0)- &
					cdfnorm((xmin-omega_m0)/(omega_s0),0.d0,1.d0,1,0))

    ! ************************************************************************************
	! * Initialize parameters
	! ************************************************************************************
	! @ MLE of beta @
	tau2=tau2_m0
	tau=dsqrt(tau2)
	tau2i=1.d0/tau2
	psi=psifixed
	omega=xmid	
	gampar=1.d0/w0
	lngampar=dlog(gampar)
	do ifun=1,nfun
	  gamvec(:,ifun)=dexp(-gampar(ifun)*kall)
	  gamvec0(1,ifun)=1.d0
	  gamvec0(2:(nbasis+1),ifun)=gamvec(:,ifun)
	  thv0(:,ifun)=tau2(ifun)*gamvec(:,ifun)
	  thv00(1,ifun)=theta0_v0
	  thv00(2:(nbasis+1),ifun)=thv0(:,ifun)
	  zeta(ifun)=dlog(tau2(ifun))-kbar*gampar(ifun) ! zeta=ln(tau2)-meanc(kall)*gamma 
	  theta(1,ifun)=0.1d0
	  theta0(ifun)=theta(1,ifun)
	  theta02(ifun)=theta0(ifun)**2.d0
	  do k=2,nbasis+1
		  theta(k,ifun)=0.1d0*dsqrt(gamvec(k-1,ifun))*rndnorm()
		  theta2(k-1,ifun)=theta(k,ifun)**2.d0
	  end do
	end do
	beta=0.d0
	wb=matmul(wdata,beta)
	alpha=0.d0
	
	! ************************************************************************************
	! * Precomputes trig functions
	! ************************************************************************************
	nfree=count(fmodel.eq.1)          ! number of functions that have unconstraints
	nfunconstraint=count(fmodel.gt.1) ! number of functions that have constraints
	allocate(phixobsfree(nobs,nbasis,nfree),phixobs(nr,nobs,nfunconstraint))
	allocate(phixgridfree(nint+1,nbasis,nfree),phixgrid(nr,nint+1,nfunconstraint))
	allocate(phi2(nbasis,nbasis,nfree),phixobsfreet(nbasis,nobs,nfree))
	ifree=1
	irest=1
	do ifun=1,nfun
	  if (fmodel(ifun).eq.1) then
	    ! ********************************************************************************
	    ! * f(x) = Phi(x)'theta
	    ! * vector f = Phi*theta
	    ! * Compute Phi matrix at observations and integration grid
	    ! ********************************************************************************
	    call CosFun(xobs(:,ifun),xmin(ifun),xrange(ifun),nobs,nbasis, &
	                phixobsfree(:,:,ifree))
	    call CosFun(xgrid(:,ifun),xmin(ifun),xrange(ifun),nint+1,nbasis, &
	                phixgridfree(:,:,ifree))
	    phixobsfreet(:,:,ifree)=transpose(phixobsfree(:,:,ifree))
	    phi2(:,:,ifree)=matmul(phixobsfreet(:,:,ifree),phixobsfree(:,:,ifree)) 
	    
	    ! @ Initialize fxobs and fxgrid @
	    call GetFreef(theta(2:(nbasis+1),ifun),phixobsfree(:,:,ifree),phixgridfree(:,:,ifree),&
	                  nbasis,nobs,nint+1,fxobs(:,ifun),fxgrid(:,ifun))
	
	    ifree=ifree+1
	  else 
	    ! ********************************************************************************
	    ! * Precomputes trig functions
	    ! *	g^c(x) = theta'Phi^c(x)*theta for c = a or b
	    ! *	Compute Phi^c, Vectorize it with vech, Use Vector multiplication
	    ! ********************************************************************************
	    
	    ! Increasing f
	    if (fmodel(ifun).eq.2) then
	      call GetPhi(xobs(:,ifun),xmin(ifun),xrange(ifun),IntCos2,IntCosCrossProd, &
	                  IntConst2,IntCos,nbasis,nobs,phixobs(:,:,irest))
		  call GetPhi(xgrid(:,ifun),xmin(ifun),xrange(ifun),IntCos2,IntCosCrossProd, &
		              IntConst2,IntCos,nbasis,nint+1,phixgrid(:,:,irest))
		              
		  ! @ Initialize fxobs and fxgrid @
	      call GetUpf(fpm(ifun),theta(:,ifun),phixobs(:,:,irest),phixgrid(:,:,irest), &
	                  quadfacts,nbasis,nr,nobs,nint+1,fxobs(:,ifun),fxgrid(:,ifun))
	                   
		end if
		
		! Increasing and convex f
	    if (fmodel(ifun).eq.3) then
	      call GetPhi(xobs(:,ifun),xmin(ifun),xrange(ifun),IntIntCos2,IntIntCrossProd, &
	                  IntIntConst2,IntIntCos,nbasis,nobs,phixobs(:,:,irest))
		  call GetPhi(xgrid(:,ifun),xmin(ifun),xrange(ifun),IntIntCos2,IntIntCrossProd, &
		              IntIntConst2,IntIntCos,nbasis,nint+1,phixgrid(:,:,irest))
		  
		  ! @ Initialize fxobs and fxgrid @
	      call GetConvexf(fpm(ifun),alpha(ifun),theta(:,ifun),xobs(:,ifun),xgrid(:,ifun),&
	                      xmid(ifun),phixobs(:,:,irest),phixgrid(:,:,irest),quadfacts, &
					      nbasis,nr,nobs,nint+1,fxobs(:,ifun),fxgrid(:,ifun))
	    end if
	    
	    ! Increasing and concave f
	    if (fmodel(ifun).eq.4) then
	      xobs2=xmin(ifun)+xmax(ifun)-xobs(:,ifun)
	      xgrid2=xmin(ifun)+xmax(ifun)-xgrid(:,ifun)
	      call GetPhi(xobs2,xmin(ifun),xrange(ifun),IntIntCos2,IntIntCrossProd, &
	                  IntIntConst2,IntIntCos,nbasis,nobs,phixobs(:,:,irest))
	      call GetPhi(xgrid2,xmin(ifun),xrange(ifun),IntIntCos2,IntIntCrossProd, &
	                  IntIntConst2,IntIntCos,nbasis,nint+1,phixgrid(:,:,irest))
	                  
	      ! @ Initialize fxobs and fxgrid @
	      call GetConcavef(fpm(ifun),alpha(ifun),theta(:,ifun),xobs(:,ifun),xgrid(:,ifun),&
	                       xmid(ifun),phixobs(:,:,irest),phixgrid(:,:,irest),quadfacts, &
					       nbasis,nr,nobs,nint+1,fxobs(:,ifun),fxgrid(:,ifun))
		end if
		
		! "S" --> Increasing Convex to Concave f
		if (fmodel(ifun).eq.5) then
		  call GetPhi(xobs(:,ifun),xmin(ifun),xrange(ifun),CosFun2,CrossProdFun,&
		              ConstFun2,ConstCosFun,nbasis,nobs,phixobs(:,:,irest))
	      call GetPhi(xgrid(:,ifun),xmin(ifun),xrange(ifun),CosFun2,CrossProdFun,&
	                  ConstFun2,ConstCosFun,nbasis,nint+1,phixgrid(:,:,irest))
				
		! @ Initialize fxobs and fxgrid @
	      call GetSf(fpm(ifun),omega(ifun),psi(ifun),alpha(ifun),theta(:,ifun),xobs(:,ifun),&
	                 xgrid(:,ifun),phixobs(:,:,irest),phixgrid(:,:,irest),xdelta(ifun),&
	                 xinxgrid(:,ifun),xidelta(:,ifun),xrange(ifun),xmid(ifun),quadfacts,&
	                 intsimpfacts,nbasis,nr,nobs,nint+1,fxobs(:,ifun),fxgrid(:,ifun))
		end if
		
		! "Rotate S" --> Increasing Concave and Convex f
		if (fmodel(ifun).eq.6) then
		  call GetPhi(xobs(:,ifun),xmin(ifun),xrange(ifun),CosFun2,CrossProdFun, &
		              ConstFun2,ConstCosFun,nbasis,nobs,phixobs(:,:,irest))
	      call GetPhi(xgrid(:,ifun),xmin(ifun),xrange(ifun),CosFun2,CrossProdFun, &
	                  ConstFun2,ConstCosFun,nbasis,nint+1,phixgrid(:,:,irest))
				
		  ! @ Initialize fxobs and fxgrid @
	      call GetRotateSf(fpm(ifun),omega(ifun),psi(ifun),alpha(ifun),theta(:,ifun), &
	                       xobs(:,ifun),xgrid(:,ifun),phixobs(:,:,irest),phixgrid(:,:,irest),&
	                       xdelta(ifun),xinxgrid(:,ifun),xidelta(:,ifun),xrange(ifun), &
	                       xmid(ifun),quadfacts,intsimpfacts,nbasis,nr,nobs,nint+1,&
	                       fxobs(:,ifun),fxgrid(:,ifun))
		end if
		
		! Increasing Concave to Decreasing or inverted "U" f
		if (fmodel(ifun).eq.7) then
		  call GetPhi(xobs(:,ifun),xmin(ifun),xrange(ifun),CosFun2,CrossProdFun, &
		              ConstFun2,ConstCosFun,nbasis,nobs,phixobs(:,:,irest))
	      call GetPhi(xgrid(:,ifun),xmin(ifun),xrange(ifun),CosFun2,CrossProdFun, &
	                  ConstFun2,ConstCosFun,nbasis,nint+1,phixgrid(:,:,irest))
				
		  ! @ Initialize fxobs and fxgrid @
	      call GetUf(fpm(ifun),omega(ifun),psi(ifun),theta(:,ifun),xobs(:,ifun), &
	                 xgrid(:,ifun),phixobs(:,:,irest),phixgrid(:,:,irest),xdelta(ifun), &
	                 xinxgrid(:,ifun),xidelta(:,ifun),xrange(ifun),xmid(ifun),quadfacts, &
	                 intsimpfacts,nbasis,nr,nobs,nint+1,fxobs(:,ifun),fxgrid(:,ifun))
		end if
		
	    irest=irest+1
	  end if
	end do
	
	! Initialize zobs (auxiliary variables)
	do iobs=1,nobs
	   if(yobs(iobs).eq.1.d0) zobs(iobs)=ltnormrnd(sum(fxobs(iobs,:)),1.d0,0.d0)
	   if(yobs(iobs).eq.0.d0) zobs(iobs)=rtnormrnd(sum(fxobs(iobs,:)),1.d0,0.d0)
	end do

	! ************************************************************************************
	! * If some f have constraints, do initial MCMC to identify
	! * good parameter values of metm for metropolis
	! ************************************************************************************
	if(maxval(fmodel).gt.1) then ! Got constraints
	  call dblepr('Initializing MCMC parameters ...',-1,1.d0,0)
	  ! **********************************************************************
	  ! * Monitor pmet = Proportion of acceptance for MCMC
	  ! * If pmet is too large, then increase metm and mets by a factor of 10
	  ! * If pmet is too small, then reduce metm and mets by a factor of 10
	  ! **********************************************************************
	  do imodmet=1,maxmodmet  ! Allow up to maxmodmet adapations of pmet
	    imodmetg=imodmet
	    
	    pmet=0
	    iflag_AM=0  ! iflag_AM = 0 for static metroplis
	    theta_mAM=0.d0
	    theta_vAM=1.d0
	    icount_AM=0
	    ! Do small "burn-in" and initial change before adaptation of IG distribution
	    do imcmc=1,nblow0
	      call rchkusr()  ! user interrupt
	      call GetMCMC()
	    end do  ! End MCMC loop
	    
	    pmet=0  ! Counts Metropolis acceptances for theta if fmodel > 1
	    iflag_AM=1  ! iflag_AM = 1 for adaptive metropolis
	    theta_mAM=0.d0
	    theta_vAM=1.d0
	    icount_AM=0
	    ! Do round of adaptive MCMC for IG distribution to get proportion of acceptances
	    do imcmc=1,nblow0
	      call rchkusr()  ! user interrupt
	      call GetMCMC()
	    end do  ! End MCMC loop
	    
	    pok=0  ! Count of functions where pmet is ok
	    do ifun=1,nfun  ! Change metm depending on pmet
	      if (fmodel(ifun).gt.1) then
	        if (dble(pmet(ifun))/dble(nblow0).gt.0.6d0) then
	        ! Too many acceptances; small MCMC steps; increase metm
	          write(acceptr,fmt='(F6.4)') dble(pmet(ifun))/dble(nblow0)
	          write(func,fmt='(I1)') ifun
			  acceptmessage='function '//func//': pmet = '//acceptr//' > 0.6. Increase metm and redo MCMC loop'	
			  call dblepr(acceptmessage,-1,1.d0,0)
			  
			  metm(ifun)=metm(ifun)*10.d0
			  met_var_all(ifun)=metm(ifun)
			  met_mAM(ifun)=metm(ifun)
			  if (met_alpha.gt.2.d0) then
			    mets(ifun)=metm(ifun)/dsqrt(met_alpha-2.d0)
			    metv(ifun)=mets(ifun)**2.d0
			  else
			    mets(ifun)=1.d0
			    metv(ifun)=1.d0
			  end if
			  met_vAM(ifun)=metv(ifun)
			  
			  ! Metropolis for omega and psi
			  met_omega_m(ifun)=10.d0*met_omega_m(ifun)
			  met_omega_beta(ifun)=(met_omega_alpha-1.d0)*met_omega_m(ifun)
			  
			  met_psi_m(ifun)=10.d0*met_psi_m(ifun)
			  met_psi_beta(ifun)=(met_psi_alpha-1.d0)*met_psi_m(ifun)
			  
			  ! Reset parameters for function k
			  theta(1,ifun)=1.d0
			  do k=2,(nbasis+1)
			    theta(k,ifun)=0.01d0*tau(ifun)*dsqrt(gamvec(k-1,ifun))*rndnorm()
			  end do
			  tau2(ifun)=tau2_m0
			  tau(ifun)=dsqrt(tau2_m0)
			  gampar(ifun)=1.d0/w0
			  lngampar(ifun)=dlog(gampar(ifun))
			  gamvec(:,ifun)=dexp(-kall*gampar(ifun))
			  psi(ifun)=psifixed
			  omega(ifun)=(xmin(ifun)+xmax(ifun))/2.d0
			  zeta(ifun)=dlog(tau2(ifun))-kbar*gampar(ifun)
			  alpha(ifun)=0.d0
			else if (dble(pmet(ifun))/dble(nblow0).lt.0.3d0) then
			  ! Too few acceptance; big MCMC steps; decrease metm
			  write(acceptr,fmt='(F6.4)') dble(pmet(ifun))/dble(nblow0)
	          write(func,fmt='(I1)') ifun
			  acceptmessage='function '//func//': pmet = '//acceptr//' < 0.3. Reduce metm and redo MCMC loop'	
			  call dblepr(acceptmessage,-1,1.d0,0)
			  
			  metm(ifun)=metm(ifun)/10.d0
			  met_var_all(ifun)=metm(ifun)
			  met_mAM(ifun)=metm(ifun)
			  if (met_alpha.gt.2.d0) then
			    mets(ifun)=metm(ifun)/dsqrt(met_alpha-2.d0)
			    metv(ifun)=mets(ifun)**2.d0
			  else
			    mets(ifun)=1.d0
			    metv(ifun)=1.d0
			  end if
			  met_vAM(ifun)=metv(ifun)
			  
			  ! Metropolis for omega and psi
			  met_omega_m(ifun)=met_omega_m(ifun)/10.d0
			  met_omega_beta(ifun)=(met_omega_alpha-1.d0)*met_omega_m(ifun)
			  
			  met_psi_m(ifun)=met_psi_m(ifun)/10.d0
			  met_psi_beta(ifun)=(met_psi_alpha-1.d0)*met_psi_m(ifun)
			  
			  ! Reset parameters for function k
			  theta(1,ifun)=1.d0
			  do k=2,(nbasis+1)
			    theta(k,ifun)=0.01d0*tau(ifun)*dsqrt(gamvec(k-1,ifun))*rndnorm()
			  end do
			  tau2(ifun)=tau2_m0
			  tau(ifun)=dsqrt(tau2_m0)
			  gampar(ifun)=1.d0/w0
			  lngampar(ifun)=dlog(gampar(ifun))
			  gamvec(:,ifun)=dexp(-kall*gampar(ifun))
			  psi(ifun)=psifixed
			  omega(ifun)=(xmin(ifun)+xmax(ifun))/2.d0
			  zeta(ifun)=dlog(tau2(ifun))-kbar*gampar(ifun)
			  alpha(ifun)=0.d0
			else
			  pok=pok+1
			end if
	      end if
	    end do
	    
	    ! All pmets look good
	    if (pok.eq.nfunconstraint) then  
	      exit  ! Stop pre-MCMC adaption loop
	    end if
	    
	    ! Still working on adaption of metm
	    if (imodmet.lt.maxmodmet) then
	      ! re-initialize parameters
	      if (met_alpha.gt.2.d0) then
	        mets=metm/dsqrt(met_alpha-2.d0)
	        metv=mets**2.d0  ! variance of inverse gamma
	      else
	        mets=1.d0
	        metv=1.d0
	      end if
	      met_beta=(met_alpha-1.d0)*metm
	      met_beta_AM=met_beta  ! updated IG parameter for adaptive metropolis
	      met_var_all=metm      ! IG values
	      met_mAM=metm
	      met_vAM=metv
	      
	      ifree=1
	      irest=1
	      do ifun=1,nfun
	        if (fmodel(ifun).eq.1) then
	          call GetFreef(theta(2:(nbasis+1),ifun),phixobsfree(:,:,ifree),phixgridfree(:,:,ifree),&
	                        nbasis,nobs,nint+1,fxobs(:,ifun),fxgrid(:,ifun))
	
	          ifree=ifree+1
	        else 
	          if (fmodel(ifun).eq.2) then
	            call GetUpf(fpm(ifun),theta(:,ifun),phixobs(:,:,irest),phixgrid(:,:,irest), &
	                        quadfacts,nbasis,nr,nobs,nint+1,fxobs(:,ifun),fxgrid(:,ifun))
	          else if (fmodel(ifun).eq.3) then
	            call GetConvexf(fpm(ifun),alpha(ifun),theta(:,ifun),xobs(:,ifun),xgrid(:,ifun),&
	                            xmid(ifun),phixobs(:,:,irest),phixgrid(:,:,irest),quadfacts, &
					            nbasis,nr,nobs,nint+1,fxobs(:,ifun),fxgrid(:,ifun))
			  else if (fmodel(ifun).eq.4) then
			    call GetConcavef(fpm(ifun),alpha(ifun),theta(:,ifun),xobs(:,ifun),xgrid(:,ifun),&
	                             xmid(ifun),phixobs(:,:,irest),phixgrid(:,:,irest),quadfacts, &
					             nbasis,nr,nobs,nint+1,fxobs(:,ifun),fxgrid(:,ifun))
			  else if (fmodel(ifun).eq.5) then
			    call GetSf(fpm(ifun),omega(ifun),psi(ifun),alpha(ifun),theta(:,ifun),xobs(:,ifun),&
	                       xgrid(:,ifun),phixobs(:,:,irest),phixgrid(:,:,irest),xdelta(ifun),&
	                       xinxgrid(:,ifun),xidelta(:,ifun),xrange(ifun),xmid(ifun),quadfacts,&
	                       intsimpfacts,nbasis,nr,nobs,nint+1,fxobs(:,ifun),fxgrid(:,ifun))
	          else if (fmodel(ifun).eq.6) then
	            call GetRotateSf(fpm(ifun),omega(ifun),psi(ifun),alpha(ifun),theta(:,ifun), &
	                             xobs(:,ifun),xgrid(:,ifun),phixobs(:,:,irest),phixgrid(:,:,irest),&
	                             xdelta(ifun),xinxgrid(:,ifun),xidelta(:,ifun),xrange(ifun), &
	                             xmid(ifun),quadfacts,intsimpfacts,nbasis,nr,nobs,nint+1,&
	                             fxobs(:,ifun),fxgrid(:,ifun))
	          else 
	            call GetUf(fpm(ifun),omega(ifun),psi(ifun),theta(:,ifun),xobs(:,ifun), &
	                       xgrid(:,ifun),phixobs(:,:,irest),phixgrid(:,:,irest),xdelta(ifun), &
	                       xinxgrid(:,ifun),xidelta(:,ifun),xrange(ifun),xmid(ifun),quadfacts, &
	                       intsimpfacts,nbasis,nr,nobs,nint+1,fxobs(:,ifun),fxgrid(:,ifun))
	          end if
	          
	          irest=irest+1  
	        end if
	      end do
	    end if
	  end do  ! End loop over monitoring pmet
	end if  ! End of determining good mcmc parameters
	
	! ************************************************************************************
	! * Do Actual MCMC
	! ************************************************************************************
	isave=1
	nmcmc=nblow+nskip*smcmc  ! total number of MCMC
	do imcmc=1,nmcmc
		if(imcmc.eq.1) then
		  call dblepr('Burnin ...',-1,1.d0,0)
		  pmet=0		! Counts Metropolis acceptances for theta if fmodel > 1
		end if
		call rchkusr()  ! user interrupt
		call GetMCMC()  ! Draw samples
		
		if(imcmc.eq.nblow) then
		  do ifun=1,nfun
		    if(fmodel(ifun).gt.1) then
		      write(acceptr,fmt='(F6.4)') dble(pmet(ifun))/dble(nblow)
	          write(func,fmt='(I1)') ifun
		      acceptmessage='function '//func//': pmet = '//acceptr	
			  call dblepr(acceptmessage,-1,1.d0,0)
		    end if
		  end do
		  call dblepr('Main iterations ...',-1,1.d0,0)
		  pmet=0	
		end if
		
		! Store MCMC iterations
		if(imcmc.gt.nblow .and. mod(imcmc,nskip).eq.0) then 
            betag(isave,:)=beta
            zetag(isave,:)=zeta
            tau2g(isave,:)=tau2
            gammag(isave,:)=gampar
            thetag(:,:,isave)=theta
            alphag(isave,:)=alpha
            psig(isave,:)=psi
            omegag(isave,:)=omega
            
            ! fx includes the constant
            wbg(isave,:)=wb
            fxobsg(:,:,isave)=fxobs
            fxgridg(:,:,isave)=fxgrid
            yest=wb+sum(fxobs,2)
            do iobs=1,nobs
               yestg(isave,iobs)=binomrnd(1,cdfnorm(yest(iobs),0.d0,1.d0,1,0))
               phatg(isave,iobs)=cdfnorm(yest(iobs),0.d0,1.d0,1,0)
            end do
            ! Compute loglikelihood
            ehat=yobs-yest
            loglikeg(isave)=0.d0
            do iobs=1,nobs
               loglikeg(isave)=loglikeg(isave)+ &
                               yobs(iobs)*cdfnorm(yest(iobs),0.d0,1.d0,1,1)+ &
                               (1.d0-yobs(iobs))*cdfnorm(yest(iobs),0.d0,1.d0,0,1)
            end do
            
            ! Compute log prior
			logpriorg(isave)=GetLogPrior()
			
			if (mod(isave,ndisp).eq.0) then
				call cpu_time(itime)					! cpu_time
				call sprint(isave,smcmc,itime-stime)	! Print
			end if
			isave=isave+1
		end if
	end do
    deallocate(phixobsfree,phixobsfreet,phixobs,phixgridfree,phixgrid,phi2)! trig function
	pmetg=dble(pmet)/dble(smcmc*nskip)
	do ifun=1,nfun
	  if(fmodel(ifun).gt.1) then
	    write(acceptr,fmt='(F6.4)') pmetg(ifun)
	    write(func,fmt='(I1)') ifun
	    acceptmessage='function '//func//': pmet = '//acceptr
	    call dblepr(acceptmessage,-1,1.d0,0)
	  end if
	end do
	call rndend()
	
!=========================================================================================	
contains

!======= Full conditionals ===============================================================
subroutine GetMCMC()
  implicit none
  
  !Internal arguments
  real*8 :: gamrnd,rndunif,rtgamrnd,normrnd,ltnormrnd,rtnormrnd,tnormrnd,ltgamrnd  ! R ftn
  real*8 :: resid(nobs),vni(nparw,nparw),vn(nparw,nparw),bn(nparw) ! for beta
  real*8 :: rk(nobs),vni1(nbasis,nbasis),vn1(nbasis,nbasis),bn1(nbasis)  ! theta of freef
  real*8 :: met_var_all_new(nfun),met_beta_new(nfun),testp  ! theta of constrains f
  real*8 :: theta0_old,theta_old(nbasis),met_var0,met_var,met_std0
  real*8 :: theta0_new,theta02_new,theta_new(nbasis),theta2_new(nbasis),thetanew(nbasis+1)
  real*8 :: theta0_new_lnp,theta0_old_lnp,fxobs_new(nobs),fxgrid_new(nint+1)
  real*8 :: met_stdS,met_varS,psi_old,psi_lnpold,psi_new,psi_lnpnew
  real*8 :: omega_old,omega_lnpold,omega_new,omega_lnpnew
  real*8 :: resid_old(nobs),resid_new(nobs),sse_old,sse_new,sold,snew
  real*8 :: fx1(nobs),fxg1(nint+1),a_vni,a_vn,a_mn   ! alpha
  real*8 :: th2gam(nbasis),sn,bup  ! tau2
  real*8 :: ck1(nbasis),ckz(nbasis),u1(nbasis),u2,bnz(nbasis),bmin ! gamma
  real*8,allocatable :: bnz1(:)
  integer :: z(nbasis)
  
  do ifun=1,nfun
    ! Partial variance for theta (without)
    gamvec0(1,ifun)=1.d0
    gamvec0(2:(nbasis+1),ifun)=gamvec(:,ifun)
    thv0(:,ifun)=tau2(ifun)*gamvec(:,ifun) ! partial variance for theta1,...,thetaK
    thv00(1,ifun)=theta0_v0                ! partial variance for theta0,theta1,...,thetaK
    thv00(2:(nbasis+1),ifun)=thv0(:,ifun)
    theta0(ifun)=theta(1,ifun)
    theta02(ifun)=theta0(ifun)**2.d0                ! theta0^2
    theta2(:,ifun)=theta(2:(nbasis+1),ifun)**2.d0  ! Other theta's squared
  end do

  !---------------------------------------------------------------------------------------
  ! * Generate beta from scale invariant prior 
  ! *	beta does not include the intercept for all models
  !---------------------------------------------------------------------------------------
  resid=zobs-sum(fxobs,2)
  vni=wtw+beta_v0i
  call inverse(vni,nparw,vn)
  bn=matmul(vn,beta_v0im0+matmul(wdatat,resid))
  
  call mvnrnd(bn,vn,nparw,beta)
  wb=matmul(wdata,beta)
  
  !---------------------------------------------------------------------------------------
  ! * Generate theta
  ! * For S-shape, also generate omega and psi
  ! * from squish functions
  !---------------------------------------------------------------------------------------
  met_var_all_new=0.d0
  met_beta_new=0.d0
  resid=zobs-wb-sum(fxobs,2)  ! Resid takes off all of the f's
  ifree=1  ! ith : no shape 
  irest=1  ! ith : shape constraints
  do ifun=1,nfun
    testp=0.d0   ! Metropolis test p
    rk=resid+fxobs(:,ifun)  ! Add back fk to residuals
    if (fmodel(ifun).eq.1) then
      ! No shape constraint: Do Gibbs. 
      ! theta ~ N(0,tau2*exp(-k*gamma)
      vni1=phi2(:,:,ifree)
      do k=1,nbasis
        vni1(k,k)=vni1(k,k)+1.d0/thv0(k,ifun)
      end do
      call inverse(vni1,nbasis,vn1)
      bn1=matmul(vn1,matmul(phixobsfreet(:,:,ifree),rk)) ! Prior mean is zero for theta
      call mvnrnd(bn1,vn1,nbasis,theta(2:(nbasis+1),ifun))
      theta(1,ifun)=0.d0
      theta02(ifun)=0.d0
      theta2(:,ifun)=theta(2:(nbasis+1),ifun)**2.d0
      call GetFreef(theta(2:(nbasis+1),ifun),phixobsfree(:,:,ifree),phixgridfree(:,:,ifree),&
	                nbasis,nobs,nint+1,fxobs(:,ifun),fxgrid(:,ifun))
	                
	  ifree=ifree+1
	else
	  ! Got shape constraints:
	  ! theta ~ N(0,tau2*exp(-k*gamma)
	  ! @ Got shape constraints: theta ~ N(0,tau2*exp(-k*gamma) ) @
	  theta0_old=theta(1,ifun)
	  theta_old=theta(2:(nbasis+1),ifun)
	
	  ! @ Random walk Metropolis for theta @
	  ! @ get variances for t-dist random walk 	@
	  met_var_all_new(ifun)=met_beta_AM(ifun)/gamrnd(met_alpha,1.d0)
	  ck=metw*met_var_all_new(ifun) + (1.d0-metw)*met_mAM(ifun)
	  met_beta_new(ifun)=(met_alpha-1.d0)*ck
	
	  ck=5.66d0
	  met_var0=ck*met_var_all_new(ifun)
	  met_var=ck*met_var_all_new(ifun)

	  met_std0=dsqrt(met_var0)
	  ! @ Random walk from normal truncated below theta > 0 @
	  theta0_new=ltnormrnd(theta0_old,met_std0,0.d0)
	  do k=1,nbasis
	    ! @ Normal random walk @
	    theta_new(k)=normrnd(theta_old(k),dsqrt(met_var*thv0(k,ifun)))
	  end do

	  theta02_new=theta0_new**2.d0
	  theta2_new=theta_new**2.d0
	
	  ! @ Normalizing constant for generating theta_new @
	  theta0_new_lnp=cdfnorm(-theta0_old/met_std0,0.d0,1.d0,0,1)		
	  ! @ Normalizing constant for generating theta_old @
	  theta0_old_lnp=cdfnorm(-theta0_new/met_std0,0.d0,1.d0,0,1)	
	  
	  thetanew(1)=theta0_new
	  thetanew(2:(nbasis+1))=theta_new
	  if (fmodel(ifun).eq.2) then
	    call  GetUpf(fpm(ifun),thetanew,phixobs(:,:,irest),phixgrid(:,:,irest),&
	                 quadfacts,nbasis,nr,nobs,nint+1,fxobs_new,fxgrid_new)
	  else if (fmodel(ifun).eq.3) then
	    call GetConvexf(fpm(ifun),alpha(ifun),thetanew,xobs(:,ifun),xgrid(:,ifun),&
	                    xmid(ifun),phixobs(:,:,irest),phixgrid(:,:,irest),quadfacts, &
					    nbasis,nr,nobs,nint+1,fxobs_new,fxgrid_new)
	  else if (fmodel(ifun).eq.4) then
	    call GetConcavef(fpm(ifun),alpha(ifun),thetanew,xobs(:,ifun),xgrid(:,ifun),&
	                     xmid(ifun),phixobs(:,:,irest),phixgrid(:,:,irest),quadfacts, &
					     nbasis,nr,nobs,nint+1,fxobs_new,fxgrid_new)
	  else
	    ! ********************************************************
	    ! * Generate omega and psi for squish function
	    ! ********************************************************
	    
	    ! ***********************************************
	    ! * Generate psi > 0 
	    ! ***********************************************
	    if (iflagpsi.eq.1) then
	      !! Begin change
	      met_varS=met_psi_beta(ifun)/gamrnd(met_psi_alpha,1.d0)
	      met_stdS=dsqrt(met_varS)
	      !! End change
		  psi_old=psi(ifun)
		  ! @ Generate from truncated normal > 0
		  psi_new=ltnormrnd(psi_old,met_stdS,0.d0)	
		  ! @ LogNormalizing constant for g(new|old)
		  psi_lnpnew=cdfnorm(-psi_old/met_stdS,0.d0,1.d0,0,1) 
		  ! @ Log constant for g(old|new) @
		  psi_lnpold=cdfnorm(-psi_new/met_stdS,0.d0,1.d0,0,1) 
		  
		  testp=testp-&
		        ((psi_new-psi_m0)**2.d0)/(2.d0*psi_v0)+&  ! prior for candidate psi_new     
			    ((psi_old-psi_m0)**2.d0)/(2.d0*psi_v0)-&  ! prior for old psi                 
			    psi_lnpold+ &							  ! constant for g(old|new)			
			    psi_lnpnew							      ! constant for g(new|old)
	    else
		  psi_old=psi(ifun)  ! Placeholder
		  psi_new=psi(ifun)  ! Placeholder
	    end if
	    
	    ! *****************************************************
	    ! * Generate xmin < omega < xmax
	    ! * It seems that omega can get stuck at endpoint
	    ! * Randomly generate it from prior distribution
	    ! *****************************************************
	    !! Begin change
	    met_varS=met_omega_beta(ifun)/gamrnd(met_omega_alpha,1.d0)
	    met_stdS=dsqrt(met_varS)
	    !! End change
	    omega_old=omega(ifun)
	    ! @ Generate omega from trunctated normal
	    omega_new=tnormrnd(omega_old,met_stdS,xmin(ifun),xmax(ifun))	
	    ! @ Normalizing constants for omega
	    ! @ Normalizing constant for g(new|old) @
	    omega_lnpnew=dlog(cdfnorm((xmax(ifun)-omega_old)/met_stdS,0.d0,1.d0,1,0)- &
		             cdfnorm((xmin(ifun)-omega_old)/met_stdS,0.d0,1.d0,1,0))
	    ! @ Normalizing constant for g(old|new) @
	    omega_lnpold=dlog(cdfnorm((xmax(ifun)-omega_new)/met_stdS,0.d0,1.d0,1,0)- &
		             cdfnorm((xmin(ifun)-omega_new)/met_stdS,0.d0,1.d0,1,0))

	    testp=testp- &  
		      ((omega_new-omega_m0)**2.d0)/(2.d0*omega_v0)+& ! prior for candidate omega_new
		      ((omega_old-omega_m0)**2.d0)/(2.d0*omega_v0)-& ! prior for old omega  
		      omega_lnpold+	&							     ! constant for g(old|new)			
		      omega_lnpnew								     ! constant for g(new|old) 
		      
		if (fmodel(ifun).eq.5) then
		  call GetSf(fpm(ifun),omega_new,psi_new,alpha(ifun),thetanew,xobs(:,ifun),&
	                 xgrid(:,ifun),phixobs(:,:,irest),phixgrid(:,:,irest),xdelta(ifun),&
	                 xinxgrid(:,ifun),xidelta(:,ifun),xrange(ifun),xmid(ifun),quadfacts,&
	                 intsimpfacts,nbasis,nr,nobs,nint+1,fxobs_new,fxgrid_new)
	    else if (fmodel(ifun).eq.6) then
	      call GetRotateSf(fpm(ifun),omega_new,psi_new,alpha(ifun),thetanew,xobs(:,ifun),&
	                       xgrid(:,ifun),phixobs(:,:,irest),phixgrid(:,:,irest),&
	                       xdelta(ifun),xinxgrid(:,ifun),xidelta(:,ifun),xrange(ifun), &
	                       xmid(ifun),quadfacts,intsimpfacts,nbasis,nr,nobs,nint+1,&
	                       fxobs_new,fxgrid_new)
	    else 
	      call GetUf(fpm(ifun),omega_new,psi_new,thetanew,xobs(:,ifun), &
	                 xgrid(:,ifun),phixobs(:,:,irest),phixgrid(:,:,irest),xdelta(ifun), &
	                 xinxgrid(:,ifun),xidelta(:,ifun),xrange(ifun),xmid(ifun),quadfacts, &
	                 intsimpfacts,nbasis,nr,nobs,nint+1,fxobs_new,fxgrid_new)
	    end if
	  end if
	  
	  ! Metropolis Test
	  resid_new=rk-fxobs_new
	  sse_new=sum(resid_new**2.d0)
	  
	  resid_old=rk-fxobs(:,ifun)
	  sse_old=sum(resid_old**2.d0)
	  
	  snew=sum(theta2_new/thv0(:,ifun))
	  sold=sum(theta2(:,ifun)/thv0(:,ifun))
	  
	  testp=testp- &
	        (sse_new-sse_old)/(2.d0)-	&	 ! Likelihood 				
	        (theta02_new)/(2.d0*theta0_v0)+ &  ! Prior for new theta0 		
	        (theta02(ifun))/(2.d0*theta0_v0)-	&	 ! Prior for old theta0 		
	        (snew-sold)/(2.d0)-	&			 ! Prior for new & old theta	
	        theta0_old_lnp+theta0_new_lnp- &		 ! Constants for truncated r-w for theta0 
	        ((met_alpha+1.d0)*dlog(met_var_all(ifun)))-&
	        (met_beta_new(ifun)/met_var_all(ifun))+ &
	        ((met_alpha+1.d0)*dlog(met_var_all_new(ifun)))+&
	        (met_beta_AM(ifun)/met_var_all_new(ifun))
	        
	  if(dlog(rndunif()).le.testp) then
	    ! Accept candidate
	    theta0(ifun)=theta0_new
		theta(1,ifun)=theta0_new
		theta(2:(nbasis+1),ifun)=theta_new
		theta2(:,ifun)=theta2_new
		theta02(ifun)=theta02_new
		if (fmodel(ifun).eq.5 .or. fmodel(ifun).eq.6 .or. fmodel(ifun).eq.7) then
		  psi(ifun)=psi_new
		  omega(ifun)=omega_new
		end if
		met_var_all(ifun)=met_var_all_new(ifun)
		fxobs(:,ifun)=fxobs_new
		fxgrid(:,ifun)=fxgrid_new
		pmet(ifun)=pmet(ifun)+1
	  end if
	  
	  if (iflag_AM.eq.1) then 
		! @ Update mean, variance, alpha, and beta for metropolis @
		icount_AM=icount_AM+1
		theta_mAM(:,ifun)=theta_mAM(:,ifun)+(theta(:,ifun)-theta_mAM(:,ifun))/dble(icount_AM)
		theta_vAM(:,ifun)=(dble(icount_AM-1)/dble(icount_AM))*theta_vAM(:,ifun)+ &
				          ((theta(:,ifun)-theta_mAM(:,ifun))**2.d0)/dble(icount_AM)
		met_mAM(ifun)=met_mAM(ifun)+(met_var_all(ifun)-met_mAM(ifun))/dble(icount_AM)
		met_vAM(ifun)=(dble(icount_AM-1)/dble(icount_AM))*met_vAM(ifun) + &
				((met_var_all(ifun)-met_mAM(ifun))**2.d0)/dble(icount_AM)
	  end if
	  ck=metw*met_var_all(ifun)+(1.d0-metw)*met_mAM(ifun)
	  met_beta_AM(ifun)=(met_alpha-1.d0)*ck
	  
	  irest=irest+1
	end if  
	resid=zobs-wb-sum(fxobs,2)  ! Resid takes off all of the f's
  end do
  
  !---------------------------------------------------------------------------------------
  ! *   Generate alpha
  ! *	Constrained convex, concave, S models have linear term
  ! *	Generate alpha from truncated normal and alpha
  !---------------------------------------------------------------------------------------
  do ifun=1,nfun
    resid=zobs-wb-sum(fxobs,2)
    if(fmodel(ifun).gt.2 .and. fmodel(ifun).lt.7) then
      ! @ Take off old alpha*(x-xmid) from fx @
	  fx1=fxobs(:,ifun)-alpha(ifun)*(xobs(:,ifun)-xmid(ifun))
	  fxg1=fxgrid(:,ifun)-alpha(ifun)*(xgrid(:,ifun)-xmid(ifun))

	  rk=resid+alpha(ifun)*(xobs(:,ifun)-xmid(ifun))  ! Add back alpha*(x-xmid)
	
	  ! @ Generate alpha ~ Truncated Normal @
	  a_vni=xtx(ifun)+alpha_v0i
	  a_vn=1.d0/a_vni
	  a_mn=a_vn*(dot_product((xobs(:,ifun)-xmid(ifun)),rk)+alpha_v0im0)
	  if (fpm(ifun).eq.1.d0) then
	    ! @ Generate Normal truncated below at alpha > 0 @
	    alpha(ifun)=ltnormrnd(a_mn,dsqrt(a_vn),0.d0)	
	  else
	    ! @ Generate Normal truncated below at alpha < 0 @
	    alpha(ifun)=rtnormrnd(a_mn,dsqrt(a_vn),0.d0)
	  end if
	  fxobs(:,ifun)=fx1+alpha(ifun)*(xobs(:,ifun)-xmid(ifun))
	  fxgrid(:,ifun)=fxg1+alpha(ifun)*(xgrid(:,ifun)-xmid(ifun))
	else
	  alpha(ifun)=0.d0
	end if
  end do
  
  !---------------------------------------------------------------------------------------
  ! *   Generate tau2 and gamma
  !---------------------------------------------------------------------------------------
  do ifun=1,nfun
    !-------------------------------------------------------------------------------------
    ! * Generate tau2
    ! *	If prior tau2 ~ IG, then generate from IG
    ! *	If prior tau2 ~ Exponential, 
    ! *		then use IG as proposal distribution, and apply Metropolis
    ! *	rn   = r0 + nbasis
    ! *	sn   = s0 + sum_k theta_k^2*exp(k*gamma)
    !-------------------------------------------------------------------------------------
    th2gam=theta2(:,ifun)/gamvec(:,ifun)
	sn=sum(th2gam)
	if (iflagprior.eq.0) then 
		! @ IG Prior is conjugate @
		tau2(ifun)=(tau2_s0+sn)/(2.d0*gamrnd(tau2_rn/2.d0,1.d0))
		tau(ifun)=dsqrt(tau2(ifun))
	else 	
		! @ Lasso prior: use slice sampling @
		bup=tau2(ifun)-dlog(rndunif())/tau2_u0
		
		! @ Gamma distribution truncated below @
		tau2i(ifun)=ltgamrnd(dble(nbasis-2)/2.d0,2.d0/sn,1.d0/bup);  
		tau2(ifun)=1.d0/tau2i(ifun)
		tau(ifun)=dsqrt(tau2(ifun))
	end if
	thv0(:,ifun)=tau2(ifun)*gamvec(:,ifun)
	
	!-------------------------------------------------------------------------------------
    ! * Generate gamma with slice sampling
    !-------------------------------------------------------------------------------------
    ! @ theta0 does not depend on gamma @
    ck1=theta2(:,ifun)/(2.d0*tau2(ifun))
	
	! @ Worry about ck1 = 0 when thetak = 0. Doesn't appear in likelihood for gamma @
	if (count(ck1.eq.0.d0).eq.nbasis) then
		gampar(ifun) = 1.d0/w0
	else 
		z=transfer(ck1.eq.0.d0, (/ (1,k=1,nbasis) /))
		ckz=dble(z)*(1.d0-ck1)+ck1 ! @ ckz = 1 when theta_k = 0 @
		u1=(/ (rndunif(),k=1,nbasis) /)
		bnz=gampar(ifun)+(dlog(ckz-dlog(u1)*gamvec(:,ifun))-dlog(ckz))/kall
		if (count(z.eq.1).gt.0) then
			allocate(bnz1(count(z.ne.1)))
			iloop=1
			do k=1,nbasis
				if (z(k).ne.1) then
					bnz1(iloop)=bnz(k)
					iloop=iloop+1
				end if
			end do
			bmin=minval(bnz1)
			deallocate(bnz1)
		else
			bmin=minval(bnz)
		end if
		u2=rndunif()
		gampar(ifun)=bmin+dlog(u2+(1.d0-u2)*dexp(-wk*bmin))/wk
	end if
	
	gamvec(:,ifun)=dexp(-gampar(ifun)*kall)
	lngampar(ifun)=dlog(gampar(ifun))
	zeta(ifun)=dlog(tau2(ifun))-kbar*gampar(ifun)
  end do
  
  !---------------------------------------------------------------------------------------
  ! *   Generate auxiliary variables
  !---------------------------------------------------------------------------------------
  do iobs=1,nobs
     if(yobs(iobs).eq.1.d0) zobs(iobs)=ltnormrnd(wb(iobs)+sum(fxobs(iobs,:)),1.d0,0.d0)
     if(yobs(iobs).eq.0.d0) zobs(iobs)=rtnormrnd(wb(iobs)+sum(fxobs(iobs,:)),1.d0,0.d0)
  end do
  
  return
end subroutine GetMCMC


!-----------------------------------------------------------------------------------------
! * GetLogPrior
! *	Compute log prior density
!-----------------------------------------------------------------------------------------
function GetLogPrior()
  implicit none	
  
  !Output argument
  real*8 :: GetLogPrior
	
  !Internal argument
  real*8 :: residb(nparw),thvar(nbasis+1),thetasq(nbasis+1,nfun)
  real*8 :: sse,theta0_lnp0,alpha_lnp0,lnpriorf
  real*8 :: cdfnorm
	
  lnpriorf=0.d0
  thetasq=theta**2.d0
  
  ! @ Normal prior for beta ~ N(bn,Bn) @
  residb=beta-beta_m0
  lnpriorf=lnpriorf-dot_product(residb,matmul(beta_v0i,residb))/(2.d0)- &
           dble(nparw)*dlog(2.d0*PI)/2.d0-beta_lnv0/2.d0

  do ifun=1,nfun
    ! @ Normal prior of theta   	
    if (fmodel(ifun).eq.1) then
      thvar(2:(nbasis+1))=tau2(ifun)*gamvec(:,ifun)
      sse=sum(thetasq(2:(nbasis+1),ifun)/thvar(2:(nbasis+1)))
      lnpriorf=lnpriorf-sse/2.d0-dble(nbasis)*dlog(2.d0*PI)/2.d0- &
               sum(dlog(thvar(2:(nbasis+1))))/2.d0
    else
      ! @ Truncated Normal Prior for theta0 ~ N(0,tau2)*I(theta0 > 0) @
      thvar(1)=theta0_v0
      thvar(2:(nbasis+1))=tau2(ifun)*gamvec(:,ifun)
      sse=sum(thetasq(:,ifun)/thvar)
      theta0_lnp0=cdfnorm(-theta0_m0/(theta0_s0),0.d0,1.d0,0,1)
      lnpriorf=lnpriorf-sse/2.d0-dble(nbasis+1)*dlog(2.d0*PI)/2.d0- &
               sum(dlog(thvar))/2.d0-theta0_lnp0
    end if

    ! @ Exponential prior for gamma ~ EXP(w0) @
    lnpriorf=lnpriorf-w0*gampar(ifun)-dlog(w0)

    ! @ Prior for tau2 @
    if (iflagprior.eq.1) then
        ! @ tau2 ~ Exp(lambda) @
        lnpriorf=lnpriorf-tau2_u0*tau2(ifun)-tau2_u0
    else 
        ! @ tau2 is IG @
        lnpriorf=lnpriorf+LogfIG(tau2(ifun),tau2_r0,tau2_s0)	
    end if
	
	if (fmodel(ifun).gt.2 .and. fmodel(ifun).lt.7) then
      ! @ Truncated normal prior for alpha~N(m0,v0)I(alpha>0) @
      alpha_lnp0=cdfnorm(-alpha_m0/alpha_s0,0.d0,1.d0,0,1)
      lnpriorf=lnpriorf-((alpha(ifun)-alpha_m0)**2.d0)/(2.d0*alpha_v0)- &
               dlog(2.d0*PI*alpha_v0)/2.d0-alpha_lnp0
	end if
	
	! Truncated normal priors for psi and omega
    if (fmodel(ifun).eq.5 .or. fmodel(ifun).eq.6 .or. fmodel(ifun).eq.7) then
      if (iflagpsi.eq.1) then
        lnpriorf=lnpriorf-((psi(ifun)-psi_m0)**2.d0)/(2.d0*psi_v0)- &
                 dlog(2.d0*PI*psi_v0)/2.d0-psi_lnp0
      end if
      lnpriorf=lnpriorf-((omega(ifun)-omega_m0)**2.d0)/(2.d0*omega_v0)- &
               dlog(2.d0*PI*omega_v0)/2.d0-omega_lnp0
    end if
  end do
  
  GetLogPrior=lnpriorf

  return
end function GetLogPrior


end subroutine bsarprobit