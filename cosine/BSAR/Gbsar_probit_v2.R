'Gbsar_probit' = function(yobs,wdata,xobs,nobs,nparw,nfun,nbasis,nint,fmodel,fpm,privals,mcvals,
                          xmin,xmax,xmid,xrange,xdelta,xgrid,xinxgrid,xidelta)
  
{
 #######################
 # Prior distributions #
 #######################
      
 # theta0 ~ N(0,v0^2) for free, theta0 ~ N(0,v0^2)I(theta0 > 0)
 theta0_m0=privals$theta0_m0
 theta0_s0=privals$theta0_s0
 theta0_v0=theta0_s0^2
 theta_m0=privals$theta_m0
 
 # IG, Exponential prior for tau
 iflagprior=privals$iflagprior
 tau2_m0=privals$tau2_m0
 tau2_v0=privals$tau2_v0
 tau2_r0=2*(2+tau2_m0^2/tau2_v0) # IG prior for tau
 tau2_s0=tau2_m0*(tau2_r0-2)
 tau2_u0=1/tau2_m0	# Exponential prior for tau
  
 # gamma ~ Exp(w0)
 w0=privals$w0
  
 # beta ~ N(bm0,bv0)
 beta_m0=privals$beta_m0
 beta_v0=privals$beta_v0
 beta_v0i=solve(beta_v0)
  
 # alpha ~ Truncated normal(m0, v0)
 alpha_m0=privals$alpha_m0
 alpha_s0=privals$alpha_s0

 # psi ~ Truncated normal is slope for squish (tanh) function in S shaped
 iflagpsi=privals$iflagpsi
 psifixed=privals$psifixed
 psi_m0=psifixed
 psi_s0=10*psifixed

 # omega ~ Truncated normal is point of inflection for "S" shaped f 
 omega_m0=privals$omega_m0
 omega_s0=privals$omega_s0
 omega_v0=omega_s0^2
		
 # MCMC
 nblow0=mcvals$nblow0
 nblow=mcvals$nblow
 smcmc=mcvals$smcmc
 nskip=mcvals$nskip
 ndisp=mcvals$ndisp
 maxmodmet=mcvals$maxmodmet
 if(max(fmodel)==1) maxmodmet=0
 nmcmc=nblow+nskip*smcmc  # total number of MCMC
  
 # Basis
 nr=(nbasis+1)*(nbasis+2)/2
 kall=seq(1,by=1,length=nbasis)
 kall0=c(0,kall)
 kbar=mean(kall0)
 kvec=w0/(w0+kall)
 wk=sum(kall)/2-w0

 # Vector used in simpson's integration
 intsimpfacts=numeric(nint+1)
 intsimpfacts[1]=1
 intsimpfacts[nint+1]=1
 intsimpfacts[seq(2,nint,2)]=4
 intsimpfacts[seq(3,nint-1,2)]=2
 
 #Index
 leftindex=vech(kronecker(matrix(1,nr=1,nc=nbasis+1),seq(1,by=1,length=nbasis+1)))
 rightindex=vech(kronecker(matrix(1,nr=nbasis+1,nc=1),t(seq(1,by=1,length=nbasis+1))))
 multfact=vech(2*matrix(1,nr=nbasis+1,nc=nbasis+1)-diag(nbasis+1))
 quadfacts=cbind(multfact,leftindex,rightindex)
 
 #########################
 # Initialize parameters #
 #########################
  
 tau2=rep(tau2_m0,nfun)
 psi=rep(psifixed,nfun)
 omega=xmid
 gampar=rep(1/w0,nfun)
 gamvec=sapply(1:nfun,function(ifun)exp(-gampar[ifun]*kall))
 zeta=sapply(1:nfun,function(ifun)log(tau2[ifun])-kbar*gampar[ifun])
 theta=sapply(1:nfun, function(ifun)0.1*sqrt(gamvec[,ifun])*rnorm(nbasis))
 theta=rbind(0.1,theta)
 alpha=rep(0,nfun)
  
 beta=numeric(nparw)
 wb=wdata%*%beta
 
 fxobs=matrix(0,nrow=nobs,ncol=nfun)
 fxgrid=matrix(0,nrow=nint+1,ncol=nfun)
 
 # Initialize uobs (auxiliary variables)
 uobs=yobs
 
 ##############################
 # Matrics for saving results #
 ##############################
  
 # Parameters
 betag=matrix(0,smcmc,nparw)
 zetag=tau2g=gammag=alphag=psig=omegag=matrix(0,smcmc,nfun)
 thetag=array(0,dim=c(nbasis+1,nfun,smcmc))
 
 # fx includes the constant
 wbg=matrix(0,smcmc,nobs)
 fxobsg=array(0,dim=c(nobs,nfun,smcmc))
 fxgridg=array(0,dim=c(nint+1,nfun,smcmc))
 muhatg= phatg=matrix(0,nr=smcmc,nc=nobs)
 
 # Compute loglikelihood
 loglikeg=logpriorg=numeric(smcmc)
 invlikeg=matrix(0,smcmc,nobs)
 
 ###############################################
 # Adaptive Metropolis parameters (fmodel > 1) #
 ###############################################
 
 metm=rep(0.01,nfun)         # Mean of inverse gamma
 met_alpha=3

 met_beta_AM=(met_alpha-1)*metm
 met_var_all=metm      # IG values
 met_mAM=metm
 pmet=rep(0,nfun)
 icount_AM=0
 
 met_m=rep(0.0001,nfun)
 met_beta=(met_alpha-1)*met_m

 ##############################
 # Precomputes trig functions #
 ##############################
  
 nfree=sum(fmodel==1)          # number of functions that have unconstraints
 nfunconstraint=sum(fmodel>1)  # number of functions that have constraints
  
 trig=pretrig(nfree,nfunconstraint,nobs,nbasis,nint,nr,nfun,fmodel,fpm,
              xobs,xgrid,kall,xmin,xmax,xmid,xrange)
 phixobsfree=trig$phixobsfree
 phixobs=trig$phixobs
 phixgridfree=trig$phixgridfree
 phixgrid=trig$phixgrid
  
 foo=prefoo(nobs,nbasis,nint,nr,nfun,fmodel,fpm,xobs,xgrid,xmid,xrange,
            phixobsfree,phixobs,phixgridfree,phixgrid,
            fxobs,fxgrid,quadfacts,theta,alpha,psi,omega,
            xinxgrid,xidelta,intsimpfacts,xdelta)
 fxobs=foo$fxobs
 fxgrid=foo$fxgrid
  
 #########################################################################
 # Initial MCMC to identify good parameter values of metm for metropolis #
 #########################################################################
  
 iflag_AM=0  # iflag_AM = 0 for static metroplis

 #input for MCMC
 param=list(beta=beta,zeta=zeta,tau2=tau2,gampar=gampar,
            theta=theta,alpha=alpha,psi=psi,omega=omega,
            wb=wb,fxobs=fxobs,fxgrid=fxgrid,uobs=uobs)
 met=list(met_var_all=met_var_all,met_mAM=met_mAM,
          met_beta_AM=met_beta_AM,pmet=pmet,icount_AM=icount_AM)
 
 if(max(fmodel) > 1) { # Got constraints
   cat("Initializing MCMC parameters ...", "\n")

    # **********************************************************************
    # * Monitor pmet = Proportion of acceptance for MCMC
    # * If pmet is too large, then increase metm and mets by a factor of 10
    # * If pmet is too small, then reduce metm and mets by a factor of 10
    # **********************************************************************    
    for (imodmet in 1:maxmodmet){  # Allow up to maxmodmet adapations of pmet

      # Do small "burn-in" and initial change before adaptation of IG distribution
      for (imcmc in 1:nblow0){
        update.param=GetMCMC.probit(yobs,xobs,xmin,xmax,xmid,xgrid,fpm,fmodel,
                                    nparw,nfun,nobs,nint,nbasis,kall,kbar,wdata,
                                    quadfacts,intsimpfacts,phixobsfree,phixgridfree,
                                    phixobs,phixgrid,iflagpsi,wk,xdelta,xinxgrid,xidelta,
                                    xrange,iflag_AM,met_alpha,beta_m0,beta_v0i,
                                    theta0_v0,tau2_s0,tau2_u0,tau2_r0,imcmc,
                                    met_beta,psi_m0,psi_s0,omega_m0,omega_v0,
                                    alpha_m0,alpha_s0,iflagprior,param,met)
        param=update.param$param.out
        met=update.param$met.out
      }  # End MCMC loop

      iflag_AM=1  # iflag_AM = 1 for adaptive metropolis
      met$pmet=rep(0,nfun)  # Counts Metropolis acceptances for theta if fmodel > 1
      met$icount_AM=0
      
      # Do round of adaptive MCMC for IG distribution to get proportion of acceptances
      for (imcmc in 1:nblow0){
        update.param=GetMCMC.probit(yobs,xobs,xmin,xmax,xmid,xgrid,fpm,fmodel,
                                    nparw,nfun,nobs,nint,nbasis,kall,kbar,wdata,
                                    quadfacts,intsimpfacts,phixobsfree,phixgridfree,
                                    phixobs,phixgrid,iflagpsi,wk,xdelta,xinxgrid,xidelta,
                                    xrange,iflag_AM,met_alpha,beta_m0,beta_v0i,
                                    theta0_v0,tau2_s0,tau2_u0,tau2_r0,imcmc,
                                    met_beta,psi_m0,psi_s0,omega_m0,omega_v0,
                                    alpha_m0,alpha_s0,iflagprior,param,met)
        param=update.param$param.out
        met=update.param$met.out
      }  # End MCMC loop
      
      pok=0  # Number of functions where pmet is ok

      for (ifun in 1:nfun){  # Change metm depending on pmet
        if (fmodel[ifun] > 1) {
          if ((met$pmet[ifun]/nblow0) > 0.6) {
            cat(paste( "function ", ifun,": pmet = ",met$pmet[ifun]/nblow0," > 0.6. Increase metm and redo MCMC loop",'\n',sep=''));
            metm[ifun]=metm[ifun]*10
            
            # Metropolis for omega and psi
            met_m[ifun]=10*met_m[ifun]
            met_beta[ifun]=(met_alpha-1)*met_m[ifun]
            
          } else if (met$pmet[ifun]/nblow0 < 0.3) {
            cat(paste( "function ", ifun,": pmet = ", met$pmet[ifun]/nblow0," < 0.3. Reduce metm and redo MCMC loop",'\n', sep=''));
            metm[ifun]=metm[ifun]/10

            # Metropolis for omega and psi
            met_m[ifun]=met_m[ifun]/10
            met_beta[ifun]=(met_alpha-1)*met_m[ifun]
            
          } else {
            pok=pok+1
          }
        } 
      }
      
      if (pok == nfunconstraint) break # All pmets look good
      if (imodmet < maxmodmet) {  # Still working on adaption of metm

        # re-initialize parameters
        param=list(beta=beta,zeta=zeta,tau2=tau2,gampar=gampar,
                   theta=theta,alpha=alpha,psi=psi,omega=omega,
                   wb=wb,fxobs=fxobs,fxgrid=fxgrid,uobs=uobs)
        
        met$met_beta_AM=(met_alpha-1)*metm
        met$met_var_all=metm      # IG values
        met$met_mAM=metm
        met$pmet=rep(0,nfun)

      }
    }	
  }
  
 ##################
 # Do Actual MCMC #
 ##################
 
 isave=1
 met$pmet=rep(1,nfun)
  
 for (imcmc in 1:nmcmc){
   if(imcmc==1) {
     cat("Burn in ...", "\n")		# Counts Metropolis acceptances for theta if fmodel > 1
   }
    update.param=GetMCMC.probit(yobs,xobs,xmin,xmax,xmid,xgrid,fpm,fmodel,
                                nparw,nfun,nobs,nint,nbasis,kall,kbar,wdata,
                                quadfacts,intsimpfacts,phixobsfree,phixgridfree,
                                phixobs,phixgrid,iflagpsi,wk,xdelta,xinxgrid,xidelta,
                                xrange,iflag_AM,met_alpha,beta_m0,beta_v0i,
                                theta0_v0,tau2_s0,tau2_u0,tau2_r0,imcmc,
                                met_beta,psi_m0,psi_s0,omega_m0,omega_v0,
                                alpha_m0,alpha_s0,iflagprior,param,met)
    param=update.param$param.out
    met=update.param$met.out

    if(imcmc==nblow) {
      for (ifun in 1:nfun){
        if (fmodel[ifun]>1){
          cat(paste( "function ", ifun,": pmet = ", met$pmet[ifun]/nblow,'\n', sep=''))
        }
      }
      cat("Main iterations ...", "\n")		# Counts Metropolis acceptances for theta if fmodel > 1
      met$pmet=rep(1,nfun)
    }
    
    # Store MCMC iterations
    if((imcmc > nblow) && (imcmc%%nskip==0)) { 
      betag[isave,]=param$beta
      zetag[isave,]=param$zeta
      tau2g[isave,]=param$tau2
      gammag[isave,]=param$gampar
      thetag[,,isave]=param$theta
      alphag[isave,]=param$alpha
      psig[isave,]=param$psi
      omegag[isave,]=param$omega
      
      # fx includes the constant
      wbg[isave,]=param$wb
      fxobsg[,,isave]=param$fxobs
      fxgridg[,,isave]=param$fxgrid
      muhat=param$wb+rowSums(param$fxobs)
      muhatg[isave,]=muhat
      phat=pnorm(muhat,0,1)
      phatg[isave,]=phat
      
      # Compute loglikelihood
      loglikeg[isave]=sum(yobs*log(phat)+(1-yobs)*log(1-phat))
      
      # Compute log prior
      logpriorg[isave]=GetLogPrior(fmodel,nparw,nfun,nbasis,kall,w0,xmin,xmax,
                                   beta,beta_m0,beta_v0,beta_v0i,
                                   theta,theta0_m0,theta0_s0,theta0_v0,
                                   gampar,tau2,tau2_r0,tau2_s0,tau2_u0,
                                   alpha,alpha_m0,alpha_s0,
                                   iflagprior,iflagpsi,psi,psi_m0,psi_s0,
                                   omega,omega_m0,omega_s0)
      
      invlikeg[isave,]=1/dbinom(yobs,1,pnorm(muhat))
      
      isave=isave+1
    }
    
    if((imcmc%%ndisp==0)) { 
      cat('iterations',imcmc,'\n')
    }
  }
  
  pmetg=met$pmet/(smcmc*nskip)
  for (ifun in 1:nfun) {
    if (fmodel[ifun]>1){
      cat(paste( "function ", ifun,": pmet = ", pmetg[ifun],'\n', sep=''))
    }
  }
  return(list(betag=betag,zetag=zetag,tau2g=tau2g,gammag=gammag,thetag=thetag,
              alphag=alphag,psig=psig,omegag=omegag,wbg=wbg,loglikeg=loglikeg,
              logpriorg=logpriorg,invlikeg=invlikeg,fxobsg=fxobsg,fxgridg=fxgridg,
              muhatg=muhatg,pmetg=pmetg,phatg=phatg))
  
}	


#---------------------------Gbsar.probit.DataAug--------------------------------
'GetMCMC.probit' = function(yobs,xobs,xmin,xmax,xmid,xgrid,fpm,fmodel,
                            nparw,nfun,nobs,nint,nbasis,kall,kbar,wdata,
                            quadfacts,intsimpfacts,phixobsfree,phixgridfree,
                            phixobs,phixgrid,iflagpsi,wk,xdelta,xinxgrid,xidelta,
                            xrange,iflag_AM,met_alpha,beta_m0,beta_v0i,
                            theta0_v0,tau2_s0,tau2_u0,tau2_r0,imcmc,
                            met_beta,psi_m0,psi_s0,omega_m0,omega_v0,
                            alpha_m0,alpha_s0,iflagprior,param,met)
  
  
{
  #extract parameters
  beta=param$beta
  zeta=param$zeta
  tau2=param$tau2
  gampar=param$gampar
  theta=param$theta
  alpha=param$alpha
  psi=param$psi
  omega=param$omega
  wb=param$wb
  fxobs=param$fxobs
  fxgrid=param$fxgrid
  uobs=param$uobs
  
  gamvec=sapply(1:nfun,function(ifun)exp(-gampar[ifun]*kall))
  thv0=sapply(1:nfun,function(ifun) tau2[ifun]*gamvec[,ifun]) # partial variance for theta1,...,thetaK
  
  ###########
  # Do MCMC #
  ###########
  #Update beta
  resid=uobs-rowSums(fxobs)
  vni=t(wdata)%*%wdata+beta_v0i
  vn=solve(vni+diag(10e-7,nparw))
  U=chol(vn)
  bn=vn%*%(beta_v0i%*%beta_m0+t(wdata)%*%resid)
  beta=bn+t(U)*as.matrix(rnorm(nparw))
  wb=wdata%*%beta
  
  #Update theta
  gtheta=gtheta.probit(yobs,nfun,wb,fxobs,fxgrid,fmodel,nbasis,nobs,nint,nr,met_alpha,
                       quadfacts,intsimpfacts,theta,thv0,theta0_v0,
                       xobs,xgrid,xmin,xmax,xmid,xrange,xdelta,xinxgrid,xidelta,
                       fpm,alpha,uobs,phixobsfree,phixgridfree,phixobs,phixgrid,
                       iflag_AM,iflagpsi,met_beta,
                       psi,psi_m0,psi_s0,omega,omega_m0,omega_v0,met)
  
  fxobs=gtheta$fxobs
  fxgrid=gtheta$fxgrid
  theta=gtheta$theta
  psi=gtheta$psi
  omega=gtheta$omega
  met=gtheta$met
  
  #Update alpha
  xtx=sapply(1:nfun,function(ifun){sum((xobs[,ifun]-xmid[ifun])^2)})
  for (ifun in 1:nfun){
    resid=uobs-wb-rowSums(fxobs)
    if(fmodel[ifun]>2 && fmodel[ifun]<7) {
      # @ Take off old alpha*(x-xmid) from fx @
      fx1=fxobs[,ifun]-alpha[ifun]*(xobs[,ifun]-xmid[ifun])
      fxg1=fxgrid[,ifun]-alpha[ifun]*(xgrid[,ifun]-xmid[ifun])
      
      rk=resid-alpha[ifun]*(xobs[,ifun]-xmid[ifun])  # take off alpha*(x-xmid)
      
      # @ Generate alpha ~ Truncated Normal @
      a_vn=1/(xtx[ifun]+solve(alpha_s0^2))
      a_mn=a_vn*(crossprod((xobs[,ifun]-xmid[ifun]),rk)+alpha_m0/alpha_s0^2)
      if (fpm[ifun] == 1) {
        # @ Generate Normal truncated below at alpha > 0 @
        alpha[ifun]=rtruncnorm(1,a=0,b=Inf,mean=a_mn,sd=sqrt(a_vn))	
      } else {
        # @ Generate Normal truncated below at alpha < 0 @
        alpha[ifun]=rtruncnorm(1,a=-Inf,b=0,mean=a_mn,sd=sqrt(a_vn))	
      }
      fxobs[,ifun]=fx1+alpha[ifun]*(xobs[,ifun]-xmid[ifun])
      fxgrid[,ifun]=fxg1+alpha[ifun]*(xgrid[,ifun]-xmid[ifun])
    } else {
      alpha[ifun]=0
    }
  }
  
  #update tau2
  tau2=gtau2(nfun,theta,gamvec,iflagprior,tau2,tau2_s0,tau2_u0,tau2_r0,nbasis)

  #update gamma
  ggamma=ggamma(nfun,theta,tau2,gampar,gamvec,w0,nbasis,wk,zeta,kbar,kall)
  gampar=ggamma$gampar
  zeta=ggamma$zeta
  
  #update uobs
  u1=rtruncnorm(nobs,a=0,b=Inf,wb+rowSums(fxobs),1)
  u2=rtruncnorm(nobs,a=-Inf,b=0,wb+rowSums(fxobs),1)
  uobs=yobs*u1+(1-yobs)*u2
  
  #save outputs
  param.out=list(beta=beta,zeta=zeta,tau2=tau2,gampar=gampar,
                 theta=theta,alpha=alpha,psi=psi,omega=omega,
                 wb=wb,fxobs=fxobs,fxgrid=fxgrid,uobs=uobs)
  
  return(list(param.out=param.out,met.out=met))
}

#**************************************************************************
# * Generate theta, omega and psi
#**************************************************************************
'gtheta.probit'=function(yobs,nfun,wb,fxobs,fxgrid,fmodel,nbasis,nobs,nint,nr,met_alpha,
                         quadfacts,intsimpfacts,theta,thv0,theta0_v0,
                         xobs,xgrid,xmin,xmax,xmid,xrange,xdelta,xinxgrid,xidelta,
                         fpm,alpha,uobs,phixobsfree,phixgridfree,phixobs,phixgrid,
                         iflag_AM,iflagpsi,met_beta,
                         psi,psi_m0,psi_s0,omega,omega_m0,omega_v0,met)
{
  #get met
  met_var_all=met$met_var_all
  met_mAM=met$met_mAM
  met_beta_AM=met$met_beta_AM
  pmet=met$pmet
  icount_AM=met$icount_AM
  
  resid=uobs-wb-rowSums(fxobs)  # Resid takes off all of the f's
  ifree=1  # ith : no shape 
  irest=1  # ith : shape constraints
  
  for (ifun in 1:nfun) {
    testp=0   # Metropolis test p
    rk=resid+fxobs[,ifun] # takes off fk
    
    if (fmodel[ifun]==1){   #No shape restriction
      vni1=crossprod(phixobsfree[,,ifree])+diag(x=1/thv0[,ifun],nr=nbasis)
      vn1=solve(vni1)
      U=chol(vn1)
      bn1=vn1%*%crossprod(phixobsfree[,,ifree],rk) # Prior mean is zero for theta
      theta[2:(nbasis+1),ifun]=bn1+t(U)%*%as.matrix(rnorm(nbasis))
      theta[1,ifun]=0

      foo=GetFreef(theta[2:(nbasis+1),ifun],phixobsfree[,,ifree],phixgridfree[,,ifree],
                   nbasis,nobs,nint+1)
      fxobs[,ifun]=foo$fxobs
      fxgrid[,ifun]=foo$fxgrid	
      ifree=ifree+1
      
    } else {
    
    theta0_old=theta[1,ifun]
    theta_old=theta[2:(nbasis+1),ifun]
    
    # Random walk Metropolis for theta 
    met_var_all_new=met_beta_AM[ifun]/rgamma(1,shape=met_alpha)
    met_std0=sqrt(5.66*met_var_all_new*theta0_v0)
    met_std=sqrt(5.66*met_var_all_new*thv0[,ifun])
    
    theta0_new=rtruncnorm(1,a=0,b=Inf,theta0_old,met_std0)
    theta_new=rnorm(nbasis,theta_old,met_std)
    thetanew=c(theta0_new,theta_new)
    
    # @ Normalizing constant for generating theta_new @
    theta0_new_lnp=pnorm(-theta0_old/met_std0,0,1,lower.tail=F,log.p=T)
    # @ Normalizing constant for generating theta_old @
    theta0_old_lnp=pnorm(-theta0_new/met_std0,0,1,lower.tail=F,log.p=T)  	
    
    ck=0.5*(met_var_all_new + met_mAM[ifun])
    met_beta_new=(met_alpha-1)*ck
    
    if (fmodel[ifun] == 2){
      foo = GetUpf(fpm[ifun],thetanew,phixobs[,,irest],phixgrid[,,irest],quadfacts,nbasis,nr,
                   nobs,nint+1)
      fxobs_new=foo$fxobs
      fxgrid_new=foo$fxgrid
      
    } else if (fmodel[ifun] == 3) {
      foo = GetConvexf(fpm[ifun],alpha[ifun],thetanew,xobs[,ifun],xgrid[,ifun],xmid[ifun],
                       phixobs[,,irest],phixgrid[,,irest],quadfacts,nbasis,nr,nobs,nint+1)
      fxobs_new=foo$fxobs
      fxgrid_new=foo$fxgrid
      
    } else if (fmodel[ifun] == 4) {
      foo=GetConcavef(fpm[ifun],alpha[ifun],thetanew,xobs[,ifun],xgrid[,ifun],xmid[ifun],
                      phixobs[,,irest],phixgrid[,,irest],quadfacts,nbasis,nr,nobs,nint+1)
      fxobs_new=foo$fxobs
      fxgrid_new=foo$fxgrid
    }  else {
    	
      ##############################################
      # Generate omega and psi for squish function #
      ##############################################
      
      if (iflagpsi == 1){  # Generate psi > 0 

        met_stdS=sqrt(met_beta[ifun]/rgamma(1,shape=met_alpha))
        
        psi_old=psi[ifun]
        psi_new=rtruncnorm(1,a=0,b=Inf,mean= psi_old,sd= met_stdS) #truncated normal > 0
        psi_lnpnew=pnorm(-psi_old/met_stdS,lower.tail=F,log.p=T)  # @ LogNormalizing constant for g(new|old)
        psi_lnpold=pnorm(-psi_new/met_stdS,lower.tail=F,log.p=T)  # @ Log constant for g(old|new) @

        testp=testp-
          ((psi_new-psi_m0)^2)/(2*psi_s0^2)+  # prior for candidate psi_new     
          ((psi_old-psi_m0)^2)/(2*psi_s0^2)-  # prior for old psi                 
          psi_lnpold+ 							  # constant for g(old|new)			
          psi_lnpnew							      # constant for g(new|old)
      } else {  #fixed psi
        psi_old=psi[ifun]  # Placeholder
        psi_new=psi[ifun]  # Placeholder
      }
      
      # Generate xmin < omega < xmax
    
      met_stdS=sqrt(met_beta[ifun]/rgamma(1,shape=met_alpha))
      
      omega_old=omega[ifun]
      omega_new=rtruncnorm(1,a=xmin[ifun],b=xmax[ifun],omega_old,met_stdS) # truncated normal
      omega_lnpnew=log(pnorm((xmax[ifun]-omega_old)/met_stdS)-pnorm((xmin[ifun]-omega_old)/met_stdS))  # @ Normalizing constants for omega for g(new|old
      omega_lnpold=log(pnorm((xmax[ifun]-omega_new)/met_stdS)-pnorm((xmin[ifun]-omega_new)/met_stdS))  # @ Normalizing constant for g(old|new) @

      testp=testp- 
        ((omega_new-omega_m0)^2)/(2*omega_v0)+ # prior for candidate omega_new
        ((omega_old-omega_m0)^2)/(2*omega_v0)- # prior for old omega  
        omega_lnpold+								     # constant for g(old|new)			
        omega_lnpnew								     # constant for g(new|old) 
      
      if (fmodel[ifun] == 5) {
        foo=GetSf(fpm[ifun],omega_new,psi_new,alpha[ifun],thetanew,xobs[,ifun],xgrid[,ifun],                 phixobs[,,irest],phixgrid[,,irest],xdelta[ifun],xinxgrid[,ifun],xidelta[,ifun],
                  xrange[ifun],xmid[ifun],quadfacts,intsimpfacts,nbasis,nr,nobs,nint+1)
        fxobs_new=foo$fxobs
        fxgrid_new=foo$fxgrid
        
      }else if (fmodel[ifun] == 6) {
        foo=GetRotateSf(fpm[ifun],omega_new,psi_new,alpha[ifun],thetanew,xobs[,ifun],xgrid[,ifun],                    phixobs[,,irest],phixgrid[,,irest],xdelta[ifun],xinxgrid[,ifun],xidelta[,ifun],
                        xrange[ifun],xmid[ifun],quadfacts,intsimpfacts,nbasis,nr,nobs,nint+1)
        fxobs_new=foo$fxobs
        fxgrid_new=foo$fxgrid
        
      }else {
        foo=GetUf(fpm[ifun],omega_new,psi_new,thetanew,xobs[,ifun],xgrid[,ifun],phixobs[,,irest],
         phixgrid[,,irest],xdelta[ifun],xinxgrid[,ifun],xidelta[,ifun],xrange[ifun],xmid[ifun],
                  quadfacts,intsimpfacts,nbasis,nr,nobs,nint+1)
        fxobs_new=foo$fxobs
        fxgrid_new=foo$fxgrid
      }
    }
    
    # Metropolis Test for theta
    
    resid_new=rk-fxobs_new
    sse_new=sum(resid_new^2)
    resid_old=rk-fxobs[,ifun]
    sse_old=sum(resid_old^2)
    
    snew=sum(theta_new^2/thv0[,ifun])
    sold=sum(theta_old^2/thv0[,ifun])
    
    testp=testp - 
      (sse_new-sse_old)/2-	          # Likelihood 				
      (theta0_new^2)/(2*theta0_v0)+    # Prior for new theta0 		
      (theta0_old^2)/(2*theta0_v0)-	#	 ! Prior for old theta0 		
      (snew-sold)/2-				          # Prior for new & old theta	
      theta0_old_lnp+theta0_new_lnp-  #	! Constants for truncated r-w for theta0 
      ((met_alpha+1)*log(met_var_all[ifun]))-
      (met_beta_new/met_var_all[ifun])+ 
      ((met_alpha+1)*log(met_var_all_new))+
      (met_beta_AM[ifun]/met_var_all_new)
    
    if(log(runif(1)) <= testp) {
      
      theta[,ifun]=thetanew
      if (fmodel[ifun] == 5 | fmodel[ifun] == 6 | fmodel[ifun] == 7) {
        psi[ifun]=psi_new
        omega[ifun]=omega_new
      }
      met_var_all[ifun]=met_var_all_new
      fxobs[,ifun]=fxobs_new
      fxgrid[,ifun]=fxgrid_new
      pmet[ifun]=pmet[ifun]+1
    }
    
    if (iflag_AM == 1) {
      icount_AM=icount_AM+1
      # @ Update mean, variance, alpha, and beta for metropolis @
      met_mAM[ifun]=met_mAM[ifun]+(met_var_all[ifun]-met_mAM[ifun])/icount_AM
    }
    ck=0.5*(met_var_all[ifun]+met_mAM[ifun])
    met_beta_AM[ifun]=(met_alpha-1)*ck
    
    irest=irest+1
    }
    resid=uobs-wb-rowSums(fxobs)  # Resid takes off all of the f's
  }
  
  met=list(met_var_all=met_var_all,met_mAM=met_mAM,
           met_beta_AM=met_beta_AM,pmet=pmet,icount_AM=icount_AM)
  
  return(list(fxobs=fxobs,fxgrid=fxgrid,theta=theta,psi=psi,
              omega=omega,met=met))
}


#*********************************************************************
# *  Compute log prior density
#*********************************************************************
'GetLogPrior' = function(fmodel,nparw,nfun,nbasis,kall,w0,xmin,xmax,
                         beta,beta_m0,beta_v0,beta_v0i,
                         theta,theta0_m0,theta0_s0,theta0_v0,
                         gampar,tau2,tau2_r0,tau2_s0,tau2_u0,
                         alpha,alpha_m0,alpha_s0,
                         iflagprior,iflagpsi,psi,psi_m0,psi_s0,
                         omega,omega_m0,omega_s0)
{
  lnpriorf=0
  thetasq=theta^2
  thvar=numeric(nbasis+1)
  
  # @ Normal prior for beta ~ N(bn,Bn) @
  residb=beta-beta_m0
  lnpriorf=lnpriorf-t(residb)%*%beta_v0i%*%residb/2- nparw*log(2*pi)/2-log(det(beta_v0))/2
  
  for (ifun in 1:nfun){
    # @ Normal prior of theta   	
    if (fmodel[ifun] == 1) {
      thvar[2:(nbasis+1)]=tau2[ifun]*exp(-gampar[ifun]*kall)
      sse=sum(thetasq[2:(nbasis+1),ifun]/thvar[2:(nbasis+1)])
      lnpriorf=lnpriorf-sse/2-nbasis*log(2*pi)/2- 
        sum(log(thvar[2:(nbasis+1)]))/2
    } else {
      # @ Truncated Normal Prior for theta0 ~ N(0,tau2)*I(theta0 > 0) @
      thvar[1]=theta0_v0
      thvar[2:(nbasis+1)]=tau2[ifun]*exp(-gampar[ifun]*kall)
      sse=sum(thetasq[,ifun]/thvar)
      theta0_lnp0=pnorm(-theta0_m0/theta0_s0,lower.tail=F,log.p=T)
      lnpriorf=lnpriorf-sse/2-(nbasis+1)*log(2*pi)/2- 
        sum(log(thvar))/2-theta0_lnp0
    }
    
    # @ Exponential prior for gamma ~ EXP(w0) @
    lnpriorf=lnpriorf-w0*gampar[ifun]-log(w0)
    
    # @ Prior for tau2 @
    if (iflagprior == 1) {
      # @ tau2 ~ Exp(lambda) @
      lnpriorf=lnpriorf-tau2_u0*tau2[ifun]-tau2_u0
    } else { 
      # @ tau2 is IG @
      lnpriorf=lnpriorf+LogfIG(tau2[ifun],tau2_r0,tau2_s0)	
    }
    
    if (fmodel[ifun] > 2 && fmodel[ifun] < 7) {
      # @ Truncated normal prior for alpha~N(m0,v0)I(alpha>0) @
      alpha_lnp0	= pnorm(-alpha_m0/alpha_s0,lower.tail=F,log.p=T)
      lnpriorf=lnpriorf-((alpha[ifun]-alpha_m0)^2)/(2*alpha_s0^2)- 
        log(2*pi*alpha_s0^2)/2-alpha_lnp0
    }
    
    # Truncated normal priors for psi and omega
    if (fmodel[ifun] == 5|fmodel[ifun] == 6|fmodel[ifun] == 7) {
      if (iflagpsi == 1) {
        psi_lnp0	=pnorm(-psi_m0/psi_s0,lower.tail=F,log.p=T)   #Log Normalizing constant 
        lnpriorf=lnpriorf-((psi[ifun]-psi_m0)^2)/(2*psi_s0^2)- 
          log(2*pi*psi_s0^2)/2-psi_lnp0
      }
      omega_lnp0=log(pnorm((xmax-omega_m0)/(omega_s0)) - pnorm( (xmin-omega_m0)/(omega_s0)))
      lnpriorf=lnpriorf-((omega[ifun]-omega_m0)^2)/(2*omega_s0^2)-
        log(2*pi*omega_s0^2)/2-omega_lnp0
    }
  }
  
  return(lnpriorf)
}
