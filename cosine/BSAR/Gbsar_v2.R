'Gbsar'=function(yobs,wdata=NULL,xobs,family,link,method,fmodel,fpm,nbasis=50,nint=200,
                 prior=list(),mcmc=list())
{
  
  # Data 
  if(!is.matrix(yobs)) yobs=as.matrix(yobs)
  nobs=nrow(yobs)
  if(is.null(wdata)) {wdata=matrix(1,nobs,1)
  } else {
    if(!is.matrix(wdata)) wdata=as.matrix(wdata)
    wdata=cbind(1,wdata)
  }
  nparw=ncol(wdata)
  if(!is.matrix(xobs)) xobs=as.matrix(xobs)
  nfun=ncol(xobs)
  
  # Sampling point for f
  xmin=apply(xobs,2,min)
  xmax=apply(xobs,2,max)
  xmid=(xmin+xmax)/2
  xrange=xmax-xmin
  xdelta=(xrange)/nint
  xgrid=sapply(1:nfun,function(i)seq(xmin[i],xmax[i],length=nint+1))
  
  xinxgrid=matrix(0,nr=nobs,nc=nfun)
  xidelta=matrix(0,nr=nobs,nc=nfun)
  for(k in 1:nfun){
    if (fmodel[k] == 5 | fmodel[k] == 6 | fmodel[k] == 7) {
      
      s = seq(1,by=1,length.out=nint+1)
      
      for (i0 in 1:nobs){ 
        i = i0
        xi = xobs[i,k]
        
        if( xi == xmin[k]){
          xinxgrid[i,k] = 1
        }else if (xi == xmax[k]){
          xinxgrid[i,k] = nint+1
        }else{
          #@ Find index of xgrid for largest xgrid < xi @
          zbot = xgrid[,k] <= xi
          szbot = max(s[which(zbot)])
          xinxgrid[i,k] = szbot
          if (xi > xgrid[szbot,k]){
            xidelta[i,k] = xi - xgrid[szbot,k]
          }
        }  
      }
    }
  }
  
  # Prior
  privals=list(iflagprior=0,theta_m0=numeric(nbasis+1),theta0_m0=0,theta0_s0=100,
               tau2_m0=1,tau2_v0=100,w0=2,beta_m0=numeric(nparw),beta_v0=diag(100,nparw),
               alpha_m0=3,alpha_s0=50,iflagpsi=1,psifixed=100,
               omega_m0=mean(xmin+xmax)/2,omega_s0=mean(xrange/8),
               kappa_m0=1,kappa_v0=100)
  privals[names(prior)]=prior
  
  # MCMC
  mcvals=list(nblow0=1000,nblow=1000,smcmc=1000,nskip=10,ndisp=1000,maxmodmet=5)
  mcvals[names(mcmc)]=mcmc
  
  #############
  # Draw MCMC #
  #############
  
  if(family=='bernoulli'){
    if(link=='probit'){
      if(method=='MH'){
        foo=Gbsar_probit_MH(yobs,wdata,xobs,nobs,nparw,nfun,nbasis,nint,fmodel,
                            fpm,theta_m0,theta0_m0,theta0_s0,tau2_m0,tau2_v0,
                            w0,beta_m0,beta_v0,alpha_m0,alpha_s0,psi_m0,psi_s0,
                            psifixed,omega_m0,omega_s0,xinxgrid,xidelta,
                            iflagprior,iflagpsi,maxmodmet,nblow0,nblow,
                            smcmc,nskip,ndisp)
        
      } else if(method=='DataAug'){
        foo=Gbsar_probit(yobs,wdata,xobs,nobs,nparw,nfun,nbasis,nint,fmodel,fpm,privals,mcvals,
                         xmin,xmax,xmid,xrange,xdelta,xgrid,xinxgrid,xidelta)
      }
      
    } else if(link=='logit'){
      if(method=='MH'){
        foo=Gbsar_logit_MH(yobs,wdata,xobs,nobs,nparw,nfun,nbasis,nint,fmodel,fpm,
                           theta_m0,theta0_m0,theta0_s0,tau2_m0,tau2_v0,w0,beta_m0,beta_v0,
                           alpha_m0,alpha_s0,psi_m0,psi_s0,psifixed,omega_m0,omega_s0,
                           xinxgrid,xidelta,iflagprior,iflagpsi,
                           maxmodmet,nblow0,nblow,smcmc,nskip,ndisp)
      }else if(method=='DataAug'){
        foo=Gbsar_logit(yobs,wdata,xobs,nobs,nparw,nfun,nbasis,nint,fmodel,fpm,privals,mcvals,
                        xmin,xmax,xmid,xrange,xdelta,xgrid,xinxgrid,xidelta)
      }
    }
  } else if(family=='poisson'){
    foo=Gbsar_pois(yobs,wdata,xobs,nobs,nparw,nfun,nbasis,nint,fmodel,fpm,privals,mcvals,
                   xmin,xmax,xmid,xrange,xdelta,xgrid,xinxgrid,xidelta)
  } else if(family=='negbin'){
    if(method=='MH'){
      # MCMC using Negative Binomial likelihood
      cat("MCMC using Negative Binomial likelihood",'\n')
      foo=Gbsar_negbinMH(yobs,wdata,xobs,nobs,nparw,nfun,nbasis,nint,fmodel,
                         fpm,theta_m0,theta0_m0,theta0_s0,tau2_m0,tau2_v0,
                         w0,beta_m0,beta_v0,alpha_m0,alpha_s0,psi_m0,psi_s0,
                         psifixed,omega_m0,omega_s0,kappa_m0,kappa_v0,
                         xinxgrid,xidelta,iflagprior,iflagpsi,
                         maxmodmet,nblow0,nblow,smcmc,nskip,ndisp)
      
    } else if(method=='PoisGam'){
      # MCMC using Poisson-Gamma mixture
      cat("MCMC using Poisson-Gamma mixture",'\n')
      foo=bsarnb2(yobs,wdata,xobs,nobs,nparw,nfun,nbasis,nint,fmodel,fpm,
                  theta_m0,theta0_m0,theta0_s0,tau2_m0,tau2_v0,w0,beta_m0,beta_v0,
                  alpha_m0,alpha_s0,psi_m0,psi_s0,psifixed,omega_m0,omega_s0,
                  xinxgrid,xidelta,iflagprior,iflagpsi,rho_a0,rho_b0,
                  maxmodmet,nblow0,nblow,smcmc,nskip,ndisp)
      
    }else if(method=='PG'){
      # MCMC using Polya-Gamma mixture
      cat("MCMC using Polya-gamma mixture",'\n')
      foo=bsarnb3(yobs,wdata,xobs,nobs,nparw,nfun,nbasis,nint,fmodel,fpm,
                  theta_m0,theta0_m0,theta0_s0,tau2_m0,tau2_v0,w0,beta_m0,beta_v0,
                  alpha_m0,alpha_s0,psi_m0,psi_s0,psifixed,omega_m0,omega_s0,
                  xinxgrid,xidelta,iflagprior,iflagpsi,rho_a0,rho_b0,
                  maxmodmet,nblow0,nblow,smcmc,nskip,ndisp)
    }
  }
  
  mcmc.draws = list()
  mcmc.draws$zeta=foo$zetag
  mcmc.draws$tau2=foo$tau2g
  mcmc.draws$alpha=foo$alphag
  mcmc.draws$psi=foo$psig
  mcmc.draws$omega=foo$omegag
  mcmc.draws$gamma=foo$gammag
  mcmc.draws$theta=foo$thetag
  mcmc.draws$beta=foo$betag
  if (family=='negbin'){
      mcmc.draws$kappag=foo$kappag
  }
  
  loglik.draws = list()
  loglik.draws$loglike = foo$loglikeg
  loglik.draws$logprior = foo$logpriorg
  loglik.draws$logjoint = foo$loglikeg + foo$logpriorg
  
  fit.draws = list()
  fit.draws$xgrid = xgrid
  fit.draws$fxobs = foo$fxobsg
  fit.draws$fxgrid = foo$fxgridg
  fit.draws$wbeta = foo$wbg
  
  post.est = list()
  muhatm = apply(foo$muhatg, 2, mean)
  muhats = apply(foo$muhatg, 2, sd)
  post.est$muhatm = muhatm
  post.est$muhats = muhats
  
  if (family=='bernoulli'){
    phatm = apply(foo$phatg, 2, mean)
    phats = apply(foo$phatg, 2, sd)
    post.est$phatm = phatm
    post.est$phats = phats
  }
  
  betam = apply(foo$betag, 2, mean)
  betas = apply(foo$betag, 2, sd)
  post.est$betam = betam
  post.est$betas = betas
  thetam = apply(foo$thetag, c(1, 2), mean)
  thetas = apply(foo$thetag, c(1, 2), sd)
  post.est$thetam = thetam
  post.est$thetas = thetas
  tau2m = colMeans(foo$tau2g)
  tau2s = apply(foo$tau2g, 2, sd)
  post.est$tau2m = tau2m
  post.est$tau2s = tau2s
  taug = sqrt(foo$tau2g)
  taum = colMeans(taug)
  taus = apply(taug, 2, sd)
  post.est$taum = taum
  post.est$taus = taus
  gammam = colMeans(foo$gammag)
  gammas = apply(foo$gammag, 2, sd)
  post.est$gammam = gammam
  post.est$gammas = gammas
  alpham = colMeans(foo$alphag)
  alphas = apply(foo$alphag, 2, sd)
  post.est$alpham = alpham
  post.est$alphas = alphas
  psim = colMeans(foo$psig)
  psis = apply(foo$psig, 2, sd)
  post.est$psim = psim
  post.est$psis = psis
  omegam = colMeans(foo$omegag)
  omegas = apply(foo$omegag, 2, sd)
  post.est$omegam = omegam
  post.est$omegas = omegas
  zetam = colMeans(foo$zetag)
  zetas = apply(foo$zetag, 2, sd)
  post.est$zetam = zetam
  post.est$zetas = zetas
  if(family=="negbin"){
    kappam=mean(foo$kappag)
    kappas=mean(foo$kappag)
    post.est$kappam = kappam
    post.est$kappas = kappas
  }
  
  # ***************************
  # * CPO - Statistics & LPML
  # ***************************
  invlikeg=foo$invlikeg
  cpo=1/apply(invlikeg,2,mean)
  lpml=sum(log(cpo))
  
  res.out=list()
  res.out$family=family
  res.out$yobs=yobs
  res.out$xobs=xobs
  res.out$wdata=wdata
  res.out$xgrid=xgrid
  res.out$fmodel=fmodel
  res.out$fpm=fpm
  res.out$prior=privals
  res.out$mcmc=mcvals
  res.out$nbasis=nbasis
  res.out$nfun=nfun
  res.out$xmin=xmin
  res.out$xmax=xmax
  res.out$nobs=nobs
  res.out$nparw=nparw
  res.out$pmet=foo$pmetg
  res.out$mcmc.draws = mcmc.draws
  res.out$fit.draws = fit.draws
  res.out$loglik.draws = loglik.draws
  res.out$post.est = post.est
  res.out$lpml = lpml
  res.out$invlikeg = invlikeg
  
  res.out
}

  

