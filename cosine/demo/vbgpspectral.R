vbgpspectral<-function(y,x,Z,T,tol,prior.parms,mupsi.q.start) {

  # Inputs:
  #
  # y - vector of responses
  # X - vector of predictors;  corresponding to the smooth term in the semiparametric model
  # Z - design matrix;  corresponding to parametric part in the semiparametric model
  # T - number of spectral basis functions
  # tol - tolerance for stopping of VB algorithm
  # prior.parms - prior parameters.  Needs to be a list with components:
  #    rsig.0, ssig.0; the prior on sigma^2 is IG(rsig.0/2,ssig.0/2)
  #    rtau.0, stau.0; the prior on tau^2 is IG(rtau.0/2,stau.0/2)
  #    w0 - the prior on gamma is Exp(w0)
  #    mubeta.0,sigbeta.0;  the prior on beta is normal, parameters mubeta.0 and sigma^2*sigbeta.0
  # 

  rsig.0<-prior.parms$rsig.0
  ssig.0<-prior.parms$ssig.0
  rtau.0<-prior.parms$rtau.0
  stau.0<-prior.parms$stau.0
  w0<-prior.parms$w0
  mubeta.0<-prior.parms$mubeta.0
  sigbeta.0<-prior.parms$sigbeta.0

  n<-length(y)
  if(is.matrix(Z)==FALSE) Z<-matrix(Z,nrow=n,ncol=1)
  p<-ncol(Z)

  mupsi.q<-mupsi.q.start
  mu2psi.q<-mupsi.q^2
  sig2psi.q<-mu2psi.q/100
  sigpsi.q<-sqrt(sig2psi.q)
  ssig.q<-ssig.0
  stau.q<-stau.0
  mubeta.q<-mubeta.0
  sigbeta.q<-sigbeta.0
  sigbeta0.inv<-solve(sigbeta.0)
  ldetsigbeta.0<-as.numeric(determinant(sigbeta.0,logarithm=TRUE)$modulus)
  sb0ibymub0<-sigbeta0.inv%*%mubeta.0
  lbold<- -Inf
  dif<-tol+1
 
  ZtZ<-t(Z)%*%Z
  lb<-c()
  const1<-sqrt(2/pi)

  Tfull<-T
  T.old<-T
  vphifull<-sqrt(2)*cos(outer(x,pi*(1:T)))
  vphitvphifull<-t(vphifull)%*%vphifull
  sigbeta.term<-solve(ZtZ+sigbeta0.inv)

  while(dif>tol | T.old!=T) {

    # Update variational distribution parameters for theta

    T.old<-T
    T<-min(floor(-15/(-0.5*mupsi.q)),Tfull)    

    bindices<-(1:T)
    bindices2<-bindices^2
    mfactor<- -T*(T+1)/(4*w0)
    rsig.q<-rsig.0+T+p+n
    rtau.q<-rtau.0+T

    # Set up spectral design matrix

    vphi<-vphifull[,1:T]
    vphitvphi<-vphitvphifull[1:T,1:T]

    # Update variational distribution parameters for theta

    term1<-exp(sig2psi.q*bindices2/2+mupsi.q*bindices)
    term2<-exp(sig2psi.q*bindices2/2-mupsi.q*bindices)
    sigpsi.q<-sqrt(sig2psi.q)
    Qvec<-term1*(1-pnorm(-mupsi.q/sigpsi.q-sigpsi.q*bindices))
    Qvec<-Qvec+term2*(1-pnorm(mupsi.q/sigpsi.q-sigpsi.q*bindices))
    sig.ratio<-rsig.q/ssig.q
    tau.ratio<-rtau.q/stau.q
    DQvec<-diag(Qvec)

    sigtheta.q<-vphitvphi+tau.ratio*DQvec
    sigtheta.q<-1/sig.ratio*solve(sigtheta.q)
    res.temp<-y-Z%*%mubeta.q
    mutheta.q<-drop(sig.ratio*sigtheta.q%*%t(vphi)%*%res.temp)
    
    # Update variational distribution parameters for sigma

    trterm1<-sum(diag((sigtheta.q+outer(mutheta.q,mutheta.q))%*%DQvec))
    trterm2<-sum(diag(ZtZ%*%sigbeta.q))
    trterm3<-sum(diag(vphitvphi%*%sigtheta.q))
    trterm4<-sum(diag(solve(sigbeta.0,sigbeta.q)))
    ssterm<-sum((y-Z%*%mubeta.q-vphi%*%mutheta.q)^2)
    ssterm2<-sum((mubeta.q-mubeta.0),solve(sigbeta.0,(mubeta.q-mubeta.0)))
    ssig.q<-ssig.0+tau.ratio*trterm1+trterm2+trterm3+trterm4+ssterm+ssterm2
    sig.ratio<-rsig.q/ssig.q
    
    # Update variational distribution parameters for tau

    stau.q<-stau.0+sig.ratio*trterm1
    tau.ratio<-rtau.q/stau.q

    # Update variational distribution parameters for beta

    sigbeta.q<-1/sig.ratio*sigbeta.term
    mubeta.q<-sig.ratio*sigbeta.q%*%(sb0ibymub0+t(Z)%*%(y-vphi%*%mutheta.q))
    mubeta.q<-as.vector(mubeta.q)

    # Calculate lower bound value before NCVMP step

    lbnew<-0
    sigpsi.q<-sqrt(sig2psi.q)
    mu2psi.q<-mupsi.q^2
    trterm2<-sum(diag(ZtZ%*%sigbeta.q))
    trterm3<-sum(diag(vphitvphi%*%sigtheta.q))
    ssterm<-sum((y-Z%*%mubeta.q-vphi%*%mutheta.q)^2)
    lbnew<-lbnew-n/2*log(2*pi)-n/2*(log(ssig.q/2)-digamma(rsig.q/2))-sig.ratio/2*(trterm2+trterm3+ssterm)
    lbnew<-lbnew-p/2*log(2*pi)-0.5*ldetsigbeta.0-p/2*(log(ssig.q/2)-digamma(rsig.q/2))-0.5*sig.ratio*(sum(diag(solve(sigbeta.0,sigbeta.q)))+sum((mubeta.q-mubeta.0)*solve(sigbeta.0,mubeta.q-mubeta.0)))
    lbnew<-lbnew-T/2*(log(2*pi)+log(ssig.q/2)-digamma(rsig.q/2)+log(stau.q/2)-digamma(rtau.q/2))
    lbnew<-lbnew+rsig.0/2*log(ssig.0/2)-lgamma(rsig.0/2)-(rsig.0/2+1)*(log(ssig.q/2)-digamma(rsig.q/2))-ssig.0/2*rsig.q/ssig.q
    lbnew<-lbnew+rtau.0/2*log(stau.0/2)-lgamma(rtau.0/2)-(rtau.0/2+1)*(log(stau.q/2)-digamma(rtau.q/2))-stau.0/2*rtau.q/stau.q
    lbnew<-lbnew-(-p/2*log(2*pi)-0.5*as.numeric(determinant(sigbeta.q,logarithm=TRUE)$modulus)-p/2)
    lbnew<-lbnew-(-T/2*log(2*pi)-0.5*as.numeric(determinant(sigtheta.q,logarithm=TRUE)$modulus)-T/2)
    lbnew<-lbnew-(rsig.q/2*log(ssig.q/2)-lgamma(rsig.q/2)-(rsig.q/2+1)*(log(ssig.q/2)-digamma(rsig.q/2))-rsig.q/2)
    lbnew<-lbnew-(rtau.q/2*log(stau.q/2)-lgamma(rtau.q/2)-(rtau.q/2+1)*(log(stau.q/2)-digamma(rtau.q/2))-rtau.q/2)
    lbnew.wopsibits<-lbnew
    lbnew<-lbnew-(-0.5*log(2*pi)-0.5*log(sig2psi.q)-0.5)
    S1<- -w0*(sigpsi.q*const1*exp(-mu2psi.q/(2*sig2psi.q))+mupsi.q*(1-2*pnorm(-mupsi.q/sigpsi.q)))
    lbnew<-lbnew+T*(T+1)/4*(sigpsi.q*const1*exp(-mu2psi.q/(2*sig2psi.q))+mupsi.q*(1-2*pnorm(-mupsi.q/sigpsi.q)))
    Qvec<-term1*(1-pnorm(-mupsi.q/sigpsi.q-sigpsi.q*bindices))
    Qvec<-Qvec+term2*(1-pnorm(mupsi.q/sigpsi.q-sigpsi.q*bindices)) 
    lbnew<-lbnew-0.5*tau.ratio*sig.ratio*sum(diag((sigtheta.q+outer(mutheta.q,mutheta.q))%*%diag(Qvec)))
    lbnew<-lbnew+log(w0/2)+S1
    lbfull<-lbnew

    # Update variational distribution parameters for psi

    sig3psi.q<-sigpsi.q^3
    sig4psi.q<-sigpsi.q^4
    mu2psi.q<-mupsi.q^2
    const1<-sqrt(2/pi)
    const2<-exp(-mu2psi.q/(2*sig2psi.q))
    pd1sig<- -w0*((1/(2*sig2psi.q)+sigpsi.q*mu2psi.q/(2*sig4psi.q))*const1*const2-mu2psi.q/sig3psi.q*dnorm(-mupsi.q/sigpsi.q))
    pd1mu<- -w0*(-mupsi.q/sigpsi.q*const1*const2+(1-2*pnorm(-mupsi.q/sigpsi.q))+2*mupsi.q/sigpsi.q*dnorm(-mupsi.q/sigpsi.q))
    
    term1<-exp(sig2psi.q*bindices2/2+mupsi.q*bindices)
    term2<-exp(sig2psi.q*bindices2/2-mupsi.q*bindices)
    term3<-mupsi.q/(2*sig3psi.q)-bindices/(2*sig2psi.q)
    term4<-mupsi.q/(2*sig3psi.q)+bindices/(2*sig2psi.q)
    term5<-dnorm(-mupsi.q/sigpsi.q-sigpsi.q*bindices)
    term6<-dnorm(mupsi.q/sigpsi.q-sigpsi.q*bindices)
    term7<-1-pnorm(-mupsi.q/sigpsi.q-sigpsi.q*bindices)
    term8<-1-pnorm(mupsi.q/sigpsi.q-sigpsi.q*bindices)
    DQmu<-term1*term5*1/sigpsi.q+bindices*term1*term7-term2*term6*1/sigpsi.q-bindices*term2*term8    
    DQsig<- -term1*term5*term3+bindices2/2*term1*term7+term2*term6*term4+bindices2/2*term2*term8    
    pd2sig<- mfactor*pd1sig-0.5*sig.ratio*tau.ratio*sum((diag(sigtheta.q)+mutheta.q^2)*DQsig)
    pd2mu<- mfactor*pd1mu-0.5*sig.ratio*tau.ratio*sum((diag(sigtheta.q)+mutheta.q^2)*DQmu)
    sig2psi.q.old<-sig2psi.q
    sig2psi.q<- -0.5/(pd1sig+pd2sig)
    mupsi.q.old<-mupsi.q
    mupsi.q<-mupsi.q+sig2psi.q*(pd1mu+pd2mu)

    # Calculate lower bound value and then check that NCVMP step increases lower bound.
    # If not, then do step halving.

    sigpsi.q<-sqrt(sig2psi.q)
    mu2psi.q<-mupsi.q^2
    term1<-exp(sig2psi.q*bindices2/2+mupsi.q*bindices)
    term2<-exp(sig2psi.q*bindices2/2-mupsi.q*bindices)

    lbnew<-lbnew.wopsibits
    lbnew<-lbnew-(-0.5*log(2*pi)-0.5*log(sig2psi.q)-0.5)
    S1<- -w0*(sigpsi.q*const1*exp(-mu2psi.q/(2*sig2psi.q))+mupsi.q*(1-2*pnorm(-mupsi.q/sigpsi.q)))
    lbnew<-lbnew+T*(T+1)/4*(sigpsi.q*const1*exp(-mu2psi.q/(2*sig2psi.q))+mupsi.q*(1-2*pnorm(-mupsi.q/sigpsi.q)))
    Qvec<-term1*(1-pnorm(-mupsi.q/sigpsi.q-sigpsi.q*bindices))
    Qvec<-Qvec+term2*(1-pnorm(mupsi.q/sigpsi.q-sigpsi.q*bindices)) 
    lbnew<-lbnew-0.5*tau.ratio*sig.ratio*sum(diag((sigtheta.q+outer(mutheta.q,mutheta.q))%*%diag(Qvec)))
    lbnew<-lbnew+log(w0/2)+S1
    
    dif<-lbnew-lbfull

    if(dif<0) {
      step<-1
      dif.try<-dif
      while(dif.try<0) {
        step<-step*0.5
        sig2psi.q.try<-1/(1/sig2psi.q.old+step*(1/sig2psi.q-1/sig2psi.q.old))
        sigpsi.q.try<-sqrt(sig2psi.q.try)
        mupsi.q.try<-sig2psi.q.try*(mupsi.q.old/sig2psi.q.old+step*(mupsi.q/sig2psi.q-mupsi.q.old/sig2psi.q.old))
        mu2psi.q.try<-mupsi.q.try^2
        term1<-exp(sig2psi.q.try*bindices2/2+mupsi.q.try*bindices)
        term2<-exp(sig2psi.q.try*bindices2/2-mupsi.q.try*bindices)
        lbnew<-lbnew.wopsibits
        lbnew<-lbnew-(-0.5*log(2*pi)-0.5*log(sig2psi.q.try)-0.5)
        S1<- -w0*(sigpsi.q.try*const1*exp(-mu2psi.q.try/(2*sig2psi.q.try))+mupsi.q.try*(1-2*pnorm(-mupsi.q.try/sigpsi.q.try)))
        lbnew<-lbnew+T*(T+1)/4*(sigpsi.q.try*const1*exp(-mu2psi.q.try/(2*sig2psi.q.try))+mupsi.q.try*(1-2*pnorm(-mupsi.q.try/sigpsi.q.try)))
        Qvec<-term1*(1-pnorm(-mupsi.q.try/sigpsi.q.try-sigpsi.q.try*bindices))
        Qvec<-Qvec+term2*(1-pnorm(mupsi.q.try/sigpsi.q.try-sigpsi.q.try*bindices)) 
        lbnew<-lbnew-0.5*tau.ratio*sig.ratio*sum(diag((sigtheta.q+outer(mutheta.q,mutheta.q))%*%diag(Qvec)))
        lbnew<-lbnew+log(w0/2)+S1
        dif.try<-lbnew-lbfull
      }

      sigpsi.q<-sigpsi.q.try
      sig2psi.q<-sig2psi.q.try
      mupsi.q<-mupsi.q.try
      mu2psi.q<-mu2psi.q.try
      
    }    

    dif<-lbnew-lbold
    dif<-dif/abs(lbnew)
    lbold<-lbnew
    lb<-c(lb,lbnew)
    
  }

  # Return results

  return(list(lb=lb,mutheta.q=mutheta.q,sigtheta.q=sigtheta.q,rsig.q=rsig.q,ssig.q=ssig.q,rtau.q=rtau.q,stau.q=stau.q,sigbeta.q=sigbeta.q,mubeta.q=mubeta.q,sig2psi.q=sig2psi.q,mupsi.q=mupsi.q))

}



