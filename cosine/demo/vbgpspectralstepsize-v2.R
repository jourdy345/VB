vbgpspectralsc.step.v2<-function(y,x,Z,T,tol,prior.parms,delta,mupsi.q.start,mutheta.q.start,n.grid,stepsize) {

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
  #    mubeta.0,sigbeta.0;  the prior on beta is normal parameters mubeta.0 and sigma^2*sigbeta.0
  #    sig02;  the prior variance on the intercept in the shape restricted term
  #      is sigma*sig02
  # delta - delta=1 for increasing function and -1 for decreasing
  # mupsi.q.start - starting value for variational posterior mean for psi
  # mutheta.q.start - starting value for variational posterior mean for theta
  # n.grid - design matrices are computed for a grid of x values with n.grid points
  #    and then interpolated if n>n.grid
  # stepsize - step size parameter for NCVMP steps for theta
  #

  maxit<-200
  rsig.0<-prior.parms$rsig.0
  ssig.0<-prior.parms$ssig.0
  rtau.0<-prior.parms$rtau.0
  stau.0<-prior.parms$stau.0
#  stau.q<-stau.0
  stau.q<-rtau.0+T
  w0<-prior.parms$w0
  mubeta.0<-prior.parms$mubeta.0
  sigbeta.0<-prior.parms$sigbeta.0
  sig02<-prior.parms$sig02

  n<-length(y)
  if(is.matrix(Z)==FALSE) Z<-matrix(Z,nrow=n,ncol=1)
  p<-ncol(Z)
  e1oversig<-10
  e1oversig2<-100
  mupsi.q<-mupsi.q.start
  mu2psi.q<-mupsi.q^2
  sig2psi.q<-mu2psi.q/(T*T)
  sigpsi.q<-sqrt(sig2psi.q)
  mubeta.q<-mean(y)
  sigbeta.q<-sigbeta.0
  sigtheta.q<-matrix(0,nrow=(T+1),ncol=(T+1))
  mutheta.q<-mutheta.q.start
    
  mutheta.q.old<-mutheta.q
  sigbeta0.inv<-solve(sigbeta.0)
  ldetsigbeta.0<-as.numeric(determinant(sigbeta.0,logarithm=TRUE)$modulus)
  sb0ibymub0<-sigbeta0.inv%*%mubeta.0
  lbold<- -Inf
  dif<-tol+1
 
  ZtZ<-t(Z)%*%Z
  lb<-c()

  Tfull<-T
  T.old<-T

  # Set up design matrices

  if(n<=n.grid) {
    x.grid<-x
    n.grid<-n
  } else {
    x.grid<-seq(from=0,to=1,length=n.grid)
    index<-round(approx(x.grid,1:length(x.grid),xout=x,rule=2,ties='ordered')$y)
    wts<-tabulate(index,n.grid) 
  }

  dmatsfull<-array(0,dim=c(n.grid,T+1,T+1))
  dmatsfull[,1,1]<-x.grid-0.5
  temparg1<-outer(x.grid,(1:T))
  temparg2<-matrix(rep((1:T),each=n.grid),nrow=n.grid,byrow=FALSE)
  dmatsfull[,1,(2:(T+1))]<-sqrt(2)/(pi*temparg2)*sin(pi*temparg1)-sqrt(2)/((pi*temparg2)^2)*(1-cos(pi*temparg2))
  dmatsfull[,(2:(T+1)),1]<-dmatsfull[,1,(2:(T+1))]
  temparg1<-outer((1:T),(1:T),'+')
  temparg2<-outer(rep(1,times=n.grid),temparg1)
  temparg1<-outer(x.grid,temparg1)
  dmatsfull[,(2:(T+1)),(2:(T+1))]<-sin(pi*temparg1)/(pi*temparg2)-(1-cos(pi*temparg2))/(pi*temparg2)^2
  temparg1<-outer((1:T),(1:T),'-')
  temparg2<-outer(rep(1,times=n.grid),temparg1)
  temparg1<-outer(x.grid,temparg1)
  dmatsfull[,(2:(T+1)),(2:(T+1))]<-dmatsfull[,(2:(T+1)),(2:(T+1))]+sin(pi*temparg1)/(pi*temparg2)-(1-cos(pi*temparg2))/(pi*temparg2)^2
  temparg1<-outer(x.grid,(1:T))
  temparg2<-outer(rep(1,times=n.grid),(1:T))
  temparg3<-outer(x.grid,rep(1,times=T))
  tempmat<-sin(2*pi*temparg1)/(2*pi*temparg2)+temparg3-0.5
  for(j in 2:(T+1)) {
    dmatsfull[,j,j]<-tempmat[,(j-1)]
  }
  
  sigbeta.term<-solve(ZtZ+sigbeta0.inv)

  matfun1<-function(psi,sigma) {
    return(4*psi%*%sigma%*%psi)
  }

  matfun2<-function(psi,mu) {
    return(4*psi%*%outer(mu,mu)%*%psi)
  }

  matfun3<-function(psi,sigma,mu) {
    return(sum(diag(sigma%*%psi))+sum(mu*(psi%*%mu)))
  }

  matfun4<-function(psi,sigma,mu) {
    return(8*psi%*%sigma%*%psi%*%mu)
  }

  matfun5<-function(psi,mu) {
    return(psi%*%mu)
  }

  matfun6<-function(psi,sigma,mu) {
    return(2*sum(diag(psi%*%sigma%*%psi%*%sigma))+4*sum(mu*(psi%*%sigma%*%psi%*%mu)))
  }

  const1<-sqrt(2/pi)
   
  while((dif>tol | T.old!=T) & length(lb)<maxit) {

    # Update variational distribution parameters for theta

    T.old<-T
    T<-min(floor(-15/(-0.5*mupsi.q)),Tfull)  
    a<-(rsig.0+n+p+(T+1)/2)/2+1 
    if(T<T.old) {
      sigtheta.q<-sigtheta.q[1:(T+1),1:(T+1)]
      mutheta.q<-mutheta.q[1:(T+1)]
      if(length(lb)>0) sigma.inv.old<-sigma.inv.old[1:(T+1),1:(T+1)]
      mutheta.q.old<-mutheta.q.old[1:(T+1)]
    }
    if(T>T.old) {
      sigtheta.q.temp<-matrix(0,nrow=(T+1),ncol=(T+1))
      sigtheta.q.temp[1:(T.old+1),1:(T.old+1)]<-sigtheta.q
      sigtheta.q<-sigtheta.q.temp
      mutheta.q.temp<-rep(0,times=(T+1))
      mutheta.q.temp[1:(T.old+1)]<-mutheta.q
      mutheta.q<-mutheta.q.temp
      if(length(lb)>0) {
        sigma.inv.old.temp<-matrix(0,nrow=(T+1),ncol=(T+1))
        sigma.inv.old.temp[1:(T.old+1),1:(T.old+1)]<-sigma.inv.old
        sigma.inv.old<-sigma.inv.old.temp
      }
      mutheta.q.old.temp<-rep(0,times=(T+1))
      mutheta.q.old.temp[1:(T.old+1)]<-mutheta.q.old
      mutheta.q.old<-mutheta.q.old.temp
    }
    dmats<-dmatsfull[,(1:(T+1)),(1:(T+1))]
   
     
    # Variational update for theta variational parameters

    bindices<-(1:T)
    bindices2<-bindices^2
    mfactor<- -T*(T+1)/(4*w0)
    rtau.q<-rtau.0+T
    tau.ratio<-rtau.q/stau.q

    term1<-exp(sig2psi.q*bindices2/2+mupsi.q*bindices)
    term2<-exp(sig2psi.q*bindices2/2-mupsi.q*bindices)
    sigpsi.q<-sqrt(sig2psi.q)
    Qvec<-term1*(1-pnorm(-mupsi.q/sigpsi.q-sigpsi.q*bindices))
    Qvec<-Qvec+term2*(1-pnorm(mupsi.q/sigpsi.q-sigpsi.q*bindices))

    e1overupsilon<-c(1/sig02,tau.ratio*Qvec)
    sigtheta.q.new<- -0.5*e1oversig*diag(e1overupsilon)
    if(length(lb)==0) {
      sigtheta.q<-1/e1oversig*diag(1/e1overupsilon)
      sigma.inv.old<-e1oversig*diag(e1overupsilon)
    }

    term1<-array(t(apply(dmats,1,matfun1,sigtheta.q)),dim=c(n.grid,T+1,T+1))
    term2<-array(t(apply(dmats,1,matfun2,mutheta.q)),dim=c(n.grid,T+1,T+1))
    term3<-apply(dmats,1,matfun3,sigtheta.q,mutheta.q)
    if(n<=n.grid) {
      data.term<-apply(term1,c(2,3),sum)+apply(term2,c(2,3),sum)
#      data.term<-data.term-apply(as.vector(2*delta*(y-Z%*%mubeta.q-delta*term3))*dmats,c(2,3),sum)
    } else {
      temp1<-aggregate(I(2*(y-Z%*%mubeta.q)) ~ index, FUN=sum)
      temp2<-rep(0,times=n.grid)
      temp2[temp1[,1]]<-temp1[,2]
      data.term<-apply(wts*term1,c(2,3),sum)+apply(wts*term2,c(2,3),sum)
#      data.term<-data.term-apply(delta*(temp2-2*wts*delta*term3)*dmats,c(2,3),sum)    
    }
    sigtheta.q.new<-sigtheta.q.new-0.5*e1oversig2*data.term
    sigma.inv.new<- -2*sigtheta.q.new
    sigtheta.q.new<-solve(sigma.inv.new)
    sigtheta.q.try<-(1-stepsize)*sigma.inv.old+stepsize*sigma.inv.new
    sigma.inv.temp<-sigtheta.q.try
    sigtheta.q.try<- solve(sigtheta.q.try)
    
    # Now update mean parameter mutheta.q

    dmatstheta<-array(t(apply(dmats,1,matfun5,mutheta.q)),dim=c(n.grid,T+1))
    term1<-array(t(apply(dmats,1,matfun4,sigtheta.q,mutheta.q)),dim=c(n.grid,T+1))
    term2<-apply(dmats,1,matfun3,sigtheta.q,mutheta.q)
    if(n<=n.grid) {
      data.term<-apply(as.vector(-4*delta*(y-Z%*%mubeta.q-delta*term2))*dmatstheta,2,sum)
    } else {
      temp1\<-aggregate(I(-4*delta*(y-Z%*%mubeta.q)) ~ index, FUN=sum)
      temp2<-rep(0,times=n.grid)
      temp2[temp1[,1]]<-temp1[,2]
      data.term<-apply((temp2+4*wts*term2)*dmatstheta,2,sum)    
    }
    mutheta.q.new<-mutheta.q+sigtheta.q.new%*%(-e1oversig*diag(e1overupsilon)%*%mutheta.q-0.5*e1oversig2*(apply(term1,2,sum)+data.term))    
    mutheta.q.new<-as.vector(mutheta.q.new)

    mutheta.q.try<-sigtheta.q.try%*%((1-stepsize)*sigma.inv.old%*%mutheta.q.old+stepsize*sigma.inv.new%*%mutheta.q.new)
    mutheta.q.try<-as.vector(mutheta.q.try)

    sigtheta.q<-sigtheta.q.try
    mutheta.q<-mutheta.q.try
    mutheta.q.old<-mutheta.q
    sigma.inv.old<-sigma.inv.temp

    # Compute expectations under MFVB update for 1/sigma and 1/sigma^2

    b<- -0.5*sum(diag((sigtheta.q+outer(mutheta.q,mutheta.q))%*%diag(e1overupsilon)))
    term2<-apply(dmats,1,matfun3,sigtheta.q,mutheta.q)
    if(n<=n.grid) {
      data.term<-sum(as.vector(y-Z%*%mubeta.q-delta*term2)^2)
    } else {
      data.term<-sum((y-Z%*%mubeta.q-delta*term2[index])^2)
    }
    cc<-0.5*(ssig.0+data.term+sum(apply(dmats,1,matfun6,sigtheta.q,mutheta.q))+sum(diag(ZtZ%*%sigbeta.q))+sum(diag(solve(sigbeta.0,sigbeta.q)))+sum((mubeta.q-mubeta.0)*solve(sigbeta.0,mubeta.q-mubeta.0)))

    mu1<-(b+sqrt(b^2+8*cc*(2*a-3)))/(4*cc)
    mu2<-(b+sqrt(b^2+8*cc*(2*a-2)))/(4*cc)
    mu3<-(b+sqrt(b^2+8*cc*(2*a-1)))/(4*cc)
    sig1<- -(2*a-3)/(mu1*mu1)-2*cc
    sig1<- 1/sqrt(-sig1)
    sig2<- -(2*a-2)/(mu2*mu2)-2*cc
    sig2<- 1/sqrt(-sig2)
    sig3<- -(2*a-1)/(mu3*mu3)-2*cc
    sig3<- 1/sqrt(-sig3)
    logI1<-log(2)-0.5*log(2*pi*sig1^2)+(2*a-3)*log(mu1)+b*mu1-cc*mu1*mu1-pnorm(mu1/sig1,log.p=TRUE)
    e1oversig<-exp(log(2)-0.5*log(2*pi*sig2^2)+(2*a-2)*log(mu2)+b*mu2-cc*mu2*mu2-pnorm(mu2/sig2,log.p=TRUE)-logI1)
    e1oversig2<-exp(log(2)-0.5*log(2*pi*sig3^2)+(2*a-1)*log(mu3)+b*mu3-cc*mu3*mu3-pnorm(mu3/sig3,log.p=TRUE)-logI1)

#    e1oversig<- -0.5*log(2*cc)+lgamma((2*a-1))-lgamma((2*a-2))
#    R1<-Rv(-b/sqrt(2*cc),(2*a-3))
#    e1oversig<-exp(e1oversig)*R1
#    e1oversig2<- -log(2*cc)+lgamma((2*a))-lgamma((2*a-2))
#    e1oversig2<- exp(e1oversig2)*R1*Rv(-b/sqrt(2*cc),(2*a-2))

    # Update variational parameters for tau

    DQvec<-diag(Qvec)
    stau.q<-stau.0+e1oversig*sum(diag((sigtheta.q[2:(T+1),2:(T+1)]+outer(mutheta.q[2:(T+1)],mutheta.q[2:(T+1)]))%*%DQvec))
    tau.ratio<-rtau.q/stau.q

    # Update variational parameters for beta

    sigbeta.q<-1/e1oversig2*sigbeta.term
    dmatstheta<-array(t(apply(dmats,1,matfun5,mutheta.q)),dim=c(n.grid,(T+1)))
    term1<-apply(dmats,1,matfun3,sigtheta.q,mutheta.q)
    if(n<=n.grid) {
       data.term<-apply(as.vector(y-delta*term1)*Z,2,sum)
    } else {
       data.term<-apply(as.vector(y-delta*term1[index])*Z,2,sum) 
    }
    mubeta.q<-as.vector(e1oversig2*sigbeta.q%*%(sb0ibymub0+data.term))

    # Update variational parameters for psi

    # Calculate lower bound before NCVMP udpate

    lbnew<-0
    sigpsi.q<-sqrt(sig2psi.q)
    mu2psi.q<-mupsi.q^2

    if(n<=n.grid) {
      data.term<-sum(as.vector(y-Z%*%mubeta.q-delta*term2)^2)
    } else {
      data.term<-sum((y-Z%*%mubeta.q-delta*term2[index])^2)
    }
    data.term<-data.term+sum(apply(dmats,1,matfun6,sigtheta.q,mutheta.q))+sum(diag(ZtZ%*%sigbeta.q))
    cc<-0.5*(ssig.0+data.term+sum((mubeta.q-mubeta.0)*solve(sigbeta.0,mubeta.q-mubeta.0))+sum(diag(solve(sigbeta.0,sigbeta.q))))
    lbnew<-lbnew-n/2*log(2*pi)-e1oversig2*data.term
  
    term1<-exp(sig2psi.q*bindices2/2+mupsi.q*bindices)
    term2<-exp(sig2psi.q*bindices2/2-mupsi.q*bindices)
    
    lbnew<-lbnew-p/2*log(2*pi)-0.5*ldetsigbeta.0-0.5*e1oversig2*(sum(diag(solve(sigbeta.0,sigbeta.q)))+sum((mubeta.q-mubeta.0)*solve(sigbeta.0,mubeta.q-mubeta.0)))
    lbnew<-lbnew-0.5*log(sig02)-T/2*(log(stau.q/2)-digamma(rtau.q/2))    
    lbnew<-lbnew+rsig.0/2*log(ssig.0/2)-lgamma(rsig.0/2)-ssig.0/2*e1oversig2
    lbnew<-lbnew+rtau.0/2*log(stau.0/2)-lgamma(rtau.0/2)-(rtau.0/2+1)*(log(stau.q/2)-digamma(rtau.q/2))-stau.0/2*rtau.q/stau.q
    lbnew<-lbnew-(-p/2*log(2*pi)-0.5*as.numeric(determinant(sigbeta.q,logarithm=TRUE)$modulus)-p/2)
    lbnew<-lbnew-(-T/2*log(2*pi)-0.5*as.numeric(determinant(sigtheta.q,logarithm=TRUE)$modulus)-T/2)
    lbnew<-lbnew-(rtau.q/2*log(stau.q/2)-lgamma(rtau.q/2)-(rtau.q/2+1)*(log(stau.q/2)-digamma(rtau.q/2))-rtau.q/2)
    lbnew.wopsibits<-lbnew
    lbnew<-lbnew-(-0.5*log(2*pi)-0.5*log(sig2psi.q)-0.5)
    S1<- -w0*(sigpsi.q*const1*exp(-mu2psi.q/(2*sig2psi.q))+mupsi.q*(1-2*pnorm(-mupsi.q/sigpsi.q)))
    lbnew<-lbnew+log(w0/2)+S1
    lbnew<-lbnew+T*(T+1)/4*(sigpsi.q*const1*exp(-mu2psi.q/(2*sig2psi.q))+mupsi.q*(1-2*pnorm(-mupsi.q/sigpsi.q)))
    Qvec<-term1*(1-pnorm(-mupsi.q/sigpsi.q-sigpsi.q*bindices))
    Qvec<-Qvec+term2*(1-pnorm(mupsi.q/sigpsi.q-sigpsi.q*bindices)) 
    e1overupsilon<-c(1/sig02,tau.ratio*Qvec)
    lbnew<-lbnew-0.5*e1oversig*sum(diag((sigtheta.q+outer(mutheta.q,mutheta.q))%*%diag(e1overupsilon)))
    b<- -0.5*sum(diag((sigtheta.q+outer(mutheta.q,mutheta.q))%*%diag(e1overupsilon)))
    lbnew<-lbnew+logI1-b*e1oversig+cc*e1oversig2    
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
    pd2sig<- mfactor*pd1sig-0.5*e1oversig*tau.ratio*sum((diag(sigtheta.q)+mutheta.q^2)*c(0,DQsig))
    pd2mu<- mfactor*pd1mu-0.5*e1oversig*tau.ratio*sum((diag(sigtheta.q)+mutheta.q^2)*c(0,DQmu))
    sig2psi.q.old<-sig2psi.q
    sig2psi.q<- -0.5/(pd1sig+pd2sig)
    #print(sig2psi.q)
    mupsi.q.old<-mupsi.q
    mupsi.q<-mupsi.q+sig2psi.q*(pd1mu+pd2mu)
    #print(mupsi.q)

    # Calculate lower bound value and then check that NCVMP step increases lower bound.
    # If not, then do step halving.

    sigpsi.q<-sqrt(sig2psi.q)
    mu2psi.q<-mupsi.q^2
    term1<-exp(sig2psi.q*bindices2/2+mupsi.q*bindices)
    term2<-exp(sig2psi.q*bindices2/2-mupsi.q*bindices)

    lbnew<-lbnew.wopsibits
    lbnew<-lbnew-(-0.5*log(2*pi)-0.5*log(sig2psi.q)-0.5)
    S1<- -w0*(sigpsi.q*const1*exp(-mu2psi.q/(2*sig2psi.q))+mupsi.q*(1-2*pnorm(-mupsi.q/sigpsi.q)))
    lbnew<-lbnew+log(w0/2)+S1
    lbnew<-lbnew+T*(T+1)/4*(sigpsi.q*const1*exp(-mu2psi.q/(2*sig2psi.q))+mupsi.q*(1-2*pnorm(-mupsi.q/sigpsi.q)))
    Qvec<-term1*(1-pnorm(-mupsi.q/sigpsi.q-sigpsi.q*bindices))
    Qvec<-Qvec+term2*(1-pnorm(mupsi.q/sigpsi.q-sigpsi.q*bindices)) 
    e1overupsilon<-c(1/sig02,tau.ratio*Qvec)
    lbnew<-lbnew-0.5*e1oversig*sum(diag((sigtheta.q+outer(mutheta.q,mutheta.q))%*%diag(e1overupsilon)))
    b<- -0.5*sum(diag((sigtheta.q+outer(mutheta.q,mutheta.q))%*%diag(e1overupsilon)))
    lbnew<-lbnew+logI1-b*e1oversig+cc*e1oversig2  
    
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
        lbnew<-lbnew+log(w0/2)+S1
        lbnew<-lbnew+T*(T+1)/4*(sigpsi.q.try*const1*exp(-mu2psi.q.try/(2*sig2psi.q.try))+mupsi.q.try*(1-2*pnorm(-mupsi.q.try/sigpsi.q.try)))
        Qvec<-term1*(1-pnorm(-mupsi.q.try/sigpsi.q.try-sigpsi.q.try*bindices))
        Qvec<-Qvec+term2*(1-pnorm(mupsi.q.try/sigpsi.q.try-sigpsi.q.try*bindices)) 
        e1overupsilon<-c(1/sig02,tau.ratio*Qvec)
        lbnew<-lbnew-0.5*e1oversig*sum(diag((sigtheta.q+outer(mutheta.q,mutheta.q))%*%diag(e1overupsilon)))
        b<- -0.5*sum(diag((sigtheta.q+outer(mutheta.q,mutheta.q))%*%diag(e1overupsilon)))
        lbnew<-lbnew+logI1-b*e1oversig+cc*e1oversig2  
        dif.try<-lbnew-lbfull
      }

      sigpsi.q<-sigpsi.q.try
      sig2psi.q<-sig2psi.q.try
      mupsi.q<-mupsi.q.try
      mu2psi.q<-mu2psi.q.try
      
    }    

    dif<-abs(lbnew-lbold)
    dif<-dif/abs(lbnew)
    lbold<-lbnew
    lb<-c(lb,lbnew)

  }

  # Return results 

  return(list(lb=lb,mutheta.q=mutheta.q,sigtheta.q=sigtheta.q,rtau.q=rtau.q,stau.q=stau.q,sigbeta.q=sigbeta.q,mubeta.q=mubeta.q,sig2psi.q=sig2psi.q,mupsi.q=mupsi.q,dmats=dmats,e1oversig=e1oversig,x.grid=x.grid))

}




