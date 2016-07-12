require('gpr') # for 'minimize'
require('stringr')

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
  if (any(x < 0 | x > 1)) {
    vphifull <- sqrt(2 * dnorm(x)) * cos(outer(pnorm(x), pi * (1:T)))
  } else {
    vphifull<-sqrt(2)*cos(outer(x,pi*(1:T)))
  }
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

minim <- function (X, f, .length, x, y)
{
    toprint = FALSE
    INT = 0.1
    EXT = 3
    MAX = 20
    RATIO = 10
    SIG = 0.1
    RHO = SIG/2
    if (is.array(.length)) {
        if (max(dim(.length)) == 2) {
            red = .length[2]
            .length = .length[1]
        }
    }
    else {
        red = 1
    }
    if (.length > 0) {
        S = "Linesearch"
    }
    else {
        S = "Function evaluation"
    }
    i.m = 0
    ls_failed = 0
    f_out = eval(call(f, X, x, y))
    f0 = c(f_out[1][[1]])
    df0 = c(f_out[2][[1]])
    fX = f0
    i.m = i.m + (.length < 0)
    s = -df0
    s = round(s * 100000)/100000
    d0 = -c(crossprod(s))
    x3 = red/(1 - d0)
    mainloop = TRUE
    while (i.m < abs(.length) && mainloop) {
        i.m = i.m + (.length > 0)
        X0 = X
        F0 = f0
        dF0 = df0
        if (.length > 0) {
            M = MAX
        }
        else {
            M = min(MAX, -.length - i.m)
        }
        whilerun = TRUE
        while (whilerun == TRUE) {
            x2 = 0
            f2 = f0
            d2 = d0
            f3 = f0
            df3 = df0
            success = FALSE
            while (success == FALSE && M > 0) {
                M = M - 1
                i.m = i.m + (.length < 0)
                options(show.error.messages = FALSE)
                f_out2 = eval(call(f, c(X + x3[1] * s), x, y))
                f3 = c(f_out2[1][[1]])
                df3 = c(f_out2[2][[1]])
                f3 = round(f3 * 100000)/100000
                df3 = round(df3 * 100000)/100000
                if (is.na(f3) || is.infinite(f3) || is.nan(f3) ||
                  any(is.nan(df3) || is.na(df3) || is.infinite(df3))) {
                  cat(" ")
                  x3 = (x2 + x3)/2
                }
                else {
                  success = TRUE
                }
                options(show.error.messages = TRUE)
            }
            if (f3 < F0) {
                X0 = X + x3[1] * s
                F0 = f3
                dF0 = df3
            }
            d3 = c(crossprod(df3, s))
            if (d3 > SIG * d0 || f3 > f0 + x3 * RHO * d0 || M ==
                0) {
                whilerun = FALSE
                break
            }
            x1 = x2
            f1 = f2
            d1 = d2
            x2 = x3
            f2 = f3
            d2 = d3
            A = 6 * (f1 - f2) + 3 * (d2 + d1) * (x2 - x1)
            B = 3 * (f2 - f1) - (2 * d1 + d2) * (x2 - x1)
            x3 = x1 - d1 * (x2 - x1)^2/(B + sqrt(abs(B * B -
                A * d1 * (x2 - x1))))
            if ((B * B - A * d1 * (x2 - x1) < 0)[1] || is.nan(x3) ||
                is.infinite(x3) || x3 < 0) {
                x3 = x2 * EXT
            }
            else if (x3 > x2 * EXT) {
                x3 = x2 * EXT
            }
            else if (x3 < x2 + INT * (x2 - x1)) {
                x3 = x2 + INT * (x2 - x1)
            }
            x3 = round(x3 * 100000)/100000
        }
        while ((abs(d3) > -SIG * d0 || f3 > f0 + x3 * RHO * d0) &&
            M > 0) {
            if (d3 > 0 || f3 > f0 + x3 * RHO * d0) {
                x4 = x3
                f4 = f3
                d4 = d3
            }
            else {
                x2 = x3
                f2 = f3
                d2 = d3
            }
            if (f4 > f0) {
                x3 = x2 - (0.5 * d2 * (x4 - x2)^2)/(f4 - f2 -
                  d2 * (x4 - x2))
            }
            else {
                A = 6 * (f2 - f4)/(x4 - x2) + 3 * (d4 + d2)
                B = 3 * (f4 - f2) - (2 * d2 + d4) * (x4 - x2)
                x3 = x2 + (sqrt(B * B - A * d2 * (x4 - x2)^2) -
                  B)/A
            }
            if (is.nan(x3) || is.infinite(x3)) {
                x3 = (x2 + x4)/2
            }
            x3 = max(min(x3, x4 - INT * (x4 - x2)), x2 + INT *
                (x4 - x2))
            f_out3 = eval(call(f, c(X + x3 * s), x, y))
            f3 = c(f_out3[1][[1]])
            df3 = c(f_out3[2][[1]])
            if (f3 < F0) {
                x3 = x3[[1]]
                X0 = X + x3 * s
                F0 = f3
                dF0 = df3
            }
            M = M - 1
            i.m = i.m + (.length < 0)
            d3 = c(crossprod(df3, s))
        }
        if (abs(d3) < -SIG * d0 && f3 < f0 + x3 * RHO * d0) {
            x3 = x3[[1]]
            X = X + x3 * s
            f0 = f3
            fX = c(fX, f0)
            cat(S, i.m, "; Value ", f0, '\n')
            s = (((c(crossprod(df3)) - c(crossprod(df0, df3)))[1])/((c(crossprod(df0)))[[1]]) * s) - df3
            df0 = df3
            d3 = d0
            d0 = c(crossprod(df0, s))
            if (d0 > 0) {
                s = -df0
                d0 = -c(crossprod(s))
            }
            x3 = x3 * min(RATIO, d3/(d0 - (2^(-1022))))
            ls_failed = 0
        }
        else {
            X = X0
            f0 = F0
            df0 = dF0
            if (ls_failed || i.m > abs(.length)) {
                mainloop = 0
                break
            }
            s = -df0
            d0 = -c(crossprod(s))
            x3 = 1/(1 - d0)
            ls_failed = 1
        }
    }
    return(list(X=X, fX=fX, i.m=i.m))
}


ssgpr <- function(optimizeparams, x_tr, y_tr, x_tst = NULL) {
  n <- dim(x_tr)[1]
  D <- dim(x_tr)[2]
  m <- (length(optimizeparams) - D - 2) / D
  ell <- exp(optimizeparams[1:D])
  sf2 <- exp(2 * optimizeparams[D+1])
  sn2 <- exp(2 * optimizeparams[D+2])
  w <- matrix(c(optimizeparams[(D+3):length(optimizeparams)]), nrow = m, ncol = D)
  # cat('dim of w: ', dim(w), '\n')
  # cat('diag(1/ell): ', diag(ell), '\n')
  # cat('ell: ', ell, '\n')
  if (length(ell) == 1) {
    w <- w / ell # dividing each row of w by ell element-wise... needs to be revised if dimension of w changes
  } else {
    w <- w %*% diag(1 / ell) # divide each row of w by ell element-wise
  }

  phi <- tcrossprod(x_tr, w)
  phi <- cbind(cos(phi), sin(phi))
  
  R <- chol((sf2/m) * crossprod(phi) + sn2 * diag(2 * m))
  PhiRi <- phi %*% solve(R)
  RtiPhit <- t(PhiRi)
  Rtiphity <- RtiPhit %*% y_tr
  
  if (is.null(x_tst)) {
    out1 <- 0.5 / sn2 * (sum(y_tr^2) - sf2 / m * sum(Rtiphity^2)) + sum(log(diag(R))) + (n / 2 - m) * log(sn2) + n / 2 + log(2 * pi)
    
    out2 <- rep(0, D+2+D*m)
    
    A <- cbind(y_tr / sn2 - PhiRi %*% ((sf2 / sn2 / m) * Rtiphity), sqrt(sf2 / sn2 / m) * PhiRi)
    diagfact <- -1 / sn2 + rowSums(A^2)
    Aphi <- crossprod(A, phi)
    B <- ((A %*% Aphi[,1:m]) * phi[, (m+1):ncol(phi)]) - ((A %*% Aphi[, (m+1):ncol(Aphi)]) * phi[,1:m])
    
    for (d in 1:D) {
      out2[d] <- -0.5 * 2 * sf2 / m * (crossprod(x_tr[,d], B) %*% w[,d])
    }
    out2[D+1] <- 0.5 * 2 * (sf2 / m) * (n * m / sn2 - sum(Aphi^2))
    out2[D+2] <- -0.5 * sum(diagfact) * 2 * sn2
    
    for (d in 1:D) {
      out2[(D+2+(d-1)*m+1):(D+2+d*m)] <- 0.5 * 2 * sf2 / m * (c(crossprod(x_tr[,d], B)) / ell[d])
    }
  } else {
    ns <- dim(x_tst)[1]
    out1 <- out2 <- rep(0, ns)
    alfa <- sf2 / m * (solve(R, Rtiphity))
    
    chunksize <- 5000
    allxstar <- x_tst
    
    for (beg_chunk in seq(from = 1, to = ns, by = chunksize)) {
      end_chunk <- min(beg_chunk + chunksize - 1, ns)
      x_tst <- allxstar[beg_chunk:end_chunk,]
      
      phistar <- tcrossprod(x_tst, w)
      phistar <- cbind(cos(phistar), sin(phistar))
      out1[beg_chunk:end_chunk] <- phistar %*% alfa
      
      out2[beg_chunk:end_chunk] <- sn2 * (1 + sf2 / m * rowSums((phistar %*% solve(R))^2))
    }
    
  }
  list(out1 = out1, out2 = out2)
}


ssgpr_ui <- function(x_tr, y_tr, x_tst, y_tst, m, iteropt = NULL, loghyper = NULL) {
  meanp <- mean(apply(y_tr, 2, mean))
  y_tr <- scale(y_tr, scale = FALSE)
  n <- dim(x_tr)[1]
  D <- dim(x_tr)[2]
  if ((!is.null(loghyper)) & (length(loghyper) == (D+2))) {
    lengthscales <- loghyper[1:D]
    covpower <- loghyper[D+1]
    noisepower <- loghyper[D+2]
    nlml <- Inf
    optimizeparams <- c(lengthscales, covpower, noisepower)
    for (k in 1:100) {
      otherparams <- rnorm(m * D)
      nlmlc <- ssgpr(c(optimizeparams, otherparams), x_tr, y_tr)$out1
      if (nlmlc < nlml) {
        w_save <- otherparams
        nlml <- nlmlc
      }
    }
    otherparams <- w_save
  } else if (is.null(loghyper)) {
    lengthscales <- log((apply(x_tr, 2, max) - apply(x_tr, 2, min)) / 2)
    lengthscales[lengthscales < -1e2] <- -1e2
    covpower <- 0.5 * log(apply(y_tr, 2, function(x) mean((x-mean(x))^2)))
    noisepower <- 0.5 * log(apply(y_tr, 2, function(x) mean((x-mean(x))^2)) / 4)
    nlml <- Inf
    optimizeparams <- c(lengthscales, covpower, noisepower)
    for (k in 1:100) {
      otherparams <- rnorm(m * D)
      nlmlc <- ssgpr(c(optimizeparams, otherparams), x_tr, y_tr)$out1
      if (nlmlc < nlml) {
        w_save <- otherparams
        nlml <- nlmlc
      }
    }
    otherparams <- w_save
  } else if ((!is.null(loghyper)) & (length(loghyper) != D+2+D*m)) {
    stop('Incorrect number of hyperparameters.')
  } else if ((!is.null(loghyper)) & (length(loghyper) == D+2+D*m)) {
    optimizeparams <- loghyper[1:(D+2)]
    otherparams <- loghyper[(D+3):(D+2+D*m)]
  }
  if (is.null(iteropt)) iteropt <- -1000
  
  optimizeparams <- c(optimizeparams, otherparams)
  temp <- minim(optimizeparams, 'ssgpr', iteropt, x_tr, y_tr)

  optimizeparams <- temp$X
  convergence <- temp$fX
  loghyper <- optimizeparams

  res <- ssgpr(optimizeparams, x_tr, y_tr, x_tst)

  mu <- res$out1
  S2 <- res$out2
  mu <- mu + meanp  
  NMSE <- mean((mu-y_tst)^2) / mean((meanp - y_tst)^2)
  
  NMLP <- -0.5 * mean((-(mu - y_tst)^2) / S2 - log(2 * pi) - log(S2))
  list(NMSE = NMSE, mu = mu, S2 = S2, NMLP = NMLP, loghyper = loghyper, convergence = convergence)
}

compareSSGPvsBSAR <- function(data = 'pendulum', fit = 'training', path = NULL, fileName_X_tr = NULL, fileName_T_tr = NULL, fileName_X_tst = NULL, fileName_T_tst = NULL) {
  ################################################################################################
  ################################################################################################
  ############     data: which data,                                                  ############
  ############       c('pendulum', 'elevators', 'kin', 'pol', 'pumadyn', 'simul')     ############
  ############                                                                        ############
  ############     fit: return training MSE or test MSE?                              ############
  ############       c('training', 'test')                                            ############
  ############                                                                        ############
  ############     path: path to the data folder                                      ############
  ############           needs to be coordinated with 'data' variable                 ############
  ############                                                                        ############
  ############     file_name_*: the name of each data file                            ############
  ############                  if fit == 'training', test files could be set to NULL ############
  ################################################################################################
  ################################################################################################
  
  # if (data != 'simul') {
  #   if (pat)
  # }



  if (data != 'simul' & any(is.null(path), is.null(fileName_X_tr), is.null(fileName_T_tr))) {
    stop("Real data fitting must be provided with 'path' and 'file_Name_*' variables.")
  }
  if (data != 'simul' & !is.null(path)) {
    if (str_sub(path, -1, -1) != '/') stop('Path variable must end with a slash indicating it is a directory.')
  }
  if (data != 'simul') {
    assign("X_tr", unname(as.matrix(read.table(paste0(path, fileName_X_tr)))), envir = .GlobalEnv)
    assign("T_tr", unname(as.matrix(read.table(paste0(path, fileName_T_tr)))), envir = .GlobalEnv)
    if (fit == 'test') {
      if (any(is.null(fileName_X_tst), missing(fileName_X_tst), is.null(fileName_T_tst), missing(fileName_T_tst))) {
        stop('You should provide the file name of the test sets.')
      } else {
        assign("X_tst", unname(as.matrix(read.table(paste0(path, fileName_X_tst)))), envir = .GlobalEnv)
        assign("T_tst", unname(as.matrix(read.table(paste0(path, fileName_T_tst)))), envir = .GlobalEnv)
      }
    }
  }


  call_SSGPR_tr  <- call('ssgpr_ui', quote(X_tr), quote(T_tr), quote(X_tr), quote(T_tr), 100, -100, quote(rep(1,ncol(X_tr)+2)))
  call_SSGPR_tst <- call('ssgpr_ui', quote(X_tr), quote(T_tr), quote(X_tst), quote(T_tst), 100, -100, quote(rep(1,ncol(X_tr)+2)))
  call_BSAR      <- call('vbgpspectral', quote(y), quote(x), quote(x_rest), quote(J), quote(tol), quote(prior.parms), 1)
  
  q_SSE_SSGP_tr <- quote( SSE_SSGP_tr <- round(mean((T_tr - res_SSGP$mu)^2), digits = 4) )
  q_SSE_BSAR_tr <- quote( SSE_BSAR_tr <- round(mean((T_tr - mu_BSAR)^2), digits = 4) )

  q_NMSE_SSGP   <- quote( NMSE_SSGP   <- round(mean((T_tst - res_SSGP$mu)^2) / mean((T_tst - mean(T_tr))^2), digits = 4) )
  q_NMSE_BSAR   <- quote( NMSE_BSAR   <- round(mean((T_tst - mu_BSAR)^2) / mean((T_tst - mean(T_tr))^2), digits = 4) )

  J      <- 20
  rsig.0 <- 0.01
  ssig.0 <- 0.01
  rtau.0 <- 0.01
  stau.0 <- 0.01
  w0     <- 1
  tol    <- 1.0e-05

  q_mubeta.0     <- quote( mubeta.0    <- rep(0, times = ncol(x_rest)) )
  q_sigbeta.0    <- quote( sigbeta.0   <- diag(1, ncol(x_rest)) )
  q_prior.params <- quote( prior.parms <- list(rsig.0=rsig.0,ssig.0=ssig.0,rtau.0=rtau.0,stau.0=stau.0,w0=w0,mubeta.0=mubeta.0,sigbeta.0=sigbeta.0) )


  switch(data,
               'pendulum' = {
                              if (fit == 'training') {
                                res_SSGP      <- eval(call_SSGPR_tr)
                                centred_SSGP  <- res_SSGP$mu - mean(res_SSGP$mu)

                                ind           <- which.min(cov(T_tr, X_tr))
                                x             <- X_tr[,ind]

                                # # squish x into [0,1]
                                # if (!all(x < 1 & x > 0)) {
                                #   x <- pnorm(x)
                                # }

                                x_rest        <- X_tr[,-ind]
                                y             <- c(T_tr)
                                y_centred     <- y - mean(y)

                                eval(q_mubeta.0)
                                eval(q_sigbeta.0)
                                eval(q_prior.params)
                                
                                res_BSAR      <- eval(call_BSAR)
                                vphi          <- sqrt(2)*cos(outer(x,pi*(1:J)))
                                mu_BSAR       <- vphi[,1:length(res_BSAR$mutheta.q)] %*% res_BSAR$mutheta.q + x_rest %*% res_BSAR$mubeta.q
                                centred_BSAR  <- mu_BSAR - mean(mu_BSAR)
                                eval(q_SSE_SSGP_tr)
                                eval(q_SSE_BSAR_tr)

                                o             <- order(x)
                                plot(x[o], y_centred[o], xlab = '', ylab = 'fitted', pch = 20, main = 'Pendulum data / training')
                                lines(x[o], centred_SSGP[o], lwd = 2, lty = 2, col = 'red')
                                lines(x[o], centred_BSAR[o], lwd = 2, lty = 3, col = 'darkgreen')
                                legend('topright', lty = c(NA, 2, 3), pch = c(20, NA, NA), col = c(1, 'red', 'darkgreen'), legend = c('true', paste0('SSGP: ', SSE_SSGP_tr, ' (SSE)'), paste0('BSAR ', SSE_BSAR_tr, ' (SSE)')))
                                return(list(res_SSGP = res_SSGP, res_BSAR = res_BSAR, mu_BSAR = mu_BSAR, mu_SSGP = res_SSGP$mu, centred_SSGP = centred_SSGP, centred_BSAR = centred_BSAR))
                              } else if (fit == 'test') {
                                res_SSGP     <- eval(call_SSGPR_tst)
                                centred_SSGP <- res_SSGP$mu - mean(res_SSGP$mu)

                                ind          <- which.min(cov(T_tr, X_tr))
                                x            <- X_tr[,ind]

                                # # squish x into [0,1]
                                # if (!all(x < 1 & x > 0)) {
                                #   x <- pnorm(x)
                                # }

                                x_rest       <- X_tr[,-ind]
                                y            <- c(T_tr)
                                y_centred    <- y - mean(y)

                                eval(q_mubeta.0)
                                eval(q_sigbeta.0)
                                eval(q_prior.params)

                                res_BSAR     <- eval(call_BSAR)
                                vphi         <- sqrt(2)*cos(outer(X_tst[,ind],pi*(1:J)))
                                mu_BSAR      <- vphi[,1:length(res_BSAR$mutheta.q)]%*%res_BSAR$mutheta.q + X_tst[,-ind] %*% res_BSAR$mubeta.q
                                centred_BSAR <- mu_BSAR - mean(mu_BSAR)
                                eval(q_NMSE_SSGP)
                                eval(q_NMSE_BSAR)

                                o            <- order(x)
                                plot(x[o], y_centred[o], xlab = '', ylab = 'fitted', pch = 20, main = 'Pendulum data / test')
                                lines(x[o], centred_SSGP[o], lwd = 2, lty = 2, col = 'red')
                                lines(x[o], centred_BSAR[o], lwd = 2, lty = 3, col = 'darkgreen')
                                legend('topright', lty = c(NA, 2, 3), pch = c(20, NA, NA), col = c(1, 'red', 'darkgreen'), legend = c('true', paste0('SSGP: ', NMSE_SSGP, ' (NMSE)'), paste0('BSAR: ', NMSE_BSAR, ' (NMSE)')))
                                return(list(res_SSGP = res_SSGP, res_BSAR = res_BSAR, mu_BSAR = mu_BSAR, mu_SSGP = res_SSGP$mu, centred_SSGP = centred_SSGP, centred_BSAR = centred_BSAR))
                              } else {
                                stop("You've inserted the wrong option. Please pick between 'training' and 'test'.")
                              }
               },
               'elevators' = {

                              if (fit == 'training') {
                                res_SSGP      <- eval(call_SSGPR_tr)
                                centred_SSGP  <- res_SSGP$mu - mean(res_SSGP$mu)

                                ind           <- which.min(cov(T_tr, X_tr))
                                x             <- X_tr[,ind]

                                # # squish x into [0,1]
                                # if (!all(x < 1 & x > 0)) {
                                #   x <- pnorm(x)
                                # }

                                x_rest        <- X_tr[,-ind]
                                y             <- c(T_tr)
                                y_centred     <- y - mean(y)
                                
                                eval(q_mubeta.0)
                                eval(q_sigbeta.0)
                                eval(q_prior.params)
                                
                                res_BSAR      <- eval(call_BSAR)
                                vphi          <- sqrt(2)*cos(outer(x,pi*(1:J)))
                                mu_BSAR       <- vphi[,1:length(res_BSAR$mutheta.q)] %*% res_BSAR$mutheta.q + x_rest %*% res_BSAR$mubeta.q
                                centred_BSAR  <- mu_BSAR - mean(mu_BSAR)
                                eval(q_SSE_SSGP_tr)
                                eval(q_SSE_BSAR_tr)

                                o             <- order(x)
                                plot(x[o], y_centred[o], xlab = '', ylab = 'fitted', pch = 20, main = 'Elevators data / training')
                                lines(x[o], centred_SSGP[o], lwd = 2, lty = 2, col = 'red')
                                lines(x[o], centred_BSAR[o], lwd = 2, lty = 3, col = 'darkgreen')
                                legend('topright', lty = c(NA, 2, 3), pch = c(20, NA, NA), col = c(1, 'red', 'darkgreen'), legend = c('true', paste0('SSGP: ', SSE_SSGP_tr, ' (SSE)'), paste0('BSAR ', SSE_BSAR_tr, ' (SSE)')))
                                return(list(res_SSGP = res_SSGP, res_BSAR = res_BSAR, mu_BSAR = mu_BSAR, mu_SSGP = res_SSGP$mu, centred_SSGP = centred_SSGP, centred_BSAR = centred_BSAR))
                              } else if (fit == 'test') {
                                res_SSGP     <- eval(call_SSGPR_tst)
                                centred_SSGP <- res_SSGP$mu - mean(res_SSGP$mu)

                                ind          <- which.min(cov(T_tr, X_tr))
                                x            <- X_tr[,ind]

                                # # squish x into [0,1]
                                # if (!all(x < 1 & x > 0)) {
                                #   x <- pnorm(x)
                                # }
                                
                                x_rest       <- X_tr[,-ind]
                                y            <- c(T_tr)
                                y_centred    <- y - mean(y)

                                eval(q_mubeta.0)
                                eval(q_sigbeta.0)
                                eval(q_prior.params)

                                res_BSAR     <- eval(call_BSAR)
                                vphi         <- sqrt(2)*cos(outer(X_tst[,ind],pi*(1:J)))
                                mu_BSAR      <- vphi[,1:length(res_BSAR$mutheta.q)]%*%res_BSAR$mutheta.q + X_tst[,-ind] %*% res_BSAR$mubeta.q
                                centred_BSAR <- mu_BSAR - mean(mu_BSAR)
                                eval(q_NMSE_SSGP)
                                eval(q_NMSE_BSAR)

                                o            <- order(x)
                                plot(x[o], y_centred[o], xlab = '', ylab = 'fitted', pch = 20, main = 'Elevators data / test')
                                lines(x[o], centred_SSGP[o], lwd = 2, lty = 2, col = 'red')
                                lines(x[o], centred_BSAR[o], lwd = 2, lty = 3, col = 'darkgreen')
                                legend('topright', lty = c(NA, 2, 3), pch = c(20, NA, NA), col = c(1, 'red', 'darkgreen'), legend = c('true', paste0('SSGP: ', NMSE_SSGP, ' (NMSE)'), paste0('BSAR: ', NMSE_BSAR, ' (NMSE)')))
                                return(list(res_SSGP = res_SSGP, res_BSAR = res_BSAR, mu_BSAR = mu_BSAR, mu_SSGP = res_SSGP$mu, centred_SSGP = centred_SSGP, centred_BSAR = centred_BSAR))
                              } else {
                                stop("You've inserted the wrong option. Please pick between 'training' and 'test'.")
                              }
               },
               'kin'       = {

                              if (fit == 'training') {
                                res_SSGP      <- eval(call_SSGPR_tr)
                                centred_SSGP  <- res_SSGP$mu - mean(res_SSGP$mu)

                                ind           <- which.min(cov(T_tr, X_tr))
                                x             <- X_tr[,ind]

                                # # squish x into [0,1]
                                # if (!all(x < 1 & x > 0)) {
                                #   x <- pnorm(x)
                                # }

                                x_rest        <- X_tr[,-ind]
                                y             <- c(T_tr)
                                y_centred     <- y - mean(y)
                                
                                eval(q_mubeta.0)
                                eval(q_sigbeta.0)
                                eval(q_prior.params)
                                
                                res_BSAR      <- eval(call_BSAR)
                                vphi          <- sqrt(2)*cos(outer(x,pi*(1:J)))
                                mu_BSAR       <- vphi[,1:length(res_BSAR$mutheta.q)] %*% res_BSAR$mutheta.q + x_rest %*% res_BSAR$mubeta.q
                                centred_BSAR  <- mu_BSAR - mean(mu_BSAR)
                                eval(q_SSE_SSGP_tr)
                                eval(q_SSE_BSAR_tr)

                                o             <- order(x)
                                plot(x[o], y_centred[o], xlab = '', ylab = 'fitted', pch = 20, main = 'Kin data / training')
                                lines(x[o], centred_SSGP[o], lwd = 2, lty = 2, col = 'red')
                                lines(x[o], centred_BSAR[o], lwd = 2, lty = 3, col = 'darkgreen')
                                legend('topright', lty = c(NA, 2, 3), pch = c(20, NA, NA), col = c(1, 'red', 'darkgreen'), legend = c('true', paste0('SSGP: ', SSE_SSGP_tr, ' (SSE)'), paste0('BSAR ', SSE_BSAR_tr, ' (SSE)')))
                                return(list(res_SSGP = res_SSGP, res_BSAR = res_BSAR, mu_BSAR = mu_BSAR, mu_SSGP = res_SSGP$mu, centred_SSGP = centred_SSGP, centred_BSAR = centred_BSAR))
                              } else if (fit == 'test') {
                                res_SSGP     <- eval(call_SSGPR_tst)
                                centred_SSGP <- res_SSGP$mu - mean(res_SSGP$mu)

                                ind          <- which.min(cov(T_tr, X_tr))
                                x            <- X_tr[,ind]

                                # # squish x into [0,1]
                                # if (!all(x < 1 & x > 0)) {
                                #   x <- pnorm(x)
                                # }
                                
                                x_rest       <- X_tr[,-ind]
                                y            <- c(T_tr)
                                y_centred    <- y - mean(y)

                                eval(q_mubeta.0)
                                eval(q_sigbeta.0)
                                eval(q_prior.params)

                                res_BSAR     <- eval(call_BSAR)
                                vphi         <- sqrt(2)*cos(outer(X_tst[,ind],pi*(1:J)))
                                mu_BSAR      <- vphi[,1:length(res_BSAR$mutheta.q)]%*%res_BSAR$mutheta.q + X_tst[,-ind] %*% res_BSAR$mubeta.q
                                centred_BSAR <- mu_BSAR - mean(mu_BSAR)
                                eval(q_NMSE_SSGP)
                                eval(q_NMSE_BSAR)

                                o            <- order(x)
                                plot(x[o], y_centred[o], xlab = '', ylab = 'fitted', pch = 20, main = 'Kin data / test')
                                lines(x[o], centred_SSGP[o], lwd = 2, lty = 2, col = 'red')
                                lines(x[o], centred_BSAR[o], lwd = 2, lty = 3, col = 'darkgreen')
                                legend('topright', lty = c(NA, 2, 3), pch = c(20, NA, NA), col = c(1, 'red', 'darkgreen'), legend = c('true', paste0('SSGP: ', NMSE_SSGP, ' (NMSE)'), paste0('BSAR: ', NMSE_BSAR, ' (NMSE)')))
                                return(list(res_SSGP = res_SSGP, res_BSAR = res_BSAR, mu_BSAR = mu_BSAR, mu_SSGP = res_SSGP$mu, centred_SSGP = centred_SSGP, centred_BSAR = centred_BSAR))
                              } else {
                                stop("You've inserted the wrong option. Please pick between 'training' and 'test'.")
                              }
               },
               'pol'       = {

                              if (fit == 'training') {
                                res_SSGP      <- eval(call_SSGPR_tr)
                                centred_SSGP  <- res_SSGP$mu - mean(res_SSGP$mu)

                                ind           <- which.min(cov(T_tr, X_tr))
                                x             <- X_tr[,ind]

                                # # squish x into [0,1]
                                # if (!all(x < 1 & x > 0)) {
                                #   x <- pnorm(x)
                                # }

                                x_rest        <- X_tr[,-ind]
                                y             <- c(T_tr)
                                y_centred     <- y - mean(y)
                                
                                eval(q_mubeta.0)
                                eval(q_sigbeta.0)
                                eval(q_prior.params)
                                
                                res_BSAR      <- eval(call_BSAR)
                                vphi          <- sqrt(2)*cos(outer(x,pi*(1:J)))
                                mu_BSAR       <- vphi[,1:length(res_BSAR$mutheta.q)] %*% res_BSAR$mutheta.q + x_rest %*% res_BSAR$mubeta.q
                                centred_BSAR  <- mu_BSAR - mean(mu_BSAR)
                                eval(q_SSE_SSGP_tr)
                                eval(q_SSE_BSAR_tr)

                                o             <- order(x)
                                plot(x[o], y_centred[o], xlab = '', ylab = 'fitted', pch = 20, main = 'Pol data / training')
                                lines(x[o], centred_SSGP[o], lwd = 2, lty = 2, col = 'red')
                                lines(x[o], centred_BSAR[o], lwd = 2, lty = 3, col = 'darkgreen')
                                legend('topright', lty = c(NA, 2, 3), pch = c(20, NA, NA), col = c(1, 'red', 'darkgreen'), legend = c('true', paste0('SSGP: ', SSE_SSGP_tr, ' (SSE)'), paste0('BSAR ', SSE_BSAR_tr, ' (SSE)')))
                                return(list(res_SSGP = res_SSGP, res_BSAR = res_BSAR, mu_BSAR = mu_BSAR, mu_SSGP = res_SSGP$mu, centred_SSGP = centred_SSGP, centred_BSAR = centred_BSAR))
                              } else if (fit == 'test') {
                                res_SSGP     <- eval(call_SSGPR_tst)
                                centred_SSGP <- res_SSGP$mu - mean(res_SSGP$mu)

                                ind          <- which.min(cov(T_tr, X_tr))
                                x            <- X_tr[,ind]

                                # # squish x into [0,1]
                                # if (!all(x < 1 & x > 0)) {
                                #   x <- pnorm(x)
                                # }
                                
                                x_rest       <- X_tr[,-ind]
                                y            <- c(T_tr)
                                y_centred    <- y - mean(y)

                                eval(q_mubeta.0)
                                eval(q_sigbeta.0)
                                eval(q_prior.params)

                                res_BSAR     <- eval(call_BSAR)
                                vphi         <- sqrt(2)*cos(outer(X_tst[,ind],pi*(1:J)))
                                mu_BSAR      <- vphi[,1:length(res_BSAR$mutheta.q)]%*%res_BSAR$mutheta.q + X_tst[,-ind] %*% res_BSAR$mubeta.q
                                centred_BSAR <- mu_BSAR - mean(mu_BSAR)
                                eval(q_NMSE_SSGP)
                                eval(q_NMSE_BSAR)

                                o            <- order(x)
                                plot(x[o], y_centred[o], xlab = '', ylab = 'fitted', pch = 20, main = 'Pol data / test')
                                lines(x[o], centred_SSGP[o], lwd = 2, lty = 2, col = 'red')
                                lines(x[o], centred_BSAR[o], lwd = 2, lty = 3, col = 'darkgreen')
                                legend('topright', lty = c(NA, 2, 3), pch = c(20, NA, NA), col = c(1, 'red', 'darkgreen'), legend = c('true', paste0('SSGP: ', NMSE_SSGP, ' (NMSE)'), paste0('BSAR: ', NMSE_BSAR, ' (NMSE)')))
                                return(list(res_SSGP = res_SSGP, res_BSAR = res_BSAR, mu_BSAR = mu_BSAR, mu_SSGP = res_SSGP$mu, centred_SSGP = centred_SSGP, centred_BSAR = centred_BSAR))
                              } else {
                                stop("You've inserted the wrong option. Please pick between 'training' and 'test'.")
                              }
               },
               'pumadyn'   = {

                              if (fit == 'training') {
                                res_SSGP      <- eval(call_SSGPR_tr)
                                centred_SSGP  <- res_SSGP$mu - mean(res_SSGP$mu)

                                ind           <- which.min(cov(T_tr, X_tr))
                                x             <- X_tr[,ind]

                                # # squish x into [0,1]
                                # if (!all(x < 1 & x > 0)) {
                                #   x <- pnorm(x)
                                # }

                                x_rest        <- X_tr[,-ind]
                                y             <- c(T_tr)
                                y_centred     <- y - mean(y)
                                
                                eval(q_mubeta.0)
                                eval(q_sigbeta.0)
                                eval(q_prior.params)
                                
                                res_BSAR      <- eval(call_BSAR)
                                vphi          <- sqrt(2)*cos(outer(x,pi*(1:J)))
                                mu_BSAR       <- vphi[,1:length(res_BSAR$mutheta.q)] %*% res_BSAR$mutheta.q + x_rest %*% res_BSAR$mubeta.q
                                centred_BSAR  <- mu_BSAR - mean(mu_BSAR)
                                eval(q_SSE_SSGP_tr)
                                eval(q_SSE_BSAR_tr)

                                o             <- order(x)
                                plot(x[o], y_centred[o], xlab = '', ylab = 'fitted', pch = 20, main = 'Pumadyn data / training')
                                lines(x[o], centred_SSGP[o], lwd = 2, lty = 2, col = 'red')
                                lines(x[o], centred_BSAR[o], lwd = 2, lty = 3, col = 'darkgreen')
                                legend('topright', lty = c(NA, 2, 3), pch = c(20, NA, NA), col = c(1, 'red', 'darkgreen'), legend = c('true', paste0('SSGP: ', SSE_SSGP_tr, ' (SSE)'), paste0('BSAR ', SSE_BSAR_tr, ' (SSE)')))
                                return(list(res_SSGP = res_SSGP, res_BSAR = res_BSAR, mu_BSAR = mu_BSAR, mu_SSGP = res_SSGP$mu, centred_SSGP = centred_SSGP, centred_BSAR = centred_BSAR))
                              } else if (fit == 'test') {
                                res_SSGP     <- eval(call_SSGPR_tst)
                                centred_SSGP <- res_SSGP$mu - mean(res_SSGP$mu)

                                ind          <- which.min(cov(T_tr, X_tr))
                                x            <- X_tr[,ind]

                                # # squish x into [0,1]
                                # if (!all(x < 1 & x > 0)) {
                                #   x <- pnorm(x)
                                # }
                                
                                x_rest       <- X_tr[,-ind]
                                y            <- c(T_tr)
                                y_centred    <- y - mean(y)

                                eval(q_mubeta.0)
                                eval(q_sigbeta.0)
                                eval(q_prior.params)

                                res_BSAR     <- eval(call_BSAR)
                                vphi         <- sqrt(2)*cos(outer(X_tst[,ind],pi*(1:J)))
                                mu_BSAR      <- vphi[,1:length(res_BSAR$mutheta.q)]%*%res_BSAR$mutheta.q + X_tst[,-ind] %*% res_BSAR$mubeta.q
                                centred_BSAR <- mu_BSAR - mean(mu_BSAR)
                                eval(q_NMSE_SSGP)
                                eval(q_NMSE_BSAR)

                                o            <- order(x)
                                plot(x[o], y_centred[o], xlab = '', ylab = 'fitted', pch = 20, col = rgb(0, 0, 0, 0.3), main = 'Pumadyn data / test')
                                lines(x[o], centred_SSGP[o], lwd = 2, lty = 2, col = 'red')
                                lines(x[o], centred_BSAR[o], lwd = 2, lty = 3, col = 'darkgreen')
                                legend('topright', lty = c(NA, 2, 3), pch = c(20, NA, NA), col = c(rgb(0, 0, 0, 0.3), 'red', 'darkgreen'), legend = c('true', paste0('SSGP: ', NMSE_SSGP, ' (NMSE)'), paste0('BSAR: ', NMSE_BSAR, ' (NMSE)')))
                                return(list(res_SSGP = res_SSGP, res_BSAR = res_BSAR, mu_BSAR = mu_BSAR, mu_SSGP = res_SSGP$mu, centred_SSGP = centred_SSGP, centred_BSAR = centred_BSAR))
                              } else {
                                stop("You've inserted the wrong option. Please pick between 'training' and 'test'.")
                              }
               },
               'simul' = {
                          e         <- readline(prompt = "Define a function to be simulated.\n")
                          f         <- eval(parse(text = e))
                          n         <- 500
                          x         <- 0.1 + 0.8 * runif(n)
                          J         <- 20
                          vphi      <- sqrt(2) * cos(outer(x, pi * (1:J)))
                          loghyper  <- rep(1, 3)
                          Z         <- rep(1, times = n)
                          yStar     <- f(x) + Z
                          y         <- yStar + 0.1 * rnorm(n)
                          y_centred <- y - mean(y)

                          X_tr      <<- X_tst <<- as.matrix(x)
                          T_tr      <<- as.matrix(y)
                          T_tst     <<- as.matrix(yStar)
                          x_rest    <<- as.matrix(Z)

                          eval(q_mubeta.0)
                          eval(q_sigbeta.0)
                          eval(q_prior.params)

                          res_SSGP     <- eval(call_SSGPR_tst)
                          res_BSAR     <- eval(call_BSAR)
                          mu_BSAR      <- vphi[,1:length(res_BSAR$mutheta.q)]%*%res_BSAR$mutheta.q + res_BSAR$mubeta.q
                          centred_SSGP <- res_SSGP$mu - mean(res_SSGP$mu)
                          centred_BSAR <- mu_BSAR - mean(mu_BSAR)
                          eval(q_SSE_SSGP_tr)
                          eval(q_SSE_BSAR_tr)

                          o            <- order(x)
                          plot(x[o], y_centred[o], xlab = '', ylab = 'fitted', pch = 20, col = rgb(0, 0, 0, 0.3), main = 'Simulated data / test')
                          lines(x[o], centred_SSGP[o], lwd = 2, lty = 2, col = 'red')
                          lines(x[o], centred_BSAR[o], lwd = 2, lty = 3, col = 'darkgreen')
                          lines(x[o], yStar[o], lty = 6, col = 'purple')
                          legend('topright', lty = c(NA, 2, 3, 6), pch = c(20, NA, NA, NA), col = c(rgb(0, 0, 0, 0.3) , 'red', 'darkgreen', 'purple'), legend = c('observed', paste0('SSGP: ', SSE_SSGP_tr, ' (SSE)'), paste0('BSAR: ', SSE_BSAR_tr, ' (SSE)'), 'True curve'))
                          return(list(res_SSGP = res_SSGP, res_BSAR = res_BSAR, mu_BSAR = mu_BSAR, mu_SSGP = res_SSGP$mu, centred_SSGP = centred_SSGP, centred_BSAR = centred_BSAR))
               })

  stop("You've inserted a nonexistent option. Please pick from:\n c('pendulum', 'elevators', 'kin', 'pol', 'pumadyn', 'simul').")
}