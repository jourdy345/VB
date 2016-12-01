SUBROUTINE quantilef(y,vphi,W,A,B,muBeta,SigmaBeta,w0,nbasis,quant,nSample,burnIn,thinIn,&
                      beta,theta,gamma,tau,n,p)
  USE ToolsRfunf
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: nbasis,n,p,nSample,burnIn,thinIn
  REAL*8, INTENT(IN) :: y(n),vphi(n,nbasis),W(n,p),A,B,muBeta(p),SigmaBeta(p,p),w0,quant
  REAL*8,INTENT(OUT) :: beta(p,nSample),theta(nbasis,nSample),gamma(nSample),tau(nSample)

  REAL*8             :: rndnorm,exprnd,rndunif,gamrnd,betatmp(p),thetatmp(nbasis)
  REAL*8             :: gammatmp,tautmp,Aq,nup,taup2,Estar(nbasis),E(nbasis)
  REAL*8             :: Jvec(nbasis),u(n),const2,sigb(p,p),mub(p),sigt(nbasis,nbasis),mut(nbasis)
  REAL*8             :: tmpb(n,p),tmpt(n,nbasis),const1,tmp(n)
  REAL*8             :: SigmaBetainv(p,p),SigmaBetainvmuBeta(p),t2(n),tp(p),tt(nbasis),gcan,rho
  INTEGER            :: i,j,t,count
  CHARACTER*20       :: process
  CHARACTER*28       :: process2
  ! CHARACTER*6        :: countm
  CALL inverse(SigmaBeta,p,SigmaBetainv)
  SigmaBetainvmuBeta = MATMUL(SigmaBetainv,muBeta)
  const1 =  DBLE(nbasis)*(DBLE(nbasis)+1.d0)/4.d0 - w0
  CALL rndstart()
  !---initialize---!
  DO i=1,p
    betatmp(i)  = rndnorm()
  END DO
  DO i=1,nbasis
    thetatmp(i) = rndnorm()
    Jvec(i)     = DBLE(i)
  END DO
  gammatmp      = 0.01d0
  tautmp        = 1.d0/exprnd(1.d0)
  Aq            = A+0.5d0*DBLE(nbasis)
  nup           = (1.d0-2.d0*quant)/(quant*(1.d0-quant))
  taup2         = 2.d0/(quant*(1.d0-quant))
  count         = 0
  const2        = 0.25d0*taup2
  DO t=1,burnIn
    CALL rchkusr()
    count = count+1
    WRITE(process,fmt='(A,I6.6)') 'Burning in... ',count
    ! process = 'Burning in... '//countm//' round'
    CALL DBLEPR(process,-1,1.d0,0)
    ! PRINT *,count," iterations (Burn-in)"
    Estar = DEXP(Jvec*gammatmp)
    CALL DBLEPR('>',-1,1.d0,0)
    E     = Estar/tautmp
    !---update u---!
    CALL DBLEPR('>>',-1,1.d0,0)
    tmp = (y-matmul(W,betatmp)-matmul(vphi,thetatmp))**2.d0/taup2
    CALL DBLEPR('>>>',-1,1.d0,0)
    DO j=1,n
      u(j) = gigrnd(0.5d0,const2,tmp(j))
    END DO
    CALL DBLEPR('>>>>',-1,1.d0,0)
    !---update beta---!
    DO j=1,nbasis
      tmpb(:,j) = W(:,j)*DSQRT(1.d0/u)
    END DO
    CALL DBLEPR('>>>>>',-1,1.d0,0)
    CALL inverse(MATMUL(transpose(tmpb),tmpb)/taup2+SigmaBetainv,p,sigb)
    CALL DBLEPR('>>>>>>',-1,1.d0,0)
    DO j=1,p
      t2 = (y-MATMUL(vphi,thetatmp)-nup*u)/u
      tmpb(:,j) = W(:,j)*t
    END DO
    CALL DBLEPR('>>>>>>>',-1,1.d0,0)
    tp = SUM(tmpb,2)/taup2+SigmaBetainvmuBeta
    mub  = MATMUL(sigb,tp)
    CALL DBLEPR('>>>>>>>>',-1,1.d0,0)
    CALL mvnrnd(mub,sigb,p,betatmp)
    !---update tau---!
    CALL DBLEPR('>>>>>>>>>',-1,1.d0,0)
    tautmp = 1.d0/gamrnd(Aq,1.d0/(0.5d0*SUM(Estar*thetatmp**2.d0)+B))
    !---update theta---!
    DO j=1,nbasis
      tmpt(:,j) = vphi(:,j)*DSQRT(1.d0/u)
    END DO
    sigt = MATMUL(TRANSPOSE(tmpt),tmpt)/taup2
    DO j=1,nbasis
      sigt(j,j)=sigt(j,j)+E(j)
    END DO
    CALL inverse(sigt,nbasis,sigt)
    tt = SUM((y-MATMUL(W,betatmp)-nup*u)/u,DIM=1)/taup2
    mut = MATMUL(sigt,tt)
    CALL mvnrnd(mut,sigt,nbasis,thetatmp)
    gcan = exprnd(1.d0)
    rho = fRho(gcan,gammatmp,const1,tautmp,thetatmp,nbasis)
    IF (rndunif().LT.rho) THEN
      gammatmp = gcan
    END IF
  END DO
  count = 0
  DO i=1,nSample
    count = count+1
    CALL rchkusr()
    ! PRINT *,count,"samples collected..."
    WRITE(process2,fmt='(A,I6.6)') 'Collecting samples... ',count
    ! process = "Collecting samples... "//countm//' samples so far'
    CALL DBLEPR(process2,-1,1.d0,0)
    DO t=1,thinIn
      Estar = DEXP(Jvec*gammatmp)
      E     = Estar/tautmp
      !---update u---!
      tmp = (y-matmul(W,betatmp)-matmul(vphi,thetatmp))**2.d0/taup2
      DO j=1,n
        u(j) = gigrnd(0.5d0,const2,tmp(j))
      END DO
      !---update beta---!
      DO j=1,nbasis
        tmpb(:,j) = W(:,j)*DSQRT(1.d0/u)
      END DO
      CALL inverse(MATMUL(transpose(tmpb),tmpb)/taup2+SigmaBetainv,p,sigb)
      DO j=1,p
        t2 = (y-MATMUL(vphi,thetatmp)-nup*u)/u
        tmpb(:,j) = W(:,j)*t
      END DO
      tp = SUM(tmpb,2)/taup2+SigmaBetainvmuBeta
      mub  = MATMUL(sigb,tp)
      CALL mvnrnd(mub,sigb,p,betatmp)
      !---update tau---!
      tautmp = 1.d0/gamrnd(Aq,1.d0/(0.5d0*SUM(Estar*thetatmp**2.d0)+B))
      !---update theta---!
      DO j=1,nbasis
        tmpt(:,j) = vphi(:,j)*DSQRT(1.d0/u)
      END DO
      sigt = MATMUL(TRANSPOSE(tmpt),tmpt)/taup2
      DO j=1,nbasis
        sigt(j,j)=sigt(j,j)+E(j)
      END DO
      CALL inverse(sigt,nbasis,sigt)
      tt = SUM((y-MATMUL(W,betatmp)-nup*u)/u,DIM=1)/taup2
      mut = MATMUL(sigt,tt)
      CALL mvnrnd(mut,sigt,nbasis,thetatmp)
      gcan = exprnd(1.d0)
      CALL SPRINT(gcan)
      rho = fRho(gcan,gammatmp,const1,tautmp,thetatmp,nbasis)
      IF (rndunif().LE.rho) THEN
        gammatmp = gcan
      END IF
    END DO
    beta(:,i) = betatmp(1:p)
    theta(:,i) = thetatmp(1:nbasis)
    gamma(i) = gammatmp
    tau(i) = tautmp
  END DO
  CALL rndend()

CONTAINS
FUNCTION fRho(Yt,x,const1,tau,theta,J)
  IMPLICIT NONE
  REAL*8,INTENT(IN) :: Yt,x,const1,tau
  INTEGER,INTENT(IN) :: J
  REAL*8,INTENT(IN) :: theta(J)
  REAL*8 :: fRho
  INTEGER :: i
  REAL*8 :: Jvec(J)

  DO i=1,J
    Jvec(i)=DBLE(i)
  END DO
  fRho = DEXP(const1*(Yt-x)-SUM((DEXP(Yt*Jvec)-DEXP(x*Jvec))*theta**2.d0)/(2.d0*tau)-x+Yt)
  RETURN
END FUNCTION fRho
END SUBROUTINE quantilef