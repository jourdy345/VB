module ToolsRfunf
implicit none
real*8,parameter :: PI=3.141592653589793238462643383279502884197d0
real*8,parameter :: TRUNC=0.64d0,TRUNC_RECIP=1.d0/0.64d0
contains

!*****************************************************************************************
! Random number gernerators
!*****************************************************************************************

! C======================================================================
real*8 function rtnorm(mu,sd,a,b,ainf,binf)
! C=======================================================================
! C     generate truncated(a,b) Normal(mu,sd**2) using the Geweke's
! C     algorithm.
! C     mu is the mean of TN distribution
! C     sd is standard deviation of TN distribution
! C     a,b  = end points of interval (ainf = binf = .false.)
! C     ainf = true, if left endpoint is infinite; otherwise = false
! C     binf = true, if right endpoint is infinite; otherwise = false
! c     A.J.V., 2006
   implicit none
   real*8 mu,sd,a,b,a1,b1
   logical ainf,binf
   a1=(a-mu)/sd
   b1=(b-mu)/sd
   rtnorm=mu+sd*rtsnorm(a1,b1,ainf,binf)
   return
end function rtnorm

! c======================================================================
real*8 function rtsnorm(a,b,la,lb)
! c======================================================================
! c     generates a N(0,1) random variable
! c     subject to the constraint that it be in an interval
! c     (a,b), where the endpoints may be finite or infinite.
! c     a, b    endpoints of interval; a < b if la = lb = .false.
! c     la      .true. if left endpoint is - infinity; in this
! c             case A is ignored.
! c     lb      .true. if right endpoint is + infinity; in this
! c             case B is ignored.
! c     A.J.V., 2006
      implicit real*8 (a-h,o-z)
      logical la,lb,lflip
      real*8 rndunif,normrnd
      data eps,t1,t2,t3,t4/2.0d0,.375d0,2.18d0,.725d0,.45d0/

      if(la.and.lb)go to 160
      lflip=.false.
      if(la.or.lb)go to 100
      if(b.le.a)go to 170
! c ******* Finite interval
      c1=a
      c2=b
      if((c1*c2).gt.0.0d0)go to 30
! c ++++ (A,B) includes 0
      if((c1.gt.-t1).and.(c2.lt.t1))go to 20
! c -- F(A) or F(B) small: full normal with rejection
   10 x=normrnd(0.d0,1.d0)
      if(x.lt.c1)go to 10
      if(x.gt.c2)go to 10
      GO TO 150
! c -- F(A) and F(B) large: uniform importance sampling
   20 cdel=c2-c1
   25 x=c1+cdel*dble(rndunif())
      if(dble(rndunif()).gt.dexpone(x))go to 25
      go to 150
! c ++++ (A,B) excludes 0
! c -- Transform to both positive
   30 if(c1.gt.0.0d0)go to 40
      c=c1
      c1=-c2
      c2=-c
      lflip=.true.
   40 f1=dexpone(c1)
      f2=dexpone(c2)
      if(f2.lt.eps)go to 60
      if((f1/f2).gt.t2)go to 60
! c  -- F(A)/F(B) not large: uniform importance sampling
      cdel=c2-c1
   55 x=c1+cdel*rndunif()
      if(dble(rndunif()).gt.(dexpone(x)/f1))go to 55
      go to 140
   60 if(c1.gt.t3)go to 80
! c -- P(X>A) and F(A)/F(B) large: half-normal with rejection
   70 x=abs(normrnd(0.d0,1.d0))
      if(x.lt.c1)go to 70
      if(x.gt.c2)go to 70
      go to 140
! c -- P(X>A) small, F(A)/F(B) large: exponential importance
! c    sampling with rejection
   80 c=c2-c1
   90 z=-log(rndunif())/c1
      if(z.gt.c)go to 90
      if(dble(rndunif()).gt.dexpone(z))GO TO 90
      x=c1+z
      go to 140
! c ****** Half-line interval
  100 c1=a
! c -- Transform to bound from below if A = -infinity
      if(lb)go to 110
      c1=-b
      lflip=.true.
  110 if(c1.gt.t4)go to 130
! c -- A not large: full normal with rejection
  120 x=normrnd(0.d0,1.d0)
      if(x.lt.c1)go to 120
      go to 140
! c -- A small: exponential importance sampling
  130 z=-log(rndunif())/c1
      if(dble(rndunif()).gt.dexpone(z))go to 130
      x=c1+z
  140 if(lflip)x=-x
  150 rtsnorm=X
      return
! c ****** Whole interval
  160 rtsnorm=normrnd(0.d0,1.d0)
      return
! c  ***** Singleton
  170 rtsnorm=A
      return
end function rtsnorm

! c=======================================================================
real*8 function dexpone(x)
! c=======================================================================
! c     evaluate a exponential function
! c     A.J.V., 2006
      implicit none
      real*8 x,expin
      expin=-.5d0*x**2
      if (expin .le. -50.0d0) then
        dexpone=0.0d0
      else
        dexpone=dexp(expin)
      end if
      return
end function dexpone


!=========================================================================================
! Adaped from MATLAB code
!	Acklam, P. J. 'MATLAB array manipulation tips and tricks' (2003)
!
! 		Random sample, with replacement
!		The function returns a weighted sample using +ve weight probs,
!		taken with replacement, from the integers 1, 2, ..., n.
!-----------------------------------------------------------------------------------------
! Input arguments :
!	n	  : population
!	probs : weights
!
! Modified by Seongil Jo.
!=========================================================================================
function discrnd(n, probs)
    implicit none

    ! Input arguments
    integer, intent(in) :: n
    real*8, intent(in) :: probs(n)

    ! Output argument
    integer :: discrnd

    ! Internal argument
    real*8 :: rndunif,cum_probs(n)
    real*8 :: u
    integer :: i

    ! cumulative sum
    cum_probs=0.d0
    cum_probs(1)=probs(1)
    do i=2,n
        cum_probs(i)=cum_probs(i-1)+probs(i)
    end do

    u=rndunif()

    discrnd=1
    do i=1,n-1
        if (u .gt. cum_probs(i)) then
            discrnd = discrnd+1
        else
            discrnd = discrnd
            exit
        end if
    end do

    return
end function discrnd


!=========================================================================================
! This code generates one draw from the generalized inverse Gaussian distribution.
! The algorithm is based on that given by Dagpunar (1989).
!
!		   (psi/chi)^(lambda/2)*X^(lambda-1)*exp{-0.5(chi/X+psi*X)}
! f(X) = ------------------------------------------------------------
!			    			2K_lambda(sqrt(chi*psi))
!
!-----------------------------------------------------------------------------------------
! Input arguments
!	- psi		: shape parameter, >= 0
!	- chi		: shape parameter, >= 0
!	- lambda    : shape parameter, real
!
! Remark.
!	lambda > 0; psi > 0,  chi >= 0 or
!	lambda = 0; psi > 0,  chi > 0  or
!	lambda < 0; psi >= 0, chi > 0
!
! Created by Seongil Jo.
!=========================================================================================
function gigrnd(lambda,psi,chi)
	implicit none

	! Input arguments
	real*8,intent(in) :: lambda,psi,chi

	! Output argument
	real*8 :: gigrnd

	! Internal arguments
	real*8 :: gamrnd,invgaussrnd
	real*8 :: param(3),a,b

	! check parameters
	if (chi.lt.0.d0) call rexit('chi must be non-negative')
	if (psi.lt.0.d0) call rexit('psi must be non-negative')
	if (chi.eq.0.d0) then
		if (lambda.le.0.d0) call rexit('lambda must be positive when chi = 0')
		if (psi.eq.0.d0) call rexit('psi and chi cannot both be 0')
	end if
	if (psi.eq.0.d0) then
		if (lambda.ge.0.d0) call rexit('lambda must be negative when psi = 0')
		if (chi.eq.0.d0) call rexit('chi and psi cannot both be 0')
	end if

	! generate random variables
	if (chi.eq.0.d0) then
		! draw r.v. from gamma dist'n with a=lambda, b=psi/2 when chi=0
		a=lambda
		b=psi/2.d0
		gigrnd=gamrnd(a,1.d0/b)
	else if (psi.eq.0.d0) then
		! draw r.v. from inv-gamma dist'n with a=-lambda, b=chi/2 when psi=0
		a=-lambda
		b=chi/2.d0
		gigrnd=1.d0/gamrnd(a,1.d0/b)
	else if (lambda.eq.-0.5d0) then
		! draw r.v. from inv-Gaussian dist'n with a=sqrt(chi/psi), b=chi when lambda=-1/2
		b=chi
		a=dsqrt(b/psi)
		gigrnd=invgaussrnd(a,b)
	else if (lambda.eq.1.d0) then
		gigrnd=rgig1(psi,chi)
	else
		gigrnd=rgig(lambda,psi,chi)
	end if

	return
end function gigrnd


!=========================================================================================
! Main version to generate random variables from a generalized inverse
! Gaussian distribution
!=========================================================================================
function rgig(lambda,psi,chi)
	implicit none

	! Input arguments
	real*8,intent(in) :: lambda,psi,chi

	! Output argument
	real*8 :: rgig

	! Internal arguments
	real*8 :: tol,param(3),alpha,beta,m,upper,yM,yP,a,b,c,output,R1,R2,Y,Yc
	real*8 :: rndunif,powerxy

	tol=epsilon(lambda)  ! machine epsilon, 2.22d-16

	alpha=dsqrt(psi/chi)
	beta=dsqrt(psi*chi)
	m=(lambda-1.d0+dsqrt(powerxy(lambda-1.d0,2.d0)+powerxy(beta,2.d0)))/beta

	param(1)=lambda
	param(2)=beta
	param(3)=m

	upper=m
	do
		if (gf(upper,param).gt.0.d0) exit
		upper=2*upper
	end do

	yM=zeroin(0.d0,m,param,tol)
	yP=zeroin(m,upper,param,tol)

	a=(yP-m)*powerxy(yP/m,0.5d0*(lambda-1.d0))*dexp(-0.25d0*beta*(yP+1.d0/yP-m-1.d0/m))
	b=(yM-m)*powerxy(YM/m,0.5d0*(lambda-1.d0))*dexp(-0.25d0*beta*(yM+1.d0/yM-m-1.d0/m))
	c=-0.25d0*beta*(m+1.d0/m)+0.5d0*(lambda-1.d0)*dlog(m)

	output=0.d0
	do
		R1=rndunif()
		R2=rndunif()
		Y=m+a*R2/R1+b*(1.d0-R2)/R1
		Yc=-0.5d0*(lambda-1.d0)*dlog(Y)+0.25d0*beta*(Y+1.d0/Y)+c
		if (Y.gt.0.d0 .and. -dlog(R1).ge.Yc) then
			output=Y
			exit
		end if
	end do

	rgig=output/alpha

	return
end function rgig


!=========================================================================================
! Modified version of gigrnd to generate random variables from a generalized inverse
! Gaussian distribution with lambda=1
!=========================================================================================
function rgig1(psi,chi)
	implicit none

	! Input arguments
	real*8,intent(in) :: psi,chi

	! Output argument
	real*8 :: rgig1

	! Internal arguments
	real*8 :: tol,param(3),alpha,beta,m,upper,yM,yP,a,b,c,output,R1,R2,Y
	real*8 :: rndunif

	tol=epsilon(psi)  ! machine epsilon, 2.22d-16

	alpha=dsqrt(psi/chi)
	beta=dsqrt(psi*chi)
	m=dabs(beta)/beta

	param(1)=1.d0
	param(2)=beta
	param(3)=m

	upper=m
	do
		if (gf(upper,param).gt.0.d0) exit
		upper=2*upper
	end do

	yM=zeroin(0.d0,m,param,tol)
	yP=zeroin(m,upper,param,tol)

	a=(yP-m)*dexp(-0.25d0*beta*(yP+1.d0/yP-m-1.d0/m))
	b=(yM-m)*dexp(-0.25d0*beta*(yM+1.d0/yM-m-1.d0/m))
	c=-0.25d0*beta*(m+1.d0/m)

	output=0.d0
	do
		R1=rndunif()
		R2=rndunif()
		Y=m+a*R2/R1+b*(1.d0-R2)/R1
		if (Y.gt.0.d0 .and. -dlog(R1).ge.(0.25d0*beta*(Y+1.d0/Y)+c)) then
			output=Y
			exit
		end if
	end do

	rgig1=output/alpha

	return
end function rgig1


!=========================================================================================
! This function is used in generalized inverse Gaussian random generator.
!=========================================================================================
function gf(x,param)
	implicit none

	! Input arguments
	real*8,intent(in) :: x,param(3)

	! Output argument
	real*8 :: gf

	! Internal arguments
	real*8 :: powerxy
	real*8 :: lambda,beta,m

	lambda=param(1)
	beta=param(2)
	m=param(3)

	gf=0.5d0*beta*powerxy(x,3.d0)-powerxy(x,2.d0)*(0.5d0*beta*m+lambda+1)+ &
	   x*((lambda-1.d0)*m-0.5d0*beta)+0.5d0*beta*m

	return
end function gf


!=========================================================================================
! Random number generator from the Polya-Gamma distribution.
! Polson, Scott, and Windle (2013)
!-----------------------------------------------------------------------------------------
! Created by Seongil Jo.
!=========================================================================================
real*8 function devroyePG(z)
   implicit none

   !input arguments
   real*8,intent(in) :: z

   !internal arguments
   logical :: goiter
   integer :: n
   real*8 :: absz,fz,X,S,Y,exprnd,rndunif

   absz=dabs(z)*0.5d0
   fz=0.125d0*PI*PI+0.5d0*absz*absz
   X=0.d0
   S=1.d0
   Y=0.d0

   do
      if(rndunif().lt.mass_texpon(absz)) then
         X=TRUNC+exprnd(1.d0)/fz
      else
         X=rtigaussPG(z)
      end if
      S=a_coef(0,X)
      Y=rndunif()*S
      n=0
      goiter=.true.

      do while(goiter)
         call rchkusr() !break infinite loop.
         n=n+1
         if(mod(n,2).eq.1) then
            S=S-a_coef(n,X)
            if(Y.le.S) then
               devroyePG=0.25d0*X
               return
            end if
         else
            S=S+a_coef(n,X)
            if(Y.gt.S) goiter=.false.
         end if
      end do
   end do
end function devroyePG


real*8 function pgrnd(n,z)
   implicit none

   !input arguments
   integer,intent(in) :: n
   real*8, intent(in) :: z

   !internal arguments
   integer :: i
   real*8 :: sumz

   sumz=0.d0
   do i=1,n
      sumz=sumz+devroyePG(z)
   end do
   pgrnd=sumz

   return
end function pgrnd

!=========================================================================================
! Random number generator from truncated inv-Gauss(1/|z|,1.d0)I(0,0.64).
! Polson, Scott, and Windle (2013)
!-----------------------------------------------------------------------------------------
! Created by Seongil Jo.
!=========================================================================================
real*8 function rtigaussPG(z)
   implicit none

   !input argument
   real*8,intent(in) :: z

   !internal arguments
   real*8 :: X,Y,t,alpha,E1,E2,mu,half_mu,mu_Y,absz
   real*8 :: exprnd,rndunif,rndnorm

   absz=dabs(z)
   t=TRUNC
   X=t+1.d0
   if(TRUNC_RECIP.gt.absz) then
      alpha=0.d0
      do while (rndunif().gt.alpha)
         E1=exprnd(1.d0)
         E2=exprnd(1.d0)
         do while (E1*E1.gt.(2.d0)*E2/t)
            E1=exprnd(1.d0)
            E2=exprnd(1.d0)
         end do
         X=1.d0+E1*t
         X=t/(X*X)
         alpha=dexp(-0.5d0*absz*absz*X)
      end do
   else
      mu=1.d0/absz
      do while (X.gt.t)
         Y=rndnorm()
         Y=Y*Y
         half_mu=0.5d0*mu
         mu_Y=mu*Y
         X=mu+half_mu*mu_Y-half_mu*dsqrt(4.d0*mu_Y+mu_Y*mu_Y)
         if (rndunif().gt.(mu/(mu+X))) then
            X=mu*mu/X
         end if
      end do
   end if
   rtigaussPG=X

   return
end function rtigaussPG

!=========================================================================================
! The following functions are used in Polya-Gamma random generator.
!=========================================================================================
real*8 function mass_texpon(z)
  implicit none

  !input argument
  real*8,intent(in) :: z

  !internal arguments
  real*8 :: x,fz,a,b,x0,xa,xb,qdivp
  real*8 :: cdfnorm

  x=TRUNC
  fz=0.125d0*PI*PI+0.5d0*z*z
  b=dsqrt(1.d0/x)*(x*z-1.d0)
  a=-1.d0*dsqrt(1.d0/x)*(x*z+1.d0)

  x0=dlog(fz)+fz*x
  xb=x0-z+cdfnorm(b,0.d0,1.d0,1,1)
  xa=x0+z+cdfnorm(a,0.d0,1.d0,1,1)

  qdivp=4.d0/PI*(dexp(xb)+dexp(xa))

  mass_texpon=1.d0/(1.d0+qdivp)

  return
end function mass_texpon


real*8 function a_coef(n,x)
   implicit none

   !input arguments
   integer,intent(in) :: n
   real*8, intent(in) :: x

   !internal arguments
   real*8 :: K,y,expnt

   y=0.d0
   K=(dble(n)+0.5d0)*PI
   if (x.gt.TRUNC) then
      y=K*dexp(-0.5d0*K*K*x)
   else if (x.gt.0.d0) then
      expnt=-1.5d0*(dlog(0.5d0*PI)+dlog(x))+dlog(K)-2.d0*(dble(n)+0.5d0)*(dble(n)+0.5d0)/x
      y=dexp(expnt)
   end if
   a_coef=y

   return
end function a_coef

!=========================================================================================
! Random number generator from the Kolmogorov-Smirnov distribution.
! Holmes and Held (2006)
!-----------------------------------------------------------------------------------------
! Created by Seongil Jo.
!=========================================================================================
real*8 function devroyeKS(r)
   implicit none

   !input argument
   real*8,intent(in) :: r

   !internal arguments
   integer :: ok
   real*8 :: r2,lamp,uacc,rndunif

   r2=r**2.d0
   do
      lamp=gigrnd(0.5d0,1.d0,r2)  ! proposal from GIG(0.5,1,r^2)
      uacc=rndunif()
      if(lamp.gt.(4.d0/3.d0)) then
         ok=rightmost_interval(uacc,lamp)
      else
         ok=leftmost_interval(uacc,lamp)
      end if
      if (ok.eq.1) then
         devroyeKS=lamp
         exit
      end if
   end do

   return
end function devroyeKS

!=========================================================================================
! The following functions are used in Kolmogorov-Smirnov random generator.
!=========================================================================================
function rightmost_interval(u,lambda)
  implicit none

  !input arguments
  real*8,intent(in) :: u,lambda

  !output
  integer :: rightmost_interval

  !internal arguments
  integer :: ok
  real*8 :: j,R,Q

  ok=0
  R=1.d0
  Q=dexp(-0.5d0*lambda)
  j=0.d0
  do
  	call rchkusr() ! check interrupt

    ! Squeezing
    j=j+1.d0
    R=R-((j+1.d0)**2.d0)*(Q**((j+1.d0)**2.d0-1.d0))
    if(R.gt.u) then
       ok=1
       exit
    end if
    j=j+1.d0
    R=R+((j+1.d0)**2.d0)*(Q**((j+1.d0)**2.d0-1.d0))
    if(R.lt.u) then
      ok=0
      exit
    end if
  end do
  rightmost_interval=ok

  return
end function rightmost_interval


function leftmost_interval(u,lambda)
  implicit none

  !input arguments
  real*8,intent(in) :: u,lambda

  !output
  integer :: leftmost_interval

  !internal arguments
  integer :: ok
  real*8 :: j,R,Q,H,lu,K

  ok=0
  H=0.5d0*dlog(2.d0)+2.5d0*dlog(PI)-2.5d0*dlog(lambda)-&
    (PI**2.d0)/(2.d0*lambda)+0.5d0*lambda
  lu=dlog(u)
  R=1.d0
  Q=dexp(-(PI**2.d0)/(2.d0*lambda))
  K=lambda/(PI**2.d0)
  j=0.d0
  do
  	call rchkusr() ! check interrupt

    ! Squeezing
    j=j+1.d0
    R=R-K*(Q**((j**2.d0)-1.d0))
    if(H+dlog(R).gt.lu) then
      ok=1
      exit
    end if
    j=j+1.d0
    R=R+((j+1.d0)**2.d0)*(Q**((j+1.d0)**2.d0-1.d0))
    if(H+dlog(R).lt.lu) then
      ok=0
      exit
    end if
  end do
  leftmost_interval=ok

  return
end function leftmost_interval


!=========================================================================================
! This code generates one draw from the multivariate normal distribution.
! The algorithm is based on Cholesky decomposition in R
!-----------------------------------------------------------------------------------------
! Input arguments:
!	- mu		: mean vector
!	- cov		: covariance matrix
!	- p			: dimension
!
! Output arguments:
!	- rn		: random draw from MV-N(mu,cov) distribution
!
! Created by Seongil Jo.
!=========================================================================================
subroutine mvnrnd(mu,cov,p,rn)

	implicit none

	! Input arguments:
	integer, intent(in) :: p
	real*8, intent(in)  :: mu(p),cov(p,p)

	! Output arguments:
	real*8, intent(out) :: rn(p)

	! R functions:
	real*8 :: rndnorm

	! Internal arguments:
	real*8 :: L(p,p),z(p)
	integer :: i,j,ok

	! Cholesky factorization of cov using Lapack: cov=LL'
	L=cov
    call dpotrf('L',p,L,p,ok)

	! standard normal distribution
	do i=1,p
		z(i)=rndnorm()
	enddo

	do i=1,p
		rn(i)=mu(i)
		do j=1,i
			rn(i)=rn(i)+L(i,j)*z(j)
		enddo
	enddo

	return
end subroutine mvnrnd


!=========================================================================================
! Wishart random number generator using Batlett's method
!
!     This subroutine generates a random number
!	  from the Wishart distribution using Bartlett's method with density
!	  proportional to |COV|**(-0.5) * |Omega|**(-0.5 * (df - k - 1)) *
!	  						exp( -trace(0.5 * COV**(-1) * Omega) )
!-----------------------------------------------------------------------------------------
!     Input :
!		df : degree of freedom, df > p-1
!		cov : p by p matrix of covariance (+ve definite, sysmetric)
!		p : dimension
!
!	  Output:
!		w : p by p Wishart random matrix
!
!    1. Cov = L'L Cholesky decomposition
!    2. Z[i][j] ~ i.i.d. N(0,1), i > j & Z[i][j] = 0, i < j
!    3. Z[i][i] ~ sqrt(c), c ~ i.i.d. chisq(df-i+1)
!    4. A = ZZ' ~ wishart(nu,I)
!    5. rn = LAL' ~ wishart(df,Cov)
!
!		Ex :
!			cov = (/1,0.3,0.3,1/)
!			To draw a 2x2 random matrix from W_2(cov, 3)
!			call wishrng(w, cov, 3, 2)
!
! Created by Seongil Jo.
!=========================================================================================
subroutine wishrnd(cov,df,p,w)
  implicit none

  !input arguments
  integer,intent(in) :: p
  real*8, intent(in) :: df
  real*8, intent(in) :: cov(p,p)

  !output argument
  real*8,intent(out) :: w(p,p)

  ! local variable
  real*8 :: Z(p,p),A(p,p),L(p,p),LA(p,p)
  real*8 :: nu,chi,chisqrnd,rndnorm
  integer :: i,j,ok

  Z=0.d0
  A=0.d0
  L=0.d0

  ! Cholesky factorization of cov using Lapack: cov=LL'
  L=cov
  call dpotrf('L',p,L,p,ok)

  ! step 2. & 3.
  do i=1,p
    nu=df-dble(i) + 1.d0
    chi=chisqrnd(nu)
    Z(i,i)=dsqrt(chi)
    do j=1,i-1
      Z(i,j)=rndnorm()
    end do
  end do

  ! compute ZZ' = A ~ Wishart(nu, I)
  A=matmul(Z,transpose(Z))

  ! compute LAL' = rn ~ Wishart(nu, Cov)
  LA=matmul(L,A)
  w=matmul(LA,transpose(L))

  return
end subroutine wishrnd


!*****************************************************************************************
! Probability density function (pdf)
!*****************************************************************************************

!=========================================================================================
! This code provides the density function from the multivariate normal distribution.
! The algorithm is based on Cholesky decomposition in R
!-----------------------------------------------------------------------------------------
! Input arguments:
!	- mu		: mean vector
!	- cov		: covariance matrix
!	- d			: dimension
!	- log_p     : logical variable
!
! Output arguments:
!	- rn		: random draw from MV-N(mu,cov) distribution
!
! Created by Seongil Jo.
!=========================================================================================
function mvnpdf(x,mu,cov,d,log_p)
	implicit none

	! Input arguments
	logical,intent(in) :: log_p
	integer,intent(in) :: d
	real*8,intent(in) :: x(d),mu(d),cov(d,d)

	! Output argument
	real*8 :: mvnpdf

	! Internal arguments
	real*8 :: detcov_half,covi(d,d),logpdf,logconst
	integer :: i,j,ok

	! Cholesky decomposition cov=U'U
	covi=cov
	call dpotrf('U',d,covi,d,ok)

	! square root of determinant of covariance
  	detcov_half=1.d0
  	do i=1,d
  		detcov_half=detcov_half*covi(i,i)
  	end do

	! inverse of covariance
	call dpotri('U',d,covi,d,ok)
	do i=1,(d-1)
   		do j=(i+1),d
       		covi(j,i) = covi(i,j)
   		end do
  	end do

  	logconst=-dble(d)*dlog(2.d0*PI)/2.d0-dlog(detcov_half)
  	logpdf=logconst-dot_product(x-mu,matmul(covi,x-mu))/2.d0

	mvnpdf=dexp(logpdf)
	if (log_p) mvnpdf=logpdf

	return
end function mvnpdf


!=========================================================================================
! This function gives the density of the mixture of normals with
!   mean vector mu and standard deviation vector sigma, evaluated at the values in x.
!
!	Input :
!		x 	    : scalar
!		mu 	    : p x 1, vector of mean values
!		sigma   : p x 1, vector of standard deviations
!		w       : p x 1, vector of weigths
!		m 	    : the number of components
!		log_p   : logical variable
!
! Created by Seongil Jo.
!========================================================================================
function dnrmmix(x,mu,sigma,w,m,log_p)
	implicit none

	! Input arguments
	integer,intent(in) :: m
	logical,intent(in) :: log_p
	real*8, intent(in) :: x,mu(m),sigma(m),w(m)

	! Output arguments
	real*8 :: dnrmmix

	! Internal arguments
	real*8 :: dnrm,wn(m)
	integer :: i

	wn=w/sum(w) ! normalization

	dnrmmix = 0.d0
	do i=1,m
		dnrmmix=dnrmmix+wn(i)*dnrm(x,mu(i),sigma(i),.false.)
	enddo

    if (log_p) dnrmmix=dlog(dnrmmix)

	return
end function dnrmmix


!*****************************************************************************************
! Linear algebra
!*****************************************************************************************

!=========================================================================================
! This code calculates an inverse of a symmetric positive matrix.
! The algorithm is based on Cholesky decomposition in R
!-----------------------------------------------------------------------------------------
! Input arguments:
!	- R		: a real symmetric positive matrix
!	- p		: dimension
!
! Output arguments:
!	- Ri	: inverse of R
!=========================================================================================
subroutine inverse(R,p,Ri)
	implicit none

	! Input arguments:
	integer, intent(in) :: p
	real*8, intent(in)  :: R(p,p)

	! Output argument
	real*8, intent(out) :: Ri(p,p)

	! Internal arguments:
	integer :: i,j,ok

	Ri=R
	call dpotrf('U',p,Ri,p,ok)	! Cholesky decomposition using Lapack
   	call dpotri('U',p,Ri,p,ok)
   	do i=1,(p-1)
   		do j=(i+1),p
       		Ri(j,i) = Ri(i,j)
   		end do
  	end do

end subroutine inverse


!=========================================================================================
! This code calculates an inverse of a symmetric positive matrix.
! The algorithm is based on LU decomposition in R's Lapack
!-----------------------------------------------------------------------------------------
! Input arguments:
!	- R		: a real symmetric positive matrix
!	- p		: dimension
!
! Output arguments:
!	- Ri	: inverse of R
!=========================================================================================
subroutine inverseLU(R,p,Ri)
implicit none

! Input arguments:
integer, intent(in) :: p
real*8, intent(in)  :: R(p,p)

! Output argument
real*8, intent(out) :: Ri(p,p)

! Internal arguments:
integer :: m,lda,lwork,info,ipiv(p)	! ipiv=pivot indices & lda=leading dimension
real*8 :: work(p)

m=p
lda=p
lwork=p

Ri=R	! matrix inversion
call dgetrf(m,m,Ri,lda,ipiv,info)	! LU decomposition using Lapack
call dgetri(m,Ri,lda,ipiv,work,lwork,info )

return
end subroutine inverseLU



!=========================================================================================
! Determinant
! Method: Based on LU decomposition
!-----------------------------------------------------------------------------------------
! Input arguments:
! A(N,N) - original matrix
! N      - dimension
!
! Created by Seongil Jo.
!=========================================================================================
function determinant(R,p)
	implicit none

	! Input arguments
	integer, intent(in) :: p
	real*8, intent(in)  :: R(p,p)

	! Output argument
	real*8 :: determinant

	! Internal arguments
	integer :: i,ipiv(p),info
	real*8 :: Rlu(p,p)

	Rlu=R
	call dgetrf(p,p,Rlu,p,ipiv,info)	! LU decomposition using Lapack in R

	determinant = 0.d0
	if (info.ne.0) then
		return
	endif
	determinant = 1.d0
	do i=1,p
		if (ipiv(i).ne.i) then
			determinant=-determinant*Rlu(i,i)
		else
			determinant=determinant*Rlu(i,i)
		endif
	end do

	return
end function determinant


!*****************************************************************************************
! Other useful functions
!*****************************************************************************************

!=========================================================================================
! Generic function whose default method centers and/or scales
! the columns of a numeric matrix.
!      input:
!            x 		: nxp matrix
!			 center : logical, If center is TRUE then centering is done by subtracting
!					  the column means of x from their corresponding columns
!			 scale  : logical, If scale is TRUE then scaling is done by dividing
!					  the (centered) columns of x by their standard deviations
!            n : number of rows
!            p : number of columns
!      output:
!            z : nxp matrix centered and/or scaled
!=========================================================================================
subroutine sweep(x,n,p,center,scale,z)
	implicit none

	! input arguments
	logical,intent(in) :: center,scale
	integer,intent(in) :: n,p
	real*8, intent(in) :: x(n,p)

	! output argument
	real*8,intent(out) :: z(n,p)

	! internal arguments
	real*8 :: xmean(p),xsd(p)
	integer :: i

	z=x
	xsd=1.d0
	xmean=0.d0

	if (center) xmean=sum(x,1)/dble(n)
	if (scale) xsd=(/ (dsqrt(sum((x(:,i)-xmean(i))**2.d0,1)/(dble(n)-1.d0)), i=1,n) /)

	do i=1,p
		z(:,i)=(x(:,i)-xmean(i))/xsd(i)
	end do
end subroutine sweep


!=========================================================================================
! Computing 100(1-alpha)% HPD and credible intervals for x
! Use Chen-Shao HPD Estimation Algorithm
! (see page 219 of Monte Carlo Methods in Bayesian Computation, Springer-Verlag, 2000)
! By Ming-hui Chen, july 23, 2001 at wpi
!      input:
!            alpha: confidence level,  0 < alpha < 1
!            n = mcmc sample size
!            y(n): a univariate vector of mcmc sample
!      output:
!            (alow(1),aupp(1)): 100(1-alpha)% HPD interval
!            (alow(2),aupp(2)): 100(1-alpha)% Bayesian credible interval
!=========================================================================================
subroutine hpd(n,alpha,y,alow,aupp)
	implicit none

	integer, intent(in) :: n
	real*8, intent(in) :: y(n),alpha

    real*8, intent(out) :: aupp(2),alow(2)

	integer :: i,j,nq1,nq2,nq
	real*8 :: whpd,aupp1,alow1,q1,q2,temp,pdiff1,pdiff2,wb,x(n)

	x=y
    whpd=0.d0
    aupp1=0.d0
    alow1=0.d0

    q1=(alpha/2.0d0)*float(n)
    q2=(1.0d0-alpha/2.0d0)*float(n)
    nq1=nint(q1)
    nq2=nint(q2)
    nq=nq2-nq1
    do i=1,n-1
    	do j=i+1,n
        	if (x(i) .gt. x(j)) then
            	temp=x(i)
                x(i)=x(j)
                x(j)=temp
            end if
 	    end do
 	end do

	do j=1,n-nq
    	pdiff1=x(j)
        pdiff2=x(j+nq)
        wb=pdiff2-pdiff1
        if (j .eq. 1) then
        	whpd=wb
            aupp1=pdiff2
            alow1=pdiff1
        else
            if (whpd .gt. wb) then
            	whpd=wb
                aupp1=pdiff2
                alow1=pdiff1
             end if
    	end if
    end do
	alow(1)=alow1
    aupp(1)=aupp1
    alow(2)=x(nq1)
    aupp(2)=x(nq2)

    return
end subroutine hpd


!=========================================================================================
! Adapted from Wikipedia
!	http://rosettacode.org/wiki/Remove_duplicate_elements
!
!	Find the unique integer values by remove the duplicate values
!-----------------------------------------------------------------------------------------
!	Input :
!		x : n x 1, duplicated values
!       n : dimension
!
!	Output :
!		uniq_x : n x 1, k unique elements and (n - k) zero values
!		k      : scalar, the number of unique elements
!
! Modified by Seongil Jo.
!=========================================================================================
subroutine find_uniquei(uniq_x, k, x, n)
    implicit none

    ! Input arguments
    integer, intent(in) :: n
    integer, intent(in) :: x(n)

    ! Output arguments
    integer, intent(out) :: uniq_x(n),k

    ! Internal arguments
    integer :: i, j, si, m

    uniq_x=0
    m=1
    do
        if (x(m) /= 0) exit
        m=m+1
    end do
    uniq_x(1)=x(m)
    si=m+1
    k=1
    outer: do i=si,n
        do j=1,k
            if (uniq_x(j) .eq. x(i)) then
                ! Found a match so start looking again
                cycle outer
            endif
        end do
        ! No match found so add it to the output
        if (x(i) /= 0) then
            k=k+1
            uniq_x(k)=x(i)
        end if
    end do outer

end subroutine find_uniquei

!=========================================================================================
! Give the TRUE indices of a logical object
!-----------------------------------------------------------------------------------------
!	Input.
!		logic : n x 1, logical variable
!       n     : dimension
!
!	Output.
!		k     : the number of .true.
!		ind   : n x 1, .true. indices (k x 1), 0 ((n-k) x 1)
!
! Created by Seongil Jo.
!=========================================================================================
subroutine which(logic,n,ind,k)
    implicit none

    ! Input arguments
    integer,intent(in) :: n
    logical,intent(in) :: logic(n)

    ! Output arguments
    integer,intent(out) :: ind(n),k

    ! Internal arguments
    integer :: i,l,x(n),y(n)

    ind=0
    k=count(logic)
    x=transfer(logic, (/ (1, i = 1, n) /))
    y=x*(/ (i, i = 1, n) /)
    l=1
    do i=1,n
        if (y(i) .ne. 0) then
            ind(l)=y(i)
            l=l+1
        endif
    enddo

end subroutine which


!=========================================================================================
! $UWHPSC/codes/fortran/zeroin/zeroin.f
! This function came from Netlib: http://www.netlib.org/go/zeroin.f
!     With some modification to remove dependence on d1mach
!     by computing eps, the machine precision, in this routine.
!
! a zero of the function  f(x)  is computed in the interval ax,bx .
! (Standard routine from netlib)
!
! This function is used in generalized inverse Gaussian random number generator.
!-----------------------------------------------------------------------------------------
! Input..
!
!  ax   : left endpoint of initial interval
!  bx   : right endpoint of initial interval
!  gf   : function subprogram which evaluates gf(x) for any x in
!         the interval  ax,bx
!  tol  : desired length of the interval of uncertainty of the
!         final result ( .ge. 0.d0)
!
!
! Output..
!
!  zeroin : abcissa approximating a zero of  f  in the interval ax,bx
!
!
! It is assumed  that   f(ax)   and   f(bx)   have  opposite  signs
!  without  a  check.  zeroin  returns a zero  x  in the given interval
!  ax,bx  to within a tolerance  4*macheps*dabs(x) + tol, where macheps
!  is the relative machine precision.
!
! This function subprogram is a slightly  modified  translation  of
!  the algol 60 procedure  zero  given in  richard brent, algorithms for
!  minimization without derivatives, prentice - hall, inc. (1973).
!
! Modified Seongil Jo.
!=========================================================================================
function zeroin(ax,bx,param,tol)
	implicit none

	! Input arguments
	real*8,intent(in) :: ax,bx,tol
	real*8,intent(inout) :: param(*)

	! Output argument
	real*8 :: zeroin

	! Internal arguments
	real*8 :: eps,tol1,a,b,c,d,e,fa,fb,fc,xm,p,q,r,s
!
!  compute eps, the relative machine precision
!
	eps = 1.d0
 10 eps = eps/2.d0
    tol1 = 1.d0 + eps
    if (tol1 .gt. 1.d0) go to 10
!
! initialization
!
    a = ax
    b = bx
    fa = gf(a,param)
    fb = gf(b,param)
!
! begin step
!
 20 c = a
    fc = fa
    d = b - a
    e = d
 30 if (dabs(fc) .ge. dabs(fb)) go to 40
    a = b
    b = c
    c = a
    fa = fb
    fb = fc
    fc = fa
!
! convergence test
!
 40 tol1 = 2.d0*eps*dabs(b) + 0.5*tol
    xm = .5*(c - b)
    if (dabs(xm) .le. tol1) go to 90
    if (fb .eq. 0.d0) go to 90
!
! is bisection necessary
!
    if (dabs(e) .lt. tol1) go to 70
    if (dabs(fa) .le. dabs(fb)) go to 70
!
! is quadratic interpolation possible
!
    if (a .ne. c) go to 50
!
! linear interpolation
!
    s = fb/fa
    p = 2.d0*xm*s
    q = 1.d0 - s
    go to 60
!
! inverse quadratic interpolation
!
 50 q = fa/fc
    r = fb/fc
    s = fb/fa
    p = s*(2.d0*xm*q*(q - r) - (b - a)*(r - 1.d0))
    q = (q - 1.d0)*(r - 1.d0)*(s - 1.d0)
!
! adjust signs
!
 60 if (p .gt. 0.d0) q = -q
    p = dabs(p)
!
! is interpolation acceptable
!
    if ((2.d0*p) .ge. (3.d0*xm*q - dabs(tol1*q))) go to 70
    if (p .ge. dabs(0.5*e*q)) go to 70
    e = d
    d = p/q
    go to 80
!
! bisection
!
 70 d = xm
    e = d
!
! complete step
!
 80 a = b
    fa = fb
    if (dabs(d) .gt. tol1) b = b + d
    if (dabs(d) .le. tol1) b = b + dsign(tol1, xm)
    fb = gf(b,param)
    if ((fb*(fc/dabs(fc))) .gt. 0.d0) go to 20
    go to 30
!
! done
!
 90 zeroin = b
    return
end function zeroin


!=========================================================================================
! Construct a real diagonal matrix
!-----------------------------------------------------------------------------------------
!	Input.
!       x : value, scalar
!       n : dimension
!
!	Output.
!		A : diagonal matrix
!
! Created by Seongil Jo.
!=========================================================================================
subroutine diag(x,n,A)
    implicit none

    !input arguments
    real*8, intent(in) :: x
    integer,intent(in) :: n

    !output arguments
    real*8,intent(out) :: A(n,n)

    !internal argument
    integer :: i

    A=0.d0
    do i=1,n
        A(i,i)=x
    end do

end subroutine diag

!=========================================================================================
! Construct an integer diagonal matrix
!-----------------------------------------------------------------------------------------
!	Input.
!       x : value, scalar
!       n : dimension
!
!	Output.
!		A : diagonal matrix
!
! Created by Seongil Jo.
!=========================================================================================
subroutine Idiag(x,n,A)
    implicit none

    !input arguments
    integer,intent(in) :: x
    integer,intent(in) :: n

    !output arguments
    integer,intent(out) :: A(n,n)

    !internal argument
    integer :: i

    A=0
    do i=1,n
        A(i,i)=x
    end do

end subroutine Idiag


!=========================================================================================
! Construct a diagonal matrix
!-----------------------------------------------------------------------------------------
!	Input.
!       x : values, vector
!       n : dimension
!
!	Output.
!		A : diagonal matrix
!
! Created by Seongil Jo.
!=========================================================================================
subroutine diagvec(x,n,A)
implicit none

!input arguments
real*8, intent(in) :: x(n)
integer,intent(in) :: n

!output arguments
real*8,intent(out) :: A(n,n)

!internal argument
integer :: i

A=0.d0
do i=1,n
A(i,i)=x(i)
end do

end subroutine diagvec

!=========================================================================================
! Extract diagonal elements of a matrix
!-----------------------------------------------------------------------------------------
!	Input.
!       A : matrix
!       n : dimension
!
!	Output.
!		x : diagonal elements
!
! Created by Seongil Jo.
!=========================================================================================
subroutine diagExtract(A,n,x)
	implicit none

	!input arguments
	integer,intent(in) :: n
	real*8, intent(in) :: A(n,n)

	!output arguments
	real*8,intent(out) :: x(n)

	!internal argument
	integer :: i

	x=0.d0
	do i=1,n
		x(i)=A(i,i)
	end do

end subroutine diagExtract


!=========================================================================================
! Performs the vectorization operation
!-----------------------------------------------------------------------------------------
! Takes a matrix and turns it into a vector
! by stacking its columns.
!
! This subroutine came from matlab code of Daniel B. Rowe
!=========================================================================================
subroutine stackvec(mat,n,p,vec)
	implicit none

	!input arguments
	integer,intent(in) :: n,p
	real*8, intent(in) :: mat(n,p)

	!output argument
	real*8,intent(out) :: vec(n*p)

	!internal arguemts
	integer :: ir,ic

	vec=0.d0
	do ic=1,p
		do ir=1,n
			vec(ir+n*(ic-1))=mat(ir,ic)
		end do
	end do

	return
end subroutine stackvec


!=========================================================================================
! Reshapes the lower triangular portion of a symmetric matrix into a column vector
!-----------------------------------------------------------------------------------------
! returns a stack of the lower triangular matrix of a real square matrix as a matrix
!  with 1 column and n * ( n + 1 ) / 2 rows
!=========================================================================================
subroutine vech(mat,nr,nc,vec)
    implicit none

    !input arguments
    integer,intent(in) :: nr,nc
    real*8, intent(in) :: mat(nr,nc)

    !output argument
    real*8,intent(out) :: vec(nr * (nc + 1)/2)

    !internal arguemts
    integer :: ir,ic,k

    vec=0.d0  ! initialize the vector

    k=1
    do ir=1,nr
        do ic=1,ir
            if(ic.le.ir) then
                vec(k)=mat(ir,ic)
                k=k+1
            end if
        end do
    end do

    return
end subroutine vech


!=========================================================================================
! Performs the vectorization operation of the lower triangular matrix of a square matrix
!-----------------------------------------------------------------------------------------
! returns a stack of the lower triangular matrix of a integer square matrix as a matrix
!  with 1 column and n * ( n + 1 ) / 2 rows
!=========================================================================================
subroutine Ivech(mat,nr,nc,vec)
    implicit none

    !input arguments
    integer,intent(in) :: nr,nc
    integer,intent(in) :: mat(nr,nc)

    !output argument
    integer,intent(out) :: vec(nr * (nc + 1)/2)

    !internal arguemts
    integer :: ir,ic,k

    vec=0  ! initialize the vector

    k=1
    do ir=1,nr
        do ic=1,ir
            if(ic.le.ir) then
                vec(k)=mat(ir,ic)
                k=k+1
            end if
        end do
    end do

    return
end subroutine Ivech


!=========================================================================================
! Kronecker product of A and B
!-----------------------------------------------------------------------------------------
! Function to compute the tensor product K of two real matrices A and B
! nra(nrb) and nca(ncb) are, respectively,
! the number of rows and columns of the Matrix A(B)
! This subroutine came from Fortran code of Jonas Maziero.
!=========================================================================================
subroutine kron(A,nra,nca,B,nrb,ncb,K)
	implicit none

	!input arguments
	integer,intent(in) :: nra,nca,nrb,ncb
	real*8, intent(in) :: A(nra,nca),B(nrb,ncb)

	!output argument
	real*8,intent(out) :: K(nra*nrb,nca*ncb)

	!internal arguments
	integer :: i,j

	K=0.d0  ! initialize the matrix K

	forall (i = 1:nra , j = 1:nca)
  		K((nrb*(i-1)+1):(nrb*i),(ncb*(j-1)+1):(ncb*j))=A(i,j)*B
	end forall

	return
end subroutine kron


!=========================================================================================
! Kronecker product of A and B
!-----------------------------------------------------------------------------------------
! Function to compute the tensor product K of two integer matrices A and B
! nra(nrb) and nca(ncb) are, respectively,
! the number of rows and columns of the Matrix A(B)
! This subroutine came from Fortran code of Jonas Maziero.
!=========================================================================================
subroutine Ikron(A,nra,nca,B,nrb,ncb,K)
    implicit none

    !input arguments
    integer,intent(in) :: nra,nca,nrb,ncb
    integer,intent(in) :: A(nra,nca),B(nrb,ncb)

    !output argument
    integer,intent(out) :: K(nra*nrb,nca*ncb)

    !internal arguments
    integer :: i,j

    K=0  ! initialize the matrix K

    forall (i = 1:nra , j = 1:nca)
        K((nrb*(i-1)+1):(nrb*i),(ncb*(j-1)+1):(ncb*j))=A(i,j)*B
    end forall

    return
end subroutine Ikron


!=========================================================================================
! Construct the correlation matrix from a covariance matrix
!-----------------------------------------------------------------------------------------
! subroutine scales a covariance matrix into
! the corresponding correlation matrix efficiently.
!=========================================================================================
subroutine cov2cor(cov,n,cor)
	implicit none

	!input arguments
	integer,intent(in) :: n
	real*8, intent(in) :: cov(n,n)

	!output argument
	real*8,intent(out) :: cor(n,n)

	!internal arguments
	integer :: i,j
	real*8 :: v(n),Is(n)

	call diagExtract(cov,n,v)  ! diagonal elements

	cor=0.d0
	do j=1,n
		do i=j,n
			if(j.eq.i) then
				cor(i,j)=1.d0
			else
				cor(i,j)=cov(i,j)/(dsqrt(v(i))*dsqrt(v(j)))
				cor(j,i)=cor(i,j)
			end if
		end do
	end do

	return
end subroutine cov2cor

!=========================================================================================
! Construct a covariance matrix
!=========================================================================================
subroutine covariance(A,n,p,cov)
	implicit none

	!input arguments
	integer,intent(in) :: n,p
	real*8, intent(in) :: A(n,p)

	!output argument
	real*8,intent(out) :: cov(p,p)

	!internal arguments
	real*8 :: z(n,p)

	call sweep(A,n,p,.true.,.false.,z)
	cov=matmul(transpose(z),z)/dble(n-1)

	return
end subroutine covariance

!=========================================================================================
! Integration using Simpson's rule
!=========================================================================================
subroutine intsim(f,delta,n,fint)
	implicit none

	! Input arguments
	integer,intent(in) :: n
	real*8, intent(in) :: f(n),delta

	! Output argument
	real*8,intent(out) :: fint

	! Internal arguments
	integer :: i
	real*8 :: t(n)

	if(n.eq.2*floor(dble(n)/2.d0)) then
		call rexit('ERROR: Even number of rows for Simpson integration')
	else if(n.eq.3) then
		! Simple Simpson integration
		t(1)=1.d0
		t(2)=4.d0
		t(3)=1.d0
		fint=sum(t*f)*delta/3.d0
	else
		! Composite Simpson integration
		t(1)=1.d0
		do i=2,n-3,2
			t(i)=4.d0
			t(i+1)=2.d0
		end do
		t(n-1)=4.d0
		t(n)=1.d0
		fint=sum(t*f)*delta/3.d0
	end if

	return
end subroutine intsim


end module ToolsRfunf
