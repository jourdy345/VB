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