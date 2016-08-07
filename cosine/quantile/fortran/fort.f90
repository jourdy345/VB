SUBROUTINE fort(X,p,res)
  USE ToolsRfunf
  IMPLICIT NONE
  INTEGER,INTENT(IN)  :: p
  REAL*8, INTENT(IN)  :: X(p,p)
  REAL*8, INTENT(OUT) :: res(p,p)
  CALL inverse(X,p,res)

! CONTAINS
! subroutine inverse(R,p,Ri)
!   implicit none

!   ! Input arguments:
!   integer, intent(in) :: p
!   real*8, intent(in)  :: R(p,p)

!   ! Output argument
!   real*8, intent(out) :: Ri(p,p)

!   ! Internal arguments:
!   integer :: i,j,ok

!   Ri=R
!   call dpotrf('U',p,Ri,p,ok)  ! Cholesky decomposition using Lapack
!     call dpotri('U',p,Ri,p,ok)
!     do i=1,(p-1)
!       do j=(i+1),p
!           Ri(j,i) = Ri(i,j)
!       end do
!     end do

! end subroutine inverse
END SUBROUTINE fort

! SUBROUTINE fort(n,x)
!   IMPLICIT NONE
!   INTEGER, INTENT(IN)  :: n
!   real*8, INTENT(OUT) :: x(n)
!   INTEGER :: i
!   REAL*8 rvnorm
!   CALL rndstart()
!   DO i=1,n
!     x(i) = rvnorm(0.d0,1.d0)
!   END DO
!   CALL rndend()
! END SUBROUTINE fort

