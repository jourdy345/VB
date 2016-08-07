MODULE aux
IMPLICIT NONE
REAL*8,PARAMETER :: PI=3.141592653589793238462643383279502884197d0
CONTAINS
SUBROUTINE sweep(x,n,p,center,scale,z)
  IMPLICIT NONE
  LOGICAL,INTENT(IN) :: center,scale
  INTEGER,INTENT(IN) :: n,p
  REAL*8, INTENT(IN) :: x(n,p)
  REAL*8,INTENT(OUT) :: z(n,p)
  REAL*8 :: xmean(p),xsd(p)
  INTEGER :: i

  z=x
  xsd=1.d0
  xmean=0.d0
  IF (center) xmean = SUM(x,1)/DBLE(n)
  IF (scale)  xsd   = (/ (DSQRT(SUM((x(:,i)-xmean(i))**2.d0,1)/(DBLE(n)-1.d0)),i=1,n) /)
  DO i=1,p
    z(:,i) = (x(:,i)-xmean(i))/xsd(i)
  END DO
END SUBROUTINE sweep

! REAL*8 FUNCTION determinant(R,p)
!   IMPLICIT NONE
!   INTEGER,INTENT(IN) :: p
!   REAL*8, INTENT(IN) :: R(p,p)
!   REAL*8 :: determinant
!   INTEGER :: i,ipiv(p),info
!   REAL*8 :: Rlu(p,p)

!   Rlu=R
!   CALL DGETRF(p,p,Rlu,p,ipiv,info)
!   determinant = 0.d0
!   IF (info .NE. 0) THEN
!     RETURN
!   ENDIF
!   determinant = 1.d0
!   DO i = 1,p
!     IF (ipiv(i) .NE. i) THEN
!       determinant = -determinant*Rlu(i,i)
!     ELSE
!       determinant = determinant*Rlu(i,i)
!     ENDIF
!   END DO
!   RETURN
! END FUNCTION determinant

SUBROUTINE inverse(R,p,Ri)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: p
  REAL*8, INTENT(IN)  :: R(p,p)
  REAL*8, INTENT(OUT) :: Ri(p,p)
  INTEGER :: i,j,ok
  Ri=R
  CALL DPOTRF('U',p,Ri,p,ok)
  CALL DPOTRI('U',p,Ri,p,ok)
  DO i = 1,(p-1)
    DO j = (i+1),p
      Ri(j,i)=Ri(i,j)
    END DO
  END DO
END SUBROUTINE inverse
END MODULE aux