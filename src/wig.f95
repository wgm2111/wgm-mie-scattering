  
! module 
!   use fact
! contains
  
subroutine wignerd_cos(x, m, n, lmax, d, ix, jx)
  ! calculate the general spherical function values
  
  use fact
  implicit none

  ! data dictionary in 
  integer, parameter :: r_=8
  integer :: ix, jx
  integer, intent(in) :: m, n, lmax
  real(kind=8), intent(in), dimension(ix, jx) :: x

  ! data dictionary out
  real(kind=8), dimension(ix, jx, 0:lmax-1), intent(out) :: d

  ! intermediate variables
  real(kind=8) :: xi, pi
  integer :: lmin, l



  ! constant pi
  pi = 3.14159265358979323

  ! define the first non-zero index (lmin)
  lmin = max(abs(m), abs(n))

  if (n>=m) then 
     xi = 1.0
  else
     xi = (-1)**(m-n)
  ENDIF

  d(:,:,:) = 0.0
  if (lmin==0) then
     ! base case
     d(:,:,0) = 1.0
     d(:,:,1) = x(:,:)
     ! loop up to lmax-1

     DO l = lmin+1, lmax-2
        ! print*, "hello from fortran, l =", l
        d(:,:,l+1) = (&
             &((2.0*l+1.0) * (l*(l+1.0)*x - m*n) * d(:,:,l) &
             &-(l+1.0)*sqrt(l**2.0-m**2)*sqrt(l**2.0-n**2) * d(:,:,l-1))&
             &/real((l * sqrt((l+1.0)**2-m**2) * sqrt((l+1.0)**2-n**2)),8))
     ENDDO
  endif
  if (lmin>=1) then  
     ! base case
     ! d(1,1,0) = xi * 2.0**(-lmin)

     d(:,:,lmin) = xi * 2.0D0**(-lmin) * (&
          & sqrt(f(real(2.0D0*lmin,r_)) / &
          & (f(real(abs(m-n),r_)) * f(real(abs(m+n), r_))))&
          & )*(1.0D0-x)**(abs(m-n)/2.0D0) * (1.0D0+x)**(abs(m+n)/2.0D0)
     ! LOOP UP TO LMAX
     DO l = lmin, lmax-2
        d(:,:,l+1) = (&
             & ((2.0D0*l+1.0D0)*(l*(l+1.0D0)*x - m*n) * d(:,:,l)-&
             &  (l+1.0D0)*sqrt(l**2.0D0-m**2.0D0)*sqrt(l**2.0D0-n**2.0D0)*d(:,:,l-1))&
             & /(l*sqrt((l+1.0D0)**2-m**2) * sqrt((l+1.0D0)**2-n**2)))
     ENDDO
  ENDIF

  ! if (mod(m-n,2)==0) then
  !    coef = (-1)**((m-n)/2.0)
  ! else
  !    gs(:,:,:) = -1.
  ! endif

end subroutine wignerd_cos





! PROGRAM MAIN
!   ! Calculate the values of wigner d functions 
  
!   USE FACT

!   IMPLICIT NONE 
  
!   ! data dictionary - declare constants 
!   INTEGER(KIND=8), PARAMETER :: ITHETA=2
!   INTEGER(KIND=8), PARAMETER :: JTHETA=3
!   INTEGER(KIND=8), PARAMETER :: LMAX=6


  
!   ! data dictionary - declare arrays 
!   INTEGER(KIND=8) :: M, N, I, J, LMIN
!   REAL(KIND=8) :: PI, XI
!   REAL(KIND=8), DIMENSION(ITHETA, JTHETA) :: THETA, X
!   REAL(KIND=8), DIMENSION(ITHETA, JTHETA, 0:LMAX-1) :: WIGNERD
  
!   PI = 3.14159265358979323
!   M = 0
!   N = 0
  
!   ! TURN THETA INTO IT'S COSINE
!   DO I = 1, ITHETA
!      DO J = 1, JTHETA
!         THETA(I,J) = PI * &
!              &((I-1)/REAL(ITHETA-1,8) + (J-1)/REAL(JTHETA-1,8))/2.0D1
!      ENDDO
!   ENDDO
  
!   X = COS(THETA)
      

!   ! DEFINE SOME CONSTANTS FOR THE ITERATION
!   LMIN = MAX(ABS(M), ABS(N))
!   IF (N>=M) THEN 
!      XI = 1.0
!   ELSE
!      XI = (-1)**(M-N)
!   ENDIF

!   call get_gensph(x, m, n, lmax, wignerd, itheta, jtheta)


!   ! WIGNERD(:,:,:) = 0.0

!   ! IF (LMIN == 0) THEN
!   !    ! DEFINE BASE CASE           
!   !    WIGNERD(:,:,0) = 1.0
!   !    WIGNERD(:,:,1) = X(:,:)
!   !    ! LOOP UP TO LMAX
!   !    WRITE(*,*) "HI"
!   !    DO L = LMIN+1, LMAX-1
!   !       WRITE(*,*) "IN LOOP ON L = ", L 
!   !       WIGNERD(:,:,L+1) = (&
!   !            &((2.0*L+1.0) * (L*(L+1.0)*X - M*N) * WIGNERD(:,:,L) &
!   !            &-(L+1.0)*SQRT(L**2.0-M**2)*SQRT(L**2.0-N**2) * WIGNERD(:,:,L-1))&
!   !            &/REAL((L * SQRT((L+1.0)**2-M**2) * SQRT((L+1.0)**2-N**2)),8))
!   !    ENDDO
     
!   ! ELSE
!   !    ! DEFINE BASE CASE
!   !    WIGNERD(:,:,LMIN) = XI * 2**(-LMIN) * (&
!   !         & SQRT(f(real(2.0*LMIN,8)) / &
!   !         & (f(real(ABS(M-N),8)) * f(real(ABS(M+N),8))))&
!   !         & ) * (1.0-X)**(ABS(M-N)/2.0) * (1.0+X)**(ABS(M+N)/2.0)
!   !    ! LOOP UP TO LMAX
!   !    DO L = LMIN, LMAX-1
!   !       WIGNERD(:,:,L+1) = (&
!   !            & ((2*L+1)*(L*(L+1)*X - M*N) * WIGNERD(:,:,L)-&
!   !            &  (L+1)*SQRT(L**2.-M**2)*SQRT(L**2.-N**2)*WIGNERD(:,:,L-1))&
!   !            & /(L*SQRT((L+1.)**2-M**2) * SQRT((L+1.)**2-N**2)))
!   !    ENDDO
!   ! ENDIF

!   WRITE(*,*) 'BUDMO', f(real(lmax,8))
!   WRITE(*,*) ITHETA, JTHETA, M, N, LMAX, PI
!   WRITE(*,*) '====='
!   ! WRITE(*,*) SHAPE(THETA)
!   DO I = 1, ITHETA
!      WRITE(*,100) X(I,:)
!   ENDDO
! 100 FORMAT(1X,5F10.4)

!   WRITE(*,*) '====='
!   WRITE(*,*) SHAPE(WIGNERD)
!   DO I = 1, ITHETA
!      DO J = 1, JTHETA
!         WRITE(*,200) WIGNERD(I, J, :)
!      ENDDO
!   ENDDO
! 200 FORMAT(1X, 8F10.4)

! END PROGRAM MAIN

     








! SUBROUTINE wignerd(THETA, M, N, LMAX, NTHETA, MTHETA, WIGNERD)
!   ! calculate and return the wigner-d functions
  
!   INTEGER :: M, N, LMAX, NTHETA, MTHETA
  


! END SUBROUTINE wignerd
