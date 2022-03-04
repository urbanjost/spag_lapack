!*==sasum.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b SASUM
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       REAL FUNCTION SASUM(N,SX,INCX)
!
!       .. Scalar Arguments ..
!       INTEGER INCX,N
!       ..
!       .. Array Arguments ..
!       REAL SX(*)
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    SASUM takes the sum of the absolute values.
!>    uses unrolled loops for increment equal to one.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>         number of elements in input vector(s)
!> \endverbatim
!>
!> \param[in] SX
!> \verbatim
!>          SX is REAL array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
!> \endverbatim
!>
!> \param[in] INCX
!> \verbatim
!>          INCX is INTEGER
!>         storage spacing between elements of SX
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \date November 2017
!
!> \ingroup single_blas_level1
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>     jack dongarra, linpack, 3/11/78.
!>     modified 3/93 to return if incx .le. 0.
!>     modified 12/3/93, array(1) declarations changed to array(*)
!> \endverbatim
!>
!  =====================================================================
      FUNCTION SASUM(N,Sx,Incx)
      IMPLICIT NONE
!*--SASUM76
      REAL :: SASUM
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(IN) , DIMENSION(*) :: Sx
      INTEGER , INTENT(IN) :: Incx
!
! Local variable declarations rewritten by SPAG
!
      INTEGER :: i , m , mp1 , nincx
      REAL :: stemp
!
! End of declarations rewritten by SPAG
!
!     ..
!     .. Array Arguments ..
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
!     ..
!     .. Intrinsic Functions ..
!     ..
      SASUM = 0.0E0
      stemp = 0.0E0
      IF ( N<=0 .OR. Incx<=0 ) RETURN
      IF ( Incx==1 ) THEN
!        code for increment equal to 1
!
!
!        clean-up loop
!
         m = MOD(N,6)
         IF ( m/=0 ) THEN
            DO i = 1 , m
               stemp = stemp + ABS(Sx(i))
            ENDDO
            IF ( N<6 ) THEN
               SASUM = stemp
               RETURN
            ENDIF
         ENDIF
         mp1 = m + 1
         DO i = mp1 , N , 6
            stemp = stemp + ABS(Sx(i)) + ABS(Sx(i+1)) + ABS(Sx(i+2))    &
     &              + ABS(Sx(i+3)) + ABS(Sx(i+4)) + ABS(Sx(i+5))
         ENDDO
      ELSE
!
!        code for increment not equal to 1
!
         nincx = N*Incx
         DO i = 1 , nincx , Incx
            stemp = stemp + ABS(Sx(i))
         ENDDO
      ENDIF
      SASUM = stemp
      END FUNCTION SASUM
