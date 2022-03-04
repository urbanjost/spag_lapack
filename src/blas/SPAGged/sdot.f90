!*==sdot.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b SDOT
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       REAL FUNCTION SDOT(N,SX,INCX,SY,INCY)
!
!       .. Scalar Arguments ..
!       INTEGER INCX,INCY,N
!       ..
!       .. Array Arguments ..
!       REAL SX(*),SY(*)
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    SDOT forms the dot product of two vectors.
!>    uses unrolled loops for increments equal to one.
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
!>
!> \param[in] SY
!> \verbatim
!>          SY is REAL array, dimension ( 1 + ( N - 1 )*abs( INCY ) )
!> \endverbatim
!>
!> \param[in] INCY
!> \verbatim
!>          INCY is INTEGER
!>         storage spacing between elements of SY
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
!>     modified 12/3/93, array(1) declarations changed to array(*)
!> \endverbatim
!>
!  =====================================================================
      FUNCTION SDOT(N,Sx,Incx,Sy,Incy)
      IMPLICIT NONE
!*--SDOT86
      REAL :: SDOT
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(IN) , DIMENSION(*) :: Sx
      INTEGER , INTENT(IN) :: Incx
      REAL , INTENT(IN) , DIMENSION(*) :: Sy
      INTEGER , INTENT(IN) :: Incy
!
! Local variable declarations rewritten by SPAG
!
      INTEGER :: i , ix , iy , m , mp1
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
      stemp = 0.0E0
      SDOT = 0.0E0
      IF ( N<=0 ) RETURN
      IF ( Incx==1 .AND. Incy==1 ) THEN
!
!        code for both increments equal to 1
!
!
!        clean-up loop
!
         m = MOD(N,5)
         IF ( m/=0 ) THEN
            DO i = 1 , m
               stemp = stemp + Sx(i)*Sy(i)
            ENDDO
            IF ( N<5 ) THEN
               SDOT = stemp
               RETURN
            ENDIF
         ENDIF
         mp1 = m + 1
         DO i = mp1 , N , 5
            stemp = stemp + Sx(i)*Sy(i) + Sx(i+1)*Sy(i+1) + Sx(i+2)     &
     &              *Sy(i+2) + Sx(i+3)*Sy(i+3) + Sx(i+4)*Sy(i+4)
         ENDDO
      ELSE
!
!        code for unequal increments or equal increments
!          not equal to 1
!
         ix = 1
         iy = 1
         IF ( Incx<0 ) ix = (-N+1)*Incx + 1
         IF ( Incy<0 ) iy = (-N+1)*Incy + 1
         DO i = 1 , N
            stemp = stemp + Sx(ix)*Sy(iy)
            ix = ix + Incx
            iy = iy + Incy
         ENDDO
      ENDIF
      SDOT = stemp
      END FUNCTION SDOT
