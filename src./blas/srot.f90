!*==srot.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
!> \brief \b SROT
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE SROT(N,SX,INCX,SY,INCY,C,S)
!
!       .. Scalar Arguments ..
!       REAL C,S
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
!>    applies a plane rotation.
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
!> \param[in,out] SX
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
!> \param[in,out] SY
!> \verbatim
!>          SY is REAL array, dimension ( 1 + ( N - 1 )*abs( INCY ) )
!> \endverbatim
!>
!> \param[in] INCY
!> \verbatim
!>          INCY is INTEGER
!>         storage spacing between elements of SY
!> \endverbatim
!>
!> \param[in] C
!> \verbatim
!>          C is REAL
!> \endverbatim
!>
!> \param[in] S
!> \verbatim
!>          S is REAL
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
      SUBROUTINE SROT(N,Sx,Incx,Sy,Incy,C,S)
      IMPLICIT NONE
!*--SROT96
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(INOUT) , DIMENSION(*) :: Sx
      INTEGER , INTENT(IN) :: Incx
      REAL , INTENT(INOUT) , DIMENSION(*) :: Sy
      INTEGER , INTENT(IN) :: Incy
      REAL , INTENT(IN) :: C
      REAL , INTENT(IN) :: S
!
! Local variable declarations rewritten by SPAG
!
      INTEGER :: i , ix , iy
      REAL :: stemp
!
! End of declarations rewritten by SPAG
!
!
! Local variable declarations rewritten by SPAG
!
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
      IF ( N<=0 ) RETURN
      IF ( Incx==1 .AND. Incy==1 ) THEN
!
!       code for both increments equal to 1
!
         DO i = 1 , N
            stemp = C*Sx(i) + S*Sy(i)
            Sy(i) = C*Sy(i) - S*Sx(i)
            Sx(i) = stemp
         ENDDO
      ELSE
!
!       code for unequal increments or equal increments not equal
!         to 1
!
         ix = 1
         iy = 1
         IF ( Incx<0 ) ix = (-N+1)*Incx + 1
         IF ( Incy<0 ) iy = (-N+1)*Incy + 1
         DO i = 1 , N
            stemp = C*Sx(ix) + S*Sy(iy)
            Sy(iy) = C*Sy(iy) - S*Sx(ix)
            Sx(ix) = stemp
            ix = ix + Incx
            iy = iy + Incy
         ENDDO
      ENDIF
      END SUBROUTINE SROT
