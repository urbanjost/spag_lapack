!*==drot.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
!> \brief \b DROT
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DROT(N,DX,INCX,DY,INCY,C,S)
!
!       .. Scalar Arguments ..
!       DOUBLE PRECISION C,S
!       INTEGER INCX,INCY,N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION DX(*),DY(*)
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    DROT applies a plane rotation.
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
!> \param[in,out] DX
!> \verbatim
!>          DX is DOUBLE PRECISION array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
!> \endverbatim
!>
!> \param[in] INCX
!> \verbatim
!>          INCX is INTEGER
!>         storage spacing between elements of DX
!> \endverbatim
!>
!> \param[in,out] DY
!> \verbatim
!>          DY is DOUBLE PRECISION array, dimension ( 1 + ( N - 1 )*abs( INCY ) )
!> \endverbatim
!>
!> \param[in] INCY
!> \verbatim
!>          INCY is INTEGER
!>         storage spacing between elements of DY
!> \endverbatim
!>
!> \param[in] C
!> \verbatim
!>          C is DOUBLE PRECISION
!> \endverbatim
!>
!> \param[in] S
!> \verbatim
!>          S is DOUBLE PRECISION
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
!> \ingroup double_blas_level1
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
      SUBROUTINE DROT(N,Dx,Incx,Dy,Incy,C,S)
      USE F77KINDS
      IMPLICIT NONE
!*--DROT97
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER , INTENT(IN) :: N
      REAL(r8kind) , INTENT(INOUT) , DIMENSION(*) :: Dx
      INTEGER , INTENT(IN) :: Incx
      REAL(r8kind) , INTENT(INOUT) , DIMENSION(*) :: Dy
      INTEGER , INTENT(IN) :: Incy
      REAL(r8kind) , INTENT(IN) :: C
      REAL(r8kind) , INTENT(IN) :: S
!
! Local variable declarations rewritten by SPAG
!
      REAL(r8kind) :: dtemp
      INTEGER :: i , ix , iy
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
            dtemp = C*Dx(i) + S*Dy(i)
            Dy(i) = C*Dy(i) - S*Dx(i)
            Dx(i) = dtemp
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
            dtemp = C*Dx(ix) + S*Dy(iy)
            Dy(iy) = C*Dy(iy) - S*Dx(ix)
            Dx(ix) = dtemp
            ix = ix + Incx
            iy = iy + Incy
         ENDDO
      ENDIF
      END SUBROUTINE DROT
