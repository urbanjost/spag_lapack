!*==zdotc.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b ZDOTC
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       COMPLEX*16 FUNCTION ZDOTC(N,ZX,INCX,ZY,INCY)
!
!       .. Scalar Arguments ..
!       INTEGER INCX,INCY,N
!       ..
!       .. Array Arguments ..
!       COMPLEX*16 ZX(*),ZY(*)
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZDOTC forms the dot product of two complex vectors
!>      ZDOTC = X^H * Y
!>
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
!> \param[in] ZX
!> \verbatim
!>          ZX is COMPLEX*16 array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
!> \endverbatim
!>
!> \param[in] INCX
!> \verbatim
!>          INCX is INTEGER
!>         storage spacing between elements of ZX
!> \endverbatim
!>
!> \param[in] ZY
!> \verbatim
!>          ZY is COMPLEX*16 array, dimension ( 1 + ( N - 1 )*abs( INCY ) )
!> \endverbatim
!>
!> \param[in] INCY
!> \verbatim
!>          INCY is INTEGER
!>         storage spacing between elements of ZY
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
!> \ingroup complex16_blas_level1
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>     jack dongarra, 3/11/78.
!>     modified 12/3/93, array(1) declarations changed to array(*)
!> \endverbatim
!>
!  =====================================================================
      FUNCTION ZDOTC(N,Zx,Incx,Zy,Incy)
      USE F77KINDS                        
      IMPLICIT NONE
!*--ZDOTC88
      COMPLEX(CX16KIND) :: ZDOTC
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(*) :: Zx
      INTEGER , INTENT(IN) :: Incx
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(*) :: Zy
      INTEGER , INTENT(IN) :: Incy
!
! Local variable declarations rewritten by SPAG
!
      INTEGER :: i , ix , iy
      COMPLEX(CX16KIND) :: ztemp
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
      ztemp = (0.0D0,0.0D0)
      ZDOTC = (0.0D0,0.0D0)
      IF ( N<=0 ) RETURN
      IF ( Incx==1 .AND. Incy==1 ) THEN
!
!        code for both increments equal to 1
!
         DO i = 1 , N
            ztemp = ztemp + DCONJG(Zx(i))*Zy(i)
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
            ztemp = ztemp + DCONJG(Zx(ix))*Zy(iy)
            ix = ix + Incx
            iy = iy + Incy
         ENDDO
      ENDIF
      ZDOTC = ztemp
      END FUNCTION ZDOTC
