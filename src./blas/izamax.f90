!*==izamax.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
!> \brief \b IZAMAX
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       INTEGER FUNCTION IZAMAX(N,ZX,INCX)
!
!       .. Scalar Arguments ..
!       INTEGER INCX,N
!       ..
!       .. Array Arguments ..
!       COMPLEX*16 ZX(*)
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    IZAMAX finds the index of the first element having maximum |Re(.)| + |Im(.)|
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
!> \ingroup aux_blas
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>     jack dongarra, 1/15/85.
!>     modified 3/93 to return if incx .le. 0.
!>     modified 12/3/93, array(1) declarations changed to array(*)
!> \endverbatim
!>
!  =====================================================================
      FUNCTION IZAMAX(N,Zx,Incx)
      USE F77KINDS
      USE S_DCABS1
      IMPLICIT NONE
!*--IZAMAX77
      INTEGER :: IZAMAX
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER , INTENT(IN) :: N
      COMPLEX(cx16kind) , INTENT(IN) , DIMENSION(*) :: Zx
      INTEGER , INTENT(IN) :: Incx
!
! Local variable declarations rewritten by SPAG
!
      REAL(r8kind) :: dmax
      INTEGER :: i , ix
!
! End of declarations rewritten by SPAG
!
!
! Dummy argument declarations rewritten by SPAG
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
!     .. External Functions ..
!     ..
      IZAMAX = 0
      IF ( N<1 .OR. Incx<=0 ) RETURN
      IZAMAX = 1
      IF ( N==1 ) RETURN
      IF ( Incx==1 ) THEN
!
!        code for increment equal to 1
!
         dmax = DCABS1(Zx(1))
         DO i = 2 , N
            IF ( DCABS1(Zx(i))>dmax ) THEN
               IZAMAX = i
               dmax = DCABS1(Zx(i))
            ENDIF
         ENDDO
      ELSE
!
!        code for increment not equal to 1
!
         ix = 1
         dmax = DCABS1(Zx(1))
         ix = ix + Incx
         DO i = 2 , N
            IF ( DCABS1(Zx(ix))>dmax ) THEN
               IZAMAX = i
               dmax = DCABS1(Zx(ix))
            ENDIF
            ix = ix + Incx
         ENDDO
      ENDIF
      END FUNCTION IZAMAX
