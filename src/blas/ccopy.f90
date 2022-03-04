!*==ccopy.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b CCOPY
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE CCOPY(N,CX,INCX,CY,INCY)
!
!       .. Scalar Arguments ..
!       INTEGER INCX,INCY,N
!       ..
!       .. Array Arguments ..
!       COMPLEX CX(*),CY(*)
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    CCOPY copies a vector x to a vector y.
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
!> \param[in] CX
!> \verbatim
!>          CX is COMPLEX array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
!> \endverbatim
!>
!> \param[in] INCX
!> \verbatim
!>          INCX is INTEGER
!>         storage spacing between elements of CX
!> \endverbatim
!>
!> \param[out] CY
!> \verbatim
!>          CY is COMPLEX array, dimension ( 1 + ( N - 1 )*abs( INCY ) )
!> \endverbatim
!>
!> \param[in] INCY
!> \verbatim
!>          INCY is INTEGER
!>         storage spacing between elements of CY
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
!> \ingroup complex_blas_level1
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
      SUBROUTINE CCOPY(N,Cx,Incx,Cy,Incy)
      IMPLICIT NONE
!*--CCOPY85
!
!  -- Reference BLAS level1 routine (version 3.8.0) --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2017
!
!     .. Scalar Arguments ..
      INTEGER Incx , Incy , N
!     ..
!     .. Array Arguments ..
      COMPLEX Cx(*) , Cy(*)
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
      INTEGER i , ix , iy
!     ..
      IF ( N<=0 ) RETURN
      IF ( Incx==1 .AND. Incy==1 ) THEN
!
!        code for both increments equal to 1
!
         DO i = 1 , N
            Cy(i) = Cx(i)
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
            Cy(iy) = Cx(ix)
            ix = ix + Incx
            iy = iy + Incy
         ENDDO
      ENDIF
      END SUBROUTINE CCOPY
