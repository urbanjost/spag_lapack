!*==cdotu.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b CDOTU
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       COMPLEX FUNCTION CDOTU(N,CX,INCX,CY,INCY)
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
!> CDOTU forms the dot product of two complex vectors
!>      CDOTU = X^T * Y
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
!> \param[in] CY
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
      FUNCTION CDOTU(N,Cx,Incx,Cy,Incy)
      IMPLICIT NONE
!*--CDOTU87
      COMPLEX :: CDOTU
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(IN) , DIMENSION(*) :: Cx
      INTEGER , INTENT(IN) :: Incx
      COMPLEX , INTENT(IN) , DIMENSION(*) :: Cy
      INTEGER , INTENT(IN) :: Incy
!
! Local variable declarations rewritten by SPAG
!
      COMPLEX :: ctemp
      INTEGER :: i , ix , iy
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
      ctemp = (0.0,0.0)
      CDOTU = (0.0,0.0)
      IF ( N<=0 ) RETURN
      IF ( Incx==1 .AND. Incy==1 ) THEN
!
!        code for both increments equal to 1
!
         DO i = 1 , N
            ctemp = ctemp + Cx(i)*Cy(i)
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
            ctemp = ctemp + Cx(ix)*Cy(iy)
            ix = ix + Incx
            iy = iy + Incy
         ENDDO
      ENDIF
      CDOTU = ctemp
      END FUNCTION CDOTU
