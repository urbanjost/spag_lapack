!*==csrot.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b CSROT
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE CSROT( N, CX, INCX, CY, INCY, C, S )
!
!       .. Scalar Arguments ..
!       INTEGER           INCX, INCY, N
!       REAL              C, S
!       ..
!       .. Array Arguments ..
!       COMPLEX           CX( * ), CY( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CSROT applies a plane rotation, where the cos and sin (c and s) are real
!> and the vectors cx and cy are complex.
!> jack dongarra, linpack, 3/11/78.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>           On entry, N specifies the order of the vectors cx and cy.
!>           N must be at least zero.
!> \endverbatim
!>
!> \param[in,out] CX
!> \verbatim
!>          CX is COMPLEX array, dimension at least
!>           ( 1 + ( N - 1 )*abs( INCX ) ).
!>           Before entry, the incremented array CX must contain the n
!>           element vector cx. On exit, CX is overwritten by the updated
!>           vector cx.
!> \endverbatim
!>
!> \param[in] INCX
!> \verbatim
!>          INCX is INTEGER
!>           On entry, INCX specifies the increment for the elements of
!>           CX. INCX must not be zero.
!> \endverbatim
!>
!> \param[in,out] CY
!> \verbatim
!>          CY is COMPLEX array, dimension at least
!>           ( 1 + ( N - 1 )*abs( INCY ) ).
!>           Before entry, the incremented array CY must contain the n
!>           element vector cy. On exit, CY is overwritten by the updated
!>           vector cy.
!> \endverbatim
!>
!> \param[in] INCY
!> \verbatim
!>          INCY is INTEGER
!>           On entry, INCY specifies the increment for the elements of
!>           CY. INCY must not be zero.
!> \endverbatim
!>
!> \param[in] C
!> \verbatim
!>          C is REAL
!>           On entry, C specifies the cosine, cos.
!> \endverbatim
!>
!> \param[in] S
!> \verbatim
!>          S is REAL
!>           On entry, S specifies the sine, sin.
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
!> \date December 2016
!
!> \ingroup complex_blas_level1
!
!  =====================================================================
      SUBROUTINE CSROT(N,Cx,Incx,Cy,Incy,C,S)
      IMPLICIT NONE
!*--CSROT102
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Cx
      INTEGER , INTENT(IN) :: Incx
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Cy
      INTEGER , INTENT(IN) :: Incy
      REAL , INTENT(IN) :: C
      REAL , INTENT(IN) :: S
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
!     .. Executable Statements ..
!
      IF ( N<=0 ) RETURN
      IF ( Incx==1 .AND. Incy==1 ) THEN
!
!        code for both increments equal to 1
!
         DO i = 1 , N
            ctemp = C*Cx(i) + S*Cy(i)
            Cy(i) = C*Cy(i) - S*Cx(i)
            Cx(i) = ctemp
         ENDDO
      ELSE
!
!        code for unequal increments or equal increments not equal
!          to 1
!
         ix = 1
         iy = 1
         IF ( Incx<0 ) ix = (-N+1)*Incx + 1
         IF ( Incy<0 ) iy = (-N+1)*Incy + 1
         DO i = 1 , N
            ctemp = C*Cx(ix) + S*Cy(iy)
            Cy(iy) = C*Cy(iy) - S*Cx(ix)
            Cx(ix) = ctemp
            ix = ix + Incx
            iy = iy + Incy
         ENDDO
      ENDIF
      END SUBROUTINE CSROT
