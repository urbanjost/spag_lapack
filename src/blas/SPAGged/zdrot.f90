!*==zdrot.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b ZDROT
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZDROT( N, ZX, INCX, ZY, INCY, C, S )
!
!       .. Scalar Arguments ..
!       INTEGER            INCX, INCY, N
!       DOUBLE PRECISION   C, S
!       ..
!       .. Array Arguments ..
!       COMPLEX*16         ZX( * ), ZY( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> Applies a plane rotation, where the cos and sin (c and s) are real
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
!> \param[in,out] ZX
!> \verbatim
!>          ZX is COMPLEX*16 array, dimension at least
!>           ( 1 + ( N - 1 )*abs( INCX ) ).
!>           Before entry, the incremented array ZX must contain the n
!>           element vector cx. On exit, ZX is overwritten by the updated
!>           vector cx.
!> \endverbatim
!>
!> \param[in] INCX
!> \verbatim
!>          INCX is INTEGER
!>           On entry, INCX specifies the increment for the elements of
!>           ZX. INCX must not be zero.
!> \endverbatim
!>
!> \param[in,out] ZY
!> \verbatim
!>          ZY is COMPLEX*16 array, dimension at least
!>           ( 1 + ( N - 1 )*abs( INCY ) ).
!>           Before entry, the incremented array ZY must contain the n
!>           element vector cy. On exit, ZY is overwritten by the updated
!>           vector cy.
!> \endverbatim
!>
!> \param[in] INCY
!> \verbatim
!>          INCY is INTEGER
!>           On entry, INCY specifies the increment for the elements of
!>           ZY. INCY must not be zero.
!> \endverbatim
!>
!> \param[in] C
!> \verbatim
!>          C is DOUBLE PRECISION
!>           On entry, C specifies the cosine, cos.
!> \endverbatim
!>
!> \param[in] S
!> \verbatim
!>          S is DOUBLE PRECISION
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
!> \ingroup complex16_blas_level1
!
!  =====================================================================
      SUBROUTINE ZDROT(N,Zx,Incx,Zy,Incy,C,S)
      USE F77KINDS                        
      IMPLICIT NONE
!*--ZDROT103
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Zx
      INTEGER , INTENT(IN) :: Incx
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Zy
      INTEGER , INTENT(IN) :: Incy
      REAL(R8KIND) , INTENT(IN) :: C
      REAL(R8KIND) , INTENT(IN) :: S
!
! Local variable declarations rewritten by SPAG
!
      COMPLEX(CX16KIND) :: ctemp
      INTEGER :: i , ix , iy
!
! End of declarations rewritten by SPAG
!
!     ..
!     .. Array Arguments ..
!     ..
!
! =====================================================================
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
            ctemp = C*Zx(i) + S*Zy(i)
            Zy(i) = C*Zy(i) - S*Zx(i)
            Zx(i) = ctemp
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
            ctemp = C*Zx(ix) + S*Zy(iy)
            Zy(iy) = C*Zy(iy) - S*Zx(ix)
            Zx(ix) = ctemp
            ix = ix + Incx
            iy = iy + Incy
         ENDDO
      ENDIF
      END SUBROUTINE ZDROT
