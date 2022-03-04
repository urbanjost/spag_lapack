!*==clarfy.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b CLARFY
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE CLARFY( UPLO, N, V, INCV, TAU, C, LDC, WORK )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            INCV, LDC, N
!       COMPLEX            TAU
!       ..
!       .. Array Arguments ..
!       COMPLEX            C( LDC, * ), V( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CLARFY applies an elementary reflector, or Householder matrix, H,
!> to an n x n Hermitian matrix C, from both the left and the right.
!>
!> H is represented in the form
!>
!>    H = I - tau * v * v'
!>
!> where  tau  is a scalar and  v  is a vector.
!>
!> If  tau  is  zero, then  H  is taken to be the unit matrix.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          Specifies whether the upper or lower triangular part of the
!>          Hermitian matrix C is stored.
!>          = 'U':  Upper triangle
!>          = 'L':  Lower triangle
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of rows and columns of the matrix C.  N >= 0.
!> \endverbatim
!>
!> \param[in] V
!> \verbatim
!>          V is COMPLEX array, dimension
!>                  (1 + (N-1)*abs(INCV))
!>          The vector v as described above.
!> \endverbatim
!>
!> \param[in] INCV
!> \verbatim
!>          INCV is INTEGER
!>          The increment between successive elements of v.  INCV must
!>          not be zero.
!> \endverbatim
!>
!> \param[in] TAU
!> \verbatim
!>          TAU is COMPLEX
!>          The value tau as described above.
!> \endverbatim
!>
!> \param[in,out] C
!> \verbatim
!>          C is COMPLEX array, dimension (LDC, N)
!>          On entry, the matrix C.
!>          On exit, C is overwritten by H * C * H'.
!> \endverbatim
!>
!> \param[in] LDC
!> \verbatim
!>          LDC is INTEGER
!>          The leading dimension of the array C.  LDC >= max( 1, N ).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX array, dimension (N)
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
!> \ingroup complexOTHERauxiliary
!
!  =====================================================================
      SUBROUTINE CLARFY(Uplo,N,V,Incv,Tau,C,Ldc,Work)
      USE S_CAXPY
      USE S_CDOTC
      USE S_CHEMV
      USE S_CHER2
      IMPLICIT NONE
!*--CLARFY116
!
! PARAMETER definitions rewritten by SPAG
!
      COMPLEX , PARAMETER  ::  ONE = (1.0E+0,0.0E+0) ,                  &
     &                         ZERO = (0.0E+0,0.0E+0) ,                 &
     &                         HALF = (0.5E+0,0.0E+0)
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX , DIMENSION(*) :: V
      INTEGER :: Incv
      COMPLEX :: Tau
      COMPLEX , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      COMPLEX , DIMENSION(*) :: Work
!
! Local variable declarations rewritten by SPAG
!
      COMPLEX :: alpha
!
! End of declarations rewritten by SPAG
!
!     ..
!     .. Array Arguments ..
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
!     ..
!     .. Local Scalars ..
!     ..
!     .. External Subroutines ..
!     ..
!     .. External Functions ..
!     ..
!     .. Executable Statements ..
!
      IF ( Tau==ZERO ) RETURN
!
!     Form  w:= C * v
!
      CALL CHEMV(Uplo,N,ONE,C,Ldc,V,Incv,ZERO,Work,1)
!
      alpha = -HALF*Tau*CDOTC(N,Work,1,V,Incv)
      CALL CAXPY(N,alpha,V,Incv,Work,1)
!
!     C := C - v * w' - w * v'
!
      CALL CHER2(Uplo,N,-Tau,V,Incv,Work,1,C,Ldc)
!
!
!     End of CLARFY
!
      END SUBROUTINE CLARFY
