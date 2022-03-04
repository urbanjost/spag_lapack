!*==sgennd.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b SGENND
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       LOGICAL FUNCTION SGENND (M, N, A, LDA)
!
!       .. Scalar Arguments ..
!       INTEGER M, N, LDA
!       ..
!       .. Array Arguments ..
!       REAL A( LDA, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    SGENND tests that its argument has a non-negative diagonal.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows in A.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns in A.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is REAL array, dimension (LDA, N)
!>          The matrix.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          Leading dimension of A.
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
!> \ingroup single_lin
!
!  =====================================================================
      LOGICAL FUNCTION SGENND(M,N,A,Lda)
      IMPLICIT NONE
!*--SGENND72
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER M , N , Lda
!     ..
!     .. Array Arguments ..
      REAL A(Lda,*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ZERO
      PARAMETER (ZERO=0.0E0)
!     ..
!     .. Local Scalars ..
      INTEGER i , k
!     ..
!     .. Intrinsics ..
      INTRINSIC MIN
!     ..
!     .. Executable Statements ..
      k = MIN(M,N)
      DO i = 1 , k
         IF ( A(i,i)<ZERO ) THEN
            SGENND = .FALSE.
            RETURN
         ENDIF
      ENDDO
      SGENND = .TRUE.
      END FUNCTION SGENND
