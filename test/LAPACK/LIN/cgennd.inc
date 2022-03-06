!*==cgennd.f90  processed by SPAG 7.51RB at 17:34 on  4 Mar 2022
!> \brief \b cgennd
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       LOGICAL FUNCTION cgennd (M, N, A, LDA)
!
!       .. Scalar Arguments ..
!       INTEGER M, N, LDA
!       ..
!       .. Array Arguments ..
!       COMPLEX A( LDA, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    cgennd tests that its argument has a real, non-negative diagonal.
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
!>          A is COMPLEX array, dimension (LDA, N)
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
!> \ingroup complex_lin
!
!  =====================================================================
      LOGICAL FUNCTION cgennd(M,N,A,Lda)
      IMPLICIT NONE
!*--CGENND72
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
      COMPLEX A(Lda,*)
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
      COMPLEX aii
!     ..
!     .. Intrinsics ..
      INTRINSIC MIN , REAL , AIMAG
!     ..
!     .. Executable Statements ..
      k = MIN(M,N)
      DO i = 1 , k
         aii = A(i,i)
         IF ( REAL(aii)<ZERO .OR. AIMAG(aii)/=ZERO ) THEN
            cgennd = .FALSE.
            RETURN
         ENDIF
      ENDDO
      cgennd = .TRUE.
      END FUNCTION cgennd
