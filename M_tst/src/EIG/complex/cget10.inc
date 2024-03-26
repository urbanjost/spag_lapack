!*==cget10.f90  processed by SPAG 7.51RB at 17:37 on  4 Mar 2022
!> \brief \b cget10
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE cget10( M, N, A, LDA, B, LDB, WORK, RWORK, RESULT )
!
!       .. Scalar Arguments ..
!       INTEGER            LDA, LDB, M, N
!       REAL               RESULT
!       ..
!       .. Array Arguments ..
!       REAL               RWORK( * )
!       COMPLEX            A( LDA, * ), B( LDB, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> cget10 compares two matrices A and B and computes the ratio
!> RESULT = norm( A - B ) / ( norm(A) * M * EPS )
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrices A and B.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrices A and B.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is COMPLEX array, dimension (LDA,N)
!>          The m by n matrix A.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,M).
!> \endverbatim
!>
!> \param[in] B
!> \verbatim
!>          B is COMPLEX array, dimension (LDB,N)
!>          The m by n matrix B.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of the array B.  LDB >= max(1,M).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX array, dimension (M)
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is COMPLEX array, dimension (M)
!> \endverbatim
!>
!> \param[out] RESULT
!> \verbatim
!>          RESULT is REAL
!>          RESULT = norm( A - B ) / ( norm(A) * M * EPS )
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
!> \ingroup complex_eig
!
!  =====================================================================
      SUBROUTINE cget10(M,N,A,Lda,B,Ldb,Work,Rwork,Result)
      IMPLICIT NONE
!*--CGET10103
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Lda , Ldb , M , N
      REAL Result
!     ..
!     .. Array Arguments ..
      REAL Rwork(*)
      COMPLEX A(Lda,*) , B(Ldb,*) , Work(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ONE , ZERO
      PARAMETER (ONE=1.0E+0,ZERO=0.0E+0)
!     ..
!     .. Local Scalars ..
      INTEGER j
      REAL anorm , eps , unfl , wnorm
!     ..
!     .. External Functions ..
      REAL SCASUM , SLAMCH , CLANGE
      EXTERNAL SCASUM , SLAMCH , CLANGE
!     ..
!     .. External Subroutines ..
      EXTERNAL CAXPY , CCOPY
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX , MIN , REAL
!     ..
!     .. Executable Statements ..
!
!     Quick return if possible
!
      IF ( M<=0 .OR. N<=0 ) THEN
         Result = ZERO
         RETURN
      ENDIF
!
      unfl = SLAMCH('Safe minimum')
      eps = SLAMCH('Precision')
!
      wnorm = ZERO
      DO j = 1 , N
         CALL CCOPY(M,A(1,j),1,Work,1)
         CALL CAXPY(M,CMPLX(-ONE),B(1,j),1,Work,1)
         wnorm = MAX(wnorm,SCASUM(N,Work,1))
      ENDDO
!
      anorm = MAX(CLANGE('1',M,N,A,Lda,Rwork),unfl)
!
      IF ( anorm>wnorm ) THEN
         Result = (wnorm/anorm)/(M*eps)
      ELSEIF ( anorm<ONE ) THEN
         Result = (MIN(wnorm,M*anorm)/anorm)/(M*eps)
      ELSE
         Result = MIN(wnorm/anorm,REAL(M))/(M*eps)
      ENDIF
!
!
!     End of cget10
!
      END SUBROUTINE cget10
