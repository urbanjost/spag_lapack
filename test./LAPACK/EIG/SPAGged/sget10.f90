!*==sget10.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b SGET10
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE SGET10( M, N, A, LDA, B, LDB, WORK, RESULT )
!
!       .. Scalar Arguments ..
!       INTEGER            LDA, LDB, M, N
!       REAL               RESULT
!       ..
!       .. Array Arguments ..
!       REAL               A( LDA, * ), B( LDB, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SGET10 compares two matrices A and B and computes the ratio
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
!>          A is REAL array, dimension (LDA,N)
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
!>          B is REAL array, dimension (LDB,N)
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
!>          WORK is REAL array, dimension (M)
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
!> \ingroup single_eig
!
!  =====================================================================
      SUBROUTINE SGET10(M,N,A,Lda,B,Ldb,Work,Result)
      IMPLICIT NONE
!*--SGET1097
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
      REAL A(Lda,*) , B(Ldb,*) , Work(*)
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
      REAL SASUM , SLAMCH , SLANGE
      EXTERNAL SASUM , SLAMCH , SLANGE
!     ..
!     .. External Subroutines ..
      EXTERNAL SAXPY , SCOPY
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
         CALL SCOPY(M,A(1,j),1,Work,1)
         CALL SAXPY(M,-ONE,B(1,j),1,Work,1)
         wnorm = MAX(wnorm,SASUM(N,Work,1))
      ENDDO
!
      anorm = MAX(SLANGE('1',M,N,A,Lda,Work),unfl)
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
!     End of SGET10
!
      END SUBROUTINE SGET10
