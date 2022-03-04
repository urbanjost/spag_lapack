!*==zget10.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b ZGET10
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZGET10( M, N, A, LDA, B, LDB, WORK, RWORK, RESULT )
!
!       .. Scalar Arguments ..
!       INTEGER            LDA, LDB, M, N
!       DOUBLE PRECISION   RESULT
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   RWORK( * )
!       COMPLEX*16         A( LDA, * ), B( LDB, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZGET10 compares two matrices A and B and computes the ratio
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
!>          A is COMPLEX*16 array, dimension (LDA,N)
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
!>          B is COMPLEX*16 array, dimension (LDB,N)
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
!>          WORK is COMPLEX*16 array, dimension (M)
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is COMPLEX*16 array, dimension (M)
!> \endverbatim
!>
!> \param[out] RESULT
!> \verbatim
!>          RESULT is DOUBLE PRECISION
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
!> \ingroup complex16_eig
!
!  =====================================================================
      SUBROUTINE ZGET10(M,N,A,Lda,B,Ldb,Work,Rwork,Result)
      IMPLICIT NONE
!*--ZGET10103
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Lda , Ldb , M , N
      DOUBLE PRECISION Result
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION Rwork(*)
      COMPLEX*16 A(Lda,*) , B(Ldb,*) , Work(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ONE , ZERO
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
!     ..
!     .. Local Scalars ..
      INTEGER j
      DOUBLE PRECISION anorm , eps , unfl , wnorm
!     ..
!     .. External Functions ..
      DOUBLE PRECISION DLAMCH , DZASUM , ZLANGE
      EXTERNAL DLAMCH , DZASUM , ZLANGE
!     ..
!     .. External Subroutines ..
      EXTERNAL ZAXPY , ZCOPY
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC DBLE , DCMPLX , MAX , MIN
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
      unfl = DLAMCH('Safe minimum')
      eps = DLAMCH('Precision')
!
      wnorm = ZERO
      DO j = 1 , N
         CALL ZCOPY(M,A(1,j),1,Work,1)
         CALL ZAXPY(M,DCMPLX(-ONE),B(1,j),1,Work,1)
         wnorm = MAX(wnorm,DZASUM(N,Work,1))
      ENDDO
!
      anorm = MAX(ZLANGE('1',M,N,A,Lda,Rwork),unfl)
!
      IF ( anorm>wnorm ) THEN
         Result = (wnorm/anorm)/(M*eps)
      ELSEIF ( anorm<ONE ) THEN
         Result = (MIN(wnorm,M*anorm)/anorm)/(M*eps)
      ELSE
         Result = MIN(wnorm/anorm,DBLE(M))/(M*eps)
      ENDIF
!
!
!     End of ZGET10
!
      END SUBROUTINE ZGET10
