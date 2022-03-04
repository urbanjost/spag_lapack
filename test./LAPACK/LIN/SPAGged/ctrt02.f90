!*==ctrt02.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b CTRT02
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE CTRT02( UPLO, TRANS, DIAG, N, NRHS, A, LDA, X, LDX, B,
!                          LDB, WORK, RWORK, RESID )
!
!       .. Scalar Arguments ..
!       CHARACTER          DIAG, TRANS, UPLO
!       INTEGER            LDA, LDB, LDX, N, NRHS
!       REAL               RESID
!       ..
!       .. Array Arguments ..
!       REAL               RWORK( * )
!       COMPLEX            A( LDA, * ), B( LDB, * ), WORK( * ),
!      $                   X( LDX, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CTRT02 computes the residual for the computed solution to a
!> triangular system of linear equations  A*x = b,  A**T *x = b,
!> or A**H *x = b.  Here A is a triangular matrix, A**T is the transpose
!> of A, A**H is the conjugate transpose of A, and x and b are N by NRHS
!> matrices.  The test ratio is the maximum over the number of right
!> hand sides of
!>    norm(b - op(A)*x) / ( norm(op(A)) * norm(x) * EPS ),
!> where op(A) denotes A, A**T, or A**H, and EPS is the machine epsilon.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          Specifies whether the matrix A is upper or lower triangular.
!>          = 'U':  Upper triangular
!>          = 'L':  Lower triangular
!> \endverbatim
!>
!> \param[in] TRANS
!> \verbatim
!>          TRANS is CHARACTER*1
!>          Specifies the operation applied to A.
!>          = 'N':  A *x = b     (No transpose)
!>          = 'T':  A**T *x = b  (Transpose)
!>          = 'C':  A**H *x = b  (Conjugate transpose)
!> \endverbatim
!>
!> \param[in] DIAG
!> \verbatim
!>          DIAG is CHARACTER*1
!>          Specifies whether or not the matrix A is unit triangular.
!>          = 'N':  Non-unit triangular
!>          = 'U':  Unit triangular
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in] NRHS
!> \verbatim
!>          NRHS is INTEGER
!>          The number of right hand sides, i.e., the number of columns
!>          of the matrices X and B.  NRHS >= 0.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is COMPLEX array, dimension (LDA,N)
!>          The triangular matrix A.  If UPLO = 'U', the leading n by n
!>          upper triangular part of the array A contains the upper
!>          triangular matrix, and the strictly lower triangular part of
!>          A is not referenced.  If UPLO = 'L', the leading n by n lower
!>          triangular part of the array A contains the lower triangular
!>          matrix, and the strictly upper triangular part of A is not
!>          referenced.  If DIAG = 'U', the diagonal elements of A are
!>          also not referenced and are assumed to be 1.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,N).
!> \endverbatim
!>
!> \param[in] X
!> \verbatim
!>          X is COMPLEX array, dimension (LDX,NRHS)
!>          The computed solution vectors for the system of linear
!>          equations.
!> \endverbatim
!>
!> \param[in] LDX
!> \verbatim
!>          LDX is INTEGER
!>          The leading dimension of the array X.  LDX >= max(1,N).
!> \endverbatim
!>
!> \param[in] B
!> \verbatim
!>          B is COMPLEX array, dimension (LDB,NRHS)
!>          The right hand side vectors for the system of linear
!>          equations.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of the array B.  LDB >= max(1,N).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX array, dimension (N)
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is REAL array, dimension (N)
!> \endverbatim
!>
!> \param[out] RESID
!> \verbatim
!>          RESID is REAL
!>          The maximum over the number of right hand sides of
!>          norm(op(A)*x - b) / ( norm(op(A)) * norm(x) * EPS ).
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
      SUBROUTINE CTRT02(Uplo,Trans,Diag,N,Nrhs,A,Lda,X,Ldx,B,Ldb,Work,  &
     &                  Rwork,Resid)
      IMPLICIT NONE
!*--CTRT02161
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Diag , Trans , Uplo
      INTEGER Lda , Ldb , Ldx , N , Nrhs
      REAL Resid
!     ..
!     .. Array Arguments ..
      REAL Rwork(*)
      COMPLEX A(Lda,*) , B(Ldb,*) , Work(*) , X(Ldx,*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ZERO , ONE
      PARAMETER (ZERO=0.0E+0,ONE=1.0E+0)
!     ..
!     .. Local Scalars ..
      INTEGER j
      REAL anorm , bnorm , eps , xnorm
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      REAL CLANTR , SCASUM , SLAMCH
      EXTERNAL LSAME , CLANTR , SCASUM , SLAMCH
!     ..
!     .. External Subroutines ..
      EXTERNAL CAXPY , CCOPY , CTRMV
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC CMPLX , MAX
!     ..
!     .. Executable Statements ..
!
!     Quick exit if N = 0 or NRHS = 0
!
      IF ( N<=0 .OR. Nrhs<=0 ) THEN
         Resid = ZERO
         RETURN
      ENDIF
!
!     Compute the 1-norm of A or A**H.
!
      IF ( LSAME(Trans,'N') ) THEN
         anorm = CLANTR('1',Uplo,Diag,N,N,A,Lda,Rwork)
      ELSE
         anorm = CLANTR('I',Uplo,Diag,N,N,A,Lda,Rwork)
      ENDIF
!
!     Exit with RESID = 1/EPS if ANORM = 0.
!
      eps = SLAMCH('Epsilon')
      IF ( anorm<=ZERO ) THEN
         Resid = ONE/eps
         RETURN
      ENDIF
!
!     Compute the maximum over the number of right hand sides of
!        norm(op(A)*x - b) / ( norm(op(A)) * norm(x) * EPS )
!
      Resid = ZERO
      DO j = 1 , Nrhs
         CALL CCOPY(N,X(1,j),1,Work,1)
         CALL CTRMV(Uplo,Trans,Diag,N,A,Lda,Work,1)
         CALL CAXPY(N,CMPLX(-ONE),B(1,j),1,Work,1)
         bnorm = SCASUM(N,Work,1)
         xnorm = SCASUM(N,X(1,j),1)
         IF ( xnorm<=ZERO ) THEN
            Resid = ONE/eps
         ELSE
            Resid = MAX(Resid,((bnorm/anorm)/xnorm)/eps)
         ENDIF
      ENDDO
!
!
!     End of CTRT02
!
      END SUBROUTINE CTRT02
