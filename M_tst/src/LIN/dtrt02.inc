!*==dtrt02.f90  processed by SPAG 7.51RB at 17:34 on  4 Mar 2022
!> \brief \b dtrt02
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE dtrt02( UPLO, TRANS, DIAG, N, NRHS, A, LDA, X, LDX, B,
!                          LDB, WORK, RESID )
!
!       .. Scalar Arguments ..
!       CHARACTER          DIAG, TRANS, UPLO
!       INTEGER            LDA, LDB, LDX, N, NRHS
!       DOUBLE PRECISION   RESID
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), WORK( * ),
!      $                   X( LDX, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> dtrt02 computes the residual for the computed solution to a
!> triangular system of linear equations  A*x = b  or  A'*x = b.
!> Here A is a triangular matrix, A' is the transpose of A, and x and b
!> are N by NRHS matrices.  The test ratio is the maximum over the
!> number of right hand sides of
!>    norm(b - op(A)*x) / ( norm(op(A)) * norm(x) * EPS ),
!> where op(A) denotes A or A' and EPS is the machine epsilon.
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
!>          = 'N':  A *x = b  (No transpose)
!>          = 'T':  A'*x = b  (Transpose)
!>          = 'C':  A'*x = b  (Conjugate transpose = Transpose)
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
!>          A is DOUBLE PRECISION array, dimension (LDA,N)
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
!>          X is DOUBLE PRECISION array, dimension (LDX,NRHS)
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
!>          B is DOUBLE PRECISION array, dimension (LDB,NRHS)
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
!>          WORK is DOUBLE PRECISION array, dimension (N)
!> \endverbatim
!>
!> \param[out] RESID
!> \verbatim
!>          RESID is DOUBLE PRECISION
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
!> \ingroup double_lin
!
!  =====================================================================
      SUBROUTINE dtrt02(Uplo,Trans,Diag,N,Nrhs,A,Lda,X,Ldx,B,Ldb,Work,  &
     &                  Resid)
      IMPLICIT NONE
!*--DTRT02154
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Diag , Trans , Uplo
      INTEGER Lda , Ldb , Ldx , N , Nrhs
      DOUBLE PRECISION Resid
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION A(Lda,*) , B(Ldb,*) , Work(*) , X(Ldx,*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO , ONE
      PARAMETER (ZERO=0.0D+0,ONE=1.0D+0)
!     ..
!     .. Local Scalars ..
      INTEGER j
      DOUBLE PRECISION anorm , bnorm , eps , xnorm
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      DOUBLE PRECISION DASUM , DLAMCH , DLANTR
EXTERNAL    LSAME,   DASUM,   DLAMCH,   DLANTR
!     ..
!     .. External Subroutines ..
EXTERNAL    DAXPY,   DCOPY,   DTRMV
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX
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
!     Compute the 1-norm of A or A'.
!
      IF ( LSAME(Trans,'N') ) THEN
         anorm = DLANTR('1',Uplo,Diag,N,N,A,Lda,Work)
      ELSE
         anorm = DLANTR('I',Uplo,Diag,N,N,A,Lda,Work)
      ENDIF
!
!     Exit with RESID = 1/EPS if ANORM = 0.
!
      eps = DLAMCH('Epsilon')
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
         CALL DCOPY(N,X(1,j),1,Work,1)
         CALL DTRMV(Uplo,Trans,Diag,N,A,Lda,Work,1)
         CALL DAXPY(N,-ONE,B(1,j),1,Work,1)
         bnorm = DASUM(N,Work,1)
         xnorm = DASUM(N,X(1,j),1)
         IF ( xnorm<=ZERO ) THEN
            Resid = ONE/eps
         ELSE
            Resid = MAX(Resid,((bnorm/anorm)/xnorm)/eps)
         ENDIF
      ENDDO
!
!
!     End of dtrt02
!
      END SUBROUTINE dtrt02
