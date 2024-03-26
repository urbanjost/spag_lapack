!*==spbt02.f90  processed by SPAG 7.51RB at 17:34 on  4 Mar 2022
!> \brief \b spbt02
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE spbt02( UPLO, N, KD, NRHS, A, LDA, X, LDX, B, LDB,
!                          RWORK, RESID )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            KD, LDA, LDB, LDX, N, NRHS
!       REAL               RESID
!       ..
!       .. Array Arguments ..
!       REAL               A( LDA, * ), B( LDB, * ), RWORK( * ),
!      $                   X( LDX, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> spbt02 computes the residual for a solution of a symmetric banded
!> system of equations  A*x = b:
!>    RESID = norm( B - A*X ) / ( norm(A) * norm(X) * EPS)
!> where EPS is the machine precision.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          Specifies whether the upper or lower triangular part of the
!>          symmetric matrix A is stored:
!>          = 'U':  Upper triangular
!>          = 'L':  Lower triangular
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of rows and columns of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in] KD
!> \verbatim
!>          KD is INTEGER
!>          The number of super-diagonals of the matrix A if UPLO = 'U',
!>          or the number of sub-diagonals if UPLO = 'L'.  KD >= 0.
!> \endverbatim
!>
!> \param[in] NRHS
!> \verbatim
!>          NRHS is INTEGER
!>          The number of right hand sides. NRHS >= 0.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is REAL array, dimension (LDA,N)
!>          The original symmetric band matrix A.  If UPLO = 'U', the
!>          upper triangular part of A is stored as a band matrix; if
!>          UPLO = 'L', the lower triangular part of A is stored.  The
!>          columns of the appropriate triangle are stored in the columns
!>          of A and the diagonals of the triangle are stored in the rows
!>          of A.  See SPBTRF for further details.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER.
!>          The leading dimension of the array A.  LDA >= max(1,KD+1).
!> \endverbatim
!>
!> \param[in] X
!> \verbatim
!>          X is REAL array, dimension (LDX,NRHS)
!>          The computed solution vectors for the system of linear
!>          equations.
!> \endverbatim
!>
!> \param[in] LDX
!> \verbatim
!>          LDX is INTEGER
!>          The leading dimension of the array X.   LDX >= max(1,N).
!> \endverbatim
!>
!> \param[in,out] B
!> \verbatim
!>          B is REAL array, dimension (LDB,NRHS)
!>          On entry, the right hand side vectors for the system of
!>          linear equations.
!>          On exit, B is overwritten with the difference B - A*X.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of the array B.  LDB >= max(1,N).
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
!>          norm(B - A*X) / ( norm(A) * norm(X) * EPS ).
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
      SUBROUTINE spbt02(Uplo,N,Kd,Nrhs,A,Lda,X,Ldx,B,Ldb,Rwork,Resid)
      IMPLICIT NONE
!*--SPBT02139
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Uplo
      INTEGER Kd , Lda , Ldb , Ldx , N , Nrhs
      REAL Resid
!     ..
!     .. Array Arguments ..
      REAL A(Lda,*) , B(Ldb,*) , Rwork(*) , X(Ldx,*)
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
      REAL SASUM , SLAMCH , SLANSB
EXTERNAL    SASUM,   SLAMCH,   SLANSB
!     ..
!     .. External Subroutines ..
EXTERNAL    SSBMV
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX
!     ..
!     .. Executable Statements ..
!
!     Quick exit if N = 0 or NRHS = 0.
!
      IF ( N<=0 .OR. Nrhs<=0 ) THEN
         Resid = ZERO
         RETURN
      ENDIF
!
!     Exit with RESID = 1/EPS if ANORM = 0.
!
      eps = SLAMCH('Epsilon')
      anorm = SLANSB('1',Uplo,N,Kd,A,Lda,Rwork)
      IF ( anorm<=ZERO ) THEN
         Resid = ONE/eps
         RETURN
      ENDIF
!
!     Compute  B - A*X
!
      DO j = 1 , Nrhs
         CALL SSBMV(Uplo,N,Kd,-ONE,A,Lda,X(1,j),1,ONE,B(1,j),1)
      ENDDO
!
!     Compute the maximum over the number of right hand sides of
!          norm( B - A*X ) / ( norm(A) * norm(X) * EPS )
!
      Resid = ZERO
      DO j = 1 , Nrhs
         bnorm = SASUM(N,B(1,j),1)
         xnorm = SASUM(N,X(1,j),1)
         IF ( xnorm<=ZERO ) THEN
            Resid = ONE/eps
         ELSE
            Resid = MAX(Resid,((bnorm/anorm)/xnorm)/eps)
         ENDIF
      ENDDO
!
!
!     End of spbt02
!
      END SUBROUTINE spbt02
