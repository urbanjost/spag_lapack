!*==sglmts.f90  processed by SPAG 7.51RB at 17:37 on  4 Mar 2022
!> \brief \b sglmts
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE sglmts( N, M, P, A, AF, LDA, B, BF, LDB, D, DF,
!                          X, U, WORK, LWORK, RWORK, RESULT )
!
!       .. Scalar Arguments ..
!       INTEGER            LDA, LDB, LWORK, M, P, N
!       REAL               RESULT
!       ..
!       .. Array Arguments ..
!       REAL               A( LDA, * ), AF( LDA, * ), B( LDB, * ),
!      $                   BF( LDB, * ), RWORK( * ), D( * ), DF( * ),
!      $                   U( * ), WORK( LWORK ), X( * )
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> sglmts tests SGGGLM - a subroutine for solving the generalized
!> linear model problem.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of rows of the matrices A and B.  N >= 0.
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of columns of the matrix A.  M >= 0.
!> \endverbatim
!>
!> \param[in] P
!> \verbatim
!>          P is INTEGER
!>          The number of columns of the matrix B.  P >= 0.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is REAL array, dimension (LDA,M)
!>          The N-by-M matrix A.
!> \endverbatim
!>
!> \param[out] AF
!> \verbatim
!>          AF is REAL array, dimension (LDA,M)
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the arrays A, AF. LDA >= max(M,N).
!> \endverbatim
!>
!> \param[in] B
!> \verbatim
!>          B is REAL array, dimension (LDB,P)
!>          The N-by-P matrix A.
!> \endverbatim
!>
!> \param[out] BF
!> \verbatim
!>          BF is REAL array, dimension (LDB,P)
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of the arrays B, BF. LDB >= max(P,N).
!> \endverbatim
!>
!> \param[in] D
!> \verbatim
!>          D is REAL array, dimension( N )
!>          On input, the left hand side of the GLM.
!> \endverbatim
!>
!> \param[out] DF
!> \verbatim
!>          DF is REAL array, dimension( N )
!> \endverbatim
!>
!> \param[out] X
!> \verbatim
!>          X is REAL array, dimension( M )
!>          solution vector X in the GLM problem.
!> \endverbatim
!>
!> \param[out] U
!> \verbatim
!>          U is REAL array, dimension( P )
!>          solution vector U in the GLM problem.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is REAL array, dimension (LWORK)
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The dimension of the array WORK.
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is REAL array, dimension (M)
!> \endverbatim
!>
!> \param[out] RESULT
!> \verbatim
!>          RESULT is REAL
!>          The test ratio:
!>                           norm( d - A*x - B*u )
!>            RESULT = -----------------------------------------
!>                     (norm(A)+norm(B))*(norm(x)+norm(u))*EPS
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
      SUBROUTINE sglmts(N,M,P,A,Af,Lda,B,Bf,Ldb,D,Df,X,U,Work,Lwork,    &
     &                  Rwork,Result)
      IMPLICIT NONE
!*--SGLMTS153
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Lda , Ldb , Lwork , M , P , N
      REAL Result
!     ..
!     .. Array Arguments ..
      REAL A(Lda,*) , Af(Lda,*) , B(Ldb,*) , Bf(Ldb,*) , Rwork(*) ,     &
     &     D(*) , Df(*) , U(*) , Work(Lwork) , X(*)
!
!  ====================================================================
!
!     .. Parameters ..
      REAL ZERO , ONE
      PARAMETER (ZERO=0.0E+0,ONE=1.0E+0)
!     ..
!     .. Local Scalars ..
      INTEGER info
      REAL anorm , bnorm , eps , xnorm , ynorm , dnorm , unfl
!     ..
!     .. External Functions ..
      REAL SASUM , SLAMCH , SLANGE
      EXTERNAL SASUM , SLAMCH , SLANGE
!     ..
!     .. External Subroutines ..
      EXTERNAL SLACPY
!
!     .. Intrinsic Functions ..
      INTRINSIC MAX
!     ..
!     .. Executable Statements ..
!
      eps = SLAMCH('Epsilon')
      unfl = SLAMCH('Safe minimum')
      anorm = MAX(SLANGE('1',N,M,A,Lda,Rwork),unfl)
      bnorm = MAX(SLANGE('1',N,P,B,Ldb,Rwork),unfl)
!
!     Copy the matrices A and B to the arrays AF and BF,
!     and the vector D the array DF.
!
      CALL SLACPY('Full',N,M,A,Lda,Af,Lda)
      CALL SLACPY('Full',N,P,B,Ldb,Bf,Ldb)
      CALL SCOPY(N,D,1,Df,1)
!
!     Solve GLM problem
!
      CALL SGGGLM(N,M,P,Af,Lda,Bf,Ldb,Df,X,U,Work,Lwork,info)
!
!     Test the residual for the solution of LSE
!
!                       norm( d - A*x - B*u )
!       RESULT = -----------------------------------------
!                (norm(A)+norm(B))*(norm(x)+norm(u))*EPS
!
      CALL SCOPY(N,D,1,Df,1)
      CALL SGEMV('No transpose',N,M,-ONE,A,Lda,X,1,ONE,Df,1)
!
      CALL SGEMV('No transpose',N,P,-ONE,B,Ldb,U,1,ONE,Df,1)
!
      dnorm = SASUM(N,Df,1)
      xnorm = SASUM(M,X,1) + SASUM(P,U,1)
      ynorm = anorm + bnorm
!
      IF ( xnorm<=ZERO ) THEN
         Result = ZERO
      ELSE
         Result = ((dnorm/ynorm)/xnorm)/eps
      ENDIF
!
!
!     End of sglmts
!
      END SUBROUTINE sglmts
