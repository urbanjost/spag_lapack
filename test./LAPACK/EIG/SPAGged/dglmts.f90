!*==dglmts.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b DGLMTS
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DGLMTS( N, M, P, A, AF, LDA, B, BF, LDB, D, DF, X, U,
!                          WORK, LWORK, RWORK, RESULT )
!
!       .. Scalar Arguments ..
!       INTEGER            LDA, LDB, LWORK, M, N, P
!       DOUBLE PRECISION   RESULT
!       ..
!       .. Array Arguments ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DGLMTS tests DGGGLM - a subroutine for solving the generalized
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
!>          A is DOUBLE PRECISION array, dimension (LDA,M)
!>          The N-by-M matrix A.
!> \endverbatim
!>
!> \param[out] AF
!> \verbatim
!>          AF is DOUBLE PRECISION array, dimension (LDA,M)
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
!>          B is DOUBLE PRECISION array, dimension (LDB,P)
!>          The N-by-P matrix A.
!> \endverbatim
!>
!> \param[out] BF
!> \verbatim
!>          BF is DOUBLE PRECISION array, dimension (LDB,P)
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
!>          D is DOUBLE PRECISION array, dimension( N )
!>          On input, the left hand side of the GLM.
!> \endverbatim
!>
!> \param[out] DF
!> \verbatim
!>          DF is DOUBLE PRECISION array, dimension( N )
!> \endverbatim
!>
!> \param[out] X
!> \verbatim
!>          X is DOUBLE PRECISION array, dimension( M )
!>          solution vector X in the GLM problem.
!> \endverbatim
!>
!> \param[out] U
!> \verbatim
!>          U is DOUBLE PRECISION array, dimension( P )
!>          solution vector U in the GLM problem.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (LWORK)
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
!>          RWORK is DOUBLE PRECISION array, dimension (M)
!> \endverbatim
!>
!> \param[out] RESULT
!> \verbatim
!>          RESULT is DOUBLE PRECISION
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
!> \ingroup double_eig
!
!  =====================================================================
      SUBROUTINE DGLMTS(N,M,P,A,Af,Lda,B,Bf,Ldb,D,Df,X,U,Work,Lwork,    &
     &                  Rwork,Result)
      IMPLICIT NONE
!*--DGLMTS150
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Lda , Ldb , Lwork , M , N , P
      DOUBLE PRECISION Result
!     ..
!     .. Array Arguments ..
!
!  ====================================================================
!
      DOUBLE PRECISION A(Lda,*) , Af(Lda,*) , B(Ldb,*) , Bf(Ldb,*) ,    &
     &                 D(*) , Df(*) , Rwork(*) , U(*) , Work(Lwork) ,   &
     &                 X(*)
!     ..
!     .. Parameters ..
      DOUBLE PRECISION ZERO , ONE
      PARAMETER (ZERO=0.0D+0,ONE=1.0D+0)
!     ..
!     .. Local Scalars ..
      INTEGER info
      DOUBLE PRECISION anorm , bnorm , dnorm , eps , unfl , xnorm ,     &
     &                 ynorm
!     ..
!     .. External Functions ..
      DOUBLE PRECISION DASUM , DLAMCH , DLANGE
      EXTERNAL DASUM , DLAMCH , DLANGE
!     ..
!     .. External Subroutines ..
!
      EXTERNAL DCOPY , DGEMV , DGGGLM , DLACPY
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX
!     ..
!     .. Executable Statements ..
!
      eps = DLAMCH('Epsilon')
      unfl = DLAMCH('Safe minimum')
      anorm = MAX(DLANGE('1',N,M,A,Lda,Rwork),unfl)
      bnorm = MAX(DLANGE('1',N,P,B,Ldb,Rwork),unfl)
!
!     Copy the matrices A and B to the arrays AF and BF,
!     and the vector D the array DF.
!
      CALL DLACPY('Full',N,M,A,Lda,Af,Lda)
      CALL DLACPY('Full',N,P,B,Ldb,Bf,Ldb)
      CALL DCOPY(N,D,1,Df,1)
!
!     Solve GLM problem
!
      CALL DGGGLM(N,M,P,Af,Lda,Bf,Ldb,Df,X,U,Work,Lwork,info)
!
!     Test the residual for the solution of LSE
!
!                       norm( d - A*x - B*u )
!       RESULT = -----------------------------------------
!                (norm(A)+norm(B))*(norm(x)+norm(u))*EPS
!
      CALL DCOPY(N,D,1,Df,1)
      CALL DGEMV('No transpose',N,M,-ONE,A,Lda,X,1,ONE,Df,1)
!
      CALL DGEMV('No transpose',N,P,-ONE,B,Ldb,U,1,ONE,Df,1)
!
      dnorm = DASUM(N,Df,1)
      xnorm = DASUM(M,X,1) + DASUM(P,U,1)
      ynorm = anorm + bnorm
!
      IF ( xnorm<=ZERO ) THEN
         Result = ZERO
      ELSE
         Result = ((dnorm/ynorm)/xnorm)/eps
      ENDIF
!
!
!     End of DGLMTS
!
      END SUBROUTINE DGLMTS
