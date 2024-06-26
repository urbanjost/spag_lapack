!*==dlsets.f90  processed by SPAG 7.51RB at 17:37 on  4 Mar 2022
!> \brief \b dlsets
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE dlsets( M, P, N, A, AF, LDA, B, BF, LDB, C, CF, D, DF,
!                          X, WORK, LWORK, RWORK, RESULT )
!
!       .. Scalar Arguments ..
!       INTEGER            LDA, LDB, LWORK, M, N, P
!       ..
!       .. Array Arguments ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> dlsets tests DGGLSE - a subroutine for solving linear equality
!> constrained least square problem (LSE).
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix A.  M >= 0.
!> \endverbatim
!>
!> \param[in] P
!> \verbatim
!>          P is INTEGER
!>          The number of rows of the matrix B.  P >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrices A and B.  N >= 0.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (LDA,N)
!>          The M-by-N matrix A.
!> \endverbatim
!>
!> \param[out] AF
!> \verbatim
!>          AF is DOUBLE PRECISION array, dimension (LDA,N)
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the arrays A, AF, Q and R.
!>          LDA >= max(M,N).
!> \endverbatim
!>
!> \param[in] B
!> \verbatim
!>          B is DOUBLE PRECISION array, dimension (LDB,N)
!>          The P-by-N matrix A.
!> \endverbatim
!>
!> \param[out] BF
!> \verbatim
!>          BF is DOUBLE PRECISION array, dimension (LDB,N)
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of the arrays B, BF, V and S.
!>          LDB >= max(P,N).
!> \endverbatim
!>
!> \param[in] C
!> \verbatim
!>          C is DOUBLE PRECISION array, dimension( M )
!>          the vector C in the LSE problem.
!> \endverbatim
!>
!> \param[out] CF
!> \verbatim
!>          CF is DOUBLE PRECISION array, dimension( M )
!> \endverbatim
!>
!> \param[in] D
!> \verbatim
!>          D is DOUBLE PRECISION array, dimension( P )
!>          the vector D in the LSE problem.
!> \endverbatim
!>
!> \param[out] DF
!> \verbatim
!>          DF is DOUBLE PRECISION array, dimension( P )
!> \endverbatim
!>
!> \param[out] X
!> \verbatim
!>          X is DOUBLE PRECISION array, dimension( N )
!>          solution vector X in the LSE problem.
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
!>          RESULT is DOUBLE PRECISION array, dimension (2)
!>          The test ratios:
!>            RESULT(1) = norm( A*x - c )/ norm(A)*norm(X)*EPS
!>            RESULT(2) = norm( B*x - d )/ norm(B)*norm(X)*EPS
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
      SUBROUTINE dlsets(M,P,N,A,Af,Lda,B,Bf,Ldb,C,Cf,D,Df,X,Work,Lwork, &
     &                  Rwork,Result)
      IMPLICIT NONE
!*--DLSETS155
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Lda , Ldb , Lwork , M , N , P
!     ..
!     .. Array Arguments ..
!
!  ====================================================================
!
      DOUBLE PRECISION A(Lda,*) , Af(Lda,*) , B(Ldb,*) , Bf(Ldb,*) ,    &
     &                 C(*) , Cf(*) , D(*) , Df(*) , Result(2) ,        &
     &                 Rwork(*) , Work(Lwork) , X(*)
!     ..
!     .. Local Scalars ..
      INTEGER info
!     ..
!     .. External Subroutines ..
      EXTERNAL DCOPY , DGGLSE , DLACPY
!     ..
!     .. Executable Statements ..
!
!     Copy the matrices A and B to the arrays AF and BF,
!     and the vectors C and D to the arrays CF and DF,
!
      CALL DLACPY('Full',M,N,A,Lda,Af,Lda)
      CALL DLACPY('Full',P,N,B,Ldb,Bf,Ldb)
      CALL DCOPY(M,C,1,Cf,1)
      CALL DCOPY(P,D,1,Df,1)
!
!     Solve LSE problem
!
      CALL DGGLSE(M,N,P,Af,Lda,Bf,Ldb,Cf,Df,X,Work,Lwork,info)
!
!     Test the residual for the solution of LSE
!
!     Compute RESULT(1) = norm( A*x - c ) / norm(A)*norm(X)*EPS
!
      CALL DCOPY(M,C,1,Cf,1)
      CALL DCOPY(P,D,1,Df,1)
      CALL dget02('No transpose',M,N,1,A,Lda,X,N,Cf,M,Rwork,Result(1))
!
!     Compute result(2) = norm( B*x - d ) / norm(B)*norm(X)*EPS
!
      CALL dget02('No transpose',P,N,1,B,Ldb,X,N,Df,P,Rwork,Result(2))
!
!
!     End of dlsets
!
      END SUBROUTINE dlsets
