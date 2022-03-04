!*==cptsvx.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief <b> CPTSVX computes the solution to system of linear equations A * X = B for PT matrices</b>
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CPTSVX + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cptsvx.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cptsvx.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cptsvx.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CPTSVX( FACT, N, NRHS, D, E, DF, EF, B, LDB, X, LDX,
!                          RCOND, FERR, BERR, WORK, RWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          FACT
!       INTEGER            INFO, LDB, LDX, N, NRHS
!       REAL               RCOND
!       ..
!       .. Array Arguments ..
!       REAL               BERR( * ), D( * ), DF( * ), FERR( * ),
!      $                   RWORK( * )
!       COMPLEX            B( LDB, * ), E( * ), EF( * ), WORK( * ),
!      $                   X( LDX, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CPTSVX uses the factorization A = L*D*L**H to compute the solution
!> to a complex system of linear equations A*X = B, where A is an
!> N-by-N Hermitian positive definite tridiagonal matrix and X and B
!> are N-by-NRHS matrices.
!>
!> Error bounds on the solution and a condition estimate are also
!> provided.
!> \endverbatim
!
!> \par Description:
!  =================
!>
!> \verbatim
!>
!> The following steps are performed:
!>
!> 1. If FACT = 'N', the matrix A is factored as A = L*D*L**H, where L
!>    is a unit lower bidiagonal matrix and D is diagonal.  The
!>    factorization can also be regarded as having the form
!>    A = U**H*D*U.
!>
!> 2. If the leading i-by-i principal minor is not positive definite,
!>    then the routine returns with INFO = i. Otherwise, the factored
!>    form of A is used to estimate the condition number of the matrix
!>    A.  If the reciprocal of the condition number is less than machine
!>    precision, INFO = N+1 is returned as a warning, but the routine
!>    still goes on to solve for X and compute error bounds as
!>    described below.
!>
!> 3. The system of equations is solved for X using the factored form
!>    of A.
!>
!> 4. Iterative refinement is applied to improve the computed solution
!>    matrix and calculate error bounds and backward error estimates
!>    for it.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] FACT
!> \verbatim
!>          FACT is CHARACTER*1
!>          Specifies whether or not the factored form of the matrix
!>          A is supplied on entry.
!>          = 'F':  On entry, DF and EF contain the factored form of A.
!>                  D, E, DF, and EF will not be modified.
!>          = 'N':  The matrix A will be copied to DF and EF and
!>                  factored.
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
!>          of the matrices B and X.  NRHS >= 0.
!> \endverbatim
!>
!> \param[in] D
!> \verbatim
!>          D is REAL array, dimension (N)
!>          The n diagonal elements of the tridiagonal matrix A.
!> \endverbatim
!>
!> \param[in] E
!> \verbatim
!>          E is COMPLEX array, dimension (N-1)
!>          The (n-1) subdiagonal elements of the tridiagonal matrix A.
!> \endverbatim
!>
!> \param[in,out] DF
!> \verbatim
!>          DF is REAL array, dimension (N)
!>          If FACT = 'F', then DF is an input argument and on entry
!>          contains the n diagonal elements of the diagonal matrix D
!>          from the L*D*L**H factorization of A.
!>          If FACT = 'N', then DF is an output argument and on exit
!>          contains the n diagonal elements of the diagonal matrix D
!>          from the L*D*L**H factorization of A.
!> \endverbatim
!>
!> \param[in,out] EF
!> \verbatim
!>          EF is COMPLEX array, dimension (N-1)
!>          If FACT = 'F', then EF is an input argument and on entry
!>          contains the (n-1) subdiagonal elements of the unit
!>          bidiagonal factor L from the L*D*L**H factorization of A.
!>          If FACT = 'N', then EF is an output argument and on exit
!>          contains the (n-1) subdiagonal elements of the unit
!>          bidiagonal factor L from the L*D*L**H factorization of A.
!> \endverbatim
!>
!> \param[in] B
!> \verbatim
!>          B is COMPLEX array, dimension (LDB,NRHS)
!>          The N-by-NRHS right hand side matrix B.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of the array B.  LDB >= max(1,N).
!> \endverbatim
!>
!> \param[out] X
!> \verbatim
!>          X is COMPLEX array, dimension (LDX,NRHS)
!>          If INFO = 0 or INFO = N+1, the N-by-NRHS solution matrix X.
!> \endverbatim
!>
!> \param[in] LDX
!> \verbatim
!>          LDX is INTEGER
!>          The leading dimension of the array X.  LDX >= max(1,N).
!> \endverbatim
!>
!> \param[out] RCOND
!> \verbatim
!>          RCOND is REAL
!>          The reciprocal condition number of the matrix A.  If RCOND
!>          is less than the machine precision (in particular, if
!>          RCOND = 0), the matrix is singular to working precision.
!>          This condition is indicated by a return code of INFO > 0.
!> \endverbatim
!>
!> \param[out] FERR
!> \verbatim
!>          FERR is REAL array, dimension (NRHS)
!>          The forward error bound for each solution vector
!>          X(j) (the j-th column of the solution matrix X).
!>          If XTRUE is the true solution corresponding to X(j), FERR(j)
!>          is an estimated upper bound for the magnitude of the largest
!>          element in (X(j) - XTRUE) divided by the magnitude of the
!>          largest element in X(j).
!> \endverbatim
!>
!> \param[out] BERR
!> \verbatim
!>          BERR is REAL array, dimension (NRHS)
!>          The componentwise relative backward error of each solution
!>          vector X(j) (i.e., the smallest relative change in any
!>          element of A or B that makes X(j) an exact solution).
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
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value
!>          > 0:  if INFO = i, and i is
!>                <= N:  the leading minor of order i of A is
!>                       not positive definite, so the factorization
!>                       could not be completed, and the solution has not
!>                       been computed. RCOND = 0 is returned.
!>                = N+1: U is nonsingular, but RCOND is less than machine
!>                       precision, meaning that the matrix is singular
!>                       to working precision.  Nevertheless, the
!>                       solution and error bounds are computed because
!>                       there are a number of situations where the
!>                       computed solution can be more accurate than the
!>                       value of RCOND would suggest.
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
!> \ingroup complexPTsolve
!
!  =====================================================================
      SUBROUTINE CPTSVX(Fact,N,Nrhs,D,E,Df,Ef,B,Ldb,X,Ldx,Rcond,Ferr,   &
     &                  Berr,Work,Rwork,Info)
      USE S_CCOPY
      USE S_CLACPY
      USE S_CLANHT
      USE S_CPTCON
      USE S_CPTRFS
      USE S_CPTTRF
      USE S_CPTTRS
      USE S_LSAME
      USE S_SCOPY
      USE S_SLAMCH
      USE S_XERBLA
      IMPLICIT NONE
!*--CPTSVX249
!
! PARAMETER definitions rewritten by SPAG
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Fact
      INTEGER :: N
      INTEGER :: Nrhs
      REAL , DIMENSION(*) :: D
      COMPLEX , DIMENSION(*) :: E
      REAL , DIMENSION(*) :: Df
      COMPLEX , DIMENSION(*) :: Ef
      COMPLEX , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX , DIMENSION(Ldx,*) :: X
      INTEGER :: Ldx
      REAL , INTENT(INOUT) :: Rcond
      REAL , DIMENSION(*) :: Ferr
      REAL , DIMENSION(*) :: Berr
      COMPLEX , DIMENSION(*) :: Work
      REAL , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      REAL :: anorm
      LOGICAL :: nofact
!
! End of declarations rewritten by SPAG
!
!     ..
!     .. Array Arguments ..
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
!     ..
!     .. Local Scalars ..
!     ..
!     .. External Functions ..
!     ..
!     .. External Subroutines ..
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      Info = 0
      nofact = LSAME(Fact,'N')
      IF ( .NOT.nofact .AND. .NOT.LSAME(Fact,'F') ) THEN
         Info = -1
      ELSEIF ( N<0 ) THEN
         Info = -2
      ELSEIF ( Nrhs<0 ) THEN
         Info = -3
      ELSEIF ( Ldb<MAX(1,N) ) THEN
         Info = -9
      ELSEIF ( Ldx<MAX(1,N) ) THEN
         Info = -11
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('CPTSVX',-Info)
         RETURN
      ENDIF
!
      IF ( nofact ) THEN
!
!        Compute the L*D*L**H (or U**H*D*U) factorization of A.
!
         CALL SCOPY(N,D,1,Df,1)
         IF ( N>1 ) CALL CCOPY(N-1,E,1,Ef,1)
         CALL CPTTRF(N,Df,Ef,Info)
!
!        Return if INFO is non-zero.
!
         IF ( Info>0 ) THEN
            Rcond = ZERO
            RETURN
         ENDIF
      ENDIF
!
!     Compute the norm of the matrix A.
!
      anorm = CLANHT('1',N,D,E)
!
!     Compute the reciprocal of the condition number of A.
!
      CALL CPTCON(N,Df,Ef,anorm,Rcond,Rwork,Info)
!
!     Compute the solution vectors X.
!
      CALL CLACPY('Full',N,Nrhs,B,Ldb,X,Ldx)
      CALL CPTTRS('Lower',N,Nrhs,Df,Ef,X,Ldx,Info)
!
!     Use iterative refinement to improve the computed solutions and
!     compute error bounds and backward error estimates for them.
!
      CALL CPTRFS('Lower',N,Nrhs,D,E,Df,Ef,B,Ldb,X,Ldx,Ferr,Berr,Work,  &
     &            Rwork,Info)
!
!     Set INFO = N+1 if the matrix is singular to working precision.
!
      IF ( Rcond<SLAMCH('Epsilon') ) Info = N + 1
!
!
!     End of CPTSVX
!
      END SUBROUTINE CPTSVX
