!*==cgtsvx.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief <b> CGTSVX computes the solution to system of linear equations A * X = B for GT matrices </b>
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CGTSVX + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgtsvx.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgtsvx.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgtsvx.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CGTSVX( FACT, TRANS, N, NRHS, DL, D, DU, DLF, DF, DUF,
!                          DU2, IPIV, B, LDB, X, LDX, RCOND, FERR, BERR,
!                          WORK, RWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          FACT, TRANS
!       INTEGER            INFO, LDB, LDX, N, NRHS
!       REAL               RCOND
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * )
!       REAL               BERR( * ), FERR( * ), RWORK( * )
!       COMPLEX            B( LDB, * ), D( * ), DF( * ), DL( * ),
!      $                   DLF( * ), DU( * ), DU2( * ), DUF( * ),
!      $                   WORK( * ), X( LDX, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CGTSVX uses the LU factorization to compute the solution to a complex
!> system of linear equations A * X = B, A**T * X = B, or A**H * X = B,
!> where A is a tridiagonal matrix of order N and X and B are N-by-NRHS
!> matrices.
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
!> 1. If FACT = 'N', the LU decomposition is used to factor the matrix A
!>    as A = L * U, where L is a product of permutation and unit lower
!>    bidiagonal matrices and U is upper triangular with nonzeros in
!>    only the main diagonal and first two superdiagonals.
!>
!> 2. If some U(i,i)=0, so that U is exactly singular, then the routine
!>    returns with INFO = i. Otherwise, the factored form of A is used
!>    to estimate the condition number of the matrix A.  If the
!>    reciprocal of the condition number is less than machine precision,
!>    INFO = N+1 is returned as a warning, but the routine still goes on
!>    to solve for X and compute error bounds as described below.
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
!>          Specifies whether or not the factored form of A has been
!>          supplied on entry.
!>          = 'F':  DLF, DF, DUF, DU2, and IPIV contain the factored form
!>                  of A; DL, D, DU, DLF, DF, DUF, DU2 and IPIV will not
!>                  be modified.
!>          = 'N':  The matrix will be copied to DLF, DF, and DUF
!>                  and factored.
!> \endverbatim
!>
!> \param[in] TRANS
!> \verbatim
!>          TRANS is CHARACTER*1
!>          Specifies the form of the system of equations:
!>          = 'N':  A * X = B     (No transpose)
!>          = 'T':  A**T * X = B  (Transpose)
!>          = 'C':  A**H * X = B  (Conjugate transpose)
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
!>          of the matrix B.  NRHS >= 0.
!> \endverbatim
!>
!> \param[in] DL
!> \verbatim
!>          DL is COMPLEX array, dimension (N-1)
!>          The (n-1) subdiagonal elements of A.
!> \endverbatim
!>
!> \param[in] D
!> \verbatim
!>          D is COMPLEX array, dimension (N)
!>          The n diagonal elements of A.
!> \endverbatim
!>
!> \param[in] DU
!> \verbatim
!>          DU is COMPLEX array, dimension (N-1)
!>          The (n-1) superdiagonal elements of A.
!> \endverbatim
!>
!> \param[in,out] DLF
!> \verbatim
!>          DLF is COMPLEX array, dimension (N-1)
!>          If FACT = 'F', then DLF is an input argument and on entry
!>          contains the (n-1) multipliers that define the matrix L from
!>          the LU factorization of A as computed by CGTTRF.
!>
!>          If FACT = 'N', then DLF is an output argument and on exit
!>          contains the (n-1) multipliers that define the matrix L from
!>          the LU factorization of A.
!> \endverbatim
!>
!> \param[in,out] DF
!> \verbatim
!>          DF is COMPLEX array, dimension (N)
!>          If FACT = 'F', then DF is an input argument and on entry
!>          contains the n diagonal elements of the upper triangular
!>          matrix U from the LU factorization of A.
!>
!>          If FACT = 'N', then DF is an output argument and on exit
!>          contains the n diagonal elements of the upper triangular
!>          matrix U from the LU factorization of A.
!> \endverbatim
!>
!> \param[in,out] DUF
!> \verbatim
!>          DUF is COMPLEX array, dimension (N-1)
!>          If FACT = 'F', then DUF is an input argument and on entry
!>          contains the (n-1) elements of the first superdiagonal of U.
!>
!>          If FACT = 'N', then DUF is an output argument and on exit
!>          contains the (n-1) elements of the first superdiagonal of U.
!> \endverbatim
!>
!> \param[in,out] DU2
!> \verbatim
!>          DU2 is COMPLEX array, dimension (N-2)
!>          If FACT = 'F', then DU2 is an input argument and on entry
!>          contains the (n-2) elements of the second superdiagonal of
!>          U.
!>
!>          If FACT = 'N', then DU2 is an output argument and on exit
!>          contains the (n-2) elements of the second superdiagonal of
!>          U.
!> \endverbatim
!>
!> \param[in,out] IPIV
!> \verbatim
!>          IPIV is INTEGER array, dimension (N)
!>          If FACT = 'F', then IPIV is an input argument and on entry
!>          contains the pivot indices from the LU factorization of A as
!>          computed by CGTTRF.
!>
!>          If FACT = 'N', then IPIV is an output argument and on exit
!>          contains the pivot indices from the LU factorization of A;
!>          row i of the matrix was interchanged with row IPIV(i).
!>          IPIV(i) will always be either i or i+1; IPIV(i) = i indicates
!>          a row interchange was not required.
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
!>          The estimate of the reciprocal condition number of the matrix
!>          A.  If RCOND is less than the machine precision (in
!>          particular, if RCOND = 0), the matrix is singular to working
!>          precision.  This condition is indicated by a return code of
!>          INFO > 0.
!> \endverbatim
!>
!> \param[out] FERR
!> \verbatim
!>          FERR is REAL array, dimension (NRHS)
!>          The estimated forward error bound for each solution vector
!>          X(j) (the j-th column of the solution matrix X).
!>          If XTRUE is the true solution corresponding to X(j), FERR(j)
!>          is an estimated upper bound for the magnitude of the largest
!>          element in (X(j) - XTRUE) divided by the magnitude of the
!>          largest element in X(j).  The estimate is as reliable as
!>          the estimate for RCOND, and is almost always a slight
!>          overestimate of the true error.
!> \endverbatim
!>
!> \param[out] BERR
!> \verbatim
!>          BERR is REAL array, dimension (NRHS)
!>          The componentwise relative backward error of each solution
!>          vector X(j) (i.e., the smallest relative change in
!>          any element of A or B that makes X(j) an exact solution).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX array, dimension (2*N)
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
!>                <= N:  U(i,i) is exactly zero.  The factorization
!>                       has not been completed unless i = N, but the
!>                       factor U is exactly singular, so the solution
!>                       and error bounds could not be computed.
!>                       RCOND = 0 is returned.
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
!> \ingroup complexGTsolve
!
!  =====================================================================
      SUBROUTINE CGTSVX(Fact,Trans,N,Nrhs,Dl,D,Du,Dlf,Df,Duf,Du2,Ipiv,B,&
     &                  Ldb,X,Ldx,Rcond,Ferr,Berr,Work,Rwork,Info)
      USE S_CCOPY
      USE S_CGTCON
      USE S_CGTRFS
      USE S_CGTTRF
      USE S_CGTTRS
      USE S_CLACPY
      USE S_CLANGT
      USE S_LSAME
      USE S_SLAMCH
      USE S_XERBLA
      IMPLICIT NONE
!*--CGTSVX307
!
! PARAMETER definitions rewritten by SPAG
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Fact
      CHARACTER :: Trans
      INTEGER :: N
      INTEGER :: Nrhs
      COMPLEX , DIMENSION(*) :: Dl
      COMPLEX , DIMENSION(*) :: D
      COMPLEX , DIMENSION(*) :: Du
      COMPLEX , DIMENSION(*) :: Dlf
      COMPLEX , DIMENSION(*) :: Df
      COMPLEX , DIMENSION(*) :: Duf
      COMPLEX , DIMENSION(*) :: Du2
      INTEGER , DIMENSION(*) :: Ipiv
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
      LOGICAL :: nofact , notran
      CHARACTER :: norm
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
      Info = 0
      nofact = LSAME(Fact,'N')
      notran = LSAME(Trans,'N')
      IF ( .NOT.nofact .AND. .NOT.LSAME(Fact,'F') ) THEN
         Info = -1
      ELSEIF ( .NOT.notran .AND. .NOT.LSAME(Trans,'T') .AND.            &
     &         .NOT.LSAME(Trans,'C') ) THEN
         Info = -2
      ELSEIF ( N<0 ) THEN
         Info = -3
      ELSEIF ( Nrhs<0 ) THEN
         Info = -4
      ELSEIF ( Ldb<MAX(1,N) ) THEN
         Info = -14
      ELSEIF ( Ldx<MAX(1,N) ) THEN
         Info = -16
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('CGTSVX',-Info)
         RETURN
      ENDIF
!
      IF ( nofact ) THEN
!
!        Compute the LU factorization of A.
!
         CALL CCOPY(N,D,1,Df,1)
         IF ( N>1 ) THEN
            CALL CCOPY(N-1,Dl,1,Dlf,1)
            CALL CCOPY(N-1,Du,1,Duf,1)
         ENDIF
         CALL CGTTRF(N,Dlf,Df,Duf,Du2,Ipiv,Info)
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
      IF ( notran ) THEN
         norm = '1'
      ELSE
         norm = 'I'
      ENDIF
      anorm = CLANGT(norm,N,Dl,D,Du)
!
!     Compute the reciprocal of the condition number of A.
!
      CALL CGTCON(norm,N,Dlf,Df,Duf,Du2,Ipiv,anorm,Rcond,Work,Info)
!
!     Compute the solution vectors X.
!
      CALL CLACPY('Full',N,Nrhs,B,Ldb,X,Ldx)
      CALL CGTTRS(Trans,N,Nrhs,Dlf,Df,Duf,Du2,Ipiv,X,Ldx,Info)
!
!     Use iterative refinement to improve the computed solutions and
!     compute error bounds and backward error estimates for them.
!
      CALL CGTRFS(Trans,N,Nrhs,Dl,D,Du,Dlf,Df,Duf,Du2,Ipiv,B,Ldb,X,Ldx, &
     &            Ferr,Berr,Work,Rwork,Info)
!
!     Set INFO = N+1 if the matrix is singular to working precision.
!
      IF ( Rcond<SLAMCH('Epsilon') ) Info = N + 1
!
!
!     End of CGTSVX
!
      END SUBROUTINE CGTSVX
