!*==zsysvxx.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief <b> ZSYSVXX computes the solution to system of linear equations A * X = B for SY matrices</b>
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZSYSVXX + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zsysvxx.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zsysvxx.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zsysvxx.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZSYSVXX( FACT, UPLO, N, NRHS, A, LDA, AF, LDAF, IPIV,
!                           EQUED, S, B, LDB, X, LDX, RCOND, RPVGRW, BERR,
!                           N_ERR_BNDS, ERR_BNDS_NORM, ERR_BNDS_COMP,
!                           NPARAMS, PARAMS, WORK, RWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          EQUED, FACT, UPLO
!       INTEGER            INFO, LDA, LDAF, LDB, LDX, N, NRHS, NPARAMS,
!      $                   N_ERR_BNDS
!       DOUBLE PRECISION   RCOND, RPVGRW
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * )
!       COMPLEX*16         A( LDA, * ), AF( LDAF, * ), B( LDB, * ),
!      $                   X( LDX, * ), WORK( * )
!       DOUBLE PRECISION   S( * ), PARAMS( * ), BERR( * ),
!      $                   ERR_BNDS_NORM( NRHS, * ),
!      $                   ERR_BNDS_COMP( NRHS, * ), RWORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    ZSYSVXX uses the diagonal pivoting factorization to compute the
!>    solution to a complex*16 system of linear equations A * X = B, where
!>    A is an N-by-N symmetric matrix and X and B are N-by-NRHS
!>    matrices.
!>
!>    If requested, both normwise and maximum componentwise error bounds
!>    are returned. ZSYSVXX will return a solution with a tiny
!>    guaranteed error (O(eps) where eps is the working machine
!>    precision) unless the matrix is very ill-conditioned, in which
!>    case a warning is returned. Relevant condition numbers also are
!>    calculated and returned.
!>
!>    ZSYSVXX accepts user-provided factorizations and equilibration
!>    factors; see the definitions of the FACT and EQUED options.
!>    Solving with refinement and using a factorization from a previous
!>    ZSYSVXX call will also produce a solution with either O(eps)
!>    errors or warnings, but we cannot make that claim for general
!>    user-provided factorizations and equilibration factors if they
!>    differ from what ZSYSVXX would itself produce.
!> \endverbatim
!
!> \par Description:
!  =================
!>
!> \verbatim
!>
!>    The following steps are performed:
!>
!>    1. If FACT = 'E', double precision scaling factors are computed to equilibrate
!>    the system:
!>
!>      diag(S)*A*diag(S)     *inv(diag(S))*X = diag(S)*B
!>
!>    Whether or not the system will be equilibrated depends on the
!>    scaling of the matrix A, but if equilibration is used, A is
!>    overwritten by diag(S)*A*diag(S) and B by diag(S)*B.
!>
!>    2. If FACT = 'N' or 'E', the LU decomposition is used to factor
!>    the matrix A (after equilibration if FACT = 'E') as
!>
!>       A = U * D * U**T,  if UPLO = 'U', or
!>       A = L * D * L**T,  if UPLO = 'L',
!>
!>    where U (or L) is a product of permutation and unit upper (lower)
!>    triangular matrices, and D is symmetric and block diagonal with
!>    1-by-1 and 2-by-2 diagonal blocks.
!>
!>    3. If some D(i,i)=0, so that D is exactly singular, then the
!>    routine returns with INFO = i. Otherwise, the factored form of A
!>    is used to estimate the condition number of the matrix A (see
!>    argument RCOND).  If the reciprocal of the condition number is
!>    less than machine precision, the routine still goes on to solve
!>    for X and compute error bounds as described below.
!>
!>    4. The system of equations is solved for X using the factored form
!>    of A.
!>
!>    5. By default (unless PARAMS(LA_LINRX_ITREF_I) is set to zero),
!>    the routine will use iterative refinement to try to get a small
!>    error and error bounds.  Refinement calculates the residual to at
!>    least twice the working precision.
!>
!>    6. If equilibration was used, the matrix X is premultiplied by
!>    diag(R) so that it solves the original system before
!>    equilibration.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \verbatim
!>     Some optional parameters are bundled in the PARAMS array.  These
!>     settings determine how refinement is performed, but often the
!>     defaults are acceptable.  If the defaults are acceptable, users
!>     can pass NPARAMS = 0 which prevents the source code from accessing
!>     the PARAMS argument.
!> \endverbatim
!>
!> \param[in] FACT
!> \verbatim
!>          FACT is CHARACTER*1
!>     Specifies whether or not the factored form of the matrix A is
!>     supplied on entry, and if not, whether the matrix A should be
!>     equilibrated before it is factored.
!>       = 'F':  On entry, AF and IPIV contain the factored form of A.
!>               If EQUED is not 'N', the matrix A has been
!>               equilibrated with scaling factors given by S.
!>               A, AF, and IPIV are not modified.
!>       = 'N':  The matrix A will be copied to AF and factored.
!>       = 'E':  The matrix A will be equilibrated if necessary, then
!>               copied to AF and factored.
!> \endverbatim
!>
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>       = 'U':  Upper triangle of A is stored;
!>       = 'L':  Lower triangle of A is stored.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>     The number of linear equations, i.e., the order of the
!>     matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in] NRHS
!> \verbatim
!>          NRHS is INTEGER
!>     The number of right hand sides, i.e., the number of columns
!>     of the matrices B and X.  NRHS >= 0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension (LDA,N)
!>     The symmetric matrix A.  If UPLO = 'U', the leading N-by-N
!>     upper triangular part of A contains the upper triangular
!>     part of the matrix A, and the strictly lower triangular
!>     part of A is not referenced.  If UPLO = 'L', the leading
!>     N-by-N lower triangular part of A contains the lower
!>     triangular part of the matrix A, and the strictly upper
!>     triangular part of A is not referenced.
!>
!>     On exit, if FACT = 'E' and EQUED = 'Y', A is overwritten by
!>     diag(S)*A*diag(S).
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>     The leading dimension of the array A.  LDA >= max(1,N).
!> \endverbatim
!>
!> \param[in,out] AF
!> \verbatim
!>          AF is COMPLEX*16 array, dimension (LDAF,N)
!>     If FACT = 'F', then AF is an input argument and on entry
!>     contains the block diagonal matrix D and the multipliers
!>     used to obtain the factor U or L from the factorization A =
!>     U*D*U**T or A = L*D*L**T as computed by DSYTRF.
!>
!>     If FACT = 'N', then AF is an output argument and on exit
!>     returns the block diagonal matrix D and the multipliers
!>     used to obtain the factor U or L from the factorization A =
!>     U*D*U**T or A = L*D*L**T.
!> \endverbatim
!>
!> \param[in] LDAF
!> \verbatim
!>          LDAF is INTEGER
!>     The leading dimension of the array AF.  LDAF >= max(1,N).
!> \endverbatim
!>
!> \param[in,out] IPIV
!> \verbatim
!>          IPIV is INTEGER array, dimension (N)
!>     If FACT = 'F', then IPIV is an input argument and on entry
!>     contains details of the interchanges and the block
!>     structure of D, as determined by DSYTRF.  If IPIV(k) > 0,
!>     then rows and columns k and IPIV(k) were interchanged and
!>     D(k,k) is a 1-by-1 diagonal block.  If UPLO = 'U' and
!>     IPIV(k) = IPIV(k-1) < 0, then rows and columns k-1 and
!>     -IPIV(k) were interchanged and D(k-1:k,k-1:k) is a 2-by-2
!>     diagonal block.  If UPLO = 'L' and IPIV(k) = IPIV(k+1) < 0,
!>     then rows and columns k+1 and -IPIV(k) were interchanged
!>     and D(k:k+1,k:k+1) is a 2-by-2 diagonal block.
!>
!>     If FACT = 'N', then IPIV is an output argument and on exit
!>     contains details of the interchanges and the block
!>     structure of D, as determined by DSYTRF.
!> \endverbatim
!>
!> \param[in,out] EQUED
!> \verbatim
!>          EQUED is CHARACTER*1
!>     Specifies the form of equilibration that was done.
!>       = 'N':  No equilibration (always true if FACT = 'N').
!>       = 'Y':  Both row and column equilibration, i.e., A has been
!>               replaced by diag(S) * A * diag(S).
!>     EQUED is an input argument if FACT = 'F'; otherwise, it is an
!>     output argument.
!> \endverbatim
!>
!> \param[in,out] S
!> \verbatim
!>          S is DOUBLE PRECISION array, dimension (N)
!>     The scale factors for A.  If EQUED = 'Y', A is multiplied on
!>     the left and right by diag(S).  S is an input argument if FACT =
!>     'F'; otherwise, S is an output argument.  If FACT = 'F' and EQUED
!>     = 'Y', each element of S must be positive.  If S is output, each
!>     element of S is a power of the radix. If S is input, each element
!>     of S should be a power of the radix to ensure a reliable solution
!>     and error estimates. Scaling by powers of the radix does not cause
!>     rounding errors unless the result underflows or overflows.
!>     Rounding errors during scaling lead to refining with a matrix that
!>     is not equivalent to the input matrix, producing error estimates
!>     that may not be reliable.
!> \endverbatim
!>
!> \param[in,out] B
!> \verbatim
!>          B is COMPLEX*16 array, dimension (LDB,NRHS)
!>     On entry, the N-by-NRHS right hand side matrix B.
!>     On exit,
!>     if EQUED = 'N', B is not modified;
!>     if EQUED = 'Y', B is overwritten by diag(S)*B;
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>     The leading dimension of the array B.  LDB >= max(1,N).
!> \endverbatim
!>
!> \param[out] X
!> \verbatim
!>          X is COMPLEX*16 array, dimension (LDX,NRHS)
!>     If INFO = 0, the N-by-NRHS solution matrix X to the original
!>     system of equations.  Note that A and B are modified on exit if
!>     EQUED .ne. 'N', and the solution to the equilibrated system is
!>     inv(diag(S))*X.
!> \endverbatim
!>
!> \param[in] LDX
!> \verbatim
!>          LDX is INTEGER
!>     The leading dimension of the array X.  LDX >= max(1,N).
!> \endverbatim
!>
!> \param[out] RCOND
!> \verbatim
!>          RCOND is DOUBLE PRECISION
!>     Reciprocal scaled condition number.  This is an estimate of the
!>     reciprocal Skeel condition number of the matrix A after
!>     equilibration (if done).  If this is less than the machine
!>     precision (in particular, if it is zero), the matrix is singular
!>     to working precision.  Note that the error may still be small even
!>     if this number is very small and the matrix appears ill-
!>     conditioned.
!> \endverbatim
!>
!> \param[out] RPVGRW
!> \verbatim
!>          RPVGRW is DOUBLE PRECISION
!>     Reciprocal pivot growth.  On exit, this contains the reciprocal
!>     pivot growth factor norm(A)/norm(U). The "max absolute element"
!>     norm is used.  If this is much less than 1, then the stability of
!>     the LU factorization of the (equilibrated) matrix A could be poor.
!>     This also means that the solution X, estimated condition numbers,
!>     and error bounds could be unreliable. If factorization fails with
!>     0<INFO<=N, then this contains the reciprocal pivot growth factor
!>     for the leading INFO columns of A.
!> \endverbatim
!>
!> \param[out] BERR
!> \verbatim
!>          BERR is DOUBLE PRECISION array, dimension (NRHS)
!>     Componentwise relative backward error.  This is the
!>     componentwise relative backward error of each solution vector X(j)
!>     (i.e., the smallest relative change in any element of A or B that
!>     makes X(j) an exact solution).
!> \endverbatim
!>
!> \param[in] N_ERR_BNDS
!> \verbatim
!>          N_ERR_BNDS is INTEGER
!>     Number of error bounds to return for each right hand side
!>     and each type (normwise or componentwise).  See ERR_BNDS_NORM and
!>     ERR_BNDS_COMP below.
!> \endverbatim
!>
!> \param[out] ERR_BNDS_NORM
!> \verbatim
!>          ERR_BNDS_NORM is DOUBLE PRECISION array, dimension (NRHS, N_ERR_BNDS)
!>     For each right-hand side, this array contains information about
!>     various error bounds and condition numbers corresponding to the
!>     normwise relative error, which is defined as follows:
!>
!>     Normwise relative error in the ith solution vector:
!>             max_j (abs(XTRUE(j,i) - X(j,i)))
!>            ------------------------------
!>                  max_j abs(X(j,i))
!>
!>     The array is indexed by the type of error information as described
!>     below. There currently are up to three pieces of information
!>     returned.
!>
!>     The first index in ERR_BNDS_NORM(i,:) corresponds to the ith
!>     right-hand side.
!>
!>     The second index in ERR_BNDS_NORM(:,err) contains the following
!>     three fields:
!>     err = 1 "Trust/don't trust" boolean. Trust the answer if the
!>              reciprocal condition number is less than the threshold
!>              sqrt(n) * dlamch('Epsilon').
!>
!>     err = 2 "Guaranteed" error bound: The estimated forward error,
!>              almost certainly within a factor of 10 of the true error
!>              so long as the next entry is greater than the threshold
!>              sqrt(n) * dlamch('Epsilon'). This error bound should only
!>              be trusted if the previous boolean is true.
!>
!>     err = 3  Reciprocal condition number: Estimated normwise
!>              reciprocal condition number.  Compared with the threshold
!>              sqrt(n) * dlamch('Epsilon') to determine if the error
!>              estimate is "guaranteed". These reciprocal condition
!>              numbers are 1 / (norm(Z^{-1},inf) * norm(Z,inf)) for some
!>              appropriately scaled matrix Z.
!>              Let Z = S*A, where S scales each row by a power of the
!>              radix so all absolute row sums of Z are approximately 1.
!>
!>     See Lapack Working Note 165 for further details and extra
!>     cautions.
!> \endverbatim
!>
!> \param[out] ERR_BNDS_COMP
!> \verbatim
!>          ERR_BNDS_COMP is DOUBLE PRECISION array, dimension (NRHS, N_ERR_BNDS)
!>     For each right-hand side, this array contains information about
!>     various error bounds and condition numbers corresponding to the
!>     componentwise relative error, which is defined as follows:
!>
!>     Componentwise relative error in the ith solution vector:
!>                    abs(XTRUE(j,i) - X(j,i))
!>             max_j ----------------------
!>                         abs(X(j,i))
!>
!>     The array is indexed by the right-hand side i (on which the
!>     componentwise relative error depends), and the type of error
!>     information as described below. There currently are up to three
!>     pieces of information returned for each right-hand side. If
!>     componentwise accuracy is not requested (PARAMS(3) = 0.0), then
!>     ERR_BNDS_COMP is not accessed.  If N_ERR_BNDS < 3, then at most
!>     the first (:,N_ERR_BNDS) entries are returned.
!>
!>     The first index in ERR_BNDS_COMP(i,:) corresponds to the ith
!>     right-hand side.
!>
!>     The second index in ERR_BNDS_COMP(:,err) contains the following
!>     three fields:
!>     err = 1 "Trust/don't trust" boolean. Trust the answer if the
!>              reciprocal condition number is less than the threshold
!>              sqrt(n) * dlamch('Epsilon').
!>
!>     err = 2 "Guaranteed" error bound: The estimated forward error,
!>              almost certainly within a factor of 10 of the true error
!>              so long as the next entry is greater than the threshold
!>              sqrt(n) * dlamch('Epsilon'). This error bound should only
!>              be trusted if the previous boolean is true.
!>
!>     err = 3  Reciprocal condition number: Estimated componentwise
!>              reciprocal condition number.  Compared with the threshold
!>              sqrt(n) * dlamch('Epsilon') to determine if the error
!>              estimate is "guaranteed". These reciprocal condition
!>              numbers are 1 / (norm(Z^{-1},inf) * norm(Z,inf)) for some
!>              appropriately scaled matrix Z.
!>              Let Z = S*(A*diag(x)), where x is the solution for the
!>              current right-hand side and S scales each row of
!>              A*diag(x) by a power of the radix so all absolute row
!>              sums of Z are approximately 1.
!>
!>     See Lapack Working Note 165 for further details and extra
!>     cautions.
!> \endverbatim
!>
!> \param[in] NPARAMS
!> \verbatim
!>          NPARAMS is INTEGER
!>     Specifies the number of parameters set in PARAMS.  If <= 0, the
!>     PARAMS array is never referenced and default values are used.
!> \endverbatim
!>
!> \param[in,out] PARAMS
!> \verbatim
!>          PARAMS is DOUBLE PRECISION array, dimension NPARAMS
!>     Specifies algorithm parameters.  If an entry is < 0.0, then
!>     that entry will be filled with default value used for that
!>     parameter.  Only positions up to NPARAMS are accessed; defaults
!>     are used for higher-numbered parameters.
!>
!>       PARAMS(LA_LINRX_ITREF_I = 1) : Whether to perform iterative
!>            refinement or not.
!>         Default: 1.0D+0
!>            = 0.0:  No refinement is performed, and no error bounds are
!>                    computed.
!>            = 1.0:  Use the extra-precise refinement algorithm.
!>              (other values are reserved for future use)
!>
!>       PARAMS(LA_LINRX_ITHRESH_I = 2) : Maximum number of residual
!>            computations allowed for refinement.
!>         Default: 10
!>         Aggressive: Set to 100 to permit convergence using approximate
!>                     factorizations or factorizations other than LU. If
!>                     the factorization uses a technique other than
!>                     Gaussian elimination, the guarantees in
!>                     err_bnds_norm and err_bnds_comp may no longer be
!>                     trustworthy.
!>
!>       PARAMS(LA_LINRX_CWISE_I = 3) : Flag determining if the code
!>            will attempt to find a solution with small componentwise
!>            relative error in the double-precision algorithm.  Positive
!>            is true, 0.0 is false.
!>         Default: 1.0 (attempt componentwise convergence)
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX*16 array, dimension (2*N)
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is DOUBLE PRECISION array, dimension (2*N)
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>       = 0:  Successful exit. The solution to every right-hand side is
!>         guaranteed.
!>       < 0:  If INFO = -i, the i-th argument had an illegal value
!>       > 0 and <= N:  U(INFO,INFO) is exactly zero.  The factorization
!>         has been completed, but the factor U is exactly singular, so
!>         the solution and error bounds could not be computed. RCOND = 0
!>         is returned.
!>       = N+J: The solution corresponding to the Jth right-hand side is
!>         not guaranteed. The solutions corresponding to other right-
!>         hand sides K with K > J may not be guaranteed as well, but
!>         only the first such right-hand side is reported. If a small
!>         componentwise error is not requested (PARAMS(3) = 0.0) then
!>         the Jth right-hand side is the first with a normwise error
!>         bound that is not guaranteed (the smallest J such
!>         that ERR_BNDS_NORM(J,1) = 0.0). By default (PARAMS(3) = 1.0)
!>         the Jth right-hand side is the first with either a normwise or
!>         componentwise error bound that is not guaranteed (the smallest
!>         J such that either ERR_BNDS_NORM(J,1) = 0.0 or
!>         ERR_BNDS_COMP(J,1) = 0.0). See the definition of
!>         ERR_BNDS_NORM(:,1) and ERR_BNDS_COMP(:,1). To get information
!>         about all of the right-hand sides check ERR_BNDS_NORM or
!>         ERR_BNDS_COMP.
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
!> \date April 2012
!
!> \ingroup complex16SYsolve
!
!  =====================================================================
      SUBROUTINE ZSYSVXX(Fact,Uplo,N,Nrhs,A,Lda,Af,Ldaf,Ipiv,Equed,S,B, &
     &                   Ldb,X,Ldx,Rcond,Rpvgrw,Berr,N_err_bnds,        &
     &                   Err_bnds_norm,Err_bnds_comp,Nparams,Params,    &
     &                   Work,Rwork,Info)
      IMPLICIT NONE
!*--ZSYSVXX510
!
!  -- LAPACK driver routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     April 2012
!
!     .. Scalar Arguments ..
      CHARACTER Equed , Fact , Uplo
      INTEGER Info , Lda , Ldaf , Ldb , Ldx , N , Nrhs , Nparams ,      &
     &        N_err_bnds
      DOUBLE PRECISION Rcond , Rpvgrw
!     ..
!     .. Array Arguments ..
      INTEGER Ipiv(*)
      COMPLEX*16 A(Lda,*) , Af(Ldaf,*) , B(Ldb,*) , X(Ldx,*) , Work(*)
      DOUBLE PRECISION S(*) , Params(*) , Berr(*) ,                     &
     &                 Err_bnds_norm(Nrhs,*) , Err_bnds_comp(Nrhs,*) ,  &
     &                 Rwork(*)
!     ..
!
!  ==================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO , ONE
      PARAMETER (ZERO=0.0D+0,ONE=1.0D+0)
!     ..
!     .. Local Scalars ..
      LOGICAL equil , nofact , rcequ
      INTEGER infequ , j
      DOUBLE PRECISION amax , bignum , smin , smax , scond , smlnum
!     ..
!     .. External Functions ..
      EXTERNAL LSAME , DLAMCH , ZLA_SYRPVGRW
      LOGICAL LSAME
      DOUBLE PRECISION DLAMCH , ZLA_SYRPVGRW
!     ..
!     .. External Subroutines ..
      EXTERNAL ZSYEQUB , ZSYTRF , ZSYTRS , ZLACPY , ZLAQSY , XERBLA ,   &
     &         ZLASCL2 , ZSYRFSX
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX , MIN
!     ..
!     .. Executable Statements ..
!
      Info = 0
      nofact = LSAME(Fact,'N')
      equil = LSAME(Fact,'E')
      smlnum = DLAMCH('Safe minimum')
      bignum = ONE/smlnum
      IF ( nofact .OR. equil ) THEN
         Equed = 'N'
         rcequ = .FALSE.
      ELSE
         rcequ = LSAME(Equed,'Y')
      ENDIF
!
!     Default is failure.  If an input parameter is wrong or
!     factorization fails, make everything look horrible.  Only the
!     pivot growth is set here, the rest is initialized in ZSYRFSX.
!
      Rpvgrw = ZERO
!
!     Test the input parameters.  PARAMS is not tested until ZSYRFSX.
!
      IF ( .NOT.nofact .AND. .NOT.equil .AND. .NOT.LSAME(Fact,'F') )    &
     &     THEN
         Info = -1
      ELSEIF ( .NOT.LSAME(Uplo,'U') .AND. .NOT.LSAME(Uplo,'L') ) THEN
         Info = -2
      ELSEIF ( N<0 ) THEN
         Info = -3
      ELSEIF ( Nrhs<0 ) THEN
         Info = -4
      ELSEIF ( Lda<MAX(1,N) ) THEN
         Info = -6
      ELSEIF ( Ldaf<MAX(1,N) ) THEN
         Info = -8
      ELSEIF ( LSAME(Fact,'F') .AND. .NOT.(rcequ .OR. LSAME(Equed,'N')) &
     &         ) THEN
         Info = -10
      ELSE
         IF ( rcequ ) THEN
            smin = bignum
            smax = ZERO
            DO j = 1 , N
               smin = MIN(smin,S(j))
               smax = MAX(smax,S(j))
            ENDDO
            IF ( smin<=ZERO ) THEN
               Info = -11
            ELSEIF ( N>0 ) THEN
               scond = MAX(smin,smlnum)/MIN(smax,bignum)
            ELSE
               scond = ONE
            ENDIF
         ENDIF
         IF ( Info==0 ) THEN
            IF ( Ldb<MAX(1,N) ) THEN
               Info = -13
            ELSEIF ( Ldx<MAX(1,N) ) THEN
               Info = -15
            ENDIF
         ENDIF
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('ZSYSVXX',-Info)
         RETURN
      ENDIF
!
      IF ( equil ) THEN
!
!     Compute row and column scalings to equilibrate the matrix A.
!
         CALL ZSYEQUB(Uplo,N,A,Lda,S,scond,amax,Work,infequ)
         IF ( infequ==0 ) THEN
!
!     Equilibrate the matrix.
!
            CALL ZLAQSY(Uplo,N,A,Lda,S,scond,amax,Equed)
            rcequ = LSAME(Equed,'Y')
         ENDIF
 
      ENDIF
!
!     Scale the right hand-side.
!
      IF ( rcequ ) CALL ZLASCL2(N,Nrhs,S,B,Ldb)
!
      IF ( nofact .OR. equil ) THEN
!
!        Compute the LDL^T or UDU^T factorization of A.
!
         CALL ZLACPY(Uplo,N,N,A,Lda,Af,Ldaf)
         CALL ZSYTRF(Uplo,N,Af,Ldaf,Ipiv,Work,5*MAX(1,N),Info)
!
!        Return if INFO is non-zero.
!
         IF ( Info>0 ) THEN
!
!           Pivot in column INFO is exactly 0
!           Compute the reciprocal pivot growth factor of the
!           leading rank-deficient INFO columns of A.
!
            IF ( N>0 ) Rpvgrw = ZLA_SYRPVGRW(Uplo,N,Info,A,Lda,Af,Ldaf, &
     &                          Ipiv,Rwork)
            RETURN
         ENDIF
      ENDIF
!
!     Compute the reciprocal pivot growth factor RPVGRW.
!
      IF ( N>0 ) Rpvgrw = ZLA_SYRPVGRW(Uplo,N,Info,A,Lda,Af,Ldaf,Ipiv,  &
     &                    Rwork)
!
!     Compute the solution matrix X.
!
      CALL ZLACPY('Full',N,Nrhs,B,Ldb,X,Ldx)
      CALL ZSYTRS(Uplo,N,Nrhs,Af,Ldaf,Ipiv,X,Ldx,Info)
!
!     Use iterative refinement to improve the computed solution and
!     compute error bounds and backward error estimates for it.
!
      CALL ZSYRFSX(Uplo,Equed,N,Nrhs,A,Lda,Af,Ldaf,Ipiv,S,B,Ldb,X,Ldx,  &
     &             Rcond,Berr,N_err_bnds,Err_bnds_norm,Err_bnds_comp,   &
     &             Nparams,Params,Work,Rwork,Info)
!
!     Scale solutions.
!
      IF ( rcequ ) CALL ZLASCL2(N,Nrhs,S,X,Ldx)
!
!
!     End of ZSYSVXX
!
      END SUBROUTINE ZSYSVXX
