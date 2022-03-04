!*==zgesvxx.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief <b> ZGESVXX computes the solution to system of linear equations A * X = B for GE matrices</b>
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZGESVXX + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zgesvxx.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zgesvxx.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zgesvxx.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZGESVXX( FACT, TRANS, N, NRHS, A, LDA, AF, LDAF, IPIV,
!                           EQUED, R, C, B, LDB, X, LDX, RCOND, RPVGRW,
!                           BERR, N_ERR_BNDS, ERR_BNDS_NORM,
!                           ERR_BNDS_COMP, NPARAMS, PARAMS, WORK, RWORK,
!                           INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          EQUED, FACT, TRANS
!       INTEGER            INFO, LDA, LDAF, LDB, LDX, N, NRHS, NPARAMS,
!      $                   N_ERR_BNDS
!       DOUBLE PRECISION   RCOND, RPVGRW
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * )
!       COMPLEX*16         A( LDA, * ), AF( LDAF, * ), B( LDB, * ),
!      $                   X( LDX , * ),WORK( * )
!       DOUBLE PRECISION   R( * ), C( * ), PARAMS( * ), BERR( * ),
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
!>    ZGESVXX uses the LU factorization to compute the solution to a
!>    complex*16 system of linear equations  A * X = B,  where A is an
!>    N-by-N matrix and X and B are N-by-NRHS matrices.
!>
!>    If requested, both normwise and maximum componentwise error bounds
!>    are returned. ZGESVXX will return a solution with a tiny
!>    guaranteed error (O(eps) where eps is the working machine
!>    precision) unless the matrix is very ill-conditioned, in which
!>    case a warning is returned. Relevant condition numbers also are
!>    calculated and returned.
!>
!>    ZGESVXX accepts user-provided factorizations and equilibration
!>    factors; see the definitions of the FACT and EQUED options.
!>    Solving with refinement and using a factorization from a previous
!>    ZGESVXX call will also produce a solution with either O(eps)
!>    errors or warnings, but we cannot make that claim for general
!>    user-provided factorizations and equilibration factors if they
!>    differ from what ZGESVXX would itself produce.
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
!>      TRANS = 'N':  diag(R)*A*diag(C)     *inv(diag(C))*X = diag(R)*B
!>      TRANS = 'T': (diag(R)*A*diag(C))**T *inv(diag(R))*X = diag(C)*B
!>      TRANS = 'C': (diag(R)*A*diag(C))**H *inv(diag(R))*X = diag(C)*B
!>
!>    Whether or not the system will be equilibrated depends on the
!>    scaling of the matrix A, but if equilibration is used, A is
!>    overwritten by diag(R)*A*diag(C) and B by diag(R)*B (if TRANS='N')
!>    or diag(C)*B (if TRANS = 'T' or 'C').
!>
!>    2. If FACT = 'N' or 'E', the LU decomposition is used to factor
!>    the matrix A (after equilibration if FACT = 'E') as
!>
!>      A = P * L * U,
!>
!>    where P is a permutation matrix, L is a unit lower triangular
!>    matrix, and U is upper triangular.
!>
!>    3. If some U(i,i)=0, so that U is exactly singular, then the
!>    routine returns with INFO = i. Otherwise, the factored form of A
!>    is used to estimate the condition number of the matrix A (see
!>    argument RCOND). If the reciprocal of the condition number is less
!>    than machine precision, the routine still goes on to solve for X
!>    and compute error bounds as described below.
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
!>    diag(C) (if TRANS = 'N') or diag(R) (if TRANS = 'T' or 'C') so
!>    that it solves the original system before equilibration.
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
!>               equilibrated with scaling factors given by R and C.
!>               A, AF, and IPIV are not modified.
!>       = 'N':  The matrix A will be copied to AF and factored.
!>       = 'E':  The matrix A will be equilibrated if necessary, then
!>               copied to AF and factored.
!> \endverbatim
!>
!> \param[in] TRANS
!> \verbatim
!>          TRANS is CHARACTER*1
!>     Specifies the form of the system of equations:
!>       = 'N':  A * X = B     (No transpose)
!>       = 'T':  A**T * X = B  (Transpose)
!>       = 'C':  A**H * X = B  (Conjugate Transpose)
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
!>     On entry, the N-by-N matrix A.  If FACT = 'F' and EQUED is
!>     not 'N', then A must have been equilibrated by the scaling
!>     factors in R and/or C.  A is not modified if FACT = 'F' or
!>     'N', or if FACT = 'E' and EQUED = 'N' on exit.
!>
!>     On exit, if EQUED .ne. 'N', A is scaled as follows:
!>     EQUED = 'R':  A := diag(R) * A
!>     EQUED = 'C':  A := A * diag(C)
!>     EQUED = 'B':  A := diag(R) * A * diag(C).
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
!>     contains the factors L and U from the factorization
!>     A = P*L*U as computed by ZGETRF.  If EQUED .ne. 'N', then
!>     AF is the factored form of the equilibrated matrix A.
!>
!>     If FACT = 'N', then AF is an output argument and on exit
!>     returns the factors L and U from the factorization A = P*L*U
!>     of the original matrix A.
!>
!>     If FACT = 'E', then AF is an output argument and on exit
!>     returns the factors L and U from the factorization A = P*L*U
!>     of the equilibrated matrix A (see the description of A for
!>     the form of the equilibrated matrix).
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
!>     contains the pivot indices from the factorization A = P*L*U
!>     as computed by ZGETRF; row i of the matrix was interchanged
!>     with row IPIV(i).
!>
!>     If FACT = 'N', then IPIV is an output argument and on exit
!>     contains the pivot indices from the factorization A = P*L*U
!>     of the original matrix A.
!>
!>     If FACT = 'E', then IPIV is an output argument and on exit
!>     contains the pivot indices from the factorization A = P*L*U
!>     of the equilibrated matrix A.
!> \endverbatim
!>
!> \param[in,out] EQUED
!> \verbatim
!>          EQUED is CHARACTER*1
!>     Specifies the form of equilibration that was done.
!>       = 'N':  No equilibration (always true if FACT = 'N').
!>       = 'R':  Row equilibration, i.e., A has been premultiplied by
!>               diag(R).
!>       = 'C':  Column equilibration, i.e., A has been postmultiplied
!>               by diag(C).
!>       = 'B':  Both row and column equilibration, i.e., A has been
!>               replaced by diag(R) * A * diag(C).
!>     EQUED is an input argument if FACT = 'F'; otherwise, it is an
!>     output argument.
!> \endverbatim
!>
!> \param[in,out] R
!> \verbatim
!>          R is DOUBLE PRECISION array, dimension (N)
!>     The row scale factors for A.  If EQUED = 'R' or 'B', A is
!>     multiplied on the left by diag(R); if EQUED = 'N' or 'C', R
!>     is not accessed.  R is an input argument if FACT = 'F';
!>     otherwise, R is an output argument.  If FACT = 'F' and
!>     EQUED = 'R' or 'B', each element of R must be positive.
!>     If R is output, each element of R is a power of the radix.
!>     If R is input, each element of R should be a power of the radix
!>     to ensure a reliable solution and error estimates. Scaling by
!>     powers of the radix does not cause rounding errors unless the
!>     result underflows or overflows. Rounding errors during scaling
!>     lead to refining with a matrix that is not equivalent to the
!>     input matrix, producing error estimates that may not be
!>     reliable.
!> \endverbatim
!>
!> \param[in,out] C
!> \verbatim
!>          C is DOUBLE PRECISION array, dimension (N)
!>     The column scale factors for A.  If EQUED = 'C' or 'B', A is
!>     multiplied on the right by diag(C); if EQUED = 'N' or 'R', C
!>     is not accessed.  C is an input argument if FACT = 'F';
!>     otherwise, C is an output argument.  If FACT = 'F' and
!>     EQUED = 'C' or 'B', each element of C must be positive.
!>     If C is output, each element of C is a power of the radix.
!>     If C is input, each element of C should be a power of the radix
!>     to ensure a reliable solution and error estimates. Scaling by
!>     powers of the radix does not cause rounding errors unless the
!>     result underflows or overflows. Rounding errors during scaling
!>     lead to refining with a matrix that is not equivalent to the
!>     input matrix, producing error estimates that may not be
!>     reliable.
!> \endverbatim
!>
!> \param[in,out] B
!> \verbatim
!>          B is COMPLEX*16 array, dimension (LDB,NRHS)
!>     On entry, the N-by-NRHS right hand side matrix B.
!>     On exit,
!>     if EQUED = 'N', B is not modified;
!>     if TRANS = 'N' and EQUED = 'R' or 'B', B is overwritten by
!>        diag(R)*B;
!>     if TRANS = 'T' or 'C' and EQUED = 'C' or 'B', B is
!>        overwritten by diag(C)*B.
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
!>     system of equations.  Note that A and B are modified on exit
!>     if EQUED .ne. 'N', and the solution to the equilibrated system is
!>     inv(diag(C))*X if TRANS = 'N' and EQUED = 'C' or 'B', or
!>     inv(diag(R))*X if TRANS = 'T' or 'C' and EQUED = 'R' or 'B'.
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
!>     for the leading INFO columns of A.  In ZGESVX, this quantity is
!>     returned in WORK(1).
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
!> \ingroup complex16GEsolve
!
!  =====================================================================
      SUBROUTINE ZGESVXX(Fact,Trans,N,Nrhs,A,Lda,Af,Ldaf,Ipiv,Equed,R,C,&
     &                   B,Ldb,X,Ldx,Rcond,Rpvgrw,Berr,N_err_bnds,      &
     &                   Err_bnds_norm,Err_bnds_comp,Nparams,Params,    &
     &                   Work,Rwork,Info)
      USE F77KINDS                        
      USE S_DLAMCH
      USE S_LSAME
      USE S_XERBLA
      USE S_ZGEEQUB
      USE S_ZGERFSX
      USE S_ZGETRF
      USE S_ZGETRS
      USE S_ZLACPY
      USE S_ZLAQGE
      USE S_ZLASCL2
      USE S_ZLA_GERPVGRW
      IMPLICIT NONE
!*--ZGESVXX555
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Fact
      CHARACTER :: Trans
      INTEGER :: N
      INTEGER :: Nrhs
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(Ldaf,*) :: Af
      INTEGER :: Ldaf
      INTEGER , DIMENSION(*) :: Ipiv
      CHARACTER :: Equed
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: R
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: C
      COMPLEX(CX16KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX(CX16KIND) , DIMENSION(Ldx,*) :: X
      INTEGER :: Ldx
      REAL(R8KIND) :: Rcond
      REAL(R8KIND) , INTENT(OUT) :: Rpvgrw
      REAL(R8KIND) , DIMENSION(*) :: Berr
      INTEGER :: N_err_bnds
      REAL(R8KIND) , DIMENSION(Nrhs,*) :: Err_bnds_norm
      REAL(R8KIND) , DIMENSION(Nrhs,*) :: Err_bnds_comp
      INTEGER :: Nparams
      REAL(R8KIND) , DIMENSION(*) :: Params
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      REAL(R8KIND) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      REAL(R8KIND) :: amax , bignum , colcnd , rcmax , rcmin , rowcnd , &
     &                smlnum
      LOGICAL :: colequ , equil , nofact , notran , rowequ
      INTEGER :: infequ , j
!
! End of declarations rewritten by SPAG
!
!     ..
!     .. Array Arguments ..
!     ..
!
!  ==================================================================
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
      equil = LSAME(Fact,'E')
      notran = LSAME(Trans,'N')
      smlnum = DLAMCH('Safe minimum')
      bignum = ONE/smlnum
      IF ( nofact .OR. equil ) THEN
         Equed = 'N'
         rowequ = .FALSE.
         colequ = .FALSE.
      ELSE
         rowequ = LSAME(Equed,'R') .OR. LSAME(Equed,'B')
         colequ = LSAME(Equed,'C') .OR. LSAME(Equed,'B')
      ENDIF
!
!     Default is failure.  If an input parameter is wrong or
!     factorization fails, make everything look horrible.  Only the
!     pivot growth is set here, the rest is initialized in ZGERFSX.
!
      Rpvgrw = ZERO
!
!     Test the input parameters.  PARAMS is not tested until ZGERFSX.
!
      IF ( .NOT.nofact .AND. .NOT.equil .AND. .NOT.LSAME(Fact,'F') )    &
     &     THEN
         Info = -1
      ELSEIF ( .NOT.notran .AND. .NOT.LSAME(Trans,'T') .AND.            &
     &         .NOT.LSAME(Trans,'C') ) THEN
         Info = -2
      ELSEIF ( N<0 ) THEN
         Info = -3
      ELSEIF ( Nrhs<0 ) THEN
         Info = -4
      ELSEIF ( Lda<MAX(1,N) ) THEN
         Info = -6
      ELSEIF ( Ldaf<MAX(1,N) ) THEN
         Info = -8
      ELSEIF ( LSAME(Fact,'F') .AND.                                    &
     &         .NOT.(rowequ .OR. colequ .OR. LSAME(Equed,'N')) ) THEN
         Info = -10
      ELSE
         IF ( rowequ ) THEN
            rcmin = bignum
            rcmax = ZERO
            DO j = 1 , N
               rcmin = MIN(rcmin,R(j))
               rcmax = MAX(rcmax,R(j))
            ENDDO
            IF ( rcmin<=ZERO ) THEN
               Info = -11
            ELSEIF ( N>0 ) THEN
               rowcnd = MAX(rcmin,smlnum)/MIN(rcmax,bignum)
            ELSE
               rowcnd = ONE
            ENDIF
         ENDIF
         IF ( colequ .AND. Info==0 ) THEN
            rcmin = bignum
            rcmax = ZERO
            DO j = 1 , N
               rcmin = MIN(rcmin,C(j))
               rcmax = MAX(rcmax,C(j))
            ENDDO
            IF ( rcmin<=ZERO ) THEN
               Info = -12
            ELSEIF ( N>0 ) THEN
               colcnd = MAX(rcmin,smlnum)/MIN(rcmax,bignum)
            ELSE
               colcnd = ONE
            ENDIF
         ENDIF
         IF ( Info==0 ) THEN
            IF ( Ldb<MAX(1,N) ) THEN
               Info = -14
            ELSEIF ( Ldx<MAX(1,N) ) THEN
               Info = -16
            ENDIF
         ENDIF
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('ZGESVXX',-Info)
         RETURN
      ENDIF
!
      IF ( equil ) THEN
!
!     Compute row and column scalings to equilibrate the matrix A.
!
         CALL ZGEEQUB(N,N,A,Lda,R,C,rowcnd,colcnd,amax,infequ)
         IF ( infequ==0 ) THEN
!
!     Equilibrate the matrix.
!
            CALL ZLAQGE(N,N,A,Lda,R,C,rowcnd,colcnd,amax,Equed)
            rowequ = LSAME(Equed,'R') .OR. LSAME(Equed,'B')
            colequ = LSAME(Equed,'C') .OR. LSAME(Equed,'B')
         ENDIF
!
!     If the scaling factors are not applied, set them to 1.0.
!
         IF ( .NOT.rowequ ) THEN
            DO j = 1 , N
               R(j) = 1.0D+0
            ENDDO
         ENDIF
         IF ( .NOT.colequ ) THEN
            DO j = 1 , N
               C(j) = 1.0D+0
            ENDDO
         ENDIF
      ENDIF
!
!     Scale the right-hand side.
!
      IF ( notran ) THEN
         IF ( rowequ ) CALL ZLASCL2(N,Nrhs,R,B,Ldb)
      ELSE
         IF ( colequ ) CALL ZLASCL2(N,Nrhs,C,B,Ldb)
      ENDIF
!
      IF ( nofact .OR. equil ) THEN
!
!        Compute the LU factorization of A.
!
         CALL ZLACPY('Full',N,N,A,Lda,Af,Ldaf)
         CALL ZGETRF(N,N,Af,Ldaf,Ipiv,Info)
!
!        Return if INFO is non-zero.
!
         IF ( Info>0 ) THEN
!
!           Pivot in column INFO is exactly 0
!           Compute the reciprocal pivot growth factor of the
!           leading rank-deficient INFO columns of A.
!
            Rpvgrw = ZLA_GERPVGRW(N,Info,A,Lda,Af,Ldaf)
            RETURN
         ENDIF
      ENDIF
!
!     Compute the reciprocal pivot growth factor RPVGRW.
!
      Rpvgrw = ZLA_GERPVGRW(N,N,A,Lda,Af,Ldaf)
!
!     Compute the solution matrix X.
!
      CALL ZLACPY('Full',N,Nrhs,B,Ldb,X,Ldx)
      CALL ZGETRS(Trans,N,Nrhs,Af,Ldaf,Ipiv,X,Ldx,Info)
!
!     Use iterative refinement to improve the computed solution and
!     compute error bounds and backward error estimates for it.
!
      CALL ZGERFSX(Trans,Equed,N,Nrhs,A,Lda,Af,Ldaf,Ipiv,R,C,B,Ldb,X,   &
     &             Ldx,Rcond,Berr,N_err_bnds,Err_bnds_norm,             &
     &             Err_bnds_comp,Nparams,Params,Work,Rwork,Info)
!
!     Scale solutions.
!
      IF ( colequ .AND. notran ) THEN
         CALL ZLASCL2(N,Nrhs,C,X,Ldx)
      ELSEIF ( rowequ .AND. .NOT.notran ) THEN
         CALL ZLASCL2(N,Nrhs,R,X,Ldx)
      ENDIF
!
!
!     End of ZGESVXX
!
      END SUBROUTINE ZGESVXX
