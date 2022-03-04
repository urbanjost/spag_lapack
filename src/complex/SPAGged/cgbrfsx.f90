!*==cgbrfsx.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b CGBRFSX
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CGBRFSX + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgbrfsx.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgbrfsx.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgbrfsx.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CGBRFSX( TRANS, EQUED, N, KL, KU, NRHS, AB, LDAB, AFB,
!                           LDAFB, IPIV, R, C, B, LDB, X, LDX, RCOND,
!                           BERR, N_ERR_BNDS, ERR_BNDS_NORM,
!                           ERR_BNDS_COMP, NPARAMS, PARAMS, WORK, RWORK,
!                           INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          TRANS, EQUED
!       INTEGER            INFO, LDAB, LDAFB, LDB, LDX, N, KL, KU, NRHS,
!      $                   NPARAMS, N_ERR_BNDS
!       REAL               RCOND
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * )
!       COMPLEX            AB( LDAB, * ), AFB( LDAFB, * ), B( LDB, * ),
!      $                   X( LDX , * ),WORK( * )
!       REAL               R( * ), C( * ), PARAMS( * ), BERR( * ),
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
!>    CGBRFSX improves the computed solution to a system of linear
!>    equations and provides error bounds and backward error estimates
!>    for the solution.  In addition to normwise error bound, the code
!>    provides maximum componentwise error bound if possible.  See
!>    comments for ERR_BNDS_NORM and ERR_BNDS_COMP for details of the
!>    error bounds.
!>
!>    The original system of linear equations may have been equilibrated
!>    before calling this routine, as described by arguments EQUED, R
!>    and C below. In this case, the solution and error bounds returned
!>    are for the original unequilibrated system.
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
!> \param[in] TRANS
!> \verbatim
!>          TRANS is CHARACTER*1
!>     Specifies the form of the system of equations:
!>       = 'N':  A * X = B     (No transpose)
!>       = 'T':  A**T * X = B  (Transpose)
!>       = 'C':  A**H * X = B  (Conjugate transpose)
!> \endverbatim
!>
!> \param[in] EQUED
!> \verbatim
!>          EQUED is CHARACTER*1
!>     Specifies the form of equilibration that was done to A
!>     before calling this routine. This is needed to compute
!>     the solution and error bounds correctly.
!>       = 'N':  No equilibration
!>       = 'R':  Row equilibration, i.e., A has been premultiplied by
!>               diag(R).
!>       = 'C':  Column equilibration, i.e., A has been postmultiplied
!>               by diag(C).
!>       = 'B':  Both row and column equilibration, i.e., A has been
!>               replaced by diag(R) * A * diag(C).
!>               The right hand side B has been changed accordingly.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>     The order of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in] KL
!> \verbatim
!>          KL is INTEGER
!>     The number of subdiagonals within the band of A.  KL >= 0.
!> \endverbatim
!>
!> \param[in] KU
!> \verbatim
!>          KU is INTEGER
!>     The number of superdiagonals within the band of A.  KU >= 0.
!> \endverbatim
!>
!> \param[in] NRHS
!> \verbatim
!>          NRHS is INTEGER
!>     The number of right hand sides, i.e., the number of columns
!>     of the matrices B and X.  NRHS >= 0.
!> \endverbatim
!>
!> \param[in] AB
!> \verbatim
!>          AB is COMPLEX array, dimension (LDAB,N)
!>     The original band matrix A, stored in rows 1 to KL+KU+1.
!>     The j-th column of A is stored in the j-th column of the
!>     array AB as follows:
!>     AB(ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(n,j+kl).
!> \endverbatim
!>
!> \param[in] LDAB
!> \verbatim
!>          LDAB is INTEGER
!>     The leading dimension of the array AB.  LDAB >= KL+KU+1.
!> \endverbatim
!>
!> \param[in] AFB
!> \verbatim
!>          AFB is COMPLEX array, dimension (LDAFB,N)
!>     Details of the LU factorization of the band matrix A, as
!>     computed by DGBTRF.  U is stored as an upper triangular band
!>     matrix with KL+KU superdiagonals in rows 1 to KL+KU+1, and
!>     the multipliers used during the factorization are stored in
!>     rows KL+KU+2 to 2*KL+KU+1.
!> \endverbatim
!>
!> \param[in] LDAFB
!> \verbatim
!>          LDAFB is INTEGER
!>     The leading dimension of the array AFB.  LDAFB >= 2*KL*KU+1.
!> \endverbatim
!>
!> \param[in] IPIV
!> \verbatim
!>          IPIV is INTEGER array, dimension (N)
!>     The pivot indices from SGETRF; for 1<=i<=N, row i of the
!>     matrix was interchanged with row IPIV(i).
!> \endverbatim
!>
!> \param[in,out] R
!> \verbatim
!>          R is REAL array, dimension (N)
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
!>          C is REAL array, dimension (N)
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
!> \param[in] B
!> \verbatim
!>          B is COMPLEX array, dimension (LDB,NRHS)
!>     The right hand side matrix B.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>     The leading dimension of the array B.  LDB >= max(1,N).
!> \endverbatim
!>
!> \param[in,out] X
!> \verbatim
!>          X is COMPLEX array, dimension (LDX,NRHS)
!>     On entry, the solution matrix X, as computed by SGETRS.
!>     On exit, the improved solution matrix X.
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
!>          RCOND is REAL
!>     Reciprocal scaled condition number.  This is an estimate of the
!>     reciprocal Skeel condition number of the matrix A after
!>     equilibration (if done).  If this is less than the machine
!>     precision (in particular, if it is zero), the matrix is singular
!>     to working precision.  Note that the error may still be small even
!>     if this number is very small and the matrix appears ill-
!>     conditioned.
!> \endverbatim
!>
!> \param[out] BERR
!> \verbatim
!>          BERR is REAL array, dimension (NRHS)
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
!>          ERR_BNDS_NORM is REAL array, dimension (NRHS, N_ERR_BNDS)
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
!>              sqrt(n) * slamch('Epsilon').
!>
!>     err = 2 "Guaranteed" error bound: The estimated forward error,
!>              almost certainly within a factor of 10 of the true error
!>              so long as the next entry is greater than the threshold
!>              sqrt(n) * slamch('Epsilon'). This error bound should only
!>              be trusted if the previous boolean is true.
!>
!>     err = 3  Reciprocal condition number: Estimated normwise
!>              reciprocal condition number.  Compared with the threshold
!>              sqrt(n) * slamch('Epsilon') to determine if the error
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
!>          ERR_BNDS_COMP is REAL array, dimension (NRHS, N_ERR_BNDS)
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
!>              sqrt(n) * slamch('Epsilon').
!>
!>     err = 2 "Guaranteed" error bound: The estimated forward error,
!>              almost certainly within a factor of 10 of the true error
!>              so long as the next entry is greater than the threshold
!>              sqrt(n) * slamch('Epsilon'). This error bound should only
!>              be trusted if the previous boolean is true.
!>
!>     err = 3  Reciprocal condition number: Estimated componentwise
!>              reciprocal condition number.  Compared with the threshold
!>              sqrt(n) * slamch('Epsilon') to determine if the error
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
!>          PARAMS is REAL array, dimension NPARAMS
!>     Specifies algorithm parameters.  If an entry is < 0.0, then
!>     that entry will be filled with default value used for that
!>     parameter.  Only positions up to NPARAMS are accessed; defaults
!>     are used for higher-numbered parameters.
!>
!>       PARAMS(LA_LINRX_ITREF_I = 1) : Whether to perform iterative
!>            refinement or not.
!>         Default: 1.0
!>            = 0.0:  No refinement is performed, and no error bounds are
!>                    computed.
!>            = 1.0:  Use the double-precision refinement algorithm,
!>                    possibly with doubled-single computations if the
!>                    compilation environment does not support DOUBLE
!>                    PRECISION.
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
!>          WORK is COMPLEX array, dimension (2*N)
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is REAL array, dimension (2*N)
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
!> \ingroup complexGBcomputational
!
!  =====================================================================
      SUBROUTINE CGBRFSX(Trans,Equed,N,Kl,Ku,Nrhs,Ab,Ldab,Afb,Ldafb,    &
     &                   Ipiv,R,C,B,Ldb,X,Ldx,Rcond,Berr,N_err_bnds,    &
     &                   Err_bnds_norm,Err_bnds_comp,Nparams,Params,    &
     &                   Work,Rwork,Info)
      USE S_CGBCON
      USE S_CLANGB
      USE S_CLA_GBRCOND_C
      USE S_CLA_GBRCOND_X
      USE S_CLA_GBRFSX_EXTENDED
      USE S_ILAPREC
      USE S_ILATRANS
      USE S_LSAME
      USE S_SLAMCH
      USE S_XERBLA
      IMPLICIT NONE
!*--CGBRFSX453
!
! PARAMETER definitions rewritten by SPAG
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ITREF_DEFAULT = 1.0 ,       &
     &                      ITHRESH_DEFAULT = 10.0 ,                    &
     &                      COMPONENTWISE_DEFAULT = 1.0 ,               &
     &                      RTHRESH_DEFAULT = 0.5 ,                     &
     &                      DZTHRESH_DEFAULT = 0.25
      INTEGER , PARAMETER  ::  LA_LINRX_ITREF_I = 1 ,                   &
     &                         LA_LINRX_ITHRESH_I = 2 ,                 &
     &                         LA_LINRX_CWISE_I = 3 ,                   &
     &                         LA_LINRX_TRUST_I = 1 ,                   &
     &                         LA_LINRX_ERR_I = 2 , LA_LINRX_RCOND_I = 3
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Trans
      CHARACTER :: Equed
      INTEGER :: N
      INTEGER :: Kl
      INTEGER :: Ku
      INTEGER :: Nrhs
      COMPLEX , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      COMPLEX , DIMENSION(Ldafb,*) :: Afb
      INTEGER :: Ldafb
      INTEGER , DIMENSION(*) :: Ipiv
      REAL , DIMENSION(*) :: R
      REAL , DIMENSION(*) :: C
      COMPLEX , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX , DIMENSION(Ldx,*) :: X
      INTEGER :: Ldx
      REAL :: Rcond
      REAL , DIMENSION(*) :: Berr
      INTEGER , INTENT(IN) :: N_err_bnds
      REAL , INTENT(INOUT) , DIMENSION(Nrhs,*) :: Err_bnds_norm
      REAL , INTENT(INOUT) , DIMENSION(Nrhs,*) :: Err_bnds_comp
      INTEGER , INTENT(IN) :: Nparams
      REAL , INTENT(INOUT) , DIMENSION(*) :: Params
      COMPLEX , DIMENSION(*) :: Work
      REAL , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      REAL :: anorm , cwise_wrong , err_lbnd , illrcond_thresh ,        &
     &        rcond_tmp , rthresh , unstable_thresh
      LOGICAL :: colequ , ignore_cwise , notran , rowequ
      INTEGER :: ithresh , j , n_norms , prec_type , ref_type ,         &
     &           trans_type
      CHARACTER(1) :: norm
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
!     .. External Subroutines ..
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. External Functions ..
!     ..
!     .. Executable Statements ..
!
!     Check the input parameters.
!
      Info = 0
      trans_type = ILATRANS(Trans)
      ref_type = INT(ITREF_DEFAULT)
      IF ( Nparams>=LA_LINRX_ITREF_I ) THEN
         IF ( Params(LA_LINRX_ITREF_I)<0.0 ) THEN
            Params(LA_LINRX_ITREF_I) = ITREF_DEFAULT
         ELSE
            ref_type = Params(LA_LINRX_ITREF_I)
         ENDIF
      ENDIF
!
!     Set default parameters.
!
      illrcond_thresh = REAL(N)*SLAMCH('Epsilon')
      ithresh = INT(ITHRESH_DEFAULT)
      rthresh = RTHRESH_DEFAULT
      unstable_thresh = DZTHRESH_DEFAULT
      ignore_cwise = COMPONENTWISE_DEFAULT==0.0
!
      IF ( Nparams>=LA_LINRX_ITHRESH_I ) THEN
         IF ( Params(LA_LINRX_ITHRESH_I)<0.0 ) THEN
            Params(LA_LINRX_ITHRESH_I) = ithresh
         ELSE
            ithresh = INT(Params(LA_LINRX_ITHRESH_I))
         ENDIF
      ENDIF
      IF ( Nparams>=LA_LINRX_CWISE_I ) THEN
         IF ( Params(LA_LINRX_CWISE_I)>=0.0 ) THEN
            ignore_cwise = Params(LA_LINRX_CWISE_I)==0.0
         ELSEIF ( ignore_cwise ) THEN
            Params(LA_LINRX_CWISE_I) = 0.0
         ELSE
            Params(LA_LINRX_CWISE_I) = 1.0
         ENDIF
      ENDIF
      IF ( ref_type==0 .OR. N_err_bnds==0 ) THEN
         n_norms = 0
      ELSEIF ( ignore_cwise ) THEN
         n_norms = 1
      ELSE
         n_norms = 2
      ENDIF
!
      notran = LSAME(Trans,'N')
      rowequ = LSAME(Equed,'R') .OR. LSAME(Equed,'B')
      colequ = LSAME(Equed,'C') .OR. LSAME(Equed,'B')
!
!     Test input parameters.
!
      IF ( trans_type==-1 ) THEN
         Info = -1
      ELSEIF ( .NOT.rowequ .AND. .NOT.colequ .AND. .NOT.LSAME(Equed,'N')&
     &         ) THEN
         Info = -2
      ELSEIF ( N<0 ) THEN
         Info = -3
      ELSEIF ( Kl<0 ) THEN
         Info = -4
      ELSEIF ( Ku<0 ) THEN
         Info = -5
      ELSEIF ( Nrhs<0 ) THEN
         Info = -6
      ELSEIF ( Ldab<Kl+Ku+1 ) THEN
         Info = -8
      ELSEIF ( Ldafb<2*Kl+Ku+1 ) THEN
         Info = -10
      ELSEIF ( Ldb<MAX(1,N) ) THEN
         Info = -13
      ELSEIF ( Ldx<MAX(1,N) ) THEN
         Info = -15
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('CGBRFSX',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible.
!
      IF ( N==0 .OR. Nrhs==0 ) THEN
         Rcond = 1.0
         DO j = 1 , Nrhs
            Berr(j) = 0.0
            IF ( N_err_bnds>=1 ) THEN
               Err_bnds_norm(j,LA_LINRX_TRUST_I) = 1.0
               Err_bnds_comp(j,LA_LINRX_TRUST_I) = 1.0
            ENDIF
            IF ( N_err_bnds>=2 ) THEN
               Err_bnds_norm(j,LA_LINRX_ERR_I) = 0.0
               Err_bnds_comp(j,LA_LINRX_ERR_I) = 0.0
            ENDIF
            IF ( N_err_bnds>=3 ) THEN
               Err_bnds_norm(j,LA_LINRX_RCOND_I) = 1.0
               Err_bnds_comp(j,LA_LINRX_RCOND_I) = 1.0
            ENDIF
         ENDDO
         RETURN
      ENDIF
!
!     Default to failure.
!
      Rcond = 0.0
      DO j = 1 , Nrhs
         Berr(j) = 1.0
         IF ( N_err_bnds>=1 ) THEN
            Err_bnds_norm(j,LA_LINRX_TRUST_I) = 1.0
            Err_bnds_comp(j,LA_LINRX_TRUST_I) = 1.0
         ENDIF
         IF ( N_err_bnds>=2 ) THEN
            Err_bnds_norm(j,LA_LINRX_ERR_I) = 1.0
            Err_bnds_comp(j,LA_LINRX_ERR_I) = 1.0
         ENDIF
         IF ( N_err_bnds>=3 ) THEN
            Err_bnds_norm(j,LA_LINRX_RCOND_I) = 0.0
            Err_bnds_comp(j,LA_LINRX_RCOND_I) = 0.0
         ENDIF
      ENDDO
!
!     Compute the norm of A and the reciprocal of the condition
!     number of A.
!
      IF ( notran ) THEN
         norm = 'I'
      ELSE
         norm = '1'
      ENDIF
      anorm = CLANGB(norm,N,Kl,Ku,Ab,Ldab,Rwork)
      CALL CGBCON(norm,N,Kl,Ku,Afb,Ldafb,Ipiv,anorm,Rcond,Work,Rwork,   &
     &            Info)
!
!     Perform refinement on each right-hand side
!
      IF ( ref_type/=0 .AND. Info==0 ) THEN
 
         prec_type = ILAPREC('D')
 
         IF ( notran ) THEN
            CALL CLA_GBRFSX_EXTENDED(prec_type,trans_type,N,Kl,Ku,Nrhs, &
     &                               Ab,Ldab,Afb,Ldafb,Ipiv,colequ,C,B, &
     &                               Ldb,X,Ldx,Berr,n_norms,            &
     &                               Err_bnds_norm,Err_bnds_comp,Work,  &
     &                               Rwork,Work(N+1),                   &
     &                               TRANSFER(Rwork(1:2*N),             &
     &                               (/(ZERO,ZERO)/),N),Rcond,ithresh,  &
     &                               rthresh,unstable_thresh,           &
     &                               ignore_cwise,Info)
         ELSE
            CALL CLA_GBRFSX_EXTENDED(prec_type,trans_type,N,Kl,Ku,Nrhs, &
     &                               Ab,Ldab,Afb,Ldafb,Ipiv,rowequ,R,B, &
     &                               Ldb,X,Ldx,Berr,n_norms,            &
     &                               Err_bnds_norm,Err_bnds_comp,Work,  &
     &                               Rwork,Work(N+1),                   &
     &                               TRANSFER(Rwork(1:2*N),             &
     &                               (/(ZERO,ZERO)/),N),Rcond,ithresh,  &
     &                               rthresh,unstable_thresh,           &
     &                               ignore_cwise,Info)
         ENDIF
      ENDIF
 
      err_lbnd = MAX(10.0,SQRT(REAL(N)))*SLAMCH('Epsilon')
      IF ( N_err_bnds>=1 .AND. n_norms>=1 ) THEN
!
!     Compute scaled normwise condition number cond(A*C).
!
         IF ( colequ .AND. notran ) THEN
            rcond_tmp = CLA_GBRCOND_C(Trans,N,Kl,Ku,Ab,Ldab,Afb,Ldafb,  &
     &                  Ipiv,C,.TRUE.,Info,Work,Rwork)
         ELSEIF ( rowequ .AND. .NOT.notran ) THEN
            rcond_tmp = CLA_GBRCOND_C(Trans,N,Kl,Ku,Ab,Ldab,Afb,Ldafb,  &
     &                  Ipiv,R,.TRUE.,Info,Work,Rwork)
         ELSE
            rcond_tmp = CLA_GBRCOND_C(Trans,N,Kl,Ku,Ab,Ldab,Afb,Ldafb,  &
     &                  Ipiv,C,.FALSE.,Info,Work,Rwork)
         ENDIF
         DO j = 1 , Nrhs
!
!     Cap the error at 1.0.
!
            IF ( N_err_bnds>=LA_LINRX_ERR_I .AND.                       &
     &           Err_bnds_norm(j,LA_LINRX_ERR_I)>1.0 )                  &
     &           Err_bnds_norm(j,LA_LINRX_ERR_I) = 1.0
!
!     Threshold the error (see LAWN).
!
            IF ( rcond_tmp<illrcond_thresh ) THEN
               Err_bnds_norm(j,LA_LINRX_ERR_I) = 1.0
               Err_bnds_norm(j,LA_LINRX_TRUST_I) = 0.0
               IF ( Info<=N ) Info = N + j
            ELSEIF ( Err_bnds_norm(j,LA_LINRX_ERR_I)<err_lbnd ) THEN
               Err_bnds_norm(j,LA_LINRX_ERR_I) = err_lbnd
               Err_bnds_norm(j,LA_LINRX_TRUST_I) = 1.0
            ENDIF
!
!     Save the condition number.
!
            IF ( N_err_bnds>=LA_LINRX_RCOND_I )                         &
     &           Err_bnds_norm(j,LA_LINRX_RCOND_I) = rcond_tmp
 
         ENDDO
      ENDIF
 
      IF ( N_err_bnds>=1 .AND. n_norms>=2 ) THEN
!
!     Compute componentwise condition number cond(A*diag(Y(:,J))) for
!     each right-hand side using the current solution as an estimate of
!     the true solution.  If the componentwise error estimate is too
!     large, then the solution is a lousy estimate of truth and the
!     estimated RCOND may be too optimistic.  To avoid misleading users,
!     the inverse condition number is set to 0.0 when the estimated
!     cwise error is at least CWISE_WRONG.
!
         cwise_wrong = SQRT(SLAMCH('Epsilon'))
         DO j = 1 , Nrhs
            IF ( Err_bnds_comp(j,LA_LINRX_ERR_I)<cwise_wrong ) THEN
               rcond_tmp = CLA_GBRCOND_X(Trans,N,Kl,Ku,Ab,Ldab,Afb,     &
     &                     Ldafb,Ipiv,X(1,j),Info,Work,Rwork)
            ELSE
               rcond_tmp = 0.0
            ENDIF
!
!     Cap the error at 1.0.
!
            IF ( N_err_bnds>=LA_LINRX_ERR_I .AND.                       &
     &           Err_bnds_comp(j,LA_LINRX_ERR_I)>1.0 )                  &
     &           Err_bnds_comp(j,LA_LINRX_ERR_I) = 1.0
!
!     Threshold the error (see LAWN).
!
            IF ( rcond_tmp<illrcond_thresh ) THEN
               Err_bnds_comp(j,LA_LINRX_ERR_I) = 1.0
               Err_bnds_comp(j,LA_LINRX_TRUST_I) = 0.0
               IF ( Params(LA_LINRX_CWISE_I)==1.0 .AND. Info<N+j )      &
     &              Info = N + j
            ELSEIF ( Err_bnds_comp(j,LA_LINRX_ERR_I)<err_lbnd ) THEN
               Err_bnds_comp(j,LA_LINRX_ERR_I) = err_lbnd
               Err_bnds_comp(j,LA_LINRX_TRUST_I) = 1.0
            ENDIF
!
!     Save the condition number.
!
            IF ( N_err_bnds>=LA_LINRX_RCOND_I )                         &
     &           Err_bnds_comp(j,LA_LINRX_RCOND_I) = rcond_tmp
 
         ENDDO
      ENDIF
!
!
!     End of CGBRFSX
!
      END SUBROUTINE CGBRFSX
