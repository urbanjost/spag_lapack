!*==zgerfsx.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b ZGERFSX
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZGERFSX + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zgerfsx.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zgerfsx.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zgerfsx.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZGERFSX( TRANS, EQUED, N, NRHS, A, LDA, AF, LDAF, IPIV,
!                           R, C, B, LDB, X, LDX, RCOND, BERR, N_ERR_BNDS,
!                           ERR_BNDS_NORM, ERR_BNDS_COMP, NPARAMS, PARAMS,
!                           WORK, RWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          TRANS, EQUED
!       INTEGER            INFO, LDA, LDAF, LDB, LDX, N, NRHS, NPARAMS,
!      $                   N_ERR_BNDS
!       DOUBLE PRECISION   RCOND
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * )
!       COMPLEX*16         A( LDA, * ), AF( LDAF, * ), B( LDB, * ),
!      $                   X( LDX , * ), WORK( * )
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
!>    ZGERFSX improves the computed solution to a system of linear
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
!> \param[in] NRHS
!> \verbatim
!>          NRHS is INTEGER
!>     The number of right hand sides, i.e., the number of columns
!>     of the matrices B and X.  NRHS >= 0.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension (LDA,N)
!>     The original N-by-N matrix A.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>     The leading dimension of the array A.  LDA >= max(1,N).
!> \endverbatim
!>
!> \param[in] AF
!> \verbatim
!>          AF is COMPLEX*16 array, dimension (LDAF,N)
!>     The factors L and U from the factorization A = P*L*U
!>     as computed by ZGETRF.
!> \endverbatim
!>
!> \param[in] LDAF
!> \verbatim
!>          LDAF is INTEGER
!>     The leading dimension of the array AF.  LDAF >= max(1,N).
!> \endverbatim
!>
!> \param[in] IPIV
!> \verbatim
!>          IPIV is INTEGER array, dimension (N)
!>     The pivot indices from ZGETRF; for 1<=i<=N, row i of the
!>     matrix was interchanged with row IPIV(i).
!> \endverbatim
!>
!> \param[in] R
!> \verbatim
!>          R is DOUBLE PRECISION array, dimension (N)
!>     The row scale factors for A.  If EQUED = 'R' or 'B', A is
!>     multiplied on the left by diag(R); if EQUED = 'N' or 'C', R
!>     is not accessed.
!>     If R is accessed, each element of R should be a power of the radix
!>     to ensure a reliable solution and error estimates. Scaling by
!>     powers of the radix does not cause rounding errors unless the
!>     result underflows or overflows. Rounding errors during scaling
!>     lead to refining with a matrix that is not equivalent to the
!>     input matrix, producing error estimates that may not be
!>     reliable.
!> \endverbatim
!>
!> \param[in] C
!> \verbatim
!>          C is DOUBLE PRECISION array, dimension (N)
!>     The column scale factors for A.  If EQUED = 'C' or 'B', A is
!>     multiplied on the right by diag(C); if EQUED = 'N' or 'R', C
!>     is not accessed.
!>     If C is accessed, each element of C should be a power of the radix
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
!>          B is COMPLEX*16 array, dimension (LDB,NRHS)
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
!>          X is COMPLEX*16 array, dimension (LDX,NRHS)
!>     On entry, the solution matrix X, as computed by ZGETRS.
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
!> \date December 2016
!
!> \ingroup complex16GEcomputational
!
!  =====================================================================
      SUBROUTINE ZGERFSX(Trans,Equed,N,Nrhs,A,Lda,Af,Ldaf,Ipiv,R,C,B,   &
     &                   Ldb,X,Ldx,Rcond,Berr,N_err_bnds,Err_bnds_norm, &
     &                   Err_bnds_comp,Nparams,Params,Work,Rwork,Info)
      IMPLICIT NONE
!*--ZGERFSX417
!
!  -- LAPACK computational routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Trans , Equed
      INTEGER Info , Lda , Ldaf , Ldb , Ldx , N , Nrhs , Nparams ,      &
     &        N_err_bnds
      DOUBLE PRECISION Rcond
!     ..
!     .. Array Arguments ..
      INTEGER Ipiv(*)
      COMPLEX*16 A(Lda,*) , Af(Ldaf,*) , B(Ldb,*) , X(Ldx,*) , Work(*)
      DOUBLE PRECISION R(*) , C(*) , Params(*) , Berr(*) ,              &
     &                 Err_bnds_norm(Nrhs,*) , Err_bnds_comp(Nrhs,*) ,  &
     &                 Rwork(*)
!     ..
!
!  ==================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D+0)
      DOUBLE PRECISION ITREF_DEFAULT , ITHRESH_DEFAULT
      DOUBLE PRECISION COMPONENTWISE_DEFAULT , RTHRESH_DEFAULT
      DOUBLE PRECISION DZTHRESH_DEFAULT
      PARAMETER (ITREF_DEFAULT=1.0D+0)
      PARAMETER (ITHRESH_DEFAULT=10.0D+0)
      PARAMETER (COMPONENTWISE_DEFAULT=1.0D+0)
      PARAMETER (RTHRESH_DEFAULT=0.5D+0)
      PARAMETER (DZTHRESH_DEFAULT=0.25D+0)
      INTEGER LA_LINRX_ITREF_I , LA_LINRX_ITHRESH_I , LA_LINRX_CWISE_I
      PARAMETER (LA_LINRX_ITREF_I=1,LA_LINRX_ITHRESH_I=2)
      PARAMETER (LA_LINRX_CWISE_I=3)
      INTEGER LA_LINRX_TRUST_I , LA_LINRX_ERR_I , LA_LINRX_RCOND_I
      PARAMETER (LA_LINRX_TRUST_I=1,LA_LINRX_ERR_I=2)
      PARAMETER (LA_LINRX_RCOND_I=3)
!     ..
!     .. Local Scalars ..
      CHARACTER(1) norm
      LOGICAL rowequ , colequ , notran
      INTEGER j , trans_type , prec_type , ref_type
      INTEGER n_norms
      DOUBLE PRECISION anorm , rcond_tmp
      DOUBLE PRECISION illrcond_thresh , err_lbnd , cwise_wrong
      LOGICAL ignore_cwise
      INTEGER ithresh
      DOUBLE PRECISION rthresh , unstable_thresh
!     ..
!     .. External Subroutines ..
      EXTERNAL XERBLA , ZGECON , ZLA_GERFSX_EXTENDED
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX , SQRT , TRANSFER
!     ..
!     .. External Functions ..
      EXTERNAL LSAME , ILATRANS , ILAPREC
      EXTERNAL DLAMCH , ZLANGE , ZLA_GERCOND_X , ZLA_GERCOND_C
      DOUBLE PRECISION DLAMCH , ZLANGE , ZLA_GERCOND_X , ZLA_GERCOND_C
      LOGICAL LSAME
      INTEGER ILATRANS , ILAPREC
!     ..
!     .. Executable Statements ..
!
!     Check the input parameters.
!
      Info = 0
      trans_type = ILATRANS(Trans)
      ref_type = INT(ITREF_DEFAULT)
      IF ( Nparams>=LA_LINRX_ITREF_I ) THEN
         IF ( Params(LA_LINRX_ITREF_I)<0.0D+0 ) THEN
            Params(LA_LINRX_ITREF_I) = ITREF_DEFAULT
         ELSE
            ref_type = Params(LA_LINRX_ITREF_I)
         ENDIF
      ENDIF
!
!     Set default parameters.
!
      illrcond_thresh = DBLE(N)*DLAMCH('Epsilon')
      ithresh = INT(ITHRESH_DEFAULT)
      rthresh = RTHRESH_DEFAULT
      unstable_thresh = DZTHRESH_DEFAULT
      ignore_cwise = COMPONENTWISE_DEFAULT==0.0D+0
!
      IF ( Nparams>=LA_LINRX_ITHRESH_I ) THEN
         IF ( Params(LA_LINRX_ITHRESH_I)<0.0D+0 ) THEN
            Params(LA_LINRX_ITHRESH_I) = ithresh
         ELSE
            ithresh = INT(Params(LA_LINRX_ITHRESH_I))
         ENDIF
      ENDIF
      IF ( Nparams>=LA_LINRX_CWISE_I ) THEN
         IF ( Params(LA_LINRX_CWISE_I)>=0.0D+0 ) THEN
            ignore_cwise = Params(LA_LINRX_CWISE_I)==0.0D+0
         ELSEIF ( ignore_cwise ) THEN
            Params(LA_LINRX_CWISE_I) = 0.0D+0
         ELSE
            Params(LA_LINRX_CWISE_I) = 1.0D+0
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
      ELSEIF ( Nrhs<0 ) THEN
         Info = -4
      ELSEIF ( Lda<MAX(1,N) ) THEN
         Info = -6
      ELSEIF ( Ldaf<MAX(1,N) ) THEN
         Info = -8
      ELSEIF ( Ldb<MAX(1,N) ) THEN
         Info = -13
      ELSEIF ( Ldx<MAX(1,N) ) THEN
         Info = -15
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('ZGERFSX',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible.
!
      IF ( N==0 .OR. Nrhs==0 ) THEN
         Rcond = 1.0D+0
         DO j = 1 , Nrhs
            Berr(j) = 0.0D+0
            IF ( N_err_bnds>=1 ) THEN
               Err_bnds_norm(j,LA_LINRX_TRUST_I) = 1.0D+0
               Err_bnds_comp(j,LA_LINRX_TRUST_I) = 1.0D+0
            ENDIF
            IF ( N_err_bnds>=2 ) THEN
               Err_bnds_norm(j,LA_LINRX_ERR_I) = 0.0D+0
               Err_bnds_comp(j,LA_LINRX_ERR_I) = 0.0D+0
            ENDIF
            IF ( N_err_bnds>=3 ) THEN
               Err_bnds_norm(j,LA_LINRX_RCOND_I) = 1.0D+0
               Err_bnds_comp(j,LA_LINRX_RCOND_I) = 1.0D+0
            ENDIF
         ENDDO
         RETURN
      ENDIF
!
!     Default to failure.
!
      Rcond = 0.0D+0
      DO j = 1 , Nrhs
         Berr(j) = 1.0D+0
         IF ( N_err_bnds>=1 ) THEN
            Err_bnds_norm(j,LA_LINRX_TRUST_I) = 1.0D+0
            Err_bnds_comp(j,LA_LINRX_TRUST_I) = 1.0D+0
         ENDIF
         IF ( N_err_bnds>=2 ) THEN
            Err_bnds_norm(j,LA_LINRX_ERR_I) = 1.0D+0
            Err_bnds_comp(j,LA_LINRX_ERR_I) = 1.0D+0
         ENDIF
         IF ( N_err_bnds>=3 ) THEN
            Err_bnds_norm(j,LA_LINRX_RCOND_I) = 0.0D+0
            Err_bnds_comp(j,LA_LINRX_RCOND_I) = 0.0D+0
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
      anorm = ZLANGE(norm,N,N,A,Lda,Rwork)
      CALL ZGECON(norm,N,Af,Ldaf,anorm,Rcond,Work,Rwork,Info)
!
!     Perform refinement on each right-hand side
!
      IF ( ref_type/=0 ) THEN
 
         prec_type = ILAPREC('E')
 
         IF ( notran ) THEN
            CALL ZLA_GERFSX_EXTENDED(prec_type,trans_type,N,Nrhs,A,Lda, &
     &                               Af,Ldaf,Ipiv,colequ,C,B,Ldb,X,Ldx, &
     &                               Berr,n_norms,Err_bnds_norm,        &
     &                               Err_bnds_comp,Work,Rwork,Work(N+1),&
     &                               TRANSFER(Rwork(1:2*N),             &
     &                               (/(ZERO,ZERO)/),N),Rcond,ithresh,  &
     &                               rthresh,unstable_thresh,           &
     &                               ignore_cwise,Info)
         ELSE
            CALL ZLA_GERFSX_EXTENDED(prec_type,trans_type,N,Nrhs,A,Lda, &
     &                               Af,Ldaf,Ipiv,rowequ,R,B,Ldb,X,Ldx, &
     &                               Berr,n_norms,Err_bnds_norm,        &
     &                               Err_bnds_comp,Work,Rwork,Work(N+1),&
     &                               TRANSFER(Rwork(1:2*N),             &
     &                               (/(ZERO,ZERO)/),N),Rcond,ithresh,  &
     &                               rthresh,unstable_thresh,           &
     &                               ignore_cwise,Info)
         ENDIF
      ENDIF
 
      err_lbnd = MAX(10.0D+0,SQRT(DBLE(N)))*DLAMCH('Epsilon')
      IF ( N_err_bnds>=1 .AND. n_norms>=1 ) THEN
!
!     Compute scaled normwise condition number cond(A*C).
!
         IF ( colequ .AND. notran ) THEN
            rcond_tmp = ZLA_GERCOND_C(Trans,N,A,Lda,Af,Ldaf,Ipiv,C,     &
     &                  .TRUE.,Info,Work,Rwork)
         ELSEIF ( rowequ .AND. .NOT.notran ) THEN
            rcond_tmp = ZLA_GERCOND_C(Trans,N,A,Lda,Af,Ldaf,Ipiv,R,     &
     &                  .TRUE.,Info,Work,Rwork)
         ELSE
            rcond_tmp = ZLA_GERCOND_C(Trans,N,A,Lda,Af,Ldaf,Ipiv,C,     &
     &                  .FALSE.,Info,Work,Rwork)
         ENDIF
         DO j = 1 , Nrhs
!
!     Cap the error at 1.0.
!
            IF ( N_err_bnds>=LA_LINRX_ERR_I .AND.                       &
     &           Err_bnds_norm(j,LA_LINRX_ERR_I)>1.0D+0 )               &
     &           Err_bnds_norm(j,LA_LINRX_ERR_I) = 1.0D+0
!
!     Threshold the error (see LAWN).
!
            IF ( rcond_tmp<illrcond_thresh ) THEN
               Err_bnds_norm(j,LA_LINRX_ERR_I) = 1.0D+0
               Err_bnds_norm(j,LA_LINRX_TRUST_I) = 0.0D+0
               IF ( Info<=N ) Info = N + j
            ELSEIF ( Err_bnds_norm(j,LA_LINRX_ERR_I)<err_lbnd ) THEN
               Err_bnds_norm(j,LA_LINRX_ERR_I) = err_lbnd
               Err_bnds_norm(j,LA_LINRX_TRUST_I) = 1.0D+0
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
         cwise_wrong = SQRT(DLAMCH('Epsilon'))
         DO j = 1 , Nrhs
            IF ( Err_bnds_comp(j,LA_LINRX_ERR_I)<cwise_wrong ) THEN
               rcond_tmp = ZLA_GERCOND_X(Trans,N,A,Lda,Af,Ldaf,Ipiv,    &
     &                     X(1,j),Info,Work,Rwork)
            ELSE
               rcond_tmp = 0.0D+0
            ENDIF
!
!     Cap the error at 1.0.
!
            IF ( N_err_bnds>=LA_LINRX_ERR_I .AND.                       &
     &           Err_bnds_comp(j,LA_LINRX_ERR_I)>1.0D+0 )               &
     &           Err_bnds_comp(j,LA_LINRX_ERR_I) = 1.0D+0
!
!     Threshold the error (see LAWN).
!
            IF ( rcond_tmp<illrcond_thresh ) THEN
               Err_bnds_comp(j,LA_LINRX_ERR_I) = 1.0D+0
               Err_bnds_comp(j,LA_LINRX_TRUST_I) = 0.0D+0
               IF ( Params(LA_LINRX_CWISE_I)==1.0D+0 .AND. Info<N+j )   &
     &              Info = N + j
            ELSEIF ( Err_bnds_comp(j,LA_LINRX_ERR_I)<err_lbnd ) THEN
               Err_bnds_comp(j,LA_LINRX_ERR_I) = err_lbnd
               Err_bnds_comp(j,LA_LINRX_TRUST_I) = 1.0D+0
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
!     End of ZGERFSX
!
      END SUBROUTINE ZGERFSX
