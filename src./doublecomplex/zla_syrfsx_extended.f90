!*==zla_syrfsx_extended.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b ZLA_SYRFSX_EXTENDED improves the computed solution to a system of linear equations for symmetric indefinite matrices by performing extra-precise iterative refinement and provides error bounds and backward error estimates for the solution.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZLA_SYRFSX_EXTENDED + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zla_syrfsx_extended.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zla_syrfsx_extended.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zla_syrfsx_extended.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZLA_SYRFSX_EXTENDED( PREC_TYPE, UPLO, N, NRHS, A, LDA,
!                                       AF, LDAF, IPIV, COLEQU, C, B, LDB,
!                                       Y, LDY, BERR_OUT, N_NORMS,
!                                       ERR_BNDS_NORM, ERR_BNDS_COMP, RES,
!                                       AYB, DY, Y_TAIL, RCOND, ITHRESH,
!                                       RTHRESH, DZ_UB, IGNORE_CWISE,
!                                       INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, LDA, LDAF, LDB, LDY, N, NRHS, PREC_TYPE,
!      $                   N_NORMS, ITHRESH
!       CHARACTER          UPLO
!       LOGICAL            COLEQU, IGNORE_CWISE
!       DOUBLE PRECISION   RTHRESH, DZ_UB
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * )
!       COMPLEX*16         A( LDA, * ), AF( LDAF, * ), B( LDB, * ),
!      $                   Y( LDY, * ), RES( * ), DY( * ), Y_TAIL( * )
!       DOUBLE PRECISION   C( * ), AYB( * ), RCOND, BERR_OUT( * ),
!      $                   ERR_BNDS_NORM( NRHS, * ),
!      $                   ERR_BNDS_COMP( NRHS, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZLA_SYRFSX_EXTENDED improves the computed solution to a system of
!> linear equations by performing extra-precise iterative refinement
!> and provides error bounds and backward error estimates for the solution.
!> This subroutine is called by ZSYRFSX to perform iterative refinement.
!> In addition to normwise error bound, the code provides maximum
!> componentwise error bound if possible. See comments for ERR_BNDS_NORM
!> and ERR_BNDS_COMP for details of the error bounds. Note that this
!> subroutine is only resonsible for setting the second fields of
!> ERR_BNDS_NORM and ERR_BNDS_COMP.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] PREC_TYPE
!> \verbatim
!>          PREC_TYPE is INTEGER
!>     Specifies the intermediate precision to be used in refinement.
!>     The value is defined by ILAPREC(P) where P is a CHARACTER and P
!>          = 'S':  Single
!>          = 'D':  Double
!>          = 'I':  Indigenous
!>          = 'X' or 'E':  Extra
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
!>     The number of right-hand-sides, i.e., the number of columns of the
!>     matrix B.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension (LDA,N)
!>     On entry, the N-by-N matrix A.
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
!>     The block diagonal matrix D and the multipliers used to
!>     obtain the factor U or L as computed by ZSYTRF.
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
!>     Details of the interchanges and the block structure of D
!>     as determined by ZSYTRF.
!> \endverbatim
!>
!> \param[in] COLEQU
!> \verbatim
!>          COLEQU is LOGICAL
!>     If .TRUE. then column equilibration was done to A before calling
!>     this routine. This is needed to compute the solution and error
!>     bounds correctly.
!> \endverbatim
!>
!> \param[in] C
!> \verbatim
!>          C is DOUBLE PRECISION array, dimension (N)
!>     The column scale factors for A. If COLEQU = .FALSE., C
!>     is not accessed. If C is input, each element of C should be a power
!>     of the radix to ensure a reliable solution and error estimates.
!>     Scaling by powers of the radix does not cause rounding errors unless
!>     the result underflows or overflows. Rounding errors during scaling
!>     lead to refining with a matrix that is not equivalent to the
!>     input matrix, producing error estimates that may not be
!>     reliable.
!> \endverbatim
!>
!> \param[in] B
!> \verbatim
!>          B is COMPLEX*16 array, dimension (LDB,NRHS)
!>     The right-hand-side matrix B.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>     The leading dimension of the array B.  LDB >= max(1,N).
!> \endverbatim
!>
!> \param[in,out] Y
!> \verbatim
!>          Y is COMPLEX*16 array, dimension (LDY,NRHS)
!>     On entry, the solution matrix X, as computed by ZSYTRS.
!>     On exit, the improved solution matrix Y.
!> \endverbatim
!>
!> \param[in] LDY
!> \verbatim
!>          LDY is INTEGER
!>     The leading dimension of the array Y.  LDY >= max(1,N).
!> \endverbatim
!>
!> \param[out] BERR_OUT
!> \verbatim
!>          BERR_OUT is DOUBLE PRECISION array, dimension (NRHS)
!>     On exit, BERR_OUT(j) contains the componentwise relative backward
!>     error for right-hand-side j from the formula
!>         max(i) ( abs(RES(i)) / ( abs(op(A_s))*abs(Y) + abs(B_s) )(i) )
!>     where abs(Z) is the componentwise absolute value of the matrix
!>     or vector Z. This is computed by ZLA_LIN_BERR.
!> \endverbatim
!>
!> \param[in] N_NORMS
!> \verbatim
!>          N_NORMS is INTEGER
!>     Determines which error bounds to return (see ERR_BNDS_NORM
!>     and ERR_BNDS_COMP).
!>     If N_NORMS >= 1 return normwise error bounds.
!>     If N_NORMS >= 2 return componentwise error bounds.
!> \endverbatim
!>
!> \param[in,out] ERR_BNDS_NORM
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
!>     This subroutine is only responsible for setting the second field
!>     above.
!>     See Lapack Working Note 165 for further details and extra
!>     cautions.
!> \endverbatim
!>
!> \param[in,out] ERR_BNDS_COMP
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
!>     This subroutine is only responsible for setting the second field
!>     above.
!>     See Lapack Working Note 165 for further details and extra
!>     cautions.
!> \endverbatim
!>
!> \param[in] RES
!> \verbatim
!>          RES is COMPLEX*16 array, dimension (N)
!>     Workspace to hold the intermediate residual.
!> \endverbatim
!>
!> \param[in] AYB
!> \verbatim
!>          AYB is DOUBLE PRECISION array, dimension (N)
!>     Workspace.
!> \endverbatim
!>
!> \param[in] DY
!> \verbatim
!>          DY is COMPLEX*16 array, dimension (N)
!>     Workspace to hold the intermediate solution.
!> \endverbatim
!>
!> \param[in] Y_TAIL
!> \verbatim
!>          Y_TAIL is COMPLEX*16 array, dimension (N)
!>     Workspace to hold the trailing bits of the intermediate solution.
!> \endverbatim
!>
!> \param[in] RCOND
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
!> \param[in] ITHRESH
!> \verbatim
!>          ITHRESH is INTEGER
!>     The maximum number of residual computations allowed for
!>     refinement. The default is 10. For 'aggressive' set to 100 to
!>     permit convergence using approximate factorizations or
!>     factorizations other than LU. If the factorization uses a
!>     technique other than Gaussian elimination, the guarantees in
!>     ERR_BNDS_NORM and ERR_BNDS_COMP may no longer be trustworthy.
!> \endverbatim
!>
!> \param[in] RTHRESH
!> \verbatim
!>          RTHRESH is DOUBLE PRECISION
!>     Determines when to stop refinement if the error estimate stops
!>     decreasing. Refinement will stop when the next solution no longer
!>     satisfies norm(dx_{i+1}) < RTHRESH * norm(dx_i) where norm(Z) is
!>     the infinity norm of Z. RTHRESH satisfies 0 < RTHRESH <= 1. The
!>     default value is 0.5. For 'aggressive' set to 0.9 to permit
!>     convergence on extremely ill-conditioned matrices. See LAWN 165
!>     for more details.
!> \endverbatim
!>
!> \param[in] DZ_UB
!> \verbatim
!>          DZ_UB is DOUBLE PRECISION
!>     Determines when to start considering componentwise convergence.
!>     Componentwise convergence is only considered after each component
!>     of the solution Y is stable, which we define as the relative
!>     change in each component being less than DZ_UB. The default value
!>     is 0.25, requiring the first bit to be stable. See LAWN 165 for
!>     more details.
!> \endverbatim
!>
!> \param[in] IGNORE_CWISE
!> \verbatim
!>          IGNORE_CWISE is LOGICAL
!>     If .TRUE. then ignore componentwise convergence. Default value
!>     is .FALSE..
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>       = 0:  Successful exit.
!>       < 0:  if INFO = -i, the ith argument to ZLA_HERFSX_EXTENDED had an illegal
!>             value
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
!> \date June 2017
!
!> \ingroup complex16SYcomputational
!
!  =====================================================================
      SUBROUTINE ZLA_SYRFSX_EXTENDED(Prec_type,Uplo,N,Nrhs,A,Lda,Af,    &
     &                               Ldaf,Ipiv,Colequ,C,B,Ldb,Y,Ldy,    &
     &                               Berr_out,N_norms,Err_bnds_norm,    &
     &                               Err_bnds_comp,Res,Ayb,Dy,Y_tail,   &
     &                               Rcond,Ithresh,Rthresh,Dz_ub,       &
     &                               Ignore_cwise,Info)
      USE F77KINDS                        
      USE S_BLAS_ZSYMV2_X
      USE S_BLAS_ZSYMV_X
      USE S_DLAMCH
      USE S_ILAUPLO
      USE S_LSAME
      USE S_XERBLA
      USE S_ZAXPY
      USE S_ZCOPY
      USE S_ZLA_LIN_BERR
      USE S_ZLA_SYAMV
      USE S_ZLA_WWADDW
      USE S_ZSYMV
      USE S_ZSYTRS
      IMPLICIT NONE
!*--ZLA_SYRFSX_EXTENDED412
!
! PARAMETER definitions rewritten by SPAG
!
      INTEGER , PARAMETER  ::  UNSTABLE_STATE = 0 , WORKING_STATE = 1 , &
     &                         CONV_STATE = 2 , NOPROG_STATE = 3 ,      &
     &                         BASE_RESIDUAL = 0 , EXTRA_RESIDUAL = 1 , &
     &                         EXTRA_Y = 2 , LA_LINRX_ERR_I = 2
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER :: Prec_type
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER , INTENT(IN) :: Nrhs
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(Ldaf,*) :: Af
      INTEGER :: Ldaf
      INTEGER , DIMENSION(*) :: Ipiv
      LOGICAL , INTENT(IN) :: Colequ
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: C
      COMPLEX(CX16KIND) , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      COMPLEX(CX16KIND) , DIMENSION(Ldy,*) :: Y
      INTEGER , INTENT(IN) :: Ldy
      REAL(R8KIND) , DIMENSION(*) :: Berr_out
      INTEGER , INTENT(IN) :: N_norms
      REAL(R8KIND) , INTENT(OUT) , DIMENSION(Nrhs,*) :: Err_bnds_norm
      REAL(R8KIND) , INTENT(OUT) , DIMENSION(Nrhs,*) :: Err_bnds_comp
      COMPLEX(CX16KIND) , DIMENSION(*) :: Res
      REAL(R8KIND) , DIMENSION(*) :: Ayb
      COMPLEX(CX16KIND) , DIMENSION(*) :: Dy
      COMPLEX(CX16KIND) , DIMENSION(*) :: Y_tail
      REAL(R8KIND) , INTENT(IN) :: Rcond
      INTEGER , INTENT(IN) :: Ithresh
      REAL(R8KIND) , INTENT(IN) :: Rthresh
      REAL(R8KIND) , INTENT(IN) :: Dz_ub
      LOGICAL , INTENT(IN) :: Ignore_cwise
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      REAL(R8KIND) :: CABS1
      INTEGER :: cnt , i , j , uplo2 , x_state , y_prec_state , z_state
      REAL(R8KIND) :: dxrat , dxratmax , dx_x , dyk , dzrat , dzratmax ,&
     &                dz_z , eps , final_dx_x , final_dz_z , hugeval ,  &
     &                incr_thresh , normdx , normx , normy ,            &
     &                prevnormdx , prev_dz_z , yk , ymin
      LOGICAL :: incr_prec , upper
      COMPLEX(CX16KIND) :: zdum
!
! End of declarations rewritten by SPAG
!
!     ..
!     .. Array Arguments ..
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
!     ..
!     .. Parameters ..
!     ..
!     .. External Functions ..
!     ..
!     .. External Subroutines ..
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. Statement Functions ..
!     ..
!     .. Statement Function Definitions ..
      CABS1(zdum) = ABS(DBLE(zdum)) + ABS(DIMAG(zdum))
!     ..
!     .. Executable Statements ..
!
      Info = 0
      upper = LSAME(Uplo,'U')
      IF ( .NOT.upper .AND. .NOT.LSAME(Uplo,'L') ) THEN
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
      ELSEIF ( Ldy<MAX(1,N) ) THEN
         Info = -15
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('ZLA_HERFSX_EXTENDED',-Info)
         RETURN
      ENDIF
      eps = DLAMCH('Epsilon')
      hugeval = DLAMCH('Overflow')
!     Force HUGEVAL to Inf
      hugeval = hugeval*hugeval
!     Using HUGEVAL may lead to spurious underflows.
      incr_thresh = DBLE(N)*eps
 
      IF ( LSAME(Uplo,'L') ) THEN
         uplo2 = ILAUPLO('L')
      ELSE
         uplo2 = ILAUPLO('U')
      ENDIF
 
      DO j = 1 , Nrhs
         y_prec_state = EXTRA_RESIDUAL
         IF ( y_prec_state==EXTRA_Y ) THEN
            DO i = 1 , N
               Y_tail(i) = 0.0D+0
            ENDDO
         ENDIF
 
         dxrat = 0.0D+0
         dxratmax = 0.0D+0
         dzrat = 0.0D+0
         dzratmax = 0.0D+0
         final_dx_x = hugeval
         final_dz_z = hugeval
         prevnormdx = hugeval
         prev_dz_z = hugeval
         dz_z = hugeval
         dx_x = hugeval
 
         x_state = WORKING_STATE
         z_state = UNSTABLE_STATE
         incr_prec = .FALSE.
 
         DO cnt = 1 , Ithresh
!
!         Compute residual RES = B_s - op(A_s) * Y,
!             op(A) = A, A**T, or A**H depending on TRANS (and type).
!
            CALL ZCOPY(N,B(1,j),1,Res,1)
            IF ( y_prec_state==BASE_RESIDUAL ) THEN
               CALL ZSYMV(Uplo,N,DCMPLX(-1.0D+0),A,Lda,Y(1,j),1,        &
     &                    DCMPLX(1.0D+0),Res,1)
            ELSEIF ( y_prec_state==EXTRA_RESIDUAL ) THEN
               CALL BLAS_ZSYMV_X(uplo2,N,DCMPLX(-1.0D+0),A,Lda,Y(1,j),1,&
     &                           DCMPLX(1.0D+0),Res,1,Prec_type)
            ELSE
               CALL BLAS_ZSYMV2_X(uplo2,N,DCMPLX(-1.0D+0),A,Lda,Y(1,j), &
     &                            Y_tail,1,DCMPLX(1.0D+0),Res,1,        &
     &                            Prec_type)
            ENDIF
 
!         XXX: RES is no longer needed.
            CALL ZCOPY(N,Res,1,Dy,1)
            CALL ZSYTRS(Uplo,N,1,Af,Ldaf,Ipiv,Dy,N,Info)
!
!         Calculate relative changes DX_X, DZ_Z and ratios DXRAT, DZRAT.
!
            normx = 0.0D+0
            normy = 0.0D+0
            normdx = 0.0D+0
            dz_z = 0.0D+0
            ymin = hugeval
 
            DO i = 1 , N
               yk = CABS1(Y(i,j))
               dyk = CABS1(Dy(i))
 
               IF ( yk/=0.0D+0 ) THEN
                  dz_z = MAX(dz_z,dyk/yk)
               ELSEIF ( dyk/=0.0D+0 ) THEN
                  dz_z = hugeval
               ENDIF
 
               ymin = MIN(ymin,yk)
 
               normy = MAX(normy,yk)
 
               IF ( Colequ ) THEN
                  normx = MAX(normx,yk*C(i))
                  normdx = MAX(normdx,dyk*C(i))
               ELSE
                  normx = normy
                  normdx = MAX(normdx,dyk)
               ENDIF
            ENDDO
 
            IF ( normx/=0.0D+0 ) THEN
               dx_x = normdx/normx
            ELSEIF ( normdx==0.0D+0 ) THEN
               dx_x = 0.0D+0
            ELSE
               dx_x = hugeval
            ENDIF
 
            dxrat = normdx/prevnormdx
            dzrat = dz_z/prev_dz_z
!
!         Check termination criteria.
!
            IF ( ymin*Rcond<incr_thresh*normy .AND.                     &
     &           y_prec_state<EXTRA_Y ) incr_prec = .TRUE.
 
            IF ( x_state==NOPROG_STATE .AND. dxrat<=Rthresh )           &
     &           x_state = WORKING_STATE
            IF ( x_state==WORKING_STATE ) THEN
               IF ( dx_x<=eps ) THEN
                  x_state = CONV_STATE
               ELSEIF ( dxrat>Rthresh ) THEN
                  IF ( y_prec_state/=EXTRA_Y ) THEN
                     incr_prec = .TRUE.
                  ELSE
                     x_state = NOPROG_STATE
                  ENDIF
               ELSE
                  IF ( dxrat>dxratmax ) dxratmax = dxrat
               ENDIF
               IF ( x_state>WORKING_STATE ) final_dx_x = dx_x
            ENDIF
 
            IF ( z_state==UNSTABLE_STATE .AND. dz_z<=Dz_ub )            &
     &           z_state = WORKING_STATE
            IF ( z_state==NOPROG_STATE .AND. dzrat<=Rthresh )           &
     &           z_state = WORKING_STATE
            IF ( z_state==WORKING_STATE ) THEN
               IF ( dz_z<=eps ) THEN
                  z_state = CONV_STATE
               ELSEIF ( dz_z>Dz_ub ) THEN
                  z_state = UNSTABLE_STATE
                  dzratmax = 0.0D+0
                  final_dz_z = hugeval
               ELSEIF ( dzrat>Rthresh ) THEN
                  IF ( y_prec_state/=EXTRA_Y ) THEN
                     incr_prec = .TRUE.
                  ELSE
                     z_state = NOPROG_STATE
                  ENDIF
               ELSE
                  IF ( dzrat>dzratmax ) dzratmax = dzrat
               ENDIF
               IF ( z_state>WORKING_STATE ) final_dz_z = dz_z
            ENDIF
 
            IF ( x_state/=WORKING_STATE .AND.                           &
     &           (Ignore_cwise .OR. z_state/=WORKING_STATE) ) EXIT
 
            IF ( incr_prec ) THEN
               incr_prec = .FALSE.
               y_prec_state = y_prec_state + 1
               DO i = 1 , N
                  Y_tail(i) = 0.0D+0
               ENDDO
            ENDIF
 
            prevnormdx = normdx
            prev_dz_z = dz_z
!
!           Update soluton.
!
            IF ( y_prec_state<EXTRA_Y ) THEN
               CALL ZAXPY(N,DCMPLX(1.0D+0),Dy,1,Y(1,j),1)
            ELSE
               CALL ZLA_WWADDW(N,Y(1,j),Y_tail,Dy)
            ENDIF
 
         ENDDO
!        Target of "IF (Z_STOP .AND. X_STOP)".  Sun's f77 won't EXIT.
!
!     Set final_* when cnt hits ithresh.
!
         IF ( x_state==WORKING_STATE ) final_dx_x = dx_x
         IF ( z_state==WORKING_STATE ) final_dz_z = dz_z
!
!     Compute error bounds.
!
         IF ( N_norms>=1 ) Err_bnds_norm(j,LA_LINRX_ERR_I)              &
     &        = final_dx_x/(1-dxratmax)
         IF ( N_norms>=2 ) Err_bnds_comp(j,LA_LINRX_ERR_I)              &
     &        = final_dz_z/(1-dzratmax)
!
!     Compute componentwise relative backward error from formula
!         max(i) ( abs(R(i)) / ( abs(op(A_s))*abs(Y) + abs(B_s) )(i) )
!     where abs(Z) is the componentwise absolute value of the matrix
!     or vector Z.
!
!        Compute residual RES = B_s - op(A_s) * Y,
!            op(A) = A, A**T, or A**H depending on TRANS (and type).
!
         CALL ZCOPY(N,B(1,j),1,Res,1)
         CALL ZSYMV(Uplo,N,DCMPLX(-1.0D+0),A,Lda,Y(1,j),1,DCMPLX(1.0D+0)&
     &              ,Res,1)
 
         DO i = 1 , N
            Ayb(i) = CABS1(B(i,j))
         ENDDO
!
!     Compute abs(op(A_s))*abs(Y) + abs(B_s).
!
         CALL ZLA_SYAMV(uplo2,N,1.0D+0,A,Lda,Y(1,j),1,1.0D+0,Ayb,1)
 
         CALL ZLA_LIN_BERR(N,N,1,Res,Ayb,Berr_out(j))
!
!     End of loop for each RHS.
!
      ENDDO
!
      END SUBROUTINE ZLA_SYRFSX_EXTENDED
