!*==cla_gerfsx_extended.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b CLA_GERFSX_EXTENDED
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CLA_GERFSX_EXTENDED + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cla_gerfsx_extended.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cla_gerfsx_extended.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cla_gerfsx_extended.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CLA_GERFSX_EXTENDED( PREC_TYPE, TRANS_TYPE, N, NRHS, A,
!                                       LDA, AF, LDAF, IPIV, COLEQU, C, B,
!                                       LDB, Y, LDY, BERR_OUT, N_NORMS,
!                                       ERRS_N, ERRS_C, RES, AYB, DY,
!                                       Y_TAIL, RCOND, ITHRESH, RTHRESH,
!                                       DZ_UB, IGNORE_CWISE, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, LDA, LDAF, LDB, LDY, N, NRHS, PREC_TYPE,
!      $                   TRANS_TYPE, N_NORMS
!       LOGICAL            COLEQU, IGNORE_CWISE
!       INTEGER            ITHRESH
!       REAL               RTHRESH, DZ_UB
!       ..
!       .. Array Arguments
!       INTEGER            IPIV( * )
!       COMPLEX            A( LDA, * ), AF( LDAF, * ), B( LDB, * ),
!      $                   Y( LDY, * ), RES( * ), DY( * ), Y_TAIL( * )
!       REAL               C( * ), AYB( * ), RCOND, BERR_OUT( * ),
!      $                   ERRS_N( NRHS, * ), ERRS_C( NRHS, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>
!> CLA_GERFSX_EXTENDED improves the computed solution to a system of
!> linear equations by performing extra-precise iterative refinement
!> and provides error bounds and backward error estimates for the solution.
!> This subroutine is called by CGERFSX to perform iterative refinement.
!> In addition to normwise error bound, the code provides maximum
!> componentwise error bound if possible. See comments for ERRS_N
!> and ERRS_C for details of the error bounds. Note that this
!> subroutine is only resonsible for setting the second fields of
!> ERRS_N and ERRS_C.
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
!> \param[in] TRANS_TYPE
!> \verbatim
!>          TRANS_TYPE is INTEGER
!>     Specifies the transposition operation on A.
!>     The value is defined by ILATRANS(T) where T is a CHARACTER and T
!>          = 'N':  No transpose
!>          = 'T':  Transpose
!>          = 'C':  Conjugate transpose
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
!>          A is COMPLEX array, dimension (LDA,N)
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
!>          AF is COMPLEX array, dimension (LDAF,N)
!>     The factors L and U from the factorization
!>     A = P*L*U as computed by CGETRF.
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
!>     The pivot indices from the factorization A = P*L*U
!>     as computed by CGETRF; row i of the matrix was interchanged
!>     with row IPIV(i).
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
!>          C is REAL array, dimension (N)
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
!>          B is COMPLEX array, dimension (LDB,NRHS)
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
!>          Y is COMPLEX array, dimension (LDY,NRHS)
!>     On entry, the solution matrix X, as computed by CGETRS.
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
!>          BERR_OUT is REAL array, dimension (NRHS)
!>     On exit, BERR_OUT(j) contains the componentwise relative backward
!>     error for right-hand-side j from the formula
!>         max(i) ( abs(RES(i)) / ( abs(op(A_s))*abs(Y) + abs(B_s) )(i) )
!>     where abs(Z) is the componentwise absolute value of the matrix
!>     or vector Z. This is computed by CLA_LIN_BERR.
!> \endverbatim
!>
!> \param[in] N_NORMS
!> \verbatim
!>          N_NORMS is INTEGER
!>     Determines which error bounds to return (see ERRS_N
!>     and ERRS_C).
!>     If N_NORMS >= 1 return normwise error bounds.
!>     If N_NORMS >= 2 return componentwise error bounds.
!> \endverbatim
!>
!> \param[in,out] ERRS_N
!> \verbatim
!>          ERRS_N is REAL array, dimension (NRHS, N_ERR_BNDS)
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
!>     The first index in ERRS_N(i,:) corresponds to the ith
!>     right-hand side.
!>
!>     The second index in ERRS_N(:,err) contains the following
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
!> \param[in,out] ERRS_C
!> \verbatim
!>          ERRS_C is REAL array, dimension (NRHS, N_ERR_BNDS)
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
!>     ERRS_C is not accessed.  If N_ERR_BNDS < 3, then at most
!>     the first (:,N_ERR_BNDS) entries are returned.
!>
!>     The first index in ERRS_C(i,:) corresponds to the ith
!>     right-hand side.
!>
!>     The second index in ERRS_C(:,err) contains the following
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
!>          RES is COMPLEX array, dimension (N)
!>     Workspace to hold the intermediate residual.
!> \endverbatim
!>
!> \param[in] AYB
!> \verbatim
!>          AYB is REAL array, dimension (N)
!>     Workspace.
!> \endverbatim
!>
!> \param[in] DY
!> \verbatim
!>          DY is COMPLEX array, dimension (N)
!>     Workspace to hold the intermediate solution.
!> \endverbatim
!>
!> \param[in] Y_TAIL
!> \verbatim
!>          Y_TAIL is COMPLEX array, dimension (N)
!>     Workspace to hold the trailing bits of the intermediate solution.
!> \endverbatim
!>
!> \param[in] RCOND
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
!> \param[in] ITHRESH
!> \verbatim
!>          ITHRESH is INTEGER
!>     The maximum number of residual computations allowed for
!>     refinement. The default is 10. For 'aggressive' set to 100 to
!>     permit convergence using approximate factorizations or
!>     factorizations other than LU. If the factorization uses a
!>     technique other than Gaussian elimination, the guarantees in
!>     ERRS_N and ERRS_C may no longer be trustworthy.
!> \endverbatim
!>
!> \param[in] RTHRESH
!> \verbatim
!>          RTHRESH is REAL
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
!>          DZ_UB is REAL
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
!>       < 0:  if INFO = -i, the ith argument to CGETRS had an illegal
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
!> \date December 2016
!
!> \ingroup complexGEcomputational
!
!  =====================================================================
      SUBROUTINE CLA_GERFSX_EXTENDED(Prec_type,Trans_type,N,Nrhs,A,Lda, &
     &                               Af,Ldaf,Ipiv,Colequ,C,B,Ldb,Y,Ldy, &
     &                               Berr_out,N_norms,Errs_n,Errs_c,Res,&
     &                               Ayb,Dy,Y_tail,Rcond,Ithresh,       &
     &                               Rthresh,Dz_ub,Ignore_cwise,Info)
      USE S_BLAS_CGEMV2_X
      USE S_BLAS_CGEMV_X
      USE S_CAXPY
      USE S_CCOPY
      USE S_CGEMV
      USE S_CGETRS
      USE S_CHLA_TRANSTYPE
      USE S_CLA_GEAMV
      USE S_CLA_LIN_BERR
      USE S_CLA_WWADDW
      USE S_SLAMCH
      IMPLICIT NONE
!*--CLA_GERFSX_EXTENDED411
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
      INTEGER :: Trans_type
      INTEGER :: N
      INTEGER , INTENT(IN) :: Nrhs
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(Ldaf,*) :: Af
      INTEGER :: Ldaf
      INTEGER , DIMENSION(*) :: Ipiv
      LOGICAL , INTENT(IN) :: Colequ
      REAL , INTENT(IN) , DIMENSION(*) :: C
      COMPLEX , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      COMPLEX , DIMENSION(Ldy,*) :: Y
      INTEGER , INTENT(IN) :: Ldy
      REAL , DIMENSION(*) :: Berr_out
      INTEGER , INTENT(IN) :: N_norms
      REAL , INTENT(OUT) , DIMENSION(Nrhs,*) :: Errs_n
      REAL , INTENT(OUT) , DIMENSION(Nrhs,*) :: Errs_c
      COMPLEX , DIMENSION(*) :: Res
      REAL , DIMENSION(*) :: Ayb
      COMPLEX , DIMENSION(*) :: Dy
      COMPLEX , DIMENSION(*) :: Y_tail
      REAL , INTENT(IN) :: Rcond
      INTEGER , INTENT(IN) :: Ithresh
      REAL , INTENT(IN) :: Rthresh
      REAL , INTENT(IN) :: Dz_ub
      LOGICAL , INTENT(IN) :: Ignore_cwise
      INTEGER :: Info
!
! Local variable declarations rewritten by SPAG
!
      REAL :: CABS1
      INTEGER :: cnt , i , j , x_state , y_prec_state , z_state
      REAL :: dxrat , dxratmax , dx_x , dyk , dzrat , dzratmax , dz_z , &
     &        eps , final_dx_x , final_dz_z , hugeval , incr_thresh ,   &
     &        normdx , normx , normy , prevnormdx , prev_dz_z , yk ,    &
     &        ymin
      LOGICAL :: incr_prec
      CHARACTER :: trans
      COMPLEX :: zdum
!
! End of declarations rewritten by SPAG
!
!     ..
!     .. Array Arguments
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
!     ..
!     .. Parameters ..
!     ..
!     .. External Subroutines ..
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. Statement Functions ..
!     ..
!     .. Statement Function Definitions ..
      CABS1(zdum) = ABS(REAL(zdum)) + ABS(AIMAG(zdum))
!     ..
!     .. Executable Statements ..
!
      IF ( Info/=0 ) RETURN
      trans = CHLA_TRANSTYPE(Trans_type)
      eps = SLAMCH('Epsilon')
      hugeval = SLAMCH('Overflow')
!     Force HUGEVAL to Inf
      hugeval = hugeval*hugeval
!     Using HUGEVAL may lead to spurious underflows.
      incr_thresh = REAL(N)*eps
!
      DO j = 1 , Nrhs
         y_prec_state = EXTRA_RESIDUAL
         IF ( y_prec_state==EXTRA_Y ) THEN
            DO i = 1 , N
               Y_tail(i) = 0.0
            ENDDO
         ENDIF
 
         dxrat = 0.0
         dxratmax = 0.0
         dzrat = 0.0
         dzratmax = 0.0
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
            CALL CCOPY(N,B(1,j),1,Res,1)
            IF ( y_prec_state==BASE_RESIDUAL ) THEN
               CALL CGEMV(trans,N,N,(-1.0E+0,0.0E+0),A,Lda,Y(1,j),1,    &
     &                    (1.0E+0,0.0E+0),Res,1)
            ELSEIF ( y_prec_state==EXTRA_RESIDUAL ) THEN
               CALL BLAS_CGEMV_X(Trans_type,N,N,(-1.0E+0,0.0E+0),A,Lda, &
     &                           Y(1,j),1,(1.0E+0,0.0E+0),Res,1,        &
     &                           Prec_type)
            ELSE
               CALL BLAS_CGEMV2_X(Trans_type,N,N,(-1.0E+0,0.0E+0),A,Lda,&
     &                            Y(1,j),Y_tail,1,(1.0E+0,0.0E+0),Res,1,&
     &                            Prec_type)
            ENDIF
 
!         XXX: RES is no longer needed.
            CALL CCOPY(N,Res,1,Dy,1)
            CALL CGETRS(trans,N,1,Af,Ldaf,Ipiv,Dy,N,Info)
!
!         Calculate relative changes DX_X, DZ_Z and ratios DXRAT, DZRAT.
!
            normx = 0.0E+0
            normy = 0.0E+0
            normdx = 0.0E+0
            dz_z = 0.0E+0
            ymin = hugeval
!
            DO i = 1 , N
               yk = CABS1(Y(i,j))
               dyk = CABS1(Dy(i))
 
               IF ( yk/=0.0E+0 ) THEN
                  dz_z = MAX(dz_z,dyk/yk)
               ELSEIF ( dyk/=0.0 ) THEN
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
 
            IF ( normx/=0.0 ) THEN
               dx_x = normdx/normx
            ELSEIF ( normdx==0.0 ) THEN
               dx_x = 0.0
            ELSE
               dx_x = hugeval
            ENDIF
 
            dxrat = normdx/prevnormdx
            dzrat = dz_z/prev_dz_z
!
!         Check termination criteria
!
            IF ( .NOT.Ignore_cwise .AND.                                &
     &           ymin*Rcond<incr_thresh*normy .AND.                     &
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
                  dzratmax = 0.0
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
!
!           Exit if both normwise and componentwise stopped working,
!           but if componentwise is unstable, let it go at least two
!           iterations.
!
            IF ( x_state/=WORKING_STATE ) THEN
               IF ( Ignore_cwise ) EXIT
               IF ( z_state==NOPROG_STATE .OR. z_state==CONV_STATE )    &
     &              EXIT
               IF ( z_state==UNSTABLE_STATE .AND. cnt>1 ) EXIT
            ENDIF
 
            IF ( incr_prec ) THEN
               incr_prec = .FALSE.
               y_prec_state = y_prec_state + 1
               DO i = 1 , N
                  Y_tail(i) = 0.0
               ENDDO
            ENDIF
 
            prevnormdx = normdx
            prev_dz_z = dz_z
!
!           Update soluton.
!
            IF ( y_prec_state<EXTRA_Y ) THEN
               CALL CAXPY(N,(1.0E+0,0.0E+0),Dy,1,Y(1,j),1)
            ELSE
               CALL CLA_WWADDW(N,Y(1,j),Y_tail,Dy)
            ENDIF
 
         ENDDO
!        Target of "IF (Z_STOP .AND. X_STOP)".  Sun's f77 won't EXIT.
!
!     Set final_* when cnt hits ithresh
!
         IF ( x_state==WORKING_STATE ) final_dx_x = dx_x
         IF ( z_state==WORKING_STATE ) final_dz_z = dz_z
!
!     Compute error bounds
!
 
         IF ( N_norms>=1 ) Errs_n(j,LA_LINRX_ERR_I)                     &
     &        = final_dx_x/(1-dxratmax)
         IF ( N_norms>=2 ) Errs_c(j,LA_LINRX_ERR_I)                     &
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
         CALL CCOPY(N,B(1,j),1,Res,1)
         CALL CGEMV(trans,N,N,(-1.0E+0,0.0E+0),A,Lda,Y(1,j),1,          &
     &              (1.0E+0,0.0E+0),Res,1)
 
         DO i = 1 , N
            Ayb(i) = CABS1(B(i,j))
         ENDDO
!
!     Compute abs(op(A_s))*abs(Y) + abs(B_s).
!
         CALL CLA_GEAMV(Trans_type,N,N,1.0E+0,A,Lda,Y(1,j),1,1.0E+0,Ayb,&
     &                  1)
 
         CALL CLA_LIN_BERR(N,N,1,Res,Ayb,Berr_out(j))
!
!     End of loop for each RHS.
!
      ENDDO
!
      END SUBROUTINE CLA_GERFSX_EXTENDED
