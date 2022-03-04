!*==ctrevc3.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b CTREVC3
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CTREVC3 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ctrevc3.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ctrevc3.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ctrevc3.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CTREVC3( SIDE, HOWMNY, SELECT, N, T, LDT, VL, LDVL, VR,
!                           LDVR, MM, M, WORK, LWORK, RWORK, LRWORK, INFO)
!
!       .. Scalar Arguments ..
!       CHARACTER          HOWMNY, SIDE
!       INTEGER            INFO, LDT, LDVL, LDVR, LWORK, M, MM, N
!       ..
!       .. Array Arguments ..
!       LOGICAL            SELECT( * )
!       REAL               RWORK( * )
!       COMPLEX            T( LDT, * ), VL( LDVL, * ), VR( LDVR, * ),
!      $                   WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CTREVC3 computes some or all of the right and/or left eigenvectors of
!> a complex upper triangular matrix T.
!> Matrices of this type are produced by the Schur factorization of
!> a complex general matrix:  A = Q*T*Q**H, as computed by CHSEQR.
!>
!> The right eigenvector x and the left eigenvector y of T corresponding
!> to an eigenvalue w are defined by:
!>
!>              T*x = w*x,     (y**H)*T = w*(y**H)
!>
!> where y**H denotes the conjugate transpose of the vector y.
!> The eigenvalues are not input to this routine, but are read directly
!> from the diagonal of T.
!>
!> This routine returns the matrices X and/or Y of right and left
!> eigenvectors of T, or the products Q*X and/or Q*Y, where Q is an
!> input matrix. If Q is the unitary factor that reduces a matrix A to
!> Schur form T, then Q*X and Q*Y are the matrices of right and left
!> eigenvectors of A.
!>
!> This uses a Level 3 BLAS version of the back transformation.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] SIDE
!> \verbatim
!>          SIDE is CHARACTER*1
!>          = 'R':  compute right eigenvectors only;
!>          = 'L':  compute left eigenvectors only;
!>          = 'B':  compute both right and left eigenvectors.
!> \endverbatim
!>
!> \param[in] HOWMNY
!> \verbatim
!>          HOWMNY is CHARACTER*1
!>          = 'A':  compute all right and/or left eigenvectors;
!>          = 'B':  compute all right and/or left eigenvectors,
!>                  backtransformed using the matrices supplied in
!>                  VR and/or VL;
!>          = 'S':  compute selected right and/or left eigenvectors,
!>                  as indicated by the logical array SELECT.
!> \endverbatim
!>
!> \param[in] SELECT
!> \verbatim
!>          SELECT is LOGICAL array, dimension (N)
!>          If HOWMNY = 'S', SELECT specifies the eigenvectors to be
!>          computed.
!>          The eigenvector corresponding to the j-th eigenvalue is
!>          computed if SELECT(j) = .TRUE..
!>          Not referenced if HOWMNY = 'A' or 'B'.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix T. N >= 0.
!> \endverbatim
!>
!> \param[in,out] T
!> \verbatim
!>          T is COMPLEX array, dimension (LDT,N)
!>          The upper triangular matrix T.  T is modified, but restored
!>          on exit.
!> \endverbatim
!>
!> \param[in] LDT
!> \verbatim
!>          LDT is INTEGER
!>          The leading dimension of the array T. LDT >= max(1,N).
!> \endverbatim
!>
!> \param[in,out] VL
!> \verbatim
!>          VL is COMPLEX array, dimension (LDVL,MM)
!>          On entry, if SIDE = 'L' or 'B' and HOWMNY = 'B', VL must
!>          contain an N-by-N matrix Q (usually the unitary matrix Q of
!>          Schur vectors returned by CHSEQR).
!>          On exit, if SIDE = 'L' or 'B', VL contains:
!>          if HOWMNY = 'A', the matrix Y of left eigenvectors of T;
!>          if HOWMNY = 'B', the matrix Q*Y;
!>          if HOWMNY = 'S', the left eigenvectors of T specified by
!>                           SELECT, stored consecutively in the columns
!>                           of VL, in the same order as their
!>                           eigenvalues.
!>          Not referenced if SIDE = 'R'.
!> \endverbatim
!>
!> \param[in] LDVL
!> \verbatim
!>          LDVL is INTEGER
!>          The leading dimension of the array VL.
!>          LDVL >= 1, and if SIDE = 'L' or 'B', LDVL >= N.
!> \endverbatim
!>
!> \param[in,out] VR
!> \verbatim
!>          VR is COMPLEX array, dimension (LDVR,MM)
!>          On entry, if SIDE = 'R' or 'B' and HOWMNY = 'B', VR must
!>          contain an N-by-N matrix Q (usually the unitary matrix Q of
!>          Schur vectors returned by CHSEQR).
!>          On exit, if SIDE = 'R' or 'B', VR contains:
!>          if HOWMNY = 'A', the matrix X of right eigenvectors of T;
!>          if HOWMNY = 'B', the matrix Q*X;
!>          if HOWMNY = 'S', the right eigenvectors of T specified by
!>                           SELECT, stored consecutively in the columns
!>                           of VR, in the same order as their
!>                           eigenvalues.
!>          Not referenced if SIDE = 'L'.
!> \endverbatim
!>
!> \param[in] LDVR
!> \verbatim
!>          LDVR is INTEGER
!>          The leading dimension of the array VR.
!>          LDVR >= 1, and if SIDE = 'R' or 'B', LDVR >= N.
!> \endverbatim
!>
!> \param[in] MM
!> \verbatim
!>          MM is INTEGER
!>          The number of columns in the arrays VL and/or VR. MM >= M.
!> \endverbatim
!>
!> \param[out] M
!> \verbatim
!>          M is INTEGER
!>          The number of columns in the arrays VL and/or VR actually
!>          used to store the eigenvectors.
!>          If HOWMNY = 'A' or 'B', M is set to N.
!>          Each selected eigenvector occupies one column.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX array, dimension (MAX(1,LWORK))
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The dimension of array WORK. LWORK >= max(1,2*N).
!>          For optimum performance, LWORK >= N + 2*N*NB, where NB is
!>          the optimal blocksize.
!>
!>          If LWORK = -1, then a workspace query is assumed; the routine
!>          only calculates the optimal size of the WORK array, returns
!>          this value as the first entry of the WORK array, and no error
!>          message related to LWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is REAL array, dimension (LRWORK)
!> \endverbatim
!>
!> \param[in] LRWORK
!> \verbatim
!>          LRWORK is INTEGER
!>          The dimension of array RWORK. LRWORK >= max(1,N).
!>
!>          If LRWORK = -1, then a workspace query is assumed; the routine
!>          only calculates the optimal size of the RWORK array, returns
!>          this value as the first entry of the RWORK array, and no error
!>          message related to LRWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value
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
!> \date November 2017
!
!  @generated from ztrevc3.f, fortran z -> c, Tue Apr 19 01:47:44 2016
!
!> \ingroup complexOTHERcomputational
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  The algorithm used in this program is basically backward (forward)
!>  substitution, with scaling to make the the code robust against
!>  possible overflow.
!>
!>  Each eigenvector is normalized so that the element of largest
!>  magnitude has magnitude 1; here the magnitude of a complex number
!>  (x,y) is taken to be |x| + |y|.
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE CTREVC3(Side,Howmny,Select,N,T,Ldt,Vl,Ldvl,Vr,Ldvr,Mm, &
     &                   M,Work,Lwork,Rwork,Lrwork,Info)
      USE S_CCOPY
      USE S_CGEMM
      USE S_CGEMV
      USE S_CLACPY
      USE S_CLASET
      USE S_CLATRS
      USE S_CSSCAL
      USE S_ICAMAX
      USE S_ILAENV
      USE S_LSAME
      USE S_SCASUM
      USE S_SLABAD
      USE S_SLAMCH
      USE S_XERBLA
      IMPLICIT NONE
!*--CTREVC3264
!
! PARAMETER definitions rewritten by SPAG
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
      COMPLEX , PARAMETER  ::  CZERO = (0.0E+0,0.0E+0) ,                &
     &                         CONE = (1.0E+0,0.0E+0)
      INTEGER , PARAMETER  ::  NBMIN = 8 , NBMAX = 128
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Side
      CHARACTER :: Howmny
      LOGICAL , INTENT(IN) , DIMENSION(*) :: Select
      INTEGER :: N
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldvl,*) :: Vl
      INTEGER :: Ldvl
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldvr,*) :: Vr
      INTEGER :: Ldvr
      INTEGER , INTENT(IN) :: Mm
      INTEGER , INTENT(INOUT) :: M
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      REAL , DIMENSION(*) :: Rwork
      INTEGER , INTENT(IN) :: Lrwork
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      LOGICAL :: allv , bothv , leftv , lquery , over , rightv , somev
      REAL :: CABS1
      COMPLEX :: cdum
      INTEGER :: i , ii , is , iv , j , k , ki , maxwrk , nb
      REAL :: ovfl , remax , scale , smin , smlnum , ulp , unfl
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
!     .. Statement Functions ..
!     ..
!     .. Statement Function definitions ..
      CABS1(cdum) = ABS(REAL(cdum)) + ABS(AIMAG(cdum))
!     ..
!     .. Executable Statements ..
!
!     Decode and test the input parameters
!
      bothv = LSAME(Side,'B')
      rightv = LSAME(Side,'R') .OR. bothv
      leftv = LSAME(Side,'L') .OR. bothv
!
      allv = LSAME(Howmny,'A')
      over = LSAME(Howmny,'B')
      somev = LSAME(Howmny,'S')
!
!     Set M to the number of columns required to store the selected
!     eigenvectors.
!
      IF ( somev ) THEN
         M = 0
         DO j = 1 , N
            IF ( Select(j) ) M = M + 1
         ENDDO
      ELSE
         M = N
      ENDIF
!
      Info = 0
      nb = ILAENV(1,'CTREVC',Side//Howmny,N,-1,-1,-1)
      maxwrk = N + 2*N*nb
      Work(1) = maxwrk
      Rwork(1) = N
      lquery = (Lwork==-1 .OR. Lrwork==-1)
      IF ( .NOT.rightv .AND. .NOT.leftv ) THEN
         Info = -1
      ELSEIF ( .NOT.allv .AND. .NOT.over .AND. .NOT.somev ) THEN
         Info = -2
      ELSEIF ( N<0 ) THEN
         Info = -4
      ELSEIF ( Ldt<MAX(1,N) ) THEN
         Info = -6
      ELSEIF ( Ldvl<1 .OR. (leftv .AND. Ldvl<N) ) THEN
         Info = -8
      ELSEIF ( Ldvr<1 .OR. (rightv .AND. Ldvr<N) ) THEN
         Info = -10
      ELSEIF ( Mm<M ) THEN
         Info = -11
      ELSEIF ( Lwork<MAX(1,2*N) .AND. .NOT.lquery ) THEN
         Info = -14
      ELSEIF ( Lrwork<MAX(1,N) .AND. .NOT.lquery ) THEN
         Info = -16
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('CTREVC3',-Info)
         RETURN
      ELSEIF ( lquery ) THEN
         RETURN
      ENDIF
!
!     Quick return if possible.
!
      IF ( N==0 ) RETURN
!
!     Use blocked version of back-transformation if sufficient workspace.
!     Zero-out the workspace to avoid potential NaN propagation.
!
      IF ( over .AND. Lwork>=N+2*N*NBMIN ) THEN
         nb = (Lwork-N)/(2*N)
         nb = MIN(nb,NBMAX)
         CALL CLASET('F',N,1+2*nb,CZERO,CZERO,Work,N)
      ELSE
         nb = 1
      ENDIF
!
!     Set the constants to control overflow.
!
      unfl = SLAMCH('Safe minimum')
      ovfl = ONE/unfl
      CALL SLABAD(unfl,ovfl)
      ulp = SLAMCH('Precision')
      smlnum = unfl*(N/ulp)
!
!     Store the diagonal elements of T in working array WORK.
!
      DO i = 1 , N
         Work(i) = T(i,i)
      ENDDO
!
!     Compute 1-norm of each column of strictly upper triangular
!     part of T to control overflow in triangular solver.
!
      Rwork(1) = ZERO
      DO j = 2 , N
         Rwork(j) = SCASUM(j-1,T(1,j),1)
      ENDDO
!
      IF ( rightv ) THEN
!
!        ============================================================
!        Compute right eigenvectors.
!
!        IV is index of column in current block.
!        Non-blocked version always uses IV=NB=1;
!        blocked     version starts with IV=NB, goes down to 1.
!        (Note the "0-th" column is used to store the original diagonal.)
         iv = nb
         is = M
         DO ki = N , 1 , -1
            IF ( somev ) THEN
               IF ( .NOT.Select(ki) ) CYCLE
            ENDIF
            smin = MAX(ulp*(CABS1(T(ki,ki))),smlnum)
!
!           --------------------------------------------------------
!           Complex right eigenvector
!
            Work(ki+iv*N) = CONE
!
!           Form right-hand side.
!
            DO k = 1 , ki - 1
               Work(k+iv*N) = -T(k,ki)
            ENDDO
!
!           Solve upper triangular system:
!           [ T(1:KI-1,1:KI-1) - T(KI,KI) ]*X = SCALE*WORK.
!
            DO k = 1 , ki - 1
               T(k,k) = T(k,k) - T(ki,ki)
               IF ( CABS1(T(k,k))<smin ) T(k,k) = smin
            ENDDO
!
            IF ( ki>1 ) THEN
               CALL CLATRS('Upper','No transpose','Non-unit','Y',ki-1,T,&
     &                     Ldt,Work(1+iv*N),scale,Rwork,Info)
               Work(ki+iv*N) = scale
            ENDIF
!
!           Copy the vector x or Q*x to VR and normalize.
!
            IF ( .NOT.over ) THEN
!              ------------------------------
!              no back-transform: copy x to VR and normalize.
               CALL CCOPY(ki,Work(1+iv*N),1,Vr(1,is),1)
!
               ii = ICAMAX(ki,Vr(1,is),1)
               remax = ONE/CABS1(Vr(ii,is))
               CALL CSSCAL(ki,remax,Vr(1,is),1)
!
               DO k = ki + 1 , N
                  Vr(k,is) = CZERO
               ENDDO
!
            ELSEIF ( nb==1 ) THEN
!              ------------------------------
!              version 1: back-transform each vector with GEMV, Q*x.
               IF ( ki>1 ) CALL CGEMV('N',N,ki-1,CONE,Vr,Ldvr,          &
     &                                Work(1+iv*N),1,CMPLX(scale),      &
     &                                Vr(1,ki),1)
!
               ii = ICAMAX(N,Vr(1,ki),1)
               remax = ONE/CABS1(Vr(ii,ki))
               CALL CSSCAL(N,remax,Vr(1,ki),1)
!
            ELSE
!              ------------------------------
!              version 2: back-transform block of vectors with GEMM
!              zero out below vector
               DO k = ki + 1 , N
                  Work(k+iv*N) = CZERO
               ENDDO
!
!              Columns IV:NB of work are valid vectors.
!              When the number of vectors stored reaches NB,
!              or if this was last vector, do the GEMM
               IF ( (iv==1) .OR. (ki==1) ) THEN
                  CALL CGEMM('N','N',N,nb-iv+1,ki+nb-iv,CONE,Vr,Ldvr,   &
     &                       Work(1+(iv)*N),N,CZERO,Work(1+(nb+iv)*N),N)
!                 normalize vectors
                  DO k = iv , nb
                     ii = ICAMAX(N,Work(1+(nb+k)*N),1)
                     remax = ONE/CABS1(Work(ii+(nb+k)*N))
                     CALL CSSCAL(N,remax,Work(1+(nb+k)*N),1)
                  ENDDO
                  CALL CLACPY('F',N,nb-iv+1,Work(1+(nb+iv)*N),N,Vr(1,ki)&
     &                        ,Ldvr)
                  iv = nb
               ELSE
                  iv = iv - 1
               ENDIF
            ENDIF
!
!           Restore the original diagonal elements of T.
!
            DO k = 1 , ki - 1
               T(k,k) = Work(k)
            ENDDO
!
            is = is - 1
         ENDDO
      ENDIF
!
      IF ( leftv ) THEN
!
!        ============================================================
!        Compute left eigenvectors.
!
!        IV is index of column in current block.
!        Non-blocked version always uses IV=1;
!        blocked     version starts with IV=1, goes up to NB.
!        (Note the "0-th" column is used to store the original diagonal.)
         iv = 1
         is = 1
         DO ki = 1 , N
!
            IF ( somev ) THEN
               IF ( .NOT.Select(ki) ) CYCLE
            ENDIF
            smin = MAX(ulp*(CABS1(T(ki,ki))),smlnum)
!
!           --------------------------------------------------------
!           Complex left eigenvector
!
            Work(ki+iv*N) = CONE
!
!           Form right-hand side.
!
            DO k = ki + 1 , N
               Work(k+iv*N) = -CONJG(T(ki,k))
            ENDDO
!
!           Solve conjugate-transposed triangular system:
!           [ T(KI+1:N,KI+1:N) - T(KI,KI) ]**H * X = SCALE*WORK.
!
            DO k = ki + 1 , N
               T(k,k) = T(k,k) - T(ki,ki)
               IF ( CABS1(T(k,k))<smin ) T(k,k) = smin
            ENDDO
!
            IF ( ki<N ) THEN
               CALL CLATRS('Upper','Conjugate transpose','Non-unit','Y',&
     &                     N-ki,T(ki+1,ki+1),Ldt,Work(ki+1+iv*N),scale, &
     &                     Rwork,Info)
               Work(ki+iv*N) = scale
            ENDIF
!
!           Copy the vector x or Q*x to VL and normalize.
!
            IF ( .NOT.over ) THEN
!              ------------------------------
!              no back-transform: copy x to VL and normalize.
               CALL CCOPY(N-ki+1,Work(ki+iv*N),1,Vl(ki,is),1)
!
               ii = ICAMAX(N-ki+1,Vl(ki,is),1) + ki - 1
               remax = ONE/CABS1(Vl(ii,is))
               CALL CSSCAL(N-ki+1,remax,Vl(ki,is),1)
!
               DO k = 1 , ki - 1
                  Vl(k,is) = CZERO
               ENDDO
!
            ELSEIF ( nb==1 ) THEN
!              ------------------------------
!              version 1: back-transform each vector with GEMV, Q*x.
               IF ( ki<N ) CALL CGEMV('N',N,N-ki,CONE,Vl(1,ki+1),Ldvl,  &
     &                                Work(ki+1+iv*N),1,CMPLX(scale),   &
     &                                Vl(1,ki),1)
!
               ii = ICAMAX(N,Vl(1,ki),1)
               remax = ONE/CABS1(Vl(ii,ki))
               CALL CSSCAL(N,remax,Vl(1,ki),1)
!
            ELSE
!              ------------------------------
!              version 2: back-transform block of vectors with GEMM
!              zero out above vector
!              could go from KI-NV+1 to KI-1
               DO k = 1 , ki - 1
                  Work(k+iv*N) = CZERO
               ENDDO
!
!              Columns 1:IV of work are valid vectors.
!              When the number of vectors stored reaches NB,
!              or if this was last vector, do the GEMM
               IF ( (iv==nb) .OR. (ki==N) ) THEN
                  CALL CGEMM('N','N',N,iv,N-ki+iv,CONE,Vl(1,ki-iv+1),   &
     &                       Ldvl,Work(ki-iv+1+(1)*N),N,CZERO,          &
     &                       Work(1+(nb+1)*N),N)
!                 normalize vectors
                  DO k = 1 , iv
                     ii = ICAMAX(N,Work(1+(nb+k)*N),1)
                     remax = ONE/CABS1(Work(ii+(nb+k)*N))
                     CALL CSSCAL(N,remax,Work(1+(nb+k)*N),1)
                  ENDDO
                  CALL CLACPY('F',N,iv,Work(1+(nb+1)*N),N,Vl(1,ki-iv+1),&
     &                        Ldvl)
                  iv = 1
               ELSE
                  iv = iv + 1
               ENDIF
            ENDIF
!
!           Restore the original diagonal elements of T.
!
            DO k = ki + 1 , N
               T(k,k) = Work(k)
            ENDDO
!
            is = is + 1
         ENDDO
      ENDIF
!
!
!     End of CTREVC3
!
      END SUBROUTINE CTREVC3
