!*==strevc.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b STREVC
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download STREVC + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/strevc.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/strevc.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/strevc.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE STREVC( SIDE, HOWMNY, SELECT, N, T, LDT, VL, LDVL, VR,
!                          LDVR, MM, M, WORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          HOWMNY, SIDE
!       INTEGER            INFO, LDT, LDVL, LDVR, M, MM, N
!       ..
!       .. Array Arguments ..
!       LOGICAL            SELECT( * )
!       REAL               T( LDT, * ), VL( LDVL, * ), VR( LDVR, * ),
!      $                   WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> STREVC computes some or all of the right and/or left eigenvectors of
!> a real upper quasi-triangular matrix T.
!> Matrices of this type are produced by the Schur factorization of
!> a real general matrix:  A = Q*T*Q**T, as computed by SHSEQR.
!>
!> The right eigenvector x and the left eigenvector y of T corresponding
!> to an eigenvalue w are defined by:
!>
!>    T*x = w*x,     (y**H)*T = w*(y**H)
!>
!> where y**H denotes the conjugate transpose of y.
!> The eigenvalues are not input to this routine, but are read directly
!> from the diagonal blocks of T.
!>
!> This routine returns the matrices X and/or Y of right and left
!> eigenvectors of T, or the products Q*X and/or Q*Y, where Q is an
!> input matrix.  If Q is the orthogonal factor that reduces a matrix
!> A to Schur form T, then Q*X and Q*Y are the matrices of right and
!> left eigenvectors of A.
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
!>                  backtransformed by the matrices in VR and/or VL;
!>          = 'S':  compute selected right and/or left eigenvectors,
!>                  as indicated by the logical array SELECT.
!> \endverbatim
!>
!> \param[in,out] SELECT
!> \verbatim
!>          SELECT is LOGICAL array, dimension (N)
!>          If HOWMNY = 'S', SELECT specifies the eigenvectors to be
!>          computed.
!>          If w(j) is a real eigenvalue, the corresponding real
!>          eigenvector is computed if SELECT(j) is .TRUE..
!>          If w(j) and w(j+1) are the real and imaginary parts of a
!>          complex eigenvalue, the corresponding complex eigenvector is
!>          computed if either SELECT(j) or SELECT(j+1) is .TRUE., and
!>          on exit SELECT(j) is set to .TRUE. and SELECT(j+1) is set to
!>          .FALSE..
!>          Not referenced if HOWMNY = 'A' or 'B'.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix T. N >= 0.
!> \endverbatim
!>
!> \param[in] T
!> \verbatim
!>          T is REAL array, dimension (LDT,N)
!>          The upper quasi-triangular matrix T in Schur canonical form.
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
!>          VL is REAL array, dimension (LDVL,MM)
!>          On entry, if SIDE = 'L' or 'B' and HOWMNY = 'B', VL must
!>          contain an N-by-N matrix Q (usually the orthogonal matrix Q
!>          of Schur vectors returned by SHSEQR).
!>          On exit, if SIDE = 'L' or 'B', VL contains:
!>          if HOWMNY = 'A', the matrix Y of left eigenvectors of T;
!>          if HOWMNY = 'B', the matrix Q*Y;
!>          if HOWMNY = 'S', the left eigenvectors of T specified by
!>                           SELECT, stored consecutively in the columns
!>                           of VL, in the same order as their
!>                           eigenvalues.
!>          A complex eigenvector corresponding to a complex eigenvalue
!>          is stored in two consecutive columns, the first holding the
!>          real part, and the second the imaginary part.
!>          Not referenced if SIDE = 'R'.
!> \endverbatim
!>
!> \param[in] LDVL
!> \verbatim
!>          LDVL is INTEGER
!>          The leading dimension of the array VL.  LDVL >= 1, and if
!>          SIDE = 'L' or 'B', LDVL >= N.
!> \endverbatim
!>
!> \param[in,out] VR
!> \verbatim
!>          VR is REAL array, dimension (LDVR,MM)
!>          On entry, if SIDE = 'R' or 'B' and HOWMNY = 'B', VR must
!>          contain an N-by-N matrix Q (usually the orthogonal matrix Q
!>          of Schur vectors returned by SHSEQR).
!>          On exit, if SIDE = 'R' or 'B', VR contains:
!>          if HOWMNY = 'A', the matrix X of right eigenvectors of T;
!>          if HOWMNY = 'B', the matrix Q*X;
!>          if HOWMNY = 'S', the right eigenvectors of T specified by
!>                           SELECT, stored consecutively in the columns
!>                           of VR, in the same order as their
!>                           eigenvalues.
!>          A complex eigenvector corresponding to a complex eigenvalue
!>          is stored in two consecutive columns, the first holding the
!>          real part and the second the imaginary part.
!>          Not referenced if SIDE = 'L'.
!> \endverbatim
!>
!> \param[in] LDVR
!> \verbatim
!>          LDVR is INTEGER
!>          The leading dimension of the array VR.  LDVR >= 1, and if
!>          SIDE = 'R' or 'B', LDVR >= N.
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
!>          Each selected real eigenvector occupies one column and each
!>          selected complex eigenvector occupies two columns.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is REAL array, dimension (3*N)
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
!> \date December 2016
!
!> \ingroup realOTHERcomputational
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
      SUBROUTINE STREVC(Side,Howmny,Select,N,T,Ldt,Vl,Ldvl,Vr,Ldvr,Mm,M,&
     &                  Work,Info)
      USE S_ISAMAX
      USE S_LSAME
      USE S_SAXPY
      USE S_SCOPY
      USE S_SDOT
      USE S_SGEMV
      USE S_SLABAD
      USE S_SLALN2
      USE S_SLAMCH
      USE S_SSCAL
      USE S_XERBLA
      IMPLICIT NONE
!*--STREVC237
!
! PARAMETER definitions rewritten by SPAG
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Side
      CHARACTER :: Howmny
      LOGICAL , INTENT(INOUT) , DIMENSION(*) :: Select
      INTEGER :: N
      REAL , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      REAL , INTENT(INOUT) , DIMENSION(Ldvl,*) :: Vl
      INTEGER :: Ldvl
      REAL , INTENT(INOUT) , DIMENSION(Ldvr,*) :: Vr
      INTEGER :: Ldvr
      INTEGER , INTENT(IN) :: Mm
      INTEGER , INTENT(INOUT) :: M
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      LOGICAL :: allv , bothv , leftv , over , pair , rightv , somev
      REAL :: beta , bignum , emax , ovfl , rec , remax , scale , smin ,&
     &        smlnum , ulp , unfl , vcrit , vmax , wi , wr , xnorm
      INTEGER :: i , ierr , ii , ip , is , j , j1 , j2 , jnxt , k , ki ,&
     &           n2
      REAL , DIMENSION(2,2) :: x
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
!     .. Local Arrays ..
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
      Info = 0
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
      ELSE
!
!        Set M to the number of columns required to store the selected
!        eigenvectors, standardize the array SELECT if necessary, and
!        test MM.
!
         IF ( somev ) THEN
            M = 0
            pair = .FALSE.
            DO j = 1 , N
               IF ( pair ) THEN
                  pair = .FALSE.
                  Select(j) = .FALSE.
               ELSEIF ( j<N ) THEN
                  IF ( T(j+1,j)==ZERO ) THEN
                     IF ( Select(j) ) M = M + 1
                  ELSE
                     pair = .TRUE.
                     IF ( Select(j) .OR. Select(j+1) ) THEN
                        Select(j) = .TRUE.
                        M = M + 2
                     ENDIF
                  ENDIF
               ELSE
                  IF ( Select(N) ) M = M + 1
               ENDIF
            ENDDO
         ELSE
            M = N
         ENDIF
!
         IF ( Mm<M ) Info = -11
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('STREVC',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible.
!
      IF ( N==0 ) RETURN
!
!     Set the constants to control overflow.
!
      unfl = SLAMCH('Safe minimum')
      ovfl = ONE/unfl
      CALL SLABAD(unfl,ovfl)
      ulp = SLAMCH('Precision')
      smlnum = unfl*(N/ulp)
      bignum = (ONE-ulp)/smlnum
!
!     Compute 1-norm of each column of strictly upper triangular
!     part of T to control overflow in triangular solver.
!
      Work(1) = ZERO
      DO j = 2 , N
         Work(j) = ZERO
         DO i = 1 , j - 1
            Work(j) = Work(j) + ABS(T(i,j))
         ENDDO
      ENDDO
!
!     Index IP is used to specify the real or complex eigenvalue:
!       IP = 0, real eigenvalue,
!            1, first of conjugate complex pair: (wr,wi)
!           -1, second of conjugate complex pair: (wr,wi)
!
      n2 = 2*N
!
      IF ( rightv ) THEN
!
!        Compute right eigenvectors.
!
         ip = 0
         is = M
         DO ki = N , 1 , -1
!
            IF ( ip/=1 ) THEN
               IF ( ki/=1 ) THEN
                  IF ( T(ki,ki-1)/=ZERO ) ip = -1
               ENDIF
!
               IF ( somev ) THEN
                  IF ( ip==0 ) THEN
                     IF ( .NOT.Select(ki) ) GOTO 20
                  ELSEIF ( .NOT.Select(ki-1) ) THEN
                     GOTO 20
                  ENDIF
               ENDIF
!
!           Compute the KI-th eigenvalue (WR,WI).
!
               wr = T(ki,ki)
               wi = ZERO
               IF ( ip/=0 ) wi = SQRT(ABS(T(ki,ki-1)))                  &
     &                           *SQRT(ABS(T(ki-1,ki)))
               smin = MAX(ulp*(ABS(wr)+ABS(wi)),smlnum)
!
               IF ( ip==0 ) THEN
!
!              Real right eigenvector
!
                  Work(ki+N) = ONE
!
!              Form right-hand side
!
                  DO k = 1 , ki - 1
                     Work(k+N) = -T(k,ki)
                  ENDDO
!
!              Solve the upper quasi-triangular system:
!                 (T(1:KI-1,1:KI-1) - WR)*X = SCALE*WORK.
!
                  jnxt = ki - 1
                  DO j = ki - 1 , 1 , -1
                     IF ( j<=jnxt ) THEN
                        j1 = j
                        j2 = j
                        jnxt = j - 1
                        IF ( j>1 ) THEN
                           IF ( T(j,j-1)/=ZERO ) THEN
                              j1 = j - 1
                              jnxt = j - 2
                           ENDIF
                        ENDIF
!
                        IF ( j1==j2 ) THEN
!
!                    1-by-1 diagonal block
!
                           CALL SLALN2(.FALSE.,1,1,smin,ONE,T(j,j),Ldt, &
     &                                 ONE,ONE,Work(j+N),N,wr,ZERO,x,2, &
     &                                 scale,xnorm,ierr)
!
!                    Scale X(1,1) to avoid overflow when updating
!                    the right-hand side.
!
                           IF ( xnorm>ONE ) THEN
                              IF ( Work(j)>bignum/xnorm ) THEN
                                 x(1,1) = x(1,1)/xnorm
                                 scale = scale/xnorm
                              ENDIF
                           ENDIF
!
!                    Scale if necessary
!
                           IF ( scale/=ONE )                            &
     &                          CALL SSCAL(ki,scale,Work(1+N),1)
                           Work(j+N) = x(1,1)
!
!                    Update right-hand side
!
                           CALL SAXPY(j-1,-x(1,1),T(1,j),1,Work(1+N),1)
!
                        ELSE
!
!                    2-by-2 diagonal block
!
                           CALL SLALN2(.FALSE.,2,1,smin,ONE,T(j-1,j-1), &
     &                                 Ldt,ONE,ONE,Work(j-1+N),N,wr,    &
     &                                 ZERO,x,2,scale,xnorm,ierr)
!
!                    Scale X(1,1) and X(2,1) to avoid overflow when
!                    updating the right-hand side.
!
                           IF ( xnorm>ONE ) THEN
                              beta = MAX(Work(j-1),Work(j))
                              IF ( beta>bignum/xnorm ) THEN
                                 x(1,1) = x(1,1)/xnorm
                                 x(2,1) = x(2,1)/xnorm
                                 scale = scale/xnorm
                              ENDIF
                           ENDIF
!
!                    Scale if necessary
!
                           IF ( scale/=ONE )                            &
     &                          CALL SSCAL(ki,scale,Work(1+N),1)
                           Work(j-1+N) = x(1,1)
                           Work(j+N) = x(2,1)
!
!                    Update right-hand side
!
                           CALL SAXPY(j-2,-x(1,1),T(1,j-1),1,Work(1+N), &
     &                                1)
                           CALL SAXPY(j-2,-x(2,1),T(1,j),1,Work(1+N),1)
                        ENDIF
                     ENDIF
                  ENDDO
!
!              Copy the vector x or Q*x to VR and normalize.
!
                  IF ( .NOT.over ) THEN
                     CALL SCOPY(ki,Work(1+N),1,Vr(1,is),1)
!
                     ii = ISAMAX(ki,Vr(1,is),1)
                     remax = ONE/ABS(Vr(ii,is))
                     CALL SSCAL(ki,remax,Vr(1,is),1)
!
                     DO k = ki + 1 , N
                        Vr(k,is) = ZERO
                     ENDDO
                  ELSE
                     IF ( ki>1 )                                        &
     &                    CALL SGEMV('N',N,ki-1,ONE,Vr,Ldvr,Work(1+N),1,&
     &                    Work(ki+N),Vr(1,ki),1)
!
                     ii = ISAMAX(N,Vr(1,ki),1)
                     remax = ONE/ABS(Vr(ii,ki))
                     CALL SSCAL(N,remax,Vr(1,ki),1)
                  ENDIF
!
               ELSE
!
!              Complex right eigenvector.
!
!              Initial solve
!                [ (T(KI-1,KI-1) T(KI-1,KI) ) - (WR + I* WI)]*X = 0.
!                [ (T(KI,KI-1)   T(KI,KI)   )               ]
!
                  IF ( ABS(T(ki-1,ki))>=ABS(T(ki,ki-1)) ) THEN
                     Work(ki-1+N) = ONE
                     Work(ki+n2) = wi/T(ki-1,ki)
                  ELSE
                     Work(ki-1+N) = -wi/T(ki,ki-1)
                     Work(ki+n2) = ONE
                  ENDIF
                  Work(ki+N) = ZERO
                  Work(ki-1+n2) = ZERO
!
!              Form right-hand side
!
                  DO k = 1 , ki - 2
                     Work(k+N) = -Work(ki-1+N)*T(k,ki-1)
                     Work(k+n2) = -Work(ki+n2)*T(k,ki)
                  ENDDO
!
!              Solve upper quasi-triangular system:
!              (T(1:KI-2,1:KI-2) - (WR+i*WI))*X = SCALE*(WORK+i*WORK2)
!
                  jnxt = ki - 2
                  DO j = ki - 2 , 1 , -1
                     IF ( j<=jnxt ) THEN
                        j1 = j
                        j2 = j
                        jnxt = j - 1
                        IF ( j>1 ) THEN
                           IF ( T(j,j-1)/=ZERO ) THEN
                              j1 = j - 1
                              jnxt = j - 2
                           ENDIF
                        ENDIF
!
                        IF ( j1==j2 ) THEN
!
!                    1-by-1 diagonal block
!
                           CALL SLALN2(.FALSE.,1,2,smin,ONE,T(j,j),Ldt, &
     &                                 ONE,ONE,Work(j+N),N,wr,wi,x,2,   &
     &                                 scale,xnorm,ierr)
!
!                    Scale X(1,1) and X(1,2) to avoid overflow when
!                    updating the right-hand side.
!
                           IF ( xnorm>ONE ) THEN
                              IF ( Work(j)>bignum/xnorm ) THEN
                                 x(1,1) = x(1,1)/xnorm
                                 x(1,2) = x(1,2)/xnorm
                                 scale = scale/xnorm
                              ENDIF
                           ENDIF
!
!                    Scale if necessary
!
                           IF ( scale/=ONE ) THEN
                              CALL SSCAL(ki,scale,Work(1+N),1)
                              CALL SSCAL(ki,scale,Work(1+n2),1)
                           ENDIF
                           Work(j+N) = x(1,1)
                           Work(j+n2) = x(1,2)
!
!                    Update the right-hand side
!
                           CALL SAXPY(j-1,-x(1,1),T(1,j),1,Work(1+N),1)
                           CALL SAXPY(j-1,-x(1,2),T(1,j),1,Work(1+n2),1)
!
                        ELSE
!
!                    2-by-2 diagonal block
!
                           CALL SLALN2(.FALSE.,2,2,smin,ONE,T(j-1,j-1), &
     &                                 Ldt,ONE,ONE,Work(j-1+N),N,wr,wi, &
     &                                 x,2,scale,xnorm,ierr)
!
!                    Scale X to avoid overflow when updating
!                    the right-hand side.
!
                           IF ( xnorm>ONE ) THEN
                              beta = MAX(Work(j-1),Work(j))
                              IF ( beta>bignum/xnorm ) THEN
                                 rec = ONE/xnorm
                                 x(1,1) = x(1,1)*rec
                                 x(1,2) = x(1,2)*rec
                                 x(2,1) = x(2,1)*rec
                                 x(2,2) = x(2,2)*rec
                                 scale = scale*rec
                              ENDIF
                           ENDIF
!
!                    Scale if necessary
!
                           IF ( scale/=ONE ) THEN
                              CALL SSCAL(ki,scale,Work(1+N),1)
                              CALL SSCAL(ki,scale,Work(1+n2),1)
                           ENDIF
                           Work(j-1+N) = x(1,1)
                           Work(j+N) = x(2,1)
                           Work(j-1+n2) = x(1,2)
                           Work(j+n2) = x(2,2)
!
!                    Update the right-hand side
!
                           CALL SAXPY(j-2,-x(1,1),T(1,j-1),1,Work(1+N), &
     &                                1)
                           CALL SAXPY(j-2,-x(2,1),T(1,j),1,Work(1+N),1)
                           CALL SAXPY(j-2,-x(1,2),T(1,j-1),1,Work(1+n2),&
     &                                1)
                           CALL SAXPY(j-2,-x(2,2),T(1,j),1,Work(1+n2),1)
                        ENDIF
                     ENDIF
                  ENDDO
!
!              Copy the vector x or Q*x to VR and normalize.
!
                  IF ( .NOT.over ) THEN
                     CALL SCOPY(ki,Work(1+N),1,Vr(1,is-1),1)
                     CALL SCOPY(ki,Work(1+n2),1,Vr(1,is),1)
!
                     emax = ZERO
                     DO k = 1 , ki
                        emax = MAX(emax,ABS(Vr(k,is-1))+ABS(Vr(k,is)))
                     ENDDO
!
                     remax = ONE/emax
                     CALL SSCAL(ki,remax,Vr(1,is-1),1)
                     CALL SSCAL(ki,remax,Vr(1,is),1)
!
                     DO k = ki + 1 , N
                        Vr(k,is-1) = ZERO
                        Vr(k,is) = ZERO
                     ENDDO
!
                  ELSE
!
                     IF ( ki>2 ) THEN
                        CALL SGEMV('N',N,ki-2,ONE,Vr,Ldvr,Work(1+N),1,  &
     &                             Work(ki-1+N),Vr(1,ki-1),1)
                        CALL SGEMV('N',N,ki-2,ONE,Vr,Ldvr,Work(1+n2),1, &
     &                             Work(ki+n2),Vr(1,ki),1)
                     ELSE
                        CALL SSCAL(N,Work(ki-1+N),Vr(1,ki-1),1)
                        CALL SSCAL(N,Work(ki+n2),Vr(1,ki),1)
                     ENDIF
!
                     emax = ZERO
                     DO k = 1 , N
                        emax = MAX(emax,ABS(Vr(k,ki-1))+ABS(Vr(k,ki)))
                     ENDDO
                     remax = ONE/emax
                     CALL SSCAL(N,remax,Vr(1,ki-1),1)
                     CALL SSCAL(N,remax,Vr(1,ki),1)
                  ENDIF
               ENDIF
!
               is = is - 1
               IF ( ip/=0 ) is = is - 1
            ENDIF
 20         IF ( ip==1 ) ip = 0
            IF ( ip==-1 ) ip = 1
         ENDDO
      ENDIF
!
      IF ( leftv ) THEN
!
!        Compute left eigenvectors.
!
         ip = 0
         is = 1
         DO ki = 1 , N
!
            IF ( ip/=-1 ) THEN
               IF ( ki/=N ) THEN
                  IF ( T(ki+1,ki)/=ZERO ) ip = 1
               ENDIF
!
               IF ( somev ) THEN
                  IF ( .NOT.Select(ki) ) GOTO 40
               ENDIF
!
!           Compute the KI-th eigenvalue (WR,WI).
!
               wr = T(ki,ki)
               wi = ZERO
               IF ( ip/=0 ) wi = SQRT(ABS(T(ki,ki+1)))                  &
     &                           *SQRT(ABS(T(ki+1,ki)))
               smin = MAX(ulp*(ABS(wr)+ABS(wi)),smlnum)
!
               IF ( ip==0 ) THEN
!
!              Real left eigenvector.
!
                  Work(ki+N) = ONE
!
!              Form right-hand side
!
                  DO k = ki + 1 , N
                     Work(k+N) = -T(ki,k)
                  ENDDO
!
!              Solve the quasi-triangular system:
!                 (T(KI+1:N,KI+1:N) - WR)**T*X = SCALE*WORK
!
                  vmax = ONE
                  vcrit = bignum
!
                  jnxt = ki + 1
                  DO j = ki + 1 , N
                     IF ( j>=jnxt ) THEN
                        j1 = j
                        j2 = j
                        jnxt = j + 1
                        IF ( j<N ) THEN
                           IF ( T(j+1,j)/=ZERO ) THEN
                              j2 = j + 1
                              jnxt = j + 2
                           ENDIF
                        ENDIF
!
                        IF ( j1==j2 ) THEN
!
!                    1-by-1 diagonal block
!
!                    Scale if necessary to avoid overflow when forming
!                    the right-hand side.
!
                           IF ( Work(j)>vcrit ) THEN
                              rec = ONE/vmax
                              CALL SSCAL(N-ki+1,rec,Work(ki+N),1)
                              vmax = ONE
                              vcrit = bignum
                           ENDIF
!
                           Work(j+N) = Work(j+N)                        &
     &                                 - SDOT(j-ki-1,T(ki+1,j),1,       &
     &                                 Work(ki+1+N),1)
!
!                    Solve (T(J,J)-WR)**T*X = WORK
!
                           CALL SLALN2(.FALSE.,1,1,smin,ONE,T(j,j),Ldt, &
     &                                 ONE,ONE,Work(j+N),N,wr,ZERO,x,2, &
     &                                 scale,xnorm,ierr)
!
!                    Scale if necessary
!
                           IF ( scale/=ONE )                            &
     &                          CALL SSCAL(N-ki+1,scale,Work(ki+N),1)
                           Work(j+N) = x(1,1)
                           vmax = MAX(ABS(Work(j+N)),vmax)
                           vcrit = bignum/vmax
!
                        ELSE
!
!                    2-by-2 diagonal block
!
!                    Scale if necessary to avoid overflow when forming
!                    the right-hand side.
!
                           beta = MAX(Work(j),Work(j+1))
                           IF ( beta>vcrit ) THEN
                              rec = ONE/vmax
                              CALL SSCAL(N-ki+1,rec,Work(ki+N),1)
                              vmax = ONE
                              vcrit = bignum
                           ENDIF
!
                           Work(j+N) = Work(j+N)                        &
     &                                 - SDOT(j-ki-1,T(ki+1,j),1,       &
     &                                 Work(ki+1+N),1)
!
                           Work(j+1+N) = Work(j+1+N)                    &
     &                        - SDOT(j-ki-1,T(ki+1,j+1),1,Work(ki+1+N), &
     &                        1)
!
!                    Solve
!                      [T(J,J)-WR   T(J,J+1)     ]**T* X = SCALE*( WORK1 )
!                      [T(J+1,J)    T(J+1,J+1)-WR]               ( WORK2 )
!
                           CALL SLALN2(.TRUE.,2,1,smin,ONE,T(j,j),Ldt,  &
     &                                 ONE,ONE,Work(j+N),N,wr,ZERO,x,2, &
     &                                 scale,xnorm,ierr)
!
!                    Scale if necessary
!
                           IF ( scale/=ONE )                            &
     &                          CALL SSCAL(N-ki+1,scale,Work(ki+N),1)
                           Work(j+N) = x(1,1)
                           Work(j+1+N) = x(2,1)
!
                           vmax = MAX(ABS(Work(j+N)),ABS(Work(j+1+N)),  &
     &                            vmax)
                           vcrit = bignum/vmax
!
                        ENDIF
                     ENDIF
                  ENDDO
!
!              Copy the vector x or Q*x to VL and normalize.
!
                  IF ( .NOT.over ) THEN
                     CALL SCOPY(N-ki+1,Work(ki+N),1,Vl(ki,is),1)
!
                     ii = ISAMAX(N-ki+1,Vl(ki,is),1) + ki - 1
                     remax = ONE/ABS(Vl(ii,is))
                     CALL SSCAL(N-ki+1,remax,Vl(ki,is),1)
!
                     DO k = 1 , ki - 1
                        Vl(k,is) = ZERO
                     ENDDO
!
                  ELSE
!
                     IF ( ki<N ) CALL SGEMV('N',N,N-ki,ONE,Vl(1,ki+1),  &
     &                    Ldvl,Work(ki+1+N),1,Work(ki+N),Vl(1,ki),1)
!
                     ii = ISAMAX(N,Vl(1,ki),1)
                     remax = ONE/ABS(Vl(ii,ki))
                     CALL SSCAL(N,remax,Vl(1,ki),1)
!
                  ENDIF
!
               ELSE
!
!              Complex left eigenvector.
!
!               Initial solve:
!                 ((T(KI,KI)    T(KI,KI+1) )**T - (WR - I* WI))*X = 0.
!                 ((T(KI+1,KI) T(KI+1,KI+1))                )
!
                  IF ( ABS(T(ki,ki+1))>=ABS(T(ki+1,ki)) ) THEN
                     Work(ki+N) = wi/T(ki,ki+1)
                     Work(ki+1+n2) = ONE
                  ELSE
                     Work(ki+N) = ONE
                     Work(ki+1+n2) = -wi/T(ki+1,ki)
                  ENDIF
                  Work(ki+1+N) = ZERO
                  Work(ki+n2) = ZERO
!
!              Form right-hand side
!
                  DO k = ki + 2 , N
                     Work(k+N) = -Work(ki+N)*T(ki,k)
                     Work(k+n2) = -Work(ki+1+n2)*T(ki+1,k)
                  ENDDO
!
!              Solve complex quasi-triangular system:
!              ( T(KI+2,N:KI+2,N) - (WR-i*WI) )*X = WORK1+i*WORK2
!
                  vmax = ONE
                  vcrit = bignum
!
                  jnxt = ki + 2
                  DO j = ki + 2 , N
                     IF ( j>=jnxt ) THEN
                        j1 = j
                        j2 = j
                        jnxt = j + 1
                        IF ( j<N ) THEN
                           IF ( T(j+1,j)/=ZERO ) THEN
                              j2 = j + 1
                              jnxt = j + 2
                           ENDIF
                        ENDIF
!
                        IF ( j1==j2 ) THEN
!
!                    1-by-1 diagonal block
!
!                    Scale if necessary to avoid overflow when
!                    forming the right-hand side elements.
!
                           IF ( Work(j)>vcrit ) THEN
                              rec = ONE/vmax
                              CALL SSCAL(N-ki+1,rec,Work(ki+N),1)
                              CALL SSCAL(N-ki+1,rec,Work(ki+n2),1)
                              vmax = ONE
                              vcrit = bignum
                           ENDIF
!
                           Work(j+N) = Work(j+N)                        &
     &                                 - SDOT(j-ki-2,T(ki+2,j),1,       &
     &                                 Work(ki+2+N),1)
                           Work(j+n2) = Work(j+n2)                      &
     &                                  - SDOT(j-ki-2,T(ki+2,j),1,      &
     &                                  Work(ki+2+n2),1)
!
!                    Solve (T(J,J)-(WR-i*WI))*(X11+i*X12)= WK+I*WK2
!
                           CALL SLALN2(.FALSE.,1,2,smin,ONE,T(j,j),Ldt, &
     &                                 ONE,ONE,Work(j+N),N,wr,-wi,x,2,  &
     &                                 scale,xnorm,ierr)
!
!                    Scale if necessary
!
                           IF ( scale/=ONE ) THEN
                              CALL SSCAL(N-ki+1,scale,Work(ki+N),1)
                              CALL SSCAL(N-ki+1,scale,Work(ki+n2),1)
                           ENDIF
                           Work(j+N) = x(1,1)
                           Work(j+n2) = x(1,2)
                           vmax = MAX(ABS(Work(j+N)),ABS(Work(j+n2)),   &
     &                            vmax)
                           vcrit = bignum/vmax
!
                        ELSE
!
!                    2-by-2 diagonal block
!
!                    Scale if necessary to avoid overflow when forming
!                    the right-hand side elements.
!
                           beta = MAX(Work(j),Work(j+1))
                           IF ( beta>vcrit ) THEN
                              rec = ONE/vmax
                              CALL SSCAL(N-ki+1,rec,Work(ki+N),1)
                              CALL SSCAL(N-ki+1,rec,Work(ki+n2),1)
                              vmax = ONE
                              vcrit = bignum
                           ENDIF
!
                           Work(j+N) = Work(j+N)                        &
     &                                 - SDOT(j-ki-2,T(ki+2,j),1,       &
     &                                 Work(ki+2+N),1)
!
                           Work(j+n2) = Work(j+n2)                      &
     &                                  - SDOT(j-ki-2,T(ki+2,j),1,      &
     &                                  Work(ki+2+n2),1)
!
                           Work(j+1+N) = Work(j+1+N)                    &
     &                        - SDOT(j-ki-2,T(ki+2,j+1),1,Work(ki+2+N), &
     &                        1)
!
                           Work(j+1+n2) = Work(j+1+n2)                  &
     &                        - SDOT(j-ki-2,T(ki+2,j+1),1,Work(ki+2+n2),&
     &                        1)
!
!                    Solve 2-by-2 complex linear equation
!                      ([T(j,j)   T(j,j+1)  ]**T-(wr-i*wi)*I)*X = SCALE*B
!                      ([T(j+1,j) T(j+1,j+1)]               )
!
                           CALL SLALN2(.TRUE.,2,2,smin,ONE,T(j,j),Ldt,  &
     &                                 ONE,ONE,Work(j+N),N,wr,-wi,x,2,  &
     &                                 scale,xnorm,ierr)
!
!                    Scale if necessary
!
                           IF ( scale/=ONE ) THEN
                              CALL SSCAL(N-ki+1,scale,Work(ki+N),1)
                              CALL SSCAL(N-ki+1,scale,Work(ki+n2),1)
                           ENDIF
                           Work(j+N) = x(1,1)
                           Work(j+n2) = x(1,2)
                           Work(j+1+N) = x(2,1)
                           Work(j+1+n2) = x(2,2)
                           vmax = MAX(ABS(x(1,1)),ABS(x(1,2)),          &
     &                            ABS(x(2,1)),ABS(x(2,2)),vmax)
                           vcrit = bignum/vmax
!
                        ENDIF
                     ENDIF
                  ENDDO
!
!              Copy the vector x or Q*x to VL and normalize.
!
                  IF ( .NOT.over ) THEN
                     CALL SCOPY(N-ki+1,Work(ki+N),1,Vl(ki,is),1)
                     CALL SCOPY(N-ki+1,Work(ki+n2),1,Vl(ki,is+1),1)
!
                     emax = ZERO
                     DO k = ki , N
                        emax = MAX(emax,ABS(Vl(k,is))+ABS(Vl(k,is+1)))
                     ENDDO
                     remax = ONE/emax
                     CALL SSCAL(N-ki+1,remax,Vl(ki,is),1)
                     CALL SSCAL(N-ki+1,remax,Vl(ki,is+1),1)
!
                     DO k = 1 , ki - 1
                        Vl(k,is) = ZERO
                        Vl(k,is+1) = ZERO
                     ENDDO
                  ELSE
                     IF ( ki<N-1 ) THEN
                        CALL SGEMV('N',N,N-ki-1,ONE,Vl(1,ki+2),Ldvl,    &
     &                             Work(ki+2+N),1,Work(ki+N),Vl(1,ki),1)
                        CALL SGEMV('N',N,N-ki-1,ONE,Vl(1,ki+2),Ldvl,    &
     &                             Work(ki+2+n2),1,Work(ki+1+n2),       &
     &                             Vl(1,ki+1),1)
                     ELSE
                        CALL SSCAL(N,Work(ki+N),Vl(1,ki),1)
                        CALL SSCAL(N,Work(ki+1+n2),Vl(1,ki+1),1)
                     ENDIF
!
                     emax = ZERO
                     DO k = 1 , N
                        emax = MAX(emax,ABS(Vl(k,ki))+ABS(Vl(k,ki+1)))
                     ENDDO
                     remax = ONE/emax
                     CALL SSCAL(N,remax,Vl(1,ki),1)
                     CALL SSCAL(N,remax,Vl(1,ki+1),1)
!
                  ENDIF
!
               ENDIF
!
               is = is + 1
               IF ( ip/=0 ) is = is + 1
            ENDIF
 40         IF ( ip==-1 ) ip = 0
            IF ( ip==1 ) ip = -1
!
         ENDDO
!
      ENDIF
!
!
!     End of STREVC
!
      END SUBROUTINE STREVC
