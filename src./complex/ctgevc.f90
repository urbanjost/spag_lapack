!*==ctgevc.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b CTGEVC
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CTGEVC + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ctgevc.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ctgevc.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ctgevc.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CTGEVC( SIDE, HOWMNY, SELECT, N, S, LDS, P, LDP, VL,
!                          LDVL, VR, LDVR, MM, M, WORK, RWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          HOWMNY, SIDE
!       INTEGER            INFO, LDP, LDS, LDVL, LDVR, M, MM, N
!       ..
!       .. Array Arguments ..
!       LOGICAL            SELECT( * )
!       REAL               RWORK( * )
!       COMPLEX            P( LDP, * ), S( LDS, * ), VL( LDVL, * ),
!      $                   VR( LDVR, * ), WORK( * )
!       ..
!
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CTGEVC computes some or all of the right and/or left eigenvectors of
!> a pair of complex matrices (S,P), where S and P are upper triangular.
!> Matrix pairs of this type are produced by the generalized Schur
!> factorization of a complex matrix pair (A,B):
!>
!>    A = Q*S*Z**H,  B = Q*P*Z**H
!>
!> as computed by CGGHRD + CHGEQZ.
!>
!> The right eigenvector x and the left eigenvector y of (S,P)
!> corresponding to an eigenvalue w are defined by:
!>
!>    S*x = w*P*x,  (y**H)*S = w*(y**H)*P,
!>
!> where y**H denotes the conjugate tranpose of y.
!> The eigenvalues are not input to this routine, but are computed
!> directly from the diagonal elements of S and P.
!>
!> This routine returns the matrices X and/or Y of right and left
!> eigenvectors of (S,P), or the products Z*X and/or Q*Y,
!> where Z and Q are input matrices.
!> If Q and Z are the unitary factors from the generalized Schur
!> factorization of a matrix pair (A,B), then Z*X and Q*Y
!> are the matrices of right and left eigenvectors of (A,B).
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] SIDE
!> \verbatim
!>          SIDE is CHARACTER*1
!>          = 'R': compute right eigenvectors only;
!>          = 'L': compute left eigenvectors only;
!>          = 'B': compute both right and left eigenvectors.
!> \endverbatim
!>
!> \param[in] HOWMNY
!> \verbatim
!>          HOWMNY is CHARACTER*1
!>          = 'A': compute all right and/or left eigenvectors;
!>          = 'B': compute all right and/or left eigenvectors,
!>                 backtransformed by the matrices in VR and/or VL;
!>          = 'S': compute selected right and/or left eigenvectors,
!>                 specified by the logical array SELECT.
!> \endverbatim
!>
!> \param[in] SELECT
!> \verbatim
!>          SELECT is LOGICAL array, dimension (N)
!>          If HOWMNY='S', SELECT specifies the eigenvectors to be
!>          computed.  The eigenvector corresponding to the j-th
!>          eigenvalue is computed if SELECT(j) = .TRUE..
!>          Not referenced if HOWMNY = 'A' or 'B'.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrices S and P.  N >= 0.
!> \endverbatim
!>
!> \param[in] S
!> \verbatim
!>          S is COMPLEX array, dimension (LDS,N)
!>          The upper triangular matrix S from a generalized Schur
!>          factorization, as computed by CHGEQZ.
!> \endverbatim
!>
!> \param[in] LDS
!> \verbatim
!>          LDS is INTEGER
!>          The leading dimension of array S.  LDS >= max(1,N).
!> \endverbatim
!>
!> \param[in] P
!> \verbatim
!>          P is COMPLEX array, dimension (LDP,N)
!>          The upper triangular matrix P from a generalized Schur
!>          factorization, as computed by CHGEQZ.  P must have real
!>          diagonal elements.
!> \endverbatim
!>
!> \param[in] LDP
!> \verbatim
!>          LDP is INTEGER
!>          The leading dimension of array P.  LDP >= max(1,N).
!> \endverbatim
!>
!> \param[in,out] VL
!> \verbatim
!>          VL is COMPLEX array, dimension (LDVL,MM)
!>          On entry, if SIDE = 'L' or 'B' and HOWMNY = 'B', VL must
!>          contain an N-by-N matrix Q (usually the unitary matrix Q
!>          of left Schur vectors returned by CHGEQZ).
!>          On exit, if SIDE = 'L' or 'B', VL contains:
!>          if HOWMNY = 'A', the matrix Y of left eigenvectors of (S,P);
!>          if HOWMNY = 'B', the matrix Q*Y;
!>          if HOWMNY = 'S', the left eigenvectors of (S,P) specified by
!>                      SELECT, stored consecutively in the columns of
!>                      VL, in the same order as their eigenvalues.
!>          Not referenced if SIDE = 'R'.
!> \endverbatim
!>
!> \param[in] LDVL
!> \verbatim
!>          LDVL is INTEGER
!>          The leading dimension of array VL.  LDVL >= 1, and if
!>          SIDE = 'L' or 'l' or 'B' or 'b', LDVL >= N.
!> \endverbatim
!>
!> \param[in,out] VR
!> \verbatim
!>          VR is COMPLEX array, dimension (LDVR,MM)
!>          On entry, if SIDE = 'R' or 'B' and HOWMNY = 'B', VR must
!>          contain an N-by-N matrix Q (usually the unitary matrix Z
!>          of right Schur vectors returned by CHGEQZ).
!>          On exit, if SIDE = 'R' or 'B', VR contains:
!>          if HOWMNY = 'A', the matrix X of right eigenvectors of (S,P);
!>          if HOWMNY = 'B', the matrix Z*X;
!>          if HOWMNY = 'S', the right eigenvectors of (S,P) specified by
!>                      SELECT, stored consecutively in the columns of
!>                      VR, in the same order as their eigenvalues.
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
!>          used to store the eigenvectors.  If HOWMNY = 'A' or 'B', M
!>          is set to N.  Each selected eigenvector occupies one column.
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
!>          = 0:  successful exit.
!>          < 0:  if INFO = -i, the i-th argument had an illegal value.
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
      SUBROUTINE CTGEVC(Side,Howmny,Select,N,S,Lds,P,Ldp,Vl,Ldvl,Vr,    &
     &                  Ldvr,Mm,M,Work,Rwork,Info)
      USE S_CGEMV
      USE S_CLADIV
      USE S_LSAME
      USE S_SLABAD
      USE S_SLAMCH
      USE S_XERBLA
      IMPLICIT NONE
!*--CTGEVC229
!
! PARAMETER definitions rewritten by SPAG
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
      COMPLEX , PARAMETER  ::  CZERO = (0.0E+0,0.0E+0) ,                &
     &                         CONE = (1.0E+0,0.0E+0)
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Side
      CHARACTER :: Howmny
      LOGICAL , INTENT(IN) , DIMENSION(*) :: Select
      INTEGER :: N
      COMPLEX , INTENT(IN) , DIMENSION(Lds,*) :: S
      INTEGER , INTENT(IN) :: Lds
      COMPLEX , INTENT(IN) , DIMENSION(Ldp,*) :: P
      INTEGER , INTENT(IN) :: Ldp
      COMPLEX , DIMENSION(Ldvl,*) :: Vl
      INTEGER :: Ldvl
      COMPLEX , DIMENSION(Ldvr,*) :: Vr
      INTEGER :: Ldvr
      INTEGER , INTENT(IN) :: Mm
      INTEGER , INTENT(OUT) :: M
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Work
      REAL , INTENT(INOUT) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      REAL :: ABS1
      REAL :: acoefa , acoeff , anorm , ascale , bcoefa , big , bignum ,&
     &        bnorm , bscale , dmin , safmin , sbeta , scale , small ,  &
     &        temp , ulp , xmax
      COMPLEX :: bcoeff , ca , cb , d , salpha , sum , suma , sumb , x
      LOGICAL :: compl , compr , ilall , ilback , ilbbad , ilcomp ,     &
     &           lsa , lsb
      INTEGER :: i , ibeg , ieig , iend , ihwmny , im , iside , isrc ,  &
     &           j , je , jr
!
! End of declarations rewritten by SPAG
!
!     ..
!     .. Array Arguments ..
!     ..
!
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
      ABS1(x) = ABS(REAL(x)) + ABS(AIMAG(x))
!     ..
!     .. Executable Statements ..
!
!     Decode and Test the input parameters
!
      IF ( LSAME(Howmny,'A') ) THEN
         ihwmny = 1
         ilall = .TRUE.
         ilback = .FALSE.
      ELSEIF ( LSAME(Howmny,'S') ) THEN
         ihwmny = 2
         ilall = .FALSE.
         ilback = .FALSE.
      ELSEIF ( LSAME(Howmny,'B') ) THEN
         ihwmny = 3
         ilall = .TRUE.
         ilback = .TRUE.
      ELSE
         ihwmny = -1
      ENDIF
!
      IF ( LSAME(Side,'R') ) THEN
         iside = 1
         compl = .FALSE.
         compr = .TRUE.
      ELSEIF ( LSAME(Side,'L') ) THEN
         iside = 2
         compl = .TRUE.
         compr = .FALSE.
      ELSEIF ( LSAME(Side,'B') ) THEN
         iside = 3
         compl = .TRUE.
         compr = .TRUE.
      ELSE
         iside = -1
      ENDIF
!
      Info = 0
      IF ( iside<0 ) THEN
         Info = -1
      ELSEIF ( ihwmny<0 ) THEN
         Info = -2
      ELSEIF ( N<0 ) THEN
         Info = -4
      ELSEIF ( Lds<MAX(1,N) ) THEN
         Info = -6
      ELSEIF ( Ldp<MAX(1,N) ) THEN
         Info = -8
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('CTGEVC',-Info)
         RETURN
      ENDIF
!
!     Count the number of eigenvectors
!
      IF ( .NOT.ilall ) THEN
         im = 0
         DO j = 1 , N
            IF ( Select(j) ) im = im + 1
         ENDDO
      ELSE
         im = N
      ENDIF
!
!     Check diagonal of B
!
      ilbbad = .FALSE.
      DO j = 1 , N
         IF ( AIMAG(P(j,j))/=ZERO ) ilbbad = .TRUE.
      ENDDO
!
      IF ( ilbbad ) THEN
         Info = -7
      ELSEIF ( compl .AND. Ldvl<N .OR. Ldvl<1 ) THEN
         Info = -10
      ELSEIF ( compr .AND. Ldvr<N .OR. Ldvr<1 ) THEN
         Info = -12
      ELSEIF ( Mm<im ) THEN
         Info = -13
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('CTGEVC',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      M = im
      IF ( N==0 ) RETURN
!
!     Machine Constants
!
      safmin = SLAMCH('Safe minimum')
      big = ONE/safmin
      CALL SLABAD(safmin,big)
      ulp = SLAMCH('Epsilon')*SLAMCH('Base')
      small = safmin*N/ulp
      big = ONE/small
      bignum = ONE/(safmin*N)
!
!     Compute the 1-norm of each column of the strictly upper triangular
!     part of A and B to check for possible overflow in the triangular
!     solver.
!
      anorm = ABS1(S(1,1))
      bnorm = ABS1(P(1,1))
      Rwork(1) = ZERO
      Rwork(N+1) = ZERO
      DO j = 2 , N
         Rwork(j) = ZERO
         Rwork(N+j) = ZERO
         DO i = 1 , j - 1
            Rwork(j) = Rwork(j) + ABS1(S(i,j))
            Rwork(N+j) = Rwork(N+j) + ABS1(P(i,j))
         ENDDO
         anorm = MAX(anorm,Rwork(j)+ABS1(S(j,j)))
         bnorm = MAX(bnorm,Rwork(N+j)+ABS1(P(j,j)))
      ENDDO
!
      ascale = ONE/MAX(anorm,safmin)
      bscale = ONE/MAX(bnorm,safmin)
!
!     Left eigenvectors
!
      IF ( compl ) THEN
         ieig = 0
!
!        Main loop over eigenvalues
!
         DO je = 1 , N
            IF ( ilall ) THEN
               ilcomp = .TRUE.
            ELSE
               ilcomp = Select(je)
            ENDIF
            IF ( ilcomp ) THEN
               ieig = ieig + 1
!
               IF ( ABS1(S(je,je))<=safmin .AND. ABS(REAL(P(je,je)))    &
     &              <=safmin ) THEN
!
!                 Singular matrix pencil -- return unit eigenvector
!
                  DO jr = 1 , N
                     Vl(jr,ieig) = CZERO
                  ENDDO
                  Vl(ieig,ieig) = CONE
                  CYCLE
               ENDIF
!
!              Non-singular eigenvalue:
!              Compute coefficients  a  and  b  in
!                   H
!                 y  ( a A - b B ) = 0
!
               temp = ONE/MAX(ABS1(S(je,je))*ascale,ABS(REAL(P(je,je))) &
     &                *bscale,safmin)
               salpha = (temp*S(je,je))*ascale
               sbeta = (temp*REAL(P(je,je)))*bscale
               acoeff = sbeta*ascale
               bcoeff = salpha*bscale
!
!              Scale to avoid underflow
!
               lsa = ABS(sbeta)>=safmin .AND. ABS(acoeff)<small
               lsb = ABS1(salpha)>=safmin .AND. ABS1(bcoeff)<small
!
               scale = ONE
               IF ( lsa ) scale = (small/ABS(sbeta))*MIN(anorm,big)
               IF ( lsb ) scale = MAX(scale,(small/ABS1(salpha))*MIN(   &
     &                            bnorm,big))
               IF ( lsa .OR. lsb ) THEN
                  scale = MIN(scale,                                    &
     &                    ONE/(safmin*MAX(ONE,ABS(acoeff),ABS1(bcoeff)))&
     &                    )
                  IF ( lsa ) THEN
                     acoeff = ascale*(scale*sbeta)
                  ELSE
                     acoeff = scale*acoeff
                  ENDIF
                  IF ( lsb ) THEN
                     bcoeff = bscale*(scale*salpha)
                  ELSE
                     bcoeff = scale*bcoeff
                  ENDIF
               ENDIF
!
               acoefa = ABS(acoeff)
               bcoefa = ABS1(bcoeff)
               xmax = ONE
               DO jr = 1 , N
                  Work(jr) = CZERO
               ENDDO
               Work(je) = CONE
               dmin = MAX(ulp*acoefa*anorm,ulp*bcoefa*bnorm,safmin)
!
!                                              H
!              Triangular solve of  (a A - b B)  y = 0
!
!                                      H
!              (rowwise in  (a A - b B) , or columnwise in a A - b B)
!
               DO j = je + 1 , N
!
!                 Compute
!                       j-1
!                 SUM = sum  conjg( a*S(k,j) - b*P(k,j) )*x(k)
!                       k=je
!                 (Scale if necessary)
!
                  temp = ONE/xmax
                  IF ( acoefa*Rwork(j)+bcoefa*Rwork(N+j)>bignum*temp )  &
     &                 THEN
                     DO jr = je , j - 1
                        Work(jr) = temp*Work(jr)
                     ENDDO
                     xmax = ONE
                  ENDIF
                  suma = CZERO
                  sumb = CZERO
!
                  DO jr = je , j - 1
                     suma = suma + CONJG(S(jr,j))*Work(jr)
                     sumb = sumb + CONJG(P(jr,j))*Work(jr)
                  ENDDO
                  sum = acoeff*suma - CONJG(bcoeff)*sumb
!
!                 Form x(j) = - SUM / conjg( a*S(j,j) - b*P(j,j) )
!
!                 with scaling and perturbation of the denominator
!
                  d = CONJG(acoeff*S(j,j)-bcoeff*P(j,j))
                  IF ( ABS1(d)<=dmin ) d = CMPLX(dmin)
!
                  IF ( ABS1(d)<ONE ) THEN
                     IF ( ABS1(sum)>=bignum*ABS1(d) ) THEN
                        temp = ONE/ABS1(sum)
                        DO jr = je , j - 1
                           Work(jr) = temp*Work(jr)
                        ENDDO
                        xmax = temp*xmax
                        sum = temp*sum
                     ENDIF
                  ENDIF
                  Work(j) = CLADIV(-sum,d)
                  xmax = MAX(xmax,ABS1(Work(j)))
               ENDDO
!
!              Back transform eigenvector if HOWMNY='B'.
!
               IF ( ilback ) THEN
                  CALL CGEMV('N',N,N+1-je,CONE,Vl(1,je),Ldvl,Work(je),1,&
     &                       CZERO,Work(N+1),1)
                  isrc = 2
                  ibeg = 1
               ELSE
                  isrc = 1
                  ibeg = je
               ENDIF
!
!              Copy and scale eigenvector into column of VL
!
               xmax = ZERO
               DO jr = ibeg , N
                  xmax = MAX(xmax,ABS1(Work((isrc-1)*N+jr)))
               ENDDO
!
               IF ( xmax>safmin ) THEN
                  temp = ONE/xmax
                  DO jr = ibeg , N
                     Vl(jr,ieig) = temp*Work((isrc-1)*N+jr)
                  ENDDO
               ELSE
                  ibeg = N + 1
               ENDIF
!
               DO jr = 1 , ibeg - 1
                  Vl(jr,ieig) = CZERO
               ENDDO
!
            ENDIF
         ENDDO
      ENDIF
!
!     Right eigenvectors
!
      IF ( compr ) THEN
         ieig = im + 1
!
!        Main loop over eigenvalues
!
         DO je = N , 1 , -1
            IF ( ilall ) THEN
               ilcomp = .TRUE.
            ELSE
               ilcomp = Select(je)
            ENDIF
            IF ( ilcomp ) THEN
               ieig = ieig - 1
!
               IF ( ABS1(S(je,je))<=safmin .AND. ABS(REAL(P(je,je)))    &
     &              <=safmin ) THEN
!
!                 Singular matrix pencil -- return unit eigenvector
!
                  DO jr = 1 , N
                     Vr(jr,ieig) = CZERO
                  ENDDO
                  Vr(ieig,ieig) = CONE
                  CYCLE
               ENDIF
!
!              Non-singular eigenvalue:
!              Compute coefficients  a  and  b  in
!
!              ( a A - b B ) x  = 0
!
               temp = ONE/MAX(ABS1(S(je,je))*ascale,ABS(REAL(P(je,je))) &
     &                *bscale,safmin)
               salpha = (temp*S(je,je))*ascale
               sbeta = (temp*REAL(P(je,je)))*bscale
               acoeff = sbeta*ascale
               bcoeff = salpha*bscale
!
!              Scale to avoid underflow
!
               lsa = ABS(sbeta)>=safmin .AND. ABS(acoeff)<small
               lsb = ABS1(salpha)>=safmin .AND. ABS1(bcoeff)<small
!
               scale = ONE
               IF ( lsa ) scale = (small/ABS(sbeta))*MIN(anorm,big)
               IF ( lsb ) scale = MAX(scale,(small/ABS1(salpha))*MIN(   &
     &                            bnorm,big))
               IF ( lsa .OR. lsb ) THEN
                  scale = MIN(scale,                                    &
     &                    ONE/(safmin*MAX(ONE,ABS(acoeff),ABS1(bcoeff)))&
     &                    )
                  IF ( lsa ) THEN
                     acoeff = ascale*(scale*sbeta)
                  ELSE
                     acoeff = scale*acoeff
                  ENDIF
                  IF ( lsb ) THEN
                     bcoeff = bscale*(scale*salpha)
                  ELSE
                     bcoeff = scale*bcoeff
                  ENDIF
               ENDIF
!
               acoefa = ABS(acoeff)
               bcoefa = ABS1(bcoeff)
               xmax = ONE
               DO jr = 1 , N
                  Work(jr) = CZERO
               ENDDO
               Work(je) = CONE
               dmin = MAX(ulp*acoefa*anorm,ulp*bcoefa*bnorm,safmin)
!
!              Triangular solve of  (a A - b B) x = 0  (columnwise)
!
!              WORK(1:j-1) contains sums w,
!              WORK(j+1:JE) contains x
!
               DO jr = 1 , je - 1
                  Work(jr) = acoeff*S(jr,je) - bcoeff*P(jr,je)
               ENDDO
               Work(je) = CONE
!
               DO j = je - 1 , 1 , -1
!
!                 Form x(j) := - w(j) / d
!                 with scaling and perturbation of the denominator
!
                  d = acoeff*S(j,j) - bcoeff*P(j,j)
                  IF ( ABS1(d)<=dmin ) d = CMPLX(dmin)
!
                  IF ( ABS1(d)<ONE ) THEN
                     IF ( ABS1(Work(j))>=bignum*ABS1(d) ) THEN
                        temp = ONE/ABS1(Work(j))
                        DO jr = 1 , je
                           Work(jr) = temp*Work(jr)
                        ENDDO
                     ENDIF
                  ENDIF
!
                  Work(j) = CLADIV(-Work(j),d)
!
                  IF ( j>1 ) THEN
!
!                    w = w + x(j)*(a S(*,j) - b P(*,j) ) with scaling
!
                     IF ( ABS1(Work(j))>ONE ) THEN
                        temp = ONE/ABS1(Work(j))
                        IF ( acoefa*Rwork(j)+bcoefa*Rwork(N+j)          &
     &                       >=bignum*temp ) THEN
                           DO jr = 1 , je
                              Work(jr) = temp*Work(jr)
                           ENDDO
                        ENDIF
                     ENDIF
!
                     ca = acoeff*Work(j)
                     cb = bcoeff*Work(j)
                     DO jr = 1 , j - 1
                        Work(jr) = Work(jr) + ca*S(jr,j) - cb*P(jr,j)
                     ENDDO
                  ENDIF
               ENDDO
!
!              Back transform eigenvector if HOWMNY='B'.
!
               IF ( ilback ) THEN
                  CALL CGEMV('N',N,je,CONE,Vr,Ldvr,Work,1,CZERO,        &
     &                       Work(N+1),1)
                  isrc = 2
                  iend = N
               ELSE
                  isrc = 1
                  iend = je
               ENDIF
!
!              Copy and scale eigenvector into column of VR
!
               xmax = ZERO
               DO jr = 1 , iend
                  xmax = MAX(xmax,ABS1(Work((isrc-1)*N+jr)))
               ENDDO
!
               IF ( xmax>safmin ) THEN
                  temp = ONE/xmax
                  DO jr = 1 , iend
                     Vr(jr,ieig) = temp*Work((isrc-1)*N+jr)
                  ENDDO
               ELSE
                  iend = 0
               ENDIF
!
               DO jr = iend + 1 , N
                  Vr(jr,ieig) = CZERO
               ENDDO
!
            ENDIF
         ENDDO
      ENDIF
!
!
!     End of CTGEVC
!
      END SUBROUTINE CTGEVC
