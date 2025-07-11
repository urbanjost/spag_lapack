!*==stgevc.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b STGEVC
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download STGEVC + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/stgevc.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/stgevc.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/stgevc.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE STGEVC( SIDE, HOWMNY, SELECT, N, S, LDS, P, LDP, VL,
!                          LDVL, VR, LDVR, MM, M, WORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          HOWMNY, SIDE
!       INTEGER            INFO, LDP, LDS, LDVL, LDVR, M, MM, N
!       ..
!       .. Array Arguments ..
!       LOGICAL            SELECT( * )
!       REAL               P( LDP, * ), S( LDS, * ), VL( LDVL, * ),
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
!> STGEVC computes some or all of the right and/or left eigenvectors of
!> a pair of real matrices (S,P), where S is a quasi-triangular matrix
!> and P is upper triangular.  Matrix pairs of this type are produced by
!> the generalized Schur factorization of a matrix pair (A,B):
!>
!>    A = Q*S*Z**T,  B = Q*P*Z**T
!>
!> as computed by SGGHRD + SHGEQZ.
!>
!> The right eigenvector x and the left eigenvector y of (S,P)
!> corresponding to an eigenvalue w are defined by:
!>
!>    S*x = w*P*x,  (y**H)*S = w*(y**H)*P,
!>
!> where y**H denotes the conjugate tranpose of y.
!> The eigenvalues are not input to this routine, but are computed
!> directly from the diagonal blocks of S and P.
!>
!> This routine returns the matrices X and/or Y of right and left
!> eigenvectors of (S,P), or the products Z*X and/or Q*Y,
!> where Z and Q are input matrices.
!> If Q and Z are the orthogonal factors from the generalized Schur
!> factorization of a matrix pair (A,B), then Z*X and Q*Y
!> are the matrices of right and left eigenvectors of (A,B).
!>
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
!>          computed.  If w(j) is a real eigenvalue, the corresponding
!>          real eigenvector is computed if SELECT(j) is .TRUE..
!>          If w(j) and w(j+1) are the real and imaginary parts of a
!>          complex eigenvalue, the corresponding complex eigenvector
!>          is computed if either SELECT(j) or SELECT(j+1) is .TRUE.,
!>          and on exit SELECT(j) is set to .TRUE. and SELECT(j+1) is
!>          set to .FALSE..
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
!>          S is REAL array, dimension (LDS,N)
!>          The upper quasi-triangular matrix S from a generalized Schur
!>          factorization, as computed by SHGEQZ.
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
!>          P is REAL array, dimension (LDP,N)
!>          The upper triangular matrix P from a generalized Schur
!>          factorization, as computed by SHGEQZ.
!>          2-by-2 diagonal blocks of P corresponding to 2-by-2 blocks
!>          of S must be in positive diagonal form.
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
!>          VL is REAL array, dimension (LDVL,MM)
!>          On entry, if SIDE = 'L' or 'B' and HOWMNY = 'B', VL must
!>          contain an N-by-N matrix Q (usually the orthogonal matrix Q
!>          of left Schur vectors returned by SHGEQZ).
!>          On exit, if SIDE = 'L' or 'B', VL contains:
!>          if HOWMNY = 'A', the matrix Y of left eigenvectors of (S,P);
!>          if HOWMNY = 'B', the matrix Q*Y;
!>          if HOWMNY = 'S', the left eigenvectors of (S,P) specified by
!>                      SELECT, stored consecutively in the columns of
!>                      VL, in the same order as their eigenvalues.
!>
!>          A complex eigenvector corresponding to a complex eigenvalue
!>          is stored in two consecutive columns, the first holding the
!>          real part, and the second the imaginary part.
!>
!>          Not referenced if SIDE = 'R'.
!> \endverbatim
!>
!> \param[in] LDVL
!> \verbatim
!>          LDVL is INTEGER
!>          The leading dimension of array VL.  LDVL >= 1, and if
!>          SIDE = 'L' or 'B', LDVL >= N.
!> \endverbatim
!>
!> \param[in,out] VR
!> \verbatim
!>          VR is REAL array, dimension (LDVR,MM)
!>          On entry, if SIDE = 'R' or 'B' and HOWMNY = 'B', VR must
!>          contain an N-by-N matrix Z (usually the orthogonal matrix Z
!>          of right Schur vectors returned by SHGEQZ).
!>
!>          On exit, if SIDE = 'R' or 'B', VR contains:
!>          if HOWMNY = 'A', the matrix X of right eigenvectors of (S,P);
!>          if HOWMNY = 'B' or 'b', the matrix Z*X;
!>          if HOWMNY = 'S' or 's', the right eigenvectors of (S,P)
!>                      specified by SELECT, stored consecutively in the
!>                      columns of VR, in the same order as their
!>                      eigenvalues.
!>
!>          A complex eigenvector corresponding to a complex eigenvalue
!>          is stored in two consecutive columns, the first holding the
!>          real part and the second the imaginary part.
!>
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
!>          is set to N.  Each selected real eigenvector occupies one
!>          column and each selected complex eigenvector occupies two
!>          columns.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is REAL array, dimension (6*N)
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit.
!>          < 0:  if INFO = -i, the i-th argument had an illegal value.
!>          > 0:  the 2-by-2 block (INFO:INFO+1) does not have a complex
!>                eigenvalue.
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
!> \ingroup realGEcomputational
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  Allocation of workspace:
!>  ---------- -- ---------
!>
!>     WORK( j ) = 1-norm of j-th column of A, above the diagonal
!>     WORK( N+j ) = 1-norm of j-th column of B, above the diagonal
!>     WORK( 2*N+1:3*N ) = real part of eigenvector
!>     WORK( 3*N+1:4*N ) = imaginary part of eigenvector
!>     WORK( 4*N+1:5*N ) = real part of back-transformed eigenvector
!>     WORK( 5*N+1:6*N ) = imaginary part of back-transformed eigenvector
!>
!>  Rowwise vs. columnwise solution methods:
!>  ------- --  ---------- -------- -------
!>
!>  Finding a generalized eigenvector consists basically of solving the
!>  singular triangular system
!>
!>   (A - w B) x = 0     (for right) or:   (A - w B)**H y = 0  (for left)
!>
!>  Consider finding the i-th right eigenvector (assume all eigenvalues
!>  are real). The equation to be solved is:
!>       n                   i
!>  0 = sum  C(j,k) v(k)  = sum  C(j,k) v(k)     for j = i,. . .,1
!>      k=j                 k=j
!>
!>  where  C = (A - w B)  (The components v(i+1:n) are 0.)
!>
!>  The "rowwise" method is:
!>
!>  (1)  v(i) := 1
!>  for j = i-1,. . .,1:
!>                          i
!>      (2) compute  s = - sum C(j,k) v(k)   and
!>                        k=j+1
!>
!>      (3) v(j) := s / C(j,j)
!>
!>  Step 2 is sometimes called the "dot product" step, since it is an
!>  inner product between the j-th row and the portion of the eigenvector
!>  that has been computed so far.
!>
!>  The "columnwise" method consists basically in doing the sums
!>  for all the rows in parallel.  As each v(j) is computed, the
!>  contribution of v(j) times the j-th column of C is added to the
!>  partial sums.  Since FORTRAN arrays are stored columnwise, this has
!>  the advantage that at each step, the elements of C that are accessed
!>  are adjacent to one another, whereas with the rowwise method, the
!>  elements accessed at a step are spaced LDS (and LDP) words apart.
!>
!>  When finding left eigenvectors, the matrix in question is the
!>  transpose of the one in storage, so the rowwise method then
!>  actually accesses columns of A and B at each step, and so is the
!>  preferred method.
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE STGEVC(Side,Howmny,Select,N,S,Lds,P,Ldp,Vl,Ldvl,Vr,    &
     &                  Ldvr,Mm,M,Work,Info)
      IMPLICIT NONE
!*--STGEVC299
!
!  -- LAPACK computational routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Howmny , Side
      INTEGER Info , Ldp , Lds , Ldvl , Ldvr , M , Mm , N
!     ..
!     .. Array Arguments ..
      LOGICAL Select(*)
      REAL P(Ldp,*) , S(Lds,*) , Vl(Ldvl,*) , Vr(Ldvr,*) , Work(*)
!     ..
!
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ZERO , ONE , SAFETY
      PARAMETER (ZERO=0.0E+0,ONE=1.0E+0,SAFETY=1.0E+2)
!     ..
!     .. Local Scalars ..
      LOGICAL compl , compr , il2by2 , ilabad , ilall , ilback ,        &
     &        ilbbad , ilcomp , ilcplx , lsa , lsb
      INTEGER i , ibeg , ieig , iend , ihwmny , iinfo , im , iside , j ,&
     &        ja , jc , je , jr , jw , na , nw
      REAL acoef , acoefa , anorm , ascale , bcoefa , bcoefi , bcoefr , &
     &     big , bignum , bnorm , bscale , cim2a , cim2b , cimaga ,     &
     &     cimagb , cre2a , cre2b , creala , crealb , dmin , safmin ,   &
     &     salfar , sbeta , scale , small , temp , temp2 , temp2i ,     &
     &     temp2r , ulp , xmax , xscale
!     ..
!     .. Local Arrays ..
      REAL bdiag(2) , sum(2,2) , sums(2,2) , sump(2,2)
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      REAL SLAMCH
      EXTERNAL LSAME , SLAMCH
!     ..
!     .. External Subroutines ..
      EXTERNAL SGEMV , SLABAD , SLACPY , SLAG2 , SLALN2 , XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , MAX , MIN
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
         ilall = .TRUE.
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
         CALL XERBLA('STGEVC',-Info)
         RETURN
      ENDIF
!
!     Count the number of eigenvectors to be computed
!
      IF ( .NOT.ilall ) THEN
         im = 0
         ilcplx = .FALSE.
         DO j = 1 , N
            IF ( ilcplx ) THEN
               ilcplx = .FALSE.
               CYCLE
            ENDIF
            IF ( j<N ) THEN
               IF ( S(j+1,j)/=ZERO ) ilcplx = .TRUE.
            ENDIF
            IF ( ilcplx ) THEN
               IF ( Select(j) .OR. Select(j+1) ) im = im + 2
            ELSE
               IF ( Select(j) ) im = im + 1
            ENDIF
         ENDDO
      ELSE
         im = N
      ENDIF
!
!     Check 2-by-2 diagonal blocks of A, B
!
      ilabad = .FALSE.
      ilbbad = .FALSE.
      DO j = 1 , N - 1
         IF ( S(j+1,j)/=ZERO ) THEN
            IF ( P(j,j)==ZERO .OR. P(j+1,j+1)==ZERO .OR. P(j,j+1)       &
     &           /=ZERO ) ilbbad = .TRUE.
            IF ( j<N-1 ) THEN
               IF ( S(j+2,j+1)/=ZERO ) ilabad = .TRUE.
            ENDIF
         ENDIF
      ENDDO
!
      IF ( ilabad ) THEN
         Info = -5
      ELSEIF ( ilbbad ) THEN
         Info = -7
      ELSEIF ( compl .AND. Ldvl<N .OR. Ldvl<1 ) THEN
         Info = -10
      ELSEIF ( compr .AND. Ldvr<N .OR. Ldvr<1 ) THEN
         Info = -12
      ELSEIF ( Mm<im ) THEN
         Info = -13
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('STGEVC',-Info)
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
!     part (i.e., excluding all elements belonging to the diagonal
!     blocks) of A and B to check for possible overflow in the
!     triangular solver.
!
      anorm = ABS(S(1,1))
      IF ( N>1 ) anorm = anorm + ABS(S(2,1))
      bnorm = ABS(P(1,1))
      Work(1) = ZERO
      Work(N+1) = ZERO
!
      DO j = 2 , N
         temp = ZERO
         temp2 = ZERO
         IF ( S(j,j-1)==ZERO ) THEN
            iend = j - 1
         ELSE
            iend = j - 2
         ENDIF
         DO i = 1 , iend
            temp = temp + ABS(S(i,j))
            temp2 = temp2 + ABS(P(i,j))
         ENDDO
         Work(j) = temp
         Work(N+j) = temp2
         DO i = iend + 1 , MIN(j+1,N)
            temp = temp + ABS(S(i,j))
            temp2 = temp2 + ABS(P(i,j))
         ENDDO
         anorm = MAX(anorm,temp)
         bnorm = MAX(bnorm,temp2)
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
         ilcplx = .FALSE.
         DO je = 1 , N
!
!           Skip this iteration if (a) HOWMNY='S' and SELECT=.FALSE., or
!           (b) this would be the second of a complex pair.
!           Check for complex eigenvalue, so as to be sure of which
!           entry(-ies) of SELECT to look at.
!
            IF ( ilcplx ) THEN
               ilcplx = .FALSE.
               CYCLE
            ENDIF
            nw = 1
            IF ( je<N ) THEN
               IF ( S(je+1,je)/=ZERO ) THEN
                  ilcplx = .TRUE.
                  nw = 2
               ENDIF
            ENDIF
            IF ( ilall ) THEN
               ilcomp = .TRUE.
            ELSEIF ( ilcplx ) THEN
               ilcomp = Select(je) .OR. Select(je+1)
            ELSE
               ilcomp = Select(je)
            ENDIF
            IF ( ilcomp ) THEN
!
!           Decide if (a) singular pencil, (b) real eigenvalue, or
!           (c) complex eigenvalue.
!
               IF ( .NOT.ilcplx ) THEN
                  IF ( ABS(S(je,je))<=safmin .AND. ABS(P(je,je))        &
     &                 <=safmin ) THEN
!
!                 Singular matrix pencil -- return unit eigenvector
!
                     ieig = ieig + 1
                     DO jr = 1 , N
                        Vl(jr,ieig) = ZERO
                     ENDDO
                     Vl(ieig,ieig) = ONE
                     CYCLE
                  ENDIF
               ENDIF
!
!           Clear vector
!
               DO jr = 1 , nw*N
                  Work(2*N+jr) = ZERO
               ENDDO
!                                                 T
!           Compute coefficients in  ( a A - b B )  y = 0
!              a  is  ACOEF
!              b  is  BCOEFR + i*BCOEFI
!
               IF ( .NOT.ilcplx ) THEN
!
!              Real eigenvalue
!
                  temp = ONE/MAX(ABS(S(je,je))*ascale,ABS(P(je,je))     &
     &                   *bscale,safmin)
                  salfar = (temp*S(je,je))*ascale
                  sbeta = (temp*P(je,je))*bscale
                  acoef = sbeta*ascale
                  bcoefr = salfar*bscale
                  bcoefi = ZERO
!
!              Scale to avoid underflow
!
                  scale = ONE
                  lsa = ABS(sbeta)>=safmin .AND. ABS(acoef)<small
                  lsb = ABS(salfar)>=safmin .AND. ABS(bcoefr)<small
                  IF ( lsa ) scale = (small/ABS(sbeta))*MIN(anorm,big)
                  IF ( lsb ) scale = MAX(scale,(small/ABS(salfar))*MIN( &
     &                               bnorm,big))
                  IF ( lsa .OR. lsb ) THEN
                     scale = MIN(scale,                                 &
     &                       ONE/(safmin*MAX(ONE,ABS(acoef),ABS(bcoefr))&
     &                       ))
                     IF ( lsa ) THEN
                        acoef = ascale*(scale*sbeta)
                     ELSE
                        acoef = scale*acoef
                     ENDIF
                     IF ( lsb ) THEN
                        bcoefr = bscale*(scale*salfar)
                     ELSE
                        bcoefr = scale*bcoefr
                     ENDIF
                  ENDIF
                  acoefa = ABS(acoef)
                  bcoefa = ABS(bcoefr)
!
!              First component is 1
!
                  Work(2*N+je) = ONE
                  xmax = ONE
               ELSE
!
!              Complex eigenvalue
!
                  CALL SLAG2(S(je,je),Lds,P(je,je),Ldp,safmin*SAFETY,   &
     &                       acoef,temp,bcoefr,temp2,bcoefi)
                  bcoefi = -bcoefi
                  IF ( bcoefi==ZERO ) THEN
                     Info = je
                     RETURN
                  ENDIF
!
!              Scale to avoid over/underflow
!
                  acoefa = ABS(acoef)
                  bcoefa = ABS(bcoefr) + ABS(bcoefi)
                  scale = ONE
                  IF ( acoefa*ulp<safmin .AND. acoefa>=safmin )         &
     &                 scale = (safmin/ulp)/acoefa
                  IF ( bcoefa*ulp<safmin .AND. bcoefa>=safmin )         &
     &                 scale = MAX(scale,(safmin/ulp)/bcoefa)
                  IF ( safmin*acoefa>ascale )                           &
     &                 scale = ascale/(safmin*acoefa)
                  IF ( safmin*bcoefa>bscale )                           &
     &                 scale = MIN(scale,bscale/(safmin*bcoefa))
                  IF ( scale/=ONE ) THEN
                     acoef = scale*acoef
                     acoefa = ABS(acoef)
                     bcoefr = scale*bcoefr
                     bcoefi = scale*bcoefi
                     bcoefa = ABS(bcoefr) + ABS(bcoefi)
                  ENDIF
!
!              Compute first two components of eigenvector
!
                  temp = acoef*S(je+1,je)
                  temp2r = acoef*S(je,je) - bcoefr*P(je,je)
                  temp2i = -bcoefi*P(je,je)
                  IF ( ABS(temp)>ABS(temp2r)+ABS(temp2i) ) THEN
                     Work(2*N+je) = ONE
                     Work(3*N+je) = ZERO
                     Work(2*N+je+1) = -temp2r/temp
                     Work(3*N+je+1) = -temp2i/temp
                  ELSE
                     Work(2*N+je+1) = ONE
                     Work(3*N+je+1) = ZERO
                     temp = acoef*S(je,je+1)
                     Work(2*N+je) = (bcoefr*P(je+1,je+1)-acoef*S(je+1,je&
     &                              +1))/temp
                     Work(3*N+je) = bcoefi*P(je+1,je+1)/temp
                  ENDIF
                  xmax = MAX(ABS(Work(2*N+je))+ABS(Work(3*N+je)),       &
     &                   ABS(Work(2*N+je+1))+ABS(Work(3*N+je+1)))
               ENDIF
!
               dmin = MAX(ulp*acoefa*anorm,ulp*bcoefa*bnorm,safmin)
!
!                                           T
!           Triangular solve of  (a A - b B)  y = 0
!
!                                   T
!           (rowwise in  (a A - b B) , or columnwise in (a A - b B) )
!
               il2by2 = .FALSE.
!
               DO j = je + nw , N
                  IF ( il2by2 ) THEN
                     il2by2 = .FALSE.
                     CYCLE
                  ENDIF
!
                  na = 1
                  bdiag(1) = P(j,j)
                  IF ( j<N ) THEN
                     IF ( S(j+1,j)/=ZERO ) THEN
                        il2by2 = .TRUE.
                        bdiag(2) = P(j+1,j+1)
                        na = 2
                     ENDIF
                  ENDIF
!
!              Check whether scaling is necessary for dot products
!
                  xscale = ONE/MAX(ONE,xmax)
                  temp = MAX(Work(j),Work(N+j),acoefa*Work(j)           &
     &                   +bcoefa*Work(N+j))
                  IF ( il2by2 ) temp = MAX(temp,Work(j+1),Work(N+j+1),  &
     &                                 acoefa*Work(j+1)                 &
     &                                 +bcoefa*Work(N+j+1))
                  IF ( temp>bignum*xscale ) THEN
                     DO jw = 0 , nw - 1
                        DO jr = je , j - 1
                           Work((jw+2)*N+jr) = xscale*Work((jw+2)*N+jr)
                        ENDDO
                     ENDDO
                     xmax = xmax*xscale
                  ENDIF
!
!              Compute dot products
!
!                    j-1
!              SUM = sum  conjg( a*S(k,j) - b*P(k,j) )*x(k)
!                    k=je
!
!              To reduce the op count, this is done as
!
!              _        j-1                  _        j-1
!              a*conjg( sum  S(k,j)*x(k) ) - b*conjg( sum  P(k,j)*x(k) )
!                       k=je                          k=je
!
!              which may cause underflow problems if A or B are close
!              to underflow.  (E.g., less than SMALL.)
!
!
                  DO jw = 1 , nw
                     DO ja = 1 , na
                        sums(ja,jw) = ZERO
                        sump(ja,jw) = ZERO
!
                        DO jr = je , j - 1
                           sums(ja,jw) = sums(ja,jw) + S(jr,j+ja-1)     &
     &                        *Work((jw+1)*N+jr)
                           sump(ja,jw) = sump(ja,jw) + P(jr,j+ja-1)     &
     &                        *Work((jw+1)*N+jr)
                        ENDDO
                     ENDDO
                  ENDDO
!
                  DO ja = 1 , na
                     IF ( ilcplx ) THEN
                        sum(ja,1) = -acoef*sums(ja,1)                   &
     &                              + bcoefr*sump(ja,1)                 &
     &                              - bcoefi*sump(ja,2)
                        sum(ja,2) = -acoef*sums(ja,2)                   &
     &                              + bcoefr*sump(ja,2)                 &
     &                              + bcoefi*sump(ja,1)
                     ELSE
                        sum(ja,1) = -acoef*sums(ja,1)                   &
     &                              + bcoefr*sump(ja,1)
                     ENDIF
                  ENDDO
!
!                                  T
!              Solve  ( a A - b B )  y = SUM(,)
!              with scaling and perturbation of the denominator
!
                  CALL SLALN2(.TRUE.,na,nw,dmin,acoef,S(j,j),Lds,       &
     &                        bdiag(1),bdiag(2),sum,2,bcoefr,bcoefi,    &
     &                        Work(2*N+j),N,scale,temp,iinfo)
                  IF ( scale<ONE ) THEN
                     DO jw = 0 , nw - 1
                        DO jr = je , j - 1
                           Work((jw+2)*N+jr) = scale*Work((jw+2)*N+jr)
                        ENDDO
                     ENDDO
                     xmax = scale*xmax
                  ENDIF
                  xmax = MAX(xmax,temp)
               ENDDO
!
!           Copy eigenvector to VL, back transforming if
!           HOWMNY='B'.
!
               ieig = ieig + 1
               IF ( ilback ) THEN
                  DO jw = 0 , nw - 1
                     CALL SGEMV('N',N,N+1-je,ONE,Vl(1,je),Ldvl,         &
     &                          Work((jw+2)*N+je),1,ZERO,               &
     &                          Work((jw+4)*N+1),1)
                  ENDDO
                  CALL SLACPY(' ',N,nw,Work(4*N+1),N,Vl(1,je),Ldvl)
                  ibeg = 1
               ELSE
                  CALL SLACPY(' ',N,nw,Work(2*N+1),N,Vl(1,ieig),Ldvl)
                  ibeg = je
               ENDIF
!
!           Scale eigenvector
!
               xmax = ZERO
               IF ( ilcplx ) THEN
                  DO j = ibeg , N
                     xmax = MAX(xmax,ABS(Vl(j,ieig))+ABS(Vl(j,ieig+1)))
                  ENDDO
               ELSE
                  DO j = ibeg , N
                     xmax = MAX(xmax,ABS(Vl(j,ieig)))
                  ENDDO
               ENDIF
!
               IF ( xmax>safmin ) THEN
                  xscale = ONE/xmax
!
                  DO jw = 0 , nw - 1
                     DO jr = ibeg , N
                        Vl(jr,ieig+jw) = xscale*Vl(jr,ieig+jw)
                     ENDDO
                  ENDDO
               ENDIF
               ieig = ieig + nw - 1
            ENDIF
!
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
         ilcplx = .FALSE.
         DO je = N , 1 , -1
!
!           Skip this iteration if (a) HOWMNY='S' and SELECT=.FALSE., or
!           (b) this would be the second of a complex pair.
!           Check for complex eigenvalue, so as to be sure of which
!           entry(-ies) of SELECT to look at -- if complex, SELECT(JE)
!           or SELECT(JE-1).
!           If this is a complex pair, the 2-by-2 diagonal block
!           corresponding to the eigenvalue is in rows/columns JE-1:JE
!
            IF ( ilcplx ) THEN
               ilcplx = .FALSE.
               CYCLE
            ENDIF
            nw = 1
            IF ( je>1 ) THEN
               IF ( S(je,je-1)/=ZERO ) THEN
                  ilcplx = .TRUE.
                  nw = 2
               ENDIF
            ENDIF
            IF ( ilall ) THEN
               ilcomp = .TRUE.
            ELSEIF ( ilcplx ) THEN
               ilcomp = Select(je) .OR. Select(je-1)
            ELSE
               ilcomp = Select(je)
            ENDIF
            IF ( ilcomp ) THEN
!
!           Decide if (a) singular pencil, (b) real eigenvalue, or
!           (c) complex eigenvalue.
!
               IF ( .NOT.ilcplx ) THEN
                  IF ( ABS(S(je,je))<=safmin .AND. ABS(P(je,je))        &
     &                 <=safmin ) THEN
!
!                 Singular matrix pencil -- unit eigenvector
!
                     ieig = ieig - 1
                     DO jr = 1 , N
                        Vr(jr,ieig) = ZERO
                     ENDDO
                     Vr(ieig,ieig) = ONE
                     CYCLE
                  ENDIF
               ENDIF
!
!           Clear vector
!
               DO jw = 0 , nw - 1
                  DO jr = 1 , N
                     Work((jw+2)*N+jr) = ZERO
                  ENDDO
               ENDDO
!
!           Compute coefficients in  ( a A - b B ) x = 0
!              a  is  ACOEF
!              b  is  BCOEFR + i*BCOEFI
!
               IF ( .NOT.ilcplx ) THEN
!
!              Real eigenvalue
!
                  temp = ONE/MAX(ABS(S(je,je))*ascale,ABS(P(je,je))     &
     &                   *bscale,safmin)
                  salfar = (temp*S(je,je))*ascale
                  sbeta = (temp*P(je,je))*bscale
                  acoef = sbeta*ascale
                  bcoefr = salfar*bscale
                  bcoefi = ZERO
!
!              Scale to avoid underflow
!
                  scale = ONE
                  lsa = ABS(sbeta)>=safmin .AND. ABS(acoef)<small
                  lsb = ABS(salfar)>=safmin .AND. ABS(bcoefr)<small
                  IF ( lsa ) scale = (small/ABS(sbeta))*MIN(anorm,big)
                  IF ( lsb ) scale = MAX(scale,(small/ABS(salfar))*MIN( &
     &                               bnorm,big))
                  IF ( lsa .OR. lsb ) THEN
                     scale = MIN(scale,                                 &
     &                       ONE/(safmin*MAX(ONE,ABS(acoef),ABS(bcoefr))&
     &                       ))
                     IF ( lsa ) THEN
                        acoef = ascale*(scale*sbeta)
                     ELSE
                        acoef = scale*acoef
                     ENDIF
                     IF ( lsb ) THEN
                        bcoefr = bscale*(scale*salfar)
                     ELSE
                        bcoefr = scale*bcoefr
                     ENDIF
                  ENDIF
                  acoefa = ABS(acoef)
                  bcoefa = ABS(bcoefr)
!
!              First component is 1
!
                  Work(2*N+je) = ONE
                  xmax = ONE
!
!              Compute contribution from column JE of A and B to sum
!              (See "Further Details", above.)
!
                  DO jr = 1 , je - 1
                     Work(2*N+jr) = bcoefr*P(jr,je) - acoef*S(jr,je)
                  ENDDO
               ELSE
!
!              Complex eigenvalue
!
                  CALL SLAG2(S(je-1,je-1),Lds,P(je-1,je-1),Ldp,         &
     &                       safmin*SAFETY,acoef,temp,bcoefr,temp2,     &
     &                       bcoefi)
                  IF ( bcoefi==ZERO ) THEN
                     Info = je - 1
                     RETURN
                  ENDIF
!
!              Scale to avoid over/underflow
!
                  acoefa = ABS(acoef)
                  bcoefa = ABS(bcoefr) + ABS(bcoefi)
                  scale = ONE
                  IF ( acoefa*ulp<safmin .AND. acoefa>=safmin )         &
     &                 scale = (safmin/ulp)/acoefa
                  IF ( bcoefa*ulp<safmin .AND. bcoefa>=safmin )         &
     &                 scale = MAX(scale,(safmin/ulp)/bcoefa)
                  IF ( safmin*acoefa>ascale )                           &
     &                 scale = ascale/(safmin*acoefa)
                  IF ( safmin*bcoefa>bscale )                           &
     &                 scale = MIN(scale,bscale/(safmin*bcoefa))
                  IF ( scale/=ONE ) THEN
                     acoef = scale*acoef
                     acoefa = ABS(acoef)
                     bcoefr = scale*bcoefr
                     bcoefi = scale*bcoefi
                     bcoefa = ABS(bcoefr) + ABS(bcoefi)
                  ENDIF
!
!              Compute first two components of eigenvector
!              and contribution to sums
!
                  temp = acoef*S(je,je-1)
                  temp2r = acoef*S(je,je) - bcoefr*P(je,je)
                  temp2i = -bcoefi*P(je,je)
                  IF ( ABS(temp)>=ABS(temp2r)+ABS(temp2i) ) THEN
                     Work(2*N+je) = ONE
                     Work(3*N+je) = ZERO
                     Work(2*N+je-1) = -temp2r/temp
                     Work(3*N+je-1) = -temp2i/temp
                  ELSE
                     Work(2*N+je-1) = ONE
                     Work(3*N+je-1) = ZERO
                     temp = acoef*S(je-1,je)
                     Work(2*N+je) = (bcoefr*P(je-1,je-1)-acoef*S(je-1,je&
     &                              -1))/temp
                     Work(3*N+je) = bcoefi*P(je-1,je-1)/temp
                  ENDIF
!
                  xmax = MAX(ABS(Work(2*N+je))+ABS(Work(3*N+je)),       &
     &                   ABS(Work(2*N+je-1))+ABS(Work(3*N+je-1)))
!
!              Compute contribution from columns JE and JE-1
!              of A and B to the sums.
!
                  creala = acoef*Work(2*N+je-1)
                  cimaga = acoef*Work(3*N+je-1)
                  crealb = bcoefr*Work(2*N+je-1) - bcoefi*Work(3*N+je-1)
                  cimagb = bcoefi*Work(2*N+je-1) + bcoefr*Work(3*N+je-1)
                  cre2a = acoef*Work(2*N+je)
                  cim2a = acoef*Work(3*N+je)
                  cre2b = bcoefr*Work(2*N+je) - bcoefi*Work(3*N+je)
                  cim2b = bcoefi*Work(2*N+je) + bcoefr*Work(3*N+je)
                  DO jr = 1 , je - 2
                     Work(2*N+jr) = -creala*S(jr,je-1)                  &
     &                              + crealb*P(jr,je-1) - cre2a*S(jr,je)&
     &                              + cre2b*P(jr,je)
                     Work(3*N+jr) = -cimaga*S(jr,je-1)                  &
     &                              + cimagb*P(jr,je-1) - cim2a*S(jr,je)&
     &                              + cim2b*P(jr,je)
                  ENDDO
               ENDIF
!
               dmin = MAX(ulp*acoefa*anorm,ulp*bcoefa*bnorm,safmin)
!
!           Columnwise triangular solve of  (a A - b B)  x = 0
!
               il2by2 = .FALSE.
               DO j = je - nw , 1 , -1
!
!              If a 2-by-2 block, is in position j-1:j, wait until
!              next iteration to process it (when it will be j:j+1)
!
                  IF ( .NOT.il2by2 .AND. j>1 ) THEN
                     IF ( S(j,j-1)/=ZERO ) THEN
                        il2by2 = .TRUE.
                        CYCLE
                     ENDIF
                  ENDIF
                  bdiag(1) = P(j,j)
                  IF ( il2by2 ) THEN
                     na = 2
                     bdiag(2) = P(j+1,j+1)
                  ELSE
                     na = 1
                  ENDIF
!
!              Compute x(j) (and x(j+1), if 2-by-2 block)
!
                  CALL SLALN2(.FALSE.,na,nw,dmin,acoef,S(j,j),Lds,      &
     &                        bdiag(1),bdiag(2),Work(2*N+j),N,bcoefr,   &
     &                        bcoefi,sum,2,scale,temp,iinfo)
                  IF ( scale<ONE ) THEN
!
                     DO jw = 0 , nw - 1
                        DO jr = 1 , je
                           Work((jw+2)*N+jr) = scale*Work((jw+2)*N+jr)
                        ENDDO
                     ENDDO
                  ENDIF
                  xmax = MAX(scale*xmax,temp)
!
                  DO jw = 1 , nw
                     DO ja = 1 , na
                        Work((jw+1)*N+j+ja-1) = sum(ja,jw)
                     ENDDO
                  ENDDO
!
!              w = w + x(j)*(a S(*,j) - b P(*,j) ) with scaling
!
                  IF ( j>1 ) THEN
!
!                 Check whether scaling is necessary for sum.
!
                     xscale = ONE/MAX(ONE,xmax)
                     temp = acoefa*Work(j) + bcoefa*Work(N+j)
                     IF ( il2by2 )                                      &
     &                    temp = MAX(temp,acoefa*Work(j+1)+bcoefa*Work  &
     &                    (N+j+1))
                     temp = MAX(temp,acoefa,bcoefa)
                     IF ( temp>bignum*xscale ) THEN
!
                        DO jw = 0 , nw - 1
                           DO jr = 1 , je
                              Work((jw+2)*N+jr)                         &
     &                           = xscale*Work((jw+2)*N+jr)
                           ENDDO
                        ENDDO
                        xmax = xmax*xscale
                     ENDIF
!
!                 Compute the contributions of the off-diagonals of
!                 column j (and j+1, if 2-by-2 block) of A and B to the
!                 sums.
!
!
                     DO ja = 1 , na
                        IF ( ilcplx ) THEN
                           creala = acoef*Work(2*N+j+ja-1)
                           cimaga = acoef*Work(3*N+j+ja-1)
                           crealb = bcoefr*Work(2*N+j+ja-1)             &
     &                              - bcoefi*Work(3*N+j+ja-1)
                           cimagb = bcoefi*Work(2*N+j+ja-1)             &
     &                              + bcoefr*Work(3*N+j+ja-1)
                           DO jr = 1 , j - 1
                              Work(2*N+jr) = Work(2*N+jr)               &
     &                           - creala*S(jr,j+ja-1)                  &
     &                           + crealb*P(jr,j+ja-1)
                              Work(3*N+jr) = Work(3*N+jr)               &
     &                           - cimaga*S(jr,j+ja-1)                  &
     &                           + cimagb*P(jr,j+ja-1)
                           ENDDO
                        ELSE
                           creala = acoef*Work(2*N+j+ja-1)
                           crealb = bcoefr*Work(2*N+j+ja-1)
                           DO jr = 1 , j - 1
                              Work(2*N+jr) = Work(2*N+jr)               &
     &                           - creala*S(jr,j+ja-1)                  &
     &                           + crealb*P(jr,j+ja-1)
                           ENDDO
                        ENDIF
                     ENDDO
                  ENDIF
!
                  il2by2 = .FALSE.
               ENDDO
!
!           Copy eigenvector to VR, back transforming if
!           HOWMNY='B'.
!
               ieig = ieig - nw
               IF ( ilback ) THEN
!
                  DO jw = 0 , nw - 1
                     DO jr = 1 , N
                        Work((jw+4)*N+jr) = Work((jw+2)*N+1)*Vr(jr,1)
                     ENDDO
!
!                 A series of compiler directives to defeat
!                 vectorization for the next loop
!
!
                     DO jc = 2 , je
                        DO jr = 1 , N
                           Work((jw+4)*N+jr) = Work((jw+4)*N+jr)        &
     &                        + Work((jw+2)*N+jc)*Vr(jr,jc)
                        ENDDO
                     ENDDO
                  ENDDO
!
                  DO jw = 0 , nw - 1
                     DO jr = 1 , N
                        Vr(jr,ieig+jw) = Work((jw+4)*N+jr)
                     ENDDO
                  ENDDO
!
                  iend = N
               ELSE
                  DO jw = 0 , nw - 1
                     DO jr = 1 , N
                        Vr(jr,ieig+jw) = Work((jw+2)*N+jr)
                     ENDDO
                  ENDDO
!
                  iend = je
               ENDIF
!
!           Scale eigenvector
!
               xmax = ZERO
               IF ( ilcplx ) THEN
                  DO j = 1 , iend
                     xmax = MAX(xmax,ABS(Vr(j,ieig))+ABS(Vr(j,ieig+1)))
                  ENDDO
               ELSE
                  DO j = 1 , iend
                     xmax = MAX(xmax,ABS(Vr(j,ieig)))
                  ENDDO
               ENDIF
!
               IF ( xmax>safmin ) THEN
                  xscale = ONE/xmax
                  DO jw = 0 , nw - 1
                     DO jr = 1 , iend
                        Vr(jr,ieig+jw) = xscale*Vr(jr,ieig+jw)
                     ENDDO
                  ENDDO
               ENDIF
            ENDIF
         ENDDO
      ENDIF
!
!
!     End of STGEVC
!
      END SUBROUTINE STGEVC
