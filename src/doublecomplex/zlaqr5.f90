!*==zlaqr5.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b ZLAQR5 performs a single small-bulge multi-shift QR sweep.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZLAQR5 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlaqr5.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlaqr5.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlaqr5.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZLAQR5( WANTT, WANTZ, KACC22, N, KTOP, KBOT, NSHFTS, S,
!                          H, LDH, ILOZ, IHIZ, Z, LDZ, V, LDV, U, LDU, NV,
!                          WV, LDWV, NH, WH, LDWH )
!
!       .. Scalar Arguments ..
!       INTEGER            IHIZ, ILOZ, KACC22, KBOT, KTOP, LDH, LDU, LDV,
!      $                   LDWH, LDWV, LDZ, N, NH, NSHFTS, NV
!       LOGICAL            WANTT, WANTZ
!       ..
!       .. Array Arguments ..
!       COMPLEX*16         H( LDH, * ), S( * ), U( LDU, * ), V( LDV, * ),
!      $                   WH( LDWH, * ), WV( LDWV, * ), Z( LDZ, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    ZLAQR5, called by ZLAQR0, performs a
!>    single small-bulge multi-shift QR sweep.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] WANTT
!> \verbatim
!>          WANTT is LOGICAL
!>             WANTT = .true. if the triangular Schur factor
!>             is being computed.  WANTT is set to .false. otherwise.
!> \endverbatim
!>
!> \param[in] WANTZ
!> \verbatim
!>          WANTZ is LOGICAL
!>             WANTZ = .true. if the unitary Schur factor is being
!>             computed.  WANTZ is set to .false. otherwise.
!> \endverbatim
!>
!> \param[in] KACC22
!> \verbatim
!>          KACC22 is INTEGER with value 0, 1, or 2.
!>             Specifies the computation mode of far-from-diagonal
!>             orthogonal updates.
!>        = 0: ZLAQR5 does not accumulate reflections and does not
!>             use matrix-matrix multiply to update far-from-diagonal
!>             matrix entries.
!>        = 1: ZLAQR5 accumulates reflections and uses matrix-matrix
!>             multiply to update the far-from-diagonal matrix entries.
!>        = 2: Same as KACC22 = 1. This option used to enable exploiting
!>             the 2-by-2 structure during matrix multiplications, but
!>             this is no longer supported.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>             N is the order of the Hessenberg matrix H upon which this
!>             subroutine operates.
!> \endverbatim
!>
!> \param[in] KTOP
!> \verbatim
!>          KTOP is INTEGER
!> \endverbatim
!>
!> \param[in] KBOT
!> \verbatim
!>          KBOT is INTEGER
!>             These are the first and last rows and columns of an
!>             isolated diagonal block upon which the QR sweep is to be
!>             applied. It is assumed without a check that
!>                       either KTOP = 1  or   H(KTOP,KTOP-1) = 0
!>             and
!>                       either KBOT = N  or   H(KBOT+1,KBOT) = 0.
!> \endverbatim
!>
!> \param[in] NSHFTS
!> \verbatim
!>          NSHFTS is INTEGER
!>             NSHFTS gives the number of simultaneous shifts.  NSHFTS
!>             must be positive and even.
!> \endverbatim
!>
!> \param[in,out] S
!> \verbatim
!>          S is COMPLEX*16 array, dimension (NSHFTS)
!>             S contains the shifts of origin that define the multi-
!>             shift QR sweep.  On output S may be reordered.
!> \endverbatim
!>
!> \param[in,out] H
!> \verbatim
!>          H is COMPLEX*16 array, dimension (LDH,N)
!>             On input H contains a Hessenberg matrix.  On output a
!>             multi-shift QR sweep with shifts SR(J)+i*SI(J) is applied
!>             to the isolated diagonal block in rows and columns KTOP
!>             through KBOT.
!> \endverbatim
!>
!> \param[in] LDH
!> \verbatim
!>          LDH is INTEGER
!>             LDH is the leading dimension of H just as declared in the
!>             calling procedure.  LDH >= MAX(1,N).
!> \endverbatim
!>
!> \param[in] ILOZ
!> \verbatim
!>          ILOZ is INTEGER
!> \endverbatim
!>
!> \param[in] IHIZ
!> \verbatim
!>          IHIZ is INTEGER
!>             Specify the rows of Z to which transformations must be
!>             applied if WANTZ is .TRUE.. 1 <= ILOZ <= IHIZ <= N
!> \endverbatim
!>
!> \param[in,out] Z
!> \verbatim
!>          Z is COMPLEX*16 array, dimension (LDZ,IHIZ)
!>             If WANTZ = .TRUE., then the QR Sweep unitary
!>             similarity transformation is accumulated into
!>             Z(ILOZ:IHIZ,ILOZ:IHIZ) from the right.
!>             If WANTZ = .FALSE., then Z is unreferenced.
!> \endverbatim
!>
!> \param[in] LDZ
!> \verbatim
!>          LDZ is INTEGER
!>             LDA is the leading dimension of Z just as declared in
!>             the calling procedure. LDZ >= N.
!> \endverbatim
!>
!> \param[out] V
!> \verbatim
!>          V is COMPLEX*16 array, dimension (LDV,NSHFTS/2)
!> \endverbatim
!>
!> \param[in] LDV
!> \verbatim
!>          LDV is INTEGER
!>             LDV is the leading dimension of V as declared in the
!>             calling procedure.  LDV >= 3.
!> \endverbatim
!>
!> \param[out] U
!> \verbatim
!>          U is COMPLEX*16 array, dimension (LDU,2*NSHFTS)
!> \endverbatim
!>
!> \param[in] LDU
!> \verbatim
!>          LDU is INTEGER
!>             LDU is the leading dimension of U just as declared in the
!>             in the calling subroutine.  LDU >= 2*NSHFTS.
!> \endverbatim
!>
!> \param[in] NV
!> \verbatim
!>          NV is INTEGER
!>             NV is the number of rows in WV agailable for workspace.
!>             NV >= 1.
!> \endverbatim
!>
!> \param[out] WV
!> \verbatim
!>          WV is COMPLEX*16 array, dimension (LDWV,2*NSHFTS)
!> \endverbatim
!>
!> \param[in] LDWV
!> \verbatim
!>          LDWV is INTEGER
!>             LDWV is the leading dimension of WV as declared in the
!>             in the calling subroutine.  LDWV >= NV.
!> \endverbatim
!
!> \param[in] NH
!> \verbatim
!>          NH is INTEGER
!>             NH is the number of columns in array WH available for
!>             workspace. NH >= 1.
!> \endverbatim
!>
!> \param[out] WH
!> \verbatim
!>          WH is COMPLEX*16 array, dimension (LDWH,NH)
!> \endverbatim
!>
!> \param[in] LDWH
!> \verbatim
!>          LDWH is INTEGER
!>             Leading dimension of WH just as declared in the
!>             calling procedure.  LDWH >= 2*NSHFTS.
!> \endverbatim
!>
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \date January 2021
!
!> \ingroup complex16OTHERauxiliary
!
!> \par Contributors:
!  ==================
!>
!>       Karen Braman and Ralph Byers, Department of Mathematics,
!>       University of Kansas, USA
!>
!>       Lars Karlsson, Daniel Kressner, and Bruno Lang
!>
!>       Thijs Steel, Department of Computer science,
!>       KU Leuven, Belgium
!
!> \par References:
!  ================
!>
!>       K. Braman, R. Byers and R. Mathias, The Multi-Shift QR
!>       Algorithm Part I: Maintaining Well Focused Shifts, and Level 3
!>       Performance, SIAM Journal of Matrix Analysis, volume 23, pages
!>       929--947, 2002.
!>
!>       Lars Karlsson, Daniel Kressner, and Bruno Lang, Optimally packed
!>       chains of bulges in multishift QR algorithms.
!>       ACM Trans. Math. Softw. 40, 2, Article 12 (February 2014).
!>
!  =====================================================================
      SUBROUTINE ZLAQR5(Wantt,Wantz,Kacc22,N,Ktop,Kbot,Nshfts,S,H,Ldh,  &
     &                  Iloz,Ihiz,Z,Ldz,V,Ldv,U,Ldu,Nv,Wv,Ldwv,Nh,Wh,   &
     &                  Ldwh)
      IMPLICIT NONE
!*--ZLAQR5261
!
!  -- LAPACK auxiliary routine (version 3.7.1) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     June 2016
!
!     .. Scalar Arguments ..
      INTEGER Ihiz , Iloz , Kacc22 , Kbot , Ktop , Ldh , Ldu , Ldv ,    &
     &        Ldwh , Ldwv , Ldz , N , Nh , Nshfts , Nv
      LOGICAL Wantt , Wantz
!     ..
!     .. Array Arguments ..
      COMPLEX*16 H(Ldh,*) , S(*) , U(Ldu,*) , V(Ldv,*) , Wh(Ldwh,*) ,   &
     &           Wv(Ldwv,*) , Z(Ldz,*)
!     ..
!
!  ================================================================
!     .. Parameters ..
      COMPLEX*16 ZERO , ONE
      PARAMETER (ZERO=(0.0D0,0.0D0),ONE=(1.0D0,0.0D0))
      DOUBLE PRECISION RZERO , RONE
      PARAMETER (RZERO=0.0D0,RONE=1.0D0)
!     ..
!     .. Local Scalars ..
      COMPLEX*16 alpha , beta , cdum , refsum
      DOUBLE PRECISION h11 , h12 , h21 , h22 , safmax , safmin , scl ,  &
     &                 smlnum , tst1 , tst2 , ulp
      INTEGER i2 , i4 , incol , j , jbot , jcol , jlen , jrow , jtop ,  &
     &        k , k1 , kdu , kms , krcol , m , m22 , mbot , mtop ,      &
     &        nbmps , ndcol , ns , nu
      LOGICAL accum , bmp22
!     ..
!     .. External Functions ..
      DOUBLE PRECISION DLAMCH
      EXTERNAL DLAMCH
!     ..
!     .. Intrinsic Functions ..
!
      INTRINSIC ABS , DBLE , DCONJG , DIMAG , MAX , MIN , MOD
!     ..
!     .. Local Arrays ..
      COMPLEX*16 vt(3)
!     ..
!     .. External Subroutines ..
      EXTERNAL DLABAD , ZGEMM , ZLACPY , ZLAQR1 , ZLARFG , ZLASET ,     &
     &         ZTRMM
!     ..
!     .. Statement Functions ..
      DOUBLE PRECISION CABS1
!     ..
!     .. Statement Function definitions ..
      CABS1(cdum) = ABS(DBLE(cdum)) + ABS(DIMAG(cdum))
!     ..
!     .. Executable Statements ..
!
!     ==== If there are no shifts, then there is nothing to do. ====
!
      IF ( Nshfts<2 ) RETURN
!
!     ==== If the active block is empty or 1-by-1, then there
!     .    is nothing to do. ====
!
      IF ( Ktop>=Kbot ) RETURN
!
!     ==== NSHFTS is supposed to be even, but if it is odd,
!     .    then simply reduce it by one.  ====
!
      ns = Nshfts - MOD(Nshfts,2)
!
!     ==== Machine constants for deflation ====
!
      safmin = DLAMCH('SAFE MINIMUM')
      safmax = RONE/safmin
      CALL DLABAD(safmin,safmax)
      ulp = DLAMCH('PRECISION')
      smlnum = safmin*(DBLE(N)/ulp)
!
!     ==== Use accumulated reflections to update far-from-diagonal
!     .    entries ? ====
!
      accum = (Kacc22==1) .OR. (Kacc22==2)
!
!     ==== clear trash ====
!
      IF ( Ktop+2<=Kbot ) H(Ktop+2,Ktop) = ZERO
!
!     ==== NBMPS = number of 2-shift bulges in the chain ====
!
      nbmps = ns/2
!
!     ==== KDU = width of slab ====
!
      kdu = 4*nbmps
!
!     ==== Create and chase chains of NBMPS bulges ====
!
      DO incol = Ktop - 2*nbmps + 1 , Kbot - 2 , 2*nbmps
!
!        JTOP = Index from which updates from the right start.
!
         IF ( accum ) THEN
            jtop = MAX(Ktop,incol)
         ELSEIF ( Wantt ) THEN
            jtop = 1
         ELSE
            jtop = Ktop
         ENDIF
!
         ndcol = incol + kdu
         IF ( accum ) CALL ZLASET('ALL',kdu,kdu,ZERO,ONE,U,Ldu)
!
!        ==== Near-the-diagonal bulge chase.  The following loop
!        .    performs the near-the-diagonal part of a small bulge
!        .    multi-shift QR sweep.  Each 4*NBMPS column diagonal
!        .    chunk extends from column INCOL to column NDCOL
!        .    (including both column INCOL and column NDCOL). The
!        .    following loop chases a 2*NBMPS+1 column long chain of
!        .    NBMPS bulges 2*NBMPS columns to the right.  (INCOL
!        .    may be less than KTOP and and NDCOL may be greater than
!        .    KBOT indicating phantom columns from which to chase
!        .    bulges before they are actually introduced or to which
!        .    to chase bulges beyond column KBOT.)  ====
!
         DO krcol = incol , MIN(incol+2*nbmps-1,Kbot-2)
!
!           ==== Bulges number MTOP to MBOT are active double implicit
!           .    shift bulges.  There may or may not also be small
!           .    2-by-2 bulge, if there is room.  The inactive bulges
!           .    (if any) must wait until the active bulges have moved
!           .    down the diagonal to make room.  The phantom matrix
!           .    paradigm described above helps keep track.  ====
!
            mtop = MAX(1,(Ktop-krcol)/2+1)
            mbot = MIN(nbmps,(Kbot-krcol-1)/2)
            m22 = mbot + 1
            bmp22 = (mbot<nbmps) .AND. (krcol+2*(m22-1))==(Kbot-2)
!
!           ==== Generate reflections to chase the chain right
!           .    one column.  (The minimum value of K is KTOP-1.) ====
!
            IF ( bmp22 ) THEN
!
!              ==== Special case: 2-by-2 reflection at bottom treated
!              .    separately ====
!
               k = krcol + 2*(m22-1)
               IF ( k==Ktop-1 ) THEN
                  CALL ZLAQR1(2,H(k+1,k+1),Ldh,S(2*m22-1),S(2*m22),     &
     &                        V(1,m22))
                  beta = V(1,m22)
                  CALL ZLARFG(2,beta,V(2,m22),1,V(1,m22))
               ELSE
                  beta = H(k+1,k)
                  V(2,m22) = H(k+2,k)
                  CALL ZLARFG(2,beta,V(2,m22),1,V(1,m22))
                  H(k+1,k) = beta
                  H(k+2,k) = ZERO
               ENDIF
 
!
!              ==== Perform update from right within
!              .    computational window. ====
!
               DO j = jtop , MIN(Kbot,k+3)
                  refsum = V(1,m22)*(H(j,k+1)+V(2,m22)*H(j,k+2))
                  H(j,k+1) = H(j,k+1) - refsum
                  H(j,k+2) = H(j,k+2) - refsum*DCONJG(V(2,m22))
               ENDDO
!
!              ==== Perform update from left within
!              .    computational window. ====
!
               IF ( accum ) THEN
                  jbot = MIN(ndcol,Kbot)
               ELSEIF ( Wantt ) THEN
                  jbot = N
               ELSE
                  jbot = Kbot
               ENDIF
               DO j = k + 1 , jbot
                  refsum = DCONJG(V(1,m22))                             &
     &                     *(H(k+1,j)+DCONJG(V(2,m22))*H(k+2,j))
                  H(k+1,j) = H(k+1,j) - refsum
                  H(k+2,j) = H(k+2,j) - refsum*V(2,m22)
               ENDDO
!
!              ==== The following convergence test requires that
!              .    the tradition small-compared-to-nearby-diagonals
!              .    criterion and the Ahues & Tisseur (LAWN 122, 1997)
!              .    criteria both be satisfied.  The latter improves
!              .    accuracy in some examples. Falling back on an
!              .    alternate convergence criterion when TST1 or TST2
!              .    is zero (as done here) is traditional but probably
!              .    unnecessary. ====
!
               IF ( k>=Ktop ) THEN
                  IF ( H(k+1,k)/=ZERO ) THEN
                     tst1 = CABS1(H(k,k)) + CABS1(H(k+1,k+1))
                     IF ( tst1==RZERO ) THEN
                        IF ( k>=Ktop+1 ) tst1 = tst1 + CABS1(H(k,k-1))
                        IF ( k>=Ktop+2 ) tst1 = tst1 + CABS1(H(k,k-2))
                        IF ( k>=Ktop+3 ) tst1 = tst1 + CABS1(H(k,k-3))
                        IF ( k<=Kbot-2 ) tst1 = tst1 + CABS1(H(k+2,k+1))
                        IF ( k<=Kbot-3 ) tst1 = tst1 + CABS1(H(k+3,k+1))
                        IF ( k<=Kbot-4 ) tst1 = tst1 + CABS1(H(k+4,k+1))
                     ENDIF
                     IF ( CABS1(H(k+1,k))<=MAX(smlnum,ulp*tst1) ) THEN
                        h12 = MAX(CABS1(H(k+1,k)),CABS1(H(k,k+1)))
                        h21 = MIN(CABS1(H(k+1,k)),CABS1(H(k,k+1)))
                        h11 = MAX(CABS1(H(k+1,k+1)),                    &
     &                        CABS1(H(k,k)-H(k+1,k+1)))
                        h22 = MIN(CABS1(H(k+1,k+1)),                    &
     &                        CABS1(H(k,k)-H(k+1,k+1)))
                        scl = h11 + h12
                        tst2 = h22*(h11/scl)
!
                        IF ( tst2==RZERO .OR. h21*(h12/scl)             &
     &                       <=MAX(smlnum,ulp*tst2) ) H(k+1,k) = ZERO
                     ENDIF
                  ENDIF
               ENDIF
!
!              ==== Accumulate orthogonal transformations. ====
!
               IF ( accum ) THEN
                  kms = k - incol
                  DO j = MAX(1,Ktop-incol) , kdu
                     refsum = V(1,m22)*(U(j,kms+1)+V(2,m22)*U(j,kms+2))
                     U(j,kms+1) = U(j,kms+1) - refsum
                     U(j,kms+2) = U(j,kms+2) - refsum*DCONJG(V(2,m22))
                  ENDDO
               ELSEIF ( Wantz ) THEN
                  DO j = Iloz , Ihiz
                     refsum = V(1,m22)*(Z(j,k+1)+V(2,m22)*Z(j,k+2))
                     Z(j,k+1) = Z(j,k+1) - refsum
                     Z(j,k+2) = Z(j,k+2) - refsum*DCONJG(V(2,m22))
                  ENDDO
               ENDIF
            ENDIF
!
!           ==== Normal case: Chain of 3-by-3 reflections ====
!
            DO m = mbot , mtop , -1
               k = krcol + 2*(m-1)
               IF ( k==Ktop-1 ) THEN
                  CALL ZLAQR1(3,H(Ktop,Ktop),Ldh,S(2*m-1),S(2*m),V(1,m))
                  alpha = V(1,m)
                  CALL ZLARFG(3,alpha,V(2,m),1,V(1,m))
               ELSE
!
!                 ==== Perform delayed transformation of row below
!                 .    Mth bulge. Exploit fact that first two elements
!                 .    of row are actually zero. ====
!
                  refsum = V(1,m)*V(3,m)*H(k+3,k+2)
                  H(k+3,k) = -refsum
                  H(k+3,k+1) = -refsum*DCONJG(V(2,m))
                  H(k+3,k+2) = H(k+3,k+2) - refsum*DCONJG(V(3,m))
!
!                 ==== Calculate reflection to move
!                 .    Mth bulge one step. ====
!
                  beta = H(k+1,k)
                  V(2,m) = H(k+2,k)
                  V(3,m) = H(k+3,k)
                  CALL ZLARFG(3,beta,V(2,m),1,V(1,m))
!
!                 ==== A Bulge may collapse because of vigilant
!                 .    deflation or destructive underflow.  In the
!                 .    underflow case, try the two-small-subdiagonals
!                 .    trick to try to reinflate the bulge.  ====
!
                  IF ( H(k+3,k)/=ZERO .OR. H(k+3,k+1)/=ZERO .OR.        &
     &                 H(k+3,k+2)==ZERO ) THEN
!
!                    ==== Typical case: not collapsed (yet). ====
!
                     H(k+1,k) = beta
                     H(k+2,k) = ZERO
                     H(k+3,k) = ZERO
                  ELSE
!
!                    ==== Atypical case: collapsed.  Attempt to
!                    .    reintroduce ignoring H(K+1,K) and H(K+2,K).
!                    .    If the fill resulting from the new
!                    .    reflector is too large, then abandon it.
!                    .    Otherwise, use the new one. ====
!
                     CALL ZLAQR1(3,H(k+1,k+1),Ldh,S(2*m-1),S(2*m),vt)
                     alpha = vt(1)
                     CALL ZLARFG(3,alpha,vt(2),1,vt(1))
                     refsum = DCONJG(vt(1))                             &
     &                        *(H(k+1,k)+DCONJG(vt(2))*H(k+2,k))
!
                     IF ( CABS1(H(k+2,k)-refsum*vt(2))                  &
     &                    +CABS1(refsum*vt(3))                          &
     &                    >ulp*(CABS1(H(k,k))+CABS1(H(k+1,k+1))         &
     &                    +CABS1(H(k+2,k+2))) ) THEN
!
!                       ==== Starting a new bulge here would
!                       .    create non-negligible fill.  Use
!                       .    the old one with trepidation. ====
!
                        H(k+1,k) = beta
                        H(k+2,k) = ZERO
                        H(k+3,k) = ZERO
                     ELSE
!
!                       ==== Starting a new bulge here would
!                       .    create only negligible fill.
!                       .    Replace the old reflector with
!                       .    the new one. ====
!
                        H(k+1,k) = H(k+1,k) - refsum
                        H(k+2,k) = ZERO
                        H(k+3,k) = ZERO
                        V(1,m) = vt(1)
                        V(2,m) = vt(2)
                        V(3,m) = vt(3)
                     ENDIF
                  ENDIF
               ENDIF
!
!              ====  Apply reflection from the right and
!              .     the first column of update from the left.
!              .     These updates are required for the vigilant
!              .     deflation check. We still delay most of the
!              .     updates from the left for efficiency. ====
!
               DO j = jtop , MIN(Kbot,k+3)
                  refsum = V(1,m)                                       &
     &                     *(H(j,k+1)+V(2,m)*H(j,k+2)+V(3,m)*H(j,k+3))
                  H(j,k+1) = H(j,k+1) - refsum
                  H(j,k+2) = H(j,k+2) - refsum*DCONJG(V(2,m))
                  H(j,k+3) = H(j,k+3) - refsum*DCONJG(V(3,m))
               ENDDO
!
!              ==== Perform update from left for subsequent
!              .    column. ====
!
               refsum = DCONJG(V(1,m))                                  &
     &                  *(H(k+1,k+1)+DCONJG(V(2,m))*H(k+2,k+1)          &
     &                  +DCONJG(V(3,m))*H(k+3,k+1))
               H(k+1,k+1) = H(k+1,k+1) - refsum
               H(k+2,k+1) = H(k+2,k+1) - refsum*V(2,m)
               H(k+3,k+1) = H(k+3,k+1) - refsum*V(3,m)
!
!              ==== The following convergence test requires that
!              .    the tradition small-compared-to-nearby-diagonals
!              .    criterion and the Ahues & Tisseur (LAWN 122, 1997)
!              .    criteria both be satisfied.  The latter improves
!              .    accuracy in some examples. Falling back on an
!              .    alternate convergence criterion when TST1 or TST2
!              .    is zero (as done here) is traditional but probably
!              .    unnecessary. ====
!
               IF ( k<Ktop ) CYCLE
               IF ( H(k+1,k)/=ZERO ) THEN
                  tst1 = CABS1(H(k,k)) + CABS1(H(k+1,k+1))
                  IF ( tst1==RZERO ) THEN
                     IF ( k>=Ktop+1 ) tst1 = tst1 + CABS1(H(k,k-1))
                     IF ( k>=Ktop+2 ) tst1 = tst1 + CABS1(H(k,k-2))
                     IF ( k>=Ktop+3 ) tst1 = tst1 + CABS1(H(k,k-3))
                     IF ( k<=Kbot-2 ) tst1 = tst1 + CABS1(H(k+2,k+1))
                     IF ( k<=Kbot-3 ) tst1 = tst1 + CABS1(H(k+3,k+1))
                     IF ( k<=Kbot-4 ) tst1 = tst1 + CABS1(H(k+4,k+1))
                  ENDIF
                  IF ( CABS1(H(k+1,k))<=MAX(smlnum,ulp*tst1) ) THEN
                     h12 = MAX(CABS1(H(k+1,k)),CABS1(H(k,k+1)))
                     h21 = MIN(CABS1(H(k+1,k)),CABS1(H(k,k+1)))
                     h11 = MAX(CABS1(H(k+1,k+1)),                       &
     &                     CABS1(H(k,k)-H(k+1,k+1)))
                     h22 = MIN(CABS1(H(k+1,k+1)),                       &
     &                     CABS1(H(k,k)-H(k+1,k+1)))
                     scl = h11 + h12
                     tst2 = h22*(h11/scl)
!
                     IF ( tst2==RZERO .OR. h21*(h12/scl)                &
     &                    <=MAX(smlnum,ulp*tst2) ) H(k+1,k) = ZERO
                  ENDIF
               ENDIF
            ENDDO
!
!           ==== Multiply H by reflections from the left ====
!
            IF ( accum ) THEN
               jbot = MIN(ndcol,Kbot)
            ELSEIF ( Wantt ) THEN
               jbot = N
            ELSE
               jbot = Kbot
            ENDIF
!
            DO m = mbot , mtop , -1
               k = krcol + 2*(m-1)
               DO j = MAX(Ktop,krcol+2*m) , jbot
                  refsum = DCONJG(V(1,m))                               &
     &                     *(H(k+1,j)+DCONJG(V(2,m))*H(k+2,j)           &
     &                     +DCONJG(V(3,m))*H(k+3,j))
                  H(k+1,j) = H(k+1,j) - refsum
                  H(k+2,j) = H(k+2,j) - refsum*V(2,m)
                  H(k+3,j) = H(k+3,j) - refsum*V(3,m)
               ENDDO
            ENDDO
!
!           ==== Accumulate orthogonal transformations. ====
!
            IF ( accum ) THEN
!
!              ==== Accumulate U. (If needed, update Z later
!              .    with an efficient matrix-matrix
!              .    multiply.) ====
!
               DO m = mbot , mtop , -1
                  k = krcol + 2*(m-1)
                  kms = k - incol
                  i2 = MAX(1,Ktop-incol)
                  i2 = MAX(i2,kms-(krcol-incol)+1)
                  i4 = MIN(kdu,krcol+2*(mbot-1)-incol+5)
                  DO j = i2 , i4
                     refsum = V(1,m)                                    &
     &                        *(U(j,kms+1)+V(2,m)*U(j,kms+2)+V(3,m)     &
     &                        *U(j,kms+3))
                     U(j,kms+1) = U(j,kms+1) - refsum
                     U(j,kms+2) = U(j,kms+2) - refsum*DCONJG(V(2,m))
                     U(j,kms+3) = U(j,kms+3) - refsum*DCONJG(V(3,m))
                  ENDDO
               ENDDO
            ELSEIF ( Wantz ) THEN
!
!              ==== U is not accumulated, so update Z
!              .    now by multiplying by reflections
!              .    from the right. ====
!
               DO m = mbot , mtop , -1
                  k = krcol + 2*(m-1)
                  DO j = Iloz , Ihiz
                     refsum = V(1,m)                                    &
     &                        *(Z(j,k+1)+V(2,m)*Z(j,k+2)+V(3,m)*Z(j,k+3)&
     &                        )
                     Z(j,k+1) = Z(j,k+1) - refsum
                     Z(j,k+2) = Z(j,k+2) - refsum*DCONJG(V(2,m))
                     Z(j,k+3) = Z(j,k+3) - refsum*DCONJG(V(3,m))
                  ENDDO
               ENDDO
            ENDIF
!
!           ==== End of near-the-diagonal bulge chase. ====
!
         ENDDO
!
!        ==== Use U (if accumulated) to update far-from-diagonal
!        .    entries in H.  If required, use U to update Z as
!        .    well. ====
!
         IF ( accum ) THEN
            IF ( Wantt ) THEN
               jtop = 1
               jbot = N
            ELSE
               jtop = Ktop
               jbot = Kbot
            ENDIF
            k1 = MAX(1,Ktop-incol)
            nu = (kdu-MAX(0,ndcol-Kbot)) - k1 + 1
!
!           ==== Horizontal Multiply ====
!
            DO jcol = MIN(ndcol,Kbot) + 1 , jbot , Nh
               jlen = MIN(Nh,jbot-jcol+1)
               CALL ZGEMM('C','N',nu,jlen,nu,ONE,U(k1,k1),Ldu,          &
     &                    H(incol+k1,jcol),Ldh,ZERO,Wh,Ldwh)
               CALL ZLACPY('ALL',nu,jlen,Wh,Ldwh,H(incol+k1,jcol),Ldh)
            ENDDO
!
!           ==== Vertical multiply ====
!
            DO jrow = jtop , MAX(Ktop,incol) - 1 , Nv
               jlen = MIN(Nv,MAX(Ktop,incol)-jrow)
               CALL ZGEMM('N','N',jlen,nu,nu,ONE,H(jrow,incol+k1),Ldh,  &
     &                    U(k1,k1),Ldu,ZERO,Wv,Ldwv)
               CALL ZLACPY('ALL',jlen,nu,Wv,Ldwv,H(jrow,incol+k1),Ldh)
            ENDDO
!
!           ==== Z multiply (also vertical) ====
!
            IF ( Wantz ) THEN
               DO jrow = Iloz , Ihiz , Nv
                  jlen = MIN(Nv,Ihiz-jrow+1)
                  CALL ZGEMM('N','N',jlen,nu,nu,ONE,Z(jrow,incol+k1),   &
     &                       Ldz,U(k1,k1),Ldu,ZERO,Wv,Ldwv)
                  CALL ZLACPY('ALL',jlen,nu,Wv,Ldwv,Z(jrow,incol+k1),   &
     &                        Ldz)
               ENDDO
            ENDIF
         ENDIF
      ENDDO
!
!     ==== End of ZLAQR5 ====
!
      END SUBROUTINE ZLAQR5
