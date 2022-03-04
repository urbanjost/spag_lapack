!*==dlaqr2.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b DLAQR2 performs the orthogonal similarity transformation of a Hessenberg matrix to detect and deflate fully converged eigenvalues from a trailing principal submatrix (aggressive early deflation).
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DLAQR2 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaqr2.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaqr2.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaqr2.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLAQR2( WANTT, WANTZ, N, KTOP, KBOT, NW, H, LDH, ILOZ,
!                          IHIZ, Z, LDZ, NS, ND, SR, SI, V, LDV, NH, T,
!                          LDT, NV, WV, LDWV, WORK, LWORK )
!
!       .. Scalar Arguments ..
!       INTEGER            IHIZ, ILOZ, KBOT, KTOP, LDH, LDT, LDV, LDWV,
!      $                   LDZ, LWORK, N, ND, NH, NS, NV, NW
!       LOGICAL            WANTT, WANTZ
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   H( LDH, * ), SI( * ), SR( * ), T( LDT, * ),
!      $                   V( LDV, * ), WORK( * ), WV( LDWV, * ),
!      $                   Z( LDZ, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    DLAQR2 is identical to DLAQR3 except that it avoids
!>    recursion by calling DLAHQR instead of DLAQR4.
!>
!>    Aggressive early deflation:
!>
!>    This subroutine accepts as input an upper Hessenberg matrix
!>    H and performs an orthogonal similarity transformation
!>    designed to detect and deflate fully converged eigenvalues from
!>    a trailing principal submatrix.  On output H has been over-
!>    written by a new Hessenberg matrix that is a perturbation of
!>    an orthogonal similarity transformation of H.  It is to be
!>    hoped that the final version of H has many zero subdiagonal
!>    entries.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] WANTT
!> \verbatim
!>          WANTT is LOGICAL
!>          If .TRUE., then the Hessenberg matrix H is fully updated
!>          so that the quasi-triangular Schur factor may be
!>          computed (in cooperation with the calling subroutine).
!>          If .FALSE., then only enough of H is updated to preserve
!>          the eigenvalues.
!> \endverbatim
!>
!> \param[in] WANTZ
!> \verbatim
!>          WANTZ is LOGICAL
!>          If .TRUE., then the orthogonal matrix Z is updated so
!>          so that the orthogonal Schur factor may be computed
!>          (in cooperation with the calling subroutine).
!>          If .FALSE., then Z is not referenced.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix H and (if WANTZ is .TRUE.) the
!>          order of the orthogonal matrix Z.
!> \endverbatim
!>
!> \param[in] KTOP
!> \verbatim
!>          KTOP is INTEGER
!>          It is assumed that either KTOP = 1 or H(KTOP,KTOP-1)=0.
!>          KBOT and KTOP together determine an isolated block
!>          along the diagonal of the Hessenberg matrix.
!> \endverbatim
!>
!> \param[in] KBOT
!> \verbatim
!>          KBOT is INTEGER
!>          It is assumed without a check that either
!>          KBOT = N or H(KBOT+1,KBOT)=0.  KBOT and KTOP together
!>          determine an isolated block along the diagonal of the
!>          Hessenberg matrix.
!> \endverbatim
!>
!> \param[in] NW
!> \verbatim
!>          NW is INTEGER
!>          Deflation window size.  1 <= NW <= (KBOT-KTOP+1).
!> \endverbatim
!>
!> \param[in,out] H
!> \verbatim
!>          H is DOUBLE PRECISION array, dimension (LDH,N)
!>          On input the initial N-by-N section of H stores the
!>          Hessenberg matrix undergoing aggressive early deflation.
!>          On output H has been transformed by an orthogonal
!>          similarity transformation, perturbed, and the returned
!>          to Hessenberg form that (it is to be hoped) has some
!>          zero subdiagonal entries.
!> \endverbatim
!>
!> \param[in] LDH
!> \verbatim
!>          LDH is INTEGER
!>          Leading dimension of H just as declared in the calling
!>          subroutine.  N <= LDH
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
!>          Specify the rows of Z to which transformations must be
!>          applied if WANTZ is .TRUE.. 1 <= ILOZ <= IHIZ <= N.
!> \endverbatim
!>
!> \param[in,out] Z
!> \verbatim
!>          Z is DOUBLE PRECISION array, dimension (LDZ,N)
!>          IF WANTZ is .TRUE., then on output, the orthogonal
!>          similarity transformation mentioned above has been
!>          accumulated into Z(ILOZ:IHIZ,ILOZ:IHIZ) from the right.
!>          If WANTZ is .FALSE., then Z is unreferenced.
!> \endverbatim
!>
!> \param[in] LDZ
!> \verbatim
!>          LDZ is INTEGER
!>          The leading dimension of Z just as declared in the
!>          calling subroutine.  1 <= LDZ.
!> \endverbatim
!>
!> \param[out] NS
!> \verbatim
!>          NS is INTEGER
!>          The number of unconverged (ie approximate) eigenvalues
!>          returned in SR and SI that may be used as shifts by the
!>          calling subroutine.
!> \endverbatim
!>
!> \param[out] ND
!> \verbatim
!>          ND is INTEGER
!>          The number of converged eigenvalues uncovered by this
!>          subroutine.
!> \endverbatim
!>
!> \param[out] SR
!> \verbatim
!>          SR is DOUBLE PRECISION array, dimension (KBOT)
!> \endverbatim
!>
!> \param[out] SI
!> \verbatim
!>          SI is DOUBLE PRECISION array, dimension (KBOT)
!>          On output, the real and imaginary parts of approximate
!>          eigenvalues that may be used for shifts are stored in
!>          SR(KBOT-ND-NS+1) through SR(KBOT-ND) and
!>          SI(KBOT-ND-NS+1) through SI(KBOT-ND), respectively.
!>          The real and imaginary parts of converged eigenvalues
!>          are stored in SR(KBOT-ND+1) through SR(KBOT) and
!>          SI(KBOT-ND+1) through SI(KBOT), respectively.
!> \endverbatim
!>
!> \param[out] V
!> \verbatim
!>          V is DOUBLE PRECISION array, dimension (LDV,NW)
!>          An NW-by-NW work array.
!> \endverbatim
!>
!> \param[in] LDV
!> \verbatim
!>          LDV is INTEGER
!>          The leading dimension of V just as declared in the
!>          calling subroutine.  NW <= LDV
!> \endverbatim
!>
!> \param[in] NH
!> \verbatim
!>          NH is INTEGER
!>          The number of columns of T.  NH >= NW.
!> \endverbatim
!>
!> \param[out] T
!> \verbatim
!>          T is DOUBLE PRECISION array, dimension (LDT,NW)
!> \endverbatim
!>
!> \param[in] LDT
!> \verbatim
!>          LDT is INTEGER
!>          The leading dimension of T just as declared in the
!>          calling subroutine.  NW <= LDT
!> \endverbatim
!>
!> \param[in] NV
!> \verbatim
!>          NV is INTEGER
!>          The number of rows of work array WV available for
!>          workspace.  NV >= NW.
!> \endverbatim
!>
!> \param[out] WV
!> \verbatim
!>          WV is DOUBLE PRECISION array, dimension (LDWV,NW)
!> \endverbatim
!>
!> \param[in] LDWV
!> \verbatim
!>          LDWV is INTEGER
!>          The leading dimension of W just as declared in the
!>          calling subroutine.  NW <= LDV
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (LWORK)
!>          On exit, WORK(1) is set to an estimate of the optimal value
!>          of LWORK for the given values of N, NW, KTOP and KBOT.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The dimension of the work array WORK.  LWORK = 2*NW
!>          suffices, but greater efficiency may result from larger
!>          values of LWORK.
!>
!>          If LWORK = -1, then a workspace query is assumed; DLAQR2
!>          only estimates the optimal workspace size for the given
!>          values of N, NW, KTOP and KBOT.  The estimate is returned
!>          in WORK(1).  No error message related to LWORK is issued
!>          by XERBLA.  Neither H nor Z are accessed.
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
!> \ingroup doubleOTHERauxiliary
!
!> \par Contributors:
!  ==================
!>
!>       Karen Braman and Ralph Byers, Department of Mathematics,
!>       University of Kansas, USA
!>
!  =====================================================================
      SUBROUTINE DLAQR2(Wantt,Wantz,N,Ktop,Kbot,Nw,H,Ldh,Iloz,Ihiz,Z,   &
     &                  Ldz,Ns,Nd,Sr,Si,V,Ldv,Nh,T,Ldt,Nv,Wv,Ldwv,Work, &
     &                  Lwork)
      USE F77KINDS                        
      USE S_DCOPY
      USE S_DGEHRD
      USE S_DGEMM
      USE S_DLABAD
      USE S_DLACPY
      USE S_DLAHQR
      USE S_DLAMCH
      USE S_DLANV2
      USE S_DLARF
      USE S_DLARFG
      USE S_DLASET
      USE S_DORMHR
      USE S_DTREXC
      IMPLICIT NONE
!*--DLAQR2296
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D0 , ONE = 1.0D0
!
! Dummy argument declarations rewritten by SPAG
!
      LOGICAL , INTENT(IN) :: Wantt
      LOGICAL , INTENT(IN) :: Wantz
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Ktop
      INTEGER , INTENT(IN) :: Kbot
      INTEGER , INTENT(IN) :: Nw
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldh,*) :: H
      INTEGER :: Ldh
      INTEGER , INTENT(IN) :: Iloz
      INTEGER , INTENT(IN) :: Ihiz
      REAL(R8KIND) , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      INTEGER , INTENT(INOUT) :: Ns
      INTEGER , INTENT(OUT) :: Nd
      REAL(R8KIND) , DIMENSION(*) :: Sr
      REAL(R8KIND) , DIMENSION(*) :: Si
      REAL(R8KIND) , DIMENSION(Ldv,*) :: V
      INTEGER :: Ldv
      INTEGER , INTENT(IN) :: Nh
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      INTEGER , INTENT(IN) :: Nv
      REAL(R8KIND) , DIMENSION(Ldwv,*) :: Wv
      INTEGER :: Ldwv
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
!
! Local variable declarations rewritten by SPAG
!
      REAL(R8KIND) :: aa , bb , beta , cc , cs , dd , evi , evk , foo , &
     &                s , safmax , safmin , smlnum , sn , tau , ulp
      LOGICAL :: bulge , sorted
      INTEGER :: i , ifst , ilst , info , infqr , j , jw , k , kcol ,   &
     &           kend , kln , krow , kwtop , ltop , lwk1 , lwk2 , lwkopt
!
! End of declarations rewritten by SPAG
!
!     ..
!     .. Array Arguments ..
!     ..
!
!  ================================================================
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
!     ==== Estimate optimal workspace. ====
!
      jw = MIN(Nw,Kbot-Ktop+1)
      IF ( jw<=2 ) THEN
         lwkopt = 1
      ELSE
!
!        ==== Workspace query call to DGEHRD ====
!
         CALL DGEHRD(jw,1,jw-1,T,Ldt,Work,Work,-1,info)
         lwk1 = INT(Work(1))
!
!        ==== Workspace query call to DORMHR ====
!
         CALL DORMHR('R','N',jw,jw,1,jw-1,T,Ldt,Work,V,Ldv,Work,-1,info)
         lwk2 = INT(Work(1))
!
!        ==== Optimal workspace ====
!
         lwkopt = jw + MAX(lwk1,lwk2)
      ENDIF
!
!     ==== Quick return in case of workspace query. ====
!
      IF ( Lwork==-1 ) THEN
         Work(1) = DBLE(lwkopt)
         RETURN
      ENDIF
!
!     ==== Nothing to do ...
!     ... for an empty active block ... ====
      Ns = 0
      Nd = 0
      Work(1) = ONE
      IF ( Ktop>Kbot ) RETURN
!     ... nor for an empty deflation window. ====
      IF ( Nw<1 ) RETURN
!
!     ==== Machine constants ====
!
      safmin = DLAMCH('SAFE MINIMUM')
      safmax = ONE/safmin
      CALL DLABAD(safmin,safmax)
      ulp = DLAMCH('PRECISION')
      smlnum = safmin*(DBLE(N)/ulp)
!
!     ==== Setup deflation window ====
!
      jw = MIN(Nw,Kbot-Ktop+1)
      kwtop = Kbot - jw + 1
      IF ( kwtop==Ktop ) THEN
         s = ZERO
      ELSE
         s = H(kwtop,kwtop-1)
      ENDIF
!
      IF ( Kbot==kwtop ) THEN
!
!        ==== 1-by-1 deflation window: not much to do ====
!
         Sr(kwtop) = H(kwtop,kwtop)
         Si(kwtop) = ZERO
         Ns = 1
         Nd = 0
         IF ( ABS(s)<=MAX(smlnum,ulp*ABS(H(kwtop,kwtop))) ) THEN
            Ns = 0
            Nd = 1
            IF ( kwtop>Ktop ) H(kwtop,kwtop-1) = ZERO
         ENDIF
         Work(1) = ONE
         RETURN
      ENDIF
!
!     ==== Convert to spike-triangular form.  (In case of a
!     .    rare QR failure, this routine continues to do
!     .    aggressive early deflation using that part of
!     .    the deflation window that converged using INFQR
!     .    here and there to keep track.) ====
!
      CALL DLACPY('U',jw,jw,H(kwtop,kwtop),Ldh,T,Ldt)
      CALL DCOPY(jw-1,H(kwtop+1,kwtop),Ldh+1,T(2,1),Ldt+1)
!
      CALL DLASET('A',jw,jw,ZERO,ONE,V,Ldv)
      CALL DLAHQR(.TRUE.,.TRUE.,jw,1,jw,T,Ldt,Sr(kwtop),Si(kwtop),1,jw, &
     &            V,Ldv,infqr)
!
!     ==== DTREXC needs a clean margin near the diagonal ====
!
      DO j = 1 , jw - 3
         T(j+2,j) = ZERO
         T(j+3,j) = ZERO
      ENDDO
      IF ( jw>2 ) T(jw,jw-2) = ZERO
!
!     ==== Deflation detection loop ====
!
      Ns = jw
      ilst = infqr + 1
      DO
         IF ( ilst<=Ns ) THEN
            IF ( Ns==1 ) THEN
               bulge = .FALSE.
            ELSE
               bulge = T(Ns,Ns-1)/=ZERO
            ENDIF
!
!        ==== Small spike tip test for deflation ====
!
            IF ( .NOT.bulge ) THEN
!
!           ==== Real eigenvalue ====
!
               foo = ABS(T(Ns,Ns))
               IF ( foo==ZERO ) foo = ABS(s)
               IF ( ABS(s*V(1,Ns))<=MAX(smlnum,ulp*foo) ) THEN
!
!              ==== Deflatable ====
!
                  Ns = Ns - 1
               ELSE
!
!              ==== Undeflatable.   Move it up out of the way.
!              .    (DTREXC can not fail in this case.) ====
!
                  ifst = Ns
                  CALL DTREXC('V',jw,T,Ldt,V,Ldv,ifst,ilst,Work,info)
                  ilst = ilst + 1
               ENDIF
            ELSE
!
!           ==== Complex conjugate pair ====
!
               foo = ABS(T(Ns,Ns)) + SQRT(ABS(T(Ns,Ns-1)))              &
     &               *SQRT(ABS(T(Ns-1,Ns)))
               IF ( foo==ZERO ) foo = ABS(s)
               IF ( MAX(ABS(s*V(1,Ns)),ABS(s*V(1,Ns-1)))                &
     &              <=MAX(smlnum,ulp*foo) ) THEN
!
!              ==== Deflatable ====
!
                  Ns = Ns - 2
               ELSE
!
!              ==== Undeflatable. Move them up out of the way.
!              .    Fortunately, DTREXC does the right thing with
!              .    ILST in case of a rare exchange failure. ====
!
                  ifst = Ns
                  CALL DTREXC('V',jw,T,Ldt,V,Ldv,ifst,ilst,Work,info)
                  ilst = ilst + 2
               ENDIF
            ENDIF
!
!        ==== End deflation detection loop ====
!
            CYCLE
         ENDIF
!
!        ==== Return to Hessenberg form ====
!
         IF ( Ns==0 ) s = ZERO
!
         IF ( Ns<jw ) THEN
!
!        ==== sorting diagonal blocks of T improves accuracy for
!        .    graded matrices.  Bubble sort deals well with
!        .    exchange failures. ====
!
            sorted = .FALSE.
            i = Ns + 1
            DO WHILE ( .NOT.(sorted) )
               sorted = .TRUE.
!
               kend = i - 1
               i = infqr + 1
               IF ( i==Ns ) THEN
                  k = i + 1
               ELSEIF ( T(i+1,i)==ZERO ) THEN
                  k = i + 1
               ELSE
                  k = i + 2
               ENDIF
               DO
                  IF ( k<=kend ) THEN
                     IF ( k==i+1 ) THEN
                        evi = ABS(T(i,i))
                     ELSE
                        evi = ABS(T(i,i)) + SQRT(ABS(T(i+1,i)))         &
     &                        *SQRT(ABS(T(i,i+1)))
                     ENDIF
!
                     IF ( k==kend ) THEN
                        evk = ABS(T(k,k))
                     ELSEIF ( T(k+1,k)==ZERO ) THEN
                        evk = ABS(T(k,k))
                     ELSE
                        evk = ABS(T(k,k)) + SQRT(ABS(T(k+1,k)))         &
     &                        *SQRT(ABS(T(k,k+1)))
                     ENDIF
!
                     IF ( evi>=evk ) THEN
                        i = k
                     ELSE
                        sorted = .FALSE.
                        ifst = i
                        ilst = k
                        CALL DTREXC('V',jw,T,Ldt,V,Ldv,ifst,ilst,Work,  &
     &                              info)
                        IF ( info==0 ) THEN
                           i = ilst
                        ELSE
                           i = k
                        ENDIF
                     ENDIF
                     IF ( i==kend ) THEN
                        k = i + 1
                     ELSEIF ( T(i+1,i)==ZERO ) THEN
                        k = i + 1
                     ELSE
                        k = i + 2
                     ENDIF
                     CYCLE
                  ENDIF
                  EXIT
               ENDDO
            ENDDO
         ENDIF
         EXIT
      ENDDO
!
!     ==== Restore shift/eigenvalue array from T ====
!
      i = jw
      DO
         IF ( i>=infqr+1 ) THEN
            IF ( i==infqr+1 ) THEN
               Sr(kwtop+i-1) = T(i,i)
               Si(kwtop+i-1) = ZERO
               i = i - 1
            ELSEIF ( T(i,i-1)==ZERO ) THEN
               Sr(kwtop+i-1) = T(i,i)
               Si(kwtop+i-1) = ZERO
               i = i - 1
            ELSE
               aa = T(i-1,i-1)
               cc = T(i,i-1)
               bb = T(i-1,i)
               dd = T(i,i)
               CALL DLANV2(aa,bb,cc,dd,Sr(kwtop+i-2),Si(kwtop+i-2),     &
     &                     Sr(kwtop+i-1),Si(kwtop+i-1),cs,sn)
               i = i - 2
            ENDIF
            CYCLE
         ENDIF
!
         IF ( Ns<jw .OR. s==ZERO ) THEN
            IF ( Ns>1 .AND. s/=ZERO ) THEN
!
!           ==== Reflect spike back into lower triangle ====
!
               CALL DCOPY(Ns,V,Ldv,Work,1)
               beta = Work(1)
               CALL DLARFG(Ns,beta,Work(2),1,tau)
               Work(1) = ONE
!
               CALL DLASET('L',jw-2,jw-2,ZERO,ZERO,T(3,1),Ldt)
!
               CALL DLARF('L',Ns,jw,Work,1,tau,T,Ldt,Work(jw+1))
               CALL DLARF('R',Ns,Ns,Work,1,tau,T,Ldt,Work(jw+1))
               CALL DLARF('R',jw,Ns,Work,1,tau,V,Ldv,Work(jw+1))
!
               CALL DGEHRD(jw,1,Ns,T,Ldt,Work,Work(jw+1),Lwork-jw,info)
            ENDIF
!
!        ==== Copy updated reduced window into place ====
!
            IF ( kwtop>1 ) H(kwtop,kwtop-1) = s*V(1,1)
            CALL DLACPY('U',jw,jw,T,Ldt,H(kwtop,kwtop),Ldh)
            CALL DCOPY(jw-1,T(2,1),Ldt+1,H(kwtop+1,kwtop),Ldh+1)
!
!        ==== Accumulate orthogonal matrix in order update
!        .    H and Z, if requested.  ====
!
            IF ( Ns>1 .AND. s/=ZERO )                                   &
     &           CALL DORMHR('R','N',jw,Ns,1,Ns,T,Ldt,Work,V,Ldv,       &
     &           Work(jw+1),Lwork-jw,info)
!
!        ==== Update vertical slab in H ====
!
            IF ( Wantt ) THEN
               ltop = 1
            ELSE
               ltop = Ktop
            ENDIF
            DO krow = ltop , kwtop - 1 , Nv
               kln = MIN(Nv,kwtop-krow)
               CALL DGEMM('N','N',kln,jw,jw,ONE,H(krow,kwtop),Ldh,V,Ldv,&
     &                    ZERO,Wv,Ldwv)
               CALL DLACPY('A',kln,jw,Wv,Ldwv,H(krow,kwtop),Ldh)
            ENDDO
!
!        ==== Update horizontal slab in H ====
!
            IF ( Wantt ) THEN
               DO kcol = Kbot + 1 , N , Nh
                  kln = MIN(Nh,N-kcol+1)
                  CALL DGEMM('C','N',jw,kln,jw,ONE,V,Ldv,H(kwtop,kcol), &
     &                       Ldh,ZERO,T,Ldt)
                  CALL DLACPY('A',jw,kln,T,Ldt,H(kwtop,kcol),Ldh)
               ENDDO
            ENDIF
!
!        ==== Update vertical slab in Z ====
!
            IF ( Wantz ) THEN
               DO krow = Iloz , Ihiz , Nv
                  kln = MIN(Nv,Ihiz-krow+1)
                  CALL DGEMM('N','N',kln,jw,jw,ONE,Z(krow,kwtop),Ldz,V, &
     &                       Ldv,ZERO,Wv,Ldwv)
                  CALL DLACPY('A',kln,jw,Wv,Ldwv,Z(krow,kwtop),Ldz)
               ENDDO
            ENDIF
         ENDIF
!
!     ==== Return the number of deflations ... ====
!
         Nd = jw - Ns
!
!     ==== ... and the number of shifts. (Subtracting
!     .    INFQR from the spike length takes care
!     .    of the case of a rare QR failure while
!     .    calculating eigenvalues of the deflation
!     .    window.)  ====
!
         Ns = Ns - infqr
!
!      ==== Return optimal workspace. ====
!
         Work(1) = DBLE(lwkopt)
         EXIT
      ENDDO
!
!     ==== End of DLAQR2 ====
!
      END SUBROUTINE DLAQR2
