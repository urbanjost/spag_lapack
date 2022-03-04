!*==zlaqr3.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b ZLAQR3 performs the unitary similarity transformation of a Hessenberg matrix to detect and deflate fully converged eigenvalues from a trailing principal submatrix (aggressive early deflation).
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZLAQR3 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlaqr3.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlaqr3.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlaqr3.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZLAQR3( WANTT, WANTZ, N, KTOP, KBOT, NW, H, LDH, ILOZ,
!                          IHIZ, Z, LDZ, NS, ND, SH, V, LDV, NH, T, LDT,
!                          NV, WV, LDWV, WORK, LWORK )
!
!       .. Scalar Arguments ..
!       INTEGER            IHIZ, ILOZ, KBOT, KTOP, LDH, LDT, LDV, LDWV,
!      $                   LDZ, LWORK, N, ND, NH, NS, NV, NW
!       LOGICAL            WANTT, WANTZ
!       ..
!       .. Array Arguments ..
!       COMPLEX*16         H( LDH, * ), SH( * ), T( LDT, * ), V( LDV, * ),
!      $                   WORK( * ), WV( LDWV, * ), Z( LDZ, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    Aggressive early deflation:
!>
!>    ZLAQR3 accepts as input an upper Hessenberg matrix
!>    H and performs an unitary similarity transformation
!>    designed to detect and deflate fully converged eigenvalues from
!>    a trailing principal submatrix.  On output H has been over-
!>    written by a new Hessenberg matrix that is a perturbation of
!>    an unitary similarity transformation of H.  It is to be
!>    hoped that the final version of H has many zero subdiagonal
!>    entries.
!>
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] WANTT
!> \verbatim
!>          WANTT is LOGICAL
!>          If .TRUE., then the Hessenberg matrix H is fully updated
!>          so that the triangular Schur factor may be
!>          computed (in cooperation with the calling subroutine).
!>          If .FALSE., then only enough of H is updated to preserve
!>          the eigenvalues.
!> \endverbatim
!>
!> \param[in] WANTZ
!> \verbatim
!>          WANTZ is LOGICAL
!>          If .TRUE., then the unitary matrix Z is updated so
!>          so that the unitary Schur factor may be computed
!>          (in cooperation with the calling subroutine).
!>          If .FALSE., then Z is not referenced.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix H and (if WANTZ is .TRUE.) the
!>          order of the unitary matrix Z.
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
!>          H is COMPLEX*16 array, dimension (LDH,N)
!>          On input the initial N-by-N section of H stores the
!>          Hessenberg matrix undergoing aggressive early deflation.
!>          On output H has been transformed by a unitary
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
!>          Z is COMPLEX*16 array, dimension (LDZ,N)
!>          IF WANTZ is .TRUE., then on output, the unitary
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
!> \param[out] SH
!> \verbatim
!>          SH is COMPLEX*16 array, dimension (KBOT)
!>          On output, approximate eigenvalues that may
!>          be used for shifts are stored in SH(KBOT-ND-NS+1)
!>          through SR(KBOT-ND).  Converged eigenvalues are
!>          stored in SH(KBOT-ND+1) through SH(KBOT).
!> \endverbatim
!>
!> \param[out] V
!> \verbatim
!>          V is COMPLEX*16 array, dimension (LDV,NW)
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
!>          T is COMPLEX*16 array, dimension (LDT,NW)
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
!>          WV is COMPLEX*16 array, dimension (LDWV,NW)
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
!>          WORK is COMPLEX*16 array, dimension (LWORK)
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
!>          If LWORK = -1, then a workspace query is assumed; ZLAQR3
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
!> \date June 2016
!
!> \ingroup complex16OTHERauxiliary
!
!> \par Contributors:
!  ==================
!>
!>       Karen Braman and Ralph Byers, Department of Mathematics,
!>       University of Kansas, USA
!>
!  =====================================================================
      SUBROUTINE ZLAQR3(Wantt,Wantz,N,Ktop,Kbot,Nw,H,Ldh,Iloz,Ihiz,Z,   &
     &                  Ldz,Ns,Nd,Sh,V,Ldv,Nh,T,Ldt,Nv,Wv,Ldwv,Work,    &
     &                  Lwork)
      USE F77KINDS                        
      USE S_DLABAD
      USE S_DLAMCH
      USE S_ILAENV
      USE S_ZCOPY
      USE S_ZGEHRD
      USE S_ZGEMM
      USE S_ZLACPY
      USE S_ZLAHQR
      USE S_ZLAQR4
      USE S_ZLARF
      USE S_ZLARFG
      USE S_ZLASET
      USE S_ZTREXC
      USE S_ZUNMHR
      IMPLICIT NONE
!*--ZLAQR3286
!
! PARAMETER definitions rewritten by SPAG
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ZERO = (0.0D0,0.0D0) ,         &
     &                 ONE = (1.0D0,0.0D0)
      REAL(R8KIND) , PARAMETER  ::  RZERO = 0.0D0 , RONE = 1.0D0
!
! Dummy argument declarations rewritten by SPAG
!
      LOGICAL , INTENT(IN) :: Wantt
      LOGICAL , INTENT(IN) :: Wantz
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: Ktop
      INTEGER , INTENT(IN) :: Kbot
      INTEGER , INTENT(IN) :: Nw
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldh,*) :: H
      INTEGER :: Ldh
      INTEGER , INTENT(IN) :: Iloz
      INTEGER , INTENT(IN) :: Ihiz
      COMPLEX(CX16KIND) , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      INTEGER , INTENT(INOUT) :: Ns
      INTEGER , INTENT(OUT) :: Nd
      COMPLEX(CX16KIND) , DIMENSION(*) :: Sh
      COMPLEX(CX16KIND) , DIMENSION(Ldv,*) :: V
      INTEGER :: Ldv
      INTEGER , INTENT(IN) :: Nh
      COMPLEX(CX16KIND) , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      INTEGER , INTENT(IN) :: Nv
      COMPLEX(CX16KIND) , DIMENSION(Ldwv,*) :: Wv
      INTEGER :: Ldwv
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER :: Lwork
!
! Local variable declarations rewritten by SPAG
!
      COMPLEX(CX16KIND) :: beta , cdum , s , tau
      REAL(R8KIND) :: CABS1
      REAL(R8KIND) :: foo , safmax , safmin , smlnum , ulp
      INTEGER :: i , ifst , ilst , info , infqr , j , jw , kcol , kln , &
     &           knt , krow , kwtop , ltop , lwk1 , lwk2 , lwk3 ,       &
     &           lwkopt , nmin
!
! End of declarations rewritten by SPAG
!
!     ..
!     .. Array Arguments ..
!     ..
!
!  ================================================================
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
      CABS1(cdum) = ABS(DBLE(cdum)) + ABS(DIMAG(cdum))
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
!        ==== Workspace query call to ZGEHRD ====
!
         CALL ZGEHRD(jw,1,jw-1,T,Ldt,Work,Work,-1,info)
         lwk1 = INT(Work(1))
!
!        ==== Workspace query call to ZUNMHR ====
!
         CALL ZUNMHR('R','N',jw,jw,1,jw-1,T,Ldt,Work,V,Ldv,Work,-1,info)
         lwk2 = INT(Work(1))
!
!        ==== Workspace query call to ZLAQR4 ====
!
         CALL ZLAQR4(.TRUE.,.TRUE.,jw,1,jw,T,Ldt,Sh,1,jw,V,Ldv,Work,-1, &
     &               infqr)
         lwk3 = INT(Work(1))
!
!        ==== Optimal workspace ====
!
         lwkopt = MAX(jw+MAX(lwk1,lwk2),lwk3)
      ENDIF
!
!     ==== Quick return in case of workspace query. ====
!
      IF ( Lwork==-1 ) THEN
         Work(1) = DCMPLX(lwkopt,0)
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
      safmax = RONE/safmin
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
         Sh(kwtop) = H(kwtop,kwtop)
         Ns = 1
         Nd = 0
         IF ( CABS1(s)<=MAX(smlnum,ulp*CABS1(H(kwtop,kwtop))) ) THEN
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
      CALL ZLACPY('U',jw,jw,H(kwtop,kwtop),Ldh,T,Ldt)
      CALL ZCOPY(jw-1,H(kwtop+1,kwtop),Ldh+1,T(2,1),Ldt+1)
!
      CALL ZLASET('A',jw,jw,ZERO,ONE,V,Ldv)
      nmin = ILAENV(12,'ZLAQR3','SV',jw,1,jw,Lwork)
      IF ( jw>nmin ) THEN
         CALL ZLAQR4(.TRUE.,.TRUE.,jw,1,jw,T,Ldt,Sh(kwtop),1,jw,V,Ldv,  &
     &               Work,Lwork,infqr)
      ELSE
         CALL ZLAHQR(.TRUE.,.TRUE.,jw,1,jw,T,Ldt,Sh(kwtop),1,jw,V,Ldv,  &
     &               infqr)
      ENDIF
!
!     ==== Deflation detection loop ====
!
      Ns = jw
      ilst = infqr + 1
      DO knt = infqr + 1 , jw
!
!        ==== Small spike tip deflation test ====
!
         foo = CABS1(T(Ns,Ns))
         IF ( foo==RZERO ) foo = CABS1(s)
         IF ( CABS1(s)*CABS1(V(1,Ns))<=MAX(smlnum,ulp*foo) ) THEN
!
!           ==== One more converged eigenvalue ====
!
            Ns = Ns - 1
         ELSE
!
!           ==== One undeflatable eigenvalue.  Move it up out of the
!           .    way.   (ZTREXC can not fail in this case.) ====
!
            ifst = Ns
            CALL ZTREXC('V',jw,T,Ldt,V,Ldv,ifst,ilst,info)
            ilst = ilst + 1
         ENDIF
      ENDDO
!
!        ==== Return to Hessenberg form ====
!
      IF ( Ns==0 ) s = ZERO
!
      IF ( Ns<jw ) THEN
!
!        ==== sorting the diagonal of T improves accuracy for
!        .    graded matrices.  ====
!
         DO i = infqr + 1 , Ns
            ifst = i
            DO j = i + 1 , Ns
               IF ( CABS1(T(j,j))>CABS1(T(ifst,ifst)) ) ifst = j
            ENDDO
            ilst = i
            IF ( ifst/=ilst ) CALL ZTREXC('V',jw,T,Ldt,V,Ldv,ifst,ilst, &
     &           info)
         ENDDO
      ENDIF
!
!     ==== Restore shift/eigenvalue array from T ====
!
      DO i = infqr + 1 , jw
         Sh(kwtop+i-1) = T(i,i)
      ENDDO
!
!
      IF ( Ns<jw .OR. s==ZERO ) THEN
         IF ( Ns>1 .AND. s/=ZERO ) THEN
!
!           ==== Reflect spike back into lower triangle ====
!
            CALL ZCOPY(Ns,V,Ldv,Work,1)
            DO i = 1 , Ns
               Work(i) = DCONJG(Work(i))
            ENDDO
            beta = Work(1)
            CALL ZLARFG(Ns,beta,Work(2),1,tau)
            Work(1) = ONE
!
            CALL ZLASET('L',jw-2,jw-2,ZERO,ZERO,T(3,1),Ldt)
!
            CALL ZLARF('L',Ns,jw,Work,1,DCONJG(tau),T,Ldt,Work(jw+1))
            CALL ZLARF('R',Ns,Ns,Work,1,tau,T,Ldt,Work(jw+1))
            CALL ZLARF('R',jw,Ns,Work,1,tau,V,Ldv,Work(jw+1))
!
            CALL ZGEHRD(jw,1,Ns,T,Ldt,Work,Work(jw+1),Lwork-jw,info)
         ENDIF
!
!        ==== Copy updated reduced window into place ====
!
         IF ( kwtop>1 ) H(kwtop,kwtop-1) = s*DCONJG(V(1,1))
         CALL ZLACPY('U',jw,jw,T,Ldt,H(kwtop,kwtop),Ldh)
         CALL ZCOPY(jw-1,T(2,1),Ldt+1,H(kwtop+1,kwtop),Ldh+1)
!
!        ==== Accumulate orthogonal matrix in order update
!        .    H and Z, if requested.  ====
!
         IF ( Ns>1 .AND. s/=ZERO ) CALL ZUNMHR('R','N',jw,Ns,1,Ns,T,Ldt,&
     &        Work,V,Ldv,Work(jw+1),Lwork-jw,info)
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
            CALL ZGEMM('N','N',kln,jw,jw,ONE,H(krow,kwtop),Ldh,V,Ldv,   &
     &                 ZERO,Wv,Ldwv)
            CALL ZLACPY('A',kln,jw,Wv,Ldwv,H(krow,kwtop),Ldh)
         ENDDO
!
!        ==== Update horizontal slab in H ====
!
         IF ( Wantt ) THEN
            DO kcol = Kbot + 1 , N , Nh
               kln = MIN(Nh,N-kcol+1)
               CALL ZGEMM('C','N',jw,kln,jw,ONE,V,Ldv,H(kwtop,kcol),Ldh,&
     &                    ZERO,T,Ldt)
               CALL ZLACPY('A',jw,kln,T,Ldt,H(kwtop,kcol),Ldh)
            ENDDO
         ENDIF
!
!        ==== Update vertical slab in Z ====
!
         IF ( Wantz ) THEN
            DO krow = Iloz , Ihiz , Nv
               kln = MIN(Nv,Ihiz-krow+1)
               CALL ZGEMM('N','N',kln,jw,jw,ONE,Z(krow,kwtop),Ldz,V,Ldv,&
     &                    ZERO,Wv,Ldwv)
               CALL ZLACPY('A',kln,jw,Wv,Ldwv,Z(krow,kwtop),Ldz)
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
      Work(1) = DCMPLX(lwkopt,0)
!
!     ==== End of ZLAQR3 ====
!
      END SUBROUTINE ZLAQR3
