!*==claqr0.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b CLAQR0 computes the eigenvalues of a Hessenberg matrix, and optionally the matrices from the Schur decomposition.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CLAQR0 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/claqr0.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/claqr0.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/claqr0.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CLAQR0( WANTT, WANTZ, N, ILO, IHI, H, LDH, W, ILOZ,
!                          IHIZ, Z, LDZ, WORK, LWORK, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            IHI, IHIZ, ILO, ILOZ, INFO, LDH, LDZ, LWORK, N
!       LOGICAL            WANTT, WANTZ
!       ..
!       .. Array Arguments ..
!       COMPLEX            H( LDH, * ), W( * ), WORK( * ), Z( LDZ, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    CLAQR0 computes the eigenvalues of a Hessenberg matrix H
!>    and, optionally, the matrices T and Z from the Schur decomposition
!>    H = Z T Z**H, where T is an upper triangular matrix (the
!>    Schur form), and Z is the unitary matrix of Schur vectors.
!>
!>    Optionally Z may be postmultiplied into an input unitary
!>    matrix Q so that this routine can give the Schur factorization
!>    of a matrix A which has been reduced to the Hessenberg form H
!>    by the unitary matrix Q:  A = Q*H*Q**H = (QZ)*H*(QZ)**H.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] WANTT
!> \verbatim
!>          WANTT is LOGICAL
!>          = .TRUE. : the full Schur form T is required;
!>          = .FALSE.: only eigenvalues are required.
!> \endverbatim
!>
!> \param[in] WANTZ
!> \verbatim
!>          WANTZ is LOGICAL
!>          = .TRUE. : the matrix of Schur vectors Z is required;
!>          = .FALSE.: Schur vectors are not required.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>           The order of the matrix H.  N >= 0.
!> \endverbatim
!>
!> \param[in] ILO
!> \verbatim
!>          ILO is INTEGER
!> \endverbatim
!>
!> \param[in] IHI
!> \verbatim
!>          IHI is INTEGER
!>           It is assumed that H is already upper triangular in rows
!>           and columns 1:ILO-1 and IHI+1:N and, if ILO > 1,
!>           H(ILO,ILO-1) is zero. ILO and IHI are normally set by a
!>           previous call to CGEBAL, and then passed to CGEHRD when the
!>           matrix output by CGEBAL is reduced to Hessenberg form.
!>           Otherwise, ILO and IHI should be set to 1 and N,
!>           respectively.  If N > 0, then 1 <= ILO <= IHI <= N.
!>           If N = 0, then ILO = 1 and IHI = 0.
!> \endverbatim
!>
!> \param[in,out] H
!> \verbatim
!>          H is COMPLEX array, dimension (LDH,N)
!>           On entry, the upper Hessenberg matrix H.
!>           On exit, if INFO = 0 and WANTT is .TRUE., then H
!>           contains the upper triangular matrix T from the Schur
!>           decomposition (the Schur form). If INFO = 0 and WANT is
!>           .FALSE., then the contents of H are unspecified on exit.
!>           (The output value of H when INFO > 0 is given under the
!>           description of INFO below.)
!>
!>           This subroutine may explicitly set H(i,j) = 0 for i > j and
!>           j = 1, 2, ... ILO-1 or j = IHI+1, IHI+2, ... N.
!> \endverbatim
!>
!> \param[in] LDH
!> \verbatim
!>          LDH is INTEGER
!>           The leading dimension of the array H. LDH >= max(1,N).
!> \endverbatim
!>
!> \param[out] W
!> \verbatim
!>          W is COMPLEX array, dimension (N)
!>           The computed eigenvalues of H(ILO:IHI,ILO:IHI) are stored
!>           in W(ILO:IHI). If WANTT is .TRUE., then the eigenvalues are
!>           stored in the same order as on the diagonal of the Schur
!>           form returned in H, with W(i) = H(i,i).
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
!>           Specify the rows of Z to which transformations must be
!>           applied if WANTZ is .TRUE..
!>           1 <= ILOZ <= ILO; IHI <= IHIZ <= N.
!> \endverbatim
!>
!> \param[in,out] Z
!> \verbatim
!>          Z is COMPLEX array, dimension (LDZ,IHI)
!>           If WANTZ is .FALSE., then Z is not referenced.
!>           If WANTZ is .TRUE., then Z(ILO:IHI,ILOZ:IHIZ) is
!>           replaced by Z(ILO:IHI,ILOZ:IHIZ)*U where U is the
!>           orthogonal Schur factor of H(ILO:IHI,ILO:IHI).
!>           (The output value of Z when INFO > 0 is given under
!>           the description of INFO below.)
!> \endverbatim
!>
!> \param[in] LDZ
!> \verbatim
!>          LDZ is INTEGER
!>           The leading dimension of the array Z.  if WANTZ is .TRUE.
!>           then LDZ >= MAX(1,IHIZ).  Otherwise, LDZ >= 1.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX array, dimension LWORK
!>           On exit, if LWORK = -1, WORK(1) returns an estimate of
!>           the optimal value for LWORK.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>           The dimension of the array WORK.  LWORK >= max(1,N)
!>           is sufficient, but LWORK typically as large as 6*N may
!>           be required for optimal performance.  A workspace query
!>           to determine the optimal workspace size is recommended.
!>
!>           If LWORK = -1, then CLAQR0 does a workspace query.
!>           In this case, CLAQR0 checks the input parameters and
!>           estimates the optimal workspace size for the given
!>           values of N, ILO and IHI.  The estimate is returned
!>           in WORK(1).  No error message related to LWORK is
!>           issued by XERBLA.  Neither H nor Z are accessed.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>             = 0:  successful exit
!>             > 0:  if INFO = i, CLAQR0 failed to compute all of
!>                the eigenvalues.  Elements 1:ilo-1 and i+1:n of WR
!>                and WI contain those eigenvalues which have been
!>                successfully computed.  (Failures are rare.)
!>
!>                If INFO > 0 and WANT is .FALSE., then on exit,
!>                the remaining unconverged eigenvalues are the eigen-
!>                values of the upper Hessenberg matrix rows and
!>                columns ILO through INFO of the final, output
!>                value of H.
!>
!>                If INFO > 0 and WANTT is .TRUE., then on exit
!>
!>           (*)  (initial value of H)*U  = U*(final value of H)
!>
!>                where U is a unitary matrix.  The final
!>                value of  H is upper Hessenberg and triangular in
!>                rows and columns INFO+1 through IHI.
!>
!>                If INFO > 0 and WANTZ is .TRUE., then on exit
!>
!>                  (final value of Z(ILO:IHI,ILOZ:IHIZ)
!>                   =  (initial value of Z(ILO:IHI,ILOZ:IHIZ)*U
!>
!>                where U is the unitary matrix in (*) (regard-
!>                less of the value of WANTT.)
!>
!>                If INFO > 0 and WANTZ is .FALSE., then Z is not
!>                accessed.
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
!> \ingroup complexOTHERauxiliary
!
!> \par Contributors:
!  ==================
!>
!>       Karen Braman and Ralph Byers, Department of Mathematics,
!>       University of Kansas, USA
!
!> \par References:
!  ================
!>
!>       K. Braman, R. Byers and R. Mathias, The Multi-Shift QR
!>       Algorithm Part I: Maintaining Well Focused Shifts, and Level 3
!>       Performance, SIAM Journal of Matrix Analysis, volume 23, pages
!>       929--947, 2002.
!> \n
!>       K. Braman, R. Byers and R. Mathias, The Multi-Shift QR
!>       Algorithm Part II: Aggressive Early Deflation, SIAM Journal
!>       of Matrix Analysis, volume 23, pages 948--973, 2002.
!>
!  =====================================================================
      SUBROUTINE CLAQR0(Wantt,Wantz,N,Ilo,Ihi,H,Ldh,W,Iloz,Ihiz,Z,Ldz,  &
     &                  Work,Lwork,Info)
      USE S_CLACPY
      USE S_CLAHQR
      USE S_CLAQR3
      USE S_CLAQR4
      USE S_CLAQR5
      USE S_ILAENV
      IMPLICIT NONE
!*--CLAQR0250
!
! PARAMETER definitions rewritten by SPAG
!
      INTEGER , PARAMETER  ::  NTINY = 15 , KEXNW = 5 , KEXSH = 6
      REAL , PARAMETER  ::  WILK1 = 0.75E0
      COMPLEX , PARAMETER  ::  ZERO = (0.0E0,0.0E0) ,                   &
     &                         ONE = (1.0E0,0.0E0)
      REAL , PARAMETER  ::  TWO = 2.0E0
!
! Dummy argument declarations rewritten by SPAG
!
      LOGICAL :: Wantt
      LOGICAL :: Wantz
      INTEGER :: N
      INTEGER :: Ilo
      INTEGER :: Ihi
      COMPLEX , DIMENSION(Ldh,*) :: H
      INTEGER :: Ldh
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: W
      INTEGER :: Iloz
      INTEGER :: Ihiz
      COMPLEX , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER :: Info
!
! Local variable declarations rewritten by SPAG
!
      COMPLEX :: aa , bb , cc , cdum , dd , det , rtdisc , swap , tr2
      REAL :: CABS1
      INTEGER :: i , inf , it , itmax , k , kacc22 , kbot , kdu , ks ,  &
     &           kt , ktop , ku , kv , kwh , kwtop , kwv , ld , ls ,    &
     &           lwkopt , ndec , ndfl , nh , nho , nibble , nmin , ns , &
     &           nsmax , nsr , nve , nw , nwmax , nwr , nwupbd
      CHARACTER(2) :: jbcmpz
      REAL :: s
      LOGICAL :: sorted
      COMPLEX , DIMENSION(1,1) :: zdum
!
! End of declarations rewritten by SPAG
!
!     ..
!     .. Array Arguments ..
!     ..
!
!  ================================================================
!     .. Parameters ..
!
!     ==== Matrices of order NTINY or smaller must be processed by
!     .    CLAHQR because of insufficient subdiagonal scratch space.
!     .    (This is a hard limit.) ====
!
!     ==== Exceptional deflation windows:  try to cure rare
!     .    slow convergence by varying the size of the
!     .    deflation window after KEXNW iterations. ====
!
!     ==== Exceptional shifts: try to cure rare slow convergence
!     .    with ad-hoc exceptional shifts every KEXSH iterations.
!     .    ====
!
!     ==== The constant WILK1 is used to form the exceptional
!     .    shifts. ====
!     ..
!     .. Local Scalars ..
!     ..
!     .. External Functions ..
!     ..
!     .. Local Arrays ..
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
      Info = 0
!
!     ==== Quick return for N = 0: nothing to do. ====
!
      IF ( N==0 ) THEN
         Work(1) = ONE
         RETURN
      ENDIF
!
      IF ( N<=NTINY ) THEN
!
!        ==== Tiny matrices must use CLAHQR. ====
!
         lwkopt = 1
         IF ( Lwork/=-1 ) CALL CLAHQR(Wantt,Wantz,N,Ilo,Ihi,H,Ldh,W,    &
     &                                Iloz,Ihiz,Z,Ldz,Info)
      ELSE
!
!        ==== Use small bulge multi-shift QR with aggressive early
!        .    deflation on larger-than-tiny matrices. ====
!
!        ==== Hope for the best. ====
!
         Info = 0
!
!        ==== Set up job flags for ILAENV. ====
!
         IF ( Wantt ) THEN
            jbcmpz(1:1) = 'S'
         ELSE
            jbcmpz(1:1) = 'E'
         ENDIF
         IF ( Wantz ) THEN
            jbcmpz(2:2) = 'V'
         ELSE
            jbcmpz(2:2) = 'N'
         ENDIF
!
!        ==== NWR = recommended deflation window size.  At this
!        .    point,  N .GT. NTINY = 15, so there is enough
!        .    subdiagonal workspace for NWR.GE.2 as required.
!        .    (In fact, there is enough subdiagonal space for
!        .    NWR.GE.4.) ====
!
         nwr = ILAENV(13,'CLAQR0',jbcmpz,N,Ilo,Ihi,Lwork)
         nwr = MAX(2,nwr)
         nwr = MIN(Ihi-Ilo+1,(N-1)/3,nwr)
!
!        ==== NSR = recommended number of simultaneous shifts.
!        .    At this point N .GT. NTINY = 15, so there is at
!        .    enough subdiagonal workspace for NSR to be even
!        .    and greater than or equal to two as required. ====
!
         nsr = ILAENV(15,'CLAQR0',jbcmpz,N,Ilo,Ihi,Lwork)
         nsr = MIN(nsr,(N-3)/6,Ihi-Ilo)
         nsr = MAX(2,nsr-MOD(nsr,2))
!
!        ==== Estimate optimal workspace ====
!
!        ==== Workspace query call to CLAQR3 ====
!
         CALL CLAQR3(Wantt,Wantz,N,Ilo,Ihi,nwr+1,H,Ldh,Iloz,Ihiz,Z,Ldz, &
     &               ls,ld,W,H,Ldh,N,H,Ldh,N,H,Ldh,Work,-1)
!
!        ==== Optimal workspace = MAX(CLAQR5, CLAQR3) ====
!
         lwkopt = MAX(3*nsr/2,INT(Work(1)))
!
!        ==== Quick return in case of workspace query. ====
!
         IF ( Lwork==-1 ) THEN
            Work(1) = CMPLX(lwkopt,0)
            RETURN
         ENDIF
!
!        ==== CLAHQR/CLAQR0 crossover point ====
!
         nmin = ILAENV(12,'CLAQR0',jbcmpz,N,Ilo,Ihi,Lwork)
         nmin = MAX(NTINY,nmin)
!
!        ==== Nibble crossover point ====
!
         nibble = ILAENV(14,'CLAQR0',jbcmpz,N,Ilo,Ihi,Lwork)
         nibble = MAX(0,nibble)
!
!        ==== Accumulate reflections during ttswp?  Use block
!        .    2-by-2 structure during matrix-matrix multiply? ====
!
         kacc22 = ILAENV(16,'CLAQR0',jbcmpz,N,Ilo,Ihi,Lwork)
         kacc22 = MAX(0,kacc22)
         kacc22 = MIN(2,kacc22)
!
!        ==== NWMAX = the largest possible deflation window for
!        .    which there is sufficient workspace. ====
!
         nwmax = MIN((N-1)/3,Lwork/2)
         nw = nwmax
!
!        ==== NSMAX = the Largest number of simultaneous shifts
!        .    for which there is sufficient workspace. ====
!
         nsmax = MIN((N-3)/6,2*Lwork/3)
         nsmax = nsmax - MOD(nsmax,2)
!
!        ==== NDFL: an iteration count restarted at deflation. ====
!
         ndfl = 1
!
!        ==== ITMAX = iteration limit ====
!
         itmax = MAX(30,2*KEXSH)*MAX(10,(Ihi-Ilo+1))
!
!        ==== Last row and column in the active block ====
!
         kbot = Ihi
!
!        ==== Main Loop ====
!
         DO it = 1 , itmax
!
!           ==== Done when KBOT falls below ILO ====
!
            IF ( kbot<Ilo ) GOTO 100
!
!           ==== Locate active block ====
!
            DO k = kbot , Ilo + 1 , -1
               IF ( H(k,k-1)==ZERO ) GOTO 20
            ENDDO
            k = Ilo
 20         ktop = k
!
!           ==== Select deflation window size:
!           .    Typical Case:
!           .      If possible and advisable, nibble the entire
!           .      active block.  If not, use size MIN(NWR,NWMAX)
!           .      or MIN(NWR+1,NWMAX) depending upon which has
!           .      the smaller corresponding subdiagonal entry
!           .      (a heuristic).
!           .
!           .    Exceptional Case:
!           .      If there have been no deflations in KEXNW or
!           .      more iterations, then vary the deflation window
!           .      size.   At first, because, larger windows are,
!           .      in general, more powerful than smaller ones,
!           .      rapidly increase the window to the maximum possible.
!           .      Then, gradually reduce the window size. ====
!
            nh = kbot - ktop + 1
            nwupbd = MIN(nh,nwmax)
            IF ( ndfl<KEXNW ) THEN
               nw = MIN(nwupbd,nwr)
            ELSE
               nw = MIN(nwupbd,2*nw)
            ENDIF
            IF ( nw<nwmax ) THEN
               IF ( nw>=nh-1 ) THEN
                  nw = nh
               ELSE
                  kwtop = kbot - nw + 1
                  IF ( CABS1(H(kwtop,kwtop-1))>CABS1(H(kwtop-1,kwtop-2))&
     &                 ) nw = nw + 1
               ENDIF
            ENDIF
            IF ( ndfl<KEXNW ) THEN
               ndec = -1
            ELSEIF ( ndec>=0 .OR. nw>=nwupbd ) THEN
               ndec = ndec + 1
               IF ( nw-ndec<2 ) ndec = 0
               nw = nw - ndec
            ENDIF
!
!           ==== Aggressive early deflation:
!           .    split workspace under the subdiagonal into
!           .      - an nw-by-nw work array V in the lower
!           .        left-hand-corner,
!           .      - an NW-by-at-least-NW-but-more-is-better
!           .        (NW-by-NHO) horizontal work array along
!           .        the bottom edge,
!           .      - an at-least-NW-but-more-is-better (NHV-by-NW)
!           .        vertical work array along the left-hand-edge.
!           .        ====
!
            kv = N - nw + 1
            kt = nw + 1
            nho = (N-nw-1) - kt + 1
            kwv = nw + 2
            nve = (N-nw) - kwv + 1
!
!           ==== Aggressive early deflation ====
!
            CALL CLAQR3(Wantt,Wantz,N,ktop,kbot,nw,H,Ldh,Iloz,Ihiz,Z,   &
     &                  Ldz,ls,ld,W,H(kv,1),Ldh,nho,H(kv,kt),Ldh,nve,   &
     &                  H(kwv,1),Ldh,Work,Lwork)
!
!           ==== Adjust KBOT accounting for new deflations. ====
!
            kbot = kbot - ld
!
!           ==== KS points to the shifts. ====
!
            ks = kbot - ls + 1
!
!           ==== Skip an expensive QR sweep if there is a (partly
!           .    heuristic) reason to expect that many eigenvalues
!           .    will deflate without it.  Here, the QR sweep is
!           .    skipped if many eigenvalues have just been deflated
!           .    or if the remaining active block is small.
!
            IF ( (ld==0) .OR.                                           &
     &           ((100*ld<=nw*nibble) .AND. (kbot-ktop+1>MIN(nmin,nwmax)&
     &           )) ) THEN
!
!              ==== NS = nominal number of simultaneous shifts.
!              .    This may be lowered (slightly) if CLAQR3
!              .    did not provide that many shifts. ====
!
               ns = MIN(nsmax,nsr,MAX(2,kbot-ktop))
               ns = ns - MOD(ns,2)
!
!              ==== If there have been no deflations
!              .    in a multiple of KEXSH iterations,
!              .    then try exceptional shifts.
!              .    Otherwise use shifts provided by
!              .    CLAQR3 above or from the eigenvalues
!              .    of a trailing principal submatrix. ====
!
               IF ( MOD(ndfl,KEXSH)==0 ) THEN
                  ks = kbot - ns + 1
                  DO i = kbot , ks + 1 , -2
                     W(i) = H(i,i) + WILK1*CABS1(H(i,i-1))
                     W(i-1) = W(i)
                  ENDDO
               ELSE
!
!                 ==== Got NS/2 or fewer shifts? Use CLAQR4 or
!                 .    CLAHQR on a trailing principal submatrix to
!                 .    get more. (Since NS.LE.NSMAX.LE.(N+6)/9,
!                 .    there is enough space below the subdiagonal
!                 .    to fit an NS-by-NS scratch array.) ====
!
                  IF ( kbot-ks+1<=ns/2 ) THEN
                     ks = kbot - ns + 1
                     kt = N - ns + 1
                     CALL CLACPY('A',ns,ns,H(ks,ks),Ldh,H(kt,1),Ldh)
                     IF ( ns>nmin ) THEN
                        CALL CLAQR4(.FALSE.,.FALSE.,ns,1,ns,H(kt,1),Ldh,&
     &                              W(ks),1,1,zdum,1,Work,Lwork,inf)
                     ELSE
                        CALL CLAHQR(.FALSE.,.FALSE.,ns,1,ns,H(kt,1),Ldh,&
     &                              W(ks),1,1,zdum,1,inf)
                     ENDIF
                     ks = ks + inf
!
!                    ==== In case of a rare QR failure use
!                    .    eigenvalues of the trailing 2-by-2
!                    .    principal submatrix.  Scale to avoid
!                    .    overflows, underflows and subnormals.
!                    .    (The scale factor S can not be zero,
!                    .    because H(KBOT,KBOT-1) is nonzero.) ====
!
                     IF ( ks>=kbot ) THEN
                        s = CABS1(H(kbot-1,kbot-1))                     &
     &                      + CABS1(H(kbot,kbot-1))                     &
     &                      + CABS1(H(kbot-1,kbot))                     &
     &                      + CABS1(H(kbot,kbot))
                        aa = H(kbot-1,kbot-1)/s
                        cc = H(kbot,kbot-1)/s
                        bb = H(kbot-1,kbot)/s
                        dd = H(kbot,kbot)/s
                        tr2 = (aa+dd)/TWO
                        det = (aa-tr2)*(dd-tr2) - bb*cc
                        rtdisc = SQRT(-det)
                        W(kbot-1) = (tr2+rtdisc)*s
                        W(kbot) = (tr2-rtdisc)*s
!
                        ks = kbot - 1
                     ENDIF
                  ENDIF
!
                  IF ( kbot-ks+1>ns ) THEN
!
!                    ==== Sort the shifts (Helps a little) ====
!
                     sorted = .FALSE.
                     DO k = kbot , ks + 1 , -1
                        IF ( sorted ) EXIT
                        sorted = .TRUE.
                        DO i = ks , k - 1
                           IF ( CABS1(W(i))<CABS1(W(i+1)) ) THEN
                              sorted = .FALSE.
                              swap = W(i)
                              W(i) = W(i+1)
                              W(i+1) = swap
                           ENDIF
                        ENDDO
                     ENDDO
                  ENDIF
               ENDIF
!
!              ==== If there are only two shifts, then use
!              .    only one.  ====
!
               IF ( kbot-ks+1==2 ) THEN
                  IF ( CABS1(W(kbot)-H(kbot,kbot))                      &
     &                 <CABS1(W(kbot-1)-H(kbot,kbot)) ) THEN
                     W(kbot-1) = W(kbot)
                  ELSE
                     W(kbot) = W(kbot-1)
                  ENDIF
               ENDIF
!
!              ==== Use up to NS of the the smallest magnitude
!              .    shifts.  If there aren't NS shifts available,
!              .    then use them all, possibly dropping one to
!              .    make the number of shifts even. ====
!
               ns = MIN(ns,kbot-ks+1)
               ns = ns - MOD(ns,2)
               ks = kbot - ns + 1
!
!              ==== Small-bulge multi-shift QR sweep:
!              .    split workspace under the subdiagonal into
!              .    - a KDU-by-KDU work array U in the lower
!              .      left-hand-corner,
!              .    - a KDU-by-at-least-KDU-but-more-is-better
!              .      (KDU-by-NHo) horizontal work array WH along
!              .      the bottom edge,
!              .    - and an at-least-KDU-but-more-is-better-by-KDU
!              .      (NVE-by-KDU) vertical work WV arrow along
!              .      the left-hand-edge. ====
!
               kdu = 2*ns
               ku = N - kdu + 1
               kwh = kdu + 1
               nho = (N-kdu+1-4) - (kdu+1) + 1
               kwv = kdu + 4
               nve = N - kdu - kwv + 1
!
!              ==== Small-bulge multi-shift QR sweep ====
!
               CALL CLAQR5(Wantt,Wantz,kacc22,N,ktop,kbot,ns,W(ks),H,   &
     &                     Ldh,Iloz,Ihiz,Z,Ldz,Work,3,H(ku,1),Ldh,nve,  &
     &                     H(kwv,1),Ldh,nho,H(ku,kwh),Ldh)
            ENDIF
!
!           ==== Note progress (or the lack of it). ====
!
            IF ( ld>0 ) THEN
               ndfl = 1
            ELSE
               ndfl = ndfl + 1
            ENDIF
!
!           ==== End of main loop ====
         ENDDO
!
!        ==== Iteration limit exceeded.  Set INFO to show where
!        .    the problem occurred and exit. ====
!
         Info = kbot
      ENDIF
!
!     ==== Return the optimal value of LWORK. ====
!
 100  Work(1) = CMPLX(lwkopt,0)
!
!     ==== End of CLAQR0 ====
!
      END SUBROUTINE CLAQR0
