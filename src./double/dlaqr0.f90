!*==dlaqr0.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b DLAQR0 computes the eigenvalues of a Hessenberg matrix, and optionally the matrices from the Schur decomposition.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DLAQR0 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaqr0.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaqr0.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaqr0.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLAQR0( WANTT, WANTZ, N, ILO, IHI, H, LDH, WR, WI,
!                          ILOZ, IHIZ, Z, LDZ, WORK, LWORK, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            IHI, IHIZ, ILO, ILOZ, INFO, LDH, LDZ, LWORK, N
!       LOGICAL            WANTT, WANTZ
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   H( LDH, * ), WI( * ), WORK( * ), WR( * ),
!      $                   Z( LDZ, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    DLAQR0 computes the eigenvalues of a Hessenberg matrix H
!>    and, optionally, the matrices T and Z from the Schur decomposition
!>    H = Z T Z**T, where T is an upper quasi-triangular matrix (the
!>    Schur form), and Z is the orthogonal matrix of Schur vectors.
!>
!>    Optionally Z may be postmultiplied into an input orthogonal
!>    matrix Q so that this routine can give the Schur factorization
!>    of a matrix A which has been reduced to the Hessenberg form H
!>    by the orthogonal matrix Q:  A = Q*H*Q**T = (QZ)*T*(QZ)**T.
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
!>           previous call to DGEBAL, and then passed to DGEHRD when the
!>           matrix output by DGEBAL is reduced to Hessenberg form.
!>           Otherwise, ILO and IHI should be set to 1 and N,
!>           respectively.  If N > 0, then 1 <= ILO <= IHI <= N.
!>           If N = 0, then ILO = 1 and IHI = 0.
!> \endverbatim
!>
!> \param[in,out] H
!> \verbatim
!>          H is DOUBLE PRECISION array, dimension (LDH,N)
!>           On entry, the upper Hessenberg matrix H.
!>           On exit, if INFO = 0 and WANTT is .TRUE., then H contains
!>           the upper quasi-triangular matrix T from the Schur
!>           decomposition (the Schur form); 2-by-2 diagonal blocks
!>           (corresponding to complex conjugate pairs of eigenvalues)
!>           are returned in standard form, with H(i,i) = H(i+1,i+1)
!>           and H(i+1,i)*H(i,i+1) < 0. If INFO = 0 and WANTT is
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
!> \param[out] WR
!> \verbatim
!>          WR is DOUBLE PRECISION array, dimension (IHI)
!> \endverbatim
!>
!> \param[out] WI
!> \verbatim
!>          WI is DOUBLE PRECISION array, dimension (IHI)
!>           The real and imaginary parts, respectively, of the computed
!>           eigenvalues of H(ILO:IHI,ILO:IHI) are stored in WR(ILO:IHI)
!>           and WI(ILO:IHI). If two eigenvalues are computed as a
!>           complex conjugate pair, they are stored in consecutive
!>           elements of WR and WI, say the i-th and (i+1)th, with
!>           WI(i) > 0 and WI(i+1) < 0. If WANTT is .TRUE., then
!>           the eigenvalues are stored in the same order as on the
!>           diagonal of the Schur form returned in H, with
!>           WR(i) = H(i,i) and, if H(i:i+1,i:i+1) is a 2-by-2 diagonal
!>           block, WI(i) = sqrt(-H(i+1,i)*H(i,i+1)) and
!>           WI(i+1) = -WI(i).
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
!>          Z is DOUBLE PRECISION array, dimension (LDZ,IHI)
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
!>          WORK is DOUBLE PRECISION array, dimension LWORK
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
!>           If LWORK = -1, then DLAQR0 does a workspace query.
!>           In this case, DLAQR0 checks the input parameters and
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
!>             > 0:  if INFO = i, DLAQR0 failed to compute all of
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
!>                where U is an orthogonal matrix.  The final
!>                value of H is upper Hessenberg and quasi-triangular
!>                in rows and columns INFO+1 through IHI.
!>
!>                If INFO > 0 and WANTZ is .TRUE., then on exit
!>
!>                  (final value of Z(ILO:IHI,ILOZ:IHIZ)
!>                   =  (initial value of Z(ILO:IHI,ILOZ:IHIZ)*U
!>
!>                where U is the orthogonal matrix in (*) (regard-
!>                less of the value of WANTT.)
!>
!>                If INFO > 0 and WANTZ is .FALSE., then Z is not
!>                accessed.
!> \endverbatim
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
!> \ingroup doubleOTHERauxiliary
!
!  =====================================================================
      SUBROUTINE DLAQR0(Wantt,Wantz,N,Ilo,Ihi,H,Ldh,Wr,Wi,Iloz,Ihiz,Z,  &
     &                  Ldz,Work,Lwork,Info)
      USE F77KINDS                        
      USE S_DLACPY
      USE S_DLAHQR
      USE S_DLANV2
      USE S_DLAQR3
      USE S_DLAQR4
      USE S_DLAQR5
      USE S_ILAENV
      IMPLICIT NONE
!*--DLAQR0268
!
! PARAMETER definitions rewritten by SPAG
!
      INTEGER , PARAMETER  ::  NTINY = 15 , KEXNW = 5 , KEXSH = 6
      REAL(R8KIND) , PARAMETER  ::  WILK1 = 0.75D0 , WILK2 = -0.4375D0 ,&
     &                              ZERO = 0.0D0 , ONE = 1.0D0
!
! Dummy argument declarations rewritten by SPAG
!
      LOGICAL :: Wantt
      LOGICAL :: Wantz
      INTEGER :: N
      INTEGER :: Ilo
      INTEGER :: Ihi
      REAL(R8KIND) , DIMENSION(Ldh,*) :: H
      INTEGER :: Ldh
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Wr
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Wi
      INTEGER :: Iloz
      INTEGER :: Ihiz
      REAL(R8KIND) , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER :: Info
!
! Local variable declarations rewritten by SPAG
!
      REAL(R8KIND) :: aa , bb , cc , cs , dd , sn , ss , swap
      INTEGER :: i , inf , it , itmax , k , kacc22 , kbot , kdu , ks ,  &
     &           kt , ktop , ku , kv , kwh , kwtop , kwv , ld , ls ,    &
     &           lwkopt , ndec , ndfl , nh , nho , nibble , nmin , ns , &
     &           nsmax , nsr , nve , nw , nwmax , nwr , nwupbd
      CHARACTER(2) :: jbcmpz
      LOGICAL :: sorted
      REAL(R8KIND) , DIMENSION(1,1) :: zdum
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
!
!     ==== Matrices of order NTINY or smaller must be processed by
!     .    DLAHQR because of insufficient subdiagonal scratch space.
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
!     ==== The constants WILK1 and WILK2 are used to form the
!     .    exceptional shifts. ====
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
!        ==== Tiny matrices must use DLAHQR. ====
!
         lwkopt = 1
         IF ( Lwork/=-1 ) CALL DLAHQR(Wantt,Wantz,N,Ilo,Ihi,H,Ldh,Wr,Wi,&
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
!        .    point,  N .GT. NTINY = 14, so there is enough
!        .    subdiagonal workspace for NWR.GE.2 as required.
!        .    (In fact, there is enough subdiagonal space for
!        .    NWR.GE.4.) ====
!
         nwr = ILAENV(13,'DLAQR0',jbcmpz,N,Ilo,Ihi,Lwork)
         nwr = MAX(2,nwr)
         nwr = MIN(Ihi-Ilo+1,(N-1)/3,nwr)
!
!        ==== NSR = recommended number of simultaneous shifts.
!        .    At this point N .GT. NTINY = 15, so there is at
!        .    enough subdiagonal workspace for NSR to be even
!        .    and greater than or equal to two as required. ====
!
         nsr = ILAENV(15,'DLAQR0',jbcmpz,N,Ilo,Ihi,Lwork)
         nsr = MIN(nsr,(N-3)/6,Ihi-Ilo)
         nsr = MAX(2,nsr-MOD(nsr,2))
!
!        ==== Estimate optimal workspace ====
!
!        ==== Workspace query call to DLAQR3 ====
!
         CALL DLAQR3(Wantt,Wantz,N,Ilo,Ihi,nwr+1,H,Ldh,Iloz,Ihiz,Z,Ldz, &
     &               ls,ld,Wr,Wi,H,Ldh,N,H,Ldh,N,H,Ldh,Work,-1)
!
!        ==== Optimal workspace = MAX(DLAQR5, DLAQR3) ====
!
         lwkopt = MAX(3*nsr/2,INT(Work(1)))
!
!        ==== Quick return in case of workspace query. ====
!
         IF ( Lwork==-1 ) THEN
            Work(1) = DBLE(lwkopt)
            RETURN
         ENDIF
!
!        ==== DLAHQR/DLAQR0 crossover point ====
!
         nmin = ILAENV(12,'DLAQR0',jbcmpz,N,Ilo,Ihi,Lwork)
         nmin = MAX(NTINY,nmin)
!
!        ==== Nibble crossover point ====
!
         nibble = ILAENV(14,'DLAQR0',jbcmpz,N,Ilo,Ihi,Lwork)
         nibble = MAX(0,nibble)
!
!        ==== Accumulate reflections during ttswp?  Use block
!        .    2-by-2 structure during matrix-matrix multiply? ====
!
         kacc22 = ILAENV(16,'DLAQR0',jbcmpz,N,Ilo,Ihi,Lwork)
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
                  IF ( ABS(H(kwtop,kwtop-1))>ABS(H(kwtop-1,kwtop-2)) )  &
     &                 nw = nw + 1
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
            CALL DLAQR3(Wantt,Wantz,N,ktop,kbot,nw,H,Ldh,Iloz,Ihiz,Z,   &
     &                  Ldz,ls,ld,Wr,Wi,H(kv,1),Ldh,nho,H(kv,kt),Ldh,   &
     &                  nve,H(kwv,1),Ldh,Work,Lwork)
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
!              .    This may be lowered (slightly) if DLAQR3
!              .    did not provide that many shifts. ====
!
               ns = MIN(nsmax,nsr,MAX(2,kbot-ktop))
               ns = ns - MOD(ns,2)
!
!              ==== If there have been no deflations
!              .    in a multiple of KEXSH iterations,
!              .    then try exceptional shifts.
!              .    Otherwise use shifts provided by
!              .    DLAQR3 above or from the eigenvalues
!              .    of a trailing principal submatrix. ====
!
               IF ( MOD(ndfl,KEXSH)==0 ) THEN
                  ks = kbot - ns + 1
                  DO i = kbot , MAX(ks+1,ktop+2) , -2
                     ss = ABS(H(i,i-1)) + ABS(H(i-1,i-2))
                     aa = WILK1*ss + H(i,i)
                     bb = ss
                     cc = WILK2*ss
                     dd = aa
                     CALL DLANV2(aa,bb,cc,dd,Wr(i-1),Wi(i-1),Wr(i),     &
     &                           Wi(i),cs,sn)
                  ENDDO
                  IF ( ks==ktop ) THEN
                     Wr(ks+1) = H(ks+1,ks+1)
                     Wi(ks+1) = ZERO
                     Wr(ks) = Wr(ks+1)
                     Wi(ks) = Wi(ks+1)
                  ENDIF
               ELSE
!
!                 ==== Got NS/2 or fewer shifts? Use DLAQR4 or
!                 .    DLAHQR on a trailing principal submatrix to
!                 .    get more. (Since NS.LE.NSMAX.LE.(N+6)/9,
!                 .    there is enough space below the subdiagonal
!                 .    to fit an NS-by-NS scratch array.) ====
!
                  IF ( kbot-ks+1<=ns/2 ) THEN
                     ks = kbot - ns + 1
                     kt = N - ns + 1
                     CALL DLACPY('A',ns,ns,H(ks,ks),Ldh,H(kt,1),Ldh)
                     IF ( ns>nmin ) THEN
                        CALL DLAQR4(.FALSE.,.FALSE.,ns,1,ns,H(kt,1),Ldh,&
     &                              Wr(ks),Wi(ks),1,1,zdum,1,Work,Lwork,&
     &                              inf)
                     ELSE
                        CALL DLAHQR(.FALSE.,.FALSE.,ns,1,ns,H(kt,1),Ldh,&
     &                              Wr(ks),Wi(ks),1,1,zdum,1,inf)
                     ENDIF
                     ks = ks + inf
!
!                    ==== In case of a rare QR failure use
!                    .    eigenvalues of the trailing 2-by-2
!                    .    principal submatrix.  ====
!
                     IF ( ks>=kbot ) THEN
                        aa = H(kbot-1,kbot-1)
                        cc = H(kbot,kbot-1)
                        bb = H(kbot-1,kbot)
                        dd = H(kbot,kbot)
                        CALL DLANV2(aa,bb,cc,dd,Wr(kbot-1),Wi(kbot-1),  &
     &                              Wr(kbot),Wi(kbot),cs,sn)
                        ks = kbot - 1
                     ENDIF
                  ENDIF
!
                  IF ( kbot-ks+1>ns ) THEN
!
!                    ==== Sort the shifts (Helps a little)
!                    .    Bubble sort keeps complex conjugate
!                    .    pairs together. ====
!
                     sorted = .FALSE.
                     DO k = kbot , ks + 1 , -1
                        IF ( sorted ) EXIT
                        sorted = .TRUE.
                        DO i = ks , k - 1
                           IF ( ABS(Wr(i))+ABS(Wi(i))<ABS(Wr(i+1))      &
     &                          +ABS(Wi(i+1)) ) THEN
                              sorted = .FALSE.
!
                              swap = Wr(i)
                              Wr(i) = Wr(i+1)
                              Wr(i+1) = swap
!
                              swap = Wi(i)
                              Wi(i) = Wi(i+1)
                              Wi(i+1) = swap
                           ENDIF
                        ENDDO
                     ENDDO
                  ENDIF
!
!                 ==== Shuffle shifts into pairs of real shifts
!                 .    and pairs of complex conjugate shifts
!                 .    assuming complex conjugate shifts are
!                 .    already adjacent to one another. (Yes,
!                 .    they are.)  ====
!
                  DO i = kbot , ks + 2 , -2
                     IF ( Wi(i)/=-Wi(i-1) ) THEN
!
                        swap = Wr(i)
                        Wr(i) = Wr(i-1)
                        Wr(i-1) = Wr(i-2)
                        Wr(i-2) = swap
!
                        swap = Wi(i)
                        Wi(i) = Wi(i-1)
                        Wi(i-1) = Wi(i-2)
                        Wi(i-2) = swap
                     ENDIF
                  ENDDO
               ENDIF
!
!              ==== If there are only two shifts and both are
!              .    real, then use only one.  ====
!
               IF ( kbot-ks+1==2 ) THEN
                  IF ( Wi(kbot)==ZERO ) THEN
                     IF ( ABS(Wr(kbot)-H(kbot,kbot))                    &
     &                    <ABS(Wr(kbot-1)-H(kbot,kbot)) ) THEN
                        Wr(kbot-1) = Wr(kbot)
                     ELSE
                        Wr(kbot) = Wr(kbot-1)
                     ENDIF
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
               CALL DLAQR5(Wantt,Wantz,kacc22,N,ktop,kbot,ns,Wr(ks),    &
     &                     Wi(ks),H,Ldh,Iloz,Ihiz,Z,Ldz,Work,3,H(ku,1), &
     &                     Ldh,nve,H(kwv,1),Ldh,nho,H(ku,kwh),Ldh)
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
 100  Work(1) = DBLE(lwkopt)
!
!     ==== End of DLAQR0 ====
!
      END SUBROUTINE DLAQR0
