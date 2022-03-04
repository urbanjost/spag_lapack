!*==dlarre.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b DLARRE given the tridiagonal matrix T, sets small off-diagonal elements to zero and for each unreduced block Ti, finds base representations and eigenvalues.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DLARRE + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlarre.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlarre.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlarre.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLARRE( RANGE, N, VL, VU, IL, IU, D, E, E2,
!                           RTOL1, RTOL2, SPLTOL, NSPLIT, ISPLIT, M,
!                           W, WERR, WGAP, IBLOCK, INDEXW, GERS, PIVMIN,
!                           WORK, IWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          RANGE
!       INTEGER            IL, INFO, IU, M, N, NSPLIT
!       DOUBLE PRECISION  PIVMIN, RTOL1, RTOL2, SPLTOL, VL, VU
!       ..
!       .. Array Arguments ..
!       INTEGER            IBLOCK( * ), ISPLIT( * ), IWORK( * ),
!      $                   INDEXW( * )
!       DOUBLE PRECISION   D( * ), E( * ), E2( * ), GERS( * ),
!      $                   W( * ),WERR( * ), WGAP( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> To find the desired eigenvalues of a given real symmetric
!> tridiagonal matrix T, DLARRE sets any "small" off-diagonal
!> elements to zero, and for each unreduced block T_i, it finds
!> (a) a suitable shift at one end of the block's spectrum,
!> (b) the base representation, T_i - sigma_i I = L_i D_i L_i^T, and
!> (c) eigenvalues of each L_i D_i L_i^T.
!> The representations and eigenvalues found are then used by
!> DSTEMR to compute the eigenvectors of T.
!> The accuracy varies depending on whether bisection is used to
!> find a few eigenvalues or the dqds algorithm (subroutine DLASQ2) to
!> conpute all and then discard any unwanted one.
!> As an added benefit, DLARRE also outputs the n
!> Gerschgorin intervals for the matrices L_i D_i L_i^T.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] RANGE
!> \verbatim
!>          RANGE is CHARACTER*1
!>          = 'A': ("All")   all eigenvalues will be found.
!>          = 'V': ("Value") all eigenvalues in the half-open interval
!>                           (VL, VU] will be found.
!>          = 'I': ("Index") the IL-th through IU-th eigenvalues (of the
!>                           entire matrix) will be found.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix. N > 0.
!> \endverbatim
!>
!> \param[in,out] VL
!> \verbatim
!>          VL is DOUBLE PRECISION
!>          If RANGE='V', the lower bound for the eigenvalues.
!>          Eigenvalues less than or equal to VL, or greater than VU,
!>          will not be returned.  VL < VU.
!>          If RANGE='I' or ='A', DLARRE computes bounds on the desired
!>          part of the spectrum.
!> \endverbatim
!>
!> \param[in,out] VU
!> \verbatim
!>          VU is DOUBLE PRECISION
!>          If RANGE='V', the upper bound for the eigenvalues.
!>          Eigenvalues less than or equal to VL, or greater than VU,
!>          will not be returned.  VL < VU.
!>          If RANGE='I' or ='A', DLARRE computes bounds on the desired
!>          part of the spectrum.
!> \endverbatim
!>
!> \param[in] IL
!> \verbatim
!>          IL is INTEGER
!>          If RANGE='I', the index of the
!>          smallest eigenvalue to be returned.
!>          1 <= IL <= IU <= N.
!> \endverbatim
!>
!> \param[in] IU
!> \verbatim
!>          IU is INTEGER
!>          If RANGE='I', the index of the
!>          largest eigenvalue to be returned.
!>          1 <= IL <= IU <= N.
!> \endverbatim
!>
!> \param[in,out] D
!> \verbatim
!>          D is DOUBLE PRECISION array, dimension (N)
!>          On entry, the N diagonal elements of the tridiagonal
!>          matrix T.
!>          On exit, the N diagonal elements of the diagonal
!>          matrices D_i.
!> \endverbatim
!>
!> \param[in,out] E
!> \verbatim
!>          E is DOUBLE PRECISION array, dimension (N)
!>          On entry, the first (N-1) entries contain the subdiagonal
!>          elements of the tridiagonal matrix T; E(N) need not be set.
!>          On exit, E contains the subdiagonal elements of the unit
!>          bidiagonal matrices L_i. The entries E( ISPLIT( I ) ),
!>          1 <= I <= NSPLIT, contain the base points sigma_i on output.
!> \endverbatim
!>
!> \param[in,out] E2
!> \verbatim
!>          E2 is DOUBLE PRECISION array, dimension (N)
!>          On entry, the first (N-1) entries contain the SQUARES of the
!>          subdiagonal elements of the tridiagonal matrix T;
!>          E2(N) need not be set.
!>          On exit, the entries E2( ISPLIT( I ) ),
!>          1 <= I <= NSPLIT, have been set to zero
!> \endverbatim
!>
!> \param[in] RTOL1
!> \verbatim
!>          RTOL1 is DOUBLE PRECISION
!> \endverbatim
!>
!> \param[in] RTOL2
!> \verbatim
!>          RTOL2 is DOUBLE PRECISION
!>           Parameters for bisection.
!>           An interval [LEFT,RIGHT] has converged if
!>           RIGHT-LEFT < MAX( RTOL1*GAP, RTOL2*MAX(|LEFT|,|RIGHT|) )
!> \endverbatim
!>
!> \param[in] SPLTOL
!> \verbatim
!>          SPLTOL is DOUBLE PRECISION
!>          The threshold for splitting.
!> \endverbatim
!>
!> \param[out] NSPLIT
!> \verbatim
!>          NSPLIT is INTEGER
!>          The number of blocks T splits into. 1 <= NSPLIT <= N.
!> \endverbatim
!>
!> \param[out] ISPLIT
!> \verbatim
!>          ISPLIT is INTEGER array, dimension (N)
!>          The splitting points, at which T breaks up into blocks.
!>          The first block consists of rows/columns 1 to ISPLIT(1),
!>          the second of rows/columns ISPLIT(1)+1 through ISPLIT(2),
!>          etc., and the NSPLIT-th consists of rows/columns
!>          ISPLIT(NSPLIT-1)+1 through ISPLIT(NSPLIT)=N.
!> \endverbatim
!>
!> \param[out] M
!> \verbatim
!>          M is INTEGER
!>          The total number of eigenvalues (of all L_i D_i L_i^T)
!>          found.
!> \endverbatim
!>
!> \param[out] W
!> \verbatim
!>          W is DOUBLE PRECISION array, dimension (N)
!>          The first M elements contain the eigenvalues. The
!>          eigenvalues of each of the blocks, L_i D_i L_i^T, are
!>          sorted in ascending order ( DLARRE may use the
!>          remaining N-M elements as workspace).
!> \endverbatim
!>
!> \param[out] WERR
!> \verbatim
!>          WERR is DOUBLE PRECISION array, dimension (N)
!>          The error bound on the corresponding eigenvalue in W.
!> \endverbatim
!>
!> \param[out] WGAP
!> \verbatim
!>          WGAP is DOUBLE PRECISION array, dimension (N)
!>          The separation from the right neighbor eigenvalue in W.
!>          The gap is only with respect to the eigenvalues of the same block
!>          as each block has its own representation tree.
!>          Exception: at the right end of a block we store the left gap
!> \endverbatim
!>
!> \param[out] IBLOCK
!> \verbatim
!>          IBLOCK is INTEGER array, dimension (N)
!>          The indices of the blocks (submatrices) associated with the
!>          corresponding eigenvalues in W; IBLOCK(i)=1 if eigenvalue
!>          W(i) belongs to the first block from the top, =2 if W(i)
!>          belongs to the second block, etc.
!> \endverbatim
!>
!> \param[out] INDEXW
!> \verbatim
!>          INDEXW is INTEGER array, dimension (N)
!>          The indices of the eigenvalues within each block (submatrix);
!>          for example, INDEXW(i)= 10 and IBLOCK(i)=2 imply that the
!>          i-th eigenvalue W(i) is the 10-th eigenvalue in block 2
!> \endverbatim
!>
!> \param[out] GERS
!> \verbatim
!>          GERS is DOUBLE PRECISION array, dimension (2*N)
!>          The N Gerschgorin intervals (the i-th Gerschgorin interval
!>          is (GERS(2*i-1), GERS(2*i)).
!> \endverbatim
!>
!> \param[out] PIVMIN
!> \verbatim
!>          PIVMIN is DOUBLE PRECISION
!>          The minimum pivot in the Sturm sequence for T.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (6*N)
!>          Workspace.
!> \endverbatim
!>
!> \param[out] IWORK
!> \verbatim
!>          IWORK is INTEGER array, dimension (5*N)
!>          Workspace.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          > 0:  A problem occurred in DLARRE.
!>          < 0:  One of the called subroutines signaled an internal problem.
!>                Needs inspection of the corresponding parameter IINFO
!>                for further information.
!>
!>          =-1:  Problem in DLARRD.
!>          = 2:  No base representation could be found in MAXTRY iterations.
!>                Increasing MAXTRY and recompilation might be a remedy.
!>          =-3:  Problem in DLARRB when computing the refined root
!>                representation for DLASQ2.
!>          =-4:  Problem in DLARRB when preforming bisection on the
!>                desired part of the spectrum.
!>          =-5:  Problem in DLASQ2.
!>          =-6:  Problem in DLASQ2.
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
!> \ingroup OTHERauxiliary
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  The base representations are required to suffer very little
!>  element growth and consequently define all their eigenvalues to
!>  high relative accuracy.
!> \endverbatim
!
!> \par Contributors:
!  ==================
!>
!>     Beresford Parlett, University of California, Berkeley, USA \n
!>     Jim Demmel, University of California, Berkeley, USA \n
!>     Inderjit Dhillon, University of Texas, Austin, USA \n
!>     Osni Marques, LBNL/NERSC, USA \n
!>     Christof Voemel, University of California, Berkeley, USA \n
!>
!  =====================================================================
      SUBROUTINE DLARRE(Range,N,Vl,Vu,Il,Iu,D,E,E2,Rtol1,Rtol2,Spltol,  &
     &                  Nsplit,Isplit,M,W,Werr,Wgap,Iblock,Indexw,Gers, &
     &                  Pivmin,Work,Iwork,Info)
      USE F77KINDS                        
      USE S_DCOPY
      USE S_DLAMCH
      USE S_DLARNV
      USE S_DLARRA
      USE S_DLARRB
      USE S_DLARRC
      USE S_DLARRD
      USE S_DLARRK
      USE S_DLASQ2
      USE S_LSAME
      IMPLICIT NONE
!*--DLARRE319
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D0 , ONE = 1.0D0 ,        &
     &                              TWO = 2.0D0 , FOUR = 4.0D0 ,        &
     &                              HNDRD = 100.0D0 , PERT = 8.0D0 ,    &
     &                              HALF = ONE/TWO , FOURTH = ONE/FOUR ,&
     &                              FAC = HALF , MAXGROWTH = 64.0D0 ,   &
     &                              FUDGE = 2.0D0
      INTEGER , PARAMETER  ::  MAXTRY = 6 , ALLRNG = 1 , INDRNG = 2 ,   &
     &                         VALRNG = 3
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Range
      INTEGER :: N
      REAL(R8KIND) , INTENT(INOUT) :: Vl
      REAL(R8KIND) , INTENT(INOUT) :: Vu
      INTEGER :: Il
      INTEGER :: Iu
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: D
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: E
      REAL(R8KIND) , DIMENSION(*) :: E2
      REAL(R8KIND) :: Rtol1
      REAL(R8KIND) :: Rtol2
      REAL(R8KIND) :: Spltol
      INTEGER :: Nsplit
      INTEGER , DIMENSION(*) :: Isplit
      INTEGER , INTENT(INOUT) :: M
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: W
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Werr
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Wgap
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iblock
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Indexw
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Gers
      REAL(R8KIND) , INTENT(INOUT) :: Pivmin
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(OUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      REAL(R8KIND) :: avgap , bsrtol , clwdth , dmax , dpivot , eabs ,  &
     &                emax , eold , eps , gl , gu , isleft , isrght ,   &
     &                rtl , rtol , s1 , s2 , safmin , sgndef , sigma ,  &
     &                spdiam , tau , tmp , tmp1
      INTEGER :: cnt , cnt1 , cnt2 , i , ibegin , idum , iend , iinfo , &
     &           in , indl , indu , irange , j , jblk , mb , mm ,       &
     &           wbegin , wend
      LOGICAL :: forceb , norep , usedqd
      INTEGER , DIMENSION(4) :: iseed
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
!     .. Local Arrays ..
!     ..
!     .. External Functions ..
 
!     ..
!     .. External Subroutines ..
!     ..
!     .. Intrinsic Functions ..
 
!     ..
!     .. Executable Statements ..
!
 
      Info = 0
!
!     Quick return if possible
!
      IF ( N<=0 ) RETURN
!
!     Decode RANGE
!
      IF ( LSAME(Range,'A') ) THEN
         irange = ALLRNG
      ELSEIF ( LSAME(Range,'V') ) THEN
         irange = VALRNG
      ELSEIF ( LSAME(Range,'I') ) THEN
         irange = INDRNG
      ENDIF
 
      M = 0
 
!     Get machine constants
      safmin = DLAMCH('S')
      eps = DLAMCH('P')
 
!     Set parameters
      rtl = SQRT(eps)
      bsrtol = SQRT(eps)
 
!     Treat case of 1x1 matrix for quick return
      IF ( N==1 ) THEN
         IF ( (irange==ALLRNG) .OR.                                     &
     &        ((irange==VALRNG) .AND. (D(1)>Vl) .AND. (D(1)<=Vu)) .OR.  &
     &        ((irange==INDRNG) .AND. (Il==1) .AND. (Iu==1)) ) THEN
            M = 1
            W(1) = D(1)
!           The computation error of the eigenvalue is zero
            Werr(1) = ZERO
            Wgap(1) = ZERO
            Iblock(1) = 1
            Indexw(1) = 1
            Gers(1) = D(1)
            Gers(2) = D(1)
         ENDIF
!        store the shift for the initial RRR, which is zero in this case
         E(1) = ZERO
         RETURN
      ENDIF
 
!     General case: tridiagonal matrix of order > 1
!
!     Init WERR, WGAP. Compute Gerschgorin intervals and spectral diameter.
!     Compute maximum off-diagonal entry and pivmin.
      gl = D(1)
      gu = D(1)
      eold = ZERO
      emax = ZERO
      E(N) = ZERO
      DO i = 1 , N
         Werr(i) = ZERO
         Wgap(i) = ZERO
         eabs = ABS(E(i))
         IF ( eabs>=emax ) emax = eabs
         tmp1 = eabs + eold
         Gers(2*i-1) = D(i) - tmp1
         gl = MIN(gl,Gers(2*i-1))
         Gers(2*i) = D(i) + tmp1
         gu = MAX(gu,Gers(2*i))
         eold = eabs
      ENDDO
!     The minimum pivot allowed in the Sturm sequence for T
      Pivmin = safmin*MAX(ONE,emax**2)
!     Compute spectral diameter. The Gerschgorin bounds give an
!     estimate that is wrong by at most a factor of SQRT(2)
      spdiam = gu - gl
 
!     Compute splitting points
      CALL DLARRA(N,D,E,E2,Spltol,spdiam,Nsplit,Isplit,iinfo)
 
!     Can force use of bisection instead of faster DQDS.
!     Option left in the code for future multisection work.
      forceb = .FALSE.
 
!     Initialize USEDQD, DQDS should be used for ALLRNG unless someone
!     explicitly wants bisection.
      usedqd = ((irange==ALLRNG) .AND. (.NOT.forceb))
 
      IF ( (irange==ALLRNG) .AND. (.NOT.forceb) ) THEN
!        Set interval [VL,VU] that contains all eigenvalues
         Vl = gl
         Vu = gu
      ELSE
!        We call DLARRD to find crude approximations to the eigenvalues
!        in the desired range. In case IRANGE = INDRNG, we also obtain the
!        interval (VL,VU] that contains all the wanted eigenvalues.
!        An interval [LEFT,RIGHT] has converged if
!        RIGHT-LEFT.LT.RTOL*MAX(ABS(LEFT),ABS(RIGHT))
!        DLARRD needs a WORK of size 4*N, IWORK of size 3*N
         CALL DLARRD(Range,'B',N,Vl,Vu,Il,Iu,Gers,bsrtol,D,E,E2,Pivmin, &
     &               Nsplit,Isplit,mm,W,Werr,Vl,Vu,Iblock,Indexw,Work,  &
     &               Iwork,iinfo)
         IF ( iinfo/=0 ) THEN
            Info = -1
            RETURN
         ENDIF
!        Make sure that the entries M+1 to N in W, WERR, IBLOCK, INDEXW are 0
         DO i = mm + 1 , N
            W(i) = ZERO
            Werr(i) = ZERO
            Iblock(i) = 0
            Indexw(i) = 0
         ENDDO
      ENDIF
 
 
!**
!     Loop over unreduced blocks
      ibegin = 1
      wbegin = 1
      DO jblk = 1 , Nsplit
         iend = Isplit(jblk)
         in = iend - ibegin + 1
 
!        1 X 1 block
         IF ( in==1 ) THEN
            IF ( (irange==ALLRNG) .OR.                                  &
     &           ((irange==VALRNG) .AND. (D(ibegin)>Vl) .AND.           &
     &           (D(ibegin)<=Vu)) .OR.                                  &
     &           ((irange==INDRNG) .AND. (Iblock(wbegin)==jblk)) ) THEN
               M = M + 1
               W(M) = D(ibegin)
               Werr(M) = ZERO
!              The gap for a single block doesn't matter for the later
!              algorithm and is assigned an arbitrary large value
               Wgap(M) = ZERO
               Iblock(M) = jblk
               Indexw(M) = 1
               wbegin = wbegin + 1
            ENDIF
!           E( IEND ) holds the shift for the initial RRR
            E(iend) = ZERO
            ibegin = iend + 1
            CYCLE
         ENDIF
!
!        Blocks of size larger than 1x1
!
!        E( IEND ) will hold the shift for the initial RRR, for now set it =0
         E(iend) = ZERO
!
!        Find local outer bounds GL,GU for the block
         gl = D(ibegin)
         gu = D(ibegin)
         DO i = ibegin , iend
            gl = MIN(Gers(2*i-1),gl)
            gu = MAX(Gers(2*i),gu)
         ENDDO
         spdiam = gu - gl
 
         IF ( .NOT.((irange==ALLRNG) .AND. (.NOT.forceb)) ) THEN
!           Count the number of eigenvalues in the current block.
            mb = 0
            DO i = wbegin , mm
               IF ( Iblock(i)/=jblk ) EXIT
               mb = mb + 1
            ENDDO
 
            IF ( mb==0 ) THEN
!              No eigenvalue in the current block lies in the desired range
!              E( IEND ) holds the shift for the initial RRR
               E(iend) = ZERO
               ibegin = iend + 1
               CYCLE
            ELSE
 
!              Decide whether dqds or bisection is more efficient
               usedqd = ((mb>FAC*in) .AND. (.NOT.forceb))
               wend = wbegin + mb - 1
!              Calculate gaps for the current block
!              In later stages, when representations for individual
!              eigenvalues are different, we use SIGMA = E( IEND ).
               sigma = ZERO
               DO i = wbegin , wend - 1
                  Wgap(i) = MAX(ZERO,W(i+1)-Werr(i+1)-(W(i)+Werr(i)))
               ENDDO
               Wgap(wend) = MAX(ZERO,Vu-sigma-(W(wend)+Werr(wend)))
!              Find local index of the first and last desired evalue.
               indl = Indexw(wbegin)
               indu = Indexw(wend)
            ENDIF
         ENDIF
         IF ( ((irange==ALLRNG) .AND. (.NOT.forceb)) .OR. usedqd ) THEN
!           Case of DQDS
!           Find approximations to the extremal eigenvalues of the block
            CALL DLARRK(in,1,gl,gu,D(ibegin),E2(ibegin),Pivmin,rtl,tmp, &
     &                  tmp1,iinfo)
            IF ( iinfo/=0 ) THEN
               Info = -1
               RETURN
            ENDIF
            isleft = MAX(gl,tmp-tmp1-HNDRD*eps*ABS(tmp-tmp1))
 
            CALL DLARRK(in,in,gl,gu,D(ibegin),E2(ibegin),Pivmin,rtl,tmp,&
     &                  tmp1,iinfo)
            IF ( iinfo/=0 ) THEN
               Info = -1
               RETURN
            ENDIF
            isrght = MIN(gu,tmp+tmp1+HNDRD*eps*ABS(tmp+tmp1))
!           Improve the estimate of the spectral diameter
            spdiam = isrght - isleft
         ELSE
!           Case of bisection
!           Find approximations to the wanted extremal eigenvalues
            isleft = MAX(gl,W(wbegin)-Werr(wbegin)                      &
     &               -HNDRD*eps*ABS(W(wbegin)-Werr(wbegin)))
            isrght = MIN(gu,W(wend)+Werr(wend)                          &
     &               +HNDRD*eps*ABS(W(wend)+Werr(wend)))
         ENDIF
 
 
!        Decide whether the base representation for the current block
!        L_JBLK D_JBLK L_JBLK^T = T_JBLK - sigma_JBLK I
!        should be on the left or the right end of the current block.
!        The strategy is to shift to the end which is "more populated"
!        Furthermore, decide whether to use DQDS for the computation of
!        the eigenvalue approximations at the end of DLARRE or bisection.
!        dqds is chosen if all eigenvalues are desired or the number of
!        eigenvalues to be computed is large compared to the blocksize.
         IF ( (irange==ALLRNG) .AND. (.NOT.forceb) ) THEN
!           If all the eigenvalues have to be computed, we use dqd
            usedqd = .TRUE.
!           INDL is the local index of the first eigenvalue to compute
            indl = 1
            indu = in
!           MB =  number of eigenvalues to compute
            mb = in
            wend = wbegin + mb - 1
!           Define 1/4 and 3/4 points of the spectrum
            s1 = isleft + FOURTH*spdiam
            s2 = isrght - FOURTH*spdiam
!           DLARRD has computed IBLOCK and INDEXW for each eigenvalue
!           approximation.
!           choose sigma
         ELSEIF ( usedqd ) THEN
            s1 = isleft + FOURTH*spdiam
            s2 = isrght - FOURTH*spdiam
         ELSE
            tmp = MIN(isrght,Vu) - MAX(isleft,Vl)
            s1 = MAX(isleft,Vl) + FOURTH*tmp
            s2 = MIN(isrght,Vu) - FOURTH*tmp
         ENDIF
 
!        Compute the negcount at the 1/4 and 3/4 points
         IF ( mb>1 ) CALL DLARRC('T',in,s1,s2,D(ibegin),E(ibegin),      &
     &                           Pivmin,cnt,cnt1,cnt2,iinfo)
 
         IF ( mb==1 ) THEN
            sigma = gl
            sgndef = ONE
         ELSEIF ( cnt1-indl>=indu-cnt2 ) THEN
            IF ( (irange==ALLRNG) .AND. (.NOT.forceb) ) THEN
               sigma = MAX(isleft,gl)
            ELSEIF ( usedqd ) THEN
!              use Gerschgorin bound as shift to get pos def matrix
!              for dqds
               sigma = isleft
            ELSE
!              use approximation of the first desired eigenvalue of the
!              block as shift
               sigma = MAX(isleft,Vl)
            ENDIF
            sgndef = ONE
         ELSE
            IF ( (irange==ALLRNG) .AND. (.NOT.forceb) ) THEN
               sigma = MIN(isrght,gu)
            ELSEIF ( usedqd ) THEN
!              use Gerschgorin bound as shift to get neg def matrix
!              for dqds
               sigma = isrght
            ELSE
!              use approximation of the first desired eigenvalue of the
!              block as shift
               sigma = MIN(isrght,Vu)
            ENDIF
            sgndef = -ONE
         ENDIF
 
 
!        An initial SIGMA has been chosen that will be used for computing
!        T - SIGMA I = L D L^T
!        Define the increment TAU of the shift in case the initial shift
!        needs to be refined to obtain a factorization with not too much
!        element growth.
         IF ( usedqd ) THEN
!           The initial SIGMA was to the outer end of the spectrum
!           the matrix is definite and we need not retreat.
            tau = spdiam*eps*N + TWO*Pivmin
            tau = MAX(tau,TWO*eps*ABS(sigma))
         ELSEIF ( mb>1 ) THEN
            clwdth = W(wend) + Werr(wend) - W(wbegin) - Werr(wbegin)
            avgap = ABS(clwdth/DBLE(wend-wbegin))
            IF ( sgndef==ONE ) THEN
               tau = HALF*MAX(Wgap(wbegin),avgap)
               tau = MAX(tau,Werr(wbegin))
            ELSE
               tau = HALF*MAX(Wgap(wend-1),avgap)
               tau = MAX(tau,Werr(wend))
            ENDIF
         ELSE
            tau = Werr(wbegin)
         ENDIF
!
         DO idum = 1 , MAXTRY
!           Compute L D L^T factorization of tridiagonal matrix T - sigma I.
!           Store D in WORK(1:IN), L in WORK(IN+1:2*IN), and reciprocals of
!           pivots in WORK(2*IN+1:3*IN)
            dpivot = D(ibegin) - sigma
            Work(1) = dpivot
            dmax = ABS(Work(1))
            j = ibegin
            DO i = 1 , in - 1
               Work(2*in+i) = ONE/Work(i)
               tmp = E(j)*Work(2*in+i)
               Work(in+i) = tmp
               dpivot = (D(j+1)-sigma) - tmp*E(j)
               Work(i+1) = dpivot
               dmax = MAX(dmax,ABS(dpivot))
               j = j + 1
            ENDDO
!           check for element growth
            IF ( dmax>MAXGROWTH*spdiam ) THEN
               norep = .TRUE.
            ELSE
               norep = .FALSE.
            ENDIF
            IF ( usedqd .AND. .NOT.norep ) THEN
!              Ensure the definiteness of the representation
!              All entries of D (of L D L^T) must have the same sign
               DO i = 1 , in
                  tmp = sgndef*Work(i)
                  IF ( tmp<ZERO ) norep = .TRUE.
               ENDDO
            ENDIF
!              an initial RRR is found
            IF ( .NOT.(norep) ) GOTO 50
!              Note that in the case of IRANGE=ALLRNG, we use the Gerschgorin
!              shift which makes the matrix definite. So we should end up
!              here really only in the case of IRANGE = VALRNG or INDRNG.
            IF ( idum/=MAXTRY-1 ) THEN
               sigma = sigma - sgndef*tau
               tau = TWO*tau
            ELSEIF ( sgndef==ONE ) THEN
!                    The fudged Gerschgorin shift should succeed
               sigma = gl - FUDGE*spdiam*eps*N - FUDGE*TWO*Pivmin
            ELSE
               sigma = gu + FUDGE*spdiam*eps*N + FUDGE*TWO*Pivmin
            ENDIF
         ENDDO
!        if the program reaches this point, no base representation could be
!        found in MAXTRY iterations.
         Info = 2
         RETURN
 
!        At this point, we have found an initial base representation
!        T - SIGMA I = L D L^T with not too much element growth.
!        Store the shift.
 50      E(iend) = sigma
!        Store D and L.
         CALL DCOPY(in,Work,1,D(ibegin),1)
         CALL DCOPY(in-1,Work(in+1),1,E(ibegin),1)
 
 
         IF ( mb>1 ) THEN
!
!           Perturb each entry of the base representation by a small
!           (but random) relative amount to overcome difficulties with
!           glued matrices.
!
            DO i = 1 , 4
               iseed(i) = 1
            ENDDO
 
            CALL DLARNV(2,iseed,2*in-1,Work(1))
            DO i = 1 , in - 1
               D(ibegin+i-1) = D(ibegin+i-1)*(ONE+eps*PERT*Work(i))
               E(ibegin+i-1) = E(ibegin+i-1)*(ONE+eps*PERT*Work(in+i))
            ENDDO
            D(iend) = D(iend)*(ONE+eps*FOUR*Work(in))
!
         ENDIF
!
!        Don't update the Gerschgorin intervals because keeping track
!        of the updates would be too much work in DLARRV.
!        We update W instead and use it to locate the proper Gerschgorin
!        intervals.
 
!        Compute the required eigenvalues of L D L' by bisection or dqds
         IF ( .NOT.usedqd ) THEN
!           If DLARRD has been used, shift the eigenvalue approximations
!           according to their representation. This is necessary for
!           a uniform DLARRV since dqds computes eigenvalues of the
!           shifted representation. In DLARRV, W will always hold the
!           UNshifted eigenvalue approximation.
            DO j = wbegin , wend
               W(j) = W(j) - sigma
               Werr(j) = Werr(j) + ABS(W(j))*eps
            ENDDO
!           call DLARRB to reduce eigenvalue error of the approximations
!           from DLARRD
            DO i = ibegin , iend - 1
               Work(i) = D(i)*E(i)**2
            ENDDO
!           use bisection to find EV from INDL to INDU
            CALL DLARRB(in,D(ibegin),Work(ibegin),indl,indu,Rtol1,Rtol2,&
     &                  indl-1,W(wbegin),Wgap(wbegin),Werr(wbegin),     &
     &                  Work(2*N+1),Iwork,Pivmin,spdiam,in,iinfo)
            IF ( iinfo/=0 ) THEN
               Info = -4
               RETURN
            ENDIF
!           DLARRB computes all gaps correctly except for the last one
!           Record distance to VU/GU
            Wgap(wend) = MAX(ZERO,(Vu-sigma)-(W(wend)+Werr(wend)))
            DO i = indl , indu
               M = M + 1
               Iblock(M) = jblk
               Indexw(M) = i
            ENDDO
         ELSE
!           Call dqds to get all eigs (and then possibly delete unwanted
!           eigenvalues).
!           Note that dqds finds the eigenvalues of the L D L^T representation
!           of T to high relative accuracy. High relative accuracy
!           might be lost when the shift of the RRR is subtracted to obtain
!           the eigenvalues of T. However, T is not guaranteed to define its
!           eigenvalues to high relative accuracy anyway.
!           Set RTOL to the order of the tolerance used in DLASQ2
!           This is an ESTIMATED error, the worst case bound is 4*N*EPS
!           which is usually too large and requires unnecessary work to be
!           done by bisection when computing the eigenvectors
            rtol = LOG(DBLE(in))*FOUR*eps
            j = ibegin
            DO i = 1 , in - 1
               Work(2*i-1) = ABS(D(j))
               Work(2*i) = E(j)*E(j)*Work(2*i-1)
               j = j + 1
            ENDDO
            Work(2*in-1) = ABS(D(iend))
            Work(2*in) = ZERO
            CALL DLASQ2(in,Work,iinfo)
            IF ( iinfo/=0 ) THEN
!              If IINFO = -5 then an index is part of a tight cluster
!              and should be changed. The index is in IWORK(1) and the
!              gap is in WORK(N+1)
               Info = -5
               RETURN
            ELSE
!              Test that all eigenvalues are positive as expected
               DO i = 1 , in
                  IF ( Work(i)<ZERO ) THEN
                     Info = -6
                     RETURN
                  ENDIF
               ENDDO
            ENDIF
            IF ( sgndef>ZERO ) THEN
               DO i = indl , indu
                  M = M + 1
                  W(M) = Work(in-i+1)
                  Iblock(M) = jblk
                  Indexw(M) = i
               ENDDO
            ELSE
               DO i = indl , indu
                  M = M + 1
                  W(M) = -Work(i)
                  Iblock(M) = jblk
                  Indexw(M) = i
               ENDDO
            ENDIF
 
            DO i = M - mb + 1 , M
!              the value of RTOL below should be the tolerance in DLASQ2
               Werr(i) = rtol*ABS(W(i))
            ENDDO
            DO i = M - mb + 1 , M - 1
!              compute the right gap between the intervals
               Wgap(i) = MAX(ZERO,W(i+1)-Werr(i+1)-(W(i)+Werr(i)))
            ENDDO
            Wgap(M) = MAX(ZERO,(Vu-sigma)-(W(M)+Werr(M)))
         ENDIF
!        proceed with next block
         ibegin = iend + 1
         wbegin = wend + 1
      ENDDO
!
 
!
!     end of DLARRE
!
      END SUBROUTINE DLARRE
