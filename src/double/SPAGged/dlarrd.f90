!*==dlarrd.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b DLARRD computes the eigenvalues of a symmetric tridiagonal matrix to suitable accuracy.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DLARRD + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlarrd.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlarrd.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlarrd.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLARRD( RANGE, ORDER, N, VL, VU, IL, IU, GERS,
!                           RELTOL, D, E, E2, PIVMIN, NSPLIT, ISPLIT,
!                           M, W, WERR, WL, WU, IBLOCK, INDEXW,
!                           WORK, IWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          ORDER, RANGE
!       INTEGER            IL, INFO, IU, M, N, NSPLIT
!       DOUBLE PRECISION    PIVMIN, RELTOL, VL, VU, WL, WU
!       ..
!       .. Array Arguments ..
!       INTEGER            IBLOCK( * ), INDEXW( * ),
!      $                   ISPLIT( * ), IWORK( * )
!       DOUBLE PRECISION   D( * ), E( * ), E2( * ),
!      $                   GERS( * ), W( * ), WERR( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLARRD computes the eigenvalues of a symmetric tridiagonal
!> matrix T to suitable accuracy. This is an auxiliary code to be
!> called from DSTEMR.
!> The user may ask for all eigenvalues, all eigenvalues
!> in the half-open interval (VL, VU], or the IL-th through IU-th
!> eigenvalues.
!>
!> To avoid overflow, the matrix must be scaled so that its
!> largest element is no greater than overflow**(1/2) * underflow**(1/4) in absolute value, and for greatest
!> accuracy, it should not be much smaller than that.
!>
!> See W. Kahan "Accurate Eigenvalues of a Symmetric Tridiagonal
!> Matrix", Report CS41, Computer Science Dept., Stanford
!> University, July 21, 1966.
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
!> \param[in] ORDER
!> \verbatim
!>          ORDER is CHARACTER*1
!>          = 'B': ("By Block") the eigenvalues will be grouped by
!>                              split-off block (see IBLOCK, ISPLIT) and
!>                              ordered from smallest to largest within
!>                              the block.
!>          = 'E': ("Entire matrix")
!>                              the eigenvalues for the entire matrix
!>                              will be ordered from smallest to
!>                              largest.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the tridiagonal matrix T.  N >= 0.
!> \endverbatim
!>
!> \param[in] VL
!> \verbatim
!>          VL is DOUBLE PRECISION
!>          If RANGE='V', the lower bound of the interval to
!>          be searched for eigenvalues.  Eigenvalues less than or equal
!>          to VL, or greater than VU, will not be returned.  VL < VU.
!>          Not referenced if RANGE = 'A' or 'I'.
!> \endverbatim
!>
!> \param[in] VU
!> \verbatim
!>          VU is DOUBLE PRECISION
!>          If RANGE='V', the upper bound of the interval to
!>          be searched for eigenvalues.  Eigenvalues less than or equal
!>          to VL, or greater than VU, will not be returned.  VL < VU.
!>          Not referenced if RANGE = 'A' or 'I'.
!> \endverbatim
!>
!> \param[in] IL
!> \verbatim
!>          IL is INTEGER
!>          If RANGE='I', the index of the
!>          smallest eigenvalue to be returned.
!>          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0.
!>          Not referenced if RANGE = 'A' or 'V'.
!> \endverbatim
!>
!> \param[in] IU
!> \verbatim
!>          IU is INTEGER
!>          If RANGE='I', the index of the
!>          largest eigenvalue to be returned.
!>          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0.
!>          Not referenced if RANGE = 'A' or 'V'.
!> \endverbatim
!>
!> \param[in] GERS
!> \verbatim
!>          GERS is DOUBLE PRECISION array, dimension (2*N)
!>          The N Gerschgorin intervals (the i-th Gerschgorin interval
!>          is (GERS(2*i-1), GERS(2*i)).
!> \endverbatim
!>
!> \param[in] RELTOL
!> \verbatim
!>          RELTOL is DOUBLE PRECISION
!>          The minimum relative width of an interval.  When an interval
!>          is narrower than RELTOL times the larger (in
!>          magnitude) endpoint, then it is considered to be
!>          sufficiently small, i.e., converged.  Note: this should
!>          always be at least radix*machine epsilon.
!> \endverbatim
!>
!> \param[in] D
!> \verbatim
!>          D is DOUBLE PRECISION array, dimension (N)
!>          The n diagonal elements of the tridiagonal matrix T.
!> \endverbatim
!>
!> \param[in] E
!> \verbatim
!>          E is DOUBLE PRECISION array, dimension (N-1)
!>          The (n-1) off-diagonal elements of the tridiagonal matrix T.
!> \endverbatim
!>
!> \param[in] E2
!> \verbatim
!>          E2 is DOUBLE PRECISION array, dimension (N-1)
!>          The (n-1) squared off-diagonal elements of the tridiagonal matrix T.
!> \endverbatim
!>
!> \param[in] PIVMIN
!> \verbatim
!>          PIVMIN is DOUBLE PRECISION
!>          The minimum pivot allowed in the Sturm sequence for T.
!> \endverbatim
!>
!> \param[in] NSPLIT
!> \verbatim
!>          NSPLIT is INTEGER
!>          The number of diagonal blocks in the matrix T.
!>          1 <= NSPLIT <= N.
!> \endverbatim
!>
!> \param[in] ISPLIT
!> \verbatim
!>          ISPLIT is INTEGER array, dimension (N)
!>          The splitting points, at which T breaks up into submatrices.
!>          The first submatrix consists of rows/columns 1 to ISPLIT(1),
!>          the second of rows/columns ISPLIT(1)+1 through ISPLIT(2),
!>          etc., and the NSPLIT-th consists of rows/columns
!>          ISPLIT(NSPLIT-1)+1 through ISPLIT(NSPLIT)=N.
!>          (Only the first NSPLIT elements will actually be used, but
!>          since the user cannot know a priori what value NSPLIT will
!>          have, N words must be reserved for ISPLIT.)
!> \endverbatim
!>
!> \param[out] M
!> \verbatim
!>          M is INTEGER
!>          The actual number of eigenvalues found. 0 <= M <= N.
!>          (See also the description of INFO=2,3.)
!> \endverbatim
!>
!> \param[out] W
!> \verbatim
!>          W is DOUBLE PRECISION array, dimension (N)
!>          On exit, the first M elements of W will contain the
!>          eigenvalue approximations. DLARRD computes an interval
!>          I_j = (a_j, b_j] that includes eigenvalue j. The eigenvalue
!>          approximation is given as the interval midpoint
!>          W(j)= ( a_j + b_j)/2. The corresponding error is bounded by
!>          WERR(j) = abs( a_j - b_j)/2
!> \endverbatim
!>
!> \param[out] WERR
!> \verbatim
!>          WERR is DOUBLE PRECISION array, dimension (N)
!>          The error bound on the corresponding eigenvalue approximation
!>          in W.
!> \endverbatim
!>
!> \param[out] WL
!> \verbatim
!>          WL is DOUBLE PRECISION
!> \endverbatim
!>
!> \param[out] WU
!> \verbatim
!>          WU is DOUBLE PRECISION
!>          The interval (WL, WU] contains all the wanted eigenvalues.
!>          If RANGE='V', then WL=VL and WU=VU.
!>          If RANGE='A', then WL and WU are the global Gerschgorin bounds
!>                        on the spectrum.
!>          If RANGE='I', then WL and WU are computed by DLAEBZ from the
!>                        index range specified.
!> \endverbatim
!>
!> \param[out] IBLOCK
!> \verbatim
!>          IBLOCK is INTEGER array, dimension (N)
!>          At each row/column j where E(j) is zero or small, the
!>          matrix T is considered to split into a block diagonal
!>          matrix.  On exit, if INFO = 0, IBLOCK(i) specifies to which
!>          block (from 1 to the number of blocks) the eigenvalue W(i)
!>          belongs.  (DLARRD may use the remaining N-M elements as
!>          workspace.)
!> \endverbatim
!>
!> \param[out] INDEXW
!> \verbatim
!>          INDEXW is INTEGER array, dimension (N)
!>          The indices of the eigenvalues within each block (submatrix);
!>          for example, INDEXW(i)= j and IBLOCK(i)=k imply that the
!>          i-th eigenvalue W(i) is the j-th eigenvalue in block k.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (4*N)
!> \endverbatim
!>
!> \param[out] IWORK
!> \verbatim
!>          IWORK is INTEGER array, dimension (3*N)
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value
!>          > 0:  some or all of the eigenvalues failed to converge or
!>                were not computed:
!>                =1 or 3: Bisection failed to converge for some
!>                        eigenvalues; these eigenvalues are flagged by a
!>                        negative block number.  The effect is that the
!>                        eigenvalues may not be as accurate as the
!>                        absolute and relative tolerances.  This is
!>                        generally caused by unexpectedly inaccurate
!>                        arithmetic.
!>                =2 or 3: RANGE='I' only: Not all of the eigenvalues
!>                        IL:IU were found.
!>                        Effect: M < IU+1-IL
!>                        Cause:  non-monotonic arithmetic, causing the
!>                                Sturm sequence to be non-monotonic.
!>                        Cure:   recalculate, using RANGE='A', and pick
!>                                out eigenvalues IL:IU.  In some cases,
!>                                increasing the PARAMETER "FUDGE" may
!>                                make things work.
!>                = 4:    RANGE='I', and the Gershgorin interval
!>                        initially used was too small.  No eigenvalues
!>                        were computed.
!>                        Probable cause: your machine has sloppy
!>                                        floating-point arithmetic.
!>                        Cure: Increase the PARAMETER "FUDGE",
!>                              recompile, and try again.
!> \endverbatim
!
!> \par Internal Parameters:
!  =========================
!>
!> \verbatim
!>  FUDGE   DOUBLE PRECISION, default = 2
!>          A "fudge factor" to widen the Gershgorin intervals.  Ideally,
!>          a value of 1 should work, but on machines with sloppy
!>          arithmetic, this needs to be larger.  The default for
!>          publicly released versions should be large enough to handle
!>          the worst machine around.  Note that this has no effect
!>          on accuracy of the solution.
!> \endverbatim
!>
!> \par Contributors:
!  ==================
!>
!>     W. Kahan, University of California, Berkeley, USA \n
!>     Beresford Parlett, University of California, Berkeley, USA \n
!>     Jim Demmel, University of California, Berkeley, USA \n
!>     Inderjit Dhillon, University of Texas, Austin, USA \n
!>     Osni Marques, LBNL/NERSC, USA \n
!>     Christof Voemel, University of California, Berkeley, USA \n
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
!  =====================================================================
      SUBROUTINE DLARRD(Range,Order,N,Vl,Vu,Il,Iu,Gers,Reltol,D,E,E2,   &
     &                  Pivmin,Nsplit,Isplit,M,W,Werr,Wl,Wu,Iblock,     &
     &                  Indexw,Work,Iwork,Info)
      USE F77KINDS                        
      USE S_DLAEBZ
      USE S_DLAMCH
      USE S_ILAENV
      USE S_LSAME
      IMPLICIT NONE
!*--DLARRD337
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D0 , ONE = 1.0D0 ,        &
     &                              TWO = 2.0D0 , HALF = ONE/TWO ,      &
     &                              FUDGE = TWO
      INTEGER , PARAMETER  ::  ALLRNG = 1 , VALRNG = 2 , INDRNG = 3
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Range
      CHARACTER :: Order
      INTEGER :: N
      REAL(R8KIND) , INTENT(IN) :: Vl
      REAL(R8KIND) , INTENT(IN) :: Vu
      INTEGER , INTENT(IN) :: Il
      INTEGER , INTENT(IN) :: Iu
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: Gers
      REAL(R8KIND) , INTENT(IN) :: Reltol
      REAL(R8KIND) , DIMENSION(*) :: D
      REAL(R8KIND) , DIMENSION(*) :: E
      REAL(R8KIND) , DIMENSION(*) :: E2
      REAL(R8KIND) :: Pivmin
      INTEGER , INTENT(IN) :: Nsplit
      INTEGER , INTENT(IN) , DIMENSION(*) :: Isplit
      INTEGER , INTENT(INOUT) :: M
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: W
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Werr
      REAL(R8KIND) , INTENT(INOUT) :: Wl
      REAL(R8KIND) , INTENT(INOUT) :: Wu
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iblock
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Indexw
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      REAL(R8KIND) :: atoli , eps , gl , gu , rtoli , tmp1 , tmp2 ,     &
     &                tnorm , uflow , wkill , wlu , wul
      INTEGER :: i , ib , ibegin , idiscl , idiscu , ie , iend , iinfo ,&
     &           im , in , ioff , iout , irange , itmax , itmp1 ,       &
     &           itmp2 , iw , iwoff , j , jblk , jdisc , je , jee , nb ,&
     &           nwl , nwu
      INTEGER , DIMENSION(1) :: idumma
      LOGICAL :: ncnvrg , toofew
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
      ELSE
         irange = 0
      ENDIF
!
!     Check for Errors
!
      IF ( irange<=0 ) THEN
         Info = -1
      ELSEIF ( .NOT.(LSAME(Order,'B') .OR. LSAME(Order,'E')) ) THEN
         Info = -2
      ELSEIF ( N<0 ) THEN
         Info = -3
      ELSEIF ( irange==VALRNG ) THEN
         IF ( Vl>=Vu ) Info = -5
      ELSEIF ( irange==INDRNG .AND. (Il<1 .OR. Il>MAX(1,N)) ) THEN
         Info = -6
      ELSEIF ( irange==INDRNG .AND. (Iu<MIN(N,Il) .OR. Iu>N) ) THEN
         Info = -7
      ENDIF
!
      IF ( Info/=0 ) RETURN
 
!     Initialize error flags
      Info = 0
      ncnvrg = .FALSE.
      toofew = .FALSE.
 
!     Quick return if possible
      M = 0
      IF ( N==0 ) RETURN
 
!     Simplification:
      IF ( irange==INDRNG .AND. Il==1 .AND. Iu==N ) irange = 1
 
!     Get machine constants
      eps = DLAMCH('P')
      uflow = DLAMCH('U')
 
 
!     Special Case when N=1
!     Treat case of 1x1 matrix for quick return
      IF ( N==1 ) THEN
         IF ( (irange==ALLRNG) .OR.                                     &
     &        ((irange==VALRNG) .AND. (D(1)>Vl) .AND. (D(1)<=Vu)) .OR.  &
     &        ((irange==INDRNG) .AND. (Il==1) .AND. (Iu==1)) ) THEN
            M = 1
            W(1) = D(1)
!           The computation error of the eigenvalue is zero
            Werr(1) = ZERO
            Iblock(1) = 1
            Indexw(1) = 1
         ENDIF
         RETURN
      ENDIF
 
!     NB is the minimum vector length for vector bisection, or 0
!     if only scalar is to be done.
      nb = ILAENV(1,'DSTEBZ',' ',N,-1,-1,-1)
      IF ( nb<=1 ) nb = 0
 
!     Find global spectral radius
      gl = D(1)
      gu = D(1)
      DO i = 1 , N
         gl = MIN(gl,Gers(2*i-1))
         gu = MAX(gu,Gers(2*i))
      ENDDO
!     Compute global Gerschgorin bounds and spectral diameter
      tnorm = MAX(ABS(gl),ABS(gu))
      gl = gl - FUDGE*tnorm*eps*N - FUDGE*TWO*Pivmin
      gu = gu + FUDGE*tnorm*eps*N + FUDGE*TWO*Pivmin
!     [JAN/28/2009] remove the line below since SPDIAM variable not use
!     SPDIAM = GU - GL
!     Input arguments for DLAEBZ:
!     The relative tolerance.  An interval (a,b] lies within
!     "relative tolerance" if  b-a < RELTOL*max(|a|,|b|),
      rtoli = Reltol
!     Set the absolute tolerance for interval convergence to zero to force
!     interval convergence based on relative size of the interval.
!     This is dangerous because intervals might not converge when RELTOL is
!     small. But at least a very small number should be selected so that for
!     strongly graded matrices, the code can get relatively accurate
!     eigenvalues.
      atoli = FUDGE*TWO*uflow + FUDGE*TWO*Pivmin
 
      IF ( irange==INDRNG ) THEN
 
!        RANGE='I': Compute an interval containing eigenvalues
!        IL through IU. The initial interval [GL,GU] from the global
!        Gerschgorin bounds GL and GU is refined by DLAEBZ.
         itmax = INT((LOG(tnorm+Pivmin)-LOG(Pivmin))/LOG(TWO)) + 2
         Work(N+1) = gl
         Work(N+2) = gl
         Work(N+3) = gu
         Work(N+4) = gu
         Work(N+5) = gl
         Work(N+6) = gu
         Iwork(1) = -1
         Iwork(2) = -1
         Iwork(3) = N + 1
         Iwork(4) = N + 1
         Iwork(5) = Il - 1
         Iwork(6) = Iu
!
         CALL DLAEBZ(3,itmax,N,2,2,nb,atoli,rtoli,Pivmin,D,E,E2,Iwork(5)&
     &               ,Work(N+1),Work(N+5),iout,Iwork,W,Iblock,iinfo)
         IF ( iinfo/=0 ) THEN
            Info = iinfo
            RETURN
         ENDIF
!        On exit, output intervals may not be ordered by ascending negcount
         IF ( Iwork(6)==Iu ) THEN
            Wl = Work(N+1)
            wlu = Work(N+3)
            nwl = Iwork(1)
            Wu = Work(N+4)
            wul = Work(N+2)
            nwu = Iwork(4)
         ELSE
            Wl = Work(N+2)
            wlu = Work(N+4)
            nwl = Iwork(2)
            Wu = Work(N+3)
            wul = Work(N+1)
            nwu = Iwork(3)
         ENDIF
!        On exit, the interval [WL, WLU] contains a value with negcount NWL,
!        and [WUL, WU] contains a value with negcount NWU.
         IF ( nwl<0 .OR. nwl>=N .OR. nwu<1 .OR. nwu>N ) THEN
            Info = 4
            RETURN
         ENDIF
 
      ELSEIF ( irange==VALRNG ) THEN
         Wl = Vl
         Wu = Vu
 
      ELSEIF ( irange==ALLRNG ) THEN
         Wl = gl
         Wu = gu
      ENDIF
 
 
 
!     Find Eigenvalues -- Loop Over blocks and recompute NWL and NWU.
!     NWL accumulates the number of eigenvalues .le. WL,
!     NWU accumulates the number of eigenvalues .le. WU
      M = 0
      iend = 0
      Info = 0
      nwl = 0
      nwu = 0
!
      DO jblk = 1 , Nsplit
         ioff = iend
         ibegin = ioff + 1
         iend = Isplit(jblk)
         in = iend - ioff
!
         IF ( in==1 ) THEN
!           1x1 block
            IF ( Wl>=D(ibegin)-Pivmin ) nwl = nwl + 1
            IF ( Wu>=D(ibegin)-Pivmin ) nwu = nwu + 1
            IF ( irange==ALLRNG .OR.                                    &
     &           (Wl<D(ibegin)-Pivmin .AND. Wu>=D(ibegin)-Pivmin) ) THEN
               M = M + 1
               W(M) = D(ibegin)
               Werr(M) = ZERO
!              The gap for a single block doesn't matter for the later
!              algorithm and is assigned an arbitrary large value
               Iblock(M) = jblk
               Indexw(M) = 1
            ENDIF
 
!        Disabled 2x2 case because of a failure on the following matrix
!        RANGE = 'I', IL = IU = 4
!          Original Tridiagonal, d = [
!           -0.150102010615740E+00
!           -0.849897989384260E+00
!           -0.128208148052635E-15
!            0.128257718286320E-15
!          ];
!          e = [
!           -0.357171383266986E+00
!           -0.180411241501588E-15
!           -0.175152352710251E-15
!          ];
!
!         ELSE IF( IN.EQ.2 ) THEN
!*           2x2 block
!            DISC = SQRT( (HALF*(D(IBEGIN)-D(IEND)))**2 + E(IBEGIN)**2 )
!            TMP1 = HALF*(D(IBEGIN)+D(IEND))
!            L1 = TMP1 - DISC
!            IF( WL.GE. L1-PIVMIN )
!     $         NWL = NWL + 1
!            IF( WU.GE. L1-PIVMIN )
!     $         NWU = NWU + 1
!            IF( IRANGE.EQ.ALLRNG .OR. ( WL.LT.L1-PIVMIN .AND. WU.GE.
!     $          L1-PIVMIN ) ) THEN
!               M = M + 1
!               W( M ) = L1
!*              The uncertainty of eigenvalues of a 2x2 matrix is very small
!               WERR( M ) = EPS * ABS( W( M ) ) * TWO
!               IBLOCK( M ) = JBLK
!               INDEXW( M ) = 1
!            ENDIF
!            L2 = TMP1 + DISC
!            IF( WL.GE. L2-PIVMIN )
!     $         NWL = NWL + 1
!            IF( WU.GE. L2-PIVMIN )
!     $         NWU = NWU + 1
!            IF( IRANGE.EQ.ALLRNG .OR. ( WL.LT.L2-PIVMIN .AND. WU.GE.
!     $          L2-PIVMIN ) ) THEN
!               M = M + 1
!               W( M ) = L2
!*              The uncertainty of eigenvalues of a 2x2 matrix is very small
!               WERR( M ) = EPS * ABS( W( M ) ) * TWO
!               IBLOCK( M ) = JBLK
!               INDEXW( M ) = 2
!            ENDIF
         ELSE
!           General Case - block of size IN >= 2
!           Compute local Gerschgorin interval and use it as the initial
!           interval for DLAEBZ
            gu = D(ibegin)
            gl = D(ibegin)
            tmp1 = ZERO
 
            DO j = ibegin , iend
               gl = MIN(gl,Gers(2*j-1))
               gu = MAX(gu,Gers(2*j))
            ENDDO
!           [JAN/28/2009]
!           change SPDIAM by TNORM in lines 2 and 3 thereafter
!           line 1: remove computation of SPDIAM (not useful anymore)
!           SPDIAM = GU - GL
!           GL = GL - FUDGE*SPDIAM*EPS*IN - FUDGE*PIVMIN
!           GU = GU + FUDGE*SPDIAM*EPS*IN + FUDGE*PIVMIN
            gl = gl - FUDGE*tnorm*eps*in - FUDGE*Pivmin
            gu = gu + FUDGE*tnorm*eps*in + FUDGE*Pivmin
!
            IF ( irange>1 ) THEN
               IF ( gu<Wl ) THEN
!                 the local block contains none of the wanted eigenvalues
                  nwl = nwl + in
                  nwu = nwu + in
                  CYCLE
               ENDIF
!              refine search interval if possible, only range (WL,WU] matters
               gl = MAX(gl,Wl)
               gu = MIN(gu,Wu)
               IF ( gl>=gu ) CYCLE
            ENDIF
 
!           Find negcount of initial interval boundaries GL and GU
            Work(N+1) = gl
            Work(N+in+1) = gu
            CALL DLAEBZ(1,0,in,in,1,nb,atoli,rtoli,Pivmin,D(ibegin),    &
     &                  E(ibegin),E2(ibegin),idumma,Work(N+1),          &
     &                  Work(N+2*in+1),im,Iwork,W(M+1),Iblock(M+1),     &
     &                  iinfo)
            IF ( iinfo/=0 ) THEN
               Info = iinfo
               RETURN
            ENDIF
!
            nwl = nwl + Iwork(1)
            nwu = nwu + Iwork(in+1)
            iwoff = M - Iwork(1)
 
!           Compute Eigenvalues
            itmax = INT((LOG(gu-gl+Pivmin)-LOG(Pivmin))/LOG(TWO)) + 2
            CALL DLAEBZ(2,itmax,in,in,1,nb,atoli,rtoli,Pivmin,D(ibegin),&
     &                  E(ibegin),E2(ibegin),idumma,Work(N+1),          &
     &                  Work(N+2*in+1),iout,Iwork,W(M+1),Iblock(M+1),   &
     &                  iinfo)
            IF ( iinfo/=0 ) THEN
               Info = iinfo
               RETURN
            ENDIF
!
!           Copy eigenvalues into W and IBLOCK
!           Use -JBLK for block number for unconverged eigenvalues.
!           Loop over the number of output intervals from DLAEBZ
            DO j = 1 , iout
!              eigenvalue approximation is middle point of interval
               tmp1 = HALF*(Work(j+N)+Work(j+in+N))
!              semi length of error interval
               tmp2 = HALF*ABS(Work(j+N)-Work(j+in+N))
               IF ( j>iout-iinfo ) THEN
!                 Flag non-convergence.
                  ncnvrg = .TRUE.
                  ib = -jblk
               ELSE
                  ib = jblk
               ENDIF
               DO je = Iwork(j) + 1 + iwoff , Iwork(j+in) + iwoff
                  W(je) = tmp1
                  Werr(je) = tmp2
                  Indexw(je) = je - iwoff
                  Iblock(je) = ib
               ENDDO
            ENDDO
!
            M = M + im
         ENDIF
      ENDDO
 
!     If RANGE='I', then (WL,WU) contains eigenvalues NWL+1,...,NWU
!     If NWL+1 < IL or NWU > IU, discard extra eigenvalues.
      IF ( irange==INDRNG ) THEN
         idiscl = Il - 1 - nwl
         idiscu = nwu - Iu
!
         IF ( idiscl>0 ) THEN
            im = 0
            DO je = 1 , M
!              Remove some of the smallest eigenvalues from the left so that
!              at the end IDISCL =0. Move all eigenvalues up to the left.
               IF ( W(je)<=wlu .AND. idiscl>0 ) THEN
                  idiscl = idiscl - 1
               ELSE
                  im = im + 1
                  W(im) = W(je)
                  Werr(im) = Werr(je)
                  Indexw(im) = Indexw(je)
                  Iblock(im) = Iblock(je)
               ENDIF
            ENDDO
            M = im
         ENDIF
         IF ( idiscu>0 ) THEN
!           Remove some of the largest eigenvalues from the right so that
!           at the end IDISCU =0. Move all eigenvalues up to the left.
            im = M + 1
            DO je = M , 1 , -1
               IF ( W(je)>=wul .AND. idiscu>0 ) THEN
                  idiscu = idiscu - 1
               ELSE
                  im = im - 1
                  W(im) = W(je)
                  Werr(im) = Werr(je)
                  Indexw(im) = Indexw(je)
                  Iblock(im) = Iblock(je)
               ENDIF
            ENDDO
            jee = 0
            DO je = im , M
               jee = jee + 1
               W(jee) = W(je)
               Werr(jee) = Werr(je)
               Indexw(jee) = Indexw(je)
               Iblock(jee) = Iblock(je)
            ENDDO
            M = M - im + 1
         ENDIF
 
         IF ( idiscl>0 .OR. idiscu>0 ) THEN
!           Code to deal with effects of bad arithmetic. (If N(w) is
!           monotone non-decreasing, this should never happen.)
!           Some low eigenvalues to be discarded are not in (WL,WLU],
!           or high eigenvalues to be discarded are not in (WUL,WU]
!           so just kill off the smallest IDISCL/largest IDISCU
!           eigenvalues, by marking the corresponding IBLOCK = 0
            IF ( idiscl>0 ) THEN
               wkill = Wu
               DO jdisc = 1 , idiscl
                  iw = 0
                  DO je = 1 , M
                     IF ( Iblock(je)/=0 .AND. (W(je)<wkill .OR. iw==0) )&
     &                    THEN
                        iw = je
                        wkill = W(je)
                     ENDIF
                  ENDDO
                  Iblock(iw) = 0
               ENDDO
            ENDIF
            IF ( idiscu>0 ) THEN
               wkill = Wl
               DO jdisc = 1 , idiscu
                  iw = 0
                  DO je = 1 , M
                     IF ( Iblock(je)/=0 .AND. (W(je)>=wkill .OR. iw==0) &
     &                    ) THEN
                        iw = je
                        wkill = W(je)
                     ENDIF
                  ENDDO
                  Iblock(iw) = 0
               ENDDO
            ENDIF
!           Now erase all eigenvalues with IBLOCK set to zero
            im = 0
            DO je = 1 , M
               IF ( Iblock(je)/=0 ) THEN
                  im = im + 1
                  W(im) = W(je)
                  Werr(im) = Werr(je)
                  Indexw(im) = Indexw(je)
                  Iblock(im) = Iblock(je)
               ENDIF
            ENDDO
            M = im
         ENDIF
         IF ( idiscl<0 .OR. idiscu<0 ) toofew = .TRUE.
      ENDIF
!
      IF ( (irange==ALLRNG .AND. M/=N) .OR.                             &
     &     (irange==INDRNG .AND. M/=Iu-Il+1) ) toofew = .TRUE.
 
!     If ORDER='B', do nothing the eigenvalues are already sorted by
!        block.
!     If ORDER='E', sort the eigenvalues from smallest to largest
 
      IF ( LSAME(Order,'E') .AND. Nsplit>1 ) THEN
         DO je = 1 , M - 1
            ie = 0
            tmp1 = W(je)
            DO j = je + 1 , M
               IF ( W(j)<tmp1 ) THEN
                  ie = j
                  tmp1 = W(j)
               ENDIF
            ENDDO
            IF ( ie/=0 ) THEN
               tmp2 = Werr(ie)
               itmp1 = Iblock(ie)
               itmp2 = Indexw(ie)
               W(ie) = W(je)
               Werr(ie) = Werr(je)
               Iblock(ie) = Iblock(je)
               Indexw(ie) = Indexw(je)
               W(je) = tmp1
               Werr(je) = tmp2
               Iblock(je) = itmp1
               Indexw(je) = itmp2
            ENDIF
         ENDDO
      ENDIF
!
      Info = 0
      IF ( ncnvrg ) Info = Info + 1
      IF ( toofew ) Info = Info + 2
!
!     End of DLARRD
!
      END SUBROUTINE DLARRD
