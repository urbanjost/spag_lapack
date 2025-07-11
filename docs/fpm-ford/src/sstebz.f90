!*==sstebz.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b SSTEBZ
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SSTEBZ + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sstebz.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sstebz.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sstebz.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SSTEBZ( RANGE, ORDER, N, VL, VU, IL, IU, ABSTOL, D, E,
!                          M, NSPLIT, W, IBLOCK, ISPLIT, WORK, IWORK,
!                          INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          ORDER, RANGE
!       INTEGER            IL, INFO, IU, M, N, NSPLIT
!       REAL               ABSTOL, VL, VU
!       ..
!       .. Array Arguments ..
!       INTEGER            IBLOCK( * ), ISPLIT( * ), IWORK( * )
!       REAL               D( * ), E( * ), W( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SSTEBZ computes the eigenvalues of a symmetric tridiagonal
!> matrix T.  The user may ask for all eigenvalues, all eigenvalues
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
!>          VL is REAL
!>
!>          If RANGE='V', the lower bound of the interval to
!>          be searched for eigenvalues.  Eigenvalues less than or equal
!>          to VL, or greater than VU, will not be returned.  VL < VU.
!>          Not referenced if RANGE = 'A' or 'I'.
!> \endverbatim
!>
!> \param[in] VU
!> \verbatim
!>          VU is REAL
!>
!>          If RANGE='V', the upper bound of the interval to
!>          be searched for eigenvalues.  Eigenvalues less than or equal
!>          to VL, or greater than VU, will not be returned.  VL < VU.
!>          Not referenced if RANGE = 'A' or 'I'.
!> \endverbatim
!>
!> \param[in] IL
!> \verbatim
!>          IL is INTEGER
!>
!>          If RANGE='I', the index of the
!>          smallest eigenvalue to be returned.
!>          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0.
!>          Not referenced if RANGE = 'A' or 'V'.
!> \endverbatim
!>
!> \param[in] IU
!> \verbatim
!>          IU is INTEGER
!>
!>          If RANGE='I', the index of the
!>          largest eigenvalue to be returned.
!>          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0.
!>          Not referenced if RANGE = 'A' or 'V'.
!> \endverbatim
!>
!> \param[in] ABSTOL
!> \verbatim
!>          ABSTOL is REAL
!>          The absolute tolerance for the eigenvalues.  An eigenvalue
!>          (or cluster) is considered to be located if it has been
!>          determined to lie in an interval whose width is ABSTOL or
!>          less.  If ABSTOL is less than or equal to zero, then ULP*|T|
!>          will be used, where |T| means the 1-norm of T.
!>
!>          Eigenvalues will be computed most accurately when ABSTOL is
!>          set to twice the underflow threshold 2*SLAMCH('S'), not zero.
!> \endverbatim
!>
!> \param[in] D
!> \verbatim
!>          D is REAL array, dimension (N)
!>          The n diagonal elements of the tridiagonal matrix T.
!> \endverbatim
!>
!> \param[in] E
!> \verbatim
!>          E is REAL array, dimension (N-1)
!>          The (n-1) off-diagonal elements of the tridiagonal matrix T.
!> \endverbatim
!>
!> \param[out] M
!> \verbatim
!>          M is INTEGER
!>          The actual number of eigenvalues found. 0 <= M <= N.
!>          (See also the description of INFO=2,3.)
!> \endverbatim
!>
!> \param[out] NSPLIT
!> \verbatim
!>          NSPLIT is INTEGER
!>          The number of diagonal blocks in the matrix T.
!>          1 <= NSPLIT <= N.
!> \endverbatim
!>
!> \param[out] W
!> \verbatim
!>          W is REAL array, dimension (N)
!>          On exit, the first M elements of W will contain the
!>          eigenvalues.  (SSTEBZ may use the remaining N-M elements as
!>          workspace.)
!> \endverbatim
!>
!> \param[out] IBLOCK
!> \verbatim
!>          IBLOCK is INTEGER array, dimension (N)
!>          At each row/column j where E(j) is zero or small, the
!>          matrix T is considered to split into a block diagonal
!>          matrix.  On exit, if INFO = 0, IBLOCK(i) specifies to which
!>          block (from 1 to the number of blocks) the eigenvalue W(i)
!>          belongs.  (SSTEBZ may use the remaining N-M elements as
!>          workspace.)
!> \endverbatim
!>
!> \param[out] ISPLIT
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
!> \param[out] WORK
!> \verbatim
!>          WORK is REAL array, dimension (4*N)
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
!>  RELFAC  REAL, default = 2.0e0
!>          The relative tolerance.  An interval (a,b] lies within
!>          "relative tolerance" if  b-a < RELFAC*ulp*max(|a|,|b|),
!>          where "ulp" is the machine precision (distance from 1 to
!>          the next larger floating point number.)
!>
!>  FUDGE   REAL, default = 2
!>          A "fudge factor" to widen the Gershgorin intervals.  Ideally,
!>          a value of 1 should work, but on machines with sloppy
!>          arithmetic, this needs to be larger.  The default for
!>          publicly released versions should be large enough to handle
!>          the worst machine around.  Note that this has no effect
!>          on accuracy of the solution.
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
!> \ingroup auxOTHERcomputational
!
!  =====================================================================
      SUBROUTINE SSTEBZ(Range,Order,N,Vl,Vu,Il,Iu,Abstol,D,E,M,Nsplit,W,&
     &                  Iblock,Isplit,Work,Iwork,Info)
      IMPLICIT NONE
!*--SSTEBZ276
!
!  -- LAPACK computational routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     June 2016
!
!     .. Scalar Arguments ..
      CHARACTER Order , Range
      INTEGER Il , Info , Iu , M , N , Nsplit
      REAL Abstol , Vl , Vu
!     ..
!     .. Array Arguments ..
      INTEGER Iblock(*) , Isplit(*) , Iwork(*)
      REAL D(*) , E(*) , W(*) , Work(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ZERO , ONE , TWO , HALF
      PARAMETER (ZERO=0.0E0,ONE=1.0E0,TWO=2.0E0,HALF=1.0E0/TWO)
      REAL FUDGE , RELFAC
      PARAMETER (FUDGE=2.1E0,RELFAC=2.0E0)
!     ..
!     .. Local Scalars ..
      LOGICAL ncnvrg , toofew
      INTEGER ib , ibegin , idiscl , idiscu , ie , iend , iinfo , im ,  &
     &        in , ioff , iorder , iout , irange , itmax , itmp1 , iw , &
     &        iwoff , j , jb , jdisc , je , nb , nwl , nwu
      REAL atoli , bnorm , gl , gu , pivmin , rtoli , safemn , tmp1 ,   &
     &     tmp2 , tnorm , ulp , wkill , wl , wlu , wu , wul
!     ..
!     .. Local Arrays ..
      INTEGER idumma(1)
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      INTEGER ILAENV
      REAL SLAMCH
      EXTERNAL LSAME , ILAENV , SLAMCH
!     ..
!     .. External Subroutines ..
      EXTERNAL SLAEBZ , XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , INT , LOG , MAX , MIN , SQRT
!     ..
!     .. Executable Statements ..
!
      Info = 0
!
!     Decode RANGE
!
      IF ( LSAME(Range,'A') ) THEN
         irange = 1
      ELSEIF ( LSAME(Range,'V') ) THEN
         irange = 2
      ELSEIF ( LSAME(Range,'I') ) THEN
         irange = 3
      ELSE
         irange = 0
      ENDIF
!
!     Decode ORDER
!
      IF ( LSAME(Order,'B') ) THEN
         iorder = 2
      ELSEIF ( LSAME(Order,'E') ) THEN
         iorder = 1
      ELSE
         iorder = 0
      ENDIF
!
!     Check for Errors
!
      IF ( irange<=0 ) THEN
         Info = -1
      ELSEIF ( iorder<=0 ) THEN
         Info = -2
      ELSEIF ( N<0 ) THEN
         Info = -3
      ELSEIF ( irange==2 ) THEN
         IF ( Vl>=Vu ) Info = -5
      ELSEIF ( irange==3 .AND. (Il<1 .OR. Il>MAX(1,N)) ) THEN
         Info = -6
      ELSEIF ( irange==3 .AND. (Iu<MIN(N,Il) .OR. Iu>N) ) THEN
         Info = -7
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('SSTEBZ',-Info)
         RETURN
      ENDIF
!
!     Initialize error flags
!
      Info = 0
      ncnvrg = .FALSE.
      toofew = .FALSE.
!
!     Quick return if possible
!
      M = 0
      IF ( N==0 ) RETURN
!
!     Simplifications:
!
      IF ( irange==3 .AND. Il==1 .AND. Iu==N ) irange = 1
!
!     Get machine constants
!     NB is the minimum vector length for vector bisection, or 0
!     if only scalar is to be done.
!
      safemn = SLAMCH('S')
      ulp = SLAMCH('P')
      rtoli = ulp*RELFAC
      nb = ILAENV(1,'SSTEBZ',' ',N,-1,-1,-1)
      IF ( nb<=1 ) nb = 0
!
!     Special Case when N=1
!
      IF ( N==1 ) THEN
         Nsplit = 1
         Isplit(1) = 1
         IF ( irange==2 .AND. (Vl>=D(1) .OR. Vu<D(1)) ) THEN
            M = 0
         ELSE
            W(1) = D(1)
            Iblock(1) = 1
            M = 1
         ENDIF
         RETURN
      ENDIF
!
!     Compute Splitting Points
!
      Nsplit = 1
      Work(N) = ZERO
      pivmin = ONE
!
      DO j = 2 , N
         tmp1 = E(j-1)**2
         IF ( ABS(D(j)*D(j-1))*ulp**2+safemn>tmp1 ) THEN
            Isplit(Nsplit) = j - 1
            Nsplit = Nsplit + 1
            Work(j-1) = ZERO
         ELSE
            Work(j-1) = tmp1
            pivmin = MAX(pivmin,tmp1)
         ENDIF
      ENDDO
      Isplit(Nsplit) = N
      pivmin = pivmin*safemn
!
!     Compute Interval and ATOLI
!
      IF ( irange==3 ) THEN
!
!        RANGE='I': Compute the interval containing eigenvalues
!                   IL through IU.
!
!        Compute Gershgorin interval for entire (split) matrix
!        and use it as the initial interval
!
         gu = D(1)
         gl = D(1)
         tmp1 = ZERO
!
         DO j = 1 , N - 1
            tmp2 = SQRT(Work(j))
            gu = MAX(gu,D(j)+tmp1+tmp2)
            gl = MIN(gl,D(j)-tmp1-tmp2)
            tmp1 = tmp2
         ENDDO
!
         gu = MAX(gu,D(N)+tmp1)
         gl = MIN(gl,D(N)-tmp1)
         tnorm = MAX(ABS(gl),ABS(gu))
         gl = gl - FUDGE*tnorm*ulp*N - FUDGE*TWO*pivmin
         gu = gu + FUDGE*tnorm*ulp*N + FUDGE*pivmin
!
!        Compute Iteration parameters
!
         itmax = INT((LOG(tnorm+pivmin)-LOG(pivmin))/LOG(TWO)) + 2
         IF ( Abstol<=ZERO ) THEN
            atoli = ulp*tnorm
         ELSE
            atoli = Abstol
         ENDIF
!
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
         CALL SLAEBZ(3,itmax,N,2,2,nb,atoli,rtoli,pivmin,D,E,Work,      &
     &               Iwork(5),Work(N+1),Work(N+5),iout,Iwork,W,Iblock,  &
     &               iinfo)
!
         IF ( Iwork(6)==Iu ) THEN
            wl = Work(N+1)
            wlu = Work(N+3)
            nwl = Iwork(1)
            wu = Work(N+4)
            wul = Work(N+2)
            nwu = Iwork(4)
         ELSE
            wl = Work(N+2)
            wlu = Work(N+4)
            nwl = Iwork(2)
            wu = Work(N+3)
            wul = Work(N+1)
            nwu = Iwork(3)
         ENDIF
!
         IF ( nwl<0 .OR. nwl>=N .OR. nwu<1 .OR. nwu>N ) THEN
            Info = 4
            RETURN
         ENDIF
      ELSE
!
!        RANGE='A' or 'V' -- Set ATOLI
!
         tnorm = MAX(ABS(D(1))+ABS(E(1)),ABS(D(N))+ABS(E(N-1)))
!
         DO j = 2 , N - 1
            tnorm = MAX(tnorm,ABS(D(j))+ABS(E(j-1))+ABS(E(j)))
         ENDDO
!
         IF ( Abstol<=ZERO ) THEN
            atoli = ulp*tnorm
         ELSE
            atoli = Abstol
         ENDIF
!
         IF ( irange==2 ) THEN
            wl = Vl
            wu = Vu
         ELSE
            wl = ZERO
            wu = ZERO
         ENDIF
      ENDIF
!
!     Find Eigenvalues -- Loop Over Blocks and recompute NWL and NWU.
!     NWL accumulates the number of eigenvalues .le. WL,
!     NWU accumulates the number of eigenvalues .le. WU
!
      M = 0
      iend = 0
      Info = 0
      nwl = 0
      nwu = 0
!
      DO jb = 1 , Nsplit
         ioff = iend
         ibegin = ioff + 1
         iend = Isplit(jb)
         in = iend - ioff
!
         IF ( in==1 ) THEN
!
!           Special Case -- IN=1
!
            IF ( irange==1 .OR. wl>=D(ibegin)-pivmin ) nwl = nwl + 1
            IF ( irange==1 .OR. wu>=D(ibegin)-pivmin ) nwu = nwu + 1
            IF ( irange==1 .OR.                                         &
     &           (wl<D(ibegin)-pivmin .AND. wu>=D(ibegin)-pivmin) ) THEN
               M = M + 1
               W(M) = D(ibegin)
               Iblock(M) = jb
            ENDIF
         ELSE
!
!           General Case -- IN > 1
!
!           Compute Gershgorin Interval
!           and use it as the initial interval
!
            gu = D(ibegin)
            gl = D(ibegin)
            tmp1 = ZERO
!
            DO j = ibegin , iend - 1
               tmp2 = ABS(E(j))
               gu = MAX(gu,D(j)+tmp1+tmp2)
               gl = MIN(gl,D(j)-tmp1-tmp2)
               tmp1 = tmp2
            ENDDO
!
            gu = MAX(gu,D(iend)+tmp1)
            gl = MIN(gl,D(iend)-tmp1)
            bnorm = MAX(ABS(gl),ABS(gu))
            gl = gl - FUDGE*bnorm*ulp*in - FUDGE*pivmin
            gu = gu + FUDGE*bnorm*ulp*in + FUDGE*pivmin
!
!           Compute ATOLI for the current submatrix
!
            IF ( Abstol<=ZERO ) THEN
               atoli = ulp*MAX(ABS(gl),ABS(gu))
            ELSE
               atoli = Abstol
            ENDIF
!
            IF ( irange>1 ) THEN
               IF ( gu<wl ) THEN
                  nwl = nwl + in
                  nwu = nwu + in
                  CYCLE
               ENDIF
               gl = MAX(gl,wl)
               gu = MIN(gu,wu)
               IF ( gl>=gu ) CYCLE
            ENDIF
!
!           Set Up Initial Interval
!
            Work(N+1) = gl
            Work(N+in+1) = gu
            CALL SLAEBZ(1,0,in,in,1,nb,atoli,rtoli,pivmin,D(ibegin),    &
     &                  E(ibegin),Work(ibegin),idumma,Work(N+1),        &
     &                  Work(N+2*in+1),im,Iwork,W(M+1),Iblock(M+1),     &
     &                  iinfo)
!
            nwl = nwl + Iwork(1)
            nwu = nwu + Iwork(in+1)
            iwoff = M - Iwork(1)
!
!           Compute Eigenvalues
!
            itmax = INT((LOG(gu-gl+pivmin)-LOG(pivmin))/LOG(TWO)) + 2
            CALL SLAEBZ(2,itmax,in,in,1,nb,atoli,rtoli,pivmin,D(ibegin),&
     &                  E(ibegin),Work(ibegin),idumma,Work(N+1),        &
     &                  Work(N+2*in+1),iout,Iwork,W(M+1),Iblock(M+1),   &
     &                  iinfo)
!
!           Copy Eigenvalues Into W and IBLOCK
!           Use -JB for block number for unconverged eigenvalues.
!
            DO j = 1 , iout
               tmp1 = HALF*(Work(j+N)+Work(j+in+N))
!
!              Flag non-convergence.
!
               IF ( j>iout-iinfo ) THEN
                  ncnvrg = .TRUE.
                  ib = -jb
               ELSE
                  ib = jb
               ENDIF
               DO je = Iwork(j) + 1 + iwoff , Iwork(j+in) + iwoff
                  W(je) = tmp1
                  Iblock(je) = ib
               ENDDO
            ENDDO
!
            M = M + im
         ENDIF
      ENDDO
!
!     If RANGE='I', then (WL,WU) contains eigenvalues NWL+1,...,NWU
!     If NWL+1 < IL or NWU > IU, discard extra eigenvalues.
!
      IF ( irange==3 ) THEN
         im = 0
         idiscl = Il - 1 - nwl
         idiscu = nwu - Iu
!
         IF ( idiscl>0 .OR. idiscu>0 ) THEN
            DO je = 1 , M
               IF ( W(je)<=wlu .AND. idiscl>0 ) THEN
                  idiscl = idiscl - 1
               ELSEIF ( W(je)>=wul .AND. idiscu>0 ) THEN
                  idiscu = idiscu - 1
               ELSE
                  im = im + 1
                  W(im) = W(je)
                  Iblock(im) = Iblock(je)
               ENDIF
            ENDDO
            M = im
         ENDIF
         IF ( idiscl>0 .OR. idiscu>0 ) THEN
!
!           Code to deal with effects of bad arithmetic:
!           Some low eigenvalues to be discarded are not in (WL,WLU],
!           or high eigenvalues to be discarded are not in (WUL,WU]
!           so just kill off the smallest IDISCL/largest IDISCU
!           eigenvalues, by simply finding the smallest/largest
!           eigenvalue(s).
!
!           (If N(w) is monotone non-decreasing, this should never
!               happen.)
!
            IF ( idiscl>0 ) THEN
               wkill = wu
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
!
               wkill = wl
               DO jdisc = 1 , idiscu
                  iw = 0
                  DO je = 1 , M
                     IF ( Iblock(je)/=0 .AND. (W(je)>wkill .OR. iw==0) )&
     &                    THEN
                        iw = je
                        wkill = W(je)
                     ENDIF
                  ENDDO
                  Iblock(iw) = 0
               ENDDO
            ENDIF
            im = 0
            DO je = 1 , M
               IF ( Iblock(je)/=0 ) THEN
                  im = im + 1
                  W(im) = W(je)
                  Iblock(im) = Iblock(je)
               ENDIF
            ENDDO
            M = im
         ENDIF
         IF ( idiscl<0 .OR. idiscu<0 ) toofew = .TRUE.
      ENDIF
!
!     If ORDER='B', do nothing -- the eigenvalues are already sorted
!        by block.
!     If ORDER='E', sort the eigenvalues from smallest to largest
!
      IF ( iorder==1 .AND. Nsplit>1 ) THEN
         DO je = 1 , M - 1
            ie = 0
            tmp1 = W(je)
            DO j = je + 1 , M
               IF ( W(j)<tmp1 ) THEN
                  ie = j
                  tmp1 = W(j)
               ENDIF
            ENDDO
!
            IF ( ie/=0 ) THEN
               itmp1 = Iblock(ie)
               W(ie) = W(je)
               Iblock(ie) = Iblock(je)
               W(je) = tmp1
               Iblock(je) = itmp1
            ENDIF
         ENDDO
      ENDIF
!
      Info = 0
      IF ( ncnvrg ) Info = Info + 1
      IF ( toofew ) Info = Info + 2
!
!     End of SSTEBZ
!
      END SUBROUTINE SSTEBZ
