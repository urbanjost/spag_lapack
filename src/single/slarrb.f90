!*==slarrb.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b SLARRB provides limited bisection to locate eigenvalues for more accuracy.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SLARRB + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slarrb.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slarrb.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slarrb.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SLARRB( N, D, LLD, IFIRST, ILAST, RTOL1,
!                          RTOL2, OFFSET, W, WGAP, WERR, WORK, IWORK,
!                          PIVMIN, SPDIAM, TWIST, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            IFIRST, ILAST, INFO, N, OFFSET, TWIST
!       REAL               PIVMIN, RTOL1, RTOL2, SPDIAM
!       ..
!       .. Array Arguments ..
!       INTEGER            IWORK( * )
!       REAL               D( * ), LLD( * ), W( * ),
!      $                   WERR( * ), WGAP( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> Given the relatively robust representation(RRR) L D L^T, SLARRB
!> does "limited" bisection to refine the eigenvalues of L D L^T,
!> W( IFIRST-OFFSET ) through W( ILAST-OFFSET ), to more accuracy. Initial
!> guesses for these eigenvalues are input in W, the corresponding estimate
!> of the error in these guesses and their gaps are input in WERR
!> and WGAP, respectively. During bisection, intervals
!> [left, right] are maintained by storing their mid-points and
!> semi-widths in the arrays W and WERR respectively.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix.
!> \endverbatim
!>
!> \param[in] D
!> \verbatim
!>          D is REAL array, dimension (N)
!>          The N diagonal elements of the diagonal matrix D.
!> \endverbatim
!>
!> \param[in] LLD
!> \verbatim
!>          LLD is REAL array, dimension (N-1)
!>          The (N-1) elements L(i)*L(i)*D(i).
!> \endverbatim
!>
!> \param[in] IFIRST
!> \verbatim
!>          IFIRST is INTEGER
!>          The index of the first eigenvalue to be computed.
!> \endverbatim
!>
!> \param[in] ILAST
!> \verbatim
!>          ILAST is INTEGER
!>          The index of the last eigenvalue to be computed.
!> \endverbatim
!>
!> \param[in] RTOL1
!> \verbatim
!>          RTOL1 is REAL
!> \endverbatim
!>
!> \param[in] RTOL2
!> \verbatim
!>          RTOL2 is REAL
!>          Tolerance for the convergence of the bisection intervals.
!>          An interval [LEFT,RIGHT] has converged if
!>          RIGHT-LEFT < MAX( RTOL1*GAP, RTOL2*MAX(|LEFT|,|RIGHT|) )
!>          where GAP is the (estimated) distance to the nearest
!>          eigenvalue.
!> \endverbatim
!>
!> \param[in] OFFSET
!> \verbatim
!>          OFFSET is INTEGER
!>          Offset for the arrays W, WGAP and WERR, i.e., the IFIRST-OFFSET
!>          through ILAST-OFFSET elements of these arrays are to be used.
!> \endverbatim
!>
!> \param[in,out] W
!> \verbatim
!>          W is REAL array, dimension (N)
!>          On input, W( IFIRST-OFFSET ) through W( ILAST-OFFSET ) are
!>          estimates of the eigenvalues of L D L^T indexed IFIRST through
!>          ILAST.
!>          On output, these estimates are refined.
!> \endverbatim
!>
!> \param[in,out] WGAP
!> \verbatim
!>          WGAP is REAL array, dimension (N-1)
!>          On input, the (estimated) gaps between consecutive
!>          eigenvalues of L D L^T, i.e., WGAP(I-OFFSET) is the gap between
!>          eigenvalues I and I+1. Note that if IFIRST = ILAST
!>          then WGAP(IFIRST-OFFSET) must be set to ZERO.
!>          On output, these gaps are refined.
!> \endverbatim
!>
!> \param[in,out] WERR
!> \verbatim
!>          WERR is REAL array, dimension (N)
!>          On input, WERR( IFIRST-OFFSET ) through WERR( ILAST-OFFSET ) are
!>          the errors in the estimates of the corresponding elements in W.
!>          On output, these errors are refined.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is REAL array, dimension (2*N)
!>          Workspace.
!> \endverbatim
!>
!> \param[out] IWORK
!> \verbatim
!>          IWORK is INTEGER array, dimension (2*N)
!>          Workspace.
!> \endverbatim
!>
!> \param[in] PIVMIN
!> \verbatim
!>          PIVMIN is REAL
!>          The minimum pivot in the Sturm sequence.
!> \endverbatim
!>
!> \param[in] SPDIAM
!> \verbatim
!>          SPDIAM is REAL
!>          The spectral diameter of the matrix.
!> \endverbatim
!>
!> \param[in] TWIST
!> \verbatim
!>          TWIST is INTEGER
!>          The twist index for the twisted factorization that is used
!>          for the negcount.
!>          TWIST = N: Compute negcount from L D L^T - LAMBDA I = L+ D+ L+^T
!>          TWIST = 1: Compute negcount from L D L^T - LAMBDA I = U- D- U-^T
!>          TWIST = R: Compute negcount from L D L^T - LAMBDA I = N(r) D(r) N(r)
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          Error flag.
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
!> \ingroup OTHERauxiliary
!
!> \par Contributors:
!  ==================
!>
!> Beresford Parlett, University of California, Berkeley, USA \n
!> Jim Demmel, University of California, Berkeley, USA \n
!> Inderjit Dhillon, University of Texas, Austin, USA \n
!> Osni Marques, LBNL/NERSC, USA \n
!> Christof Voemel, University of California, Berkeley, USA
!
!  =====================================================================
      SUBROUTINE SLARRB(N,D,Lld,Ifirst,Ilast,Rtol1,Rtol2,Offset,W,Wgap, &
     &                  Werr,Work,Iwork,Pivmin,Spdiam,Twist,Info)
      IMPLICIT NONE
!*--SLARRB199
!
!  -- LAPACK auxiliary routine (version 3.7.1) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     June 2017
!
!     .. Scalar Arguments ..
      INTEGER Ifirst , Ilast , Info , N , Offset , Twist
      REAL Pivmin , Rtol1 , Rtol2 , Spdiam
!     ..
!     .. Array Arguments ..
      INTEGER Iwork(*)
      REAL D(*) , Lld(*) , W(*) , Werr(*) , Wgap(*) , Work(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ZERO , TWO , HALF
      PARAMETER (ZERO=0.0E0,TWO=2.0E0,HALF=0.5E0)
      INTEGER maxitr
!     ..
!     .. Local Scalars ..
      INTEGER i , i1 , ii , ip , iter , k , negcnt , next , nint ,      &
     &        olnint , prev , r
      REAL back , cvrgd , gap , left , lgap , mid , mnwdth , rgap ,     &
     &     right , tmp , width
!     ..
!     .. External Functions ..
      INTEGER SLANEG
      EXTERNAL SLANEG
!
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , MAX , MIN
!     ..
!     .. Executable Statements ..
!
      Info = 0
!
!     Quick return if possible
!
      IF ( N<=0 ) RETURN
!
      maxitr = INT((LOG(Spdiam+Pivmin)-LOG(Pivmin))/LOG(TWO)) + 2
      mnwdth = TWO*Pivmin
!
      r = Twist
      IF ( (r<1) .OR. (r>N) ) r = N
!
!     Initialize unconverged intervals in [ WORK(2*I-1), WORK(2*I) ].
!     The Sturm Count, Count( WORK(2*I-1) ) is arranged to be I-1, while
!     Count( WORK(2*I) ) is stored in IWORK( 2*I ). The integer IWORK( 2*I-1 )
!     for an unconverged interval is set to the index of the next unconverged
!     interval, and is -1 or 0 for a converged interval. Thus a linked
!     list of unconverged intervals is set up.
!
      i1 = Ifirst
!     The number of unconverged intervals
      nint = 0
!     The last unconverged interval found
      prev = 0
 
      rgap = Wgap(i1-Offset)
      DO i = i1 , Ilast
         k = 2*i
         ii = i - Offset
         left = W(ii) - Werr(ii)
         right = W(ii) + Werr(ii)
         lgap = rgap
         rgap = Wgap(ii)
         gap = MIN(lgap,rgap)
 
!        Make sure that [LEFT,RIGHT] contains the desired eigenvalue
!        Compute negcount from dstqds facto L+D+L+^T = L D L^T - LEFT
!
!        Do while( NEGCNT(LEFT).GT.I-1 )
!
         back = Werr(ii)
         DO
            negcnt = SLANEG(N,D,Lld,left,Pivmin,r)
            IF ( negcnt>i-1 ) THEN
               left = left - back
               back = TWO*back
               CYCLE
            ENDIF
!
!        Do while( NEGCNT(RIGHT).LT.I )
!        Compute negcount from dstqds facto L+D+L+^T = L D L^T - RIGHT
!
            back = Werr(ii)
            EXIT
         ENDDO
         DO
 
            negcnt = SLANEG(N,D,Lld,right,Pivmin,r)
            IF ( negcnt<i ) THEN
               right = right + back
               back = TWO*back
               CYCLE
            ENDIF
            width = HALF*ABS(left-right)
            tmp = MAX(ABS(left),ABS(right))
            cvrgd = MAX(Rtol1*gap,Rtol2*tmp)
            IF ( width<=cvrgd .OR. width<=mnwdth ) THEN
!           This interval has already converged and does not need refinement.
!           (Note that the gaps might change through refining the
!            eigenvalues, however, they can only get bigger.)
!           Remove it from the list.
               Iwork(k-1) = -1
!           Make sure that I1 always points to the first unconverged interval
               IF ( (i==i1) .AND. (i<Ilast) ) i1 = i + 1
               IF ( (prev>=i1) .AND. (i<=Ilast) ) Iwork(2*prev-1)       &
     &              = i + 1
            ELSE
!           unconverged interval found
               prev = i
               nint = nint + 1
               Iwork(k-1) = i + 1
               Iwork(k) = negcnt
            ENDIF
            Work(k-1) = left
            Work(k) = right
            EXIT
         ENDDO
      ENDDO
 
!
!     Do while( NINT.GT.0 ), i.e. there are still unconverged intervals
!     and while (ITER.LT.MAXITR)
!
      iter = 0
      DO
         prev = i1 - 1
         i = i1
         olnint = nint
 
         DO ip = 1 , olnint
            k = 2*i
            ii = i - Offset
            rgap = Wgap(ii)
            lgap = rgap
            IF ( ii>1 ) lgap = Wgap(ii-1)
            gap = MIN(lgap,rgap)
            next = Iwork(k-1)
            left = Work(k-1)
            right = Work(k)
            mid = HALF*(left+right)
 
!        semiwidth of interval
            width = right - mid
            tmp = MAX(ABS(left),ABS(right))
            cvrgd = MAX(Rtol1*gap,Rtol2*tmp)
            IF ( (width<=cvrgd) .OR. (width<=mnwdth) .OR. (iter==maxitr)&
     &           ) THEN
!           reduce number of unconverged intervals
               nint = nint - 1
!           Mark interval as converged.
               Iwork(k-1) = 0
               IF ( i1==i ) THEN
                  i1 = next
               ELSE
!              Prev holds the last unconverged interval previously examined
                  IF ( prev>=i1 ) Iwork(2*prev-1) = next
               ENDIF
               i = next
               CYCLE
            ENDIF
            prev = i
!
!        Perform one bisection step
!
            negcnt = SLANEG(N,D,Lld,mid,Pivmin,r)
            IF ( negcnt<=i-1 ) THEN
               Work(k-1) = mid
            ELSE
               Work(k) = mid
            ENDIF
            i = next
         ENDDO
         iter = iter + 1
!     do another loop if there are still unconverged intervals
!     However, in the last iteration, all intervals are accepted
!     since this is the best we can do.
         IF ( (nint<=0) .OR. (iter>maxitr) ) THEN
!
!
!     At this point, all the intervals have converged
            DO i = Ifirst , Ilast
               k = 2*i
               ii = i - Offset
!        All intervals marked by '0' have been refined.
               IF ( Iwork(k-1)==0 ) THEN
                  W(ii) = HALF*(Work(k-1)+Work(k))
                  Werr(ii) = Work(k) - W(ii)
               ENDIF
            ENDDO
!
            DO i = Ifirst + 1 , Ilast
               k = 2*i
               ii = i - Offset
               Wgap(ii-1) = MAX(ZERO,W(ii)-Werr(ii)-W(ii-1)-Werr(ii-1))
            ENDDO
            EXIT
         ENDIF
      ENDDO
 
!
!     End of SLARRB
!
      END SUBROUTINE SLARRB
