!*==dlarrj.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b DLARRJ performs refinement of the initial estimates of the eigenvalues of the matrix T.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DLARRJ + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlarrj.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlarrj.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlarrj.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLARRJ( N, D, E2, IFIRST, ILAST,
!                          RTOL, OFFSET, W, WERR, WORK, IWORK,
!                          PIVMIN, SPDIAM, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            IFIRST, ILAST, INFO, N, OFFSET
!       DOUBLE PRECISION   PIVMIN, RTOL, SPDIAM
!       ..
!       .. Array Arguments ..
!       INTEGER            IWORK( * )
!       DOUBLE PRECISION   D( * ), E2( * ), W( * ),
!      $                   WERR( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> Given the initial eigenvalue approximations of T, DLARRJ
!> does  bisection to refine the eigenvalues of T,
!> W( IFIRST-OFFSET ) through W( ILAST-OFFSET ), to more accuracy. Initial
!> guesses for these eigenvalues are input in W, the corresponding estimate
!> of the error in these guesses in WERR. During bisection, intervals
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
!>          D is DOUBLE PRECISION array, dimension (N)
!>          The N diagonal elements of T.
!> \endverbatim
!>
!> \param[in] E2
!> \verbatim
!>          E2 is DOUBLE PRECISION array, dimension (N-1)
!>          The Squares of the (N-1) subdiagonal elements of T.
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
!> \param[in] RTOL
!> \verbatim
!>          RTOL is DOUBLE PRECISION
!>          Tolerance for the convergence of the bisection intervals.
!>          An interval [LEFT,RIGHT] has converged if
!>          RIGHT-LEFT < RTOL*MAX(|LEFT|,|RIGHT|).
!> \endverbatim
!>
!> \param[in] OFFSET
!> \verbatim
!>          OFFSET is INTEGER
!>          Offset for the arrays W and WERR, i.e., the IFIRST-OFFSET
!>          through ILAST-OFFSET elements of these arrays are to be used.
!> \endverbatim
!>
!> \param[in,out] W
!> \verbatim
!>          W is DOUBLE PRECISION array, dimension (N)
!>          On input, W( IFIRST-OFFSET ) through W( ILAST-OFFSET ) are
!>          estimates of the eigenvalues of L D L^T indexed IFIRST through
!>          ILAST.
!>          On output, these estimates are refined.
!> \endverbatim
!>
!> \param[in,out] WERR
!> \verbatim
!>          WERR is DOUBLE PRECISION array, dimension (N)
!>          On input, WERR( IFIRST-OFFSET ) through WERR( ILAST-OFFSET ) are
!>          the errors in the estimates of the corresponding elements in W.
!>          On output, these errors are refined.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (2*N)
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
!>          PIVMIN is DOUBLE PRECISION
!>          The minimum pivot in the Sturm sequence for T.
!> \endverbatim
!>
!> \param[in] SPDIAM
!> \verbatim
!>          SPDIAM is DOUBLE PRECISION
!>          The spectral diameter of T.
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
      SUBROUTINE DLARRJ(N,D,E2,Ifirst,Ilast,Rtol,Offset,W,Werr,Work,    &
     &                  Iwork,Pivmin,Spdiam,Info)
      USE F77KINDS                        
      IMPLICIT NONE
!*--DLARRJ172
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D0 , ONE = 1.0D0 ,        &
     &                              TWO = 2.0D0 , HALF = 0.5D0
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: D
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: E2
      INTEGER , INTENT(IN) :: Ifirst
      INTEGER , INTENT(IN) :: Ilast
      REAL(R8KIND) , INTENT(IN) :: Rtol
      INTEGER , INTENT(IN) :: Offset
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: W
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Werr
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iwork
      REAL(R8KIND) , INTENT(IN) :: Pivmin
      REAL(R8KIND) , INTENT(IN) :: Spdiam
      INTEGER , INTENT(OUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      INTEGER :: cnt , i , i1 , i2 , ii , iter , j , k , maxitr , next ,&
     &           nint , olnint , p , prev , savi1
      REAL(R8KIND) :: dplus , fac , left , mid , right , s , tmp , width
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
!
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
      maxitr = INT((LOG(Spdiam+Pivmin)-LOG(Pivmin))/LOG(TWO)) + 2
!
!     Initialize unconverged intervals in [ WORK(2*I-1), WORK(2*I) ].
!     The Sturm Count, Count( WORK(2*I-1) ) is arranged to be I-1, while
!     Count( WORK(2*I) ) is stored in IWORK( 2*I ). The integer IWORK( 2*I-1 )
!     for an unconverged interval is set to the index of the next unconverged
!     interval, and is -1 or 0 for a converged interval. Thus a linked
!     list of unconverged intervals is set up.
!
 
      i1 = Ifirst
      i2 = Ilast
!     The number of unconverged intervals
      nint = 0
!     The last unconverged interval found
      prev = 0
      DO i = i1 , i2
         k = 2*i
         ii = i - Offset
         left = W(ii) - Werr(ii)
         mid = W(ii)
         right = W(ii) + Werr(ii)
         width = right - mid
         tmp = MAX(ABS(left),ABS(right))
 
!        The following test prevents the test of converged intervals
         IF ( width<Rtol*tmp ) THEN
!           This interval has already converged and does not need refinement.
!           (Note that the gaps might change through refining the
!            eigenvalues, however, they can only get bigger.)
!           Remove it from the list.
            Iwork(k-1) = -1
!           Make sure that I1 always points to the first unconverged interval
            IF ( (i==i1) .AND. (i<i2) ) i1 = i + 1
            IF ( (prev>=i1) .AND. (i<=i2) ) Iwork(2*prev-1) = i + 1
         ELSE
!           unconverged interval found
            prev = i
!           Make sure that [LEFT,RIGHT] contains the desired eigenvalue
!
!           Do while( CNT(LEFT).GT.I-1 )
!
            fac = ONE
            DO
               cnt = 0
               s = left
               dplus = D(1) - s
               IF ( dplus<ZERO ) cnt = cnt + 1
               DO j = 2 , N
                  dplus = D(j) - s - E2(j-1)/dplus
                  IF ( dplus<ZERO ) cnt = cnt + 1
               ENDDO
               IF ( cnt>i-1 ) THEN
                  left = left - Werr(ii)*fac
                  fac = TWO*fac
                  CYCLE
               ENDIF
!
!           Do while( CNT(RIGHT).LT.I )
!
               fac = ONE
               EXIT
            ENDDO
            DO
               cnt = 0
               s = right
               dplus = D(1) - s
               IF ( dplus<ZERO ) cnt = cnt + 1
               DO j = 2 , N
                  dplus = D(j) - s - E2(j-1)/dplus
                  IF ( dplus<ZERO ) cnt = cnt + 1
               ENDDO
               IF ( cnt<i ) THEN
                  right = right + Werr(ii)*fac
                  fac = TWO*fac
                  CYCLE
               ENDIF
               nint = nint + 1
               Iwork(k-1) = i + 1
               Iwork(k) = cnt
               EXIT
            ENDDO
         ENDIF
         Work(k-1) = left
         Work(k) = right
      ENDDO
 
 
      savi1 = i1
!
!     Do while( NINT.GT.0 ), i.e. there are still unconverged intervals
!     and while (ITER.LT.MAXITR)
!
      iter = 0
      DO
         prev = i1 - 1
         i = i1
         olnint = nint
 
         DO p = 1 , olnint
            k = 2*i
            ii = i - Offset
            next = Iwork(k-1)
            left = Work(k-1)
            right = Work(k)
            mid = HALF*(left+right)
 
!        semiwidth of interval
            width = right - mid
            tmp = MAX(ABS(left),ABS(right))
 
            IF ( (width<Rtol*tmp) .OR. (iter==maxitr) ) THEN
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
            cnt = 0
            s = mid
            dplus = D(1) - s
            IF ( dplus<ZERO ) cnt = cnt + 1
            DO j = 2 , N
               dplus = D(j) - s - E2(j-1)/dplus
               IF ( dplus<ZERO ) cnt = cnt + 1
            ENDDO
            IF ( cnt<=i-1 ) THEN
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
            DO i = savi1 , Ilast
               k = 2*i
               ii = i - Offset
!        All intervals marked by '0' have been refined.
               IF ( Iwork(k-1)==0 ) THEN
                  W(ii) = HALF*(Work(k-1)+Work(k))
                  Werr(ii) = Work(k) - W(ii)
               ENDIF
            ENDDO
            EXIT
         ENDIF
      ENDDO
!
 
!
!     End of DLARRJ
!
      END SUBROUTINE DLARRJ
