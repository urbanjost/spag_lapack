!*==ssvdch.f90  processed by SPAG 7.51RB at 17:37 on  4 Mar 2022
!> \brief \b SSVDCH
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE SSVDCH( N, S, E, SVD, TOL, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, N
!       REAL               TOL
!       ..
!       .. Array Arguments ..
!       REAL               E( * ), S( * ), SVD( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SSVDCH checks to see if SVD(1) ,..., SVD(N) are accurate singular
!> values of the bidiagonal matrix B with diagonal entries
!> S(1) ,..., S(N) and superdiagonal entries E(1) ,..., E(N-1)).
!> It does this by expanding each SVD(I) into an interval
!> [SVD(I) * (1-EPS) , SVD(I) * (1+EPS)], merging overlapping intervals
!> if any, and using Sturm sequences to count and verify whether each
!> resulting interval has the correct number of singular values (using
!> SSVDCT). Here EPS=TOL*MAX(N/10,1)*MACHEP, where MACHEP is the
!> machine precision. The routine assumes the singular values are sorted
!> with SVD(1) the largest and SVD(N) smallest.  If each interval
!> contains the correct number of singular values, INFO = 0 is returned,
!> otherwise INFO is the index of the first singular value in the first
!> bad interval.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The dimension of the bidiagonal matrix B.
!> \endverbatim
!>
!> \param[in] S
!> \verbatim
!>          S is REAL array, dimension (N)
!>          The diagonal entries of the bidiagonal matrix B.
!> \endverbatim
!>
!> \param[in] E
!> \verbatim
!>          E is REAL array, dimension (N-1)
!>          The superdiagonal entries of the bidiagonal matrix B.
!> \endverbatim
!>
!> \param[in] SVD
!> \verbatim
!>          SVD is REAL array, dimension (N)
!>          The computed singular values to be checked.
!> \endverbatim
!>
!> \param[in] TOL
!> \verbatim
!>          TOL is REAL
!>          Error tolerance for checking, a multiplier of the
!>          machine precision.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          =0 if the singular values are all correct (to within
!>             1 +- TOL*MACHEPS)
!>          >0 if the interval containing the INFO-th singular value
!>             contains the incorrect number of singular values.
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
!> \ingroup single_eig
!
!  =====================================================================
      SUBROUTINE SSVDCH(N,S,E,Svd,Tol,Info)
      IMPLICIT NONE
!*--SSVDCH101
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Info , N
      REAL Tol
!     ..
!     .. Array Arguments ..
      REAL E(*) , S(*) , Svd(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ONE
      PARAMETER (ONE=1.0E0)
      REAL ZERO
      PARAMETER (ZERO=0.0E0)
!     ..
!     .. Local Scalars ..
      INTEGER bpnt , count , numl , numu , tpnt
      REAL eps , lower , ovfl , tuppr , unfl , unflep , upper
!     ..
!     .. External Functions ..
      REAL SLAMCH
      EXTERNAL SLAMCH
!     ..
!     .. External Subroutines ..
      EXTERNAL SSVDCT
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX , SQRT
!     ..
!     .. Executable Statements ..
!
!     Get machine constants
!
      Info = 0
      IF ( N<=0 ) RETURN
      unfl = SLAMCH('Safe minimum')
      ovfl = SLAMCH('Overflow')
      eps = SLAMCH('Epsilon')*SLAMCH('Base')
!
!     UNFLEP is chosen so that when an eigenvalue is multiplied by the
!     scale factor sqrt(OVFL)*sqrt(sqrt(UNFL))/MX in SSVDCT, it exceeds
!     sqrt(UNFL), which is the lower limit for SSVDCT.
!
      unflep = (SQRT(SQRT(unfl))/SQRT(ovfl))*Svd(1) + unfl/eps
!
!     The value of EPS works best when TOL .GE. 10.
!
      eps = Tol*MAX(N/10,1)*eps
!
!     TPNT points to singular value at right endpoint of interval
!     BPNT points to singular value at left  endpoint of interval
!
      tpnt = 1
      bpnt = 1
!
!     Begin loop over all intervals
!
 100  upper = (ONE+eps)*Svd(tpnt) + unflep
      lower = (ONE-eps)*Svd(bpnt) - unflep
      IF ( lower<=unflep ) lower = -upper
!
!     Begin loop merging overlapping intervals
!
      DO WHILE ( bpnt/=N )
         tuppr = (ONE+eps)*Svd(bpnt+1) + unflep
         IF ( tuppr<lower ) EXIT
!
!     Merge
!
         bpnt = bpnt + 1
         lower = (ONE-eps)*Svd(bpnt) - unflep
         IF ( lower<=unflep ) lower = -upper
      ENDDO
!
!     Count singular values in interval [ LOWER, UPPER ]
!
      CALL SSVDCT(N,S,E,lower,numl)
      CALL SSVDCT(N,S,E,upper,numu)
      count = numu - numl
      IF ( lower<ZERO ) count = count/2
      IF ( count/=bpnt-tpnt+1 ) THEN
!
!        Wrong number of singular values in interval
!
         Info = tpnt
         GOTO 99999
      ENDIF
      tpnt = bpnt + 1
      bpnt = tpnt
      IF ( tpnt<=N ) GOTO 100
!
!     End of SSVDCH
!
99999 END SUBROUTINE SSVDCH
