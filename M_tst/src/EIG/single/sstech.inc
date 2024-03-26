!*==sstech.f90  processed by SPAG 7.51RB at 17:37 on  4 Mar 2022
!> \brief \b sstech
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE sstech( N, A, B, EIG, TOL, WORK, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, N
!       REAL               TOL
!       ..
!       .. Array Arguments ..
!       REAL               A( * ), B( * ), EIG( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    Let T be the tridiagonal matrix with diagonal entries A(1) ,...,
!>    A(N) and offdiagonal entries B(1) ,..., B(N-1)).  sstech checks to
!>    see if EIG(1) ,..., EIG(N) are indeed accurate eigenvalues of T.
!>    It does this by expanding each EIG(I) into an interval
!>    [SVD(I) - EPS, SVD(I) + EPS], merging overlapping intervals if
!>    any, and using Sturm sequences to count and verify whether each
!>    resulting interval has the correct number of eigenvalues (using
!>    sstect).  Here EPS = TOL*MACHEPS*MAXEIG, where MACHEPS is the
!>    machine precision and MAXEIG is the absolute value of the largest
!>    eigenvalue. If each interval contains the correct number of
!>    eigenvalues, INFO = 0 is returned, otherwise INFO is the index of
!>    the first eigenvalue in the first bad interval.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The dimension of the tridiagonal matrix T.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is REAL array, dimension (N)
!>          The diagonal entries of the tridiagonal matrix T.
!> \endverbatim
!>
!> \param[in] B
!> \verbatim
!>          B is REAL array, dimension (N-1)
!>          The offdiagonal entries of the tridiagonal matrix T.
!> \endverbatim
!>
!> \param[in] EIG
!> \verbatim
!>          EIG is REAL array, dimension (N)
!>          The purported eigenvalues to be checked.
!> \endverbatim
!>
!> \param[in] TOL
!> \verbatim
!>          TOL is REAL
!>          Error tolerance for checking, a multiple of the
!>          machine precision.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is REAL array, dimension (N)
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          0  if the eigenvalues are all correct (to within
!>             1 +- TOL*MACHEPS*MAXEIG)
!>          >0 if the interval containing the INFO-th eigenvalue
!>             contains the incorrect number of eigenvalues.
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
      SUBROUTINE sstech(N,A,B,Eig,Tol,Work,Info)
      IMPLICIT NONE
!*--SSTECH105
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
      REAL A(*) , B(*) , Eig(*) , Work(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ZERO
      PARAMETER (ZERO=0.0E0)
!     ..
!     .. Local Scalars ..
      INTEGER bpnt , count , i , isub , j , numl , numu , tpnt
      REAL emin , eps , lower , mx , tuppr , unflep , upper
!     ..
!     .. External Functions ..
      REAL SLAMCH
      EXTERNAL SLAMCH
!     ..
!     .. External Subroutines ..
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , MAX
!     ..
!     .. Executable Statements ..
!
!     Check input parameters
!
      Info = 0
      IF ( N==0 ) RETURN
      IF ( N<0 ) THEN
         Info = -1
         RETURN
      ENDIF
      IF ( Tol<ZERO ) THEN
         Info = -5
         RETURN
      ENDIF
!
!     Get machine constants
!
      eps = SLAMCH('Epsilon')*SLAMCH('Base')
      unflep = SLAMCH('Safe minimum')/eps
      eps = Tol*eps
!
!     Compute maximum absolute eigenvalue, error tolerance
!
      mx = ABS(Eig(1))
      DO i = 2 , N
         mx = MAX(mx,ABS(Eig(i)))
      ENDDO
      eps = MAX(eps*mx,unflep)
!
!     Sort eigenvalues from EIG into WORK
!
      DO i = 1 , N
         Work(i) = Eig(i)
      ENDDO
      DO i = 1 , N - 1
         isub = 1
         emin = Work(1)
         DO j = 2 , N + 1 - i
            IF ( Work(j)<emin ) THEN
               isub = j
               emin = Work(j)
            ENDIF
         ENDDO
         IF ( isub/=N+1-i ) THEN
            Work(isub) = Work(N+1-i)
            Work(N+1-i) = emin
         ENDIF
      ENDDO
!
!     TPNT points to singular value at right endpoint of interval
!     BPNT points to singular value at left  endpoint of interval
!
      tpnt = 1
      bpnt = 1
!
!     Begin loop over all intervals
!
 100  upper = Work(tpnt) + eps
      lower = Work(bpnt) - eps
!
!     Begin loop merging overlapping intervals
!
      DO WHILE ( bpnt/=N )
         tuppr = Work(bpnt+1) + eps
         IF ( tuppr<lower ) EXIT
!
!     Merge
!
         bpnt = bpnt + 1
         lower = Work(bpnt) - eps
      ENDDO
!
!     Count singular values in interval [ LOWER, UPPER ]
!
      CALL sstect(N,A,B,lower,numl)
      CALL sstect(N,A,B,upper,numu)
      count = numu - numl
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
!     End of sstech
!
99999 END SUBROUTINE sstech
