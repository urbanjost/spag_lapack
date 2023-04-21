!*==dlarrk.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b DLARRK computes one eigenvalue of a symmetric tridiagonal matrix T to suitable accuracy.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DLARRK + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlarrk.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlarrk.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlarrk.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLARRK( N, IW, GL, GU,
!                           D, E2, PIVMIN, RELTOL, W, WERR, INFO)
!
!       .. Scalar Arguments ..
!       INTEGER   INFO, IW, N
!       DOUBLE PRECISION    PIVMIN, RELTOL, GL, GU, W, WERR
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   D( * ), E2( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLARRK computes one eigenvalue of a symmetric tridiagonal
!> matrix T to suitable accuracy. This is an auxiliary code to be
!> called from DSTEMR.
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
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the tridiagonal matrix T.  N >= 0.
!> \endverbatim
!>
!> \param[in] IW
!> \verbatim
!>          IW is INTEGER
!>          The index of the eigenvalues to be returned.
!> \endverbatim
!>
!> \param[in] GL
!> \verbatim
!>          GL is DOUBLE PRECISION
!> \endverbatim
!>
!> \param[in] GU
!> \verbatim
!>          GU is DOUBLE PRECISION
!>          An upper and a lower bound on the eigenvalue.
!> \endverbatim
!>
!> \param[in] D
!> \verbatim
!>          D is DOUBLE PRECISION array, dimension (N)
!>          The n diagonal elements of the tridiagonal matrix T.
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
!> \param[out] W
!> \verbatim
!>          W is DOUBLE PRECISION
!> \endverbatim
!>
!> \param[out] WERR
!> \verbatim
!>          WERR is DOUBLE PRECISION
!>          The error bound on the corresponding eigenvalue approximation
!>          in W.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:       Eigenvalue converged
!>          = -1:      Eigenvalue did NOT converge
!> \endverbatim
!
!> \par Internal Parameters:
!  =========================
!>
!> \verbatim
!>  FUDGE   DOUBLE PRECISION, default = 2
!>          A "fudge factor" to widen the Gershgorin intervals.
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
!  =====================================================================
      SUBROUTINE DLARRK(N,Iw,Gl,Gu,D,E2,Pivmin,Reltol,W,Werr,Info)
      IMPLICIT NONE
!*--DLARRK148
!
!  -- LAPACK auxiliary routine (version 3.7.1) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     June 2017
!
!     .. Scalar Arguments ..
      INTEGER Info , Iw , N
      DOUBLE PRECISION Pivmin , Reltol , Gl , Gu , W , Werr
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION D(*) , E2(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION FUDGE , HALF , TWO , ZERO
      PARAMETER (HALF=0.5D0,TWO=2.0D0,FUDGE=TWO,ZERO=0.0D0)
!     ..
!     .. Local Scalars ..
      INTEGER i , it , itmax , negcnt
      DOUBLE PRECISION atoli , eps , left , mid , right , rtoli , tmp1 ,&
     &                 tmp2 , tnorm
!     ..
!     .. External Functions ..
      DOUBLE PRECISION DLAMCH
      EXTERNAL DLAMCH
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , INT , LOG , MAX
!     ..
!     .. Executable Statements ..
!
!     Quick return if possible
!
      IF ( N<=0 ) THEN
         Info = 0
         RETURN
      ENDIF
!
!     Get machine constants
      eps = DLAMCH('P')
 
      tnorm = MAX(ABS(Gl),ABS(Gu))
      rtoli = Reltol
      atoli = FUDGE*TWO*Pivmin
 
      itmax = INT((LOG(tnorm+Pivmin)-LOG(Pivmin))/LOG(TWO)) + 2
 
      Info = -1
 
      left = Gl - FUDGE*tnorm*eps*N - FUDGE*TWO*Pivmin
      right = Gu + FUDGE*tnorm*eps*N + FUDGE*TWO*Pivmin
      it = 0
      DO
 
!
!     Check if interval converged or maximum number of iterations reached
!
         tmp1 = ABS(right-left)
         tmp2 = MAX(ABS(right),ABS(left))
         IF ( tmp1<MAX(atoli,Pivmin,rtoli*tmp2) ) THEN
            Info = 0
            EXIT
         ENDIF
         IF ( it>itmax ) EXIT
 
!
!     Count number of negative pivots for mid-point
!
         it = it + 1
         mid = HALF*(left+right)
         negcnt = 0
         tmp1 = D(1) - mid
         IF ( ABS(tmp1)<Pivmin ) tmp1 = -Pivmin
         IF ( tmp1<=ZERO ) negcnt = negcnt + 1
!
         DO i = 2 , N
            tmp1 = D(i) - E2(i-1)/tmp1 - mid
            IF ( ABS(tmp1)<Pivmin ) tmp1 = -Pivmin
            IF ( tmp1<=ZERO ) negcnt = negcnt + 1
         ENDDO
 
         IF ( negcnt>=Iw ) THEN
            right = mid
         ELSE
            left = mid
         ENDIF
      ENDDO
 
!
!     Converged or maximum number of iterations reached
!
      W = HALF*(left+right)
      Werr = HALF*ABS(right-left)
 
!
!     End of DLARRK
!
      END SUBROUTINE DLARRK
