!*==sget06.f90  processed by SPAG 7.51RB at 17:34 on  4 Mar 2022
!> \brief \b sget06
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       REAL             FUNCTION sget06( RCOND, RCONDC )
!
!       .. Scalar Arguments ..
!       REAL               RCOND, RCONDC
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> sget06 computes a test ratio to compare two values for RCOND.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] RCOND
!> \verbatim
!>          RCOND is REAL
!>          The estimate of the reciprocal of the condition number of A,
!>          as computed by SGECON.
!> \endverbatim
!>
!> \param[in] RCONDC
!> \verbatim
!>          RCONDC is REAL
!>          The reciprocal of the condition number of A, computed as
!>          ( 1/norm(A) ) / norm(inv(A)).
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
!> \ingroup single_lin
!
!  =====================================================================
      REAL FUNCTION sget06(Rcond,Rcondc)
      IMPLICIT NONE
!*--SGET0659
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      REAL Rcond , Rcondc
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ZERO , ONE
      PARAMETER (ZERO=0.0E+0,ONE=1.0E+0)
!     ..
!     .. Local Scalars ..
      REAL eps , rat
!     ..
!     .. External Functions ..
      REAL SLAMCH
EXTERNAL    SLAMCH
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX , MIN
!     ..
!     .. Executable Statements ..
!
      eps = SLAMCH('Epsilon')
      IF ( Rcond>ZERO ) THEN
         IF ( Rcondc>ZERO ) THEN
            rat = MAX(Rcond,Rcondc)/MIN(Rcond,Rcondc) - (ONE-eps)
         ELSE
            rat = Rcond/eps
         ENDIF
      ELSEIF ( Rcondc>ZERO ) THEN
         rat = Rcondc/eps
      ELSE
         rat = ZERO
      ENDIF
      sget06 = rat
!
!     End of sget06
!
      END FUNCTION sget06
