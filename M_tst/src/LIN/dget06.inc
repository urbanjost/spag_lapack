!*==dget06.f90  processed by SPAG 7.51RB at 17:34 on  4 Mar 2022
!> \brief \b dget06
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       DOUBLE PRECISION FUNCTION dget06( RCOND, RCONDC )
!
!       .. Scalar Arguments ..
!       DOUBLE PRECISION   RCOND, RCONDC
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> dget06 computes a test ratio to compare two values for RCOND.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] RCOND
!> \verbatim
!>          RCOND is DOUBLE PRECISION
!>          The estimate of the reciprocal of the condition number of A,
!>          as computed by DGECON.
!> \endverbatim
!>
!> \param[in] RCONDC
!> \verbatim
!>          RCONDC is DOUBLE PRECISION
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
!> \ingroup double_lin
!
!  =====================================================================
      DOUBLE PRECISION FUNCTION dget06(Rcond,Rcondc)
      IMPLICIT NONE
!*--DGET0659
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      DOUBLE PRECISION Rcond , Rcondc
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO , ONE
      PARAMETER (ZERO=0.0D+0,ONE=1.0D+0)
!     ..
!     .. Local Scalars ..
      DOUBLE PRECISION eps , rat
!     ..
!     .. External Functions ..
      DOUBLE PRECISION DLAMCH
EXTERNAL    DLAMCH
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX , MIN
!     ..
!     .. Executable Statements ..
!
      eps = DLAMCH('Epsilon')
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
      dget06 = rat
!
!     End of dget06
!
      END FUNCTION dget06
