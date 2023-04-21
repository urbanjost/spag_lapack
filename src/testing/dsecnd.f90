!*==dsecnd.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b DSECND Using INTERNAL function CPU_TIME.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!      DOUBLE PRECISION FUNCTION DSECND( )
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>  DSECND returns the user time for a process in seconds.
!>  This version gets the time from the INTERNAL function CPU_TIME.
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
!> \ingroup auxOTHERauxiliary
!
!  =====================================================================
      DOUBLE PRECISION FUNCTION DSECND()
      IMPLICIT NONE
!*--DSECND39
!
!  -- LAPACK auxiliary routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
! =====================================================================
!
!     .. Local Scalars ..
!
      REAL t
!
! .. Intrinsic Functions ..
!
      INTRINSIC CPU_TIME
!
! .. Executable Statements .. *
!
      CALL CPU_TIME(t)
      DSECND = t
!
!     End of DSECND
!
      END FUNCTION DSECND
