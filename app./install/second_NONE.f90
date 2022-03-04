!*==second.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b SECOND returns nothing
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!      REAL FUNCTION SECOND( )
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>  SECOND returns nothing instead of returning the user time for a process in seconds.
!>  If you are using that routine, it means that neither EXTERNAL ETIME,
!>  EXTERNAL ETIME_, INTERNAL ETIME, INTERNAL CPU_TIME is available  on
!>  your machine.
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
      REAL FUNCTION SECOND()
      IMPLICIT NONE
!*--SECOND41
!
!  -- LAPACK auxiliary routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
! =====================================================================
!
      SECOND = 0.0E+0
!
!     End of SECOND
!
      END FUNCTION SECOND
