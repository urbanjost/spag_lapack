!*==second.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b SECOND  Using ETIME_
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
!>  SECOND returns the user time for a process in seconds.
!>  This version gets the time from the system function ETIME_.
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
!*--SECOND39
!
!  -- LAPACK auxiliary routine (version 3.7.0) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     February 2007
! =====================================================================
!
!     .. Local Scalars ..
      REAL t1
!     ..
!     .. Local Arrays ..
      REAL tarray(2)
!     ..
!     .. External Functions ..
      REAL ETIME_
      EXTERNAL ETIME_
!     ..
!     .. Executable Statements ..
!
      t1 = ETIME_(tarray)
      SECOND = tarray(1)
!
!     End of SECOND
!
      END FUNCTION SECOND
