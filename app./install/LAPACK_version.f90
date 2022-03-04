!*==lapack_version.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b LAPACK_VERSION
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!      PROGRAM LAPACK_VERSION
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \date November 2017
!
!> \ingroup auxOTHERauxiliary
!
!  =====================================================================
      PROGRAM LAPACK_VERSION
      IMPLICIT NONE
!*--LAPACK_VERSION29
!
!  -- LAPACK auxiliary routine (version 3.8.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2017
!
      INTEGER major , minor , patch
!     ..
!     .. External Subroutines ..
      EXTERNAL ILAVER
!
      CALL ILAVER(major,minor,patch)
      WRITE (*,*) "LAPACK " , major , "." , minor , "." , patch
!
      END PROGRAM LAPACK_VERSION
