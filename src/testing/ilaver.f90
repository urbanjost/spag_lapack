!*==ilaver.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b ILAVER returns the LAPACK version.
!*
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!     SUBROUTINE ILAVER( VERS_MAJOR, VERS_MINOR, VERS_PATCH )
!
!     INTEGER VERS_MAJOR, VERS_MINOR, VERS_PATCH
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>  This subroutine returns the LAPACK version.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!>  \param[out] VERS_MAJOR
!>      VERS_MAJOR is INTEGER
!>      return the lapack major version
!>
!>  \param[out] VERS_MINOR
!>      VERS_MINOR is INTEGER
!>      return the lapack minor version from the major version
!>
!>  \param[out] VERS_PATCH
!>      VERS_PATCH is INTEGER
!>      return the lapack patch version from the minor version
!
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \date November 2019
!
!> \ingroup auxOTHERauxiliary
!
!  =====================================================================
      SUBROUTINE ILAVER(Vers_major,Vers_minor,Vers_patch)
      IMPLICIT NONE
!*--ILAVER55
!
!  -- LAPACK computational routine --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!
!  =====================================================================
!
      INTEGER Vers_major , Vers_minor , Vers_patch
!  =====================================================================
      Vers_major = 3
      Vers_minor = 9
      Vers_patch = 0
!  =====================================================================
!
      END SUBROUTINE ILAVER
