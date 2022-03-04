!*==slabad.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b SLABAD
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SLABAD + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slabad.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slabad.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slabad.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SLABAD( SMALL, LARGE )
!
!       .. Scalar Arguments ..
!       REAL               LARGE, SMALL
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SLABAD takes as input the values computed by SLAMCH for underflow and
!> overflow, and returns the square root of each of these values if the
!> log of LARGE is sufficiently large.  This subroutine is intended to
!> identify machines with a large exponent range, such as the Crays, and
!> redefine the underflow and overflow limits to be the square roots of
!> the values computed by SLAMCH.  This subroutine is needed because
!> SLAMCH does not compensate for poor arithmetic in the upper half of
!> the exponent range, as is found on a Cray.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in,out] SMALL
!> \verbatim
!>          SMALL is REAL
!>          On entry, the underflow threshold as computed by SLAMCH.
!>          On exit, if LOG10(LARGE) is sufficiently large, the square
!>          root of SMALL, otherwise unchanged.
!> \endverbatim
!>
!> \param[in,out] LARGE
!> \verbatim
!>          LARGE is REAL
!>          On entry, the overflow threshold as computed by SLAMCH.
!>          On exit, if LOG10(LARGE) is sufficiently large, the square
!>          root of LARGE, otherwise unchanged.
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
!> \ingroup OTHERauxiliary
!
!  =====================================================================
      SUBROUTINE SLABAD(Small,Large)
      IMPLICIT NONE
!*--SLABAD78
!
! Dummy argument declarations rewritten by SPAG
!
      REAL , INTENT(INOUT) :: Small
      REAL , INTENT(INOUT) :: Large
!
! End of declarations rewritten by SPAG
!
!     ..
!
!  =====================================================================
!
!     .. Intrinsic Functions ..
!     ..
!     .. Executable Statements ..
!
!     If it looks like we're on a Cray, take the square root of
!     SMALL and LARGE to avoid overflow and underflow problems.
!
      IF ( LOG10(Large)>2000. ) THEN
         Small = SQRT(Small)
         Large = SQRT(Large)
      ENDIF
!
!
!     End of SLABAD
!
      END SUBROUTINE SLABAD
