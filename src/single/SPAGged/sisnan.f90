!*==sisnan.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b SISNAN tests input for NaN.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SISNAN + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sisnan.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sisnan.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sisnan.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       LOGICAL FUNCTION SISNAN( SIN )
!
!       .. Scalar Arguments ..
!       REAL, INTENT(IN) :: SIN
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SISNAN returns .TRUE. if its argument is NaN, and .FALSE.
!> otherwise.  To be replaced by the Fortran 2003 intrinsic in the
!> future.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] SIN
!> \verbatim
!>          SIN is REAL
!>          Input to test for NaN.
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
      FUNCTION SISNAN(Sin)
      USE S_SLAISNAN
      IMPLICIT NONE
!*--SISNAN64
      LOGICAL :: SISNAN
!
! Dummy argument declarations rewritten by SPAG
!
      REAL , INTENT(IN) :: Sin
!
! End of declarations rewritten by SPAG
!
!     ..
!
!  =====================================================================
!
!  .. External Functions ..
!  ..
!  .. Executable Statements ..
      SISNAN = SLAISNAN(Sin,Sin)
      END FUNCTION SISNAN
