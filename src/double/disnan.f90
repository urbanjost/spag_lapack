!*==disnan.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b DISNAN tests input for NaN.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DISNAN + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/disnan.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/disnan.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/disnan.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       LOGICAL FUNCTION DISNAN( DIN )
!
!       .. Scalar Arguments ..
!       DOUBLE PRECISION, INTENT(IN) :: DIN
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DISNAN returns .TRUE. if its argument is NaN, and .FALSE.
!> otherwise.  To be replaced by the Fortran 2003 intrinsic in the
!> future.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] DIN
!> \verbatim
!>          DIN is DOUBLE PRECISION
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
      LOGICAL FUNCTION DISNAN(Din)
      IMPLICIT NONE
!*--DISNAN63
!
!  -- LAPACK auxiliary routine (version 3.7.1) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     June 2017
!
!     .. Scalar Arguments ..
      DOUBLE PRECISION , INTENT(IN)  ::  Din
!     ..
!
!  =====================================================================
!
!  .. External Functions ..
      LOGICAL DLAISNAN
      EXTERNAL DLAISNAN
!  ..
!  .. Executable Statements ..
      DISNAN = DLAISNAN(Din,Din)
      END FUNCTION DISNAN
