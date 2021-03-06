!*==slapy3.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b SLAPY3 returns sqrt(x2+y2+z2).
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SLAPY3 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slapy3.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slapy3.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slapy3.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       REAL             FUNCTION SLAPY3( X, Y, Z )
!
!       .. Scalar Arguments ..
!       REAL               X, Y, Z
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SLAPY3 returns sqrt(x**2+y**2+z**2), taking care not to cause
!> unnecessary overflow.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] X
!> \verbatim
!>          X is REAL
!> \endverbatim
!>
!> \param[in] Y
!> \verbatim
!>          Y is REAL
!> \endverbatim
!>
!> \param[in] Z
!> \verbatim
!>          Z is REAL
!>          X, Y and Z specify the values x, y and z.
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
      REAL FUNCTION SLAPY3(X,Y,Z)
      IMPLICIT NONE
!*--SLAPY372
!
!  -- LAPACK auxiliary routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      REAL X , Y , Z
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ZERO
      PARAMETER (ZERO=0.0E0)
!     ..
!     .. Local Scalars ..
      REAL w , xabs , yabs , zabs
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , MAX , SQRT
!     ..
!     .. Executable Statements ..
!
      xabs = ABS(X)
      yabs = ABS(Y)
      zabs = ABS(Z)
      w = MAX(xabs,yabs,zabs)
      IF ( w==ZERO ) THEN
!     W can be zero for max(0,nan,0)
!     adding all three entries together will make sure
!     NaN will not disappear.
         SLAPY3 = xabs + yabs + zabs
      ELSE
         SLAPY3 = w*SQRT((xabs/w)**2+(yabs/w)**2+(zabs/w)**2)
      ENDIF
!
!     End of SLAPY3
!
      END FUNCTION SLAPY3
