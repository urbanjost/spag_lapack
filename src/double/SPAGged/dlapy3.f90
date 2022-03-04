!*==dlapy3.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b DLAPY3 returns sqrt(x2+y2+z2).
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DLAPY3 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlapy3.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlapy3.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlapy3.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       DOUBLE PRECISION FUNCTION DLAPY3( X, Y, Z )
!
!       .. Scalar Arguments ..
!       DOUBLE PRECISION   X, Y, Z
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLAPY3 returns sqrt(x**2+y**2+z**2), taking care not to cause
!> unnecessary overflow.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] X
!> \verbatim
!>          X is DOUBLE PRECISION
!> \endverbatim
!>
!> \param[in] Y
!> \verbatim
!>          Y is DOUBLE PRECISION
!> \endverbatim
!>
!> \param[in] Z
!> \verbatim
!>          Z is DOUBLE PRECISION
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
      FUNCTION DLAPY3(X,Y,Z)
      USE F77KINDS                        
      IMPLICIT NONE
!*--DLAPY373
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D0
      REAL(R8KIND) :: DLAPY3
!
! Dummy argument declarations rewritten by SPAG
!
      REAL(R8KIND) , INTENT(IN) :: X
      REAL(R8KIND) , INTENT(IN) :: Y
      REAL(R8KIND) , INTENT(IN) :: Z
!
! Local variable declarations rewritten by SPAG
!
      REAL(R8KIND) :: w , xabs , yabs , zabs
!
! End of declarations rewritten by SPAG
!
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
!     ..
!     .. Local Scalars ..
!     ..
!     .. Intrinsic Functions ..
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
         DLAPY3 = xabs + yabs + zabs
      ELSE
         DLAPY3 = w*SQRT((xabs/w)**2+(yabs/w)**2+(zabs/w)**2)
      ENDIF
!
!     End of DLAPY3
!
      END FUNCTION DLAPY3
