!*==slapy2.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b SLAPY2 returns sqrt(x2+y2).
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SLAPY2 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slapy2.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slapy2.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slapy2.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       REAL             FUNCTION SLAPY2( X, Y )
!
!       .. Scalar Arguments ..
!       REAL               X, Y
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SLAPY2 returns sqrt(x**2+y**2), taking care not to cause unnecessary
!> overflow.
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
!>          X and Y specify the values x and y.
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
      REAL FUNCTION SLAPY2(X,Y)
      IMPLICIT NONE
!*--SLAPY267
!
!  -- LAPACK auxiliary routine (version 3.7.1) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     June 2017
!
!     .. Scalar Arguments ..
      REAL X , Y
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ZERO
      PARAMETER (ZERO=0.0E0)
      REAL ONE
      PARAMETER (ONE=1.0E0)
!     ..
!     .. Local Scalars ..
      REAL w , xabs , yabs , z
      LOGICAL x_is_nan , y_is_nan
!     ..
!     .. External Functions ..
      LOGICAL SISNAN
      EXTERNAL SISNAN
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , MAX , MIN , SQRT
!     ..
!     .. Executable Statements ..
!
!     ..
!     .. Executable Statements ..
!
      x_is_nan = SISNAN(X)
      y_is_nan = SISNAN(Y)
      IF ( x_is_nan ) SLAPY2 = X
      IF ( y_is_nan ) SLAPY2 = Y
!
      IF ( .NOT.(x_is_nan .OR. y_is_nan) ) THEN
         xabs = ABS(X)
         yabs = ABS(Y)
         w = MAX(xabs,yabs)
         z = MIN(xabs,yabs)
         IF ( z==ZERO ) THEN
            SLAPY2 = w
         ELSE
            SLAPY2 = w*SQRT(ONE+(z/w)**2)
         ENDIF
      ENDIF
!
!     End of SLAPY2
!
      END FUNCTION SLAPY2
