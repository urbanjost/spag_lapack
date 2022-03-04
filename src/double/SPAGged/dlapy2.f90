!*==dlapy2.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b DLAPY2 returns sqrt(x2+y2).
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DLAPY2 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlapy2.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlapy2.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlapy2.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       DOUBLE PRECISION FUNCTION DLAPY2( X, Y )
!
!       .. Scalar Arguments ..
!       DOUBLE PRECISION   X, Y
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLAPY2 returns sqrt(x**2+y**2), taking care not to cause unnecessary
!> overflow.
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
      FUNCTION DLAPY2(X,Y)
      USE F77KINDS                        
      USE S_DISNAN
      IMPLICIT NONE
!*--DLAPY269
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D0 , ONE = 1.0D0
      REAL(R8KIND) :: DLAPY2
!
! Dummy argument declarations rewritten by SPAG
!
      REAL(R8KIND) :: X
      REAL(R8KIND) :: Y
!
! Local variable declarations rewritten by SPAG
!
      REAL(R8KIND) :: w , xabs , yabs , z
      LOGICAL :: x_is_nan , y_is_nan
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
!     .. External Functions ..
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. Executable Statements ..
!
      x_is_nan = DISNAN(X)
      y_is_nan = DISNAN(Y)
      IF ( x_is_nan ) DLAPY2 = X
      IF ( y_is_nan ) DLAPY2 = Y
!
      IF ( .NOT.(x_is_nan .OR. y_is_nan) ) THEN
         xabs = ABS(X)
         yabs = ABS(Y)
         w = MAX(xabs,yabs)
         z = MIN(xabs,yabs)
         IF ( z==ZERO ) THEN
            DLAPY2 = w
         ELSE
            DLAPY2 = w*SQRT(ONE+(z/w)**2)
         ENDIF
      ENDIF
!
!     End of DLAPY2
!
      END FUNCTION DLAPY2
