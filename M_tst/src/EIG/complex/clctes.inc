!*==clctes.f90  processed by SPAG 7.51RB at 17:37 on  4 Mar 2022
!> \brief \b clctes
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       LOGICAL          FUNCTION clctes( Z, D )
!
!       .. Scalar Arguments ..
!       COMPLEX            D, Z
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> clctes returns .TRUE. if the eigenvalue Z/D is to be selected
!> (specifically, in this subroutine, if the real part of the
!> eigenvalue is negative), and otherwise it returns .FALSE..
!>
!> It is used by the test routine cdrges to test whether the driver
!> routine CGGES successfully sorts eigenvalues.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] Z
!> \verbatim
!>          Z is COMPLEX
!>          The numerator part of a complex eigenvalue Z/D.
!> \endverbatim
!>
!> \param[in] D
!> \verbatim
!>          D is COMPLEX
!>          The denominator part of a complex eigenvalue Z/D.
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
!> \date June 2016
!
!> \ingroup complex_eig
!
!  =====================================================================
      LOGICAL FUNCTION clctes(Z,D)
      IMPLICIT NONE
!*--CLCTES62
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     June 2016
!
!     .. Scalar Arguments ..
      COMPLEX D , Z
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
!
      REAL ZERO , ONE
      PARAMETER (ZERO=0.0E+0,ONE=1.0E+0)
      COMPLEX CZERO
      PARAMETER (CZERO=(0.0E+0,0.0E+0))
!     ..
!     .. Local Scalars ..
      REAL zmax
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , AIMAG , MAX , REAL , SIGN
!     ..
!     .. Executable Statements ..
!
      IF ( D==CZERO ) THEN
         clctes = (REAL(Z)<ZERO)
      ELSEIF ( REAL(Z)==ZERO .OR. REAL(D)==ZERO ) THEN
         clctes = (SIGN(ONE,AIMAG(Z))/=SIGN(ONE,AIMAG(D)))
      ELSEIF ( AIMAG(Z)==ZERO .OR. AIMAG(D)==ZERO ) THEN
         clctes = (SIGN(ONE,REAL(Z))/=SIGN(ONE,REAL(D)))
      ELSE
         zmax = MAX(ABS(REAL(Z)),ABS(AIMAG(Z)))
         clctes = ((REAL(Z)/zmax)*REAL(D)+(AIMAG(Z)/zmax)*AIMAG(D)<ZERO)
      ENDIF
!
!
!     End of clctes
!
      END FUNCTION clctes
