!*==zlctes.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b ZLCTES
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       LOGICAL          FUNCTION ZLCTES( Z, D )
!
!       .. Scalar Arguments ..
!       COMPLEX*16         D, Z
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZLCTES returns .TRUE. if the eigenvalue Z/D is to be selected
!> (specifically, in this subroutine, if the real part of the
!> eigenvalue is negative), and otherwise it returns .FALSE..
!>
!> It is used by the test routine ZDRGES to test whether the driver
!> routine ZGGES successfully sorts eigenvalues.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] Z
!> \verbatim
!>          Z is COMPLEX*16
!>          The numerator part of a complex eigenvalue Z/D.
!> \endverbatim
!>
!> \param[in] D
!> \verbatim
!>          D is COMPLEX*16
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
!> \ingroup complex16_eig
!
!  =====================================================================
      LOGICAL FUNCTION ZLCTES(Z,D)
      IMPLICIT NONE
!*--ZLCTES62
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     June 2016
!
!     .. Scalar Arguments ..
      COMPLEX*16 D , Z
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
!
      DOUBLE PRECISION ZERO , ONE
      PARAMETER (ZERO=0.0D+0,ONE=1.0D+0)
      COMPLEX*16 CZERO
      PARAMETER (CZERO=(0.0D+0,0.0D+0))
!     ..
!     .. Local Scalars ..
      DOUBLE PRECISION zmax
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , DBLE , DIMAG , MAX , SIGN
!     ..
!     .. Executable Statements ..
!
      IF ( D==CZERO ) THEN
         ZLCTES = (DBLE(Z)<ZERO)
      ELSEIF ( DBLE(Z)==ZERO .OR. DBLE(D)==ZERO ) THEN
         ZLCTES = (SIGN(ONE,DIMAG(Z))/=SIGN(ONE,DIMAG(D)))
      ELSEIF ( DIMAG(Z)==ZERO .OR. DIMAG(D)==ZERO ) THEN
         ZLCTES = (SIGN(ONE,DBLE(Z))/=SIGN(ONE,DBLE(D)))
      ELSE
         zmax = MAX(ABS(DBLE(Z)),ABS(DIMAG(Z)))
         ZLCTES = ((DBLE(Z)/zmax)*DBLE(D)+(DIMAG(Z)/zmax)*DIMAG(D)<ZERO)
      ENDIF
!
!
!     End of ZLCTES
!
      END FUNCTION ZLCTES
