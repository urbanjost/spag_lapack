!*==slctes.f90  processed by SPAG 7.51RB at 17:37 on  4 Mar 2022
!> \brief \b slctes
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       LOGICAL          FUNCTION slctes( ZR, ZI, D )
!
!       .. Scalar Arguments ..
!       REAL               D, ZI, ZR
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> slctes returns .TRUE. if the eigenvalue (ZR/D) + sqrt(-1)*(ZI/D)
!> is to be selected (specifically, in this subroutine, if the real
!> part of the eigenvalue is negative), and otherwise it returns
!> .FALSE..
!>
!> It is used by the test routine sdrges to test whether the driver
!> routine SGGES successfully sorts eigenvalues.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] ZR
!> \verbatim
!>          ZR is REAL
!>          The numerator of the real part of a complex eigenvalue
!>          (ZR/D) + i*(ZI/D).
!> \endverbatim
!>
!> \param[in] ZI
!> \verbatim
!>          ZI is REAL
!>          The numerator of the imaginary part of a complex eigenvalue
!>          (ZR/D) + i*(ZI).
!> \endverbatim
!>
!> \param[in] D
!> \verbatim
!>          D is REAL
!>          The denominator part of a complex eigenvalue
!>          (ZR/D) + i*(ZI/D).
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
!> \ingroup single_eig
!
!  =====================================================================
      LOGICAL FUNCTION slctes(Zr,Zi,D)
      IMPLICIT NONE
!*--SLCTES72
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     June 2016
!
!     .. Scalar Arguments ..
      REAL D , Zi , Zr
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ZERO , ONE
      PARAMETER (ZERO=0.0E+0,ONE=1.0E+0)
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC SIGN
!     ..
!     .. Executable Statements ..
!
      IF ( D==ZERO ) THEN
         slctes = (Zr<ZERO)
      ELSE
         slctes = (SIGN(ONE,Zr)/=SIGN(ONE,D))
      ENDIF
!
!
!     End of slctes
!
      END FUNCTION slctes
