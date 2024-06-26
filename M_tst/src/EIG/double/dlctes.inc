!*==dlctes.f90  processed by SPAG 7.51RB at 17:37 on  4 Mar 2022
!> \brief \b dlctes
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       LOGICAL          FUNCTION dlctes( ZR, ZI, D )
!
!       .. Scalar Arguments ..
!       DOUBLE PRECISION   D, ZI, ZR
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> dlctes returns .TRUE. if the eigenvalue (ZR/D) + sqrt(-1)*(ZI/D)
!> is to be selected (specifically, in this subroutine, if the real
!> part of the eigenvalue is negative), and otherwise it returns
!> .FALSE..
!>
!> It is used by the test routine ddrges to test whether the driver
!> routine DGGES successfully sorts eigenvalues.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] ZR
!> \verbatim
!>          ZR is DOUBLE PRECISION
!>          The numerator of the real part of a complex eigenvalue
!>          (ZR/D) + i*(ZI/D).
!> \endverbatim
!>
!> \param[in] ZI
!> \verbatim
!>          ZI is DOUBLE PRECISION
!>          The numerator of the imaginary part of a complex eigenvalue
!>          (ZR/D) + i*(ZI).
!> \endverbatim
!>
!> \param[in] D
!> \verbatim
!>          D is DOUBLE PRECISION
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
!> \ingroup double_eig
!
!  =====================================================================
      LOGICAL FUNCTION dlctes(Zr,Zi,D)
      IMPLICIT NONE
!*--DLCTES72
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     June 2016
!
!     .. Scalar Arguments ..
      DOUBLE PRECISION D , Zi , Zr
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO , ONE
      PARAMETER (ZERO=0.0D+0,ONE=1.0D+0)
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC SIGN
!     ..
!     .. Executable Statements ..
!
      IF ( D==ZERO ) THEN
         dlctes = (Zr<ZERO)
      ELSE
         dlctes = (SIGN(ONE,Zr)/=SIGN(ONE,D))
      ENDIF
!
!
!     End of dlctes
!
      END FUNCTION dlctes
