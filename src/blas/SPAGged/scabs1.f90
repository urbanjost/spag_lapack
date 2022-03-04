!*==scabs1.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b SCABS1
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       REAL FUNCTION SCABS1(Z)
!
!       .. Scalar Arguments ..
!       COMPLEX Z
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SCABS1 computes |Re(.)| + |Im(.)| of a complex number
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] Z
!> \verbatim
!>          Z is COMPLEX
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
!> \date November 2017
!
!> \ingroup single_blas_level1
!
!  =====================================================================
      FUNCTION SCABS1(Z)
      IMPLICIT NONE
!*--SCABS150
      REAL :: SCABS1
!
! Dummy argument declarations rewritten by SPAG
!
      COMPLEX , INTENT(IN) :: Z
!
! End of declarations rewritten by SPAG
!
!     ..
!
!  =====================================================================
!
!     .. Intrinsic Functions ..
!     ..
      SCABS1 = ABS(REAL(Z)) + ABS(AIMAG(Z))
      END FUNCTION SCABS1
