!*==dcabs1.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b DCABS1
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       DOUBLE PRECISION FUNCTION DCABS1(Z)
!
!       .. Scalar Arguments ..
!       COMPLEX*16 Z
!       ..
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DCABS1 computes |Re(.)| + |Im(.)| of a double complex number
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] Z
!> \verbatim
!>          Z is COMPLEX*16
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
!> \ingroup double_blas_level1
!
!  =====================================================================
      FUNCTION DCABS1(Z)
      USE F77KINDS                        
      IMPLICIT NONE
!*--DCABS152
      REAL(R8KIND) :: DCABS1
!
! Dummy argument declarations rewritten by SPAG
!
      COMPLEX(CX16KIND) , INTENT(IN) :: Z
!
! End of declarations rewritten by SPAG
!
!     ..
!     ..
!  =====================================================================
!
!     .. Intrinsic Functions ..
!
      DCABS1 = ABS(DBLE(Z)) + ABS(DIMAG(Z))
      END FUNCTION DCABS1
