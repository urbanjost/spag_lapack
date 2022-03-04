!*==dcombssq.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b DCOMBSSQ adds two scaled sum of squares quantities.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!
!  Definition:
!  ===========
!
!       SUBROUTINE DCOMBSSQ( V1, V2 )
!
!       .. Array Arguments ..
!       DOUBLE PRECISION   V1( 2 ), V2( 2 )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DCOMBSSQ adds two scaled sum of squares quantities, V1 := V1 + V2.
!> That is,
!>
!>    V1_scale**2 * V1_sumsq := V1_scale**2 * V1_sumsq
!>                            + V2_scale**2 * V2_sumsq
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in,out] V1
!> \verbatim
!>          V1 is DOUBLE PRECISION array, dimension (2).
!>          The first scaled sum.
!>          V1(1) = V1_scale, V1(2) = V1_sumsq.
!> \endverbatim
!>
!> \param[in] V2
!> \verbatim
!>          V2 is DOUBLE PRECISION array, dimension (2).
!>          The second scaled sum.
!>          V2(1) = V2_scale, V2(2) = V2_sumsq.
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
!> \date November 2018
!
!> \ingroup OTHERauxiliary
!
!  =====================================================================
      SUBROUTINE DCOMBSSQ(V1,V2)
      USE F77KINDS                        
      IMPLICIT NONE
!*--DCOMBSSQ65
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0
!
! Dummy argument declarations rewritten by SPAG
!
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(2) :: V1
      REAL(R8KIND) , INTENT(IN) , DIMENSION(2) :: V2
!
! End of declarations rewritten by SPAG
!
!     ..
!
! =====================================================================
!
!     .. Parameters ..
!     ..
!     .. Executable Statements ..
!
      IF ( V1(1)<V2(1) ) THEN
         V1(2) = V2(2) + (V1(1)/V2(1))**2*V1(2)
         V1(1) = V2(1)
      ELSEIF ( V1(1)/=ZERO ) THEN
         V1(2) = V1(2) + (V2(1)/V1(1))**2*V2(2)
      ELSE
         V1(2) = V1(2) + V2(2)
      ENDIF
!
!     End of DCOMBSSQ
!
      END SUBROUTINE DCOMBSSQ
