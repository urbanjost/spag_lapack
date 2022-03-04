!*==srotg.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b SROTG
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE SROTG(SA,SB,C,S)
!
!       .. Scalar Arguments ..
!       REAL C,S,SA,SB
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    SROTG construct givens plane rotation.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in,out] SA
!> \verbatim
!>          SA is REAL
!> \endverbatim
!>
!> \param[in,out] SB
!> \verbatim
!>          SB is REAL
!> \endverbatim
!>
!> \param[out] C
!> \verbatim
!>          C is REAL
!> \endverbatim
!>
!> \param[out] S
!> \verbatim
!>          S is REAL
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
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>     jack dongarra, linpack, 3/11/78.
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE SROTG(Sa,Sb,C,S)
      IMPLICIT NONE
!*--SROTG73
!
! Dummy argument declarations rewritten by SPAG
!
      REAL , INTENT(INOUT) :: Sa
      REAL , INTENT(INOUT) :: Sb
      REAL , INTENT(INOUT) :: C
      REAL , INTENT(INOUT) :: S
!
! Local variable declarations rewritten by SPAG
!
      REAL :: r , roe , scale , z
!
! End of declarations rewritten by SPAG
!
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
!     ..
!     .. Intrinsic Functions ..
!     ..
      scale = ABS(Sa) + ABS(Sb)
      IF ( scale==0.0 ) THEN
         C = 1.0
         S = 0.0
         r = 0.0
         z = 0.0
      ELSE
         roe = Sb
         IF ( ABS(Sa)>ABS(Sb) ) roe = Sa
         r = scale*SQRT((Sa/scale)**2+(Sb/scale)**2)
         r = SIGN(1.0,roe)*r
         C = Sa/r
         S = Sb/r
         z = 1.0
         IF ( ABS(Sa)>ABS(Sb) ) z = S
         IF ( ABS(Sb)>=ABS(Sa) .AND. C/=0.0 ) z = 1.0/C
      ENDIF
      Sa = r
      Sb = z
      END SUBROUTINE SROTG
