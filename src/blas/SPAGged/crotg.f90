!*==crotg.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b CROTG
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE CROTG(CA,CB,C,S)
!
!       .. Scalar Arguments ..
!       COMPLEX CA,CB,S
!       REAL C
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CROTG determines a complex Givens rotation.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in,out] CA
!> \verbatim
!>          CA is COMPLEX
!> \endverbatim
!>
!> \param[in] CB
!> \verbatim
!>          CB is COMPLEX
!> \endverbatim
!>
!> \param[out] C
!> \verbatim
!>          C is REAL
!> \endverbatim
!>
!> \param[out] S
!> \verbatim
!>          S is COMPLEX
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
!> \ingroup complex_blas_level1
!
!  =====================================================================
      SUBROUTINE CROTG(Ca,Cb,C,S)
      IMPLICIT NONE
!*--CROTG66
!
! Dummy argument declarations rewritten by SPAG
!
      COMPLEX , INTENT(INOUT) :: Ca
      COMPLEX , INTENT(IN) :: Cb
      REAL , INTENT(OUT) :: C
      COMPLEX , INTENT(OUT) :: S
!
! Local variable declarations rewritten by SPAG
!
      COMPLEX :: alpha
      REAL :: norm , scale
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
      IF ( CABS(Ca)==0. ) THEN
         C = 0.
         S = (1.,0.)
         Ca = Cb
      ELSE
         scale = CABS(Ca) + CABS(Cb)
         norm = scale*SQRT((CABS(Ca/scale))**2+(CABS(Cb/scale))**2)
         alpha = Ca/CABS(Ca)
         C = CABS(Ca)/norm
         S = alpha*CONJG(Cb)/norm
         Ca = alpha*norm
      ENDIF
      END SUBROUTINE CROTG
