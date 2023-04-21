!*==crotg.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
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
!  -- Reference BLAS level1 routine (version 3.8.0) --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2017
!
!     .. Scalar Arguments ..
      COMPLEX Ca , Cb , S
      REAL C
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
      COMPLEX alpha
      REAL norm , scale
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC CABS , CONJG , SQRT
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
