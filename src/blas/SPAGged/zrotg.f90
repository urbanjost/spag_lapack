!*==zrotg.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b ZROTG
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZROTG(CA,CB,C,S)
!
!       .. Scalar Arguments ..
!       COMPLEX*16 CA,CB,S
!       DOUBLE PRECISION C
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    ZROTG determines a double complex Givens rotation.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in,out] CA
!> \verbatim
!>          CA is COMPLEX*16
!> \endverbatim
!>
!> \param[in] CB
!> \verbatim
!>          CB is COMPLEX*16
!> \endverbatim
!>
!> \param[out] C
!> \verbatim
!>          C is DOUBLE PRECISION
!> \endverbatim
!>
!> \param[out] S
!> \verbatim
!>          S is COMPLEX*16
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
!> \ingroup complex16_blas_level1
!
!  =====================================================================
      SUBROUTINE ZROTG(Ca,Cb,C,S)
      USE F77KINDS                        
      IMPLICIT NONE
!*--ZROTG67
!
! Dummy argument declarations rewritten by SPAG
!
      COMPLEX(CX16KIND) , INTENT(INOUT) :: Ca
      COMPLEX(CX16KIND) , INTENT(IN) :: Cb
      REAL(R8KIND) , INTENT(OUT) :: C
      COMPLEX(CX16KIND) , INTENT(OUT) :: S
!
! Local variable declarations rewritten by SPAG
!
      COMPLEX(CX16KIND) :: alpha
      REAL(R8KIND) :: norm , scale
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
      IF ( CDABS(Ca)==0.0D0 ) THEN
         C = 0.0D0
         S = (1.0D0,0.0D0)
         Ca = Cb
      ELSE
         scale = CDABS(Ca) + CDABS(Cb)
         norm = scale*DSQRT((CDABS(Ca/DCMPLX(scale,0.0D0)))             &
     &          **2+(CDABS(Cb/DCMPLX(scale,0.0D0)))**2)
         alpha = Ca/CDABS(Ca)
         C = CDABS(Ca)/norm
         S = alpha*DCONJG(Cb)/norm
         Ca = alpha*norm
      ENDIF
      END SUBROUTINE ZROTG
