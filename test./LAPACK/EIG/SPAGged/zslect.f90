!*==zslect.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b ZSLECT
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       LOGICAL          FUNCTION ZSLECT( Z )
!
!       .. Scalar Arguments ..
!       COMPLEX*16         Z
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZSLECT returns .TRUE. if the eigenvalue Z is to be selected,
!> otherwise it returns .FALSE.
!> It is used by ZCHK41 to test if ZGEES successfully sorts eigenvalues,
!> and by ZCHK43 to test if ZGEESX successfully sorts eigenvalues.
!>
!> The common block /SSLCT/ controls how eigenvalues are selected.
!> If SELOPT = 0, then ZSLECT return .TRUE. when real(Z) is less than
!> zero, and .FALSE. otherwise.
!> If SELOPT is at least 1, ZSLECT returns SELVAL(SELOPT) and adds 1
!> to SELOPT, cycling back to 1 at SELMAX.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] Z
!> \verbatim
!>          Z is COMPLEX*16
!>          The eigenvalue Z.
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
      LOGICAL FUNCTION ZSLECT(Z)
      IMPLICIT NONE
!*--ZSLECT60
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     June 2016
!
!     .. Scalar Arguments ..
      COMPLEX*16 Z
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
!     ..
!     .. Local Scalars ..
      INTEGER i
      DOUBLE PRECISION rmin , x
!     ..
!     .. Scalars in Common ..
      INTEGER SELdim , SELopt
!     ..
!     .. Arrays in Common ..
      LOGICAL SELval(20)
      DOUBLE PRECISION SELwi(20) , SELwr(20)
!     ..
!     .. Common blocks ..
      COMMON /SSLCT / SELopt , SELdim , SELval , SELwr , SELwi
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , DBLE , DCMPLX
!     ..
!     .. Executable Statements ..
!
      IF ( SELopt==0 ) THEN
         ZSLECT = (DBLE(Z)<ZERO)
      ELSE
         rmin = ABS(Z-DCMPLX(SELwr(1),SELwi(1)))
         ZSLECT = SELval(1)
         DO i = 2 , SELdim
            x = ABS(Z-DCMPLX(SELwr(i),SELwi(i)))
            IF ( x<=rmin ) THEN
               rmin = x
               ZSLECT = SELval(i)
            ENDIF
         ENDDO
      ENDIF
!
!     End of ZSLECT
!
      END FUNCTION ZSLECT
