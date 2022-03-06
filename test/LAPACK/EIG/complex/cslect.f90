!*==cslect.f90  processed by SPAG 7.51RB at 17:37 on  4 Mar 2022
!> \brief \b CSLECT
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       LOGICAL          FUNCTION CSLECT( Z )
!
!       .. Scalar Arguments ..
!       COMPLEX            Z
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CSLECT returns .TRUE. if the eigenvalue Z is to be selected,
!> otherwise it returns .FALSE.
!> It is used by CCHK41 to test if CGEES successfully sorts eigenvalues,
!> and by CCHK43 to test if CGEESX successfully sorts eigenvalues.
!>
!> The common block /SSLCT/ controls how eigenvalues are selected.
!> If SELOPT = 0, then CSLECT return .TRUE. when real(Z) is less than
!> zero, and .FALSE. otherwise.
!> If SELOPT is at least 1, CSLECT returns SELVAL(SELOPT) and adds 1
!> to SELOPT, cycling back to 1 at SELMAX.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] Z
!> \verbatim
!>          Z is COMPLEX
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
!> \ingroup complex_eig
!
!  =====================================================================
      LOGICAL FUNCTION CSLECT(Z)
      IMPLICIT NONE
!*--CSLECT60
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     June 2016
!
!     .. Scalar Arguments ..
      COMPLEX Z
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ZERO
      PARAMETER (ZERO=0.0E0)
!     ..
!     .. Local Scalars ..
      INTEGER i
      REAL rmin , x
!     ..
!     .. Scalars in Common ..
      INTEGER SELdim , SELopt
!     ..
!     .. Arrays in Common ..
      LOGICAL SELval(20)
      REAL SELwi(20) , SELwr(20)
!     ..
!     .. Common blocks ..
      COMMON /SSLCT / SELopt , SELdim , SELval , SELwr , SELwi
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , CMPLX , REAL
!     ..
!     .. Executable Statements ..
!
      IF ( SELopt==0 ) THEN
         CSLECT = (REAL(Z)<ZERO)
      ELSE
         rmin = ABS(Z-CMPLX(SELwr(1),SELwi(1)))
         CSLECT = SELval(1)
         DO i = 2 , SELdim
            x = ABS(Z-CMPLX(SELwr(i),SELwi(i)))
            IF ( x<=rmin ) THEN
               rmin = x
               CSLECT = SELval(i)
            ENDIF
         ENDDO
      ENDIF
!
!     End of CSLECT
!
      END FUNCTION CSLECT
