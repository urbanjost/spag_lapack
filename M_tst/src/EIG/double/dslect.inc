!*==dslect.f90  processed by SPAG 7.51RB at 17:37 on  4 Mar 2022
!> \brief \b dslect
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       LOGICAL          FUNCTION dslect( ZR, ZI )
!
!       .. Scalar Arguments ..
!       DOUBLE PRECISION   ZI, ZR
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> dslect returns .TRUE. if the eigenvalue ZR+sqrt(-1)*ZI is to be
!> selected, and otherwise it returns .FALSE.
!> It is used by DCHK41 to test if DGEES successfully sorts eigenvalues,
!> and by DCHK43 to test if DGEESX successfully sorts eigenvalues.
!>
!> The common block /SSLCT/ controls how eigenvalues are selected.
!> If SELOPT = 0, then dslect return .TRUE. when ZR is less than zero,
!> and .FALSE. otherwise.
!> If SELOPT is at least 1, dslect returns SELVAL(SELOPT) and adds 1
!> to SELOPT, cycling back to 1 at SELMAX.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] ZR
!> \verbatim
!>          ZR is DOUBLE PRECISION
!>          The real part of a complex eigenvalue ZR + i*ZI.
!> \endverbatim
!>
!> \param[in] ZI
!> \verbatim
!>          ZI is DOUBLE PRECISION
!>          The imaginary part of a complex eigenvalue ZR + i*ZI.
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
      LOGICAL FUNCTION dslect(Zr,Zi)
      IMPLICIT NONE
!*--DSLECT66
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     June 2016
!
!     .. Scalar Arguments ..
      DOUBLE PRECISION Zi , Zr
!     ..
!
!  =====================================================================
!
!     .. Arrays in Common ..
      LOGICAL SELval(20)
      DOUBLE PRECISION SELwi(20) , SELwr(20)
!     ..
!     .. Scalars in Common ..
      INTEGER SELdim , SELopt
!     ..
!     .. Common blocks ..
      COMMON /SSLCT / SELopt , SELdim , SELval , SELwr , SELwi
!     ..
!     .. Local Scalars ..
      INTEGER i
      DOUBLE PRECISION rmin , x
!     ..
!     .. Parameters ..
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
!     ..
!     .. External Functions ..
      DOUBLE PRECISION DLAPY2
      EXTERNAL DLAPY2
!     ..
!     .. Executable Statements ..
!
      IF ( SELopt==0 ) THEN
         dslect = (Zr<ZERO)
      ELSE
         rmin = DLAPY2(Zr-SELwr(1),Zi-SELwi(1))
         dslect = SELval(1)
         DO i = 2 , SELdim
            x = DLAPY2(Zr-SELwr(i),Zi-SELwi(i))
            IF ( x<=rmin ) THEN
               rmin = x
               dslect = SELval(i)
            ENDIF
         ENDDO
      ENDIF
!
!     End of dslect
!
      END FUNCTION dslect
