!*==sslect.f90  processed by SPAG 7.51RB at 17:37 on  4 Mar 2022
!> \brief \b SSLECT
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       LOGICAL          FUNCTION SSLECT( ZR, ZI )
!
!       .. Scalar Arguments ..
!       REAL               ZI, ZR
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SSLECT returns .TRUE. if the eigenvalue ZR+sqrt(-1)*ZI is to be
!> selected, and otherwise it returns .FALSE.
!> It is used by SCHK41 to test if SGEES successfully sorts eigenvalues,
!> and by SCHK43 to test if SGEESX successfully sorts eigenvalues.
!>
!> The common block /SSLCT/ controls how eigenvalues are selected.
!> If SELOPT = 0, then SSLECT return .TRUE. when ZR is less than zero,
!> and .FALSE. otherwise.
!> If SELOPT is at least 1, SSLECT returns SELVAL(SELOPT) and adds 1
!> to SELOPT, cycling back to 1 at SELMAX.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] ZR
!> \verbatim
!>          ZR is REAL
!>          The real part of a complex eigenvalue ZR + i*ZI.
!> \endverbatim
!>
!> \param[in] ZI
!> \verbatim
!>          ZI is REAL
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
!> \ingroup single_eig
!
!  =====================================================================
      LOGICAL FUNCTION SSLECT(Zr,Zi)
      IMPLICIT NONE
!*--SSLECT66
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     June 2016
!
!     .. Scalar Arguments ..
      REAL Zi , Zr
!     ..
!
!  =====================================================================
!
!     .. Arrays in Common ..
      LOGICAL SELval(20)
      REAL SELwi(20) , SELwr(20)
!     ..
!     .. Scalars in Common ..
      INTEGER SELdim , SELopt
!     ..
!     .. Common blocks ..
      COMMON /SSLCT / SELopt , SELdim , SELval , SELwr , SELwi
!     ..
!     .. Local Scalars ..
      INTEGER i
      REAL rmin , x
!     ..
!     .. Parameters ..
      REAL ZERO
      PARAMETER (ZERO=0.0E0)
!     ..
!     .. External Functions ..
      REAL SLAPY2
      EXTERNAL SLAPY2
!     ..
!     .. Executable Statements ..
!
      IF ( SELopt==0 ) THEN
         SSLECT = (Zr<ZERO)
      ELSE
         rmin = SLAPY2(Zr-SELwr(1),Zi-SELwi(1))
         SSLECT = SELval(1)
         DO i = 2 , SELdim
            x = SLAPY2(Zr-SELwr(i),Zi-SELwi(i))
            IF ( x<=rmin ) THEN
               rmin = x
               SSLECT = SELval(i)
            ENDIF
         ENDDO
      ENDIF
!
!     End of SSLECT
!
      END FUNCTION SSLECT
