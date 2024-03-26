!*==clctsx.f90  processed by SPAG 7.51RB at 17:37 on  4 Mar 2022
!> \brief \b clctsx
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       LOGICAL          FUNCTION clctsx( ALPHA, BETA )
!
!       .. Scalar Arguments ..
!       COMPLEX            ALPHA, BETA
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> This function is used to determine what eigenvalues will be
!> selected.  If this is part of the test driver cdrgsx, do not
!> change the code UNLESS you are testing input examples and not
!> using the built-in examples.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] ALPHA
!> \verbatim
!>          ALPHA is COMPLEX
!> \endverbatim
!>
!> \param[in] BETA
!> \verbatim
!>          BETA is COMPLEX
!>
!>          parameters to decide whether the pair (ALPHA, BETA) is
!>          selected.
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
!> \date December 2016
!
!> \ingroup complex_eig
!
!  =====================================================================
      LOGICAL FUNCTION clctsx(Alpha,Beta)
      IMPLICIT NONE
!*--CLCTSX61
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      COMPLEX Alpha , Beta
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
!     REAL               ZERO
!     PARAMETER          ( ZERO = 0.0E+0 )
!     COMPLEX            CZERO
!     PARAMETER          ( CZERO = ( 0.0E+0, 0.0E+0 ) )
!     ..
!     .. Scalars in Common ..
      LOGICAL FS
      INTEGER I , M , MPLusn , N
!     ..
!     .. Common blocks ..
      COMMON /MN    / M , N , MPLusn , I , FS
!     ..
!     .. Save statement ..
      SAVE 
!     ..
!     .. Executable Statements ..
!
      IF ( FS ) THEN
         I = I + 1
         IF ( I<=M ) THEN
            clctsx = .FALSE.
         ELSE
            clctsx = .TRUE.
         ENDIF
         IF ( I==MPLusn ) THEN
            FS = .FALSE.
            I = 0
         ENDIF
      ELSE
         I = I + 1
         IF ( I<=N ) THEN
            clctsx = .TRUE.
         ELSE
            clctsx = .FALSE.
         ENDIF
         IF ( I==MPLusn ) THEN
            FS = .TRUE.
            I = 0
         ENDIF
      ENDIF
!
!      IF( BETA.EQ.CZERO ) THEN
!         clctsx = ( REAL( ALPHA ).GT.ZERO )
!      ELSE
!         clctsx = ( REAL( ALPHA/BETA ).GT.ZERO )
!      END IF
!
!
!     End of clctsx
!
      END FUNCTION clctsx
