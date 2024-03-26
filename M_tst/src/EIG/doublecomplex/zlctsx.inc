!*==zlctsx.f90  processed by SPAG 7.51RB at 17:37 on  4 Mar 2022
!> \brief \b zlctsx
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       LOGICAL          FUNCTION zlctsx( ALPHA, BETA )
!
!       .. Scalar Arguments ..
!       COMPLEX*16         ALPHA, BETA
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> This function is used to determine what eigenvalues will be
!> selected.  If this is part of the test driver zdrgsx, do not
!> change the code UNLESS you are testing input examples and not
!> using the built-in examples.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] ALPHA
!> \verbatim
!>          ALPHA is COMPLEX*16
!> \endverbatim
!>
!> \param[in] BETA
!> \verbatim
!>          BETA is COMPLEX*16
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
!> \ingroup complex16_eig
!
!  =====================================================================
      LOGICAL FUNCTION zlctsx(Alpha,Beta)
      IMPLICIT NONE
!*--ZLCTSX61
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      COMPLEX*16 Alpha , Beta
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
!     DOUBLE PRECISION               ZERO
!     PARAMETER          ( ZERO = 0.0E+0 )
!     COMPLEX*16            CZERO
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
            zlctsx = .FALSE.
         ELSE
            zlctsx = .TRUE.
         ENDIF
         IF ( I==MPLusn ) THEN
            FS = .FALSE.
            I = 0
         ENDIF
      ELSE
         I = I + 1
         IF ( I<=N ) THEN
            zlctsx = .TRUE.
         ELSE
            zlctsx = .FALSE.
         ENDIF
         IF ( I==MPLusn ) THEN
            FS = .TRUE.
            I = 0
         ENDIF
      ENDIF
!
!      IF( BETA.EQ.CZERO ) THEN
!         zlctsx = ( DBLE( ALPHA ).GT.ZERO )
!      ELSE
!         zlctsx = ( DBLE( ALPHA/BETA ).GT.ZERO )
!      END IF
!
!
!     End of zlctsx
!
      END FUNCTION zlctsx
