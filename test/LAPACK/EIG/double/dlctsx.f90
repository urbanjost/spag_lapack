!*==dlctsx.f90  processed by SPAG 7.51RB at 17:37 on  4 Mar 2022
!> \brief \b DLCTSX
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       LOGICAL          FUNCTION DLCTSX( AR, AI, BETA )
!
!       .. Scalar Arguments ..
!       DOUBLE PRECISION   AI, AR, BETA
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> This function is used to determine what eigenvalues will be
!> selected.  If this is part of the test driver DDRGSX, do not
!> change the code UNLESS you are testing input examples and not
!> using the built-in examples.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] AR
!> \verbatim
!>          AR is DOUBLE PRECISION
!>          The numerator of the real part of a complex eigenvalue
!>          (AR/BETA) + i*(AI/BETA).
!> \endverbatim
!>
!> \param[in] AI
!> \verbatim
!>          AI is DOUBLE PRECISION
!>          The numerator of the imaginary part of a complex eigenvalue
!>          (AR/BETA) + i*(AI).
!> \endverbatim
!>
!> \param[in] BETA
!> \verbatim
!>          BETA is DOUBLE PRECISION
!>          The denominator part of a complex eigenvalue
!>          (AR/BETA) + i*(AI/BETA).
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
!> \ingroup double_eig
!
!  =====================================================================
      LOGICAL FUNCTION DLCTSX(Ar,Ai,Beta)
      IMPLICIT NONE
!*--DLCTSX69
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      DOUBLE PRECISION Ai , Ar , Beta
!     ..
!
!  =====================================================================
!
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
            DLCTSX = .FALSE.
         ELSE
            DLCTSX = .TRUE.
         ENDIF
         IF ( I==MPLusn ) THEN
            FS = .FALSE.
            I = 0
         ENDIF
      ELSE
         I = I + 1
         IF ( I<=N ) THEN
            DLCTSX = .TRUE.
         ELSE
            DLCTSX = .FALSE.
         ENDIF
         IF ( I==MPLusn ) THEN
            FS = .TRUE.
            I = 0
         ENDIF
      ENDIF
!
!       IF( AR/BETA.GT.0.0 )THEN
!          DLCTSX = .TRUE.
!       ELSE
!          DLCTSX = .FALSE.
!       END IF
!
!
!     End of DLCTSX
!
      END FUNCTION DLCTSX
