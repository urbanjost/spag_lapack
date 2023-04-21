!*==ieeeck.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b IEEECK
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download IEEECK + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ieeeck.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ieeeck.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ieeeck.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       INTEGER          FUNCTION IEEECK( ISPEC, ZERO, ONE )
!
!       .. Scalar Arguments ..
!       INTEGER            ISPEC
!       REAL               ONE, ZERO
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> IEEECK is called from the ILAENV to verify that Infinity and
!> possibly NaN arithmetic is safe (i.e. will not trap).
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] ISPEC
!> \verbatim
!>          ISPEC is INTEGER
!>          Specifies whether to test just for inifinity arithmetic
!>          or whether to test for infinity and NaN arithmetic.
!>          = 0: Verify infinity arithmetic only.
!>          = 1: Verify infinity and NaN arithmetic.
!> \endverbatim
!>
!> \param[in] ZERO
!> \verbatim
!>          ZERO is REAL
!>          Must contain the value 0.0
!>          This is passed to prevent the compiler from optimizing
!>          away this code.
!> \endverbatim
!>
!> \param[in] ONE
!> \verbatim
!>          ONE is REAL
!>          Must contain the value 1.0
!>          This is passed to prevent the compiler from optimizing
!>          away this code.
!>
!>  RETURN VALUE:  INTEGER
!>          = 0:  Arithmetic failed to produce the correct answers
!>          = 1:  Arithmetic produced the correct answers
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
!> \ingroup OTHERauxiliary
!
!  =====================================================================
      INTEGER FUNCTION IEEECK(Ispec,Zero,One)
      IMPLICIT NONE
!*--IEEECK86
!
!  -- LAPACK auxiliary routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Ispec
      REAL One , Zero
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
      REAL nan1 , nan2 , nan3 , nan4 , nan5 , nan6 , neginf , negzro ,  &
     &     newzro , posinf
!     ..
!     .. Executable Statements ..
      IEEECK = 1
!
      posinf = One/Zero
      IF ( posinf<=One ) THEN
         IEEECK = 0
         RETURN
      ENDIF
!
      neginf = -One/Zero
      IF ( neginf>=Zero ) THEN
         IEEECK = 0
         RETURN
      ENDIF
!
      negzro = One/(neginf+One)
      IF ( negzro/=Zero ) THEN
         IEEECK = 0
         RETURN
      ENDIF
!
      neginf = One/negzro
      IF ( neginf>=Zero ) THEN
         IEEECK = 0
         RETURN
      ENDIF
!
      newzro = negzro + Zero
      IF ( newzro/=Zero ) THEN
         IEEECK = 0
         RETURN
      ENDIF
!
      posinf = One/newzro
      IF ( posinf<=One ) THEN
         IEEECK = 0
         RETURN
      ENDIF
!
      neginf = neginf*posinf
      IF ( neginf>=Zero ) THEN
         IEEECK = 0
         RETURN
      ENDIF
!
      posinf = posinf*posinf
      IF ( posinf<=One ) THEN
         IEEECK = 0
         RETURN
      ENDIF
!
!
!
!
!     Return if we were only asked to check infinity arithmetic
!
      IF ( Ispec==0 ) RETURN
!
      nan1 = posinf + neginf
!
      nan2 = posinf/neginf
!
      nan3 = posinf/posinf
!
      nan4 = posinf*Zero
!
      nan5 = neginf*negzro
!
      nan6 = nan5*Zero
!
      IF ( nan1==nan1 ) THEN
         IEEECK = 0
         RETURN
      ENDIF
!
      IF ( nan2==nan2 ) THEN
         IEEECK = 0
         RETURN
      ENDIF
!
      IF ( nan3==nan3 ) THEN
         IEEECK = 0
         RETURN
      ENDIF
!
      IF ( nan4==nan4 ) THEN
         IEEECK = 0
         RETURN
      ENDIF
!
      IF ( nan5==nan5 ) THEN
         IEEECK = 0
         RETURN
      ENDIF
!
      IF ( nan6==nan6 ) THEN
         IEEECK = 0
         RETURN
      ENDIF
!
      END FUNCTION IEEECK
