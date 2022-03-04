!*==ilatrans.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b ILATRANS
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ILATRANS + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ilatrans.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ilatrans.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ilatrans.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       INTEGER FUNCTION ILATRANS( TRANS )
!
!       .. Scalar Arguments ..
!       CHARACTER          TRANS
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> This subroutine translates from a character string specifying a
!> transposition operation to the relevant BLAST-specified integer
!> constant.
!>
!> ILATRANS returns an INTEGER.  If ILATRANS < 0, then the input is not
!> a character indicating a transposition operator.  Otherwise ILATRANS
!> returns the constant value corresponding to TRANS.
!> \endverbatim
!
!  Arguments:
!  ==========
!
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
!> \ingroup auxOTHERcomputational
!
!  =====================================================================
      FUNCTION ILATRANS(Trans)
      USE S_LSAME
      IMPLICIT NONE
!*--ILATRANS63
!
! PARAMETER definitions rewritten by SPAG
!
      INTEGER , PARAMETER  ::  BLAS_NO_TRANS = 111 , BLAS_TRANS = 112 , &
     &                         BLAS_CONJ_TRANS = 113
      INTEGER :: ILATRANS
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Trans
!
! End of declarations rewritten by SPAG
!
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
!     ..
!     .. External Functions ..
!     ..
!     .. Executable Statements ..
      IF ( LSAME(Trans,'N') ) THEN
         ILATRANS = BLAS_NO_TRANS
      ELSEIF ( LSAME(Trans,'T') ) THEN
         ILATRANS = BLAS_TRANS
      ELSEIF ( LSAME(Trans,'C') ) THEN
         ILATRANS = BLAS_CONJ_TRANS
      ELSE
         ILATRANS = -1
      ENDIF
!
!     End of ILATRANS
!
      END FUNCTION ILATRANS
