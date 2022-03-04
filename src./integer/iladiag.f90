!*==iladiag.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b ILADIAG
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ILADIAG + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/iladiag.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/iladiag.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/iladiag.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       INTEGER FUNCTION ILADIAG( DIAG )
!
!       .. Scalar Arguments ..
!       CHARACTER          DIAG
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> This subroutine translated from a character string specifying if a
!> matrix has unit diagonal or not to the relevant BLAST-specified
!> integer constant.
!>
!> ILADIAG returns an INTEGER.  If ILADIAG < 0, then the input is not a
!> character indicating a unit or non-unit diagonal.  Otherwise ILADIAG
!> returns the constant value corresponding to DIAG.
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
      FUNCTION ILADIAG(Diag)
      USE S_LSAME
      IMPLICIT NONE
!*--ILADIAG63
!
! PARAMETER definitions rewritten by SPAG
!
      INTEGER , PARAMETER  ::  BLAS_NON_UNIT_DIAG = 131 ,               &
     &                         BLAS_UNIT_DIAG = 132
      INTEGER :: ILADIAG
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Diag
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
      IF ( LSAME(Diag,'N') ) THEN
         ILADIAG = BLAS_NON_UNIT_DIAG
      ELSEIF ( LSAME(Diag,'U') ) THEN
         ILADIAG = BLAS_UNIT_DIAG
      ELSE
         ILADIAG = -1
      ENDIF
!
!     End of ILADIAG
!
      END FUNCTION ILADIAG
