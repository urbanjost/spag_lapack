!*==ilauplo.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b ILAUPLO
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ILAUPLO + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ilauplo.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ilauplo.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ilauplo.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       INTEGER FUNCTION ILAUPLO( UPLO )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> This subroutine translated from a character string specifying a
!> upper- or lower-triangular matrix to the relevant BLAST-specified
!> integer constant.
!>
!> ILAUPLO returns an INTEGER.  If ILAUPLO < 0, then the input is not
!> a character indicating an upper- or lower-triangular matrix.
!> Otherwise ILAUPLO returns the constant value corresponding to UPLO.
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
      FUNCTION ILAUPLO(Uplo)
      USE S_LSAME
      IMPLICIT NONE
!*--ILAUPLO63
!
! PARAMETER definitions rewritten by SPAG
!
      INTEGER , PARAMETER  ::  BLAS_UPPER = 121 , BLAS_LOWER = 122
      INTEGER :: ILAUPLO
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Uplo
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
      IF ( LSAME(Uplo,'U') ) THEN
         ILAUPLO = BLAS_UPPER
      ELSEIF ( LSAME(Uplo,'L') ) THEN
         ILAUPLO = BLAS_LOWER
      ELSE
         ILAUPLO = -1
      ENDIF
!
!     End of ILAUPLO
!
      END FUNCTION ILAUPLO
