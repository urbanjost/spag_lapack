!*==cladiv.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b CLADIV performs complex division in real arithmetic, avoiding unnecessary overflow.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CLADIV + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cladiv.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cladiv.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cladiv.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       COMPLEX FUNCTION CLADIV( X, Y )
!
!       .. Scalar Arguments ..
!       COMPLEX            X, Y
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CLADIV := X / Y, where X and Y are complex.  The computation of X / Y
!> will not overflow on an intermediary step unless the results
!> overflows.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] X
!> \verbatim
!>          X is COMPLEX
!> \endverbatim
!>
!> \param[in] Y
!> \verbatim
!>          Y is COMPLEX
!>          The complex scalars X and Y.
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
!> \ingroup complexOTHERauxiliary
!
!  =====================================================================
      FUNCTION CLADIV(X,Y)
      USE S_SLADIV
      IMPLICIT NONE
!*--CLADIV69
      COMPLEX :: CLADIV
!
! Dummy argument declarations rewritten by SPAG
!
      COMPLEX , INTENT(IN) :: X
      COMPLEX , INTENT(IN) :: Y
!
! Local variable declarations rewritten by SPAG
!
      REAL :: zi , zr
!
! End of declarations rewritten by SPAG
!
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
!     ..
!     .. External Subroutines ..
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. Executable Statements ..
!
      CALL SLADIV(REAL(X),AIMAG(X),REAL(Y),AIMAG(Y),zr,zi)
      CLADIV = CMPLX(zr,zi)
!
!
!     End of CLADIV
!
      END FUNCTION CLADIV
