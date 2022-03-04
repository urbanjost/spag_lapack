!*==zladiv.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b ZLADIV performs complex division in real arithmetic, avoiding unnecessary overflow.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZLADIV + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zladiv.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zladiv.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zladiv.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       COMPLEX*16     FUNCTION ZLADIV( X, Y )
!
!       .. Scalar Arguments ..
!       COMPLEX*16         X, Y
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZLADIV := X / Y, where X and Y are complex.  The computation of X / Y
!> will not overflow on an intermediary step unless the results
!> overflows.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] X
!> \verbatim
!>          X is COMPLEX*16
!> \endverbatim
!>
!> \param[in] Y
!> \verbatim
!>          Y is COMPLEX*16
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
!> \ingroup complex16OTHERauxiliary
!
!  =====================================================================
      FUNCTION ZLADIV(X,Y)
      USE F77KINDS                        
      USE S_DLADIV
      IMPLICIT NONE
!*--ZLADIV70
      COMPLEX(CX16KIND) :: ZLADIV
!
! Dummy argument declarations rewritten by SPAG
!
      COMPLEX(CX16KIND) , INTENT(IN) :: X
      COMPLEX(CX16KIND) , INTENT(IN) :: Y
!
! Local variable declarations rewritten by SPAG
!
      REAL(R8KIND) :: zi , zr
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
      CALL DLADIV(DBLE(X),DIMAG(X),DBLE(Y),DIMAG(Y),zr,zi)
      ZLADIV = DCMPLX(zr,zi)
!
!
!     End of ZLADIV
!
      END FUNCTION ZLADIV
