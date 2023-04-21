!*==sceil.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
 
!> \brief \b SCEIL
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       REAL FUNCTION SCEIL( A )
!
!       .. Scalar Arguments ..
!       REAL A
!       ..
!
!    =====================================================================
!
!       .. Intrinsic Functions ..
! 	      INTRINSIC          INT
!       ..
!       .. Executable Statements ..*
!
!       IF (A-INT(A).EQ.0) THEN
!           SCEIL = A
!       ELSE IF (A.GT.0) THEN
!           SCEIL = INT(A)+1;
!       ELSE
!           SCEIL = INT(A)
!       END IF
!
!       RETURN
!
!       END
!  Purpose
!  =======
!
!>\details \b Purpose:
!>\verbatim
!>\endverbatim
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
!> \ingroup variantsOTHERcomputational
!
!  =====================================================================
      REAL FUNCTION SCEIL(A)
      IMPLICIT NONE
!*--SCEIL63
!
!  -- LAPACK computational routine (version 3.1) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..*
      REAL A
!     ..
!
!  =====================================================================
!
!     .. Intrinsic Functions ..
      INTRINSIC INT
!     ..
!     .. Executable Statements ..*
!
      IF ( A==INT(A) ) THEN
         SCEIL = A
      ELSEIF ( A>0 ) THEN
         SCEIL = INT(A) + 1
      ELSE
         SCEIL = INT(A)
      ENDIF
 
!
      END FUNCTION SCEIL
