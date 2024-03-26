!*==dlasum.f90  processed by SPAG 7.51RB at 17:37 on  4 Mar 2022
!> \brief \b dlasum
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE dlasum( TYPE, IOUNIT, IE, NRUN )
!
!       .. Scalar Arguments ..
!       CHARACTER*3        TYPE
!       INTEGER            IE, IOUNIT, NRUN
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> dlasum prints a summary of the results from one of the test routines.
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
!> \ingroup double_eig
!
!  =====================================================================
      SUBROUTINE dlasum(Type,Iounit,Ie,Nrun)
      IMPLICIT NONE
!*--DLASUM47
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER*3 Type
      INTEGER Ie , Iounit , Nrun
!     ..
!
!  =====================================================================
!
!     .. Executable Statements ..
!
      IF ( Ie>0 ) THEN
         WRITE (Iounit,FMT=99001) Type , ': ' , Ie , ' out of ' , Nrun ,&
     &                            ' tests failed to pass the threshold'
      ELSE
         WRITE (Iounit,FMT=99002) 'All tests for ' , Type ,             &
     &                            ' passed the threshold ( ' , Nrun ,   &
     &                            ' tests run)'
      ENDIF
99001 FORMAT (1X,A3,A2,I4,A8,I5,A35)
99002 FORMAT (/1X,A14,A3,A24,I5,A11)
!
!     End of dlasum
!
      END SUBROUTINE dlasum
