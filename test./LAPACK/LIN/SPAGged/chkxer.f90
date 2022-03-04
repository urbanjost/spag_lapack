!*==chkxer.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b CHKXER
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE CHKXER( SRNAMT, INFOT, NOUT, LERR, OK )
!
!       .. Scalar Arguments ..
!       LOGICAL            LERR, OK
!       CHARACTER*(*)      SRNAMT
!       INTEGER            INFOT, NOUT
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
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
!> \date June 2017
!
!> \ingroup complex_lin
!
!  =====================================================================
      SUBROUTINE CHKXER(Srnamt,Infot,Nout,Lerr,Ok)
      IMPLICIT NONE
!*--CHKXER45
!
!  -- LAPACK test routine (version 3.7.1) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     June 2017
!
!     .. Scalar Arguments ..
      LOGICAL Lerr , Ok
      CHARACTER*(*) Srnamt
      INTEGER Infot , Nout
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC LEN_TRIM
!     ..
!     .. Executable Statements ..
      IF ( .NOT.Lerr ) THEN
         WRITE (Nout,FMT=99001) Infot , Srnamt(1:LEN_TRIM(Srnamt))
         Ok = .FALSE.
      ENDIF
      Lerr = .FALSE.
      RETURN
!
99001 FORMAT (' *** Illegal value of parameter number ',I2,             &
     &        ' not detected by ',A6,' ***')
!
!     End of CHKXER.
!
      END SUBROUTINE CHKXER
