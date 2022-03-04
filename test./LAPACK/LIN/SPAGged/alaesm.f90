!*==alaesm.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b ALAESM
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ALAESM( PATH, OK, NOUT )
!
!       .. Scalar Arguments ..
!       LOGICAL            OK
!       CHARACTER*3        PATH
!       INTEGER            NOUT
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ALAESM prints a summary of results from one of the -ERR- routines.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] PATH
!> \verbatim
!>          PATH is CHARACTER*3
!>          The LAPACK path name.
!> \endverbatim
!>
!> \param[in] OK
!> \verbatim
!>          OK is LOGICAL
!>          The flag from CHKXER that indicates whether or not the tests
!>          of error exits passed.
!> \endverbatim
!>
!> \param[in] NOUT
!> \verbatim
!>          NOUT is INTEGER
!>          The unit number on which results are to be printed.
!>          NOUT >= 0.
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
!> \ingroup aux_lin
!
!  =====================================================================
      SUBROUTINE ALAESM(Path,Ok,Nout)
      IMPLICIT NONE
!*--ALAESM67
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      LOGICAL Ok
      CHARACTER*3 Path
      INTEGER Nout
!     ..
!
!  =====================================================================
!
!     .. Executable Statements ..
!
      IF ( Ok ) THEN
         WRITE (Nout,FMT=99001) Path
      ELSE
         WRITE (Nout,FMT=99002) Path
      ENDIF
!
99001 FORMAT (1X,A3,' routines passed the tests of the error exits')
99002 FORMAT (' *** ',A3,' routines failed the tests of the error ',    &
     &        'exits ***')
!
!     End of ALAESM
!
      END SUBROUTINE ALAESM
