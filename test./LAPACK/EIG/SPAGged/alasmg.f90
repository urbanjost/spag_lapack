!*==alasmg.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b ALASMG
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ALASMG( TYPE, NOUT, NFAIL, NRUN, NERRS )
!
!       .. Scalar Arguments ..
!       CHARACTER*3        TYPE
!       INTEGER            NFAIL, NOUT, NRUN, NERRS
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ALASMG prints a summary of results from one of the -CHK- routines.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] TYPE
!> \verbatim
!>          TYPE is CHARACTER*3
!>          The LAPACK path name.
!> \endverbatim
!>
!> \param[in] NOUT
!> \verbatim
!>          NOUT is INTEGER
!>          The unit number on which results are to be printed.
!>          NOUT >= 0.
!> \endverbatim
!>
!> \param[in] NFAIL
!> \verbatim
!>          NFAIL is INTEGER
!>          The number of tests which did not pass the threshold ratio.
!> \endverbatim
!>
!> \param[in] NRUN
!> \verbatim
!>          NRUN is INTEGER
!>          The total number of tests.
!> \endverbatim
!>
!> \param[in] NERRS
!> \verbatim
!>          NERRS is INTEGER
!>          The number of error messages recorded.
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
!> \ingroup aux_eig
!
!  =====================================================================
      SUBROUTINE ALASMG(Type,Nout,Nfail,Nrun,Nerrs)
      IMPLICIT NONE
!*--ALASMG77
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER*3 Type
      INTEGER Nfail , Nout , Nrun , Nerrs
!     ..
!
! ======================================================================
!
!     .. Executable Statements ..
!
      IF ( Nfail>0 ) THEN
         WRITE (Nout,FMT=99001) Type , Nfail , Nrun
      ELSE
         WRITE (Nout,FMT=99002) Type , Nrun
      ENDIF
      IF ( Nerrs>0 ) WRITE (Nout,FMT=99003) Nerrs
!
99001 FORMAT (1X,A3,': ',I6,' out of ',I6,                              &
     &        ' tests failed to pass the threshold')
99002 FORMAT (/1X,'All tests for ',A3,                                  &
     &        ' routines passed the threshold ( ',I6,' tests run)')
99003 FORMAT (6X,I6,' error messages recorded')
!
!     End of ALASMG
!
      END SUBROUTINE ALASMG
