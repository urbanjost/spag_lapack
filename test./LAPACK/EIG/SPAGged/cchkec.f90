!*==cchkec.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b CCHKEC
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE CCHKEC( THRESH, TSTERR, NIN, NOUT )
!
!       .. Scalar Arguments ..
!       LOGICAL            TSTERR
!       INTEGER            NIN, NOUT
!       REAL               THRESH
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CCHKEC tests eigen- condition estimation routines
!>        CTRSYL, CTREXC, CTRSNA, CTRSEN
!>
!> In all cases, the routine runs through a fixed set of numerical
!> examples, subjects them to various tests, and compares the test
!> results to a threshold THRESH. In addition, CTRSNA and CTRSEN are
!> tested by reading in precomputed examples from a file (on input unit
!> NIN).  Output is written to output unit NOUT.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] THRESH
!> \verbatim
!>          THRESH is REAL
!>          Threshold for residual tests.  A computed test ratio passes
!>          the threshold if it is less than THRESH.
!> \endverbatim
!>
!> \param[in] TSTERR
!> \verbatim
!>          TSTERR is LOGICAL
!>          Flag that indicates whether error exits are to be tested.
!> \endverbatim
!>
!> \param[in] NIN
!> \verbatim
!>          NIN is INTEGER
!>          The logical unit number for input.
!> \endverbatim
!>
!> \param[in] NOUT
!> \verbatim
!>          NOUT is INTEGER
!>          The logical unit number for output.
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
!> \ingroup complex_eig
!
!  =====================================================================
      SUBROUTINE CCHKEC(Thresh,Tsterr,Nin,Nout)
      IMPLICIT NONE
!*--CCHKEC79
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      LOGICAL Tsterr
      INTEGER Nin , Nout
      REAL Thresh
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
      LOGICAL ok
      CHARACTER*3 path
      INTEGER ktrexc , ktrsen , ktrsna , ktrsyl , ltrexc , ltrsyl ,     &
     &        ntests , ntrexc , ntrsyl
      REAL eps , rtrexc , rtrsyl , sfmin
!     ..
!     .. Local Arrays ..
      INTEGER ltrsen(3) , ltrsna(3) , ntrsen(3) , ntrsna(3)
      REAL rtrsen(3) , rtrsna(3)
!     ..
!     .. External Subroutines ..
      EXTERNAL CERREC , CGET35 , CGET36 , CGET37 , CGET38
!     ..
!     .. External Functions ..
      REAL SLAMCH
      EXTERNAL SLAMCH
!     ..
!     .. Executable Statements ..
!
      path(1:1) = 'Complex precision'
      path(2:3) = 'EC'
      eps = SLAMCH('P')
      sfmin = SLAMCH('S')
      WRITE (Nout,FMT=99006)
      WRITE (Nout,FMT=99007) eps , sfmin
      WRITE (Nout,FMT=99008) Thresh
!
!     Test error exits if TSTERR is .TRUE.
!
      IF ( Tsterr ) CALL CERREC(path,Nout)
!
      ok = .TRUE.
      CALL CGET35(rtrsyl,ltrsyl,ntrsyl,ktrsyl,Nin)
      IF ( rtrsyl>Thresh ) THEN
         ok = .FALSE.
         WRITE (Nout,FMT=99001) rtrsyl , ltrsyl , ntrsyl , ktrsyl
      ENDIF
!
      CALL CGET36(rtrexc,ltrexc,ntrexc,ktrexc,Nin)
      IF ( rtrexc>Thresh .OR. ntrexc>0 ) THEN
         ok = .FALSE.
         WRITE (Nout,FMT=99002) rtrexc , ltrexc , ntrexc , ktrexc
      ENDIF
!
      CALL CGET37(rtrsna,ltrsna,ntrsna,ktrsna,Nin)
      IF ( rtrsna(1)>Thresh .OR. rtrsna(2)>Thresh .OR. ntrsna(1)/=0 .OR.&
     &     ntrsna(2)/=0 .OR. ntrsna(3)/=0 ) THEN
         ok = .FALSE.
         WRITE (Nout,FMT=99003) rtrsna , ltrsna , ntrsna , ktrsna
      ENDIF
!
      CALL CGET38(rtrsen,ltrsen,ntrsen,ktrsen,Nin)
      IF ( rtrsen(1)>Thresh .OR. rtrsen(2)>Thresh .OR. ntrsen(1)/=0 .OR.&
     &     ntrsen(2)/=0 .OR. ntrsen(3)/=0 ) THEN
         ok = .FALSE.
         WRITE (Nout,FMT=99004) rtrsen , ltrsen , ntrsen , ktrsen
      ENDIF
!
      ntests = ktrsyl + ktrexc + ktrsna + ktrsen
      IF ( ok ) WRITE (Nout,FMT=99005) path , ntests
!
99001 FORMAT (' Error in CTRSYL: RMAX =',E12.3,/' LMAX = ',I8,' NINFO=',&
     &        I8,' KNT=',I8)
99002 FORMAT (' Error in CTREXC: RMAX =',E12.3,/' LMAX = ',I8,' NINFO=',&
     &        I8,' KNT=',I8)
99003 FORMAT (' Error in CTRSNA: RMAX =',3E12.3,/' LMAX = ',3I8,        &
     &        ' NINFO=',3I8,' KNT=',I8)
99004 FORMAT (' Error in CTRSEN: RMAX =',3E12.3,/' LMAX = ',3I8,        &
     &        ' NINFO=',3I8,' KNT=',I8)
99005 FORMAT (/1X,'All tests for ',A3,                                  &
     &        ' routines passed the threshold ( ',I6,' tests run)')
99006 FORMAT (' Tests of the Nonsymmetric eigenproblem condition',      &
     &        ' estimation routines',/' CTRSYL, CTREXC, CTRSNA, CTRSEN',&
     &        /)
99007 FORMAT (' Relative machine precision (EPS) = ',E16.6,             &
     &        /' Safe minimum (SFMIN)             = ',E16.6,/)
99008 FORMAT (' Routines pass computational tests if test ratio is ',   &
     &        'less than',F8.2,//)
!
!     End of CCHKEC
!
      END SUBROUTINE CCHKEC
