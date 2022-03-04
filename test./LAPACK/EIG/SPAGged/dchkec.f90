!*==dchkec.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b DCHKEC
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DCHKEC( THRESH, TSTERR, NIN, NOUT )
!
!       .. Scalar Arguments ..
!       LOGICAL            TSTERR
!       INTEGER            NIN, NOUT
!       DOUBLE PRECISION   THRESH
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DCHKEC tests eigen- condition estimation routines
!>        DLALN2, DLASY2, DLANV2, DLAQTR, DLAEXC,
!>        DTRSYL, DTREXC, DTRSNA, DTRSEN, DTGEXC
!>
!> In all cases, the routine runs through a fixed set of numerical
!> examples, subjects them to various tests, and compares the test
!> results to a threshold THRESH. In addition, DTREXC, DTRSNA and DTRSEN
!> are tested by reading in precomputed examples from a file (on input
!> unit NIN).  Output is written to output unit NOUT.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] THRESH
!> \verbatim
!>          THRESH is DOUBLE PRECISION
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
!> \ingroup double_eig
!
!  =====================================================================
      SUBROUTINE DCHKEC(Thresh,Tsterr,Nin,Nout)
      IMPLICIT NONE
!*--DCHKEC80
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      LOGICAL Tsterr
      INTEGER Nin , Nout
      DOUBLE PRECISION Thresh
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
      LOGICAL ok
      CHARACTER*3 path
      INTEGER klaexc , klaln2 , klanv2 , klaqtr , klasy2 , ktrexc ,     &
     &        ktrsen , ktrsna , ktrsyl , llaexc , llaln2 , llanv2 ,     &
     &        llaqtr , llasy2 , ltrexc , ltrsyl , nlanv2 , nlaqtr ,     &
     &        nlasy2 , ntests , ntrsyl , ktgexc , ntgexc , ltgexc
      DOUBLE PRECISION eps , rlaexc , rlaln2 , rlanv2 , rlaqtr ,        &
     &                 rlasy2 , rtrexc , rtrsyl , sfmin , rtgexc
!     ..
!     .. Local Arrays ..
      INTEGER ltrsen(3) , ltrsna(3) , nlaexc(2) , nlaln2(2) , ntrexc(3) &
     &        , ntrsen(3) , ntrsna(3)
      DOUBLE PRECISION rtrsen(3) , rtrsna(3)
!     ..
!     .. External Subroutines ..
      EXTERNAL DERREC , DGET31 , DGET32 , DGET33 , DGET34 , DGET35 ,    &
     &         DGET36 , DGET37 , DGET38 , DGET39 , DGET40
!     ..
!     .. External Functions ..
      DOUBLE PRECISION DLAMCH
      EXTERNAL DLAMCH
!     ..
!     .. Executable Statements ..
!
      path(1:1) = 'Double precision'
      path(2:3) = 'EC'
      eps = DLAMCH('P')
      sfmin = DLAMCH('S')
!
!     Print header information
!
      WRITE (Nout,FMT=99011)
      WRITE (Nout,FMT=99012) eps , sfmin
      WRITE (Nout,FMT=99013) Thresh
!
!     Test error exits if TSTERR is .TRUE.
!
      IF ( Tsterr ) CALL DERREC(path,Nout)
!
      ok = .TRUE.
      CALL DGET31(rlaln2,llaln2,nlaln2,klaln2)
      IF ( rlaln2>Thresh .OR. nlaln2(1)/=0 ) THEN
         ok = .FALSE.
         WRITE (Nout,FMT=99001) rlaln2 , llaln2 , nlaln2 , klaln2
      ENDIF
!
      CALL DGET32(rlasy2,llasy2,nlasy2,klasy2)
      IF ( rlasy2>Thresh ) THEN
         ok = .FALSE.
         WRITE (Nout,FMT=99002) rlasy2 , llasy2 , nlasy2 , klasy2
      ENDIF
!
      CALL DGET33(rlanv2,llanv2,nlanv2,klanv2)
      IF ( rlanv2>Thresh .OR. nlanv2/=0 ) THEN
         ok = .FALSE.
         WRITE (Nout,FMT=99003) rlanv2 , llanv2 , nlanv2 , klanv2
      ENDIF
!
      CALL DGET34(rlaexc,llaexc,nlaexc,klaexc)
      IF ( rlaexc>Thresh .OR. nlaexc(2)/=0 ) THEN
         ok = .FALSE.
         WRITE (Nout,FMT=99004) rlaexc , llaexc , nlaexc , klaexc
      ENDIF
!
      CALL DGET35(rtrsyl,ltrsyl,ntrsyl,ktrsyl)
      IF ( rtrsyl>Thresh ) THEN
         ok = .FALSE.
         WRITE (Nout,FMT=99005) rtrsyl , ltrsyl , ntrsyl , ktrsyl
      ENDIF
!
      CALL DGET36(rtrexc,ltrexc,ntrexc,ktrexc,Nin)
      IF ( rtrexc>Thresh .OR. ntrexc(3)>0 ) THEN
         ok = .FALSE.
         WRITE (Nout,FMT=99006) rtrexc , ltrexc , ntrexc , ktrexc
      ENDIF
!
      CALL DGET37(rtrsna,ltrsna,ntrsna,ktrsna,Nin)
      IF ( rtrsna(1)>Thresh .OR. rtrsna(2)>Thresh .OR. ntrsna(1)/=0 .OR.&
     &     ntrsna(2)/=0 .OR. ntrsna(3)/=0 ) THEN
         ok = .FALSE.
         WRITE (Nout,FMT=99007) rtrsna , ltrsna , ntrsna , ktrsna
      ENDIF
!
      CALL DGET38(rtrsen,ltrsen,ntrsen,ktrsen,Nin)
      IF ( rtrsen(1)>Thresh .OR. rtrsen(2)>Thresh .OR. ntrsen(1)/=0 .OR.&
     &     ntrsen(2)/=0 .OR. ntrsen(3)/=0 ) THEN
         ok = .FALSE.
         WRITE (Nout,FMT=99008) rtrsen , ltrsen , ntrsen , ktrsen
      ENDIF
!
      CALL DGET39(rlaqtr,llaqtr,nlaqtr,klaqtr)
      IF ( rlaqtr>Thresh ) THEN
         ok = .FALSE.
         WRITE (Nout,FMT=99009) rlaqtr , llaqtr , nlaqtr , klaqtr
      ENDIF
!
      CALL DGET40(rtgexc,ltgexc,ntgexc,ktgexc,Nin)
      IF ( rtgexc>Thresh ) THEN
         ok = .FALSE.
         WRITE (Nout,FMT=99014) rtgexc , ltgexc , ntgexc , ktgexc
      ENDIF
!
      ntests = klaln2 + klasy2 + klanv2 + klaexc + ktrsyl + ktrexc +    &
     &         ktrsna + ktrsen + klaqtr + ktgexc
      IF ( ok ) WRITE (Nout,FMT=99010) path , ntests
!
      RETURN
99001 FORMAT (' Error in DLALN2: RMAX =',D12.3,/' LMAX = ',I8,' N',     &
     &        'INFO=',2I8,' KNT=',I8)
99002 FORMAT (' Error in DLASY2: RMAX =',D12.3,/' LMAX = ',I8,' N',     &
     &        'INFO=',I8,' KNT=',I8)
99003 FORMAT (' Error in DLANV2: RMAX =',D12.3,/' LMAX = ',I8,' N',     &
     &        'INFO=',I8,' KNT=',I8)
99004 FORMAT (' Error in DLAEXC: RMAX =',D12.3,/' LMAX = ',I8,' N',     &
     &        'INFO=',2I8,' KNT=',I8)
99005 FORMAT (' Error in DTRSYL: RMAX =',D12.3,/' LMAX = ',I8,' N',     &
     &        'INFO=',I8,' KNT=',I8)
99006 FORMAT (' Error in DTREXC: RMAX =',D12.3,/' LMAX = ',I8,' N',     &
     &        'INFO=',3I8,' KNT=',I8)
99007 FORMAT (' Error in DTRSNA: RMAX =',3D12.3,/' LMAX = ',3I8,        &
     &        ' NINFO=',3I8,' KNT=',I8)
99008 FORMAT (' Error in DTRSEN: RMAX =',3D12.3,/' LMAX = ',3I8,        &
     &        ' NINFO=',3I8,' KNT=',I8)
99009 FORMAT (' Error in DLAQTR: RMAX =',D12.3,/' LMAX = ',I8,' N',     &
     &        'INFO=',I8,' KNT=',I8)
99010 FORMAT (/1X,'All tests for ',A3,' routines passed the thresh',    &
     &        'old ( ',I6,' tests run)')
99011 FORMAT (' Tests of the Nonsymmetric eigenproblem condition estim',&
     &        'ation routines',/' DLALN2, DLASY2, DLANV2, DLAEXC, DTRS',&
     &        'YL, DTREXC, DTRSNA, DTRSEN, DLAQTR, DTGEXC',/)
99012 FORMAT (' Relative machine precision (EPS) = ',D16.6,/' Safe ',   &
     &        'minimum (SFMIN)             = ',D16.6,/)
99013 FORMAT (' Routines pass computational tests if test ratio is les',&
     &        's than',F8.2,//)
99014 FORMAT (' Error in DTGEXC: RMAX =',D12.3,/' LMAX = ',I8,' N',     &
     &        'INFO=',I8,' KNT=',I8)
!
!     End of DCHKEC
!
      END SUBROUTINE DCHKEC
