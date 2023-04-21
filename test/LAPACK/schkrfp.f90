!*==schkrfp.f90 processed by SPAG 7.51RB at 17:34 on  4 Mar 2022
!> \brief \b SCHKRFP
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       PROGRAM SCHKRFP
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SCHKRFP is the main test program for the REAL linear
!> equation routines with RFP storage format
!>
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \verbatim
!>  MAXIN   INTEGER
!>          The number of different values that can be used for each of
!>          M, N, or NB
!>
!>  MAXRHS  INTEGER
!>          The maximum number of right hand sides
!>
!>  NTYPES  INTEGER
!>
!>  NMAX    INTEGER
!>          The maximum allowable value for N.
!>
!>  NIN     INTEGER
!>          The unit number for input
!>
!>  NOUT    INTEGER
!>          The unit number for output
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
!> \date April 2012
!
!> \ingroup single_lin
!
!  =====================================================================
      PROGRAM SCHKRFP
use M_tst_lin, only: aladhd, alahd, cgennd, cqpt01, cqrt11, cqrt12, cqrt14, cqrt17
use M_tst_lin, only: crzt01, crzt02, dgennd, dget06, dqpt01, dqrt11, dqrt12, dqrt14
use M_tst_lin, only: dqrt17, drzt01, drzt02, alaerh, alasvm, sdrvrf2, sdrvrf3, sdrvrf4
use M_tst_lin, only: slarhs, spot01, spot02, spot03, sdrvrf1, sdrvrfp, serrrfp, sgennd
use M_tst_lin, only: sqpt01, sqrt11, sqrt12, sqrt14, sqrt17, srzt01, srzt02, zgennd
use M_tst_lin, only: zqpt01, zqrt11, zqrt12, zqrt14, zqrt17, zrzt01, zrzt02
use M_tst_matgen, only : clagge, claghe, clagsy, clakf2, clarge, clarnd, claror, clarot
use M_tst_matgen, only : clatm1, clatm2, clatm3, clatms, clatmt, dlagge, dlagsy, dlakf2
use M_tst_matgen, only : dlaran, dlarge, dlarnd, dlaror, dlarot, dlatm1, dlatm2, dlatm3
use M_tst_matgen, only : dlatm7, dlatms, dlatmt, slagge, slagsy, slakf2, slaran, slarge
use M_tst_matgen, only : slarnd, slaror, slarot, slatm1, slatm2, slatm3, slatm7, slatms
use M_tst_matgen, only : slatmt, zlagge, zlaghe, zlagsy, zlakf2, zlarge, zlarnd, zlaror
use M_tst_matgen, only : zlarot, zlatm1, zlatm2, zlatm3, zlatms, zlatmt
      IMPLICIT NONE
!*--SCHKRFP63
!
!  -- LAPACK test routine (version 3.7.1) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     April 2012
!
!  =====================================================================
!
!     .. Parameters ..
      INTEGER MAXIN
      PARAMETER (MAXIN=12)
      INTEGER NMAX
      PARAMETER (NMAX=50)
      INTEGER MAXRHS
      PARAMETER (MAXRHS=16)
      INTEGER NTYPES
      PARAMETER (NTYPES=9)
      INTEGER NIN , NOUT
      PARAMETER (NIN=5,NOUT=6)
!     ..
!     .. Local Scalars ..
      LOGICAL fatal , tsterr
      INTEGER vers_major , vers_minor , vers_patch
      INTEGER i , nn , nns , nnt
      REAL eps , s1 , s2 , thresh
!     ..
!     .. Local Arrays ..
      INTEGER nval(MAXIN) , nsval(MAXIN) , ntval(NTYPES)
      REAL worka(NMAX,NMAX)
      REAL workasav(NMAX,NMAX)
      REAL workb(NMAX,MAXRHS)
      REAL workxact(NMAX,MAXRHS)
      REAL workbsav(NMAX,MAXRHS)
      REAL workx(NMAX,MAXRHS)
      REAL workafac(NMAX,NMAX)
      REAL workainv(NMAX,NMAX)
      REAL workarf((NMAX*(NMAX+1))/2)
      REAL workap((NMAX*(NMAX+1))/2)
      REAL workarfinv((NMAX*(NMAX+1))/2)
      REAL s_work_slatms(3*NMAX)
      REAL s_work_spot01(NMAX)
      REAL s_temp_spot02(NMAX,MAXRHS)
      REAL s_temp_spot03(NMAX,NMAX)
      REAL s_work_slansy(NMAX)
      REAL s_work_spot02(NMAX)
      REAL s_work_spot03(NMAX)
!     ..
!     .. External Functions ..
      REAL SLAMCH , SECOND
      EXTERNAL SLAMCH , SECOND
!     ..
!     .. External Subroutines ..
      EXTERNAL ILAVER
!     ..
!     .. Executable Statements ..
!
      s1 = SECOND()
      fatal = .FALSE.
!
!     Read a dummy line.
!
      READ (NIN,FMT=*)
!
!     Report LAPACK version tag (e.g. LAPACK-3.2.0)
!
      CALL ILAVER(vers_major,vers_minor,vers_patch)
      WRITE (NOUT,FMT=99006) vers_major , vers_minor , vers_patch
!
!     Read the values of N
!
      READ (NIN,FMT=*) nn
      IF ( nn<1 ) THEN
         WRITE (NOUT,FMT=99004) ' NN ' , nn , 1
         nn = 0
         fatal = .TRUE.
      ELSEIF ( nn>MAXIN ) THEN
         WRITE (NOUT,FMT=99005) ' NN ' , nn , MAXIN
         nn = 0
         fatal = .TRUE.
      ENDIF
      READ (NIN,FMT=*) (nval(i),i=1,nn)
      DO i = 1 , nn
         IF ( nval(i)<0 ) THEN
            WRITE (NOUT,FMT=99004) ' M  ' , nval(i) , 0
            fatal = .TRUE.
         ELSEIF ( nval(i)>NMAX ) THEN
            WRITE (NOUT,FMT=99005) ' M  ' , nval(i) , NMAX
            fatal = .TRUE.
         ENDIF
      ENDDO
      IF ( nn>0 ) WRITE (NOUT,FMT=99007) 'N   ' , (nval(i),i=1,nn)
!
!     Read the values of NRHS
!
      READ (NIN,FMT=*) nns
      IF ( nns<1 ) THEN
         WRITE (NOUT,FMT=99004) ' NNS' , nns , 1
         nns = 0
         fatal = .TRUE.
      ELSEIF ( nns>MAXIN ) THEN
         WRITE (NOUT,FMT=99005) ' NNS' , nns , MAXIN
         nns = 0
         fatal = .TRUE.
      ENDIF
      READ (NIN,FMT=*) (nsval(i),i=1,nns)
      DO i = 1 , nns
         IF ( nsval(i)<0 ) THEN
            WRITE (NOUT,FMT=99004) 'NRHS' , nsval(i) , 0
            fatal = .TRUE.
         ELSEIF ( nsval(i)>MAXRHS ) THEN
            WRITE (NOUT,FMT=99005) 'NRHS' , nsval(i) , MAXRHS
            fatal = .TRUE.
         ENDIF
      ENDDO
      IF ( nns>0 ) WRITE (NOUT,FMT=99007) 'NRHS' , (nsval(i),i=1,nns)
!
!     Read the matrix types
!
      READ (NIN,FMT=*) nnt
      IF ( nnt<1 ) THEN
         WRITE (NOUT,FMT=99004) ' NMA' , nnt , 1
         nnt = 0
         fatal = .TRUE.
      ELSEIF ( nnt>NTYPES ) THEN
         WRITE (NOUT,FMT=99005) ' NMA' , nnt , NTYPES
         nnt = 0
         fatal = .TRUE.
      ENDIF
      READ (NIN,FMT=*) (ntval(i),i=1,nnt)
      DO i = 1 , nnt
         IF ( ntval(i)<0 ) THEN
            WRITE (NOUT,FMT=99004) 'TYPE' , ntval(i) , 0
            fatal = .TRUE.
         ELSEIF ( ntval(i)>NTYPES ) THEN
            WRITE (NOUT,FMT=99005) 'TYPE' , ntval(i) , NTYPES
            fatal = .TRUE.
         ENDIF
      ENDDO
      IF ( nnt>0 ) WRITE (NOUT,FMT=99007) 'TYPE' , (ntval(i),i=1,nnt)
!
!     Read the threshold value for the test ratios.
!
      READ (NIN,FMT=*) thresh
      WRITE (NOUT,FMT=99008) thresh
!
!     Read the flag that indicates whether to test the error exits.
!
      READ (NIN,FMT=*) tsterr
!
      IF ( fatal ) THEN
         WRITE (NOUT,FMT=99001)
         STOP
      ENDIF
!
!     Calculate and print the machine dependent constants.
!
      eps = SLAMCH('Underflow threshold')
      WRITE (NOUT,FMT=99009) 'underflow' , eps
      eps = SLAMCH('Overflow threshold')
      WRITE (NOUT,FMT=99009) 'overflow ' , eps
      eps = SLAMCH('Epsilon')
      WRITE (NOUT,FMT=99009) 'precision' , eps
      WRITE (NOUT,FMT=*)
!
!     Test the error exit of:
!
      IF ( tsterr ) CALL SERRRFP(NOUT)
!
!     Test the routines: spftrf, spftri, spftrs (as in SDRVPO).
!     This also tests the routines: stfsm, stftri, stfttr, strttf.
!
      CALL SDRVRFP(NOUT,nn,nval,nns,nsval,nnt,ntval,thresh,worka,       &
     &             workasav,workafac,workainv,workb,workbsav,workxact,  &
     &             workx,workarf,workarfinv,s_work_slatms,s_work_spot01,&
     &             s_temp_spot02,s_temp_spot03,s_work_slansy,           &
     &             s_work_spot02,s_work_spot03)
!
!     Test the routine: slansf
!
      CALL SDRVRF1(NOUT,nn,nval,thresh,worka,NMAX,workarf,s_work_slansy)
!
!     Test the conversion routines:
!       stfttp, stpttf, stfttr, strttf, strttp and stpttr.
!
      CALL SDRVRF2(NOUT,nn,nval,worka,NMAX,workarf,workap,workasav)
!
!     Test the routine: stfsm
!
      CALL SDRVRF3(NOUT,nn,nval,thresh,worka,NMAX,workarf,workainv,     &
     &             workafac,s_work_slansy,s_work_spot03,s_work_spot01)
!
!
!     Test the routine: ssfrk
!
      CALL SDRVRF4(NOUT,nn,nval,thresh,worka,workafac,NMAX,workarf,     &
     &             workainv,NMAX,s_work_slansy)
!
      CLOSE (NIN)
      s2 = SECOND()
      WRITE (NOUT,FMT=99002)
      WRITE (NOUT,FMT=99003) s2 - s1
!
99001 FORMAT (/' Execution not attempted due to input errors')
99002 FORMAT (/' End of tests')
99003 FORMAT (' Total time used = ',F12.2,' seconds',/)
99004 FORMAT (' !! Invalid input value: ',A4,'=',I6,'; must be >=',I6)
99005 FORMAT (' !! Invalid input value: ',A4,'=',I6,'; must be <=',I6)
99006 FORMAT (/' Tests of the REAL LAPACK RFP routines ',               &
     &        /' LAPACK VERSION ',I1,'.',I1,'.',I1,                     &
     &        //' The following parameter values will be used:')
99007 FORMAT (4X,A4,':  ',10I6,/11X,10I6)
99008 FORMAT (/' Routines pass computational tests if test ratio is ',  &
     &        'less than',F8.2,/)
99009 FORMAT (' Relative machine ',A,' is taken to be',D16.6)
!
!     End of SCHKRFP
!
      END PROGRAM SCHKRFP
