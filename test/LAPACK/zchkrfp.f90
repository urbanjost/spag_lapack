!*==zchkrfp.f90 processed by SPAG 7.51RB at 17:35 on  4 Mar 2022
!> \brief \b ZCHKRFP
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       PROGRAM ZCHKRFP
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZCHKRFP is the main test program for the COMPLEX*16 linear equation
!> routines with RFP storage format
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
!> \ingroup complex16_lin
!
!  =====================================================================
      PROGRAM ZCHKRFP
      use M_tst__lin, only : ZDRVRFP , ZDRVRF1 , ZDRVRF2 , ZDRVRF3 , ZDRVRF4, ZERRRFP
      IMPLICIT NONE
!*--ZCHKRFP63
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
      DOUBLE PRECISION eps , s1 , s2 , thresh

!     ..
!     .. Local Arrays ..
      INTEGER nval(MAXIN) , nsval(MAXIN) , ntval(NTYPES)
      COMPLEX*16 worka(NMAX,NMAX)
      COMPLEX*16 workasav(NMAX,NMAX)
      COMPLEX*16 workb(NMAX,MAXRHS)
      COMPLEX*16 workxact(NMAX,MAXRHS)
      COMPLEX*16 workbsav(NMAX,MAXRHS)
      COMPLEX*16 workx(NMAX,MAXRHS)
      COMPLEX*16 workafac(NMAX,NMAX)
      COMPLEX*16 workainv(NMAX,NMAX)
      COMPLEX*16 workarf((NMAX*(NMAX+1))/2)
      COMPLEX*16 workap((NMAX*(NMAX+1))/2)
      COMPLEX*16 workarfinv((NMAX*(NMAX+1))/2)
      COMPLEX*16 z_work_zlatms(3*NMAX)
      COMPLEX*16 z_work_zpot02(NMAX,MAXRHS)
      COMPLEX*16 z_work_zpot03(NMAX,NMAX)
      DOUBLE PRECISION d_work_zlatms(NMAX)
      DOUBLE PRECISION d_work_zlanhe(NMAX)
      DOUBLE PRECISION d_work_zpot01(NMAX)
      DOUBLE PRECISION d_work_zpot02(NMAX)
      DOUBLE PRECISION d_work_zpot03(NMAX)
!     ..
!     .. External Functions ..
      DOUBLE PRECISION DLAMCH , DSECND
      EXTERNAL DLAMCH , DSECND
!     ..
!     .. External Subroutines ..
      EXTERNAL ILAVER
!     ..
!     .. Executable Statements ..
!
      s1 = DSECND()
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
      eps = DLAMCH('Underflow threshold')
      WRITE (NOUT,FMT=99009) 'underflow' , eps
      eps = DLAMCH('Overflow threshold')
      WRITE (NOUT,FMT=99009) 'overflow ' , eps
      eps = DLAMCH('Epsilon')
      WRITE (NOUT,FMT=99009) 'precision' , eps
      WRITE (NOUT,FMT=*)
!
!     Test the error exit of:
!
      IF ( tsterr ) CALL ZERRRFP(NOUT)
!
!    Test the routines: zpftrf, zpftri, zpftrs (as in ZDRVPO).
!    This also tests the routines: ztfsm, ztftri, ztfttr, ztrttf.
!
      CALL ZDRVRFP(NOUT,nn,nval,nns,nsval,nnt,ntval,thresh,worka,       &
     &             workasav,workafac,workainv,workb,workbsav,workxact,  &
     &             workx,workarf,workarfinv,z_work_zlatms,z_work_zpot02,&
     &             z_work_zpot03,d_work_zlatms,d_work_zlanhe,           &
     &             d_work_zpot01,d_work_zpot02,d_work_zpot03)
!
!    Test the routine: zlanhf
!
      CALL ZDRVRF1(NOUT,nn,nval,thresh,worka,NMAX,workarf,d_work_zlanhe)
!
!    Test the conversion routines:
!       zhfttp, ztpthf, ztfttr, ztrttf, ztrttp and ztpttr.
!
      CALL ZDRVRF2(NOUT,nn,nval,worka,NMAX,workarf,workap,workasav)
!
!    Test the routine: ztfsm
!
      CALL ZDRVRF3(NOUT,nn,nval,thresh,worka,NMAX,workarf,workainv,     &
     &             workafac,d_work_zlanhe,z_work_zpot03,z_work_zpot02)

!
!    Test the routine: zhfrk
!
      CALL ZDRVRF4(NOUT,nn,nval,thresh,worka,workafac,NMAX,workarf,     &
     &             workainv,NMAX,d_work_zlanhe)
!
      CLOSE (NIN)
      s2 = DSECND()
      WRITE (NOUT,FMT=99002)
      WRITE (NOUT,FMT=99003) s2 - s1
!
99001 FORMAT (/' Execution not attempted due to input errors')
99002 FORMAT (/' End of tests')
99003 FORMAT (' Total time used = ',F12.2,' seconds',/)
99004 FORMAT (' !! Invalid input value: ',A4,'=',I6,'; must be >=',I6)
99005 FORMAT (' !! Invalid input value: ',A4,'=',I6,'; must be <=',I6)
99006 FORMAT (/' Tests of the COMPLEX*16 LAPACK RFP routines ',         &
     &        /' LAPACK VERSION ',I1,'.',I1,'.',I1,                     &
     &        //' The following parameter values will be used:')
99007 FORMAT (4X,A4,':  ',10I6,/11X,10I6)
99008 FORMAT (/' Routines pass computational tests if test ratio is ',  &
     &        'less than',F8.2,/)
99009 FORMAT (' Relative machine ',A,' is taken to be',D16.6)
!
!     End of ZCHKRFP
!
      END PROGRAM ZCHKRFP
