!*==zchkab.f90  processed by SPAG 7.51RB at 17:34 on  4 Mar 2022
!> \brief \b ZCHKAB
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       PROGRAM ZCHKAB
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZCHKAB is the test program for the COMPLEX*16 LAPACK
!> ZCGESV/ZCPOSV routine
!>
!> The program must be driven by a short data file. The first 5 records
!> specify problem dimensions and program options using list-directed
!> input. The remaining lines specify the LAPACK test paths and the
!> number of matrix types to use in testing.  An annotated example of a
!> data file can be obtained by deleting the first 3 characters from the
!> following 9 lines:
!> Data file for testing COMPLEX*16 LAPACK ZCGESV
!> 7                      Number of values of M
!> 0 1 2 3 5 10 16        Values of M (row dimension)
!> 1                      Number of values of NRHS
!> 2                      Values of NRHS (number of right hand sides)
!> 20.0                   Threshold value of test ratio
!> T                      Put T to test the LAPACK routine
!> T                      Put T to test the error exits
!> DGE    11              List types on next line if 0 < NTYPES < 11
!> DPO    9               List types on next line if 0 < NTYPES <  9
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \verbatim
!>  NMAX    INTEGER
!>          The maximum allowable value for N
!>
!>  MAXIN   INTEGER
!>          The number of different values that can be used for each of
!>          M, N, NRHS, NB, and NX
!>
!>  MAXRHS  INTEGER
!>          The maximum number of right hand sides
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
      PROGRAM ZCHKAB
      use M_tst_lin, only : ALAREQ , ZDRVAB , ZDRVAC , ZERRAB , ZERRAC
      IMPLICIT NONE
!*--ZCHKAB77
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     April 2012
!
!  =====================================================================
!
!     .. Parameters ..
      INTEGER NMAX
      PARAMETER (NMAX=132)
      INTEGER MAXIN
      PARAMETER (MAXIN=12)
      INTEGER MAXRHS
      PARAMETER (MAXRHS=16)
      INTEGER MATMAX
      PARAMETER (MATMAX=30)
      INTEGER NIN , NOUT
      PARAMETER (NIN=5,NOUT=6)
      INTEGER LDAMAX
      PARAMETER (LDAMAX=NMAX)
!     ..
!     .. Local Scalars ..
      LOGICAL fatal , tstdrv , tsterr
      CHARACTER c1
      CHARACTER*2 c2
      CHARACTER*3 path
      CHARACTER*10 intstr
      CHARACTER*72 aline
      INTEGER i , ic , k , lda , nm , nmats , nns , nrhs , ntypes ,     &
     &        vers_major , vers_minor , vers_patch
      DOUBLE PRECISION eps , s1 , s2 , thresh
      REAL seps
!     ..
!     .. Local Arrays ..
      LOGICAL dotype(MATMAX)
      INTEGER iwork(NMAX) , mval(MAXIN) , nsval(MAXIN)
      DOUBLE PRECISION rwork(NMAX)
      COMPLEX*16 a(LDAMAX*NMAX,2) , b(NMAX*MAXRHS,2) ,                  &
     &           work(NMAX*MAXRHS*2)
      COMPLEX swork(NMAX*(NMAX+MAXRHS))
!     ..
!     .. External Functions ..
      DOUBLE PRECISION DLAMCH , DSECND
      LOGICAL LSAME , LSAMEN
      REAL SLAMCH
      EXTERNAL DLAMCH , DSECND , LSAME , LSAMEN , SLAMCH
!     ..
!     .. External Subroutines ..
      EXTERNAL ILAVER
!     ..
!     .. Scalars in Common ..
      LOGICAL LERr , OK
      CHARACTER*32 SRNamt
      INTEGER INFot , NUNit
!     ..
!     .. Common blocks ..
      COMMON /INFOC / INFot , NUNit , OK , LERr
      COMMON /SRNAMC/ SRNamt
!
!     .. Data statements ..
      DATA intstr/'0123456789'/
!     ..
!     .. Executable Statements ..
!
      s1 = DSECND()
      lda = NMAX
      fatal = .FALSE.
!
!     Read a dummy line.
!
      READ (NIN,FMT=*)
!
!     Report values of parameters.
!
      CALL ILAVER(vers_major,vers_minor,vers_patch)
      WRITE (NOUT,FMT=99006) vers_major , vers_minor , vers_patch
!
!     Read the values of M
!
      READ (NIN,FMT=*) nm
      IF ( nm<1 ) THEN
         WRITE (NOUT,FMT=99004) ' NM ' , nm , 1
         nm = 0
         fatal = .TRUE.
      ELSEIF ( nm>MAXIN ) THEN
         WRITE (NOUT,FMT=99005) ' NM ' , nm , MAXIN
         nm = 0
         fatal = .TRUE.
      ENDIF
      READ (NIN,FMT=*) (mval(i),i=1,nm)
      DO i = 1 , nm
         IF ( mval(i)<0 ) THEN
            WRITE (NOUT,FMT=99004) ' M  ' , mval(i) , 0
            fatal = .TRUE.
         ELSEIF ( mval(i)>NMAX ) THEN
            WRITE (NOUT,FMT=99005) ' M  ' , mval(i) , NMAX
            fatal = .TRUE.
         ENDIF
      ENDDO
      IF ( nm>0 ) WRITE (NOUT,FMT=99007) 'M   ' , (mval(i),i=1,nm)
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
!     Read the threshold value for the test ratios.
!
      READ (NIN,FMT=*) thresh
      WRITE (NOUT,FMT=99008) thresh
!
!     Read the flag that indicates whether to test the driver routine.
!
      READ (NIN,FMT=*) tstdrv
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
      seps = SLAMCH('Underflow threshold')
      WRITE (NOUT,FMT=99009) '(single precision) underflow' , seps
      seps = SLAMCH('Overflow threshold')
      WRITE (NOUT,FMT=99009) '(single precision) overflow ' , seps
      seps = SLAMCH('Epsilon')
      WRITE (NOUT,FMT=99009) '(single precision) precision' , seps
      WRITE (NOUT,FMT=*)
!
      eps = DLAMCH('Underflow threshold')
      WRITE (NOUT,FMT=99009) '(double precision) underflow' , eps
      eps = DLAMCH('Overflow threshold')
      WRITE (NOUT,FMT=99009) '(double precision) overflow ' , eps
      eps = DLAMCH('Epsilon')
      WRITE (NOUT,FMT=99009) '(double precision) precision' , eps
      WRITE (NOUT,FMT=*)
!
!
!     Read a test path and the number of matrix types to use.
!
 100  READ (NIN,FMT='(A72)',END=500) aline
      path = aline(1:3)
      nmats = MATMAX
      i = 3
      DO
         i = i + 1
         IF ( i>72 ) THEN
            nmats = MATMAX
            GOTO 400
         ENDIF
         IF ( aline(i:i)/=' ' ) THEN
            nmats = 0
            EXIT
         ENDIF
      ENDDO
 200  c1 = aline(i:i)
      DO k = 1 , 10
         IF ( c1==intstr(k:k) ) THEN
            ic = k - 1
            GOTO 300
         ENDIF
      ENDDO
      GOTO 400
 300  nmats = nmats*10 + ic
      i = i + 1
      IF ( i<=72 ) GOTO 200
 400  c1 = path(1:1)
      c2 = path(2:3)
      nrhs = nsval(1)
      nrhs = nsval(1)
!
!     Check first character for correct precision.
!
      IF ( .NOT.LSAME(c1,'Zomplex precision') ) THEN
         WRITE (NOUT,FMT=99010) path
!
      ELSEIF ( nmats<=0 ) THEN
!
!        Check for a positive number of tests requested.
!
         WRITE (NOUT,FMT=99010) 'ZCGESV'
         GOTO 500
!
      ELSEIF ( LSAMEN(2,c2,'GE') ) THEN
!
!        GE:  general matrices
!
         ntypes = 11
         CALL ALAREQ('ZGE',nmats,dotype,ntypes,NIN,NOUT)
!
!        Test the error exits
!
         IF ( tsterr ) CALL ZERRAB(NOUT)
!
         IF ( tstdrv ) THEN
            CALL ZDRVAB(dotype,nm,mval,nns,nsval,thresh,lda,a(1,1),     &
     &                  a(1,2),b(1,1),b(1,2),work,rwork,swork,iwork,    &
     &                  NOUT)
         ELSE
            WRITE (NOUT,FMT=99011) 'ZCGESV'
         ENDIF
!
      ELSEIF ( LSAMEN(2,c2,'PO') ) THEN
!
!        PO:  positive definite matrices
!
         ntypes = 9
         CALL ALAREQ('DPO',nmats,dotype,ntypes,NIN,NOUT)
!
         IF ( tsterr ) CALL ZERRAC(NOUT)
!
!
         IF ( tstdrv ) THEN
            CALL ZDRVAC(dotype,nm,mval,nns,nsval,thresh,lda,a(1,1),     &
     &                  a(1,2),b(1,1),b(1,2),work,rwork,swork,NOUT)
         ELSE
            WRITE (NOUT,FMT=99011) 'ZCPOSV'
         ENDIF
!
!
      ENDIF
!
!     Go back to get another input line.
!
      GOTO 100
!
!     Branch to this line when the last record is read.
!
 500  CLOSE (NIN)
      s2 = DSECND()
      WRITE (NOUT,FMT=99002)
      WRITE (NOUT,FMT=99003) s2 - s1
!
99001 FORMAT (/' Execution not attempted due to input errors')
99002 FORMAT (/' End of tests')
99003 FORMAT (' Total time used = ',F12.2,' seconds',/)
99004 FORMAT (' Invalid input value: ',A4,'=',I6,'; must be >=',I6)
99005 FORMAT (' Invalid input value: ',A4,'=',I6,'; must be <=',I6)
99006 FORMAT (' Tests of the COMPLEX*16 LAPACK ZCGESV/ZCPOSV routines ',&
     &        /' LAPACK VERSION ',I1,'.',I1,'.',I1,                     &
     &        //' The following parameter values will be used:')
99007 FORMAT (4X,A4,':  ',10I6,/11X,10I6)
99008 FORMAT (/' Routines pass computational tests if test ratio is ',  &
     &        'less than',F8.2,/)
99009 FORMAT (' Relative machine ',A,' is taken to be',D16.6)
99010 FORMAT (/1X,A6,' routines were not tested')
99011 FORMAT (/1X,A6,' driver routines were not tested')
!
!     End of ZCHKAB
!
      END PROGRAM ZCHKAB
