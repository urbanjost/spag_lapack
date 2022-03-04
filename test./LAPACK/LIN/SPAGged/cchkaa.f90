!*==cchkaa.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b CCHKAA
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       PROGRAM CCHKAA
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CCHKAA is the main test program for the COMPLEX linear equation
!> routines.
!>
!> The program must be driven by a short data file. The first 15 records
!> (not including the first comment  line) specify problem dimensions
!> and program options using list-directed input. The remaining lines
!> specify the LAPACK test paths and the number of matrix types to use
!> in testing.  An annotated example of a data file can be obtained by
!> deleting the first 3 characters from the following 42 lines:
!> Data file for testing COMPLEX LAPACK linear equation routines
!> 7                      Number of values of M
!> 0 1 2 3 5 10 16        Values of M (row dimension)
!> 7                      Number of values of N
!> 0 1 2 3 5 10 16        Values of N (column dimension)
!> 1                      Number of values of NRHS
!> 2                      Values of NRHS (number of right hand sides)
!> 5                      Number of values of NB
!> 1 3 3 3 20             Values of NB (the blocksize)
!> 1 0 5 9 1              Values of NX (crossover point)
!> 3                      Number of values of RANK
!> 30 50 90               Values of rank (as a % of N)
!> 30.0                   Threshold value of test ratio
!> T                      Put T to test the LAPACK routines
!> T                      Put T to test the driver routines
!> T                      Put T to test the error exits
!> CGE   11               List types on next line if 0 < NTYPES < 11
!> CGB    8               List types on next line if 0 < NTYPES <  8
!> CGT   12               List types on next line if 0 < NTYPES < 12
!> CPO    9               List types on next line if 0 < NTYPES <  9
!> CPO    9               List types on next line if 0 < NTYPES <  9
!> CPP    9               List types on next line if 0 < NTYPES <  9
!> CPB    8               List types on next line if 0 < NTYPES <  8
!> CPT   12               List types on next line if 0 < NTYPES < 12
!> CHE   10               List types on next line if 0 < NTYPES < 10
!> CHR   10               List types on next line if 0 < NTYPES < 10
!> CHK   10               List types on next line if 0 < NTYPES < 10
!> CHA   10               List types on next line if 0 < NTYPES < 10
!> CH2   10               List types on next line if 0 < NTYPES < 10
!> CSA   11               List types on next line if 0 < NTYPES < 10
!> CS2   11               List types on next line if 0 < NTYPES < 10
!> CHP   10               List types on next line if 0 < NTYPES < 10
!> CSY   11               List types on next line if 0 < NTYPES < 11
!> CSK   11               List types on next line if 0 < NTYPES < 11
!> CSR   11               List types on next line if 0 < NTYPES < 11
!> CSP   11               List types on next line if 0 < NTYPES < 11
!> CTR   18               List types on next line if 0 < NTYPES < 18
!> CTP   18               List types on next line if 0 < NTYPES < 18
!> CTB   17               List types on next line if 0 < NTYPES < 17
!> CQR    8               List types on next line if 0 < NTYPES <  8
!> CRQ    8               List types on next line if 0 < NTYPES <  8
!> CLQ    8               List types on next line if 0 < NTYPES <  8
!> CQL    8               List types on next line if 0 < NTYPES <  8
!> CQP    6               List types on next line if 0 < NTYPES <  6
!> CTZ    3               List types on next line if 0 < NTYPES <  3
!> CLS    6               List types on next line if 0 < NTYPES <  6
!> CEQ
!> CQT
!> CQX
!> CTS
!> CHH
!> \endverbatim
!
!  Parameters:
!  ==========
!
!> \verbatim
!>  NMAX    INTEGER
!>          The maximum allowable value for M and N.
!>
!>  MAXIN   INTEGER
!>          The number of different values that can be used for each of
!>          M, N, NRHS, NB, NX and RANK
!>
!>  MAXRHS  INTEGER
!>          The maximum number of right hand sides
!>
!>  MATMAX  INTEGER
!>          The maximum number of matrix types to use for testing
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
!> \date November 2019
!
!> \ingroup complex_lin
!
!  =====================================================================
      PROGRAM CCHKAA
      IMPLICIT NONE
!*--CCHKAA121
!
!  -- LAPACK test routine (version 3.9.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2017
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
      INTEGER KDMAX
      PARAMETER (KDMAX=NMAX+(NMAX+1)/4)
!     ..
!     .. Local Scalars ..
      LOGICAL fatal , tstchk , tstdrv , tsterr
      CHARACTER c1
      CHARACTER*2 c2
      CHARACTER*3 path
      CHARACTER*10 intstr
      CHARACTER*72 aline
      INTEGER i , ic , j , k , la , lafac , lda , nb , nm , nmats , nn ,&
     &        nnb , nnb2 , nns , nrhs , ntypes , nrank , vers_major ,   &
     &        vers_minor , vers_patch
      REAL eps , s1 , s2 , threq , thresh
!     ..
!     .. Local Arrays ..
      LOGICAL dotype(MATMAX)
      INTEGER iwork(25*NMAX) , mval(MAXIN) , nbval(MAXIN) ,             &
     &        nbval2(MAXIN) , nsval(MAXIN) , nval(MAXIN) , nxval(MAXIN) &
     &        , rankval(MAXIN) , piv(NMAX)
      REAL rwork(150*NMAX+2*MAXRHS) , s(2*NMAX)
      COMPLEX a((KDMAX+1)*NMAX,7) , b(NMAX*MAXRHS,4) , e(NMAX) ,        &
     &        work(NMAX,NMAX+MAXRHS+10)
!     ..
!     .. External Functions ..
      LOGICAL LSAME , LSAMEN
      REAL SECOND , SLAMCH
      EXTERNAL LSAME , LSAMEN , SECOND , SLAMCH
!     ..
!     .. External Subroutines ..
      EXTERNAL ALAREQ , CCHKEQ , CCHKGB , CCHKGE , CCHKGT , CCHKHE ,    &
     &         CCHKHE_ROOK , CCHKHE_RK , CCHKHE_AA , CCHKHP , CCHKLQ ,  &
     &         CCHKUNHR_COL , CCHKPB , CCHKPO , CCHKPS , CCHKPP ,       &
     &         CCHKPT , CCHKQ3 , CCHKQL , CCHKQR , CCHKRQ , CCHKSP ,    &
     &         CCHKSY , CCHKSY_ROOK , CCHKSY_RK , CCHKSY_AA , CCHKTB ,  &
     &         CCHKTP , CCHKTR , CCHKTZ , CDRVGB , CDRVGE , CDRVGT ,    &
     &         CDRVHE , CDRVHE_ROOK , CDRVHE_RK , CDRVHE_AA , CDRVHP ,  &
     &         CDRVLS , CDRVPB , CDRVPO , CDRVPP , CDRVPT , CDRVSP ,    &
     &         CDRVSY , CDRVSY_ROOK , CDRVSY_RK , CDRVSY_AA , ILAVER ,  &
     &         CCHKQRT , CCHKQRTP
!     ..
!     .. Scalars in Common ..
      LOGICAL LERr , OK
      CHARACTER*32 SRNamt
      INTEGER INFot , NUNit
!     ..
!     .. Arrays in Common ..
      INTEGER IPArms(100)
!     ..
!     .. Common blocks ..
      COMMON /CLAENV/ IPArms
      COMMON /INFOC / INFot , NUNit , OK , LERr
      COMMON /SRNAMC/ SRNamt
!     ..
!     .. Data statements ..
      DATA threq/2.0/ , intstr/'0123456789'/
!     ..
!     .. Executable Statements ..
!
      s1 = SECOND()
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
            WRITE (NOUT,FMT=99004) ' N  ' , nval(i) , 0
            fatal = .TRUE.
         ELSEIF ( nval(i)>NMAX ) THEN
            WRITE (NOUT,FMT=99005) ' N  ' , nval(i) , NMAX
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
!     Read the values of NB
!
      READ (NIN,FMT=*) nnb
      IF ( nnb<1 ) THEN
         WRITE (NOUT,FMT=99004) 'NNB ' , nnb , 1
         nnb = 0
         fatal = .TRUE.
      ELSEIF ( nnb>MAXIN ) THEN
         WRITE (NOUT,FMT=99005) 'NNB ' , nnb , MAXIN
         nnb = 0
         fatal = .TRUE.
      ENDIF
      READ (NIN,FMT=*) (nbval(i),i=1,nnb)
      DO i = 1 , nnb
         IF ( nbval(i)<0 ) THEN
            WRITE (NOUT,FMT=99004) ' NB ' , nbval(i) , 0
            fatal = .TRUE.
         ENDIF
      ENDDO
      IF ( nnb>0 ) WRITE (NOUT,FMT=99007) 'NB  ' , (nbval(i),i=1,nnb)
!
!     Set NBVAL2 to be the set of unique values of NB
!
      nnb2 = 0
      DO i = 1 , nnb
         nb = nbval(i)
         DO j = 1 , nnb2
            IF ( nb==nbval2(j) ) GOTO 100
         ENDDO
         nnb2 = nnb2 + 1
         nbval2(nnb2) = nb
 100  ENDDO
!
!     Read the values of NX
!
      READ (NIN,FMT=*) (nxval(i),i=1,nnb)
      DO i = 1 , nnb
         IF ( nxval(i)<0 ) THEN
            WRITE (NOUT,FMT=99004) ' NX ' , nxval(i) , 0
            fatal = .TRUE.
         ENDIF
      ENDDO
      IF ( nnb>0 ) WRITE (NOUT,FMT=99007) 'NX  ' , (nxval(i),i=1,nnb)
!
!     Read the values of RANKVAL
!
      READ (NIN,FMT=*) nrank
      IF ( nn<1 ) THEN
         WRITE (NOUT,FMT=99004) ' NRANK ' , nrank , 1
         nrank = 0
         fatal = .TRUE.
      ELSEIF ( nn>MAXIN ) THEN
         WRITE (NOUT,FMT=99005) ' NRANK ' , nrank , MAXIN
         nrank = 0
         fatal = .TRUE.
      ENDIF
      READ (NIN,FMT=*) (rankval(i),i=1,nrank)
      DO i = 1 , nrank
         IF ( rankval(i)<0 ) THEN
            WRITE (NOUT,FMT=99004) ' RANK  ' , rankval(i) , 0
            fatal = .TRUE.
         ELSEIF ( rankval(i)>100 ) THEN
            WRITE (NOUT,FMT=99005) ' RANK  ' , rankval(i) , 100
            fatal = .TRUE.
         ENDIF
      ENDDO
      IF ( nrank>0 ) WRITE (NOUT,FMT=99007) 'RANK % OF N' ,             &
     &                      (rankval(i),i=1,nrank)
!
!     Read the threshold value for the test ratios.
!
      READ (NIN,FMT=*) thresh
      WRITE (NOUT,FMT=99008) thresh
!
!     Read the flag that indicates whether to test the LAPACK routines.
!
      READ (NIN,FMT=*) tstchk
!
!     Read the flag that indicates whether to test the driver routines.
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
      eps = SLAMCH('Underflow threshold')
      WRITE (NOUT,FMT=99009) 'underflow' , eps
      eps = SLAMCH('Overflow threshold')
      WRITE (NOUT,FMT=99009) 'overflow ' , eps
      eps = SLAMCH('Epsilon')
      WRITE (NOUT,FMT=99009) 'precision' , eps
      WRITE (NOUT,FMT=*)
      nrhs = nsval(1)
!
!
!     Read a test path and the number of matrix types to use.
!
 200  READ (NIN,FMT='(A72)',END=600) aline
      path = aline(1:3)
      nmats = MATMAX
      i = 3
      DO
         i = i + 1
         IF ( i>72 ) GOTO 500
         IF ( aline(i:i)/=' ' ) THEN
            nmats = 0
            EXIT
         ENDIF
      ENDDO
 300  c1 = aline(i:i)
      DO k = 1 , 10
         IF ( c1==intstr(k:k) ) THEN
            ic = k - 1
            GOTO 400
         ENDIF
      ENDDO
      GOTO 500
 400  nmats = nmats*10 + ic
      i = i + 1
      IF ( i<=72 ) GOTO 300
 500  c1 = path(1:1)
      c2 = path(2:3)
!
!     Check first character for correct precision.
!
      IF ( .NOT.LSAME(c1,'Complex precision') ) THEN
         WRITE (NOUT,FMT=99010) path
!
      ELSEIF ( nmats<=0 ) THEN
!
!        Check for a positive number of tests requested.
!
         WRITE (NOUT,FMT=99011) path
!
      ELSEIF ( LSAMEN(2,c2,'GE') ) THEN
!
!        GE:  general matrices
!
         ntypes = 11
         CALL ALAREQ(path,nmats,dotype,ntypes,NIN,NOUT)
!
         IF ( tstchk ) THEN
            CALL CCHKGE(dotype,nm,mval,nn,nval,nnb2,nbval2,nns,nsval,   &
     &                  thresh,tsterr,lda,a(1,1),a(1,2),a(1,3),b(1,1),  &
     &                  b(1,2),b(1,3),work,rwork,iwork,NOUT)
         ELSE
            WRITE (NOUT,FMT=99011) path
         ENDIF
!
         IF ( tstdrv ) THEN
            CALL CDRVGE(dotype,nn,nval,nrhs,thresh,tsterr,lda,a(1,1),   &
     &                  a(1,2),a(1,3),b(1,1),b(1,2),b(1,3),b(1,4),s,    &
     &                  work,rwork,iwork,NOUT)
         ELSE
            WRITE (NOUT,FMT=99012) path
         ENDIF
!
      ELSEIF ( LSAMEN(2,c2,'GB') ) THEN
!
!        GB:  general banded matrices
!
         la = (2*KDMAX+1)*NMAX
         lafac = (3*KDMAX+1)*NMAX
         ntypes = 8
         CALL ALAREQ(path,nmats,dotype,ntypes,NIN,NOUT)
!
         IF ( tstchk ) THEN
            CALL CCHKGB(dotype,nm,mval,nn,nval,nnb2,nbval2,nns,nsval,   &
     &                  thresh,tsterr,a(1,1),la,a(1,3),lafac,b(1,1),    &
     &                  b(1,2),b(1,3),work,rwork,iwork,NOUT)
         ELSE
            WRITE (NOUT,FMT=99011) path
         ENDIF
!
         IF ( tstdrv ) THEN
            CALL CDRVGB(dotype,nn,nval,nrhs,thresh,tsterr,a(1,1),la,    &
     &                  a(1,3),lafac,a(1,6),b(1,1),b(1,2),b(1,3),b(1,4),&
     &                  s,work,rwork,iwork,NOUT)
         ELSE
            WRITE (NOUT,FMT=99012) path
         ENDIF
!
      ELSEIF ( LSAMEN(2,c2,'GT') ) THEN
!
!        GT:  general tridiagonal matrices
!
         ntypes = 12
         CALL ALAREQ(path,nmats,dotype,ntypes,NIN,NOUT)
!
         IF ( tstchk ) THEN
            CALL CCHKGT(dotype,nn,nval,nns,nsval,thresh,tsterr,a(1,1),  &
     &                  a(1,2),b(1,1),b(1,2),b(1,3),work,rwork,iwork,   &
     &                  NOUT)
         ELSE
            WRITE (NOUT,FMT=99011) path
         ENDIF
!
         IF ( tstdrv ) THEN
            CALL CDRVGT(dotype,nn,nval,nrhs,thresh,tsterr,a(1,1),a(1,2),&
     &                  b(1,1),b(1,2),b(1,3),work,rwork,iwork,NOUT)
         ELSE
            WRITE (NOUT,FMT=99012) path
         ENDIF
!
      ELSEIF ( LSAMEN(2,c2,'PO') ) THEN
!
!        PO:  positive definite matrices
!
         ntypes = 9
         CALL ALAREQ(path,nmats,dotype,ntypes,NIN,NOUT)
!
         IF ( tstchk ) THEN
            CALL CCHKPO(dotype,nn,nval,nnb2,nbval2,nns,nsval,thresh,    &
     &                  tsterr,lda,a(1,1),a(1,2),a(1,3),b(1,1),b(1,2),  &
     &                  b(1,3),work,rwork,NOUT)
         ELSE
            WRITE (NOUT,FMT=99011) path
         ENDIF
!
         IF ( tstdrv ) THEN
            CALL CDRVPO(dotype,nn,nval,nrhs,thresh,tsterr,lda,a(1,1),   &
     &                  a(1,2),a(1,3),b(1,1),b(1,2),b(1,3),b(1,4),s,    &
     &                  work,rwork,NOUT)
         ELSE
            WRITE (NOUT,FMT=99012) path
         ENDIF
!
      ELSEIF ( LSAMEN(2,c2,'PS') ) THEN
!
!        PS:  positive semi-definite matrices
!
         ntypes = 9
!
         CALL ALAREQ(path,nmats,dotype,ntypes,NIN,NOUT)
!
         IF ( tstchk ) THEN
            CALL CCHKPS(dotype,nn,nval,nnb2,nbval2,nrank,rankval,thresh,&
     &                  tsterr,lda,a(1,1),a(1,2),a(1,3),piv,work,rwork, &
     &                  NOUT)
         ELSE
            WRITE (NOUT,FMT=99011) path
         ENDIF
!
      ELSEIF ( LSAMEN(2,c2,'PP') ) THEN
!
!        PP:  positive definite packed matrices
!
         ntypes = 9
         CALL ALAREQ(path,nmats,dotype,ntypes,NIN,NOUT)
!
         IF ( tstchk ) THEN
            CALL CCHKPP(dotype,nn,nval,nns,nsval,thresh,tsterr,lda,     &
     &                  a(1,1),a(1,2),a(1,3),b(1,1),b(1,2),b(1,3),work, &
     &                  rwork,NOUT)
         ELSE
            WRITE (NOUT,FMT=99011) path
         ENDIF
!
         IF ( tstdrv ) THEN
            CALL CDRVPP(dotype,nn,nval,nrhs,thresh,tsterr,lda,a(1,1),   &
     &                  a(1,2),a(1,3),b(1,1),b(1,2),b(1,3),b(1,4),s,    &
     &                  work,rwork,NOUT)
         ELSE
            WRITE (NOUT,FMT=99012) path
         ENDIF
!
      ELSEIF ( LSAMEN(2,c2,'PB') ) THEN
!
!        PB:  positive definite banded matrices
!
         ntypes = 8
         CALL ALAREQ(path,nmats,dotype,ntypes,NIN,NOUT)
!
         IF ( tstchk ) THEN
            CALL CCHKPB(dotype,nn,nval,nnb2,nbval2,nns,nsval,thresh,    &
     &                  tsterr,lda,a(1,1),a(1,2),a(1,3),b(1,1),b(1,2),  &
     &                  b(1,3),work,rwork,NOUT)
         ELSE
            WRITE (NOUT,FMT=99011) path
         ENDIF
!
         IF ( tstdrv ) THEN
            CALL CDRVPB(dotype,nn,nval,nrhs,thresh,tsterr,lda,a(1,1),   &
     &                  a(1,2),a(1,3),b(1,1),b(1,2),b(1,3),b(1,4),s,    &
     &                  work,rwork,NOUT)
         ELSE
            WRITE (NOUT,FMT=99012) path
         ENDIF
!
      ELSEIF ( LSAMEN(2,c2,'PT') ) THEN
!
!        PT:  positive definite tridiagonal matrices
!
         ntypes = 12
         CALL ALAREQ(path,nmats,dotype,ntypes,NIN,NOUT)
!
         IF ( tstchk ) THEN
            CALL CCHKPT(dotype,nn,nval,nns,nsval,thresh,tsterr,a(1,1),s,&
     &                  a(1,2),b(1,1),b(1,2),b(1,3),work,rwork,NOUT)
         ELSE
            WRITE (NOUT,FMT=99011) path
         ENDIF
!
         IF ( tstdrv ) THEN
            CALL CDRVPT(dotype,nn,nval,nrhs,thresh,tsterr,a(1,1),s,     &
     &                  a(1,2),b(1,1),b(1,2),b(1,3),work,rwork,NOUT)
         ELSE
            WRITE (NOUT,FMT=99012) path
         ENDIF
!
      ELSEIF ( LSAMEN(2,c2,'HE') ) THEN
!
!        HE:  Hermitian indefinite matrices,
!             with partial (Bunch-Kaufman) pivoting algorithm
!
         ntypes = 10
         CALL ALAREQ(path,nmats,dotype,ntypes,NIN,NOUT)
!
         IF ( tstchk ) THEN
            CALL CCHKHE(dotype,nn,nval,nnb2,nbval2,nns,nsval,thresh,    &
     &                  tsterr,lda,a(1,1),a(1,2),a(1,3),b(1,1),b(1,2),  &
     &                  b(1,3),work,rwork,iwork,NOUT)
         ELSE
            WRITE (NOUT,FMT=99011) path
         ENDIF
!
         IF ( tstdrv ) THEN
            CALL CDRVHE(dotype,nn,nval,nrhs,thresh,tsterr,lda,a(1,1),   &
     &                  a(1,2),a(1,3),b(1,1),b(1,2),b(1,3),work,rwork,  &
     &                  iwork,NOUT)
         ELSE
            WRITE (NOUT,FMT=99012) path
         ENDIF
!
      ELSEIF ( LSAMEN(2,c2,'HR') ) THEN
!
!        HR:  Hermitian indefinite matrices,
!             with bounded Bunch-Kaufman (rook) pivoting algorithm
!
         ntypes = 10
         CALL ALAREQ(path,nmats,dotype,ntypes,NIN,NOUT)
!
         IF ( tstchk ) THEN
            CALL CCHKHE_ROOK(dotype,nn,nval,nnb2,nbval2,nns,nsval,      &
     &                       thresh,tsterr,lda,a(1,1),a(1,2),a(1,3),    &
     &                       b(1,1),b(1,2),b(1,3),work,rwork,iwork,NOUT)
         ELSE
            WRITE (NOUT,FMT=99011) path
         ENDIF
!
         IF ( tstdrv ) THEN
            CALL CDRVHE_ROOK(dotype,nn,nval,nrhs,thresh,tsterr,lda,     &
     &                       a(1,1),a(1,2),a(1,3),b(1,1),b(1,2),b(1,3), &
     &                       work,rwork,iwork,NOUT)
         ELSE
            WRITE (NOUT,FMT=99012) path
         ENDIF
!
      ELSEIF ( LSAMEN(2,c2,'HK') ) THEN
!
!        HK:  Hermitian indefinite matrices,
!             with bounded Bunch-Kaufman (rook) pivoting algorithm,
!             different matrix storage format than HR path version.
!
         ntypes = 10
         CALL ALAREQ(path,nmats,dotype,ntypes,NIN,NOUT)
!
         IF ( tstchk ) THEN
            CALL CCHKHE_RK(dotype,nn,nval,nnb2,nbval2,nns,nsval,thresh, &
     &                     tsterr,lda,a(1,1),a(1,2),e,a(1,3),b(1,1),    &
     &                     b(1,2),b(1,3),work,rwork,iwork,NOUT)
         ELSE
            WRITE (NOUT,FMT=99011) path
         ENDIF
!
         IF ( tstdrv ) THEN
            CALL CDRVHE_RK(dotype,nn,nval,nrhs,thresh,tsterr,lda,a(1,1),&
     &                     a(1,2),e,a(1,3),b(1,1),b(1,2),b(1,3),work,   &
     &                     rwork,iwork,NOUT)
         ELSE
            WRITE (NOUT,FMT=99012) path
         ENDIF
!
      ELSEIF ( LSAMEN(2,c2,'HA') ) THEN
!
!        HA:  Hermitian matrices,
!             Aasen Algorithm
!
         ntypes = 10
         CALL ALAREQ(path,nmats,dotype,ntypes,NIN,NOUT)
!
         IF ( tstchk ) THEN
            CALL CCHKHE_AA(dotype,nn,nval,nnb2,nbval2,nns,nsval,thresh, &
     &                     tsterr,lda,a(1,1),a(1,2),a(1,3),b(1,1),b(1,2)&
     &                     ,b(1,3),work,rwork,iwork,NOUT)
         ELSE
            WRITE (NOUT,FMT=99011) path
         ENDIF
!
         IF ( tstdrv ) THEN
            CALL CDRVHE_AA(dotype,nn,nval,nrhs,thresh,tsterr,lda,a(1,1),&
     &                     a(1,2),a(1,3),b(1,1),b(1,2),b(1,3),work,     &
     &                     rwork,iwork,NOUT)
         ELSE
            WRITE (NOUT,FMT=99012) path
         ENDIF
!
      ELSEIF ( LSAMEN(2,c2,'H2') ) THEN
!
!        H2:  Hermitian matrices,
!             with partial (Aasen's) pivoting algorithm
!
         ntypes = 10
         CALL ALAREQ(path,nmats,dotype,ntypes,NIN,NOUT)
!
         IF ( tstchk ) THEN
            CALL CCHKHE_AA_2STAGE(dotype,nn,nval,nnb2,nbval2,nns,nsval, &
     &                            thresh,tsterr,lda,a(1,1),a(1,2),a(1,3)&
     &                            ,b(1,1),b(1,2),b(1,3),work,rwork,     &
     &                            iwork,NOUT)
         ELSE
            WRITE (NOUT,FMT=99011) path
         ENDIF
!
         IF ( tstdrv ) THEN
            CALL CDRVHE_AA_2STAGE(dotype,nn,nval,nrhs,thresh,tsterr,lda,&
     &                            a(1,1),a(1,2),a(1,3),b(1,1),b(1,2),   &
     &                            b(1,3),work,rwork,iwork,NOUT)
         ELSE
            WRITE (NOUT,FMT=99012) path
         ENDIF
!
      ELSEIF ( LSAMEN(2,c2,'HP') ) THEN
!
!        HP:  Hermitian indefinite packed matrices,
!             with partial (Bunch-Kaufman) pivoting algorithm
!
         ntypes = 10
         CALL ALAREQ(path,nmats,dotype,ntypes,NIN,NOUT)
!
         IF ( tstchk ) THEN
            CALL CCHKHP(dotype,nn,nval,nns,nsval,thresh,tsterr,lda,     &
     &                  a(1,1),a(1,2),a(1,3),b(1,1),b(1,2),b(1,3),work, &
     &                  rwork,iwork,NOUT)
         ELSE
            WRITE (NOUT,FMT=99011) path
         ENDIF
!
         IF ( tstdrv ) THEN
            CALL CDRVHP(dotype,nn,nval,nrhs,thresh,tsterr,lda,a(1,1),   &
     &                  a(1,2),a(1,3),b(1,1),b(1,2),b(1,3),work,rwork,  &
     &                  iwork,NOUT)
         ELSE
            WRITE (NOUT,FMT=99012) path
         ENDIF
!
      ELSEIF ( LSAMEN(2,c2,'SY') ) THEN
!
!        SY:  symmetric indefinite matrices,
!             with partial (Bunch-Kaufman) pivoting algorithm
!
         ntypes = 11
         CALL ALAREQ(path,nmats,dotype,ntypes,NIN,NOUT)
!
         IF ( tstchk ) THEN
            CALL CCHKSY(dotype,nn,nval,nnb2,nbval2,nns,nsval,thresh,    &
     &                  tsterr,lda,a(1,1),a(1,2),a(1,3),b(1,1),b(1,2),  &
     &                  b(1,3),work,rwork,iwork,NOUT)
         ELSE
            WRITE (NOUT,FMT=99011) path
         ENDIF
!
         IF ( tstdrv ) THEN
            CALL CDRVSY(dotype,nn,nval,nrhs,thresh,tsterr,lda,a(1,1),   &
     &                  a(1,2),a(1,3),b(1,1),b(1,2),b(1,3),work,rwork,  &
     &                  iwork,NOUT)
         ELSE
            WRITE (NOUT,FMT=99012) path
         ENDIF
!
      ELSEIF ( LSAMEN(2,c2,'SR') ) THEN
!
!        SR:  symmetric indefinite matrices,
!             with bounded Bunch-Kaufman (rook) pivoting algorithm
!
         ntypes = 11
         CALL ALAREQ(path,nmats,dotype,ntypes,NIN,NOUT)
!
         IF ( tstchk ) THEN
            CALL CCHKSY_ROOK(dotype,nn,nval,nnb2,nbval2,nns,nsval,      &
     &                       thresh,tsterr,lda,a(1,1),a(1,2),a(1,3),    &
     &                       b(1,1),b(1,2),b(1,3),work,rwork,iwork,NOUT)
         ELSE
            WRITE (NOUT,FMT=99011) path
         ENDIF
!
         IF ( tstdrv ) THEN
            CALL CDRVSY_ROOK(dotype,nn,nval,nrhs,thresh,tsterr,lda,     &
     &                       a(1,1),a(1,2),a(1,3),b(1,1),b(1,2),b(1,3), &
     &                       work,rwork,iwork,NOUT)
         ELSE
            WRITE (NOUT,FMT=99012) path
         ENDIF
!
      ELSEIF ( LSAMEN(2,c2,'SK') ) THEN
!
!        SK:  symmetric indefinite matrices,
!             with bounded Bunch-Kaufman (rook) pivoting algorithm,
!             different matrix storage format than SR path version.
!
         ntypes = 11
         CALL ALAREQ(path,nmats,dotype,ntypes,NIN,NOUT)
!
         IF ( tstchk ) THEN
            CALL CCHKSY_RK(dotype,nn,nval,nnb2,nbval2,nns,nsval,thresh, &
     &                     tsterr,lda,a(1,1),a(1,2),e,a(1,3),b(1,1),    &
     &                     b(1,2),b(1,3),work,rwork,iwork,NOUT)
         ELSE
            WRITE (NOUT,FMT=99011) path
         ENDIF
!
         IF ( tstdrv ) THEN
            CALL CDRVSY_RK(dotype,nn,nval,nrhs,thresh,tsterr,lda,a(1,1),&
     &                     a(1,2),e,a(1,3),b(1,1),b(1,2),b(1,3),work,   &
     &                     rwork,iwork,NOUT)
         ELSE
            WRITE (NOUT,FMT=99012) path
         ENDIF
!
      ELSEIF ( LSAMEN(2,c2,'SA') ) THEN
!
!        SA:  symmetric indefinite matrices with Aasen's algorithm,
!
         ntypes = 11
         CALL ALAREQ(path,nmats,dotype,ntypes,NIN,NOUT)
!
         IF ( tstchk ) THEN
            CALL CCHKSY_AA(dotype,nn,nval,nnb2,nbval2,nns,nsval,thresh, &
     &                     tsterr,lda,a(1,1),a(1,2),a(1,3),b(1,1),b(1,2)&
     &                     ,b(1,3),work,rwork,iwork,NOUT)
         ELSE
            WRITE (NOUT,FMT=99011) path
         ENDIF
!
         IF ( tstdrv ) THEN
            CALL CDRVSY_AA(dotype,nn,nval,nrhs,thresh,tsterr,lda,a(1,1),&
     &                     a(1,2),a(1,3),b(1,1),b(1,2),b(1,3),work,     &
     &                     rwork,iwork,NOUT)
         ELSE
            WRITE (NOUT,FMT=99012) path
         ENDIF
!
      ELSEIF ( LSAMEN(2,c2,'S2') ) THEN
!
!        S2:  symmetric indefinite matrices with Aasen's algorithm
!             2 stage
!
         ntypes = 11
         CALL ALAREQ(path,nmats,dotype,ntypes,NIN,NOUT)
!
         IF ( tstchk ) THEN
            CALL CCHKSY_AA_2STAGE(dotype,nn,nval,nnb2,nbval2,nns,nsval, &
     &                            thresh,tsterr,lda,a(1,1),a(1,2),a(1,3)&
     &                            ,b(1,1),b(1,2),b(1,3),work,rwork,     &
     &                            iwork,NOUT)
         ELSE
            WRITE (NOUT,FMT=99011) path
         ENDIF
!
         IF ( tstdrv ) THEN
            CALL CDRVSY_AA_2STAGE(dotype,nn,nval,nrhs,thresh,tsterr,lda,&
     &                            a(1,1),a(1,2),a(1,3),b(1,1),b(1,2),   &
     &                            b(1,3),work,rwork,iwork,NOUT)
         ELSE
            WRITE (NOUT,FMT=99012) path
         ENDIF
!
      ELSEIF ( LSAMEN(2,c2,'SP') ) THEN
!
!        SP:  symmetric indefinite packed matrices,
!             with partial (Bunch-Kaufman) pivoting algorithm
!
         ntypes = 11
         CALL ALAREQ(path,nmats,dotype,ntypes,NIN,NOUT)
!
         IF ( tstchk ) THEN
            CALL CCHKSP(dotype,nn,nval,nns,nsval,thresh,tsterr,lda,     &
     &                  a(1,1),a(1,2),a(1,3),b(1,1),b(1,2),b(1,3),work, &
     &                  rwork,iwork,NOUT)
         ELSE
            WRITE (NOUT,FMT=99011) path
         ENDIF
!
         IF ( tstdrv ) THEN
            CALL CDRVSP(dotype,nn,nval,nrhs,thresh,tsterr,lda,a(1,1),   &
     &                  a(1,2),a(1,3),b(1,1),b(1,2),b(1,3),work,rwork,  &
     &                  iwork,NOUT)
         ELSE
            WRITE (NOUT,FMT=99012) path
         ENDIF
!
      ELSEIF ( LSAMEN(2,c2,'TR') ) THEN
!
!        TR:  triangular matrices
!
         ntypes = 18
         CALL ALAREQ(path,nmats,dotype,ntypes,NIN,NOUT)
!
         IF ( tstchk ) THEN
            CALL CCHKTR(dotype,nn,nval,nnb2,nbval2,nns,nsval,thresh,    &
     &                  tsterr,lda,a(1,1),a(1,2),b(1,1),b(1,2),b(1,3),  &
     &                  work,rwork,NOUT)
         ELSE
            WRITE (NOUT,FMT=99011) path
         ENDIF
!
      ELSEIF ( LSAMEN(2,c2,'TP') ) THEN
!
!        TP:  triangular packed matrices
!
         ntypes = 18
         CALL ALAREQ(path,nmats,dotype,ntypes,NIN,NOUT)
!
         IF ( tstchk ) THEN
            CALL CCHKTP(dotype,nn,nval,nns,nsval,thresh,tsterr,lda,     &
     &                  a(1,1),a(1,2),b(1,1),b(1,2),b(1,3),work,rwork,  &
     &                  NOUT)
         ELSE
            WRITE (NOUT,FMT=99011) path
         ENDIF
!
      ELSEIF ( LSAMEN(2,c2,'TB') ) THEN
!
!        TB:  triangular banded matrices
!
         ntypes = 17
         CALL ALAREQ(path,nmats,dotype,ntypes,NIN,NOUT)
!
         IF ( tstchk ) THEN
            CALL CCHKTB(dotype,nn,nval,nns,nsval,thresh,tsterr,lda,     &
     &                  a(1,1),a(1,2),b(1,1),b(1,2),b(1,3),work,rwork,  &
     &                  NOUT)
         ELSE
            WRITE (NOUT,FMT=99011) path
         ENDIF
!
      ELSEIF ( LSAMEN(2,c2,'QR') ) THEN
!
!        QR:  QR factorization
!
         ntypes = 8
         CALL ALAREQ(path,nmats,dotype,ntypes,NIN,NOUT)
!
         IF ( tstchk ) THEN
            CALL CCHKQR(dotype,nm,mval,nn,nval,nnb,nbval,nxval,nrhs,    &
     &                  thresh,tsterr,NMAX,a(1,1),a(1,2),a(1,3),a(1,4), &
     &                  a(1,5),b(1,1),b(1,2),b(1,3),b(1,4),work,rwork,  &
     &                  iwork,NOUT)
         ELSE
            WRITE (NOUT,FMT=99011) path
         ENDIF
!
      ELSEIF ( LSAMEN(2,c2,'LQ') ) THEN
!
!        LQ:  LQ factorization
!
         ntypes = 8
         CALL ALAREQ(path,nmats,dotype,ntypes,NIN,NOUT)
!
         IF ( tstchk ) THEN
            CALL CCHKLQ(dotype,nm,mval,nn,nval,nnb,nbval,nxval,nrhs,    &
     &                  thresh,tsterr,NMAX,a(1,1),a(1,2),a(1,3),a(1,4), &
     &                  a(1,5),b(1,1),b(1,2),b(1,3),b(1,4),work,rwork,  &
     &                  NOUT)
         ELSE
            WRITE (NOUT,FMT=99011) path
         ENDIF
!
      ELSEIF ( LSAMEN(2,c2,'QL') ) THEN
!
!        QL:  QL factorization
!
         ntypes = 8
         CALL ALAREQ(path,nmats,dotype,ntypes,NIN,NOUT)
!
         IF ( tstchk ) THEN
            CALL CCHKQL(dotype,nm,mval,nn,nval,nnb,nbval,nxval,nrhs,    &
     &                  thresh,tsterr,NMAX,a(1,1),a(1,2),a(1,3),a(1,4), &
     &                  a(1,5),b(1,1),b(1,2),b(1,3),b(1,4),work,rwork,  &
     &                  NOUT)
         ELSE
            WRITE (NOUT,FMT=99011) path
         ENDIF
!
      ELSEIF ( LSAMEN(2,c2,'RQ') ) THEN
!
!        RQ:  RQ factorization
!
         ntypes = 8
         CALL ALAREQ(path,nmats,dotype,ntypes,NIN,NOUT)
!
         IF ( tstchk ) THEN
            CALL CCHKRQ(dotype,nm,mval,nn,nval,nnb,nbval,nxval,nrhs,    &
     &                  thresh,tsterr,NMAX,a(1,1),a(1,2),a(1,3),a(1,4), &
     &                  a(1,5),b(1,1),b(1,2),b(1,3),b(1,4),work,rwork,  &
     &                  iwork,NOUT)
         ELSE
            WRITE (NOUT,FMT=99011) path
         ENDIF
!
      ELSEIF ( LSAMEN(2,c2,'EQ') ) THEN
!
!        EQ:  Equilibration routines for general and positive definite
!             matrices (THREQ should be between 2 and 10)
!
         IF ( tstchk ) THEN
            CALL CCHKEQ(threq,NOUT)
         ELSE
            WRITE (NOUT,FMT=99011) path
         ENDIF
!
      ELSEIF ( LSAMEN(2,c2,'TZ') ) THEN
!
!        TZ:  Trapezoidal matrix
!
         ntypes = 3
         CALL ALAREQ(path,nmats,dotype,ntypes,NIN,NOUT)
!
         IF ( tstchk ) THEN
            CALL CCHKTZ(dotype,nm,mval,nn,nval,thresh,tsterr,a(1,1),    &
     &                  a(1,2),s(1),b(1,1),work,rwork,NOUT)
         ELSE
            WRITE (NOUT,FMT=99011) path
         ENDIF
!
      ELSEIF ( LSAMEN(2,c2,'QP') ) THEN
!
!        QP:  QR factorization with pivoting
!
         ntypes = 6
         CALL ALAREQ(path,nmats,dotype,ntypes,NIN,NOUT)
!
         IF ( tstchk ) THEN
            CALL CCHKQ3(dotype,nm,mval,nn,nval,nnb,nbval,nxval,thresh,  &
     &                  a(1,1),a(1,2),s(1),b(1,1),work,rwork,iwork,NOUT)
         ELSE
            WRITE (NOUT,FMT=99011) path
         ENDIF
!
      ELSEIF ( LSAMEN(2,c2,'LS') ) THEN
!
!        LS:  Least squares drivers
!
         ntypes = 6
         CALL ALAREQ(path,nmats,dotype,ntypes,NIN,NOUT)
!
         IF ( tstdrv ) THEN
            CALL CDRVLS(dotype,nm,mval,nn,nval,nns,nsval,nnb,nbval,     &
     &                  nxval,thresh,tsterr,a(1,1),a(1,2),a(1,3),a(1,4),&
     &                  a(1,5),s(1),s(NMAX+1),NOUT)
         ELSE
            WRITE (NOUT,FMT=99011) path
         ENDIF
!
      ELSEIF ( LSAMEN(2,c2,'QT') ) THEN
!
!        QT:  QRT routines for general matrices
!
         IF ( tstchk ) THEN
            CALL CCHKQRT(thresh,tsterr,nm,mval,nn,nval,nnb,nbval,NOUT)
         ELSE
            WRITE (NOUT,FMT=99011) path
         ENDIF
!
      ELSEIF ( LSAMEN(2,c2,'QX') ) THEN
!
!        QX:  QRT routines for triangular-pentagonal matrices
!
         IF ( tstchk ) THEN
            CALL CCHKQRTP(thresh,tsterr,nm,mval,nn,nval,nnb,nbval,NOUT)
         ELSE
            WRITE (NOUT,FMT=99011) path
         ENDIF
!
      ELSEIF ( LSAMEN(2,c2,'TQ') ) THEN
!
!        TQ:  LQT routines for general matrices
!
         IF ( tstchk ) THEN
            CALL CCHKLQT(thresh,tsterr,nm,mval,nn,nval,nnb,nbval,NOUT)
         ELSE
            WRITE (NOUT,FMT=99011) path
         ENDIF
!
      ELSEIF ( LSAMEN(2,c2,'XQ') ) THEN
!
!        XQ:  LQT routines for triangular-pentagonal matrices
!
         IF ( tstchk ) THEN
            CALL CCHKLQTP(thresh,tsterr,nm,mval,nn,nval,nnb,nbval,NOUT)
         ELSE
            WRITE (NOUT,FMT=99011) path
         ENDIF
!
      ELSEIF ( LSAMEN(2,c2,'TS') ) THEN
!
!        TS:  QR routines for tall-skinny matrices
!
         IF ( tstchk ) THEN
            CALL CCHKTSQR(thresh,tsterr,nm,mval,nn,nval,nnb,nbval,NOUT)
         ELSE
            WRITE (NOUT,FMT=99011) path
         ENDIF
!
      ELSEIF ( LSAMEN(2,c2,'HH') ) THEN
!
!        HH:  Householder reconstruction for tall-skinny matrices
!
         IF ( tstchk ) THEN
            CALL CCHKUNHR_COL(thresh,tsterr,nm,mval,nn,nval,nnb,nbval,  &
     &                        NOUT)
         ELSE
            WRITE (NOUT,FMT=99011) path
         ENDIF
!
      ELSE
!
         WRITE (NOUT,FMT=99010) path
      ENDIF
!
!     Go back to get another input line.
!
      GOTO 200
!
!     Branch to this line when the last record is read.
!
 600  CLOSE (NIN)
      s2 = SECOND()
      WRITE (NOUT,FMT=99002)
      WRITE (NOUT,FMT=99003) s2 - s1
!
99001 FORMAT (/' Execution not attempted due to input errors')
99002 FORMAT (/' End of tests')
99003 FORMAT (' Total time used = ',F12.2,' seconds',/)
99004 FORMAT (' Invalid input value: ',A4,'=',I6,'; must be >=',I6)
99005 FORMAT (' Invalid input value: ',A4,'=',I6,'; must be <=',I6)
99006 FORMAT (' Tests of the COMPLEX LAPACK routines ',                 &
     &        /' LAPACK VERSION ',I1,'.',I1,'.',I1,                     &
     &        //' The following parameter values will be used:')
99007 FORMAT (4X,A4,':  ',10I6,/11X,10I6)
99008 FORMAT (/' Routines pass computational tests if test ratio is ',  &
     &        'less than',F8.2,/)
99009 FORMAT (' Relative machine ',A,' is taken to be',E16.6)
99010 FORMAT (/1X,A3,':  Unrecognized path name')
99011 FORMAT (/1X,A3,' routines were not tested')
99012 FORMAT (/1X,A3,' driver routines were not tested')
!
!     End of CCHKAA
!
      END PROGRAM CCHKAA
