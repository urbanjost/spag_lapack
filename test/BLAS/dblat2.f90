!*==dblat2.f90  processed by SPAG 7.51RB at 17:52 on  4 Mar 2022
!> \brief \b DBLAT2
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       PROGRAM DBLAT2
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> Test program for the DOUBLE PRECISION Level 2 Blas.
!>
!> The program must be driven by a short data file. The first 18 records
!> of the file are read using list-directed input, the last 16 records
!> are read using the format ( A6, L2 ). An annotated example of a data
!> file can be obtained by deleting the first 3 characters from the
!> following 34 lines:
!> 'dblat2.out'      NAME OF SUMMARY OUTPUT FILE
!> 6                 UNIT NUMBER OF SUMMARY FILE
!> 'DBLAT2.SNAP'     NAME OF SNAPSHOT OUTPUT FILE
!> -1                UNIT NUMBER OF SNAPSHOT FILE (NOT USED IF .LT. 0)
!> F        LOGICAL FLAG, T TO REWIND SNAPSHOT FILE AFTER EACH RECORD.
!> F        LOGICAL FLAG, T TO STOP ON FAILURES.
!> T        LOGICAL FLAG, T TO TEST ERROR EXITS.
!> 16.0     THRESHOLD VALUE OF TEST RATIO
!> 6                 NUMBER OF VALUES OF N
!> 0 1 2 3 5 9       VALUES OF N
!> 4                 NUMBER OF VALUES OF K
!> 0 1 2 4           VALUES OF K
!> 4                 NUMBER OF VALUES OF INCX AND INCY
!> 1 2 -1 -2         VALUES OF INCX AND INCY
!> 3                 NUMBER OF VALUES OF ALPHA
!> 0.0 1.0 0.7       VALUES OF ALPHA
!> 3                 NUMBER OF VALUES OF BETA
!> 0.0 1.0 0.9       VALUES OF BETAC
!> DGEMV  T PUT F FOR NO TEST. SAME COLUMNS.
!> DGBMV  T PUT F FOR NO TEST. SAME COLUMNS.
!> DSYMV  T PUT F FOR NO TEST. SAME COLUMNS.
!> DSBMV  T PUT F FOR NO TEST. SAME COLUMNS.
!> DSPMV  T PUT F FOR NO TEST. SAME COLUMNS.
!> DTRMV  T PUT F FOR NO TEST. SAME COLUMNS.
!> DTBMV  T PUT F FOR NO TEST. SAME COLUMNS.
!> DTPMV  T PUT F FOR NO TEST. SAME COLUMNS.
!> DTRSV  T PUT F FOR NO TEST. SAME COLUMNS.
!> DTBSV  T PUT F FOR NO TEST. SAME COLUMNS.
!> DTPSV  T PUT F FOR NO TEST. SAME COLUMNS.
!> DGER   T PUT F FOR NO TEST. SAME COLUMNS.
!> DSYR   T PUT F FOR NO TEST. SAME COLUMNS.
!> DSPR   T PUT F FOR NO TEST. SAME COLUMNS.
!> DSYR2  T PUT F FOR NO TEST. SAME COLUMNS.
!> DSPR2  T PUT F FOR NO TEST. SAME COLUMNS.
!>
!> Further Details
!> ===============
!>
!>    See:
!>
!>       Dongarra J. J., Du Croz J. J., Hammarling S.  and Hanson R. J..
!>       An  extended  set of Fortran  Basic Linear Algebra Subprograms.
!>
!>       Technical  Memoranda  Nos. 41 (revision 3) and 81,  Mathematics
!>       and  Computer Science  Division,  Argonne  National Laboratory,
!>       9700 South Cass Avenue, Argonne, Illinois 60439, US.
!>
!>       Or
!>
!>       NAG  Technical Reports TR3/87 and TR4/87,  Numerical Algorithms
!>       Group  Ltd.,  NAG  Central  Office,  256  Banbury  Road, Oxford
!>       OX2 7DE, UK,  and  Numerical Algorithms Group Inc.,  1101  31st
!>       Street,  Suite 100,  Downers Grove,  Illinois 60515-1263,  USA.
!>
!>
!> -- Written on 10-August-1987.
!>    Richard Hanson, Sandia National Labs.
!>    Jeremy Du Croz, NAG Central Office.
!>
!>    10-9-00:  Change STATUS='NEW' to 'UNKNOWN' so that the testers
!>              can be run multiple times without deleting generated
!>              output files (susan)
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
!> \ingroup double_blas_testing
!
!  =====================================================================
      PROGRAM DBLAT2
      IMPLICIT NONE
!*--DBLAT2106
!
!  -- Reference BLAS test routine (version 3.7.0) --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     April 2012
!
!  =====================================================================
!
!     .. Parameters ..
      INTEGER NIN
      PARAMETER (NIN=5)
      INTEGER NSUBS
      PARAMETER (NSUBS=16)
      DOUBLE PRECISION ZERO , ONE
      PARAMETER (ZERO=0.0D0,ONE=1.0D0)
      INTEGER NMAX , INCMAX
      PARAMETER (NMAX=65,INCMAX=2)
      INTEGER NINMAX , NIDMAX , NKBMAX , NALMAX , NBEMAX
      PARAMETER (NINMAX=7,NIDMAX=9,NKBMAX=7,NALMAX=7,NBEMAX=7)
!     .. Local Scalars ..
      DOUBLE PRECISION eps , err , thresh
      INTEGER i , isnum , j , n , nalf , nbet , nidim , ninc , nkb ,    &
     &        nout , ntra
      LOGICAL fatal , ltestt , rewi , same , sfatal , trace , tsterr
      CHARACTER*1 trans
      CHARACTER*6 snamet
      CHARACTER*32 snaps , summry
!     .. Local Arrays ..
      DOUBLE PRECISION a(NMAX,NMAX) , aa(NMAX*NMAX) , alf(NALMAX) ,     &
     &                 as(NMAX*NMAX) , bet(NBEMAX) , g(NMAX) , x(NMAX) ,&
     &                 xs(NMAX*INCMAX) , xx(NMAX*INCMAX) , y(NMAX) ,    &
     &                 ys(NMAX*INCMAX) , yt(NMAX) , yy(NMAX*INCMAX) ,   &
     &                 z(2*NMAX)
      INTEGER idim(NIDMAX) , inc(NINMAX) , kb(NKBMAX)
      LOGICAL ltest(NSUBS)
      CHARACTER*6 snames(NSUBS)
!     .. External Functions ..
      DOUBLE PRECISION DDIFF
      LOGICAL LDE
      EXTERNAL DDIFF , LDE
!     .. External Subroutines ..
      EXTERNAL DCHK1 , DCHK2 , DCHK3 , DCHK4 , DCHK5 , DCHK6 , DCHKE ,  &
     &         DMVCH
!     .. Intrinsic Functions ..
      INTRINSIC ABS , MAX , MIN
!     .. Scalars in Common ..
      INTEGER INFot , NOUtc
      LOGICAL LERr , OK
      CHARACTER*6 SRNamt
!     .. Common blocks ..
      COMMON /INFOC / INFot , NOUtc , OK , LERr
      COMMON /SRNAMC/ SRNamt
!     .. Data statements ..
      DATA snames/'DGEMV ' , 'DGBMV ' , 'DSYMV ' , 'DSBMV ' , 'DSPMV ' ,&
     &     'DTRMV ' , 'DTBMV ' , 'DTPMV ' , 'DTRSV ' , 'DTBSV ' ,       &
     &     'DTPSV ' , 'DGER  ' , 'DSYR  ' , 'DSPR  ' , 'DSYR2 ' ,       &
     &     'DSPR2 '/
!     .. Executable Statements ..
!
!     Read name and unit number for summary output file and open file.
!
      READ (NIN,FMT=*) summry
      READ (NIN,FMT=*) nout
      OPEN (nout,FILE=summry,STATUS='UNKNOWN')
      NOUtc = nout
!
!     Read name and unit number for snapshot output file and open file.
!
      READ (NIN,FMT=*) snaps
      READ (NIN,FMT=*) ntra
      trace = ntra>=0
      IF ( trace ) OPEN (ntra,FILE=snaps,STATUS='UNKNOWN')
!     Read the flag that directs rewinding of the snapshot file.
      READ (NIN,FMT=*) rewi
      rewi = rewi .AND. trace
!     Read the flag that directs stopping on any failure.
      READ (NIN,FMT=*) sfatal
!     Read the flag that indicates whether error exits are to be tested.
      READ (NIN,FMT=*) tsterr
!     Read the threshold value of the test ratio
      READ (NIN,FMT=*) thresh
!
!     Read and check the parameter values for the tests.
!
!     Values of N
      READ (NIN,FMT=*) nidim
      IF ( nidim<1 .OR. nidim>NIDMAX ) THEN
         WRITE (nout,FMT=99003) 'N' , NIDMAX
         GOTO 500
      ENDIF
      READ (NIN,FMT=*) (idim(i),i=1,nidim)
      DO i = 1 , nidim
         IF ( idim(i)<0 .OR. idim(i)>NMAX ) THEN
            WRITE (nout,FMT=99004) NMAX
            GOTO 500
         ENDIF
      ENDDO
!     Values of K
      READ (NIN,FMT=*) nkb
      IF ( nkb<1 .OR. nkb>NKBMAX ) THEN
         WRITE (nout,FMT=99003) 'K' , NKBMAX
         GOTO 500
      ENDIF
      READ (NIN,FMT=*) (kb(i),i=1,nkb)
      DO i = 1 , nkb
         IF ( kb(i)<0 ) THEN
            WRITE (nout,FMT=99005)
            GOTO 500
         ENDIF
      ENDDO
!     Values of INCX and INCY
      READ (NIN,FMT=*) ninc
      IF ( ninc<1 .OR. ninc>NINMAX ) THEN
         WRITE (nout,FMT=99003) 'INCX AND INCY' , NINMAX
         GOTO 500
      ENDIF
      READ (NIN,FMT=*) (inc(i),i=1,ninc)
      DO i = 1 , ninc
         IF ( inc(i)==0 .OR. ABS(inc(i))>INCMAX ) THEN
            WRITE (nout,FMT=99006) INCMAX
            GOTO 500
         ENDIF
      ENDDO
!     Values of ALPHA
      READ (NIN,FMT=*) nalf
      IF ( nalf<1 .OR. nalf>NALMAX ) THEN
         WRITE (nout,FMT=99003) 'ALPHA' , NALMAX
         GOTO 500
      ENDIF
      READ (NIN,FMT=*) (alf(i),i=1,nalf)
!     Values of BETA
      READ (NIN,FMT=*) nbet
      IF ( nbet<1 .OR. nbet>NBEMAX ) THEN
         WRITE (nout,FMT=99003) 'BETA' , NBEMAX
         GOTO 500
      ENDIF
      READ (NIN,FMT=*) (bet(i),i=1,nbet)
!
!     Report values of parameters.
!
      WRITE (nout,FMT=99007)
      WRITE (nout,FMT=99008) (idim(i),i=1,nidim)
      WRITE (nout,FMT=99009) (kb(i),i=1,nkb)
      WRITE (nout,FMT=99010) (inc(i),i=1,ninc)
      WRITE (nout,FMT=99011) (alf(i),i=1,nalf)
      WRITE (nout,FMT=99012) (bet(i),i=1,nbet)
      IF ( .NOT.tsterr ) THEN
         WRITE (nout,FMT=*)
         WRITE (nout,FMT=99020)
      ENDIF
      WRITE (nout,FMT=*)
      WRITE (nout,FMT=99001) thresh
      WRITE (nout,FMT=*)
!
!     Read names of subroutines and flags which indicate
!     whether they are to be tested.
!
      DO i = 1 , NSUBS
         ltest(i) = .FALSE.
      ENDDO
 100  READ (NIN,FMT=99016,END=300) snamet , ltestt
      DO i = 1 , NSUBS
         IF ( snamet==snames(i) ) GOTO 200
      ENDDO
      WRITE (nout,FMT=99014) snamet
      STOP
 200  ltest(i) = ltestt
      GOTO 100
!
 300  CLOSE (NIN)
!
!     Compute EPS (the machine precision).
!
      eps = EPSILON(ZERO)
      WRITE (nout,FMT=99002) eps
!
!     Check the reliability of DMVCH using exact data.
!
      n = MIN(32,NMAX)
      DO j = 1 , n
         DO i = 1 , n
            a(i,j) = MAX(i-j+1,0)
         ENDDO
         x(j) = j
         y(j) = ZERO
      ENDDO
      DO j = 1 , n
         yy(j) = j*((j+1)*j)/2 - ((j+1)*j*(j-1))/3
      ENDDO
!     YY holds the exact result. On exit from DMVCH YT holds
!     the result computed by DMVCH.
      trans = 'N'
      CALL DMVCH(trans,n,n,ONE,a,NMAX,x,1,ZERO,y,1,yt,g,yy,eps,err,     &
     &           fatal,nout,.TRUE.)
      same = LDE(yy,yt,n)
      IF ( .NOT.same .OR. err/=ZERO ) THEN
         WRITE (nout,FMT=99015) trans , same , err
         STOP
      ENDIF
      trans = 'T'
      CALL DMVCH(trans,n,n,ONE,a,NMAX,x,-1,ZERO,y,-1,yt,g,yy,eps,err,   &
     &           fatal,nout,.TRUE.)
      same = LDE(yy,yt,n)
      IF ( .NOT.same .OR. err/=ZERO ) THEN
         WRITE (nout,FMT=99015) trans , same , err
         STOP
      ENDIF
!
!     Test each subroutine in turn.
!
      DO isnum = 1 , NSUBS
         WRITE (nout,FMT=*)
         IF ( .NOT.ltest(isnum) ) THEN
!           Subprogram is not to be tested.
            WRITE (nout,FMT=99017) snames(isnum)
         ELSE
            SRNamt = snames(isnum)
!           Test error exits.
            IF ( tsterr ) THEN
               CALL DCHKE(isnum,snames(isnum),nout)
               WRITE (nout,FMT=*)
            ENDIF
!           Test computations.
            INFot = 0
            OK = .TRUE.
            fatal = .FALSE.
            IF ( isnum==3 .OR. isnum==4 .OR. isnum==5 ) THEN
!           Test DSYMV, 03, DSBMV, 04, and DSPMV, 05.
               CALL DCHK2(snames(isnum),eps,thresh,nout,ntra,trace,rewi,&
     &                    fatal,nidim,idim,nkb,kb,nalf,alf,nbet,bet,    &
     &                    ninc,inc,NMAX,INCMAX,a,aa,as,x,xx,xs,y,yy,ys, &
     &                    yt,g)
            ELSEIF ( isnum==6 .OR. isnum==7 .OR. isnum==8 .OR.          &
     &               isnum==9 .OR. isnum==10 .OR. isnum==11 ) THEN
!           Test DTRMV, 06, DTBMV, 07, DTPMV, 08,
!           DTRSV, 09, DTBSV, 10, and DTPSV, 11.
               CALL DCHK3(snames(isnum),eps,thresh,nout,ntra,trace,rewi,&
     &                    fatal,nidim,idim,nkb,kb,ninc,inc,NMAX,INCMAX, &
     &                    a,aa,as,y,yy,ys,yt,g,z)
            ELSEIF ( isnum==12 ) THEN
!           Test DGER, 12.
               CALL DCHK4(snames(isnum),eps,thresh,nout,ntra,trace,rewi,&
     &                    fatal,nidim,idim,nalf,alf,ninc,inc,NMAX,      &
     &                    INCMAX,a,aa,as,x,xx,xs,y,yy,ys,yt,g,z)
            ELSEIF ( isnum==13 .OR. isnum==14 ) THEN
!           Test DSYR, 13, and DSPR, 14.
               CALL DCHK5(snames(isnum),eps,thresh,nout,ntra,trace,rewi,&
     &                    fatal,nidim,idim,nalf,alf,ninc,inc,NMAX,      &
     &                    INCMAX,a,aa,as,x,xx,xs,y,yy,ys,yt,g,z)
            ELSEIF ( isnum==15 .OR. isnum==16 ) THEN
!           Test DSYR2, 15, and DSPR2, 16.
               CALL DCHK6(snames(isnum),eps,thresh,nout,ntra,trace,rewi,&
     &                    fatal,nidim,idim,nalf,alf,ninc,inc,NMAX,      &
     &                    INCMAX,a,aa,as,x,xx,xs,y,yy,ys,yt,g,z)
            ELSE
!           Test DGEMV, 01, and DGBMV, 02.
               CALL DCHK1(snames(isnum),eps,thresh,nout,ntra,trace,rewi,&
     &                    fatal,nidim,idim,nkb,kb,nalf,alf,nbet,bet,    &
     &                    ninc,inc,NMAX,INCMAX,a,aa,as,x,xx,xs,y,yy,ys, &
     &                    yt,g)
            ENDIF
!
            IF ( fatal .AND. sfatal ) GOTO 400
         ENDIF
      ENDDO
      WRITE (nout,FMT=99018)
      GOTO 600
!
 400  WRITE (nout,FMT=99019)
      GOTO 600
!
 500  WRITE (nout,FMT=99013)
!
 600  IF ( trace ) CLOSE (ntra)
      CLOSE (nout)
      STOP
!
99001 FORMAT (' ROUTINES PASS COMPUTATIONAL TESTS IF TEST RATIO IS LES',&
     &        'S THAN',F8.2)
99002 FORMAT (' RELATIVE MACHINE PRECISION IS TAKEN TO BE',1P,D9.1)
99003 FORMAT (' NUMBER OF VALUES OF ',A,' IS LESS THAN 1 OR GREATER ',  &
     &        'THAN ',I2)
99004 FORMAT (' VALUE OF N IS LESS THAN 0 OR GREATER THAN ',I2)
99005 FORMAT (' VALUE OF K IS LESS THAN 0')
99006 FORMAT (' ABSOLUTE VALUE OF INCX OR INCY IS 0 OR GREATER THAN ',  &
     &        I2)
99007 FORMAT (' TESTS OF THE DOUBLE PRECISION LEVEL 2 BLAS',//' THE F', &
     &        'OLLOWING PARAMETER VALUES WILL BE USED:')
99008 FORMAT ('   FOR N              ',9I6)
99009 FORMAT ('   FOR K              ',7I6)
99010 FORMAT ('   FOR INCX AND INCY  ',7I6)
99011 FORMAT ('   FOR ALPHA          ',7F6.1)
99012 FORMAT ('   FOR BETA           ',7F6.1)
99013 FORMAT (' AMEND DATA FILE OR INCREASE ARRAY SIZES IN PROGRAM',    &
     &        /' ******* TESTS ABANDONED *******')
99014 FORMAT (' SUBPROGRAM NAME ',A6,' NOT RECOGNIZED',/' ******* T',   &
     &        'ESTS ABANDONED *******')
99015 FORMAT (' ERROR IN DMVCH -  IN-LINE DOT PRODUCTS ARE BEING EVALU',&
     &        'ATED WRONGLY.',/' DMVCH WAS CALLED WITH TRANS = ',A1,    &
     &        ' AND RETURNED SAME = ',L1,' AND ERR = ',F12.3,'.',/      &
     &   ' THIS MAY BE DUE TO FAULTS IN THE ARITHMETIC OR THE COMPILER.'&
     &   ,/' ******* TESTS ABANDONED *******')
99016 FORMAT (A6,L2)
99017 FORMAT (1X,A6,' WAS NOT TESTED')
99018 FORMAT (/' END OF TESTS')
99019 FORMAT (/' ******* FATAL ERROR - TESTS ABANDONED *******')
99020 FORMAT (' ERROR-EXITS WILL NOT BE TESTED')
!
!     End of DBLAT2.
!
      END PROGRAM DBLAT2
!*==dchk1.f90  processed by SPAG 7.51RB at 17:52 on  4 Mar 2022
      SUBROUTINE DCHK1(Sname,Eps,Thresh,Nout,Ntra,Trace,Rewi,Fatal,     &
     &                 Nidim,Idim,Nkb,Kb,Nalf,Alf,Nbet,Bet,Ninc,Inc,    &
     &                 Nmax,Incmax,A,Aa,As,X,Xx,Xs,Y,Yy,Ys,Yt,G)
      IMPLICIT NONE
!*--DCHK1423
!
!  Tests DGEMV and DGBMV.
!
!  Auxiliary routine for test program for Level 2 Blas.
!
!  -- Written on 10-August-1987.
!     Richard Hanson, Sandia National Labs.
!     Jeremy Du Croz, NAG Central Office.
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO , HALF
      PARAMETER (ZERO=0.0D0,HALF=0.5D0)
!     .. Scalar Arguments ..
      DOUBLE PRECISION Eps , Thresh
      INTEGER Incmax , Nalf , Nbet , Nidim , Ninc , Nkb , Nmax , Nout , &
     &        Ntra
      LOGICAL Fatal , Rewi , Trace
      CHARACTER*6 Sname
!     .. Array Arguments ..
      DOUBLE PRECISION A(Nmax,Nmax) , Aa(Nmax*Nmax) , Alf(Nalf) ,       &
     &                 As(Nmax*Nmax) , Bet(Nbet) , G(Nmax) , X(Nmax) ,  &
     &                 Xs(Nmax*Incmax) , Xx(Nmax*Incmax) , Y(Nmax) ,    &
     &                 Ys(Nmax*Incmax) , Yt(Nmax) , Yy(Nmax*Incmax)
      INTEGER Idim(Nidim) , Inc(Ninc) , Kb(Nkb)
!     .. Local Scalars ..
      DOUBLE PRECISION alpha , als , beta , bls , err , errmax , transl
      INTEGER i , ia , ib , ic , iku , im , in , incx , incxs , incy ,  &
     &        incys , ix , iy , kl , kls , ku , kus , laa , lda , ldas ,&
     &        lx , ly , m , ml , ms , n , nargs , nc , nd , nk , nl , ns
      LOGICAL banded , full , null , reset , same , tran
      CHARACTER*1 trans , transs
      CHARACTER*3 ich
!     .. Local Arrays ..
      LOGICAL isame(13)
!     .. External Functions ..
      LOGICAL LDE , LDERES
      EXTERNAL LDE , LDERES
!     .. External Subroutines ..
      EXTERNAL DGBMV , DGEMV , DMAKE , DMVCH
!     .. Intrinsic Functions ..
      INTRINSIC ABS , MAX , MIN
!     .. Scalars in Common ..
      INTEGER INFot , NOUtc
      LOGICAL LERr , OK
!     .. Common blocks ..
      COMMON /INFOC / INFot , NOUtc , OK , LERr
!     .. Data statements ..
      DATA ich/'NTC'/
!     .. Executable Statements ..
      full = Sname(3:3)=='E'
      banded = Sname(3:3)=='B'
!     Define the number of arguments.
      IF ( full ) THEN
         nargs = 11
      ELSEIF ( banded ) THEN
         nargs = 13
      ENDIF
!
      nc = 0
      reset = .TRUE.
      errmax = ZERO
!
      DO in = 1 , Nidim
         n = Idim(in)
         nd = n/2 + 1
!
         DO im = 1 , 2
            IF ( im==1 ) m = MAX(n-nd,0)
            IF ( im==2 ) m = MIN(n+nd,Nmax)
!
            IF ( banded ) THEN
               nk = Nkb
            ELSE
               nk = 1
            ENDIF
            DO iku = 1 , nk
               IF ( banded ) THEN
                  ku = Kb(iku)
                  kl = MAX(ku-1,0)
               ELSE
                  ku = n - 1
                  kl = m - 1
               ENDIF
!              Set LDA to 1 more than minimum value if room.
               IF ( banded ) THEN
                  lda = kl + ku + 1
               ELSE
                  lda = m
               ENDIF
               IF ( lda<Nmax ) lda = lda + 1
!              Skip tests if not enough room.
               IF ( lda<=Nmax ) THEN
                  laa = lda*n
                  null = n<=0 .OR. m<=0
!
!              Generate the matrix A.
!
                  transl = ZERO
                  CALL DMAKE(Sname(2:3),' ',' ',m,n,A,Nmax,Aa,lda,kl,ku,&
     &                       reset,transl)
!
                  DO ic = 1 , 3
                     trans = ich(ic:ic)
                     tran = trans=='T' .OR. trans=='C'
!
                     IF ( tran ) THEN
                        ml = n
                        nl = m
                     ELSE
                        ml = m
                        nl = n
                     ENDIF
!
                     DO ix = 1 , Ninc
                        incx = Inc(ix)
                        lx = ABS(incx)*nl
!
!                    Generate the vector X.
!
                        transl = HALF
                        CALL DMAKE('GE',' ',' ',1,nl,X,1,Xx,ABS(incx),0,&
     &                             nl-1,reset,transl)
                        IF ( nl>1 ) THEN
                           X(nl/2) = ZERO
                           Xx(1+ABS(incx)*(nl/2-1)) = ZERO
                        ENDIF
!
                        DO iy = 1 , Ninc
                           incy = Inc(iy)
                           ly = ABS(incy)*ml
!
                           DO ia = 1 , Nalf
                              alpha = Alf(ia)
!
                              DO ib = 1 , Nbet
                                 beta = Bet(ib)
!
!                             Generate the vector Y.
!
                                 transl = ZERO
                                 CALL DMAKE('GE',' ',' ',1,ml,Y,1,Yy,   &
     &                              ABS(incy),0,ml-1,reset,transl)
!
                                 nc = nc + 1
!
!                             Save every datum before calling the
!                             subroutine.
!
                                 transs = trans
                                 ms = m
                                 ns = n
                                 kls = kl
                                 kus = ku
                                 als = alpha
                                 DO i = 1 , laa
                                    As(i) = Aa(i)
                                 ENDDO
                                 ldas = lda
                                 DO i = 1 , lx
                                    Xs(i) = Xx(i)
                                 ENDDO
                                 incxs = incx
                                 bls = beta
                                 DO i = 1 , ly
                                    Ys(i) = Yy(i)
                                 ENDDO
                                 incys = incy
!
!                             Call the subroutine.
!
                                 IF ( full ) THEN
                                    IF ( Trace ) WRITE (Ntra,FMT=99006) &
     &                                 nc , Sname , trans , m , n ,     &
     &                                 alpha , lda , incx , beta , incy
                                    IF ( Rewi ) REWIND Ntra
                                    CALL DGEMV(trans,m,n,alpha,Aa,lda,  &
     &                                 Xx,incx,beta,Yy,incy)
                                 ELSEIF ( banded ) THEN
                                    IF ( Trace ) WRITE (Ntra,FMT=99005) &
     &                                 nc , Sname , trans , m , n , kl ,&
     &                                 ku , alpha , lda , incx , beta , &
     &                                 incy
                                    IF ( Rewi ) REWIND Ntra
                                    CALL DGBMV(trans,m,n,kl,ku,alpha,Aa,&
     &                                 lda,Xx,incx,beta,Yy,incy)
                                 ENDIF
!
!                             Check if error-exit was taken incorrectly.
!
                                 IF ( .NOT.OK ) THEN
                                    WRITE (Nout,FMT=99007)
                                    Fatal = .TRUE.
                                    GOTO 100
                                 ENDIF
!
!                             See what data changed inside subroutines.
!
                                 isame(1) = trans==transs
                                 isame(2) = ms==m
                                 isame(3) = ns==n
                                 IF ( full ) THEN
                                    isame(4) = als==alpha
                                    isame(5) = LDE(As,Aa,laa)
                                    isame(6) = ldas==lda
                                    isame(7) = LDE(Xs,Xx,lx)
                                    isame(8) = incxs==incx
                                    isame(9) = bls==beta
                                    IF ( null ) THEN
                                       isame(10) = LDE(Ys,Yy,ly)
                                    ELSE
                                       isame(10)                        &
     &                                    = LDERES('GE',' ',1,ml,Ys,Yy, &
     &                                    ABS(incy))
                                    ENDIF
                                    isame(11) = incys==incy
                                 ELSEIF ( banded ) THEN
                                    isame(4) = kls==kl
                                    isame(5) = kus==ku
                                    isame(6) = als==alpha
                                    isame(7) = LDE(As,Aa,laa)
                                    isame(8) = ldas==lda
                                    isame(9) = LDE(Xs,Xx,lx)
                                    isame(10) = incxs==incx
                                    isame(11) = bls==beta
                                    IF ( null ) THEN
                                       isame(12) = LDE(Ys,Yy,ly)
                                    ELSE
                                       isame(12)                        &
     &                                    = LDERES('GE',' ',1,ml,Ys,Yy, &
     &                                    ABS(incy))
                                    ENDIF
                                    isame(13) = incys==incy
                                 ENDIF
!
!                             If data was incorrectly changed, report
!                             and return.
!
                                 same = .TRUE.
                                 DO i = 1 , nargs
                                    same = same .AND. isame(i)
                                    IF ( .NOT.isame(i) )                &
     &                                 WRITE (Nout,FMT=99002) i
                                 ENDDO
                                 IF ( .NOT.same ) THEN
                                    Fatal = .TRUE.
                                    GOTO 100
                                 ENDIF
!
!                                Avoid repeating tests with M.le.0 or
!                                N.le.0.
                                 IF ( null ) GOTO 50
!
!                                Check the result.
!
                                 CALL DMVCH(trans,m,n,alpha,A,Nmax,X,   &
     &                              incx,beta,Y,incy,Yt,G,Yy,Eps,err,   &
     &                              Fatal,Nout,.TRUE.)
                                 errmax = MAX(errmax,err)
!                                If got really bad answer, report and
!                                return.
                                 IF ( Fatal ) GOTO 100
!
                              ENDDO
!
                           ENDDO
!
                        ENDDO
!
                     ENDDO
!
                  ENDDO
               ENDIF
!
            ENDDO
!
 50      ENDDO
!
      ENDDO
!
!     Report result.
!
      IF ( errmax<Thresh ) THEN
         WRITE (Nout,FMT=99001) Sname , nc
      ELSE
         WRITE (Nout,FMT=99003) Sname , nc , errmax
      ENDIF
      GOTO 200
!
 100  WRITE (Nout,FMT=99004) Sname
      IF ( full ) THEN
         WRITE (Nout,FMT=99006) nc , Sname , trans , m , n , alpha ,    &
     &                          lda , incx , beta , incy
      ELSEIF ( banded ) THEN
         WRITE (Nout,FMT=99005) nc , Sname , trans , m , n , kl , ku ,  &
     &                          alpha , lda , incx , beta , incy
      ENDIF
!
 200  RETURN
!
99001 FORMAT (' ',A6,' PASSED THE COMPUTATIONAL TESTS (',I6,' CALL',    &
     &        'S)')
99002 FORMAT (' ******* FATAL ERROR - PARAMETER NUMBER ',I2,' WAS CH',  &
     &        'ANGED INCORRECTLY *******')
99003 FORMAT (' ',A6,' COMPLETED THE COMPUTATIONAL TESTS (',I6,' C',    &
     &        'ALLS)',/' ******* BUT WITH MAXIMUM TEST RATIO',F8.2,     &
     &        ' - SUSPECT *******')
99004 FORMAT (' ******* ',A6,' FAILED ON CALL NUMBER:')
99005 FORMAT (1X,I6,': ',A6,'(''',A1,''',',4(I3,','),F4.1,', A,',I3,    &
     &        ', X,',I2,',',F4.1,', Y,',I2,') .')
99006 FORMAT (1X,I6,': ',A6,'(''',A1,''',',2(I3,','),F4.1,', A,',I3,    &
     &        ', X,',I2,',',F4.1,', Y,',I2,')         .')
99007 FORMAT (' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *',&
     &        '******')
!
!     End of DCHK1.
!
      END SUBROUTINE DCHK1
!*==dchk2.f90  processed by SPAG 7.51RB at 17:52 on  4 Mar 2022
      SUBROUTINE DCHK2(Sname,Eps,Thresh,Nout,Ntra,Trace,Rewi,Fatal,     &
     &                 Nidim,Idim,Nkb,Kb,Nalf,Alf,Nbet,Bet,Ninc,Inc,    &
     &                 Nmax,Incmax,A,Aa,As,X,Xx,Xs,Y,Yy,Ys,Yt,G)
      IMPLICIT NONE
!*--DCHK2746
!
!  Tests DSYMV, DSBMV and DSPMV.
!
!  Auxiliary routine for test program for Level 2 Blas.
!
!  -- Written on 10-August-1987.
!     Richard Hanson, Sandia National Labs.
!     Jeremy Du Croz, NAG Central Office.
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO , HALF
      PARAMETER (ZERO=0.0D0,HALF=0.5D0)
!     .. Scalar Arguments ..
      DOUBLE PRECISION Eps , Thresh
      INTEGER Incmax , Nalf , Nbet , Nidim , Ninc , Nkb , Nmax , Nout , &
     &        Ntra
      LOGICAL Fatal , Rewi , Trace
      CHARACTER*6 Sname
!     .. Array Arguments ..
      DOUBLE PRECISION A(Nmax,Nmax) , Aa(Nmax*Nmax) , Alf(Nalf) ,       &
     &                 As(Nmax*Nmax) , Bet(Nbet) , G(Nmax) , X(Nmax) ,  &
     &                 Xs(Nmax*Incmax) , Xx(Nmax*Incmax) , Y(Nmax) ,    &
     &                 Ys(Nmax*Incmax) , Yt(Nmax) , Yy(Nmax*Incmax)
      INTEGER Idim(Nidim) , Inc(Ninc) , Kb(Nkb)
!     .. Local Scalars ..
      DOUBLE PRECISION alpha , als , beta , bls , err , errmax , transl
      INTEGER i , ia , ib , ic , ik , in , incx , incxs , incy , incys ,&
     &        ix , iy , k , ks , laa , lda , ldas , lx , ly , n ,       &
     &        nargs , nc , nk , ns
      LOGICAL banded , full , null , packed , reset , same
      CHARACTER*1 uplo , uplos
      CHARACTER*2 ich
!     .. Local Arrays ..
      LOGICAL isame(13)
!     .. External Functions ..
      LOGICAL LDE , LDERES
      EXTERNAL LDE , LDERES
!     .. External Subroutines ..
      EXTERNAL DMAKE , DMVCH , DSBMV , DSPMV , DSYMV
!     .. Intrinsic Functions ..
      INTRINSIC ABS , MAX
!     .. Scalars in Common ..
      INTEGER INFot , NOUtc
      LOGICAL LERr , OK
!     .. Common blocks ..
      COMMON /INFOC / INFot , NOUtc , OK , LERr
!     .. Data statements ..
      DATA ich/'UL'/
!     .. Executable Statements ..
      full = Sname(3:3)=='Y'
      banded = Sname(3:3)=='B'
      packed = Sname(3:3)=='P'
!     Define the number of arguments.
      IF ( full ) THEN
         nargs = 10
      ELSEIF ( banded ) THEN
         nargs = 11
      ELSEIF ( packed ) THEN
         nargs = 9
      ENDIF
!
      nc = 0
      reset = .TRUE.
      errmax = ZERO
!
      DO in = 1 , Nidim
         n = Idim(in)
!
         IF ( banded ) THEN
            nk = Nkb
         ELSE
            nk = 1
         ENDIF
         DO ik = 1 , nk
            IF ( banded ) THEN
               k = Kb(ik)
            ELSE
               k = n - 1
            ENDIF
!           Set LDA to 1 more than minimum value if room.
            IF ( banded ) THEN
               lda = k + 1
            ELSE
               lda = n
            ENDIF
            IF ( lda<Nmax ) lda = lda + 1
!           Skip tests if not enough room.
            IF ( lda<=Nmax ) THEN
               IF ( packed ) THEN
                  laa = (n*(n+1))/2
               ELSE
                  laa = lda*n
               ENDIF
               null = n<=0
!
               DO ic = 1 , 2
                  uplo = ich(ic:ic)
!
!              Generate the matrix A.
!
                  transl = ZERO
                  CALL DMAKE(Sname(2:3),uplo,' ',n,n,A,Nmax,Aa,lda,k,k, &
     &                       reset,transl)
!
                  DO ix = 1 , Ninc
                     incx = Inc(ix)
                     lx = ABS(incx)*n
!
!                 Generate the vector X.
!
                     transl = HALF
                     CALL DMAKE('GE',' ',' ',1,n,X,1,Xx,ABS(incx),0,n-1,&
     &                          reset,transl)
                     IF ( n>1 ) THEN
                        X(n/2) = ZERO
                        Xx(1+ABS(incx)*(n/2-1)) = ZERO
                     ENDIF
!
                     DO iy = 1 , Ninc
                        incy = Inc(iy)
                        ly = ABS(incy)*n
!
                        DO ia = 1 , Nalf
                           alpha = Alf(ia)
!
                           DO ib = 1 , Nbet
                              beta = Bet(ib)
!
!                          Generate the vector Y.
!
                              transl = ZERO
                              CALL DMAKE('GE',' ',' ',1,n,Y,1,Yy,       &
     &                           ABS(incy),0,n-1,reset,transl)
!
                              nc = nc + 1
!
!                          Save every datum before calling the
!                          subroutine.
!
                              uplos = uplo
                              ns = n
                              ks = k
                              als = alpha
                              DO i = 1 , laa
                                 As(i) = Aa(i)
                              ENDDO
                              ldas = lda
                              DO i = 1 , lx
                                 Xs(i) = Xx(i)
                              ENDDO
                              incxs = incx
                              bls = beta
                              DO i = 1 , ly
                                 Ys(i) = Yy(i)
                              ENDDO
                              incys = incy
!
!                          Call the subroutine.
!
                              IF ( full ) THEN
                                 IF ( Trace ) WRITE (Ntra,FMT=99007)    &
     &                                nc , Sname , uplo , n , alpha ,   &
     &                                lda , incx , beta , incy
                                 IF ( Rewi ) REWIND Ntra
                                 CALL DSYMV(uplo,n,alpha,Aa,lda,Xx,incx,&
     &                              beta,Yy,incy)
                              ELSEIF ( banded ) THEN
                                 IF ( Trace ) WRITE (Ntra,FMT=99006)    &
     &                                nc , Sname , uplo , n , k ,       &
     &                                alpha , lda , incx , beta , incy
                                 IF ( Rewi ) REWIND Ntra
                                 CALL DSBMV(uplo,n,k,alpha,Aa,lda,Xx,   &
     &                              incx,beta,Yy,incy)
                              ELSEIF ( packed ) THEN
                                 IF ( Trace ) WRITE (Ntra,FMT=99005)    &
     &                                nc , Sname , uplo , n , alpha ,   &
     &                                incx , beta , incy
                                 IF ( Rewi ) REWIND Ntra
                                 CALL DSPMV(uplo,n,alpha,Aa,Xx,incx,    &
     &                              beta,Yy,incy)
                              ENDIF
!
!                          Check if error-exit was taken incorrectly.
!
                              IF ( .NOT.OK ) THEN
                                 WRITE (Nout,FMT=99008)
                                 Fatal = .TRUE.
                                 GOTO 200
                              ENDIF
!
!                          See what data changed inside subroutines.
!
                              isame(1) = uplo==uplos
                              isame(2) = ns==n
                              IF ( full ) THEN
                                 isame(3) = als==alpha
                                 isame(4) = LDE(As,Aa,laa)
                                 isame(5) = ldas==lda
                                 isame(6) = LDE(Xs,Xx,lx)
                                 isame(7) = incxs==incx
                                 isame(8) = bls==beta
                                 IF ( null ) THEN
                                    isame(9) = LDE(Ys,Yy,ly)
                                 ELSE
                                    isame(9)                            &
     &                                 = LDERES('GE',' ',1,n,Ys,Yy,     &
     &                                 ABS(incy))
                                 ENDIF
                                 isame(10) = incys==incy
                              ELSEIF ( banded ) THEN
                                 isame(3) = ks==k
                                 isame(4) = als==alpha
                                 isame(5) = LDE(As,Aa,laa)
                                 isame(6) = ldas==lda
                                 isame(7) = LDE(Xs,Xx,lx)
                                 isame(8) = incxs==incx
                                 isame(9) = bls==beta
                                 IF ( null ) THEN
                                    isame(10) = LDE(Ys,Yy,ly)
                                 ELSE
                                    isame(10)                           &
     &                                 = LDERES('GE',' ',1,n,Ys,Yy,     &
     &                                 ABS(incy))
                                 ENDIF
                                 isame(11) = incys==incy
                              ELSEIF ( packed ) THEN
                                 isame(3) = als==alpha
                                 isame(4) = LDE(As,Aa,laa)
                                 isame(5) = LDE(Xs,Xx,lx)
                                 isame(6) = incxs==incx
                                 isame(7) = bls==beta
                                 IF ( null ) THEN
                                    isame(8) = LDE(Ys,Yy,ly)
                                 ELSE
                                    isame(8)                            &
     &                                 = LDERES('GE',' ',1,n,Ys,Yy,     &
     &                                 ABS(incy))
                                 ENDIF
                                 isame(9) = incys==incy
                              ENDIF
!
!                          If data was incorrectly changed, report and
!                          return.
!
                              same = .TRUE.
                              DO i = 1 , nargs
                                 same = same .AND. isame(i)
                                 IF ( .NOT.isame(i) )                   &
     &                                WRITE (Nout,FMT=99002) i
                              ENDDO
                              IF ( .NOT.same ) THEN
                                 Fatal = .TRUE.
                                 GOTO 200
                              ENDIF
!
!                             Avoid repeating tests with N.le.0
                              IF ( null ) GOTO 100
!
!                             Check the result.
!
                              CALL DMVCH('N',n,n,alpha,A,Nmax,X,incx,   &
     &                           beta,Y,incy,Yt,G,Yy,Eps,err,Fatal,Nout,&
     &                           .TRUE.)
                              errmax = MAX(errmax,err)
!                             If got really bad answer, report and
!                             return.
                              IF ( Fatal ) GOTO 200
!
                           ENDDO
!
                        ENDDO
!
                     ENDDO
!
                  ENDDO
!
               ENDDO
            ENDIF
!
         ENDDO
!
 100  ENDDO
!
!     Report result.
!
      IF ( errmax<Thresh ) THEN
         WRITE (Nout,FMT=99001) Sname , nc
      ELSE
         WRITE (Nout,FMT=99003) Sname , nc , errmax
      ENDIF
      GOTO 300
!
 200  WRITE (Nout,FMT=99004) Sname
      IF ( full ) THEN
         WRITE (Nout,FMT=99007) nc , Sname , uplo , n , alpha , lda ,   &
     &                          incx , beta , incy
      ELSEIF ( banded ) THEN
         WRITE (Nout,FMT=99006) nc , Sname , uplo , n , k , alpha ,     &
     &                          lda , incx , beta , incy
      ELSEIF ( packed ) THEN
         WRITE (Nout,FMT=99005) nc , Sname , uplo , n , alpha , incx ,  &
     &                          beta , incy
      ENDIF
!
 300  RETURN
!
99001 FORMAT (' ',A6,' PASSED THE COMPUTATIONAL TESTS (',I6,' CALL',    &
     &        'S)')
99002 FORMAT (' ******* FATAL ERROR - PARAMETER NUMBER ',I2,' WAS CH',  &
     &        'ANGED INCORRECTLY *******')
99003 FORMAT (' ',A6,' COMPLETED THE COMPUTATIONAL TESTS (',I6,' C',    &
     &        'ALLS)',/' ******* BUT WITH MAXIMUM TEST RATIO',F8.2,     &
     &        ' - SUSPECT *******')
99004 FORMAT (' ******* ',A6,' FAILED ON CALL NUMBER:')
99005 FORMAT (1X,I6,': ',A6,'(''',A1,''',',I3,',',F4.1,', AP',', X,',I2,&
     &        ',',F4.1,', Y,',I2,')                .')
99006 FORMAT (1X,I6,': ',A6,'(''',A1,''',',2(I3,','),F4.1,', A,',I3,    &
     &        ', X,',I2,',',F4.1,', Y,',I2,')         .')
99007 FORMAT (1X,I6,': ',A6,'(''',A1,''',',I3,',',F4.1,', A,',I3,', X,',&
     &        I2,',',F4.1,', Y,',I2,')             .')
99008 FORMAT (' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *',&
     &        '******')
!
!     End of DCHK2.
!
      END SUBROUTINE DCHK2
!*==dchk3.f90  processed by SPAG 7.51RB at 17:52 on  4 Mar 2022
      SUBROUTINE DCHK3(Sname,Eps,Thresh,Nout,Ntra,Trace,Rewi,Fatal,     &
     &                 Nidim,Idim,Nkb,Kb,Ninc,Inc,Nmax,Incmax,A,Aa,As,X,&
     &                 Xx,Xs,Xt,G,Z)
      IMPLICIT NONE
!*--DCHK31078
!
!  Tests DTRMV, DTBMV, DTPMV, DTRSV, DTBSV and DTPSV.
!
!  Auxiliary routine for test program for Level 2 Blas.
!
!  -- Written on 10-August-1987.
!     Richard Hanson, Sandia National Labs.
!     Jeremy Du Croz, NAG Central Office.
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO , HALF , ONE
      PARAMETER (ZERO=0.0D0,HALF=0.5D0,ONE=1.0D0)
!     .. Scalar Arguments ..
      DOUBLE PRECISION Eps , Thresh
      INTEGER Incmax , Nidim , Ninc , Nkb , Nmax , Nout , Ntra
      LOGICAL Fatal , Rewi , Trace
      CHARACTER*6 Sname
!     .. Array Arguments ..
      DOUBLE PRECISION A(Nmax,Nmax) , Aa(Nmax*Nmax) , As(Nmax*Nmax) ,   &
     &                 G(Nmax) , X(Nmax) , Xs(Nmax*Incmax) , Xt(Nmax) , &
     &                 Xx(Nmax*Incmax) , Z(Nmax)
      INTEGER Idim(Nidim) , Inc(Ninc) , Kb(Nkb)
!     .. Local Scalars ..
      DOUBLE PRECISION err , errmax , transl
      INTEGER i , icd , ict , icu , ik , in , incx , incxs , ix , k ,   &
     &        ks , laa , lda , ldas , lx , n , nargs , nc , nk , ns
      LOGICAL banded , full , null , packed , reset , same
      CHARACTER*1 diag , diags , trans , transs , uplo , uplos
      CHARACTER*2 ichd , ichu
      CHARACTER*3 icht
!     .. Local Arrays ..
      LOGICAL isame(13)
!     .. External Functions ..
      LOGICAL LDE , LDERES
      EXTERNAL LDE , LDERES
!     .. External Subroutines ..
      EXTERNAL DMAKE , DMVCH , DTBMV , DTBSV , DTPMV , DTPSV , DTRMV ,  &
     &         DTRSV
!     .. Intrinsic Functions ..
      INTRINSIC ABS , MAX
!     .. Scalars in Common ..
      INTEGER INFot , NOUtc
      LOGICAL LERr , OK
!     .. Common blocks ..
      COMMON /INFOC / INFot , NOUtc , OK , LERr
!     .. Data statements ..
      DATA ichu/'UL'/ , icht/'NTC'/ , ichd/'UN'/
!     .. Executable Statements ..
      full = Sname(3:3)=='R'
      banded = Sname(3:3)=='B'
      packed = Sname(3:3)=='P'
!     Define the number of arguments.
      IF ( full ) THEN
         nargs = 8
      ELSEIF ( banded ) THEN
         nargs = 9
      ELSEIF ( packed ) THEN
         nargs = 7
      ENDIF
!
      nc = 0
      reset = .TRUE.
      errmax = ZERO
!     Set up zero vector for DMVCH.
      DO i = 1 , Nmax
         Z(i) = ZERO
      ENDDO
!
      DO in = 1 , Nidim
         n = Idim(in)
!
         IF ( banded ) THEN
            nk = Nkb
         ELSE
            nk = 1
         ENDIF
         DO ik = 1 , nk
            IF ( banded ) THEN
               k = Kb(ik)
            ELSE
               k = n - 1
            ENDIF
!           Set LDA to 1 more than minimum value if room.
            IF ( banded ) THEN
               lda = k + 1
            ELSE
               lda = n
            ENDIF
            IF ( lda<Nmax ) lda = lda + 1
!           Skip tests if not enough room.
            IF ( lda<=Nmax ) THEN
               IF ( packed ) THEN
                  laa = (n*(n+1))/2
               ELSE
                  laa = lda*n
               ENDIF
               null = n<=0
!
               DO icu = 1 , 2
                  uplo = ichu(icu:icu)
!
                  DO ict = 1 , 3
                     trans = icht(ict:ict)
!
                     DO icd = 1 , 2
                        diag = ichd(icd:icd)
!
!                    Generate the matrix A.
!
                        transl = ZERO
                        CALL DMAKE(Sname(2:3),uplo,diag,n,n,A,Nmax,Aa,  &
     &                             lda,k,k,reset,transl)
!
                        DO ix = 1 , Ninc
                           incx = Inc(ix)
                           lx = ABS(incx)*n
!
!                       Generate the vector X.
!
                           transl = HALF
                           CALL DMAKE('GE',' ',' ',1,n,X,1,Xx,ABS(incx),&
     &                                0,n-1,reset,transl)
                           IF ( n>1 ) THEN
                              X(n/2) = ZERO
                              Xx(1+ABS(incx)*(n/2-1)) = ZERO
                           ENDIF
!
                           nc = nc + 1
!
!                       Save every datum before calling the subroutine.
!
                           uplos = uplo
                           transs = trans
                           diags = diag
                           ns = n
                           ks = k
                           DO i = 1 , laa
                              As(i) = Aa(i)
                           ENDDO
                           ldas = lda
                           DO i = 1 , lx
                              Xs(i) = Xx(i)
                           ENDDO
                           incxs = incx
!
!                       Call the subroutine.
!
                           IF ( Sname(4:5)=='MV' ) THEN
                              IF ( full ) THEN
                                 IF ( Trace ) WRITE (Ntra,FMT=99007)    &
     &                                nc , Sname , uplo , trans , diag ,&
     &                                n , lda , incx
                                 IF ( Rewi ) REWIND Ntra
                                 CALL DTRMV(uplo,trans,diag,n,Aa,lda,Xx,&
     &                              incx)
                              ELSEIF ( banded ) THEN
                                 IF ( Trace ) WRITE (Ntra,FMT=99006)    &
     &                                nc , Sname , uplo , trans , diag ,&
     &                                n , k , lda , incx
                                 IF ( Rewi ) REWIND Ntra
                                 CALL DTBMV(uplo,trans,diag,n,k,Aa,lda, &
     &                              Xx,incx)
                              ELSEIF ( packed ) THEN
                                 IF ( Trace ) WRITE (Ntra,FMT=99005)    &
     &                                nc , Sname , uplo , trans , diag ,&
     &                                n , incx
                                 IF ( Rewi ) REWIND Ntra
                                 CALL DTPMV(uplo,trans,diag,n,Aa,Xx,    &
     &                              incx)
                              ENDIF
                           ELSEIF ( Sname(4:5)=='SV' ) THEN
                              IF ( full ) THEN
                                 IF ( Trace ) WRITE (Ntra,FMT=99007)    &
     &                                nc , Sname , uplo , trans , diag ,&
     &                                n , lda , incx
                                 IF ( Rewi ) REWIND Ntra
                                 CALL DTRSV(uplo,trans,diag,n,Aa,lda,Xx,&
     &                              incx)
                              ELSEIF ( banded ) THEN
                                 IF ( Trace ) WRITE (Ntra,FMT=99006)    &
     &                                nc , Sname , uplo , trans , diag ,&
     &                                n , k , lda , incx
                                 IF ( Rewi ) REWIND Ntra
                                 CALL DTBSV(uplo,trans,diag,n,k,Aa,lda, &
     &                              Xx,incx)
                              ELSEIF ( packed ) THEN
                                 IF ( Trace ) WRITE (Ntra,FMT=99005)    &
     &                                nc , Sname , uplo , trans , diag ,&
     &                                n , incx
                                 IF ( Rewi ) REWIND Ntra
                                 CALL DTPSV(uplo,trans,diag,n,Aa,Xx,    &
     &                              incx)
                              ENDIF
                           ENDIF
!
!                       Check if error-exit was taken incorrectly.
!
                           IF ( .NOT.OK ) THEN
                              WRITE (Nout,FMT=99008)
                              Fatal = .TRUE.
                              GOTO 200
                           ENDIF
!
!                       See what data changed inside subroutines.
!
                           isame(1) = uplo==uplos
                           isame(2) = trans==transs
                           isame(3) = diag==diags
                           isame(4) = ns==n
                           IF ( full ) THEN
                              isame(5) = LDE(As,Aa,laa)
                              isame(6) = ldas==lda
                              IF ( null ) THEN
                                 isame(7) = LDE(Xs,Xx,lx)
                              ELSE
                                 isame(7)                               &
     &                              = LDERES('GE',' ',1,n,Xs,Xx,ABS     &
     &                              (incx))
                              ENDIF
                              isame(8) = incxs==incx
                           ELSEIF ( banded ) THEN
                              isame(5) = ks==k
                              isame(6) = LDE(As,Aa,laa)
                              isame(7) = ldas==lda
                              IF ( null ) THEN
                                 isame(8) = LDE(Xs,Xx,lx)
                              ELSE
                                 isame(8)                               &
     &                              = LDERES('GE',' ',1,n,Xs,Xx,ABS     &
     &                              (incx))
                              ENDIF
                              isame(9) = incxs==incx
                           ELSEIF ( packed ) THEN
                              isame(5) = LDE(As,Aa,laa)
                              IF ( null ) THEN
                                 isame(6) = LDE(Xs,Xx,lx)
                              ELSE
                                 isame(6)                               &
     &                              = LDERES('GE',' ',1,n,Xs,Xx,ABS     &
     &                              (incx))
                              ENDIF
                              isame(7) = incxs==incx
                           ENDIF
!
!                       If data was incorrectly changed, report and
!                       return.
!
                           same = .TRUE.
                           DO i = 1 , nargs
                              same = same .AND. isame(i)
                              IF ( .NOT.isame(i) )                      &
     &                             WRITE (Nout,FMT=99002) i
                           ENDDO
                           IF ( .NOT.same ) THEN
                              Fatal = .TRUE.
                              GOTO 200
                           ENDIF
!
!                          Avoid repeating tests with N.le.0.
                           IF ( null ) GOTO 100
                           IF ( Sname(4:5)=='MV' ) THEN
!
!                             Check the result.
!
                              CALL DMVCH(trans,n,n,ONE,A,Nmax,X,incx,   &
     &                           ZERO,Z,incx,Xt,G,Xx,Eps,err,Fatal,Nout,&
     &                           .TRUE.)
                           ELSEIF ( Sname(4:5)=='SV' ) THEN
!
!                             Compute approximation to original vector.
!
                              DO i = 1 , n
                                 Z(i) = Xx(1+(i-1)*ABS(incx))
                                 Xx(1+(i-1)*ABS(incx)) = X(i)
                              ENDDO
                              CALL DMVCH(trans,n,n,ONE,A,Nmax,Z,incx,   &
     &                           ZERO,X,incx,Xt,G,Xx,Eps,err,Fatal,Nout,&
     &                           .FALSE.)
                           ENDIF
                           errmax = MAX(errmax,err)
!                          If got really bad answer, report and return.
                           IF ( Fatal ) GOTO 200
!
                        ENDDO
!
                     ENDDO
!
                  ENDDO
!
               ENDDO
            ENDIF
!
         ENDDO
!
 100  ENDDO
!
!     Report result.
!
      IF ( errmax<Thresh ) THEN
         WRITE (Nout,FMT=99001) Sname , nc
      ELSE
         WRITE (Nout,FMT=99003) Sname , nc , errmax
      ENDIF
      GOTO 300
!
 200  WRITE (Nout,FMT=99004) Sname
      IF ( full ) THEN
         WRITE (Nout,FMT=99007) nc , Sname , uplo , trans , diag , n ,  &
     &                          lda , incx
      ELSEIF ( banded ) THEN
         WRITE (Nout,FMT=99006) nc , Sname , uplo , trans , diag , n ,  &
     &                          k , lda , incx
      ELSEIF ( packed ) THEN
         WRITE (Nout,FMT=99005) nc , Sname , uplo , trans , diag , n ,  &
     &                          incx
      ENDIF
!
 300  RETURN
!
99001 FORMAT (' ',A6,' PASSED THE COMPUTATIONAL TESTS (',I6,' CALL',    &
     &        'S)')
99002 FORMAT (' ******* FATAL ERROR - PARAMETER NUMBER ',I2,' WAS CH',  &
     &        'ANGED INCORRECTLY *******')
99003 FORMAT (' ',A6,' COMPLETED THE COMPUTATIONAL TESTS (',I6,' C',    &
     &        'ALLS)',/' ******* BUT WITH MAXIMUM TEST RATIO',F8.2,     &
     &        ' - SUSPECT *******')
99004 FORMAT (' ******* ',A6,' FAILED ON CALL NUMBER:')
99005 FORMAT (1X,I6,': ',A6,'(',3('''',A1,''','),I3,', AP, ','X,',I2,   &
     &        ')                        .')
99006 FORMAT (1X,I6,': ',A6,'(',3('''',A1,''','),2(I3,','),' A,',I3,    &
     &        ', X,',I2,')                 .')
99007 FORMAT (1X,I6,': ',A6,'(',3('''',A1,''','),I3,', A,',I3,', X,',I2,&
     &        ')                     .')
99008 FORMAT (' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *',&
     &        '******')
!
!     End of DCHK3.
!
      END SUBROUTINE DCHK3
!*==dchk4.f90  processed by SPAG 7.51RB at 17:52 on  4 Mar 2022
      SUBROUTINE DCHK4(Sname,Eps,Thresh,Nout,Ntra,Trace,Rewi,Fatal,     &
     &                 Nidim,Idim,Nalf,Alf,Ninc,Inc,Nmax,Incmax,A,Aa,As,&
     &                 X,Xx,Xs,Y,Yy,Ys,Yt,G,Z)
      IMPLICIT NONE
!*--DCHK41423
!
!  Tests DGER.
!
!  Auxiliary routine for test program for Level 2 Blas.
!
!  -- Written on 10-August-1987.
!     Richard Hanson, Sandia National Labs.
!     Jeremy Du Croz, NAG Central Office.
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO , HALF , ONE
      PARAMETER (ZERO=0.0D0,HALF=0.5D0,ONE=1.0D0)
!     .. Scalar Arguments ..
      DOUBLE PRECISION Eps , Thresh
      INTEGER Incmax , Nalf , Nidim , Ninc , Nmax , Nout , Ntra
      LOGICAL Fatal , Rewi , Trace
      CHARACTER*6 Sname
!     .. Array Arguments ..
      DOUBLE PRECISION A(Nmax,Nmax) , Aa(Nmax*Nmax) , Alf(Nalf) ,       &
     &                 As(Nmax*Nmax) , G(Nmax) , X(Nmax) ,              &
     &                 Xs(Nmax*Incmax) , Xx(Nmax*Incmax) , Y(Nmax) ,    &
     &                 Ys(Nmax*Incmax) , Yt(Nmax) , Yy(Nmax*Incmax) ,   &
     &                 Z(Nmax)
      INTEGER Idim(Nidim) , Inc(Ninc)
!     .. Local Scalars ..
      DOUBLE PRECISION alpha , als , err , errmax , transl
      INTEGER i , ia , im , in , incx , incxs , incy , incys , ix , iy ,&
     &        j , laa , lda , ldas , lx , ly , m , ms , n , nargs , nc ,&
     &        nd , ns
      LOGICAL null , reset , same
!     .. Local Arrays ..
      DOUBLE PRECISION w(1)
      LOGICAL isame(13)
!     .. External Functions ..
      LOGICAL LDE , LDERES
      EXTERNAL LDE , LDERES
!     .. External Subroutines ..
      EXTERNAL DGER , DMAKE , DMVCH
!     .. Intrinsic Functions ..
      INTRINSIC ABS , MAX , MIN
!     .. Scalars in Common ..
      INTEGER INFot , NOUtc
      LOGICAL LERr , OK
!     .. Common blocks ..
      COMMON /INFOC / INFot , NOUtc , OK , LERr
!     .. Executable Statements ..
!     Define the number of arguments.
      nargs = 9
!
      nc = 0
      reset = .TRUE.
      errmax = ZERO
!
      DO in = 1 , Nidim
         n = Idim(in)
         nd = n/2 + 1
!
         DO im = 1 , 2
            IF ( im==1 ) m = MAX(n-nd,0)
            IF ( im==2 ) m = MIN(n+nd,Nmax)
!
!           Set LDA to 1 more than minimum value if room.
            lda = m
            IF ( lda<Nmax ) lda = lda + 1
!           Skip tests if not enough room.
            IF ( lda<=Nmax ) THEN
               laa = lda*n
               null = n<=0 .OR. m<=0
!
               DO ix = 1 , Ninc
                  incx = Inc(ix)
                  lx = ABS(incx)*m
!
!              Generate the vector X.
!
                  transl = HALF
                  CALL DMAKE('GE',' ',' ',1,m,X,1,Xx,ABS(incx),0,m-1,   &
     &                       reset,transl)
                  IF ( m>1 ) THEN
                     X(m/2) = ZERO
                     Xx(1+ABS(incx)*(m/2-1)) = ZERO
                  ENDIF
!
                  DO iy = 1 , Ninc
                     incy = Inc(iy)
                     ly = ABS(incy)*n
!
!                 Generate the vector Y.
!
                     transl = ZERO
                     CALL DMAKE('GE',' ',' ',1,n,Y,1,Yy,ABS(incy),0,n-1,&
     &                          reset,transl)
                     IF ( n>1 ) THEN
                        Y(n/2) = ZERO
                        Yy(1+ABS(incy)*(n/2-1)) = ZERO
                     ENDIF
!
                     DO ia = 1 , Nalf
                        alpha = Alf(ia)
!
!                    Generate the matrix A.
!
                        transl = ZERO
                        CALL DMAKE(Sname(2:3),' ',' ',m,n,A,Nmax,Aa,lda,&
     &                             m-1,n-1,reset,transl)
!
                        nc = nc + 1
!
!                    Save every datum before calling the subroutine.
!
                        ms = m
                        ns = n
                        als = alpha
                        DO i = 1 , laa
                           As(i) = Aa(i)
                        ENDDO
                        ldas = lda
                        DO i = 1 , lx
                           Xs(i) = Xx(i)
                        ENDDO
                        incxs = incx
                        DO i = 1 , ly
                           Ys(i) = Yy(i)
                        ENDDO
                        incys = incy
!
!                    Call the subroutine.
!
                        IF ( Trace ) WRITE (Ntra,FMT=99006) nc , Sname ,&
     &                       m , n , alpha , incx , incy , lda
                        IF ( Rewi ) REWIND Ntra
                        CALL DGER(m,n,alpha,Xx,incx,Yy,incy,Aa,lda)
!
!                    Check if error-exit was taken incorrectly.
!
                        IF ( .NOT.OK ) THEN
                           WRITE (Nout,FMT=99007)
                           Fatal = .TRUE.
                           GOTO 200
                        ENDIF
!
!                    See what data changed inside subroutine.
!
                        isame(1) = ms==m
                        isame(2) = ns==n
                        isame(3) = als==alpha
                        isame(4) = LDE(Xs,Xx,lx)
                        isame(5) = incxs==incx
                        isame(6) = LDE(Ys,Yy,ly)
                        isame(7) = incys==incy
                        IF ( null ) THEN
                           isame(8) = LDE(As,Aa,laa)
                        ELSE
                           isame(8) = LDERES('GE',' ',m,n,As,Aa,lda)
                        ENDIF
                        isame(9) = ldas==lda
!
!                    If data was incorrectly changed, report and return.
!
                        same = .TRUE.
                        DO i = 1 , nargs
                           same = same .AND. isame(i)
                           IF ( .NOT.isame(i) ) WRITE (Nout,FMT=99002) i
                        ENDDO
                        IF ( .NOT.same ) THEN
                           Fatal = .TRUE.
                           GOTO 200
                        ENDIF
!
!                       Avoid repeating tests with M.le.0 or N.le.0.
                        IF ( null ) GOTO 50
!
!                       Check the result column by column.
!
                        IF ( incx>0 ) THEN
                           DO i = 1 , m
                              Z(i) = X(i)
                           ENDDO
                        ELSE
                           DO i = 1 , m
                              Z(i) = X(m-i+1)
                           ENDDO
                        ENDIF
                        DO j = 1 , n
                           IF ( incy>0 ) THEN
                              w(1) = Y(j)
                           ELSE
                              w(1) = Y(n-j+1)
                           ENDIF
                           CALL DMVCH('N',m,1,alpha,Z,Nmax,w,1,ONE,     &
     &                                A(1,j),1,Yt,G,Aa(1+(j-1)*lda),Eps,&
     &                                err,Fatal,Nout,.TRUE.)
                           errmax = MAX(errmax,err)
!                          If got really bad answer, report and return.
                           IF ( Fatal ) GOTO 100
                        ENDDO
!
                     ENDDO
!
                  ENDDO
!
               ENDDO
            ENDIF
!
 50      ENDDO
!
      ENDDO
!
!     Report result.
!
      IF ( errmax<Thresh ) THEN
         WRITE (Nout,FMT=99001) Sname , nc
      ELSE
         WRITE (Nout,FMT=99003) Sname , nc , errmax
      ENDIF
      GOTO 300
!
 100  WRITE (Nout,FMT=99005) j
!
 200  WRITE (Nout,FMT=99004) Sname
      WRITE (Nout,FMT=99006) nc , Sname , m , n , alpha , incx , incy , &
     &                       lda
!
 300  RETURN
!
99001 FORMAT (' ',A6,' PASSED THE COMPUTATIONAL TESTS (',I6,' CALL',    &
     &        'S)')
99002 FORMAT (' ******* FATAL ERROR - PARAMETER NUMBER ',I2,' WAS CH',  &
     &        'ANGED INCORRECTLY *******')
99003 FORMAT (' ',A6,' COMPLETED THE COMPUTATIONAL TESTS (',I6,' C',    &
     &        'ALLS)',/' ******* BUT WITH MAXIMUM TEST RATIO',F8.2,     &
     &        ' - SUSPECT *******')
99004 FORMAT (' ******* ',A6,' FAILED ON CALL NUMBER:')
99005 FORMAT ('      THESE ARE THE RESULTS FOR COLUMN ',I3)
99006 FORMAT (1X,I6,': ',A6,'(',2(I3,','),F4.1,', X,',I2,', Y,',I2,     &
     &        ', A,',I3,')                  .')
99007 FORMAT (' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *',&
     &        '******')
!
!     End of DCHK4.
!
      END SUBROUTINE DCHK4
!*==dchk5.f90  processed by SPAG 7.51RB at 17:52 on  4 Mar 2022
      SUBROUTINE DCHK5(Sname,Eps,Thresh,Nout,Ntra,Trace,Rewi,Fatal,     &
     &                 Nidim,Idim,Nalf,Alf,Ninc,Inc,Nmax,Incmax,A,Aa,As,&
     &                 X,Xx,Xs,Y,Yy,Ys,Yt,G,Z)
      IMPLICIT NONE
!*--DCHK51671
!
!  Tests DSYR and DSPR.
!
!  Auxiliary routine for test program for Level 2 Blas.
!
!  -- Written on 10-August-1987.
!     Richard Hanson, Sandia National Labs.
!     Jeremy Du Croz, NAG Central Office.
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO , HALF , ONE
      PARAMETER (ZERO=0.0D0,HALF=0.5D0,ONE=1.0D0)
!     .. Scalar Arguments ..
      DOUBLE PRECISION Eps , Thresh
      INTEGER Incmax , Nalf , Nidim , Ninc , Nmax , Nout , Ntra
      LOGICAL Fatal , Rewi , Trace
      CHARACTER*6 Sname
!     .. Array Arguments ..
      DOUBLE PRECISION A(Nmax,Nmax) , Aa(Nmax*Nmax) , Alf(Nalf) ,       &
     &                 As(Nmax*Nmax) , G(Nmax) , X(Nmax) ,              &
     &                 Xs(Nmax*Incmax) , Xx(Nmax*Incmax) , Y(Nmax) ,    &
     &                 Ys(Nmax*Incmax) , Yt(Nmax) , Yy(Nmax*Incmax) ,   &
     &                 Z(Nmax)
      INTEGER Idim(Nidim) , Inc(Ninc)
!     .. Local Scalars ..
      DOUBLE PRECISION alpha , als , err , errmax , transl
      INTEGER i , ia , ic , in , incx , incxs , ix , j , ja , jj , laa ,&
     &        lda , ldas , lj , lx , n , nargs , nc , ns
      LOGICAL full , null , packed , reset , same , upper
      CHARACTER*1 uplo , uplos
      CHARACTER*2 ich
!     .. Local Arrays ..
      DOUBLE PRECISION w(1)
      LOGICAL isame(13)
!     .. External Functions ..
      LOGICAL LDE , LDERES
      EXTERNAL LDE , LDERES
!     .. External Subroutines ..
      EXTERNAL DMAKE , DMVCH , DSPR , DSYR
!     .. Intrinsic Functions ..
      INTRINSIC ABS , MAX
!     .. Scalars in Common ..
      INTEGER INFot , NOUtc
      LOGICAL LERr , OK
!     .. Common blocks ..
      COMMON /INFOC / INFot , NOUtc , OK , LERr
!     .. Data statements ..
      DATA ich/'UL'/
!     .. Executable Statements ..
      full = Sname(3:3)=='Y'
      packed = Sname(3:3)=='P'
!     Define the number of arguments.
      IF ( full ) THEN
         nargs = 7
      ELSEIF ( packed ) THEN
         nargs = 6
      ENDIF
!
      nc = 0
      reset = .TRUE.
      errmax = ZERO
!
      DO in = 1 , Nidim
         n = Idim(in)
!        Set LDA to 1 more than minimum value if room.
         lda = n
         IF ( lda<Nmax ) lda = lda + 1
!        Skip tests if not enough room.
         IF ( lda<=Nmax ) THEN
            IF ( packed ) THEN
               laa = (n*(n+1))/2
            ELSE
               laa = lda*n
            ENDIF
!
            DO ic = 1 , 2
               uplo = ich(ic:ic)
               upper = uplo=='U'
!
               DO ix = 1 , Ninc
                  incx = Inc(ix)
                  lx = ABS(incx)*n
!
!              Generate the vector X.
!
                  transl = HALF
                  CALL DMAKE('GE',' ',' ',1,n,X,1,Xx,ABS(incx),0,n-1,   &
     &                       reset,transl)
                  IF ( n>1 ) THEN
                     X(n/2) = ZERO
                     Xx(1+ABS(incx)*(n/2-1)) = ZERO
                  ENDIF
!
                  DO ia = 1 , Nalf
                     alpha = Alf(ia)
                     null = n<=0 .OR. alpha==ZERO
!
!                 Generate the matrix A.
!
                     transl = ZERO
                     CALL DMAKE(Sname(2:3),uplo,' ',n,n,A,Nmax,Aa,lda,  &
     &                          n-1,n-1,reset,transl)
!
                     nc = nc + 1
!
!                 Save every datum before calling the subroutine.
!
                     uplos = uplo
                     ns = n
                     als = alpha
                     DO i = 1 , laa
                        As(i) = Aa(i)
                     ENDDO
                     ldas = lda
                     DO i = 1 , lx
                        Xs(i) = Xx(i)
                     ENDDO
                     incxs = incx
!
!                 Call the subroutine.
!
                     IF ( full ) THEN
                        IF ( Trace ) WRITE (Ntra,FMT=99007) nc , Sname ,&
     &                       uplo , n , alpha , incx , lda
                        IF ( Rewi ) REWIND Ntra
                        CALL DSYR(uplo,n,alpha,Xx,incx,Aa,lda)
                     ELSEIF ( packed ) THEN
                        IF ( Trace ) WRITE (Ntra,FMT=99006) nc , Sname ,&
     &                       uplo , n , alpha , incx
                        IF ( Rewi ) REWIND Ntra
                        CALL DSPR(uplo,n,alpha,Xx,incx,Aa)
                     ENDIF
!
!                 Check if error-exit was taken incorrectly.
!
                     IF ( .NOT.OK ) THEN
                        WRITE (Nout,FMT=99008)
                        Fatal = .TRUE.
                        GOTO 300
                     ENDIF
!
!                 See what data changed inside subroutines.
!
                     isame(1) = uplo==uplos
                     isame(2) = ns==n
                     isame(3) = als==alpha
                     isame(4) = LDE(Xs,Xx,lx)
                     isame(5) = incxs==incx
                     IF ( null ) THEN
                        isame(6) = LDE(As,Aa,laa)
                     ELSE
                        isame(6) = LDERES(Sname(2:3),uplo,n,n,As,Aa,lda)
                     ENDIF
                     IF ( .NOT.packed ) isame(7) = ldas==lda
!
!                 If data was incorrectly changed, report and return.
!
                     same = .TRUE.
                     DO i = 1 , nargs
                        same = same .AND. isame(i)
                        IF ( .NOT.isame(i) ) WRITE (Nout,FMT=99002) i
                     ENDDO
                     IF ( .NOT.same ) THEN
                        Fatal = .TRUE.
                        GOTO 300
                     ENDIF
!
                     IF ( .NOT.null ) THEN
!
!                    Check the result column by column.
!
                        IF ( incx>0 ) THEN
                           DO i = 1 , n
                              Z(i) = X(i)
                           ENDDO
                        ELSE
                           DO i = 1 , n
                              Z(i) = X(n-i+1)
                           ENDDO
                        ENDIF
                        ja = 1
                        DO j = 1 , n
                           w(1) = Z(j)
                           IF ( upper ) THEN
                              jj = 1
                              lj = j
                           ELSE
                              jj = j
                              lj = n - j + 1
                           ENDIF
                           CALL DMVCH('N',lj,1,alpha,Z(jj),lj,w,1,ONE,  &
     &                                A(jj,j),1,Yt,G,Aa(ja),Eps,err,    &
     &                                Fatal,Nout,.TRUE.)
                           IF ( .NOT.(full) ) THEN
                              ja = ja + lj
                           ELSEIF ( upper ) THEN
                              ja = ja + lda
                           ELSE
                              ja = ja + lda + 1
                           ENDIF
                           errmax = MAX(errmax,err)
!                       If got really bad answer, report and return.
                           IF ( Fatal ) GOTO 200
                        ENDDO
!                    Avoid repeating tests if N.le.0.
                     ELSEIF ( n<=0 ) THEN
                        GOTO 100
                     ENDIF
!
                  ENDDO
!
               ENDDO
!
            ENDDO
         ENDIF
!
 100  ENDDO
!
!     Report result.
!
      IF ( errmax<Thresh ) THEN
         WRITE (Nout,FMT=99001) Sname , nc
      ELSE
         WRITE (Nout,FMT=99003) Sname , nc , errmax
      ENDIF
      GOTO 400
!
 200  WRITE (Nout,FMT=99005) j
!
 300  WRITE (Nout,FMT=99004) Sname
      IF ( full ) THEN
         WRITE (Nout,FMT=99007) nc , Sname , uplo , n , alpha , incx ,  &
     &                          lda
      ELSEIF ( packed ) THEN
         WRITE (Nout,FMT=99006) nc , Sname , uplo , n , alpha , incx
      ENDIF
!
 400  RETURN
!
99001 FORMAT (' ',A6,' PASSED THE COMPUTATIONAL TESTS (',I6,' CALL',    &
     &        'S)')
99002 FORMAT (' ******* FATAL ERROR - PARAMETER NUMBER ',I2,' WAS CH',  &
     &        'ANGED INCORRECTLY *******')
99003 FORMAT (' ',A6,' COMPLETED THE COMPUTATIONAL TESTS (',I6,' C',    &
     &        'ALLS)',/' ******* BUT WITH MAXIMUM TEST RATIO',F8.2,     &
     &        ' - SUSPECT *******')
99004 FORMAT (' ******* ',A6,' FAILED ON CALL NUMBER:')
99005 FORMAT ('      THESE ARE THE RESULTS FOR COLUMN ',I3)
99006 FORMAT (1X,I6,': ',A6,'(''',A1,''',',I3,',',F4.1,', X,',I2,       &
     &        ', AP)                           .')
99007 FORMAT (1X,I6,': ',A6,'(''',A1,''',',I3,',',F4.1,', X,',I2,', A,',&
     &        I3,')                        .')
99008 FORMAT (' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *',&
     &        '******')
!
!     End of DCHK5.
!
      END SUBROUTINE DCHK5
!*==dchk6.f90  processed by SPAG 7.51RB at 17:52 on  4 Mar 2022
      SUBROUTINE DCHK6(Sname,Eps,Thresh,Nout,Ntra,Trace,Rewi,Fatal,     &
     &                 Nidim,Idim,Nalf,Alf,Ninc,Inc,Nmax,Incmax,A,Aa,As,&
     &                 X,Xx,Xs,Y,Yy,Ys,Yt,G,Z)
      IMPLICIT NONE
!*--DCHK61935
!
!  Tests DSYR2 and DSPR2.
!
!  Auxiliary routine for test program for Level 2 Blas.
!
!  -- Written on 10-August-1987.
!     Richard Hanson, Sandia National Labs.
!     Jeremy Du Croz, NAG Central Office.
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO , HALF , ONE
      PARAMETER (ZERO=0.0D0,HALF=0.5D0,ONE=1.0D0)
!     .. Scalar Arguments ..
      DOUBLE PRECISION Eps , Thresh
      INTEGER Incmax , Nalf , Nidim , Ninc , Nmax , Nout , Ntra
      LOGICAL Fatal , Rewi , Trace
      CHARACTER*6 Sname
!     .. Array Arguments ..
      DOUBLE PRECISION A(Nmax,Nmax) , Aa(Nmax*Nmax) , Alf(Nalf) ,       &
     &                 As(Nmax*Nmax) , G(Nmax) , X(Nmax) ,              &
     &                 Xs(Nmax*Incmax) , Xx(Nmax*Incmax) , Y(Nmax) ,    &
     &                 Ys(Nmax*Incmax) , Yt(Nmax) , Yy(Nmax*Incmax) ,   &
     &                 Z(Nmax,2)
      INTEGER Idim(Nidim) , Inc(Ninc)
!     .. Local Scalars ..
      DOUBLE PRECISION alpha , als , err , errmax , transl
      INTEGER i , ia , ic , in , incx , incxs , incy , incys , ix , iy ,&
     &        j , ja , jj , laa , lda , ldas , lj , lx , ly , n ,       &
     &        nargs , nc , ns
      LOGICAL full , null , packed , reset , same , upper
      CHARACTER*1 uplo , uplos
      CHARACTER*2 ich
!     .. Local Arrays ..
      DOUBLE PRECISION w(2)
      LOGICAL isame(13)
!     .. External Functions ..
      LOGICAL LDE , LDERES
      EXTERNAL LDE , LDERES
!     .. External Subroutines ..
      EXTERNAL DMAKE , DMVCH , DSPR2 , DSYR2
!     .. Intrinsic Functions ..
      INTRINSIC ABS , MAX
!     .. Scalars in Common ..
      INTEGER INFot , NOUtc
      LOGICAL LERr , OK
!     .. Common blocks ..
      COMMON /INFOC / INFot , NOUtc , OK , LERr
!     .. Data statements ..
      DATA ich/'UL'/
!     .. Executable Statements ..
      full = Sname(3:3)=='Y'
      packed = Sname(3:3)=='P'
!     Define the number of arguments.
      IF ( full ) THEN
         nargs = 9
      ELSEIF ( packed ) THEN
         nargs = 8
      ENDIF
!
      nc = 0
      reset = .TRUE.
      errmax = ZERO
!
      DO in = 1 , Nidim
         n = Idim(in)
!        Set LDA to 1 more than minimum value if room.
         lda = n
         IF ( lda<Nmax ) lda = lda + 1
!        Skip tests if not enough room.
         IF ( lda<=Nmax ) THEN
            IF ( packed ) THEN
               laa = (n*(n+1))/2
            ELSE
               laa = lda*n
            ENDIF
!
            DO ic = 1 , 2
               uplo = ich(ic:ic)
               upper = uplo=='U'
!
               DO ix = 1 , Ninc
                  incx = Inc(ix)
                  lx = ABS(incx)*n
!
!              Generate the vector X.
!
                  transl = HALF
                  CALL DMAKE('GE',' ',' ',1,n,X,1,Xx,ABS(incx),0,n-1,   &
     &                       reset,transl)
                  IF ( n>1 ) THEN
                     X(n/2) = ZERO
                     Xx(1+ABS(incx)*(n/2-1)) = ZERO
                  ENDIF
!
                  DO iy = 1 , Ninc
                     incy = Inc(iy)
                     ly = ABS(incy)*n
!
!                 Generate the vector Y.
!
                     transl = ZERO
                     CALL DMAKE('GE',' ',' ',1,n,Y,1,Yy,ABS(incy),0,n-1,&
     &                          reset,transl)
                     IF ( n>1 ) THEN
                        Y(n/2) = ZERO
                        Yy(1+ABS(incy)*(n/2-1)) = ZERO
                     ENDIF
!
                     DO ia = 1 , Nalf
                        alpha = Alf(ia)
                        null = n<=0 .OR. alpha==ZERO
!
!                    Generate the matrix A.
!
                        transl = ZERO
                        CALL DMAKE(Sname(2:3),uplo,' ',n,n,A,Nmax,Aa,   &
     &                             lda,n-1,n-1,reset,transl)
!
                        nc = nc + 1
!
!                    Save every datum before calling the subroutine.
!
                        uplos = uplo
                        ns = n
                        als = alpha
                        DO i = 1 , laa
                           As(i) = Aa(i)
                        ENDDO
                        ldas = lda
                        DO i = 1 , lx
                           Xs(i) = Xx(i)
                        ENDDO
                        incxs = incx
                        DO i = 1 , ly
                           Ys(i) = Yy(i)
                        ENDDO
                        incys = incy
!
!                    Call the subroutine.
!
                        IF ( full ) THEN
                           IF ( Trace ) WRITE (Ntra,FMT=99007) nc ,     &
     &                          Sname , uplo , n , alpha , incx , incy ,&
     &                          lda
                           IF ( Rewi ) REWIND Ntra
                           CALL DSYR2(uplo,n,alpha,Xx,incx,Yy,incy,Aa,  &
     &                                lda)
                        ELSEIF ( packed ) THEN
                           IF ( Trace ) WRITE (Ntra,FMT=99006) nc ,     &
     &                          Sname , uplo , n , alpha , incx , incy
                           IF ( Rewi ) REWIND Ntra
                           CALL DSPR2(uplo,n,alpha,Xx,incx,Yy,incy,Aa)
                        ENDIF
!
!                    Check if error-exit was taken incorrectly.
!
                        IF ( .NOT.OK ) THEN
                           WRITE (Nout,FMT=99008)
                           Fatal = .TRUE.
                           GOTO 300
                        ENDIF
!
!                    See what data changed inside subroutines.
!
                        isame(1) = uplo==uplos
                        isame(2) = ns==n
                        isame(3) = als==alpha
                        isame(4) = LDE(Xs,Xx,lx)
                        isame(5) = incxs==incx
                        isame(6) = LDE(Ys,Yy,ly)
                        isame(7) = incys==incy
                        IF ( null ) THEN
                           isame(8) = LDE(As,Aa,laa)
                        ELSE
                           isame(8) = LDERES(Sname(2:3),uplo,n,n,As,Aa, &
     &                                lda)
                        ENDIF
                        IF ( .NOT.packed ) isame(9) = ldas==lda
!
!                    If data was incorrectly changed, report and return.
!
                        same = .TRUE.
                        DO i = 1 , nargs
                           same = same .AND. isame(i)
                           IF ( .NOT.isame(i) ) WRITE (Nout,FMT=99002) i
                        ENDDO
                        IF ( .NOT.same ) THEN
                           Fatal = .TRUE.
                           GOTO 300
                        ENDIF
!
                        IF ( .NOT.null ) THEN
!
!                       Check the result column by column.
!
                           IF ( incx>0 ) THEN
                              DO i = 1 , n
                                 Z(i,1) = X(i)
                              ENDDO
                           ELSE
                              DO i = 1 , n
                                 Z(i,1) = X(n-i+1)
                              ENDDO
                           ENDIF
                           IF ( incy>0 ) THEN
                              DO i = 1 , n
                                 Z(i,2) = Y(i)
                              ENDDO
                           ELSE
                              DO i = 1 , n
                                 Z(i,2) = Y(n-i+1)
                              ENDDO
                           ENDIF
                           ja = 1
                           DO j = 1 , n
                              w(1) = Z(j,2)
                              w(2) = Z(j,1)
                              IF ( upper ) THEN
                                 jj = 1
                                 lj = j
                              ELSE
                                 jj = j
                                 lj = n - j + 1
                              ENDIF
                              CALL DMVCH('N',lj,2,alpha,Z(jj,1),Nmax,w, &
     &                           1,ONE,A(jj,j),1,Yt,G,Aa(ja),Eps,err,   &
     &                           Fatal,Nout,.TRUE.)
                              IF ( .NOT.(full) ) THEN
                                 ja = ja + lj
                              ELSEIF ( upper ) THEN
                                 ja = ja + lda
                              ELSE
                                 ja = ja + lda + 1
                              ENDIF
                              errmax = MAX(errmax,err)
!                          If got really bad answer, report and return.
                              IF ( Fatal ) GOTO 200
                           ENDDO
!                       Avoid repeating tests with N.le.0.
                        ELSEIF ( n<=0 ) THEN
                           GOTO 100
                        ENDIF
!
                     ENDDO
!
                  ENDDO
!
               ENDDO
!
            ENDDO
         ENDIF
!
 100  ENDDO
!
!     Report result.
!
      IF ( errmax<Thresh ) THEN
         WRITE (Nout,FMT=99001) Sname , nc
      ELSE
         WRITE (Nout,FMT=99003) Sname , nc , errmax
      ENDIF
      GOTO 400
!
 200  WRITE (Nout,FMT=99005) j
!
 300  WRITE (Nout,FMT=99004) Sname
      IF ( full ) THEN
         WRITE (Nout,FMT=99007) nc , Sname , uplo , n , alpha , incx ,  &
     &                          incy , lda
      ELSEIF ( packed ) THEN
         WRITE (Nout,FMT=99006) nc , Sname , uplo , n , alpha , incx ,  &
     &                          incy
      ENDIF
!
 400  RETURN
!
99001 FORMAT (' ',A6,' PASSED THE COMPUTATIONAL TESTS (',I6,' CALL',    &
     &        'S)')
99002 FORMAT (' ******* FATAL ERROR - PARAMETER NUMBER ',I2,' WAS CH',  &
     &        'ANGED INCORRECTLY *******')
99003 FORMAT (' ',A6,' COMPLETED THE COMPUTATIONAL TESTS (',I6,' C',    &
     &        'ALLS)',/' ******* BUT WITH MAXIMUM TEST RATIO',F8.2,     &
     &        ' - SUSPECT *******')
99004 FORMAT (' ******* ',A6,' FAILED ON CALL NUMBER:')
99005 FORMAT ('      THESE ARE THE RESULTS FOR COLUMN ',I3)
99006 FORMAT (1X,I6,': ',A6,'(''',A1,''',',I3,',',F4.1,', X,',I2,', Y,',&
     &        I2,', AP)                     .')
99007 FORMAT (1X,I6,': ',A6,'(''',A1,''',',I3,',',F4.1,', X,',I2,', Y,',&
     &        I2,', A,',I3,')                  .')
99008 FORMAT (' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *',&
     &        '******')
!
!     End of DCHK6.
!
      END SUBROUTINE DCHK6
!*==dchke.f90  processed by SPAG 7.51RB at 17:52 on  4 Mar 2022
      SUBROUTINE DCHKE(Isnum,Srnamt,Nout)
      IMPLICIT NONE
!*--DCHKE2234
!
!  Tests the error exits from the Level 2 Blas.
!  Requires a special version of the error-handling routine XERBLA.
!  ALPHA, BETA, A, X and Y should not need to be defined.
!
!  Auxiliary routine for test program for Level 2 Blas.
!
!  -- Written on 10-August-1987.
!     Richard Hanson, Sandia National Labs.
!     Jeremy Du Croz, NAG Central Office.
!
!     .. Scalar Arguments ..
      INTEGER Isnum , Nout
      CHARACTER*6 Srnamt
!     .. Scalars in Common ..
      INTEGER INFot , NOUtc
      LOGICAL LERr , OK
!     .. Local Scalars ..
      DOUBLE PRECISION alpha , beta
!     .. Local Arrays ..
      DOUBLE PRECISION a(1,1) , x(1) , y(1)
!     .. External Subroutines ..
      EXTERNAL CHKXER , DGBMV , DGEMV , DGER , DSBMV , DSPMV , DSPR ,   &
     &         DSPR2 , DSYMV , DSYR , DSYR2 , DTBMV , DTBSV , DTPMV ,   &
     &         DTPSV , DTRMV , DTRSV
!     .. Common blocks ..
      COMMON /INFOC / INFot , NOUtc , OK , LERr
!     .. Executable Statements ..
!     OK is set to .FALSE. by the special version of XERBLA or by CHKXER
!     if anything is wrong.
      OK = .TRUE.
!     LERR is set to .TRUE. by the special version of XERBLA each time
!     it is called, and is then tested and re-set by CHKXER.
      LERr = .FALSE.
      IF ( Isnum==2 ) THEN
         INFot = 1
         CALL DGBMV('/',0,0,0,0,alpha,a,1,x,1,beta,y,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 2
         CALL DGBMV('N',-1,0,0,0,alpha,a,1,x,1,beta,y,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 3
         CALL DGBMV('N',0,-1,0,0,alpha,a,1,x,1,beta,y,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 4
         CALL DGBMV('N',0,0,-1,0,alpha,a,1,x,1,beta,y,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 5
         CALL DGBMV('N',2,0,0,-1,alpha,a,1,x,1,beta,y,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 8
         CALL DGBMV('N',0,0,1,0,alpha,a,1,x,1,beta,y,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 10
         CALL DGBMV('N',0,0,0,0,alpha,a,1,x,0,beta,y,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 13
         CALL DGBMV('N',0,0,0,0,alpha,a,1,x,1,beta,y,0)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
      ELSEIF ( Isnum==3 ) THEN
         INFot = 1
         CALL DSYMV('/',0,alpha,a,1,x,1,beta,y,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 2
         CALL DSYMV('U',-1,alpha,a,1,x,1,beta,y,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 5
         CALL DSYMV('U',2,alpha,a,1,x,1,beta,y,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 7
         CALL DSYMV('U',0,alpha,a,1,x,0,beta,y,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 10
         CALL DSYMV('U',0,alpha,a,1,x,1,beta,y,0)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
      ELSEIF ( Isnum==4 ) THEN
         INFot = 1
         CALL DSBMV('/',0,0,alpha,a,1,x,1,beta,y,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 2
         CALL DSBMV('U',-1,0,alpha,a,1,x,1,beta,y,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 3
         CALL DSBMV('U',0,-1,alpha,a,1,x,1,beta,y,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 6
         CALL DSBMV('U',0,1,alpha,a,1,x,1,beta,y,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 8
         CALL DSBMV('U',0,0,alpha,a,1,x,0,beta,y,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 11
         CALL DSBMV('U',0,0,alpha,a,1,x,1,beta,y,0)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
      ELSEIF ( Isnum==5 ) THEN
         INFot = 1
         CALL DSPMV('/',0,alpha,a,x,1,beta,y,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 2
         CALL DSPMV('U',-1,alpha,a,x,1,beta,y,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 6
         CALL DSPMV('U',0,alpha,a,x,0,beta,y,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 9
         CALL DSPMV('U',0,alpha,a,x,1,beta,y,0)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
      ELSEIF ( Isnum==6 ) THEN
         INFot = 1
         CALL DTRMV('/','N','N',0,a,1,x,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 2
         CALL DTRMV('U','/','N',0,a,1,x,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 3
         CALL DTRMV('U','N','/',0,a,1,x,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 4
         CALL DTRMV('U','N','N',-1,a,1,x,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 6
         CALL DTRMV('U','N','N',2,a,1,x,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 8
         CALL DTRMV('U','N','N',0,a,1,x,0)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
      ELSEIF ( Isnum==7 ) THEN
         INFot = 1
         CALL DTBMV('/','N','N',0,0,a,1,x,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 2
         CALL DTBMV('U','/','N',0,0,a,1,x,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 3
         CALL DTBMV('U','N','/',0,0,a,1,x,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 4
         CALL DTBMV('U','N','N',-1,0,a,1,x,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 5
         CALL DTBMV('U','N','N',0,-1,a,1,x,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 7
         CALL DTBMV('U','N','N',0,1,a,1,x,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 9
         CALL DTBMV('U','N','N',0,0,a,1,x,0)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
      ELSEIF ( Isnum==8 ) THEN
         INFot = 1
         CALL DTPMV('/','N','N',0,a,x,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 2
         CALL DTPMV('U','/','N',0,a,x,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 3
         CALL DTPMV('U','N','/',0,a,x,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 4
         CALL DTPMV('U','N','N',-1,a,x,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 7
         CALL DTPMV('U','N','N',0,a,x,0)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
      ELSEIF ( Isnum==9 ) THEN
         INFot = 1
         CALL DTRSV('/','N','N',0,a,1,x,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 2
         CALL DTRSV('U','/','N',0,a,1,x,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 3
         CALL DTRSV('U','N','/',0,a,1,x,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 4
         CALL DTRSV('U','N','N',-1,a,1,x,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 6
         CALL DTRSV('U','N','N',2,a,1,x,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 8
         CALL DTRSV('U','N','N',0,a,1,x,0)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
      ELSEIF ( Isnum==10 ) THEN
         INFot = 1
         CALL DTBSV('/','N','N',0,0,a,1,x,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 2
         CALL DTBSV('U','/','N',0,0,a,1,x,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 3
         CALL DTBSV('U','N','/',0,0,a,1,x,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 4
         CALL DTBSV('U','N','N',-1,0,a,1,x,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 5
         CALL DTBSV('U','N','N',0,-1,a,1,x,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 7
         CALL DTBSV('U','N','N',0,1,a,1,x,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 9
         CALL DTBSV('U','N','N',0,0,a,1,x,0)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
      ELSEIF ( Isnum==11 ) THEN
         INFot = 1
         CALL DTPSV('/','N','N',0,a,x,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 2
         CALL DTPSV('U','/','N',0,a,x,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 3
         CALL DTPSV('U','N','/',0,a,x,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 4
         CALL DTPSV('U','N','N',-1,a,x,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 7
         CALL DTPSV('U','N','N',0,a,x,0)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
      ELSEIF ( Isnum==12 ) THEN
         INFot = 1
         CALL DGER(-1,0,alpha,x,1,y,1,a,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 2
         CALL DGER(0,-1,alpha,x,1,y,1,a,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 5
         CALL DGER(0,0,alpha,x,0,y,1,a,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 7
         CALL DGER(0,0,alpha,x,1,y,0,a,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 9
         CALL DGER(2,0,alpha,x,1,y,1,a,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
      ELSEIF ( Isnum==13 ) THEN
         INFot = 1
         CALL DSYR('/',0,alpha,x,1,a,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 2
         CALL DSYR('U',-1,alpha,x,1,a,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 5
         CALL DSYR('U',0,alpha,x,0,a,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 7
         CALL DSYR('U',2,alpha,x,1,a,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
      ELSEIF ( Isnum==14 ) THEN
         INFot = 1
         CALL DSPR('/',0,alpha,x,1,a)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 2
         CALL DSPR('U',-1,alpha,x,1,a)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 5
         CALL DSPR('U',0,alpha,x,0,a)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
      ELSEIF ( Isnum==15 ) THEN
         INFot = 1
         CALL DSYR2('/',0,alpha,x,1,y,1,a,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 2
         CALL DSYR2('U',-1,alpha,x,1,y,1,a,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 5
         CALL DSYR2('U',0,alpha,x,0,y,1,a,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 7
         CALL DSYR2('U',0,alpha,x,1,y,0,a,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 9
         CALL DSYR2('U',2,alpha,x,1,y,1,a,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
      ELSEIF ( Isnum==16 ) THEN
         INFot = 1
         CALL DSPR2('/',0,alpha,x,1,y,1,a)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 2
         CALL DSPR2('U',-1,alpha,x,1,y,1,a)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 5
         CALL DSPR2('U',0,alpha,x,0,y,1,a)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 7
         CALL DSPR2('U',0,alpha,x,1,y,0,a)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
      ELSE
         INFot = 1
         CALL DGEMV('/',0,0,alpha,a,1,x,1,beta,y,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 2
         CALL DGEMV('N',-1,0,alpha,a,1,x,1,beta,y,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 3
         CALL DGEMV('N',0,-1,alpha,a,1,x,1,beta,y,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 6
         CALL DGEMV('N',2,0,alpha,a,1,x,1,beta,y,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 8
         CALL DGEMV('N',0,0,alpha,a,1,x,0,beta,y,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 11
         CALL DGEMV('N',0,0,alpha,a,1,x,1,beta,y,0)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
      ENDIF
!
      IF ( OK ) THEN
         WRITE (Nout,FMT=99001) Srnamt
      ELSE
         WRITE (Nout,FMT=99002) Srnamt
      ENDIF
      RETURN
!
99001 FORMAT (' ',A6,' PASSED THE TESTS OF ERROR-EXITS')
99002 FORMAT (' ******* ',A6,' FAILED THE TESTS OF ERROR-EXITS *****',  &
     &        '**')
!
!     End of DCHKE.
!
      END SUBROUTINE DCHKE
!*==dmake.f90  processed by SPAG 7.51RB at 17:52 on  4 Mar 2022
      SUBROUTINE DMAKE(Type,Uplo,Diag,M,N,A,Nmax,Aa,Lda,Kl,Ku,Reset,    &
     &                 Transl)
      IMPLICIT NONE
!*--DMAKE2563
!
!  Generates values for an M by N matrix A within the bandwidth
!  defined by KL and KU.
!  Stores the values in the array AA in the data structure required
!  by the routine, with unwanted elements set to rogue value.
!
!  TYPE is 'GE', 'GB', 'SY', 'SB', 'SP', 'TR', 'TB' OR 'TP'.
!
!  Auxiliary routine for test program for Level 2 Blas.
!
!  -- Written on 10-August-1987.
!     Richard Hanson, Sandia National Labs.
!     Jeremy Du Croz, NAG Central Office.
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO , ONE
      PARAMETER (ZERO=0.0D0,ONE=1.0D0)
      DOUBLE PRECISION ROGUE
      PARAMETER (ROGUE=-1.0D10)
!     .. Scalar Arguments ..
      DOUBLE PRECISION Transl
      INTEGER Kl , Ku , Lda , M , N , Nmax
      LOGICAL Reset
      CHARACTER*1 Diag , Uplo
      CHARACTER*2 Type
!     .. Array Arguments ..
      DOUBLE PRECISION A(Nmax,*) , Aa(*)
!     .. Local Scalars ..
      INTEGER i , i1 , i2 , i3 , ibeg , iend , ioff , j , kk
      LOGICAL gen , lower , sym , tri , unit , upper
!     .. External Functions ..
      DOUBLE PRECISION DBEG
      EXTERNAL DBEG
!     .. Intrinsic Functions ..
      INTRINSIC MAX , MIN
!     .. Executable Statements ..
      gen = Type(1:1)=='G'
      sym = Type(1:1)=='S'
      tri = Type(1:1)=='T'
      upper = (sym .OR. tri) .AND. Uplo=='U'
      lower = (sym .OR. tri) .AND. Uplo=='L'
      unit = tri .AND. Diag=='U'
!
!     Generate data in array A.
!
      DO j = 1 , N
         DO i = 1 , M
            IF ( gen .OR. (upper .AND. i<=j) .OR. (lower .AND. i>=j) )  &
     &           THEN
               IF ( (i<=j .AND. j-i<=Ku) .OR. (i>=j .AND. i-j<=Kl) )    &
     &              THEN
                  A(i,j) = DBEG(Reset) + Transl
               ELSE
                  A(i,j) = ZERO
               ENDIF
               IF ( i/=j ) THEN
                  IF ( sym ) THEN
                     A(j,i) = A(i,j)
                  ELSEIF ( tri ) THEN
                     A(j,i) = ZERO
                  ENDIF
               ENDIF
            ENDIF
         ENDDO
         IF ( tri ) A(j,j) = A(j,j) + ONE
         IF ( unit ) A(j,j) = ONE
      ENDDO
!
!     Store elements in array AS in data structure required by routine.
!
      IF ( Type=='GE' ) THEN
         DO j = 1 , N
            DO i = 1 , M
               Aa(i+(j-1)*Lda) = A(i,j)
            ENDDO
            DO i = M + 1 , Lda
               Aa(i+(j-1)*Lda) = ROGUE
            ENDDO
         ENDDO
      ELSEIF ( Type=='GB' ) THEN
         DO j = 1 , N
            DO i1 = 1 , Ku + 1 - j
               Aa(i1+(j-1)*Lda) = ROGUE
            ENDDO
            DO i2 = i1 , MIN(Kl+Ku+1,Ku+1+M-j)
               Aa(i2+(j-1)*Lda) = A(i2+j-Ku-1,j)
            ENDDO
            DO i3 = i2 , Lda
               Aa(i3+(j-1)*Lda) = ROGUE
            ENDDO
         ENDDO
      ELSEIF ( Type=='SY' .OR. Type=='TR' ) THEN
         DO j = 1 , N
            IF ( upper ) THEN
               ibeg = 1
               IF ( unit ) THEN
                  iend = j - 1
               ELSE
                  iend = j
               ENDIF
            ELSE
               IF ( unit ) THEN
                  ibeg = j + 1
               ELSE
                  ibeg = j
               ENDIF
               iend = N
            ENDIF
            DO i = 1 , ibeg - 1
               Aa(i+(j-1)*Lda) = ROGUE
            ENDDO
            DO i = ibeg , iend
               Aa(i+(j-1)*Lda) = A(i,j)
            ENDDO
            DO i = iend + 1 , Lda
               Aa(i+(j-1)*Lda) = ROGUE
            ENDDO
         ENDDO
      ELSEIF ( Type=='SB' .OR. Type=='TB' ) THEN
         DO j = 1 , N
            IF ( upper ) THEN
               kk = Kl + 1
               ibeg = MAX(1,Kl+2-j)
               IF ( unit ) THEN
                  iend = Kl
               ELSE
                  iend = Kl + 1
               ENDIF
            ELSE
               kk = 1
               IF ( unit ) THEN
                  ibeg = 2
               ELSE
                  ibeg = 1
               ENDIF
               iend = MIN(Kl+1,1+M-j)
            ENDIF
            DO i = 1 , ibeg - 1
               Aa(i+(j-1)*Lda) = ROGUE
            ENDDO
            DO i = ibeg , iend
               Aa(i+(j-1)*Lda) = A(i+j-kk,j)
            ENDDO
            DO i = iend + 1 , Lda
               Aa(i+(j-1)*Lda) = ROGUE
            ENDDO
         ENDDO
      ELSEIF ( Type=='SP' .OR. Type=='TP' ) THEN
         ioff = 0
         DO j = 1 , N
            IF ( upper ) THEN
               ibeg = 1
               iend = j
            ELSE
               ibeg = j
               iend = N
            ENDIF
            DO i = ibeg , iend
               ioff = ioff + 1
               Aa(ioff) = A(i,j)
               IF ( i==j ) THEN
                  IF ( unit ) Aa(ioff) = ROGUE
               ENDIF
            ENDDO
         ENDDO
      ENDIF
!
!     End of DMAKE.
!
      END SUBROUTINE DMAKE
!*==dmvch.f90  processed by SPAG 7.51RB at 17:52 on  4 Mar 2022
      SUBROUTINE DMVCH(Trans,M,N,Alpha,A,Nmax,X,Incx,Beta,Y,Incy,Yt,G,  &
     &                 Yy,Eps,Err,Fatal,Nout,Mv)
      IMPLICIT NONE
!*--DMVCH2738
!
!  Checks the results of the computational tests.
!
!  Auxiliary routine for test program for Level 2 Blas.
!
!  -- Written on 10-August-1987.
!     Richard Hanson, Sandia National Labs.
!     Jeremy Du Croz, NAG Central Office.
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO , ONE
      PARAMETER (ZERO=0.0D0,ONE=1.0D0)
!     .. Scalar Arguments ..
      DOUBLE PRECISION Alpha , Beta , Eps , Err
      INTEGER Incx , Incy , M , N , Nmax , Nout
      LOGICAL Fatal , Mv
      CHARACTER*1 Trans
!     .. Array Arguments ..
      DOUBLE PRECISION A(Nmax,*) , G(*) , X(*) , Y(*) , Yt(*) , Yy(*)
!     .. Local Scalars ..
      DOUBLE PRECISION erri
      INTEGER i , incxl , incyl , iy , j , jx , kx , ky , ml , nl
      LOGICAL tran
!     .. Intrinsic Functions ..
      INTRINSIC ABS , MAX , SQRT
!     .. Executable Statements ..
      tran = Trans=='T' .OR. Trans=='C'
      IF ( tran ) THEN
         ml = N
         nl = M
      ELSE
         ml = M
         nl = N
      ENDIF
      IF ( Incx<0 ) THEN
         kx = nl
         incxl = -1
      ELSE
         kx = 1
         incxl = 1
      ENDIF
      IF ( Incy<0 ) THEN
         ky = ml
         incyl = -1
      ELSE
         ky = 1
         incyl = 1
      ENDIF
!
!     Compute expected result in YT using data in A, X and Y.
!     Compute gauges in G.
!
      iy = ky
      DO i = 1 , ml
         Yt(iy) = ZERO
         G(iy) = ZERO
         jx = kx
         IF ( tran ) THEN
            DO j = 1 , nl
               Yt(iy) = Yt(iy) + A(j,i)*X(jx)
               G(iy) = G(iy) + ABS(A(j,i)*X(jx))
               jx = jx + incxl
            ENDDO
         ELSE
            DO j = 1 , nl
               Yt(iy) = Yt(iy) + A(i,j)*X(jx)
               G(iy) = G(iy) + ABS(A(i,j)*X(jx))
               jx = jx + incxl
            ENDDO
         ENDIF
         Yt(iy) = Alpha*Yt(iy) + Beta*Y(iy)
         G(iy) = ABS(Alpha)*G(iy) + ABS(Beta*Y(iy))
         iy = iy + incyl
      ENDDO
!
!     Compute the error ratio for this result.
!
      Err = ZERO
      DO i = 1 , ml
         erri = ABS(Yt(i)-Yy(1+(i-1)*ABS(Incy)))/Eps
         IF ( G(i)/=ZERO ) erri = erri/G(i)
         Err = MAX(Err,erri)
         IF ( Err*SQRT(Eps)>=ONE ) GOTO 100
      ENDDO
!     If the loop completes, all results are at least half accurate.
      GOTO 200
!
!     Report fatal error.
!
 100  Fatal = .TRUE.
      WRITE (Nout,FMT=99001)
      DO i = 1 , ml
         IF ( Mv ) THEN
            WRITE (Nout,FMT=99002) i , Yt(i) , Yy(1+(i-1)*ABS(Incy))
         ELSE
            WRITE (Nout,FMT=99002) i , Yy(1+(i-1)*ABS(Incy)) , Yt(i)
         ENDIF
      ENDDO
!
 200  RETURN
!
99001 FORMAT (' ******* FATAL ERROR - COMPUTED RESULT IS LESS THAN HAL',&
     &        'F ACCURATE *******',                                     &
     &        /'           EXPECTED RESULT   COMPU','TED RESULT')
99002 FORMAT (1X,I7,2G18.6)
!
!     End of DMVCH.
!
      END SUBROUTINE DMVCH
!*==lde.f90  processed by SPAG 7.51RB at 17:52 on  4 Mar 2022
      LOGICAL FUNCTION LDE(Ri,Rj,Lr)
      IMPLICIT NONE
!*--LDE2851
!
!  Tests if two arrays are identical.
!
!  Auxiliary routine for test program for Level 2 Blas.
!
!  -- Written on 10-August-1987.
!     Richard Hanson, Sandia National Labs.
!     Jeremy Du Croz, NAG Central Office.
!
!     .. Scalar Arguments ..
      INTEGER Lr
!     .. Array Arguments ..
      DOUBLE PRECISION Ri(*) , Rj(*)
!     .. Local Scalars ..
      INTEGER i
!     .. Executable Statements ..
      DO i = 1 , Lr
         IF ( Ri(i)/=Rj(i) ) GOTO 100
      ENDDO
      LDE = .TRUE.
      GOTO 99999
 100  LDE = .FALSE.
!
!     End of LDE.
!
99999 END FUNCTION LDE
!*==lderes.f90  processed by SPAG 7.51RB at 17:52 on  4 Mar 2022
      LOGICAL FUNCTION LDERES(Type,Uplo,M,N,Aa,As,Lda)
      IMPLICIT NONE
!*--LDERES2881
!
!  Tests if selected elements in two arrays are equal.
!
!  TYPE is 'GE', 'SY' or 'SP'.
!
!  Auxiliary routine for test program for Level 2 Blas.
!
!  -- Written on 10-August-1987.
!     Richard Hanson, Sandia National Labs.
!     Jeremy Du Croz, NAG Central Office.
!
!     .. Scalar Arguments ..
      INTEGER Lda , M , N
      CHARACTER*1 Uplo
      CHARACTER*2 Type
!     .. Array Arguments ..
      DOUBLE PRECISION Aa(Lda,*) , As(Lda,*)
!     .. Local Scalars ..
      INTEGER i , ibeg , iend , j
      LOGICAL upper
!     .. Executable Statements ..
      upper = Uplo=='U'
      IF ( Type=='GE' ) THEN
         DO j = 1 , N
            DO i = M + 1 , Lda
               IF ( Aa(i,j)/=As(i,j) ) GOTO 100
            ENDDO
         ENDDO
      ELSEIF ( Type=='SY' ) THEN
         DO j = 1 , N
            IF ( upper ) THEN
               ibeg = 1
               iend = j
            ELSE
               ibeg = j
               iend = N
            ENDIF
            DO i = 1 , ibeg - 1
               IF ( Aa(i,j)/=As(i,j) ) GOTO 100
            ENDDO
            DO i = iend + 1 , Lda
               IF ( Aa(i,j)/=As(i,j) ) GOTO 100
            ENDDO
         ENDDO
      ENDIF
!
      LDERES = .TRUE.
      GOTO 99999
 100  LDERES = .FALSE.
!
!     End of LDERES.
!
99999 END FUNCTION LDERES
!*==dbeg.f90  processed by SPAG 7.51RB at 17:52 on  4 Mar 2022
      DOUBLE PRECISION FUNCTION DBEG(Reset)
      IMPLICIT NONE
!*--DBEG2938
!
!  Generates random numbers uniformly distributed between -0.5 and 0.5.
!
!  Auxiliary routine for test program for Level 2 Blas.
!
!  -- Written on 10-August-1987.
!     Richard Hanson, Sandia National Labs.
!     Jeremy Du Croz, NAG Central Office.
!
!     .. Scalar Arguments ..
      LOGICAL Reset
!     .. Local Scalars ..
      INTEGER i , ic , mi
!     .. Save statement ..
      SAVE i , ic , mi
!     .. Intrinsic Functions ..
      INTRINSIC DBLE
!     .. Executable Statements ..
      IF ( Reset ) THEN
!        Initialize local variables.
         mi = 891
         i = 7
         ic = 0
         Reset = .FALSE.
      ENDIF
!
!     The sequence of values of I is bounded between 1 and 999.
!     If initial I = 1,2,3,6,7 or 9, the period will be 50.
!     If initial I = 4 or 8, the period will be 25.
!     If initial I = 5, the period will be 10.
!     IC is used to break up the period by skipping 1 value of I in 6.
!
      ic = ic + 1
      DO
         i = i*mi
         i = i - 1000*(i/1000)
         IF ( ic>=5 ) THEN
            ic = 0
            CYCLE
         ENDIF
         DBEG = DBLE(i-500)/1001.0D0
         EXIT
      ENDDO
!
!     End of DBEG.
!
      END FUNCTION DBEG
!*==ddiff.f90  processed by SPAG 7.51RB at 17:52 on  4 Mar 2022
      DOUBLE PRECISION FUNCTION DDIFF(X,Y)
      IMPLICIT NONE
!*--DDIFF2989
!
!  Auxiliary routine for test program for Level 2 Blas.
!
!  -- Written on 10-August-1987.
!     Richard Hanson, Sandia National Labs.
!
!     .. Scalar Arguments ..
      DOUBLE PRECISION X , Y
!     .. Executable Statements ..
      DDIFF = X - Y
!
!     End of DDIFF.
!
      END FUNCTION DDIFF
!*==chkxer.f90  processed by SPAG 7.51RB at 17:52 on  4 Mar 2022
      SUBROUTINE CHKXER(Srnamt,Infot,Nout,Lerr,Ok)
      IMPLICIT NONE
!*--CHKXER3007
!
!  Tests whether XERBLA has detected an error when it should.
!
!  Auxiliary routine for test program for Level 2 Blas.
!
!  -- Written on 10-August-1987.
!     Richard Hanson, Sandia National Labs.
!     Jeremy Du Croz, NAG Central Office.
!
!     .. Scalar Arguments ..
      INTEGER Infot , Nout
      LOGICAL Lerr , Ok
      CHARACTER*6 Srnamt
!     .. Executable Statements ..
      IF ( .NOT.Lerr ) THEN
         WRITE (Nout,FMT=99001) Infot , Srnamt
         Ok = .FALSE.
      ENDIF
      Lerr = .FALSE.
      RETURN
!
99001 FORMAT (' ***** ILLEGAL VALUE OF PARAMETER NUMBER ',I2,' NOT D',  &
     &        'ETECTED BY ',A6,' *****')
!
!     End of CHKXER.
!
      END SUBROUTINE CHKXER
!*==xerbla.f90  processed by SPAG 7.51RB at 17:52 on  4 Mar 2022
      SUBROUTINE XERBLA(Srname,Info)
      IMPLICIT NONE
!*--XERBLA3038
!
!  This is a special version of XERBLA to be used only as part of
!  the test program for testing error exits from the Level 2 BLAS
!  routines.
!
!  XERBLA  is an error handler for the Level 2 BLAS routines.
!
!  It is called by the Level 2 BLAS routines if an input parameter is
!  invalid.
!
!  Auxiliary routine for test program for Level 2 Blas.
!
!  -- Written on 10-August-1987.
!     Richard Hanson, Sandia National Labs.
!     Jeremy Du Croz, NAG Central Office.
!
!     .. Scalar Arguments ..
      INTEGER Info
      CHARACTER*6 Srname
!     .. Scalars in Common ..
      INTEGER INFot , NOUt
      LOGICAL LERr , OK
      CHARACTER*6 SRNamt
!     .. Common blocks ..
      COMMON /INFOC / INFot , NOUt , OK , LERr
      COMMON /SRNAMC/ SRNamt
!     .. Executable Statements ..
      LERr = .TRUE.
      IF ( Info/=INFot ) THEN
         IF ( INFot/=0 ) THEN
            WRITE (NOUt,FMT=99001) Info , INFot
         ELSE
            WRITE (NOUt,FMT=99003) Info
         ENDIF
         OK = .FALSE.
      ENDIF
      IF ( Srname/=SRNamt ) THEN
         WRITE (NOUt,FMT=99002) Srname , SRNamt
         OK = .FALSE.
      ENDIF
      RETURN
!
99001 FORMAT (' ******* XERBLA WAS CALLED WITH INFO = ',I6,' INSTEAD',  &
     &        ' OF ',I2,' *******')
99002 FORMAT (' ******* XERBLA WAS CALLED WITH SRNAME = ',A6,' INSTE',  &
     &        'AD OF ',A6,' *******')
99003 FORMAT (' ******* XERBLA WAS CALLED WITH INFO = ',I6,' *******')
!
!     End of XERBLA
!
      END SUBROUTINE XERBLA
