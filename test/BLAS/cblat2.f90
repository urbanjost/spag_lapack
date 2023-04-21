!*==cblat2.f90  processed by SPAG 7.51RB at 17:52 on  4 Mar 2022
!> \brief \b CBLAT2
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       PROGRAM CBLAT2
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> Test program for the COMPLEX          Level 2 Blas.
!>
!> The program must be driven by a short data file. The first 18 records
!> of the file are read using list-directed input, the last 17 records
!> are read using the format ( A6, L2 ). An annotated example of a data
!> file can be obtained by deleting the first 3 characters from the
!> following 35 lines:
!> 'cblat2.out'      NAME OF SUMMARY OUTPUT FILE
!> 6                 UNIT NUMBER OF SUMMARY FILE
!> 'CBLA2T.SNAP'     NAME OF SNAPSHOT OUTPUT FILE
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
!> (0.0,0.0) (1.0,0.0) (0.7,-0.9)       VALUES OF ALPHA
!> 3                 NUMBER OF VALUES OF BETA
!> (0.0,0.0) (1.0,0.0) (1.3,-1.1)       VALUES OF BETA
!> CGEMV  T PUT F FOR NO TEST. SAME COLUMNS.
!> CGBMV  T PUT F FOR NO TEST. SAME COLUMNS.
!> CHEMV  T PUT F FOR NO TEST. SAME COLUMNS.
!> CHBMV  T PUT F FOR NO TEST. SAME COLUMNS.
!> CHPMV  T PUT F FOR NO TEST. SAME COLUMNS.
!> CTRMV  T PUT F FOR NO TEST. SAME COLUMNS.
!> CTBMV  T PUT F FOR NO TEST. SAME COLUMNS.
!> CTPMV  T PUT F FOR NO TEST. SAME COLUMNS.
!> CTRSV  T PUT F FOR NO TEST. SAME COLUMNS.
!> CTBSV  T PUT F FOR NO TEST. SAME COLUMNS.
!> CTPSV  T PUT F FOR NO TEST. SAME COLUMNS.
!> CGERC  T PUT F FOR NO TEST. SAME COLUMNS.
!> CGERU  T PUT F FOR NO TEST. SAME COLUMNS.
!> CHER   T PUT F FOR NO TEST. SAME COLUMNS.
!> CHPR   T PUT F FOR NO TEST. SAME COLUMNS.
!> CHER2  T PUT F FOR NO TEST. SAME COLUMNS.
!> CHPR2  T PUT F FOR NO TEST. SAME COLUMNS.
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
!> \ingroup complex_blas_testing
!
!  =====================================================================
      PROGRAM CBLAT2
      IMPLICIT NONE
!*--CBLAT2107
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
      PARAMETER (NSUBS=17)
      COMPLEX ZERO , ONE
      PARAMETER (ZERO=(0.0,0.0),ONE=(1.0,0.0))
      REAL RZERO
      PARAMETER (RZERO=0.0)
      INTEGER NMAX , INCMAX
      PARAMETER (NMAX=65,INCMAX=2)
      INTEGER NINMAX , NIDMAX , NKBMAX , NALMAX , NBEMAX
      PARAMETER (NINMAX=7,NIDMAX=9,NKBMAX=7,NALMAX=7,NBEMAX=7)
!     .. Local Scalars ..
      REAL eps , err , thresh
      INTEGER i , isnum , j , n , nalf , nbet , nidim , ninc , nkb ,    &
     &        nout , ntra
      LOGICAL fatal , ltestt , rewi , same , sfatal , trace , tsterr
      CHARACTER*1 trans
      CHARACTER*6 snamet
      CHARACTER*32 snaps , summry
!     .. Local Arrays ..
      COMPLEX a(NMAX,NMAX) , aa(NMAX*NMAX) , alf(NALMAX) , as(NMAX*NMAX)&
     &        , bet(NBEMAX) , x(NMAX) , xs(NMAX*INCMAX) ,               &
     &        xx(NMAX*INCMAX) , y(NMAX) , ys(NMAX*INCMAX) , yt(NMAX) ,  &
     &        yy(NMAX*INCMAX) , z(2*NMAX)
      REAL g(NMAX)
      INTEGER idim(NIDMAX) , inc(NINMAX) , kb(NKBMAX)
      LOGICAL ltest(NSUBS)
      CHARACTER*6 snames(NSUBS)
!     .. External Functions ..
      REAL SDIFF
      LOGICAL LCE
      EXTERNAL SDIFF , LCE
!     .. External Subroutines ..
      EXTERNAL CCHK1 , CCHK2 , CCHK3 , CCHK4 , CCHK5 , CCHK6 , CCHKE ,  &
     &         CMVCH
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
      DATA snames/'CGEMV ' , 'CGBMV ' , 'CHEMV ' , 'CHBMV ' , 'CHPMV ' ,&
     &     'CTRMV ' , 'CTBMV ' , 'CTPMV ' , 'CTRSV ' , 'CTBSV ' ,       &
     &     'CTPSV ' , 'CGERC ' , 'CGERU ' , 'CHER  ' , 'CHPR  ' ,       &
     &     'CHER2 ' , 'CHPR2 '/
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
      eps = EPSILON(RZERO)
      WRITE (nout,FMT=99002) eps
!
!     Check the reliability of CMVCH using exact data.
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
!     YY holds the exact result. On exit from CMVCH YT holds
!     the result computed by CMVCH.
      trans = 'N'
      CALL CMVCH(trans,n,n,ONE,a,NMAX,x,1,ZERO,y,1,yt,g,yy,eps,err,     &
     &           fatal,nout,.TRUE.)
      same = LCE(yy,yt,n)
      IF ( .NOT.same .OR. err/=RZERO ) THEN
         WRITE (nout,FMT=99015) trans , same , err
         STOP
      ENDIF
      trans = 'T'
      CALL CMVCH(trans,n,n,ONE,a,NMAX,x,-1,ZERO,y,-1,yt,g,yy,eps,err,   &
     &           fatal,nout,.TRUE.)
      same = LCE(yy,yt,n)
      IF ( .NOT.same .OR. err/=RZERO ) THEN
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
               CALL CCHKE(isnum,snames(isnum),nout)
               WRITE (nout,FMT=*)
            ENDIF
!           Test computations.
            INFot = 0
            OK = .TRUE.
            fatal = .FALSE.
            IF ( isnum==3 .OR. isnum==4 .OR. isnum==5 ) THEN
!           Test CHEMV, 03, CHBMV, 04, and CHPMV, 05.
               CALL CCHK2(snames(isnum),eps,thresh,nout,ntra,trace,rewi,&
     &                    fatal,nidim,idim,nkb,kb,nalf,alf,nbet,bet,    &
     &                    ninc,inc,NMAX,INCMAX,a,aa,as,x,xx,xs,y,yy,ys, &
     &                    yt,g)
            ELSEIF ( isnum==6 .OR. isnum==7 .OR. isnum==8 .OR.          &
     &               isnum==9 .OR. isnum==10 .OR. isnum==11 ) THEN
!           Test CTRMV, 06, CTBMV, 07, CTPMV, 08,
!           CTRSV, 09, CTBSV, 10, and CTPSV, 11.
               CALL CCHK3(snames(isnum),eps,thresh,nout,ntra,trace,rewi,&
     &                    fatal,nidim,idim,nkb,kb,ninc,inc,NMAX,INCMAX, &
     &                    a,aa,as,y,yy,ys,yt,g,z)
            ELSEIF ( isnum==12 .OR. isnum==13 ) THEN
!           Test CGERC, 12, CGERU, 13.
               CALL CCHK4(snames(isnum),eps,thresh,nout,ntra,trace,rewi,&
     &                    fatal,nidim,idim,nalf,alf,ninc,inc,NMAX,      &
     &                    INCMAX,a,aa,as,x,xx,xs,y,yy,ys,yt,g,z)
            ELSEIF ( isnum==14 .OR. isnum==15 ) THEN
!           Test CHER, 14, and CHPR, 15.
               CALL CCHK5(snames(isnum),eps,thresh,nout,ntra,trace,rewi,&
     &                    fatal,nidim,idim,nalf,alf,ninc,inc,NMAX,      &
     &                    INCMAX,a,aa,as,x,xx,xs,y,yy,ys,yt,g,z)
            ELSEIF ( isnum==16 .OR. isnum==17 ) THEN
!           Test CHER2, 16, and CHPR2, 17.
               CALL CCHK6(snames(isnum),eps,thresh,nout,ntra,trace,rewi,&
     &                    fatal,nidim,idim,nalf,alf,ninc,inc,NMAX,      &
     &                    INCMAX,a,aa,as,x,xx,xs,y,yy,ys,yt,g,z)
            ELSE
!           Test CGEMV, 01, and CGBMV, 02.
               CALL CCHK1(snames(isnum),eps,thresh,nout,ntra,trace,rewi,&
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
99002 FORMAT (' RELATIVE MACHINE PRECISION IS TAKEN TO BE',1P,E9.1)
99003 FORMAT (' NUMBER OF VALUES OF ',A,' IS LESS THAN 1 OR GREATER ',  &
     &        'THAN ',I2)
99004 FORMAT (' VALUE OF N IS LESS THAN 0 OR GREATER THAN ',I2)
99005 FORMAT (' VALUE OF K IS LESS THAN 0')
99006 FORMAT (' ABSOLUTE VALUE OF INCX OR INCY IS 0 OR GREATER THAN ',  &
     &        I2)
99007 FORMAT (' TESTS OF THE COMPLEX          LEVEL 2 BLAS',//' THE F', &
     &        'OLLOWING PARAMETER VALUES WILL BE USED:')
99008 FORMAT ('   FOR N              ',9I6)
99009 FORMAT ('   FOR K              ',7I6)
99010 FORMAT ('   FOR INCX AND INCY  ',7I6)
99011 FORMAT ('   FOR ALPHA          ',7('(',F4.1,',',F4.1,')  ',:))
99012 FORMAT ('   FOR BETA           ',7('(',F4.1,',',F4.1,')  ',:))
99013 FORMAT (' AMEND DATA FILE OR INCREASE ARRAY SIZES IN PROGRAM',    &
     &        /' ******* TESTS ABANDONED *******')
99014 FORMAT (' SUBPROGRAM NAME ',A6,' NOT RECOGNIZED',/' ******* T',   &
     &        'ESTS ABANDONED *******')
99015 FORMAT (' ERROR IN CMVCH -  IN-LINE DOT PRODUCTS ARE BEING EVALU',&
     &        'ATED WRONGLY.',/' CMVCH WAS CALLED WITH TRANS = ',A1,    &
     &        ' AND RETURNED SAME = ',L1,' AND ERR = ',F12.3,'.',/      &
     &   ' THIS MAY BE DUE TO FAULTS IN THE ARITHMETIC OR THE COMPILER.'&
     &   ,/' ******* TESTS ABANDONED *******')
99016 FORMAT (A6,L2)
99017 FORMAT (1X,A6,' WAS NOT TESTED')
99018 FORMAT (/' END OF TESTS')
99019 FORMAT (/' ******* FATAL ERROR - TESTS ABANDONED *******')
99020 FORMAT (' ERROR-EXITS WILL NOT BE TESTED')
!
!     End of CBLAT2.
!
      END PROGRAM CBLAT2
!*==cchk1.f90  processed by SPAG 7.51RB at 17:52 on  4 Mar 2022
      SUBROUTINE CCHK1(Sname,Eps,Thresh,Nout,Ntra,Trace,Rewi,Fatal,     &
     &                 Nidim,Idim,Nkb,Kb,Nalf,Alf,Nbet,Bet,Ninc,Inc,    &
     &                 Nmax,Incmax,A,Aa,As,X,Xx,Xs,Y,Yy,Ys,Yt,G)
      IMPLICIT NONE
!*--CCHK1426
!
!  Tests CGEMV and CGBMV.
!
!  Auxiliary routine for test program for Level 2 Blas.
!
!  -- Written on 10-August-1987.
!     Richard Hanson, Sandia National Labs.
!     Jeremy Du Croz, NAG Central Office.
!
!     .. Parameters ..
      COMPLEX ZERO , HALF
      PARAMETER (ZERO=(0.0,0.0),HALF=(0.5,0.0))
      REAL RZERO
      PARAMETER (RZERO=0.0)
!     .. Scalar Arguments ..
      REAL Eps , Thresh
      INTEGER Incmax , Nalf , Nbet , Nidim , Ninc , Nkb , Nmax , Nout , &
     &        Ntra
      LOGICAL Fatal , Rewi , Trace
      CHARACTER*6 Sname
!     .. Array Arguments ..
      COMPLEX A(Nmax,Nmax) , Aa(Nmax*Nmax) , Alf(Nalf) , As(Nmax*Nmax) ,&
     &        Bet(Nbet) , X(Nmax) , Xs(Nmax*Incmax) , Xx(Nmax*Incmax) , &
     &        Y(Nmax) , Ys(Nmax*Incmax) , Yt(Nmax) , Yy(Nmax*Incmax)
      REAL G(Nmax)
      INTEGER Idim(Nidim) , Inc(Ninc) , Kb(Nkb)
!     .. Local Scalars ..
      COMPLEX alpha , als , beta , bls , transl
      REAL err , errmax
      INTEGER i , ia , ib , ic , iku , im , in , incx , incxs , incy ,  &
     &        incys , ix , iy , kl , kls , ku , kus , laa , lda , ldas ,&
     &        lx , ly , m , ml , ms , n , nargs , nc , nd , nk , nl , ns
      LOGICAL banded , full , null , reset , same , tran
      CHARACTER*1 trans , transs
      CHARACTER*3 ich
!     .. Local Arrays ..
      LOGICAL isame(13)
!     .. External Functions ..
      LOGICAL LCE , LCERES
      EXTERNAL LCE , LCERES
!     .. External Subroutines ..
      EXTERNAL CGBMV , CGEMV , CMAKE , CMVCH
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
      errmax = RZERO
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
                  CALL CMAKE(Sname(2:3),' ',' ',m,n,A,Nmax,Aa,lda,kl,ku,&
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
                        CALL CMAKE('GE',' ',' ',1,nl,X,1,Xx,ABS(incx),0,&
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
                                 CALL CMAKE('GE',' ',' ',1,ml,Y,1,Yy,   &
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
                                    CALL CGEMV(trans,m,n,alpha,Aa,lda,  &
     &                                 Xx,incx,beta,Yy,incy)
                                 ELSEIF ( banded ) THEN
                                    IF ( Trace ) WRITE (Ntra,FMT=99005) &
     &                                 nc , Sname , trans , m , n , kl ,&
     &                                 ku , alpha , lda , incx , beta , &
     &                                 incy
                                    IF ( Rewi ) REWIND Ntra
                                    CALL CGBMV(trans,m,n,kl,ku,alpha,Aa,&
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
                                    isame(5) = LCE(As,Aa,laa)
                                    isame(6) = ldas==lda
                                    isame(7) = LCE(Xs,Xx,lx)
                                    isame(8) = incxs==incx
                                    isame(9) = bls==beta
                                    IF ( null ) THEN
                                       isame(10) = LCE(Ys,Yy,ly)
                                    ELSE
                                       isame(10)                        &
     &                                    = LCERES('GE',' ',1,ml,Ys,Yy, &
     &                                    ABS(incy))
                                    ENDIF
                                    isame(11) = incys==incy
                                 ELSEIF ( banded ) THEN
                                    isame(4) = kls==kl
                                    isame(5) = kus==ku
                                    isame(6) = als==alpha
                                    isame(7) = LCE(As,Aa,laa)
                                    isame(8) = ldas==lda
                                    isame(9) = LCE(Xs,Xx,lx)
                                    isame(10) = incxs==incx
                                    isame(11) = bls==beta
                                    IF ( null ) THEN
                                       isame(12) = LCE(Ys,Yy,ly)
                                    ELSE
                                       isame(12)                        &
     &                                    = LCERES('GE',' ',1,ml,Ys,Yy, &
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
                                 CALL CMVCH(trans,m,n,alpha,A,Nmax,X,   &
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
99005 FORMAT (1X,I6,': ',A6,'(''',A1,''',',4(I3,','),'(',F4.1,',',F4.1, &
     &        '), A,',I3,', X,',I2,',(',F4.1,',',F4.1,'), Y,',I2,') .')
99006 FORMAT (1X,I6,': ',A6,'(''',A1,''',',2(I3,','),'(',F4.1,',',F4.1, &
     &        '), A,',I3,', X,',I2,',(',F4.1,',',F4.1,'), Y,',I2,       &
     &        ')         .')
99007 FORMAT (' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *',&
     &        '******')
!
!     End of CCHK1.
!
      END SUBROUTINE CCHK1
!*==cchk2.f90  processed by SPAG 7.51RB at 17:52 on  4 Mar 2022
      SUBROUTINE CCHK2(Sname,Eps,Thresh,Nout,Ntra,Trace,Rewi,Fatal,     &
     &                 Nidim,Idim,Nkb,Kb,Nalf,Alf,Nbet,Bet,Ninc,Inc,    &
     &                 Nmax,Incmax,A,Aa,As,X,Xx,Xs,Y,Yy,Ys,Yt,G)
      IMPLICIT NONE
!*--CCHK2753
!
!  Tests CHEMV, CHBMV and CHPMV.
!
!  Auxiliary routine for test program for Level 2 Blas.
!
!  -- Written on 10-August-1987.
!     Richard Hanson, Sandia National Labs.
!     Jeremy Du Croz, NAG Central Office.
!
!     .. Parameters ..
      COMPLEX ZERO , HALF
      PARAMETER (ZERO=(0.0,0.0),HALF=(0.5,0.0))
      REAL RZERO
      PARAMETER (RZERO=0.0)
!     .. Scalar Arguments ..
      REAL Eps , Thresh
      INTEGER Incmax , Nalf , Nbet , Nidim , Ninc , Nkb , Nmax , Nout , &
     &        Ntra
      LOGICAL Fatal , Rewi , Trace
      CHARACTER*6 Sname
!     .. Array Arguments ..
      COMPLEX A(Nmax,Nmax) , Aa(Nmax*Nmax) , Alf(Nalf) , As(Nmax*Nmax) ,&
     &        Bet(Nbet) , X(Nmax) , Xs(Nmax*Incmax) , Xx(Nmax*Incmax) , &
     &        Y(Nmax) , Ys(Nmax*Incmax) , Yt(Nmax) , Yy(Nmax*Incmax)
      REAL G(Nmax)
      INTEGER Idim(Nidim) , Inc(Ninc) , Kb(Nkb)
!     .. Local Scalars ..
      COMPLEX alpha , als , beta , bls , transl
      REAL err , errmax
      INTEGER i , ia , ib , ic , ik , in , incx , incxs , incy , incys ,&
     &        ix , iy , k , ks , laa , lda , ldas , lx , ly , n ,       &
     &        nargs , nc , nk , ns
      LOGICAL banded , full , null , packed , reset , same
      CHARACTER*1 uplo , uplos
      CHARACTER*2 ich
!     .. Local Arrays ..
      LOGICAL isame(13)
!     .. External Functions ..
      LOGICAL LCE , LCERES
      EXTERNAL LCE , LCERES
!     .. External Subroutines ..
      EXTERNAL CHBMV , CHEMV , CHPMV , CMAKE , CMVCH
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
      full = Sname(3:3)=='E'
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
      errmax = RZERO
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
                  CALL CMAKE(Sname(2:3),uplo,' ',n,n,A,Nmax,Aa,lda,k,k, &
     &                       reset,transl)
!
                  DO ix = 1 , Ninc
                     incx = Inc(ix)
                     lx = ABS(incx)*n
!
!                 Generate the vector X.
!
                     transl = HALF
                     CALL CMAKE('GE',' ',' ',1,n,X,1,Xx,ABS(incx),0,n-1,&
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
                              CALL CMAKE('GE',' ',' ',1,n,Y,1,Yy,       &
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
                                 CALL CHEMV(uplo,n,alpha,Aa,lda,Xx,incx,&
     &                              beta,Yy,incy)
                              ELSEIF ( banded ) THEN
                                 IF ( Trace ) WRITE (Ntra,FMT=99006)    &
     &                                nc , Sname , uplo , n , k ,       &
     &                                alpha , lda , incx , beta , incy
                                 IF ( Rewi ) REWIND Ntra
                                 CALL CHBMV(uplo,n,k,alpha,Aa,lda,Xx,   &
     &                              incx,beta,Yy,incy)
                              ELSEIF ( packed ) THEN
                                 IF ( Trace ) WRITE (Ntra,FMT=99005)    &
     &                                nc , Sname , uplo , n , alpha ,   &
     &                                incx , beta , incy
                                 IF ( Rewi ) REWIND Ntra
                                 CALL CHPMV(uplo,n,alpha,Aa,Xx,incx,    &
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
                                 isame(4) = LCE(As,Aa,laa)
                                 isame(5) = ldas==lda
                                 isame(6) = LCE(Xs,Xx,lx)
                                 isame(7) = incxs==incx
                                 isame(8) = bls==beta
                                 IF ( null ) THEN
                                    isame(9) = LCE(Ys,Yy,ly)
                                 ELSE
                                    isame(9)                            &
     &                                 = LCERES('GE',' ',1,n,Ys,Yy,     &
     &                                 ABS(incy))
                                 ENDIF
                                 isame(10) = incys==incy
                              ELSEIF ( banded ) THEN
                                 isame(3) = ks==k
                                 isame(4) = als==alpha
                                 isame(5) = LCE(As,Aa,laa)
                                 isame(6) = ldas==lda
                                 isame(7) = LCE(Xs,Xx,lx)
                                 isame(8) = incxs==incx
                                 isame(9) = bls==beta
                                 IF ( null ) THEN
                                    isame(10) = LCE(Ys,Yy,ly)
                                 ELSE
                                    isame(10)                           &
     &                                 = LCERES('GE',' ',1,n,Ys,Yy,     &
     &                                 ABS(incy))
                                 ENDIF
                                 isame(11) = incys==incy
                              ELSEIF ( packed ) THEN
                                 isame(3) = als==alpha
                                 isame(4) = LCE(As,Aa,laa)
                                 isame(5) = LCE(Xs,Xx,lx)
                                 isame(6) = incxs==incx
                                 isame(7) = bls==beta
                                 IF ( null ) THEN
                                    isame(8) = LCE(Ys,Yy,ly)
                                 ELSE
                                    isame(8)                            &
     &                                 = LCERES('GE',' ',1,n,Ys,Yy,     &
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
                              CALL CMVCH('N',n,n,alpha,A,Nmax,X,incx,   &
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
99005 FORMAT (1X,I6,': ',A6,'(''',A1,''',',I3,',(',F4.1,',',F4.1,       &
     &        '), AP, X,',I2,',(',F4.1,',',F4.1,'), Y,',I2,             &
     &        ')                .')
99006 FORMAT (1X,I6,': ',A6,'(''',A1,''',',2(I3,','),'(',F4.1,',',F4.1, &
     &        '), A,',I3,', X,',I2,',(',F4.1,',',F4.1,'), Y,',I2,       &
     &        ')         .')
99007 FORMAT (1X,I6,': ',A6,'(''',A1,''',',I3,',(',F4.1,',',F4.1,       &
     &        '), A,',I3,', X,',I2,',(',F4.1,',',F4.1,'), ','Y,',I2,    &
     &        ')             .')
99008 FORMAT (' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *',&
     &        '******')
!
!     End of CCHK2.
!
      END SUBROUTINE CCHK2
!*==cchk3.f90  processed by SPAG 7.51RB at 17:52 on  4 Mar 2022
      SUBROUTINE CCHK3(Sname,Eps,Thresh,Nout,Ntra,Trace,Rewi,Fatal,     &
     &                 Nidim,Idim,Nkb,Kb,Ninc,Inc,Nmax,Incmax,A,Aa,As,X,&
     &                 Xx,Xs,Xt,G,Z)
      IMPLICIT NONE
!*--CCHK31091
!
!  Tests CTRMV, CTBMV, CTPMV, CTRSV, CTBSV and CTPSV.
!
!  Auxiliary routine for test program for Level 2 Blas.
!
!  -- Written on 10-August-1987.
!     Richard Hanson, Sandia National Labs.
!     Jeremy Du Croz, NAG Central Office.
!
!     .. Parameters ..
      COMPLEX ZERO , HALF , ONE
      PARAMETER (ZERO=(0.0,0.0),HALF=(0.5,0.0),ONE=(1.0,0.0))
      REAL RZERO
      PARAMETER (RZERO=0.0)
!     .. Scalar Arguments ..
      REAL Eps , Thresh
      INTEGER Incmax , Nidim , Ninc , Nkb , Nmax , Nout , Ntra
      LOGICAL Fatal , Rewi , Trace
      CHARACTER*6 Sname
!     .. Array Arguments ..
      COMPLEX A(Nmax,Nmax) , Aa(Nmax*Nmax) , As(Nmax*Nmax) , X(Nmax) ,  &
     &        Xs(Nmax*Incmax) , Xt(Nmax) , Xx(Nmax*Incmax) , Z(Nmax)
      REAL G(Nmax)
      INTEGER Idim(Nidim) , Inc(Ninc) , Kb(Nkb)
!     .. Local Scalars ..
      COMPLEX transl
      REAL err , errmax
      INTEGER i , icd , ict , icu , ik , in , incx , incxs , ix , k ,   &
     &        ks , laa , lda , ldas , lx , n , nargs , nc , nk , ns
      LOGICAL banded , full , null , packed , reset , same
      CHARACTER*1 diag , diags , trans , transs , uplo , uplos
      CHARACTER*2 ichd , ichu
      CHARACTER*3 icht
!     .. Local Arrays ..
      LOGICAL isame(13)
!     .. External Functions ..
      LOGICAL LCE , LCERES
      EXTERNAL LCE , LCERES
!     .. External Subroutines ..
      EXTERNAL CMAKE , CMVCH , CTBMV , CTBSV , CTPMV , CTPSV , CTRMV ,  &
     &         CTRSV
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
      errmax = RZERO
!     Set up zero vector for CMVCH.
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
                        CALL CMAKE(Sname(2:3),uplo,diag,n,n,A,Nmax,Aa,  &
     &                             lda,k,k,reset,transl)
!
                        DO ix = 1 , Ninc
                           incx = Inc(ix)
                           lx = ABS(incx)*n
!
!                       Generate the vector X.
!
                           transl = HALF
                           CALL CMAKE('GE',' ',' ',1,n,X,1,Xx,ABS(incx),&
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
                                 CALL CTRMV(uplo,trans,diag,n,Aa,lda,Xx,&
     &                              incx)
                              ELSEIF ( banded ) THEN
                                 IF ( Trace ) WRITE (Ntra,FMT=99006)    &
     &                                nc , Sname , uplo , trans , diag ,&
     &                                n , k , lda , incx
                                 IF ( Rewi ) REWIND Ntra
                                 CALL CTBMV(uplo,trans,diag,n,k,Aa,lda, &
     &                              Xx,incx)
                              ELSEIF ( packed ) THEN
                                 IF ( Trace ) WRITE (Ntra,FMT=99005)    &
     &                                nc , Sname , uplo , trans , diag ,&
     &                                n , incx
                                 IF ( Rewi ) REWIND Ntra
                                 CALL CTPMV(uplo,trans,diag,n,Aa,Xx,    &
     &                              incx)
                              ENDIF
                           ELSEIF ( Sname(4:5)=='SV' ) THEN
                              IF ( full ) THEN
                                 IF ( Trace ) WRITE (Ntra,FMT=99007)    &
     &                                nc , Sname , uplo , trans , diag ,&
     &                                n , lda , incx
                                 IF ( Rewi ) REWIND Ntra
                                 CALL CTRSV(uplo,trans,diag,n,Aa,lda,Xx,&
     &                              incx)
                              ELSEIF ( banded ) THEN
                                 IF ( Trace ) WRITE (Ntra,FMT=99006)    &
     &                                nc , Sname , uplo , trans , diag ,&
     &                                n , k , lda , incx
                                 IF ( Rewi ) REWIND Ntra
                                 CALL CTBSV(uplo,trans,diag,n,k,Aa,lda, &
     &                              Xx,incx)
                              ELSEIF ( packed ) THEN
                                 IF ( Trace ) WRITE (Ntra,FMT=99005)    &
     &                                nc , Sname , uplo , trans , diag ,&
     &                                n , incx
                                 IF ( Rewi ) REWIND Ntra
                                 CALL CTPSV(uplo,trans,diag,n,Aa,Xx,    &
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
                              isame(5) = LCE(As,Aa,laa)
                              isame(6) = ldas==lda
                              IF ( null ) THEN
                                 isame(7) = LCE(Xs,Xx,lx)
                              ELSE
                                 isame(7)                               &
     &                              = LCERES('GE',' ',1,n,Xs,Xx,ABS     &
     &                              (incx))
                              ENDIF
                              isame(8) = incxs==incx
                           ELSEIF ( banded ) THEN
                              isame(5) = ks==k
                              isame(6) = LCE(As,Aa,laa)
                              isame(7) = ldas==lda
                              IF ( null ) THEN
                                 isame(8) = LCE(Xs,Xx,lx)
                              ELSE
                                 isame(8)                               &
     &                              = LCERES('GE',' ',1,n,Xs,Xx,ABS     &
     &                              (incx))
                              ENDIF
                              isame(9) = incxs==incx
                           ELSEIF ( packed ) THEN
                              isame(5) = LCE(As,Aa,laa)
                              IF ( null ) THEN
                                 isame(6) = LCE(Xs,Xx,lx)
                              ELSE
                                 isame(6)                               &
     &                              = LCERES('GE',' ',1,n,Xs,Xx,ABS     &
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
                              CALL CMVCH(trans,n,n,ONE,A,Nmax,X,incx,   &
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
                              CALL CMVCH(trans,n,n,ONE,A,Nmax,Z,incx,   &
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
     &        ')                                      .')
99006 FORMAT (1X,I6,': ',A6,'(',3('''',A1,''','),2(I3,','),' A,',I3,    &
     &        ', X,',I2,')                               .')
99007 FORMAT (1X,I6,': ',A6,'(',3('''',A1,''','),I3,', A,',I3,', X,',I2,&
     &        ')                                   .')
99008 FORMAT (' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *',&
     &        '******')
!
!     End of CCHK3.
!
      END SUBROUTINE CCHK3
!*==cchk4.f90  processed by SPAG 7.51RB at 17:52 on  4 Mar 2022
      SUBROUTINE CCHK4(Sname,Eps,Thresh,Nout,Ntra,Trace,Rewi,Fatal,     &
     &                 Nidim,Idim,Nalf,Alf,Ninc,Inc,Nmax,Incmax,A,Aa,As,&
     &                 X,Xx,Xs,Y,Yy,Ys,Yt,G,Z)
      IMPLICIT NONE
!*--CCHK41439
!
!  Tests CGERC and CGERU.
!
!  Auxiliary routine for test program for Level 2 Blas.
!
!  -- Written on 10-August-1987.
!     Richard Hanson, Sandia National Labs.
!     Jeremy Du Croz, NAG Central Office.
!
!     .. Parameters ..
      COMPLEX ZERO , HALF , ONE
      PARAMETER (ZERO=(0.0,0.0),HALF=(0.5,0.0),ONE=(1.0,0.0))
      REAL RZERO
      PARAMETER (RZERO=0.0)
!     .. Scalar Arguments ..
      REAL Eps , Thresh
      INTEGER Incmax , Nalf , Nidim , Ninc , Nmax , Nout , Ntra
      LOGICAL Fatal , Rewi , Trace
      CHARACTER*6 Sname
!     .. Array Arguments ..
      COMPLEX A(Nmax,Nmax) , Aa(Nmax*Nmax) , Alf(Nalf) , As(Nmax*Nmax) ,&
     &        X(Nmax) , Xs(Nmax*Incmax) , Xx(Nmax*Incmax) , Y(Nmax) ,   &
     &        Ys(Nmax*Incmax) , Yt(Nmax) , Yy(Nmax*Incmax) , Z(Nmax)
      REAL G(Nmax)
      INTEGER Idim(Nidim) , Inc(Ninc)
!     .. Local Scalars ..
      COMPLEX alpha , als , transl
      REAL err , errmax
      INTEGER i , ia , im , in , incx , incxs , incy , incys , ix , iy ,&
     &        j , laa , lda , ldas , lx , ly , m , ms , n , nargs , nc ,&
     &        nd , ns
      LOGICAL conj , null , reset , same
!     .. Local Arrays ..
      COMPLEX w(1)
      LOGICAL isame(13)
!     .. External Functions ..
      LOGICAL LCE , LCERES
      EXTERNAL LCE , LCERES
!     .. External Subroutines ..
      EXTERNAL CGERC , CGERU , CMAKE , CMVCH
!     .. Intrinsic Functions ..
      INTRINSIC ABS , CONJG , MAX , MIN
!     .. Scalars in Common ..
      INTEGER INFot , NOUtc
      LOGICAL LERr , OK
!     .. Common blocks ..
      COMMON /INFOC / INFot , NOUtc , OK , LERr
!     .. Executable Statements ..
      conj = Sname(5:5)=='C'
!     Define the number of arguments.
      nargs = 9
!
      nc = 0
      reset = .TRUE.
      errmax = RZERO
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
                  CALL CMAKE('GE',' ',' ',1,m,X,1,Xx,ABS(incx),0,m-1,   &
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
                     CALL CMAKE('GE',' ',' ',1,n,Y,1,Yy,ABS(incy),0,n-1,&
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
                        CALL CMAKE(Sname(2:3),' ',' ',m,n,A,Nmax,Aa,lda,&
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
                        IF ( conj ) THEN
                           IF ( Rewi ) REWIND Ntra
                           CALL CGERC(m,n,alpha,Xx,incx,Yy,incy,Aa,lda)
                        ELSE
                           IF ( Rewi ) REWIND Ntra
                           CALL CGERU(m,n,alpha,Xx,incx,Yy,incy,Aa,lda)
                        ENDIF
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
                        isame(4) = LCE(Xs,Xx,lx)
                        isame(5) = incxs==incx
                        isame(6) = LCE(Ys,Yy,ly)
                        isame(7) = incys==incy
                        IF ( null ) THEN
                           isame(8) = LCE(As,Aa,laa)
                        ELSE
                           isame(8) = LCERES('GE',' ',m,n,As,Aa,lda)
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
                           IF ( conj ) w(1) = CONJG(w(1))
                           CALL CMVCH('N',m,1,alpha,Z,Nmax,w,1,ONE,     &
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
99006 FORMAT (1X,I6,': ',A6,'(',2(I3,','),'(',F4.1,',',F4.1,'), X,',I2, &
     &        ', Y,',I2,', A,',I3,')                   ','      .')
99007 FORMAT (' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *',&
     &        '******')
!
!     End of CCHK4.
!
      END SUBROUTINE CCHK4
!*==cchk5.f90  processed by SPAG 7.51RB at 17:52 on  4 Mar 2022
      SUBROUTINE CCHK5(Sname,Eps,Thresh,Nout,Ntra,Trace,Rewi,Fatal,     &
     &                 Nidim,Idim,Nalf,Alf,Ninc,Inc,Nmax,Incmax,A,Aa,As,&
     &                 X,Xx,Xs,Y,Yy,Ys,Yt,G,Z)
      IMPLICIT NONE
!*--CCHK51696
!
!  Tests CHER and CHPR.
!
!  Auxiliary routine for test program for Level 2 Blas.
!
!  -- Written on 10-August-1987.
!     Richard Hanson, Sandia National Labs.
!     Jeremy Du Croz, NAG Central Office.
!
!     .. Parameters ..
      COMPLEX ZERO , HALF , ONE
      PARAMETER (ZERO=(0.0,0.0),HALF=(0.5,0.0),ONE=(1.0,0.0))
      REAL RZERO
      PARAMETER (RZERO=0.0)
!     .. Scalar Arguments ..
      REAL Eps , Thresh
      INTEGER Incmax , Nalf , Nidim , Ninc , Nmax , Nout , Ntra
      LOGICAL Fatal , Rewi , Trace
      CHARACTER*6 Sname
!     .. Array Arguments ..
      COMPLEX A(Nmax,Nmax) , Aa(Nmax*Nmax) , Alf(Nalf) , As(Nmax*Nmax) ,&
     &        X(Nmax) , Xs(Nmax*Incmax) , Xx(Nmax*Incmax) , Y(Nmax) ,   &
     &        Ys(Nmax*Incmax) , Yt(Nmax) , Yy(Nmax*Incmax) , Z(Nmax)
      REAL G(Nmax)
      INTEGER Idim(Nidim) , Inc(Ninc)
!     .. Local Scalars ..
      COMPLEX alpha , transl
      REAL err , errmax , ralpha , rals
      INTEGER i , ia , ic , in , incx , incxs , ix , j , ja , jj , laa ,&
     &        lda , ldas , lj , lx , n , nargs , nc , ns
      LOGICAL full , null , packed , reset , same , upper
      CHARACTER*1 uplo , uplos
      CHARACTER*2 ich
!     .. Local Arrays ..
      COMPLEX w(1)
      LOGICAL isame(13)
!     .. External Functions ..
      LOGICAL LCE , LCERES
      EXTERNAL LCE , LCERES
!     .. External Subroutines ..
      EXTERNAL CHER , CHPR , CMAKE , CMVCH
!     .. Intrinsic Functions ..
      INTRINSIC ABS , CMPLX , CONJG , MAX , REAL
!     .. Scalars in Common ..
      INTEGER INFot , NOUtc
      LOGICAL LERr , OK
!     .. Common blocks ..
      COMMON /INFOC / INFot , NOUtc , OK , LERr
!     .. Data statements ..
      DATA ich/'UL'/
!     .. Executable Statements ..
      full = Sname(3:3)=='E'
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
      errmax = RZERO
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
                  CALL CMAKE('GE',' ',' ',1,n,X,1,Xx,ABS(incx),0,n-1,   &
     &                       reset,transl)
                  IF ( n>1 ) THEN
                     X(n/2) = ZERO
                     Xx(1+ABS(incx)*(n/2-1)) = ZERO
                  ENDIF
!
                  DO ia = 1 , Nalf
                     ralpha = REAL(Alf(ia))
                     alpha = CMPLX(ralpha,RZERO)
                     null = n<=0 .OR. ralpha==RZERO
!
!                 Generate the matrix A.
!
                     transl = ZERO
                     CALL CMAKE(Sname(2:3),uplo,' ',n,n,A,Nmax,Aa,lda,  &
     &                          n-1,n-1,reset,transl)
!
                     nc = nc + 1
!
!                 Save every datum before calling the subroutine.
!
                     uplos = uplo
                     ns = n
                     rals = ralpha
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
     &                       uplo , n , ralpha , incx , lda
                        IF ( Rewi ) REWIND Ntra
                        CALL CHER(uplo,n,ralpha,Xx,incx,Aa,lda)
                     ELSEIF ( packed ) THEN
                        IF ( Trace ) WRITE (Ntra,FMT=99006) nc , Sname ,&
     &                       uplo , n , ralpha , incx
                        IF ( Rewi ) REWIND Ntra
                        CALL CHPR(uplo,n,ralpha,Xx,incx,Aa)
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
                     isame(3) = rals==ralpha
                     isame(4) = LCE(Xs,Xx,lx)
                     isame(5) = incxs==incx
                     IF ( null ) THEN
                        isame(6) = LCE(As,Aa,laa)
                     ELSE
                        isame(6) = LCERES(Sname(2:3),uplo,n,n,As,Aa,lda)
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
                           w(1) = CONJG(Z(j))
                           IF ( upper ) THEN
                              jj = 1
                              lj = j
                           ELSE
                              jj = j
                              lj = n - j + 1
                           ENDIF
                           CALL CMVCH('N',lj,1,alpha,Z(jj),lj,w,1,ONE,  &
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
         WRITE (Nout,FMT=99007) nc , Sname , uplo , n , ralpha , incx , &
     &                          lda
      ELSEIF ( packed ) THEN
         WRITE (Nout,FMT=99006) nc , Sname , uplo , n , ralpha , incx
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
     &        ', AP)                                         .')
99007 FORMAT (1X,I6,': ',A6,'(''',A1,''',',I3,',',F4.1,', X,',I2,', A,',&
     &        I3,')                                      .')
99008 FORMAT (' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *',&
     &        '******')
!
!     End of CCHK5.
!
      END SUBROUTINE CCHK5
!*==cchk6.f90  processed by SPAG 7.51RB at 17:52 on  4 Mar 2022
      SUBROUTINE CCHK6(Sname,Eps,Thresh,Nout,Ntra,Trace,Rewi,Fatal,     &
     &                 Nidim,Idim,Nalf,Alf,Ninc,Inc,Nmax,Incmax,A,Aa,As,&
     &                 X,Xx,Xs,Y,Yy,Ys,Yt,G,Z)
      IMPLICIT NONE
!*--CCHK61963
!
!  Tests CHER2 and CHPR2.
!
!  Auxiliary routine for test program for Level 2 Blas.
!
!  -- Written on 10-August-1987.
!     Richard Hanson, Sandia National Labs.
!     Jeremy Du Croz, NAG Central Office.
!
!     .. Parameters ..
      COMPLEX ZERO , HALF , ONE
      PARAMETER (ZERO=(0.0,0.0),HALF=(0.5,0.0),ONE=(1.0,0.0))
      REAL RZERO
      PARAMETER (RZERO=0.0)
!     .. Scalar Arguments ..
      REAL Eps , Thresh
      INTEGER Incmax , Nalf , Nidim , Ninc , Nmax , Nout , Ntra
      LOGICAL Fatal , Rewi , Trace
      CHARACTER*6 Sname
!     .. Array Arguments ..
      COMPLEX A(Nmax,Nmax) , Aa(Nmax*Nmax) , Alf(Nalf) , As(Nmax*Nmax) ,&
     &        X(Nmax) , Xs(Nmax*Incmax) , Xx(Nmax*Incmax) , Y(Nmax) ,   &
     &        Ys(Nmax*Incmax) , Yt(Nmax) , Yy(Nmax*Incmax) , Z(Nmax,2)
      REAL G(Nmax)
      INTEGER Idim(Nidim) , Inc(Ninc)
!     .. Local Scalars ..
      COMPLEX alpha , als , transl
      REAL err , errmax
      INTEGER i , ia , ic , in , incx , incxs , incy , incys , ix , iy ,&
     &        j , ja , jj , laa , lda , ldas , lj , lx , ly , n ,       &
     &        nargs , nc , ns
      LOGICAL full , null , packed , reset , same , upper
      CHARACTER*1 uplo , uplos
      CHARACTER*2 ich
!     .. Local Arrays ..
      COMPLEX w(2)
      LOGICAL isame(13)
!     .. External Functions ..
      LOGICAL LCE , LCERES
      EXTERNAL LCE , LCERES
!     .. External Subroutines ..
      EXTERNAL CHER2 , CHPR2 , CMAKE , CMVCH
!     .. Intrinsic Functions ..
      INTRINSIC ABS , CONJG , MAX
!     .. Scalars in Common ..
      INTEGER INFot , NOUtc
      LOGICAL LERr , OK
!     .. Common blocks ..
      COMMON /INFOC / INFot , NOUtc , OK , LERr
!     .. Data statements ..
      DATA ich/'UL'/
!     .. Executable Statements ..
      full = Sname(3:3)=='E'
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
      errmax = RZERO
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
                  CALL CMAKE('GE',' ',' ',1,n,X,1,Xx,ABS(incx),0,n-1,   &
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
                     CALL CMAKE('GE',' ',' ',1,n,Y,1,Yy,ABS(incy),0,n-1,&
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
                        CALL CMAKE(Sname(2:3),uplo,' ',n,n,A,Nmax,Aa,   &
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
                           CALL CHER2(uplo,n,alpha,Xx,incx,Yy,incy,Aa,  &
     &                                lda)
                        ELSEIF ( packed ) THEN
                           IF ( Trace ) WRITE (Ntra,FMT=99006) nc ,     &
     &                          Sname , uplo , n , alpha , incx , incy
                           IF ( Rewi ) REWIND Ntra
                           CALL CHPR2(uplo,n,alpha,Xx,incx,Yy,incy,Aa)
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
                        isame(4) = LCE(Xs,Xx,lx)
                        isame(5) = incxs==incx
                        isame(6) = LCE(Ys,Yy,ly)
                        isame(7) = incys==incy
                        IF ( null ) THEN
                           isame(8) = LCE(As,Aa,laa)
                        ELSE
                           isame(8) = LCERES(Sname(2:3),uplo,n,n,As,Aa, &
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
                              w(1) = alpha*CONJG(Z(j,2))
                              w(2) = CONJG(alpha)*CONJG(Z(j,1))
                              IF ( upper ) THEN
                                 jj = 1
                                 lj = j
                              ELSE
                                 jj = j
                                 lj = n - j + 1
                              ENDIF
                              CALL CMVCH('N',lj,2,ONE,Z(jj,1),Nmax,w,1, &
     &                           ONE,A(jj,j),1,Yt,G,Aa(ja),Eps,err,     &
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
99006 FORMAT (1X,I6,': ',A6,'(''',A1,''',',I3,',(',F4.1,',',F4.1,       &
     &        '), X,',I2,', Y,',I2,', AP)                     ',        &
     &        '       .')
99007 FORMAT (1X,I6,': ',A6,'(''',A1,''',',I3,',(',F4.1,',',F4.1,       &
     &        '), X,',I2,', Y,',I2,', A,',I3,')             ',          &
     &        '            .')
99008 FORMAT (' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *',&
     &        '******')
!
!     End of CCHK6.
!
      END SUBROUTINE CCHK6
!*==cchke.f90  processed by SPAG 7.51RB at 17:52 on  4 Mar 2022
      SUBROUTINE CCHKE(Isnum,Srnamt,Nout)
      IMPLICIT NONE
!*--CCHKE2266
!
!  Tests the error exits from the Level 2 Blas.
!  Requires a special version of the error-handling routine XERBLA.
!  ALPHA, RALPHA, BETA, A, X and Y should not need to be defined.
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
      COMPLEX alpha , beta
      REAL ralpha
!     .. Local Arrays ..
      COMPLEX a(1,1) , x(1) , y(1)
!     .. External Subroutines ..
      EXTERNAL CGBMV , CGEMV , CGERC , CGERU , CHBMV , CHEMV , CHER ,   &
     &         CHER2 , CHKXER , CHPMV , CHPR , CHPR2 , CTBMV , CTBSV ,  &
     &         CTPMV , CTPSV , CTRMV , CTRSV
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
         CALL CGBMV('/',0,0,0,0,alpha,a,1,x,1,beta,y,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 2
         CALL CGBMV('N',-1,0,0,0,alpha,a,1,x,1,beta,y,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 3
         CALL CGBMV('N',0,-1,0,0,alpha,a,1,x,1,beta,y,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 4
         CALL CGBMV('N',0,0,-1,0,alpha,a,1,x,1,beta,y,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 5
         CALL CGBMV('N',2,0,0,-1,alpha,a,1,x,1,beta,y,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 8
         CALL CGBMV('N',0,0,1,0,alpha,a,1,x,1,beta,y,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 10
         CALL CGBMV('N',0,0,0,0,alpha,a,1,x,0,beta,y,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 13
         CALL CGBMV('N',0,0,0,0,alpha,a,1,x,1,beta,y,0)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
      ELSEIF ( Isnum==3 ) THEN
         INFot = 1
         CALL CHEMV('/',0,alpha,a,1,x,1,beta,y,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 2
         CALL CHEMV('U',-1,alpha,a,1,x,1,beta,y,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 5
         CALL CHEMV('U',2,alpha,a,1,x,1,beta,y,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 7
         CALL CHEMV('U',0,alpha,a,1,x,0,beta,y,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 10
         CALL CHEMV('U',0,alpha,a,1,x,1,beta,y,0)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
      ELSEIF ( Isnum==4 ) THEN
         INFot = 1
         CALL CHBMV('/',0,0,alpha,a,1,x,1,beta,y,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 2
         CALL CHBMV('U',-1,0,alpha,a,1,x,1,beta,y,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 3
         CALL CHBMV('U',0,-1,alpha,a,1,x,1,beta,y,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 6
         CALL CHBMV('U',0,1,alpha,a,1,x,1,beta,y,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 8
         CALL CHBMV('U',0,0,alpha,a,1,x,0,beta,y,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 11
         CALL CHBMV('U',0,0,alpha,a,1,x,1,beta,y,0)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
      ELSEIF ( Isnum==5 ) THEN
         INFot = 1
         CALL CHPMV('/',0,alpha,a,x,1,beta,y,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 2
         CALL CHPMV('U',-1,alpha,a,x,1,beta,y,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 6
         CALL CHPMV('U',0,alpha,a,x,0,beta,y,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 9
         CALL CHPMV('U',0,alpha,a,x,1,beta,y,0)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
      ELSEIF ( Isnum==6 ) THEN
         INFot = 1
         CALL CTRMV('/','N','N',0,a,1,x,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 2
         CALL CTRMV('U','/','N',0,a,1,x,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 3
         CALL CTRMV('U','N','/',0,a,1,x,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 4
         CALL CTRMV('U','N','N',-1,a,1,x,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 6
         CALL CTRMV('U','N','N',2,a,1,x,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 8
         CALL CTRMV('U','N','N',0,a,1,x,0)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
      ELSEIF ( Isnum==7 ) THEN
         INFot = 1
         CALL CTBMV('/','N','N',0,0,a,1,x,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 2
         CALL CTBMV('U','/','N',0,0,a,1,x,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 3
         CALL CTBMV('U','N','/',0,0,a,1,x,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 4
         CALL CTBMV('U','N','N',-1,0,a,1,x,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 5
         CALL CTBMV('U','N','N',0,-1,a,1,x,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 7
         CALL CTBMV('U','N','N',0,1,a,1,x,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 9
         CALL CTBMV('U','N','N',0,0,a,1,x,0)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
      ELSEIF ( Isnum==8 ) THEN
         INFot = 1
         CALL CTPMV('/','N','N',0,a,x,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 2
         CALL CTPMV('U','/','N',0,a,x,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 3
         CALL CTPMV('U','N','/',0,a,x,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 4
         CALL CTPMV('U','N','N',-1,a,x,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 7
         CALL CTPMV('U','N','N',0,a,x,0)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
      ELSEIF ( Isnum==9 ) THEN
         INFot = 1
         CALL CTRSV('/','N','N',0,a,1,x,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 2
         CALL CTRSV('U','/','N',0,a,1,x,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 3
         CALL CTRSV('U','N','/',0,a,1,x,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 4
         CALL CTRSV('U','N','N',-1,a,1,x,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 6
         CALL CTRSV('U','N','N',2,a,1,x,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 8
         CALL CTRSV('U','N','N',0,a,1,x,0)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
      ELSEIF ( Isnum==10 ) THEN
         INFot = 1
         CALL CTBSV('/','N','N',0,0,a,1,x,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 2
         CALL CTBSV('U','/','N',0,0,a,1,x,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 3
         CALL CTBSV('U','N','/',0,0,a,1,x,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 4
         CALL CTBSV('U','N','N',-1,0,a,1,x,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 5
         CALL CTBSV('U','N','N',0,-1,a,1,x,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 7
         CALL CTBSV('U','N','N',0,1,a,1,x,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 9
         CALL CTBSV('U','N','N',0,0,a,1,x,0)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
      ELSEIF ( Isnum==11 ) THEN
         INFot = 1
         CALL CTPSV('/','N','N',0,a,x,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 2
         CALL CTPSV('U','/','N',0,a,x,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 3
         CALL CTPSV('U','N','/',0,a,x,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 4
         CALL CTPSV('U','N','N',-1,a,x,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 7
         CALL CTPSV('U','N','N',0,a,x,0)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
      ELSEIF ( Isnum==12 ) THEN
         INFot = 1
         CALL CGERC(-1,0,alpha,x,1,y,1,a,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 2
         CALL CGERC(0,-1,alpha,x,1,y,1,a,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 5
         CALL CGERC(0,0,alpha,x,0,y,1,a,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 7
         CALL CGERC(0,0,alpha,x,1,y,0,a,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 9
         CALL CGERC(2,0,alpha,x,1,y,1,a,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
      ELSEIF ( Isnum==13 ) THEN
         INFot = 1
         CALL CGERU(-1,0,alpha,x,1,y,1,a,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 2
         CALL CGERU(0,-1,alpha,x,1,y,1,a,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 5
         CALL CGERU(0,0,alpha,x,0,y,1,a,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 7
         CALL CGERU(0,0,alpha,x,1,y,0,a,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 9
         CALL CGERU(2,0,alpha,x,1,y,1,a,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
      ELSEIF ( Isnum==14 ) THEN
         INFot = 1
         CALL CHER('/',0,ralpha,x,1,a,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 2
         CALL CHER('U',-1,ralpha,x,1,a,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 5
         CALL CHER('U',0,ralpha,x,0,a,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 7
         CALL CHER('U',2,ralpha,x,1,a,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
      ELSEIF ( Isnum==15 ) THEN
         INFot = 1
         CALL CHPR('/',0,ralpha,x,1,a)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 2
         CALL CHPR('U',-1,ralpha,x,1,a)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 5
         CALL CHPR('U',0,ralpha,x,0,a)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
      ELSEIF ( Isnum==16 ) THEN
         INFot = 1
         CALL CHER2('/',0,alpha,x,1,y,1,a,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 2
         CALL CHER2('U',-1,alpha,x,1,y,1,a,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 5
         CALL CHER2('U',0,alpha,x,0,y,1,a,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 7
         CALL CHER2('U',0,alpha,x,1,y,0,a,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 9
         CALL CHER2('U',2,alpha,x,1,y,1,a,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
      ELSEIF ( Isnum==17 ) THEN
         INFot = 1
         CALL CHPR2('/',0,alpha,x,1,y,1,a)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 2
         CALL CHPR2('U',-1,alpha,x,1,y,1,a)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 5
         CALL CHPR2('U',0,alpha,x,0,y,1,a)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 7
         CALL CHPR2('U',0,alpha,x,1,y,0,a)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
      ELSE
         INFot = 1
         CALL CGEMV('/',0,0,alpha,a,1,x,1,beta,y,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 2
         CALL CGEMV('N',-1,0,alpha,a,1,x,1,beta,y,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 3
         CALL CGEMV('N',0,-1,alpha,a,1,x,1,beta,y,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 6
         CALL CGEMV('N',2,0,alpha,a,1,x,1,beta,y,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 8
         CALL CGEMV('N',0,0,alpha,a,1,x,0,beta,y,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 11
         CALL CGEMV('N',0,0,alpha,a,1,x,1,beta,y,0)
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
!     End of CCHKE.
!
      END SUBROUTINE CCHKE
!*==cmake.f90  processed by SPAG 7.51RB at 17:52 on  4 Mar 2022
      SUBROUTINE CMAKE(Type,Uplo,Diag,M,N,A,Nmax,Aa,Lda,Kl,Ku,Reset,    &
     &                 Transl)
      IMPLICIT NONE
!*--CMAKE2612
!
!  Generates values for an M by N matrix A within the bandwidth
!  defined by KL and KU.
!  Stores the values in the array AA in the data structure required
!  by the routine, with unwanted elements set to rogue value.
!
!  TYPE is 'GE', 'GB', 'HE', 'HB', 'HP', 'TR', 'TB' OR 'TP'.
!
!  Auxiliary routine for test program for Level 2 Blas.
!
!  -- Written on 10-August-1987.
!     Richard Hanson, Sandia National Labs.
!     Jeremy Du Croz, NAG Central Office.
!
!     .. Parameters ..
      COMPLEX ZERO , ONE
      PARAMETER (ZERO=(0.0,0.0),ONE=(1.0,0.0))
      COMPLEX ROGUE
      PARAMETER (ROGUE=(-1.0E10,1.0E10))
      REAL RZERO
      PARAMETER (RZERO=0.0)
      REAL RROGUE
      PARAMETER (RROGUE=-1.0E10)
!     .. Scalar Arguments ..
      COMPLEX Transl
      INTEGER Kl , Ku , Lda , M , N , Nmax
      LOGICAL Reset
      CHARACTER*1 Diag , Uplo
      CHARACTER*2 Type
!     .. Array Arguments ..
      COMPLEX A(Nmax,*) , Aa(*)
!     .. Local Scalars ..
      INTEGER i , i1 , i2 , i3 , ibeg , iend , ioff , j , jj , kk
      LOGICAL gen , lower , sym , tri , unit , upper
!     .. External Functions ..
      COMPLEX CBEG
      EXTERNAL CBEG
!     .. Intrinsic Functions ..
      INTRINSIC CMPLX , CONJG , MAX , MIN , REAL
!     .. Executable Statements ..
      gen = Type(1:1)=='G'
      sym = Type(1:1)=='H'
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
                  A(i,j) = CBEG(Reset) + Transl
               ELSE
                  A(i,j) = ZERO
               ENDIF
               IF ( i/=j ) THEN
                  IF ( sym ) THEN
                     A(j,i) = CONJG(A(i,j))
                  ELSEIF ( tri ) THEN
                     A(j,i) = ZERO
                  ENDIF
               ENDIF
            ENDIF
         ENDDO
         IF ( sym ) A(j,j) = CMPLX(REAL(A(j,j)),RZERO)
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
      ELSEIF ( Type=='HE' .OR. Type=='TR' ) THEN
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
            IF ( sym ) THEN
               jj = j + (j-1)*Lda
               Aa(jj) = CMPLX(REAL(Aa(jj)),RROGUE)
            ENDIF
         ENDDO
      ELSEIF ( Type=='HB' .OR. Type=='TB' ) THEN
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
            IF ( sym ) THEN
               jj = kk + (j-1)*Lda
               Aa(jj) = CMPLX(REAL(Aa(jj)),RROGUE)
            ENDIF
         ENDDO
      ELSEIF ( Type=='HP' .OR. Type=='TP' ) THEN
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
                  IF ( sym ) Aa(ioff) = CMPLX(REAL(Aa(ioff)),RROGUE)
               ENDIF
            ENDDO
         ENDDO
      ENDIF
!
!     End of CMAKE.
!
      END SUBROUTINE CMAKE
!*==cmvch.f90  processed by SPAG 7.51RB at 17:52 on  4 Mar 2022
      SUBROUTINE CMVCH(Trans,M,N,Alpha,A,Nmax,X,Incx,Beta,Y,Incy,Yt,G,  &
     &                 Yy,Eps,Err,Fatal,Nout,Mv)
      IMPLICIT NONE
!*--CMVCH2801
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
      COMPLEX ZERO
      PARAMETER (ZERO=(0.0,0.0))
      REAL RZERO , RONE
      PARAMETER (RZERO=0.0,RONE=1.0)
!     .. Scalar Arguments ..
      COMPLEX Alpha , Beta
      REAL Eps , Err
      INTEGER Incx , Incy , M , N , Nmax , Nout
      LOGICAL Fatal , Mv
      CHARACTER*1 Trans
!     .. Array Arguments ..
      COMPLEX A(Nmax,*) , X(*) , Y(*) , Yt(*) , Yy(*)
      REAL G(*)
!     .. Local Scalars ..
      COMPLEX c
      REAL erri
      INTEGER i , incxl , incyl , iy , j , jx , kx , ky , ml , nl
      LOGICAL ctran , tran
!     .. Intrinsic Functions ..
      INTRINSIC ABS , AIMAG , CONJG , MAX , REAL , SQRT
!     .. Statement Functions ..
      REAL ABS1
!     .. Statement Function definitions ..
      ABS1(c) = ABS(REAL(c)) + ABS(AIMAG(c))
!     .. Executable Statements ..
      tran = Trans=='T'
      ctran = Trans=='C'
      IF ( tran .OR. ctran ) THEN
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
         G(iy) = RZERO
         jx = kx
         IF ( tran ) THEN
            DO j = 1 , nl
               Yt(iy) = Yt(iy) + A(j,i)*X(jx)
               G(iy) = G(iy) + ABS1(A(j,i))*ABS1(X(jx))
               jx = jx + incxl
            ENDDO
         ELSEIF ( ctran ) THEN
            DO j = 1 , nl
               Yt(iy) = Yt(iy) + CONJG(A(j,i))*X(jx)
               G(iy) = G(iy) + ABS1(A(j,i))*ABS1(X(jx))
               jx = jx + incxl
            ENDDO
         ELSE
            DO j = 1 , nl
               Yt(iy) = Yt(iy) + A(i,j)*X(jx)
               G(iy) = G(iy) + ABS1(A(i,j))*ABS1(X(jx))
               jx = jx + incxl
            ENDDO
         ENDIF
         Yt(iy) = Alpha*Yt(iy) + Beta*Y(iy)
         G(iy) = ABS1(Alpha)*G(iy) + ABS1(Beta)*ABS1(Y(iy))
         iy = iy + incyl
      ENDDO
!
!     Compute the error ratio for this result.
!
      Err = ZERO
      DO i = 1 , ml
         erri = ABS(Yt(i)-Yy(1+(i-1)*ABS(Incy)))/Eps
         IF ( G(i)/=RZERO ) erri = erri/G(i)
         Err = MAX(Err,erri)
         IF ( Err*SQRT(Eps)>=RONE ) GOTO 100
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
     &        /'                       EXPECTED RE',                    &
     &        'SULT                    COMPUTED RESULT')
99002 FORMAT (1X,I7,2('  (',G15.6,',',G15.6,')'))
!
!     End of CMVCH.
!
      END SUBROUTINE CMVCH
!*==lce.f90  processed by SPAG 7.51RB at 17:52 on  4 Mar 2022
      LOGICAL FUNCTION LCE(Ri,Rj,Lr)
      IMPLICIT NONE
!*--LCE2931
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
      COMPLEX Ri(*) , Rj(*)
!     .. Local Scalars ..
      INTEGER i
!     .. Executable Statements ..
      DO i = 1 , Lr
         IF ( Ri(i)/=Rj(i) ) GOTO 100
      ENDDO
      LCE = .TRUE.
      GOTO 99999
 100  LCE = .FALSE.
!
!     End of LCE.
!
99999 END FUNCTION LCE
!*==lceres.f90  processed by SPAG 7.51RB at 17:52 on  4 Mar 2022
      LOGICAL FUNCTION LCERES(Type,Uplo,M,N,Aa,As,Lda)
      IMPLICIT NONE
!*--LCERES2961
!
!  Tests if selected elements in two arrays are equal.
!
!  TYPE is 'GE', 'HE' or 'HP'.
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
      COMPLEX Aa(Lda,*) , As(Lda,*)
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
      ELSEIF ( Type=='HE' ) THEN
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
      LCERES = .TRUE.
      GOTO 99999
 100  LCERES = .FALSE.
!
!     End of LCERES.
!
99999 END FUNCTION LCERES
!*==cbeg.f90  processed by SPAG 7.51RB at 17:52 on  4 Mar 2022
      COMPLEX FUNCTION CBEG(Reset)
      IMPLICIT NONE
!*--CBEG3018
!
!  Generates complex numbers as pairs of random numbers uniformly
!  distributed between -0.5 and 0.5.
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
      INTEGER i , ic , j , mi , mj
!     .. Save statement ..
      SAVE i , ic , j , mi , mj
!     .. Intrinsic Functions ..
      INTRINSIC CMPLX
!     .. Executable Statements ..
      IF ( Reset ) THEN
!        Initialize local variables.
         mi = 891
         mj = 457
         i = 7
         j = 7
         ic = 0
         Reset = .FALSE.
      ENDIF
!
!     The sequence of values of I or J is bounded between 1 and 999.
!     If initial I or J = 1,2,3,6,7 or 9, the period will be 50.
!     If initial I or J = 4 or 8, the period will be 25.
!     If initial I or J = 5, the period will be 10.
!     IC is used to break up the period by skipping 1 value of I or J
!     in 6.
!
      ic = ic + 1
      DO
         i = i*mi
         j = j*mj
         i = i - 1000*(i/1000)
         j = j - 1000*(j/1000)
         IF ( ic>=5 ) THEN
            ic = 0
            CYCLE
         ENDIF
         CBEG = CMPLX((i-500)/1001.0,(j-500)/1001.0)
         EXIT
      ENDDO
!
!     End of CBEG.
!
      END FUNCTION CBEG
!*==sdiff.f90  processed by SPAG 7.51RB at 17:52 on  4 Mar 2022
      REAL FUNCTION SDIFF(X,Y)
      IMPLICIT NONE
!*--SDIFF3075
!
!  Auxiliary routine for test program for Level 2 Blas.
!
!  -- Written on 10-August-1987.
!     Richard Hanson, Sandia National Labs.
!
!     .. Scalar Arguments ..
      REAL X , Y
!     .. Executable Statements ..
      SDIFF = X - Y
!
!     End of SDIFF.
!
      END FUNCTION SDIFF
!*==chkxer.f90  processed by SPAG 7.51RB at 17:52 on  4 Mar 2022
      SUBROUTINE CHKXER(Srnamt,Infot,Nout,Lerr,Ok)
      IMPLICIT NONE
!*--CHKXER3093
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
!*--XERBLA3124
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
