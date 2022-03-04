!*==cblat3.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
 
!> \brief \b CBLAT3
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       PROGRAM CBLAT3
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> Test program for the COMPLEX          Level 3 Blas.
!>
!> The program must be driven by a short data file. The first 14 records
!> of the file are read using list-directed input, the last 9 records
!> are read using the format ( A6, L2 ). An annotated example of a data
!> file can be obtained by deleting the first 3 characters from the
!> following 23 lines:
!> 'cblat3.out'      NAME OF SUMMARY OUTPUT FILE
!> 6                 UNIT NUMBER OF SUMMARY FILE
!> 'CBLAT3.SNAP'     NAME OF SNAPSHOT OUTPUT FILE
!> -1                UNIT NUMBER OF SNAPSHOT FILE (NOT USED IF .LT. 0)
!> F        LOGICAL FLAG, T TO REWIND SNAPSHOT FILE AFTER EACH RECORD.
!> F        LOGICAL FLAG, T TO STOP ON FAILURES.
!> T        LOGICAL FLAG, T TO TEST ERROR EXITS.
!> 16.0     THRESHOLD VALUE OF TEST RATIO
!> 6                 NUMBER OF VALUES OF N
!> 0 1 2 3 5 9       VALUES OF N
!> 3                 NUMBER OF VALUES OF ALPHA
!> (0.0,0.0) (1.0,0.0) (0.7,-0.9)       VALUES OF ALPHA
!> 3                 NUMBER OF VALUES OF BETA
!> (0.0,0.0) (1.0,0.0) (1.3,-1.1)       VALUES OF BETA
!> CGEMM  T PUT F FOR NO TEST. SAME COLUMNS.
!> CHEMM  T PUT F FOR NO TEST. SAME COLUMNS.
!> CSYMM  T PUT F FOR NO TEST. SAME COLUMNS.
!> CTRMM  T PUT F FOR NO TEST. SAME COLUMNS.
!> CTRSM  T PUT F FOR NO TEST. SAME COLUMNS.
!> CHERK  T PUT F FOR NO TEST. SAME COLUMNS.
!> CSYRK  T PUT F FOR NO TEST. SAME COLUMNS.
!> CHER2K T PUT F FOR NO TEST. SAME COLUMNS.
!> CSYR2K T PUT F FOR NO TEST. SAME COLUMNS.
!>
!> Further Details
!> ===============
!>
!> See:
!>
!>    Dongarra J. J., Du Croz J. J., Duff I. S. and Hammarling S.
!>    A Set of Level 3 Basic Linear Algebra Subprograms.
!>
!>    Technical Memorandum No.88 (Revision 1), Mathematics and
!>    Computer Science Division, Argonne National Laboratory, 9700
!>    South Cass Avenue, Argonne, Illinois 60439, US.
!>
!> -- Written on 8-February-1989.
!>    Jack Dongarra, Argonne National Laboratory.
!>    Iain Duff, AERE Harwell.
!>    Jeremy Du Croz, Numerical Algorithms Group Ltd.
!>    Sven Hammarling, Numerical Algorithms Group Ltd.
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
      PROGRAM CBLAT3
      IMPLICIT NONE
!*--CBLAT390
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
      PARAMETER (NSUBS=9)
      COMPLEX ZERO , ONE
      PARAMETER (ZERO=(0.0,0.0),ONE=(1.0,0.0))
      REAL RZERO
      PARAMETER (RZERO=0.0)
      INTEGER NMAX
      PARAMETER (NMAX=65)
      INTEGER NIDMAX , NALMAX , NBEMAX
      PARAMETER (NIDMAX=9,NALMAX=7,NBEMAX=7)
!     .. Local Scalars ..
      REAL eps , err , thresh
      INTEGER i , isnum , j , n , nalf , nbet , nidim , nout , ntra
      LOGICAL fatal , ltestt , rewi , same , sfatal , trace , tsterr
      CHARACTER*1 transa , transb
      CHARACTER*6 snamet
      CHARACTER*32 snaps , summry
!     .. Local Arrays ..
      COMPLEX aa(NMAX*NMAX) , ab(NMAX,2*NMAX) , alf(NALMAX) ,           &
     &        as(NMAX*NMAX) , bb(NMAX*NMAX) , bet(NBEMAX) ,             &
     &        bs(NMAX*NMAX) , c(NMAX,NMAX) , cc(NMAX*NMAX) ,            &
     &        cs(NMAX*NMAX) , ct(NMAX) , w(2*NMAX)
      REAL g(NMAX)
      INTEGER idim(NIDMAX)
      LOGICAL ltest(NSUBS)
      CHARACTER*6 snames(NSUBS)
!     .. External Functions ..
      REAL SDIFF
      LOGICAL LCE
      EXTERNAL SDIFF , LCE
!     .. External Subroutines ..
      EXTERNAL CCHK1 , CCHK2 , CCHK3 , CCHK4 , CCHK5 , CCHKE , CMMCH
!     .. Intrinsic Functions ..
      INTRINSIC MAX , MIN
!     .. Scalars in Common ..
      INTEGER INFot , NOUtc
      LOGICAL LERr , OK
      CHARACTER*6 SRNamt
!     .. Common blocks ..
      COMMON /INFOC / INFot , NOUtc , OK , LERr
      COMMON /SRNAMC/ SRNamt
!     .. Data statements ..
      DATA snames/'CGEMM ' , 'CHEMM ' , 'CSYMM ' , 'CTRMM ' , 'CTRSM ' ,&
     &     'CHERK ' , 'CSYRK ' , 'CHER2K' , 'CSYR2K'/
!     .. Executable Statements ..
!
!     Read name and unit number for summary output file and open file.
!
      READ (NIN,FMT=*) summry
      READ (NIN,FMT=*) nout
      OPEN (nout,FILE=summry)
      NOUtc = nout
!
!     Read name and unit number for snapshot output file and open file.
!
      READ (NIN,FMT=*) snaps
      READ (NIN,FMT=*) ntra
      trace = ntra>=0
      IF ( trace ) OPEN (ntra,FILE=snaps)
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
      WRITE (nout,FMT=99005)
      WRITE (nout,FMT=99006) (idim(i),i=1,nidim)
      WRITE (nout,FMT=99007) (alf(i),i=1,nalf)
      WRITE (nout,FMT=99008) (bet(i),i=1,nbet)
      IF ( .NOT.tsterr ) THEN
         WRITE (nout,FMT=*)
         WRITE (nout,FMT=99016)
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
 100  READ (NIN,FMT=99012,END=300) snamet , ltestt
      DO i = 1 , NSUBS
         IF ( snamet==snames(i) ) GOTO 200
      ENDDO
      WRITE (nout,FMT=99010) snamet
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
!     Check the reliability of CMMCH using exact data.
!
      n = MIN(32,NMAX)
      DO j = 1 , n
         DO i = 1 , n
            ab(i,j) = MAX(i-j+1,0)
         ENDDO
         ab(j,NMAX+1) = j
         ab(1,NMAX+j) = j
         c(j,1) = ZERO
      ENDDO
      DO j = 1 , n
         cc(j) = j*((j+1)*j)/2 - ((j+1)*j*(j-1))/3
      ENDDO
!     CC holds the exact result. On exit from CMMCH CT holds
!     the result computed by CMMCH.
      transa = 'N'
      transb = 'N'
      CALL CMMCH(transa,transb,n,1,n,ONE,ab,NMAX,ab(1,NMAX+1),NMAX,ZERO,&
     &           c,NMAX,ct,g,cc,NMAX,eps,err,fatal,nout,.TRUE.)
      same = LCE(cc,ct,n)
      IF ( .NOT.same .OR. err/=RZERO ) THEN
         WRITE (nout,FMT=99011) transa , transb , same , err
         STOP
      ENDIF
      transb = 'C'
      CALL CMMCH(transa,transb,n,1,n,ONE,ab,NMAX,ab(1,NMAX+1),NMAX,ZERO,&
     &           c,NMAX,ct,g,cc,NMAX,eps,err,fatal,nout,.TRUE.)
      same = LCE(cc,ct,n)
      IF ( .NOT.same .OR. err/=RZERO ) THEN
         WRITE (nout,FMT=99011) transa , transb , same , err
         STOP
      ENDIF
      DO j = 1 , n
         ab(j,NMAX+1) = n - j + 1
         ab(1,NMAX+j) = n - j + 1
      ENDDO
      DO j = 1 , n
         cc(n-j+1) = j*((j+1)*j)/2 - ((j+1)*j*(j-1))/3
      ENDDO
      transa = 'C'
      transb = 'N'
      CALL CMMCH(transa,transb,n,1,n,ONE,ab,NMAX,ab(1,NMAX+1),NMAX,ZERO,&
     &           c,NMAX,ct,g,cc,NMAX,eps,err,fatal,nout,.TRUE.)
      same = LCE(cc,ct,n)
      IF ( .NOT.same .OR. err/=RZERO ) THEN
         WRITE (nout,FMT=99011) transa , transb , same , err
         STOP
      ENDIF
      transb = 'C'
      CALL CMMCH(transa,transb,n,1,n,ONE,ab,NMAX,ab(1,NMAX+1),NMAX,ZERO,&
     &           c,NMAX,ct,g,cc,NMAX,eps,err,fatal,nout,.TRUE.)
      same = LCE(cc,ct,n)
      IF ( .NOT.same .OR. err/=RZERO ) THEN
         WRITE (nout,FMT=99011) transa , transb , same , err
         STOP
      ENDIF
!
!     Test each subroutine in turn.
!
      DO isnum = 1 , NSUBS
         WRITE (nout,FMT=*)
         IF ( .NOT.ltest(isnum) ) THEN
!           Subprogram is not to be tested.
            WRITE (nout,FMT=99013) snames(isnum)
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
            IF ( isnum==2 .OR. isnum==3 ) THEN
!           Test CHEMM, 02, CSYMM, 03.
               CALL CCHK2(snames(isnum),eps,thresh,nout,ntra,trace,rewi,&
     &                    fatal,nidim,idim,nalf,alf,nbet,bet,NMAX,ab,aa,&
     &                    as,ab(1,NMAX+1),bb,bs,c,cc,cs,ct,g)
            ELSEIF ( isnum==4 .OR. isnum==5 ) THEN
!           Test CTRMM, 04, CTRSM, 05.
               CALL CCHK3(snames(isnum),eps,thresh,nout,ntra,trace,rewi,&
     &                    fatal,nidim,idim,nalf,alf,NMAX,ab,aa,as,      &
     &                    ab(1,NMAX+1),bb,bs,ct,g,c)
            ELSEIF ( isnum==6 .OR. isnum==7 ) THEN
!           Test CHERK, 06, CSYRK, 07.
               CALL CCHK4(snames(isnum),eps,thresh,nout,ntra,trace,rewi,&
     &                    fatal,nidim,idim,nalf,alf,nbet,bet,NMAX,ab,aa,&
     &                    as,ab(1,NMAX+1),bb,bs,c,cc,cs,ct,g)
            ELSEIF ( isnum==8 .OR. isnum==9 ) THEN
!           Test CHER2K, 08, CSYR2K, 09.
               CALL CCHK5(snames(isnum),eps,thresh,nout,ntra,trace,rewi,&
     &                    fatal,nidim,idim,nalf,alf,nbet,bet,NMAX,ab,aa,&
     &                    as,bb,bs,c,cc,cs,ct,g,w)
            ELSE
!           Test CGEMM, 01.
               CALL CCHK1(snames(isnum),eps,thresh,nout,ntra,trace,rewi,&
     &                    fatal,nidim,idim,nalf,alf,nbet,bet,NMAX,ab,aa,&
     &                    as,ab(1,NMAX+1),bb,bs,c,cc,cs,ct,g)
            ENDIF
!
            IF ( fatal .AND. sfatal ) GOTO 400
         ENDIF
      ENDDO
      WRITE (nout,FMT=99014)
      GOTO 600
!
 400  WRITE (nout,FMT=99015)
      GOTO 600
!
 500  WRITE (nout,FMT=99009)
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
99005 FORMAT (' TESTS OF THE COMPLEX          LEVEL 3 BLAS',//' THE F', &
     &        'OLLOWING PARAMETER VALUES WILL BE USED:')
99006 FORMAT ('   FOR N              ',9I6)
99007 FORMAT ('   FOR ALPHA          ',7('(',F4.1,',',F4.1,')  ',:))
99008 FORMAT ('   FOR BETA           ',7('(',F4.1,',',F4.1,')  ',:))
99009 FORMAT (' AMEND DATA FILE OR INCREASE ARRAY SIZES IN PROGRAM',    &
     &        /' ******* TESTS ABANDONED *******')
99010 FORMAT (' SUBPROGRAM NAME ',A6,' NOT RECOGNIZED',/' ******* T',   &
     &        'ESTS ABANDONED *******')
99011 FORMAT (' ERROR IN CMMCH -  IN-LINE DOT PRODUCTS ARE BEING EVALU',&
     &        'ATED WRONGLY.',/' CMMCH WAS CALLED WITH TRANSA = ',A1,   &
     &        ' AND TRANSB = ',A1,/' AND RETURNED SAME = ',L1,' AND ',  &
     &        'ERR = ',F12.3,'.',/' THIS MAY BE DUE TO FAULTS IN THE ', &
     &        'ARITHMETIC OR THE COMPILER.',                            &
     &        /' ******* TESTS ABANDONED ','*******')
99012 FORMAT (A6,L2)
99013 FORMAT (1X,A6,' WAS NOT TESTED')
99014 FORMAT (/' END OF TESTS')
99015 FORMAT (/' ******* FATAL ERROR - TESTS ABANDONED *******')
99016 FORMAT (' ERROR-EXITS WILL NOT BE TESTED')
!
!     End of CBLAT3.
!
      END PROGRAM CBLAT3
!*==cchk1.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
      SUBROUTINE CCHK1(Sname,Eps,Thresh,Nout,Ntra,Trace,Rewi,Fatal,     &
     &                 Nidim,Idim,Nalf,Alf,Nbet,Bet,Nmax,A,Aa,As,B,Bb,  &
     &                 Bs,C,Cc,Cs,Ct,G)
      IMPLICIT NONE
!*--CCHK1390
!
!  Tests CGEMM.
!
!  Auxiliary routine for test program for Level 3 Blas.
!
!  -- Written on 8-February-1989.
!     Jack Dongarra, Argonne National Laboratory.
!     Iain Duff, AERE Harwell.
!     Jeremy Du Croz, Numerical Algorithms Group Ltd.
!     Sven Hammarling, Numerical Algorithms Group Ltd.
!
!     .. Parameters ..
      COMPLEX ZERO
      PARAMETER (ZERO=(0.0,0.0))
      REAL RZERO
      PARAMETER (RZERO=0.0)
!     .. Scalar Arguments ..
      REAL Eps , Thresh
      INTEGER Nalf , Nbet , Nidim , Nmax , Nout , Ntra
      LOGICAL Fatal , Rewi , Trace
      CHARACTER*6 Sname
!     .. Array Arguments ..
      COMPLEX A(Nmax,Nmax) , Aa(Nmax*Nmax) , Alf(Nalf) , As(Nmax*Nmax) ,&
     &        B(Nmax,Nmax) , Bb(Nmax*Nmax) , Bet(Nbet) , Bs(Nmax*Nmax) ,&
     &        C(Nmax,Nmax) , Cc(Nmax*Nmax) , Cs(Nmax*Nmax) , Ct(Nmax)
      REAL G(Nmax)
      INTEGER Idim(Nidim)
!     .. Local Scalars ..
      COMPLEX alpha , als , beta , bls
      REAL err , errmax
      INTEGER i , ia , ib , ica , icb , ik , im , in , k , ks , laa ,   &
     &        lbb , lcc , lda , ldas , ldb , ldbs , ldc , ldcs , m ,    &
     &        ma , mb , ms , n , na , nargs , nb , nc , ns
      LOGICAL null , reset , same , trana , tranb
      CHARACTER*1 tranas , tranbs , transa , transb
      CHARACTER*3 ich
!     .. Local Arrays ..
      LOGICAL isame(13)
!     .. External Functions ..
      LOGICAL LCE , LCERES
      EXTERNAL LCE , LCERES
!     .. External Subroutines ..
      EXTERNAL CGEMM , CMAKE , CMMCH
!     .. Intrinsic Functions ..
      INTRINSIC MAX
!     .. Scalars in Common ..
      INTEGER INFot , NOUtc
      LOGICAL LERr , OK
!     .. Common blocks ..
      COMMON /INFOC / INFot , NOUtc , OK , LERr
!     .. Data statements ..
      DATA ich/'NTC'/
!     .. Executable Statements ..
!
      nargs = 13
      nc = 0
      reset = .TRUE.
      errmax = RZERO
!
      DO im = 1 , Nidim
         m = Idim(im)
!
         DO in = 1 , Nidim
            n = Idim(in)
!           Set LDC to 1 more than minimum value if room.
            ldc = m
            IF ( ldc<Nmax ) ldc = ldc + 1
!           Skip tests if not enough room.
            IF ( ldc<=Nmax ) THEN
               lcc = ldc*n
               null = n<=0 .OR. m<=0
!
               DO ik = 1 , Nidim
                  k = Idim(ik)
!
                  DO ica = 1 , 3
                     transa = ich(ica:ica)
                     trana = transa=='T' .OR. transa=='C'
!
                     IF ( trana ) THEN
                        ma = k
                        na = m
                     ELSE
                        ma = m
                        na = k
                     ENDIF
!                 Set LDA to 1 more than minimum value if room.
                     lda = ma
                     IF ( lda<Nmax ) lda = lda + 1
!                 Skip tests if not enough room.
                     IF ( lda<=Nmax ) THEN
                        laa = lda*na
!
!                 Generate the matrix A.
!
                        CALL CMAKE('GE',' ',' ',ma,na,A,Nmax,Aa,lda,    &
     &                             reset,ZERO)
!
                        DO icb = 1 , 3
                           transb = ich(icb:icb)
                           tranb = transb=='T' .OR. transb=='C'
!
                           IF ( tranb ) THEN
                              mb = n
                              nb = k
                           ELSE
                              mb = k
                              nb = n
                           ENDIF
!                    Set LDB to 1 more than minimum value if room.
                           ldb = mb
                           IF ( ldb<Nmax ) ldb = ldb + 1
!                    Skip tests if not enough room.
                           IF ( ldb<=Nmax ) THEN
                              lbb = ldb*nb
!
!                    Generate the matrix B.
!
                              CALL CMAKE('GE',' ',' ',mb,nb,B,Nmax,Bb,  &
     &                           ldb,reset,ZERO)
!
                              DO ia = 1 , Nalf
                                 alpha = Alf(ia)
!
                                 DO ib = 1 , Nbet
                                    beta = Bet(ib)
!
!                          Generate the matrix C.
!
                                    CALL CMAKE('GE',' ',' ',m,n,C,Nmax, &
     &                                 Cc,ldc,reset,ZERO)
!
                                    nc = nc + 1
!
!                          Save every datum before calling the
!                          subroutine.
!
                                    tranas = transa
                                    tranbs = transb
                                    ms = m
                                    ns = n
                                    ks = k
                                    als = alpha
                                    DO i = 1 , laa
                                       As(i) = Aa(i)
                                    ENDDO
                                    ldas = lda
                                    DO i = 1 , lbb
                                       Bs(i) = Bb(i)
                                    ENDDO
                                    ldbs = ldb
                                    bls = beta
                                    DO i = 1 , lcc
                                       Cs(i) = Cc(i)
                                    ENDDO
                                    ldcs = ldc
!
!                          Call the subroutine.
!
                                    IF ( Trace ) WRITE (Ntra,FMT=99005) &
     &                                 nc , Sname , transa , transb ,   &
     &                                 m , n , k , alpha , lda , ldb ,  &
     &                                 beta , ldc
                                    IF ( Rewi ) REWIND Ntra
                                    CALL CGEMM(transa,transb,m,n,k,     &
     &                                 alpha,Aa,lda,Bb,ldb,beta,Cc,ldc)
!
!                          Check if error-exit was taken incorrectly.
!
                                    IF ( .NOT.OK ) THEN
                                       WRITE (Nout,FMT=99006)
                                       Fatal = .TRUE.
                                       GOTO 100
                                    ENDIF
!
!                          See what data changed inside subroutines.
!
                                    isame(1) = transa==tranas
                                    isame(2) = transb==tranbs
                                    isame(3) = ms==m
                                    isame(4) = ns==n
                                    isame(5) = ks==k
                                    isame(6) = als==alpha
                                    isame(7) = LCE(As,Aa,laa)
                                    isame(8) = ldas==lda
                                    isame(9) = LCE(Bs,Bb,lbb)
                                    isame(10) = ldbs==ldb
                                    isame(11) = bls==beta
                                    IF ( null ) THEN
                                       isame(12) = LCE(Cs,Cc,lcc)
                                    ELSE
                                       isame(12)                        &
     &                                    = LCERES('GE',' ',m,n,Cs,Cc,  &
     &                                    ldc)
                                    ENDIF
                                    isame(13) = ldcs==ldc
!
!                          If data was incorrectly changed, report
!                          and return.
!
                                    same = .TRUE.
                                    DO i = 1 , nargs
                                       same = same .AND. isame(i)
                                       IF ( .NOT.isame(i) )             &
     &                                    WRITE (Nout,FMT=99002) i
                                    ENDDO
                                    IF ( .NOT.same ) THEN
                                       Fatal = .TRUE.
                                       GOTO 100
                                    ENDIF
!
                                    IF ( .NOT.null ) THEN
!
!                             Check the result.
!
                                       CALL CMMCH(transa,transb,m,n,k,  &
     &                                    alpha,A,Nmax,B,Nmax,beta,C,   &
     &                                    Nmax,Ct,G,Cc,ldc,Eps,err,     &
     &                                    Fatal,Nout,.TRUE.)
                                       errmax = MAX(errmax,err)
!                             If got really bad answer, report and
!                             return.
                                       IF ( Fatal ) GOTO 100
                                    ENDIF
!
                                 ENDDO
!
                              ENDDO
                           ENDIF
!
                        ENDDO
                     ENDIF
!
                  ENDDO
!
               ENDDO
            ENDIF
!
         ENDDO
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
      WRITE (Nout,FMT=99005) nc , Sname , transa , transb , m , n , k , &
     &                       alpha , lda , ldb , beta , ldc
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
99005 FORMAT (1X,I6,': ',A6,'(''',A1,''',''',A1,''',',3(I3,','),'(',    &
     &        F4.1,',',F4.1,'), A,',I3,', B,',I3,',(',F4.1,',',F4.1,    &
     &        '), C,',I3,').')
99006 FORMAT (' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *',&
     &        '******')
!
!     End of CCHK1.
!
      END SUBROUTINE CCHK1
!*==cchk2.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
      SUBROUTINE CCHK2(Sname,Eps,Thresh,Nout,Ntra,Trace,Rewi,Fatal,     &
     &                 Nidim,Idim,Nalf,Alf,Nbet,Bet,Nmax,A,Aa,As,B,Bb,  &
     &                 Bs,C,Cc,Cs,Ct,G)
      IMPLICIT NONE
!*--CCHK2670
!
!  Tests CHEMM and CSYMM.
!
!  Auxiliary routine for test program for Level 3 Blas.
!
!  -- Written on 8-February-1989.
!     Jack Dongarra, Argonne National Laboratory.
!     Iain Duff, AERE Harwell.
!     Jeremy Du Croz, Numerical Algorithms Group Ltd.
!     Sven Hammarling, Numerical Algorithms Group Ltd.
!
!     .. Parameters ..
      COMPLEX ZERO
      PARAMETER (ZERO=(0.0,0.0))
      REAL RZERO
      PARAMETER (RZERO=0.0)
!     .. Scalar Arguments ..
      REAL Eps , Thresh
      INTEGER Nalf , Nbet , Nidim , Nmax , Nout , Ntra
      LOGICAL Fatal , Rewi , Trace
      CHARACTER*6 Sname
!     .. Array Arguments ..
      COMPLEX A(Nmax,Nmax) , Aa(Nmax*Nmax) , Alf(Nalf) , As(Nmax*Nmax) ,&
     &        B(Nmax,Nmax) , Bb(Nmax*Nmax) , Bet(Nbet) , Bs(Nmax*Nmax) ,&
     &        C(Nmax,Nmax) , Cc(Nmax*Nmax) , Cs(Nmax*Nmax) , Ct(Nmax)
      REAL G(Nmax)
      INTEGER Idim(Nidim)
!     .. Local Scalars ..
      COMPLEX alpha , als , beta , bls
      REAL err , errmax
      INTEGER i , ia , ib , ics , icu , im , in , laa , lbb , lcc ,     &
     &        lda , ldas , ldb , ldbs , ldc , ldcs , m , ms , n , na ,  &
     &        nargs , nc , ns
      LOGICAL conj , left , null , reset , same
      CHARACTER*1 side , sides , uplo , uplos
      CHARACTER*2 ichs , ichu
!     .. Local Arrays ..
      LOGICAL isame(13)
!     .. External Functions ..
      LOGICAL LCE , LCERES
      EXTERNAL LCE , LCERES
!     .. External Subroutines ..
      EXTERNAL CHEMM , CMAKE , CMMCH , CSYMM
!     .. Intrinsic Functions ..
      INTRINSIC MAX
!     .. Scalars in Common ..
      INTEGER INFot , NOUtc
      LOGICAL LERr , OK
!     .. Common blocks ..
      COMMON /INFOC / INFot , NOUtc , OK , LERr
!     .. Data statements ..
      DATA ichs/'LR'/ , ichu/'UL'/
!     .. Executable Statements ..
      conj = Sname(2:3)=='HE'
!
      nargs = 12
      nc = 0
      reset = .TRUE.
      errmax = RZERO
!
      DO im = 1 , Nidim
         m = Idim(im)
!
         DO in = 1 , Nidim
            n = Idim(in)
!           Set LDC to 1 more than minimum value if room.
            ldc = m
            IF ( ldc<Nmax ) ldc = ldc + 1
!           Skip tests if not enough room.
            IF ( ldc<=Nmax ) THEN
               lcc = ldc*n
               null = n<=0 .OR. m<=0
!           Set LDB to 1 more than minimum value if room.
               ldb = m
               IF ( ldb<Nmax ) ldb = ldb + 1
!           Skip tests if not enough room.
               IF ( ldb<=Nmax ) THEN
                  lbb = ldb*n
!
!           Generate the matrix B.
!
                  CALL CMAKE('GE',' ',' ',m,n,B,Nmax,Bb,ldb,reset,ZERO)
!
                  DO ics = 1 , 2
                     side = ichs(ics:ics)
                     left = side=='L'
!
                     IF ( left ) THEN
                        na = m
                     ELSE
                        na = n
                     ENDIF
!              Set LDA to 1 more than minimum value if room.
                     lda = na
                     IF ( lda<Nmax ) lda = lda + 1
!              Skip tests if not enough room.
                     IF ( lda<=Nmax ) THEN
                        laa = lda*na
!
                        DO icu = 1 , 2
                           uplo = ichu(icu:icu)
!
!                 Generate the hermitian or symmetric matrix A.
!
                           CALL CMAKE(Sname(2:3),uplo,' ',na,na,A,Nmax, &
     &                                Aa,lda,reset,ZERO)
!
                           DO ia = 1 , Nalf
                              alpha = Alf(ia)
!
                              DO ib = 1 , Nbet
                                 beta = Bet(ib)
!
!                       Generate the matrix C.
!
                                 CALL CMAKE('GE',' ',' ',m,n,C,Nmax,Cc, &
     &                              ldc,reset,ZERO)
!
                                 nc = nc + 1
!
!                       Save every datum before calling the
!                       subroutine.
!
                                 sides = side
                                 uplos = uplo
                                 ms = m
                                 ns = n
                                 als = alpha
                                 DO i = 1 , laa
                                    As(i) = Aa(i)
                                 ENDDO
                                 ldas = lda
                                 DO i = 1 , lbb
                                    Bs(i) = Bb(i)
                                 ENDDO
                                 ldbs = ldb
                                 bls = beta
                                 DO i = 1 , lcc
                                    Cs(i) = Cc(i)
                                 ENDDO
                                 ldcs = ldc
!
!                       Call the subroutine.
!
                                 IF ( Trace ) WRITE (Ntra,FMT=99005)    &
     &                                nc , Sname , side , uplo , m , n ,&
     &                                alpha , lda , ldb , beta , ldc
                                 IF ( Rewi ) REWIND Ntra
                                 IF ( conj ) THEN
                                    CALL CHEMM(side,uplo,m,n,alpha,Aa,  &
     &                                 lda,Bb,ldb,beta,Cc,ldc)
                                 ELSE
                                    CALL CSYMM(side,uplo,m,n,alpha,Aa,  &
     &                                 lda,Bb,ldb,beta,Cc,ldc)
                                 ENDIF
!
!                       Check if error-exit was taken incorrectly.
!
                                 IF ( .NOT.OK ) THEN
                                    WRITE (Nout,FMT=99006)
                                    Fatal = .TRUE.
                                    GOTO 100
                                 ENDIF
!
!                       See what data changed inside subroutines.
!
                                 isame(1) = sides==side
                                 isame(2) = uplos==uplo
                                 isame(3) = ms==m
                                 isame(4) = ns==n
                                 isame(5) = als==alpha
                                 isame(6) = LCE(As,Aa,laa)
                                 isame(7) = ldas==lda
                                 isame(8) = LCE(Bs,Bb,lbb)
                                 isame(9) = ldbs==ldb
                                 isame(10) = bls==beta
                                 IF ( null ) THEN
                                    isame(11) = LCE(Cs,Cc,lcc)
                                 ELSE
                                    isame(11)                           &
     &                                 = LCERES('GE',' ',m,n,Cs,Cc,ldc)
                                 ENDIF
                                 isame(12) = ldcs==ldc
!
!                       If data was incorrectly changed, report and
!                       return.
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
                                 IF ( .NOT.null ) THEN
!
!                          Check the result.
!
                                    IF ( left ) THEN
                                       CALL CMMCH('N','N',m,n,m,alpha,A,&
     &                                    Nmax,B,Nmax,beta,C,Nmax,Ct,G, &
     &                                    Cc,ldc,Eps,err,Fatal,Nout,    &
     &                                    .TRUE.)
                                    ELSE
                                       CALL CMMCH('N','N',m,n,n,alpha,B,&
     &                                    Nmax,A,Nmax,beta,C,Nmax,Ct,G, &
     &                                    Cc,ldc,Eps,err,Fatal,Nout,    &
     &                                    .TRUE.)
                                    ENDIF
                                    errmax = MAX(errmax,err)
!                          If got really bad answer, report and
!                          return.
                                    IF ( Fatal ) GOTO 100
                                 ENDIF
!
                              ENDDO
!
                           ENDDO
!
                        ENDDO
                     ENDIF
!
                  ENDDO
               ENDIF
            ENDIF
!
         ENDDO
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
      WRITE (Nout,FMT=99005) nc , Sname , side , uplo , m , n , alpha , &
     &                       lda , ldb , beta , ldc
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
99005 FORMAT (1X,I6,': ',A6,'(',2('''',A1,''','),2(I3,','),'(',F4.1,',',&
     &        F4.1,'), A,',I3,', B,',I3,',(',F4.1,',',F4.1,'), C,',I3,  &
     &        ')    .')
99006 FORMAT (' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *',&
     &        '******')
!
!     End of CCHK2.
!
      END SUBROUTINE CCHK2
!*==cchk3.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
      SUBROUTINE CCHK3(Sname,Eps,Thresh,Nout,Ntra,Trace,Rewi,Fatal,     &
     &                 Nidim,Idim,Nalf,Alf,Nmax,A,Aa,As,B,Bb,Bs,Ct,G,C)
      IMPLICIT NONE
!*--CCHK3941
!
!  Tests CTRMM and CTRSM.
!
!  Auxiliary routine for test program for Level 3 Blas.
!
!  -- Written on 8-February-1989.
!     Jack Dongarra, Argonne National Laboratory.
!     Iain Duff, AERE Harwell.
!     Jeremy Du Croz, Numerical Algorithms Group Ltd.
!     Sven Hammarling, Numerical Algorithms Group Ltd.
!
!     .. Parameters ..
      COMPLEX ZERO , ONE
      PARAMETER (ZERO=(0.0,0.0),ONE=(1.0,0.0))
      REAL RZERO
      PARAMETER (RZERO=0.0)
!     .. Scalar Arguments ..
      REAL Eps , Thresh
      INTEGER Nalf , Nidim , Nmax , Nout , Ntra
      LOGICAL Fatal , Rewi , Trace
      CHARACTER*6 Sname
!     .. Array Arguments ..
      COMPLEX A(Nmax,Nmax) , Aa(Nmax*Nmax) , Alf(Nalf) , As(Nmax*Nmax) ,&
     &        B(Nmax,Nmax) , Bb(Nmax*Nmax) , Bs(Nmax*Nmax) ,            &
     &        C(Nmax,Nmax) , Ct(Nmax)
      REAL G(Nmax)
      INTEGER Idim(Nidim)
!     .. Local Scalars ..
      COMPLEX alpha , als
      REAL err , errmax
      INTEGER i , ia , icd , ics , ict , icu , im , in , j , laa , lbb ,&
     &        lda , ldas , ldb , ldbs , m , ms , n , na , nargs , nc ,  &
     &        ns
      LOGICAL left , null , reset , same
      CHARACTER*1 diag , diags , side , sides , tranas , transa , uplo ,&
     &            uplos
      CHARACTER*2 ichd , ichs , ichu
      CHARACTER*3 icht
!     .. Local Arrays ..
      LOGICAL isame(13)
!     .. External Functions ..
      LOGICAL LCE , LCERES
      EXTERNAL LCE , LCERES
!     .. External Subroutines ..
      EXTERNAL CMAKE , CMMCH , CTRMM , CTRSM
!     .. Intrinsic Functions ..
      INTRINSIC MAX
!     .. Scalars in Common ..
      INTEGER INFot , NOUtc
      LOGICAL LERr , OK
!     .. Common blocks ..
      COMMON /INFOC / INFot , NOUtc , OK , LERr
!     .. Data statements ..
      DATA ichu/'UL'/ , icht/'NTC'/ , ichd/'UN'/ , ichs/'LR'/
!     .. Executable Statements ..
!
      nargs = 11
      nc = 0
      reset = .TRUE.
      errmax = RZERO
!     Set up zero matrix for CMMCH.
      DO j = 1 , Nmax
         DO i = 1 , Nmax
            C(i,j) = ZERO
         ENDDO
      ENDDO
!
      DO im = 1 , Nidim
         m = Idim(im)
!
         DO in = 1 , Nidim
            n = Idim(in)
!           Set LDB to 1 more than minimum value if room.
            ldb = m
            IF ( ldb<Nmax ) ldb = ldb + 1
!           Skip tests if not enough room.
            IF ( ldb<=Nmax ) THEN
               lbb = ldb*n
               null = m<=0 .OR. n<=0
!
               DO ics = 1 , 2
                  side = ichs(ics:ics)
                  left = side=='L'
                  IF ( left ) THEN
                     na = m
                  ELSE
                     na = n
                  ENDIF
!              Set LDA to 1 more than minimum value if room.
                  lda = na
                  IF ( lda<Nmax ) lda = lda + 1
!              Skip tests if not enough room.
                  IF ( lda>Nmax ) EXIT
                  laa = lda*na
!
                  DO icu = 1 , 2
                     uplo = ichu(icu:icu)
!
                     DO ict = 1 , 3
                        transa = icht(ict:ict)
!
                        DO icd = 1 , 2
                           diag = ichd(icd:icd)
!
                           DO ia = 1 , Nalf
                              alpha = Alf(ia)
!
!                          Generate the matrix A.
!
                              CALL CMAKE('TR',uplo,diag,na,na,A,Nmax,Aa,&
     &                           lda,reset,ZERO)
!
!                          Generate the matrix B.
!
                              CALL CMAKE('GE',' ',' ',m,n,B,Nmax,Bb,ldb,&
     &                           reset,ZERO)
!
                              nc = nc + 1
!
!                          Save every datum before calling the
!                          subroutine.
!
                              sides = side
                              uplos = uplo
                              tranas = transa
                              diags = diag
                              ms = m
                              ns = n
                              als = alpha
                              DO i = 1 , laa
                                 As(i) = Aa(i)
                              ENDDO
                              ldas = lda
                              DO i = 1 , lbb
                                 Bs(i) = Bb(i)
                              ENDDO
                              ldbs = ldb
!
!                          Call the subroutine.
!
                              IF ( Sname(4:5)=='MM' ) THEN
                                 IF ( Trace ) WRITE (Ntra,FMT=99005)    &
     &                                nc , Sname , side , uplo ,        &
     &                                transa , diag , m , n , alpha ,   &
     &                                lda , ldb
                                 IF ( Rewi ) REWIND Ntra
                                 CALL CTRMM(side,uplo,transa,diag,m,n,  &
     &                              alpha,Aa,lda,Bb,ldb)
                              ELSEIF ( Sname(4:5)=='SM' ) THEN
                                 IF ( Trace ) WRITE (Ntra,FMT=99005)    &
     &                                nc , Sname , side , uplo ,        &
     &                                transa , diag , m , n , alpha ,   &
     &                                lda , ldb
                                 IF ( Rewi ) REWIND Ntra
                                 CALL CTRSM(side,uplo,transa,diag,m,n,  &
     &                              alpha,Aa,lda,Bb,ldb)
                              ENDIF
!
!                          Check if error-exit was taken incorrectly.
!
                              IF ( .NOT.OK ) THEN
                                 WRITE (Nout,FMT=99006)
                                 Fatal = .TRUE.
                                 GOTO 100
                              ENDIF
!
!                          See what data changed inside subroutines.
!
                              isame(1) = sides==side
                              isame(2) = uplos==uplo
                              isame(3) = tranas==transa
                              isame(4) = diags==diag
                              isame(5) = ms==m
                              isame(6) = ns==n
                              isame(7) = als==alpha
                              isame(8) = LCE(As,Aa,laa)
                              isame(9) = ldas==lda
                              IF ( null ) THEN
                                 isame(10) = LCE(Bs,Bb,lbb)
                              ELSE
                                 isame(10)                              &
     &                              = LCERES('GE',' ',m,n,Bs,Bb,ldb)
                              ENDIF
                              isame(11) = ldbs==ldb
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
                                 GOTO 100
                              ENDIF
!
                              IF ( .NOT.null ) THEN
                                 IF ( Sname(4:5)=='MM' ) THEN
!
!                                Check the result.
!
                                    IF ( left ) THEN
                                       CALL CMMCH(transa,'N',m,n,m,     &
     &                                    alpha,A,Nmax,B,Nmax,ZERO,C,   &
     &                                    Nmax,Ct,G,Bb,ldb,Eps,err,     &
     &                                    Fatal,Nout,.TRUE.)
                                    ELSE
                                       CALL CMMCH('N',transa,m,n,n,     &
     &                                    alpha,B,Nmax,A,Nmax,ZERO,C,   &
     &                                    Nmax,Ct,G,Bb,ldb,Eps,err,     &
     &                                    Fatal,Nout,.TRUE.)
                                    ENDIF
                                 ELSEIF ( Sname(4:5)=='SM' ) THEN
!
!                                Compute approximation to original
!                                matrix.
!
                                    DO j = 1 , n
                                       DO i = 1 , m
                                         C(i,j) = Bb(i+(j-1)*ldb)
                                         Bb(i+(j-1)*ldb) = alpha*B(i,j)
                                       ENDDO
                                    ENDDO
!
                                    IF ( left ) THEN
                                       CALL CMMCH(transa,'N',m,n,m,ONE, &
     &                                    A,Nmax,C,Nmax,ZERO,B,Nmax,Ct, &
     &                                    G,Bb,ldb,Eps,err,Fatal,Nout,  &
     &                                    .FALSE.)
                                    ELSE
                                       CALL CMMCH('N',transa,m,n,n,ONE, &
     &                                    C,Nmax,A,Nmax,ZERO,B,Nmax,Ct, &
     &                                    G,Bb,ldb,Eps,err,Fatal,Nout,  &
     &                                    .FALSE.)
                                    ENDIF
                                 ENDIF
                                 errmax = MAX(errmax,err)
!                             If got really bad answer, report and
!                             return.
                                 IF ( Fatal ) GOTO 100
                              ENDIF
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
      WRITE (Nout,FMT=99005) nc , Sname , side , uplo , transa , diag , &
     &                       m , n , alpha , lda , ldb
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
99005 FORMAT (1X,I6,': ',A6,'(',4('''',A1,''','),2(I3,','),'(',F4.1,',',&
     &        F4.1,'), A,',I3,', B,',I3,')         ','      .')
99006 FORMAT (' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *',&
     &        '******')
!
!     End of CCHK3.
!
      END SUBROUTINE CCHK3
!*==cchk4.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
      SUBROUTINE CCHK4(Sname,Eps,Thresh,Nout,Ntra,Trace,Rewi,Fatal,     &
     &                 Nidim,Idim,Nalf,Alf,Nbet,Bet,Nmax,A,Aa,As,B,Bb,  &
     &                 Bs,C,Cc,Cs,Ct,G)
      IMPLICIT NONE
!*--CCHK41238
!
!  Tests CHERK and CSYRK.
!
!  Auxiliary routine for test program for Level 3 Blas.
!
!  -- Written on 8-February-1989.
!     Jack Dongarra, Argonne National Laboratory.
!     Iain Duff, AERE Harwell.
!     Jeremy Du Croz, Numerical Algorithms Group Ltd.
!     Sven Hammarling, Numerical Algorithms Group Ltd.
!
!     .. Parameters ..
      COMPLEX ZERO
      PARAMETER (ZERO=(0.0,0.0))
      REAL RONE , RZERO
      PARAMETER (RONE=1.0,RZERO=0.0)
!     .. Scalar Arguments ..
      REAL Eps , Thresh
      INTEGER Nalf , Nbet , Nidim , Nmax , Nout , Ntra
      LOGICAL Fatal , Rewi , Trace
      CHARACTER*6 Sname
!     .. Array Arguments ..
      COMPLEX A(Nmax,Nmax) , Aa(Nmax*Nmax) , Alf(Nalf) , As(Nmax*Nmax) ,&
     &        B(Nmax,Nmax) , Bb(Nmax*Nmax) , Bet(Nbet) , Bs(Nmax*Nmax) ,&
     &        C(Nmax,Nmax) , Cc(Nmax*Nmax) , Cs(Nmax*Nmax) , Ct(Nmax)
      REAL G(Nmax)
      INTEGER Idim(Nidim)
!     .. Local Scalars ..
      COMPLEX alpha , als , beta , bets
      REAL err , errmax , ralpha , rals , rbeta , rbets
      INTEGER i , ia , ib , ict , icu , ik , in , j , jc , jj , k , ks ,&
     &        laa , lcc , lda , ldas , ldc , ldcs , lj , ma , n , na ,  &
     &        nargs , nc , ns
      LOGICAL conj , null , reset , same , tran , upper
      CHARACTER*1 trans , transs , transt , uplo , uplos
      CHARACTER*2 icht , ichu
!     .. Local Arrays ..
      LOGICAL isame(13)
!     .. External Functions ..
      LOGICAL LCE , LCERES
      EXTERNAL LCE , LCERES
!     .. External Subroutines ..
      EXTERNAL CHERK , CMAKE , CMMCH , CSYRK
!     .. Intrinsic Functions ..
      INTRINSIC CMPLX , MAX , REAL
!     .. Scalars in Common ..
      INTEGER INFot , NOUtc
      LOGICAL LERr , OK
!     .. Common blocks ..
      COMMON /INFOC / INFot , NOUtc , OK , LERr
!     .. Data statements ..
      DATA icht/'NC'/ , ichu/'UL'/
!     .. Executable Statements ..
      conj = Sname(2:3)=='HE'
!
      nargs = 10
      nc = 0
      reset = .TRUE.
      errmax = RZERO
!
      DO in = 1 , Nidim
         n = Idim(in)
!        Set LDC to 1 more than minimum value if room.
         ldc = n
         IF ( ldc<Nmax ) ldc = ldc + 1
!        Skip tests if not enough room.
         IF ( ldc<=Nmax ) THEN
            lcc = ldc*n
!
            DO ik = 1 , Nidim
               k = Idim(ik)
!
               DO ict = 1 , 2
                  trans = icht(ict:ict)
                  tran = trans=='C'
                  IF ( tran .AND. .NOT.conj ) trans = 'T'
                  IF ( tran ) THEN
                     ma = k
                     na = n
                  ELSE
                     ma = n
                     na = k
                  ENDIF
!              Set LDA to 1 more than minimum value if room.
                  lda = ma
                  IF ( lda<Nmax ) lda = lda + 1
!              Skip tests if not enough room.
                  IF ( lda<=Nmax ) THEN
                     laa = lda*na
!
!              Generate the matrix A.
!
                     CALL CMAKE('GE',' ',' ',ma,na,A,Nmax,Aa,lda,reset, &
     &                          ZERO)
!
                     DO icu = 1 , 2
                        uplo = ichu(icu:icu)
                        upper = uplo=='U'
!
                        DO ia = 1 , Nalf
                           alpha = Alf(ia)
                           IF ( conj ) THEN
                              ralpha = REAL(alpha)
                              alpha = CMPLX(ralpha,RZERO)
                           ENDIF
!
                           DO ib = 1 , Nbet
                              beta = Bet(ib)
                              IF ( conj ) THEN
                                 rbeta = REAL(beta)
                                 beta = CMPLX(rbeta,RZERO)
                              ENDIF
                              null = n<=0
                              IF ( conj ) null = null .OR.              &
     &                             ((k<=0 .OR. ralpha==RZERO) .AND.     &
     &                             rbeta==RONE)
!
!                       Generate the matrix C.
!
                              CALL CMAKE(Sname(2:3),uplo,' ',n,n,C,Nmax,&
     &                           Cc,ldc,reset,ZERO)
!
                              nc = nc + 1
!
!                       Save every datum before calling the subroutine.
!
                              uplos = uplo
                              transs = trans
                              ns = n
                              ks = k
                              IF ( conj ) THEN
                                 rals = ralpha
                              ELSE
                                 als = alpha
                              ENDIF
                              DO i = 1 , laa
                                 As(i) = Aa(i)
                              ENDDO
                              ldas = lda
                              IF ( conj ) THEN
                                 rbets = rbeta
                              ELSE
                                 bets = beta
                              ENDIF
                              DO i = 1 , lcc
                                 Cs(i) = Cc(i)
                              ENDDO
                              ldcs = ldc
!
!                       Call the subroutine.
!
                              IF ( conj ) THEN
                                 IF ( Trace ) WRITE (Ntra,FMT=99006)    &
     &                                nc , Sname , uplo , trans , n ,   &
     &                                k , ralpha , lda , rbeta , ldc
                                 IF ( Rewi ) REWIND Ntra
                                 CALL CHERK(uplo,trans,n,k,ralpha,Aa,   &
     &                              lda,rbeta,Cc,ldc)
                              ELSE
                                 IF ( Trace ) WRITE (Ntra,FMT=99007)    &
     &                                nc , Sname , uplo , trans , n ,   &
     &                                k , alpha , lda , beta , ldc
                                 IF ( Rewi ) REWIND Ntra
                                 CALL CSYRK(uplo,trans,n,k,alpha,Aa,lda,&
     &                              beta,Cc,ldc)
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
                              isame(1) = uplos==uplo
                              isame(2) = transs==trans
                              isame(3) = ns==n
                              isame(4) = ks==k
                              IF ( conj ) THEN
                                 isame(5) = rals==ralpha
                              ELSE
                                 isame(5) = als==alpha
                              ENDIF
                              isame(6) = LCE(As,Aa,laa)
                              isame(7) = ldas==lda
                              IF ( conj ) THEN
                                 isame(8) = rbets==rbeta
                              ELSE
                                 isame(8) = bets==beta
                              ENDIF
                              IF ( null ) THEN
                                 isame(9) = LCE(Cs,Cc,lcc)
                              ELSE
                                 isame(9)                               &
     &                              = LCERES(Sname(2:3),uplo,n,n,Cs,Cc, &
     &                              ldc)
                              ENDIF
                              isame(10) = ldcs==ldc
!
!                       If data was incorrectly changed, report and
!                       return.
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
                              IF ( .NOT.null ) THEN
!
!                          Check the result column by column.
!
                                 IF ( conj ) THEN
                                    transt = 'C'
                                 ELSE
                                    transt = 'T'
                                 ENDIF
                                 jc = 1
                                 DO j = 1 , n
                                    IF ( upper ) THEN
                                       jj = 1
                                       lj = j
                                    ELSE
                                       jj = j
                                       lj = n - j + 1
                                    ENDIF
                                    IF ( tran ) THEN
                                       CALL CMMCH(transt,'N',lj,1,k,    &
     &                                    alpha,A(1,jj),Nmax,A(1,j),    &
     &                                    Nmax,beta,C(jj,j),Nmax,Ct,G,  &
     &                                    Cc(jc),ldc,Eps,err,Fatal,Nout,&
     &                                    .TRUE.)
                                    ELSE
                                       CALL CMMCH('N',transt,lj,1,k,    &
     &                                    alpha,A(jj,1),Nmax,A(j,1),    &
     &                                    Nmax,beta,C(jj,j),Nmax,Ct,G,  &
     &                                    Cc(jc),ldc,Eps,err,Fatal,Nout,&
     &                                    .TRUE.)
                                    ENDIF
                                    IF ( upper ) THEN
                                       jc = jc + ldc
                                    ELSE
                                       jc = jc + ldc + 1
                                    ENDIF
                                    errmax = MAX(errmax,err)
!                             If got really bad answer, report and
!                             return.
                                    IF ( Fatal ) GOTO 100
                                 ENDDO
                              ENDIF
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
            ENDDO
         ENDIF
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
 100  IF ( n>1 ) WRITE (Nout,FMT=99005) j
!
 200  WRITE (Nout,FMT=99004) Sname
      IF ( conj ) THEN
         WRITE (Nout,FMT=99006) nc , Sname , uplo , trans , n , k ,     &
     &                          ralpha , lda , rbeta , ldc
      ELSE
         WRITE (Nout,FMT=99007) nc , Sname , uplo , trans , n , k ,     &
     &                          alpha , lda , beta , ldc
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
99005 FORMAT ('      THESE ARE THE RESULTS FOR COLUMN ',I3)
99006 FORMAT (1X,I6,': ',A6,'(',2('''',A1,''','),2(I3,','),F4.1,', A,', &
     &        I3,',',F4.1,', C,',I3,')               ','          .')
99007 FORMAT (1X,I6,': ',A6,'(',2('''',A1,''','),2(I3,','),'(',F4.1,',',&
     &        F4.1,') , A,',I3,',(',F4.1,',',F4.1,'), C,',I3,           &
     &        ')          .')
99008 FORMAT (' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *',&
     &        '******')
!
!     End of CCHK4.
!
      END SUBROUTINE CCHK4
!*==cchk5.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
      SUBROUTINE CCHK5(Sname,Eps,Thresh,Nout,Ntra,Trace,Rewi,Fatal,     &
     &                 Nidim,Idim,Nalf,Alf,Nbet,Bet,Nmax,Ab,Aa,As,Bb,Bs,&
     &                 C,Cc,Cs,Ct,G,W)
      IMPLICIT NONE
!*--CCHK51559
!
!  Tests CHER2K and CSYR2K.
!
!  Auxiliary routine for test program for Level 3 Blas.
!
!  -- Written on 8-February-1989.
!     Jack Dongarra, Argonne National Laboratory.
!     Iain Duff, AERE Harwell.
!     Jeremy Du Croz, Numerical Algorithms Group Ltd.
!     Sven Hammarling, Numerical Algorithms Group Ltd.
!
!     .. Parameters ..
      COMPLEX ZERO , ONE
      PARAMETER (ZERO=(0.0,0.0),ONE=(1.0,0.0))
      REAL RONE , RZERO
      PARAMETER (RONE=1.0,RZERO=0.0)
!     .. Scalar Arguments ..
      REAL Eps , Thresh
      INTEGER Nalf , Nbet , Nidim , Nmax , Nout , Ntra
      LOGICAL Fatal , Rewi , Trace
      CHARACTER*6 Sname
!     .. Array Arguments ..
      COMPLEX Aa(Nmax*Nmax) , Ab(2*Nmax*Nmax) , Alf(Nalf) ,             &
     &        As(Nmax*Nmax) , Bb(Nmax*Nmax) , Bet(Nbet) , Bs(Nmax*Nmax) &
     &        , C(Nmax,Nmax) , Cc(Nmax*Nmax) , Cs(Nmax*Nmax) , Ct(Nmax) &
     &        , W(2*Nmax)
      REAL G(Nmax)
      INTEGER Idim(Nidim)
!     .. Local Scalars ..
      COMPLEX alpha , als , beta , bets
      REAL err , errmax , rbeta , rbets
      INTEGER i , ia , ib , ict , icu , ik , in , j , jc , jj , jjab ,  &
     &        k , ks , laa , lbb , lcc , lda , ldas , ldb , ldbs , ldc ,&
     &        ldcs , lj , ma , n , na , nargs , nc , ns
      LOGICAL conj , null , reset , same , tran , upper
      CHARACTER*1 trans , transs , transt , uplo , uplos
      CHARACTER*2 icht , ichu
!     .. Local Arrays ..
      LOGICAL isame(13)
!     .. External Functions ..
      LOGICAL LCE , LCERES
      EXTERNAL LCE , LCERES
!     .. External Subroutines ..
      EXTERNAL CHER2K , CMAKE , CMMCH , CSYR2K
!     .. Intrinsic Functions ..
      INTRINSIC CMPLX , CONJG , MAX , REAL
!     .. Scalars in Common ..
      INTEGER INFot , NOUtc
      LOGICAL LERr , OK
!     .. Common blocks ..
      COMMON /INFOC / INFot , NOUtc , OK , LERr
!     .. Data statements ..
      DATA icht/'NC'/ , ichu/'UL'/
!     .. Executable Statements ..
      conj = Sname(2:3)=='HE'
!
      nargs = 12
      nc = 0
      reset = .TRUE.
      errmax = RZERO
!
      DO in = 1 , Nidim
         n = Idim(in)
!        Set LDC to 1 more than minimum value if room.
         ldc = n
         IF ( ldc<Nmax ) ldc = ldc + 1
!        Skip tests if not enough room.
         IF ( ldc<=Nmax ) THEN
            lcc = ldc*n
!
            DO ik = 1 , Nidim
               k = Idim(ik)
!
               DO ict = 1 , 2
                  trans = icht(ict:ict)
                  tran = trans=='C'
                  IF ( tran .AND. .NOT.conj ) trans = 'T'
                  IF ( tran ) THEN
                     ma = k
                     na = n
                  ELSE
                     ma = n
                     na = k
                  ENDIF
!              Set LDA to 1 more than minimum value if room.
                  lda = ma
                  IF ( lda<Nmax ) lda = lda + 1
!              Skip tests if not enough room.
                  IF ( lda<=Nmax ) THEN
                     laa = lda*na
!
!              Generate the matrix A.
!
                     IF ( tran ) THEN
                        CALL CMAKE('GE',' ',' ',ma,na,Ab,2*Nmax,Aa,lda, &
     &                             reset,ZERO)
                     ELSE
                        CALL CMAKE('GE',' ',' ',ma,na,Ab,Nmax,Aa,lda,   &
     &                             reset,ZERO)
                     ENDIF
!
!              Generate the matrix B.
!
                     ldb = lda
                     lbb = laa
                     IF ( tran ) THEN
                        CALL CMAKE('GE',' ',' ',ma,na,Ab(k+1),2*Nmax,Bb,&
     &                             ldb,reset,ZERO)
                     ELSE
                        CALL CMAKE('GE',' ',' ',ma,na,Ab(k*Nmax+1),Nmax,&
     &                             Bb,ldb,reset,ZERO)
                     ENDIF
!
                     DO icu = 1 , 2
                        uplo = ichu(icu:icu)
                        upper = uplo=='U'
!
                        DO ia = 1 , Nalf
                           alpha = Alf(ia)
!
                           DO ib = 1 , Nbet
                              beta = Bet(ib)
                              IF ( conj ) THEN
                                 rbeta = REAL(beta)
                                 beta = CMPLX(rbeta,RZERO)
                              ENDIF
                              null = n<=0
                              IF ( conj ) null = null .OR.              &
     &                             ((k<=0 .OR. alpha==ZERO) .AND.       &
     &                             rbeta==RONE)
!
!                       Generate the matrix C.
!
                              CALL CMAKE(Sname(2:3),uplo,' ',n,n,C,Nmax,&
     &                           Cc,ldc,reset,ZERO)
!
                              nc = nc + 1
!
!                       Save every datum before calling the subroutine.
!
                              uplos = uplo
                              transs = trans
                              ns = n
                              ks = k
                              als = alpha
                              DO i = 1 , laa
                                 As(i) = Aa(i)
                              ENDDO
                              ldas = lda
                              DO i = 1 , lbb
                                 Bs(i) = Bb(i)
                              ENDDO
                              ldbs = ldb
                              IF ( conj ) THEN
                                 rbets = rbeta
                              ELSE
                                 bets = beta
                              ENDIF
                              DO i = 1 , lcc
                                 Cs(i) = Cc(i)
                              ENDDO
                              ldcs = ldc
!
!                       Call the subroutine.
!
                              IF ( conj ) THEN
                                 IF ( Trace ) WRITE (Ntra,FMT=99006)    &
     &                                nc , Sname , uplo , trans , n ,   &
     &                                k , alpha , lda , ldb , rbeta ,   &
     &                                ldc
                                 IF ( Rewi ) REWIND Ntra
                                 CALL CHER2K(uplo,trans,n,k,alpha,Aa,   &
     &                              lda,Bb,ldb,rbeta,Cc,ldc)
                              ELSE
                                 IF ( Trace ) WRITE (Ntra,FMT=99007)    &
     &                                nc , Sname , uplo , trans , n ,   &
     &                                k , alpha , lda , ldb , beta , ldc
                                 IF ( Rewi ) REWIND Ntra
                                 CALL CSYR2K(uplo,trans,n,k,alpha,Aa,   &
     &                              lda,Bb,ldb,beta,Cc,ldc)
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
                              isame(1) = uplos==uplo
                              isame(2) = transs==trans
                              isame(3) = ns==n
                              isame(4) = ks==k
                              isame(5) = als==alpha
                              isame(6) = LCE(As,Aa,laa)
                              isame(7) = ldas==lda
                              isame(8) = LCE(Bs,Bb,lbb)
                              isame(9) = ldbs==ldb
                              IF ( conj ) THEN
                                 isame(10) = rbets==rbeta
                              ELSE
                                 isame(10) = bets==beta
                              ENDIF
                              IF ( null ) THEN
                                 isame(11) = LCE(Cs,Cc,lcc)
                              ELSE
                                 isame(11)                              &
     &                              = LCERES('HE',uplo,n,n,Cs,Cc,ldc)
                              ENDIF
                              isame(12) = ldcs==ldc
!
!                       If data was incorrectly changed, report and
!                       return.
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
                              IF ( .NOT.null ) THEN
!
!                          Check the result column by column.
!
                                 IF ( conj ) THEN
                                    transt = 'C'
                                 ELSE
                                    transt = 'T'
                                 ENDIF
                                 jjab = 1
                                 jc = 1
                                 DO j = 1 , n
                                    IF ( upper ) THEN
                                       jj = 1
                                       lj = j
                                    ELSE
                                       jj = j
                                       lj = n - j + 1
                                    ENDIF
                                    IF ( tran ) THEN
                                       DO i = 1 , k
                                         W(i)                           &
     &                                      = alpha*Ab((j-1)*2*Nmax+k+i)
                                         IF ( conj ) THEN
                                         W(k+i) = CONJG(alpha)          &
     &                                      *Ab((j-1)*2*Nmax+i)
                                         ELSE
                                         W(k+i)                         &
     &                                      = alpha*Ab((j-1)*2*Nmax+i)
                                         ENDIF
                                       ENDDO
                                       CALL CMMCH(transt,'N',lj,1,2*k,  &
     &                                    ONE,Ab(jjab),2*Nmax,W,2*Nmax, &
     &                                    beta,C(jj,j),Nmax,Ct,G,Cc(jc),&
     &                                    ldc,Eps,err,Fatal,Nout,.TRUE.)
                                    ELSE
                                       DO i = 1 , k
                                         IF ( conj ) THEN
                                         W(i) = alpha*CONJG(Ab((k+i-1)* &
     &                                      Nmax+j))
                                         W(k+i)                         &
     &                                      = CONJG(alpha*Ab((i-1)*Nmax+&
     &                                      j))
                                         ELSE
                                         W(i) = alpha*Ab((k+i-1)*Nmax+j)
                                         W(k+i) = alpha*Ab((i-1)*Nmax+j)
                                         ENDIF
                                       ENDDO
                                       CALL CMMCH('N','N',lj,1,2*k,ONE, &
     &                                    Ab(jj),Nmax,W,2*Nmax,beta,    &
     &                                    C(jj,j),Nmax,Ct,G,Cc(jc),ldc, &
     &                                    Eps,err,Fatal,Nout,.TRUE.)
                                    ENDIF
                                    IF ( upper ) THEN
                                       jc = jc + ldc
                                    ELSE
                                       jc = jc + ldc + 1
                                       IF ( tran ) jjab = jjab + 2*Nmax
                                    ENDIF
                                    errmax = MAX(errmax,err)
!                             If got really bad answer, report and
!                             return.
                                    IF ( Fatal ) GOTO 100
                                 ENDDO
                              ENDIF
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
            ENDDO
         ENDIF
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
 100  IF ( n>1 ) WRITE (Nout,FMT=99005) j
!
 200  WRITE (Nout,FMT=99004) Sname
      IF ( conj ) THEN
         WRITE (Nout,FMT=99006) nc , Sname , uplo , trans , n , k ,     &
     &                          alpha , lda , ldb , rbeta , ldc
      ELSE
         WRITE (Nout,FMT=99007) nc , Sname , uplo , trans , n , k ,     &
     &                          alpha , lda , ldb , beta , ldc
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
99005 FORMAT ('      THESE ARE THE RESULTS FOR COLUMN ',I3)
99006 FORMAT (1X,I6,': ',A6,'(',2('''',A1,''','),2(I3,','),'(',F4.1,',',&
     &        F4.1,'), A,',I3,', B,',I3,',',F4.1,', C,',I3,             &
     &        ')           .')
99007 FORMAT (1X,I6,': ',A6,'(',2('''',A1,''','),2(I3,','),'(',F4.1,',',&
     &        F4.1,'), A,',I3,', B,',I3,',(',F4.1,',',F4.1,'), C,',I3,  &
     &        ')    .')
99008 FORMAT (' ******* FATAL ERROR - ERROR-EXIT TAKEN ON VALID CALL *',&
     &        '******')
!
!     End of CCHK5.
!
      END SUBROUTINE CCHK5
!*==cchke.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
      SUBROUTINE CCHKE(Isnum,Srnamt,Nout)
      IMPLICIT NONE
!*--CCHKE1914
!
!  Tests the error exits from the Level 3 Blas.
!  Requires a special version of the error-handling routine XERBLA.
!  A, B and C should not need to be defined.
!
!  Auxiliary routine for test program for Level 3 Blas.
!
!  -- Written on 8-February-1989.
!     Jack Dongarra, Argonne National Laboratory.
!     Iain Duff, AERE Harwell.
!     Jeremy Du Croz, Numerical Algorithms Group Ltd.
!     Sven Hammarling, Numerical Algorithms Group Ltd.
!
!  3-19-92:  Initialize ALPHA, BETA, RALPHA, and RBETA  (eca)
!  3-19-92:  Fix argument 12 in calls to CSYMM and CHEMM
!            with INFOT = 9  (eca)
!
!     .. Scalar Arguments ..
      INTEGER Isnum , Nout
      CHARACTER*6 Srnamt
!     .. Scalars in Common ..
      INTEGER INFot , NOUtc
      LOGICAL LERr , OK
!     .. Parameters ..
      REAL ONE , TWO
      PARAMETER (ONE=1.0E0,TWO=2.0E0)
!     .. Local Scalars ..
      COMPLEX alpha , beta
      REAL ralpha , rbeta
!     .. Local Arrays ..
      COMPLEX a(2,1) , b(2,1) , c(2,1)
!     .. External Subroutines ..
      EXTERNAL CGEMM , CHEMM , CHER2K , CHERK , CHKXER , CSYMM ,        &
     &         CSYR2K , CSYRK , CTRMM , CTRSM
!     .. Common blocks ..
      COMMON /INFOC / INFot , NOUtc , OK , LERr
!     .. Executable Statements ..
!     OK is set to .FALSE. by the special version of XERBLA or by CHKXER
!     if anything is wrong.
      OK = .TRUE.
!     LERR is set to .TRUE. by the special version of XERBLA each time
!     it is called, and is then tested and re-set by CHKXER.
      LERr = .FALSE.
!
!     Initialize ALPHA, BETA, RALPHA, and RBETA.
!
      alpha = CMPLX(ONE,-ONE)
      beta = CMPLX(TWO,-TWO)
      ralpha = ONE
      rbeta = TWO
!
      IF ( Isnum==2 ) THEN
         INFot = 1
         CALL CHEMM('/','U',0,0,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 2
         CALL CHEMM('L','/',0,0,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 3
         CALL CHEMM('L','U',-1,0,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 3
         CALL CHEMM('R','U',-1,0,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 3
         CALL CHEMM('L','L',-1,0,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 3
         CALL CHEMM('R','L',-1,0,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 4
         CALL CHEMM('L','U',0,-1,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 4
         CALL CHEMM('R','U',0,-1,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 4
         CALL CHEMM('L','L',0,-1,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 4
         CALL CHEMM('R','L',0,-1,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 7
         CALL CHEMM('L','U',2,0,alpha,a,1,b,2,beta,c,2)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 7
         CALL CHEMM('R','U',0,2,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 7
         CALL CHEMM('L','L',2,0,alpha,a,1,b,2,beta,c,2)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 7
         CALL CHEMM('R','L',0,2,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 9
         CALL CHEMM('L','U',2,0,alpha,a,2,b,1,beta,c,2)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 9
         CALL CHEMM('R','U',2,0,alpha,a,1,b,1,beta,c,2)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 9
         CALL CHEMM('L','L',2,0,alpha,a,2,b,1,beta,c,2)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 9
         CALL CHEMM('R','L',2,0,alpha,a,1,b,1,beta,c,2)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 12
         CALL CHEMM('L','U',2,0,alpha,a,2,b,2,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 12
         CALL CHEMM('R','U',2,0,alpha,a,1,b,2,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 12
         CALL CHEMM('L','L',2,0,alpha,a,2,b,2,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 12
         CALL CHEMM('R','L',2,0,alpha,a,1,b,2,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
      ELSEIF ( Isnum==3 ) THEN
         INFot = 1
         CALL CSYMM('/','U',0,0,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 2
         CALL CSYMM('L','/',0,0,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 3
         CALL CSYMM('L','U',-1,0,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 3
         CALL CSYMM('R','U',-1,0,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 3
         CALL CSYMM('L','L',-1,0,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 3
         CALL CSYMM('R','L',-1,0,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 4
         CALL CSYMM('L','U',0,-1,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 4
         CALL CSYMM('R','U',0,-1,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 4
         CALL CSYMM('L','L',0,-1,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 4
         CALL CSYMM('R','L',0,-1,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 7
         CALL CSYMM('L','U',2,0,alpha,a,1,b,2,beta,c,2)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 7
         CALL CSYMM('R','U',0,2,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 7
         CALL CSYMM('L','L',2,0,alpha,a,1,b,2,beta,c,2)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 7
         CALL CSYMM('R','L',0,2,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 9
         CALL CSYMM('L','U',2,0,alpha,a,2,b,1,beta,c,2)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 9
         CALL CSYMM('R','U',2,0,alpha,a,1,b,1,beta,c,2)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 9
         CALL CSYMM('L','L',2,0,alpha,a,2,b,1,beta,c,2)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 9
         CALL CSYMM('R','L',2,0,alpha,a,1,b,1,beta,c,2)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 12
         CALL CSYMM('L','U',2,0,alpha,a,2,b,2,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 12
         CALL CSYMM('R','U',2,0,alpha,a,1,b,2,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 12
         CALL CSYMM('L','L',2,0,alpha,a,2,b,2,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 12
         CALL CSYMM('R','L',2,0,alpha,a,1,b,2,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
      ELSEIF ( Isnum==4 ) THEN
         INFot = 1
         CALL CTRMM('/','U','N','N',0,0,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 2
         CALL CTRMM('L','/','N','N',0,0,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 3
         CALL CTRMM('L','U','/','N',0,0,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 4
         CALL CTRMM('L','U','N','/',0,0,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 5
         CALL CTRMM('L','U','N','N',-1,0,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 5
         CALL CTRMM('L','U','C','N',-1,0,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 5
         CALL CTRMM('L','U','T','N',-1,0,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 5
         CALL CTRMM('R','U','N','N',-1,0,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 5
         CALL CTRMM('R','U','C','N',-1,0,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 5
         CALL CTRMM('R','U','T','N',-1,0,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 5
         CALL CTRMM('L','L','N','N',-1,0,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 5
         CALL CTRMM('L','L','C','N',-1,0,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 5
         CALL CTRMM('L','L','T','N',-1,0,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 5
         CALL CTRMM('R','L','N','N',-1,0,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 5
         CALL CTRMM('R','L','C','N',-1,0,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 5
         CALL CTRMM('R','L','T','N',-1,0,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 6
         CALL CTRMM('L','U','N','N',0,-1,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 6
         CALL CTRMM('L','U','C','N',0,-1,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 6
         CALL CTRMM('L','U','T','N',0,-1,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 6
         CALL CTRMM('R','U','N','N',0,-1,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 6
         CALL CTRMM('R','U','C','N',0,-1,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 6
         CALL CTRMM('R','U','T','N',0,-1,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 6
         CALL CTRMM('L','L','N','N',0,-1,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 6
         CALL CTRMM('L','L','C','N',0,-1,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 6
         CALL CTRMM('L','L','T','N',0,-1,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 6
         CALL CTRMM('R','L','N','N',0,-1,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 6
         CALL CTRMM('R','L','C','N',0,-1,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 6
         CALL CTRMM('R','L','T','N',0,-1,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 9
         CALL CTRMM('L','U','N','N',2,0,alpha,a,1,b,2)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 9
         CALL CTRMM('L','U','C','N',2,0,alpha,a,1,b,2)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 9
         CALL CTRMM('L','U','T','N',2,0,alpha,a,1,b,2)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 9
         CALL CTRMM('R','U','N','N',0,2,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 9
         CALL CTRMM('R','U','C','N',0,2,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 9
         CALL CTRMM('R','U','T','N',0,2,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 9
         CALL CTRMM('L','L','N','N',2,0,alpha,a,1,b,2)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 9
         CALL CTRMM('L','L','C','N',2,0,alpha,a,1,b,2)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 9
         CALL CTRMM('L','L','T','N',2,0,alpha,a,1,b,2)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 9
         CALL CTRMM('R','L','N','N',0,2,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 9
         CALL CTRMM('R','L','C','N',0,2,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 9
         CALL CTRMM('R','L','T','N',0,2,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 11
         CALL CTRMM('L','U','N','N',2,0,alpha,a,2,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 11
         CALL CTRMM('L','U','C','N',2,0,alpha,a,2,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 11
         CALL CTRMM('L','U','T','N',2,0,alpha,a,2,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 11
         CALL CTRMM('R','U','N','N',2,0,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 11
         CALL CTRMM('R','U','C','N',2,0,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 11
         CALL CTRMM('R','U','T','N',2,0,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 11
         CALL CTRMM('L','L','N','N',2,0,alpha,a,2,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 11
         CALL CTRMM('L','L','C','N',2,0,alpha,a,2,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 11
         CALL CTRMM('L','L','T','N',2,0,alpha,a,2,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 11
         CALL CTRMM('R','L','N','N',2,0,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 11
         CALL CTRMM('R','L','C','N',2,0,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 11
         CALL CTRMM('R','L','T','N',2,0,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
      ELSEIF ( Isnum==5 ) THEN
         INFot = 1
         CALL CTRSM('/','U','N','N',0,0,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 2
         CALL CTRSM('L','/','N','N',0,0,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 3
         CALL CTRSM('L','U','/','N',0,0,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 4
         CALL CTRSM('L','U','N','/',0,0,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 5
         CALL CTRSM('L','U','N','N',-1,0,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 5
         CALL CTRSM('L','U','C','N',-1,0,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 5
         CALL CTRSM('L','U','T','N',-1,0,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 5
         CALL CTRSM('R','U','N','N',-1,0,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 5
         CALL CTRSM('R','U','C','N',-1,0,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 5
         CALL CTRSM('R','U','T','N',-1,0,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 5
         CALL CTRSM('L','L','N','N',-1,0,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 5
         CALL CTRSM('L','L','C','N',-1,0,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 5
         CALL CTRSM('L','L','T','N',-1,0,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 5
         CALL CTRSM('R','L','N','N',-1,0,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 5
         CALL CTRSM('R','L','C','N',-1,0,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 5
         CALL CTRSM('R','L','T','N',-1,0,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 6
         CALL CTRSM('L','U','N','N',0,-1,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 6
         CALL CTRSM('L','U','C','N',0,-1,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 6
         CALL CTRSM('L','U','T','N',0,-1,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 6
         CALL CTRSM('R','U','N','N',0,-1,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 6
         CALL CTRSM('R','U','C','N',0,-1,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 6
         CALL CTRSM('R','U','T','N',0,-1,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 6
         CALL CTRSM('L','L','N','N',0,-1,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 6
         CALL CTRSM('L','L','C','N',0,-1,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 6
         CALL CTRSM('L','L','T','N',0,-1,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 6
         CALL CTRSM('R','L','N','N',0,-1,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 6
         CALL CTRSM('R','L','C','N',0,-1,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 6
         CALL CTRSM('R','L','T','N',0,-1,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 9
         CALL CTRSM('L','U','N','N',2,0,alpha,a,1,b,2)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 9
         CALL CTRSM('L','U','C','N',2,0,alpha,a,1,b,2)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 9
         CALL CTRSM('L','U','T','N',2,0,alpha,a,1,b,2)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 9
         CALL CTRSM('R','U','N','N',0,2,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 9
         CALL CTRSM('R','U','C','N',0,2,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 9
         CALL CTRSM('R','U','T','N',0,2,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 9
         CALL CTRSM('L','L','N','N',2,0,alpha,a,1,b,2)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 9
         CALL CTRSM('L','L','C','N',2,0,alpha,a,1,b,2)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 9
         CALL CTRSM('L','L','T','N',2,0,alpha,a,1,b,2)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 9
         CALL CTRSM('R','L','N','N',0,2,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 9
         CALL CTRSM('R','L','C','N',0,2,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 9
         CALL CTRSM('R','L','T','N',0,2,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 11
         CALL CTRSM('L','U','N','N',2,0,alpha,a,2,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 11
         CALL CTRSM('L','U','C','N',2,0,alpha,a,2,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 11
         CALL CTRSM('L','U','T','N',2,0,alpha,a,2,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 11
         CALL CTRSM('R','U','N','N',2,0,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 11
         CALL CTRSM('R','U','C','N',2,0,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 11
         CALL CTRSM('R','U','T','N',2,0,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 11
         CALL CTRSM('L','L','N','N',2,0,alpha,a,2,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 11
         CALL CTRSM('L','L','C','N',2,0,alpha,a,2,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 11
         CALL CTRSM('L','L','T','N',2,0,alpha,a,2,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 11
         CALL CTRSM('R','L','N','N',2,0,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 11
         CALL CTRSM('R','L','C','N',2,0,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 11
         CALL CTRSM('R','L','T','N',2,0,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
      ELSEIF ( Isnum==6 ) THEN
         INFot = 1
         CALL CHERK('/','N',0,0,ralpha,a,1,rbeta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 2
         CALL CHERK('U','T',0,0,ralpha,a,1,rbeta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 3
         CALL CHERK('U','N',-1,0,ralpha,a,1,rbeta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 3
         CALL CHERK('U','C',-1,0,ralpha,a,1,rbeta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 3
         CALL CHERK('L','N',-1,0,ralpha,a,1,rbeta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 3
         CALL CHERK('L','C',-1,0,ralpha,a,1,rbeta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 4
         CALL CHERK('U','N',0,-1,ralpha,a,1,rbeta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 4
         CALL CHERK('U','C',0,-1,ralpha,a,1,rbeta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 4
         CALL CHERK('L','N',0,-1,ralpha,a,1,rbeta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 4
         CALL CHERK('L','C',0,-1,ralpha,a,1,rbeta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 7
         CALL CHERK('U','N',2,0,ralpha,a,1,rbeta,c,2)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 7
         CALL CHERK('U','C',0,2,ralpha,a,1,rbeta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 7
         CALL CHERK('L','N',2,0,ralpha,a,1,rbeta,c,2)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 7
         CALL CHERK('L','C',0,2,ralpha,a,1,rbeta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 10
         CALL CHERK('U','N',2,0,ralpha,a,2,rbeta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 10
         CALL CHERK('U','C',2,0,ralpha,a,1,rbeta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 10
         CALL CHERK('L','N',2,0,ralpha,a,2,rbeta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 10
         CALL CHERK('L','C',2,0,ralpha,a,1,rbeta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
      ELSEIF ( Isnum==7 ) THEN
         INFot = 1
         CALL CSYRK('/','N',0,0,alpha,a,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 2
         CALL CSYRK('U','C',0,0,alpha,a,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 3
         CALL CSYRK('U','N',-1,0,alpha,a,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 3
         CALL CSYRK('U','T',-1,0,alpha,a,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 3
         CALL CSYRK('L','N',-1,0,alpha,a,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 3
         CALL CSYRK('L','T',-1,0,alpha,a,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 4
         CALL CSYRK('U','N',0,-1,alpha,a,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 4
         CALL CSYRK('U','T',0,-1,alpha,a,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 4
         CALL CSYRK('L','N',0,-1,alpha,a,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 4
         CALL CSYRK('L','T',0,-1,alpha,a,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 7
         CALL CSYRK('U','N',2,0,alpha,a,1,beta,c,2)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 7
         CALL CSYRK('U','T',0,2,alpha,a,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 7
         CALL CSYRK('L','N',2,0,alpha,a,1,beta,c,2)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 7
         CALL CSYRK('L','T',0,2,alpha,a,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 10
         CALL CSYRK('U','N',2,0,alpha,a,2,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 10
         CALL CSYRK('U','T',2,0,alpha,a,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 10
         CALL CSYRK('L','N',2,0,alpha,a,2,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 10
         CALL CSYRK('L','T',2,0,alpha,a,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
      ELSEIF ( Isnum==8 ) THEN
         INFot = 1
         CALL CHER2K('/','N',0,0,alpha,a,1,b,1,rbeta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 2
         CALL CHER2K('U','T',0,0,alpha,a,1,b,1,rbeta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 3
         CALL CHER2K('U','N',-1,0,alpha,a,1,b,1,rbeta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 3
         CALL CHER2K('U','C',-1,0,alpha,a,1,b,1,rbeta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 3
         CALL CHER2K('L','N',-1,0,alpha,a,1,b,1,rbeta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 3
         CALL CHER2K('L','C',-1,0,alpha,a,1,b,1,rbeta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 4
         CALL CHER2K('U','N',0,-1,alpha,a,1,b,1,rbeta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 4
         CALL CHER2K('U','C',0,-1,alpha,a,1,b,1,rbeta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 4
         CALL CHER2K('L','N',0,-1,alpha,a,1,b,1,rbeta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 4
         CALL CHER2K('L','C',0,-1,alpha,a,1,b,1,rbeta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 7
         CALL CHER2K('U','N',2,0,alpha,a,1,b,1,rbeta,c,2)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 7
         CALL CHER2K('U','C',0,2,alpha,a,1,b,1,rbeta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 7
         CALL CHER2K('L','N',2,0,alpha,a,1,b,1,rbeta,c,2)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 7
         CALL CHER2K('L','C',0,2,alpha,a,1,b,1,rbeta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 9
         CALL CHER2K('U','N',2,0,alpha,a,2,b,1,rbeta,c,2)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 9
         CALL CHER2K('U','C',0,2,alpha,a,2,b,1,rbeta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 9
         CALL CHER2K('L','N',2,0,alpha,a,2,b,1,rbeta,c,2)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 9
         CALL CHER2K('L','C',0,2,alpha,a,2,b,1,rbeta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 12
         CALL CHER2K('U','N',2,0,alpha,a,2,b,2,rbeta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 12
         CALL CHER2K('U','C',2,0,alpha,a,1,b,1,rbeta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 12
         CALL CHER2K('L','N',2,0,alpha,a,2,b,2,rbeta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 12
         CALL CHER2K('L','C',2,0,alpha,a,1,b,1,rbeta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
      ELSEIF ( Isnum==9 ) THEN
         INFot = 1
         CALL CSYR2K('/','N',0,0,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 2
         CALL CSYR2K('U','C',0,0,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 3
         CALL CSYR2K('U','N',-1,0,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 3
         CALL CSYR2K('U','T',-1,0,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 3
         CALL CSYR2K('L','N',-1,0,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 3
         CALL CSYR2K('L','T',-1,0,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 4
         CALL CSYR2K('U','N',0,-1,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 4
         CALL CSYR2K('U','T',0,-1,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 4
         CALL CSYR2K('L','N',0,-1,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 4
         CALL CSYR2K('L','T',0,-1,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 7
         CALL CSYR2K('U','N',2,0,alpha,a,1,b,1,beta,c,2)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 7
         CALL CSYR2K('U','T',0,2,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 7
         CALL CSYR2K('L','N',2,0,alpha,a,1,b,1,beta,c,2)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 7
         CALL CSYR2K('L','T',0,2,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 9
         CALL CSYR2K('U','N',2,0,alpha,a,2,b,1,beta,c,2)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 9
         CALL CSYR2K('U','T',0,2,alpha,a,2,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 9
         CALL CSYR2K('L','N',2,0,alpha,a,2,b,1,beta,c,2)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 9
         CALL CSYR2K('L','T',0,2,alpha,a,2,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 12
         CALL CSYR2K('U','N',2,0,alpha,a,2,b,2,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 12
         CALL CSYR2K('U','T',2,0,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 12
         CALL CSYR2K('L','N',2,0,alpha,a,2,b,2,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 12
         CALL CSYR2K('L','T',2,0,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
      ELSE
         INFot = 1
         CALL CGEMM('/','N',0,0,0,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 1
         CALL CGEMM('/','C',0,0,0,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 1
         CALL CGEMM('/','T',0,0,0,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 2
         CALL CGEMM('N','/',0,0,0,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 2
         CALL CGEMM('C','/',0,0,0,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 2
         CALL CGEMM('T','/',0,0,0,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 3
         CALL CGEMM('N','N',-1,0,0,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 3
         CALL CGEMM('N','C',-1,0,0,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 3
         CALL CGEMM('N','T',-1,0,0,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 3
         CALL CGEMM('C','N',-1,0,0,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 3
         CALL CGEMM('C','C',-1,0,0,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 3
         CALL CGEMM('C','T',-1,0,0,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 3
         CALL CGEMM('T','N',-1,0,0,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 3
         CALL CGEMM('T','C',-1,0,0,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 3
         CALL CGEMM('T','T',-1,0,0,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 4
         CALL CGEMM('N','N',0,-1,0,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 4
         CALL CGEMM('N','C',0,-1,0,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 4
         CALL CGEMM('N','T',0,-1,0,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 4
         CALL CGEMM('C','N',0,-1,0,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 4
         CALL CGEMM('C','C',0,-1,0,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 4
         CALL CGEMM('C','T',0,-1,0,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 4
         CALL CGEMM('T','N',0,-1,0,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 4
         CALL CGEMM('T','C',0,-1,0,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 4
         CALL CGEMM('T','T',0,-1,0,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 5
         CALL CGEMM('N','N',0,0,-1,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 5
         CALL CGEMM('N','C',0,0,-1,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 5
         CALL CGEMM('N','T',0,0,-1,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 5
         CALL CGEMM('C','N',0,0,-1,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 5
         CALL CGEMM('C','C',0,0,-1,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 5
         CALL CGEMM('C','T',0,0,-1,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 5
         CALL CGEMM('T','N',0,0,-1,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 5
         CALL CGEMM('T','C',0,0,-1,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 5
         CALL CGEMM('T','T',0,0,-1,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 8
         CALL CGEMM('N','N',2,0,0,alpha,a,1,b,1,beta,c,2)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 8
         CALL CGEMM('N','C',2,0,0,alpha,a,1,b,1,beta,c,2)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 8
         CALL CGEMM('N','T',2,0,0,alpha,a,1,b,1,beta,c,2)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 8
         CALL CGEMM('C','N',0,0,2,alpha,a,1,b,2,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 8
         CALL CGEMM('C','C',0,0,2,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 8
         CALL CGEMM('C','T',0,0,2,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 8
         CALL CGEMM('T','N',0,0,2,alpha,a,1,b,2,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 8
         CALL CGEMM('T','C',0,0,2,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 8
         CALL CGEMM('T','T',0,0,2,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 10
         CALL CGEMM('N','N',0,0,2,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 10
         CALL CGEMM('C','N',0,0,2,alpha,a,2,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 10
         CALL CGEMM('T','N',0,0,2,alpha,a,2,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 10
         CALL CGEMM('N','C',0,2,0,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 10
         CALL CGEMM('C','C',0,2,0,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 10
         CALL CGEMM('T','C',0,2,0,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 10
         CALL CGEMM('N','T',0,2,0,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 10
         CALL CGEMM('C','T',0,2,0,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 10
         CALL CGEMM('T','T',0,2,0,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 13
         CALL CGEMM('N','N',2,0,0,alpha,a,2,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 13
         CALL CGEMM('N','C',2,0,0,alpha,a,2,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 13
         CALL CGEMM('N','T',2,0,0,alpha,a,2,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 13
         CALL CGEMM('C','N',2,0,0,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 13
         CALL CGEMM('C','C',2,0,0,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 13
         CALL CGEMM('C','T',2,0,0,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 13
         CALL CGEMM('T','N',2,0,0,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 13
         CALL CGEMM('T','C',2,0,0,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 13
         CALL CGEMM('T','T',2,0,0,alpha,a,1,b,1,beta,c,1)
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
!*==cmake.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
      SUBROUTINE CMAKE(Type,Uplo,Diag,M,N,A,Nmax,Aa,Lda,Reset,Transl)
      IMPLICIT NONE
!*--CMAKE2858
!
!  Generates values for an M by N matrix A.
!  Stores the values in the array AA in the data structure required
!  by the routine, with unwanted elements set to rogue value.
!
!  TYPE is 'GE', 'HE', 'SY' or 'TR'.
!
!  Auxiliary routine for test program for Level 3 Blas.
!
!  -- Written on 8-February-1989.
!     Jack Dongarra, Argonne National Laboratory.
!     Iain Duff, AERE Harwell.
!     Jeremy Du Croz, Numerical Algorithms Group Ltd.
!     Sven Hammarling, Numerical Algorithms Group Ltd.
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
      INTEGER Lda , M , N , Nmax
      LOGICAL Reset
      CHARACTER*1 Diag , Uplo
      CHARACTER*2 Type
!     .. Array Arguments ..
      COMPLEX A(Nmax,*) , Aa(*)
!     .. Local Scalars ..
      INTEGER i , ibeg , iend , j , jj
      LOGICAL gen , her , lower , sym , tri , unit , upper
!     .. External Functions ..
      COMPLEX CBEG
      EXTERNAL CBEG
!     .. Intrinsic Functions ..
      INTRINSIC CMPLX , CONJG , REAL
!     .. Executable Statements ..
      gen = Type=='GE'
      her = Type=='HE'
      sym = Type=='SY'
      tri = Type=='TR'
      upper = (her .OR. sym .OR. tri) .AND. Uplo=='U'
      lower = (her .OR. sym .OR. tri) .AND. Uplo=='L'
      unit = tri .AND. Diag=='U'
!
!     Generate data in array A.
!
      DO j = 1 , N
         DO i = 1 , M
            IF ( gen .OR. (upper .AND. i<=j) .OR. (lower .AND. i>=j) )  &
     &           THEN
               A(i,j) = CBEG(Reset) + Transl
               IF ( i/=j ) THEN
!                 Set some elements to zero
                  IF ( N>3 .AND. j==N/2 ) A(i,j) = ZERO
                  IF ( her ) THEN
                     A(j,i) = CONJG(A(i,j))
                  ELSEIF ( sym ) THEN
                     A(j,i) = A(i,j)
                  ELSEIF ( tri ) THEN
                     A(j,i) = ZERO
                  ENDIF
               ENDIF
            ENDIF
         ENDDO
         IF ( her ) A(j,j) = CMPLX(REAL(A(j,j)),RZERO)
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
      ELSEIF ( Type=='HE' .OR. Type=='SY' .OR. Type=='TR' ) THEN
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
            IF ( her ) THEN
               jj = j + (j-1)*Lda
               Aa(jj) = CMPLX(REAL(Aa(jj)),RROGUE)
            ENDIF
         ENDDO
      ENDIF
!
!     End of CMAKE.
!
      END SUBROUTINE CMAKE
!*==cmmch.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
      SUBROUTINE CMMCH(Transa,Transb,M,N,Kk,Alpha,A,Lda,B,Ldb,Beta,C,   &
     &                 Ldc,Ct,G,Cc,Ldcc,Eps,Err,Fatal,Nout,Mv)
      IMPLICIT NONE
!*--CMMCH2984
!
!  Checks the results of the computational tests.
!
!  Auxiliary routine for test program for Level 3 Blas.
!
!  -- Written on 8-February-1989.
!     Jack Dongarra, Argonne National Laboratory.
!     Iain Duff, AERE Harwell.
!     Jeremy Du Croz, Numerical Algorithms Group Ltd.
!     Sven Hammarling, Numerical Algorithms Group Ltd.
!
!     .. Parameters ..
      COMPLEX ZERO
      PARAMETER (ZERO=(0.0,0.0))
      REAL RZERO , RONE
      PARAMETER (RZERO=0.0,RONE=1.0)
!     .. Scalar Arguments ..
      COMPLEX Alpha , Beta
      REAL Eps , Err
      INTEGER Kk , Lda , Ldb , Ldc , Ldcc , M , N , Nout
      LOGICAL Fatal , Mv
      CHARACTER*1 Transa , Transb
!     .. Array Arguments ..
      COMPLEX A(Lda,*) , B(Ldb,*) , C(Ldc,*) , Cc(Ldcc,*) , Ct(*)
      REAL G(*)
!     .. Local Scalars ..
      COMPLEX cl
      REAL erri
      INTEGER i , j , k
      LOGICAL ctrana , ctranb , trana , tranb
!     .. Intrinsic Functions ..
      INTRINSIC ABS , AIMAG , CONJG , MAX , REAL , SQRT
!     .. Statement Functions ..
      REAL ABS1
!     .. Statement Function definitions ..
      ABS1(cl) = ABS(REAL(cl)) + ABS(AIMAG(cl))
!     .. Executable Statements ..
      trana = Transa=='T' .OR. Transa=='C'
      tranb = Transb=='T' .OR. Transb=='C'
      ctrana = Transa=='C'
      ctranb = Transb=='C'
!
!     Compute expected result, one column at a time, in CT using data
!     in A, B and C.
!     Compute gauges in G.
!
      DO j = 1 , N
!
         DO i = 1 , M
            Ct(i) = ZERO
            G(i) = RZERO
         ENDDO
         IF ( .NOT.trana .AND. .NOT.tranb ) THEN
            DO k = 1 , Kk
               DO i = 1 , M
                  Ct(i) = Ct(i) + A(i,k)*B(k,j)
                  G(i) = G(i) + ABS1(A(i,k))*ABS1(B(k,j))
               ENDDO
            ENDDO
         ELSEIF ( trana .AND. .NOT.tranb ) THEN
            IF ( ctrana ) THEN
               DO k = 1 , Kk
                  DO i = 1 , M
                     Ct(i) = Ct(i) + CONJG(A(k,i))*B(k,j)
                     G(i) = G(i) + ABS1(A(k,i))*ABS1(B(k,j))
                  ENDDO
               ENDDO
            ELSE
               DO k = 1 , Kk
                  DO i = 1 , M
                     Ct(i) = Ct(i) + A(k,i)*B(k,j)
                     G(i) = G(i) + ABS1(A(k,i))*ABS1(B(k,j))
                  ENDDO
               ENDDO
            ENDIF
         ELSEIF ( .NOT.trana .AND. tranb ) THEN
            IF ( ctranb ) THEN
               DO k = 1 , Kk
                  DO i = 1 , M
                     Ct(i) = Ct(i) + A(i,k)*CONJG(B(j,k))
                     G(i) = G(i) + ABS1(A(i,k))*ABS1(B(j,k))
                  ENDDO
               ENDDO
            ELSE
               DO k = 1 , Kk
                  DO i = 1 , M
                     Ct(i) = Ct(i) + A(i,k)*B(j,k)
                     G(i) = G(i) + ABS1(A(i,k))*ABS1(B(j,k))
                  ENDDO
               ENDDO
            ENDIF
         ELSEIF ( trana .AND. tranb ) THEN
            IF ( ctrana ) THEN
               IF ( ctranb ) THEN
                  DO k = 1 , Kk
                     DO i = 1 , M
                        Ct(i) = Ct(i) + CONJG(A(k,i))*CONJG(B(j,k))
                        G(i) = G(i) + ABS1(A(k,i))*ABS1(B(j,k))
                     ENDDO
                  ENDDO
               ELSE
                  DO k = 1 , Kk
                     DO i = 1 , M
                        Ct(i) = Ct(i) + CONJG(A(k,i))*B(j,k)
                        G(i) = G(i) + ABS1(A(k,i))*ABS1(B(j,k))
                     ENDDO
                  ENDDO
               ENDIF
            ELSEIF ( ctranb ) THEN
               DO k = 1 , Kk
                  DO i = 1 , M
                     Ct(i) = Ct(i) + A(k,i)*CONJG(B(j,k))
                     G(i) = G(i) + ABS1(A(k,i))*ABS1(B(j,k))
                  ENDDO
               ENDDO
            ELSE
               DO k = 1 , Kk
                  DO i = 1 , M
                     Ct(i) = Ct(i) + A(k,i)*B(j,k)
                     G(i) = G(i) + ABS1(A(k,i))*ABS1(B(j,k))
                  ENDDO
               ENDDO
            ENDIF
         ENDIF
         DO i = 1 , M
            Ct(i) = Alpha*Ct(i) + Beta*C(i,j)
            G(i) = ABS1(Alpha)*G(i) + ABS1(Beta)*ABS1(C(i,j))
         ENDDO
!
!        Compute the error ratio for this result.
!
         Err = ZERO
         DO i = 1 , M
            erri = ABS1(Ct(i)-Cc(i,j))/Eps
            IF ( G(i)/=RZERO ) erri = erri/G(i)
            Err = MAX(Err,erri)
            IF ( Err*SQRT(Eps)>=RONE ) GOTO 100
         ENDDO
!
      ENDDO
!
!     If the loop completes, all results are at least half accurate.
      GOTO 200
!
!     Report fatal error.
!
 100  Fatal = .TRUE.
      WRITE (Nout,FMT=99001)
      DO i = 1 , M
         IF ( Mv ) THEN
            WRITE (Nout,FMT=99002) i , Ct(i) , Cc(i,j)
         ELSE
            WRITE (Nout,FMT=99002) i , Cc(i,j) , Ct(i)
         ENDIF
      ENDDO
      IF ( N>1 ) WRITE (Nout,FMT=99003) j
!
 200  RETURN
!
99001 FORMAT (' ******* FATAL ERROR - COMPUTED RESULT IS LESS THAN HAL',&
     &        'F ACCURATE *******',                                     &
     &        /'                       EXPECTED RE',                    &
     &        'SULT                    COMPUTED RESULT')
99002 FORMAT (1X,I7,2('  (',G15.6,',',G15.6,')'))
99003 FORMAT ('      THESE ARE THE RESULTS FOR COLUMN ',I3)
!
!     End of CMMCH.
!
      END SUBROUTINE CMMCH
!*==lce.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
      LOGICAL FUNCTION LCE(Ri,Rj,Lr)
      IMPLICIT NONE
!*--LCE3157
!
!  Tests if two arrays are identical.
!
!  Auxiliary routine for test program for Level 3 Blas.
!
!  -- Written on 8-February-1989.
!     Jack Dongarra, Argonne National Laboratory.
!     Iain Duff, AERE Harwell.
!     Jeremy Du Croz, Numerical Algorithms Group Ltd.
!     Sven Hammarling, Numerical Algorithms Group Ltd.
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
!*==lceres.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
      LOGICAL FUNCTION LCERES(Type,Uplo,M,N,Aa,As,Lda)
      IMPLICIT NONE
!*--LCERES3189
!
!  Tests if selected elements in two arrays are equal.
!
!  TYPE is 'GE' or 'HE' or 'SY'.
!
!  Auxiliary routine for test program for Level 3 Blas.
!
!  -- Written on 8-February-1989.
!     Jack Dongarra, Argonne National Laboratory.
!     Iain Duff, AERE Harwell.
!     Jeremy Du Croz, Numerical Algorithms Group Ltd.
!     Sven Hammarling, Numerical Algorithms Group Ltd.
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
      ELSEIF ( Type=='HE' .OR. Type=='SY' ) THEN
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
!*==cbeg.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
      COMPLEX FUNCTION CBEG(Reset)
      IMPLICIT NONE
!*--CBEG3248
!
!  Generates complex numbers as pairs of random numbers uniformly
!  distributed between -0.5 and 0.5.
!
!  Auxiliary routine for test program for Level 3 Blas.
!
!  -- Written on 8-February-1989.
!     Jack Dongarra, Argonne National Laboratory.
!     Iain Duff, AERE Harwell.
!     Jeremy Du Croz, Numerical Algorithms Group Ltd.
!     Sven Hammarling, Numerical Algorithms Group Ltd.
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
!*==sdiff.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
      REAL FUNCTION SDIFF(X,Y)
      IMPLICIT NONE
!*--SDIFF3307
!
!  Auxiliary routine for test program for Level 3 Blas.
!
!  -- Written on 8-February-1989.
!     Jack Dongarra, Argonne National Laboratory.
!     Iain Duff, AERE Harwell.
!     Jeremy Du Croz, Numerical Algorithms Group Ltd.
!     Sven Hammarling, Numerical Algorithms Group Ltd.
!
!     .. Scalar Arguments ..
      REAL X , Y
!     .. Executable Statements ..
      SDIFF = X - Y
!
!     End of SDIFF.
!
      END FUNCTION SDIFF
!*==chkxer.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
      SUBROUTINE CHKXER(Srnamt,Infot,Nout,Lerr,Ok)
      IMPLICIT NONE
!*--CHKXER3328
!
!  Tests whether XERBLA has detected an error when it should.
!
!  Auxiliary routine for test program for Level 3 Blas.
!
!  -- Written on 8-February-1989.
!     Jack Dongarra, Argonne National Laboratory.
!     Iain Duff, AERE Harwell.
!     Jeremy Du Croz, Numerical Algorithms Group Ltd.
!     Sven Hammarling, Numerical Algorithms Group Ltd.
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
!*==xerbla.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
      SUBROUTINE XERBLA(Srname,Info)
      IMPLICIT NONE
!*--XERBLA3361
!
!  This is a special version of XERBLA to be used only as part of
!  the test program for testing error exits from the Level 3 BLAS
!  routines.
!
!  XERBLA  is an error handler for the Level 3 BLAS routines.
!
!  It is called by the Level 3 BLAS routines if an input parameter is
!  invalid.
!
!  Auxiliary routine for test program for Level 3 Blas.
!
!  -- Written on 8-February-1989.
!     Jack Dongarra, Argonne National Laboratory.
!     Iain Duff, AERE Harwell.
!     Jeremy Du Croz, Numerical Algorithms Group Ltd.
!     Sven Hammarling, Numerical Algorithms Group Ltd.
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
