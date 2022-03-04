!*==zblat3.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
 
!> \brief \b ZBLAT3
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       PROGRAM ZBLAT3
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> Test program for the COMPLEX*16       Level 3 Blas.
!>
!> The program must be driven by a short data file. The first 14 records
!> of the file are read using list-directed input, the last 9 records
!> are read using the format ( A6, L2 ). An annotated example of a data
!> file can be obtained by deleting the first 3 characters from the
!> following 23 lines:
!> 'zblat3.out'      NAME OF SUMMARY OUTPUT FILE
!> 6                 UNIT NUMBER OF SUMMARY FILE
!> 'ZBLAT3.SNAP'     NAME OF SNAPSHOT OUTPUT FILE
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
!> ZGEMM  T PUT F FOR NO TEST. SAME COLUMNS.
!> ZHEMM  T PUT F FOR NO TEST. SAME COLUMNS.
!> ZSYMM  T PUT F FOR NO TEST. SAME COLUMNS.
!> ZTRMM  T PUT F FOR NO TEST. SAME COLUMNS.
!> ZTRSM  T PUT F FOR NO TEST. SAME COLUMNS.
!> ZHERK  T PUT F FOR NO TEST. SAME COLUMNS.
!> ZSYRK  T PUT F FOR NO TEST. SAME COLUMNS.
!> ZHER2K T PUT F FOR NO TEST. SAME COLUMNS.
!> ZSYR2K T PUT F FOR NO TEST. SAME COLUMNS.
!>
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
!> \ingroup complex16_blas_testing
!
!  =====================================================================
      PROGRAM ZBLAT3
      IMPLICIT NONE
!*--ZBLAT391
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
      COMPLEX*16 ZERO , ONE
      PARAMETER (ZERO=(0.0D0,0.0D0),ONE=(1.0D0,0.0D0))
      DOUBLE PRECISION RZERO
      PARAMETER (RZERO=0.0D0)
      INTEGER NMAX
      PARAMETER (NMAX=65)
      INTEGER NIDMAX , NALMAX , NBEMAX
      PARAMETER (NIDMAX=9,NALMAX=7,NBEMAX=7)
!     .. Local Scalars ..
      DOUBLE PRECISION eps , err , thresh
      INTEGER i , isnum , j , n , nalf , nbet , nidim , nout , ntra
      LOGICAL fatal , ltestt , rewi , same , sfatal , trace , tsterr
      CHARACTER*1 transa , transb
      CHARACTER*6 snamet
      CHARACTER*32 snaps , summry
!     .. Local Arrays ..
      COMPLEX*16 aa(NMAX*NMAX) , ab(NMAX,2*NMAX) , alf(NALMAX) ,        &
     &           as(NMAX*NMAX) , bb(NMAX*NMAX) , bet(NBEMAX) ,          &
     &           bs(NMAX*NMAX) , c(NMAX,NMAX) , cc(NMAX*NMAX) ,         &
     &           cs(NMAX*NMAX) , ct(NMAX) , w(2*NMAX)
      DOUBLE PRECISION g(NMAX)
      INTEGER idim(NIDMAX)
      LOGICAL ltest(NSUBS)
      CHARACTER*6 snames(NSUBS)
!     .. External Functions ..
      DOUBLE PRECISION DDIFF
      LOGICAL LZE
      EXTERNAL DDIFF , LZE
!     .. External Subroutines ..
      EXTERNAL ZCHK1 , ZCHK2 , ZCHK3 , ZCHK4 , ZCHK5 , ZCHKE , ZMMCH
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
      DATA snames/'ZGEMM ' , 'ZHEMM ' , 'ZSYMM ' , 'ZTRMM ' , 'ZTRSM ' ,&
     &     'ZHERK ' , 'ZSYRK ' , 'ZHER2K' , 'ZSYR2K'/
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
!     Check the reliability of ZMMCH using exact data.
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
!     CC holds the exact result. On exit from ZMMCH CT holds
!     the result computed by ZMMCH.
      transa = 'N'
      transb = 'N'
      CALL ZMMCH(transa,transb,n,1,n,ONE,ab,NMAX,ab(1,NMAX+1),NMAX,ZERO,&
     &           c,NMAX,ct,g,cc,NMAX,eps,err,fatal,nout,.TRUE.)
      same = LZE(cc,ct,n)
      IF ( .NOT.same .OR. err/=RZERO ) THEN
         WRITE (nout,FMT=99011) transa , transb , same , err
         STOP
      ENDIF
      transb = 'C'
      CALL ZMMCH(transa,transb,n,1,n,ONE,ab,NMAX,ab(1,NMAX+1),NMAX,ZERO,&
     &           c,NMAX,ct,g,cc,NMAX,eps,err,fatal,nout,.TRUE.)
      same = LZE(cc,ct,n)
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
      CALL ZMMCH(transa,transb,n,1,n,ONE,ab,NMAX,ab(1,NMAX+1),NMAX,ZERO,&
     &           c,NMAX,ct,g,cc,NMAX,eps,err,fatal,nout,.TRUE.)
      same = LZE(cc,ct,n)
      IF ( .NOT.same .OR. err/=RZERO ) THEN
         WRITE (nout,FMT=99011) transa , transb , same , err
         STOP
      ENDIF
      transb = 'C'
      CALL ZMMCH(transa,transb,n,1,n,ONE,ab,NMAX,ab(1,NMAX+1),NMAX,ZERO,&
     &           c,NMAX,ct,g,cc,NMAX,eps,err,fatal,nout,.TRUE.)
      same = LZE(cc,ct,n)
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
               CALL ZCHKE(isnum,snames(isnum),nout)
               WRITE (nout,FMT=*)
            ENDIF
!           Test computations.
            INFot = 0
            OK = .TRUE.
            fatal = .FALSE.
            IF ( isnum==2 .OR. isnum==3 ) THEN
!           Test ZHEMM, 02, ZSYMM, 03.
               CALL ZCHK2(snames(isnum),eps,thresh,nout,ntra,trace,rewi,&
     &                    fatal,nidim,idim,nalf,alf,nbet,bet,NMAX,ab,aa,&
     &                    as,ab(1,NMAX+1),bb,bs,c,cc,cs,ct,g)
            ELSEIF ( isnum==4 .OR. isnum==5 ) THEN
!           Test ZTRMM, 04, ZTRSM, 05.
               CALL ZCHK3(snames(isnum),eps,thresh,nout,ntra,trace,rewi,&
     &                    fatal,nidim,idim,nalf,alf,NMAX,ab,aa,as,      &
     &                    ab(1,NMAX+1),bb,bs,ct,g,c)
            ELSEIF ( isnum==6 .OR. isnum==7 ) THEN
!           Test ZHERK, 06, ZSYRK, 07.
               CALL ZCHK4(snames(isnum),eps,thresh,nout,ntra,trace,rewi,&
     &                    fatal,nidim,idim,nalf,alf,nbet,bet,NMAX,ab,aa,&
     &                    as,ab(1,NMAX+1),bb,bs,c,cc,cs,ct,g)
            ELSEIF ( isnum==8 .OR. isnum==9 ) THEN
!           Test ZHER2K, 08, ZSYR2K, 09.
               CALL ZCHK5(snames(isnum),eps,thresh,nout,ntra,trace,rewi,&
     &                    fatal,nidim,idim,nalf,alf,nbet,bet,NMAX,ab,aa,&
     &                    as,bb,bs,c,cc,cs,ct,g,w)
            ELSE
!           Test ZGEMM, 01.
               CALL ZCHK1(snames(isnum),eps,thresh,nout,ntra,trace,rewi,&
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
99002 FORMAT (' RELATIVE MACHINE PRECISION IS TAKEN TO BE',1P,D9.1)
99003 FORMAT (' NUMBER OF VALUES OF ',A,' IS LESS THAN 1 OR GREATER ',  &
     &        'THAN ',I2)
99004 FORMAT (' VALUE OF N IS LESS THAN 0 OR GREATER THAN ',I2)
99005 FORMAT (' TESTS OF THE COMPLEX*16       LEVEL 3 BLAS',//' THE F', &
     &        'OLLOWING PARAMETER VALUES WILL BE USED:')
99006 FORMAT ('   FOR N              ',9I6)
99007 FORMAT ('   FOR ALPHA          ',7('(',F4.1,',',F4.1,')  ',:))
99008 FORMAT ('   FOR BETA           ',7('(',F4.1,',',F4.1,')  ',:))
99009 FORMAT (' AMEND DATA FILE OR INCREASE ARRAY SIZES IN PROGRAM',    &
     &        /' ******* TESTS ABANDONED *******')
99010 FORMAT (' SUBPROGRAM NAME ',A6,' NOT RECOGNIZED',/' ******* T',   &
     &        'ESTS ABANDONED *******')
99011 FORMAT (' ERROR IN ZMMCH -  IN-LINE DOT PRODUCTS ARE BEING EVALU',&
     &        'ATED WRONGLY.',/' ZMMCH WAS CALLED WITH TRANSA = ',A1,   &
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
!     End of ZBLAT3.
!
      END PROGRAM ZBLAT3
!*==zchk1.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
      SUBROUTINE ZCHK1(Sname,Eps,Thresh,Nout,Ntra,Trace,Rewi,Fatal,     &
     &                 Nidim,Idim,Nalf,Alf,Nbet,Bet,Nmax,A,Aa,As,B,Bb,  &
     &                 Bs,C,Cc,Cs,Ct,G)
      IMPLICIT NONE
!*--ZCHK1391
!
!  Tests ZGEMM.
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
      COMPLEX*16 ZERO
      PARAMETER (ZERO=(0.0D0,0.0D0))
      DOUBLE PRECISION RZERO
      PARAMETER (RZERO=0.0D0)
!     .. Scalar Arguments ..
      DOUBLE PRECISION Eps , Thresh
      INTEGER Nalf , Nbet , Nidim , Nmax , Nout , Ntra
      LOGICAL Fatal , Rewi , Trace
      CHARACTER*6 Sname
!     .. Array Arguments ..
      COMPLEX*16 A(Nmax,Nmax) , Aa(Nmax*Nmax) , Alf(Nalf) ,             &
     &           As(Nmax*Nmax) , B(Nmax,Nmax) , Bb(Nmax*Nmax) ,         &
     &           Bet(Nbet) , Bs(Nmax*Nmax) , C(Nmax,Nmax) ,             &
     &           Cc(Nmax*Nmax) , Cs(Nmax*Nmax) , Ct(Nmax)
      DOUBLE PRECISION G(Nmax)
      INTEGER Idim(Nidim)
!     .. Local Scalars ..
      COMPLEX*16 alpha , als , beta , bls
      DOUBLE PRECISION err , errmax
      INTEGER i , ia , ib , ica , icb , ik , im , in , k , ks , laa ,   &
     &        lbb , lcc , lda , ldas , ldb , ldbs , ldc , ldcs , m ,    &
     &        ma , mb , ms , n , na , nargs , nb , nc , ns
      LOGICAL null , reset , same , trana , tranb
      CHARACTER*1 tranas , tranbs , transa , transb
      CHARACTER*3 ich
!     .. Local Arrays ..
      LOGICAL isame(13)
!     .. External Functions ..
      LOGICAL LZE , LZERES
      EXTERNAL LZE , LZERES
!     .. External Subroutines ..
      EXTERNAL ZGEMM , ZMAKE , ZMMCH
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
                        CALL ZMAKE('GE',' ',' ',ma,na,A,Nmax,Aa,lda,    &
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
                              CALL ZMAKE('GE',' ',' ',mb,nb,B,Nmax,Bb,  &
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
                                    CALL ZMAKE('GE',' ',' ',m,n,C,Nmax, &
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
                                    CALL ZGEMM(transa,transb,m,n,k,     &
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
                                    isame(7) = LZE(As,Aa,laa)
                                    isame(8) = ldas==lda
                                    isame(9) = LZE(Bs,Bb,lbb)
                                    isame(10) = ldbs==ldb
                                    isame(11) = bls==beta
                                    IF ( null ) THEN
                                       isame(12) = LZE(Cs,Cc,lcc)
                                    ELSE
                                       isame(12)                        &
     &                                    = LZERES('GE',' ',m,n,Cs,Cc,  &
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
                                       CALL ZMMCH(transa,transb,m,n,k,  &
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
!     End of ZCHK1.
!
      END SUBROUTINE ZCHK1
!*==zchk2.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
      SUBROUTINE ZCHK2(Sname,Eps,Thresh,Nout,Ntra,Trace,Rewi,Fatal,     &
     &                 Nidim,Idim,Nalf,Alf,Nbet,Bet,Nmax,A,Aa,As,B,Bb,  &
     &                 Bs,C,Cc,Cs,Ct,G)
      IMPLICIT NONE
!*--ZCHK2672
!
!  Tests ZHEMM and ZSYMM.
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
      COMPLEX*16 ZERO
      PARAMETER (ZERO=(0.0D0,0.0D0))
      DOUBLE PRECISION RZERO
      PARAMETER (RZERO=0.0D0)
!     .. Scalar Arguments ..
      DOUBLE PRECISION Eps , Thresh
      INTEGER Nalf , Nbet , Nidim , Nmax , Nout , Ntra
      LOGICAL Fatal , Rewi , Trace
      CHARACTER*6 Sname
!     .. Array Arguments ..
      COMPLEX*16 A(Nmax,Nmax) , Aa(Nmax*Nmax) , Alf(Nalf) ,             &
     &           As(Nmax*Nmax) , B(Nmax,Nmax) , Bb(Nmax*Nmax) ,         &
     &           Bet(Nbet) , Bs(Nmax*Nmax) , C(Nmax,Nmax) ,             &
     &           Cc(Nmax*Nmax) , Cs(Nmax*Nmax) , Ct(Nmax)
      DOUBLE PRECISION G(Nmax)
      INTEGER Idim(Nidim)
!     .. Local Scalars ..
      COMPLEX*16 alpha , als , beta , bls
      DOUBLE PRECISION err , errmax
      INTEGER i , ia , ib , ics , icu , im , in , laa , lbb , lcc ,     &
     &        lda , ldas , ldb , ldbs , ldc , ldcs , m , ms , n , na ,  &
     &        nargs , nc , ns
      LOGICAL conj , left , null , reset , same
      CHARACTER*1 side , sides , uplo , uplos
      CHARACTER*2 ichs , ichu
!     .. Local Arrays ..
      LOGICAL isame(13)
!     .. External Functions ..
      LOGICAL LZE , LZERES
      EXTERNAL LZE , LZERES
!     .. External Subroutines ..
      EXTERNAL ZHEMM , ZMAKE , ZMMCH , ZSYMM
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
                  CALL ZMAKE('GE',' ',' ',m,n,B,Nmax,Bb,ldb,reset,ZERO)
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
                           CALL ZMAKE(Sname(2:3),uplo,' ',na,na,A,Nmax, &
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
                                 CALL ZMAKE('GE',' ',' ',m,n,C,Nmax,Cc, &
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
                                    CALL ZHEMM(side,uplo,m,n,alpha,Aa,  &
     &                                 lda,Bb,ldb,beta,Cc,ldc)
                                 ELSE
                                    CALL ZSYMM(side,uplo,m,n,alpha,Aa,  &
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
                                 isame(6) = LZE(As,Aa,laa)
                                 isame(7) = ldas==lda
                                 isame(8) = LZE(Bs,Bb,lbb)
                                 isame(9) = ldbs==ldb
                                 isame(10) = bls==beta
                                 IF ( null ) THEN
                                    isame(11) = LZE(Cs,Cc,lcc)
                                 ELSE
                                    isame(11)                           &
     &                                 = LZERES('GE',' ',m,n,Cs,Cc,ldc)
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
                                       CALL ZMMCH('N','N',m,n,m,alpha,A,&
     &                                    Nmax,B,Nmax,beta,C,Nmax,Ct,G, &
     &                                    Cc,ldc,Eps,err,Fatal,Nout,    &
     &                                    .TRUE.)
                                    ELSE
                                       CALL ZMMCH('N','N',m,n,n,alpha,B,&
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
!     End of ZCHK2.
!
      END SUBROUTINE ZCHK2
!*==zchk3.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
      SUBROUTINE ZCHK3(Sname,Eps,Thresh,Nout,Ntra,Trace,Rewi,Fatal,     &
     &                 Nidim,Idim,Nalf,Alf,Nmax,A,Aa,As,B,Bb,Bs,Ct,G,C)
      IMPLICIT NONE
!*--ZCHK3944
!
!  Tests ZTRMM and ZTRSM.
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
      COMPLEX*16 ZERO , ONE
      PARAMETER (ZERO=(0.0D0,0.0D0),ONE=(1.0D0,0.0D0))
      DOUBLE PRECISION RZERO
      PARAMETER (RZERO=0.0D0)
!     .. Scalar Arguments ..
      DOUBLE PRECISION Eps , Thresh
      INTEGER Nalf , Nidim , Nmax , Nout , Ntra
      LOGICAL Fatal , Rewi , Trace
      CHARACTER*6 Sname
!     .. Array Arguments ..
      COMPLEX*16 A(Nmax,Nmax) , Aa(Nmax*Nmax) , Alf(Nalf) ,             &
     &           As(Nmax*Nmax) , B(Nmax,Nmax) , Bb(Nmax*Nmax) ,         &
     &           Bs(Nmax*Nmax) , C(Nmax,Nmax) , Ct(Nmax)
      DOUBLE PRECISION G(Nmax)
      INTEGER Idim(Nidim)
!     .. Local Scalars ..
      COMPLEX*16 alpha , als
      DOUBLE PRECISION err , errmax
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
      LOGICAL LZE , LZERES
      EXTERNAL LZE , LZERES
!     .. External Subroutines ..
      EXTERNAL ZMAKE , ZMMCH , ZTRMM , ZTRSM
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
!     Set up zero matrix for ZMMCH.
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
                              CALL ZMAKE('TR',uplo,diag,na,na,A,Nmax,Aa,&
     &                           lda,reset,ZERO)
!
!                          Generate the matrix B.
!
                              CALL ZMAKE('GE',' ',' ',m,n,B,Nmax,Bb,ldb,&
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
                                 CALL ZTRMM(side,uplo,transa,diag,m,n,  &
     &                              alpha,Aa,lda,Bb,ldb)
                              ELSEIF ( Sname(4:5)=='SM' ) THEN
                                 IF ( Trace ) WRITE (Ntra,FMT=99005)    &
     &                                nc , Sname , side , uplo ,        &
     &                                transa , diag , m , n , alpha ,   &
     &                                lda , ldb
                                 IF ( Rewi ) REWIND Ntra
                                 CALL ZTRSM(side,uplo,transa,diag,m,n,  &
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
                              isame(8) = LZE(As,Aa,laa)
                              isame(9) = ldas==lda
                              IF ( null ) THEN
                                 isame(10) = LZE(Bs,Bb,lbb)
                              ELSE
                                 isame(10)                              &
     &                              = LZERES('GE',' ',m,n,Bs,Bb,ldb)
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
                                       CALL ZMMCH(transa,'N',m,n,m,     &
     &                                    alpha,A,Nmax,B,Nmax,ZERO,C,   &
     &                                    Nmax,Ct,G,Bb,ldb,Eps,err,     &
     &                                    Fatal,Nout,.TRUE.)
                                    ELSE
                                       CALL ZMMCH('N',transa,m,n,n,     &
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
                                       CALL ZMMCH(transa,'N',m,n,m,ONE, &
     &                                    A,Nmax,C,Nmax,ZERO,B,Nmax,Ct, &
     &                                    G,Bb,ldb,Eps,err,Fatal,Nout,  &
     &                                    .FALSE.)
                                    ELSE
                                       CALL ZMMCH('N',transa,m,n,n,ONE, &
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
!     End of ZCHK3.
!
      END SUBROUTINE ZCHK3
!*==zchk4.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
      SUBROUTINE ZCHK4(Sname,Eps,Thresh,Nout,Ntra,Trace,Rewi,Fatal,     &
     &                 Nidim,Idim,Nalf,Alf,Nbet,Bet,Nmax,A,Aa,As,B,Bb,  &
     &                 Bs,C,Cc,Cs,Ct,G)
      IMPLICIT NONE
!*--ZCHK41241
!
!  Tests ZHERK and ZSYRK.
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
      COMPLEX*16 ZERO
      PARAMETER (ZERO=(0.0D0,0.0D0))
      DOUBLE PRECISION RONE , RZERO
      PARAMETER (RONE=1.0D0,RZERO=0.0D0)
!     .. Scalar Arguments ..
      DOUBLE PRECISION Eps , Thresh
      INTEGER Nalf , Nbet , Nidim , Nmax , Nout , Ntra
      LOGICAL Fatal , Rewi , Trace
      CHARACTER*6 Sname
!     .. Array Arguments ..
      COMPLEX*16 A(Nmax,Nmax) , Aa(Nmax*Nmax) , Alf(Nalf) ,             &
     &           As(Nmax*Nmax) , B(Nmax,Nmax) , Bb(Nmax*Nmax) ,         &
     &           Bet(Nbet) , Bs(Nmax*Nmax) , C(Nmax,Nmax) ,             &
     &           Cc(Nmax*Nmax) , Cs(Nmax*Nmax) , Ct(Nmax)
      DOUBLE PRECISION G(Nmax)
      INTEGER Idim(Nidim)
!     .. Local Scalars ..
      COMPLEX*16 alpha , als , beta , bets
      DOUBLE PRECISION err , errmax , ralpha , rals , rbeta , rbets
      INTEGER i , ia , ib , ict , icu , ik , in , j , jc , jj , k , ks ,&
     &        laa , lcc , lda , ldas , ldc , ldcs , lj , ma , n , na ,  &
     &        nargs , nc , ns
      LOGICAL conj , null , reset , same , tran , upper
      CHARACTER*1 trans , transs , transt , uplo , uplos
      CHARACTER*2 icht , ichu
!     .. Local Arrays ..
      LOGICAL isame(13)
!     .. External Functions ..
      LOGICAL LZE , LZERES
      EXTERNAL LZE , LZERES
!     .. External Subroutines ..
      EXTERNAL ZHERK , ZMAKE , ZMMCH , ZSYRK
!     .. Intrinsic Functions ..
      INTRINSIC DCMPLX , MAX , DBLE
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
                     CALL ZMAKE('GE',' ',' ',ma,na,A,Nmax,Aa,lda,reset, &
     &                          ZERO)
!
                     DO icu = 1 , 2
                        uplo = ichu(icu:icu)
                        upper = uplo=='U'
!
                        DO ia = 1 , Nalf
                           alpha = Alf(ia)
                           IF ( conj ) THEN
                              ralpha = DBLE(alpha)
                              alpha = DCMPLX(ralpha,RZERO)
                           ENDIF
!
                           DO ib = 1 , Nbet
                              beta = Bet(ib)
                              IF ( conj ) THEN
                                 rbeta = DBLE(beta)
                                 beta = DCMPLX(rbeta,RZERO)
                              ENDIF
                              null = n<=0
                              IF ( conj ) null = null .OR.              &
     &                             ((k<=0 .OR. ralpha==RZERO) .AND.     &
     &                             rbeta==RONE)
!
!                       Generate the matrix C.
!
                              CALL ZMAKE(Sname(2:3),uplo,' ',n,n,C,Nmax,&
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
                                 CALL ZHERK(uplo,trans,n,k,ralpha,Aa,   &
     &                              lda,rbeta,Cc,ldc)
                              ELSE
                                 IF ( Trace ) WRITE (Ntra,FMT=99007)    &
     &                                nc , Sname , uplo , trans , n ,   &
     &                                k , alpha , lda , beta , ldc
                                 IF ( Rewi ) REWIND Ntra
                                 CALL ZSYRK(uplo,trans,n,k,alpha,Aa,lda,&
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
                              isame(6) = LZE(As,Aa,laa)
                              isame(7) = ldas==lda
                              IF ( conj ) THEN
                                 isame(8) = rbets==rbeta
                              ELSE
                                 isame(8) = bets==beta
                              ENDIF
                              IF ( null ) THEN
                                 isame(9) = LZE(Cs,Cc,lcc)
                              ELSE
                                 isame(9)                               &
     &                              = LZERES(Sname(2:3),uplo,n,n,Cs,Cc, &
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
                                       CALL ZMMCH(transt,'N',lj,1,k,    &
     &                                    alpha,A(1,jj),Nmax,A(1,j),    &
     &                                    Nmax,beta,C(jj,j),Nmax,Ct,G,  &
     &                                    Cc(jc),ldc,Eps,err,Fatal,Nout,&
     &                                    .TRUE.)
                                    ELSE
                                       CALL ZMMCH('N',transt,lj,1,k,    &
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
!     End of ZCHK4.
!
      END SUBROUTINE ZCHK4
!*==zchk5.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
      SUBROUTINE ZCHK5(Sname,Eps,Thresh,Nout,Ntra,Trace,Rewi,Fatal,     &
     &                 Nidim,Idim,Nalf,Alf,Nbet,Bet,Nmax,Ab,Aa,As,Bb,Bs,&
     &                 C,Cc,Cs,Ct,G,W)
      IMPLICIT NONE
!*--ZCHK51563
!
!  Tests ZHER2K and ZSYR2K.
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
      COMPLEX*16 ZERO , ONE
      PARAMETER (ZERO=(0.0D0,0.0D0),ONE=(1.0D0,0.0D0))
      DOUBLE PRECISION RONE , RZERO
      PARAMETER (RONE=1.0D0,RZERO=0.0D0)
!     .. Scalar Arguments ..
      DOUBLE PRECISION Eps , Thresh
      INTEGER Nalf , Nbet , Nidim , Nmax , Nout , Ntra
      LOGICAL Fatal , Rewi , Trace
      CHARACTER*6 Sname
!     .. Array Arguments ..
      COMPLEX*16 Aa(Nmax*Nmax) , Ab(2*Nmax*Nmax) , Alf(Nalf) ,          &
     &           As(Nmax*Nmax) , Bb(Nmax*Nmax) , Bet(Nbet) ,            &
     &           Bs(Nmax*Nmax) , C(Nmax,Nmax) , Cc(Nmax*Nmax) ,         &
     &           Cs(Nmax*Nmax) , Ct(Nmax) , W(2*Nmax)
      DOUBLE PRECISION G(Nmax)
      INTEGER Idim(Nidim)
!     .. Local Scalars ..
      COMPLEX*16 alpha , als , beta , bets
      DOUBLE PRECISION err , errmax , rbeta , rbets
      INTEGER i , ia , ib , ict , icu , ik , in , j , jc , jj , jjab ,  &
     &        k , ks , laa , lbb , lcc , lda , ldas , ldb , ldbs , ldc ,&
     &        ldcs , lj , ma , n , na , nargs , nc , ns
      LOGICAL conj , null , reset , same , tran , upper
      CHARACTER*1 trans , transs , transt , uplo , uplos
      CHARACTER*2 icht , ichu
!     .. Local Arrays ..
      LOGICAL isame(13)
!     .. External Functions ..
      LOGICAL LZE , LZERES
      EXTERNAL LZE , LZERES
!     .. External Subroutines ..
      EXTERNAL ZHER2K , ZMAKE , ZMMCH , ZSYR2K
!     .. Intrinsic Functions ..
      INTRINSIC DCMPLX , DCONJG , MAX , DBLE
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
                        CALL ZMAKE('GE',' ',' ',ma,na,Ab,2*Nmax,Aa,lda, &
     &                             reset,ZERO)
                     ELSE
                        CALL ZMAKE('GE',' ',' ',ma,na,Ab,Nmax,Aa,lda,   &
     &                             reset,ZERO)
                     ENDIF
!
!              Generate the matrix B.
!
                     ldb = lda
                     lbb = laa
                     IF ( tran ) THEN
                        CALL ZMAKE('GE',' ',' ',ma,na,Ab(k+1),2*Nmax,Bb,&
     &                             ldb,reset,ZERO)
                     ELSE
                        CALL ZMAKE('GE',' ',' ',ma,na,Ab(k*Nmax+1),Nmax,&
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
                                 rbeta = DBLE(beta)
                                 beta = DCMPLX(rbeta,RZERO)
                              ENDIF
                              null = n<=0
                              IF ( conj ) null = null .OR.              &
     &                             ((k<=0 .OR. alpha==ZERO) .AND.       &
     &                             rbeta==RONE)
!
!                       Generate the matrix C.
!
                              CALL ZMAKE(Sname(2:3),uplo,' ',n,n,C,Nmax,&
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
                                 CALL ZHER2K(uplo,trans,n,k,alpha,Aa,   &
     &                              lda,Bb,ldb,rbeta,Cc,ldc)
                              ELSE
                                 IF ( Trace ) WRITE (Ntra,FMT=99007)    &
     &                                nc , Sname , uplo , trans , n ,   &
     &                                k , alpha , lda , ldb , beta , ldc
                                 IF ( Rewi ) REWIND Ntra
                                 CALL ZSYR2K(uplo,trans,n,k,alpha,Aa,   &
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
                              isame(6) = LZE(As,Aa,laa)
                              isame(7) = ldas==lda
                              isame(8) = LZE(Bs,Bb,lbb)
                              isame(9) = ldbs==ldb
                              IF ( conj ) THEN
                                 isame(10) = rbets==rbeta
                              ELSE
                                 isame(10) = bets==beta
                              ENDIF
                              IF ( null ) THEN
                                 isame(11) = LZE(Cs,Cc,lcc)
                              ELSE
                                 isame(11)                              &
     &                              = LZERES('HE',uplo,n,n,Cs,Cc,ldc)
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
                                         W(k+i) = DCONJG(alpha)         &
     &                                      *Ab((j-1)*2*Nmax+i)
                                         ELSE
                                         W(k+i)                         &
     &                                      = alpha*Ab((j-1)*2*Nmax+i)
                                         ENDIF
                                       ENDDO
                                       CALL ZMMCH(transt,'N',lj,1,2*k,  &
     &                                    ONE,Ab(jjab),2*Nmax,W,2*Nmax, &
     &                                    beta,C(jj,j),Nmax,Ct,G,Cc(jc),&
     &                                    ldc,Eps,err,Fatal,Nout,.TRUE.)
                                    ELSE
                                       DO i = 1 , k
                                         IF ( conj ) THEN
                                         W(i) = alpha*DCONJG(Ab((k+i-1)*&
     &                                      Nmax+j))
                                         W(k+i) = DCONJG(alpha*Ab((i-1)*&
     &                                      Nmax+j))
                                         ELSE
                                         W(i) = alpha*Ab((k+i-1)*Nmax+j)
                                         W(k+i) = alpha*Ab((i-1)*Nmax+j)
                                         ENDIF
                                       ENDDO
                                       CALL ZMMCH('N','N',lj,1,2*k,ONE, &
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
!     End of ZCHK5.
!
      END SUBROUTINE ZCHK5
!*==zchke.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
      SUBROUTINE ZCHKE(Isnum,Srnamt,Nout)
      IMPLICIT NONE
!*--ZCHKE1917
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
!  3-19-92:  Fix argument 12 in calls to ZSYMM and ZHEMM
!            with INFOT = 9  (eca)
!  10-9-00:  Declared INTRINSIC DCMPLX (susan)
!
!     .. Scalar Arguments ..
      INTEGER Isnum , Nout
      CHARACTER*6 Srnamt
!     .. Scalars in Common ..
      INTEGER INFot , NOUtc
      LOGICAL LERr , OK
!     .. Parameters ..
      REAL ONE , TWO
      PARAMETER (ONE=1.0D0,TWO=2.0D0)
!     .. Local Scalars ..
      COMPLEX*16 alpha , beta
      DOUBLE PRECISION ralpha , rbeta
!     .. Local Arrays ..
      COMPLEX*16 a(2,1) , b(2,1) , c(2,1)
!     .. External Subroutines ..
      EXTERNAL ZGEMM , ZHEMM , ZHER2K , ZHERK , CHKXER , ZSYMM ,        &
     &         ZSYR2K , ZSYRK , ZTRMM , ZTRSM
!     .. Intrinsic Functions ..
      INTRINSIC DCMPLX
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
      alpha = DCMPLX(ONE,-ONE)
      beta = DCMPLX(TWO,-TWO)
      ralpha = ONE
      rbeta = TWO
!
      IF ( Isnum==2 ) THEN
         INFot = 1
         CALL ZHEMM('/','U',0,0,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 2
         CALL ZHEMM('L','/',0,0,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 3
         CALL ZHEMM('L','U',-1,0,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 3
         CALL ZHEMM('R','U',-1,0,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 3
         CALL ZHEMM('L','L',-1,0,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 3
         CALL ZHEMM('R','L',-1,0,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 4
         CALL ZHEMM('L','U',0,-1,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 4
         CALL ZHEMM('R','U',0,-1,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 4
         CALL ZHEMM('L','L',0,-1,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 4
         CALL ZHEMM('R','L',0,-1,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 7
         CALL ZHEMM('L','U',2,0,alpha,a,1,b,2,beta,c,2)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 7
         CALL ZHEMM('R','U',0,2,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 7
         CALL ZHEMM('L','L',2,0,alpha,a,1,b,2,beta,c,2)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 7
         CALL ZHEMM('R','L',0,2,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 9
         CALL ZHEMM('L','U',2,0,alpha,a,2,b,1,beta,c,2)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 9
         CALL ZHEMM('R','U',2,0,alpha,a,1,b,1,beta,c,2)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 9
         CALL ZHEMM('L','L',2,0,alpha,a,2,b,1,beta,c,2)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 9
         CALL ZHEMM('R','L',2,0,alpha,a,1,b,1,beta,c,2)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 12
         CALL ZHEMM('L','U',2,0,alpha,a,2,b,2,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 12
         CALL ZHEMM('R','U',2,0,alpha,a,1,b,2,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 12
         CALL ZHEMM('L','L',2,0,alpha,a,2,b,2,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 12
         CALL ZHEMM('R','L',2,0,alpha,a,1,b,2,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
      ELSEIF ( Isnum==3 ) THEN
         INFot = 1
         CALL ZSYMM('/','U',0,0,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 2
         CALL ZSYMM('L','/',0,0,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 3
         CALL ZSYMM('L','U',-1,0,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 3
         CALL ZSYMM('R','U',-1,0,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 3
         CALL ZSYMM('L','L',-1,0,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 3
         CALL ZSYMM('R','L',-1,0,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 4
         CALL ZSYMM('L','U',0,-1,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 4
         CALL ZSYMM('R','U',0,-1,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 4
         CALL ZSYMM('L','L',0,-1,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 4
         CALL ZSYMM('R','L',0,-1,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 7
         CALL ZSYMM('L','U',2,0,alpha,a,1,b,2,beta,c,2)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 7
         CALL ZSYMM('R','U',0,2,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 7
         CALL ZSYMM('L','L',2,0,alpha,a,1,b,2,beta,c,2)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 7
         CALL ZSYMM('R','L',0,2,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 9
         CALL ZSYMM('L','U',2,0,alpha,a,2,b,1,beta,c,2)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 9
         CALL ZSYMM('R','U',2,0,alpha,a,1,b,1,beta,c,2)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 9
         CALL ZSYMM('L','L',2,0,alpha,a,2,b,1,beta,c,2)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 9
         CALL ZSYMM('R','L',2,0,alpha,a,1,b,1,beta,c,2)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 12
         CALL ZSYMM('L','U',2,0,alpha,a,2,b,2,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 12
         CALL ZSYMM('R','U',2,0,alpha,a,1,b,2,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 12
         CALL ZSYMM('L','L',2,0,alpha,a,2,b,2,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 12
         CALL ZSYMM('R','L',2,0,alpha,a,1,b,2,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
      ELSEIF ( Isnum==4 ) THEN
         INFot = 1
         CALL ZTRMM('/','U','N','N',0,0,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 2
         CALL ZTRMM('L','/','N','N',0,0,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 3
         CALL ZTRMM('L','U','/','N',0,0,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 4
         CALL ZTRMM('L','U','N','/',0,0,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 5
         CALL ZTRMM('L','U','N','N',-1,0,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 5
         CALL ZTRMM('L','U','C','N',-1,0,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 5
         CALL ZTRMM('L','U','T','N',-1,0,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 5
         CALL ZTRMM('R','U','N','N',-1,0,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 5
         CALL ZTRMM('R','U','C','N',-1,0,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 5
         CALL ZTRMM('R','U','T','N',-1,0,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 5
         CALL ZTRMM('L','L','N','N',-1,0,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 5
         CALL ZTRMM('L','L','C','N',-1,0,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 5
         CALL ZTRMM('L','L','T','N',-1,0,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 5
         CALL ZTRMM('R','L','N','N',-1,0,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 5
         CALL ZTRMM('R','L','C','N',-1,0,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 5
         CALL ZTRMM('R','L','T','N',-1,0,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 6
         CALL ZTRMM('L','U','N','N',0,-1,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 6
         CALL ZTRMM('L','U','C','N',0,-1,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 6
         CALL ZTRMM('L','U','T','N',0,-1,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 6
         CALL ZTRMM('R','U','N','N',0,-1,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 6
         CALL ZTRMM('R','U','C','N',0,-1,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 6
         CALL ZTRMM('R','U','T','N',0,-1,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 6
         CALL ZTRMM('L','L','N','N',0,-1,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 6
         CALL ZTRMM('L','L','C','N',0,-1,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 6
         CALL ZTRMM('L','L','T','N',0,-1,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 6
         CALL ZTRMM('R','L','N','N',0,-1,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 6
         CALL ZTRMM('R','L','C','N',0,-1,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 6
         CALL ZTRMM('R','L','T','N',0,-1,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 9
         CALL ZTRMM('L','U','N','N',2,0,alpha,a,1,b,2)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 9
         CALL ZTRMM('L','U','C','N',2,0,alpha,a,1,b,2)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 9
         CALL ZTRMM('L','U','T','N',2,0,alpha,a,1,b,2)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 9
         CALL ZTRMM('R','U','N','N',0,2,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 9
         CALL ZTRMM('R','U','C','N',0,2,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 9
         CALL ZTRMM('R','U','T','N',0,2,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 9
         CALL ZTRMM('L','L','N','N',2,0,alpha,a,1,b,2)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 9
         CALL ZTRMM('L','L','C','N',2,0,alpha,a,1,b,2)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 9
         CALL ZTRMM('L','L','T','N',2,0,alpha,a,1,b,2)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 9
         CALL ZTRMM('R','L','N','N',0,2,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 9
         CALL ZTRMM('R','L','C','N',0,2,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 9
         CALL ZTRMM('R','L','T','N',0,2,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 11
         CALL ZTRMM('L','U','N','N',2,0,alpha,a,2,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 11
         CALL ZTRMM('L','U','C','N',2,0,alpha,a,2,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 11
         CALL ZTRMM('L','U','T','N',2,0,alpha,a,2,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 11
         CALL ZTRMM('R','U','N','N',2,0,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 11
         CALL ZTRMM('R','U','C','N',2,0,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 11
         CALL ZTRMM('R','U','T','N',2,0,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 11
         CALL ZTRMM('L','L','N','N',2,0,alpha,a,2,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 11
         CALL ZTRMM('L','L','C','N',2,0,alpha,a,2,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 11
         CALL ZTRMM('L','L','T','N',2,0,alpha,a,2,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 11
         CALL ZTRMM('R','L','N','N',2,0,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 11
         CALL ZTRMM('R','L','C','N',2,0,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 11
         CALL ZTRMM('R','L','T','N',2,0,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
      ELSEIF ( Isnum==5 ) THEN
         INFot = 1
         CALL ZTRSM('/','U','N','N',0,0,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 2
         CALL ZTRSM('L','/','N','N',0,0,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 3
         CALL ZTRSM('L','U','/','N',0,0,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 4
         CALL ZTRSM('L','U','N','/',0,0,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 5
         CALL ZTRSM('L','U','N','N',-1,0,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 5
         CALL ZTRSM('L','U','C','N',-1,0,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 5
         CALL ZTRSM('L','U','T','N',-1,0,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 5
         CALL ZTRSM('R','U','N','N',-1,0,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 5
         CALL ZTRSM('R','U','C','N',-1,0,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 5
         CALL ZTRSM('R','U','T','N',-1,0,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 5
         CALL ZTRSM('L','L','N','N',-1,0,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 5
         CALL ZTRSM('L','L','C','N',-1,0,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 5
         CALL ZTRSM('L','L','T','N',-1,0,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 5
         CALL ZTRSM('R','L','N','N',-1,0,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 5
         CALL ZTRSM('R','L','C','N',-1,0,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 5
         CALL ZTRSM('R','L','T','N',-1,0,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 6
         CALL ZTRSM('L','U','N','N',0,-1,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 6
         CALL ZTRSM('L','U','C','N',0,-1,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 6
         CALL ZTRSM('L','U','T','N',0,-1,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 6
         CALL ZTRSM('R','U','N','N',0,-1,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 6
         CALL ZTRSM('R','U','C','N',0,-1,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 6
         CALL ZTRSM('R','U','T','N',0,-1,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 6
         CALL ZTRSM('L','L','N','N',0,-1,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 6
         CALL ZTRSM('L','L','C','N',0,-1,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 6
         CALL ZTRSM('L','L','T','N',0,-1,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 6
         CALL ZTRSM('R','L','N','N',0,-1,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 6
         CALL ZTRSM('R','L','C','N',0,-1,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 6
         CALL ZTRSM('R','L','T','N',0,-1,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 9
         CALL ZTRSM('L','U','N','N',2,0,alpha,a,1,b,2)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 9
         CALL ZTRSM('L','U','C','N',2,0,alpha,a,1,b,2)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 9
         CALL ZTRSM('L','U','T','N',2,0,alpha,a,1,b,2)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 9
         CALL ZTRSM('R','U','N','N',0,2,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 9
         CALL ZTRSM('R','U','C','N',0,2,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 9
         CALL ZTRSM('R','U','T','N',0,2,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 9
         CALL ZTRSM('L','L','N','N',2,0,alpha,a,1,b,2)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 9
         CALL ZTRSM('L','L','C','N',2,0,alpha,a,1,b,2)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 9
         CALL ZTRSM('L','L','T','N',2,0,alpha,a,1,b,2)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 9
         CALL ZTRSM('R','L','N','N',0,2,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 9
         CALL ZTRSM('R','L','C','N',0,2,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 9
         CALL ZTRSM('R','L','T','N',0,2,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 11
         CALL ZTRSM('L','U','N','N',2,0,alpha,a,2,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 11
         CALL ZTRSM('L','U','C','N',2,0,alpha,a,2,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 11
         CALL ZTRSM('L','U','T','N',2,0,alpha,a,2,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 11
         CALL ZTRSM('R','U','N','N',2,0,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 11
         CALL ZTRSM('R','U','C','N',2,0,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 11
         CALL ZTRSM('R','U','T','N',2,0,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 11
         CALL ZTRSM('L','L','N','N',2,0,alpha,a,2,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 11
         CALL ZTRSM('L','L','C','N',2,0,alpha,a,2,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 11
         CALL ZTRSM('L','L','T','N',2,0,alpha,a,2,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 11
         CALL ZTRSM('R','L','N','N',2,0,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 11
         CALL ZTRSM('R','L','C','N',2,0,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 11
         CALL ZTRSM('R','L','T','N',2,0,alpha,a,1,b,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
      ELSEIF ( Isnum==6 ) THEN
         INFot = 1
         CALL ZHERK('/','N',0,0,ralpha,a,1,rbeta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 2
         CALL ZHERK('U','T',0,0,ralpha,a,1,rbeta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 3
         CALL ZHERK('U','N',-1,0,ralpha,a,1,rbeta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 3
         CALL ZHERK('U','C',-1,0,ralpha,a,1,rbeta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 3
         CALL ZHERK('L','N',-1,0,ralpha,a,1,rbeta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 3
         CALL ZHERK('L','C',-1,0,ralpha,a,1,rbeta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 4
         CALL ZHERK('U','N',0,-1,ralpha,a,1,rbeta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 4
         CALL ZHERK('U','C',0,-1,ralpha,a,1,rbeta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 4
         CALL ZHERK('L','N',0,-1,ralpha,a,1,rbeta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 4
         CALL ZHERK('L','C',0,-1,ralpha,a,1,rbeta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 7
         CALL ZHERK('U','N',2,0,ralpha,a,1,rbeta,c,2)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 7
         CALL ZHERK('U','C',0,2,ralpha,a,1,rbeta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 7
         CALL ZHERK('L','N',2,0,ralpha,a,1,rbeta,c,2)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 7
         CALL ZHERK('L','C',0,2,ralpha,a,1,rbeta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 10
         CALL ZHERK('U','N',2,0,ralpha,a,2,rbeta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 10
         CALL ZHERK('U','C',2,0,ralpha,a,1,rbeta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 10
         CALL ZHERK('L','N',2,0,ralpha,a,2,rbeta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 10
         CALL ZHERK('L','C',2,0,ralpha,a,1,rbeta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
      ELSEIF ( Isnum==7 ) THEN
         INFot = 1
         CALL ZSYRK('/','N',0,0,alpha,a,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 2
         CALL ZSYRK('U','C',0,0,alpha,a,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 3
         CALL ZSYRK('U','N',-1,0,alpha,a,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 3
         CALL ZSYRK('U','T',-1,0,alpha,a,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 3
         CALL ZSYRK('L','N',-1,0,alpha,a,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 3
         CALL ZSYRK('L','T',-1,0,alpha,a,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 4
         CALL ZSYRK('U','N',0,-1,alpha,a,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 4
         CALL ZSYRK('U','T',0,-1,alpha,a,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 4
         CALL ZSYRK('L','N',0,-1,alpha,a,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 4
         CALL ZSYRK('L','T',0,-1,alpha,a,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 7
         CALL ZSYRK('U','N',2,0,alpha,a,1,beta,c,2)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 7
         CALL ZSYRK('U','T',0,2,alpha,a,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 7
         CALL ZSYRK('L','N',2,0,alpha,a,1,beta,c,2)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 7
         CALL ZSYRK('L','T',0,2,alpha,a,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 10
         CALL ZSYRK('U','N',2,0,alpha,a,2,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 10
         CALL ZSYRK('U','T',2,0,alpha,a,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 10
         CALL ZSYRK('L','N',2,0,alpha,a,2,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 10
         CALL ZSYRK('L','T',2,0,alpha,a,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
      ELSEIF ( Isnum==8 ) THEN
         INFot = 1
         CALL ZHER2K('/','N',0,0,alpha,a,1,b,1,rbeta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 2
         CALL ZHER2K('U','T',0,0,alpha,a,1,b,1,rbeta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 3
         CALL ZHER2K('U','N',-1,0,alpha,a,1,b,1,rbeta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 3
         CALL ZHER2K('U','C',-1,0,alpha,a,1,b,1,rbeta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 3
         CALL ZHER2K('L','N',-1,0,alpha,a,1,b,1,rbeta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 3
         CALL ZHER2K('L','C',-1,0,alpha,a,1,b,1,rbeta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 4
         CALL ZHER2K('U','N',0,-1,alpha,a,1,b,1,rbeta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 4
         CALL ZHER2K('U','C',0,-1,alpha,a,1,b,1,rbeta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 4
         CALL ZHER2K('L','N',0,-1,alpha,a,1,b,1,rbeta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 4
         CALL ZHER2K('L','C',0,-1,alpha,a,1,b,1,rbeta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 7
         CALL ZHER2K('U','N',2,0,alpha,a,1,b,1,rbeta,c,2)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 7
         CALL ZHER2K('U','C',0,2,alpha,a,1,b,1,rbeta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 7
         CALL ZHER2K('L','N',2,0,alpha,a,1,b,1,rbeta,c,2)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 7
         CALL ZHER2K('L','C',0,2,alpha,a,1,b,1,rbeta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 9
         CALL ZHER2K('U','N',2,0,alpha,a,2,b,1,rbeta,c,2)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 9
         CALL ZHER2K('U','C',0,2,alpha,a,2,b,1,rbeta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 9
         CALL ZHER2K('L','N',2,0,alpha,a,2,b,1,rbeta,c,2)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 9
         CALL ZHER2K('L','C',0,2,alpha,a,2,b,1,rbeta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 12
         CALL ZHER2K('U','N',2,0,alpha,a,2,b,2,rbeta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 12
         CALL ZHER2K('U','C',2,0,alpha,a,1,b,1,rbeta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 12
         CALL ZHER2K('L','N',2,0,alpha,a,2,b,2,rbeta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 12
         CALL ZHER2K('L','C',2,0,alpha,a,1,b,1,rbeta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
      ELSEIF ( Isnum==9 ) THEN
         INFot = 1
         CALL ZSYR2K('/','N',0,0,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 2
         CALL ZSYR2K('U','C',0,0,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 3
         CALL ZSYR2K('U','N',-1,0,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 3
         CALL ZSYR2K('U','T',-1,0,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 3
         CALL ZSYR2K('L','N',-1,0,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 3
         CALL ZSYR2K('L','T',-1,0,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 4
         CALL ZSYR2K('U','N',0,-1,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 4
         CALL ZSYR2K('U','T',0,-1,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 4
         CALL ZSYR2K('L','N',0,-1,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 4
         CALL ZSYR2K('L','T',0,-1,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 7
         CALL ZSYR2K('U','N',2,0,alpha,a,1,b,1,beta,c,2)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 7
         CALL ZSYR2K('U','T',0,2,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 7
         CALL ZSYR2K('L','N',2,0,alpha,a,1,b,1,beta,c,2)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 7
         CALL ZSYR2K('L','T',0,2,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 9
         CALL ZSYR2K('U','N',2,0,alpha,a,2,b,1,beta,c,2)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 9
         CALL ZSYR2K('U','T',0,2,alpha,a,2,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 9
         CALL ZSYR2K('L','N',2,0,alpha,a,2,b,1,beta,c,2)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 9
         CALL ZSYR2K('L','T',0,2,alpha,a,2,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 12
         CALL ZSYR2K('U','N',2,0,alpha,a,2,b,2,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 12
         CALL ZSYR2K('U','T',2,0,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 12
         CALL ZSYR2K('L','N',2,0,alpha,a,2,b,2,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 12
         CALL ZSYR2K('L','T',2,0,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
      ELSE
         INFot = 1
         CALL ZGEMM('/','N',0,0,0,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 1
         CALL ZGEMM('/','C',0,0,0,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 1
         CALL ZGEMM('/','T',0,0,0,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 2
         CALL ZGEMM('N','/',0,0,0,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 2
         CALL ZGEMM('C','/',0,0,0,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 2
         CALL ZGEMM('T','/',0,0,0,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 3
         CALL ZGEMM('N','N',-1,0,0,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 3
         CALL ZGEMM('N','C',-1,0,0,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 3
         CALL ZGEMM('N','T',-1,0,0,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 3
         CALL ZGEMM('C','N',-1,0,0,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 3
         CALL ZGEMM('C','C',-1,0,0,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 3
         CALL ZGEMM('C','T',-1,0,0,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 3
         CALL ZGEMM('T','N',-1,0,0,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 3
         CALL ZGEMM('T','C',-1,0,0,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 3
         CALL ZGEMM('T','T',-1,0,0,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 4
         CALL ZGEMM('N','N',0,-1,0,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 4
         CALL ZGEMM('N','C',0,-1,0,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 4
         CALL ZGEMM('N','T',0,-1,0,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 4
         CALL ZGEMM('C','N',0,-1,0,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 4
         CALL ZGEMM('C','C',0,-1,0,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 4
         CALL ZGEMM('C','T',0,-1,0,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 4
         CALL ZGEMM('T','N',0,-1,0,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 4
         CALL ZGEMM('T','C',0,-1,0,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 4
         CALL ZGEMM('T','T',0,-1,0,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 5
         CALL ZGEMM('N','N',0,0,-1,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 5
         CALL ZGEMM('N','C',0,0,-1,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 5
         CALL ZGEMM('N','T',0,0,-1,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 5
         CALL ZGEMM('C','N',0,0,-1,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 5
         CALL ZGEMM('C','C',0,0,-1,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 5
         CALL ZGEMM('C','T',0,0,-1,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 5
         CALL ZGEMM('T','N',0,0,-1,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 5
         CALL ZGEMM('T','C',0,0,-1,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 5
         CALL ZGEMM('T','T',0,0,-1,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 8
         CALL ZGEMM('N','N',2,0,0,alpha,a,1,b,1,beta,c,2)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 8
         CALL ZGEMM('N','C',2,0,0,alpha,a,1,b,1,beta,c,2)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 8
         CALL ZGEMM('N','T',2,0,0,alpha,a,1,b,1,beta,c,2)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 8
         CALL ZGEMM('C','N',0,0,2,alpha,a,1,b,2,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 8
         CALL ZGEMM('C','C',0,0,2,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 8
         CALL ZGEMM('C','T',0,0,2,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 8
         CALL ZGEMM('T','N',0,0,2,alpha,a,1,b,2,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 8
         CALL ZGEMM('T','C',0,0,2,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 8
         CALL ZGEMM('T','T',0,0,2,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 10
         CALL ZGEMM('N','N',0,0,2,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 10
         CALL ZGEMM('C','N',0,0,2,alpha,a,2,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 10
         CALL ZGEMM('T','N',0,0,2,alpha,a,2,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 10
         CALL ZGEMM('N','C',0,2,0,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 10
         CALL ZGEMM('C','C',0,2,0,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 10
         CALL ZGEMM('T','C',0,2,0,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 10
         CALL ZGEMM('N','T',0,2,0,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 10
         CALL ZGEMM('C','T',0,2,0,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 10
         CALL ZGEMM('T','T',0,2,0,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 13
         CALL ZGEMM('N','N',2,0,0,alpha,a,2,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 13
         CALL ZGEMM('N','C',2,0,0,alpha,a,2,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 13
         CALL ZGEMM('N','T',2,0,0,alpha,a,2,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 13
         CALL ZGEMM('C','N',2,0,0,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 13
         CALL ZGEMM('C','C',2,0,0,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 13
         CALL ZGEMM('C','T',2,0,0,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 13
         CALL ZGEMM('T','N',2,0,0,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 13
         CALL ZGEMM('T','C',2,0,0,alpha,a,1,b,1,beta,c,1)
         CALL CHKXER(Srnamt,INFot,Nout,LERr,OK)
         INFot = 13
         CALL ZGEMM('T','T',2,0,0,alpha,a,1,b,1,beta,c,1)
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
!     End of ZCHKE.
!
      END SUBROUTINE ZCHKE
!*==zmake.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
      SUBROUTINE ZMAKE(Type,Uplo,Diag,M,N,A,Nmax,Aa,Lda,Reset,Transl)
      IMPLICIT NONE
!*--ZMAKE2864
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
      COMPLEX*16 ZERO , ONE
      PARAMETER (ZERO=(0.0D0,0.0D0),ONE=(1.0D0,0.0D0))
      COMPLEX*16 ROGUE
      PARAMETER (ROGUE=(-1.0D10,1.0D10))
      DOUBLE PRECISION RZERO
      PARAMETER (RZERO=0.0D0)
      DOUBLE PRECISION RROGUE
      PARAMETER (RROGUE=-1.0D10)
!     .. Scalar Arguments ..
      COMPLEX*16 Transl
      INTEGER Lda , M , N , Nmax
      LOGICAL Reset
      CHARACTER*1 Diag , Uplo
      CHARACTER*2 Type
!     .. Array Arguments ..
      COMPLEX*16 A(Nmax,*) , Aa(*)
!     .. Local Scalars ..
      INTEGER i , ibeg , iend , j , jj
      LOGICAL gen , her , lower , sym , tri , unit , upper
!     .. External Functions ..
      COMPLEX*16 ZBEG
      EXTERNAL ZBEG
!     .. Intrinsic Functions ..
      INTRINSIC DCMPLX , DCONJG , DBLE
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
               A(i,j) = ZBEG(Reset) + Transl
               IF ( i/=j ) THEN
!                 Set some elements to zero
                  IF ( N>3 .AND. j==N/2 ) A(i,j) = ZERO
                  IF ( her ) THEN
                     A(j,i) = DCONJG(A(i,j))
                  ELSEIF ( sym ) THEN
                     A(j,i) = A(i,j)
                  ELSEIF ( tri ) THEN
                     A(j,i) = ZERO
                  ENDIF
               ENDIF
            ENDIF
         ENDDO
         IF ( her ) A(j,j) = DCMPLX(DBLE(A(j,j)),RZERO)
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
               Aa(jj) = DCMPLX(DBLE(Aa(jj)),RROGUE)
            ENDIF
         ENDDO
      ENDIF
!
!     End of ZMAKE.
!
      END SUBROUTINE ZMAKE
!*==zmmch.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
      SUBROUTINE ZMMCH(Transa,Transb,M,N,Kk,Alpha,A,Lda,B,Ldb,Beta,C,   &
     &                 Ldc,Ct,G,Cc,Ldcc,Eps,Err,Fatal,Nout,Mv)
      IMPLICIT NONE
!*--ZMMCH2990
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
      COMPLEX*16 ZERO
      PARAMETER (ZERO=(0.0D0,0.0D0))
      DOUBLE PRECISION RZERO , RONE
      PARAMETER (RZERO=0.0D0,RONE=1.0D0)
!     .. Scalar Arguments ..
      COMPLEX*16 Alpha , Beta
      DOUBLE PRECISION Eps , Err
      INTEGER Kk , Lda , Ldb , Ldc , Ldcc , M , N , Nout
      LOGICAL Fatal , Mv
      CHARACTER*1 Transa , Transb
!     .. Array Arguments ..
      COMPLEX*16 A(Lda,*) , B(Ldb,*) , C(Ldc,*) , Cc(Ldcc,*) , Ct(*)
      DOUBLE PRECISION G(*)
!     .. Local Scalars ..
      COMPLEX*16 cl
      DOUBLE PRECISION erri
      INTEGER i , j , k
      LOGICAL ctrana , ctranb , trana , tranb
!     .. Intrinsic Functions ..
      INTRINSIC ABS , DIMAG , DCONJG , MAX , DBLE , SQRT
!     .. Statement Functions ..
      DOUBLE PRECISION ABS1
!     .. Statement Function definitions ..
      ABS1(cl) = ABS(DBLE(cl)) + ABS(DIMAG(cl))
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
                     Ct(i) = Ct(i) + DCONJG(A(k,i))*B(k,j)
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
                     Ct(i) = Ct(i) + A(i,k)*DCONJG(B(j,k))
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
                        Ct(i) = Ct(i) + DCONJG(A(k,i))*DCONJG(B(j,k))
                        G(i) = G(i) + ABS1(A(k,i))*ABS1(B(j,k))
                     ENDDO
                  ENDDO
               ELSE
                  DO k = 1 , Kk
                     DO i = 1 , M
                        Ct(i) = Ct(i) + DCONJG(A(k,i))*B(j,k)
                        G(i) = G(i) + ABS1(A(k,i))*ABS1(B(j,k))
                     ENDDO
                  ENDDO
               ENDIF
            ELSEIF ( ctranb ) THEN
               DO k = 1 , Kk
                  DO i = 1 , M
                     Ct(i) = Ct(i) + A(k,i)*DCONJG(B(j,k))
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
!     End of ZMMCH.
!
      END SUBROUTINE ZMMCH
!*==lze.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
      LOGICAL FUNCTION LZE(Ri,Rj,Lr)
      IMPLICIT NONE
!*--LZE3163
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
      COMPLEX*16 Ri(*) , Rj(*)
!     .. Local Scalars ..
      INTEGER i
!     .. Executable Statements ..
      DO i = 1 , Lr
         IF ( Ri(i)/=Rj(i) ) GOTO 100
      ENDDO
      LZE = .TRUE.
      GOTO 99999
 100  LZE = .FALSE.
!
!     End of LZE.
!
99999 END FUNCTION LZE
!*==lzeres.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
      LOGICAL FUNCTION LZERES(Type,Uplo,M,N,Aa,As,Lda)
      IMPLICIT NONE
!*--LZERES3195
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
      COMPLEX*16 Aa(Lda,*) , As(Lda,*)
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
      LZERES = .TRUE.
      GOTO 99999
 100  LZERES = .FALSE.
!
!     End of LZERES.
!
99999 END FUNCTION LZERES
!*==zbeg.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
      COMPLEX*16 FUNCTION ZBEG(Reset)
      IMPLICIT NONE
!*--ZBEG3254
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
      INTRINSIC DCMPLX
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
         ZBEG = DCMPLX((i-500)/1001.0D0,(j-500)/1001.0D0)
         EXIT
      ENDDO
!
!     End of ZBEG.
!
      END FUNCTION ZBEG
!*==ddiff.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
      DOUBLE PRECISION FUNCTION DDIFF(X,Y)
      IMPLICIT NONE
!*--DDIFF3313
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
      DOUBLE PRECISION X , Y
!     .. Executable Statements ..
      DDIFF = X - Y
!
!     End of DDIFF.
!
      END FUNCTION DDIFF
!*==chkxer.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
      SUBROUTINE CHKXER(Srnamt,Infot,Nout,Lerr,Ok)
      IMPLICIT NONE
!*--CHKXER3334
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
!*--XERBLA3367
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
 
