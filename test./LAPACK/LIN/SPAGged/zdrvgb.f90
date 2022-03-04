!*==zdrvgb.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b ZDRVGB
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZDRVGB( DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, A, LA,
!                          AFB, LAFB, ASAV, B, BSAV, X, XACT, S, WORK,
!                          RWORK, IWORK, NOUT )
!
!       .. Scalar Arguments ..
!       LOGICAL            TSTERR
!       INTEGER            LA, LAFB, NN, NOUT, NRHS
!       DOUBLE PRECISION   THRESH
!       ..
!       .. Array Arguments ..
!       LOGICAL            DOTYPE( * )
!       INTEGER            IWORK( * ), NVAL( * )
!       DOUBLE PRECISION   RWORK( * ), S( * )
!       COMPLEX*16         A( * ), AFB( * ), ASAV( * ), B( * ), BSAV( * ),
!      $                   WORK( * ), X( * ), XACT( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZDRVGB tests the driver routines ZGBSV and -SVX.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] DOTYPE
!> \verbatim
!>          DOTYPE is LOGICAL array, dimension (NTYPES)
!>          The matrix types to be used for testing.  Matrices of type j
!>          (for 1 <= j <= NTYPES) are used for testing if DOTYPE(j) =
!>          .TRUE.; if DOTYPE(j) = .FALSE., then type j is not used.
!> \endverbatim
!>
!> \param[in] NN
!> \verbatim
!>          NN is INTEGER
!>          The number of values of N contained in the vector NVAL.
!> \endverbatim
!>
!> \param[in] NVAL
!> \verbatim
!>          NVAL is INTEGER array, dimension (NN)
!>          The values of the matrix column dimension N.
!> \endverbatim
!>
!> \param[in] NRHS
!> \verbatim
!>          NRHS is INTEGER
!>          The number of right hand side vectors to be generated for
!>          each linear system.
!> \endverbatim
!>
!> \param[in] THRESH
!> \verbatim
!>          THRESH is DOUBLE PRECISION
!>          The threshold value for the test ratios.  A result is
!>          included in the output file if RESULT >= THRESH.  To have
!>          every test ratio printed, use THRESH = 0.
!> \endverbatim
!>
!> \param[in] TSTERR
!> \verbatim
!>          TSTERR is LOGICAL
!>          Flag that indicates whether error exits are to be tested.
!> \endverbatim
!>
!> \param[out] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension (LA)
!> \endverbatim
!>
!> \param[in] LA
!> \verbatim
!>          LA is INTEGER
!>          The length of the array A.  LA >= (2*NMAX-1)*NMAX
!>          where NMAX is the largest entry in NVAL.
!> \endverbatim
!>
!> \param[out] AFB
!> \verbatim
!>          AFB is COMPLEX*16 array, dimension (LAFB)
!> \endverbatim
!>
!> \param[in] LAFB
!> \verbatim
!>          LAFB is INTEGER
!>          The length of the array AFB.  LAFB >= (3*NMAX-2)*NMAX
!>          where NMAX is the largest entry in NVAL.
!> \endverbatim
!>
!> \param[out] ASAV
!> \verbatim
!>          ASAV is COMPLEX*16 array, dimension (LA)
!> \endverbatim
!>
!> \param[out] B
!> \verbatim
!>          B is COMPLEX*16 array, dimension (NMAX*NRHS)
!> \endverbatim
!>
!> \param[out] BSAV
!> \verbatim
!>          BSAV is COMPLEX*16 array, dimension (NMAX*NRHS)
!> \endverbatim
!>
!> \param[out] X
!> \verbatim
!>          X is COMPLEX*16 array, dimension (NMAX*NRHS)
!> \endverbatim
!>
!> \param[out] XACT
!> \verbatim
!>          XACT is COMPLEX*16 array, dimension (NMAX*NRHS)
!> \endverbatim
!>
!> \param[out] S
!> \verbatim
!>          S is DOUBLE PRECISION array, dimension (2*NMAX)
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX*16 array, dimension
!>                      (NMAX*max(3,NRHS,NMAX))
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is DOUBLE PRECISION array, dimension
!>                      (max(NMAX,2*NRHS))
!> \endverbatim
!>
!> \param[out] IWORK
!> \verbatim
!>          IWORK is INTEGER array, dimension (NMAX)
!> \endverbatim
!>
!> \param[in] NOUT
!> \verbatim
!>          NOUT is INTEGER
!>          The unit number for output.
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
!> \ingroup complex16_lin
!
!  =====================================================================
      SUBROUTINE ZDRVGB(Dotype,Nn,Nval,Nrhs,Thresh,Tsterr,A,La,Afb,Lafb,&
     &                  Asav,B,Bsav,X,Xact,S,Work,Rwork,Iwork,Nout)
      IMPLICIT NONE
!*--ZDRVGB175
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      LOGICAL Tsterr
      INTEGER La , Lafb , Nn , Nout , Nrhs
      DOUBLE PRECISION Thresh
!     ..
!     .. Array Arguments ..
      LOGICAL Dotype(*)
      INTEGER Iwork(*) , Nval(*)
      DOUBLE PRECISION Rwork(*) , S(*)
      COMPLEX*16 A(*) , Afb(*) , Asav(*) , B(*) , Bsav(*) , Work(*) ,   &
     &           X(*) , Xact(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ONE , ZERO
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
      INTEGER NTYPES
      PARAMETER (NTYPES=8)
      INTEGER NTESTS
      PARAMETER (NTESTS=7)
      INTEGER NTRAN
      PARAMETER (NTRAN=3)
!     ..
!     .. Local Scalars ..
      LOGICAL equil , nofact , prefac , trfcon , zerot
      CHARACTER dist , equed , fact , trans , type , xtype
      CHARACTER*3 path
      INTEGER i , i1 , i2 , iequed , ifact , ikl , iku , imat , in ,    &
     &        info , ioff , itran , izero , j , k , k1 , kl , ku , lda ,&
     &        ldafb , ldb , mode , n , nb , nbmin , nerrs , nfact ,     &
     &        nfail , nimat , nkl , nku , nrun , nt
      DOUBLE PRECISION ainvnm , amax , anorm , anormi , anormo ,        &
     &                 anrmpv , cndnum , colcnd , rcond , rcondc ,      &
     &                 rcondi , rcondo , roldc , roldi , roldo ,        &
     &                 rowcnd , rpvgrw
!     ..
!     .. Local Arrays ..
      CHARACTER equeds(4) , facts(3) , transs(NTRAN)
      INTEGER iseed(4) , iseedy(4)
      DOUBLE PRECISION rdum(1) , result(NTESTS)
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      DOUBLE PRECISION DGET06 , DLAMCH , ZLANGB , ZLANGE , ZLANTB
      EXTERNAL LSAME , DGET06 , DLAMCH , ZLANGB , ZLANGE , ZLANTB
!     ..
!     .. External Subroutines ..
      EXTERNAL ALADHD , ALAERH , ALASVM , XLAENV , ZERRVX , ZGBEQU ,    &
     &         ZGBSV , ZGBSVX , ZGBT01 , ZGBT02 , ZGBT05 , ZGBTRF ,     &
     &         ZGBTRS , ZGET04 , ZLACPY , ZLAQGB , ZLARHS , ZLASET ,    &
     &         ZLATB4 , ZLATMS
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , DCMPLX , MAX , MIN
!     ..
!     .. Scalars in Common ..
      LOGICAL LERr , OK
      CHARACTER*32 SRNamt
      INTEGER INFot , NUNit
!     ..
!     .. Common blocks ..
      COMMON /INFOC / INFot , NUNit , OK , LERr
      COMMON /SRNAMC/ SRNamt
!     ..
!     .. Data statements ..
      DATA iseedy/1988 , 1989 , 1990 , 1991/
      DATA transs/'N' , 'T' , 'C'/
      DATA facts/'F' , 'N' , 'E'/
      DATA equeds/'N' , 'R' , 'C' , 'B'/
!     ..
!     .. Executable Statements ..
!
!     Initialize constants and the random number seed.
!
      path(1:1) = 'Zomplex precision'
      path(2:3) = 'GB'
      nrun = 0
      nfail = 0
      nerrs = 0
      DO i = 1 , 4
         iseed(i) = iseedy(i)
      ENDDO
!
!     Test the error exits
!
      IF ( Tsterr ) CALL ZERRVX(path,Nout)
      INFot = 0
!
!     Set the block size and minimum block size for testing.
!
      nb = 1
      nbmin = 2
      CALL XLAENV(1,nb)
      CALL XLAENV(2,nbmin)
!
!     Do for each value of N in NVAL
!
      DO in = 1 , Nn
         n = Nval(in)
         ldb = MAX(n,1)
         xtype = 'N'
!
!        Set limits on the number of loop iterations.
!
         nkl = MAX(1,MIN(n,4))
         IF ( n==0 ) nkl = 1
         nku = nkl
         nimat = NTYPES
         IF ( n<=0 ) nimat = 1
!
         DO ikl = 1 , nkl
!
!           Do for KL = 0, N-1, (3N-1)/4, and (N+1)/4. This order makes
!           it easier to skip redundant values for small values of N.
!
            IF ( ikl==1 ) THEN
               kl = 0
            ELSEIF ( ikl==2 ) THEN
               kl = MAX(n-1,0)
            ELSEIF ( ikl==3 ) THEN
               kl = (3*n-1)/4
            ELSEIF ( ikl==4 ) THEN
               kl = (n+1)/4
            ENDIF
            DO iku = 1 , nku
!
!              Do for KU = 0, N-1, (3N-1)/4, and (N+1)/4. This order
!              makes it easier to skip redundant values for small
!              values of N.
!
               IF ( iku==1 ) THEN
                  ku = 0
               ELSEIF ( iku==2 ) THEN
                  ku = MAX(n-1,0)
               ELSEIF ( iku==3 ) THEN
                  ku = (3*n-1)/4
               ELSEIF ( iku==4 ) THEN
                  ku = (n+1)/4
               ENDIF
!
!              Check that A and AFB are big enough to generate this
!              matrix.
!
               lda = kl + ku + 1
               ldafb = 2*kl + ku + 1
               IF ( lda*n>La .OR. ldafb*n>Lafb ) THEN
                  IF ( nfail==0 .AND. nerrs==0 ) CALL ALADHD(Nout,path)
                  IF ( lda*n>La ) THEN
                     WRITE (Nout,FMT=99001) La , n , kl , ku ,          &
     &                      n*(kl+ku+1)
                     nerrs = nerrs + 1
                  ENDIF
                  IF ( ldafb*n>Lafb ) THEN
                     WRITE (Nout,FMT=99002) Lafb , n , kl , ku ,        &
     &                      n*(2*kl+ku+1)
                     nerrs = nerrs + 1
                  ENDIF
                  CYCLE
               ENDIF
!
               DO imat = 1 , nimat
!
!                 Do the tests only if DOTYPE( IMAT ) is true.
!
                  IF ( Dotype(imat) ) THEN
!
!                 Skip types 2, 3, or 4 if the matrix is too small.
!
                     zerot = imat>=2 .AND. imat<=4
                     IF ( .NOT.(zerot .AND. n<imat-1) ) THEN
!
!                 Set up parameters with ZLATB4 and generate a
!                 test matrix with ZLATMS.
!
                        CALL ZLATB4(path,imat,n,n,type,kl,ku,anorm,mode,&
     &                              cndnum,dist)
                        rcondc = ONE/cndnum
!
                        SRNamt = 'ZLATMS'
                        CALL ZLATMS(n,n,dist,iseed,type,Rwork,mode,     &
     &                              cndnum,anorm,kl,ku,'Z',A,lda,Work,  &
     &                              info)
!
!                 Check the error code from ZLATMS.
!
                        IF ( info/=0 ) THEN
                           CALL ALAERH(path,'ZLATMS',info,0,' ',n,n,kl, &
     &                                 ku,-1,imat,nfail,nerrs,Nout)
                           CYCLE
                        ENDIF
!
!                 For types 2, 3, and 4, zero one or more columns of
!                 the matrix to test that INFO is returned correctly.
!
                        izero = 0
                        IF ( zerot ) THEN
                           IF ( imat==2 ) THEN
                              izero = 1
                           ELSEIF ( imat==3 ) THEN
                              izero = n
                           ELSE
                              izero = n/2 + 1
                           ENDIF
                           ioff = (izero-1)*lda
                           IF ( imat<4 ) THEN
                              i1 = MAX(1,ku+2-izero)
                              i2 = MIN(kl+ku+1,ku+1+(n-izero))
                              DO i = i1 , i2
                                 A(ioff+i) = ZERO
                              ENDDO
                           ELSE
                              DO j = izero , n
                                 DO i = MAX(1,ku+2-j) ,                 &
     &                              MIN(kl+ku+1,ku+1+(n-j))
                                    A(ioff+i) = ZERO
                                 ENDDO
                                 ioff = ioff + lda
                              ENDDO
                           ENDIF
                        ENDIF
!
!                 Save a copy of the matrix A in ASAV.
!
                        CALL ZLACPY('Full',kl+ku+1,n,A,lda,Asav,lda)
!
                        DO iequed = 1 , 4
                           equed = equeds(iequed)
                           IF ( iequed==1 ) THEN
                              nfact = 3
                           ELSE
                              nfact = 1
                           ENDIF
!
                           DO ifact = 1 , nfact
                              fact = facts(ifact)
                              prefac = LSAME(fact,'F')
                              nofact = LSAME(fact,'N')
                              equil = LSAME(fact,'E')
!
                              IF ( zerot ) THEN
                                 IF ( prefac ) CYCLE
                                 rcondo = ZERO
                                 rcondi = ZERO
!
                              ELSEIF ( .NOT.nofact ) THEN
!
!                          Compute the condition number for comparison
!                          with the value returned by DGESVX (FACT =
!                          'N' reuses the condition number from the
!                          previous iteration with FACT = 'F').
!
                                 CALL ZLACPY('Full',kl+ku+1,n,Asav,lda, &
     &                              Afb(kl+1),ldafb)
                                 IF ( equil .OR. iequed>1 ) THEN
!
!                             Compute row and column scale factors to
!                             equilibrate the matrix A.
!
                                    CALL ZGBEQU(n,n,kl,ku,Afb(kl+1),    &
     &                                 ldafb,S,S(n+1),rowcnd,colcnd,    &
     &                                 amax,info)
                                    IF ( info==0 .AND. n>0 ) THEN
                                       IF ( LSAME(equed,'R') ) THEN
                                         rowcnd = ZERO
                                         colcnd = ONE
                                       ELSEIF ( LSAME(equed,'C') ) THEN
                                         rowcnd = ONE
                                         colcnd = ZERO
                                       ELSEIF ( LSAME(equed,'B') ) THEN
                                         rowcnd = ZERO
                                         colcnd = ZERO
                                       ENDIF
!
!                                Equilibrate the matrix.
!
                                       CALL ZLAQGB(n,n,kl,ku,Afb(kl+1), &
     &                                    ldafb,S,S(n+1),rowcnd,colcnd, &
     &                                    amax,equed)
                                    ENDIF
                                 ENDIF
!
!                          Save the condition number of the
!                          non-equilibrated system for use in ZGET04.
!
                                 IF ( equil ) THEN
                                    roldo = rcondo
                                    roldi = rcondi
                                 ENDIF
!
!                          Compute the 1-norm and infinity-norm of A.
!
                                 anormo = ZLANGB('1',n,kl,ku,Afb(kl+1), &
     &                              ldafb,Rwork)
                                 anormi = ZLANGB('I',n,kl,ku,Afb(kl+1), &
     &                              ldafb,Rwork)
!
!                          Factor the matrix A.
!
                                 CALL ZGBTRF(n,n,kl,ku,Afb,ldafb,Iwork, &
     &                              info)
!
!                          Form the inverse of A.
!
                                 CALL ZLASET('Full',n,n,DCMPLX(ZERO),   &
     &                              DCMPLX(ONE),Work,ldb)
                                 SRNamt = 'ZGBTRS'
                                 CALL ZGBTRS('No transpose',n,kl,ku,n,  &
     &                              Afb,ldafb,Iwork,Work,ldb,info)
!
!                          Compute the 1-norm condition number of A.
!
                                 ainvnm = ZLANGE('1',n,n,Work,ldb,Rwork)
                                 IF ( anormo<=ZERO .OR. ainvnm<=ZERO )  &
     &                                THEN
                                    rcondo = ONE
                                 ELSE
                                    rcondo = (ONE/anormo)/ainvnm
                                 ENDIF
!
!                          Compute the infinity-norm condition number
!                          of A.
!
                                 ainvnm = ZLANGE('I',n,n,Work,ldb,Rwork)
                                 IF ( anormi<=ZERO .OR. ainvnm<=ZERO )  &
     &                                THEN
                                    rcondi = ONE
                                 ELSE
                                    rcondi = (ONE/anormi)/ainvnm
                                 ENDIF
                              ENDIF
!
                              DO itran = 1 , NTRAN
!
!                          Do for each value of TRANS.
!
                                 trans = transs(itran)
                                 IF ( itran==1 ) THEN
                                    rcondc = rcondo
                                 ELSE
                                    rcondc = rcondi
                                 ENDIF
!
!                          Restore the matrix A.
!
                                 CALL ZLACPY('Full',kl+ku+1,n,Asav,lda, &
     &                              A,lda)
!
!                          Form an exact solution and set the right hand
!                          side.
!
                                 SRNamt = 'ZLARHS'
                                 CALL ZLARHS(path,xtype,'Full',trans,n, &
     &                              n,kl,ku,Nrhs,A,lda,Xact,ldb,B,ldb,  &
     &                              iseed,info)
                                 xtype = 'C'
                                 CALL ZLACPY('Full',n,Nrhs,B,ldb,Bsav,  &
     &                              ldb)
!
                                 IF ( nofact .AND. itran==1 ) THEN
!
!                             --- Test ZGBSV  ---
!
!                             Compute the LU factorization of the matrix
!                             and solve the system.
!
                                    CALL ZLACPY('Full',kl+ku+1,n,A,lda, &
     &                                 Afb(kl+1),ldafb)
                                    CALL ZLACPY('Full',n,Nrhs,B,ldb,X,  &
     &                                 ldb)
!
                                    SRNamt = 'ZGBSV '
                                    CALL ZGBSV(n,kl,ku,Nrhs,Afb,ldafb,  &
     &                                 Iwork,X,ldb,info)
!
!                             Check error code from ZGBSV .
!
                                    IF ( info/=izero )                  &
     &                                  CALL ALAERH(path,'ZGBSV ',info, &
     &                                 izero,' ',n,n,kl,ku,Nrhs,imat,   &
     &                                 nfail,nerrs,Nout)
!
!                             Reconstruct matrix from factors and
!                             compute residual.
!
                                    CALL ZGBT01(n,n,kl,ku,A,lda,Afb,    &
     &                                 ldafb,Iwork,Work,result(1))
                                    nt = 1
                                    IF ( izero==0 ) THEN
!
!                                Compute residual of the computed
!                                solution.
!
                                       CALL ZLACPY('Full',n,Nrhs,B,ldb, &
     &                                    Work,ldb)
                                       CALL ZGBT02('No transpose',n,n,  &
     &                                    kl,ku,Nrhs,A,lda,X,ldb,Work,  &
     &                                    ldb,result(2))
!
!                                Check solution from generated exact
!                                solution.
!
                                       CALL ZGET04(n,Nrhs,X,ldb,Xact,   &
     &                                    ldb,rcondc,result(3))
                                       nt = 3
                                    ENDIF
!
!                             Print information about the tests that did
!                             not pass the threshold.
!
                                    DO k = 1 , nt
                                       IF ( result(k)>=Thresh ) THEN
                                         IF ( nfail==0 .AND. nerrs==0 ) &
     &                                      CALL ALADHD(Nout,path)
                                         WRITE (Nout,FMT=99003)         &
     &                                      'ZGBSV ' , n , kl , ku ,    &
     &                                      imat , k , result(k)
                                         nfail = nfail + 1
                                       ENDIF
                                    ENDDO
                                    nrun = nrun + nt
                                 ENDIF
!
!                          --- Test ZGBSVX ---
!
                                 IF ( .NOT.prefac )                     &
     &                                CALL ZLASET('Full',2*kl+ku+1,n,   &
     &                                DCMPLX(ZERO),DCMPLX(ZERO),Afb,    &
     &                                ldafb)
                                 CALL ZLASET('Full',n,Nrhs,DCMPLX(ZERO),&
     &                              DCMPLX(ZERO),X,ldb)
!
!                             Equilibrate the matrix if FACT = 'F' and
!                             EQUED = 'R', 'C', or 'B'.
!
                                 IF ( iequed>1 .AND. n>0 )              &
     &                                CALL ZLAQGB(n,n,kl,ku,A,lda,S,    &
     &                                S(n+1),rowcnd,colcnd,amax,equed)
!
!                          Solve the system and compute the condition
!                          number and error bounds using ZGBSVX.
!
                                 SRNamt = 'ZGBSVX'
                                 CALL ZGBSVX(fact,trans,n,kl,ku,Nrhs,A, &
     &                              lda,Afb,ldafb,Iwork,equed,S,S(ldb+1)&
     &                              ,B,ldb,X,ldb,rcond,Rwork,           &
     &                              Rwork(Nrhs+1),Work,Rwork(2*Nrhs+1), &
     &                              info)
!
!                          Check the error code from ZGBSVX.
!
                                 IF ( info/=izero )                     &
     &                                 CALL ALAERH(path,'ZGBSVX',info,  &
     &                                izero,fact//trans,n,n,kl,ku,Nrhs, &
     &                                imat,nfail,nerrs,Nout)
!                          Compare RWORK(2*NRHS+1) from ZGBSVX with the
!                          computed reciprocal pivot growth RPVGRW
!
                                 IF ( info/=0 .AND. info<=n ) THEN
                                    anrmpv = ZERO
                                    DO j = 1 , info
                                       DO i = MAX(ku+2-j,1) ,           &
     &                                    MIN(n+ku+1-j,kl+ku+1)
                                         anrmpv = MAX(anrmpv,           &
     &                                      ABS(A(i+(j-1)*lda)))
                                       ENDDO
                                    ENDDO
                                    rpvgrw = ZLANTB('M','U','N',info,   &
     &                                 MIN(info-1,kl+ku),               &
     &                                 Afb(MAX(1,kl+ku+2-info)),ldafb,  &
     &                                 rdum)
                                    IF ( rpvgrw==ZERO ) THEN
                                       rpvgrw = ONE
                                    ELSE
                                       rpvgrw = anrmpv/rpvgrw
                                    ENDIF
                                 ELSE
                                    rpvgrw = ZLANTB('M','U','N',n,kl+ku,&
     &                                 Afb,ldafb,rdum)
                                    IF ( rpvgrw==ZERO ) THEN
                                       rpvgrw = ONE
                                    ELSE
                                       rpvgrw = ZLANGB('M',n,kl,ku,A,   &
     &                                    lda,rdum)/rpvgrw
                                    ENDIF
                                 ENDIF
                                 result(7) = ABS(rpvgrw-Rwork(2*Nrhs+1))&
     &                              /MAX(Rwork(2*Nrhs+1),rpvgrw)        &
     &                              /DLAMCH('E')
!
                                 IF ( .NOT.prefac ) THEN
!
!                             Reconstruct matrix from factors and
!                             compute residual.
!
                                    CALL ZGBT01(n,n,kl,ku,A,lda,Afb,    &
     &                                 ldafb,Iwork,Work,result(1))
                                    k1 = 1
                                 ELSE
                                    k1 = 2
                                 ENDIF
!
                                 IF ( info==0 ) THEN
                                    trfcon = .FALSE.
!
!                             Compute residual of the computed solution.
!
                                    CALL ZLACPY('Full',n,Nrhs,Bsav,ldb, &
     &                                 Work,ldb)
                                    CALL ZGBT02(trans,n,n,kl,ku,Nrhs,   &
     &                                 Asav,lda,X,ldb,Work,ldb,result(2)&
     &                                 )
!
!                             Check solution from generated exact
!                             solution.
!
                                    IF ( nofact .OR.                    &
     &                                 (prefac .AND. LSAME(equed,'N')) )&
     &                                 THEN
                                       CALL ZGET04(n,Nrhs,X,ldb,Xact,   &
     &                                    ldb,rcondc,result(3))
                                    ELSE
                                       IF ( itran==1 ) THEN
                                         roldc = roldo
                                       ELSE
                                         roldc = roldi
                                       ENDIF
                                       CALL ZGET04(n,Nrhs,X,ldb,Xact,   &
     &                                    ldb,roldc,result(3))
                                    ENDIF
!
!                             Check the error bounds from iterative
!                             refinement.
!
                                    CALL ZGBT05(trans,n,kl,ku,Nrhs,Asav,&
     &                                 lda,Bsav,ldb,X,ldb,Xact,ldb,     &
     &                                 Rwork,Rwork(Nrhs+1),result(4))
                                 ELSE
                                    trfcon = .TRUE.
                                 ENDIF
!
!                          Compare RCOND from ZGBSVX with the computed
!                          value in RCONDC.
!
                                 result(6) = DGET06(rcond,rcondc)
!
!                          Print information about the tests that did
!                          not pass the threshold.
!
                                 IF ( .NOT.trfcon ) THEN
                                    DO k = k1 , NTESTS
                                       IF ( result(k)>=Thresh ) THEN
                                         IF ( nfail==0 .AND. nerrs==0 ) &
     &                                      CALL ALADHD(Nout,path)
                                         IF ( prefac ) THEN
                                         WRITE (Nout,FMT=99005)         &
     &                                      'ZGBSVX' , fact , trans ,   &
     &                                      n , kl , ku , equed , imat ,&
     &                                      k , result(k)
                                         ELSE
                                         WRITE (Nout,FMT=99004)         &
     &                                      'ZGBSVX' , fact , trans ,   &
     &                                      n , kl , ku , imat , k ,    &
     &                                      result(k)
                                         ENDIF
                                         nfail = nfail + 1
                                       ENDIF
                                    ENDDO
                                    nrun = nrun + NTESTS - k1 + 1
                                 ELSE
                                    IF ( result(1)>=Thresh .AND.        &
     &                                 .NOT.prefac ) THEN
                                       IF ( nfail==0 .AND. nerrs==0 )   &
     &                                    CALL ALADHD(Nout,path)
                                       IF ( prefac ) THEN
                                         WRITE (Nout,FMT=99005)         &
     &                                      'ZGBSVX' , fact , trans ,   &
     &                                      n , kl , ku , equed , imat ,&
     &                                      1 , result(1)
                                       ELSE
                                         WRITE (Nout,FMT=99004)         &
     &                                      'ZGBSVX' , fact , trans ,   &
     &                                      n , kl , ku , imat , 1 ,    &
     &                                      result(1)
                                       ENDIF
                                       nfail = nfail + 1
                                       nrun = nrun + 1
                                    ENDIF
                                    IF ( result(6)>=Thresh ) THEN
                                       IF ( nfail==0 .AND. nerrs==0 )   &
     &                                    CALL ALADHD(Nout,path)
                                       IF ( prefac ) THEN
                                         WRITE (Nout,FMT=99005)         &
     &                                      'ZGBSVX' , fact , trans ,   &
     &                                      n , kl , ku , equed , imat ,&
     &                                      6 , result(6)
                                       ELSE
                                         WRITE (Nout,FMT=99004)         &
     &                                      'ZGBSVX' , fact , trans ,   &
     &                                      n , kl , ku , imat , 6 ,    &
     &                                      result(6)
                                       ENDIF
                                       nfail = nfail + 1
                                       nrun = nrun + 1
                                    ENDIF
                                    IF ( result(7)>=Thresh ) THEN
                                       IF ( nfail==0 .AND. nerrs==0 )   &
     &                                    CALL ALADHD(Nout,path)
                                       IF ( prefac ) THEN
                                         WRITE (Nout,FMT=99005)         &
     &                                      'ZGBSVX' , fact , trans ,   &
     &                                      n , kl , ku , equed , imat ,&
     &                                      7 , result(7)
                                       ELSE
                                         WRITE (Nout,FMT=99004)         &
     &                                      'ZGBSVX' , fact , trans ,   &
     &                                      n , kl , ku , imat , 7 ,    &
     &                                      result(7)
                                       ENDIF
                                       nfail = nfail + 1
                                       nrun = nrun + 1
                                    ENDIF
                                 ENDIF
                              ENDDO
                           ENDDO
                        ENDDO
                     ENDIF
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
      ENDDO
!
!     Print a summary of the results.
!
      CALL ALASVM(path,Nout,nfail,nrun,nerrs)
!
99001 FORMAT (' *** In ZDRVGB, LA=',I5,' is too small for N=',I5,       &
     &        ', KU=',I5,', KL=',I5,/' ==> Increase LA to at least ',I5)
99002 FORMAT (' *** In ZDRVGB, LAFB=',I5,' is too small for N=',I5,     &
     &        ', KU=',I5,', KL=',I5,/' ==> Increase LAFB to at least ', &
     &        I5)
99003 FORMAT (1X,A,', N=',I5,', KL=',I5,', KU=',I5,', type ',I1,        &
     &        ', test(',I1,')=',G12.5)
99004 FORMAT (1X,A,'( ''',A1,''',''',A1,''',',I5,',',I5,',',I5,         &
     &        ',...), type ',I1,', test(',I1,')=',G12.5)
99005 FORMAT (1X,A,'( ''',A1,''',''',A1,''',',I5,',',I5,',',I5,         &
     &        ',...), EQUED=''',A1,''', type ',I1,', test(',I1,')=',    &
     &        G12.5)
!
!
!     End of ZDRVGB
!
      END SUBROUTINE ZDRVGB
