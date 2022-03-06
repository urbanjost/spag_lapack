!*==zdrvge.f90  processed by SPAG 7.51RB at 17:35 on  4 Mar 2022
!> \brief \b zdrvge
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE zdrvge( DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, NMAX,
!                          A, AFAC, ASAV, B, BSAV, X, XACT, S, WORK,
!                          RWORK, IWORK, NOUT )
!
!       .. Scalar Arguments ..
!       LOGICAL            TSTERR
!       INTEGER            NMAX, NN, NOUT, NRHS
!       DOUBLE PRECISION   THRESH
!       ..
!       .. Array Arguments ..
!       LOGICAL            DOTYPE( * )
!       INTEGER            IWORK( * ), NVAL( * )
!       DOUBLE PRECISION   RWORK( * ), S( * )
!       COMPLEX*16         A( * ), AFAC( * ), ASAV( * ), B( * ),
!      $                   BSAV( * ), WORK( * ), X( * ), XACT( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> zdrvge tests the driver routines ZGESV and -SVX.
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
!> \param[in] NMAX
!> \verbatim
!>          NMAX is INTEGER
!>          The maximum value permitted for N, used in dimensioning the
!>          work arrays.
!> \endverbatim
!>
!> \param[out] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension (NMAX*NMAX)
!> \endverbatim
!>
!> \param[out] AFAC
!> \verbatim
!>          AFAC is COMPLEX*16 array, dimension (NMAX*NMAX)
!> \endverbatim
!>
!> \param[out] ASAV
!> \verbatim
!>          ASAV is COMPLEX*16 array, dimension (NMAX*NMAX)
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
!>                      (NMAX*max(3,NRHS))
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is DOUBLE PRECISION array, dimension (2*NRHS+NMAX)
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
      SUBROUTINE zdrvge(Dotype,Nn,Nval,Nrhs,Thresh,Tsterr,Nmax,A,Afac,  &
     &                  Asav,B,Bsav,X,Xact,S,Work,Rwork,Iwork,Nout)
      IMPLICIT NONE
!*--ZDRVGE167
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      LOGICAL Tsterr
      INTEGER Nmax , Nn , Nout , Nrhs
      DOUBLE PRECISION Thresh
!     ..
!     .. Array Arguments ..
      LOGICAL Dotype(*)
      INTEGER Iwork(*) , Nval(*)
      DOUBLE PRECISION Rwork(*) , S(*)
      COMPLEX*16 A(*) , Afac(*) , Asav(*) , B(*) , Bsav(*) , Work(*) ,  &
     &           X(*) , Xact(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ONE , ZERO
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
      INTEGER NTYPES
      PARAMETER (NTYPES=11)
      INTEGER NTESTS
      PARAMETER (NTESTS=7)
      INTEGER NTRAN
      PARAMETER (NTRAN=3)
!     ..
!     .. Local Scalars ..
      LOGICAL equil , nofact , prefac , trfcon , zerot
      CHARACTER dist , equed , fact , trans , type , xtype
      CHARACTER*3 path
      INTEGER i , iequed , ifact , imat , in , info , ioff , itran ,    &
     &        izero , k , k1 , kl , ku , lda , lwork , mode , n , nb ,  &
     &        nbmin , nerrs , nfact , nfail , nimat , nrun , nt
      DOUBLE PRECISION ainvnm , amax , anorm , anormi , anormo ,        &
     &                 cndnum , colcnd , rcond , rcondc , rcondi ,      &
     &                 rcondo , roldc , roldi , roldo , rowcnd , rpvgrw
!     ..
!     .. Local Arrays ..
      CHARACTER equeds(4) , facts(3) , transs(NTRAN)
      INTEGER iseed(4) , iseedy(4)
      DOUBLE PRECISION rdum(1) , result(NTESTS)
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      DOUBLE PRECISION DLAMCH , ZLANGE , ZLANTR
EXTERNAL    LSAME,   DLAMCH,   ZLANGE,   ZLANTR
!     ..
!     .. External Subroutines ..
EXTERNAL    ZGEEQU
EXTERNAL    ZGESV,   ZGESVX
EXTERNAL    ZGETRF,   ZGETRI,   ZLACPY,   ZLAQGE,   ZLASET
EXTERNAL    ZLATMS
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , DCMPLX , MAX
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
      path(2:3) = 'GE'
      nrun = 0
      nfail = 0
      nerrs = 0
      DO i = 1 , 4
         iseed(i) = iseedy(i)
      ENDDO
!
!     Test the error exits
!
      IF ( Tsterr ) CALL zerrvx(path,Nout)
      INFot = 0
!
!     Set the block size and minimum block size for testing.
!
      nb = 1
      nbmin = 2
      CALL xlaenv(1,nb)
      CALL xlaenv(2,nbmin)
!
!     Do for each value of N in NVAL
!
      DO in = 1 , Nn
         n = Nval(in)
         lda = MAX(n,1)
         xtype = 'N'
         nimat = NTYPES
         IF ( n<=0 ) nimat = 1
!
         DO imat = 1 , nimat
!
!           Do the tests only if DOTYPE( IMAT ) is true.
!
            IF ( Dotype(imat) ) THEN
!
!           Skip types 5, 6, or 7 if the matrix size is too small.
!
               zerot = imat>=5 .AND. imat<=7
               IF ( .NOT.(zerot .AND. n<imat-4) ) THEN
!
!           Set up parameters with zlatb4 and generate a test matrix
!           with ZLATMS.
!
                  CALL zlatb4(path,imat,n,n,type,kl,ku,anorm,mode,      &
     &                        cndnum,dist)
                  rcondc = ONE/cndnum
!
                  SRNamt = 'ZLATMS'
                  CALL ZLATMS(n,n,dist,iseed,type,Rwork,mode,cndnum,    &
     &                        anorm,kl,ku,'No packing',A,lda,Work,info)
!
!           Check error code from ZLATMS.
!
                  IF ( info/=0 ) THEN
                     CALL alaerh(path,'ZLATMS',info,0,' ',n,n,-1,-1,-1, &
     &                           imat,nfail,nerrs,Nout)
                     CYCLE
                  ENDIF
!
!           For types 5-7, zero one or more columns of the matrix to
!           test that INFO is returned correctly.
!
                  IF ( zerot ) THEN
                     IF ( imat==5 ) THEN
                        izero = 1
                     ELSEIF ( imat==6 ) THEN
                        izero = n
                     ELSE
                        izero = n/2 + 1
                     ENDIF
                     ioff = (izero-1)*lda
                     IF ( imat<7 ) THEN
                        DO i = 1 , n
                           A(ioff+i) = ZERO
                        ENDDO
                     ELSE
                        CALL ZLASET('Full',n,n-izero+1,DCMPLX(ZERO),    &
     &                              DCMPLX(ZERO),A(ioff+1),lda)
                     ENDIF
                  ELSE
                     izero = 0
                  ENDIF
!
!           Save a copy of the matrix A in ASAV.
!
                  CALL ZLACPY('Full',n,n,A,lda,Asav,lda)
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
!                    Compute the condition number for comparison with
!                    the value returned by ZGESVX (FACT = 'N' reuses
!                    the condition number from the previous iteration
!                    with FACT = 'F').
!
                           CALL ZLACPY('Full',n,n,Asav,lda,Afac,lda)
                           IF ( equil .OR. iequed>1 ) THEN
!
!                       Compute row and column scale factors to
!                       equilibrate the matrix A.
!
                              CALL ZGEEQU(n,n,Afac,lda,S,S(n+1),rowcnd, &
     &                           colcnd,amax,info)
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
!                          Equilibrate the matrix.
!
                                 CALL ZLAQGE(n,n,Afac,lda,S,S(n+1),     &
     &                              rowcnd,colcnd,amax,equed)
                              ENDIF
                           ENDIF
!
!                    Save the condition number of the non-equilibrated
!                    system for use in zget04.
!
                           IF ( equil ) THEN
                              roldo = rcondo
                              roldi = rcondi
                           ENDIF
!
!                    Compute the 1-norm and infinity-norm of A.
!
                           anormo = ZLANGE('1',n,n,Afac,lda,Rwork)
                           anormi = ZLANGE('I',n,n,Afac,lda,Rwork)
!
!                    Factor the matrix A.
!
                           SRNamt = 'ZGETRF'
                           CALL ZGETRF(n,n,Afac,lda,Iwork,info)
!
!                    Form the inverse of A.
!
                           CALL ZLACPY('Full',n,n,Afac,lda,A,lda)
                           lwork = Nmax*MAX(3,Nrhs)
                           SRNamt = 'ZGETRI'
                           CALL ZGETRI(n,A,lda,Iwork,Work,lwork,info)
!
!                    Compute the 1-norm condition number of A.
!
                           ainvnm = ZLANGE('1',n,n,A,lda,Rwork)
                           IF ( anormo<=ZERO .OR. ainvnm<=ZERO ) THEN
                              rcondo = ONE
                           ELSE
                              rcondo = (ONE/anormo)/ainvnm
                           ENDIF
!
!                    Compute the infinity-norm condition number of A.
!
                           ainvnm = ZLANGE('I',n,n,A,lda,Rwork)
                           IF ( anormi<=ZERO .OR. ainvnm<=ZERO ) THEN
                              rcondi = ONE
                           ELSE
                              rcondi = (ONE/anormi)/ainvnm
                           ENDIF
                        ENDIF
!
                        DO itran = 1 , NTRAN
!
!                    Do for each value of TRANS.
!
                           trans = transs(itran)
                           IF ( itran==1 ) THEN
                              rcondc = rcondo
                           ELSE
                              rcondc = rcondi
                           ENDIF
!
!                    Restore the matrix A.
!
                           CALL ZLACPY('Full',n,n,Asav,lda,A,lda)
!
!                    Form an exact solution and set the right hand side.
!
                           SRNamt = 'zlarhs'
                           CALL zlarhs(path,xtype,'Full',trans,n,n,kl,  &
     &                                 ku,Nrhs,A,lda,Xact,lda,B,lda,    &
     &                                 iseed,info)
                           xtype = 'C'
                           CALL ZLACPY('Full',n,Nrhs,B,lda,Bsav,lda)
!
                           IF ( nofact .AND. itran==1 ) THEN
!
!                       --- Test ZGESV  ---
!
!                       Compute the LU factorization of the matrix and
!                       solve the system.
!
                              CALL ZLACPY('Full',n,n,A,lda,Afac,lda)
                              CALL ZLACPY('Full',n,Nrhs,B,lda,X,lda)
!
                              SRNamt = 'ZGESV '
                              CALL ZGESV(n,Nrhs,Afac,lda,Iwork,X,lda,   &
     &                           info)
!
!                       Check error code from ZGESV .
!
                              IF ( info/=izero )                        &
     &                              CALL alaerh(path,'ZGESV ',info,     &
     &                             izero,' ',n,n,-1,-1,Nrhs,imat,nfail, &
     &                             nerrs,Nout)
!
!                       Reconstruct matrix from factors and compute
!                       residual.
!
                              CALL zget01(n,n,A,lda,Afac,lda,Iwork,     &
     &                           Rwork,result(1))
                              nt = 1
                              IF ( izero==0 ) THEN
!
!                          Compute residual of the computed solution.
!
                                 CALL ZLACPY('Full',n,Nrhs,B,lda,Work,  &
     &                              lda)
                                 CALL zget02('No transpose',n,n,Nrhs,A, &
     &                              lda,X,lda,Work,lda,Rwork,result(2))
!
!                          Check solution from generated exact solution.
!
                                 CALL zget04(n,Nrhs,X,lda,Xact,lda,     &
     &                              rcondc,result(3))
                                 nt = 3
                              ENDIF
!
!                       Print information about the tests that did not
!                       pass the threshold.
!
                              DO k = 1 , nt
                                 IF ( result(k)>=Thresh ) THEN
                                    IF ( nfail==0 .AND. nerrs==0 )      &
     &                                 CALL aladhd(Nout,path)
                                    WRITE (Nout,FMT=99001) 'ZGESV ' ,   &
     &                                 n , imat , k , result(k)
                                    nfail = nfail + 1
                                 ENDIF
                              ENDDO
                              nrun = nrun + nt
                           ENDIF
!
!                    --- Test ZGESVX ---
!
                           IF ( .NOT.prefac )                           &
     &                          CALL ZLASET('Full',n,n,DCMPLX(ZERO),    &
     &                          DCMPLX(ZERO),Afac,lda)
                           CALL ZLASET('Full',n,Nrhs,DCMPLX(ZERO),      &
     &                                 DCMPLX(ZERO),X,lda)
!
!                       Equilibrate the matrix if FACT = 'F' and
!                       EQUED = 'R', 'C', or 'B'.
!
                           IF ( iequed>1 .AND. n>0 )                    &
     &                          CALL ZLAQGE(n,n,A,lda,S,S(n+1),rowcnd,  &
     &                          colcnd,amax,equed)
!
!                    Solve the system and compute the condition number
!                    and error bounds using ZGESVX.
!
                           SRNamt = 'ZGESVX'
                           CALL ZGESVX(fact,trans,n,Nrhs,A,lda,Afac,lda,&
     &                                 Iwork,equed,S,S(n+1),B,lda,X,lda,&
     &                                 rcond,Rwork,Rwork(Nrhs+1),Work,  &
     &                                 Rwork(2*Nrhs+1),info)
!
!                    Check the error code from ZGESVX.
!
                           IF ( info/=izero )                           &
     &                          CALL alaerh(path,'ZGESVX',info,izero,   &
     &                          fact//trans,n,n,-1,-1,Nrhs,imat,nfail,  &
     &                          nerrs,Nout)
!
!                    Compare RWORK(2*NRHS+1) from ZGESVX with the
!                    computed reciprocal pivot growth factor RPVGRW
!
                           IF ( info/=0 .AND. info<=n ) THEN
                              rpvgrw = ZLANTR('M','U','N',info,info,    &
     &                                 Afac,lda,rdum)
                              IF ( rpvgrw==ZERO ) THEN
                                 rpvgrw = ONE
                              ELSE
                                 rpvgrw = ZLANGE('M',n,info,A,lda,rdum) &
     &                              /rpvgrw
                              ENDIF
                           ELSE
                              rpvgrw = ZLANTR('M','U','N',n,n,Afac,lda, &
     &                                 rdum)
                              IF ( rpvgrw==ZERO ) THEN
                                 rpvgrw = ONE
                              ELSE
                                 rpvgrw = ZLANGE('M',n,n,A,lda,rdum)    &
     &                              /rpvgrw
                              ENDIF
                           ENDIF
                           result(7) = ABS(rpvgrw-Rwork(2*Nrhs+1))      &
     &                                 /MAX(Rwork(2*Nrhs+1),rpvgrw)     &
     &                                 /DLAMCH('E')
!
                           IF ( .NOT.prefac ) THEN
!
!                       Reconstruct matrix from factors and compute
!                       residual.
!
                              CALL zget01(n,n,A,lda,Afac,lda,Iwork,     &
     &                           Rwork(2*Nrhs+1),result(1))
                              k1 = 1
                           ELSE
                              k1 = 2
                           ENDIF
!
                           IF ( info==0 ) THEN
                              trfcon = .FALSE.
!
!                       Compute residual of the computed solution.
!
                              CALL ZLACPY('Full',n,Nrhs,Bsav,lda,Work,  &
     &                           lda)
                              CALL zget02(trans,n,n,Nrhs,Asav,lda,X,lda,&
     &                           Work,lda,Rwork(2*Nrhs+1),result(2))
!
!                       Check solution from generated exact solution.
!
                              IF ( nofact .OR.                          &
     &                             (prefac .AND. LSAME(equed,'N')) )    &
     &                             THEN
                                 CALL zget04(n,Nrhs,X,lda,Xact,lda,     &
     &                              rcondc,result(3))
                              ELSE
                                 IF ( itran==1 ) THEN
                                    roldc = roldo
                                 ELSE
                                    roldc = roldi
                                 ENDIF
                                 CALL zget04(n,Nrhs,X,lda,Xact,lda,     &
     &                              roldc,result(3))
                              ENDIF
!
!                       Check the error bounds from iterative
!                       refinement.
!
                              CALL zget07(trans,n,Nrhs,Asav,lda,B,lda,X,&
     &                           lda,Xact,lda,Rwork,.TRUE.,Rwork(Nrhs+1)&
     &                           ,result(4))
                           ELSE
                              trfcon = .TRUE.
                           ENDIF
!
!                    Compare RCOND from ZGESVX with the computed value
!                    in RCONDC.
!
                           result(6) = dget06(rcond,rcondc)
!
!                    Print information about the tests that did not pass
!                    the threshold.
!
                           IF ( .NOT.trfcon ) THEN
                              DO k = k1 , NTESTS
                                 IF ( result(k)>=Thresh ) THEN
                                    IF ( nfail==0 .AND. nerrs==0 )      &
     &                                 CALL aladhd(Nout,path)
                                    IF ( prefac ) THEN
                                       WRITE (Nout,FMT=99003) 'ZGESVX' ,&
     &                                    fact , trans , n , equed ,    &
     &                                    imat , k , result(k)
                                    ELSE
                                       WRITE (Nout,FMT=99002) 'ZGESVX' ,&
     &                                    fact , trans , n , imat , k , &
     &                                    result(k)
                                    ENDIF
                                    nfail = nfail + 1
                                 ENDIF
                              ENDDO
                              nrun = nrun + NTESTS - k1 + 1
                           ELSE
                              IF ( result(1)>=Thresh .AND. .NOT.prefac )&
     &                             THEN
                                 IF ( nfail==0 .AND. nerrs==0 )         &
     &                                CALL aladhd(Nout,path)
                                 IF ( prefac ) THEN
                                    WRITE (Nout,FMT=99003) 'ZGESVX' ,   &
     &                                 fact , trans , n , equed , imat ,&
     &                                 1 , result(1)
                                 ELSE
                                    WRITE (Nout,FMT=99002) 'ZGESVX' ,   &
     &                                 fact , trans , n , imat , 1 ,    &
     &                                 result(1)
                                 ENDIF
                                 nfail = nfail + 1
                                 nrun = nrun + 1
                              ENDIF
                              IF ( result(6)>=Thresh ) THEN
                                 IF ( nfail==0 .AND. nerrs==0 )         &
     &                                CALL aladhd(Nout,path)
                                 IF ( prefac ) THEN
                                    WRITE (Nout,FMT=99003) 'ZGESVX' ,   &
     &                                 fact , trans , n , equed , imat ,&
     &                                 6 , result(6)
                                 ELSE
                                    WRITE (Nout,FMT=99002) 'ZGESVX' ,   &
     &                                 fact , trans , n , imat , 6 ,    &
     &                                 result(6)
                                 ENDIF
                                 nfail = nfail + 1
                                 nrun = nrun + 1
                              ENDIF
                              IF ( result(7)>=Thresh ) THEN
                                 IF ( nfail==0 .AND. nerrs==0 )         &
     &                                CALL aladhd(Nout,path)
                                 IF ( prefac ) THEN
                                    WRITE (Nout,FMT=99003) 'ZGESVX' ,   &
     &                                 fact , trans , n , equed , imat ,&
     &                                 7 , result(7)
                                 ELSE
                                    WRITE (Nout,FMT=99002) 'ZGESVX' ,   &
     &                                 fact , trans , n , imat , 7 ,    &
     &                                 result(7)
                                 ENDIF
                                 nfail = nfail + 1
                                 nrun = nrun + 1
                              ENDIF
!
                           ENDIF
!
                        ENDDO
                     ENDDO
                  ENDDO
               ENDIF
            ENDIF
         ENDDO
      ENDDO
!
!     Print a summary of the results.
!
      CALL alasvm(path,Nout,nfail,nrun,nerrs)
!
99001 FORMAT (1X,A,', N =',I5,', type ',I2,', test(',I2,') =',G12.5)
99002 FORMAT (1X,A,', FACT=''',A1,''', TRANS=''',A1,''', N=',I5,        &
     &        ', type ',I2,', test(',I1,')=',G12.5)
99003 FORMAT (1X,A,', FACT=''',A1,''', TRANS=''',A1,''', N=',I5,        &
     &        ', EQUED=''',A1,''', type ',I2,', test(',I1,')=',G12.5)
!
!     End of zdrvge
!
      END SUBROUTINE zdrvge
