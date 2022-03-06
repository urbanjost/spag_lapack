!*==zchkgb.f90  processed by SPAG 7.51RB at 17:34 on  4 Mar 2022
!> \brief \b zchkgb
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE zchkgb( DOTYPE, NM, MVAL, NN, NVAL, NNB, NBVAL, NNS,
!                          NSVAL, THRESH, TSTERR, A, LA, AFAC, LAFAC, B,
!                          X, XACT, WORK, RWORK, IWORK, NOUT )
!
!       .. Scalar Arguments ..
!       LOGICAL            TSTERR
!       INTEGER            LA, LAFAC, NM, NN, NNB, NNS, NOUT
!       DOUBLE PRECISION   THRESH
!       ..
!       .. Array Arguments ..
!       LOGICAL            DOTYPE( * )
!       INTEGER            IWORK( * ), MVAL( * ), NBVAL( * ), NSVAL( * ),
!      $                   NVAL( * )
!       DOUBLE PRECISION   RWORK( * )
!       COMPLEX*16         A( * ), AFAC( * ), B( * ), WORK( * ), X( * ),
!      $                   XACT( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> zchkgb tests ZGBTRF, -TRS, -RFS, and -CON
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
!> \param[in] NM
!> \verbatim
!>          NM is INTEGER
!>          The number of values of M contained in the vector MVAL.
!> \endverbatim
!>
!> \param[in] MVAL
!> \verbatim
!>          MVAL is INTEGER array, dimension (NM)
!>          The values of the matrix row dimension M.
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
!> \param[in] NNB
!> \verbatim
!>          NNB is INTEGER
!>          The number of values of NB contained in the vector NBVAL.
!> \endverbatim
!>
!> \param[in] NBVAL
!> \verbatim
!>          NBVAL is INTEGER array, dimension (NBVAL)
!>          The values of the blocksize NB.
!> \endverbatim
!>
!> \param[in] NNS
!> \verbatim
!>          NNS is INTEGER
!>          The number of values of NRHS contained in the vector NSVAL.
!> \endverbatim
!>
!> \param[in] NSVAL
!> \verbatim
!>          NSVAL is INTEGER array, dimension (NNS)
!>          The values of the number of right hand sides NRHS.
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
!>          The length of the array A.  LA >= (KLMAX+KUMAX+1)*NMAX
!>          where KLMAX is the largest entry in the local array KLVAL,
!>                KUMAX is the largest entry in the local array KUVAL and
!>                NMAX is the largest entry in the input array NVAL.
!> \endverbatim
!>
!> \param[out] AFAC
!> \verbatim
!>          AFAC is COMPLEX*16 array, dimension (LAFAC)
!> \endverbatim
!>
!> \param[in] LAFAC
!> \verbatim
!>          LAFAC is INTEGER
!>          The length of the array AFAC. LAFAC >= (2*KLMAX+KUMAX+1)*NMAX
!>          where KLMAX is the largest entry in the local array KLVAL,
!>                KUMAX is the largest entry in the local array KUVAL and
!>                NMAX is the largest entry in the input array NVAL.
!> \endverbatim
!>
!> \param[out] B
!> \verbatim
!>          B is COMPLEX*16 array, dimension (NMAX*NSMAX)
!> \endverbatim
!>
!> \param[out] X
!> \verbatim
!>          X is COMPLEX*16 array, dimension (NMAX*NSMAX)
!> \endverbatim
!>
!> \param[out] XACT
!> \verbatim
!>          XACT is COMPLEX*16 array, dimension (NMAX*NSMAX)
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX*16 array, dimension
!>                      (NMAX*max(3,NSMAX,NMAX))
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is DOUBLE PRECISION array, dimension
!>                      (max(NMAX,2*NSMAX))
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
      SUBROUTINE zchkgb(Dotype,Nm,Mval,Nn,Nval,Nnb,Nbval,Nns,Nsval,     &
     &                  Thresh,Tsterr,A,La,Afac,Lafac,B,X,Xact,Work,    &
     &                  Rwork,Iwork,Nout)
      IMPLICIT NONE
!*--ZCHKGB195
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      LOGICAL Tsterr
      INTEGER La , Lafac , Nm , Nn , Nnb , Nns , Nout
      DOUBLE PRECISION Thresh
!     ..
!     .. Array Arguments ..
      LOGICAL Dotype(*)
      INTEGER Iwork(*) , Mval(*) , Nbval(*) , Nsval(*) , Nval(*)
      DOUBLE PRECISION Rwork(*)
      COMPLEX*16 A(*) , Afac(*) , B(*) , Work(*) , X(*) , Xact(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ONE , ZERO
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
      INTEGER NTYPES , NTESTS
      PARAMETER (NTYPES=8,NTESTS=7)
      INTEGER NBW , NTRAN
      PARAMETER (NBW=4,NTRAN=3)
!     ..
!     .. Local Scalars ..
      LOGICAL trfcon , zerot
      CHARACTER dist , norm , trans , type , xtype
      CHARACTER*3 path
      INTEGER i , i1 , i2 , ikl , iku , im , imat , in , inb , info ,   &
     &        ioff , irhs , itran , izero , j , k , kl , koff , ku ,    &
     &        lda , ldafac , ldb , m , mode , n , nb , nerrs , nfail ,  &
     &        nimat , nkl , nku , nrhs , nrun
      DOUBLE PRECISION ainvnm , anorm , anormi , anormo , cndnum ,      &
     &                 rcond , rcondc , rcondi , rcondo
!     ..
!     .. Local Arrays ..
      CHARACTER transs(NTRAN)
      INTEGER iseed(4) , iseedy(4) , klval(NBW) , kuval(NBW)
      DOUBLE PRECISION result(NTESTS)
!     ..
!     .. External Functions ..
      DOUBLE PRECISION ZLANGB , ZLANGE
EXTERNAL    ZLANGB,   ZLANGE
!     ..
!     .. External Subroutines ..
EXTERNAL    ZCOPY
EXTERNAL    ZGBCON,   ZGBRFS,   ZGBTRF
EXTERNAL    ZGBTRS,   ZLACPY,   ZLASET
EXTERNAL    ZLATMS
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC DCMPLX , MAX , MIN
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
      DATA iseedy/1988 , 1989 , 1990 , 1991/ , transs/'N' , 'T' , 'C'/
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
      IF ( Tsterr ) CALL zerrge(path,Nout)
      INFot = 0
!
!     Initialize the first value for the lower and upper bandwidths.
!
      klval(1) = 0
      kuval(1) = 0
!
!     Do for each value of M in MVAL
!
      DO im = 1 , Nm
         m = Mval(im)
!
!        Set values to use for the lower bandwidth.
!
         klval(2) = m + (m+1)/4
!
!        KLVAL( 2 ) = MAX( M-1, 0 )
!
         klval(3) = (3*m-1)/4
         klval(4) = (m+1)/4
!
!        Do for each value of N in NVAL
!
         DO in = 1 , Nn
            n = Nval(in)
            xtype = 'N'
!
!           Set values to use for the upper bandwidth.
!
            kuval(2) = n + (n+1)/4
!
!           KUVAL( 2 ) = MAX( N-1, 0 )
!
            kuval(3) = (3*n-1)/4
            kuval(4) = (n+1)/4
!
!           Set limits on the number of loop iterations.
!
            nkl = MIN(m+1,4)
            IF ( n==0 ) nkl = 2
            nku = MIN(n+1,4)
            IF ( m==0 ) nku = 2
            nimat = NTYPES
            IF ( m<=0 .OR. n<=0 ) nimat = 1
!
            DO ikl = 1 , nkl
!
!              Do for KL = 0, (5*M+1)/4, (3M-1)/4, and (M+1)/4. This
!              order makes it easier to skip redundant values for small
!              values of M.
!
               kl = klval(ikl)
               DO iku = 1 , nku
!
!                 Do for KU = 0, (5*N+1)/4, (3N-1)/4, and (N+1)/4. This
!                 order makes it easier to skip redundant values for
!                 small values of N.
!
                  ku = kuval(iku)
!
!                 Check that A and AFAC are big enough to generate this
!                 matrix.
!
                  lda = kl + ku + 1
                  ldafac = 2*kl + ku + 1
                  IF ( (lda*n)>La .OR. (ldafac*n)>Lafac ) THEN
                     IF ( nfail==0 .AND. nerrs==0 )                     &
     &                    CALL alahd(Nout,path)
                     IF ( n*(kl+ku+1)>La ) THEN
                        WRITE (Nout,FMT=99001) La , m , n , kl , ku ,   &
     &                         n*(kl+ku+1)
                        nerrs = nerrs + 1
                     ENDIF
                     IF ( n*(2*kl+ku+1)>Lafac ) THEN
                        WRITE (Nout,FMT=99002) Lafac , m , n , kl , ku ,&
     &                         n*(2*kl+ku+1)
                        nerrs = nerrs + 1
                     ENDIF
                     CYCLE
                  ENDIF
!
                  DO imat = 1 , nimat
!
!                    Do the tests only if DOTYPE( IMAT ) is true.
!
                     IF ( Dotype(imat) ) THEN
!
!                    Skip types 2, 3, or 4 if the matrix size is too
!                    small.
!
                        zerot = imat>=2 .AND. imat<=4
                        IF ( .NOT.(zerot .AND. n<imat-1) ) THEN
!
                           IF ( .NOT.zerot .OR. .NOT.Dotype(1) ) THEN
!
!                       Set up parameters with zlatb4 and generate a
!                       test matrix with ZLATMS.
!
                              CALL zlatb4(path,imat,m,n,type,kl,ku,     &
     &                           anorm,mode,cndnum,dist)
!
                              koff = MAX(1,ku+2-n)
                              DO i = 1 , koff - 1
                                 A(i) = ZERO
                              ENDDO
                              SRNamt = 'ZLATMS'
                              CALL ZLATMS(m,n,dist,iseed,type,Rwork,    &
     &                           mode,cndnum,anorm,kl,ku,'Z',A(koff),   &
     &                           lda,Work,info)
!
!                       Check the error code from ZLATMS.
!
                              IF ( info/=0 ) THEN
                                 CALL alaerh(path,'ZLATMS',info,0,' ',m,&
     &                              n,kl,ku,-1,imat,nfail,nerrs,Nout)
                                 CYCLE
                              ENDIF
                           ELSEIF ( izero>0 ) THEN
!
!                       Use the same matrix for types 3 and 4 as for
!                       type 2 by copying back the zeroed out column.
!
                              CALL ZCOPY(i2-i1+1,B,1,A(ioff+i1),1)
                           ENDIF
!
!                    For types 2, 3, and 4, zero one or more columns of
!                    the matrix to test that INFO is returned correctly.
!
                           izero = 0
                           IF ( zerot ) THEN
                              IF ( imat==2 ) THEN
                                 izero = 1
                              ELSEIF ( imat==3 ) THEN
                                 izero = MIN(m,n)
                              ELSE
                                 izero = MIN(m,n)/2 + 1
                              ENDIF
                              ioff = (izero-1)*lda
                              IF ( imat<4 ) THEN
!
!                          Store the column to be zeroed out in B.
!
                                 i1 = MAX(1,ku+2-izero)
                                 i2 = MIN(kl+ku+1,ku+1+(m-izero))
                                 CALL ZCOPY(i2-i1+1,A(ioff+i1),1,B,1)
!
                                 DO i = i1 , i2
                                    A(ioff+i) = ZERO
                                 ENDDO
                              ELSE
                                 DO j = izero , n
                                    DO i = MAX(1,ku+2-j) ,              &
     &                                 MIN(kl+ku+1,ku+1+(m-j))
                                       A(ioff+i) = ZERO
                                    ENDDO
                                    ioff = ioff + lda
                                 ENDDO
                              ENDIF
                           ENDIF
!
!                    These lines, if used in place of the calls in the
!                    loop over INB, cause the code to bomb on a Sun
!                    SPARCstation.
!
!                     ANORMO = ZLANGB( 'O', N, KL, KU, A, LDA, RWORK )
!                     ANORMI = ZLANGB( 'I', N, KL, KU, A, LDA, RWORK )
!
!                    Do for each blocksize in NBVAL
!
                           DO inb = 1 , Nnb
                              nb = Nbval(inb)
                              CALL xlaenv(1,nb)
!
!                       Compute the LU factorization of the band matrix.
!
                              IF ( m>0 .AND. n>0 )                      &
     &                             CALL ZLACPY('Full',kl+ku+1,n,A,lda,  &
     &                             Afac(kl+1),ldafac)
                              SRNamt = 'ZGBTRF'
                              CALL ZGBTRF(m,n,kl,ku,Afac,ldafac,Iwork,  &
     &                           info)
!
!                       Check error code from ZGBTRF.
!
                              IF ( info/=izero )                        &
     &                              CALL alaerh(path,'ZGBTRF',info,     &
     &                             izero,' ',m,n,kl,ku,nb,imat,nfail,   &
     &                             nerrs,Nout)
                              trfcon = .FALSE.
!
!+    TEST 1
!                       Reconstruct matrix from factors and compute
!                       residual.
!
                              CALL zgbt01(m,n,kl,ku,A,lda,Afac,ldafac,  &
     &                           Iwork,Work,result(1))
!
!                       Print information about the tests so far that
!                       did not pass the threshold.
!
                              IF ( result(1)>=Thresh ) THEN
                                 IF ( nfail==0 .AND. nerrs==0 )         &
     &                                CALL alahd(Nout,path)
                                 WRITE (Nout,FMT=99003) m , n , kl ,    &
     &                                  ku , nb , imat , 1 , result(1)
                                 nfail = nfail + 1
                              ENDIF
                              nrun = nrun + 1
!
!                       Skip the remaining tests if this is not the
!                       first block size or if M .ne. N.
!
                              IF ( inb<=1 .AND. m==n ) THEN
!
                                 anormo = ZLANGB('O',n,kl,ku,A,lda,     &
     &                              Rwork)
                                 anormi = ZLANGB('I',n,kl,ku,A,lda,     &
     &                              Rwork)
!
                                 IF ( info==0 ) THEN
!
!                          Form the inverse of A so we can get a good
!                          estimate of CNDNUM = norm(A) * norm(inv(A)).
!
                                    ldb = MAX(1,n)
                                    CALL ZLASET('Full',n,n,DCMPLX(ZERO),&
     &                                 DCMPLX(ONE),Work,ldb)
                                    SRNamt = 'ZGBTRS'
                                    CALL ZGBTRS('No transpose',n,kl,ku, &
     &                                 n,Afac,ldafac,Iwork,Work,ldb,    &
     &                                 info)
!
!                          Compute the 1-norm condition number of A.
!
                                    ainvnm = ZLANGE('O',n,n,Work,ldb,   &
     &                                 Rwork)
                                    IF ( anormo<=ZERO .OR.              &
     &                                 ainvnm<=ZERO ) THEN
                                       rcondo = ONE
                                    ELSE
                                       rcondo = (ONE/anormo)/ainvnm
                                    ENDIF
!
!                          Compute the infinity-norm condition number of
!                          A.
!
                                    ainvnm = ZLANGE('I',n,n,Work,ldb,   &
     &                                 Rwork)
                                    IF ( anormi<=ZERO .OR.              &
     &                                 ainvnm<=ZERO ) THEN
                                       rcondi = ONE
                                    ELSE
                                       rcondi = (ONE/anormi)/ainvnm
                                    ENDIF
                                 ELSE
!
!                          Do only the condition estimate if INFO.NE.0.
!
                                    trfcon = .TRUE.
                                    rcondo = ZERO
                                    rcondi = ZERO
                                 ENDIF
!
!                       Skip the solve tests if the matrix is singular.
!
                                 IF ( .NOT.(trfcon) ) THEN
!
                                    DO irhs = 1 , Nns
                                       nrhs = Nsval(irhs)
                                       xtype = 'N'
!
                                       DO itran = 1 , NTRAN
                                         trans = transs(itran)
                                         IF ( itran==1 ) THEN
                                         rcondc = rcondo
                                         norm = 'O'
                                         ELSE
                                         rcondc = rcondi
                                         norm = 'I'
                                         ENDIF
!
!+    TEST 2:
!                             Solve and compute residual for A * X = B.
!
                                         SRNamt = 'zlarhs'
                                         CALL zlarhs(path,xtype,' ',    &
     &                                      trans,n,n,kl,ku,nrhs,A,lda, &
     &                                      Xact,ldb,B,ldb,iseed,info)
                                         xtype = 'C'
                                         CALL ZLACPY('Full',n,nrhs,B,   &
     &                                      ldb,X,ldb)
!
                                         SRNamt = 'ZGBTRS'
                                         CALL ZGBTRS(trans,n,kl,ku,nrhs,&
     &                                      Afac,ldafac,Iwork,X,ldb,    &
     &                                      info)
!
!                             Check error code from ZGBTRS.
!
                                         IF ( info/=0 )                 &
     &                                      CALL alaerh(path,'ZGBTRS',  &
     &                                      info,0,trans,n,n,kl,ku,-1,  &
     &                                      imat,nfail,nerrs,Nout)
!
                                         CALL ZLACPY('Full',n,nrhs,B,   &
     &                                      ldb,Work,ldb)
                                         CALL zgbt02(trans,m,n,kl,ku,   &
     &                                      nrhs,A,lda,X,ldb,Work,ldb,  &
     &                                      result(2))
!
!+    TEST 3:
!                             Check solution from generated exact
!                             solution.
!
                                         CALL zget04(n,nrhs,X,ldb,Xact, &
     &                                      ldb,rcondc,result(3))
!
!+    TESTS 4, 5, 6:
!                             Use iterative refinement to improve the
!                             solution.
!
                                         SRNamt = 'ZGBRFS'
                                         CALL ZGBRFS(trans,n,kl,ku,nrhs,&
     &                                      A,lda,Afac,ldafac,Iwork,B,  &
     &                                      ldb,X,ldb,Rwork,            &
     &                                      Rwork(nrhs+1),Work,         &
     &                                      Rwork(2*nrhs+1),info)
!
!                             Check error code from ZGBRFS.
!
                                         IF ( info/=0 )                 &
     &                                      CALL alaerh(path,'ZGBRFS',  &
     &                                      info,0,trans,n,n,kl,ku,nrhs,&
     &                                      imat,nfail,nerrs,Nout)
!
                                         CALL zget04(n,nrhs,X,ldb,Xact, &
     &                                      ldb,rcondc,result(4))
                                         CALL zgbt05(trans,n,kl,ku,nrhs,&
     &                                      A,lda,B,ldb,X,ldb,Xact,ldb, &
     &                                      Rwork,Rwork(nrhs+1),        &
     &                                      result(5))
!
!                             Print information about the tests that did
!                             not pass the threshold.
!
                                         DO k = 2 , 6
                                         IF ( result(k)>=Thresh ) THEN
                                         IF ( nfail==0 .AND. nerrs==0 ) &
     &                                      CALL alahd(Nout,path)
                                         WRITE (Nout,FMT=99004) trans , &
     &                                      n , kl , ku , nrhs , imat , &
     &                                      k , result(k)
                                         nfail = nfail + 1
                                         ENDIF
                                         ENDDO
                                         nrun = nrun + 5
                                       ENDDO
                                    ENDDO
                                 ENDIF
!
!+    TEST 7:
!                          Get an estimate of RCOND = 1/CNDNUM.
!
                                 DO itran = 1 , 2
                                    IF ( itran==1 ) THEN
                                       anorm = anormo
                                       rcondc = rcondo
                                       norm = 'O'
                                    ELSE
                                       anorm = anormi
                                       rcondc = rcondi
                                       norm = 'I'
                                    ENDIF
                                    SRNamt = 'ZGBCON'
                                    CALL ZGBCON(norm,n,kl,ku,Afac,      &
     &                                 ldafac,Iwork,anorm,rcond,Work,   &
     &                                 Rwork,info)
!
!                             Check error code from ZGBCON.
!
                                    IF ( info/=0 )                      &
     &                                  CALL alaerh(path,'ZGBCON',info, &
     &                                 0,norm,n,n,kl,ku,-1,imat,nfail,  &
     &                                 nerrs,Nout)
!
                                    result(7) = dget06(rcond,rcondc)
!
!                          Print information about the tests that did
!                          not pass the threshold.
!
                                    IF ( result(7)>=Thresh ) THEN
                                       IF ( nfail==0 .AND. nerrs==0 )   &
     &                                    CALL alahd(Nout,path)
                                       WRITE (Nout,FMT=99005) norm , n ,&
     &                                    kl , ku , imat , 7 , result(7)
                                       nfail = nfail + 1
                                    ENDIF
                                    nrun = nrun + 1
                                 ENDDO
                              ENDIF
                           ENDDO
                        ENDIF
                     ENDIF
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
!
!     Print a summary of the results.
!
      CALL alasum(path,Nout,nfail,nrun,nerrs)
!
99001 FORMAT (' *** In zchkgb, LA=',I5,' is too small for M=',I5,', N=',&
     &        I5,', KL=',I4,', KU=',I4,/' ==> Increase LA to at least ',&
     &        I5)
99002 FORMAT (' *** In zchkgb, LAFAC=',I5,' is too small for M=',I5,    &
     &        ', N=',I5,', KL=',I4,', KU=',I4,                          &
     &        /' ==> Increase LAFAC to at least ',I5)
99003 FORMAT (' M =',I5,', N =',I5,', KL=',I5,', KU=',I5,', NB =',I4,   &
     &        ', type ',I1,', test(',I1,')=',G12.5)
99004 FORMAT (' TRANS=''',A1,''', N=',I5,', KL=',I5,', KU=',I5,         &
     &        ', NRHS=',I3,', type ',I1,', test(',I1,')=',G12.5)
99005 FORMAT (' NORM =''',A1,''', N=',I5,', KL=',I5,', KU=',I5,',',10X, &
     &        ' type ',I1,', test(',I1,')=',G12.5)
!
!
!     End of zchkgb
!
      END SUBROUTINE zchkgb
