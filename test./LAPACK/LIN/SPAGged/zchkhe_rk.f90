!*==zchkhe_rk.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b ZCHKHE_RK
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZCHKHE_RK( DOTYPE, NN, NVAL, NNB, NBVAL, NNS, NSVAL,
!                             THRESH, TSTERR, NMAX, A, AFAC, AINV, B, X,
!                             XACT, WORK, RWORK, IWORK, NOUT )
!
!       .. Scalar Arguments ..
!       LOGICAL            TSTERR
!       INTEGER            NMAX, NN, NNB, NNS, NOUT
!       DOUBLE PRECISION   THRESH
!       ..
!       .. Array Arguments ..
!       LOGICAL            DOTYPE( * )
!       INTEGER            IWORK( * ), NBVAL( * ), NSVAL( * ), NVAL( * )
!       DOUBLE PRECISION   RWORK( * )
!       COMPLEX*16         A( * ), AFAC( * ), AINV( * ), B( * ), E( * ),
!      $                   WORK( * ), X( * ), XACT( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZCHKHE_RK tests ZHETRF_RK, -TRI_3, -TRS_3,
!> and -CON_3.
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
!>          The values of the matrix dimension N.
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
!> \param[in] NMAX
!> \verbatim
!>          NMAX is INTEGER
!>          The maximum value permitted for N, used in dimensioning the
!>          work arrays.
!> \endverbatim
!>
!> \param[out] A
!> \verbatim
!>          A is CCOMPLEX*16 array, dimension (NMAX*NMAX)
!> \endverbatim
!>
!> \param[out] AFAC
!> \verbatim
!>          AFAC is COMPLEX*16 array, dimension (NMAX*NMAX)
!> \endverbatim
!>
!> \param[out] E
!> \verbatim
!>          E is COMPLEX*16 array, dimension (NMAX)
!> \endverbatim
!>
!> \param[out] AINV
!> \verbatim
!>          AINV is COMPLEX*16 array, dimension (NMAX*NMAX)
!> \endverbatim
!>
!> \param[out] B
!> \verbatim
!>          B is CCOMPLEX*16 array, dimension (NMAX*NSMAX)
!>          where NSMAX is the largest entry in NSVAL.
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
!>          WORK is COMPLEX*16 array, dimension (NMAX*max(3,NSMAX))
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is DOUBLE PRECISION array, dimension (max(NMAX,2*NSMAX)
!> \endverbatim
!>
!> \param[out] IWORK
!> \verbatim
!>          IWORK is INTEGER array, dimension (2*NMAX)
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
      SUBROUTINE ZCHKHE_RK(Dotype,Nn,Nval,Nnb,Nbval,Nns,Nsval,Thresh,   &
     &                     Tsterr,Nmax,A,Afac,E,Ainv,B,X,Xact,Work,     &
     &                     Rwork,Iwork,Nout)
      IMPLICIT NONE
!*--ZCHKHE_RK181
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      LOGICAL Tsterr
      INTEGER Nmax , Nn , Nnb , Nns , Nout
      DOUBLE PRECISION Thresh
!     ..
!     .. Array Arguments ..
      LOGICAL Dotype(*)
      INTEGER Iwork(*) , Nbval(*) , Nsval(*) , Nval(*)
      DOUBLE PRECISION Rwork(*)
      COMPLEX*16 A(*) , Afac(*) , Ainv(*) , B(*) , E(*) , Work(*) ,     &
     &           X(*) , Xact(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO , ONE
      PARAMETER (ZERO=0.0D+0,ONE=1.0D+0)
      DOUBLE PRECISION ONEHALF
      PARAMETER (ONEHALF=0.5D+0)
      DOUBLE PRECISION EIGHT , SEVTEN
      PARAMETER (EIGHT=8.0D+0,SEVTEN=17.0D+0)
      COMPLEX*16 CZERO
      PARAMETER (CZERO=(0.0D+0,0.0D+0))
      INTEGER NTYPES
      PARAMETER (NTYPES=10)
      INTEGER NTESTS
      PARAMETER (NTESTS=7)
!     ..
!     .. Local Scalars ..
      LOGICAL trfcon , zerot
      CHARACTER dist , type , uplo , xtype
      CHARACTER*3 path , matpath
      INTEGER i , i1 , i2 , imat , in , inb , info , ioff , irhs ,      &
     &        iuplo , izero , j , k , kl , ku , lda , lwork , mode , n ,&
     &        nb , nerrs , nfail , nimat , nrhs , nrun , nt
      DOUBLE PRECISION alpha , anorm , cndnum , const , sing_max ,      &
     &                 sing_min , rcond , rcondc , dtemp
!     ..
!     .. Local Arrays ..
      CHARACTER uplos(2)
      INTEGER iseed(4) , iseedy(4)
      DOUBLE PRECISION result(NTESTS)
      COMPLEX*16 block(2,2) , zdummy(1)
!     ..
!     .. External Functions ..
      DOUBLE PRECISION DGET06 , ZLANGE , ZLANHE
      EXTERNAL DGET06 , ZLANGE , ZLANHE
!     ..
!     .. External Subroutines ..
      EXTERNAL ALAERH , ALAHD , ALASUM , ZERRHE , ZGESVD , ZGET04 ,     &
     &         ZLACPY , ZLARHS , ZLATB4 , ZLATMS , ZPOT02 , ZPOT03 ,    &
     &         ZHECON_3 , ZHET01_3 , ZHETRF_RK , ZHETRI_3 , ZHETRS_3 ,  &
     &         XLAENV
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC DCONJG , MAX , MIN , SQRT
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
      DATA uplos/'U' , 'L'/
!     ..
!     .. Executable Statements ..
!
!     Initialize constants and the random number seed.
!
      alpha = (ONE+SQRT(SEVTEN))/EIGHT
!
!     Test path
!
      path(1:1) = 'Zomplex precision'
      path(2:3) = 'HK'
!
!     Path to generate matrices
!
      matpath(1:1) = 'Zomplex precision'
      matpath(2:3) = 'HE'
!
      nrun = 0
      nfail = 0
      nerrs = 0
      DO i = 1 , 4
         iseed(i) = iseedy(i)
      ENDDO
!
!     Test the error exits
!
      IF ( Tsterr ) CALL ZERRHE(path,Nout)
      INFot = 0
!
!     Set the minimum block size for which the block routine should
!     be used, which will be later returned by ILAENV
!
      CALL XLAENV(2,2)
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
         izero = 0
!
!        Do for each value of matrix type IMAT
!
         DO imat = 1 , nimat
!
!           Do the tests only if DOTYPE( IMAT ) is true.
!
            IF ( Dotype(imat) ) THEN
!
!           Skip types 3, 4, 5, or 6 if the matrix size is too small.
!
               zerot = imat>=3 .AND. imat<=6
               IF ( .NOT.(zerot .AND. n<imat-2) ) THEN
!
!           Do first for UPLO = 'U', then for UPLO = 'L'
!
                  DO iuplo = 1 , 2
                     uplo = uplos(iuplo)
!
!                 Begin generate the test matrix A.
!
!                 Set up parameters with ZLATB4 for the matrix generator
!                 based on the type of matrix to be generated.
!
                     CALL ZLATB4(matpath,imat,n,n,type,kl,ku,anorm,mode,&
     &                           cndnum,dist)
!
!                 Generate a matrix with ZLATMS.
!
                     SRNamt = 'ZLATMS'
                     CALL ZLATMS(n,n,dist,iseed,type,Rwork,mode,cndnum, &
     &                           anorm,kl,ku,uplo,A,lda,Work,info)
!
!                 Check error code from ZLATMS and handle error.
!
                     IF ( info/=0 ) THEN
                        CALL ALAERH(path,'ZLATMS',info,0,uplo,n,n,-1,-1,&
     &                              -1,imat,nfail,nerrs,Nout)
!
!                    Skip all tests for this generated matrix
!
                        CYCLE
                     ENDIF
!
!                 For matrix types 3-6, zero one or more rows and
!                 columns of the matrix to test that INFO is returned
!                 correctly.
!
                     IF ( zerot ) THEN
                        IF ( imat==3 ) THEN
                           izero = 1
                        ELSEIF ( imat==4 ) THEN
                           izero = n
                        ELSE
                           izero = n/2 + 1
                        ENDIF
!
                        IF ( imat<6 ) THEN
!
!                       Set row and column IZERO to zero.
!
                           IF ( iuplo==1 ) THEN
                              ioff = (izero-1)*lda
                              DO i = 1 , izero - 1
                                 A(ioff+i) = CZERO
                              ENDDO
                              ioff = ioff + izero
                              DO i = izero , n
                                 A(ioff) = CZERO
                                 ioff = ioff + lda
                              ENDDO
                           ELSE
                              ioff = izero
                              DO i = 1 , izero - 1
                                 A(ioff) = CZERO
                                 ioff = ioff + lda
                              ENDDO
                              ioff = ioff - izero
                              DO i = izero , n
                                 A(ioff+i) = CZERO
                              ENDDO
                           ENDIF
                        ELSEIF ( iuplo==1 ) THEN
!
!                          Set the first IZERO rows and columns to zero.
!
                           ioff = 0
                           DO j = 1 , n
                              i2 = MIN(j,izero)
                              DO i = 1 , i2
                                 A(ioff+i) = CZERO
                              ENDDO
                              ioff = ioff + lda
                           ENDDO
                        ELSE
!
!                          Set the last IZERO rows and columns to zero.
!
                           ioff = 0
                           DO j = 1 , n
                              i1 = MAX(j,izero)
                              DO i = i1 , n
                                 A(ioff+i) = CZERO
                              ENDDO
                              ioff = ioff + lda
                           ENDDO
                        ENDIF
                     ELSE
                        izero = 0
                     ENDIF
!
!                 End generate the test matrix A.
!
!
!              Do for each value of NB in NBVAL
!
                     DO inb = 1 , Nnb
!
!                 Set the optimal blocksize, which will be later
!                 returned by ILAENV.
!
                        nb = Nbval(inb)
                        CALL XLAENV(1,nb)
!
!                 Copy the test matrix A into matrix AFAC which
!                 will be factorized in place. This is needed to
!                 preserve the test matrix A for subsequent tests.
!
                        CALL ZLACPY(uplo,n,n,A,lda,Afac,lda)
!
!                 Compute the L*D*L**T or U*D*U**T factorization of the
!                 matrix. IWORK stores details of the interchanges and
!                 the block structure of D. AINV is a work array for
!                 block factorization, LWORK is the length of AINV.
!
                        lwork = MAX(2,nb)*lda
                        SRNamt = 'ZHETRF_RK'
                        CALL ZHETRF_RK(uplo,n,Afac,lda,E,Iwork,Ainv,    &
     &                                 lwork,info)
!
!                 Adjust the expected value of INFO to account for
!                 pivoting.
!
                        k = izero
                        IF ( k>0 ) THEN
                           DO
                              IF ( Iwork(k)<0 ) THEN
                                 IF ( Iwork(k)/=-k ) THEN
                                    k = -Iwork(k)
                                    CYCLE
                                 ENDIF
                              ELSEIF ( Iwork(k)/=k ) THEN
                                 k = Iwork(k)
                                 CYCLE
                              ENDIF
                              EXIT
                           ENDDO
                        ENDIF
!
!                 Check error code from ZHETRF_RK and handle error.
!
                        IF ( info/=k )                                  &
     &                       CALL ALAERH(path,'ZHETRF_RK',info,k,uplo,n,&
     &                       n,-1,-1,nb,imat,nfail,nerrs,Nout)
!
!                 Set the condition estimate flag if the INFO is not 0.
!
                        IF ( info/=0 ) THEN
                           trfcon = .TRUE.
                        ELSE
                           trfcon = .FALSE.
                        ENDIF
!
!+    TEST 1
!                 Reconstruct matrix from factors and compute residual.
!
                        CALL ZHET01_3(uplo,n,A,lda,Afac,lda,E,Iwork,    &
     &                                Ainv,lda,Rwork,result(1))
                        nt = 1
!
!+    TEST 2
!                 Form the inverse and compute the residual,
!                 if the factorization was competed without INFO > 0
!                 (i.e. there is no zero rows and columns).
!                 Do it only for the first block size.
!
                        IF ( inb==1 .AND. .NOT.trfcon ) THEN
                           CALL ZLACPY(uplo,n,n,Afac,lda,Ainv,lda)
                           SRNamt = 'ZHETRI_3'
!
!                    Another reason that we need to compute the inverse
!                    is that ZPOT03 produces RCONDC which is used later
!                    in TEST6 and TEST7.
!
                           lwork = (n+nb+1)*(nb+3)
                           CALL ZHETRI_3(uplo,n,Ainv,lda,E,Iwork,Work,  &
     &                        lwork,info)
!
!                    Check error code from ZHETRI_3 and handle error.
!
                           IF ( info/=0 )                               &
     &                          CALL ALAERH(path,'ZHETRI_3',info,-1,    &
     &                          uplo,n,n,-1,-1,-1,imat,nfail,nerrs,Nout)
!
!                    Compute the residual for a Hermitian matrix times
!                    its inverse.
!
                           CALL ZPOT03(uplo,n,A,lda,Ainv,lda,Work,lda,  &
     &                                 Rwork,rcondc,result(2))
                           nt = 2
                        ENDIF
!
!                 Print information about the tests that did not pass
!                 the threshold.
!
                        DO k = 1 , nt
                           IF ( result(k)>=Thresh ) THEN
                              IF ( nfail==0 .AND. nerrs==0 )            &
     &                             CALL ALAHD(Nout,path)
                              WRITE (Nout,FMT=99001) uplo , n , nb ,    &
     &                               imat , k , result(k)
                              nfail = nfail + 1
                           ENDIF
                        ENDDO
                        nrun = nrun + nt
!
!+    TEST 3
!                 Compute largest element in U or L
!
                        result(3) = ZERO
                        dtemp = ZERO
!
                        const = ((alpha**2-ONE)/(alpha**2-ONEHALF))     &
     &                          /(ONE-alpha)
!
                        IF ( iuplo==1 ) THEN
!
!                 Compute largest element in U
!
                           k = n
                           DO WHILE ( k>1 )
!
                              IF ( Iwork(k)>ZERO ) THEN
!
!                       Get max absolute value from elements
!                       in column k in U
!
                                 dtemp = ZLANGE('M',k-1,1,              &
     &                              Afac((k-1)*lda+1),lda,Rwork)
                              ELSE
!
!                       Get max absolute value from elements
!                       in columns k and k-1 in U
!
                                 dtemp = ZLANGE('M',k-2,2,              &
     &                              Afac((k-2)*lda+1),lda,Rwork)
                                 k = k - 1
!
                              ENDIF
!
!                    DTEMP should be bounded by CONST
!
                              dtemp = dtemp - const + Thresh
                              IF ( dtemp>result(3) ) result(3) = dtemp
!
!
                              k = k - 1
                           ENDDO
!
                        ELSE
!
!                 Compute largest element in L
!
                           k = 1
                           DO WHILE ( k<n )
!
                              IF ( Iwork(k)>ZERO ) THEN
!
!                       Get max absolute value from elements
!                       in column k in L
!
                                 dtemp = ZLANGE('M',n-k,1,              &
     &                              Afac((k-1)*lda+k+1),lda,Rwork)
                              ELSE
!
!                       Get max absolute value from elements
!                       in columns k and k+1 in L
!
                                 dtemp = ZLANGE('M',n-k-1,2,            &
     &                              Afac((k-1)*lda+k+2),lda,Rwork)
                                 k = k + 1
!
                              ENDIF
!
!                    DTEMP should be bounded by CONST
!
                              dtemp = dtemp - const + Thresh
                              IF ( dtemp>result(3) ) result(3) = dtemp
!
!
                              k = k + 1
                           ENDDO
                        ENDIF
!
!
!+    TEST 4
!                 Compute largest 2-Norm (condition number)
!                 of 2-by-2 diag blocks
!
                        result(4) = ZERO
                        dtemp = ZERO
!
                        const = ((alpha**2-ONE)/(alpha**2-ONEHALF))     &
     &                          *((ONE+alpha)/(ONE-alpha))
                        CALL ZLACPY(uplo,n,n,Afac,lda,Ainv,lda)
!
                        IF ( iuplo==1 ) THEN
!
!                    Loop backward for UPLO = 'U'
!
                           k = n
                           DO WHILE ( k>1 )
!
                              IF ( Iwork(k)<ZERO ) THEN
!
!                       Get the two singular values
!                       (real and non-negative) of a 2-by-2 block,
!                       store them in RWORK array
!
                                 block(1,1) = Afac((k-2)*lda+k-1)
                                 block(1,2) = E(k)
                                 block(2,1) = DCONJG(block(1,2))
                                 block(2,2) = Afac((k-1)*lda+k)
!
                                 CALL ZGESVD('N','N',2,2,block,2,Rwork, &
     &                              zdummy,1,zdummy,1,Work,6,Rwork(3),  &
     &                              info)
!
!
                                 sing_max = Rwork(1)
                                 sing_min = Rwork(2)
!
                                 dtemp = sing_max/sing_min
!
!                       DTEMP should be bounded by CONST
!
                                 dtemp = dtemp - const + Thresh
                                 IF ( dtemp>result(4) ) result(4)       &
     &                                = dtemp
                                 k = k - 1
!
                              ENDIF
!
!
                              k = k - 1
                           ENDDO
!
                        ELSE
!
!                    Loop forward for UPLO = 'L'
!
                           k = 1
                           DO WHILE ( k<n )
!
                              IF ( Iwork(k)<ZERO ) THEN
!
!                       Get the two singular values
!                       (real and non-negative) of a 2-by-2 block,
!                       store them in RWORK array
!
                                 block(1,1) = Afac((k-1)*lda+k)
                                 block(2,1) = E(k)
                                 block(1,2) = DCONJG(block(2,1))
                                 block(2,2) = Afac(k*lda+k+1)
!
                                 CALL ZGESVD('N','N',2,2,block,2,Rwork, &
     &                              zdummy,1,zdummy,1,Work,6,Rwork(3),  &
     &                              info)
!
                                 sing_max = Rwork(1)
                                 sing_min = Rwork(2)
!
                                 dtemp = sing_max/sing_min
!
!                       DTEMP should be bounded by CONST
!
                                 dtemp = dtemp - const + Thresh
                                 IF ( dtemp>result(4) ) result(4)       &
     &                                = dtemp
                                 k = k + 1
!
                              ENDIF
!
!
                              k = k + 1
                           ENDDO
                        ENDIF
!
!                 Print information about the tests that did not pass
!                 the threshold.
!
                        DO k = 3 , 4
                           IF ( result(k)>=Thresh ) THEN
                              IF ( nfail==0 .AND. nerrs==0 )            &
     &                             CALL ALAHD(Nout,path)
                              WRITE (Nout,FMT=99001) uplo , n , nb ,    &
     &                               imat , k , result(k)
                              nfail = nfail + 1
                           ENDIF
                        ENDDO
                        nrun = nrun + 2
!
!                 Skip the other tests if this is not the first block
!                 size.
!
                        IF ( inb<=1 ) THEN
!
!                 Do only the condition estimate if INFO is not 0.
!
                           IF ( trfcon ) THEN
                              rcondc = ZERO
                              GOTO 2
                           ENDIF
!
!                 Do for each value of NRHS in NSVAL.
!
                           DO irhs = 1 , Nns
                              nrhs = Nsval(irhs)
!
!                    Begin loop over NRHS values
!
!
!+    TEST 5 ( Using TRS_3)
!                 Solve and compute residual for  A * X = B.
!
!                    Choose a set of NRHS random solution vectors
!                    stored in XACT and set up the right hand side B
!
                              SRNamt = 'ZLARHS'
                              CALL ZLARHS(matpath,xtype,uplo,' ',n,n,kl,&
     &                           ku,nrhs,A,lda,Xact,lda,B,lda,iseed,    &
     &                           info)
                              CALL ZLACPY('Full',n,nrhs,B,lda,X,lda)
!
                              SRNamt = 'ZHETRS_3'
                              CALL ZHETRS_3(uplo,n,nrhs,Afac,lda,E,     &
     &                           Iwork,X,lda,info)
!
!                    Check error code from ZHETRS_3 and handle error.
!
                              IF ( info/=0 )                            &
     &                             CALL ALAERH(path,'ZHETRS_3',info,0,  &
     &                             uplo,n,n,-1,-1,nrhs,imat,nfail,nerrs,&
     &                             Nout)
!
                              CALL ZLACPY('Full',n,nrhs,B,lda,Work,lda)
!
!                    Compute the residual for the solution
!
                              CALL ZPOT02(uplo,n,nrhs,A,lda,X,lda,Work, &
     &                           lda,Rwork,result(5))
!
!+    TEST 6
!                 Check solution from generated exact solution.
!
                              CALL ZGET04(n,nrhs,X,lda,Xact,lda,rcondc, &
     &                           result(6))
!
!                    Print information about the tests that did not pass
!                    the threshold.
!
                              DO k = 5 , 6
                                 IF ( result(k)>=Thresh ) THEN
                                    IF ( nfail==0 .AND. nerrs==0 )      &
     &                                 CALL ALAHD(Nout,path)
                                    WRITE (Nout,FMT=99002) uplo , n ,   &
     &                                 nrhs , imat , k , result(k)
                                    nfail = nfail + 1
                                 ENDIF
                              ENDDO
                              nrun = nrun + 2
!
!                 End do for each value of NRHS in NSVAL.
!
                           ENDDO
!
!+    TEST 7
!                 Get an estimate of RCOND = 1/CNDNUM.
!
 2                         anorm = ZLANHE('1',uplo,n,A,lda,Rwork)
                           SRNamt = 'ZHECON_3'
                           CALL ZHECON_3(uplo,n,Afac,lda,E,Iwork,anorm, &
     &                        rcond,Work,info)
!
!                 Check error code from ZHECON_3 and handle error.
!
                           IF ( info/=0 )                               &
     &                          CALL ALAERH(path,'ZHECON_3',info,0,uplo,&
     &                          n,n,-1,-1,-1,imat,nfail,nerrs,Nout)
!
!                 Compute the test ratio to compare values of RCOND
!
                           result(7) = DGET06(rcond,rcondc)
!
!                 Print information about the tests that did not pass
!                 the threshold.
!
                           IF ( result(7)>=Thresh ) THEN
                              IF ( nfail==0 .AND. nerrs==0 )            &
     &                             CALL ALAHD(Nout,path)
                              WRITE (Nout,FMT=99003) uplo , n , imat ,  &
     &                               7 , result(7)
                              nfail = nfail + 1
                           ENDIF
                           nrun = nrun + 1
                        ENDIF
                     ENDDO
!
                  ENDDO
               ENDIF
            ENDIF
         ENDDO
      ENDDO
!
!     Print a summary of the results.
!
      CALL ALASUM(path,Nout,nfail,nrun,nerrs)
!
99001 FORMAT (' UPLO = ''',A1,''', N =',I5,', NB =',I4,', type ',I2,    &
     &        ', test ',I2,', ratio =',G12.5)
99002 FORMAT (' UPLO = ''',A1,''', N =',I5,', NRHS=',I3,', type ',I2,   &
     &        ', test ',I2,', ratio =',G12.5)
99003 FORMAT (' UPLO = ''',A1,''', N =',I5,',',10X,' type ',I2,         &
     &        ', test ',I2,', ratio =',G12.5)
!
!     End of ZCHKHE_RK
!
      END SUBROUTINE ZCHKHE_RK
