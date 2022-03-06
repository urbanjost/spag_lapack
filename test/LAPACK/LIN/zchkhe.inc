!*==zchkhe.f90  processed by SPAG 7.51RB at 17:35 on  4 Mar 2022
!> \brief \b zchkhe
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE zchkhe( DOTYPE, NN, NVAL, NNB, NBVAL, NNS, NSVAL,
!                          THRESH, TSTERR, NMAX, A, AFAC, AINV, B, X,
!                          XACT, WORK, RWORK, IWORK, NOUT )
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
!       COMPLEX*16         A( * ), AFAC( * ), AINV( * ), B( * ),
!      $                   WORK( * ), X( * ), XACT( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> zchkhe tests ZHETRF, -TRI2, -TRS, -TRS2, -RFS, and -CON.
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
!>          A is COMPLEX*16 array, dimension (NMAX*NMAX)
!> \endverbatim
!>
!> \param[out] AFAC
!> \verbatim
!>          AFAC is COMPLEX*16 array, dimension (NMAX*NMAX)
!> \endverbatim
!>
!> \param[out] AINV
!> \verbatim
!>          AINV is COMPLEX*16 array, dimension (NMAX*NMAX)
!> \endverbatim
!>
!> \param[out] B
!> \verbatim
!>          B is COMPLEX*16 array, dimension (NMAX*NSMAX)
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
!>          RWORK is DOUBLE PRECISION array, dimension (max(NMAX,2*NSMAX))
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
!> \date November 2013
!
!> \ingroup complex16_lin
!
!  =====================================================================
      SUBROUTINE zchkhe(Dotype,Nn,Nval,Nnb,Nbval,Nns,Nsval,Thresh,      &
     &                  Tsterr,Nmax,A,Afac,Ainv,B,X,Xact,Work,Rwork,    &
     &                  Iwork,Nout)
      IMPLICIT NONE
!*--ZCHKHE175
!
!  -- LAPACK test routine (version 3.5.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2013
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
      COMPLEX*16 A(*) , Afac(*) , Ainv(*) , B(*) , Work(*) , X(*) ,     &
     &           Xact(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D+0)
      COMPLEX*16 CZERO
      PARAMETER (CZERO=(0.0D+0,0.0D+0))
      INTEGER NTYPES
      PARAMETER (NTYPES=10)
      INTEGER NTESTS
      PARAMETER (NTESTS=9)
!     ..
!     .. Local Scalars ..
      LOGICAL trfcon , zerot
      CHARACTER dist , type , uplo , xtype
      CHARACTER*3 path
      INTEGER i , i1 , i2 , imat , in , inb , info , ioff , irhs ,      &
     &        iuplo , izero , j , k , kl , ku , lda , lwork , mode , n ,&
     &        nb , nerrs , nfail , nimat , nrhs , nrun , nt
      DOUBLE PRECISION anorm , cndnum , rcond , rcondc
!     ..
!     .. Local Arrays ..
      CHARACTER uplos(2)
      INTEGER iseed(4) , iseedy(4)
      DOUBLE PRECISION result(NTESTS)
!     ..
!     .. External Functions ..
      DOUBLE PRECISION ZLANHE
EXTERNAL    ZLANHE
!     ..
!     .. External Subroutines ..
EXTERNAL    ZHECON,   ZHERFS,   ZHETRF,   ZHETRI2,   ZHETRS
EXTERNAL    ZLACPY,   ZLATMS
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX , MIN
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
      path(1:1) = 'Zomplex precision'
      path(2:3) = 'HE'
      nrun = 0
      nfail = 0
      nerrs = 0
      DO i = 1 , 4
         iseed(i) = iseedy(i)
      ENDDO
!
!     Test the error exits
!
      IF ( Tsterr ) CALL zerrhe(path,Nout)
      INFot = 0
!
!     Set the minimum block size for which the block routine should
!     be used, which will be later returned by ilaenv
!
      CALL xlaenv(2,2)
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
!              Set up parameters with zlatb4 for the matrix generator
!              based on the type of matrix to be generated.
!
                     CALL zlatb4(path,imat,n,n,type,kl,ku,anorm,mode,   &
     &                           cndnum,dist)
!
!              Generate a matrix with ZLATMS.
!
                     SRNamt = 'ZLATMS'
                     CALL ZLATMS(n,n,dist,iseed,type,Rwork,mode,cndnum, &
     &                           anorm,kl,ku,uplo,A,lda,Work,info)
!
!              Check error code from ZLATMS and handle error.
!
                     IF ( info/=0 ) THEN
                        CALL alaerh(path,'ZLATMS',info,0,uplo,n,n,-1,-1,&
     &                              -1,imat,nfail,nerrs,Nout)
!
!                 Skip all tests for this generated matrix
!
                        CYCLE
                     ENDIF
!
!              For types 3-6, zero one or more rows and columns of
!              the matrix to test that INFO is returned correctly.
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
!                    Set row and column IZERO to zero.
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
!                       Set the first IZERO rows and columns to zero.
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
!                       Set the last IZERO rows and columns to zero.
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
!              End generate test matrix A.
!
!
!              Set the imaginary part of the diagonals.
!
                     CALL zlaipd(n,A,lda+1,0)
!
!              Do for each value of NB in NBVAL
!
                     DO inb = 1 , Nnb
!
!                 Set the optimal blocksize, which will be later
!                 returned by ilaenv.
!
                        nb = Nbval(inb)
                        CALL xlaenv(1,nb)
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
                        SRNamt = 'ZHETRF'
                        CALL ZHETRF(uplo,n,Afac,lda,Iwork,Ainv,lwork,   &
     &                              info)
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
!                 Check error code from ZHETRF and handle error.
!
                        IF ( info/=k )                                  &
     &                       CALL alaerh(path,'ZHETRF',info,k,uplo,n,n, &
     &                       -1,-1,nb,imat,nfail,nerrs,Nout)
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
                        CALL zhet01(uplo,n,A,lda,Afac,lda,Iwork,Ainv,   &
     &                              lda,Rwork,result(1))
                        nt = 1
!
!+    TEST 2
!                 Form the inverse and compute the residual.
!
                        IF ( inb==1 .AND. .NOT.trfcon ) THEN
                           CALL ZLACPY(uplo,n,n,Afac,lda,Ainv,lda)
                           SRNamt = 'ZHETRI2'
                           lwork = (n+nb+1)*(nb+3)
                           CALL ZHETRI2(uplo,n,Ainv,lda,Iwork,Work,     &
     &                                  lwork,info)
!
!                    Check error code from ZHETRI and handle error.
!
                           IF ( info/=0 )                               &
     &                          CALL alaerh(path,'ZHETRI',info,-1,uplo, &
     &                          n,n,-1,-1,-1,imat,nfail,nerrs,Nout)
!
!                    Compute the residual for a symmetric matrix times
!                    its inverse.
!
                           CALL zpot03(uplo,n,A,lda,Ainv,lda,Work,lda,  &
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
     &                             CALL alahd(Nout,path)
                              WRITE (Nout,FMT=99001) uplo , n , nb ,    &
     &                               imat , k , result(k)
                              nfail = nfail + 1
                           ENDIF
                        ENDDO
                        nrun = nrun + nt
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
!+    TEST 3 (Using TRS)
!                 Solve and compute residual for  A * X = B.
!
!                    Choose a set of NRHS random solution vectors
!                    stored in XACT and set up the right hand side B
!
                              SRNamt = 'zlarhs'
                              CALL zlarhs(path,xtype,uplo,' ',n,n,kl,ku,&
     &                           nrhs,A,lda,Xact,lda,B,lda,iseed,info)
                              CALL ZLACPY('Full',n,nrhs,B,lda,X,lda)
!
                              SRNamt = 'ZHETRS'
                              CALL ZHETRS(uplo,n,nrhs,Afac,lda,Iwork,X, &
     &                           lda,info)
!
!                    Check error code from ZHETRS and handle error.
!
                              IF ( info/=0 )                            &
     &                             CALL alaerh(path,'ZHETRS',info,0,    &
     &                             uplo,n,n,-1,-1,nrhs,imat,nfail,nerrs,&
     &                             Nout)
!
                              CALL ZLACPY('Full',n,nrhs,B,lda,Work,lda)
!
!                    Compute the residual for the solution
!
                              CALL zpot02(uplo,n,nrhs,A,lda,X,lda,Work, &
     &                           lda,Rwork,result(3))
!
!+    TEST 4 (Using TRS2)
!                 Solve and compute residual for  A * X = B.
!
!                    Choose a set of NRHS random solution vectors
!                    stored in XACT and set up the right hand side B
!
                              SRNamt = 'zlarhs'
                              CALL zlarhs(path,xtype,uplo,' ',n,n,kl,ku,&
     &                           nrhs,A,lda,Xact,lda,B,lda,iseed,info)
                              CALL ZLACPY('Full',n,nrhs,B,lda,X,lda)
!
                              SRNamt = 'ZHETRS2'
                              CALL ZHETRS2(uplo,n,nrhs,Afac,lda,Iwork,X,&
     &                           lda,Work,info)
!
!                    Check error code from ZHETRS2 and handle error.
!
                              IF ( info/=0 )                            &
     &                             CALL alaerh(path,'ZHETRS2',info,0,   &
     &                             uplo,n,n,-1,-1,nrhs,imat,nfail,nerrs,&
     &                             Nout)
!
                              CALL ZLACPY('Full',n,nrhs,B,lda,Work,lda)
!
!                    Compute the residual for the solution
!
                              CALL zpot02(uplo,n,nrhs,A,lda,X,lda,Work, &
     &                           lda,Rwork,result(4))
!
!+    TEST 5
!                 Check solution from generated exact solution.
!
                              CALL zget04(n,nrhs,X,lda,Xact,lda,rcondc, &
     &                           result(5))
!
!+    TESTS 6, 7, and 8
!                 Use iterative refinement to improve the solution.
!
                              SRNamt = 'ZHERFS'
                              CALL ZHERFS(uplo,n,nrhs,A,lda,Afac,lda,   &
     &                           Iwork,B,lda,X,lda,Rwork,Rwork(nrhs+1), &
     &                           Work,Rwork(2*nrhs+1),info)
!
!                    Check error code from ZHERFS.
!
                              IF ( info/=0 )                            &
     &                             CALL alaerh(path,'ZHERFS',info,0,    &
     &                             uplo,n,n,-1,-1,nrhs,imat,nfail,nerrs,&
     &                             Nout)
!
                              CALL zget04(n,nrhs,X,lda,Xact,lda,rcondc, &
     &                           result(6))
                              CALL zpot05(uplo,n,nrhs,A,lda,B,lda,X,lda,&
     &                           Xact,lda,Rwork,Rwork(nrhs+1),result(7))
!
!                    Print information about the tests that did not pass
!                    the threshold.
!
                              DO k = 3 , 8
                                 IF ( result(k)>=Thresh ) THEN
                                    IF ( nfail==0 .AND. nerrs==0 )      &
     &                                 CALL alahd(Nout,path)
                                    WRITE (Nout,FMT=99002) uplo , n ,   &
     &                                 nrhs , imat , k , result(k)
                                    nfail = nfail + 1
                                 ENDIF
                              ENDDO
                              nrun = nrun + 6
!
!                 End do for each value of NRHS in NSVAL.
!
                           ENDDO
!
!+    TEST 9
!                 Get an estimate of RCOND = 1/CNDNUM.
!
 2                         anorm = ZLANHE('1',uplo,n,A,lda,Rwork)
                           SRNamt = 'ZHECON'
                           CALL ZHECON(uplo,n,Afac,lda,Iwork,anorm,     &
     &                                 rcond,Work,info)
!
!                 Check error code from ZHECON and handle error.
!
                           IF ( info/=0 )                               &
     &                          CALL alaerh(path,'ZHECON',info,0,uplo,n,&
     &                          n,-1,-1,-1,imat,nfail,nerrs,Nout)
!
                           result(9) = dget06(rcond,rcondc)
!
!                 Print information about the tests that did not pass
!                 the threshold.
!
                           IF ( result(9)>=Thresh ) THEN
                              IF ( nfail==0 .AND. nerrs==0 )            &
     &                             CALL alahd(Nout,path)
                              WRITE (Nout,FMT=99003) uplo , n , imat ,  &
     &                               9 , result(9)
                              nfail = nfail + 1
                           ENDIF
                           nrun = nrun + 1
                        ENDIF
                     ENDDO
                  ENDDO
               ENDIF
            ENDIF
         ENDDO
      ENDDO
!
!     Print a summary of the results.
!
      CALL alasum(path,Nout,nfail,nrun,nerrs)
!
99001 FORMAT (' UPLO = ''',A1,''', N =',I5,', NB =',I4,', type ',I2,    &
     &        ', test ',I2,', ratio =',G12.5)
99002 FORMAT (' UPLO = ''',A1,''', N =',I5,', NRHS=',I3,', type ',I2,   &
     &        ', test(',I2,') =',G12.5)
99003 FORMAT (' UPLO = ''',A1,''', N =',I5,',',10X,' type ',I2,         &
     &        ', test(',I2,') =',G12.5)
!
!     End of zchkhe
!
      END SUBROUTINE zchkhe
