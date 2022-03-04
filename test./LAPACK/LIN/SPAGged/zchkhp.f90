!*==zchkhp.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b ZCHKHP
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZCHKHP( DOTYPE, NN, NVAL, NNS, NSVAL, THRESH, TSTERR,
!                          NMAX, A, AFAC, AINV, B, X, XACT, WORK, RWORK,
!                          IWORK, NOUT )
!
!       .. Scalar Arguments ..
!       LOGICAL            TSTERR
!       INTEGER            NMAX, NN, NNS, NOUT
!       DOUBLE PRECISION   THRESH
!       ..
!       .. Array Arguments ..
!       LOGICAL            DOTYPE( * )
!       INTEGER            IWORK( * ), NSVAL( * ), NVAL( * )
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
!> ZCHKHP tests ZHPTRF, -TRI, -TRS, -RFS, and -CON
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
!>          A is COMPLEX*16 array, dimension
!>                      (NMAX*(NMAX+1)/2)
!> \endverbatim
!>
!> \param[out] AFAC
!> \verbatim
!>          AFAC is COMPLEX*16 array, dimension
!>                      (NMAX*(NMAX+1)/2)
!> \endverbatim
!>
!> \param[out] AINV
!> \verbatim
!>          AINV is COMPLEX*16 array, dimension
!>                      (NMAX*(NMAX+1)/2)
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
!>          WORK is COMPLEX*16 array, dimension
!>                      (NMAX*max(2,NSMAX))
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is DOUBLE PRECISION array,
!>                                 dimension (NMAX+2*NSMAX)
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
      SUBROUTINE ZCHKHP(Dotype,Nn,Nval,Nns,Nsval,Thresh,Tsterr,Nmax,A,  &
     &                  Afac,Ainv,B,X,Xact,Work,Rwork,Iwork,Nout)
      IMPLICIT NONE
!*--ZCHKHP167
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      LOGICAL Tsterr
      INTEGER Nmax , Nn , Nns , Nout
      DOUBLE PRECISION Thresh
!     ..
!     .. Array Arguments ..
      LOGICAL Dotype(*)
      INTEGER Iwork(*) , Nsval(*) , Nval(*)
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
      INTEGER NTYPES
      PARAMETER (NTYPES=10)
      INTEGER NTESTS
      PARAMETER (NTESTS=8)
!     ..
!     .. Local Scalars ..
      LOGICAL trfcon , zerot
      CHARACTER dist , packit , type , uplo , xtype
      CHARACTER*3 path
      INTEGER i , i1 , i2 , imat , in , info , ioff , irhs , iuplo ,    &
     &        izero , j , k , kl , ku , lda , mode , n , nerrs , nfail ,&
     &        nimat , npp , nrhs , nrun , nt
      DOUBLE PRECISION anorm , cndnum , rcond , rcondc
!     ..
!     .. Local Arrays ..
      CHARACTER uplos(2)
      INTEGER iseed(4) , iseedy(4)
      DOUBLE PRECISION result(NTESTS)
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      DOUBLE PRECISION DGET06 , ZLANHP
      EXTERNAL LSAME , DGET06 , ZLANHP
!     ..
!     .. External Subroutines ..
      EXTERNAL ALAERH , ALAHD , ALASUM , ZCOPY , ZERRSY , ZGET04 ,      &
     &         ZHPCON , ZHPRFS , ZHPT01 , ZHPTRF , ZHPTRI , ZHPTRS ,    &
     &         ZLACPY , ZLAIPD , ZLARHS , ZLATB4 , ZLATMS , ZPPT02 ,    &
     &         ZPPT03 , ZPPT05
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
      path(2:3) = 'HP'
      nrun = 0
      nfail = 0
      nerrs = 0
      DO i = 1 , 4
         iseed(i) = iseedy(i)
      ENDDO
!
!     Test the error exits
!
      IF ( Tsterr ) CALL ZERRSY(path,Nout)
      INFot = 0
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
                     IF ( LSAME(uplo,'U') ) THEN
                        packit = 'C'
                     ELSE
                        packit = 'R'
                     ENDIF
!
!              Set up parameters with ZLATB4 and generate a test matrix
!              with ZLATMS.
!
                     CALL ZLATB4(path,imat,n,n,type,kl,ku,anorm,mode,   &
     &                           cndnum,dist)
!
                     SRNamt = 'ZLATMS'
                     CALL ZLATMS(n,n,dist,iseed,type,Rwork,mode,cndnum, &
     &                           anorm,kl,ku,packit,A,lda,Work,info)
!
!              Check error code from ZLATMS.
!
                     IF ( info/=0 ) THEN
                        CALL ALAERH(path,'ZLATMS',info,0,uplo,n,n,-1,-1,&
     &                              -1,imat,nfail,nerrs,Nout)
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
                        IF ( imat>=6 ) THEN
                           ioff = 0
                           IF ( iuplo==1 ) THEN
!
!                       Set the first IZERO rows and columns to zero.
!
                              DO j = 1 , n
                                 i2 = MIN(j,izero)
                                 DO i = 1 , i2
                                    A(ioff+i) = ZERO
                                 ENDDO
                                 ioff = ioff + j
                              ENDDO
                           ELSE
!
!                       Set the last IZERO rows and columns to zero.
!
                              DO j = 1 , n
                                 i1 = MAX(j,izero)
                                 DO i = i1 , n
                                    A(ioff+i) = ZERO
                                 ENDDO
                                 ioff = ioff + n - j
                              ENDDO
                           ENDIF
!
!                    Set row and column IZERO to zero.
!
                        ELSEIF ( iuplo==1 ) THEN
                           ioff = (izero-1)*izero/2
                           DO i = 1 , izero - 1
                              A(ioff+i) = ZERO
                           ENDDO
                           ioff = ioff + izero
                           DO i = izero , n
                              A(ioff) = ZERO
                              ioff = ioff + i
                           ENDDO
                        ELSE
                           ioff = izero
                           DO i = 1 , izero - 1
                              A(ioff) = ZERO
                              ioff = ioff + n - i
                           ENDDO
                           ioff = ioff - izero
                           DO i = izero , n
                              A(ioff+i) = ZERO
                           ENDDO
                        ENDIF
                     ELSE
                        izero = 0
                     ENDIF
!
!              Set the imaginary part of the diagonals.
!
                     IF ( iuplo==1 ) THEN
                        CALL ZLAIPD(n,A,2,1)
                     ELSE
                        CALL ZLAIPD(n,A,n,-1)
                     ENDIF
!
!              Compute the L*D*L' or U*D*U' factorization of the matrix.
!
                     npp = n*(n+1)/2
                     CALL ZCOPY(npp,A,1,Afac,1)
                     SRNamt = 'ZHPTRF'
                     CALL ZHPTRF(uplo,n,Afac,Iwork,info)
!
!              Adjust the expected value of INFO to account for
!              pivoting.
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
!              Check error code from ZHPTRF.
!
                     IF ( info/=k )                                     &
     &                    CALL ALAERH(path,'ZHPTRF',info,k,uplo,n,n,-1, &
     &                    -1,-1,imat,nfail,nerrs,Nout)
                     IF ( info/=0 ) THEN
                        trfcon = .TRUE.
                     ELSE
                        trfcon = .FALSE.
                     ENDIF
!
!+    TEST 1
!              Reconstruct matrix from factors and compute residual.
!
                     CALL ZHPT01(uplo,n,A,Afac,Iwork,Ainv,lda,Rwork,    &
     &                           result(1))
                     nt = 1
!
!+    TEST 2
!              Form the inverse and compute the residual.
!
                     IF ( .NOT.trfcon ) THEN
                        CALL ZCOPY(npp,Afac,1,Ainv,1)
                        SRNamt = 'ZHPTRI'
                        CALL ZHPTRI(uplo,n,Ainv,Iwork,Work,info)
!
!              Check error code from ZHPTRI.
!
                        IF ( info/=0 )                                  &
     &                       CALL ALAERH(path,'ZHPTRI',info,0,uplo,n,n, &
     &                       -1,-1,-1,imat,nfail,nerrs,Nout)
!
                        CALL ZPPT03(uplo,n,A,Ainv,Work,lda,Rwork,rcondc,&
     &                              result(2))
                        nt = 2
                     ENDIF
!
!              Print information about the tests that did not pass
!              the threshold.
!
                     DO k = 1 , nt
                        IF ( result(k)>=Thresh ) THEN
                           IF ( nfail==0 .AND. nerrs==0 )               &
     &                          CALL ALAHD(Nout,path)
                           WRITE (Nout,FMT=99001) uplo , n , imat , k , &
     &                            result(k)
                           nfail = nfail + 1
                        ENDIF
                     ENDDO
                     nrun = nrun + nt
!
!              Do only the condition estimate if INFO is not 0.
!
                     IF ( trfcon ) THEN
                        rcondc = ZERO
                        GOTO 2
                     ENDIF
!
                     DO irhs = 1 , Nns
                        nrhs = Nsval(irhs)
!
!+    TEST 3
!              Solve and compute residual for  A * X = B.
!
                        SRNamt = 'ZLARHS'
                        CALL ZLARHS(path,xtype,uplo,' ',n,n,kl,ku,nrhs, &
     &                              A,lda,Xact,lda,B,lda,iseed,info)
                        xtype = 'C'
                        CALL ZLACPY('Full',n,nrhs,B,lda,X,lda)
!
                        SRNamt = 'ZHPTRS'
                        CALL ZHPTRS(uplo,n,nrhs,Afac,Iwork,X,lda,info)
!
!              Check error code from ZHPTRS.
!
                        IF ( info/=0 )                                  &
     &                       CALL ALAERH(path,'ZHPTRS',info,0,uplo,n,n, &
     &                       -1,-1,nrhs,imat,nfail,nerrs,Nout)
!
                        CALL ZLACPY('Full',n,nrhs,B,lda,Work,lda)
                        CALL ZPPT02(uplo,n,nrhs,A,X,lda,Work,lda,Rwork, &
     &                              result(3))
!
!+    TEST 4
!              Check solution from generated exact solution.
!
                        CALL ZGET04(n,nrhs,X,lda,Xact,lda,rcondc,       &
     &                              result(4))
!
!+    TESTS 5, 6, and 7
!              Use iterative refinement to improve the solution.
!
                        SRNamt = 'ZHPRFS'
                        CALL ZHPRFS(uplo,n,nrhs,A,Afac,Iwork,B,lda,X,   &
     &                              lda,Rwork,Rwork(nrhs+1),Work,       &
     &                              Rwork(2*nrhs+1),info)
!
!              Check error code from ZHPRFS.
!
                        IF ( info/=0 )                                  &
     &                       CALL ALAERH(path,'ZHPRFS',info,0,uplo,n,n, &
     &                       -1,-1,nrhs,imat,nfail,nerrs,Nout)
!
                        CALL ZGET04(n,nrhs,X,lda,Xact,lda,rcondc,       &
     &                              result(5))
                        CALL ZPPT05(uplo,n,nrhs,A,B,lda,X,lda,Xact,lda, &
     &                              Rwork,Rwork(nrhs+1),result(6))
!
!                 Print information about the tests that did not pass
!                 the threshold.
!
                        DO k = 3 , 7
                           IF ( result(k)>=Thresh ) THEN
                              IF ( nfail==0 .AND. nerrs==0 )            &
     &                             CALL ALAHD(Nout,path)
                              WRITE (Nout,FMT=99002) uplo , n , nrhs ,  &
     &                               imat , k , result(k)
                              nfail = nfail + 1
                           ENDIF
                        ENDDO
                        nrun = nrun + 5
                     ENDDO
!
!+    TEST 8
!              Get an estimate of RCOND = 1/CNDNUM.
!
 2                   anorm = ZLANHP('1',uplo,n,A,Rwork)
                     SRNamt = 'ZHPCON'
                     CALL ZHPCON(uplo,n,Afac,Iwork,anorm,rcond,Work,    &
     &                           info)
!
!              Check error code from ZHPCON.
!
                     IF ( info/=0 )                                     &
     &                    CALL ALAERH(path,'ZHPCON',info,0,uplo,n,n,-1, &
     &                    -1,-1,imat,nfail,nerrs,Nout)
!
                     result(8) = DGET06(rcond,rcondc)
!
!              Print the test ratio if it is .GE. THRESH.
!
                     IF ( result(8)>=Thresh ) THEN
                        IF ( nfail==0 .AND. nerrs==0 )                  &
     &                       CALL ALAHD(Nout,path)
                        WRITE (Nout,FMT=99001) uplo , n , imat , 8 ,    &
     &                         result(8)
                        nfail = nfail + 1
                     ENDIF
                     nrun = nrun + 1
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
99001 FORMAT (' UPLO = ''',A1,''', N =',I5,', type ',I2,', test ',I2,   &
     &        ', ratio =',G12.5)
99002 FORMAT (' UPLO = ''',A1,''', N =',I5,', NRHS=',I3,', type ',I2,   &
     &        ', test(',I2,') =',G12.5)
!
!     End of ZCHKHP
!
      END SUBROUTINE ZCHKHP