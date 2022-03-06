!*==cchksp.f90  processed by SPAG 7.51RB at 17:34 on  4 Mar 2022
!> \brief \b cchksp
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE cchksp( DOTYPE, NN, NVAL, NNS, NSVAL, THRESH, TSTERR,
!                          NMAX, A, AFAC, AINV, B, X, XACT, WORK, RWORK,
!                          IWORK, NOUT )
!
!       .. Scalar Arguments ..
!       LOGICAL            TSTERR
!       INTEGER            NMAX, NN, NNS, NOUT
!       REAL               THRESH
!       ..
!       .. Array Arguments ..
!       LOGICAL            DOTYPE( * )
!       INTEGER            IWORK( * ), NSVAL( * ), NVAL( * )
!       REAL               RWORK( * )
!       COMPLEX            A( * ), AFAC( * ), AINV( * ), B( * ),
!      $                   WORK( * ), X( * ), XACT( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> cchksp tests CSPTRF, -TRI, -TRS, -RFS, and -CON
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
!>          THRESH is REAL
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
!>          A is COMPLEX array, dimension
!>                      (NMAX*(NMAX+1)/2)
!> \endverbatim
!>
!> \param[out] AFAC
!> \verbatim
!>          AFAC is COMPLEX array, dimension
!>                      (NMAX*(NMAX+1)/2)
!> \endverbatim
!>
!> \param[out] AINV
!> \verbatim
!>          AINV is COMPLEX array, dimension
!>                      (NMAX*(NMAX+1)/2)
!> \endverbatim
!>
!> \param[out] B
!> \verbatim
!>          B is COMPLEX array, dimension (NMAX*NSMAX)
!>          where NSMAX is the largest entry in NSVAL.
!> \endverbatim
!>
!> \param[out] X
!> \verbatim
!>          X is COMPLEX array, dimension (NMAX*NSMAX)
!> \endverbatim
!>
!> \param[out] XACT
!> \verbatim
!>          XACT is COMPLEX array, dimension (NMAX*NSMAX)
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX array, dimension
!>                      (NMAX*max(2,NSMAX))
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is REAL array,
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
!> \ingroup complex_lin
!
!  =====================================================================
      SUBROUTINE cchksp(Dotype,Nn,Nval,Nns,Nsval,Thresh,Tsterr,Nmax,A,  &
     &                  Afac,Ainv,B,X,Xact,Work,Rwork,Iwork,Nout)
      IMPLICIT NONE
!*--CCHKSP167
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      LOGICAL Tsterr
      INTEGER Nmax , Nn , Nns , Nout
      REAL Thresh
!     ..
!     .. Array Arguments ..
      LOGICAL Dotype(*)
      INTEGER Iwork(*) , Nsval(*) , Nval(*)
      REAL Rwork(*)
      COMPLEX A(*) , Afac(*) , Ainv(*) , B(*) , Work(*) , X(*) , Xact(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ZERO
      PARAMETER (ZERO=0.0E+0)
      INTEGER NTYPES
      PARAMETER (NTYPES=11)
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
      REAL anorm , cndnum , rcond , rcondc
!     ..
!     .. Local Arrays ..
      CHARACTER uplos(2)
      INTEGER iseed(4) , iseedy(4)
      REAL result(NTESTS)
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      REAL CLANSP
EXTERNAL    LSAME,   CLANSP
!     ..
!     .. External Subroutines ..
EXTERNAL    CCOPY
EXTERNAL    CLACPY,   CLATMS
EXTERNAL    CSPCON,   CSPRFS,   CSPTRF
EXTERNAL    CSPTRI,   CSPTRS
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
      path(1:1) = 'Complex precision'
      path(2:3) = 'SP'
      nrun = 0
      nfail = 0
      nerrs = 0
      DO i = 1 , 4
         iseed(i) = iseedy(i)
      ENDDO
!
!     Test the error exits
!
      IF ( Tsterr ) CALL cerrsy(path,Nout)
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
                     IF ( imat/=NTYPES ) THEN
!
!                 Set up parameters with clatb4 and generate a test
!                 matrix with CLATMS.
!
                        CALL clatb4(path,imat,n,n,type,kl,ku,anorm,mode,&
     &                              cndnum,dist)
!
                        SRNamt = 'CLATMS'
                        CALL CLATMS(n,n,dist,iseed,type,Rwork,mode,     &
     &                              cndnum,anorm,kl,ku,packit,A,lda,    &
     &                              Work,info)
!
!                 Check error code from CLATMS.
!
                        IF ( info/=0 ) THEN
                           CALL alaerh(path,'CLATMS',info,0,uplo,n,n,-1,&
     &                                 -1,-1,imat,nfail,nerrs,Nout)
                           CYCLE
                        ENDIF
!
!                 For types 3-6, zero one or more rows and columns of
!                 the matrix to test that INFO is returned correctly.
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
                           ELSEIF ( iuplo==1 ) THEN
!
!                          Set the first IZERO rows and columns to zero.
!
                              ioff = 0
                              DO j = 1 , n
                                 i2 = MIN(j,izero)
                                 DO i = 1 , i2
                                    A(ioff+i) = ZERO
                                 ENDDO
                                 ioff = ioff + j
                              ENDDO
                           ELSE
!
!                          Set the last IZERO rows and columns to zero.
!
                              ioff = 0
                              DO j = 1 , n
                                 i1 = MAX(j,izero)
                                 DO i = i1 , n
                                    A(ioff+i) = ZERO
                                 ENDDO
                                 ioff = ioff + n - j
                              ENDDO
                           ENDIF
                        ELSE
                           izero = 0
                        ENDIF
                     ELSE
!
!                 Use a special block diagonal matrix to test alternate
!                 code for the 2 x 2 blocks.
!
                        CALL clatsp(uplo,n,A,iseed)
                     ENDIF
!
!              Compute the L*D*L' or U*D*U' factorization of the matrix.
!
                     npp = n*(n+1)/2
                     CALL CCOPY(npp,A,1,Afac,1)
                     SRNamt = 'CSPTRF'
                     CALL CSPTRF(uplo,n,Afac,Iwork,info)
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
!              Check error code from CSPTRF.
!
                     IF ( info/=k )                                     &
     &                    CALL alaerh(path,'CSPTRF',info,k,uplo,n,n,-1, &
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
                     CALL cspt01(uplo,n,A,Afac,Iwork,Ainv,lda,Rwork,    &
     &                           result(1))
                     nt = 1
!
!+    TEST 2
!              Form the inverse and compute the residual.
!
                     IF ( .NOT.trfcon ) THEN
                        CALL CCOPY(npp,Afac,1,Ainv,1)
                        SRNamt = 'CSPTRI'
                        CALL CSPTRI(uplo,n,Ainv,Iwork,Work,info)
!
!              Check error code from CSPTRI.
!
                        IF ( info/=0 )                                  &
     &                       CALL alaerh(path,'CSPTRI',info,0,uplo,n,n, &
     &                       -1,-1,-1,imat,nfail,nerrs,Nout)
!
                        CALL cspt03(uplo,n,A,Ainv,Work,lda,Rwork,rcondc,&
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
     &                          CALL alahd(Nout,path)
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
                        SRNamt = 'clarhs'
                        CALL clarhs(path,xtype,uplo,' ',n,n,kl,ku,nrhs, &
     &                              A,lda,Xact,lda,B,lda,iseed,info)
                        CALL CLACPY('Full',n,nrhs,B,lda,X,lda)
!
                        SRNamt = 'CSPTRS'
                        CALL CSPTRS(uplo,n,nrhs,Afac,Iwork,X,lda,info)
!
!              Check error code from CSPTRS.
!
                        IF ( info/=0 )                                  &
     &                       CALL alaerh(path,'CSPTRS',info,0,uplo,n,n, &
     &                       -1,-1,nrhs,imat,nfail,nerrs,Nout)
!
                        CALL CLACPY('Full',n,nrhs,B,lda,Work,lda)
                        CALL cspt02(uplo,n,nrhs,A,X,lda,Work,lda,Rwork, &
     &                              result(3))
!
!+    TEST 4
!              Check solution from generated exact solution.
!
                        CALL cget04(n,nrhs,X,lda,Xact,lda,rcondc,       &
     &                              result(4))
!
!+    TESTS 5, 6, and 7
!              Use iterative refinement to improve the solution.
!
                        SRNamt = 'CSPRFS'
                        CALL CSPRFS(uplo,n,nrhs,A,Afac,Iwork,B,lda,X,   &
     &                              lda,Rwork,Rwork(nrhs+1),Work,       &
     &                              Rwork(2*nrhs+1),info)
!
!              Check error code from CSPRFS.
!
                        IF ( info/=0 )                                  &
     &                       CALL alaerh(path,'CSPRFS',info,0,uplo,n,n, &
     &                       -1,-1,nrhs,imat,nfail,nerrs,Nout)
!
                        CALL cget04(n,nrhs,X,lda,Xact,lda,rcondc,       &
     &                              result(5))
                        CALL cppt05(uplo,n,nrhs,A,B,lda,X,lda,Xact,lda, &
     &                              Rwork,Rwork(nrhs+1),result(6))
!
!                 Print information about the tests that did not pass
!                 the threshold.
!
                        DO k = 3 , 7
                           IF ( result(k)>=Thresh ) THEN
                              IF ( nfail==0 .AND. nerrs==0 )            &
     &                             CALL alahd(Nout,path)
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
 2                   anorm = CLANSP('1',uplo,n,A,Rwork)
                     SRNamt = 'CSPCON'
                     CALL CSPCON(uplo,n,Afac,Iwork,anorm,rcond,Work,    &
     &                           info)
!
!              Check error code from CSPCON.
!
                     IF ( info/=0 )                                     &
     &                    CALL alaerh(path,'CSPCON',info,0,uplo,n,n,-1, &
     &                    -1,-1,imat,nfail,nerrs,Nout)
!
                     result(8) = sget06(rcond,rcondc)
!
!              Print the test ratio if it is .GE. THRESH.
!
                     IF ( result(8)>=Thresh ) THEN
                        IF ( nfail==0 .AND. nerrs==0 )                  &
     &                       CALL alahd(Nout,path)
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
      CALL alasum(path,Nout,nfail,nrun,nerrs)
!
99001 FORMAT (' UPLO = ''',A1,''', N =',I5,', type ',I2,', test ',I2,   &
     &        ', ratio =',G12.5)
99002 FORMAT (' UPLO = ''',A1,''', N =',I5,', NRHS=',I3,', type ',I2,   &
     &        ', test(',I2,') =',G12.5)
!
!     End of cchksp
!
      END SUBROUTINE cchksp
