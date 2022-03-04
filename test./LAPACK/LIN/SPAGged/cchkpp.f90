!*==cchkpp.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b CCHKPP
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE CCHKPP( DOTYPE, NN, NVAL, NNS, NSVAL, THRESH, TSTERR,
!                          NMAX, A, AFAC, AINV, B, X, XACT, WORK, RWORK,
!                          NOUT )
!
!       .. Scalar Arguments ..
!       LOGICAL            TSTERR
!       INTEGER            NMAX, NN, NNS, NOUT
!       REAL               THRESH
!       ..
!       .. Array Arguments ..
!       LOGICAL            DOTYPE( * )
!       INTEGER            NSVAL( * ), NVAL( * )
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
!> CCHKPP tests CPPTRF, -TRI, -TRS, -RFS, and -CON
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
!>                      (NMAX*max(3,NSMAX))
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is REAL array, dimension
!>                      (max(NMAX,2*NSMAX))
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
      SUBROUTINE CCHKPP(Dotype,Nn,Nval,Nns,Nsval,Thresh,Tsterr,Nmax,A,  &
     &                  Afac,Ainv,B,X,Xact,Work,Rwork,Nout)
      IMPLICIT NONE
!*--CCHKPP162
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
      INTEGER Nsval(*) , Nval(*)
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
      PARAMETER (NTYPES=9)
      INTEGER NTESTS
      PARAMETER (NTESTS=8)
!     ..
!     .. Local Scalars ..
      LOGICAL zerot
      CHARACTER dist , packit , type , uplo , xtype
      CHARACTER*3 path
      INTEGER i , imat , in , info , ioff , irhs , iuplo , izero , k ,  &
     &        kl , ku , lda , mode , n , nerrs , nfail , nimat , npp ,  &
     &        nrhs , nrun
      REAL anorm , cndnum , rcond , rcondc
!     ..
!     .. Local Arrays ..
      CHARACTER packs(2) , uplos(2)
      INTEGER iseed(4) , iseedy(4)
      REAL result(NTESTS)
!     ..
!     .. External Functions ..
      REAL CLANHP , SGET06
      EXTERNAL CLANHP , SGET06
!     ..
!     .. External Subroutines ..
      EXTERNAL ALAERH , ALAHD , ALASUM , CCOPY , CERRPO , CGET04 ,      &
     &         CLACPY , CLAIPD , CLARHS , CLATB4 , CLATMS , CPPCON ,    &
     &         CPPRFS , CPPT01 , CPPT02 , CPPT03 , CPPT05 , CPPTRF ,    &
     &         CPPTRI , CPPTRS
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
!     .. Intrinsic Functions ..
      INTRINSIC MAX
!     ..
!     .. Data statements ..
      DATA iseedy/1988 , 1989 , 1990 , 1991/
      DATA uplos/'U' , 'L'/ , packs/'C' , 'R'/
!     ..
!     .. Executable Statements ..
!
!     Initialize constants and the random number seed.
!
      path(1:1) = 'Complex precision'
      path(2:3) = 'PP'
      nrun = 0
      nfail = 0
      nerrs = 0
      DO i = 1 , 4
         iseed(i) = iseedy(i)
      ENDDO
!
!     Test the error exits
!
      IF ( Tsterr ) CALL CERRPO(path,Nout)
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
!           Skip types 3, 4, or 5 if the matrix size is too small.
!
               zerot = imat>=3 .AND. imat<=5
               IF ( .NOT.(zerot .AND. n<imat-2) ) THEN
!
!           Do first for UPLO = 'U', then for UPLO = 'L'
!
                  DO iuplo = 1 , 2
                     uplo = uplos(iuplo)
                     packit = packs(iuplo)
!
!              Set up parameters with CLATB4 and generate a test matrix
!              with CLATMS.
!
                     CALL CLATB4(path,imat,n,n,type,kl,ku,anorm,mode,   &
     &                           cndnum,dist)
!
                     SRNamt = 'CLATMS'
                     CALL CLATMS(n,n,dist,iseed,type,Rwork,mode,cndnum, &
     &                           anorm,kl,ku,packit,A,lda,Work,info)
!
!              Check error code from CLATMS.
!
                     IF ( info/=0 ) THEN
                        CALL ALAERH(path,'CLATMS',info,0,uplo,n,n,-1,-1,&
     &                              -1,imat,nfail,nerrs,Nout)
                        CYCLE
                     ENDIF
!
!              For types 3-5, zero one row and column of the matrix to
!              test that INFO is returned correctly.
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
!                 Set row and column IZERO of A to 0.
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
                     ELSE
                        izero = 0
                     ENDIF
!
!              Set the imaginary part of the diagonals.
!
                     IF ( iuplo==1 ) THEN
                        CALL CLAIPD(n,A,2,1)
                     ELSE
                        CALL CLAIPD(n,A,n,-1)
                     ENDIF
!
!              Compute the L*L' or U'*U factorization of the matrix.
!
                     npp = n*(n+1)/2
                     CALL CCOPY(npp,A,1,Afac,1)
                     SRNamt = 'CPPTRF'
                     CALL CPPTRF(uplo,n,Afac,info)
!
!              Check error code from CPPTRF.
!
                     IF ( info/=izero ) THEN
                        CALL ALAERH(path,'CPPTRF',info,izero,uplo,n,n,  &
     &                              -1,-1,-1,imat,nfail,nerrs,Nout)
                        CYCLE
                     ENDIF
!
!              Skip the tests if INFO is not 0.
!
                     IF ( info==0 ) THEN
!
!+    TEST 1
!              Reconstruct matrix from factors and compute residual.
!
                        CALL CCOPY(npp,Afac,1,Ainv,1)
                        CALL CPPT01(uplo,n,A,Ainv,Rwork,result(1))
!
!+    TEST 2
!              Form the inverse and compute the residual.
!
                        CALL CCOPY(npp,Afac,1,Ainv,1)
                        SRNamt = 'CPPTRI'
                        CALL CPPTRI(uplo,n,Ainv,info)
!
!              Check error code from CPPTRI.
!
                        IF ( info/=0 )                                  &
     &                       CALL ALAERH(path,'CPPTRI',info,0,uplo,n,n, &
     &                       -1,-1,-1,imat,nfail,nerrs,Nout)
!
                        CALL CPPT03(uplo,n,A,Ainv,Work,lda,Rwork,rcondc,&
     &                              result(2))
!
!              Print information about the tests that did not pass
!              the threshold.
!
                        DO k = 1 , 2
                           IF ( result(k)>=Thresh ) THEN
                              IF ( nfail==0 .AND. nerrs==0 )            &
     &                             CALL ALAHD(Nout,path)
                              WRITE (Nout,FMT=99001) uplo , n , imat ,  &
     &                               k , result(k)
                              nfail = nfail + 1
                           ENDIF
                        ENDDO
                        nrun = nrun + 2
!
                        DO irhs = 1 , Nns
                           nrhs = Nsval(irhs)
!
!+    TEST 3
!              Solve and compute residual for  A * X = B.
!
                           SRNamt = 'CLARHS'
                           CALL CLARHS(path,xtype,uplo,' ',n,n,kl,ku,   &
     &                                 nrhs,A,lda,Xact,lda,B,lda,iseed, &
     &                                 info)
                           CALL CLACPY('Full',n,nrhs,B,lda,X,lda)
!
                           SRNamt = 'CPPTRS'
                           CALL CPPTRS(uplo,n,nrhs,Afac,X,lda,info)
!
!              Check error code from CPPTRS.
!
                           IF ( info/=0 )                               &
     &                          CALL ALAERH(path,'CPPTRS',info,0,uplo,n,&
     &                          n,-1,-1,nrhs,imat,nfail,nerrs,Nout)
!
                           CALL CLACPY('Full',n,nrhs,B,lda,Work,lda)
                           CALL CPPT02(uplo,n,nrhs,A,X,lda,Work,lda,    &
     &                                 Rwork,result(3))
!
!+    TEST 4
!              Check solution from generated exact solution.
!
                           CALL CGET04(n,nrhs,X,lda,Xact,lda,rcondc,    &
     &                                 result(4))
!
!+    TESTS 5, 6, and 7
!              Use iterative refinement to improve the solution.
!
                           SRNamt = 'CPPRFS'
                           CALL CPPRFS(uplo,n,nrhs,A,Afac,B,lda,X,lda,  &
     &                                 Rwork,Rwork(nrhs+1),Work,        &
     &                                 Rwork(2*nrhs+1),info)
!
!              Check error code from CPPRFS.
!
                           IF ( info/=0 )                               &
     &                          CALL ALAERH(path,'CPPRFS',info,0,uplo,n,&
     &                          n,-1,-1,nrhs,imat,nfail,nerrs,Nout)
!
                           CALL CGET04(n,nrhs,X,lda,Xact,lda,rcondc,    &
     &                                 result(5))
                           CALL CPPT05(uplo,n,nrhs,A,B,lda,X,lda,Xact,  &
     &                                 lda,Rwork,Rwork(nrhs+1),result(6)&
     &                                 )
!
!                 Print information about the tests that did not pass
!                 the threshold.
!
                           DO k = 3 , 7
                              IF ( result(k)>=Thresh ) THEN
                                 IF ( nfail==0 .AND. nerrs==0 )         &
     &                                CALL ALAHD(Nout,path)
                                 WRITE (Nout,FMT=99002) uplo , n ,      &
     &                                  nrhs , imat , k , result(k)
                                 nfail = nfail + 1
                              ENDIF
                           ENDDO
                           nrun = nrun + 5
                        ENDDO
!
!+    TEST 8
!              Get an estimate of RCOND = 1/CNDNUM.
!
                        anorm = CLANHP('1',uplo,n,A,Rwork)
                        SRNamt = 'CPPCON'
                        CALL CPPCON(uplo,n,Afac,anorm,rcond,Work,Rwork, &
     &                              info)
!
!              Check error code from CPPCON.
!
                        IF ( info/=0 )                                  &
     &                       CALL ALAERH(path,'CPPCON',info,0,uplo,n,n, &
     &                       -1,-1,-1,imat,nfail,nerrs,Nout)
!
                        result(8) = SGET06(rcond,rcondc)
!
!              Print the test ratio if greater than or equal to THRESH.
!
                        IF ( result(8)>=Thresh ) THEN
                           IF ( nfail==0 .AND. nerrs==0 )               &
     &                          CALL ALAHD(Nout,path)
                           WRITE (Nout,FMT=99001) uplo , n , imat , 8 , &
     &                            result(8)
                           nfail = nfail + 1
                        ENDIF
                        nrun = nrun + 1
                     ENDIF
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
99001 FORMAT (' UPLO = ''',A1,''', N =',I5,', type ',I2,', test ',I2,   &
     &        ', ratio =',G12.5)
99002 FORMAT (' UPLO = ''',A1,''', N =',I5,', NRHS=',I3,', type ',I2,   &
     &        ', test(',I2,') =',G12.5)
!
!     End of CCHKPP
!
      END SUBROUTINE CCHKPP
