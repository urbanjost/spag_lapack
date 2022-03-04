!*==schkpo.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b SCHKPO
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE SCHKPO( DOTYPE, NN, NVAL, NNB, NBVAL, NNS, NSVAL,
!                          THRESH, TSTERR, NMAX, A, AFAC, AINV, B, X,
!                          XACT, WORK, RWORK, IWORK, NOUT )
!
!       .. Scalar Arguments ..
!       LOGICAL            TSTERR
!       INTEGER            NMAX, NN, NNB, NNS, NOUT
!       REAL               THRESH
!       ..
!       .. Array Arguments ..
!       LOGICAL            DOTYPE( * )
!       INTEGER            IWORK( * ), NBVAL( * ), NSVAL( * ), NVAL( * )
!       REAL               A( * ), AFAC( * ), AINV( * ), B( * ),
!      $                   RWORK( * ), WORK( * ), X( * ), XACT( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SCHKPO tests SPOTRF, -TRI, -TRS, -RFS, and -CON
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
!>          A is REAL array, dimension (NMAX*NMAX)
!> \endverbatim
!>
!> \param[out] AFAC
!> \verbatim
!>          AFAC is REAL array, dimension (NMAX*NMAX)
!> \endverbatim
!>
!> \param[out] AINV
!> \verbatim
!>          AINV is REAL array, dimension (NMAX*NMAX)
!> \endverbatim
!>
!> \param[out] B
!> \verbatim
!>          B is REAL array, dimension (NMAX*NSMAX)
!>          where NSMAX is the largest entry in NSVAL.
!> \endverbatim
!>
!> \param[out] X
!> \verbatim
!>          X is REAL array, dimension (NMAX*NSMAX)
!> \endverbatim
!>
!> \param[out] XACT
!> \verbatim
!>          XACT is REAL array, dimension (NMAX*NSMAX)
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is REAL array, dimension
!>                      (NMAX*max(3,NSMAX))
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is REAL array, dimension
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
!> \ingroup single_lin
!
!  =====================================================================
      SUBROUTINE SCHKPO(Dotype,Nn,Nval,Nnb,Nbval,Nns,Nsval,Thresh,      &
     &                  Tsterr,Nmax,A,Afac,Ainv,B,X,Xact,Work,Rwork,    &
     &                  Iwork,Nout)
      IMPLICIT NONE
!*--SCHKPO176
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      LOGICAL Tsterr
      INTEGER Nmax , Nn , Nnb , Nns , Nout
      REAL Thresh
!     ..
!     .. Array Arguments ..
      LOGICAL Dotype(*)
      INTEGER Iwork(*) , Nbval(*) , Nsval(*) , Nval(*)
      REAL A(*) , Afac(*) , Ainv(*) , B(*) , Rwork(*) , Work(*) , X(*) ,&
     &     Xact(*)
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
      CHARACTER dist , type , uplo , xtype
      CHARACTER*3 path
      INTEGER i , imat , in , inb , info , ioff , irhs , iuplo , izero ,&
     &        k , kl , ku , lda , mode , n , nb , nerrs , nfail ,       &
     &        nimat , nrhs , nrun
      REAL anorm , cndnum , rcond , rcondc
!     ..
!     .. Local Arrays ..
      CHARACTER uplos(2)
      INTEGER iseed(4) , iseedy(4)
      REAL result(NTESTS)
!     ..
!     .. External Functions ..
      REAL SGET06 , SLANSY
      EXTERNAL SGET06 , SLANSY
!     ..
!     .. External Subroutines ..
      EXTERNAL ALAERH , ALAHD , ALASUM , SERRPO , SGET04 , SLACPY ,     &
     &         SLARHS , SLATB4 , SLATMS , SPOCON , SPORFS , SPOT01 ,    &
     &         SPOT02 , SPOT03 , SPOT05 , SPOTRF , SPOTRI , SPOTRS ,    &
     &         XLAENV
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
      DATA uplos/'U' , 'L'/
!     ..
!     .. Executable Statements ..
!
!     Initialize constants and the random number seed.
!
      path(1:1) = 'Single precision'
      path(2:3) = 'PO'
      nrun = 0
      nfail = 0
      nerrs = 0
      DO i = 1 , 4
         iseed(i) = iseedy(i)
      ENDDO
!
!     Test the error exits
!
      IF ( Tsterr ) CALL SERRPO(path,Nout)
      INFot = 0
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
!
!              Set up parameters with SLATB4 and generate a test matrix
!              with SLATMS.
!
                     CALL SLATB4(path,imat,n,n,type,kl,ku,anorm,mode,   &
     &                           cndnum,dist)
!
                     SRNamt = 'SLATMS'
                     CALL SLATMS(n,n,dist,iseed,type,Rwork,mode,cndnum, &
     &                           anorm,kl,ku,uplo,A,lda,Work,info)
!
!              Check error code from SLATMS.
!
                     IF ( info/=0 ) THEN
                        CALL ALAERH(path,'SLATMS',info,0,uplo,n,n,-1,-1,&
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
                        ioff = (izero-1)*lda
!
!                 Set row and column IZERO of A to 0.
!
                        IF ( iuplo==1 ) THEN
                           DO i = 1 , izero - 1
                              A(ioff+i) = ZERO
                           ENDDO
                           ioff = ioff + izero
                           DO i = izero , n
                              A(ioff) = ZERO
                              ioff = ioff + lda
                           ENDDO
                        ELSE
                           ioff = izero
                           DO i = 1 , izero - 1
                              A(ioff) = ZERO
                              ioff = ioff + lda
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
!              Do for each value of NB in NBVAL
!
                     DO inb = 1 , Nnb
                        nb = Nbval(inb)
                        CALL XLAENV(1,nb)
!
!                 Compute the L*L' or U'*U factorization of the matrix.
!
                        CALL SLACPY(uplo,n,n,A,lda,Afac,lda)
                        SRNamt = 'SPOTRF'
                        CALL SPOTRF(uplo,n,Afac,lda,info)
!
!                 Check error code from SPOTRF.
!
                        IF ( info/=izero ) THEN
                           CALL ALAERH(path,'SPOTRF',info,izero,uplo,n, &
     &                                 n,-1,-1,nb,imat,nfail,nerrs,Nout)
                           CYCLE
                        ENDIF
!
!                 Skip the tests if INFO is not 0.
!
                        IF ( info==0 ) THEN
!
!+    TEST 1
!                 Reconstruct matrix from factors and compute residual.
!
                           CALL SLACPY(uplo,n,n,Afac,lda,Ainv,lda)
                           CALL SPOT01(uplo,n,A,lda,Ainv,lda,Rwork,     &
     &                                 result(1))
!
!+    TEST 2
!                 Form the inverse and compute the residual.
!
                           CALL SLACPY(uplo,n,n,Afac,lda,Ainv,lda)
                           SRNamt = 'SPOTRI'
                           CALL SPOTRI(uplo,n,Ainv,lda,info)
!
!                 Check error code from SPOTRI.
!
                           IF ( info/=0 )                               &
     &                          CALL ALAERH(path,'SPOTRI',info,0,uplo,n,&
     &                          n,-1,-1,-1,imat,nfail,nerrs,Nout)
!
                           CALL SPOT03(uplo,n,A,lda,Ainv,lda,Work,lda,  &
     &                                 Rwork,rcondc,result(2))
!
!                 Print information about the tests that did not pass
!                 the threshold.
!
                           DO k = 1 , 2
                              IF ( result(k)>=Thresh ) THEN
                                 IF ( nfail==0 .AND. nerrs==0 )         &
     &                                CALL ALAHD(Nout,path)
                                 WRITE (Nout,FMT=99001) uplo , n , nb , &
     &                                  imat , k , result(k)
                                 nfail = nfail + 1
                              ENDIF
                           ENDDO
                           nrun = nrun + 2
!
!                 Skip the rest of the tests unless this is the first
!                 blocksize.
!
                           IF ( inb==1 ) THEN
!
                              DO irhs = 1 , Nns
                                 nrhs = Nsval(irhs)
!
!+    TEST 3
!                 Solve and compute residual for A * X = B .
!
                                 SRNamt = 'SLARHS'
                                 CALL SLARHS(path,xtype,uplo,' ',n,n,kl,&
     &                              ku,nrhs,A,lda,Xact,lda,B,lda,iseed, &
     &                              info)
                                 CALL SLACPY('Full',n,nrhs,B,lda,X,lda)
!
                                 SRNamt = 'SPOTRS'
                                 CALL SPOTRS(uplo,n,nrhs,Afac,lda,X,lda,&
     &                              info)
!
!                 Check error code from SPOTRS.
!
                                 IF ( info/=0 )                         &
     &                                 CALL ALAERH(path,'SPOTRS',info,0,&
     &                                uplo,n,n,-1,-1,nrhs,imat,nfail,   &
     &                                nerrs,Nout)
!
                                 CALL SLACPY('Full',n,nrhs,B,lda,Work,  &
     &                              lda)
                                 CALL SPOT02(uplo,n,nrhs,A,lda,X,lda,   &
     &                              Work,lda,Rwork,result(3))
!
!+    TEST 4
!                 Check solution from generated exact solution.
!
                                 CALL SGET04(n,nrhs,X,lda,Xact,lda,     &
     &                              rcondc,result(4))
!
!+    TESTS 5, 6, and 7
!                 Use iterative refinement to improve the solution.
!
                                 SRNamt = 'SPORFS'
                                 CALL SPORFS(uplo,n,nrhs,A,lda,Afac,lda,&
     &                              B,lda,X,lda,Rwork,Rwork(nrhs+1),    &
     &                              Work,Iwork,info)
!
!                 Check error code from SPORFS.
!
                                 IF ( info/=0 )                         &
     &                                 CALL ALAERH(path,'SPORFS',info,0,&
     &                                uplo,n,n,-1,-1,nrhs,imat,nfail,   &
     &                                nerrs,Nout)
!
                                 CALL SGET04(n,nrhs,X,lda,Xact,lda,     &
     &                              rcondc,result(5))
                                 CALL SPOT05(uplo,n,nrhs,A,lda,B,lda,X, &
     &                              lda,Xact,lda,Rwork,Rwork(nrhs+1),   &
     &                              result(6))
!
!                    Print information about the tests that did not pass
!                    the threshold.
!
                                 DO k = 3 , 7
                                    IF ( result(k)>=Thresh ) THEN
                                       IF ( nfail==0 .AND. nerrs==0 )   &
     &                                    CALL ALAHD(Nout,path)
                                       WRITE (Nout,FMT=99002) uplo , n ,&
     &                                    nrhs , imat , k , result(k)
                                       nfail = nfail + 1
                                    ENDIF
                                 ENDDO
                                 nrun = nrun + 5
                              ENDDO
!
!+    TEST 8
!                 Get an estimate of RCOND = 1/CNDNUM.
!
                              anorm = SLANSY('1',uplo,n,A,lda,Rwork)
                              SRNamt = 'SPOCON'
                              CALL SPOCON(uplo,n,Afac,lda,anorm,rcond,  &
     &                           Work,Iwork,info)
!
!                 Check error code from SPOCON.
!
                              IF ( info/=0 )                            &
     &                             CALL ALAERH(path,'SPOCON',info,0,    &
     &                             uplo,n,n,-1,-1,-1,imat,nfail,nerrs,  &
     &                             Nout)
!
                              result(8) = SGET06(rcond,rcondc)
!
!                 Print the test ratio if it is .GE. THRESH.
!
                              IF ( result(8)>=Thresh ) THEN
                                 IF ( nfail==0 .AND. nerrs==0 )         &
     &                                CALL ALAHD(Nout,path)
                                 WRITE (Nout,FMT=99003) uplo , n ,      &
     &                                  imat , 8 , result(8)
                                 nfail = nfail + 1
                              ENDIF
                              nrun = nrun + 1
                           ENDIF
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
      CALL ALASUM(path,Nout,nfail,nrun,nerrs)
!
99001 FORMAT (' UPLO = ''',A1,''', N =',I5,', NB =',I4,', type ',I2,    &
     &        ', test ',I2,', ratio =',G12.5)
99002 FORMAT (' UPLO = ''',A1,''', N =',I5,', NRHS=',I3,', type ',I2,   &
     &        ', test(',I2,') =',G12.5)
99003 FORMAT (' UPLO = ''',A1,''', N =',I5,',',10X,' type ',I2,         &
     &        ', test(',I2,') =',G12.5)
!
!     End of SCHKPO
!
      END SUBROUTINE SCHKPO
