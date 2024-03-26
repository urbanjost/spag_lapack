!*==dchksy_aa.f90  processed by SPAG 7.51RB at 17:34 on  4 Mar 2022
!> \brief \b dchksy_aa
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE dchksy_aa( DOTYPE, NN, NVAL, NNB, NBVAL, NNS, NSVAL,
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
!       DOUBLE PRECISION   A( * ), AFAC( * ), AINV( * ), B( * ),
!      $                   RWORK( * ), WORK( * ), X( * ), XACT( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> dchksy_aa tests DSYTRF_AA, -TRS_AA.
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
!>          A is DOUBLE PRECISION array, dimension (NMAX*NMAX)
!> \endverbatim
!>
!> \param[out] AFAC
!> \verbatim
!>          AFAC is DOUBLE PRECISION array, dimension (NMAX*NMAX)
!> \endverbatim
!>
!> \param[out] AINV
!> \verbatim
!>          AINV is DOUBLE PRECISION array, dimension (NMAX*NMAX)
!> \endverbatim
!>
!> \param[out] B
!> \verbatim
!>          B is DOUBLE PRECISION array, dimension (NMAX*NSMAX)
!>          where NSMAX is the largest entry in NSVAL.
!> \endverbatim
!>
!> \param[out] X
!> \verbatim
!>          X is DOUBLE PRECISION array, dimension (NMAX*NSMAX)
!> \endverbatim
!>
!> \param[out] XACT
!> \verbatim
!>          XACT is DOUBLE PRECISION array, dimension (NMAX*NSMAX)
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (NMAX*max(3,NSMAX))
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is DOUBLE PRECISION array, dimension (max(NMAX,2*NSMAX))
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
!> \date November 2017
!
!  @precisions fortran d -> z c
!
!> \ingroup double_lin
!
!  =====================================================================
      SUBROUTINE dchksy_aa(Dotype,Nn,Nval,Nnb,Nbval,Nns,Nsval,Thresh,   &
     &                     Tsterr,Nmax,A,Afac,Ainv,B,X,Xact,Work,Rwork, &
     &                     Iwork,Nout)
!
!  -- LAPACK test routine (version 3.8.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2017
!
      IMPLICIT NONE
!*--DCHKSY_AA182
!
!     .. Scalar Arguments ..
      LOGICAL Tsterr
      INTEGER Nn , Nnb , Nns , Nmax , Nout
      DOUBLE PRECISION Thresh
!     ..
!     .. Array Arguments ..
      LOGICAL Dotype(*)
      INTEGER Iwork(*) , Nbval(*) , Nsval(*) , Nval(*)
      DOUBLE PRECISION A(*) , Afac(*) , Ainv(*) , B(*) , Rwork(*) ,     &
     &                 Work(*) , X(*) , Xact(*)
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
      PARAMETER (NTESTS=9)
!     ..
!     .. Local Scalars ..
      LOGICAL zerot
      CHARACTER dist , type , uplo , xtype
      CHARACTER*3 path , matpath
      INTEGER i , i1 , i2 , imat , in , inb , info , ioff , irhs ,      &
     &        iuplo , izero , j , k , kl , ku , lda , lwork , mode , n ,&
     &        nb , nerrs , nfail , nimat , nrhs , nrun , nt
      DOUBLE PRECISION anorm , cndnum
!     ..
!     .. Local Arrays ..
      CHARACTER uplos(2)
      INTEGER iseed(4) , iseedy(4)
      DOUBLE PRECISION result(NTESTS)
!     ..
!     .. External Subroutines ..
EXTERNAL    DLACPY
EXTERNAL    DSYTRF_AA
EXTERNAL    DSYTRS_AA
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
!     Test path
!
      path(1:1) = 'Double precision'
      path(2:3) = 'SA'
!
!     Path to generate matrices
!
      matpath(1:1) = 'Double precision'
      matpath(2:3) = 'SY'
      nrun = 0
      nfail = 0
      nerrs = 0
      DO i = 1 , 4
         iseed(i) = iseedy(i)
      ENDDO
!
!     Test the error exits
!
      IF ( Tsterr ) CALL derrsy(path,Nout)
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
         IF ( n>Nmax ) THEN
            nfail = nfail + 1
            WRITE (Nout,99001) 'M ' , n , Nmax
99001       FORMAT (' Invalid input value: ',A4,'=',I6,'; must be <=',  &
     &              I6)
            CYCLE
         ENDIF
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
!              Begin generate the test matrix A.
!
!
!              Set up parameters with dlatb4 for the matrix generator
!              based on the type of matrix to be generated.
!
                     CALL dlatb4(matpath,imat,n,n,type,kl,ku,anorm,mode,&
     &                           cndnum,dist)
!
!              Generate a matrix with DLATMS.
!
                     SRNamt = 'DLATMS'
                     CALL DLATMS(n,n,dist,iseed,type,Rwork,mode,cndnum, &
     &                           anorm,kl,ku,uplo,A,lda,Work,info)
!
!              Check error code from DLATMS and handle error.
!
                     IF ( info/=0 ) THEN
                        CALL alaerh(path,'DLATMS',info,0,uplo,n,n,-1,-1,&
     &                              -1,imat,nfail,nerrs,Nout)
!
!                    Skip all tests for this generated matrix
!
                        CYCLE
                     ENDIF
!
!              For matrix types 3-6, zero one or more rows and
!              columns of the matrix to test that INFO is returned
!              correctly.
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
                        ELSEIF ( iuplo==1 ) THEN
!
!                       Set the first IZERO rows and columns to zero.
!
                           ioff = 0
                           DO j = 1 , n
                              i2 = MIN(j,izero)
                              DO i = 1 , i2
                                 A(ioff+i) = ZERO
                              ENDDO
                              ioff = ioff + lda
                           ENDDO
                           izero = 1
                        ELSE
!
!                       Set the last IZERO rows and columns to zero.
!
                           ioff = 0
                           DO j = 1 , n
                              i1 = MAX(j,izero)
                              DO i = i1 , n
                                 A(ioff+i) = ZERO
                              ENDDO
                              ioff = ioff + lda
                           ENDDO
                        ENDIF
                     ELSE
                        izero = 0
                     ENDIF
!
!              End generate the test matrix A.
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
                        CALL DLACPY(uplo,n,n,A,lda,Afac,lda)
!
!                 Compute the L*D*L**T or U*D*U**T factorization of the
!                 matrix. IWORK stores details of the interchanges and
!                 the block structure of D. AINV is a work array for
!                 block factorization, LWORK is the length of AINV.
!
                        SRNamt = 'DSYTRF_AA'
                        lwork = MAX(1,n*nb+n)
                        CALL DSYTRF_AA(uplo,n,Afac,lda,Iwork,Ainv,lwork,&
     &                                 info)
!
!                 Adjust the expected value of INFO to account for
!                 pivoting.
!
!                  IF( IZERO.GT.0 ) THEN
!                     J = 1
!                     K = IZERO
!  100                CONTINUE
!                     IF( J.EQ.K ) THEN
!                        K = IWORK( J )
!                     ELSE IF( IWORK( J ).EQ.K ) THEN
!                        K = J
!                     END IF
!                     IF( J.LT.K ) THEN
!                        J = J + 1
!                        GO TO 100
!                     END IF
!                  ELSE
                        k = 0
!                  END IF
!
!                 Check error code from DSYTRF and handle error.
!
                        IF ( info/=k )                                  &
     &                       CALL alaerh(path,'DSYTRF_AA',info,k,uplo,n,&
     &                       n,-1,-1,nb,imat,nfail,nerrs,Nout)
!
!+    TEST 1
!                 Reconstruct matrix from factors and compute residual.
!
                        CALL dsyt01_aa(uplo,n,A,lda,Afac,lda,Iwork,Ainv,&
     &                                 lda,Rwork,result(1))
                        nt = 1
!
!
!                 Print information about the tests that did not pass
!                 the threshold.
!
                        DO k = 1 , nt
                           IF ( result(k)>=Thresh ) THEN
                              IF ( nfail==0 .AND. nerrs==0 )            &
     &                             CALL alahd(Nout,path)
                              WRITE (Nout,FMT=99002) uplo , n , nb ,    &
     &                               imat , k , result(k)
                              nfail = nfail + 1
                           ENDIF
                        ENDDO
                        nrun = nrun + nt
!
!                 Skip solver test if INFO is not 0.
!
                        IF ( info/=0 ) CYCLE
!
!                 Do for each value of NRHS in NSVAL.
!
                        DO irhs = 1 , Nns
                           nrhs = Nsval(irhs)
!
!+    TEST 2 (Using TRS)
!                 Solve and compute residual for  A * X = B.
!
!                    Choose a set of NRHS random solution vectors
!                    stored in XACT and set up the right hand side B
!
                           SRNamt = 'dlarhs'
                           CALL dlarhs(matpath,xtype,uplo,' ',n,n,kl,ku,&
     &                                 nrhs,A,lda,Xact,lda,B,lda,iseed, &
     &                                 info)
                           CALL DLACPY('Full',n,nrhs,B,lda,X,lda)
!
                           SRNamt = 'DSYTRS_AA'
                           lwork = MAX(1,3*n-2)
                           CALL DSYTRS_AA(uplo,n,nrhs,Afac,lda,Iwork,X, &
     &                        lda,Work,lwork,info)
!
!                    Check error code from DSYTRS and handle error.
!
                           IF ( info==0 ) THEN
                              CALL DLACPY('Full',n,nrhs,B,lda,Work,lda)
!
!                       Compute the residual for the solution
!
                              CALL dpot02(uplo,n,nrhs,A,lda,X,lda,Work, &
     &                           lda,Rwork,result(2))
!
!
!                       Print information about the tests that did not pass
!                       the threshold.
!
                              DO k = 2 , 2
                                 IF ( result(k)>=Thresh ) THEN
                                    IF ( nfail==0 .AND. nerrs==0 )      &
     &                                 CALL alahd(Nout,path)
                                    WRITE (Nout,FMT=99003) uplo , n ,   &
     &                                 nrhs , imat , k , result(k)
                                    nfail = nfail + 1
                                 ENDIF
                              ENDDO
                           ELSEIF ( izero==0 ) THEN
                              CALL alaerh(path,'DSYTRS_AA',info,0,uplo, &
     &                           n,n,-1,-1,nrhs,imat,nfail,nerrs,Nout)
                           ENDIF
                           nrun = nrun + 1
!
!                 End do for each value of NRHS in NSVAL.
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
      CALL alasum(path,Nout,nfail,nrun,nerrs)
!
99002 FORMAT (' UPLO = ''',A1,''', N =',I5,', NB =',I4,', type ',I2,    &
     &        ', test ',I2,', ratio =',G12.5)
99003 FORMAT (' UPLO = ''',A1,''', N =',I5,', NRHS=',I3,', type ',I2,   &
     &        ', test(',I2,') =',G12.5)
!
!     End of dchksy_aa
!
      END SUBROUTINE dchksy_aa
