!*==ddrvsy_rk.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b DDRVSY_RK
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!     SUBROUTINE DDRVSY_RK( DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR,
!    $                      NMAX, A, AFAC, E, AINV, B, X, XACT, WORK,
!    $                      RWORK, IWORK, NOUT )
!
!       .. Scalar Arguments ..
!       LOGICAL            TSTERR
!       INTEGER            NMAX, NN, NOUT, NRHS
!       DOUBLE PRECISION   THRESH
!       ..
!       .. Array Arguments ..
!       LOGICAL            DOTYPE( * )
!       INTEGER            IWORK( * ), NVAL( * )
!       DOUBLE PRECISION   A( * ), AFAC( * ), AINV( * ), B( * ), E( * ),
!      $                   RWORK( * ), WORK( * ), X( * ), XACT( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!> DDRVSY_RK tests the driver routines DSYSV_RK.
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
!>          A is DOUBLE PRECISION array, dimension (NMAX*NMAX)
!> \endverbatim
!>
!> \param[out] AFAC
!> \verbatim
!>          AFAC is DOUBLE PRECISION array, dimension (NMAX*NMAX)
!> \endverbatim
!>
!> \param[out] E
!> \verbatim
!>          E is DOUBLE PRECISION array, dimension (NMAX)
!> \endverbatim
!>
!> \param[out] AINV
!> \verbatim
!>          AINV is DOUBLE PRECISION array, dimension (NMAX*NMAX)
!> \endverbatim
!>
!> \param[out] B
!> \verbatim
!>          B is DOUBLE PRECISION array, dimension (NMAX*NRHS)
!> \endverbatim
!>
!> \param[out] X
!> \verbatim
!>          X is DOUBLE PRECISION array, dimension (NMAX*NRHS)
!> \endverbatim
!>
!> \param[out] XACT
!> \verbatim
!>          XACT is DOUBLE PRECISION array, dimension (NMAX*NRHS)
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (NMAX*max(2,NRHS))
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is DOUBLE PRECISION array, dimension (NMAX+2*NRHS)
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
!> \ingroup double_lin
!
!  =====================================================================
      SUBROUTINE DDRVSY_RK(Dotype,Nn,Nval,Nrhs,Thresh,Tsterr,Nmax,A,    &
     &                     Afac,E,Ainv,B,X,Xact,Work,Rwork,Iwork,Nout)
      IMPLICIT NONE
!*--DDRVSY_RK159
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
      DOUBLE PRECISION A(*) , Afac(*) , Ainv(*) , B(*) , E(*) , Rwork(*)&
     &                 , Work(*) , X(*) , Xact(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ONE , ZERO
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
      INTEGER NTYPES , NTESTS
      PARAMETER (NTYPES=10,NTESTS=3)
      INTEGER NFACT
      PARAMETER (NFACT=2)
!     ..
!     .. Local Scalars ..
      LOGICAL zerot
      CHARACTER dist , fact , type , uplo , xtype
      CHARACTER*3 path , matpath
      INTEGER i , i1 , i2 , ifact , imat , in , info , ioff , iuplo ,   &
     &        izero , j , k , kl , ku , lda , lwork , mode , n , nb ,   &
     &        nbmin , nerrs , nfail , nimat , nrun , nt
      DOUBLE PRECISION ainvnm , anorm , cndnum , rcondc
!     ..
!     .. Local Arrays ..
      CHARACTER facts(NFACT) , uplos(2)
      INTEGER iseed(4) , iseedy(4)
      DOUBLE PRECISION result(NTESTS)
!     ..
!     .. External Functions ..
      DOUBLE PRECISION DLANSY
      EXTERNAL DLANSY
!     ..
!     .. External Subroutines ..
      EXTERNAL ALADHD , ALAERH , ALASVM , DERRVX , DGET04 , DLACPY ,    &
     &         DLARHS , DLATB4 , DLATMS , DPOT02 , DSYSV_RK , DSYT01_3 ,&
     &         DSYTRF_RK , DSYTRI_3 , XLAENV
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
      INTRINSIC MAX , MIN
!     ..
!     .. Data statements ..
      DATA iseedy/1988 , 1989 , 1990 , 1991/
      DATA uplos/'U' , 'L'/ , facts/'F' , 'N'/
!     ..
!     .. Executable Statements ..
!
!     Initialize constants and the random number seed.
!
!     Test path
!
      path(1:1) = 'Double precision'
      path(2:3) = 'SK'
!
!     Path to generate matrices
!
      matpath(1:1) = 'Double precision'
      matpath(2:3) = 'SY'
!
      nrun = 0
      nfail = 0
      nerrs = 0
      DO i = 1 , 4
         iseed(i) = iseedy(i)
      ENDDO
      lwork = MAX(2*Nmax,Nmax*Nrhs)
!
!     Test the error exits
!
      IF ( Tsterr ) CALL DERRVX(path,Nout)
      INFot = 0
!
!     Set the block size and minimum block size for which the block
!     routine should be used, which will be later returned by ILAENV.
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
!
!              Begin generate the test matrix A.
!
!              Set up parameters with DLATB4 for the matrix generator
!              based on the type of matrix to be generated.
!
                     CALL DLATB4(matpath,imat,n,n,type,kl,ku,anorm,mode,&
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
                        CALL ALAERH(path,'DLATMS',info,0,uplo,n,n,-1,-1,&
     &                              -1,imat,nfail,nerrs,Nout)
!
!                 Skip all tests for this generated matrix
!
                        CYCLE
                     ENDIF
!
!              For types 3-6, zero one or more rows and columns of the
!              matrix to test that INFO is returned correctly.
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
                                 ioff = ioff + lda
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
                                 ioff = ioff + lda
                              ENDDO
                           ENDIF
!
!                    Set row and column IZERO to zero.
!
                        ELSEIF ( iuplo==1 ) THEN
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
                     ELSE
                        izero = 0
                     ENDIF
!
!              End generate the test matrix A.
!
                     DO ifact = 1 , NFACT
!
!                 Do first for FACT = 'F', then for other values.
!
                        fact = facts(ifact)
!
!                 Compute the condition number
!
                        IF ( zerot ) THEN
                           IF ( ifact==1 ) CYCLE
                           rcondc = ZERO
!
                        ELSEIF ( ifact==1 ) THEN
!
!                    Compute the 1-norm of A.
!
                           anorm = DLANSY('1',uplo,n,A,lda,Rwork)
!
!                    Factor the matrix A.
!
                           CALL DLACPY(uplo,n,n,A,lda,Afac,lda)
                           CALL DSYTRF_RK(uplo,n,Afac,lda,E,Iwork,Work, &
     &                        lwork,info)
!
!                    Compute inv(A) and take its norm.
!
                           CALL DLACPY(uplo,n,n,Afac,lda,Ainv,lda)
                           lwork = (n+nb+1)*(nb+3)
!
!                    We need to compute the inverse to compute
!                    RCONDC that is used later in TEST3.
!
                           CALL DSYTRI_3(uplo,n,Ainv,lda,E,Iwork,Work,  &
     &                        lwork,info)
                           ainvnm = DLANSY('1',uplo,n,Ainv,lda,Rwork)
!
!                    Compute the 1-norm condition number of A.
!
                           IF ( anorm<=ZERO .OR. ainvnm<=ZERO ) THEN
                              rcondc = ONE
                           ELSE
                              rcondc = (ONE/anorm)/ainvnm
                           ENDIF
                        ENDIF
!
!                 Form an exact solution and set the right hand side.
!
                        SRNamt = 'DLARHS'
                        CALL DLARHS(matpath,xtype,uplo,' ',n,n,kl,ku,   &
     &                              Nrhs,A,lda,Xact,lda,B,lda,iseed,    &
     &                              info)
                        xtype = 'C'
!
!                 --- Test DSYSV_RK  ---
!
                        IF ( ifact==2 ) THEN
                           CALL DLACPY(uplo,n,n,A,lda,Afac,lda)
                           CALL DLACPY('Full',n,Nrhs,B,lda,X,lda)
!
!                    Factor the matrix and solve the system using
!                    DSYSV_RK.
!
                           SRNamt = 'DSYSV_RK'
                           CALL DSYSV_RK(uplo,n,Nrhs,Afac,lda,E,Iwork,X,&
     &                        lda,Work,lwork,info)
!
!                    Adjust the expected value of INFO to account for
!                    pivoting.
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
!                    Check error code from DSYSV_RK and handle error.
!
                           IF ( info/=k ) THEN
                              CALL ALAERH(path,'DSYSV_RK',info,k,uplo,n,&
     &                           n,-1,-1,Nrhs,imat,nfail,nerrs,Nout)
                              CYCLE
                           ELSEIF ( info/=0 ) THEN
                              CYCLE
                           ENDIF
!
!+    TEST 1      Reconstruct matrix from factors and compute
!                 residual.
!
                           CALL DSYT01_3(uplo,n,A,lda,Afac,lda,E,Iwork, &
     &                        Ainv,lda,Rwork,result(1))
!
!+    TEST 2      Compute residual of the computed solution.
!
                           CALL DLACPY('Full',n,Nrhs,B,lda,Work,lda)
                           CALL DPOT02(uplo,n,Nrhs,A,lda,X,lda,Work,lda,&
     &                                 Rwork,result(2))
!
!+    TEST 3
!                 Check solution from generated exact solution.
!
                           CALL DGET04(n,Nrhs,X,lda,Xact,lda,rcondc,    &
     &                                 result(3))
                           nt = 3
!
!                    Print information about the tests that did not pass
!                    the threshold.
!
                           DO k = 1 , nt
                              IF ( result(k)>=Thresh ) THEN
                                 IF ( nfail==0 .AND. nerrs==0 )         &
     &                                CALL ALADHD(Nout,path)
                                 WRITE (Nout,FMT=99001) 'DSYSV_RK' ,    &
     &                                  uplo , n , imat , k , result(k)
                                 nfail = nfail + 1
                              ENDIF
                           ENDDO
                           nrun = nrun + nt
                        ENDIF
!
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
      CALL ALASVM(path,Nout,nfail,nrun,nerrs)
!
99001 FORMAT (1X,A,', UPLO=''',A1,''', N =',I5,', type ',I2,', test ',  &
     &        I2,', ratio =',G12.5)
!
!     End of DDRVSY_RK
!
      END SUBROUTINE DDRVSY_RK
