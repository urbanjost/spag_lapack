!*==ddrvsp.f90  processed by SPAG 7.51RB at 17:34 on  4 Mar 2022
!> \brief \b ddrvsp
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ddrvsp( DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, NMAX,
!                          A, AFAC, AINV, B, X, XACT, WORK, RWORK, IWORK,
!                          NOUT )
!
!       .. Scalar Arguments ..
!       LOGICAL            TSTERR
!       INTEGER            NMAX, NN, NOUT, NRHS
!       DOUBLE PRECISION   THRESH
!       ..
!       .. Array Arguments ..
!       LOGICAL            DOTYPE( * )
!       INTEGER            IWORK( * ), NVAL( * )
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
!> ddrvsp tests the driver routines DSPSV and -SVX.
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
!>          A is DOUBLE PRECISION array, dimension
!>                      (NMAX*(NMAX+1)/2)
!> \endverbatim
!>
!> \param[out] AFAC
!> \verbatim
!>          AFAC is DOUBLE PRECISION array, dimension
!>                      (NMAX*(NMAX+1)/2)
!> \endverbatim
!>
!> \param[out] AINV
!> \verbatim
!>          AINV is DOUBLE PRECISION array, dimension
!>                      (NMAX*(NMAX+1)/2)
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
!>          WORK is DOUBLE PRECISION array, dimension
!>                      (NMAX*max(2,NRHS))
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
      SUBROUTINE ddrvsp(Dotype,Nn,Nval,Nrhs,Thresh,Tsterr,Nmax,A,Afac,  &
     &                  Ainv,B,X,Xact,Work,Rwork,Iwork,Nout)
      IMPLICIT NONE
!*--DDRVSP159
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
      DOUBLE PRECISION A(*) , Afac(*) , Ainv(*) , B(*) , Rwork(*) ,     &
     &                 Work(*) , X(*) , Xact(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ONE , ZERO
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
      INTEGER NTYPES , NTESTS
      PARAMETER (NTYPES=10,NTESTS=6)
      INTEGER NFACT
      PARAMETER (NFACT=2)
!     ..
!     .. Local Scalars ..
      LOGICAL zerot
      CHARACTER dist , fact , packit , type , uplo , xtype
      CHARACTER*3 path
      INTEGER i , i1 , i2 , ifact , imat , in , info , ioff , iuplo ,   &
     &        izero , j , k , k1 , kl , ku , lda , lwork , mode , n ,   &
     &        nerrs , nfail , nimat , npp , nrun , nt
      DOUBLE PRECISION ainvnm , anorm , cndnum , rcond , rcondc
!     ..
!     .. Local Arrays ..
      CHARACTER facts(NFACT)
      INTEGER iseed(4) , iseedy(4)
      DOUBLE PRECISION result(NTESTS)
!     ..
!     .. External Functions ..
      DOUBLE PRECISION DLANSP
EXTERNAL    DLANSP
!     ..
!     .. External Subroutines ..
EXTERNAL    DCOPY
EXTERNAL    DLACPY,   DLASET
EXTERNAL    DSPSV,   DSPSVX,   DSPTRF,   DSPTRI
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
      DATA facts/'F' , 'N'/
!     ..
!     .. Executable Statements ..
!
!     Initialize constants and the random number seed.
!
      path(1:1) = 'Double precision'
      path(2:3) = 'SP'
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
      IF ( Tsterr ) CALL derrvx(path,Nout)
      INFot = 0
!
!     Do for each value of N in NVAL
!
      DO in = 1 , Nn
         n = Nval(in)
         lda = MAX(n,1)
         npp = n*(n+1)/2
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
                     IF ( iuplo==1 ) THEN
                        uplo = 'U'
                        packit = 'C'
                     ELSE
                        uplo = 'L'
                        packit = 'R'
                     ENDIF
!
!              Set up parameters with dlatb4 and generate a test matrix
!              with DLATMS.
!
                     CALL dlatb4(path,imat,n,n,type,kl,ku,anorm,mode,   &
     &                           cndnum,dist)
!
                     SRNamt = 'DLATMS'
                     CALL DLATMS(n,n,dist,iseed,type,Rwork,mode,cndnum, &
     &                           anorm,kl,ku,packit,A,lda,Work,info)
!
!              Check error code from DLATMS.
!
                     IF ( info/=0 ) THEN
                        CALL alaerh(path,'DLATMS',info,0,uplo,n,n,-1,-1,&
     &                              -1,imat,nfail,nerrs,Nout)
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
                     DO ifact = 1 , NFACT
!
!                 Do first for FACT = 'F', then for other values.
!
                        fact = facts(ifact)
!
!                 Compute the condition number for comparison with
!                 the value returned by DSPSVX.
!
                        IF ( zerot ) THEN
                           IF ( ifact==1 ) CYCLE
                           rcondc = ZERO
!
                        ELSEIF ( ifact==1 ) THEN
!
!                    Compute the 1-norm of A.
!
                           anorm = DLANSP('1',uplo,n,A,Rwork)
!
!                    Factor the matrix A.
!
                           CALL DCOPY(npp,A,1,Afac,1)
                           CALL DSPTRF(uplo,n,Afac,Iwork,info)
!
!                    Compute inv(A) and take its norm.
!
                           CALL DCOPY(npp,Afac,1,Ainv,1)
                           CALL DSPTRI(uplo,n,Ainv,Iwork,Work,info)
                           ainvnm = DLANSP('1',uplo,n,Ainv,Rwork)
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
                        SRNamt = 'dlarhs'
                        CALL dlarhs(path,xtype,uplo,' ',n,n,kl,ku,Nrhs, &
     &                              A,lda,Xact,lda,B,lda,iseed,info)
                        xtype = 'C'
!
!                 --- Test DSPSV  ---
!
                        IF ( ifact==2 ) THEN
                           CALL DCOPY(npp,A,1,Afac,1)
                           CALL DLACPY('Full',n,Nrhs,B,lda,X,lda)
!
!                    Factor the matrix and solve the system using DSPSV.
!
                           SRNamt = 'DSPSV '
                           CALL DSPSV(uplo,n,Nrhs,Afac,Iwork,X,lda,info)
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
!                    Check error code from DSPSV .
!
                           IF ( info/=k ) THEN
                              CALL alaerh(path,'DSPSV ',info,k,uplo,n,n,&
     &                           -1,-1,Nrhs,imat,nfail,nerrs,Nout)
                              GOTO 2
                           ELSEIF ( info/=0 ) THEN
                              GOTO 2
                           ENDIF
!
!                    Reconstruct matrix from factors and compute
!                    residual.
!
                           CALL dspt01(uplo,n,A,Afac,Iwork,Ainv,lda,    &
     &                                 Rwork,result(1))
!
!                    Compute residual of the computed solution.
!
                           CALL DLACPY('Full',n,Nrhs,B,lda,Work,lda)
                           CALL dppt02(uplo,n,Nrhs,A,X,lda,Work,lda,    &
     &                                 Rwork,result(2))
!
!                    Check solution from generated exact solution.
!
                           CALL dget04(n,Nrhs,X,lda,Xact,lda,rcondc,    &
     &                                 result(3))
                           nt = 3
!
!                    Print information about the tests that did not pass
!                    the threshold.
!
                           DO k = 1 , nt
                              IF ( result(k)>=Thresh ) THEN
                                 IF ( nfail==0 .AND. nerrs==0 )         &
     &                                CALL aladhd(Nout,path)
                                 WRITE (Nout,FMT=99001) 'DSPSV ' ,      &
     &                                  uplo , n , imat , k , result(k)
                                 nfail = nfail + 1
                              ENDIF
                           ENDDO
                           nrun = nrun + nt
                        ENDIF
!
!                 --- Test DSPSVX ---
!
 2                      IF ( ifact==2 .AND. npp>0 )                     &
     &                       CALL DLASET('Full',npp,1,ZERO,ZERO,Afac,   &
     &                       npp)
                        CALL DLASET('Full',n,Nrhs,ZERO,ZERO,X,lda)
!
!                 Solve the system and compute the condition number and
!                 error bounds using DSPSVX.
!
                        SRNamt = 'DSPSVX'
                        CALL DSPSVX(fact,uplo,n,Nrhs,A,Afac,Iwork,B,lda,&
     &                              X,lda,rcond,Rwork,Rwork(Nrhs+1),    &
     &                              Work,Iwork(n+1),info)
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
!                 Check the error code from DSPSVX.
!
                        IF ( info/=k ) THEN
                           CALL alaerh(path,'DSPSVX',info,k,fact//uplo, &
     &                                 n,n,-1,-1,Nrhs,imat,nfail,nerrs, &
     &                                 Nout)
                           CYCLE
                        ENDIF
!
                        IF ( info==0 ) THEN
                           IF ( ifact>=2 ) THEN
!
!                       Reconstruct matrix from factors and compute
!                       residual.
!
                              CALL dspt01(uplo,n,A,Afac,Iwork,Ainv,lda, &
     &                           Rwork(2*Nrhs+1),result(1))
                              k1 = 1
                           ELSE
                              k1 = 2
                           ENDIF
!
!                    Compute residual of the computed solution.
!
                           CALL DLACPY('Full',n,Nrhs,B,lda,Work,lda)
                           CALL dppt02(uplo,n,Nrhs,A,X,lda,Work,lda,    &
     &                                 Rwork(2*Nrhs+1),result(2))
!
!                    Check solution from generated exact solution.
!
                           CALL dget04(n,Nrhs,X,lda,Xact,lda,rcondc,    &
     &                                 result(3))
!
!                    Check the error bounds from iterative refinement.
!
                           CALL dppt05(uplo,n,Nrhs,A,B,lda,X,lda,Xact,  &
     &                                 lda,Rwork,Rwork(Nrhs+1),result(4)&
     &                                 )
                        ELSE
                           k1 = 6
                        ENDIF
!
!                 Compare RCOND from DSPSVX with the computed value
!                 in RCONDC.
!
                        result(6) = dget06(rcond,rcondc)
!
!                 Print information about the tests that did not pass
!                 the threshold.
!
                        DO k = k1 , 6
                           IF ( result(k)>=Thresh ) THEN
                              IF ( nfail==0 .AND. nerrs==0 )            &
     &                             CALL aladhd(Nout,path)
                              WRITE (Nout,FMT=99002) 'DSPSVX' , fact ,  &
     &                               uplo , n , imat , k , result(k)
                              nfail = nfail + 1
                           ENDIF
                        ENDDO
                        nrun = nrun + 7 - k1
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
      CALL alasvm(path,Nout,nfail,nrun,nerrs)
!
99001 FORMAT (1X,A,', UPLO=''',A1,''', N =',I5,', type ',I2,', test ',  &
     &        I2,', ratio =',G12.5)
99002 FORMAT (1X,A,', FACT=''',A1,''', UPLO=''',A1,''', N =',I5,        &
     &        ', type ',I2,', test ',I2,', ratio =',G12.5)
!
!     End of ddrvsp
!
      END SUBROUTINE ddrvsp
