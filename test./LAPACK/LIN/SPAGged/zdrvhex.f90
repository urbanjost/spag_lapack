!*==zdrvhe.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b ZDRVHEX
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZDRVHE( DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, NMAX,
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
!> ZDRVHE tests the driver routines ZHESV, -SVX, and -SVXX.
!>
!> Note that this file is used only when the XBLAS are available,
!> otherwise zdrvhe.f defines this subroutine.
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
!>          B is COMPLEX*16 array, dimension (NMAX*NRHS)
!> \endverbatim
!>
!> \param[out] X
!> \verbatim
!>          X is COMPLEX*16 array, dimension (NMAX*NRHS)
!> \endverbatim
!>
!> \param[out] XACT
!> \verbatim
!>          XACT is COMPLEX*16 array, dimension (NMAX*NRHS)
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX*16 array, dimension
!>                      (NMAX*max(2,NRHS))
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is DOUBLE PRECISION array, dimension (2*NMAX+2*NRHS)
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
!> \date April 2012
!
!> \ingroup complex16_lin
!
!  =====================================================================
      SUBROUTINE ZDRVHE(Dotype,Nn,Nval,Nrhs,Thresh,Tsterr,Nmax,A,Afac,  &
     &                  Ainv,B,X,Xact,Work,Rwork,Iwork,Nout)
      IMPLICIT NONE
!*--ZDRVHE160
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     April 2012
!
!     .. Scalar Arguments ..
      LOGICAL Tsterr
      INTEGER Nmax , Nn , Nout , Nrhs
      DOUBLE PRECISION Thresh
!     ..
!     .. Array Arguments ..
      LOGICAL Dotype(*)
      INTEGER Iwork(*) , Nval(*)
      DOUBLE PRECISION Rwork(*)
      COMPLEX*16 A(*) , Afac(*) , Ainv(*) , B(*) , Work(*) , X(*) ,     &
     &           Xact(*)
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
      CHARACTER dist , equed , fact , type , uplo , xtype
      CHARACTER*3 path
      INTEGER i , i1 , i2 , ifact , imat , in , info , ioff , iuplo ,   &
     &        izero , j , k , k1 , kl , ku , lda , lwork , mode , n ,   &
     &        nb , nbmin , nerrs , nfail , nimat , nrun , nt ,          &
     &        n_err_bnds
      DOUBLE PRECISION ainvnm , anorm , cndnum , rcond , rcondc ,       &
     &                 rpvgrw_svxx
!     ..
!     .. Local Arrays ..
      CHARACTER facts(NFACT) , uplos(2)
      INTEGER iseed(4) , iseedy(4)
      DOUBLE PRECISION result(NTESTS) , berr(Nrhs) , errbnds_n(Nrhs,3) ,&
     &                 errbnds_c(Nrhs,3)
!     ..
!     .. External Functions ..
      DOUBLE PRECISION DGET06 , ZLANHE
      EXTERNAL DGET06 , ZLANHE
!     ..
!     .. External Subroutines ..
      EXTERNAL ALADHD , ALAERH , ALASVM , XLAENV , ZERRVX , ZGET04 ,    &
     &         ZHESV , ZHESVX , ZHET01 , ZHETRF , ZHETRI2 , ZLACPY ,    &
     &         ZLAIPD , ZLARHS , ZLASET , ZLATB4 , ZLATMS , ZPOT02 ,    &
     &         ZPOT05 , ZHESVXX
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
      INTRINSIC DCMPLX , MAX , MIN
!     ..
!     .. Data statements ..
      DATA iseedy/1988 , 1989 , 1990 , 1991/
      DATA uplos/'U' , 'L'/ , facts/'F' , 'N'/
!     ..
!     .. Executable Statements ..
!
!     Initialize constants and the random number seed.
!
      path(1:1) = 'Z'
      path(2:3) = 'HE'
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
      IF ( Tsterr ) CALL ZERRVX(path,Nout)
      INFot = 0
!
!     Set the block size and minimum block size for testing.
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
!              Set up parameters with ZLATB4 and generate a test matrix
!              with ZLATMS.
!
                     CALL ZLATB4(path,imat,n,n,type,kl,ku,anorm,mode,   &
     &                           cndnum,dist)
!
                     SRNamt = 'ZLATMS'
                     CALL ZLATMS(n,n,dist,iseed,type,Rwork,mode,cndnum, &
     &                           anorm,kl,ku,uplo,A,lda,Work,info)
!
!              Check error code from ZLATMS.
!
                     IF ( info/=0 ) THEN
                        CALL ALAERH(path,'ZLATMS',info,0,uplo,n,n,-1,-1,&
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
!              Set the imaginary part of the diagonals.
!
                     CALL ZLAIPD(n,A,lda+1,0)
!
                     DO ifact = 1 , NFACT
!
!                 Do first for FACT = 'F', then for other values.
!
                        fact = facts(ifact)
!
!                 Compute the condition number for comparison with
!                 the value returned by ZHESVX.
!
                        IF ( zerot ) THEN
                           IF ( ifact==1 ) CYCLE
                           rcondc = ZERO
!
                        ELSEIF ( ifact==1 ) THEN
!
!                    Compute the 1-norm of A.
!
                           anorm = ZLANHE('1',uplo,n,A,lda,Rwork)
!
!                    Factor the matrix A.
!
                           CALL ZLACPY(uplo,n,n,A,lda,Afac,lda)
                           CALL ZHETRF(uplo,n,Afac,lda,Iwork,Work,lwork,&
     &                                 info)
!
!                    Compute inv(A) and take its norm.
!
                           CALL ZLACPY(uplo,n,n,Afac,lda,Ainv,lda)
                           lwork = (n+nb+1)*(nb+3)
                           CALL ZHETRI2(uplo,n,Ainv,lda,Iwork,Work,     &
     &                                  lwork,info)
                           ainvnm = ZLANHE('1',uplo,n,Ainv,lda,Rwork)
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
                        SRNamt = 'ZLARHS'
                        CALL ZLARHS(path,xtype,uplo,' ',n,n,kl,ku,Nrhs, &
     &                              A,lda,Xact,lda,B,lda,iseed,info)
                        xtype = 'C'
!
!                 --- Test ZHESV  ---
!
                        IF ( ifact==2 ) THEN
                           CALL ZLACPY(uplo,n,n,A,lda,Afac,lda)
                           CALL ZLACPY('Full',n,Nrhs,B,lda,X,lda)
!
!                    Factor the matrix and solve the system using ZHESV.
!
                           SRNamt = 'ZHESV '
                           CALL ZHESV(uplo,n,Nrhs,Afac,lda,Iwork,X,lda, &
     &                                Work,lwork,info)
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
!                    Check error code from ZHESV .
!
                           IF ( info/=k ) THEN
                              CALL ALAERH(path,'ZHESV ',info,k,uplo,n,n,&
     &                           -1,-1,Nrhs,imat,nfail,nerrs,Nout)
                              GOTO 2
                           ELSEIF ( info/=0 ) THEN
                              GOTO 2
                           ENDIF
!
!                    Reconstruct matrix from factors and compute
!                    residual.
!
                           CALL ZHET01(uplo,n,A,lda,Afac,lda,Iwork,Ainv,&
     &                                 lda,Rwork,result(1))
!
!                    Compute residual of the computed solution.
!
                           CALL ZLACPY('Full',n,Nrhs,B,lda,Work,lda)
                           CALL ZPOT02(uplo,n,Nrhs,A,lda,X,lda,Work,lda,&
     &                                 Rwork,result(2))
!
!                    Check solution from generated exact solution.
!
                           CALL ZGET04(n,Nrhs,X,lda,Xact,lda,rcondc,    &
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
                                 WRITE (Nout,FMT=99001) 'ZHESV ' ,      &
     &                                  uplo , n , imat , k , result(k)
                                 nfail = nfail + 1
                              ENDIF
                           ENDDO
                           nrun = nrun + nt
                        ENDIF
!
!                 --- Test ZHESVX ---
!
 2                      IF ( ifact==2 )                                 &
     &                       CALL ZLASET(uplo,n,n,DCMPLX(ZERO),         &
     &                       DCMPLX(ZERO),Afac,lda)
                        CALL ZLASET('Full',n,Nrhs,DCMPLX(ZERO),         &
     &                              DCMPLX(ZERO),X,lda)
!
!                 Solve the system and compute the condition number and
!                 error bounds using ZHESVX.
!
                        SRNamt = 'ZHESVX'
                        CALL ZHESVX(fact,uplo,n,Nrhs,A,lda,Afac,lda,    &
     &                              Iwork,B,lda,X,lda,rcond,Rwork,      &
     &                              Rwork(Nrhs+1),Work,lwork,           &
     &                              Rwork(2*Nrhs+1),info)
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
!                 Check the error code from ZHESVX.
!
                        IF ( info/=k ) THEN
                           CALL ALAERH(path,'ZHESVX',info,k,fact//uplo, &
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
                              CALL ZHET01(uplo,n,A,lda,Afac,lda,Iwork,  &
     &                           Ainv,lda,Rwork(2*Nrhs+1),result(1))
                              k1 = 1
                           ELSE
                              k1 = 2
                           ENDIF
!
!                    Compute residual of the computed solution.
!
                           CALL ZLACPY('Full',n,Nrhs,B,lda,Work,lda)
                           CALL ZPOT02(uplo,n,Nrhs,A,lda,X,lda,Work,lda,&
     &                                 Rwork(2*Nrhs+1),result(2))
!
!                    Check solution from generated exact solution.
!
                           CALL ZGET04(n,Nrhs,X,lda,Xact,lda,rcondc,    &
     &                                 result(3))
!
!                    Check the error bounds from iterative refinement.
!
                           CALL ZPOT05(uplo,n,Nrhs,A,lda,B,lda,X,lda,   &
     &                                 Xact,lda,Rwork,Rwork(Nrhs+1),    &
     &                                 result(4))
                        ELSE
                           k1 = 6
                        ENDIF
!
!                 Compare RCOND from ZHESVX with the computed value
!                 in RCONDC.
!
                        result(6) = DGET06(rcond,rcondc)
!
!                 Print information about the tests that did not pass
!                 the threshold.
!
                        DO k = k1 , 6
                           IF ( result(k)>=Thresh ) THEN
                              IF ( nfail==0 .AND. nerrs==0 )            &
     &                             CALL ALADHD(Nout,path)
                              WRITE (Nout,FMT=99002) 'ZHESVX' , fact ,  &
     &                               uplo , n , imat , k , result(k)
                              nfail = nfail + 1
                           ENDIF
                        ENDDO
                        nrun = nrun + 7 - k1
!
!                 --- Test ZHESVXX ---
!
!                 Restore the matrices A and B.
!
                        IF ( ifact==2 )                                 &
     &                       CALL ZLASET(uplo,n,n,CMPLX(ZERO),          &
     &                       CMPLX(ZERO),Afac,lda)
                        CALL ZLASET('Full',n,Nrhs,CMPLX(ZERO),          &
     &                              CMPLX(ZERO),X,lda)
!
!                 Solve the system and compute the condition number
!                 and error bounds using ZHESVXX.
!
                        SRNamt = 'ZHESVXX'
                        n_err_bnds = 3
                        equed = 'N'
                        CALL ZHESVXX(fact,uplo,n,Nrhs,A,lda,Afac,lda,   &
     &                               Iwork,equed,Work(n+1),B,lda,X,lda, &
     &                               rcond,rpvgrw_svxx,berr,n_err_bnds, &
     &                               errbnds_n,errbnds_c,0,ZERO,Work,   &
     &                               Rwork(2*Nrhs+1),info)
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
!                 Check the error code from ZHESVXX.
!
                        IF ( info/=k .AND. info<=n ) THEN
                           CALL ALAERH(path,'ZHESVXX',info,k,fact//uplo,&
     &                                 n,n,-1,-1,Nrhs,imat,nfail,nerrs, &
     &                                 Nout)
                           CYCLE
                        ENDIF
!
                        IF ( info==0 ) THEN
                           IF ( ifact>=2 ) THEN
!
!                 Reconstruct matrix from factors and compute
!                 residual.
!
                              CALL ZHET01(uplo,n,A,lda,Afac,lda,Iwork,  &
     &                           Ainv,lda,Rwork(2*Nrhs+1),result(1))
                              k1 = 1
                           ELSE
                              k1 = 2
                           ENDIF
!
!                 Compute residual of the computed solution.
!
                           CALL ZLACPY('Full',n,Nrhs,B,lda,Work,lda)
                           CALL ZPOT02(uplo,n,Nrhs,A,lda,X,lda,Work,lda,&
     &                                 Rwork(2*Nrhs+1),result(2))
                           result(2) = 0.0
!
!                 Check solution from generated exact solution.
!
                           CALL ZGET04(n,Nrhs,X,lda,Xact,lda,rcondc,    &
     &                                 result(3))
!
!                 Check the error bounds from iterative refinement.
!
                           CALL ZPOT05(uplo,n,Nrhs,A,lda,B,lda,X,lda,   &
     &                                 Xact,lda,Rwork,Rwork(Nrhs+1),    &
     &                                 result(4))
                        ELSE
                           k1 = 6
                        ENDIF
!
!                 Compare RCOND from ZHESVXX with the computed value
!                 in RCONDC.
!
                        result(6) = DGET06(rcond,rcondc)
!
!                 Print information about the tests that did not pass
!                 the threshold.
!
                        DO k = k1 , 6
                           IF ( result(k)>=Thresh ) THEN
                              IF ( nfail==0 .AND. nerrs==0 )            &
     &                             CALL ALADHD(Nout,path)
                              WRITE (Nout,FMT=99002) 'ZHESVXX' , fact , &
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
      CALL ALASVM(path,Nout,nfail,nrun,nerrs)
!
 
!     Test Error Bounds from ZHESVXX
 
      CALL ZEBCHVXX(Thresh,path)
 
99001 FORMAT (1X,A,', UPLO=''',A1,''', N =',I5,', type ',I2,', test ',  &
     &        I2,', ratio =',G12.5)
99002 FORMAT (1X,A,', FACT=''',A1,''', UPLO=''',A1,''', N =',I5,        &
     &        ', type ',I2,', test ',I2,', ratio =',G12.5)
!
!     End of ZDRVHE
!
      END SUBROUTINE ZDRVHE
