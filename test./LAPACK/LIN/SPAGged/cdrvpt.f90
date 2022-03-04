!*==cdrvpt.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b CDRVPT
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE CDRVPT( DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, A, D,
!                          E, B, X, XACT, WORK, RWORK, NOUT )
!
!       .. Scalar Arguments ..
!       LOGICAL            TSTERR
!       INTEGER            NN, NOUT, NRHS
!       REAL               THRESH
!       ..
!       .. Array Arguments ..
!       LOGICAL            DOTYPE( * )
!       INTEGER            NVAL( * )
!       REAL               D( * ), RWORK( * )
!       COMPLEX            A( * ), B( * ), E( * ), WORK( * ), X( * ),
!      $                   XACT( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CDRVPT tests CPTSV and -SVX.
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
!> \param[out] A
!> \verbatim
!>          A is COMPLEX array, dimension (NMAX*2)
!> \endverbatim
!>
!> \param[out] D
!> \verbatim
!>          D is REAL array, dimension (NMAX*2)
!> \endverbatim
!>
!> \param[out] E
!> \verbatim
!>          E is COMPLEX array, dimension (NMAX*2)
!> \endverbatim
!>
!> \param[out] B
!> \verbatim
!>          B is COMPLEX array, dimension (NMAX*NRHS)
!> \endverbatim
!>
!> \param[out] X
!> \verbatim
!>          X is COMPLEX array, dimension (NMAX*NRHS)
!> \endverbatim
!>
!> \param[out] XACT
!> \verbatim
!>          XACT is COMPLEX array, dimension (NMAX*NRHS)
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX array, dimension
!>                      (NMAX*max(3,NRHS))
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is REAL array, dimension (NMAX+2*NRHS)
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
      SUBROUTINE CDRVPT(Dotype,Nn,Nval,Nrhs,Thresh,Tsterr,A,D,E,B,X,    &
     &                  Xact,Work,Rwork,Nout)
      IMPLICIT NONE
!*--CDRVPT144
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      LOGICAL Tsterr
      INTEGER Nn , Nout , Nrhs
      REAL Thresh
!     ..
!     .. Array Arguments ..
      LOGICAL Dotype(*)
      INTEGER Nval(*)
      REAL D(*) , Rwork(*)
      COMPLEX A(*) , B(*) , E(*) , Work(*) , X(*) , Xact(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ONE , ZERO
      PARAMETER (ONE=1.0E+0,ZERO=0.0E+0)
      INTEGER NTYPES
      PARAMETER (NTYPES=12)
      INTEGER NTESTS
      PARAMETER (NTESTS=6)
!     ..
!     .. Local Scalars ..
      LOGICAL zerot
      CHARACTER dist , fact , type
      CHARACTER*3 path
      INTEGER i , ia , ifact , imat , in , info , ix , izero , j , k ,  &
     &        k1 , kl , ku , lda , mode , n , nerrs , nfail , nimat ,   &
     &        nrun , nt
      REAL ainvnm , anorm , cond , dmax , rcond , rcondc
!     ..
!     .. Local Arrays ..
      INTEGER iseed(4) , iseedy(4)
      REAL result(NTESTS) , z(3)
!     ..
!     .. External Functions ..
      INTEGER ISAMAX
      REAL CLANHT , SCASUM , SGET06
      EXTERNAL ISAMAX , CLANHT , SCASUM , SGET06
!     ..
!     .. External Subroutines ..
      EXTERNAL ALADHD , ALAERH , ALASVM , CCOPY , CERRVX , CGET04 ,     &
     &         CLACPY , CLAPTM , CLARNV , CLASET , CLATB4 , CLATMS ,    &
     &         CPTSV , CPTSVX , CPTT01 , CPTT02 , CPTT05 , CPTTRF ,     &
     &         CPTTRS , CSSCAL , SCOPY , SLARNV , SSCAL
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , CMPLX , MAX
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
      DATA iseedy/0 , 0 , 0 , 1/
!     ..
!     .. Executable Statements ..
!
      path(1:1) = 'Complex precision'
      path(2:3) = 'PT'
      nrun = 0
      nfail = 0
      nerrs = 0
      DO i = 1 , 4
         iseed(i) = iseedy(i)
      ENDDO
!
!     Test the error exits
!
      IF ( Tsterr ) CALL CERRVX(path,Nout)
      INFot = 0
!
      DO in = 1 , Nn
!
!        Do for each value of N in NVAL.
!
         n = Nval(in)
         lda = MAX(1,n)
         nimat = NTYPES
         IF ( n<=0 ) nimat = 1
!
         DO imat = 1 , nimat
!
!           Do the tests only if DOTYPE( IMAT ) is true.
!
            IF ( .NOT.(n>0 .AND. .NOT.Dotype(imat)) ) THEN
!
!           Set up parameters with CLATB4.
!
               CALL CLATB4(path,imat,n,n,type,kl,ku,anorm,mode,cond,    &
     &                     dist)
!
               zerot = imat>=8 .AND. imat<=10
               IF ( imat<=6 ) THEN
!
!              Type 1-6:  generate a symmetric tridiagonal matrix of
!              known condition number in lower triangular band storage.
!
                  SRNamt = 'CLATMS'
                  CALL CLATMS(n,n,dist,iseed,type,Rwork,mode,cond,anorm,&
     &                        kl,ku,'B',A,2,Work,info)
!
!              Check the error code from CLATMS.
!
                  IF ( info/=0 ) THEN
                     CALL ALAERH(path,'CLATMS',info,0,' ',n,n,kl,ku,-1, &
     &                           imat,nfail,nerrs,Nout)
                     CYCLE
                  ENDIF
                  izero = 0
!
!              Copy the matrix to D and E.
!
                  ia = 1
                  DO i = 1 , n - 1
                     D(i) = A(ia)
                     E(i) = A(ia+1)
                     ia = ia + 2
                  ENDDO
                  IF ( n>0 ) D(n) = A(ia)
               ELSE
!
!              Type 7-12:  generate a diagonally dominant matrix with
!              unknown condition number in the vectors D and E.
!
                  IF ( .NOT.zerot .OR. .NOT.Dotype(7) ) THEN
!
!                 Let D and E have values from [-1,1].
!
                     CALL SLARNV(2,iseed,n,D)
                     CALL CLARNV(2,iseed,n-1,E)
!
!                 Make the tridiagonal matrix diagonally dominant.
!
                     IF ( n==1 ) THEN
                        D(1) = ABS(D(1))
                     ELSE
                        D(1) = ABS(D(1)) + ABS(E(1))
                        D(n) = ABS(D(n)) + ABS(E(n-1))
                        DO i = 2 , n - 1
                           D(i) = ABS(D(i)) + ABS(E(i)) + ABS(E(i-1))
                        ENDDO
                     ENDIF
!
!                 Scale D and E so the maximum element is ANORM.
!
                     ix = ISAMAX(n,D,1)
                     dmax = D(ix)
                     CALL SSCAL(n,anorm/dmax,D,1)
                     IF ( n>1 ) CALL CSSCAL(n-1,anorm/dmax,E,1)
!
                  ELSEIF ( izero>0 ) THEN
!
!                 Reuse the last matrix by copying back the zeroed out
!                 elements.
!
                     IF ( izero==1 ) THEN
                        D(1) = z(2)
                        IF ( n>1 ) E(1) = z(3)
                     ELSEIF ( izero==n ) THEN
                        E(n-1) = z(1)
                        D(n) = z(2)
                     ELSE
                        E(izero-1) = z(1)
                        D(izero) = z(2)
                        E(izero) = z(3)
                     ENDIF
                  ENDIF
!
!              For types 8-10, set one row and column of the matrix to
!              zero.
!
                  izero = 0
                  IF ( imat==8 ) THEN
                     izero = 1
                     z(2) = D(1)
                     D(1) = ZERO
                     IF ( n>1 ) THEN
                        z(3) = E(1)
                        E(1) = ZERO
                     ENDIF
                  ELSEIF ( imat==9 ) THEN
                     izero = n
                     IF ( n>1 ) THEN
                        z(1) = E(n-1)
                        E(n-1) = ZERO
                     ENDIF
                     z(2) = D(n)
                     D(n) = ZERO
                  ELSEIF ( imat==10 ) THEN
                     izero = (n+1)/2
                     IF ( izero>1 ) THEN
                        z(1) = E(izero-1)
                        E(izero-1) = ZERO
                        z(3) = E(izero)
                        E(izero) = ZERO
                     ENDIF
                     z(2) = D(izero)
                     D(izero) = ZERO
                  ENDIF
               ENDIF
!
!           Generate NRHS random solution vectors.
!
               ix = 1
               DO j = 1 , Nrhs
                  CALL CLARNV(2,iseed,n,Xact(ix))
                  ix = ix + lda
               ENDDO
!
!           Set the right hand side.
!
               CALL CLAPTM('Lower',n,Nrhs,ONE,D,E,Xact,lda,ZERO,B,lda)
!
               DO ifact = 1 , 2
                  IF ( ifact==1 ) THEN
                     fact = 'F'
                  ELSE
                     fact = 'N'
                  ENDIF
!
!              Compute the condition number for comparison with
!              the value returned by CPTSVX.
!
                  IF ( zerot ) THEN
                     IF ( ifact==1 ) CYCLE
                     rcondc = ZERO
!
                  ELSEIF ( ifact==1 ) THEN
!
!                 Compute the 1-norm of A.
!
                     anorm = CLANHT('1',n,D,E)
!
                     CALL SCOPY(n,D,1,D(n+1),1)
                     IF ( n>1 ) CALL CCOPY(n-1,E,1,E(n+1),1)
!
!                 Factor the matrix A.
!
                     CALL CPTTRF(n,D(n+1),E(n+1),info)
!
!                 Use CPTTRS to solve for one column at a time of
!                 inv(A), computing the maximum column sum as we go.
!
                     ainvnm = ZERO
                     DO i = 1 , n
                        DO j = 1 , n
                           X(j) = ZERO
                        ENDDO
                        X(i) = ONE
                        CALL CPTTRS('Lower',n,1,D(n+1),E(n+1),X,lda,    &
     &                              info)
                        ainvnm = MAX(ainvnm,SCASUM(n,X,1))
                     ENDDO
!
!                 Compute the 1-norm condition number of A.
!
                     IF ( anorm<=ZERO .OR. ainvnm<=ZERO ) THEN
                        rcondc = ONE
                     ELSE
                        rcondc = (ONE/anorm)/ainvnm
                     ENDIF
                  ENDIF
!
                  IF ( ifact==2 ) THEN
!
!                 --- Test CPTSV --
!
                     CALL SCOPY(n,D,1,D(n+1),1)
                     IF ( n>1 ) CALL CCOPY(n-1,E,1,E(n+1),1)
                     CALL CLACPY('Full',n,Nrhs,B,lda,X,lda)
!
!                 Factor A as L*D*L' and solve the system A*X = B.
!
                     SRNamt = 'CPTSV '
                     CALL CPTSV(n,Nrhs,D(n+1),E(n+1),X,lda,info)
!
!                 Check error code from CPTSV .
!
                     IF ( info/=izero )                                 &
     &                    CALL ALAERH(path,'CPTSV ',info,izero,' ',n,n, &
     &                    1,1,Nrhs,imat,nfail,nerrs,Nout)
                     nt = 0
                     IF ( izero==0 ) THEN
!
!                    Check the factorization by computing the ratio
!                       norm(L*D*L' - A) / (n * norm(A) * EPS )
!
                        CALL CPTT01(n,D,E,D(n+1),E(n+1),Work,result(1))
!
!                    Compute the residual in the solution.
!
                        CALL CLACPY('Full',n,Nrhs,B,lda,Work,lda)
                        CALL CPTT02('Lower',n,Nrhs,D,E,X,lda,Work,lda,  &
     &                              result(2))
!
!                    Check solution from generated exact solution.
!
                        CALL CGET04(n,Nrhs,X,lda,Xact,lda,rcondc,       &
     &                              result(3))
                        nt = 3
                     ENDIF
!
!                 Print information about the tests that did not pass
!                 the threshold.
!
                     DO k = 1 , nt
                        IF ( result(k)>=Thresh ) THEN
                           IF ( nfail==0 .AND. nerrs==0 )               &
     &                          CALL ALADHD(Nout,path)
                           WRITE (Nout,FMT=99001) 'CPTSV ' , n , imat , &
     &                            k , result(k)
                           nfail = nfail + 1
                        ENDIF
                     ENDDO
                     nrun = nrun + nt
                  ENDIF
!
!              --- Test CPTSVX ---
!
                  IF ( ifact>1 ) THEN
!
!                 Initialize D( N+1:2*N ) and E( N+1:2*N ) to zero.
!
                     DO i = 1 , n - 1
                        D(n+i) = ZERO
                        E(n+i) = ZERO
                     ENDDO
                     IF ( n>0 ) D(n+n) = ZERO
                  ENDIF
!
                  CALL CLASET('Full',n,Nrhs,CMPLX(ZERO),CMPLX(ZERO),X,  &
     &                        lda)
!
!              Solve the system and compute the condition number and
!              error bounds using CPTSVX.
!
                  SRNamt = 'CPTSVX'
                  CALL CPTSVX(fact,n,Nrhs,D,E,D(n+1),E(n+1),B,lda,X,lda,&
     &                        rcond,Rwork,Rwork(Nrhs+1),Work,           &
     &                        Rwork(2*Nrhs+1),info)
!
!              Check the error code from CPTSVX.
!
                  IF ( info/=izero )                                    &
     &                 CALL ALAERH(path,'CPTSVX',info,izero,fact,n,n,1, &
     &                 1,Nrhs,imat,nfail,nerrs,Nout)
                  IF ( izero==0 ) THEN
                     IF ( ifact==2 ) THEN
!
!                    Check the factorization by computing the ratio
!                       norm(L*D*L' - A) / (n * norm(A) * EPS )
!
                        k1 = 1
                        CALL CPTT01(n,D,E,D(n+1),E(n+1),Work,result(1))
                     ELSE
                        k1 = 2
                     ENDIF
!
!                 Compute the residual in the solution.
!
                     CALL CLACPY('Full',n,Nrhs,B,lda,Work,lda)
                     CALL CPTT02('Lower',n,Nrhs,D,E,X,lda,Work,lda,     &
     &                           result(2))
!
!                 Check solution from generated exact solution.
!
                     CALL CGET04(n,Nrhs,X,lda,Xact,lda,rcondc,result(3))
!
!                 Check error bounds from iterative refinement.
!
                     CALL CPTT05(n,Nrhs,D,E,B,lda,X,lda,Xact,lda,Rwork, &
     &                           Rwork(Nrhs+1),result(4))
                  ELSE
                     k1 = 6
                  ENDIF
!
!              Check the reciprocal of the condition number.
!
                  result(6) = SGET06(rcond,rcondc)
!
!              Print information about the tests that did not pass
!              the threshold.
!
                  DO k = k1 , 6
                     IF ( result(k)>=Thresh ) THEN
                        IF ( nfail==0 .AND. nerrs==0 )                  &
     &                       CALL ALADHD(Nout,path)
                        WRITE (Nout,FMT=99002) 'CPTSVX' , fact , n ,    &
     &                         imat , k , result(k)
                        nfail = nfail + 1
                     ENDIF
                  ENDDO
                  nrun = nrun + 7 - k1
               ENDDO
            ENDIF
         ENDDO
      ENDDO
!
!     Print a summary of the results.
!
      CALL ALASVM(path,Nout,nfail,nrun,nerrs)
!
99001 FORMAT (1X,A,', N =',I5,', type ',I2,', test ',I2,', ratio = ',   &
     &        G12.5)
99002 FORMAT (1X,A,', FACT=''',A1,''', N =',I5,', type ',I2,', test ',  &
     &        I2,', ratio = ',G12.5)
!
!     End of CDRVPT
!
      END SUBROUTINE CDRVPT
