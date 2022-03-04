!*==dchkpt.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b DCHKPT
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DCHKPT( DOTYPE, NN, NVAL, NNS, NSVAL, THRESH, TSTERR,
!                          A, D, E, B, X, XACT, WORK, RWORK, NOUT )
!
!       .. Scalar Arguments ..
!       LOGICAL            TSTERR
!       INTEGER            NN, NNS, NOUT
!       DOUBLE PRECISION   THRESH
!       ..
!       .. Array Arguments ..
!       LOGICAL            DOTYPE( * )
!       INTEGER            NSVAL( * ), NVAL( * )
!       DOUBLE PRECISION   A( * ), B( * ), D( * ), E( * ), RWORK( * ),
!      $                   WORK( * ), X( * ), XACT( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DCHKPT tests DPTTRF, -TRS, -RFS, and -CON
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
!> \param[out] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (NMAX*2)
!> \endverbatim
!>
!> \param[out] D
!> \verbatim
!>          D is DOUBLE PRECISION array, dimension (NMAX*2)
!> \endverbatim
!>
!> \param[out] E
!> \verbatim
!>          E is DOUBLE PRECISION array, dimension (NMAX*2)
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
!>          WORK is DOUBLE PRECISION array, dimension
!>                      (NMAX*max(3,NSMAX))
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is DOUBLE PRECISION array, dimension
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
!> \ingroup double_lin
!
!  =====================================================================
      SUBROUTINE DCHKPT(Dotype,Nn,Nval,Nns,Nsval,Thresh,Tsterr,A,D,E,B, &
     &                  X,Xact,Work,Rwork,Nout)
      IMPLICIT NONE
!*--DCHKPT150
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      LOGICAL Tsterr
      INTEGER Nn , Nns , Nout
      DOUBLE PRECISION Thresh
!     ..
!     .. Array Arguments ..
      LOGICAL Dotype(*)
      INTEGER Nsval(*) , Nval(*)
      DOUBLE PRECISION A(*) , B(*) , D(*) , E(*) , Rwork(*) , Work(*) , &
     &                 X(*) , Xact(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ONE , ZERO
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
      INTEGER NTYPES
      PARAMETER (NTYPES=12)
      INTEGER NTESTS
      PARAMETER (NTESTS=7)
!     ..
!     .. Local Scalars ..
      LOGICAL zerot
      CHARACTER dist , type
      CHARACTER*3 path
      INTEGER i , ia , imat , in , info , irhs , ix , izero , j , k ,   &
     &        kl , ku , lda , mode , n , nerrs , nfail , nimat , nrhs , &
     &        nrun
      DOUBLE PRECISION ainvnm , anorm , cond , dmax , rcond , rcondc
!     ..
!     .. Local Arrays ..
      INTEGER iseed(4) , iseedy(4)
      DOUBLE PRECISION result(NTESTS) , z(3)
!     ..
!     .. External Functions ..
      INTEGER IDAMAX
      DOUBLE PRECISION DASUM , DGET06 , DLANST
      EXTERNAL IDAMAX , DASUM , DGET06 , DLANST
!     ..
!     .. External Subroutines ..
      EXTERNAL ALAERH , ALAHD , ALASUM , DCOPY , DERRGT , DGET04 ,      &
     &         DLACPY , DLAPTM , DLARNV , DLATB4 , DLATMS , DPTCON ,    &
     &         DPTRFS , DPTT01 , DPTT02 , DPTT05 , DPTTRF , DPTTRS ,    &
     &         DSCAL
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , MAX
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
      path(1:1) = 'Double precision'
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
      IF ( Tsterr ) CALL DERRGT(path,Nout)
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
!           Set up parameters with DLATB4.
!
               CALL DLATB4(path,imat,n,n,type,kl,ku,anorm,mode,cond,    &
     &                     dist)
!
               zerot = imat>=8 .AND. imat<=10
               IF ( imat<=6 ) THEN
!
!              Type 1-6:  generate a symmetric tridiagonal matrix of
!              known condition number in lower triangular band storage.
!
                  SRNamt = 'DLATMS'
                  CALL DLATMS(n,n,dist,iseed,type,Rwork,mode,cond,anorm,&
     &                        kl,ku,'B',A,2,Work,info)
!
!              Check the error code from DLATMS.
!
                  IF ( info/=0 ) THEN
                     CALL ALAERH(path,'DLATMS',info,0,' ',n,n,kl,ku,-1, &
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
                     CALL DLARNV(2,iseed,n,D)
                     CALL DLARNV(2,iseed,n-1,E)
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
                     ix = IDAMAX(n,D,1)
                     dmax = D(ix)
                     CALL DSCAL(n,anorm/dmax,D,1)
                     CALL DSCAL(n-1,anorm/dmax,E,1)
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
               CALL DCOPY(n,D,1,D(n+1),1)
               IF ( n>1 ) CALL DCOPY(n-1,E,1,E(n+1),1)
!
!+    TEST 1
!           Factor A as L*D*L' and compute the ratio
!              norm(L*D*L' - A) / (n * norm(A) * EPS )
!
               CALL DPTTRF(n,D(n+1),E(n+1),info)
!
!           Check error code from DPTTRF.
!
               IF ( info/=izero ) THEN
                  CALL ALAERH(path,'DPTTRF',info,izero,' ',n,n,-1,-1,-1,&
     &                        imat,nfail,nerrs,Nout)
                  CYCLE
               ENDIF
!
               IF ( info>0 ) THEN
                  rcondc = ZERO
                  GOTO 10
               ENDIF
!
               CALL DPTT01(n,D,E,D(n+1),E(n+1),Work,result(1))
!
!           Print the test ratio if greater than or equal to THRESH.
!
               IF ( result(1)>=Thresh ) THEN
                  IF ( nfail==0 .AND. nerrs==0 ) CALL ALAHD(Nout,path)
                  WRITE (Nout,FMT=99001) n , imat , 1 , result(1)
                  nfail = nfail + 1
               ENDIF
               nrun = nrun + 1
!
!           Compute RCONDC = 1 / (norm(A) * norm(inv(A))
!
!           Compute norm(A).
!
               anorm = DLANST('1',n,D,E)
!
!           Use DPTTRS to solve for one column at a time of inv(A),
!           computing the maximum column sum as we go.
!
               ainvnm = ZERO
               DO i = 1 , n
                  DO j = 1 , n
                     X(j) = ZERO
                  ENDDO
                  X(i) = ONE
                  CALL DPTTRS(n,1,D(n+1),E(n+1),X,lda,info)
                  ainvnm = MAX(ainvnm,DASUM(n,X,1))
               ENDDO
               rcondc = ONE/MAX(ONE,anorm*ainvnm)
!
               DO irhs = 1 , Nns
                  nrhs = Nsval(irhs)
!
!           Generate NRHS random solution vectors.
!
                  ix = 1
                  DO j = 1 , nrhs
                     CALL DLARNV(2,iseed,n,Xact(ix))
                     ix = ix + lda
                  ENDDO
!
!           Set the right hand side.
!
                  CALL DLAPTM(n,nrhs,ONE,D,E,Xact,lda,ZERO,B,lda)
!
!+    TEST 2
!           Solve A*x = b and compute the residual.
!
                  CALL DLACPY('Full',n,nrhs,B,lda,X,lda)
                  CALL DPTTRS(n,nrhs,D(n+1),E(n+1),X,lda,info)
!
!           Check error code from DPTTRS.
!
                  IF ( info/=0 ) CALL ALAERH(path,'DPTTRS',info,0,' ',n,&
     &                 n,-1,-1,nrhs,imat,nfail,nerrs,Nout)
!
                  CALL DLACPY('Full',n,nrhs,B,lda,Work,lda)
                  CALL DPTT02(n,nrhs,D,E,X,lda,Work,lda,result(2))
!
!+    TEST 3
!           Check solution from generated exact solution.
!
                  CALL DGET04(n,nrhs,X,lda,Xact,lda,rcondc,result(3))
!
!+    TESTS 4, 5, and 6
!           Use iterative refinement to improve the solution.
!
                  SRNamt = 'DPTRFS'
                  CALL DPTRFS(n,nrhs,D,E,D(n+1),E(n+1),B,lda,X,lda,     &
     &                        Rwork,Rwork(nrhs+1),Work,info)
!
!           Check error code from DPTRFS.
!
                  IF ( info/=0 ) CALL ALAERH(path,'DPTRFS',info,0,' ',n,&
     &                 n,-1,-1,nrhs,imat,nfail,nerrs,Nout)
!
                  CALL DGET04(n,nrhs,X,lda,Xact,lda,rcondc,result(4))
                  CALL DPTT05(n,nrhs,D,E,B,lda,X,lda,Xact,lda,Rwork,    &
     &                        Rwork(nrhs+1),result(5))
!
!           Print information about the tests that did not pass the
!           threshold.
!
                  DO k = 2 , 6
                     IF ( result(k)>=Thresh ) THEN
                        IF ( nfail==0 .AND. nerrs==0 )                  &
     &                       CALL ALAHD(Nout,path)
                        WRITE (Nout,FMT=99002) n , nrhs , imat , k ,    &
     &                         result(k)
                        nfail = nfail + 1
                     ENDIF
                  ENDDO
                  nrun = nrun + 5
               ENDDO
!
!+    TEST 7
!           Estimate the reciprocal of the condition number of the
!           matrix.
!
 10            SRNamt = 'DPTCON'
               CALL DPTCON(n,D(n+1),E(n+1),anorm,rcond,Rwork,info)
!
!           Check error code from DPTCON.
!
               IF ( info/=0 ) CALL ALAERH(path,'DPTCON',info,0,' ',n,n, &
     &              -1,-1,-1,imat,nfail,nerrs,Nout)
!
               result(7) = DGET06(rcond,rcondc)
!
!           Print the test ratio if greater than or equal to THRESH.
!
               IF ( result(7)>=Thresh ) THEN
                  IF ( nfail==0 .AND. nerrs==0 ) CALL ALAHD(Nout,path)
                  WRITE (Nout,FMT=99001) n , imat , 7 , result(7)
                  nfail = nfail + 1
               ENDIF
               nrun = nrun + 1
            ENDIF
         ENDDO
      ENDDO
!
!     Print a summary of the results.
!
      CALL ALASUM(path,Nout,nfail,nrun,nerrs)
!
99001 FORMAT (' N =',I5,', type ',I2,', test ',I2,', ratio = ',G12.5)
99002 FORMAT (' N =',I5,', NRHS=',I3,', type ',I2,', test(',I2,') = ',  &
     &        G12.5)
!
!     End of DCHKPT
!
      END SUBROUTINE DCHKPT
