!*==cchkgt.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b CCHKGT
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE CCHKGT( DOTYPE, NN, NVAL, NNS, NSVAL, THRESH, TSTERR,
!                          A, AF, B, X, XACT, WORK, RWORK, IWORK, NOUT )
!
!       .. Scalar Arguments ..
!       LOGICAL            TSTERR
!       INTEGER            NN, NNS, NOUT
!       REAL               THRESH
!       ..
!       .. Array Arguments ..
!       LOGICAL            DOTYPE( * )
!       INTEGER            IWORK( * ), NSVAL( * ), NVAL( * )
!       REAL               RWORK( * )
!       COMPLEX            A( * ), AF( * ), B( * ), WORK( * ), X( * ),
!      $                   XACT( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CCHKGT tests CGTTRF, -TRS, -RFS, and -CON
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
!> \param[out] A
!> \verbatim
!>          A is COMPLEX array, dimension (NMAX*4)
!> \endverbatim
!>
!> \param[out] AF
!> \verbatim
!>          AF is COMPLEX array, dimension (NMAX*4)
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
!>                      (max(NMAX)+2*NSMAX)
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
      SUBROUTINE CCHKGT(Dotype,Nn,Nval,Nns,Nsval,Thresh,Tsterr,A,Af,B,X,&
     &                  Xact,Work,Rwork,Iwork,Nout)
      IMPLICIT NONE
!*--CCHKGT151
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      LOGICAL Tsterr
      INTEGER Nn , Nns , Nout
      REAL Thresh
!     ..
!     .. Array Arguments ..
      LOGICAL Dotype(*)
      INTEGER Iwork(*) , Nsval(*) , Nval(*)
      REAL Rwork(*)
      COMPLEX A(*) , Af(*) , B(*) , Work(*) , X(*) , Xact(*)
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
      PARAMETER (NTESTS=7)
!     ..
!     .. Local Scalars ..
      LOGICAL trfcon , zerot
      CHARACTER dist , norm , trans , type
      CHARACTER*3 path
      INTEGER i , imat , in , info , irhs , itran , ix , izero , j , k ,&
     &        kl , koff , ku , lda , m , mode , n , nerrs , nfail ,     &
     &        nimat , nrhs , nrun
      REAL ainvnm , anorm , cond , rcond , rcondc , rcondi , rcondo
!     ..
!     .. Local Arrays ..
      CHARACTER transs(3)
      INTEGER iseed(4) , iseedy(4)
      REAL result(NTESTS)
      COMPLEX z(3)
!     ..
!     .. External Functions ..
      REAL CLANGT , SCASUM , SGET06
      EXTERNAL CLANGT , SCASUM , SGET06
!     ..
!     .. External Subroutines ..
      EXTERNAL ALAERH , ALAHD , ALASUM , CCOPY , CERRGE , CGET04 ,      &
     &         CGTCON , CGTRFS , CGTT01 , CGTT02 , CGTT05 , CGTTRF ,    &
     &         CGTTRS , CLACPY , CLAGTM , CLARNV , CLATB4 , CLATMS ,    &
     &         CSSCAL
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX
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
      DATA iseedy/0 , 0 , 0 , 1/ , transs/'N' , 'T' , 'C'/
!     ..
!     .. Executable Statements ..
!
      path(1:1) = 'Complex precision'
      path(2:3) = 'GT'
      nrun = 0
      nfail = 0
      nerrs = 0
      DO i = 1 , 4
         iseed(i) = iseedy(i)
      ENDDO
!
!     Test the error exits
!
      IF ( Tsterr ) CALL CERRGE(path,Nout)
      INFot = 0
!
      DO in = 1 , Nn
!
!        Do for each value of N in NVAL.
!
         n = Nval(in)
         m = MAX(n-1,0)
         lda = MAX(1,n)
         nimat = NTYPES
         IF ( n<=0 ) nimat = 1
!
         DO imat = 1 , nimat
!
!           Do the tests only if DOTYPE( IMAT ) is true.
!
            IF ( Dotype(imat) ) THEN
!
!           Set up parameters with CLATB4.
!
               CALL CLATB4(path,imat,n,n,type,kl,ku,anorm,mode,cond,    &
     &                     dist)
!
               zerot = imat>=8 .AND. imat<=10
               IF ( imat<=6 ) THEN
!
!              Types 1-6:  generate matrices of known condition number.
!
                  koff = MAX(2-ku,3-MAX(1,n))
                  SRNamt = 'CLATMS'
                  CALL CLATMS(n,n,dist,iseed,type,Rwork,mode,cond,anorm,&
     &                        kl,ku,'Z',Af(koff),3,Work,info)
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
                  IF ( n>1 ) THEN
                     CALL CCOPY(n-1,Af(4),3,A,1)
                     CALL CCOPY(n-1,Af(3),3,A(n+m+1),1)
                  ENDIF
                  CALL CCOPY(n,Af(2),3,A(m+1),1)
               ELSE
!
!              Types 7-12:  generate tridiagonal matrices with
!              unknown condition numbers.
!
                  IF ( .NOT.zerot .OR. .NOT.Dotype(7) ) THEN
!
!                 Generate a matrix with elements whose real and
!                 imaginary parts are from [-1,1].
!
                     CALL CLARNV(2,iseed,n+2*m,A)
                     IF ( anorm/=ONE ) CALL CSSCAL(n+2*m,anorm,A,1)
                  ELSEIF ( izero>0 ) THEN
!
!                 Reuse the last matrix by copying back the zeroed out
!                 elements.
!
                     IF ( izero==1 ) THEN
                        A(n) = z(2)
                        IF ( n>1 ) A(1) = z(3)
                     ELSEIF ( izero==n ) THEN
                        A(3*n-2) = z(1)
                        A(2*n-1) = z(2)
                     ELSE
                        A(2*n-2+izero) = z(1)
                        A(n-1+izero) = z(2)
                        A(izero) = z(3)
                     ENDIF
                  ENDIF
!
!              If IMAT > 7, set one column of the matrix to 0.
!
                  IF ( .NOT.zerot ) THEN
                     izero = 0
                  ELSEIF ( imat==8 ) THEN
                     izero = 1
                     z(2) = A(n)
                     A(n) = ZERO
                     IF ( n>1 ) THEN
                        z(3) = A(1)
                        A(1) = ZERO
                     ENDIF
                  ELSEIF ( imat==9 ) THEN
                     izero = n
                     z(1) = A(3*n-2)
                     z(2) = A(2*n-1)
                     A(3*n-2) = ZERO
                     A(2*n-1) = ZERO
                  ELSE
                     izero = (n+1)/2
                     DO i = izero , n - 1
                        A(2*n-2+i) = ZERO
                        A(n-1+i) = ZERO
                        A(i) = ZERO
                     ENDDO
                     A(3*n-2) = ZERO
                     A(2*n-1) = ZERO
                  ENDIF
               ENDIF
!
!+    TEST 1
!           Factor A as L*U and compute the ratio
!              norm(L*U - A) / (n * norm(A) * EPS )
!
               CALL CCOPY(n+2*m,A,1,Af,1)
               SRNamt = 'CGTTRF'
               CALL CGTTRF(n,Af,Af(m+1),Af(n+m+1),Af(n+2*m+1),Iwork,    &
     &                     info)
!
!           Check error code from CGTTRF.
!
               IF ( info/=izero ) CALL ALAERH(path,'CGTTRF',info,izero, &
     &              ' ',n,n,1,1,-1,imat,nfail,nerrs,Nout)
               trfcon = info/=0
!
               CALL CGTT01(n,A,A(m+1),A(n+m+1),Af,Af(m+1),Af(n+m+1),    &
     &                     Af(n+2*m+1),Iwork,Work,lda,Rwork,result(1))
!
!           Print the test ratio if it is .GE. THRESH.
!
               IF ( result(1)>=Thresh ) THEN
                  IF ( nfail==0 .AND. nerrs==0 ) CALL ALAHD(Nout,path)
                  WRITE (Nout,FMT=99001) n , imat , 1 , result(1)
                  nfail = nfail + 1
               ENDIF
               nrun = nrun + 1
!
               DO itran = 1 , 2
                  trans = transs(itran)
                  IF ( itran==1 ) THEN
                     norm = 'O'
                  ELSE
                     norm = 'I'
                  ENDIF
                  anorm = CLANGT(norm,n,A,A(m+1),A(n+m+1))
!
                  IF ( .NOT.trfcon ) THEN
!
!                 Use CGTTRS to solve for one column at a time of
!                 inv(A), computing the maximum column sum as we go.
!
                     ainvnm = ZERO
                     DO i = 1 , n
                        DO j = 1 , n
                           X(j) = ZERO
                        ENDDO
                        X(i) = ONE
                        CALL CGTTRS(trans,n,1,Af,Af(m+1),Af(n+m+1),     &
     &                              Af(n+2*m+1),Iwork,X,lda,info)
                        ainvnm = MAX(ainvnm,SCASUM(n,X,1))
                     ENDDO
!
!                 Compute RCONDC = 1 / (norm(A) * norm(inv(A))
!
                     IF ( anorm<=ZERO .OR. ainvnm<=ZERO ) THEN
                        rcondc = ONE
                     ELSE
                        rcondc = (ONE/anorm)/ainvnm
                     ENDIF
                     IF ( itran==1 ) THEN
                        rcondo = rcondc
                     ELSE
                        rcondi = rcondc
                     ENDIF
                  ELSE
                     rcondc = ZERO
                  ENDIF
!
!+    TEST 7
!              Estimate the reciprocal of the condition number of the
!              matrix.
!
                  SRNamt = 'CGTCON'
                  CALL CGTCON(norm,n,Af,Af(m+1),Af(n+m+1),Af(n+2*m+1),  &
     &                        Iwork,anorm,rcond,Work,info)
!
!              Check error code from CGTCON.
!
                  IF ( info/=0 ) CALL ALAERH(path,'CGTCON',info,0,norm, &
     &                 n,n,-1,-1,-1,imat,nfail,nerrs,Nout)
!
                  result(7) = SGET06(rcond,rcondc)
!
!              Print the test ratio if it is .GE. THRESH.
!
                  IF ( result(7)>=Thresh ) THEN
                     IF ( nfail==0 .AND. nerrs==0 )                     &
     &                    CALL ALAHD(Nout,path)
                     WRITE (Nout,FMT=99003) norm , n , imat , 7 ,       &
     &                      result(7)
                     nfail = nfail + 1
                  ENDIF
                  nrun = nrun + 1
               ENDDO
!
!           Skip the remaining tests if the matrix is singular.
!
               IF ( .NOT.(trfcon) ) THEN
!
                  DO irhs = 1 , Nns
                     nrhs = Nsval(irhs)
!
!              Generate NRHS random solution vectors.
!
                     ix = 1
                     DO j = 1 , nrhs
                        CALL CLARNV(2,iseed,n,Xact(ix))
                        ix = ix + lda
                     ENDDO
!
                     DO itran = 1 , 3
                        trans = transs(itran)
                        IF ( itran==1 ) THEN
                           rcondc = rcondo
                        ELSE
                           rcondc = rcondi
                        ENDIF
!
!                 Set the right hand side.
!
                        CALL CLAGTM(trans,n,nrhs,ONE,A,A(m+1),A(n+m+1), &
     &                              Xact,lda,ZERO,B,lda)
!
!+    TEST 2
!              Solve op(A) * X = B and compute the residual.
!
                        CALL CLACPY('Full',n,nrhs,B,lda,X,lda)
                        SRNamt = 'CGTTRS'
                        CALL CGTTRS(trans,n,nrhs,Af,Af(m+1),Af(n+m+1),  &
     &                              Af(n+2*m+1),Iwork,X,lda,info)
!
!              Check error code from CGTTRS.
!
                        IF ( info/=0 )                                  &
     &                       CALL ALAERH(path,'CGTTRS',info,0,trans,n,n,&
     &                       -1,-1,nrhs,imat,nfail,nerrs,Nout)
!
                        CALL CLACPY('Full',n,nrhs,B,lda,Work,lda)
                        CALL CGTT02(trans,n,nrhs,A,A(m+1),A(n+m+1),X,   &
     &                              lda,Work,lda,result(2))
!
!+    TEST 3
!              Check solution from generated exact solution.
!
                        CALL CGET04(n,nrhs,X,lda,Xact,lda,rcondc,       &
     &                              result(3))
!
!+    TESTS 4, 5, and 6
!              Use iterative refinement to improve the solution.
!
                        SRNamt = 'CGTRFS'
                        CALL CGTRFS(trans,n,nrhs,A,A(m+1),A(n+m+1),Af,  &
     &                              Af(m+1),Af(n+m+1),Af(n+2*m+1),Iwork,&
     &                              B,lda,X,lda,Rwork,Rwork(nrhs+1),    &
     &                              Work,Rwork(2*nrhs+1),info)
!
!              Check error code from CGTRFS.
!
                        IF ( info/=0 )                                  &
     &                       CALL ALAERH(path,'CGTRFS',info,0,trans,n,n,&
     &                       -1,-1,nrhs,imat,nfail,nerrs,Nout)
!
                        CALL CGET04(n,nrhs,X,lda,Xact,lda,rcondc,       &
     &                              result(4))
                        CALL CGTT05(trans,n,nrhs,A,A(m+1),A(n+m+1),B,   &
     &                              lda,X,lda,Xact,lda,Rwork,           &
     &                              Rwork(nrhs+1),result(5))
!
!              Print information about the tests that did not pass the
!              threshold.
!
                        DO k = 2 , 6
                           IF ( result(k)>=Thresh ) THEN
                              IF ( nfail==0 .AND. nerrs==0 )            &
     &                             CALL ALAHD(Nout,path)
                              WRITE (Nout,FMT=99002) trans , n , nrhs , &
     &                               imat , k , result(k)
                              nfail = nfail + 1
                           ENDIF
                        ENDDO
                        nrun = nrun + 5
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
99001 FORMAT (12X,'N =',I5,',',10X,' type ',I2,', test(',I2,') = ',     &
     &        G12.5)
99002 FORMAT (' TRANS=''',A1,''', N =',I5,', NRHS=',I3,', type ',I2,    &
     &        ', test(',I2,') = ',G12.5)
99003 FORMAT (' NORM =''',A1,''', N =',I5,',',10X,' type ',I2,', test(',&
     &        I2,') = ',G12.5)
!
!     End of CCHKGT
!
      END SUBROUTINE CCHKGT
