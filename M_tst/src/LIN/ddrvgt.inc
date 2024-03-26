!*==ddrvgt.f90  processed by SPAG 7.51RB at 17:34 on  4 Mar 2022
!> \brief \b ddrvgt
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ddrvgt( DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, A, AF,
!                          B, X, XACT, WORK, RWORK, IWORK, NOUT )
!
!       .. Scalar Arguments ..
!       LOGICAL            TSTERR
!       INTEGER            NN, NOUT, NRHS
!       DOUBLE PRECISION   THRESH
!       ..
!       .. Array Arguments ..
!       LOGICAL            DOTYPE( * )
!       INTEGER            IWORK( * ), NVAL( * )
!       DOUBLE PRECISION   A( * ), AF( * ), B( * ), RWORK( * ), WORK( * ),
!      $                   X( * ), XACT( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ddrvgt tests DGTSV and -SVX.
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
!>          The number of right hand sides, NRHS >= 0.
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
!>          A is DOUBLE PRECISION array, dimension (NMAX*4)
!> \endverbatim
!>
!> \param[out] AF
!> \verbatim
!>          AF is DOUBLE PRECISION array, dimension (NMAX*4)
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
!>                      (NMAX*max(3,NRHS))
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is DOUBLE PRECISION array, dimension
!>                      (max(NMAX,2*NRHS))
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
      SUBROUTINE ddrvgt(Dotype,Nn,Nval,Nrhs,Thresh,Tsterr,A,Af,B,X,Xact,&
     &                  Work,Rwork,Iwork,Nout)
      IMPLICIT NONE
!*--DDRVGT143
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      LOGICAL Tsterr
      INTEGER Nn , Nout , Nrhs
      DOUBLE PRECISION Thresh
!     ..
!     .. Array Arguments ..
      LOGICAL Dotype(*)
      INTEGER Iwork(*) , Nval(*)
      DOUBLE PRECISION A(*) , Af(*) , B(*) , Rwork(*) , Work(*) , X(*) ,&
     &                 Xact(*)
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
      PARAMETER (NTESTS=6)
!     ..
!     .. Local Scalars ..
      LOGICAL trfcon , zerot
      CHARACTER dist , fact , trans , type
      CHARACTER*3 path
      INTEGER i , ifact , imat , in , info , itran , ix , izero , j ,   &
     &        k , k1 , kl , koff , ku , lda , m , mode , n , nerrs ,    &
     &        nfail , nimat , nrun , nt
      DOUBLE PRECISION ainvnm , anorm , anormi , anormo , cond , rcond ,&
     &                 rcondc , rcondi , rcondo
!     ..
!     .. Local Arrays ..
      CHARACTER transs(3)
      INTEGER iseed(4) , iseedy(4)
      DOUBLE PRECISION result(NTESTS) , z(3)
!     ..
!     .. External Functions ..
      DOUBLE PRECISION DASUM , DLANGT
EXTERNAL    DASUM,   DLANGT
!     ..
!     .. External Subroutines ..
EXTERNAL    DCOPY
EXTERNAL    DGTSV,   DGTSVX,   DGTTRF
EXTERNAL    DGTTRS,   DLACPY,   DLAGTM,   DLARNV,   DLASET
EXTERNAL    DSCAL
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
      path(1:1) = 'Double precision'
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
      IF ( Tsterr ) CALL derrvx(path,Nout)
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
!           Set up parameters with dlatb4.
!
               CALL dlatb4(path,imat,n,n,type,kl,ku,anorm,mode,cond,    &
     &                     dist)
!
               zerot = imat>=8 .AND. imat<=10
               IF ( imat<=6 ) THEN
!
!              Types 1-6:  generate matrices of known condition number.
!
                  koff = MAX(2-ku,3-MAX(1,n))
                  SRNamt = 'DLATMS'
                  CALL DLATMS(n,n,dist,iseed,type,Rwork,mode,cond,anorm,&
     &                        kl,ku,'Z',Af(koff),3,Work,info)
!
!              Check the error code from DLATMS.
!
                  IF ( info/=0 ) THEN
                     CALL alaerh(path,'DLATMS',info,0,' ',n,n,kl,ku,-1, &
     &                           imat,nfail,nerrs,Nout)
                     CYCLE
                  ENDIF
                  izero = 0
!
                  IF ( n>1 ) THEN
                     CALL DCOPY(n-1,Af(4),3,A,1)
                     CALL DCOPY(n-1,Af(3),3,A(n+m+1),1)
                  ENDIF
                  CALL DCOPY(n,Af(2),3,A(m+1),1)
               ELSE
!
!              Types 7-12:  generate tridiagonal matrices with
!              unknown condition numbers.
!
                  IF ( .NOT.zerot .OR. .NOT.Dotype(7) ) THEN
!
!                 Generate a matrix with elements from [-1,1].
!
                     CALL DLARNV(2,iseed,n+2*m,A)
                     IF ( anorm/=ONE ) CALL DSCAL(n+2*m,anorm,A,1)
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
               DO ifact = 1 , 2
                  IF ( ifact==1 ) THEN
                     fact = 'F'
                  ELSE
                     fact = 'N'
                  ENDIF
!
!              Compute the condition number for comparison with
!              the value returned by DGTSVX.
!
                  IF ( zerot ) THEN
                     IF ( ifact==1 ) CYCLE
                     rcondo = ZERO
                     rcondi = ZERO
!
                  ELSEIF ( ifact==1 ) THEN
                     CALL DCOPY(n+2*m,A,1,Af,1)
!
!                 Compute the 1-norm and infinity-norm of A.
!
                     anormo = DLANGT('1',n,A,A(m+1),A(n+m+1))
                     anormi = DLANGT('I',n,A,A(m+1),A(n+m+1))
!
!                 Factor the matrix A.
!
                     CALL DGTTRF(n,Af,Af(m+1),Af(n+m+1),Af(n+2*m+1),    &
     &                           Iwork,info)
!
!                 Use DGTTRS to solve for one column at a time of
!                 inv(A), computing the maximum column sum as we go.
!
                     ainvnm = ZERO
                     DO i = 1 , n
                        DO j = 1 , n
                           X(j) = ZERO
                        ENDDO
                        X(i) = ONE
                        CALL DGTTRS('No transpose',n,1,Af,Af(m+1),      &
     &                              Af(n+m+1),Af(n+2*m+1),Iwork,X,lda,  &
     &                              info)
                        ainvnm = MAX(ainvnm,DASUM(n,X,1))
                     ENDDO
!
!                 Compute the 1-norm condition number of A.
!
                     IF ( anormo<=ZERO .OR. ainvnm<=ZERO ) THEN
                        rcondo = ONE
                     ELSE
                        rcondo = (ONE/anormo)/ainvnm
                     ENDIF
!
!                 Use DGTTRS to solve for one column at a time of
!                 inv(A'), computing the maximum column sum as we go.
!
                     ainvnm = ZERO
                     DO i = 1 , n
                        DO j = 1 , n
                           X(j) = ZERO
                        ENDDO
                        X(i) = ONE
                        CALL DGTTRS('Transpose',n,1,Af,Af(m+1),Af(n+m+1)&
     &                              ,Af(n+2*m+1),Iwork,X,lda,info)
                        ainvnm = MAX(ainvnm,DASUM(n,X,1))
                     ENDDO
!
!                 Compute the infinity-norm condition number of A.
!
                     IF ( anormi<=ZERO .OR. ainvnm<=ZERO ) THEN
                        rcondi = ONE
                     ELSE
                        rcondi = (ONE/anormi)/ainvnm
                     ENDIF
                  ENDIF
!
                  DO itran = 1 , 3
                     trans = transs(itran)
                     IF ( itran==1 ) THEN
                        rcondc = rcondo
                     ELSE
                        rcondc = rcondi
                     ENDIF
!
!                 Generate NRHS random solution vectors.
!
                     ix = 1
                     DO j = 1 , Nrhs
                        CALL DLARNV(2,iseed,n,Xact(ix))
                        ix = ix + lda
                     ENDDO
!
!                 Set the right hand side.
!
                     CALL DLAGTM(trans,n,Nrhs,ONE,A,A(m+1),A(n+m+1),    &
     &                           Xact,lda,ZERO,B,lda)
!
                     IF ( ifact==2 .AND. itran==1 ) THEN
!
!                    --- Test DGTSV  ---
!
!                    Solve the system using Gaussian elimination with
!                    partial pivoting.
!
                        CALL DCOPY(n+2*m,A,1,Af,1)
                        CALL DLACPY('Full',n,Nrhs,B,lda,X,lda)
!
                        SRNamt = 'DGTSV '
                        CALL DGTSV(n,Nrhs,Af,Af(m+1),Af(n+m+1),X,lda,   &
     &                             info)
!
!                    Check error code from DGTSV .
!
                        IF ( info/=izero )                              &
     &                       CALL alaerh(path,'DGTSV ',info,izero,' ',n,&
     &                       n,1,1,Nrhs,imat,nfail,nerrs,Nout)
                        nt = 1
                        IF ( izero==0 ) THEN
!
!                       Check residual of computed solution.
!
                           CALL DLACPY('Full',n,Nrhs,B,lda,Work,lda)
                           CALL dgtt02(trans,n,Nrhs,A,A(m+1),A(n+m+1),X,&
     &                                 lda,Work,lda,result(2))
!
!                       Check solution from generated exact solution.
!
                           CALL dget04(n,Nrhs,X,lda,Xact,lda,rcondc,    &
     &                                 result(3))
                           nt = 3
                        ENDIF
!
!                    Print information about the tests that did not pass
!                    the threshold.
!
                        DO k = 2 , nt
                           IF ( result(k)>=Thresh ) THEN
                              IF ( nfail==0 .AND. nerrs==0 )            &
     &                             CALL aladhd(Nout,path)
                              WRITE (Nout,FMT=99001) 'DGTSV ' , n ,     &
     &                               imat , k , result(k)
                              nfail = nfail + 1
                           ENDIF
                        ENDDO
                        nrun = nrun + nt - 1
                     ENDIF
!
!                 --- Test DGTSVX ---
!
                     IF ( ifact>1 ) THEN
!
!                    Initialize AF to zero.
!
                        DO i = 1 , 3*n - 2
                           Af(i) = ZERO
                        ENDDO
                     ENDIF
                     CALL DLASET('Full',n,Nrhs,ZERO,ZERO,X,lda)
!
!                 Solve the system and compute the condition number and
!                 error bounds using DGTSVX.
!
                     SRNamt = 'DGTSVX'
                     CALL DGTSVX(fact,trans,n,Nrhs,A,A(m+1),A(n+m+1),Af,&
     &                           Af(m+1),Af(n+m+1),Af(n+2*m+1),Iwork,B, &
     &                           lda,X,lda,rcond,Rwork,Rwork(Nrhs+1),   &
     &                           Work,Iwork(n+1),info)
!
!                 Check the error code from DGTSVX.
!
                     IF ( info/=izero )                                 &
     &                    CALL alaerh(path,'DGTSVX',info,izero,         &
     &                    fact//trans,n,n,1,1,Nrhs,imat,nfail,nerrs,    &
     &                    Nout)
!
                     IF ( ifact>=2 ) THEN
!
!                    Reconstruct matrix from factors and compute
!                    residual.
!
                        CALL dgtt01(n,A,A(m+1),A(n+m+1),Af,Af(m+1),     &
     &                              Af(n+m+1),Af(n+2*m+1),Iwork,Work,   &
     &                              lda,Rwork,result(1))
                        k1 = 1
                     ELSE
                        k1 = 2
                     ENDIF
!
                     IF ( info==0 ) THEN
                        trfcon = .FALSE.
!
!                    Check residual of computed solution.
!
                        CALL DLACPY('Full',n,Nrhs,B,lda,Work,lda)
                        CALL dgtt02(trans,n,Nrhs,A,A(m+1),A(n+m+1),X,   &
     &                              lda,Work,lda,result(2))
!
!                    Check solution from generated exact solution.
!
                        CALL dget04(n,Nrhs,X,lda,Xact,lda,rcondc,       &
     &                              result(3))
!
!                    Check the error bounds from iterative refinement.
!
                        CALL dgtt05(trans,n,Nrhs,A,A(m+1),A(n+m+1),B,   &
     &                              lda,X,lda,Xact,lda,Rwork,           &
     &                              Rwork(Nrhs+1),result(4))
                        nt = 5
                     ENDIF
!
!                 Print information about the tests that did not pass
!                 the threshold.
!
                     DO k = k1 , nt
                        IF ( result(k)>=Thresh ) THEN
                           IF ( nfail==0 .AND. nerrs==0 )               &
     &                          CALL aladhd(Nout,path)
                           WRITE (Nout,FMT=99002) 'DGTSVX' , fact ,     &
     &                            trans , n , imat , k , result(k)
                           nfail = nfail + 1
                        ENDIF
                     ENDDO
!
!                 Check the reciprocal of the condition number.
!
                     result(6) = dget06(rcond,rcondc)
                     IF ( result(6)>=Thresh ) THEN
                        IF ( nfail==0 .AND. nerrs==0 )                  &
     &                       CALL aladhd(Nout,path)
                        WRITE (Nout,FMT=99002) 'DGTSVX' , fact , trans ,&
     &                         n , imat , k , result(k)
                        nfail = nfail + 1
                     ENDIF
                     nrun = nrun + nt - k1 + 2
!
                  ENDDO
               ENDDO
            ENDIF
         ENDDO
      ENDDO
!
!     Print a summary of the results.
!
      CALL alasvm(path,Nout,nfail,nrun,nerrs)
!
99001 FORMAT (1X,A,', N =',I5,', type ',I2,', test ',I2,', ratio = ',   &
     &        G12.5)
99002 FORMAT (1X,A,', FACT=''',A1,''', TRANS=''',A1,''', N =',I5,       &
     &        ', type ',I2,', test ',I2,', ratio = ',G12.5)
!
!     End of ddrvgt
!
      END SUBROUTINE ddrvgt
