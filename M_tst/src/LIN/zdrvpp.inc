!*==zdrvpp.f90  processed by SPAG 7.51RB at 17:35 on  4 Mar 2022
!> \brief \b zdrvpp
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE zdrvpp( DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, NMAX,
!                          A, AFAC, ASAV, B, BSAV, X, XACT, S, WORK,
!                          RWORK, NOUT )
!
!       .. Scalar Arguments ..
!       LOGICAL            TSTERR
!       INTEGER            NMAX, NN, NOUT, NRHS
!       DOUBLE PRECISION   THRESH
!       ..
!       .. Array Arguments ..
!       LOGICAL            DOTYPE( * )
!       INTEGER            NVAL( * )
!       DOUBLE PRECISION   RWORK( * ), S( * )
!       COMPLEX*16         A( * ), AFAC( * ), ASAV( * ), B( * ),
!      $                   BSAV( * ), WORK( * ), X( * ), XACT( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> zdrvpp tests the driver routines ZPPSV and -SVX.
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
!>          A is COMPLEX*16 array, dimension (NMAX*(NMAX+1)/2)
!> \endverbatim
!>
!> \param[out] AFAC
!> \verbatim
!>          AFAC is COMPLEX*16 array, dimension (NMAX*(NMAX+1)/2)
!> \endverbatim
!>
!> \param[out] ASAV
!> \verbatim
!>          ASAV is COMPLEX*16 array, dimension (NMAX*(NMAX+1)/2)
!> \endverbatim
!>
!> \param[out] B
!> \verbatim
!>          B is COMPLEX*16 array, dimension (NMAX*NRHS)
!> \endverbatim
!>
!> \param[out] BSAV
!> \verbatim
!>          BSAV is COMPLEX*16 array, dimension (NMAX*NRHS)
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
!> \param[out] S
!> \verbatim
!>          S is DOUBLE PRECISION array, dimension (NMAX)
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX*16 array, dimension
!>                      (NMAX*max(3,NRHS))
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is DOUBLE PRECISION array, dimension (NMAX+2*NRHS)
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
!> \ingroup complex16_lin
!
!  =====================================================================
      SUBROUTINE zdrvpp(Dotype,Nn,Nval,Nrhs,Thresh,Tsterr,Nmax,A,Afac,  &
     &                  Asav,B,Bsav,X,Xact,S,Work,Rwork,Nout)
      IMPLICIT NONE
!*--ZDRVPP162
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
      INTEGER Nval(*)
      DOUBLE PRECISION Rwork(*) , S(*)
      COMPLEX*16 A(*) , Afac(*) , Asav(*) , B(*) , Bsav(*) , Work(*) ,  &
     &           X(*) , Xact(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ONE , ZERO
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
      INTEGER NTYPES
      PARAMETER (NTYPES=9)
      INTEGER NTESTS
      PARAMETER (NTESTS=6)
!     ..
!     .. Local Scalars ..
      LOGICAL equil , nofact , prefac , zerot
      CHARACTER dist , equed , fact , packit , type , uplo , xtype
      CHARACTER*3 path
      INTEGER i , iequed , ifact , imat , in , info , ioff , iuplo ,    &
     &        izero , k , k1 , kl , ku , lda , mode , n , nerrs ,       &
     &        nfact , nfail , nimat , npp , nrun , nt
      DOUBLE PRECISION ainvnm , amax , anorm , cndnum , rcond , rcondc ,&
     &                 roldc , scond
!     ..
!     .. Local Arrays ..
      CHARACTER equeds(2) , facts(3) , packs(2) , uplos(2)
      INTEGER iseed(4) , iseedy(4)
      DOUBLE PRECISION result(NTESTS)
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      DOUBLE PRECISION ZLANHP
EXTERNAL    LSAME,   ZLANHP
!     ..
!     .. External Subroutines ..
EXTERNAL    ZCOPY
EXTERNAL    ZLACPY,   ZLAQHP,   ZLASET
EXTERNAL    ZPPEQU,   ZPPSV,   ZPPSVX
EXTERNAL    ZPPTRF,   ZPPTRI
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
      INTRINSIC DCMPLX , MAX
!     ..
!     .. Data statements ..
      DATA iseedy/1988 , 1989 , 1990 , 1991/
      DATA uplos/'U' , 'L'/ , facts/'F' , 'N' , 'E'/ , packs/'C' ,      &
     &     'R'/ , equeds/'N' , 'Y'/
!     ..
!     .. Executable Statements ..
!
!     Initialize constants and the random number seed.
!
      path(1:1) = 'Zomplex precision'
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
      IF ( Tsterr ) CALL zerrvx(path,Nout)
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
!              Set up parameters with zlatb4 and generate a test matrix
!              with ZLATMS.
!
                     CALL zlatb4(path,imat,n,n,type,kl,ku,anorm,mode,   &
     &                           cndnum,dist)
                     rcondc = ONE/cndnum
!
                     SRNamt = 'ZLATMS'
                     CALL ZLATMS(n,n,dist,iseed,type,Rwork,mode,cndnum, &
     &                           anorm,kl,ku,packit,A,lda,Work,info)
!
!              Check error code from ZLATMS.
!
                     IF ( info/=0 ) THEN
                        CALL alaerh(path,'ZLATMS',info,0,uplo,n,n,-1,-1,&
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
                        CALL zlaipd(n,A,2,1)
                     ELSE
                        CALL zlaipd(n,A,n,-1)
                     ENDIF
!
!              Save a copy of the matrix A in ASAV.
!
                     CALL ZCOPY(npp,A,1,Asav,1)
!
                     DO iequed = 1 , 2
                        equed = equeds(iequed)
                        IF ( iequed==1 ) THEN
                           nfact = 3
                        ELSE
                           nfact = 1
                        ENDIF
!
                        DO ifact = 1 , nfact
                           fact = facts(ifact)
                           prefac = LSAME(fact,'F')
                           nofact = LSAME(fact,'N')
                           equil = LSAME(fact,'E')
!
                           IF ( zerot ) THEN
                              IF ( prefac ) CYCLE
                              rcondc = ZERO
!
                           ELSEIF ( .NOT.LSAME(fact,'N') ) THEN
!
!                       Compute the condition number for comparison with
!                       the value returned by ZPPSVX (FACT = 'N' reuses
!                       the condition number from the previous iteration
!                          with FACT = 'F').
!
                              CALL ZCOPY(npp,Asav,1,Afac,1)
                              IF ( equil .OR. iequed>1 ) THEN
!
!                          Compute row and column scale factors to
!                          equilibrate the matrix A.
!
                                 CALL ZPPEQU(uplo,n,Afac,S,scond,amax,  &
     &                              info)
                                 IF ( info==0 .AND. n>0 ) THEN
                                    IF ( iequed>1 ) scond = ZERO
!
!                             Equilibrate the matrix.
!
                                    CALL ZLAQHP(uplo,n,Afac,S,scond,    &
     &                                 amax,equed)
                                 ENDIF
                              ENDIF
!
!                       Save the condition number of the
!                       non-equilibrated system for use in zget04.
!
                              IF ( equil ) roldc = rcondc
!
!                       Compute the 1-norm of A.
!
                              anorm = ZLANHP('1',uplo,n,Afac,Rwork)
!
!                       Factor the matrix A.
!
                              CALL ZPPTRF(uplo,n,Afac,info)
!
!                       Form the inverse of A.
!
                              CALL ZCOPY(npp,Afac,1,A,1)
                              CALL ZPPTRI(uplo,n,A,info)
!
!                       Compute the 1-norm condition number of A.
!
                              ainvnm = ZLANHP('1',uplo,n,A,Rwork)
                              IF ( anorm<=ZERO .OR. ainvnm<=ZERO ) THEN
                                 rcondc = ONE
                              ELSE
                                 rcondc = (ONE/anorm)/ainvnm
                              ENDIF
                           ENDIF
!
!                    Restore the matrix A.
!
                           CALL ZCOPY(npp,Asav,1,A,1)
!
!                    Form an exact solution and set the right hand side.
!
                           SRNamt = 'zlarhs'
                           CALL zlarhs(path,xtype,uplo,' ',n,n,kl,ku,   &
     &                                 Nrhs,A,lda,Xact,lda,B,lda,iseed, &
     &                                 info)
                           xtype = 'C'
                           CALL ZLACPY('Full',n,Nrhs,B,lda,Bsav,lda)
!
                           IF ( nofact ) THEN
!
!                       --- Test ZPPSV  ---
!
!                       Compute the L*L' or U'*U factorization of the
!                       matrix and solve the system.
!
                              CALL ZCOPY(npp,A,1,Afac,1)
                              CALL ZLACPY('Full',n,Nrhs,B,lda,X,lda)
!
                              SRNamt = 'ZPPSV '
                              CALL ZPPSV(uplo,n,Nrhs,Afac,X,lda,info)
!
!                       Check error code from ZPPSV .
!
                              IF ( info/=izero ) THEN
                                 CALL alaerh(path,'ZPPSV ',info,izero,  &
     &                              uplo,n,n,-1,-1,Nrhs,imat,nfail,     &
     &                              nerrs,Nout)
                                 GOTO 2
                              ELSEIF ( info/=0 ) THEN
                                 GOTO 2
                              ENDIF
!
!                       Reconstruct matrix from factors and compute
!                       residual.
!
                              CALL zppt01(uplo,n,A,Afac,Rwork,result(1))
!
!                       Compute residual of the computed solution.
!
                              CALL ZLACPY('Full',n,Nrhs,B,lda,Work,lda)
                              CALL zppt02(uplo,n,Nrhs,A,X,lda,Work,lda, &
     &                           Rwork,result(2))
!
!                       Check solution from generated exact solution.
!
                              CALL zget04(n,Nrhs,X,lda,Xact,lda,rcondc, &
     &                           result(3))
                              nt = 3
!
!                       Print information about the tests that did not
!                       pass the threshold.
!
                              DO k = 1 , nt
                                 IF ( result(k)>=Thresh ) THEN
                                    IF ( nfail==0 .AND. nerrs==0 )      &
     &                                 CALL aladhd(Nout,path)
                                    WRITE (Nout,FMT=99001) 'ZPPSV ' ,   &
     &                                 uplo , n , imat , k , result(k)
                                    nfail = nfail + 1
                                 ENDIF
                              ENDDO
                              nrun = nrun + nt
                           ENDIF
!
!                    --- Test ZPPSVX ---
!
 2                         IF ( .NOT.prefac .AND. npp>0 )               &
     &                           CALL ZLASET('Full',npp,1,DCMPLX(ZERO), &
     &                          DCMPLX(ZERO),Afac,npp)
                           CALL ZLASET('Full',n,Nrhs,DCMPLX(ZERO),      &
     &                                 DCMPLX(ZERO),X,lda)
!
!                       Equilibrate the matrix if FACT='F' and
!                       EQUED='Y'.
!
                           IF ( iequed>1 .AND. n>0 )                    &
     &                          CALL ZLAQHP(uplo,n,A,S,scond,amax,equed)
!
!                    Solve the system and compute the condition number
!                    and error bounds using ZPPSVX.
!
                           SRNamt = 'ZPPSVX'
                           CALL ZPPSVX(fact,uplo,n,Nrhs,A,Afac,equed,S, &
     &                                 B,lda,X,lda,rcond,Rwork,         &
     &                                 Rwork(Nrhs+1),Work,              &
     &                                 Rwork(2*Nrhs+1),info)
!
!                    Check the error code from ZPPSVX.
!
                           IF ( info/=izero ) THEN
                              CALL alaerh(path,'ZPPSVX',info,izero,     &
     &                           fact//uplo,n,n,-1,-1,Nrhs,imat,nfail,  &
     &                           nerrs,Nout)
                              CYCLE
                           ENDIF
!
                           IF ( info==0 ) THEN
                              IF ( .NOT.prefac ) THEN
!
!                          Reconstruct matrix from factors and compute
!                          residual.
!
                                 CALL zppt01(uplo,n,A,Afac,             &
     &                              Rwork(2*Nrhs+1),result(1))
                                 k1 = 1
                              ELSE
                                 k1 = 2
                              ENDIF
!
!                       Compute residual of the computed solution.
!
                              CALL ZLACPY('Full',n,Nrhs,Bsav,lda,Work,  &
     &                           lda)
                              CALL zppt02(uplo,n,Nrhs,Asav,X,lda,Work,  &
     &                           lda,Rwork(2*Nrhs+1),result(2))
!
!                       Check solution from generated exact solution.
!
                              IF ( nofact .OR.                          &
     &                             (prefac .AND. LSAME(equed,'N')) )    &
     &                             THEN
                                 CALL zget04(n,Nrhs,X,lda,Xact,lda,     &
     &                              rcondc,result(3))
                              ELSE
                                 CALL zget04(n,Nrhs,X,lda,Xact,lda,     &
     &                              roldc,result(3))
                              ENDIF
!
!                       Check the error bounds from iterative
!                       refinement.
!
                              CALL zppt05(uplo,n,Nrhs,Asav,B,lda,X,lda, &
     &                           Xact,lda,Rwork,Rwork(Nrhs+1),result(4))
                           ELSE
                              k1 = 6
                           ENDIF
!
!                    Compare RCOND from ZPPSVX with the computed value
!                    in RCONDC.
!
                           result(6) = dget06(rcond,rcondc)
!
!                    Print information about the tests that did not pass
!                    the threshold.
!
                           DO k = k1 , 6
                              IF ( result(k)>=Thresh ) THEN
                                 IF ( nfail==0 .AND. nerrs==0 )         &
     &                                CALL aladhd(Nout,path)
                                 IF ( prefac ) THEN
                                    WRITE (Nout,FMT=99003) 'ZPPSVX' ,   &
     &                                 fact , uplo , n , equed , imat , &
     &                                 k , result(k)
                                 ELSE
                                    WRITE (Nout,FMT=99002) 'ZPPSVX' ,   &
     &                                 fact , uplo , n , imat , k ,     &
     &                                 result(k)
                                 ENDIF
                                 nfail = nfail + 1
                              ENDIF
                           ENDDO
                           nrun = nrun + 7 - k1
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
      CALL alasvm(path,Nout,nfail,nrun,nerrs)
!
99001 FORMAT (1X,A,', UPLO=''',A1,''', N =',I5,', type ',I1,', test(',  &
     &        I1,')=',G12.5)
99002 FORMAT (1X,A,', FACT=''',A1,''', UPLO=''',A1,''', N=',I5,         &
     &        ', type ',I1,', test(',I1,')=',G12.5)
99003 FORMAT (1X,A,', FACT=''',A1,''', UPLO=''',A1,''', N=',I5,         &
     &        ', EQUED=''',A1,''', type ',I1,', test(',I1,')=',G12.5)
!
!     End of zdrvpp
!
      END SUBROUTINE zdrvpp
