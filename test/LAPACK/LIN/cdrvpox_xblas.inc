!*==cdrvpo.f90  processed by SPAG 7.51RB at 17:34 on  4 Mar 2022
!> \brief \b CDRVPOX
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE cdrvpo( DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, NMAX,
!                          A, AFAC, ASAV, B, BSAV, X, XACT, S, WORK,
!                          RWORK, NOUT )
!
!       .. Scalar Arguments ..
!       LOGICAL            TSTERR
!       INTEGER            NMAX, NN, NOUT, NRHS
!       REAL               THRESH
!       ..
!       .. Array Arguments ..
!       LOGICAL            DOTYPE( * )
!       INTEGER            NVAL( * )
!       REAL               RWORK( * ), S( * )
!       COMPLEX            A( * ), AFAC( * ), ASAV( * ), B( * ),
!      $                   BSAV( * ), WORK( * ), X( * ), XACT( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> cdrvpo tests the driver routines CPOSV, -SVX, and -SVXX.
!>
!> Note that this file is used only when the XBLAS are available,
!> otherwise cdrvpo.f defines this subroutine.
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
!> \param[in] NMAX
!> \verbatim
!>          NMAX is INTEGER
!>          The maximum value permitted for N, used in dimensioning the
!>          work arrays.
!> \endverbatim
!>
!> \param[out] A
!> \verbatim
!>          A is COMPLEX array, dimension (NMAX*NMAX)
!> \endverbatim
!>
!> \param[out] AFAC
!> \verbatim
!>          AFAC is COMPLEX array, dimension (NMAX*NMAX)
!> \endverbatim
!>
!> \param[out] ASAV
!> \verbatim
!>          ASAV is COMPLEX array, dimension (NMAX*NMAX)
!> \endverbatim
!>
!> \param[out] B
!> \verbatim
!>          B is COMPLEX array, dimension (NMAX*NRHS)
!> \endverbatim
!>
!> \param[out] BSAV
!> \verbatim
!>          BSAV is COMPLEX array, dimension (NMAX*NRHS)
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
!> \param[out] S
!> \verbatim
!>          S is REAL array, dimension (NMAX)
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
      SUBROUTINE cdrvpo(Dotype,Nn,Nval,Nrhs,Thresh,Tsterr,Nmax,A,Afac,  &
     &                  Asav,B,Bsav,X,Xact,S,Work,Rwork,Nout)
      IMPLICIT NONE
!*--CDRVPO165
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      LOGICAL Tsterr
      INTEGER Nmax , Nn , Nout , Nrhs
      REAL Thresh
!     ..
!     .. Array Arguments ..
      LOGICAL Dotype(*)
      INTEGER Nval(*)
      REAL Rwork(*) , S(*)
      COMPLEX A(*) , Afac(*) , Asav(*) , B(*) , Bsav(*) , Work(*) ,     &
     &        X(*) , Xact(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ONE , ZERO
      PARAMETER (ONE=1.0E+0,ZERO=0.0E+0)
      INTEGER NTYPES
      PARAMETER (NTYPES=9)
      INTEGER NTESTS
      PARAMETER (NTESTS=6)
!     ..
!     .. Local Scalars ..
      LOGICAL equil , nofact , prefac , zerot
      CHARACTER dist , equed , fact , type , uplo , xtype
      CHARACTER*3 path
      INTEGER i , iequed , ifact , imat , in , info , ioff , iuplo ,    &
     &        izero , k , k1 , kl , ku , lda , mode , n , nb , nbmin ,  &
     &        nerrs , nfact , nfail , nimat , nrun , nt , n_err_bnds
      REAL ainvnm , amax , anorm , cndnum , rcond , rcondc , roldc ,    &
     &     scond , rpvgrw_svxx
!     ..
!     .. Local Arrays ..
      CHARACTER equeds(2) , facts(3) , uplos(2)
      INTEGER iseed(4) , iseedy(4)
      REAL result(NTESTS) , berr(Nrhs) , errbnds_n(Nrhs,3) ,            &
     &     errbnds_c(Nrhs,3)
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      REAL CLANHE  
EXTERNAL    LSAME,   CLANHE
!     ..
!     .. External Subroutines ..
EXTERNAL    CLACPY
EXTERNAL    CLAQHE,   CLASET,   CLATMS
EXTERNAL    CPOEQU,   CPOSV,   CPOSVX
EXTERNAL    CPOTRF,   CPOTRI,   CPOSVXX
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
      INTRINSIC CMPLX , MAX
!     ..
!     .. Data statements ..
      DATA iseedy/1988 , 1989 , 1990 , 1991/
      DATA uplos/'U' , 'L'/
      DATA facts/'F' , 'N' , 'E'/
      DATA equeds/'N' , 'Y'/
!     ..
!     .. Executable Statements ..
!
!     Initialize constants and the random number seed.
!
      path(1:1) = 'Complex precision'
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
      IF ( Tsterr ) CALL cerrvx(path,Nout)
      INFot = 0
!
!     Set the block size and minimum block size for testing.
!
      nb = 1
      nbmin = 2
      CALL xlaenv(1,nb)
      CALL xlaenv(2,nbmin)
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
!              Set up parameters with clatb4 and generate a test matrix
!              with CLATMS.
!
                     CALL clatb4(path,imat,n,n,type,kl,ku,anorm,mode,   &
     &                           cndnum,dist)
!
                     SRNamt = 'CLATMS'
                     CALL CLATMS(n,n,dist,iseed,type,Rwork,mode,cndnum, &
     &                           anorm,kl,ku,uplo,A,lda,Work,info)
!
!              Check error code from CLATMS.
!
                     IF ( info/=0 ) THEN
                        CALL alaerh(path,'CLATMS',info,0,uplo,n,n,-1,-1,&
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
!              Set the imaginary part of the diagonals.
!
                     CALL claipd(n,A,lda+1,0)
!
!              Save a copy of the matrix A in ASAV.
!
                     CALL CLACPY(uplo,n,n,A,lda,Asav,lda)
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
!                       the value returned by CPOSVX (FACT = 'N' reuses
!                       the condition number from the previous iteration
!                       with FACT = 'F').
!
                              CALL CLACPY(uplo,n,n,Asav,lda,Afac,lda)
                              IF ( equil .OR. iequed>1 ) THEN
!
!                          Compute row and column scale factors to
!                          equilibrate the matrix A.
!
                                 CALL CPOEQU(n,Afac,lda,S,scond,amax,   &
     &                              info)
                                 IF ( info==0 .AND. n>0 ) THEN
                                    IF ( iequed>1 ) scond = ZERO
!
!                             Equilibrate the matrix.
!
                                    CALL CLAQHE(uplo,n,Afac,lda,S,scond,&
     &                                 amax,equed)
                                 ENDIF
                              ENDIF
!
!                       Save the condition number of the
!                       non-equilibrated system for use in cget04.
!
                              IF ( equil ) roldc = rcondc
!
!                       Compute the 1-norm of A.
!
                              anorm = CLANHE('1',uplo,n,Afac,lda,Rwork)
!
!                       Factor the matrix A.
!
                              CALL CPOTRF(uplo,n,Afac,lda,info)
!
!                       Form the inverse of A.
!
                              CALL CLACPY(uplo,n,n,Afac,lda,A,lda)
                              CALL CPOTRI(uplo,n,A,lda,info)
!
!                       Compute the 1-norm condition number of A.
!
                              ainvnm = CLANHE('1',uplo,n,A,lda,Rwork)
                              IF ( anorm<=ZERO .OR. ainvnm<=ZERO ) THEN
                                 rcondc = ONE
                              ELSE
                                 rcondc = (ONE/anorm)/ainvnm
                              ENDIF
                           ENDIF
!
!                    Restore the matrix A.
!
                           CALL CLACPY(uplo,n,n,Asav,lda,A,lda)
!
!                    Form an exact solution and set the right hand side.
!
                           SRNamt = 'clarhs'
                           CALL clarhs(path,xtype,uplo,' ',n,n,kl,ku,   &
     &                                 Nrhs,A,lda,Xact,lda,B,lda,iseed, &
     &                                 info)
                           xtype = 'C'
                           CALL CLACPY('Full',n,Nrhs,B,lda,Bsav,lda)
!
                           IF ( nofact ) THEN
!
!                       --- Test CPOSV  ---
!
!                       Compute the L*L' or U'*U factorization of the
!                       matrix and solve the system.
!
                              CALL CLACPY(uplo,n,n,A,lda,Afac,lda)
                              CALL CLACPY('Full',n,Nrhs,B,lda,X,lda)
!
                              SRNamt = 'CPOSV '
                              CALL CPOSV(uplo,n,Nrhs,Afac,lda,X,lda,    &
     &                           info)
!
!                       Check error code from CPOSV .
!
                              IF ( info/=izero ) THEN
                                 CALL alaerh(path,'CPOSV ',info,izero,  &
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
                              CALL cpot01(uplo,n,A,lda,Afac,lda,Rwork,  &
     &                           result(1))
!
!                       Compute residual of the computed solution.
!
                              CALL CLACPY('Full',n,Nrhs,B,lda,Work,lda)
                              CALL cpot02(uplo,n,Nrhs,A,lda,X,lda,Work, &
     &                           lda,Rwork,result(2))
!
!                       Check solution from generated exact solution.
!
                              CALL cget04(n,Nrhs,X,lda,Xact,lda,rcondc, &
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
                                    WRITE (Nout,FMT=99001) 'CPOSV ' ,   &
     &                                 uplo , n , imat , k , result(k)
                                    nfail = nfail + 1
                                 ENDIF
                              ENDDO
                              nrun = nrun + nt
                           ENDIF
!
!                    --- Test CPOSVX ---
!
 2                         IF ( .NOT.prefac )                           &
     &                          CALL CLASET(uplo,n,n,CMPLX(ZERO),       &
     &                          CMPLX(ZERO),Afac,lda)
                           CALL CLASET('Full',n,Nrhs,CMPLX(ZERO),       &
     &                                 CMPLX(ZERO),X,lda)
!
!                       Equilibrate the matrix if FACT='F' and
!                       EQUED='Y'.
!
                           IF ( iequed>1 .AND. n>0 )                    &
     &                          CALL CLAQHE(uplo,n,A,lda,S,scond,amax,  &
     &                          equed)
!
!                    Solve the system and compute the condition number
!                    and error bounds using CPOSVX.
!
                           SRNamt = 'CPOSVX'
                           CALL CPOSVX(fact,uplo,n,Nrhs,A,lda,Afac,lda, &
     &                                 equed,S,B,lda,X,lda,rcond,Rwork, &
     &                                 Rwork(Nrhs+1),Work,              &
     &                                 Rwork(2*Nrhs+1),info)
!
!                    Check the error code from CPOSVX.
!
                           IF ( info/=izero ) THEN
                              CALL alaerh(path,'CPOSVX',info,izero,     &
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
                                 CALL cpot01(uplo,n,A,lda,Afac,lda,     &
     &                              Rwork(2*Nrhs+1),result(1))
                                 k1 = 1
                              ELSE
                                 k1 = 2
                              ENDIF
!
!                       Compute residual of the computed solution.
!
                              CALL CLACPY('Full',n,Nrhs,Bsav,lda,Work,  &
     &                           lda)
                              CALL cpot02(uplo,n,Nrhs,Asav,lda,X,lda,   &
     &                           Work,lda,Rwork(2*Nrhs+1),result(2))
!
!                       Check solution from generated exact solution.
!
                              IF ( nofact .OR.                          &
     &                             (prefac .AND. LSAME(equed,'N')) )    &
     &                             THEN
                                 CALL cget04(n,Nrhs,X,lda,Xact,lda,     &
     &                              rcondc,result(3))
                              ELSE
                                 CALL cget04(n,Nrhs,X,lda,Xact,lda,     &
     &                              roldc,result(3))
                              ENDIF
!
!                       Check the error bounds from iterative
!                       refinement.
!
                              CALL cpot05(uplo,n,Nrhs,Asav,lda,B,lda,X, &
     &                           lda,Xact,lda,Rwork,Rwork(Nrhs+1),      &
     &                           result(4))
                           ELSE
                              k1 = 6
                           ENDIF
!
!                    Compare RCOND from CPOSVX with the computed value
!                    in RCONDC.
!
                           result(6) = sget06(rcond,rcondc)
!
!                    Print information about the tests that did not pass
!                    the threshold.
!
                           DO k = k1 , 6
                              IF ( result(k)>=Thresh ) THEN
                                 IF ( nfail==0 .AND. nerrs==0 )         &
     &                                CALL aladhd(Nout,path)
                                 IF ( prefac ) THEN
                                    WRITE (Nout,FMT=99003) 'CPOSVX' ,   &
     &                                 fact , uplo , n , equed , imat , &
     &                                 k , result(k)
                                 ELSE
                                    WRITE (Nout,FMT=99002) 'CPOSVX' ,   &
     &                                 fact , uplo , n , imat , k ,     &
     &                                 result(k)
                                 ENDIF
                                 nfail = nfail + 1
                              ENDIF
                           ENDDO
                           nrun = nrun + 7 - k1
!
!                    --- Test CPOSVXX ---
!
!                    Restore the matrices A and B.
!
                           CALL CLACPY('Full',n,n,Asav,lda,A,lda)
                           CALL CLACPY('Full',n,Nrhs,Bsav,lda,B,lda)
 
                           IF ( .NOT.prefac )                           &
     &                          CALL CLASET(uplo,n,n,CMPLX(ZERO),       &
     &                          CMPLX(ZERO),Afac,lda)
                           CALL CLASET('Full',n,Nrhs,CMPLX(ZERO),       &
     &                                 CMPLX(ZERO),X,lda)
!
!                       Equilibrate the matrix if FACT='F' and
!                       EQUED='Y'.
!
                           IF ( iequed>1 .AND. n>0 )                    &
     &                          CALL CLAQHE(uplo,n,A,lda,S,scond,amax,  &
     &                          equed)
!
!                    Solve the system and compute the condition number
!                    and error bounds using CPOSVXX.
!
                           SRNamt = 'CPOSVXX'
                           n_err_bnds = 3
                           CALL CPOSVXX(fact,uplo,n,Nrhs,A,lda,Afac,lda,&
     &                                  equed,S,B,lda,X,lda,rcond,      &
     &                                  rpvgrw_svxx,berr,n_err_bnds,    &
     &                                  errbnds_n,errbnds_c,0,ZERO,Work,&
     &                                  Rwork(2*Nrhs+1),info)
!
!                    Check the error code from CPOSVXX.
!
                           IF ( info/=n+1 ) THEN
                              IF ( info/=izero ) THEN
                                 CALL alaerh(path,'CPOSVXX',info,izero, &
     &                              fact//uplo,n,n,-1,-1,Nrhs,imat,     &
     &                              nfail,nerrs,Nout)
                                 CYCLE
                              ENDIF
!
                              IF ( info==0 ) THEN
                                 IF ( .NOT.prefac ) THEN
!
!                          Reconstruct matrix from factors and compute
!                          residual.
!
                                    CALL cpot01(uplo,n,A,lda,Afac,lda,  &
     &                                 Rwork(2*Nrhs+1),result(1))
                                    k1 = 1
                                 ELSE
                                    k1 = 2
                                 ENDIF
!
!                       Compute residual of the computed solution.
!
                                 CALL CLACPY('Full',n,Nrhs,Bsav,lda,    &
     &                              Work,lda)
                                 CALL cpot02(uplo,n,Nrhs,Asav,lda,X,lda,&
     &                              Work,lda,Rwork(2*Nrhs+1),result(2))
!
!                       Check solution from generated exact solution.
!
                                 IF ( nofact .OR.                       &
     &                                (prefac .AND. LSAME(equed,'N')) ) &
     &                                THEN
                                    CALL cget04(n,Nrhs,X,lda,Xact,lda,  &
     &                                 rcondc,result(3))
                                 ELSE
                                    CALL cget04(n,Nrhs,X,lda,Xact,lda,  &
     &                                 roldc,result(3))
                                 ENDIF
!
!                       Check the error bounds from iterative
!                       refinement.
!
                                 CALL cpot05(uplo,n,Nrhs,Asav,lda,B,lda,&
     &                              X,lda,Xact,lda,Rwork,Rwork(Nrhs+1), &
     &                              result(4))
                              ELSE
                                 k1 = 6
                              ENDIF
!
!                    Compare RCOND from CPOSVXX with the computed value
!                    in RCONDC.
!
                              result(6) = sget06(rcond,rcondc)
!
!                    Print information about the tests that did not pass
!                    the threshold.
!
                              DO k = k1 , 6
                                 IF ( result(k)>=Thresh ) THEN
                                    IF ( nfail==0 .AND. nerrs==0 )      &
     &                                 CALL aladhd(Nout,path)
                                    IF ( prefac ) THEN
                                       WRITE (Nout,FMT=99003)           &
     &                                    'CPOSVXX' , fact , uplo , n , &
     &                                    equed , imat , k , result(k)
                                    ELSE
                                       WRITE (Nout,FMT=99002)           &
     &                                    'CPOSVXX' , fact , uplo , n , &
     &                                    imat , k , result(k)
                                    ENDIF
                                    nfail = nfail + 1
                                 ENDIF
                              ENDDO
                              nrun = nrun + 7 - k1
                           ENDIF
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
 
!     Test Error Bounds for CGESVXX
 
      CALL CEBCHVXX(Thresh,path)
 
99001 FORMAT (1X,A,', UPLO=''',A1,''', N =',I5,', type ',I1,', test(',  &
     &        I1,')=',G12.5)
99002 FORMAT (1X,A,', FACT=''',A1,''', UPLO=''',A1,''', N=',I5,         &
     &        ', type ',I1,', test(',I1,')=',G12.5)
99003 FORMAT (1X,A,', FACT=''',A1,''', UPLO=''',A1,''', N=',I5,         &
     &        ', EQUED=''',A1,''', type ',I1,', test(',I1,') =',G12.5)
!
!     End of cdrvpo
!
      END SUBROUTINE cdrvpo
