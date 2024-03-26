!*==sdrvpb.f90  processed by SPAG 7.51RB at 17:34 on  4 Mar 2022
!> \brief \b sdrvpb
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE sdrvpb( DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, NMAX,
!                          A, AFAC, ASAV, B, BSAV, X, XACT, S, WORK,
!                          RWORK, IWORK, NOUT )
!
!       .. Scalar Arguments ..
!       LOGICAL            TSTERR
!       INTEGER            NMAX, NN, NOUT, NRHS
!       REAL               THRESH
!       ..
!       .. Array Arguments ..
!       LOGICAL            DOTYPE( * )
!       INTEGER            IWORK( * ), NVAL( * )
!       REAL               A( * ), AFAC( * ), ASAV( * ), B( * ),
!      $                   BSAV( * ), RWORK( * ), S( * ), WORK( * ),
!      $                   X( * ), XACT( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> sdrvpb tests the driver routines SPBSV and -SVX.
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
!>          A is REAL array, dimension (NMAX*NMAX)
!> \endverbatim
!>
!> \param[out] AFAC
!> \verbatim
!>          AFAC is REAL array, dimension (NMAX*NMAX)
!> \endverbatim
!>
!> \param[out] ASAV
!> \verbatim
!>          ASAV is REAL array, dimension (NMAX*NMAX)
!> \endverbatim
!>
!> \param[out] B
!> \verbatim
!>          B is REAL array, dimension (NMAX*NRHS)
!> \endverbatim
!>
!> \param[out] BSAV
!> \verbatim
!>          BSAV is REAL array, dimension (NMAX*NRHS)
!> \endverbatim
!>
!> \param[out] X
!> \verbatim
!>          X is REAL array, dimension (NMAX*NRHS)
!> \endverbatim
!>
!> \param[out] XACT
!> \verbatim
!>          XACT is REAL array, dimension (NMAX*NRHS)
!> \endverbatim
!>
!> \param[out] S
!> \verbatim
!>          S is REAL array, dimension (NMAX)
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is REAL array, dimension
!>                      (NMAX*max(3,NRHS))
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is REAL array, dimension (NMAX+2*NRHS)
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
!> \ingroup single_lin
!
!  =====================================================================
      SUBROUTINE sdrvpb(Dotype,Nn,Nval,Nrhs,Thresh,Tsterr,Nmax,A,Afac,  &
     &                  Asav,B,Bsav,X,Xact,S,Work,Rwork,Iwork,Nout)
      IMPLICIT NONE
!*--SDRVPB167
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
      INTEGER Iwork(*) , Nval(*)
      REAL A(*) , Afac(*) , Asav(*) , B(*) , Bsav(*) , Rwork(*) , S(*) ,&
     &     Work(*) , X(*) , Xact(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ONE , ZERO
      PARAMETER (ONE=1.0E+0,ZERO=0.0E+0)
      INTEGER NTYPES , NTESTS
      PARAMETER (NTYPES=8,NTESTS=6)
      INTEGER NBW
      PARAMETER (NBW=4)
!     ..
!     .. Local Scalars ..
      LOGICAL equil , nofact , prefac , zerot
      CHARACTER dist , equed , fact , packit , type , uplo , xtype
      CHARACTER*3 path
      INTEGER i , i1 , i2 , iequed , ifact , ikd , imat , in , info ,   &
     &        ioff , iuplo , iw , izero , k , k1 , kd , kl , koff , ku ,&
     &        lda , ldab , mode , n , nb , nbmin , nerrs , nfact ,      &
     &        nfail , nimat , nkd , nrun , nt
      REAL ainvnm , amax , anorm , cndnum , rcond , rcondc , roldc ,    &
     &     scond
!     ..
!     .. Local Arrays ..
      CHARACTER equeds(2) , facts(3)
      INTEGER iseed(4) , iseedy(4) , kdval(NBW)
      REAL result(NTESTS)
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      REAL SLANGE , SLANSB
EXTERNAL    LSAME,   SLANGE,   SLANSB
!     ..
!     .. External Subroutines ..
EXTERNAL    SCOPY
EXTERNAL    SLACPY,   SLAQSB,   SLASET
EXTERNAL    SPBEQU,   SPBSV,   SPBSVX
EXTERNAL    SPBTRF,   SPBTRS,   SSWAP
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
      DATA facts/'F' , 'N' , 'E'/
      DATA equeds/'N' , 'Y'/
!     ..
!     .. Executable Statements ..
!
!     Initialize constants and the random number seed.
!
      path(1:1) = 'Single precision'
      path(2:3) = 'PB'
      nrun = 0
      nfail = 0
      nerrs = 0
      DO i = 1 , 4
         iseed(i) = iseedy(i)
      ENDDO
!
!     Test the error exits
!
      IF ( Tsterr ) CALL serrvx(path,Nout)
      INFot = 0
      kdval(1) = 0
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
!
!        Set limits on the number of loop iterations.
!
         nkd = MAX(1,MIN(n,4))
         nimat = NTYPES
         IF ( n==0 ) nimat = 1
!
         kdval(2) = n + (n+1)/4
         kdval(3) = (3*n-1)/4
         kdval(4) = (n+1)/4
!
         DO ikd = 1 , nkd
!
!           Do for KD = 0, (5*N+1)/4, (3N-1)/4, and (N+1)/4. This order
!           makes it easier to skip redundant values for small values
!           of N.
!
            kd = kdval(ikd)
            ldab = kd + 1
!
!           Do first for UPLO = 'U', then for UPLO = 'L'
!
            DO iuplo = 1 , 2
               koff = 1
               IF ( iuplo==1 ) THEN
                  uplo = 'U'
                  packit = 'Q'
                  koff = MAX(1,kd+2-n)
               ELSE
                  uplo = 'L'
                  packit = 'B'
               ENDIF
!
               DO imat = 1 , nimat
!
!                 Do the tests only if DOTYPE( IMAT ) is true.
!
                  IF ( Dotype(imat) ) THEN
!
!                 Skip types 2, 3, or 4 if the matrix size is too small.
!
                     zerot = imat>=2 .AND. imat<=4
                     IF ( .NOT.(zerot .AND. n<imat-1) ) THEN
!
                        IF ( .NOT.zerot .OR. .NOT.Dotype(1) ) THEN
!
!                    Set up parameters with slatb4 and generate a test
!                    matrix with SLATMS.
!
                           CALL slatb4(path,imat,n,n,type,kl,ku,anorm,  &
     &                                 mode,cndnum,dist)
!
                           SRNamt = 'SLATMS'
                           CALL SLATMS(n,n,dist,iseed,type,Rwork,mode,  &
     &                                 cndnum,anorm,kd,kd,packit,A(koff)&
     &                                 ,ldab,Work,info)
!
!                    Check error code from SLATMS.
!
                           IF ( info/=0 ) THEN
                              CALL alaerh(path,'SLATMS',info,0,uplo,n,n,&
     &                           -1,-1,-1,imat,nfail,nerrs,Nout)
                              CYCLE
                           ENDIF
                        ELSEIF ( izero>0 ) THEN
!
!                    Use the same matrix for types 3 and 4 as for type
!                    2 by copying back the zeroed out column,
!
                           iw = 2*lda + 1
                           IF ( iuplo==1 ) THEN
                              ioff = (izero-1)*ldab + kd + 1
                              CALL SCOPY(izero-i1,Work(iw),1,           &
     &                           A(ioff-izero+i1),1)
                              iw = iw + izero - i1
                              CALL SCOPY(i2-izero+1,Work(iw),1,A(ioff), &
     &                           MAX(ldab-1,1))
                           ELSE
                              ioff = (i1-1)*ldab + 1
                              CALL SCOPY(izero-i1,Work(iw),1,           &
     &                           A(ioff+izero-i1),MAX(ldab-1,1))
                              ioff = (izero-1)*ldab + 1
                              iw = iw + izero - i1
                              CALL SCOPY(i2-izero+1,Work(iw),1,A(ioff), &
     &                           1)
                           ENDIF
                        ENDIF
!
!                 For types 2-4, zero one row and column of the matrix
!                 to test that INFO is returned correctly.
!
                        izero = 0
                        IF ( zerot ) THEN
                           IF ( imat==2 ) THEN
                              izero = 1
                           ELSEIF ( imat==3 ) THEN
                              izero = n
                           ELSE
                              izero = n/2 + 1
                           ENDIF
!
!                    Save the zeroed out row and column in WORK(*,3)
!
                           iw = 2*lda
                           DO i = 1 , MIN(2*kd+1,n)
                              Work(iw+i) = ZERO
                           ENDDO
                           iw = iw + 1
                           i1 = MAX(izero-kd,1)
                           i2 = MIN(izero+kd,n)
!
                           IF ( iuplo==1 ) THEN
                              ioff = (izero-1)*ldab + kd + 1
                              CALL SSWAP(izero-i1,A(ioff-izero+i1),1,   &
     &                           Work(iw),1)
                              iw = iw + izero - i1
                              CALL SSWAP(i2-izero+1,A(ioff),            &
     &                           MAX(ldab-1,1),Work(iw),1)
                           ELSE
                              ioff = (i1-1)*ldab + 1
                              CALL SSWAP(izero-i1,A(ioff+izero-i1),     &
     &                           MAX(ldab-1,1),Work(iw),1)
                              ioff = (izero-1)*ldab + 1
                              iw = iw + izero - i1
                              CALL SSWAP(i2-izero+1,A(ioff),1,Work(iw), &
     &                           1)
                           ENDIF
                        ENDIF
!
!                 Save a copy of the matrix A in ASAV.
!
                        CALL SLACPY('Full',kd+1,n,A,ldab,Asav,ldab)
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
!                          Compute the condition number for comparison
!                          with the value returned by SPBSVX (FACT =
!                          'N' reuses the condition number from the
!                          previous iteration with FACT = 'F').
!
                                 CALL SLACPY('Full',kd+1,n,Asav,ldab,   &
     &                              Afac,ldab)
                                 IF ( equil .OR. iequed>1 ) THEN
!
!                             Compute row and column scale factors to
!                             equilibrate the matrix A.
!
                                    CALL SPBEQU(uplo,n,kd,Afac,ldab,S,  &
     &                                 scond,amax,info)
                                    IF ( info==0 .AND. n>0 ) THEN
                                       IF ( iequed>1 ) scond = ZERO
!
!                                Equilibrate the matrix.
!
                                       CALL SLAQSB(uplo,n,kd,Afac,ldab, &
     &                                    S,scond,amax,equed)
                                    ENDIF
                                 ENDIF
!
!                          Save the condition number of the
!                          non-equilibrated system for use in sget04.
!
                                 IF ( equil ) roldc = rcondc
!
!                          Compute the 1-norm of A.
!
                                 anorm = SLANSB('1',uplo,n,kd,Afac,ldab,&
     &                              Rwork)
!
!                          Factor the matrix A.
!
                                 CALL SPBTRF(uplo,n,kd,Afac,ldab,info)
!
!                          Form the inverse of A.
!
                                 CALL SLASET('Full',n,n,ZERO,ONE,A,lda)
                                 SRNamt = 'SPBTRS'
                                 CALL SPBTRS(uplo,n,kd,n,Afac,ldab,A,   &
     &                              lda,info)
!
!                          Compute the 1-norm condition number of A.
!
                                 ainvnm = SLANGE('1',n,n,A,lda,Rwork)
                                 IF ( anorm<=ZERO .OR. ainvnm<=ZERO )   &
     &                                THEN
                                    rcondc = ONE
                                 ELSE
                                    rcondc = (ONE/anorm)/ainvnm
                                 ENDIF
                              ENDIF
!
!                       Restore the matrix A.
!
                              CALL SLACPY('Full',kd+1,n,Asav,ldab,A,    &
     &                           ldab)
!
!                       Form an exact solution and set the right hand
!                       side.
!
                              SRNamt = 'slarhs'
                              CALL slarhs(path,xtype,uplo,' ',n,n,kd,kd,&
     &                           Nrhs,A,ldab,Xact,lda,B,lda,iseed,info)
                              xtype = 'C'
                              CALL SLACPY('Full',n,Nrhs,B,lda,Bsav,lda)
!
                              IF ( nofact ) THEN
!
!                          --- Test SPBSV  ---
!
!                          Compute the L*L' or U'*U factorization of the
!                          matrix and solve the system.
!
                                 CALL SLACPY('Full',kd+1,n,A,ldab,Afac, &
     &                              ldab)
                                 CALL SLACPY('Full',n,Nrhs,B,lda,X,lda)
!
                                 SRNamt = 'SPBSV '
                                 CALL SPBSV(uplo,n,kd,Nrhs,Afac,ldab,X, &
     &                              lda,info)
!
!                          Check error code from SPBSV .
!
                                 IF ( info/=izero ) THEN
                                    CALL alaerh(path,'SPBSV ',info,     &
     &                                 izero,uplo,n,n,kd,kd,Nrhs,imat,  &
     &                                 nfail,nerrs,Nout)
                                    GOTO 2
                                 ELSEIF ( info/=0 ) THEN
                                    GOTO 2
                                 ENDIF
!
!                          Reconstruct matrix from factors and compute
!                          residual.
!
                                 CALL spbt01(uplo,n,kd,A,ldab,Afac,ldab,&
     &                              Rwork,result(1))
!
!                          Compute residual of the computed solution.
!
                                 CALL SLACPY('Full',n,Nrhs,B,lda,Work,  &
     &                              lda)
                                 CALL spbt02(uplo,n,kd,Nrhs,A,ldab,X,   &
     &                              lda,Work,lda,Rwork,result(2))
!
!                          Check solution from generated exact solution.
!
                                 CALL sget04(n,Nrhs,X,lda,Xact,lda,     &
     &                              rcondc,result(3))
                                 nt = 3
!
!                          Print information about the tests that did
!                          not pass the threshold.
!
                                 DO k = 1 , nt
                                    IF ( result(k)>=Thresh ) THEN
                                       IF ( nfail==0 .AND. nerrs==0 )   &
     &                                    CALL aladhd(Nout,path)
                                       WRITE (Nout,FMT=99001) 'SPBSV ' ,&
     &                                    uplo , n , kd , imat , k ,    &
     &                                    result(k)
                                       nfail = nfail + 1
                                    ENDIF
                                 ENDDO
                                 nrun = nrun + nt
                              ENDIF
!
!                       --- Test SPBSVX ---
!
 2                            IF ( .NOT.prefac )                        &
     &                             CALL SLASET('Full',kd+1,n,ZERO,ZERO, &
     &                             Afac,ldab)
                              CALL SLASET('Full',n,Nrhs,ZERO,ZERO,X,lda)
!
!                          Equilibrate the matrix if FACT='F' and
!                          EQUED='Y'
!
                              IF ( iequed>1 .AND. n>0 )                 &
     &                             CALL SLAQSB(uplo,n,kd,A,ldab,S,scond,&
     &                             amax,equed)
!
!                       Solve the system and compute the condition
!                       number and error bounds using SPBSVX.
!
                              SRNamt = 'SPBSVX'
                              CALL SPBSVX(fact,uplo,n,kd,Nrhs,A,ldab,   &
     &                           Afac,ldab,equed,S,B,lda,X,lda,rcond,   &
     &                           Rwork,Rwork(Nrhs+1),Work,Iwork,info)
!
!                       Check the error code from SPBSVX.
!
                              IF ( info/=izero ) THEN
                                 CALL alaerh(path,'SPBSVX',info,izero,  &
     &                              fact//uplo,n,n,kd,kd,Nrhs,imat,     &
     &                              nfail,nerrs,Nout)
                                 CYCLE
                              ENDIF
!
                              IF ( info==0 ) THEN
                                 IF ( .NOT.prefac ) THEN
!
!                             Reconstruct matrix from factors and
!                             compute residual.
!
                                    CALL spbt01(uplo,n,kd,A,ldab,Afac,  &
     &                                 ldab,Rwork(2*Nrhs+1),result(1))
                                    k1 = 1
                                 ELSE
                                    k1 = 2
                                 ENDIF
!
!                          Compute residual of the computed solution.
!
                                 CALL SLACPY('Full',n,Nrhs,Bsav,lda,    &
     &                              Work,lda)
                                 CALL spbt02(uplo,n,kd,Nrhs,Asav,ldab,X,&
     &                              lda,Work,lda,Rwork(2*Nrhs+1),       &
     &                              result(2))
!
!                          Check solution from generated exact solution.
!
                                 IF ( nofact .OR.                       &
     &                                (prefac .AND. LSAME(equed,'N')) ) &
     &                                THEN
                                    CALL sget04(n,Nrhs,X,lda,Xact,lda,  &
     &                                 rcondc,result(3))
                                 ELSE
                                    CALL sget04(n,Nrhs,X,lda,Xact,lda,  &
     &                                 roldc,result(3))
                                 ENDIF
!
!                          Check the error bounds from iterative
!                          refinement.
!
                                 CALL spbt05(uplo,n,kd,Nrhs,Asav,ldab,B,&
     &                              lda,X,lda,Xact,lda,Rwork,           &
     &                              Rwork(Nrhs+1),result(4))
                              ELSE
                                 k1 = 6
                              ENDIF
!
!                       Compare RCOND from SPBSVX with the computed
!                       value in RCONDC.
!
                              result(6) = sget06(rcond,rcondc)
!
!                       Print information about the tests that did not
!                       pass the threshold.
!
                              DO k = k1 , 6
                                 IF ( result(k)>=Thresh ) THEN
                                    IF ( nfail==0 .AND. nerrs==0 )      &
     &                                 CALL aladhd(Nout,path)
                                    IF ( prefac ) THEN
                                       WRITE (Nout,FMT=99003) 'SPBSVX' ,&
     &                                    fact , uplo , n , kd , equed ,&
     &                                    imat , k , result(k)
                                    ELSE
                                       WRITE (Nout,FMT=99002) 'SPBSVX' ,&
     &                                    fact , uplo , n , kd , imat , &
     &                                    k , result(k)
                                    ENDIF
                                    nfail = nfail + 1
                                 ENDIF
                              ENDDO
                              nrun = nrun + 7 - k1
                           ENDDO
                        ENDDO
                     ENDIF
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
      ENDDO
!
!     Print a summary of the results.
!
      CALL alasvm(path,Nout,nfail,nrun,nerrs)
!
99001 FORMAT (1X,A,', UPLO=''',A1,''', N =',I5,', KD =',I5,', type ',I1,&
     &        ', test(',I1,')=',G12.5)
99002 FORMAT (1X,A,'( ''',A1,''', ''',A1,''', ',I5,', ',I5,             &
     &        ', ... ), type ',I1,', test(',I1,')=',G12.5)
99003 FORMAT (1X,A,'( ''',A1,''', ''',A1,''', ',I5,', ',I5,             &
     &        ', ... ), EQUED=''',A1,''', type ',I1,', test(',I1,')=',  &
     &        G12.5)
!
!     End of sdrvpb
!
      END SUBROUTINE sdrvpb
