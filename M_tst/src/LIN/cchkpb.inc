!*==cchkpb.f90  processed by SPAG 7.51RB at 17:34 on  4 Mar 2022
!> \brief \b cchkpb
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE cchkpb( DOTYPE, NN, NVAL, NNB, NBVAL, NNS, NSVAL,
!                          THRESH, TSTERR, NMAX, A, AFAC, AINV, B, X,
!                          XACT, WORK, RWORK, NOUT )
!
!       .. Scalar Arguments ..
!       LOGICAL            TSTERR
!       INTEGER            NMAX, NN, NNB, NNS, NOUT
!       REAL               THRESH
!       ..
!       .. Array Arguments ..
!       LOGICAL            DOTYPE( * )
!       INTEGER            NBVAL( * ), NSVAL( * ), NVAL( * )
!       REAL               RWORK( * )
!       COMPLEX            A( * ), AFAC( * ), AINV( * ), B( * ),
!      $                   WORK( * ), X( * ), XACT( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> cchkpb tests CPBTRF, -TRS, -RFS, and -CON.
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
!> \param[in] NNB
!> \verbatim
!>          NNB is INTEGER
!>          The number of values of NB contained in the vector NBVAL.
!> \endverbatim
!>
!> \param[in] NBVAL
!> \verbatim
!>          NBVAL is INTEGER array, dimension (NBVAL)
!>          The values of the blocksize NB.
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
!> \param[out] AINV
!> \verbatim
!>          AINV is REAL array, dimension (NMAX*NMAX)
!> \endverbatim
!>
!> \param[out] B
!> \verbatim
!>          B is REAL array, dimension (NMAX*NSMAX)
!>          where NSMAX is the largest entry in NSVAL.
!> \endverbatim
!>
!> \param[out] X
!> \verbatim
!>          X is REAL array, dimension (NMAX*NSMAX)
!> \endverbatim
!>
!> \param[out] XACT
!> \verbatim
!>          XACT is REAL array, dimension (NMAX*NSMAX)
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is REAL array, dimension
!>                      (NMAX*max(3,NSMAX))
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is REAL array, dimension
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
!> \ingroup complex_lin
!
!  =====================================================================
      SUBROUTINE cchkpb(Dotype,Nn,Nval,Nnb,Nbval,Nns,Nsval,Thresh,      &
     &                  Tsterr,Nmax,A,Afac,Ainv,B,X,Xact,Work,Rwork,    &
     &                  Nout)
      IMPLICIT NONE
!*--CCHKPB172
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      LOGICAL Tsterr
      INTEGER Nmax , Nn , Nnb , Nns , Nout
      REAL Thresh
!     ..
!     .. Array Arguments ..
      LOGICAL Dotype(*)
      INTEGER Nbval(*) , Nsval(*) , Nval(*)
      REAL Rwork(*)
      COMPLEX A(*) , Afac(*) , Ainv(*) , B(*) , Work(*) , X(*) , Xact(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ONE , ZERO
      PARAMETER (ONE=1.0E+0,ZERO=0.0E+0)
      INTEGER NTYPES , NTESTS
      PARAMETER (NTYPES=8,NTESTS=7)
      INTEGER NBW
      PARAMETER (NBW=4)
!     ..
!     .. Local Scalars ..
      LOGICAL zerot
      CHARACTER dist , packit , type , uplo , xtype
      CHARACTER*3 path
      INTEGER i , i1 , i2 , ikd , imat , in , inb , info , ioff , irhs ,&
     &        iuplo , iw , izero , k , kd , kl , koff , ku , lda ,      &
     &        ldab , mode , n , nb , nerrs , nfail , nimat , nkd ,      &
     &        nrhs , nrun
      REAL ainvnm , anorm , cndnum , rcond , rcondc
!     ..
!     .. Local Arrays ..
      INTEGER iseed(4) , iseedy(4) , kdval(NBW)
      REAL result(NTESTS)
!     ..
!     .. External Functions ..
      REAL CLANGE , CLANHB  
EXTERNAL    CLANGE,   CLANHB
!     ..
!     .. External Subroutines ..
EXTERNAL    CCOPY
EXTERNAL    CLACPY,   CLASET
EXTERNAL    CPBCON,   CPBRFS,   CPBTRF
EXTERNAL    CPBTRS,   CSWAP
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC CMPLX , MAX , MIN
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
!     ..
!     .. Executable Statements ..
!
!     Initialize constants and the random number seed.
!
      path(1:1) = 'Complex precision'
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
      IF ( Tsterr ) CALL cerrpo(path,Nout)
      INFot = 0
      kdval(1) = 0
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
                  koff = MAX(1,kd+2-n)
                  packit = 'Q'
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
!                    Set up parameters with clatb4 and generate a test
!                    matrix with CLATMS.
!
                           CALL clatb4(path,imat,n,n,type,kl,ku,anorm,  &
     &                                 mode,cndnum,dist)
!
                           SRNamt = 'CLATMS'
                           CALL CLATMS(n,n,dist,iseed,type,Rwork,mode,  &
     &                                 cndnum,anorm,kd,kd,packit,A(koff)&
     &                                 ,ldab,Work,info)
!
!                    Check error code from CLATMS.
!
                           IF ( info/=0 ) THEN
                              CALL alaerh(path,'CLATMS',info,0,uplo,n,n,&
     &                           kd,kd,-1,imat,nfail,nerrs,Nout)
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
                              CALL CCOPY(izero-i1,Work(iw),1,           &
     &                           A(ioff-izero+i1),1)
                              iw = iw + izero - i1
                              CALL CCOPY(i2-izero+1,Work(iw),1,A(ioff), &
     &                           MAX(ldab-1,1))
                           ELSE
                              ioff = (i1-1)*ldab + 1
                              CALL CCOPY(izero-i1,Work(iw),1,           &
     &                           A(ioff+izero-i1),MAX(ldab-1,1))
                              ioff = (izero-1)*ldab + 1
                              iw = iw + izero - i1
                              CALL CCOPY(i2-izero+1,Work(iw),1,A(ioff), &
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
                              CALL CSWAP(izero-i1,A(ioff-izero+i1),1,   &
     &                           Work(iw),1)
                              iw = iw + izero - i1
                              CALL CSWAP(i2-izero+1,A(ioff),            &
     &                           MAX(ldab-1,1),Work(iw),1)
                           ELSE
                              ioff = (i1-1)*ldab + 1
                              CALL CSWAP(izero-i1,A(ioff+izero-i1),     &
     &                           MAX(ldab-1,1),Work(iw),1)
                              ioff = (izero-1)*ldab + 1
                              iw = iw + izero - i1
                              CALL CSWAP(i2-izero+1,A(ioff),1,Work(iw), &
     &                           1)
                           ENDIF
                        ENDIF
!
!                 Set the imaginary part of the diagonals.
!
                        IF ( iuplo==1 ) THEN
                           CALL claipd(n,A(kd+1),ldab,0)
                        ELSE
                           CALL claipd(n,A(1),ldab,0)
                        ENDIF
!
!                 Do for each value of NB in NBVAL
!
                        DO inb = 1 , Nnb
                           nb = Nbval(inb)
                           CALL xlaenv(1,nb)
!
!                    Compute the L*L' or U'*U factorization of the band
!                    matrix.
!
                           CALL CLACPY('Full',kd+1,n,A,ldab,Afac,ldab)
                           SRNamt = 'CPBTRF'
                           CALL CPBTRF(uplo,n,kd,Afac,ldab,info)
!
!                    Check error code from CPBTRF.
!
                           IF ( info/=izero ) THEN
                              CALL alaerh(path,'CPBTRF',info,izero,uplo,&
     &                           n,n,kd,kd,nb,imat,nfail,nerrs,Nout)
                              CYCLE
                           ENDIF
!
!                    Skip the tests if INFO is not 0.
!
                           IF ( info==0 ) THEN
!
!+    TEST 1
!                    Reconstruct matrix from factors and compute
!                    residual.
!
                              CALL CLACPY('Full',kd+1,n,Afac,ldab,Ainv, &
     &                           ldab)
                              CALL cpbt01(uplo,n,kd,A,ldab,Ainv,ldab,   &
     &                           Rwork,result(1))
!
!                    Print the test ratio if it is .GE. THRESH.
!
                              IF ( result(1)>=Thresh ) THEN
                                 IF ( nfail==0 .AND. nerrs==0 )         &
     &                                CALL alahd(Nout,path)
                                 WRITE (Nout,FMT=99001) uplo , n , kd , &
     &                                  nb , imat , 1 , result(1)
                                 nfail = nfail + 1
                              ENDIF
                              nrun = nrun + 1
!
!                    Only do other tests if this is the first blocksize.
!
                              IF ( inb<=1 ) THEN
!
!                    Form the inverse of A so we can get a good estimate
!                    of RCONDC = 1/(norm(A) * norm(inv(A))).
!
                                 CALL CLASET('Full',n,n,CMPLX(ZERO),    &
     &                              CMPLX(ONE),Ainv,lda)
                                 SRNamt = 'CPBTRS'
                                 CALL CPBTRS(uplo,n,kd,n,Afac,ldab,Ainv,&
     &                              lda,info)
!
!                    Compute RCONDC = 1/(norm(A) * norm(inv(A))).
!
                                 anorm = CLANHB('1',uplo,n,kd,A,ldab,   &
     &                              Rwork)
                                 ainvnm = CLANGE('1',n,n,Ainv,lda,Rwork)
                                 IF ( anorm<=ZERO .OR. ainvnm<=ZERO )   &
     &                                THEN
                                    rcondc = ONE
                                 ELSE
                                    rcondc = (ONE/anorm)/ainvnm
                                 ENDIF
!
                                 DO irhs = 1 , Nns
                                    nrhs = Nsval(irhs)
!
!+    TEST 2
!                    Solve and compute residual for A * X = B.
!
                                    SRNamt = 'clarhs'
                                    CALL clarhs(path,xtype,uplo,' ',n,n,&
     &                                 kd,kd,nrhs,A,ldab,Xact,lda,B,lda,&
     &                                 iseed,info)
                                    CALL CLACPY('Full',n,nrhs,B,lda,X,  &
     &                                 lda)
!
                                    SRNamt = 'CPBTRS'
                                    CALL CPBTRS(uplo,n,kd,nrhs,Afac,    &
     &                                 ldab,X,lda,info)
!
!                    Check error code from CPBTRS.
!
                                    IF ( info/=0 )                      &
     &                                  CALL alaerh(path,'CPBTRS',info, &
     &                                 0,uplo,n,n,kd,kd,nrhs,imat,nfail,&
     &                                 nerrs,Nout)
!
                                    CALL CLACPY('Full',n,nrhs,B,lda,    &
     &                                 Work,lda)
                                    CALL cpbt02(uplo,n,kd,nrhs,A,ldab,X,&
     &                                 lda,Work,lda,Rwork,result(2))
!
!+    TEST 3
!                    Check solution from generated exact solution.
!
                                    CALL cget04(n,nrhs,X,lda,Xact,lda,  &
     &                                 rcondc,result(3))
!
!+    TESTS 4, 5, and 6
!                    Use iterative refinement to improve the solution.
!
                                    SRNamt = 'CPBRFS'
                                    CALL CPBRFS(uplo,n,kd,nrhs,A,ldab,  &
     &                                 Afac,ldab,B,lda,X,lda,Rwork,     &
     &                                 Rwork(nrhs+1),Work,              &
     &                                 Rwork(2*nrhs+1),info)
!
!                    Check error code from CPBRFS.
!
                                    IF ( info/=0 )                      &
     &                                  CALL alaerh(path,'CPBRFS',info, &
     &                                 0,uplo,n,n,kd,kd,nrhs,imat,nfail,&
     &                                 nerrs,Nout)
!
                                    CALL cget04(n,nrhs,X,lda,Xact,lda,  &
     &                                 rcondc,result(4))
                                    CALL cpbt05(uplo,n,kd,nrhs,A,ldab,B,&
     &                                 lda,X,lda,Xact,lda,Rwork,        &
     &                                 Rwork(nrhs+1),result(5))
!
!                       Print information about the tests that did not
!                       pass the threshold.
!
                                    DO k = 2 , 6
                                       IF ( result(k)>=Thresh ) THEN
                                         IF ( nfail==0 .AND. nerrs==0 ) &
     &                                      CALL alahd(Nout,path)
                                         WRITE (Nout,FMT=99002) uplo ,  &
     &                                      n , kd , nrhs , imat , k ,  &
     &                                      result(k)
                                         nfail = nfail + 1
                                       ENDIF
                                    ENDDO
                                    nrun = nrun + 5
                                 ENDDO
!
!+    TEST 7
!                    Get an estimate of RCOND = 1/CNDNUM.
!
                                 SRNamt = 'CPBCON'
                                 CALL CPBCON(uplo,n,kd,Afac,ldab,anorm, &
     &                              rcond,Work,Rwork,info)
!
!                    Check error code from CPBCON.
!
                                 IF ( info/=0 )                         &
     &                                 CALL alaerh(path,'CPBCON',info,0,&
     &                                uplo,n,n,kd,kd,-1,imat,nfail,     &
     &                                nerrs,Nout)
!
                                 result(7) = sget06(rcond,rcondc)
!
!                    Print the test ratio if it is .GE. THRESH.
!
                                 IF ( result(7)>=Thresh ) THEN
                                    IF ( nfail==0 .AND. nerrs==0 )      &
     &                                 CALL alahd(Nout,path)
                                    WRITE (Nout,FMT=99003) uplo , n ,   &
     &                                 kd , imat , 7 , result(7)
                                    nfail = nfail + 1
                                 ENDIF
                                 nrun = nrun + 1
                              ENDIF
                           ENDIF
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
      CALL alasum(path,Nout,nfail,nrun,nerrs)
!
99001 FORMAT (' UPLO=''',A1,''', N=',I5,', KD=',I5,', NB=',I4,', type ',&
     &        I2,', test ',I2,', ratio= ',G12.5)
99002 FORMAT (' UPLO=''',A1,''', N=',I5,', KD=',I5,', NRHS=',I3,        &
     &        ', type ',I2,', test(',I2,') = ',G12.5)
99003 FORMAT (' UPLO=''',A1,''', N=',I5,', KD=',I5,',',10X,' type ',I2, &
     &        ', test(',I2,') = ',G12.5)
!
!     End of cchkpb
!
      END SUBROUTINE cchkpb
