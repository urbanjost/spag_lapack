!*==schktb.f90  processed by SPAG 7.51RB at 17:34 on  4 Mar 2022
!> \brief \b schktb
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE schktb( DOTYPE, NN, NVAL, NNS, NSVAL, THRESH, TSTERR,
!                          NMAX, AB, AINV, B, X, XACT, WORK, RWORK, IWORK,
!                          NOUT )
!
!       .. Scalar Arguments ..
!       LOGICAL            TSTERR
!       INTEGER            NMAX, NN, NNS, NOUT
!       REAL               THRESH
!       ..
!       .. Array Arguments ..
!       LOGICAL            DOTYPE( * )
!       INTEGER            IWORK( * ), NSVAL( * ), NVAL( * )
!       REAL               AB( * ), AINV( * ), B( * ), RWORK( * ),
!      $                   WORK( * ), X( * ), XACT( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> schktb tests STBTRS, -RFS, and -CON, and SLATBS.
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
!>          The values of the matrix column dimension N.
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
!>          The leading dimension of the work arrays.
!>          NMAX >= the maximum value of N in NVAL.
!> \endverbatim
!>
!> \param[out] AB
!> \verbatim
!>          AB is REAL array, dimension (NMAX*NMAX)
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
      SUBROUTINE schktb(Dotype,Nn,Nval,Nns,Nsval,Thresh,Tsterr,Nmax,Ab, &
     &                  Ainv,B,X,Xact,Work,Rwork,Iwork,Nout)
      IMPLICIT NONE
!*--SCHKTB158
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      LOGICAL Tsterr
      INTEGER Nmax , Nn , Nns , Nout
      REAL Thresh
!     ..
!     .. Array Arguments ..
      LOGICAL Dotype(*)
      INTEGER Iwork(*) , Nsval(*) , Nval(*)
      REAL Ab(*) , Ainv(*) , B(*) , Rwork(*) , Work(*) , X(*) , Xact(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      INTEGER NTYPE1 , NTYPES
      PARAMETER (NTYPE1=9,NTYPES=17)
      INTEGER NTESTS
      PARAMETER (NTESTS=8)
      INTEGER NTRAN
      PARAMETER (NTRAN=3)
      REAL ONE , ZERO
      PARAMETER (ONE=1.0E+0,ZERO=0.0E+0)
!     ..
!     .. Local Scalars ..
      CHARACTER diag , norm , trans , uplo , xtype
      CHARACTER*3 path
      INTEGER i , idiag , ik , imat , in , info , irhs , itran , iuplo ,&
     &        j , k , kd , lda , ldab , n , nerrs , nfail , nimat ,     &
     &        nimat2 , nk , nrhs , nrun
      REAL ainvnm , anorm , rcond , rcondc , rcondi , rcondo , scale
!     ..
!     .. Local Arrays ..
      CHARACTER transs(NTRAN) , uplos(2)
      INTEGER iseed(4) , iseedy(4)
      REAL result(NTESTS)
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      REAL SLANTB , SLANTR
EXTERNAL    LSAME,   SLANTB,   SLANTR
!     ..
!     .. External Subroutines ..
EXTERNAL    SCOPY
EXTERNAL    SLACPY,   SLASET,   SLATBS,   STBCON
EXTERNAL    STBRFS,   STBSV
EXTERNAL    STBTRS
!     ..
!     .. Scalars in Common ..
      LOGICAL LERr , OK
      CHARACTER*32 SRNamt
      INTEGER INFot , IOUnit
!     ..
!     .. Common blocks ..
      COMMON /INFOC / INFot , IOUnit , OK , LERr
      COMMON /SRNAMC/ SRNamt
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX , MIN
!     ..
!     .. Data statements ..
      DATA iseedy/1988 , 1989 , 1990 , 1991/
      DATA uplos/'U' , 'L'/ , transs/'N' , 'T' , 'C'/
!     ..
!     .. Executable Statements ..
!
!     Initialize constants and the random number seed.
!
      path(1:1) = 'Single precision'
      path(2:3) = 'TB'
      nrun = 0
      nfail = 0
      nerrs = 0
      DO i = 1 , 4
         iseed(i) = iseedy(i)
      ENDDO
!
!     Test the error exits
!
      IF ( Tsterr ) CALL serrtr(path,Nout)
      INFot = 0
!
      DO in = 1 , Nn
!
!        Do for each value of N in NVAL
!
         n = Nval(in)
         lda = MAX(1,n)
         xtype = 'N'
         nimat = NTYPE1
         nimat2 = NTYPES
         IF ( n<=0 ) THEN
            nimat = 1
            nimat2 = NTYPE1 + 1
         ENDIF
!
         nk = MIN(n+1,4)
         DO ik = 1 , nk
!
!           Do for KD = 0, N, (3N-1)/4, and (N+1)/4. This order makes
!           it easier to skip redundant values for small values of N.
!
            IF ( ik==1 ) THEN
               kd = 0
            ELSEIF ( ik==2 ) THEN
               kd = MAX(n,0)
            ELSEIF ( ik==3 ) THEN
               kd = (3*n-1)/4
            ELSEIF ( ik==4 ) THEN
               kd = (n+1)/4
            ENDIF
            ldab = kd + 1
!
            DO imat = 1 , nimat
!
!              Do the tests only if DOTYPE( IMAT ) is true.
!
               IF ( Dotype(imat) ) THEN
!
                  DO iuplo = 1 , 2
!
!                 Do first for UPLO = 'U', then for UPLO = 'L'
!
                     uplo = uplos(iuplo)
!
!                 Call slattb to generate a triangular test matrix.
!
                     SRNamt = 'slattb'
                     CALL slattb(imat,uplo,'No transpose',diag,iseed,n, &
     &                           kd,Ab,ldab,X,Work,info)
!
!                 Set IDIAG = 1 for non-unit matrices, 2 for unit.
!
                     IF ( LSAME(diag,'N') ) THEN
                        idiag = 1
                     ELSE
                        idiag = 2
                     ENDIF
!
!                 Form the inverse of A so we can get a good estimate
!                 of RCONDC = 1/(norm(A) * norm(inv(A))).
!
                     CALL SLASET('Full',n,n,ZERO,ONE,Ainv,lda)
                     IF ( LSAME(uplo,'U') ) THEN
                        DO j = 1 , n
                           CALL STBSV(uplo,'No transpose',diag,j,kd,Ab, &
     &                                ldab,Ainv((j-1)*lda+1),1)
                        ENDDO
                     ELSE
                        DO j = 1 , n
                           CALL STBSV(uplo,'No transpose',diag,n-j+1,kd,&
     &                                Ab((j-1)*ldab+1),ldab,            &
     &                                Ainv((j-1)*lda+j),1)
                        ENDDO
                     ENDIF
!
!                 Compute the 1-norm condition number of A.
!
                     anorm = SLANTB('1',uplo,diag,n,kd,Ab,ldab,Rwork)
                     ainvnm = SLANTR('1',uplo,diag,n,n,Ainv,lda,Rwork)
                     IF ( anorm<=ZERO .OR. ainvnm<=ZERO ) THEN
                        rcondo = ONE
                     ELSE
                        rcondo = (ONE/anorm)/ainvnm
                     ENDIF
!
!                 Compute the infinity-norm condition number of A.
!
                     anorm = SLANTB('I',uplo,diag,n,kd,Ab,ldab,Rwork)
                     ainvnm = SLANTR('I',uplo,diag,n,n,Ainv,lda,Rwork)
                     IF ( anorm<=ZERO .OR. ainvnm<=ZERO ) THEN
                        rcondi = ONE
                     ELSE
                        rcondi = (ONE/anorm)/ainvnm
                     ENDIF
!
                     DO irhs = 1 , Nns
                        nrhs = Nsval(irhs)
                        xtype = 'N'
!
                        DO itran = 1 , NTRAN
!
!                    Do for op(A) = A, A**T, or A**H.
!
                           trans = transs(itran)
                           IF ( itran==1 ) THEN
                              norm = 'O'
                              rcondc = rcondo
                           ELSE
                              norm = 'I'
                              rcondc = rcondi
                           ENDIF
!
!+    TEST 1
!                    Solve and compute residual for op(A)*x = b.
!
                           SRNamt = 'slarhs'
                           CALL slarhs(path,xtype,uplo,trans,n,n,kd,    &
     &                                 idiag,nrhs,Ab,ldab,Xact,lda,B,   &
     &                                 lda,iseed,info)
                           xtype = 'C'
                           CALL SLACPY('Full',n,nrhs,B,lda,X,lda)
!
                           SRNamt = 'STBTRS'
                           CALL STBTRS(uplo,trans,diag,n,kd,nrhs,Ab,    &
     &                                 ldab,X,lda,info)
!
!                    Check error code from STBTRS.
!
                           IF ( info/=0 )                               &
     &                          CALL alaerh(path,'STBTRS',info,0,       &
     &                          uplo//trans//diag,n,n,kd,kd,nrhs,imat,  &
     &                          nfail,nerrs,Nout)
!
                           CALL stbt02(uplo,trans,diag,n,kd,nrhs,Ab,    &
     &                                 ldab,X,lda,B,lda,Work,result(1))
!
!+    TEST 2
!                    Check solution from generated exact solution.
!
                           CALL sget04(n,nrhs,X,lda,Xact,lda,rcondc,    &
     &                                 result(2))
!
!+    TESTS 3, 4, and 5
!                    Use iterative refinement to improve the solution
!                    and compute error bounds.
!
                           SRNamt = 'STBRFS'
                           CALL STBRFS(uplo,trans,diag,n,kd,nrhs,Ab,    &
     &                                 ldab,B,lda,X,lda,Rwork,          &
     &                                 Rwork(nrhs+1),Work,Iwork,info)
!
!                    Check error code from STBRFS.
!
                           IF ( info/=0 )                               &
     &                          CALL alaerh(path,'STBRFS',info,0,       &
     &                          uplo//trans//diag,n,n,kd,kd,nrhs,imat,  &
     &                          nfail,nerrs,Nout)
!
                           CALL sget04(n,nrhs,X,lda,Xact,lda,rcondc,    &
     &                                 result(3))
                           CALL stbt05(uplo,trans,diag,n,kd,nrhs,Ab,    &
     &                                 ldab,B,lda,X,lda,Xact,lda,Rwork, &
     &                                 Rwork(nrhs+1),result(4))
!
!                       Print information about the tests that did not
!                       pass the threshold.
!
                           DO k = 1 , 5
                              IF ( result(k)>=Thresh ) THEN
                                 IF ( nfail==0 .AND. nerrs==0 )         &
     &                                CALL alahd(Nout,path)
                                 WRITE (Nout,FMT=99001) uplo , trans ,  &
     &                                  diag , n , kd , nrhs , imat ,   &
     &                                  k , result(k)
                                 nfail = nfail + 1
                              ENDIF
                           ENDDO
                           nrun = nrun + 5
                        ENDDO
                     ENDDO
!
!+    TEST 6
!                    Get an estimate of RCOND = 1/CNDNUM.
!
                     DO itran = 1 , 2
                        IF ( itran==1 ) THEN
                           norm = 'O'
                           rcondc = rcondo
                        ELSE
                           norm = 'I'
                           rcondc = rcondi
                        ENDIF
                        SRNamt = 'STBCON'
                        CALL STBCON(norm,uplo,diag,n,kd,Ab,ldab,rcond,  &
     &                              Work,Iwork,info)
!
!                    Check error code from STBCON.
!
                        IF ( info/=0 )                                  &
     &                       CALL alaerh(path,'STBCON',info,0,norm//    &
     &                       uplo//diag,n,n,kd,kd,-1,imat,nfail,nerrs,  &
     &                       Nout)
!
                        CALL stbt06(rcond,rcondc,uplo,diag,n,kd,Ab,ldab,&
     &                              Rwork,result(6))
!
!                    Print information about the tests that did not pass
!                    the threshold.
!
                        IF ( result(6)>=Thresh ) THEN
                           IF ( nfail==0 .AND. nerrs==0 )               &
     &                          CALL alahd(Nout,path)
                           WRITE (Nout,FMT=99002) 'STBCON' , norm ,     &
     &                            uplo , diag , n , kd , imat , 6 ,     &
     &                            result(6)
                           nfail = nfail + 1
                        ENDIF
                        nrun = nrun + 1
                     ENDDO
                  ENDDO
               ENDIF
            ENDDO
!
!           Use pathological test matrices to test SLATBS.
!
            DO imat = NTYPE1 + 1 , nimat2
!
!              Do the tests only if DOTYPE( IMAT ) is true.
!
               IF ( Dotype(imat) ) THEN
!
                  DO iuplo = 1 , 2
!
!                 Do first for UPLO = 'U', then for UPLO = 'L'
!
                     uplo = uplos(iuplo)
                     DO itran = 1 , NTRAN
!
!                    Do for op(A) = A, A**T, and A**H.
!
                        trans = transs(itran)
!
!                    Call slattb to generate a triangular test matrix.
!
                        SRNamt = 'slattb'
                        CALL slattb(imat,uplo,trans,diag,iseed,n,kd,Ab, &
     &                              ldab,X,Work,info)
!
!+    TEST 7
!                    Solve the system op(A)*x = b
!
                        SRNamt = 'SLATBS'
                        CALL SCOPY(n,X,1,B,1)
                        CALL SLATBS(uplo,trans,diag,'N',n,kd,Ab,ldab,B, &
     &                              scale,Rwork,info)
!
!                    Check error code from SLATBS.
!
                        IF ( info/=0 )                                  &
     &                       CALL alaerh(path,'SLATBS',info,0,uplo//    &
     &                       trans//diag//'N',n,n,kd,kd,-1,imat,nfail,  &
     &                       nerrs,Nout)
!
                        CALL stbt03(uplo,trans,diag,n,kd,1,Ab,ldab,     &
     &                              scale,Rwork,ONE,B,lda,X,lda,Work,   &
     &                              result(7))
!
!+    TEST 8
!                    Solve op(A)*x = b again with NORMIN = 'Y'.
!
                        CALL SCOPY(n,X,1,B,1)
                        CALL SLATBS(uplo,trans,diag,'Y',n,kd,Ab,ldab,B, &
     &                              scale,Rwork,info)
!
!                    Check error code from SLATBS.
!
                        IF ( info/=0 )                                  &
     &                       CALL alaerh(path,'SLATBS',info,0,uplo//    &
     &                       trans//diag//'Y',n,n,kd,kd,-1,imat,nfail,  &
     &                       nerrs,Nout)
!
                        CALL stbt03(uplo,trans,diag,n,kd,1,Ab,ldab,     &
     &                              scale,Rwork,ONE,B,lda,X,lda,Work,   &
     &                              result(8))
!
!                    Print information about the tests that did not pass
!                    the threshold.
!
                        IF ( result(7)>=Thresh ) THEN
                           IF ( nfail==0 .AND. nerrs==0 )               &
     &                          CALL alahd(Nout,path)
                           WRITE (Nout,FMT=99003) 'SLATBS' , uplo ,     &
     &                            trans , diag , 'N' , n , kd , imat ,  &
     &                            7 , result(7)
                           nfail = nfail + 1
                        ENDIF
                        IF ( result(8)>=Thresh ) THEN
                           IF ( nfail==0 .AND. nerrs==0 )               &
     &                          CALL alahd(Nout,path)
                           WRITE (Nout,FMT=99003) 'SLATBS' , uplo ,     &
     &                            trans , diag , 'Y' , n , kd , imat ,  &
     &                            8 , result(8)
                           nfail = nfail + 1
                        ENDIF
                        nrun = nrun + 2
                     ENDDO
                  ENDDO
               ENDIF
            ENDDO
         ENDDO
      ENDDO
!
!     Print a summary of the results.
!
      CALL alasum(path,Nout,nfail,nrun,nerrs)
!
99001 FORMAT (' UPLO=''',A1,''', TRANS=''',A1,                          &
     &''',                                                              &
     &                                                                  &
     &                                                                  &
     &        DIAG=''',A1,''', N=',I5,', KD=',I5,', NRHS=',I5,', type ',&
     &I2,', test(',I2,')=',G12.5)
99002 FORMAT (1X,A,'( ''',A1,''', ''',A1,''', ''',A1,''',',I5,',',I5,   &
     &        ',  ... ), type ',I2,', test(',I2,')=',G12.5)
99003 FORMAT (1X,A,'( ''',A1,''', ''',A1,''', ''',A1,''', ''',A1,''',', &
     &        I5,',',I5,', ...  ),  type ',I2,', test(',I1,')=',G12.5)
!
!     End of schktb
!
      END SUBROUTINE schktb
