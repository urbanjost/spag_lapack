!*==dchktr.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b DCHKTR
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DCHKTR( DOTYPE, NN, NVAL, NNB, NBVAL, NNS, NSVAL,
!                          THRESH, TSTERR, NMAX, A, AINV, B, X, XACT,
!                          WORK, RWORK, IWORK, NOUT )
!
!       .. Scalar Arguments ..
!       LOGICAL            TSTERR
!       INTEGER            NMAX, NN, NNB, NNS, NOUT
!       DOUBLE PRECISION   THRESH
!       ..
!       .. Array Arguments ..
!       LOGICAL            DOTYPE( * )
!       INTEGER            IWORK( * ), NBVAL( * ), NSVAL( * ), NVAL( * )
!       DOUBLE PRECISION   A( * ), AINV( * ), B( * ), RWORK( * ),
!      $                   WORK( * ), X( * ), XACT( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DCHKTR tests DTRTRI, -TRS, -RFS, and -CON, and DLATRS
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
!> \param[in] NNB
!> \verbatim
!>          NNB is INTEGER
!>          The number of values of NB contained in the vector NBVAL.
!> \endverbatim
!>
!> \param[in] NBVAL
!> \verbatim
!>          NBVAL is INTEGER array, dimension (NNB)
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
!>          The leading dimension of the work arrays.
!>          NMAX >= the maximum value of N in NVAL.
!> \endverbatim
!>
!> \param[out] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (NMAX*NMAX)
!> \endverbatim
!>
!> \param[out] AINV
!> \verbatim
!>          AINV is DOUBLE PRECISION array, dimension (NMAX*NMAX)
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
!> \ingroup double_lin
!
!  =====================================================================
      SUBROUTINE DCHKTR(Dotype,Nn,Nval,Nnb,Nbval,Nns,Nsval,Thresh,      &
     &                  Tsterr,Nmax,A,Ainv,B,X,Xact,Work,Rwork,Iwork,   &
     &                  Nout)
      IMPLICIT NONE
!*--DCHKTR171
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      LOGICAL Tsterr
      INTEGER Nmax , Nn , Nnb , Nns , Nout
      DOUBLE PRECISION Thresh
!     ..
!     .. Array Arguments ..
      LOGICAL Dotype(*)
      INTEGER Iwork(*) , Nbval(*) , Nsval(*) , Nval(*)
      DOUBLE PRECISION A(*) , Ainv(*) , B(*) , Rwork(*) , Work(*) ,     &
     &                 X(*) , Xact(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      INTEGER NTYPE1 , NTYPES
      PARAMETER (NTYPE1=10,NTYPES=18)
      INTEGER NTESTS
      PARAMETER (NTESTS=9)
      INTEGER NTRAN
      PARAMETER (NTRAN=3)
      DOUBLE PRECISION ONE , ZERO
      PARAMETER (ONE=1.0D0,ZERO=0.0D0)
!     ..
!     .. Local Scalars ..
      CHARACTER diag , norm , trans , uplo , xtype
      CHARACTER*3 path
      INTEGER i , idiag , imat , in , inb , info , irhs , itran ,       &
     &        iuplo , k , lda , n , nb , nerrs , nfail , nrhs , nrun
      DOUBLE PRECISION ainvnm , anorm , dummy , rcond , rcondc ,        &
     &                 rcondi , rcondo , scale
!     ..
!     .. Local Arrays ..
      CHARACTER transs(NTRAN) , uplos(2)
      INTEGER iseed(4) , iseedy(4)
      DOUBLE PRECISION result(NTESTS)
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      DOUBLE PRECISION DLANTR
      EXTERNAL LSAME , DLANTR
!     ..
!     .. External Subroutines ..
      EXTERNAL ALAERH , ALAHD , ALASUM , DCOPY , DERRTR , DGET04 ,      &
     &         DLACPY , DLARHS , DLATRS , DLATTR , DTRCON , DTRRFS ,    &
     &         DTRT01 , DTRT02 , DTRT03 , DTRT05 , DTRT06 , DTRTRI ,    &
     &         DTRTRS , XLAENV
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
      INTRINSIC MAX
!     ..
!     .. Data statements ..
      DATA iseedy/1988 , 1989 , 1990 , 1991/
      DATA uplos/'U' , 'L'/ , transs/'N' , 'T' , 'C'/
!     ..
!     .. Executable Statements ..
!
!     Initialize constants and the random number seed.
!
      path(1:1) = 'Double precision'
      path(2:3) = 'TR'
      nrun = 0
      nfail = 0
      nerrs = 0
      DO i = 1 , 4
         iseed(i) = iseedy(i)
      ENDDO
!
!     Test the error exits
!
      IF ( Tsterr ) CALL DERRTR(path,Nout)
      INFot = 0
      CALL XLAENV(2,2)
!
      DO in = 1 , Nn
!
!        Do for each value of N in NVAL
!
         n = Nval(in)
         lda = MAX(1,n)
         xtype = 'N'
!
         DO imat = 1 , NTYPE1
!
!           Do the tests only if DOTYPE( IMAT ) is true.
!
            IF ( Dotype(imat) ) THEN
!
               DO iuplo = 1 , 2
!
!              Do first for UPLO = 'U', then for UPLO = 'L'
!
                  uplo = uplos(iuplo)
!
!              Call DLATTR to generate a triangular test matrix.
!
                  SRNamt = 'DLATTR'
                  CALL DLATTR(imat,uplo,'No transpose',diag,iseed,n,A,  &
     &                        lda,X,Work,info)
!
!              Set IDIAG = 1 for non-unit matrices, 2 for unit.
!
                  IF ( LSAME(diag,'N') ) THEN
                     idiag = 1
                  ELSE
                     idiag = 2
                  ENDIF
!
                  DO inb = 1 , Nnb
!
!                 Do for each blocksize in NBVAL
!
                     nb = Nbval(inb)
                     CALL XLAENV(1,nb)
!
!+    TEST 1
!                 Form the inverse of A.
!
                     CALL DLACPY(uplo,n,n,A,lda,Ainv,lda)
                     SRNamt = 'DTRTRI'
                     CALL DTRTRI(uplo,diag,n,Ainv,lda,info)
!
!                 Check error code from DTRTRI.
!
                     IF ( info/=0 )                                     &
     &                    CALL ALAERH(path,'DTRTRI',info,0,uplo//diag,n,&
     &                    n,-1,-1,nb,imat,nfail,nerrs,Nout)
!
!                 Compute the infinity-norm condition number of A.
!
                     anorm = DLANTR('I',uplo,diag,n,n,A,lda,Rwork)
                     ainvnm = DLANTR('I',uplo,diag,n,n,Ainv,lda,Rwork)
                     IF ( anorm<=ZERO .OR. ainvnm<=ZERO ) THEN
                        rcondi = ONE
                     ELSE
                        rcondi = (ONE/anorm)/ainvnm
                     ENDIF
!
!                 Compute the residual for the triangular matrix times
!                 its inverse.  Also compute the 1-norm condition number
!                 of A.
!
                     CALL DTRT01(uplo,diag,n,A,lda,Ainv,lda,rcondo,     &
     &                           Rwork,result(1))
!
!                 Print the test ratio if it is .GE. THRESH.
!
                     IF ( result(1)>=Thresh ) THEN
                        IF ( nfail==0 .AND. nerrs==0 )                  &
     &                       CALL ALAHD(Nout,path)
                        WRITE (Nout,FMT=99001) uplo , diag , n , nb ,   &
     &                         imat , 1 , result(1)
                        nfail = nfail + 1
                     ENDIF
                     nrun = nrun + 1
!
!                 Skip remaining tests if not the first block size.
!
                     IF ( inb==1 ) THEN
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
!+    TEST 2
!                       Solve and compute residual for op(A)*x = b.
!
                              SRNamt = 'DLARHS'
                              CALL DLARHS(path,xtype,uplo,trans,n,n,0,  &
     &                           idiag,nrhs,A,lda,Xact,lda,B,lda,iseed, &
     &                           info)
                              xtype = 'C'
                              CALL DLACPY('Full',n,nrhs,B,lda,X,lda)
!
                              SRNamt = 'DTRTRS'
                              CALL DTRTRS(uplo,trans,diag,n,nrhs,A,lda, &
     &                           X,lda,info)
!
!                       Check error code from DTRTRS.
!
                              IF ( info/=0 )                            &
     &                             CALL ALAERH(path,'DTRTRS',info,0,    &
     &                             uplo//trans//diag,n,n,-1,-1,nrhs,    &
     &                             imat,nfail,nerrs,Nout)
!
!                       This line is needed on a Sun SPARCstation.
!
                              IF ( n>0 ) dummy = A(1)
!
                              CALL DTRT02(uplo,trans,diag,n,nrhs,A,lda, &
     &                           X,lda,B,lda,Work,result(2))
!
!+    TEST 3
!                       Check solution from generated exact solution.
!
                              CALL DGET04(n,nrhs,X,lda,Xact,lda,rcondc, &
     &                           result(3))
!
!+    TESTS 4, 5, and 6
!                       Use iterative refinement to improve the solution
!                       and compute error bounds.
!
                              SRNamt = 'DTRRFS'
                              CALL DTRRFS(uplo,trans,diag,n,nrhs,A,lda, &
     &                           B,lda,X,lda,Rwork,Rwork(nrhs+1),Work,  &
     &                           Iwork,info)
!
!                       Check error code from DTRRFS.
!
                              IF ( info/=0 )                            &
     &                             CALL ALAERH(path,'DTRRFS',info,0,    &
     &                             uplo//trans//diag,n,n,-1,-1,nrhs,    &
     &                             imat,nfail,nerrs,Nout)
!
                              CALL DGET04(n,nrhs,X,lda,Xact,lda,rcondc, &
     &                           result(4))
                              CALL DTRT05(uplo,trans,diag,n,nrhs,A,lda, &
     &                           B,lda,X,lda,Xact,lda,Rwork,            &
     &                           Rwork(nrhs+1),result(5))
!
!                       Print information about the tests that did not
!                       pass the threshold.
!
                              DO k = 2 , 6
                                 IF ( result(k)>=Thresh ) THEN
                                    IF ( nfail==0 .AND. nerrs==0 )      &
     &                                 CALL ALAHD(Nout,path)
                                    WRITE (Nout,FMT=99002) uplo ,       &
     &                                 trans , diag , n , nrhs , imat , &
     &                                 k , result(k)
                                    nfail = nfail + 1
                                 ENDIF
                              ENDDO
                              nrun = nrun + 5
                           ENDDO
                        ENDDO
!
!+    TEST 7
!                       Get an estimate of RCOND = 1/CNDNUM.
!
                        DO itran = 1 , 2
                           IF ( itran==1 ) THEN
                              norm = 'O'
                              rcondc = rcondo
                           ELSE
                              norm = 'I'
                              rcondc = rcondi
                           ENDIF
                           SRNamt = 'DTRCON'
                           CALL DTRCON(norm,uplo,diag,n,A,lda,rcond,    &
     &                                 Work,Iwork,info)
!
!                       Check error code from DTRCON.
!
                           IF ( info/=0 )                               &
     &                          CALL ALAERH(path,'DTRCON',info,0,       &
     &                          norm//uplo//diag,n,n,-1,-1,-1,imat,     &
     &                          nfail,nerrs,Nout)
!
                           CALL DTRT06(rcond,rcondc,uplo,diag,n,A,lda,  &
     &                                 Rwork,result(7))
!
!                    Print the test ratio if it is .GE. THRESH.
!
                           IF ( result(7)>=Thresh ) THEN
                              IF ( nfail==0 .AND. nerrs==0 )            &
     &                             CALL ALAHD(Nout,path)
                              WRITE (Nout,FMT=99003) norm , uplo , n ,  &
     &                               imat , 7 , result(7)
                              nfail = nfail + 1
                           ENDIF
                           nrun = nrun + 1
                        ENDDO
                     ENDIF
                  ENDDO
               ENDDO
            ENDIF
         ENDDO
!
!        Use pathological test matrices to test DLATRS.
!
         DO imat = NTYPE1 + 1 , NTYPES
!
!           Do the tests only if DOTYPE( IMAT ) is true.
!
            IF ( Dotype(imat) ) THEN
!
               DO iuplo = 1 , 2
!
!              Do first for UPLO = 'U', then for UPLO = 'L'
!
                  uplo = uplos(iuplo)
                  DO itran = 1 , NTRAN
!
!                 Do for op(A) = A, A**T, and A**H.
!
                     trans = transs(itran)
!
!                 Call DLATTR to generate a triangular test matrix.
!
                     SRNamt = 'DLATTR'
                     CALL DLATTR(imat,uplo,trans,diag,iseed,n,A,lda,X,  &
     &                           Work,info)
!
!+    TEST 8
!                 Solve the system op(A)*x = b.
!
                     SRNamt = 'DLATRS'
                     CALL DCOPY(n,X,1,B,1)
                     CALL DLATRS(uplo,trans,diag,'N',n,A,lda,B,scale,   &
     &                           Rwork,info)
!
!                 Check error code from DLATRS.
!
                     IF ( info/=0 )                                     &
     &                    CALL ALAERH(path,'DLATRS',info,0,uplo//trans//&
     &                    diag//'N',n,n,-1,-1,-1,imat,nfail,nerrs,Nout)
!
                     CALL DTRT03(uplo,trans,diag,n,1,A,lda,scale,Rwork, &
     &                           ONE,B,lda,X,lda,Work,result(8))
!
!+    TEST 9
!                 Solve op(A)*X = b again with NORMIN = 'Y'.
!
                     CALL DCOPY(n,X,1,B(n+1),1)
                     CALL DLATRS(uplo,trans,diag,'Y',n,A,lda,B(n+1),    &
     &                           scale,Rwork,info)
!
!                 Check error code from DLATRS.
!
                     IF ( info/=0 )                                     &
     &                    CALL ALAERH(path,'DLATRS',info,0,uplo//trans//&
     &                    diag//'Y',n,n,-1,-1,-1,imat,nfail,nerrs,Nout)
!
                     CALL DTRT03(uplo,trans,diag,n,1,A,lda,scale,Rwork, &
     &                           ONE,B(n+1),lda,X,lda,Work,result(9))
!
!                 Print information about the tests that did not pass
!                 the threshold.
!
                     IF ( result(8)>=Thresh ) THEN
                        IF ( nfail==0 .AND. nerrs==0 )                  &
     &                       CALL ALAHD(Nout,path)
                        WRITE (Nout,FMT=99004) 'DLATRS' , uplo , trans ,&
     &                         diag , 'N' , n , imat , 8 , result(8)
                        nfail = nfail + 1
                     ENDIF
                     IF ( result(9)>=Thresh ) THEN
                        IF ( nfail==0 .AND. nerrs==0 )                  &
     &                       CALL ALAHD(Nout,path)
                        WRITE (Nout,FMT=99004) 'DLATRS' , uplo , trans ,&
     &                         diag , 'Y' , n , imat , 9 , result(9)
                        nfail = nfail + 1
                     ENDIF
                     nrun = nrun + 2
                  ENDDO
               ENDDO
            ENDIF
         ENDDO
      ENDDO
!
!     Print a summary of the results.
!
      CALL ALASUM(path,Nout,nfail,nrun,nerrs)
!
99001 FORMAT (' UPLO=''',A1,''', DIAG=''',A1,''', N=',I5,', NB=',I4,    &
     &        ', type ',I2,', test(',I2,')= ',G12.5)
99002 FORMAT (' UPLO=''',A1,''', TRANS=''',A1,''', DIAG=''',A1,''', N=',&
     &        I5,', NB=',I4,', type ',I2,                               &
     &',                                                                &
     &                                                                  &
     &                                                                  &
     &    test(',I2,')= ',G12.5)
99003 FORMAT (' NORM=''',A1,''', UPLO =''',A1,''', N=',I5,',',11X,      &
     &        ' type ',I2,', test(',I2,')=',G12.5)
99004 FORMAT (1X,A,'( ''',A1,''', ''',A1,''', ''',A1,''', ''',A1,''',', &
     &        I5,', ... ), type ',I2,', test(',I2,')=',G12.5)
!
!     End of DCHKTR
!
      END SUBROUTINE DCHKTR
