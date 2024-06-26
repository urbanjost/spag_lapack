!*==zchktp.f90  processed by SPAG 7.51RB at 17:35 on  4 Mar 2022
!> \brief \b zchktp
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE zchktp( DOTYPE, NN, NVAL, NNS, NSVAL, THRESH, TSTERR,
!                          NMAX, AP, AINVP, B, X, XACT, WORK, RWORK,
!                          NOUT )
!
!       .. Scalar Arguments ..
!       LOGICAL            TSTERR
!       INTEGER            NMAX, NN, NNS, NOUT
!       DOUBLE PRECISION   THRESH
!       ..
!       .. Array Arguments ..
!       LOGICAL            DOTYPE( * )
!       INTEGER            NSVAL( * ), NVAL( * )
!       DOUBLE PRECISION   RWORK( * )
!       COMPLEX*16         AINVP( * ), AP( * ), B( * ), WORK( * ), X( * ),
!      $                   XACT( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> zchktp tests ZTPTRI, -TRS, -RFS, and -CON, and ZLATPS
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
!>          The leading dimension of the work arrays.  NMAX >= the
!>          maximumm value of N in NVAL.
!> \endverbatim
!>
!> \param[out] AP
!> \verbatim
!>          AP is COMPLEX*16 array, dimension (NMAX*(NMAX+1)/2)
!> \endverbatim
!>
!> \param[out] AINVP
!> \verbatim
!>          AINVP is COMPLEX*16 array, dimension (NMAX*(NMAX+1)/2)
!> \endverbatim
!>
!> \param[out] B
!> \verbatim
!>          B is COMPLEX*16 array, dimension (NMAX*NSMAX)
!>          where NSMAX is the largest entry in NSVAL.
!> \endverbatim
!>
!> \param[out] X
!> \verbatim
!>          X is COMPLEX*16 array, dimension (NMAX*NSMAX)
!> \endverbatim
!>
!> \param[out] XACT
!> \verbatim
!>          XACT is COMPLEX*16 array, dimension (NMAX*NSMAX)
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX*16 array, dimension
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
!> \ingroup complex16_lin
!
!  =====================================================================
      SUBROUTINE zchktp(Dotype,Nn,Nval,Nns,Nsval,Thresh,Tsterr,Nmax,Ap, &
     &                  Ainvp,B,X,Xact,Work,Rwork,Nout)
      IMPLICIT NONE
!*--ZCHKTP154
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      LOGICAL Tsterr
      INTEGER Nmax , Nn , Nns , Nout
      DOUBLE PRECISION Thresh
!     ..
!     .. Array Arguments ..
      LOGICAL Dotype(*)
      INTEGER Nsval(*) , Nval(*)
      DOUBLE PRECISION Rwork(*)
      COMPLEX*16 Ainvp(*) , Ap(*) , B(*) , Work(*) , X(*) , Xact(*)
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
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
!     ..
!     .. Local Scalars ..
      CHARACTER diag , norm , trans , uplo , xtype
      CHARACTER*3 path
      INTEGER i , idiag , imat , in , info , irhs , itran , iuplo , k , &
     &        lap , lda , n , nerrs , nfail , nrhs , nrun
      DOUBLE PRECISION ainvnm , anorm , rcond , rcondc , rcondi ,       &
     &                 rcondo , scale
!     ..
!     .. Local Arrays ..
      CHARACTER transs(NTRAN) , uplos(2)
      INTEGER iseed(4) , iseedy(4)
      DOUBLE PRECISION result(NTESTS)
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      DOUBLE PRECISION ZLANTP
EXTERNAL    LSAME,   ZLANTP
!     ..
!     .. External Subroutines ..
EXTERNAL    ZCOPY
EXTERNAL    ZLACPY,   ZLATPS,   ZTPCON,   ZTPRFS
EXTERNAL    ZTPTRI
EXTERNAL    ZTPTRS
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
      path(1:1) = 'Zomplex precision'
      path(2:3) = 'TP'
      nrun = 0
      nfail = 0
      nerrs = 0
      DO i = 1 , 4
         iseed(i) = iseedy(i)
      ENDDO
!
!     Test the error exits
!
      IF ( Tsterr ) CALL zerrtr(path,Nout)
      INFot = 0
!
      DO in = 1 , Nn
!
!        Do for each value of N in NVAL
!
         n = Nval(in)
         lda = MAX(1,n)
         lap = lda*(lda+1)/2
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
!              Call zlattp to generate a triangular test matrix.
!
                  SRNamt = 'zlattp'
                  CALL zlattp(imat,uplo,'No transpose',diag,iseed,n,Ap, &
     &                        X,Work,Rwork,info)
!
!              Set IDIAG = 1 for non-unit matrices, 2 for unit.
!
                  IF ( LSAME(diag,'N') ) THEN
                     idiag = 1
                  ELSE
                     idiag = 2
                  ENDIF
!
!+    TEST 1
!              Form the inverse of A.
!
                  IF ( n>0 ) CALL ZCOPY(lap,Ap,1,Ainvp,1)
                  SRNamt = 'ZTPTRI'
                  CALL ZTPTRI(uplo,diag,n,Ainvp,info)
!
!              Check error code from ZTPTRI.
!
                  IF ( info/=0 )                                        &
     &                 CALL alaerh(path,'ZTPTRI',info,0,uplo//diag,n,n, &
     &                 -1,-1,-1,imat,nfail,nerrs,Nout)
!
!              Compute the infinity-norm condition number of A.
!
                  anorm = ZLANTP('I',uplo,diag,n,Ap,Rwork)
                  ainvnm = ZLANTP('I',uplo,diag,n,Ainvp,Rwork)
                  IF ( anorm<=ZERO .OR. ainvnm<=ZERO ) THEN
                     rcondi = ONE
                  ELSE
                     rcondi = (ONE/anorm)/ainvnm
                  ENDIF
!
!              Compute the residual for the triangular matrix times its
!              inverse.  Also compute the 1-norm condition number of A.
!
                  CALL ztpt01(uplo,diag,n,Ap,Ainvp,rcondo,Rwork,        &
     &                        result(1))
!
!              Print the test ratio if it is .GE. THRESH.
!
                  IF ( result(1)>=Thresh ) THEN
                     IF ( nfail==0 .AND. nerrs==0 )                     &
     &                    CALL alahd(Nout,path)
                     WRITE (Nout,FMT=99001) uplo , diag , n , imat , 1 ,&
     &                      result(1)
                     nfail = nfail + 1
                  ENDIF
                  nrun = nrun + 1
!
                  DO irhs = 1 , Nns
                     nrhs = Nsval(irhs)
                     xtype = 'N'
!
                     DO itran = 1 , NTRAN
!
!                 Do for op(A) = A, A**T, or A**H.
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
!                 Solve and compute residual for op(A)*x = b.
!
                        SRNamt = 'zlarhs'
                        CALL zlarhs(path,xtype,uplo,trans,n,n,0,idiag,  &
     &                              nrhs,Ap,lap,Xact,lda,B,lda,iseed,   &
     &                              info)
                        xtype = 'C'
                        CALL ZLACPY('Full',n,nrhs,B,lda,X,lda)
!
                        SRNamt = 'ZTPTRS'
                        CALL ZTPTRS(uplo,trans,diag,n,nrhs,Ap,X,lda,    &
     &                              info)
!
!                 Check error code from ZTPTRS.
!
                        IF ( info/=0 )                                  &
     &                       CALL alaerh(path,'ZTPTRS',info,0,uplo//    &
     &                       trans//diag,n,n,-1,-1,-1,imat,nfail,nerrs, &
     &                       Nout)
!
                        CALL ztpt02(uplo,trans,diag,n,nrhs,Ap,X,lda,B,  &
     &                              lda,Work,Rwork,result(2))
!
!+    TEST 3
!                 Check solution from generated exact solution.
!
                        CALL zget04(n,nrhs,X,lda,Xact,lda,rcondc,       &
     &                              result(3))
!
!+    TESTS 4, 5, and 6
!                 Use iterative refinement to improve the solution and
!                 compute error bounds.
!
                        SRNamt = 'ZTPRFS'
                        CALL ZTPRFS(uplo,trans,diag,n,nrhs,Ap,B,lda,X,  &
     &                              lda,Rwork,Rwork(nrhs+1),Work,       &
     &                              Rwork(2*nrhs+1),info)
!
!                 Check error code from ZTPRFS.
!
                        IF ( info/=0 )                                  &
     &                       CALL alaerh(path,'ZTPRFS',info,0,uplo//    &
     &                       trans//diag,n,n,-1,-1,nrhs,imat,nfail,     &
     &                       nerrs,Nout)
!
                        CALL zget04(n,nrhs,X,lda,Xact,lda,rcondc,       &
     &                              result(4))
                        CALL ztpt05(uplo,trans,diag,n,nrhs,Ap,B,lda,X,  &
     &                              lda,Xact,lda,Rwork,Rwork(nrhs+1),   &
     &                              result(5))
!
!                    Print information about the tests that did not pass
!                    the threshold.
!
                        DO k = 2 , 6
                           IF ( result(k)>=Thresh ) THEN
                              IF ( nfail==0 .AND. nerrs==0 )            &
     &                             CALL alahd(Nout,path)
                              WRITE (Nout,FMT=99002) uplo , trans ,     &
     &                               diag , n , nrhs , imat , k ,       &
     &                               result(k)
                              nfail = nfail + 1
                           ENDIF
                        ENDDO
                        nrun = nrun + 5
                     ENDDO
                  ENDDO
!
!+    TEST 7
!                 Get an estimate of RCOND = 1/CNDNUM.
!
                  DO itran = 1 , 2
                     IF ( itran==1 ) THEN
                        norm = 'O'
                        rcondc = rcondo
                     ELSE
                        norm = 'I'
                        rcondc = rcondi
                     ENDIF
                     SRNamt = 'ZTPCON'
                     CALL ZTPCON(norm,uplo,diag,n,Ap,rcond,Work,Rwork,  &
     &                           info)
!
!                 Check error code from ZTPCON.
!
                     IF ( info/=0 )                                     &
     &                    CALL alaerh(path,'ZTPCON',info,0,norm//uplo// &
     &                    diag,n,n,-1,-1,-1,imat,nfail,nerrs,Nout)
!
                     CALL ztpt06(rcond,rcondc,uplo,diag,n,Ap,Rwork,     &
     &                           result(7))
!
!                 Print the test ratio if it is .GE. THRESH.
!
                     IF ( result(7)>=Thresh ) THEN
                        IF ( nfail==0 .AND. nerrs==0 )                  &
     &                       CALL alahd(Nout,path)
                        WRITE (Nout,FMT=99003) 'ZTPCON' , norm , uplo , &
     &                         diag , n , imat , 7 , result(7)
                        nfail = nfail + 1
                     ENDIF
                     nrun = nrun + 1
                  ENDDO
               ENDDO
            ENDIF
         ENDDO
!
!        Use pathological test matrices to test ZLATPS.
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
!                 Do for op(A) = A, A**T, or A**H.
!
                     trans = transs(itran)
!
!                 Call zlattp to generate a triangular test matrix.
!
                     SRNamt = 'zlattp'
                     CALL zlattp(imat,uplo,trans,diag,iseed,n,Ap,X,Work,&
     &                           Rwork,info)
!
!+    TEST 8
!                 Solve the system op(A)*x = b.
!
                     SRNamt = 'ZLATPS'
                     CALL ZCOPY(n,X,1,B,1)
                     CALL ZLATPS(uplo,trans,diag,'N',n,Ap,B,scale,Rwork,&
     &                           info)
!
!                 Check error code from ZLATPS.
!
                     IF ( info/=0 )                                     &
     &                    CALL alaerh(path,'ZLATPS',info,0,uplo//trans//&
     &                    diag//'N',n,n,-1,-1,-1,imat,nfail,nerrs,Nout)
!
                     CALL ztpt03(uplo,trans,diag,n,1,Ap,scale,Rwork,ONE,&
     &                           B,lda,X,lda,Work,result(8))
!
!+    TEST 9
!                 Solve op(A)*x = b again with NORMIN = 'Y'.
!
                     CALL ZCOPY(n,X,1,B(n+1),1)
                     CALL ZLATPS(uplo,trans,diag,'Y',n,Ap,B(n+1),scale, &
     &                           Rwork,info)
!
!                 Check error code from ZLATPS.
!
                     IF ( info/=0 )                                     &
     &                    CALL alaerh(path,'ZLATPS',info,0,uplo//trans//&
     &                    diag//'Y',n,n,-1,-1,-1,imat,nfail,nerrs,Nout)
!
                     CALL ztpt03(uplo,trans,diag,n,1,Ap,scale,Rwork,ONE,&
     &                           B(n+1),lda,X,lda,Work,result(9))
!
!                 Print information about the tests that did not pass
!                 the threshold.
!
                     IF ( result(8)>=Thresh ) THEN
                        IF ( nfail==0 .AND. nerrs==0 )                  &
     &                       CALL alahd(Nout,path)
                        WRITE (Nout,FMT=99004) 'ZLATPS' , uplo , trans ,&
     &                         diag , 'N' , n , imat , 8 , result(8)
                        nfail = nfail + 1
                     ENDIF
                     IF ( result(9)>=Thresh ) THEN
                        IF ( nfail==0 .AND. nerrs==0 )                  &
     &                       CALL alahd(Nout,path)
                        WRITE (Nout,FMT=99004) 'ZLATPS' , uplo , trans ,&
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
      CALL alasum(path,Nout,nfail,nrun,nerrs)
!
99001 FORMAT (' UPLO=''',A1,''', DIAG=''',A1,''', N=',I5,', type ',I2,  &
     &        ', test(',I2,')= ',G12.5)
99002 FORMAT (' UPLO=''',A1,''', TRANS=''',A1,''', DIAG=''',A1,''', N=',&
     &        I5,''', NRHS=',I5,', type ',I2,', test(',I2,')= ',G12.5)
99003 FORMAT (1X,A,'( ''',A1,''', ''',A1,''', ''',A1,''',',I5,          &
     &        ', ... ), type ',I2,', test(',I2,')=',G12.5)
99004 FORMAT (1X,A,'( ''',A1,''', ''',A1,''', ''',A1,''', ''',A1,''',', &
     &        I5,', ... ), type ',I2,', test(',I2,')=',G12.5)
!
!     End of zchktp
!
      END SUBROUTINE zchktp
