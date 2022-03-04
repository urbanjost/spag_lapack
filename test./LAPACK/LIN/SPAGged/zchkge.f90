!*==zchkge.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b ZCHKGE
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZCHKGE( DOTYPE, NM, MVAL, NN, NVAL, NNB, NBVAL, NNS,
!                          NSVAL, THRESH, TSTERR, NMAX, A, AFAC, AINV, B,
!                          X, XACT, WORK, RWORK, IWORK, NOUT )
!
!       .. Scalar Arguments ..
!       LOGICAL            TSTERR
!       INTEGER            NM, NMAX, NN, NNB, NNS, NOUT
!       DOUBLE PRECISION   THRESH
!       ..
!       .. Array Arguments ..
!       LOGICAL            DOTYPE( * )
!       INTEGER            IWORK( * ), MVAL( * ), NBVAL( * ), NSVAL( * ),
!      $                   NVAL( * )
!       DOUBLE PRECISION   RWORK( * )
!       COMPLEX*16         A( * ), AFAC( * ), AINV( * ), B( * ),
!      $                   WORK( * ), X( * ), XACT( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZCHKGE tests ZGETRF, -TRI, -TRS, -RFS, and -CON.
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
!> \param[in] NM
!> \verbatim
!>          NM is INTEGER
!>          The number of values of M contained in the vector MVAL.
!> \endverbatim
!>
!> \param[in] MVAL
!> \verbatim
!>          MVAL is INTEGER array, dimension (NM)
!>          The values of the matrix row dimension M.
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
!>          The maximum value permitted for M or N, used in dimensioning
!>          the work arrays.
!> \endverbatim
!>
!> \param[out] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension (NMAX*NMAX)
!> \endverbatim
!>
!> \param[out] AFAC
!> \verbatim
!>          AFAC is COMPLEX*16 array, dimension (NMAX*NMAX)
!> \endverbatim
!>
!> \param[out] AINV
!> \verbatim
!>          AINV is COMPLEX*16 array, dimension (NMAX*NMAX)
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
!>                      (max(2*NMAX,2*NSMAX+NWORK))
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
!> \ingroup complex16_lin
!
!  =====================================================================
      SUBROUTINE ZCHKGE(Dotype,Nm,Mval,Nn,Nval,Nnb,Nbval,Nns,Nsval,     &
     &                  Thresh,Tsterr,Nmax,A,Afac,Ainv,B,X,Xact,Work,   &
     &                  Rwork,Iwork,Nout)
      IMPLICIT NONE
!*--ZCHKGE190
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      LOGICAL Tsterr
      INTEGER Nm , Nmax , Nn , Nnb , Nns , Nout
      DOUBLE PRECISION Thresh
!     ..
!     .. Array Arguments ..
      LOGICAL Dotype(*)
      INTEGER Iwork(*) , Mval(*) , Nbval(*) , Nsval(*) , Nval(*)
      DOUBLE PRECISION Rwork(*)
      COMPLEX*16 A(*) , Afac(*) , Ainv(*) , B(*) , Work(*) , X(*) ,     &
     &           Xact(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ONE , ZERO
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
      INTEGER NTYPES
      PARAMETER (NTYPES=11)
      INTEGER NTESTS
      PARAMETER (NTESTS=8)
      INTEGER NTRAN
      PARAMETER (NTRAN=3)
!     ..
!     .. Local Scalars ..
      LOGICAL trfcon , zerot
      CHARACTER dist , norm , trans , type , xtype
      CHARACTER*3 path
      INTEGER i , im , imat , in , inb , info , ioff , irhs , itran ,   &
     &        izero , k , kl , ku , lda , lwork , m , mode , n , nb ,   &
     &        nerrs , nfail , nimat , nrhs , nrun , nt
      DOUBLE PRECISION ainvnm , anorm , anormi , anormo , cndnum ,      &
     &                 dummy , rcond , rcondc , rcondi , rcondo
!     ..
!     .. Local Arrays ..
      CHARACTER transs(NTRAN)
      INTEGER iseed(4) , iseedy(4)
      DOUBLE PRECISION result(NTESTS)
!     ..
!     .. External Functions ..
      DOUBLE PRECISION DGET06 , ZLANGE
      EXTERNAL DGET06 , ZLANGE
!     ..
!     .. External Subroutines ..
      EXTERNAL ALAERH , ALAHD , ALASUM , XLAENV , ZERRGE , ZGECON ,     &
     &         ZGERFS , ZGET01 , ZGET02 , ZGET03 , ZGET04 , ZGET07 ,    &
     &         ZGETRF , ZGETRI , ZGETRS , ZLACPY , ZLARHS , ZLASET ,    &
     &         ZLATB4 , ZLATMS
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC DCMPLX , MAX , MIN
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
      DATA iseedy/1988 , 1989 , 1990 , 1991/ , transs/'N' , 'T' , 'C'/
!     ..
!     .. Executable Statements ..
!
!     Initialize constants and the random number seed.
!
      path(1:1) = 'Zomplex precision'
      path(2:3) = 'GE'
      nrun = 0
      nfail = 0
      nerrs = 0
      DO i = 1 , 4
         iseed(i) = iseedy(i)
      ENDDO
!
!     Test the error exits
!
      CALL XLAENV(1,1)
      IF ( Tsterr ) CALL ZERRGE(path,Nout)
      INFot = 0
      CALL XLAENV(2,2)
!
!     Do for each value of M in MVAL
!
      DO im = 1 , Nm
         m = Mval(im)
         lda = MAX(1,m)
!
!        Do for each value of N in NVAL
!
         DO in = 1 , Nn
            n = Nval(in)
            xtype = 'N'
            nimat = NTYPES
            IF ( m<=0 .OR. n<=0 ) nimat = 1
!
            DO imat = 1 , nimat
!
!              Do the tests only if DOTYPE( IMAT ) is true.
!
               IF ( Dotype(imat) ) THEN
!
!              Skip types 5, 6, or 7 if the matrix size is too small.
!
                  zerot = imat>=5 .AND. imat<=7
                  IF ( .NOT.(zerot .AND. n<imat-4) ) THEN
!
!              Set up parameters with ZLATB4 and generate a test matrix
!              with ZLATMS.
!
                     CALL ZLATB4(path,imat,m,n,type,kl,ku,anorm,mode,   &
     &                           cndnum,dist)
!
                     SRNamt = 'ZLATMS'
                     CALL ZLATMS(m,n,dist,iseed,type,Rwork,mode,cndnum, &
     &                           anorm,kl,ku,'No packing',A,lda,Work,   &
     &                           info)
!
!              Check error code from ZLATMS.
!
                     IF ( info/=0 ) THEN
                        CALL ALAERH(path,'ZLATMS',info,0,' ',m,n,-1,-1, &
     &                              -1,imat,nfail,nerrs,Nout)
                        CYCLE
                     ENDIF
!
!              For types 5-7, zero one or more columns of the matrix to
!              test that INFO is returned correctly.
!
                     IF ( zerot ) THEN
                        IF ( imat==5 ) THEN
                           izero = 1
                        ELSEIF ( imat==6 ) THEN
                           izero = MIN(m,n)
                        ELSE
                           izero = MIN(m,n)/2 + 1
                        ENDIF
                        ioff = (izero-1)*lda
                        IF ( imat<7 ) THEN
                           DO i = 1 , m
                              A(ioff+i) = ZERO
                           ENDDO
                        ELSE
                           CALL ZLASET('Full',m,n-izero+1,DCMPLX(ZERO), &
     &                                 DCMPLX(ZERO),A(ioff+1),lda)
                        ENDIF
                     ELSE
                        izero = 0
                     ENDIF
!
!              These lines, if used in place of the calls in the DO 60
!              loop, cause the code to bomb on a Sun SPARCstation.
!
!               ANORMO = ZLANGE( 'O', M, N, A, LDA, RWORK )
!               ANORMI = ZLANGE( 'I', M, N, A, LDA, RWORK )
!
!              Do for each blocksize in NBVAL
!
                     DO inb = 1 , Nnb
                        nb = Nbval(inb)
                        CALL XLAENV(1,nb)
!
!                 Compute the LU factorization of the matrix.
!
                        CALL ZLACPY('Full',m,n,A,lda,Afac,lda)
                        SRNamt = 'ZGETRF'
                        CALL ZGETRF(m,n,Afac,lda,Iwork,info)
!
!                 Check error code from ZGETRF.
!
                        IF ( info/=izero )                              &
     &                       CALL ALAERH(path,'ZGETRF',info,izero,' ',m,&
     &                       n,-1,-1,nb,imat,nfail,nerrs,Nout)
                        trfcon = .FALSE.
!
!+    TEST 1
!                 Reconstruct matrix from factors and compute residual.
!
                        CALL ZLACPY('Full',m,n,Afac,lda,Ainv,lda)
                        CALL ZGET01(m,n,A,lda,Ainv,lda,Iwork,Rwork,     &
     &                              result(1))
                        nt = 1
!
!+    TEST 2
!                 Form the inverse if the factorization was successful
!                 and compute the residual.
!
                        IF ( m==n .AND. info==0 ) THEN
                           CALL ZLACPY('Full',n,n,Afac,lda,Ainv,lda)
                           SRNamt = 'ZGETRI'
                           nrhs = Nsval(1)
                           lwork = Nmax*MAX(3,nrhs)
                           CALL ZGETRI(n,Ainv,lda,Iwork,Work,lwork,info)
!
!                    Check error code from ZGETRI.
!
                           IF ( info/=0 )                               &
     &                          CALL ALAERH(path,'ZGETRI',info,0,' ',n, &
     &                          n,-1,-1,nb,imat,nfail,nerrs,Nout)
!
!                    Compute the residual for the matrix times its
!                    inverse.  Also compute the 1-norm condition number
!                    of A.
!
                           CALL ZGET03(n,A,lda,Ainv,lda,Work,lda,Rwork, &
     &                                 rcondo,result(2))
                           anormo = ZLANGE('O',m,n,A,lda,Rwork)
!
!                    Compute the infinity-norm condition number of A.
!
                           anormi = ZLANGE('I',m,n,A,lda,Rwork)
                           ainvnm = ZLANGE('I',n,n,Ainv,lda,Rwork)
                           IF ( anormi<=ZERO .OR. ainvnm<=ZERO ) THEN
                              rcondi = ONE
                           ELSE
                              rcondi = (ONE/anormi)/ainvnm
                           ENDIF
                           nt = 2
                        ELSE
!
!                    Do only the condition estimate if INFO > 0.
!
                           trfcon = .TRUE.
                           anormo = ZLANGE('O',m,n,A,lda,Rwork)
                           anormi = ZLANGE('I',m,n,A,lda,Rwork)
                           rcondo = ZERO
                           rcondi = ZERO
                        ENDIF
!
!                 Print information about the tests so far that did not
!                 pass the threshold.
!
                        DO k = 1 , nt
                           IF ( result(k)>=Thresh ) THEN
                              IF ( nfail==0 .AND. nerrs==0 )            &
     &                             CALL ALAHD(Nout,path)
                              WRITE (Nout,FMT=99001) m , n , nb , imat ,&
     &                               k , result(k)
                              nfail = nfail + 1
                           ENDIF
                        ENDDO
                        nrun = nrun + nt
!
!                 Skip the remaining tests if this is not the first
!                 block size or if M .ne. N.  Skip the solve tests if
!                 the matrix is singular.
!
                        IF ( inb<=1 .AND. m==n ) THEN
                           IF ( .NOT.(trfcon) ) THEN
!
                              DO irhs = 1 , Nns
                                 nrhs = Nsval(irhs)
                                 xtype = 'N'
!
                                 DO itran = 1 , NTRAN
                                    trans = transs(itran)
                                    IF ( itran==1 ) THEN
                                       rcondc = rcondo
                                    ELSE
                                       rcondc = rcondi
                                    ENDIF
!
!+    TEST 3
!                       Solve and compute residual for A * X = B.
!
                                    SRNamt = 'ZLARHS'
                                    CALL ZLARHS(path,xtype,' ',trans,n, &
     &                                 n,kl,ku,nrhs,A,lda,Xact,lda,B,   &
     &                                 lda,iseed,info)
                                    xtype = 'C'
!
                                    CALL ZLACPY('Full',n,nrhs,B,lda,X,  &
     &                                 lda)
                                    SRNamt = 'ZGETRS'
                                    CALL ZGETRS(trans,n,nrhs,Afac,lda,  &
     &                                 Iwork,X,lda,info)
!
!                       Check error code from ZGETRS.
!
                                    IF ( info/=0 )                      &
     &                                  CALL ALAERH(path,'ZGETRS',info, &
     &                                 0,trans,n,n,-1,-1,nrhs,imat,     &
     &                                 nfail,nerrs,Nout)
!
                                    CALL ZLACPY('Full',n,nrhs,B,lda,    &
     &                                 Work,lda)
                                    CALL ZGET02(trans,n,n,nrhs,A,lda,X, &
     &                                 lda,Work,lda,Rwork,result(3))
!
!+    TEST 4
!                       Check solution from generated exact solution.
!
                                    CALL ZGET04(n,nrhs,X,lda,Xact,lda,  &
     &                                 rcondc,result(4))
!
!+    TESTS 5, 6, and 7
!                       Use iterative refinement to improve the
!                       solution.
!
                                    SRNamt = 'ZGERFS'
                                    CALL ZGERFS(trans,n,nrhs,A,lda,Afac,&
     &                                 lda,Iwork,B,lda,X,lda,Rwork,     &
     &                                 Rwork(nrhs+1),Work,              &
     &                                 Rwork(2*nrhs+1),info)
!
!                       Check error code from ZGERFS.
!
                                    IF ( info/=0 )                      &
     &                                  CALL ALAERH(path,'ZGERFS',info, &
     &                                 0,trans,n,n,-1,-1,nrhs,imat,     &
     &                                 nfail,nerrs,Nout)
!
                                    CALL ZGET04(n,nrhs,X,lda,Xact,lda,  &
     &                                 rcondc,result(5))
                                    CALL ZGET07(trans,n,nrhs,A,lda,B,   &
     &                                 lda,X,lda,Xact,lda,Rwork,.TRUE., &
     &                                 Rwork(nrhs+1),result(6))
!
!                       Print information about the tests that did not
!                       pass the threshold.
!
                                    DO k = 3 , 7
                                       IF ( result(k)>=Thresh ) THEN
                                         IF ( nfail==0 .AND. nerrs==0 ) &
     &                                      CALL ALAHD(Nout,path)
                                         WRITE (Nout,FMT=99002) trans , &
     &                                      n , nrhs , imat , k ,       &
     &                                      result(k)
                                         nfail = nfail + 1
                                       ENDIF
                                    ENDDO
                                    nrun = nrun + 5
                                 ENDDO
                              ENDDO
                           ENDIF
!
!+    TEST 8
!                    Get an estimate of RCOND = 1/CNDNUM.
!
                           DO itran = 1 , 2
                              IF ( itran==1 ) THEN
                                 anorm = anormo
                                 rcondc = rcondo
                                 norm = 'O'
                              ELSE
                                 anorm = anormi
                                 rcondc = rcondi
                                 norm = 'I'
                              ENDIF
                              SRNamt = 'ZGECON'
                              CALL ZGECON(norm,n,Afac,lda,anorm,rcond,  &
     &                           Work,Rwork,info)
!
!                       Check error code from ZGECON.
!
                              IF ( info/=0 )                            &
     &                             CALL ALAERH(path,'ZGECON',info,0,    &
     &                             norm,n,n,-1,-1,-1,imat,nfail,nerrs,  &
     &                             Nout)
!
!                       This line is needed on a Sun SPARCstation.
!
                              dummy = rcond
!
                              result(8) = DGET06(rcond,rcondc)
!
!                    Print information about the tests that did not pass
!                    the threshold.
!
                              IF ( result(8)>=Thresh ) THEN
                                 IF ( nfail==0 .AND. nerrs==0 )         &
     &                                CALL ALAHD(Nout,path)
                                 WRITE (Nout,FMT=99003) norm , n ,      &
     &                                  imat , 8 , result(8)
                                 nfail = nfail + 1
                              ENDIF
                              nrun = nrun + 1
                           ENDDO
                        ENDIF
                     ENDDO
                  ENDIF
               ENDIF
            ENDDO
!
         ENDDO
      ENDDO
!
!     Print a summary of the results.
!
      CALL ALASUM(path,Nout,nfail,nrun,nerrs)
!
99001 FORMAT (' M = ',I5,', N =',I5,', NB =',I4,', type ',I2,', test(', &
     &        I2,') =',G12.5)
99002 FORMAT (' TRANS=''',A1,''', N =',I5,', NRHS=',I3,', type ',I2,    &
     &        ', test(',I2,') =',G12.5)
99003 FORMAT (' NORM =''',A1,''', N =',I5,',',10X,' type ',I2,', test(',&
     &        I2,') =',G12.5)
!
!     End of ZCHKGE
!
      END SUBROUTINE ZCHKGE
