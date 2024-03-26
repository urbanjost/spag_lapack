!*==zdrvab.f90  processed by SPAG 7.51RB at 17:35 on  4 Mar 2022
!> \brief \b zdrvab
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE zdrvab( DOTYPE, NM, MVAL, NNS,
!                          NSVAL, THRESH, NMAX, A, AFAC, B,
!                          X, WORK, RWORK, SWORK, IWORK, NOUT )
!
!       .. Scalar Arguments ..
!       INTEGER            NM, NMAX, NNS, NOUT
!       DOUBLE PRECISION   THRESH
!       ..
!       .. Array Arguments ..
!       LOGICAL            DOTYPE( * )
!       INTEGER            MVAL( * ), NSVAL( * ), IWORK( * )
!       DOUBLE PRECISION   RWORK( * )
!       COMPLEX            SWORK( * )
!       COMPLEX*16         A( * ), AFAC( * ), B( * ),
!      $                   WORK( * ), X( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> zdrvab tests ZCGESV
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
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX*16 array, dimension
!>                      (NMAX*max(3,NSMAX*2))
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is DOUBLE PRECISION array, dimension
!>                      NMAX
!> \endverbatim
!>
!> \param[out] SWORK
!> \verbatim
!>          SWORK is COMPLEX array, dimension
!>                      (NMAX*(NSMAX+NMAX))
!> \endverbatim
!>
!> \param[out] IWORK
!> \verbatim
!>          IWORK is INTEGER array, dimension
!>                      NMAX
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
      SUBROUTINE zdrvab(Dotype,Nm,Mval,Nns,Nsval,Thresh,Nmax,A,Afac,B,X,&
     &                  Work,Rwork,Swork,Iwork,Nout)
      IMPLICIT NONE
!*--ZDRVAB155
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Nm , Nmax , Nns , Nout
      DOUBLE PRECISION Thresh
!     ..
!     .. Array Arguments ..
      LOGICAL Dotype(*)
      INTEGER Mval(*) , Nsval(*) , Iwork(*)
      DOUBLE PRECISION Rwork(*)
      COMPLEX Swork(*)
      COMPLEX*16 A(*) , Afac(*) , B(*) , Work(*) , X(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D+0)
      INTEGER NTYPES
      PARAMETER (NTYPES=11)
      INTEGER NTESTS
      PARAMETER (NTESTS=1)
!     ..
!     .. Local Scalars ..
      LOGICAL zerot
      CHARACTER dist , trans , type , xtype
      CHARACTER*3 path
      INTEGER i , im , imat , info , ioff , irhs , izero , kl , ku ,    &
     &        lda , m , mode , n , nerrs , nfail , nimat , nrhs , nrun
      DOUBLE PRECISION anorm , cndnum
!     ..
!     .. Local Arrays ..
      INTEGER iseed(4) , iseedy(4)
      DOUBLE PRECISION result(NTESTS)
!     ..
!     .. Local Variables ..
      INTEGER iter , kase
!     ..
!     .. External Subroutines ..
EXTERNAL    ZLACPY,   ZLASET
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC DCMPLX , DBLE , MAX , MIN , SQRT
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
      DATA iseedy/2006 , 2007 , 2008 , 2009/
!     ..
!     .. Executable Statements ..
!
!     Initialize constants and the random number seed.
!
      kase = 0
      path(1:1) = 'Zomplex precision'
      path(2:3) = 'GE'
      nrun = 0
      nfail = 0
      nerrs = 0
      DO i = 1 , 4
         iseed(i) = iseedy(i)
      ENDDO
!
      INFot = 0
!
!     Do for each value of M in MVAL
!
      DO im = 1 , Nm
         m = Mval(im)
         lda = MAX(1,m)
!
         n = m
         nimat = NTYPES
         IF ( m<=0 .OR. n<=0 ) nimat = 1
!
         DO imat = 1 , nimat
!
!           Do the tests only if DOTYPE( IMAT ) is true.
!
            IF ( Dotype(imat) ) THEN
!
!           Skip types 5, 6, or 7 if the matrix size is too small.
!
               zerot = imat>=5 .AND. imat<=7
               IF ( .NOT.(zerot .AND. n<imat-4) ) THEN
!
!           Set up parameters with zlatb4 and generate a test matrix
!           with ZLATMS.
!
                  CALL zlatb4(path,imat,m,n,type,kl,ku,anorm,mode,      &
     &                        cndnum,dist)
!
                  SRNamt = 'ZLATMS'
                  CALL ZLATMS(m,n,dist,iseed,type,Rwork,mode,cndnum,    &
     &                        anorm,kl,ku,'No packing',A,lda,Work,info)
!
!           Check error code from ZLATMS.
!
                  IF ( info/=0 ) THEN
                     CALL alaerh(path,'ZLATMS',info,0,' ',m,n,-1,-1,-1, &
     &                           imat,nfail,nerrs,Nout)
                     CYCLE
                  ENDIF
!
!           For types 5-7, zero one or more columns of the matrix to
!           test that INFO is returned correctly.
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
                        CALL ZLASET('Full',m,n-izero+1,DCMPLX(ZERO),    &
     &                              DCMPLX(ZERO),A(ioff+1),lda)
                     ENDIF
                  ELSE
                     izero = 0
                  ENDIF
!
                  DO irhs = 1 , Nns
                     nrhs = Nsval(irhs)
                     xtype = 'N'
                     trans = 'N'
!
                     SRNamt = 'zlarhs'
                     CALL zlarhs(path,xtype,' ',trans,n,n,kl,ku,nrhs,A, &
     &                           lda,X,lda,B,lda,iseed,info)
!
                     SRNamt = 'ZCGESV'
!
                     kase = kase + 1
!
                     CALL ZLACPY('Full',m,n,A,lda,Afac,lda)
!
                     CALL ZCGESV(n,nrhs,A,lda,Iwork,B,lda,X,lda,Work,   &
     &                           Swork,Rwork,iter,info)
!
                     IF ( iter<0 )                                      &
     &                    CALL ZLACPY('Full',m,n,Afac,lda,A,lda)
!
!              Check error code from ZCGESV. This should be the same as
!              the one of DGETRF.
!
                     IF ( info/=izero ) THEN
!
                        IF ( nfail==0 .AND. nerrs==0 )                  &
     &                       CALL alahd(Nout,path)
                        nerrs = nerrs + 1
!
                        IF ( info/=izero .AND. izero/=0 ) THEN
                           WRITE (Nout,FMT=99005) 'ZCGESV' , info ,     &
     &                            izero , m , imat
                        ELSE
                           WRITE (Nout,FMT=99006) 'ZCGESV' , info , m , &
     &                            imat
                        ENDIF
                     ENDIF
!
!              Skip the remaining test if the matrix is singular.
!
                     IF ( info/=0 ) EXIT
!
!              Check the quality of the solution
!
                     CALL ZLACPY('Full',n,nrhs,B,lda,Work,lda)
!
                     CALL zget08(trans,n,n,nrhs,A,lda,X,lda,Work,lda,   &
     &                           Rwork,result(1))
!
!              Check if the test passes the tesing.
!              Print information about the tests that did not
!              pass the testing.
!
!              If iterative refinement has been used and claimed to
!              be successful (ITER>0), we want
!                NORMI(B - A*X)/(NORMI(A)*NORMI(X)*EPS*SRQT(N)) < 1
!
!              If double precision has been used (ITER<0), we want
!                NORMI(B - A*X)/(NORMI(A)*NORMI(X)*EPS) < THRES
!              (Cf. the linear solver testing routines)
!
                     IF ( (Thresh<=0.0E+00) .OR.                        &
     &                    ((iter>=0) .AND. (n>0) .AND.                  &
     &                    (result(1)>=SQRT(DBLE(n)))) .OR.              &
     &                    ((iter<0) .AND. (result(1)>=Thresh)) ) THEN
!
                        IF ( nfail==0 .AND. nerrs==0 ) THEN
                           WRITE (Nout,FMT=99007) 'DGE'
                           WRITE (Nout,FMT='( '' Matrix types:'' )')
                           WRITE (Nout,FMT=99008)
                           WRITE (Nout,FMT='( '' Test ratios:'' )')
                           WRITE (Nout,FMT=99009) 1
                           WRITE (Nout,FMT='( '' Messages:'' )')
                        ENDIF
!
                        WRITE (Nout,FMT=99001) trans , n , nrhs , imat ,&
     &                         1 , result(1)
                        nfail = nfail + 1
                     ENDIF
                     nrun = nrun + 1
                  ENDDO
               ENDIF
            ENDIF
         ENDDO
      ENDDO
!
!     Print a summary of the results.
!
      IF ( nfail>0 ) THEN
         WRITE (Nout,FMT=99002) 'ZCGESV' , nfail , nrun
      ELSE
         WRITE (Nout,FMT=99003) 'ZCGESV' , nrun
      ENDIF
      IF ( nerrs>0 ) WRITE (Nout,FMT=99004) nerrs
!
99001 FORMAT (' TRANS=''',A1,''', N =',I5,', NRHS=',I3,', type ',I2,    &
     &        ', test(',I2,') =',G12.5)
99002 FORMAT (1X,A6,': ',I6,' out of ',I6,                              &
     &        ' tests failed to pass the threshold')
99003 FORMAT (/1X,'All tests for ',A6,                                  &
     &        ' routines passed the threshold ( ',I6,' tests run)')
99004 FORMAT (6X,I6,' error messages recorded')
!
!     SUBNAM, INFO, INFOE, M, IMAT
!
99005 FORMAT (' *** ',A6,' returned with INFO =',I5,' instead of ',I5,  &
     &        /' ==> M =',I5,', type ',I2)
!
!     SUBNAM, INFO, M, IMAT
!
99006 FORMAT (' *** Error code from ',A6,'=',I5,' for M=',I5,', type ', &
     &        I2)
99007 FORMAT (/1X,A3,':  General dense matrices')
99008 FORMAT (4X,'1. Diagonal',24X,'7. Last n/2 columns zero',/4X,      &
     &        '2. Upper triangular',16X,                                &
     &        '8. Random, CNDNUM = sqrt(0.1/EPS)',/4X,                  &
     &        '3. Lower triangular',16X,'9. Random, CNDNUM = 0.1/EPS',  &
     &        /4X,'4. Random, CNDNUM = 2',13X,                          &
     &        '10. Scaled near underflow',/4X,'5. First column zero',   &
     &        14X,'11. Scaled near overflow',/4X,'6. Last column zero')
99009 FORMAT (3X,I2,': norm_1( B - A * X )  / ',                        &
     &        '( norm_1(A) * norm_1(X) * EPS * SQRT(N) ) > 1 if ITERREF'&
     &        ,/4x,'or norm_1( B - A * X )  / ',                        &
     &        '( norm_1(A) * norm_1(X) * EPS ) > THRES if DGETRF')
!
!     End of zdrvab
!
      END SUBROUTINE zdrvab
