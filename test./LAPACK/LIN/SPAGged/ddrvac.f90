!*==ddrvac.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b DDRVAC
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DDRVAC( DOTYPE, NM, MVAL, NNS, NSVAL, THRESH, NMAX,
!                          A, AFAC, B, X, WORK,
!                          RWORK, SWORK, NOUT )
!
!       .. Scalar Arguments ..
!       INTEGER            NMAX, NM, NNS, NOUT
!       DOUBLE PRECISION   THRESH
!       ..
!       .. Array Arguments ..
!       LOGICAL            DOTYPE( * )
!       INTEGER            MVAL( * ), NSVAL( * )
!       REAL               SWORK(*)
!       DOUBLE PRECISION   A( * ), AFAC( * ), B( * ),
!      $                   RWORK( * ), WORK( * ), X( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DDRVAC tests DSPOSV.
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
!>          The number of values of N contained in the vector MVAL.
!> \endverbatim
!>
!> \param[in] MVAL
!> \verbatim
!>          MVAL is INTEGER array, dimension (NM)
!>          The values of the matrix dimension N.
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
!>          The maximum value permitted for N, used in dimensioning the
!>          work arrays.
!> \endverbatim
!>
!> \param[out] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (NMAX*NMAX)
!> \endverbatim
!>
!> \param[out] AFAC
!> \verbatim
!>          AFAC is DOUBLE PRECISION array, dimension (NMAX*NMAX)
!> \endverbatim
!>
!> \param[out] B
!> \verbatim
!>          B is DOUBLE PRECISION array, dimension (NMAX*NSMAX)
!> \endverbatim
!>
!> \param[out] X
!> \verbatim
!>          X is DOUBLE PRECISION array, dimension (NMAX*NSMAX)
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
!>                      (max(2*NMAX,2*NSMAX+NWORK))
!> \endverbatim
!>
!> \param[out] SWORK
!> \verbatim
!>          SWORK is REAL array, dimension
!>                      (NMAX*(NSMAX+NMAX))
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
      SUBROUTINE DDRVAC(Dotype,Nm,Mval,Nns,Nsval,Thresh,Nmax,A,Afac,B,X,&
     &                  Work,Rwork,Swork,Nout)
      IMPLICIT NONE
!*--DDRVAC147
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Nmax , Nm , Nns , Nout
      DOUBLE PRECISION Thresh
!     ..
!     .. Array Arguments ..
      LOGICAL Dotype(*)
      INTEGER Mval(*) , Nsval(*)
      REAL Swork(*)
      DOUBLE PRECISION A(*) , Afac(*) , B(*) , Rwork(*) , Work(*) , X(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D+0)
      INTEGER NTYPES
      PARAMETER (NTYPES=9)
      INTEGER NTESTS
      PARAMETER (NTESTS=1)
!     ..
!     .. Local Scalars ..
      LOGICAL zerot
      CHARACTER dist , type , uplo , xtype
      CHARACTER*3 path
      INTEGER i , im , imat , info , ioff , irhs , iuplo , izero , kl , &
     &        ku , lda , mode , n , nerrs , nfail , nimat , nrhs , nrun
      DOUBLE PRECISION anorm , cndnum
!     ..
!     .. Local Arrays ..
      CHARACTER uplos(2)
      INTEGER iseed(4) , iseedy(4)
      DOUBLE PRECISION result(NTESTS)
!     ..
!     .. Local Variables ..
      INTEGER iter , kase
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL ALAERH , DLACPY , DLARHS , DLASET , DLATB4 , DLATMS ,    &
     &         DPOT06 , DSPOSV
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC DBLE , MAX , SQRT
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
      DATA uplos/'U' , 'L'/
!     ..
!     .. Executable Statements ..
!
!     Initialize constants and the random number seed.
!
      kase = 0
      path(1:1) = 'Double precision'
      path(2:3) = 'PO'
      nrun = 0
      nfail = 0
      nerrs = 0
      DO i = 1 , 4
         iseed(i) = iseedy(i)
      ENDDO
!
      INFot = 0
!
!     Do for each value of N in MVAL
!
      DO im = 1 , Nm
         n = Mval(im)
         lda = MAX(n,1)
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
!              Set up parameters with DLATB4 and generate a test matrix
!              with DLATMS.
!
                     CALL DLATB4(path,imat,n,n,type,kl,ku,anorm,mode,   &
     &                           cndnum,dist)
!
                     SRNamt = 'DLATMS'
                     CALL DLATMS(n,n,dist,iseed,type,Rwork,mode,cndnum, &
     &                           anorm,kl,ku,uplo,A,lda,Work,info)
!
!              Check error code from DLATMS.
!
                     IF ( info/=0 ) THEN
                        CALL ALAERH(path,'DLATMS',info,0,uplo,n,n,-1,-1,&
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
                     DO irhs = 1 , Nns
                        nrhs = Nsval(irhs)
                        xtype = 'N'
!
!                 Form an exact solution and set the right hand side.
!
                        SRNamt = 'DLARHS'
                        CALL DLARHS(path,xtype,uplo,' ',n,n,kl,ku,nrhs, &
     &                              A,lda,X,lda,B,lda,iseed,info)
!
!                 Compute the L*L' or U'*U factorization of the
!                 matrix and solve the system.
!
                        SRNamt = 'DSPOSV '
                        kase = kase + 1
!
                        CALL DLACPY('All',n,n,A,lda,Afac,lda)
!
                        CALL DSPOSV(uplo,n,nrhs,Afac,lda,B,lda,X,lda,   &
     &                              Work,Swork,iter,info)
 
                        IF ( iter<0 )                                   &
     &                       CALL DLACPY('All',n,n,A,lda,Afac,lda)
!
!                 Check error code from DSPOSV .
!
                        IF ( info/=izero ) THEN
!
                           IF ( nfail==0 .AND. nerrs==0 )               &
     &                          CALL ALAHD(Nout,path)
                           nerrs = nerrs + 1
!
                           IF ( info/=izero .AND. izero/=0 ) THEN
                              WRITE (Nout,FMT=99005) 'DSPOSV' , info ,  &
     &                               izero , n , imat
                           ELSE
                              WRITE (Nout,FMT=99006) 'DSPOSV' , info ,  &
     &                               n , imat
                           ENDIF
                        ENDIF
!
!                 Skip the remaining test if the matrix is singular.
!
                        IF ( info/=0 ) GOTO 50
!
!                 Check the quality of the solution
!
                        CALL DLACPY('All',n,nrhs,B,lda,Work,lda)
!
                        CALL DPOT06(uplo,n,nrhs,A,lda,X,lda,Work,lda,   &
     &                              Rwork,result(1))
!
!                 Check if the test passes the tesing.
!                 Print information about the tests that did not
!                 pass the testing.
!
!                 If iterative refinement has been used and claimed to
!                 be successful (ITER>0), we want
!                 NORM1(B - A*X)/(NORM1(A)*NORM1(X)*EPS*SRQT(N)) < 1
!
!                 If double precision has been used (ITER<0), we want
!                 NORM1(B - A*X)/(NORM1(A)*NORM1(X)*EPS) < THRES
!                 (Cf. the linear solver testing routines)
!
                        IF ( (Thresh<=0.0E+00) .OR.                     &
     &                       ((iter>=0) .AND. (n>0) .AND.               &
     &                       (result(1)>=SQRT(DBLE(n)))) .OR.           &
     &                       ((iter<0) .AND. (result(1)>=Thresh)) ) THEN
!
                           IF ( nfail==0 .AND. nerrs==0 ) THEN
                              WRITE (Nout,FMT=99007) 'DPO'
                              WRITE (Nout,FMT='( '' Matrix types:'' )')
                              WRITE (Nout,FMT=99008)
                              WRITE (Nout,FMT='( '' Test ratios:'' )')
                              WRITE (Nout,FMT=99009) 1
                              WRITE (Nout,FMT='( '' Messages:'' )')
                           ENDIF
!
                           WRITE (Nout,FMT=99001) uplo , n , nrhs ,     &
     &                            imat , 1 , result(1)
!
                           nfail = nfail + 1
!
                        ENDIF
!
                        nrun = nrun + 1
!
                     ENDDO
                  ENDDO
               ENDIF
            ENDIF
 50      ENDDO
      ENDDO
!
!     Print a summary of the results.
!
      IF ( nfail>0 ) THEN
         WRITE (Nout,FMT=99002) 'DSPOSV' , nfail , nrun
      ELSE
         WRITE (Nout,FMT=99003) 'DSPOSV' , nrun
      ENDIF
      IF ( nerrs>0 ) WRITE (Nout,FMT=99004) nerrs
!
99001 FORMAT (' UPLO=''',A1,''', N =',I5,', NRHS=',I3,', type ',I2,     &
     &        ', test(',I2,') =',G12.5)
99002 FORMAT (1X,A6,': ',I6,' out of ',I6,                              &
     &        ' tests failed to pass the threshold')
99003 FORMAT (/1X,'All tests for ',A6,                                  &
     &        ' routines passed the threshold ( ',I6,' tests run)')
99004 FORMAT (6X,I6,' error messages recorded')
!
!     SUBNAM, INFO, INFOE, N, IMAT
!
99005 FORMAT (' *** ',A6,' returned with INFO =',I5,' instead of ',I5,  &
     &        /' ==> N =',I5,', type ',I2)
!
!     SUBNAM, INFO, N, IMAT
!
99006 FORMAT (' *** Error code from ',A6,'=',I5,' for M=',I5,', type ', &
     &        I2)
99007 FORMAT (/1X,A3,':  positive definite dense matrices')
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
     &        '( norm_1(A) * norm_1(X) * EPS ) > THRES if DPOTRF')
 
!
!     End of DDRVAC
!
      END SUBROUTINE DDRVAC
