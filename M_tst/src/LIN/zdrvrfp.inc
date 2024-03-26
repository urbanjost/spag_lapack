!*==zdrvrfp.f90  processed by SPAG 7.51RB at 17:35 on  4 Mar 2022
!> \brief \b zdrvrfp
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE zdrvrfp( NOUT, NN, NVAL, NNS, NSVAL, NNT, NTVAL,
!      +              THRESH, A, ASAV, AFAC, AINV, B,
!      +              BSAV, XACT, X, ARF, ARFINV,
!      +              Z_WORK_ZLATMS, Z_WORK_ZPOT02,
!      +              Z_WORK_ZPOT03, D_WORK_ZLATMS, D_WORK_ZLANHE,
!      +              D_WORK_ZPOT01, D_WORK_ZPOT02, D_WORK_ZPOT03 )
!
!       .. Scalar Arguments ..
!       INTEGER            NN, NNS, NNT, NOUT
!       DOUBLE PRECISION   THRESH
!       ..
!       .. Array Arguments ..
!       INTEGER            NVAL( NN ), NSVAL( NNS ), NTVAL( NNT )
!       COMPLEX*16         A( * )
!       COMPLEX*16         AINV( * )
!       COMPLEX*16         ASAV( * )
!       COMPLEX*16         B( * )
!       COMPLEX*16         BSAV( * )
!       COMPLEX*16         AFAC( * )
!       COMPLEX*16         ARF( * )
!       COMPLEX*16         ARFINV( * )
!       COMPLEX*16         XACT( * )
!       COMPLEX*16         X( * )
!       COMPLEX*16         Z_WORK_ZLATMS( * )
!       COMPLEX*16         Z_WORK_ZPOT02( * )
!       COMPLEX*16         Z_WORK_ZPOT03( * )
!       DOUBLE PRECISION   D_WORK_ZLATMS( * )
!       DOUBLE PRECISION   D_WORK_ZLANHE( * )
!       DOUBLE PRECISION   D_WORK_ZPOT01( * )
!       DOUBLE PRECISION   D_WORK_ZPOT02( * )
!       DOUBLE PRECISION   D_WORK_ZPOT03( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> zdrvrfp tests the LAPACK RFP routines:
!>     ZPFTRF, ZPFTRS, and ZPFTRI.
!>
!> This testing routine follow the same tests as zdrvpo (test for the full
!> format Symmetric Positive Definite solver).
!>
!> The tests are performed in Full Format, conversion back and forth from
!> full format to RFP format are performed using the routines ZTRTTF and
!> ZTFTTR.
!>
!> First, a specific matrix A of size N is created. There is nine types of
!> different matrixes possible.
!>  1. Diagonal                        6. Random, CNDNUM = sqrt(0.1/EPS)
!>  2. Random, CNDNUM = 2              7. Random, CNDNUM = 0.1/EPS
!> *3. First row and column zero       8. Scaled near underflow
!> *4. Last row and column zero        9. Scaled near overflow
!> *5. Middle row and column zero
!> (* - tests error exits from ZPFTRF, no test ratios are computed)
!> A solution XACT of size N-by-NRHS is created and the associated right
!> hand side B as well. Then ZPFTRF is called to compute L (or U), the
!> Cholesky factor of A. Then L (or U) is used to solve the linear system
!> of equations AX = B. This gives X. Then L (or U) is used to compute the
!> inverse of A, AINV. The following four tests are then performed:
!> (1) norm( L*L' - A ) / ( N * norm(A) * EPS ) or
!>     norm( U'*U - A ) / ( N * norm(A) * EPS ),
!> (2) norm(B - A*X) / ( norm(A) * norm(X) * EPS ),
!> (3) norm( I - A*AINV ) / ( N * norm(A) * norm(AINV) * EPS ),
!> (4) ( norm(X-XACT) * RCOND ) / ( norm(XACT) * EPS ),
!> where EPS is the machine precision, RCOND the condition number of A, and
!> norm( . ) the 1-norm for (1,2,3) and the inf-norm for (4).
!> Errors occur when INFO parameter is not as expected. Failures occur when
!> a test ratios is greater than THRES.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] NOUT
!> \verbatim
!>          NOUT is INTEGER
!>                The unit number for output.
!> \endverbatim
!>
!> \param[in] NN
!> \verbatim
!>          NN is INTEGER
!>                The number of values of N contained in the vector NVAL.
!> \endverbatim
!>
!> \param[in] NVAL
!> \verbatim
!>          NVAL is INTEGER array, dimension (NN)
!>                The values of the matrix dimension N.
!> \endverbatim
!>
!> \param[in] NNS
!> \verbatim
!>          NNS is INTEGER
!>                The number of values of NRHS contained in the vector NSVAL.
!> \endverbatim
!>
!> \param[in] NSVAL
!> \verbatim
!>          NSVAL is INTEGER array, dimension (NNS)
!>                The values of the number of right-hand sides NRHS.
!> \endverbatim
!>
!> \param[in] NNT
!> \verbatim
!>          NNT is INTEGER
!>                The number of values of MATRIX TYPE contained in the vector NTVAL.
!> \endverbatim
!>
!> \param[in] NTVAL
!> \verbatim
!>          NTVAL is INTEGER array, dimension (NNT)
!>                The values of matrix type (between 0 and 9 for PO/PP/PF matrices).
!> \endverbatim
!>
!> \param[in] THRESH
!> \verbatim
!>          THRESH is DOUBLE PRECISION
!>                The threshold value for the test ratios.  A result is
!>                included in the output file if RESULT >= THRESH.  To have
!>                every test ratio printed, use THRESH = 0.
!> \endverbatim
!>
!> \param[out] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension (NMAX*NMAX)
!> \endverbatim
!>
!> \param[out] ASAV
!> \verbatim
!>          ASAV is COMPLEX*16 array, dimension (NMAX*NMAX)
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
!>          B is COMPLEX*16 array, dimension (NMAX*MAXRHS)
!> \endverbatim
!>
!> \param[out] BSAV
!> \verbatim
!>          BSAV is COMPLEX*16 array, dimension (NMAX*MAXRHS)
!> \endverbatim
!>
!> \param[out] XACT
!> \verbatim
!>          XACT is COMPLEX*16 array, dimension (NMAX*MAXRHS)
!> \endverbatim
!>
!> \param[out] X
!> \verbatim
!>          X is COMPLEX*16 array, dimension (NMAX*MAXRHS)
!> \endverbatim
!>
!> \param[out] ARF
!> \verbatim
!>          ARF is COMPLEX*16 array, dimension ((NMAX*(NMAX+1))/2)
!> \endverbatim
!>
!> \param[out] ARFINV
!> \verbatim
!>          ARFINV is COMPLEX*16 array, dimension ((NMAX*(NMAX+1))/2)
!> \endverbatim
!>
!> \param[out] Z_WORK_ZLATMS
!> \verbatim
!>          Z_WORK_ZLATMS is COMPLEX*16 array, dimension ( 3*NMAX )
!> \endverbatim
!>
!> \param[out] Z_WORK_ZPOT02
!> \verbatim
!>          Z_WORK_ZPOT02 is COMPLEX*16 array, dimension ( NMAX*MAXRHS )
!> \endverbatim
!>
!> \param[out] Z_WORK_ZPOT03
!> \verbatim
!>          Z_WORK_ZPOT03 is COMPLEX*16 array, dimension ( NMAX*NMAX )
!> \endverbatim
!>
!> \param[out] D_WORK_ZLATMS
!> \verbatim
!>          D_WORK_ZLATMS is DOUBLE PRECISION array, dimension ( NMAX )
!> \endverbatim
!>
!> \param[out] D_WORK_ZLANHE
!> \verbatim
!>          D_WORK_ZLANHE is DOUBLE PRECISION array, dimension ( NMAX )
!> \endverbatim
!>
!> \param[out] D_WORK_ZPOT01
!> \verbatim
!>          D_WORK_ZPOT01 is DOUBLE PRECISION array, dimension ( NMAX )
!> \endverbatim
!>
!> \param[out] D_WORK_ZPOT02
!> \verbatim
!>          D_WORK_ZPOT02 is DOUBLE PRECISION array, dimension ( NMAX )
!> \endverbatim
!>
!> \param[out] D_WORK_ZPOT03
!> \verbatim
!>          D_WORK_ZPOT03 is DOUBLE PRECISION array, dimension ( NMAX )
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
      SUBROUTINE zdrvrfp(Nout,Nn,Nval,Nns,Nsval,Nnt,Ntval,Thresh,A,Asav,&
     &                   Afac,Ainv,B,Bsav,Xact,X,Arf,Arfinv,            &
     &                   Z_work_zlatms,Z_work_zpot02,Z_work_zpot03,     &
     &                   D_work_zlatms,D_work_zlanhe,D_work_zpot01,     &
     &                   D_work_zpot02,D_work_zpot03)
      IMPLICIT NONE
!*--ZDRVRFP247
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Nn , Nns , Nnt , Nout
      DOUBLE PRECISION Thresh
!     ..
!     .. Array Arguments ..
      INTEGER Nval(Nn) , Nsval(Nns) , Ntval(Nnt)
      COMPLEX*16 A(*)
      COMPLEX*16 Ainv(*)
      COMPLEX*16 Asav(*)
      COMPLEX*16 B(*)
      COMPLEX*16 Bsav(*)
      COMPLEX*16 Afac(*)
      COMPLEX*16 Arf(*)
      COMPLEX*16 Arfinv(*)
      COMPLEX*16 Xact(*)
      COMPLEX*16 X(*)
      COMPLEX*16 Z_work_zlatms(*)
      COMPLEX*16 Z_work_zpot02(*)
      COMPLEX*16 Z_work_zpot03(*)
      DOUBLE PRECISION D_work_zlatms(*)
      DOUBLE PRECISION D_work_zlanhe(*)
      DOUBLE PRECISION D_work_zpot01(*)
      DOUBLE PRECISION D_work_zpot02(*)
      DOUBLE PRECISION D_work_zpot03(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ONE , ZERO
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
      INTEGER NTESTS
      PARAMETER (NTESTS=4)
!     ..
!     .. Local Scalars ..
      LOGICAL zerot
      INTEGER i , info , iuplo , lda , ldb , imat , nerrs , nfail ,     &
     &        nrhs , nrun , izero , ioff , k , nt , n , iform , iin ,   &
     &        iit , iis
      CHARACTER dist , ctype , uplo , cform
      INTEGER kl , ku , mode
      DOUBLE PRECISION anorm , ainvnm , cndnum , rcondc
!     ..
!     .. Local Arrays ..
      CHARACTER uplos(2) , forms(2)
      INTEGER iseed(4) , iseedy(4)
      DOUBLE PRECISION result(NTESTS)
!     ..
!     .. External Functions ..
      DOUBLE PRECISION ZLANHE
EXTERNAL    ZLANHE
!     ..
!     .. External Subroutines ..
EXTERNAL    ZTFTTR,   ZLACPY
EXTERNAL    ZPFTRI,   ZPFTRF
EXTERNAL    ZPFTRS,   ZPOTRI,   ZPOTRF
EXTERNAL    ZTRTTF
!     ..
!     .. Scalars in Common ..
      CHARACTER*32 SRNamt
!     ..
!     .. Common blocks ..
      COMMON /SRNAMC/ SRNamt
!     ..
!     .. Data statements ..
      DATA iseedy/1988 , 1989 , 1990 , 1991/
      DATA uplos/'U' , 'L'/
      DATA forms/'N' , 'C'/
!     ..
!     .. Executable Statements ..
!
!     Initialize constants and the random number seed.
!
      nrun = 0
      nfail = 0
      nerrs = 0
      DO i = 1 , 4
         iseed(i) = iseedy(i)
      ENDDO
!
      DO iin = 1 , Nn
!
         n = Nval(iin)
         lda = MAX(n,1)
         ldb = MAX(n,1)
!
         DO iis = 1 , Nns
!
            nrhs = Nsval(iis)
!
            DO iit = 1 , Nnt
!
               imat = Ntval(iit)
!
!              If N.EQ.0, only consider the first type
!
               IF ( n/=0 .OR. iit<1 ) THEN
!
!              Skip types 3, 4, or 5 if the matrix size is too small.
!
                  IF ( imat/=4 .OR. n>1 ) THEN
                     IF ( imat/=5 .OR. n>2 ) THEN
!
!              Do first for UPLO = 'U', then for UPLO = 'L'
!
                        DO iuplo = 1 , 2
                           uplo = uplos(iuplo)
!
!                 Do first for CFORM = 'N', then for CFORM = 'C'
!
                           DO iform = 1 , 2
                              cform = forms(iform)
!
!                    Set up parameters with zlatb4 and generate a test
!                    matrix with ZLATMS.
!
                              CALL zlatb4('ZPO',imat,n,n,ctype,kl,ku,   &
     &                           anorm,mode,cndnum,dist)
!
                              SRNamt = 'ZLATMS'
                              CALL ZLATMS(n,n,dist,iseed,ctype,         &
     &                           D_work_zlatms,mode,cndnum,anorm,kl,ku, &
     &                           uplo,A,lda,Z_work_zlatms,info)
!
!                    Check error code from ZLATMS.
!
                              IF ( info/=0 ) THEN
                                 CALL alaerh('ZPF','ZLATMS',info,0,uplo,&
     &                              n,n,-1,-1,-1,iit,nfail,nerrs,Nout)
                                 CYCLE
                              ENDIF
!
!                    For types 3-5, zero one row and column of the matrix to
!                    test that INFO is returned correctly.
!
                              zerot = imat>=3 .AND. imat<=5
                              IF ( zerot ) THEN
                                 IF ( iit==3 ) THEN
                                    izero = 1
                                 ELSEIF ( iit==4 ) THEN
                                    izero = n
                                 ELSE
                                    izero = n/2 + 1
                                 ENDIF
                                 ioff = (izero-1)*lda
!
!                       Set row and column IZERO of A to 0.
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
!                    Set the imaginary part of the diagonals.
!
                              CALL zlaipd(n,A,lda+1,0)
!
!                    Save a copy of the matrix A in ASAV.
!
                              CALL ZLACPY(uplo,n,n,A,lda,Asav,lda)
!
!                    Compute the condition number of A (RCONDC).
!
                              IF ( zerot ) THEN
                                 rcondc = ZERO
                              ELSE
!
!                       Compute the 1-norm of A.
!
                                 anorm = ZLANHE('1',uplo,n,A,lda,       &
     &                              D_work_zlanhe)
!
!                       Factor the matrix A.
!
                                 CALL ZPOTRF(uplo,n,A,lda,info)
!
!                       Form the inverse of A.
!
                                 CALL ZPOTRI(uplo,n,A,lda,info)
 
                                 IF ( n/=0 ) THEN
!
!                          Compute the 1-norm condition number of A.
!
                                    ainvnm = ZLANHE('1',uplo,n,A,lda,   &
     &                                 D_work_zlanhe)
                                    rcondc = (ONE/anorm)/ainvnm
!
!                          Restore the matrix A.
!
                                    CALL ZLACPY(uplo,n,n,Asav,lda,A,lda)
                                 ENDIF
!
                              ENDIF
!
!                    Form an exact solution and set the right hand side.
!
                              SRNamt = 'zlarhs'
                              CALL zlarhs('ZPO','N',uplo,' ',n,n,kl,ku, &
     &                           nrhs,A,lda,Xact,lda,B,lda,iseed,info)
                              CALL ZLACPY('Full',n,nrhs,B,lda,Bsav,lda)
!
!                    Compute the L*L' or U'*U factorization of the
!                    matrix and solve the system.
!
                              CALL ZLACPY(uplo,n,n,A,lda,Afac,lda)
                              CALL ZLACPY('Full',n,nrhs,B,ldb,X,ldb)
!
                              SRNamt = 'ZTRTTF'
                              CALL ZTRTTF(cform,uplo,n,Afac,lda,Arf,    &
     &                           info)
                              SRNamt = 'ZPFTRF'
                              CALL ZPFTRF(cform,uplo,n,Arf,info)
!
!                    Check error code from ZPFTRF.
!
                              IF ( info/=izero ) THEN
!
!                       LANGOU: there is a small hick here: IZERO should
!                       always be INFO however if INFO is ZERO, alaerh does not
!                       complain.
!
                                 CALL alaerh('ZPF','ZPFSV ',info,izero, &
     &                              uplo,n,n,-1,-1,nrhs,iit,nfail,nerrs,&
     &                              Nout)
                                 CYCLE
                              ENDIF
!
!                     Skip the tests if INFO is not 0.
!
                              IF ( info/=0 ) CYCLE
!
                              SRNamt = 'ZPFTRS'
                              CALL ZPFTRS(cform,uplo,n,nrhs,Arf,X,ldb,  &
     &                           info)
!
                              SRNamt = 'ZTFTTR'
                              CALL ZTFTTR(cform,uplo,n,Arf,Afac,lda,    &
     &                           info)
!
!                    Reconstruct matrix from factors and compute
!                    residual.
!
                              CALL ZLACPY(uplo,n,n,Afac,lda,Asav,lda)
                              CALL zpot01(uplo,n,A,lda,Afac,lda,        &
     &                           D_work_zpot01,result(1))
                              CALL ZLACPY(uplo,n,n,Asav,lda,Afac,lda)
!
!                    Form the inverse and compute the residual.
!
                              IF ( MOD(n,2)==0 ) THEN
                                 CALL ZLACPY('A',n+1,n/2,Arf,n+1,Arfinv,&
     &                              n+1)
                              ELSE
                                 CALL ZLACPY('A',n,(n+1)/2,Arf,n,Arfinv,&
     &                              n)
                              ENDIF
!
                              SRNamt = 'ZPFTRI'
                              CALL ZPFTRI(cform,uplo,n,Arfinv,info)
!
                              SRNamt = 'ZTFTTR'
                              CALL ZTFTTR(cform,uplo,n,Arfinv,Ainv,lda, &
     &                           info)
!
!                    Check error code from ZPFTRI.
!
                              IF ( info/=0 )                            &
     &                             CALL alaerh('ZPO','ZPFTRI',info,0,   &
     &                             uplo,n,n,-1,-1,-1,imat,nfail,nerrs,  &
     &                             Nout)
!
                              CALL zpot03(uplo,n,A,lda,Ainv,lda,        &
     &                           Z_work_zpot03,lda,D_work_zpot03,rcondc,&
     &                           result(2))
!
!                    Compute residual of the computed solution.
!
                              CALL ZLACPY('Full',n,nrhs,B,lda,          &
     &                           Z_work_zpot02,lda)
                              CALL zpot02(uplo,n,nrhs,A,lda,X,lda,      &
     &                           Z_work_zpot02,lda,D_work_zpot02,       &
     &                           result(3))
!
!                    Check solution from generated exact solution.
!
                              CALL zget04(n,nrhs,X,lda,Xact,lda,rcondc, &
     &                           result(4))
                              nt = 4
!
!                    Print information about the tests that did not
!                    pass the threshold.
!
                              DO k = 1 , nt
                                 IF ( result(k)>=Thresh ) THEN
                                    IF ( nfail==0 .AND. nerrs==0 )      &
     &                                 CALL aladhd(Nout,'ZPF')
                                    WRITE (Nout,FMT=99001) 'ZPFSV ' ,   &
     &                                 uplo , n , iit , k , result(k)
                                    nfail = nfail + 1
                                 ENDIF
                              ENDDO
                              nrun = nrun + nt
                           ENDDO
                        ENDDO
                     ENDIF
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
      ENDDO
!
!     Print a summary of the results.
!
      CALL alasvm('ZPF',Nout,nfail,nrun,nerrs)
!
99001 FORMAT (1X,A6,', UPLO=''',A1,''', N =',I5,', type ',I1,', test(', &
     &        I1,')=',G12.5)
!
!
!     End of zdrvrfp
!
      END SUBROUTINE zdrvrfp
