!*==ddrvrfp.f90  processed by SPAG 7.51RB at 17:34 on  4 Mar 2022
!> \brief \b ddrvrfp
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ddrvrfp( NOUT, NN, NVAL, NNS, NSVAL, NNT, NTVAL,
!      +              THRESH, A, ASAV, AFAC, AINV, B,
!      +              BSAV, XACT, X, ARF, ARFINV,
!      +              D_WORK_DLATMS, D_WORK_DPOT01, D_TEMP_DPOT02,
!      +              D_TEMP_DPOT03, D_WORK_DLANSY,
!      +              D_WORK_DPOT02, D_WORK_DPOT03 )
!
!       .. Scalar Arguments ..
!       INTEGER            NN, NNS, NNT, NOUT
!       DOUBLE PRECISION   THRESH
!       ..
!       .. Array Arguments ..
!       INTEGER            NVAL( NN ), NSVAL( NNS ), NTVAL( NNT )
!       DOUBLE PRECISION   A( * )
!       DOUBLE PRECISION   AINV( * )
!       DOUBLE PRECISION   ASAV( * )
!       DOUBLE PRECISION   B( * )
!       DOUBLE PRECISION   BSAV( * )
!       DOUBLE PRECISION   AFAC( * )
!       DOUBLE PRECISION   ARF( * )
!       DOUBLE PRECISION   ARFINV( * )
!       DOUBLE PRECISION   XACT( * )
!       DOUBLE PRECISION   X( * )
!       DOUBLE PRECISION   D_WORK_DLATMS( * )
!       DOUBLE PRECISION   D_WORK_DPOT01( * )
!       DOUBLE PRECISION   D_TEMP_DPOT02( * )
!       DOUBLE PRECISION   D_TEMP_DPOT03( * )
!       DOUBLE PRECISION   D_WORK_DLANSY( * )
!       DOUBLE PRECISION   D_WORK_DPOT02( * )
!       DOUBLE PRECISION   D_WORK_DPOT03( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ddrvrfp tests the LAPACK RFP routines:
!>     DPFTRF, DPFTRS, and DPFTRI.
!>
!> This testing routine follow the same tests as ddrvpo (test for the full
!> format Symmetric Positive Definite solver).
!>
!> The tests are performed in Full Format, conversion back and forth from
!> full format to RFP format are performed using the routines DTRTTF and
!> DTFTTR.
!>
!> First, a specific matrix A of size N is created. There is nine types of
!> different matrixes possible.
!>  1. Diagonal                        6. Random, CNDNUM = sqrt(0.1/EPS)
!>  2. Random, CNDNUM = 2              7. Random, CNDNUM = 0.1/EPS
!> *3. First row and column zero       8. Scaled near underflow
!> *4. Last row and column zero        9. Scaled near overflow
!> *5. Middle row and column zero
!> (* - tests error exits from DPFTRF, no test ratios are computed)
!> A solution XACT of size N-by-NRHS is created and the associated right
!> hand side B as well. Then DPFTRF is called to compute L (or U), the
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
!>          A is DOUBLE PRECISION array, dimension (NMAX*NMAX)
!> \endverbatim
!>
!> \param[out] ASAV
!> \verbatim
!>          ASAV is DOUBLE PRECISION array, dimension (NMAX*NMAX)
!> \endverbatim
!>
!> \param[out] AFAC
!> \verbatim
!>          AFAC is DOUBLE PRECISION array, dimension (NMAX*NMAX)
!> \endverbatim
!>
!> \param[out] AINV
!> \verbatim
!>          AINV is DOUBLE PRECISION array, dimension (NMAX*NMAX)
!> \endverbatim
!>
!> \param[out] B
!> \verbatim
!>          B is DOUBLE PRECISION array, dimension (NMAX*MAXRHS)
!> \endverbatim
!>
!> \param[out] BSAV
!> \verbatim
!>          BSAV is DOUBLE PRECISION array, dimension (NMAX*MAXRHS)
!> \endverbatim
!>
!> \param[out] XACT
!> \verbatim
!>          XACT is DOUBLE PRECISION array, dimension (NMAX*MAXRHS)
!> \endverbatim
!>
!> \param[out] X
!> \verbatim
!>          X is DOUBLE PRECISION array, dimension (NMAX*MAXRHS)
!> \endverbatim
!>
!> \param[out] ARF
!> \verbatim
!>          ARF is DOUBLE PRECISION array, dimension ((NMAX*(NMAX+1))/2)
!> \endverbatim
!>
!> \param[out] ARFINV
!> \verbatim
!>          ARFINV is DOUBLE PRECISION array, dimension ((NMAX*(NMAX+1))/2)
!> \endverbatim
!>
!> \param[out] D_WORK_DLATMS
!> \verbatim
!>          D_WORK_DLATMS is DOUBLE PRECISION array, dimension ( 3*NMAX )
!> \endverbatim
!>
!> \param[out] D_WORK_DPOT01
!> \verbatim
!>          D_WORK_DPOT01 is DOUBLE PRECISION array, dimension ( NMAX )
!> \endverbatim
!>
!> \param[out] D_TEMP_DPOT02
!> \verbatim
!>          D_TEMP_DPOT02 is DOUBLE PRECISION array, dimension ( NMAX*MAXRHS )
!> \endverbatim
!>
!> \param[out] D_TEMP_DPOT03
!> \verbatim
!>          D_TEMP_DPOT03 is DOUBLE PRECISION array, dimension ( NMAX*NMAX )
!> \endverbatim
!>
!> \param[out] D_WORK_DLANSY
!> \verbatim
!>          D_WORK_DLANSY is DOUBLE PRECISION array, dimension ( NMAX )
!> \endverbatim
!>
!> \param[out] D_WORK_DPOT02
!> \verbatim
!>          D_WORK_DPOT02 is DOUBLE PRECISION array, dimension ( NMAX )
!> \endverbatim
!>
!> \param[out] D_WORK_DPOT03
!> \verbatim
!>          D_WORK_DPOT03 is DOUBLE PRECISION array, dimension ( NMAX )
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
      SUBROUTINE ddrvrfp(Nout,Nn,Nval,Nns,Nsval,Nnt,Ntval,Thresh,A,Asav,&
     &                   Afac,Ainv,B,Bsav,Xact,X,Arf,Arfinv,            &
     &                   D_work_dlatms,D_work_dpot01,D_temp_dpot02,     &
     &                   D_temp_dpot03,D_work_dlansy,D_work_dpot02,     &
     &                   D_work_dpot03)
      IMPLICIT NONE
!*--DDRVRFP241
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
      DOUBLE PRECISION A(*)
      DOUBLE PRECISION Ainv(*)
      DOUBLE PRECISION Asav(*)
      DOUBLE PRECISION B(*)
      DOUBLE PRECISION Bsav(*)
      DOUBLE PRECISION Afac(*)
      DOUBLE PRECISION Arf(*)
      DOUBLE PRECISION Arfinv(*)
      DOUBLE PRECISION Xact(*)
      DOUBLE PRECISION X(*)
      DOUBLE PRECISION D_work_dlatms(*)
      DOUBLE PRECISION D_work_dpot01(*)
      DOUBLE PRECISION D_temp_dpot02(*)
      DOUBLE PRECISION D_temp_dpot03(*)
      DOUBLE PRECISION D_work_dlansy(*)
      DOUBLE PRECISION D_work_dpot02(*)
      DOUBLE PRECISION D_work_dpot03(*)
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
      DOUBLE PRECISION DLANSY
EXTERNAL    DLANSY
!     ..
!     .. External Subroutines ..
EXTERNAL    DTFTTR,   DLACPY
EXTERNAL    DPFTRI,   DPFTRF,   DPFTRS
EXTERNAL    DPOTRI,   DPOTRF,   DTRTTF
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
      DATA forms/'N' , 'T'/
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
!                    Set up parameters with dlatb4 and generate a test
!                    matrix with DLATMS.
!
                              CALL dlatb4('DPO',imat,n,n,ctype,kl,ku,   &
     &                           anorm,mode,cndnum,dist)
!
                              SRNamt = 'DLATMS'
                              CALL DLATMS(n,n,dist,iseed,ctype,         &
     &                           D_work_dlatms,mode,cndnum,anorm,kl,ku, &
     &                           uplo,A,lda,D_work_dlatms,info)
!
!                    Check error code from DLATMS.
!
                              IF ( info/=0 ) THEN
                                 CALL alaerh('DPF','DLATMS',info,0,uplo,&
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
!                    Save a copy of the matrix A in ASAV.
!
                              CALL DLACPY(uplo,n,n,A,lda,Asav,lda)
!
!                    Compute the condition number of A (RCONDC).
!
                              IF ( zerot ) THEN
                                 rcondc = ZERO
                              ELSE
!
!                       Compute the 1-norm of A.
!
                                 anorm = DLANSY('1',uplo,n,A,lda,       &
     &                              D_work_dlansy)
!
!                       Factor the matrix A.
!
                                 CALL DPOTRF(uplo,n,A,lda,info)
!
!                       Form the inverse of A.
!
                                 CALL DPOTRI(uplo,n,A,lda,info)
 
                                 IF ( n/=0 ) THEN
 
!
!                          Compute the 1-norm condition number of A.
!
                                    ainvnm = DLANSY('1',uplo,n,A,lda,   &
     &                                 D_work_dlansy)
                                    rcondc = (ONE/anorm)/ainvnm
!
!                          Restore the matrix A.
!
                                    CALL DLACPY(uplo,n,n,Asav,lda,A,lda)
                                 ENDIF
!
                              ENDIF
!
!                    Form an exact solution and set the right hand side.
!
                              SRNamt = 'dlarhs'
                              CALL dlarhs('DPO','N',uplo,' ',n,n,kl,ku, &
     &                           nrhs,A,lda,Xact,lda,B,lda,iseed,info)
                              CALL DLACPY('Full',n,nrhs,B,lda,Bsav,lda)
!
!                    Compute the L*L' or U'*U factorization of the
!                    matrix and solve the system.
!
                              CALL DLACPY(uplo,n,n,A,lda,Afac,lda)
                              CALL DLACPY('Full',n,nrhs,B,ldb,X,ldb)
!
                              SRNamt = 'DTRTTF'
                              CALL DTRTTF(cform,uplo,n,Afac,lda,Arf,    &
     &                           info)
                              SRNamt = 'DPFTRF'
                              CALL DPFTRF(cform,uplo,n,Arf,info)
!
!                    Check error code from DPFTRF.
!
                              IF ( info/=izero ) THEN
!
!                       LANGOU: there is a small hick here: IZERO should
!                       always be INFO however if INFO is ZERO, alaerh does not
!                       complain.
!
                                 CALL alaerh('DPF','DPFSV ',info,izero, &
     &                              uplo,n,n,-1,-1,nrhs,iit,nfail,nerrs,&
     &                              Nout)
                                 CYCLE
                              ENDIF
!
!                    Skip the tests if INFO is not 0.
!
                              IF ( info/=0 ) CYCLE
!
                              SRNamt = 'DPFTRS'
                              CALL DPFTRS(cform,uplo,n,nrhs,Arf,X,ldb,  &
     &                           info)
!
                              SRNamt = 'DTFTTR'
                              CALL DTFTTR(cform,uplo,n,Arf,Afac,lda,    &
     &                           info)
!
!                    Reconstruct matrix from factors and compute
!                    residual.
!
                              CALL DLACPY(uplo,n,n,Afac,lda,Asav,lda)
                              CALL dpot01(uplo,n,A,lda,Afac,lda,        &
     &                           D_work_dpot01,result(1))
                              CALL DLACPY(uplo,n,n,Asav,lda,Afac,lda)
!
!                    Form the inverse and compute the residual.
!
                              IF ( MOD(n,2)==0 ) THEN
                                 CALL DLACPY('A',n+1,n/2,Arf,n+1,Arfinv,&
     &                              n+1)
                              ELSE
                                 CALL DLACPY('A',n,(n+1)/2,Arf,n,Arfinv,&
     &                              n)
                              ENDIF
!
                              SRNamt = 'DPFTRI'
                              CALL DPFTRI(cform,uplo,n,Arfinv,info)
!
                              SRNamt = 'DTFTTR'
                              CALL DTFTTR(cform,uplo,n,Arfinv,Ainv,lda, &
     &                           info)
!
!                    Check error code from DPFTRI.
!
                              IF ( info/=0 )                            &
     &                             CALL alaerh('DPO','DPFTRI',info,0,   &
     &                             uplo,n,n,-1,-1,-1,imat,nfail,nerrs,  &
     &                             Nout)
!
                              CALL dpot03(uplo,n,A,lda,Ainv,lda,        &
     &                           D_temp_dpot03,lda,D_work_dpot03,rcondc,&
     &                           result(2))
!
!                    Compute residual of the computed solution.
!
                              CALL DLACPY('Full',n,nrhs,B,lda,          &
     &                           D_temp_dpot02,lda)
                              CALL dpot02(uplo,n,nrhs,A,lda,X,lda,      &
     &                           D_temp_dpot02,lda,D_work_dpot02,       &
     &                           result(3))
!
!                    Check solution from generated exact solution.
 
                              CALL dget04(n,nrhs,X,lda,Xact,lda,rcondc, &
     &                           result(4))
                              nt = 4
!
!                    Print information about the tests that did not
!                    pass the threshold.
!
                              DO k = 1 , nt
                                 IF ( result(k)>=Thresh ) THEN
                                    IF ( nfail==0 .AND. nerrs==0 )      &
     &                                 CALL aladhd(Nout,'DPF')
                                    WRITE (Nout,FMT=99001) 'DPFSV ' ,   &
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
      CALL alasvm('DPF',Nout,nfail,nrun,nerrs)
!
99001 FORMAT (1X,A6,', UPLO=''',A1,''', N =',I5,', type ',I1,', test(', &
     &        I1,')=',G12.5)
!
!
!     End of ddrvrfp
!
      END SUBROUTINE ddrvrfp
