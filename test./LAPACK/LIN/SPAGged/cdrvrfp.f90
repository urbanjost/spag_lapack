!*==cdrvrfp.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b CDRVRFP
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE CDRVRFP( NOUT, NN, NVAL, NNS, NSVAL, NNT, NTVAL,
!      +              THRESH, A, ASAV, AFAC, AINV, B,
!      +              BSAV, XACT, X, ARF, ARFINV,
!      +              C_WORK_CLATMS, C_WORK_CPOT02,
!      +              C_WORK_CPOT03, S_WORK_CLATMS, S_WORK_CLANHE,
!      +              S_WORK_CPOT01, S_WORK_CPOT02, S_WORK_CPOT03 )
!
!       .. Scalar Arguments ..
!       INTEGER            NN, NNS, NNT, NOUT
!       REAL               THRESH
!       ..
!       .. Array Arguments ..
!       INTEGER            NVAL( NN ), NSVAL( NNS ), NTVAL( NNT )
!       COMPLEX            A( * )
!       COMPLEX            AINV( * )
!       COMPLEX            ASAV( * )
!       COMPLEX            B( * )
!       COMPLEX            BSAV( * )
!       COMPLEX            AFAC( * )
!       COMPLEX            ARF( * )
!       COMPLEX            ARFINV( * )
!       COMPLEX            XACT( * )
!       COMPLEX            X( * )
!       COMPLEX            C_WORK_CLATMS( * )
!       COMPLEX            C_WORK_CPOT02( * )
!       COMPLEX            C_WORK_CPOT03( * )
!       REAL               S_WORK_CLATMS( * )
!       REAL               S_WORK_CLANHE( * )
!       REAL               S_WORK_CPOT01( * )
!       REAL               S_WORK_CPOT02( * )
!       REAL               S_WORK_CPOT03( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CDRVRFP tests the LAPACK RFP routines:
!>     CPFTRF, CPFTRS, and CPFTRI.
!>
!> This testing routine follow the same tests as CDRVPO (test for the full
!> format Symmetric Positive Definite solver).
!>
!> The tests are performed in Full Format, conversion back and forth from
!> full format to RFP format are performed using the routines CTRTTF and
!> CTFTTR.
!>
!> First, a specific matrix A of size N is created. There is nine types of
!> different matrixes possible.
!>  1. Diagonal                        6. Random, CNDNUM = sqrt(0.1/EPS)
!>  2. Random, CNDNUM = 2              7. Random, CNDNUM = 0.1/EPS
!> *3. First row and column zero       8. Scaled near underflow
!> *4. Last row and column zero        9. Scaled near overflow
!> *5. Middle row and column zero
!> (* - tests error exits from CPFTRF, no test ratios are computed)
!> A solution XACT of size N-by-NRHS is created and the associated right
!> hand side B as well. Then CPFTRF is called to compute L (or U), the
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
!>          THRESH is REAL
!>                The threshold value for the test ratios.  A result is
!>                included in the output file if RESULT >= THRESH.  To have
!>                every test ratio printed, use THRESH = 0.
!> \endverbatim
!>
!> \param[out] A
!> \verbatim
!>          A is COMPLEX array, dimension (NMAX*NMAX)
!> \endverbatim
!>
!> \param[out] ASAV
!> \verbatim
!>          ASAV is COMPLEX array, dimension (NMAX*NMAX)
!> \endverbatim
!>
!> \param[out] AFAC
!> \verbatim
!>          AFAC is COMPLEX array, dimension (NMAX*NMAX)
!> \endverbatim
!>
!> \param[out] AINV
!> \verbatim
!>          AINV is COMPLEX array, dimension (NMAX*NMAX)
!> \endverbatim
!>
!> \param[out] B
!> \verbatim
!>          B is COMPLEX array, dimension (NMAX*MAXRHS)
!> \endverbatim
!>
!> \param[out] BSAV
!> \verbatim
!>          BSAV is COMPLEX array, dimension (NMAX*MAXRHS)
!> \endverbatim
!>
!> \param[out] XACT
!> \verbatim
!>          XACT is COMPLEX array, dimension (NMAX*MAXRHS)
!> \endverbatim
!>
!> \param[out] X
!> \verbatim
!>          X is COMPLEX array, dimension (NMAX*MAXRHS)
!> \endverbatim
!>
!> \param[out] ARF
!> \verbatim
!>          ARF is COMPLEX array, dimension ((NMAX*(NMAX+1))/2)
!> \endverbatim
!>
!> \param[out] ARFINV
!> \verbatim
!>          ARFINV is COMPLEX array, dimension ((NMAX*(NMAX+1))/2)
!> \endverbatim
!>
!> \param[out] C_WORK_CLATMS
!> \verbatim
!>          C_WORK_CLATMS is COMPLEX array, dimension ( 3*NMAX )
!> \endverbatim
!>
!> \param[out] C_WORK_CPOT02
!> \verbatim
!>          C_WORK_CPOT02 is COMPLEX array, dimension ( NMAX*MAXRHS )
!> \endverbatim
!>
!> \param[out] C_WORK_CPOT03
!> \verbatim
!>          C_WORK_CPOT03 is COMPLEX array, dimension ( NMAX*NMAX )
!> \endverbatim
!>
!> \param[out] S_WORK_CLATMS
!> \verbatim
!>          S_WORK_CLATMS is REAL array, dimension ( NMAX )
!> \endverbatim
!>
!> \param[out] S_WORK_CLANHE
!> \verbatim
!>          S_WORK_CLANHE is REAL array, dimension ( NMAX )
!> \endverbatim
!>
!> \param[out] S_WORK_CPOT01
!> \verbatim
!>          S_WORK_CPOT01 is REAL array, dimension ( NMAX )
!> \endverbatim
!>
!> \param[out] S_WORK_CPOT02
!> \verbatim
!>          S_WORK_CPOT02 is REAL array, dimension ( NMAX )
!> \endverbatim
!>
!> \param[out] S_WORK_CPOT03
!> \verbatim
!>          S_WORK_CPOT03 is REAL array, dimension ( NMAX )
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
      SUBROUTINE CDRVRFP(Nout,Nn,Nval,Nns,Nsval,Nnt,Ntval,Thresh,A,Asav,&
     &                   Afac,Ainv,B,Bsav,Xact,X,Arf,Arfinv,            &
     &                   C_work_clatms,C_work_cpot02,C_work_cpot03,     &
     &                   S_work_clatms,S_work_clanhe,S_work_cpot01,     &
     &                   S_work_cpot02,S_work_cpot03)
      IMPLICIT NONE
!*--CDRVRFP247
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Nn , Nns , Nnt , Nout
      REAL Thresh
!     ..
!     .. Array Arguments ..
      INTEGER Nval(Nn) , Nsval(Nns) , Ntval(Nnt)
      COMPLEX A(*)
      COMPLEX Ainv(*)
      COMPLEX Asav(*)
      COMPLEX B(*)
      COMPLEX Bsav(*)
      COMPLEX Afac(*)
      COMPLEX Arf(*)
      COMPLEX Arfinv(*)
      COMPLEX Xact(*)
      COMPLEX X(*)
      COMPLEX C_work_clatms(*)
      COMPLEX C_work_cpot02(*)
      COMPLEX C_work_cpot03(*)
      REAL S_work_clatms(*)
      REAL S_work_clanhe(*)
      REAL S_work_cpot01(*)
      REAL S_work_cpot02(*)
      REAL S_work_cpot03(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ONE , ZERO
      PARAMETER (ONE=1.0E+0,ZERO=0.0E+0)
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
      REAL anorm , ainvnm , cndnum , rcondc
!     ..
!     .. Local Arrays ..
      CHARACTER uplos(2) , forms(2)
      INTEGER iseed(4) , iseedy(4)
      REAL result(NTESTS)
!     ..
!     .. External Functions ..
      REAL CLANHE
      EXTERNAL CLANHE
!     ..
!     .. External Subroutines ..
      EXTERNAL ALADHD , ALAERH , ALASVM , CGET04 , CTFTTR , CLACPY ,    &
     &         CLAIPD , CLARHS , CLATB4 , CLATMS , CPFTRI , CPFTRF ,    &
     &         CPFTRS , CPOT01 , CPOT02 , CPOT03 , CPOTRI , CPOTRF ,    &
     &         CTRTTF
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
!                    Set up parameters with CLATB4 and generate a test
!                    matrix with CLATMS.
!
                              CALL CLATB4('CPO',imat,n,n,ctype,kl,ku,   &
     &                           anorm,mode,cndnum,dist)
!
                              SRNamt = 'CLATMS'
                              CALL CLATMS(n,n,dist,iseed,ctype,         &
     &                           S_work_clatms,mode,cndnum,anorm,kl,ku, &
     &                           uplo,A,lda,C_work_clatms,info)
!
!                    Check error code from CLATMS.
!
                              IF ( info/=0 ) THEN
                                 CALL ALAERH('CPF','CLATMS',info,0,uplo,&
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
                              CALL CLAIPD(n,A,lda+1,0)
!
!                    Save a copy of the matrix A in ASAV.
!
                              CALL CLACPY(uplo,n,n,A,lda,Asav,lda)
!
!                    Compute the condition number of A (RCONDC).
!
                              IF ( zerot ) THEN
                                 rcondc = ZERO
                              ELSE
!
!                       Compute the 1-norm of A.
!
                                 anorm = CLANHE('1',uplo,n,A,lda,       &
     &                              S_work_clanhe)
!
!                       Factor the matrix A.
!
                                 CALL CPOTRF(uplo,n,A,lda,info)
!
!                       Form the inverse of A.
!
                                 CALL CPOTRI(uplo,n,A,lda,info)
!
!                       Compute the 1-norm condition number of A.
!
                                 IF ( n/=0 ) THEN
                                    ainvnm = CLANHE('1',uplo,n,A,lda,   &
     &                                 S_work_clanhe)
                                    rcondc = (ONE/anorm)/ainvnm
!
!                          Restore the matrix A.
!
                                    CALL CLACPY(uplo,n,n,Asav,lda,A,lda)
                                 ENDIF
 
!
                              ENDIF
!
!                    Form an exact solution and set the right hand side.
!
                              SRNamt = 'CLARHS'
                              CALL CLARHS('CPO','N',uplo,' ',n,n,kl,ku, &
     &                           nrhs,A,lda,Xact,lda,B,lda,iseed,info)
                              CALL CLACPY('Full',n,nrhs,B,lda,Bsav,lda)
!
!                    Compute the L*L' or U'*U factorization of the
!                    matrix and solve the system.
!
                              CALL CLACPY(uplo,n,n,A,lda,Afac,lda)
                              CALL CLACPY('Full',n,nrhs,B,ldb,X,ldb)
!
                              SRNamt = 'CTRTTF'
                              CALL CTRTTF(cform,uplo,n,Afac,lda,Arf,    &
     &                           info)
                              SRNamt = 'CPFTRF'
                              CALL CPFTRF(cform,uplo,n,Arf,info)
!
!                    Check error code from CPFTRF.
!
                              IF ( info/=izero ) THEN
!
!                       LANGOU: there is a small hick here: IZERO should
!                       always be INFO however if INFO is ZERO, ALAERH does not
!                       complain.
!
                                 CALL ALAERH('CPF','CPFSV ',info,izero, &
     &                              uplo,n,n,-1,-1,nrhs,iit,nfail,nerrs,&
     &                              Nout)
                                 CYCLE
                              ENDIF
!
!                     Skip the tests if INFO is not 0.
!
                              IF ( info/=0 ) CYCLE
!
                              SRNamt = 'CPFTRS'
                              CALL CPFTRS(cform,uplo,n,nrhs,Arf,X,ldb,  &
     &                           info)
!
                              SRNamt = 'CTFTTR'
                              CALL CTFTTR(cform,uplo,n,Arf,Afac,lda,    &
     &                           info)
!
!                    Reconstruct matrix from factors and compute
!                    residual.
!
                              CALL CLACPY(uplo,n,n,Afac,lda,Asav,lda)
                              CALL CPOT01(uplo,n,A,lda,Afac,lda,        &
     &                           S_work_cpot01,result(1))
                              CALL CLACPY(uplo,n,n,Asav,lda,Afac,lda)
!
!                    Form the inverse and compute the residual.
!
                              IF ( MOD(n,2)==0 ) THEN
                                 CALL CLACPY('A',n+1,n/2,Arf,n+1,Arfinv,&
     &                              n+1)
                              ELSE
                                 CALL CLACPY('A',n,(n+1)/2,Arf,n,Arfinv,&
     &                              n)
                              ENDIF
!
                              SRNamt = 'CPFTRI'
                              CALL CPFTRI(cform,uplo,n,Arfinv,info)
!
                              SRNamt = 'CTFTTR'
                              CALL CTFTTR(cform,uplo,n,Arfinv,Ainv,lda, &
     &                           info)
!
!                    Check error code from CPFTRI.
!
                              IF ( info/=0 )                            &
     &                             CALL ALAERH('CPO','CPFTRI',info,0,   &
     &                             uplo,n,n,-1,-1,-1,imat,nfail,nerrs,  &
     &                             Nout)
!
                              CALL CPOT03(uplo,n,A,lda,Ainv,lda,        &
     &                           C_work_cpot03,lda,S_work_cpot03,rcondc,&
     &                           result(2))
!
!                    Compute residual of the computed solution.
!
                              CALL CLACPY('Full',n,nrhs,B,lda,          &
     &                           C_work_cpot02,lda)
                              CALL CPOT02(uplo,n,nrhs,A,lda,X,lda,      &
     &                           C_work_cpot02,lda,S_work_cpot02,       &
     &                           result(3))
!
!                    Check solution from generated exact solution.
!
                              CALL CGET04(n,nrhs,X,lda,Xact,lda,rcondc, &
     &                           result(4))
                              nt = 4
!
!                    Print information about the tests that did not
!                    pass the threshold.
!
                              DO k = 1 , nt
                                 IF ( result(k)>=Thresh ) THEN
                                    IF ( nfail==0 .AND. nerrs==0 )      &
     &                                 CALL ALADHD(Nout,'CPF')
                                    WRITE (Nout,FMT=99001) 'CPFSV ' ,   &
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
      CALL ALASVM('CPF',Nout,nfail,nrun,nerrs)
!
99001 FORMAT (1X,A6,', UPLO=''',A1,''', N =',I5,', type ',I1,', test(', &
     &        I1,')=',G12.5)
!
!
!     End of CDRVRFP
!
      END SUBROUTINE CDRVRFP
