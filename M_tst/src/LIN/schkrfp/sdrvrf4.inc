!*==sdrvrf4.f90  processed by SPAG 7.51RB at 17:34 on  4 Mar 2022
!> \brief \b sdrvrf4
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE sdrvrf4( NOUT, NN, NVAL, THRESH, C1, C2, LDC, CRF, A,
!      +                    LDA, S_WORK_SLANGE )
!
!       .. Scalar Arguments ..
!       INTEGER            LDA, LDC, NN, NOUT
!       REAL               THRESH
!       ..
!       .. Array Arguments ..
!       INTEGER            NVAL( NN )
!       REAL               A( LDA, * ), C1( LDC, * ), C2( LDC, *),
!      +                   CRF( * ), S_WORK_SLANGE( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> sdrvrf4 tests the LAPACK RFP routines:
!>     SSFRK
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
!> \param[in] THRESH
!> \verbatim
!>          THRESH is REAL
!>                The threshold value for the test ratios.  A result is
!>                included in the output file if RESULT >= THRESH.  To
!>                have every test ratio printed, use THRESH = 0.
!> \endverbatim
!>
!> \param[out] C1
!> \verbatim
!>          C1 is REAL array,
!>                dimension (LDC,NMAX)
!> \endverbatim
!>
!> \param[out] C2
!> \verbatim
!>          C2 is REAL array,
!>                dimension (LDC,NMAX)
!> \endverbatim
!>
!> \param[in] LDC
!> \verbatim
!>          LDC is INTEGER
!>                The leading dimension of the array A.
!>                LDA >= max(1,NMAX).
!> \endverbatim
!>
!> \param[out] CRF
!> \verbatim
!>          CRF is REAL array,
!>                dimension ((NMAX*(NMAX+1))/2).
!> \endverbatim
!>
!> \param[out] A
!> \verbatim
!>          A is REAL array,
!>                dimension (LDA,NMAX)
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>                The leading dimension of the array A.  LDA >= max(1,NMAX).
!> \endverbatim
!>
!> \param[out] S_WORK_SLANGE
!> \verbatim
!>          S_WORK_SLANGE is REAL array, dimension (NMAX)
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
      SUBROUTINE sdrvrf4(Nout,Nn,Nval,Thresh,C1,C2,Ldc,Crf,A,Lda,       &
     &                   S_work_slange)
      IMPLICIT NONE
!*--SDRVRF4122
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Lda , Ldc , Nn , Nout
      REAL Thresh
!     ..
!     .. Array Arguments ..
      INTEGER Nval(Nn)
      REAL A(Lda,*) , C1(Ldc,*) , C2(Ldc,*) , Crf(*) , S_work_slange(*)
!     ..
!
!  =====================================================================
!     ..
!     .. Parameters ..
      REAL ZERO , ONE
      PARAMETER (ZERO=0.0E+0,ONE=1.0E+0)
      INTEGER NTESTS
      PARAMETER (NTESTS=1)
!     ..
!     .. Local Scalars ..
      CHARACTER uplo , cform , trans
      INTEGER i , iform , iik , iin , info , iuplo , j , k , n , nfail ,&
     &        nrun , ialpha , itrans
      REAL alpha , beta , eps , norma , normc
!     ..
!     .. Local Arrays ..
      CHARACTER uplos(2) , forms(2) , transs(2)
      INTEGER iseed(4) , iseedy(4)
      REAL result(NTESTS)
!     ..
!     .. External Functions ..
      REAL SLAMCH , SLANGE
      EXTERNAL SLAMCH , SLANGE
!     ..
!     .. External Subroutines ..
      EXTERNAL SSYRK , SSFRK , STFTTR , STRTTF
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , MAX
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
      DATA transs/'N' , 'T'/
!     ..
!     .. Executable Statements ..
!
!     Initialize constants and the random number seed.
!
      nrun = 0
      nfail = 0
      info = 0
      DO i = 1 , 4
         iseed(i) = iseedy(i)
      ENDDO
      eps = SLAMCH('Precision')
!
      DO iin = 1 , Nn
!
         n = Nval(iin)
!
         DO iik = 1 , Nn
!
            k = Nval(iin)
!
            DO iform = 1 , 2
!
               cform = forms(iform)
!
               DO iuplo = 1 , 2
!
                  uplo = uplos(iuplo)
!
                  DO itrans = 1 , 2
!
                     trans = transs(itrans)
!
                     DO ialpha = 1 , 4
!
                        IF ( ialpha==1 ) THEN
                           alpha = ZERO
                           beta = ZERO
                        ELSEIF ( ialpha==2 ) THEN
                           alpha = ONE
                           beta = ZERO
                        ELSEIF ( ialpha==3 ) THEN
                           alpha = ZERO
                           beta = ONE
                        ELSE
                           alpha = SLARND(2,iseed)
                           beta = SLARND(2,iseed)
                        ENDIF
!
!                       All the parameters are set:
!                          CFORM, UPLO, TRANS, M, N,
!                          ALPHA, and BETA
!                       READY TO TEST!
!
                        nrun = nrun + 1
!
                        IF ( itrans==1 ) THEN
!
!                          In this case we are NOTRANS, so A is N-by-K
!
                           DO j = 1 , k
                              DO i = 1 , n
                                 A(i,j) = SLARND(2,iseed)
                              ENDDO
                           ENDDO
!
                           norma = SLANGE('I',n,k,A,Lda,S_work_slange)
!
 
                        ELSE
!
!                          In this case we are TRANS, so A is K-by-N
!
                           DO j = 1 , n
                              DO i = 1 , k
                                 A(i,j) = SLARND(2,iseed)
                              ENDDO
                           ENDDO
!
                           norma = SLANGE('I',k,n,A,Lda,S_work_slange)
!
                        ENDIF
!
!                       Generate C1 our N--by--N symmetric matrix.
!                       Make sure C2 has the same upper/lower part,
!                       (the one that we do not touch), so
!                       copy the initial C1 in C2 in it.
!
                        DO j = 1 , n
                           DO i = 1 , n
                              C1(i,j) = SLARND(2,iseed)
                              C2(i,j) = C1(i,j)
                           ENDDO
                        ENDDO
!
!                       (See comment later on for why we use SLANGE and
!                       not SLANSY for C1.)
!
                        normc = SLANGE('I',n,n,C1,Ldc,S_work_slange)
!
                        SRNamt = 'STRTTF'
                        CALL STRTTF(cform,uplo,n,C1,Ldc,Crf,info)
!
!                       call ssyrk the BLAS routine -> gives C1
!
                        SRNamt = 'SSYRK '
                        CALL SSYRK(uplo,trans,n,k,alpha,A,Lda,beta,C1,  &
     &                             Ldc)
!
!                       call ssfrk the RFP routine -> gives CRF
!
                        SRNamt = 'SSFRK '
                        CALL SSFRK(cform,uplo,trans,n,k,alpha,A,Lda,    &
     &                             beta,Crf)
!
!                       convert CRF in full format -> gives C2
!
                        SRNamt = 'STFTTR'
                        CALL STFTTR(cform,uplo,n,Crf,C2,Ldc,info)
!
!                       compare C1 and C2
!
                        DO j = 1 , n
                           DO i = 1 , n
                              C1(i,j) = C1(i,j) - C2(i,j)
                           ENDDO
                        ENDDO
!
!                       Yes, C1 is symmetric so we could call SLANSY,
!                       but we want to check the upper part that is
!                       supposed to be unchanged and the diagonal that
!                       is supposed to be real -> SLANGE
!
                        result(1) = SLANGE('I',n,n,C1,Ldc,S_work_slange)
                        result(1) = result(1)                           &
     &                              /MAX(ABS(alpha)*norma+ABS(beta),ONE)&
     &                              /MAX(n,1)/eps
!
                        IF ( result(1)>=Thresh ) THEN
                           IF ( nfail==0 ) THEN
                              WRITE (Nout,*)
                              WRITE (Nout,FMT=99001)
                           ENDIF
                           WRITE (Nout,FMT=99002) 'SSFRK' , cform ,     &
     &                            uplo , trans , n , k , result(1)
                           nfail = nfail + 1
                        ENDIF
!
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
!
!     Print a summary of the results.
!
      IF ( nfail==0 ) THEN
         WRITE (Nout,FMT=99003) 'SSFRK' , nrun
      ELSE
         WRITE (Nout,FMT=99004) 'SSFRK' , nfail , nrun
      ENDIF
!
99001 FORMAT (1X,                                                       &
     &' *** Error(s) or Failure(s) while testing SSFRK                  &
     &                                                                  &
     &                                                                  &
     &                                           ***')
99002 FORMAT (1X,'     Failure in ',A5,', CFORM=''',A1,''',',' UPLO=''',&
     &        A1,''',',' TRANS=''',A1,''',',' N=',I3,', K =',I3,        &
     &        ', test=',G12.5)
99003 FORMAT (1X,'All tests for ',A5,' auxiliary routine passed the ',  &
     &        'threshold ( ',I5,' tests run)')
99004 FORMAT (1X,A6,' auxiliary routine: ',I5,' out of ',I5,            &
     &        ' tests failed to pass the threshold')
!
!
!     End of sdrvrf4
!
      END SUBROUTINE sdrvrf4
