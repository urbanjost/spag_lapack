!*==cdrvrf4.f90  processed by SPAG 7.51RB at 17:34 on  4 Mar 2022
!> \brief \b cdrvrf4
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE cdrvrf4( NOUT, NN, NVAL, THRESH, C1, C2, LDC, CRF, A,
!      +                    LDA, S_WORK_CLANGE )
!
!       .. Scalar Arguments ..
!       INTEGER            LDA, LDC, NN, NOUT
!       REAL               THRESH
!       ..
!       .. Array Arguments ..
!       INTEGER            NVAL( NN )
!       REAL               S_WORK_CLANGE( * )
!       COMPLEX            A( LDA, * ), C1( LDC, * ), C2( LDC, *),
!      +                   CRF( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> cdrvrf4 tests the LAPACK RFP routines:
!>     CHFRK
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
!>                included in the output file if RESULT >= THRESH.  To have
!>                every test ratio printed, use THRESH = 0.
!> \endverbatim
!>
!> \param[out] C1
!> \verbatim
!>          C1 is COMPLEX array, dimension (LDC,NMAX)
!> \endverbatim
!>
!> \param[out] C2
!> \verbatim
!>          C2 is COMPLEX array, dimension (LDC,NMAX)
!> \endverbatim
!>
!> \param[in] LDC
!> \verbatim
!>          LDC is INTEGER
!>                The leading dimension of the array A.  LDA >= max(1,NMAX).
!> \endverbatim
!>
!> \param[out] CRF
!> \verbatim
!>          CRF is COMPLEX array, dimension ((NMAX*(NMAX+1))/2).
!> \endverbatim
!>
!> \param[out] A
!> \verbatim
!>          A is COMPLEX array, dimension (LDA,NMAX)
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>                The leading dimension of the array A.  LDA >= max(1,NMAX).
!> \endverbatim
!>
!> \param[out] S_WORK_CLANGE
!> \verbatim
!>          S_WORK_CLANGE is REAL array, dimension (NMAX)
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
!> \date June 2017
!
!> \ingroup complex_lin
!
!  =====================================================================
      SUBROUTINE cdrvrf4(Nout,Nn,Nval,Thresh,C1,C2,Ldc,Crf,A,Lda,       &
     &                   S_work_clange)
      IMPLICIT NONE
!*--CDRVRF4118
!
!  -- LAPACK test routine (version 3.7.1) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     June 2017
!
!     .. Scalar Arguments ..
      INTEGER Lda , Ldc , Nn , Nout
      REAL Thresh
!     ..
!     .. Array Arguments ..
      INTEGER Nval(Nn)
      REAL S_work_clange(*)
      COMPLEX A(Lda,*) , C1(Ldc,*) , C2(Ldc,*) , Crf(*)
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
      REAL SLAMCH , CLANGE
EXTERNAL    SLAMCH,   CLANGE   
!     ..
!     .. External Subroutines ..
EXTERNAL    CHERK,   CHFRK,   CTFTTR,   CTRTTF
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
      DATA forms/'N' , 'C'/
      DATA transs/'N' , 'C'/
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
                                 A(i,j) = CLARND(4,iseed)
                              ENDDO
                           ENDDO
!
                           norma = CLANGE('I',n,k,A,Lda,S_work_clange)
!
                        ELSE
!
!                          In this case we are TRANS, so A is K-by-N
!
                           DO j = 1 , n
                              DO i = 1 , k
                                 A(i,j) = CLARND(4,iseed)
                              ENDDO
                           ENDDO
!
                           norma = CLANGE('I',k,n,A,Lda,S_work_clange)
!
                        ENDIF
!
!
!                       Generate C1 our N--by--N Hermitian matrix.
!                       Make sure C2 has the same upper/lower part,
!                       (the one that we do not touch), so
!                       copy the initial C1 in C2 in it.
!
                        DO j = 1 , n
                           DO i = 1 , n
                              C1(i,j) = CLARND(4,iseed)
                              C2(i,j) = C1(i,j)
                           ENDDO
                        ENDDO
!
!                       (See comment later on for why we use CLANGE and
!                       not CLANHE for C1.)
!
                        normc = CLANGE('I',n,n,C1,Ldc,S_work_clange)
!
                        SRNamt = 'CTRTTF'
                        CALL CTRTTF(cform,uplo,n,C1,Ldc,Crf,info)
!
!                       call zherk the BLAS routine -> gives C1
!
                        SRNamt = 'CHERK '
                        CALL CHERK(uplo,trans,n,k,alpha,A,Lda,beta,C1,  &
     &                             Ldc)
!
!                       call zhfrk the RFP routine -> gives CRF
!
                        SRNamt = 'CHFRK '
                        CALL CHFRK(cform,uplo,trans,n,k,alpha,A,Lda,    &
     &                             beta,Crf)
!
!                       convert CRF in full format -> gives C2
!
                        SRNamt = 'CTFTTR'
                        CALL CTFTTR(cform,uplo,n,Crf,C2,Ldc,info)
!
!                       compare C1 and C2
!
                        DO j = 1 , n
                           DO i = 1 , n
                              C1(i,j) = C1(i,j) - C2(i,j)
                           ENDDO
                        ENDDO
!
!                       Yes, C1 is Hermitian so we could call CLANHE,
!                       but we want to check the upper part that is
!                       supposed to be unchanged and the diagonal that
!                       is supposed to be real -> CLANGE
!
                        result(1) = CLANGE('I',n,n,C1,Ldc,S_work_clange)
                        result(1) = result(1)                           &
     &                              /MAX(ABS(alpha)*norma*norma+        &
     &                              ABS(beta)*normc,ONE)/MAX(n,1)/eps
!
                        IF ( result(1)>=Thresh ) THEN
                           IF ( nfail==0 ) THEN
                              WRITE (Nout,*)
                              WRITE (Nout,FMT=99001)
                           ENDIF
                           WRITE (Nout,FMT=99002) 'CHFRK' , cform ,     &
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
         WRITE (Nout,FMT=99003) 'CHFRK' , nrun
      ELSE
         WRITE (Nout,FMT=99004) 'CHFRK' , nfail , nrun
      ENDIF
!
99001 FORMAT (1X,                                                       &
     &' *** Error(s) or Failure(s) while testing CHFRK                  &
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
!     End of cdrvrf4
!
      END SUBROUTINE cdrvrf4
