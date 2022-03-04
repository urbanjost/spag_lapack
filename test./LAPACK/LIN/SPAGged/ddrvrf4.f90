!*==ddrvrf4.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b DDRVRF4
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DDRVRF4( NOUT, NN, NVAL, THRESH, C1, C2, LDC, CRF, A,
!      +                    LDA, D_WORK_DLANGE )
!
!       .. Scalar Arguments ..
!       INTEGER            LDA, LDC, NN, NOUT
!       DOUBLE PRECISION   THRESH
!       ..
!       .. Array Arguments ..
!       INTEGER            NVAL( NN )
!       DOUBLE PRECISION   A( LDA, * ), C1( LDC, * ), C2( LDC, *),
!      +                   CRF( * ), D_WORK_DLANGE( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DDRVRF4 tests the LAPACK RFP routines:
!>     DSFRK
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
!>          THRESH is DOUBLE PRECISION
!>                The threshold value for the test ratios.  A result is
!>                included in the output file if RESULT >= THRESH.  To
!>                have every test ratio printed, use THRESH = 0.
!> \endverbatim
!>
!> \param[out] C1
!> \verbatim
!>          C1 is DOUBLE PRECISION array,
!>                dimension (LDC,NMAX)
!> \endverbatim
!>
!> \param[out] C2
!> \verbatim
!>          C2 is DOUBLE PRECISION array,
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
!>          CRF is DOUBLE PRECISION array,
!>                dimension ((NMAX*(NMAX+1))/2).
!> \endverbatim
!>
!> \param[out] A
!> \verbatim
!>          A is DOUBLE PRECISION array,
!>                dimension (LDA,NMAX)
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>                The leading dimension of the array A.  LDA >= max(1,NMAX).
!> \endverbatim
!>
!> \param[out] D_WORK_DLANGE
!> \verbatim
!>          D_WORK_DLANGE is DOUBLE PRECISION array, dimension (NMAX)
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
      SUBROUTINE DDRVRF4(Nout,Nn,Nval,Thresh,C1,C2,Ldc,Crf,A,Lda,       &
     &                   D_work_dlange)
      IMPLICIT NONE
!*--DDRVRF4122
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Lda , Ldc , Nn , Nout
      DOUBLE PRECISION Thresh
!     ..
!     .. Array Arguments ..
      INTEGER Nval(Nn)
      DOUBLE PRECISION A(Lda,*) , C1(Ldc,*) , C2(Ldc,*) , Crf(*) ,      &
     &                 D_work_dlange(*)
!     ..
!
!  =====================================================================
!     ..
!     .. Parameters ..
      DOUBLE PRECISION ZERO , ONE
      PARAMETER (ZERO=0.0D+0,ONE=1.0D+0)
      INTEGER NTESTS
      PARAMETER (NTESTS=1)
!     ..
!     .. Local Scalars ..
      CHARACTER uplo , cform , trans
      INTEGER i , iform , iik , iin , info , iuplo , j , k , n , nfail ,&
     &        nrun , ialpha , itrans
      DOUBLE PRECISION alpha , beta , eps , norma , normc
!     ..
!     .. Local Arrays ..
      CHARACTER uplos(2) , forms(2) , transs(2)
      INTEGER iseed(4) , iseedy(4)
      DOUBLE PRECISION result(NTESTS)
!     ..
!     .. External Functions ..
      DOUBLE PRECISION DLAMCH , DLARND , DLANGE
      EXTERNAL DLAMCH , DLARND , DLANGE
!     ..
!     .. External Subroutines ..
      EXTERNAL DSYRK , DSFRK , DTFTTR , DTRTTF
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
      eps = DLAMCH('Precision')
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
                           alpha = DLARND(2,iseed)
                           beta = DLARND(2,iseed)
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
                                 A(i,j) = DLARND(2,iseed)
                              ENDDO
                           ENDDO
!
                           norma = DLANGE('I',n,k,A,Lda,D_work_dlange)
!
 
                        ELSE
!
!                          In this case we are TRANS, so A is K-by-N
!
                           DO j = 1 , n
                              DO i = 1 , k
                                 A(i,j) = DLARND(2,iseed)
                              ENDDO
                           ENDDO
!
                           norma = DLANGE('I',k,n,A,Lda,D_work_dlange)
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
                              C1(i,j) = DLARND(2,iseed)
                              C2(i,j) = C1(i,j)
                           ENDDO
                        ENDDO
!
!                       (See comment later on for why we use DLANGE and
!                       not DLANSY for C1.)
!
                        normc = DLANGE('I',n,n,C1,Ldc,D_work_dlange)
!
                        SRNamt = 'DTRTTF'
                        CALL DTRTTF(cform,uplo,n,C1,Ldc,Crf,info)
!
!                       call dsyrk the BLAS routine -> gives C1
!
                        SRNamt = 'DSYRK '
                        CALL DSYRK(uplo,trans,n,k,alpha,A,Lda,beta,C1,  &
     &                             Ldc)
!
!                       call dsfrk the RFP routine -> gives CRF
!
                        SRNamt = 'DSFRK '
                        CALL DSFRK(cform,uplo,trans,n,k,alpha,A,Lda,    &
     &                             beta,Crf)
!
!                       convert CRF in full format -> gives C2
!
                        SRNamt = 'DTFTTR'
                        CALL DTFTTR(cform,uplo,n,Crf,C2,Ldc,info)
!
!                       compare C1 and C2
!
                        DO j = 1 , n
                           DO i = 1 , n
                              C1(i,j) = C1(i,j) - C2(i,j)
                           ENDDO
                        ENDDO
!
!                       Yes, C1 is symmetric so we could call DLANSY,
!                       but we want to check the upper part that is
!                       supposed to be unchanged and the diagonal that
!                       is supposed to be real -> DLANGE
!
                        result(1) = DLANGE('I',n,n,C1,Ldc,D_work_dlange)
                        result(1) = result(1)                           &
     &                              /MAX(ABS(alpha)*norma+ABS(beta),ONE)&
     &                              /MAX(n,1)/eps
!
                        IF ( result(1)>=Thresh ) THEN
                           IF ( nfail==0 ) THEN
                              WRITE (Nout,*)
                              WRITE (Nout,FMT=99001)
                           ENDIF
                           WRITE (Nout,FMT=99002) 'DSFRK' , cform ,     &
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
         WRITE (Nout,FMT=99003) 'DSFRK' , nrun
      ELSE
         WRITE (Nout,FMT=99004) 'DSFRK' , nfail , nrun
      ENDIF
!
99001 FORMAT (1X,                                                       &
     &' *** Error(s) or Failure(s) while testing DSFRK                  &
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
!     End of DDRVRF4
!
      END SUBROUTINE DDRVRF4
