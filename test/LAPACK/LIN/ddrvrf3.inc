!*==ddrvrf3.f90  processed by SPAG 7.51RB at 17:34 on  4 Mar 2022
!> \brief \b ddrvrf3
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ddrvrf3( NOUT, NN, NVAL, THRESH, A, LDA, ARF, B1, B2,
!      +                    D_WORK_DLANGE, D_WORK_DGEQRF, TAU )
!
!       .. Scalar Arguments ..
!       INTEGER            LDA, NN, NOUT
!       DOUBLE PRECISION   THRESH
!       ..
!       .. Array Arguments ..
!       INTEGER            NVAL( NN )
!       DOUBLE PRECISION   A( LDA, * ), ARF( * ), B1( LDA, * ),
!      +                   B2( LDA, * ), D_WORK_DGEQRF( * ),
!      +                   D_WORK_DLANGE( * ), TAU( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ddrvrf3 tests the LAPACK RFP routines:
!>     DTFSM
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
!>                included in the output file if RESULT >= THRESH.  To have
!>                every test ratio printed, use THRESH = 0.
!> \endverbatim
!>
!> \param[out] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (LDA,NMAX)
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>                The leading dimension of the array A.  LDA >= max(1,NMAX).
!> \endverbatim
!>
!> \param[out] ARF
!> \verbatim
!>          ARF is DOUBLE PRECISION array, dimension ((NMAX*(NMAX+1))/2).
!> \endverbatim
!>
!> \param[out] B1
!> \verbatim
!>          B1 is DOUBLE PRECISION array, dimension (LDA,NMAX)
!> \endverbatim
!>
!> \param[out] B2
!> \verbatim
!>          B2 is DOUBLE PRECISION array, dimension (LDA,NMAX)
!> \endverbatim
!>
!> \param[out] D_WORK_DLANGE
!> \verbatim
!>          D_WORK_DLANGE is DOUBLE PRECISION array, dimension (NMAX)
!> \endverbatim
!>
!> \param[out] D_WORK_DGEQRF
!> \verbatim
!>          D_WORK_DGEQRF is DOUBLE PRECISION array, dimension (NMAX)
!> \endverbatim
!>
!> \param[out] TAU
!> \verbatim
!>          TAU is DOUBLE PRECISION array, dimension (NMAX)
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
!> \ingroup double_lin
!
!  =====================================================================
      SUBROUTINE ddrvrf3(Nout,Nn,Nval,Thresh,A,Lda,Arf,B1,B2,           &
     &                   D_work_dlange,D_work_dgeqrf,Tau)
      IMPLICIT NONE
!*--DDRVRF3122
!
!  -- LAPACK test routine (version 3.7.1) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     June 2017
!
!     .. Scalar Arguments ..
      INTEGER Lda , Nn , Nout
      DOUBLE PRECISION Thresh
!     ..
!     .. Array Arguments ..
      INTEGER Nval(Nn)
      DOUBLE PRECISION A(Lda,*) , Arf(*) , B1(Lda,*) , B2(Lda,*) ,      &
     &                 D_work_dgeqrf(*) , D_work_dlange(*) , Tau(*)
!     ..
!
!  =====================================================================
!     ..
!     .. Parameters ..
      DOUBLE PRECISION ZERO , ONE
      PARAMETER (ZERO=(0.0D+0,0.0D+0),ONE=(1.0D+0,0.0D+0))
      INTEGER NTESTS
      PARAMETER (NTESTS=1)
!     ..
!     .. Local Scalars ..
      CHARACTER uplo , cform , diag , trans , side
      INTEGER i , iform , iim , iin , info , iuplo , j , m , n , na ,   &
     &        nfail , nrun , iside , idiag , ialpha , itrans
      DOUBLE PRECISION eps , alpha
!     ..
!     .. Local Arrays ..
      CHARACTER uplos(2) , forms(2) , transs(2) , diags(2) , sides(2)
      INTEGER iseed(4) , iseedy(4)
      DOUBLE PRECISION result(NTESTS)
!     ..
!     .. External Functions ..
      DOUBLE PRECISION DLAMCH , DLANGE , DLARND
EXTERNAL    DLAMCH,   DLANGE,   DLARND
!     ..
!     .. External Subroutines ..
EXTERNAL    DTRTTF,   DGEQRF,   DGEQLF,   DTFSM,   DTRSM
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX , SQRT
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
      DATA sides/'L' , 'R'/
      DATA transs/'N' , 'T'/
      DATA diags/'N' , 'U'/
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
      DO iim = 1 , Nn
!
         m = Nval(iim)
!
         DO iin = 1 , Nn
!
            n = Nval(iin)
!
            DO iform = 1 , 2
!
               cform = forms(iform)
!
               DO iuplo = 1 , 2
!
                  uplo = uplos(iuplo)
!
                  DO iside = 1 , 2
!
                     side = sides(iside)
!
                     DO itrans = 1 , 2
!
                        trans = transs(itrans)
!
                        DO idiag = 1 , 2
!
                           diag = diags(idiag)
!
                           DO ialpha = 1 , 3
!
                              IF ( ialpha==1 ) THEN
                                 alpha = ZERO
                              ELSEIF ( ialpha==2 ) THEN
                                 alpha = ONE
                              ELSE
                                 alpha = DLARND(2,iseed)
                              ENDIF
!
!                             All the parameters are set:
!                                CFORM, SIDE, UPLO, TRANS, DIAG, M, N,
!                                and ALPHA
!                             READY TO TEST!
!
                              nrun = nrun + 1
!
                              IF ( iside==1 ) THEN
!
!                                The case ISIDE.EQ.1 is when SIDE.EQ.'L'
!                                -> A is M-by-M ( B is M-by-N )
!
                                 na = m
!
                              ELSE
!
!                                The case ISIDE.EQ.2 is when SIDE.EQ.'R'
!                                -> A is N-by-N ( B is M-by-N )
!
                                 na = n
!
                              ENDIF
!
!                             Generate A our NA--by--NA triangular
!                             matrix.
!                             Our test is based on forward error so we
!                             do want A to be well conditioned! To get
!                             a well-conditioned triangular matrix, we
!                             take the R factor of the QR/LQ factorization
!                             of a random matrix.
!
                              DO j = 1 , na
                                 DO i = 1 , na
                                    A(i,j) = DLARND(2,iseed)
                                 ENDDO
                              ENDDO
!
                              IF ( iuplo==1 ) THEN
!
!                                The case IUPLO.EQ.1 is when SIDE.EQ.'U'
!                                -> QR factorization.
!
                                 SRNamt = 'DGEQRF'
                                 CALL DGEQRF(na,na,A,Lda,Tau,           &
     &                              D_work_dgeqrf,Lda,info)
                              ELSE
!
!                                The case IUPLO.EQ.2 is when SIDE.EQ.'L'
!                                -> QL factorization.
!
                                 SRNamt = 'DGELQF'
                                 CALL DGELQF(na,na,A,Lda,Tau,           &
     &                              D_work_dgeqrf,Lda,info)
                              ENDIF
!
!                             Store a copy of A in RFP format (in ARF).
!
                              SRNamt = 'DTRTTF'
                              CALL DTRTTF(cform,uplo,na,A,Lda,Arf,info)
!
!                             Generate B1 our M--by--N right-hand side
!                             and store a copy in B2.
!
                              DO j = 1 , n
                                 DO i = 1 , m
                                    B1(i,j) = DLARND(2,iseed)
                                    B2(i,j) = B1(i,j)
                                 ENDDO
                              ENDDO
!
!                             Solve op( A ) X = B or X op( A ) = B
!                             with DTRSM
!
                              SRNamt = 'DTRSM'
                              CALL DTRSM(side,uplo,trans,diag,m,n,alpha,&
     &                           A,Lda,B1,Lda)
!
!                             Solve op( A ) X = B or X op( A ) = B
!                             with DTFSM
!
                              SRNamt = 'DTFSM'
                              CALL DTFSM(cform,side,uplo,trans,diag,m,n,&
     &                           alpha,Arf,B2,Lda)
!
!                             Check that the result agrees.
!
                              DO j = 1 , n
                                 DO i = 1 , m
                                    B1(i,j) = B2(i,j) - B1(i,j)
                                 ENDDO
                              ENDDO
!
                              result(1)                                 &
     &                           = DLANGE('I',m,n,B1,Lda,D_work_dlange)
!
                              result(1) = result(1)/SQRT(eps)           &
     &                           /MAX(MAX(m,n),1)
!
                              IF ( result(1)>=Thresh ) THEN
                                 IF ( nfail==0 ) THEN
                                    WRITE (Nout,*)
                                    WRITE (Nout,FMT=99001)
                                 ENDIF
                                 WRITE (Nout,FMT=99002) 'DTFSM' ,       &
     &                                  cform , side , uplo , trans ,   &
     &                                  diag , m , n , result(1)
                                 nfail = nfail + 1
                              ENDIF
!
                           ENDDO
                        ENDDO
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
         WRITE (Nout,FMT=99003) 'DTFSM' , nrun
      ELSE
         WRITE (Nout,FMT=99004) 'DTFSM' , nfail , nrun
      ENDIF
!
99001 FORMAT (1X,                                                       &
     &' *** Error(s) or Failure(s) while testing DTFSM                  &
     &                                                                  &
     &                                                                  &
     &                                           ***')
99002 FORMAT (1X,'     Failure in ',A5,', CFORM=''',A1,''',',' SIDE=''',&
     &        A1,''',',' UPLO=''',A1,''',',' TRANS=''',A1,''',',        &
     &        ' DIAG=''',A1,''',',' M=',I3,', N =',I3,', test=',G12.5)
99003 FORMAT (1X,'All tests for ',A5,' auxiliary routine passed the ',  &
     &        'threshold ( ',I5,' tests run)')
99004 FORMAT (1X,A6,' auxiliary routine: ',I5,' out of ',I5,            &
     &        ' tests failed to pass the threshold')
!
!
!     End of ddrvrf3
!
      END SUBROUTINE ddrvrf3
