!*==ddrvrf1.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b DDRVRF1
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DDRVRF1( NOUT, NN, NVAL, THRESH, A, LDA, ARF, WORK )
!
!       .. Scalar Arguments ..
!       INTEGER            LDA, NN, NOUT
!       DOUBLE PRECISION   THRESH
!       ..
!       .. Array Arguments ..
!       INTEGER            NVAL( NN )
!       DOUBLE PRECISION   A( LDA, * ), ARF( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DDRVRF1 tests the LAPACK RFP routines:
!>     DLANSF
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
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension ( NMAX )
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
      SUBROUTINE DDRVRF1(Nout,Nn,Nval,Thresh,A,Lda,Arf,Work)
      IMPLICIT NONE
!*--DDRVRF198
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Lda , Nn , Nout
      DOUBLE PRECISION Thresh
!     ..
!     .. Array Arguments ..
      INTEGER Nval(Nn)
      DOUBLE PRECISION A(Lda,*) , Arf(*) , Work(*)
!     ..
!
!  =====================================================================
!     ..
!     .. Parameters ..
      DOUBLE PRECISION ONE
      PARAMETER (ONE=1.0D+0)
      INTEGER NTESTS
      PARAMETER (NTESTS=1)
!     ..
!     .. Local Scalars ..
      CHARACTER uplo , cform , norm
      INTEGER i , iform , iin , iit , info , inorm , iuplo , j , n ,    &
     &        nerrs , nfail , nrun
      DOUBLE PRECISION eps , large , norma , normarf , small
!     ..
!     .. Local Arrays ..
      CHARACTER uplos(2) , forms(2) , norms(4)
      INTEGER iseed(4) , iseedy(4)
      DOUBLE PRECISION result(NTESTS)
!     ..
!     .. External Functions ..
      DOUBLE PRECISION DLAMCH , DLANSY , DLANSF , DLARND
      EXTERNAL DLAMCH , DLANSY , DLANSF , DLARND
!     ..
!     .. External Subroutines ..
      EXTERNAL DTRTTF
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
      DATA norms/'M' , '1' , 'I' , 'F'/
!     ..
!     .. Executable Statements ..
!
!     Initialize constants and the random number seed.
!
      nrun = 0
      nfail = 0
      nerrs = 0
      info = 0
      DO i = 1 , 4
         iseed(i) = iseedy(i)
      ENDDO
!
      eps = DLAMCH('Precision')
      small = DLAMCH('Safe minimum')
      large = ONE/small
      small = small*Lda*Lda
      large = large/Lda/Lda
!
      DO iin = 1 , Nn
!
         n = Nval(iin)
!
         DO iit = 1 , 3
!           Nothing to do for N=0
            IF ( n==0 ) EXIT
!
!           IIT = 1 : random matrix
!           IIT = 2 : random matrix scaled near underflow
!           IIT = 3 : random matrix scaled near overflow
!
            DO j = 1 , n
               DO i = 1 , n
                  A(i,j) = DLARND(2,iseed)
               ENDDO
            ENDDO
!
            IF ( iit==2 ) THEN
               DO j = 1 , n
                  DO i = 1 , n
                     A(i,j) = A(i,j)*large
                  ENDDO
               ENDDO
            ENDIF
!
            IF ( iit==3 ) THEN
               DO j = 1 , n
                  DO i = 1 , n
                     A(i,j) = A(i,j)*small
                  ENDDO
               ENDDO
            ENDIF
!
!           Do first for UPLO = 'U', then for UPLO = 'L'
!
            DO iuplo = 1 , 2
!
               uplo = uplos(iuplo)
!
!              Do first for CFORM = 'N', then for CFORM = 'C'
!
               DO iform = 1 , 2
!
                  cform = forms(iform)
!
                  SRNamt = 'DTRTTF'
                  CALL DTRTTF(cform,uplo,n,A,Lda,Arf,info)
!
!                 Check error code from DTRTTF
!
                  IF ( info/=0 ) THEN
                     IF ( nfail==0 .AND. nerrs==0 ) THEN
                        WRITE (Nout,*)
                        WRITE (Nout,FMT=99001)
                     ENDIF
                     WRITE (Nout,FMT=99002) SRNamt , uplo , cform , n
                     nerrs = nerrs + 1
                     CYCLE
                  ENDIF
!
                  DO inorm = 1 , 4
!
!                    Check all four norms: 'M', '1', 'I', 'F'
!
                     norm = norms(inorm)
                     normarf = DLANSF(norm,cform,uplo,n,Arf,Work)
                     norma = DLANSY(norm,uplo,n,A,Lda,Work)
!
                     result(1) = (norma-normarf)/norma/eps
                     nrun = nrun + 1
!
                     IF ( result(1)>=Thresh ) THEN
                        IF ( nfail==0 .AND. nerrs==0 ) THEN
                           WRITE (Nout,*)
                           WRITE (Nout,FMT=99001)
                        ENDIF
                        WRITE (Nout,FMT=99003) 'DLANSF' , n , iit ,     &
     &                         uplo , cform , norm , result(1)
                        nfail = nfail + 1
                     ENDIF
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO
!
!     Print a summary of the results.
!
      IF ( nfail==0 ) THEN
         WRITE (Nout,FMT=99004) 'DLANSF' , nrun
      ELSE
         WRITE (Nout,FMT=99005) 'DLANSF' , nfail , nrun
      ENDIF
      IF ( nerrs/=0 ) WRITE (Nout,FMT=99006) nerrs , 'DLANSF'
!
99001 FORMAT (1X,                                                       &
     &' *** Error(s) or Failure(s) while testing DLANSF                 &
     &                                                                  &
     &                                                                  &
     &                                           ***')
99002 FORMAT (1X,'     Error in ',A6,' with UPLO=''',A1,''', FORM=''',  &
     &        A1,''', N=',I5)
99003 FORMAT (1X,'     Failure in ',A6,' N=',I5,' TYPE=',I5,' UPLO=''', &
     &        A1,''', FORM =''',A1,''', NORM=''',A1,''', test=',G12.5)
99004 FORMAT (1X,'All tests for ',A6,' auxiliary routine passed the ',  &
     &        'threshold ( ',I5,' tests run)')
99005 FORMAT (1X,A6,' auxiliary routine: ',I5,' out of ',I5,            &
     &        ' tests failed to pass the threshold')
99006 FORMAT (26X,I5,' error message recorded (',A6,')')
!
!
!     End of DDRVRF1
!
      END SUBROUTINE DDRVRF1
