!*==cdrvrf2.f90  processed by SPAG 7.51RB at 17:34 on  4 Mar 2022
!> \brief \b cdrvrf2
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE cdrvrf2( NOUT, NN, NVAL, A, LDA, ARF, AP, ASAV  )
!
!       .. Scalar Arguments ..
!       INTEGER            LDA, NN, NOUT
!       ..
!       .. Array Arguments ..
!       INTEGER            NVAL( NN )
!       COMPLEX            A( LDA, * ), ARF( * ), AP(*), ASAV( LDA, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> cdrvrf2 tests the LAPACK RFP conversion routines.
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
!> \param[out] ARF
!> \verbatim
!>          ARF is COMPLEX array, dimension ((NMAX*(NMAX+1))/2).
!> \endverbatim
!>
!> \param[out] AP
!> \verbatim
!>          AP is COMPLEX array, dimension ((NMAX*(NMAX+1))/2).
!> \endverbatim
!>
!> \param[out] ASAV
!> \verbatim
!>          ASAV is COMPLEX6 array, dimension (LDA,NMAX)
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
      SUBROUTINE cdrvrf2(Nout,Nn,Nval,A,Lda,Arf,Ap,Asav)
      IMPLICIT NONE
!*--CDRVRF293
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Lda , Nn , Nout
!     ..
!     .. Array Arguments ..
      INTEGER Nval(Nn)
      COMPLEX A(Lda,*) , Arf(*) , Ap(*) , Asav(Lda,*)
!     ..
!
!  =====================================================================
!     ..
!     .. Local Scalars ..
      LOGICAL lower , ok1 , ok2
      CHARACTER uplo , cform
      INTEGER i , iform , iin , info , iuplo , j , n , nerrs , nrun
!     ..
!     .. Local Arrays ..
      CHARACTER uplos(2) , forms(2)
      INTEGER iseed(4) , iseedy(4)
!     ..
!     .. External Functions ..
!     ..
!     .. External Subroutines ..
EXTERNAL    CTFTTR,   CTFTTP,   CTRTTF,   CTRTTP,   CTPTTR,   CTPTTF
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
      nerrs = 0
      info = 0
      DO i = 1 , 4
         iseed(i) = iseedy(i)
      ENDDO
!
      DO iin = 1 , Nn
!
         n = Nval(iin)
!
!        Do first for UPLO = 'U', then for UPLO = 'L'
!
         DO iuplo = 1 , 2
!
            uplo = uplos(iuplo)
            lower = .TRUE.
            IF ( iuplo==1 ) lower = .FALSE.
!
!           Do first for CFORM = 'N', then for CFORM = 'C'
!
            DO iform = 1 , 2
!
               cform = forms(iform)
!
               nrun = nrun + 1
!
               DO j = 1 , n
                  DO i = 1 , n
                     A(i,j) = CLARND(4,iseed)
                  ENDDO
               ENDDO
!
               SRNamt = 'CTRTTF'
               CALL CTRTTF(cform,uplo,n,A,Lda,Arf,info)
!
               SRNamt = 'CTFTTP'
               CALL CTFTTP(cform,uplo,n,Arf,Ap,info)
!
               SRNamt = 'CTPTTR'
               CALL CTPTTR(uplo,n,Ap,Asav,Lda,info)
!
               ok1 = .TRUE.
               IF ( lower ) THEN
                  DO j = 1 , n
                     DO i = j , n
                        IF ( A(i,j)/=Asav(i,j) ) ok1 = .FALSE.
                     ENDDO
                  ENDDO
               ELSE
                  DO j = 1 , n
                     DO i = 1 , j
                        IF ( A(i,j)/=Asav(i,j) ) ok1 = .FALSE.
                     ENDDO
                  ENDDO
               ENDIF
!
               nrun = nrun + 1
!
               SRNamt = 'CTRTTP'
               CALL CTRTTP(uplo,n,A,Lda,Ap,info)
!
               SRNamt = 'CTPTTF'
               CALL CTPTTF(cform,uplo,n,Ap,Arf,info)
!
               SRNamt = 'CTFTTR'
               CALL CTFTTR(cform,uplo,n,Arf,Asav,Lda,info)
!
               ok2 = .TRUE.
               IF ( lower ) THEN
                  DO j = 1 , n
                     DO i = j , n
                        IF ( A(i,j)/=Asav(i,j) ) ok2 = .FALSE.
                     ENDDO
                  ENDDO
               ELSE
                  DO j = 1 , n
                     DO i = 1 , j
                        IF ( A(i,j)/=Asav(i,j) ) ok2 = .FALSE.
                     ENDDO
                  ENDDO
               ENDIF
!
               IF ( (.NOT.ok1) .OR. (.NOT.ok2) ) THEN
                  IF ( nerrs==0 ) THEN
                     WRITE (Nout,*)
                     WRITE (Nout,FMT=99001)
                  ENDIF
                  WRITE (Nout,FMT=99002) n , uplo , cform
                  nerrs = nerrs + 1
               ENDIF
!
            ENDDO
         ENDDO
      ENDDO
!
!     Print a summary of the results.
!
      IF ( nerrs==0 ) THEN
         WRITE (Nout,FMT=99003) nrun
      ELSE
         WRITE (Nout,FMT=99004) nerrs , nrun
      ENDIF
!
99001 FORMAT (1X,' *** Error(s) while testing the RFP conversion',      &
     &        ' routines ***')
99002 FORMAT (1X,'     Error in RFP,conversion routines N=',I5,         &
     &        ' UPLO=''',A1,''', FORM =''',A1,'''')
99003 FORMAT (1X,'All tests for the RFP conversion routines passed ( ', &
     &        I5,' tests run)')
99004 FORMAT (1X,'RFP conversion routines: ',I5,' out of ',I5,          &
     &        ' error message recorded')
!
!
!     End of cdrvrf2
!
      END SUBROUTINE cdrvrf2
