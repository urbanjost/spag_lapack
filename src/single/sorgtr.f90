!*==sorgtr.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b SORGTR
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SORGTR + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sorgtr.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sorgtr.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sorgtr.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SORGTR( UPLO, N, A, LDA, TAU, WORK, LWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            INFO, LDA, LWORK, N
!       ..
!       .. Array Arguments ..
!       REAL               A( LDA, * ), TAU( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SORGTR generates a real orthogonal matrix Q which is defined as the
!> product of n-1 elementary reflectors of order N, as returned by
!> SSYTRD:
!>
!> if UPLO = 'U', Q = H(n-1) . . . H(2) H(1),
!>
!> if UPLO = 'L', Q = H(1) H(2) . . . H(n-1).
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          = 'U': Upper triangle of A contains elementary reflectors
!>                 from SSYTRD;
!>          = 'L': Lower triangle of A contains elementary reflectors
!>                 from SSYTRD.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix Q. N >= 0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is REAL array, dimension (LDA,N)
!>          On entry, the vectors which define the elementary reflectors,
!>          as returned by SSYTRD.
!>          On exit, the N-by-N orthogonal matrix Q.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A. LDA >= max(1,N).
!> \endverbatim
!>
!> \param[in] TAU
!> \verbatim
!>          TAU is REAL array, dimension (N-1)
!>          TAU(i) must contain the scalar factor of the elementary
!>          reflector H(i), as returned by SSYTRD.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is REAL array, dimension (MAX(1,LWORK))
!>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The dimension of the array WORK. LWORK >= max(1,N-1).
!>          For optimum performance LWORK >= (N-1)*NB, where NB is
!>          the optimal blocksize.
!>
!>          If LWORK = -1, then a workspace query is assumed; the routine
!>          only calculates the optimal size of the WORK array, returns
!>          this value as the first entry of the WORK array, and no error
!>          message related to LWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value
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
!> \ingroup realOTHERcomputational
!
!  =====================================================================
      SUBROUTINE SORGTR(Uplo,N,A,Lda,Tau,Work,Lwork,Info)
      IMPLICIT NONE
!*--SORGTR127
!
!  -- LAPACK computational routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Uplo
      INTEGER Info , Lda , Lwork , N
!     ..
!     .. Array Arguments ..
      REAL A(Lda,*) , Tau(*) , Work(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ZERO , ONE
      PARAMETER (ZERO=0.0E+0,ONE=1.0E+0)
!     ..
!     .. Local Scalars ..
      LOGICAL lquery , upper
      INTEGER i , iinfo , j , lwkopt , nb
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      INTEGER ILAENV
      EXTERNAL ILAENV , LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL SORGQL , SORGQR , XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      Info = 0
      lquery = (Lwork==-1)
      upper = LSAME(Uplo,'U')
      IF ( .NOT.upper .AND. .NOT.LSAME(Uplo,'L') ) THEN
         Info = -1
      ELSEIF ( N<0 ) THEN
         Info = -2
      ELSEIF ( Lda<MAX(1,N) ) THEN
         Info = -4
      ELSEIF ( Lwork<MAX(1,N-1) .AND. .NOT.lquery ) THEN
         Info = -7
      ENDIF
!
      IF ( Info==0 ) THEN
         IF ( upper ) THEN
            nb = ILAENV(1,'SORGQL',' ',N-1,N-1,N-1,-1)
         ELSE
            nb = ILAENV(1,'SORGQR',' ',N-1,N-1,N-1,-1)
         ENDIF
         lwkopt = MAX(1,N-1)*nb
         Work(1) = lwkopt
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('SORGTR',-Info)
         RETURN
      ELSEIF ( lquery ) THEN
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( N==0 ) THEN
         Work(1) = 1
         RETURN
      ENDIF
!
      IF ( upper ) THEN
!
!        Q was determined by a call to SSYTRD with UPLO = 'U'
!
!        Shift the vectors which define the elementary reflectors one
!        column to the left, and set the last row and column of Q to
!        those of the unit matrix
!
         DO j = 1 , N - 1
            DO i = 1 , j - 1
               A(i,j) = A(i,j+1)
            ENDDO
            A(N,j) = ZERO
         ENDDO
         DO i = 1 , N - 1
            A(i,N) = ZERO
         ENDDO
         A(N,N) = ONE
!
!        Generate Q(1:n-1,1:n-1)
!
         CALL SORGQL(N-1,N-1,N-1,A,Lda,Tau,Work,Lwork,iinfo)
!
      ELSE
!
!        Q was determined by a call to SSYTRD with UPLO = 'L'.
!
!        Shift the vectors which define the elementary reflectors one
!        column to the right, and set the first row and column of Q to
!        those of the unit matrix
!
         DO j = N , 2 , -1
            A(1,j) = ZERO
            DO i = j + 1 , N
               A(i,j) = A(i,j-1)
            ENDDO
         ENDDO
         A(1,1) = ONE
         DO i = 2 , N
            A(i,1) = ZERO
         ENDDO
!
!           Generate Q(2:n,2:n)
!
         IF ( N>1 ) CALL SORGQR(N-1,N-1,N-1,A(2,2),Lda,Tau,Work,Lwork,  &
     &                          iinfo)
      ENDIF
      Work(1) = lwkopt
!
!     End of SORGTR
!
      END SUBROUTINE SORGTR
