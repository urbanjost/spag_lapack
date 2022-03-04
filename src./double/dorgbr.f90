!*==dorgbr.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b DORGBR
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DORGBR + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dorgbr.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dorgbr.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dorgbr.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DORGBR( VECT, M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          VECT
!       INTEGER            INFO, K, LDA, LWORK, M, N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DORGBR generates one of the real orthogonal matrices Q or P**T
!> determined by DGEBRD when reducing a real matrix A to bidiagonal
!> form: A = Q * B * P**T.  Q and P**T are defined as products of
!> elementary reflectors H(i) or G(i) respectively.
!>
!> If VECT = 'Q', A is assumed to have been an M-by-K matrix, and Q
!> is of order M:
!> if m >= k, Q = H(1) H(2) . . . H(k) and DORGBR returns the first n
!> columns of Q, where m >= n >= k;
!> if m < k, Q = H(1) H(2) . . . H(m-1) and DORGBR returns Q as an
!> M-by-M matrix.
!>
!> If VECT = 'P', A is assumed to have been a K-by-N matrix, and P**T
!> is of order N:
!> if k < n, P**T = G(k) . . . G(2) G(1) and DORGBR returns the first m
!> rows of P**T, where n >= m >= k;
!> if k >= n, P**T = G(n-1) . . . G(2) G(1) and DORGBR returns P**T as
!> an N-by-N matrix.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] VECT
!> \verbatim
!>          VECT is CHARACTER*1
!>          Specifies whether the matrix Q or the matrix P**T is
!>          required, as defined in the transformation applied by DGEBRD:
!>          = 'Q':  generate Q;
!>          = 'P':  generate P**T.
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix Q or P**T to be returned.
!>          M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix Q or P**T to be returned.
!>          N >= 0.
!>          If VECT = 'Q', M >= N >= min(M,K);
!>          if VECT = 'P', N >= M >= min(N,K).
!> \endverbatim
!>
!> \param[in] K
!> \verbatim
!>          K is INTEGER
!>          If VECT = 'Q', the number of columns in the original M-by-K
!>          matrix reduced by DGEBRD.
!>          If VECT = 'P', the number of rows in the original K-by-N
!>          matrix reduced by DGEBRD.
!>          K >= 0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (LDA,N)
!>          On entry, the vectors which define the elementary reflectors,
!>          as returned by DGEBRD.
!>          On exit, the M-by-N matrix Q or P**T.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A. LDA >= max(1,M).
!> \endverbatim
!>
!> \param[in] TAU
!> \verbatim
!>          TAU is DOUBLE PRECISION array, dimension
!>                                (min(M,K)) if VECT = 'Q'
!>                                (min(N,K)) if VECT = 'P'
!>          TAU(i) must contain the scalar factor of the elementary
!>          reflector H(i) or G(i), which determines Q or P**T, as
!>          returned by DGEBRD in its array argument TAUQ or TAUP.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))
!>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The dimension of the array WORK. LWORK >= max(1,min(M,N)).
!>          For optimum performance LWORK >= min(M,N)*NB, where NB
!>          is the optimal blocksize.
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
!> \date April 2012
!
!> \ingroup doubleGBcomputational
!
!  =====================================================================
      SUBROUTINE DORGBR(Vect,M,N,K,A,Lda,Tau,Work,Lwork,Info)
      USE F77KINDS                        
      USE S_DORGLQ
      USE S_DORGQR
      USE S_LSAME
      USE S_XERBLA
      IMPLICIT NONE
!*--DORGBR166
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Vect
      INTEGER :: M
      INTEGER :: N
      INTEGER :: K
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(*) :: Tau
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      INTEGER :: i , iinfo , j , lwkopt , mn
      LOGICAL :: lquery , wantq
!
! End of declarations rewritten by SPAG
!
!     ..
!     .. Array Arguments ..
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
!     ..
!     .. Local Scalars ..
!     ..
!     .. External Functions ..
!     ..
!     .. External Subroutines ..
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      Info = 0
      wantq = LSAME(Vect,'Q')
      mn = MIN(M,N)
      lquery = (Lwork==-1)
      IF ( .NOT.wantq .AND. .NOT.LSAME(Vect,'P') ) THEN
         Info = -1
      ELSEIF ( M<0 ) THEN
         Info = -2
      ELSEIF ( N<0 .OR. (wantq .AND. (N>M .OR. N<MIN(M,K))) .OR.        &
     &         (.NOT.wantq .AND. (M>N .OR. M<MIN(N,K))) ) THEN
         Info = -3
      ELSEIF ( K<0 ) THEN
         Info = -4
      ELSEIF ( Lda<MAX(1,M) ) THEN
         Info = -6
      ELSEIF ( Lwork<MAX(1,mn) .AND. .NOT.lquery ) THEN
         Info = -9
      ENDIF
!
      IF ( Info==0 ) THEN
         Work(1) = 1
         IF ( wantq ) THEN
            IF ( M>=K ) THEN
               CALL DORGQR(M,N,K,A,Lda,Tau,Work,-1,iinfo)
            ELSEIF ( M>1 ) THEN
               CALL DORGQR(M-1,M-1,M-1,A(2,2),Lda,Tau,Work,-1,iinfo)
            ENDIF
         ELSEIF ( K<N ) THEN
            CALL DORGLQ(M,N,K,A,Lda,Tau,Work,-1,iinfo)
         ELSEIF ( N>1 ) THEN
            CALL DORGLQ(N-1,N-1,N-1,A(2,2),Lda,Tau,Work,-1,iinfo)
         ENDIF
         lwkopt = Work(1)
         lwkopt = MAX(lwkopt,mn)
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('DORGBR',-Info)
         RETURN
      ELSEIF ( lquery ) THEN
         Work(1) = lwkopt
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( M==0 .OR. N==0 ) THEN
         Work(1) = 1
         RETURN
      ENDIF
!
      IF ( wantq ) THEN
!
!        Form Q, determined by a call to DGEBRD to reduce an m-by-k
!        matrix
!
         IF ( M>=K ) THEN
!
!           If m >= k, assume m >= n >= k
!
            CALL DORGQR(M,N,K,A,Lda,Tau,Work,Lwork,iinfo)
!
         ELSE
!
!           If m < k, assume m = n
!
!           Shift the vectors which define the elementary reflectors one
!           column to the right, and set the first row and column of Q
!           to those of the unit matrix
!
            DO j = M , 2 , -1
               A(1,j) = ZERO
               DO i = j + 1 , M
                  A(i,j) = A(i,j-1)
               ENDDO
            ENDDO
            A(1,1) = ONE
            DO i = 2 , M
               A(i,1) = ZERO
            ENDDO
!
!              Form Q(2:m,2:m)
!
            IF ( M>1 ) CALL DORGQR(M-1,M-1,M-1,A(2,2),Lda,Tau,Work,     &
     &                             Lwork,iinfo)
         ENDIF
!
!        Form P**T, determined by a call to DGEBRD to reduce a k-by-n
!        matrix
!
      ELSEIF ( K<N ) THEN
!
!           If k < n, assume k <= m <= n
!
         CALL DORGLQ(M,N,K,A,Lda,Tau,Work,Lwork,iinfo)
!
      ELSE
!
!           If k >= n, assume m = n
!
!           Shift the vectors which define the elementary reflectors one
!           row downward, and set the first row and column of P**T to
!           those of the unit matrix
!
         A(1,1) = ONE
         DO i = 2 , N
            A(i,1) = ZERO
         ENDDO
         DO j = 2 , N
            DO i = j - 1 , 2 , -1
               A(i,j) = A(i-1,j)
            ENDDO
            A(1,j) = ZERO
         ENDDO
!
!              Form P**T(2:n,2:n)
!
         IF ( N>1 ) CALL DORGLQ(N-1,N-1,N-1,A(2,2),Lda,Tau,Work,Lwork,  &
     &                          iinfo)
      ENDIF
      Work(1) = lwkopt
!
!     End of DORGBR
!
      END SUBROUTINE DORGBR