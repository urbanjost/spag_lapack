!*==dormbr.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b DORMBR
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DORMBR + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dormbr.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dormbr.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dormbr.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DORMBR( VECT, SIDE, TRANS, M, N, K, A, LDA, TAU, C,
!                          LDC, WORK, LWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          SIDE, TRANS, VECT
!       INTEGER            INFO, K, LDA, LDC, LWORK, M, N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> If VECT = 'Q', DORMBR overwrites the general real M-by-N matrix C
!> with
!>                 SIDE = 'L'     SIDE = 'R'
!> TRANS = 'N':      Q * C          C * Q
!> TRANS = 'T':      Q**T * C       C * Q**T
!>
!> If VECT = 'P', DORMBR overwrites the general real M-by-N matrix C
!> with
!>                 SIDE = 'L'     SIDE = 'R'
!> TRANS = 'N':      P * C          C * P
!> TRANS = 'T':      P**T * C       C * P**T
!>
!> Here Q and P**T are the orthogonal matrices determined by DGEBRD when
!> reducing a real matrix A to bidiagonal form: A = Q * B * P**T. Q and
!> P**T are defined as products of elementary reflectors H(i) and G(i)
!> respectively.
!>
!> Let nq = m if SIDE = 'L' and nq = n if SIDE = 'R'. Thus nq is the
!> order of the orthogonal matrix Q or P**T that is applied.
!>
!> If VECT = 'Q', A is assumed to have been an NQ-by-K matrix:
!> if nq >= k, Q = H(1) H(2) . . . H(k);
!> if nq < k, Q = H(1) H(2) . . . H(nq-1).
!>
!> If VECT = 'P', A is assumed to have been a K-by-NQ matrix:
!> if k < nq, P = G(1) G(2) . . . G(k);
!> if k >= nq, P = G(1) G(2) . . . G(nq-1).
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] VECT
!> \verbatim
!>          VECT is CHARACTER*1
!>          = 'Q': apply Q or Q**T;
!>          = 'P': apply P or P**T.
!> \endverbatim
!>
!> \param[in] SIDE
!> \verbatim
!>          SIDE is CHARACTER*1
!>          = 'L': apply Q, Q**T, P or P**T from the Left;
!>          = 'R': apply Q, Q**T, P or P**T from the Right.
!> \endverbatim
!>
!> \param[in] TRANS
!> \verbatim
!>          TRANS is CHARACTER*1
!>          = 'N':  No transpose, apply Q  or P;
!>          = 'T':  Transpose, apply Q**T or P**T.
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix C. M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix C. N >= 0.
!> \endverbatim
!>
!> \param[in] K
!> \verbatim
!>          K is INTEGER
!>          If VECT = 'Q', the number of columns in the original
!>          matrix reduced by DGEBRD.
!>          If VECT = 'P', the number of rows in the original
!>          matrix reduced by DGEBRD.
!>          K >= 0.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension
!>                                (LDA,min(nq,K)) if VECT = 'Q'
!>                                (LDA,nq)        if VECT = 'P'
!>          The vectors which define the elementary reflectors H(i) and
!>          G(i), whose products determine the matrices Q and P, as
!>          returned by DGEBRD.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.
!>          If VECT = 'Q', LDA >= max(1,nq);
!>          if VECT = 'P', LDA >= max(1,min(nq,K)).
!> \endverbatim
!>
!> \param[in] TAU
!> \verbatim
!>          TAU is DOUBLE PRECISION array, dimension (min(nq,K))
!>          TAU(i) must contain the scalar factor of the elementary
!>          reflector H(i) or G(i) which determines Q or P, as returned
!>          by DGEBRD in the array argument TAUQ or TAUP.
!> \endverbatim
!>
!> \param[in,out] C
!> \verbatim
!>          C is DOUBLE PRECISION array, dimension (LDC,N)
!>          On entry, the M-by-N matrix C.
!>          On exit, C is overwritten by Q*C or Q**T*C or C*Q**T or C*Q
!>          or P*C or P**T*C or C*P or C*P**T.
!> \endverbatim
!>
!> \param[in] LDC
!> \verbatim
!>          LDC is INTEGER
!>          The leading dimension of the array C. LDC >= max(1,M).
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
!>          The dimension of the array WORK.
!>          If SIDE = 'L', LWORK >= max(1,N);
!>          if SIDE = 'R', LWORK >= max(1,M).
!>          For optimum performance LWORK >= N*NB if SIDE = 'L', and
!>          LWORK >= M*NB if SIDE = 'R', where NB is the optimal
!>          blocksize.
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
!> \ingroup doubleOTHERcomputational
!
!  =====================================================================
      SUBROUTINE DORMBR(Vect,Side,Trans,M,N,K,A,Lda,Tau,C,Ldc,Work,     &
     &                  Lwork,Info)
      USE F77KINDS                        
      USE S_DORMLQ
      USE S_DORMQR
      USE S_ILAENV
      USE S_LSAME
      USE S_XERBLA
      IMPLICIT NONE
!*--DORMBR205
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Vect
      CHARACTER :: Side
      CHARACTER :: Trans
      INTEGER :: M
      INTEGER :: N
      INTEGER :: K
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(*) :: Tau
      REAL(R8KIND) , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      LOGICAL :: applyq , left , lquery , notran
      INTEGER :: i1 , i2 , iinfo , lwkopt , mi , nb , ni , nq , nw
      CHARACTER :: transt
!
! End of declarations rewritten by SPAG
!
!     ..
!     .. Array Arguments ..
!     ..
!
!  =====================================================================
!
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
      applyq = LSAME(Vect,'Q')
      left = LSAME(Side,'L')
      notran = LSAME(Trans,'N')
      lquery = (Lwork==-1)
!
!     NQ is the order of Q or P and NW is the minimum dimension of WORK
!
      IF ( left ) THEN
         nq = M
         nw = N
      ELSE
         nq = N
         nw = M
      ENDIF
      IF ( .NOT.applyq .AND. .NOT.LSAME(Vect,'P') ) THEN
         Info = -1
      ELSEIF ( .NOT.left .AND. .NOT.LSAME(Side,'R') ) THEN
         Info = -2
      ELSEIF ( .NOT.notran .AND. .NOT.LSAME(Trans,'T') ) THEN
         Info = -3
      ELSEIF ( M<0 ) THEN
         Info = -4
      ELSEIF ( N<0 ) THEN
         Info = -5
      ELSEIF ( K<0 ) THEN
         Info = -6
      ELSEIF ( (applyq .AND. Lda<MAX(1,nq)) .OR.                        &
     &         (.NOT.applyq .AND. Lda<MAX(1,MIN(nq,K))) ) THEN
         Info = -8
      ELSEIF ( Ldc<MAX(1,M) ) THEN
         Info = -11
      ELSEIF ( Lwork<MAX(1,nw) .AND. .NOT.lquery ) THEN
         Info = -13
      ENDIF
!
      IF ( Info==0 ) THEN
         IF ( applyq ) THEN
            IF ( left ) THEN
               nb = ILAENV(1,'DORMQR',Side//Trans,M-1,N,M-1,-1)
            ELSE
               nb = ILAENV(1,'DORMQR',Side//Trans,M,N-1,N-1,-1)
            ENDIF
         ELSEIF ( left ) THEN
            nb = ILAENV(1,'DORMLQ',Side//Trans,M-1,N,M-1,-1)
         ELSE
            nb = ILAENV(1,'DORMLQ',Side//Trans,M,N-1,N-1,-1)
         ENDIF
         lwkopt = MAX(1,nw)*nb
         Work(1) = lwkopt
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('DORMBR',-Info)
         RETURN
      ELSEIF ( lquery ) THEN
         RETURN
      ENDIF
!
!     Quick return if possible
!
      Work(1) = 1
      IF ( M==0 .OR. N==0 ) RETURN
!
      IF ( .NOT.(applyq) ) THEN
!
!        Apply P
!
         IF ( notran ) THEN
            transt = 'T'
         ELSE
            transt = 'N'
         ENDIF
         IF ( nq>K ) THEN
!
!           P was determined by a call to DGEBRD with nq > k
!
            CALL DORMLQ(Side,transt,M,N,K,A,Lda,Tau,C,Ldc,Work,Lwork,   &
     &                  iinfo)
         ELSEIF ( nq>1 ) THEN
!
!           P was determined by a call to DGEBRD with nq <= k
!
            IF ( left ) THEN
               mi = M - 1
               ni = N
               i1 = 2
               i2 = 1
            ELSE
               mi = M
               ni = N - 1
               i1 = 1
               i2 = 2
            ENDIF
            CALL DORMLQ(Side,transt,mi,ni,nq-1,A(1,2),Lda,Tau,C(i1,i2), &
     &                  Ldc,Work,Lwork,iinfo)
         ENDIF
!
!        Apply Q
!
      ELSEIF ( nq>=K ) THEN
!
!           Q was determined by a call to DGEBRD with nq >= k
!
         CALL DORMQR(Side,Trans,M,N,K,A,Lda,Tau,C,Ldc,Work,Lwork,iinfo)
      ELSEIF ( nq>1 ) THEN
!
!           Q was determined by a call to DGEBRD with nq < k
!
         IF ( left ) THEN
            mi = M - 1
            ni = N
            i1 = 2
            i2 = 1
         ELSE
            mi = M
            ni = N - 1
            i1 = 1
            i2 = 2
         ENDIF
         CALL DORMQR(Side,Trans,mi,ni,nq-1,A(2,1),Lda,Tau,C(i1,i2),Ldc, &
     &               Work,Lwork,iinfo)
      ENDIF
      Work(1) = lwkopt
!
!     End of DORMBR
!
      END SUBROUTINE DORMBR
