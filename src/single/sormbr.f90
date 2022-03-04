!*==sormbr.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b SORMBR
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SORMBR + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sormbr.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sormbr.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sormbr.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SORMBR( VECT, SIDE, TRANS, M, N, K, A, LDA, TAU, C,
!                          LDC, WORK, LWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          SIDE, TRANS, VECT
!       INTEGER            INFO, K, LDA, LDC, LWORK, M, N
!       ..
!       .. Array Arguments ..
!       REAL               A( LDA, * ), C( LDC, * ), TAU( * ),
!      $                   WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> If VECT = 'Q', SORMBR overwrites the general real M-by-N matrix C
!> with
!>                 SIDE = 'L'     SIDE = 'R'
!> TRANS = 'N':      Q * C          C * Q
!> TRANS = 'T':      Q**T * C       C * Q**T
!>
!> If VECT = 'P', SORMBR overwrites the general real M-by-N matrix C
!> with
!>                 SIDE = 'L'     SIDE = 'R'
!> TRANS = 'N':      P * C          C * P
!> TRANS = 'T':      P**T * C       C * P**T
!>
!> Here Q and P**T are the orthogonal matrices determined by SGEBRD when
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
!>          matrix reduced by SGEBRD.
!>          If VECT = 'P', the number of rows in the original
!>          matrix reduced by SGEBRD.
!>          K >= 0.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is REAL array, dimension
!>                                (LDA,min(nq,K)) if VECT = 'Q'
!>                                (LDA,nq)        if VECT = 'P'
!>          The vectors which define the elementary reflectors H(i) and
!>          G(i), whose products determine the matrices Q and P, as
!>          returned by SGEBRD.
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
!>          TAU is REAL array, dimension (min(nq,K))
!>          TAU(i) must contain the scalar factor of the elementary
!>          reflector H(i) or G(i) which determines Q or P, as returned
!>          by SGEBRD in the array argument TAUQ or TAUP.
!> \endverbatim
!>
!> \param[in,out] C
!> \verbatim
!>          C is REAL array, dimension (LDC,N)
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
!>          WORK is REAL array, dimension (MAX(1,LWORK))
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
!> \ingroup realOTHERcomputational
!
!  =====================================================================
      SUBROUTINE SORMBR(Vect,Side,Trans,M,N,K,A,Lda,Tau,C,Ldc,Work,     &
     &                  Lwork,Info)
      IMPLICIT NONE
!*--SORMBR200
!
!  -- LAPACK computational routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Side , Trans , Vect
      INTEGER Info , K , Lda , Ldc , Lwork , M , N
!     ..
!     .. Array Arguments ..
      REAL A(Lda,*) , C(Ldc,*) , Tau(*) , Work(*)
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
      LOGICAL applyq , left , lquery , notran
      CHARACTER transt
      INTEGER i1 , i2 , iinfo , lwkopt , mi , nb , ni , nq , nw
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      INTEGER ILAENV
      EXTERNAL ILAENV , LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL SORMLQ , SORMQR , XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX , MIN
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
               nb = ILAENV(1,'SORMQR',Side//Trans,M-1,N,M-1,-1)
            ELSE
               nb = ILAENV(1,'SORMQR',Side//Trans,M,N-1,N-1,-1)
            ENDIF
         ELSEIF ( left ) THEN
            nb = ILAENV(1,'SORMLQ',Side//Trans,M-1,N,M-1,-1)
         ELSE
            nb = ILAENV(1,'SORMLQ',Side//Trans,M,N-1,N-1,-1)
         ENDIF
         lwkopt = MAX(1,nw)*nb
         Work(1) = lwkopt
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('SORMBR',-Info)
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
!           P was determined by a call to SGEBRD with nq > k
!
            CALL SORMLQ(Side,transt,M,N,K,A,Lda,Tau,C,Ldc,Work,Lwork,   &
     &                  iinfo)
         ELSEIF ( nq>1 ) THEN
!
!           P was determined by a call to SGEBRD with nq <= k
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
            CALL SORMLQ(Side,transt,mi,ni,nq-1,A(1,2),Lda,Tau,C(i1,i2), &
     &                  Ldc,Work,Lwork,iinfo)
         ENDIF
!
!        Apply Q
!
      ELSEIF ( nq>=K ) THEN
!
!           Q was determined by a call to SGEBRD with nq >= k
!
         CALL SORMQR(Side,Trans,M,N,K,A,Lda,Tau,C,Ldc,Work,Lwork,iinfo)
      ELSEIF ( nq>1 ) THEN
!
!           Q was determined by a call to SGEBRD with nq < k
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
         CALL SORMQR(Side,Trans,mi,ni,nq-1,A(2,1),Lda,Tau,C(i1,i2),Ldc, &
     &               Work,Lwork,iinfo)
      ENDIF
      Work(1) = lwkopt
!
!     End of SORMBR
!
      END SUBROUTINE SORMBR
