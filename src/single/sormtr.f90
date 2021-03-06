!*==sormtr.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b SORMTR
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SORMTR + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sormtr.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sormtr.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sormtr.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SORMTR( SIDE, UPLO, TRANS, M, N, A, LDA, TAU, C, LDC,
!                          WORK, LWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          SIDE, TRANS, UPLO
!       INTEGER            INFO, LDA, LDC, LWORK, M, N
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
!> SORMTR overwrites the general real M-by-N matrix C with
!>
!>                 SIDE = 'L'     SIDE = 'R'
!> TRANS = 'N':      Q * C          C * Q
!> TRANS = 'T':      Q**T * C       C * Q**T
!>
!> where Q is a real orthogonal matrix of order nq, with nq = m if
!> SIDE = 'L' and nq = n if SIDE = 'R'. Q is defined as the product of
!> nq-1 elementary reflectors, as returned by SSYTRD:
!>
!> if UPLO = 'U', Q = H(nq-1) . . . H(2) H(1);
!>
!> if UPLO = 'L', Q = H(1) H(2) . . . H(nq-1).
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] SIDE
!> \verbatim
!>          SIDE is CHARACTER*1
!>          = 'L': apply Q or Q**T from the Left;
!>          = 'R': apply Q or Q**T from the Right.
!> \endverbatim
!>
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          = 'U': Upper triangle of A contains elementary reflectors
!>                 from SSYTRD;
!>          = 'L': Lower triangle of A contains elementary reflectors
!>                 from SSYTRD.
!> \endverbatim
!>
!> \param[in] TRANS
!> \verbatim
!>          TRANS is CHARACTER*1
!>          = 'N':  No transpose, apply Q;
!>          = 'T':  Transpose, apply Q**T.
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
!> \param[in] A
!> \verbatim
!>          A is REAL array, dimension
!>                               (LDA,M) if SIDE = 'L'
!>                               (LDA,N) if SIDE = 'R'
!>          The vectors which define the elementary reflectors, as
!>          returned by SSYTRD.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.
!>          LDA >= max(1,M) if SIDE = 'L'; LDA >= max(1,N) if SIDE = 'R'.
!> \endverbatim
!>
!> \param[in] TAU
!> \verbatim
!>          TAU is REAL array, dimension
!>                               (M-1) if SIDE = 'L'
!>                               (N-1) if SIDE = 'R'
!>          TAU(i) must contain the scalar factor of the elementary
!>          reflector H(i), as returned by SSYTRD.
!> \endverbatim
!>
!> \param[in,out] C
!> \verbatim
!>          C is REAL array, dimension (LDC,N)
!>          On entry, the M-by-N matrix C.
!>          On exit, C is overwritten by Q*C or Q**T*C or C*Q**T or C*Q.
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
      SUBROUTINE SORMTR(Side,Uplo,Trans,M,N,A,Lda,Tau,C,Ldc,Work,Lwork, &
     &                  Info)
      IMPLICIT NONE
!*--SORMTR176
!
!  -- LAPACK computational routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Side , Trans , Uplo
      INTEGER Info , Lda , Ldc , Lwork , M , N
!     ..
!     .. Array Arguments ..
      REAL A(Lda,*) , C(Ldc,*) , Tau(*) , Work(*)
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
      LOGICAL left , lquery , upper
      INTEGER i1 , i2 , iinfo , lwkopt , mi , ni , nb , nq , nw
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      INTEGER ILAENV
      EXTERNAL ILAENV , LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL SORMQL , SORMQR , XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      Info = 0
      left = LSAME(Side,'L')
      upper = LSAME(Uplo,'U')
      lquery = (Lwork==-1)
!
!     NQ is the order of Q and NW is the minimum dimension of WORK
!
      IF ( left ) THEN
         nq = M
         nw = N
      ELSE
         nq = N
         nw = M
      ENDIF
      IF ( .NOT.left .AND. .NOT.LSAME(Side,'R') ) THEN
         Info = -1
      ELSEIF ( .NOT.upper .AND. .NOT.LSAME(Uplo,'L') ) THEN
         Info = -2
      ELSEIF ( .NOT.LSAME(Trans,'N') .AND. .NOT.LSAME(Trans,'T') ) THEN
         Info = -3
      ELSEIF ( M<0 ) THEN
         Info = -4
      ELSEIF ( N<0 ) THEN
         Info = -5
      ELSEIF ( Lda<MAX(1,nq) ) THEN
         Info = -7
      ELSEIF ( Ldc<MAX(1,M) ) THEN
         Info = -10
      ELSEIF ( Lwork<MAX(1,nw) .AND. .NOT.lquery ) THEN
         Info = -12
      ENDIF
!
      IF ( Info==0 ) THEN
         IF ( upper ) THEN
            IF ( left ) THEN
               nb = ILAENV(1,'SORMQL',Side//Trans,M-1,N,M-1,-1)
            ELSE
               nb = ILAENV(1,'SORMQL',Side//Trans,M,N-1,N-1,-1)
            ENDIF
         ELSEIF ( left ) THEN
            nb = ILAENV(1,'SORMQR',Side//Trans,M-1,N,M-1,-1)
         ELSE
            nb = ILAENV(1,'SORMQR',Side//Trans,M,N-1,N-1,-1)
         ENDIF
         lwkopt = MAX(1,nw)*nb
         Work(1) = lwkopt
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('SORMTR',-Info)
         RETURN
      ELSEIF ( lquery ) THEN
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( M==0 .OR. N==0 .OR. nq==1 ) THEN
         Work(1) = 1
         RETURN
      ENDIF
!
      IF ( left ) THEN
         mi = M - 1
         ni = N
      ELSE
         mi = M
         ni = N - 1
      ENDIF
!
      IF ( upper ) THEN
!
!        Q was determined by a call to SSYTRD with UPLO = 'U'
!
         CALL SORMQL(Side,Trans,mi,ni,nq-1,A(1,2),Lda,Tau,C,Ldc,Work,   &
     &               Lwork,iinfo)
      ELSE
!
!        Q was determined by a call to SSYTRD with UPLO = 'L'
!
         IF ( left ) THEN
            i1 = 2
            i2 = 1
         ELSE
            i1 = 1
            i2 = 2
         ENDIF
         CALL SORMQR(Side,Trans,mi,ni,nq-1,A(2,1),Lda,Tau,C(i1,i2),Ldc, &
     &               Work,Lwork,iinfo)
      ENDIF
      Work(1) = lwkopt
!
!     End of SORMTR
!
      END SUBROUTINE SORMTR
