!*==dorm2r.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b DORM2R multiplies a general matrix by the orthogonal matrix from a QR factorization determined by sgeqrf (unblocked algorithm).
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DORM2R + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dorm2r.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dorm2r.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dorm2r.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DORM2R( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC,
!                          WORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          SIDE, TRANS
!       INTEGER            INFO, K, LDA, LDC, M, N
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
!> DORM2R overwrites the general real m by n matrix C with
!>
!>       Q * C  if SIDE = 'L' and TRANS = 'N', or
!>
!>       Q**T* C  if SIDE = 'L' and TRANS = 'T', or
!>
!>       C * Q  if SIDE = 'R' and TRANS = 'N', or
!>
!>       C * Q**T if SIDE = 'R' and TRANS = 'T',
!>
!> where Q is a real orthogonal matrix defined as the product of k
!> elementary reflectors
!>
!>       Q = H(1) H(2) . . . H(k)
!>
!> as returned by DGEQRF. Q is of order m if SIDE = 'L' and of order n
!> if SIDE = 'R'.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] SIDE
!> \verbatim
!>          SIDE is CHARACTER*1
!>          = 'L': apply Q or Q**T from the Left
!>          = 'R': apply Q or Q**T from the Right
!> \endverbatim
!>
!> \param[in] TRANS
!> \verbatim
!>          TRANS is CHARACTER*1
!>          = 'N': apply Q  (No transpose)
!>          = 'T': apply Q**T (Transpose)
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
!>          The number of elementary reflectors whose product defines
!>          the matrix Q.
!>          If SIDE = 'L', M >= K >= 0;
!>          if SIDE = 'R', N >= K >= 0.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (LDA,K)
!>          The i-th column must contain the vector which defines the
!>          elementary reflector H(i), for i = 1,2,...,k, as returned by
!>          DGEQRF in the first k columns of its array argument A.
!>          A is modified by the routine but restored on exit.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.
!>          If SIDE = 'L', LDA >= max(1,M);
!>          if SIDE = 'R', LDA >= max(1,N).
!> \endverbatim
!>
!> \param[in] TAU
!> \verbatim
!>          TAU is DOUBLE PRECISION array, dimension (K)
!>          TAU(i) must contain the scalar factor of the elementary
!>          reflector H(i), as returned by DGEQRF.
!> \endverbatim
!>
!> \param[in,out] C
!> \verbatim
!>          C is DOUBLE PRECISION array, dimension (LDC,N)
!>          On entry, the m by n matrix C.
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
!>          WORK is DOUBLE PRECISION array, dimension
!>                                   (N) if SIDE = 'L',
!>                                   (M) if SIDE = 'R'
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0: successful exit
!>          < 0: if INFO = -i, the i-th argument had an illegal value
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
      SUBROUTINE DORM2R(Side,Trans,M,N,K,A,Lda,Tau,C,Ldc,Work,Info)
      IMPLICIT NONE
!*--DORM2R162
!
!  -- LAPACK computational routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Side , Trans
      INTEGER Info , K , Lda , Ldc , M , N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION A(Lda,*) , C(Ldc,*) , Tau(*) , Work(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ONE
      PARAMETER (ONE=1.0D+0)
!     ..
!     .. Local Scalars ..
      LOGICAL left , notran
      INTEGER i , i1 , i2 , i3 , ic , jc , mi , ni , nq
      DOUBLE PRECISION aii
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL DLARF , XERBLA
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
      notran = LSAME(Trans,'N')
!
!     NQ is the order of Q
!
      IF ( left ) THEN
         nq = M
      ELSE
         nq = N
      ENDIF
      IF ( .NOT.left .AND. .NOT.LSAME(Side,'R') ) THEN
         Info = -1
      ELSEIF ( .NOT.notran .AND. .NOT.LSAME(Trans,'T') ) THEN
         Info = -2
      ELSEIF ( M<0 ) THEN
         Info = -3
      ELSEIF ( N<0 ) THEN
         Info = -4
      ELSEIF ( K<0 .OR. K>nq ) THEN
         Info = -5
      ELSEIF ( Lda<MAX(1,nq) ) THEN
         Info = -7
      ELSEIF ( Ldc<MAX(1,M) ) THEN
         Info = -10
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('DORM2R',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( M==0 .OR. N==0 .OR. K==0 ) RETURN
!
      IF ( (left .AND. .NOT.notran) .OR. (.NOT.left .AND. notran) ) THEN
         i1 = 1
         i2 = K
         i3 = 1
      ELSE
         i1 = K
         i2 = 1
         i3 = -1
      ENDIF
!
      IF ( left ) THEN
         ni = N
         jc = 1
      ELSE
         mi = M
         ic = 1
      ENDIF
!
      DO i = i1 , i2 , i3
         IF ( left ) THEN
!
!           H(i) is applied to C(i:m,1:n)
!
            mi = M - i + 1
            ic = i
         ELSE
!
!           H(i) is applied to C(1:m,i:n)
!
            ni = N - i + 1
            jc = i
         ENDIF
!
!        Apply H(i)
!
         aii = A(i,i)
         A(i,i) = ONE
         CALL DLARF(Side,mi,ni,A(i,i),1,Tau(i),C(ic,jc),Ldc,Work)
         A(i,i) = aii
      ENDDO
!
!     End of DORM2R
!
      END SUBROUTINE DORM2R
