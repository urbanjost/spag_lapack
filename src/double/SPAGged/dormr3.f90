!*==dormr3.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b DORMR3 multiplies a general matrix by the orthogonal matrix from a RZ factorization determined by stzrzf (unblocked algorithm).
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DORMR3 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dormr3.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dormr3.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dormr3.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DORMR3( SIDE, TRANS, M, N, K, L, A, LDA, TAU, C, LDC,
!                          WORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          SIDE, TRANS
!       INTEGER            INFO, K, L, LDA, LDC, M, N
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
!> DORMR3 overwrites the general real m by n matrix C with
!>
!>       Q * C  if SIDE = 'L' and TRANS = 'N', or
!>
!>       Q**T* C  if SIDE = 'L' and TRANS = 'C', or
!>
!>       C * Q  if SIDE = 'R' and TRANS = 'N', or
!>
!>       C * Q**T if SIDE = 'R' and TRANS = 'C',
!>
!> where Q is a real orthogonal matrix defined as the product of k
!> elementary reflectors
!>
!>       Q = H(1) H(2) . . . H(k)
!>
!> as returned by DTZRZF. Q is of order m if SIDE = 'L' and of order n
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
!> \param[in] L
!> \verbatim
!>          L is INTEGER
!>          The number of columns of the matrix A containing
!>          the meaningful part of the Householder reflectors.
!>          If SIDE = 'L', M >= L >= 0, if SIDE = 'R', N >= L >= 0.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension
!>                               (LDA,M) if SIDE = 'L',
!>                               (LDA,N) if SIDE = 'R'
!>          The i-th row must contain the vector which defines the
!>          elementary reflector H(i), for i = 1,2,...,k, as returned by
!>          DTZRZF in the last k rows of its array argument A.
!>          A is modified by the routine but restored on exit.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A. LDA >= max(1,K).
!> \endverbatim
!>
!> \param[in] TAU
!> \verbatim
!>          TAU is DOUBLE PRECISION array, dimension (K)
!>          TAU(i) must contain the scalar factor of the elementary
!>          reflector H(i), as returned by DTZRZF.
!> \endverbatim
!>
!> \param[in,out] C
!> \verbatim
!>          C is DOUBLE PRECISION array, dimension (LDC,N)
!>          On entry, the m-by-n matrix C.
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
!> \par Contributors:
!  ==================
!>
!>    A. Petitet, Computer Science Dept., Univ. of Tenn., Knoxville, USA
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE DORMR3(Side,Trans,M,N,K,L,A,Lda,Tau,C,Ldc,Work,Info)
      USE F77KINDS                        
      USE S_DLARZ
      USE S_LSAME
      USE S_XERBLA
      IMPLICIT NONE
!*--DORMR3185
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Side
      CHARACTER :: Trans
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: K
      INTEGER :: L
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(*) :: Tau
      REAL(R8KIND) , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      INTEGER :: i , i1 , i2 , i3 , ic , ja , jc , mi , ni , nq
      LOGICAL :: left , notran
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
      ELSEIF ( L<0 .OR. (left .AND. (L>M)) .OR. (.NOT.left .AND. (L>N)) &
     &         ) THEN
         Info = -6
      ELSEIF ( Lda<MAX(1,K) ) THEN
         Info = -8
      ELSEIF ( Ldc<MAX(1,M) ) THEN
         Info = -11
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('DORMR3',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( M==0 .OR. N==0 .OR. K==0 ) RETURN
!
      IF ( (left .AND. .NOT.notran .OR. .NOT.left .AND. notran) ) THEN
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
         ja = M - L + 1
         jc = 1
      ELSE
         mi = M
         ja = N - L + 1
         ic = 1
      ENDIF
!
      DO i = i1 , i2 , i3
         IF ( left ) THEN
!
!           H(i) or H(i)**T is applied to C(i:m,1:n)
!
            mi = M - i + 1
            ic = i
         ELSE
!
!           H(i) or H(i)**T is applied to C(1:m,i:n)
!
            ni = N - i + 1
            jc = i
         ENDIF
!
!        Apply H(i) or H(i)**T
!
         CALL DLARZ(Side,mi,ni,L,A(i,ja),Lda,Tau(i),C(ic,jc),Ldc,Work)
!
      ENDDO
!
!
!     End of DORMR3
!
      END SUBROUTINE DORMR3
