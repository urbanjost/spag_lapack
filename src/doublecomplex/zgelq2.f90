!*==zgelq2.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!
!> \brief \b ZGELQ2 computes the LQ factorization of a general rectangular matrix using an unblocked algorithm.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZGELQ2 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zgelq2.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zgelq2.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zgelq2.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZGELQ2( M, N, A, LDA, TAU, WORK, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, LDA, M, N
!       ..
!       .. Array Arguments ..
!       COMPLEX*16         A( LDA, * ), TAU( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZGELQ2 computes an LQ factorization of a complex m-by-n matrix A:
!>
!>    A = ( L 0 ) *  Q
!>
!> where:
!>
!>    Q is a n-by-n orthogonal matrix;
!>    L is a lower-triangular m-by-m matrix;
!>    0 is a m-by-(n-m) zero matrix, if m < n.
!>
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix A.  M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension (LDA,N)
!>          On entry, the m by n matrix A.
!>          On exit, the elements on and below the diagonal of the array
!>          contain the m by min(m,n) lower trapezoidal matrix L (L is
!>          lower triangular if m <= n); the elements above the diagonal,
!>          with the array TAU, represent the unitary matrix Q as a
!>          product of elementary reflectors (see Further Details).
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,M).
!> \endverbatim
!>
!> \param[out] TAU
!> \verbatim
!>          TAU is COMPLEX*16 array, dimension (min(M,N))
!>          The scalar factors of the elementary reflectors (see Further
!>          Details).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX*16 array, dimension (M)
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
!> \date November 2019
!
!> \ingroup complex16GEcomputational
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  The matrix Q is represented as a product of elementary reflectors
!>
!>     Q = H(k)**H . . . H(2)**H H(1)**H, where k = min(m,n).
!>
!>  Each H(i) has the form
!>
!>     H(i) = I - tau * v * v**H
!>
!>  where tau is a complex scalar, and v is a complex vector with
!>  v(1:i-1) = 0 and v(i) = 1; conjg(v(i+1:n)) is stored on exit in
!>  A(i,i+1:n), and tau in TAU(i).
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE ZGELQ2(M,N,A,Lda,Tau,Work,Info)
      IMPLICIT NONE
!*--ZGELQ2134
!
!  -- LAPACK computational routine (version 3.9.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2019
!
!     .. Scalar Arguments ..
      INTEGER Info , Lda , M , N
!     ..
!     .. Array Arguments ..
      COMPLEX*16 A(Lda,*) , Tau(*) , Work(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      COMPLEX*16 ONE
      PARAMETER (ONE=(1.0D+0,0.0D+0))
!     ..
!     .. Local Scalars ..
      INTEGER i , k
      COMPLEX*16 alpha
!     ..
!     .. External Subroutines ..
      EXTERNAL XERBLA , ZLACGV , ZLARF , ZLARFG
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX , MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      Info = 0
      IF ( M<0 ) THEN
         Info = -1
      ELSEIF ( N<0 ) THEN
         Info = -2
      ELSEIF ( Lda<MAX(1,M) ) THEN
         Info = -4
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('ZGELQ2',-Info)
         RETURN
      ENDIF
!
      k = MIN(M,N)
!
      DO i = 1 , k
!
!        Generate elementary reflector H(i) to annihilate A(i,i+1:n)
!
         CALL ZLACGV(N-i+1,A(i,i),Lda)
         alpha = A(i,i)
         CALL ZLARFG(N-i+1,alpha,A(i,MIN(i+1,N)),Lda,Tau(i))
         IF ( i<M ) THEN
!
!           Apply H(i) to A(i+1:m,i:n) from the right
!
            A(i,i) = ONE
            CALL ZLARF('Right',M-i,N-i+1,A(i,i),Lda,Tau(i),A(i+1,i),Lda,&
     &                 Work)
         ENDIF
         A(i,i) = alpha
         CALL ZLACGV(N-i+1,A(i,i),Lda)
      ENDDO
!
!     End of ZGELQ2
!
      END SUBROUTINE ZGELQ2
