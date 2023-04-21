!*==zgerq2.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b ZGERQ2 computes the RQ factorization of a general rectangular matrix using an unblocked algorithm.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZGERQ2 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zgerq2.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zgerq2.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zgerq2.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZGERQ2( M, N, A, LDA, TAU, WORK, INFO )
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
!> ZGERQ2 computes an RQ factorization of a complex m by n matrix A:
!> A = R * Q.
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
!>          On exit, if m <= n, the upper triangle of the subarray
!>          A(1:m,n-m+1:n) contains the m by m upper triangular matrix R;
!>          if m >= n, the elements on and above the (m-n)-th subdiagonal
!>          contain the m by n upper trapezoidal matrix R; the remaining
!>          elements, with the array TAU, represent the unitary matrix
!>          Q as a product of elementary reflectors (see Further
!>          Details).
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
!> \date December 2016
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
!>     Q = H(1)**H H(2)**H . . . H(k)**H, where k = min(m,n).
!>
!>  Each H(i) has the form
!>
!>     H(i) = I - tau * v * v**H
!>
!>  where tau is a complex scalar, and v is a complex vector with
!>  v(n-k+i+1:n) = 0 and v(n-k+i) = 1; conjg(v(1:n-k+i-1)) is stored on
!>  exit in A(m-k+i,1:n-k+i-1), and tau in TAU(i).
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE ZGERQ2(M,N,A,Lda,Tau,Work,Info)
      IMPLICIT NONE
!*--ZGERQ2127
!
!  -- LAPACK computational routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
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
         CALL XERBLA('ZGERQ2',-Info)
         RETURN
      ENDIF
!
      k = MIN(M,N)
!
      DO i = k , 1 , -1
!
!        Generate elementary reflector H(i) to annihilate
!        A(m-k+i,1:n-k+i-1)
!
         CALL ZLACGV(N-k+i,A(M-k+i,1),Lda)
         alpha = A(M-k+i,N-k+i)
         CALL ZLARFG(N-k+i,alpha,A(M-k+i,1),Lda,Tau(i))
!
!        Apply H(i) to A(1:m-k+i-1,1:n-k+i) from the right
!
         A(M-k+i,N-k+i) = ONE
         CALL ZLARF('Right',M-k+i-1,N-k+i,A(M-k+i,1),Lda,Tau(i),A,Lda,  &
     &              Work)
         A(M-k+i,N-k+i) = alpha
         CALL ZLACGV(N-k+i-1,A(M-k+i,1),Lda)
      ENDDO
!
!     End of ZGERQ2
!
      END SUBROUTINE ZGERQ2
