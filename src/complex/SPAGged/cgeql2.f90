!*==cgeql2.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b CGEQL2 computes the QL factorization of a general rectangular matrix using an unblocked algorithm.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CGEQL2 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgeql2.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgeql2.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgeql2.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CGEQL2( M, N, A, LDA, TAU, WORK, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, LDA, M, N
!       ..
!       .. Array Arguments ..
!       COMPLEX            A( LDA, * ), TAU( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CGEQL2 computes a QL factorization of a complex m by n matrix A:
!> A = Q * L.
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
!>          A is COMPLEX array, dimension (LDA,N)
!>          On entry, the m by n matrix A.
!>          On exit, if m >= n, the lower triangle of the subarray
!>          A(m-n+1:m,1:n) contains the n by n lower triangular matrix L;
!>          if m <= n, the elements on and below the (n-m)-th
!>          superdiagonal contain the m by n lower trapezoidal matrix L;
!>          the remaining elements, with the array TAU, represent the
!>          unitary matrix Q as a product of elementary reflectors
!>          (see Further Details).
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
!>          TAU is COMPLEX array, dimension (min(M,N))
!>          The scalar factors of the elementary reflectors (see Further
!>          Details).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX array, dimension (N)
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
!> \ingroup complexGEcomputational
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  The matrix Q is represented as a product of elementary reflectors
!>
!>     Q = H(k) . . . H(2) H(1), where k = min(m,n).
!>
!>  Each H(i) has the form
!>
!>     H(i) = I - tau * v * v**H
!>
!>  where tau is a complex scalar, and v is a complex vector with
!>  v(m-k+i+1:m) = 0 and v(m-k+i) = 1; v(1:m-k+i-1) is stored on exit in
!>  A(1:m-k+i-1,n-k+i), and tau in TAU(i).
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE CGEQL2(M,N,A,Lda,Tau,Work,Info)
      USE S_CLARF
      USE S_CLARFG
      USE S_XERBLA
      IMPLICIT NONE
!*--CGEQL2130
!
! PARAMETER definitions rewritten by SPAG
!
      COMPLEX , PARAMETER  ::  ONE = (1.0E+0,0.0E+0)
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(*) :: Tau
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      COMPLEX :: alpha
      INTEGER :: i , k
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
!     .. External Subroutines ..
!     ..
!     .. Intrinsic Functions ..
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
         CALL XERBLA('CGEQL2',-Info)
         RETURN
      ENDIF
!
      k = MIN(M,N)
!
      DO i = k , 1 , -1
!
!        Generate elementary reflector H(i) to annihilate
!        A(1:m-k+i-1,n-k+i)
!
         alpha = A(M-k+i,N-k+i)
         CALL CLARFG(M-k+i,alpha,A(1,N-k+i),1,Tau(i))
!
!        Apply H(i)**H to A(1:m-k+i,1:n-k+i-1) from the left
!
         A(M-k+i,N-k+i) = ONE
         CALL CLARF('Left',M-k+i,N-k+i-1,A(1,N-k+i),1,CONJG(Tau(i)),A,  &
     &              Lda,Work)
         A(M-k+i,N-k+i) = alpha
      ENDDO
!
!     End of CGEQL2
!
      END SUBROUTINE CGEQL2
