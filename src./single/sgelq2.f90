!*==sgelq2.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!
!> \brief \b SGELQ2 computes the LQ factorization of a general rectangular matrix using an unblocked algorithm.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SGELQ2 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgelq2.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgelq2.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgelq2.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SGELQ2( M, N, A, LDA, TAU, WORK, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, LDA, M, N
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
!> SGELQ2 computes an LQ factorization of a real m-by-n matrix A:
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
!>          A is REAL array, dimension (LDA,N)
!>          On entry, the m by n matrix A.
!>          On exit, the elements on and below the diagonal of the array
!>          contain the m by min(m,n) lower trapezoidal matrix L (L is
!>          lower triangular if m <= n); the elements above the diagonal,
!>          with the array TAU, represent the orthogonal matrix Q as a
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
!>          TAU is REAL array, dimension (min(M,N))
!>          The scalar factors of the elementary reflectors (see Further
!>          Details).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is REAL array, dimension (M)
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
!> \ingroup realGEcomputational
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
!>     H(i) = I - tau * v * v**T
!>
!>  where tau is a real scalar, and v is a real vector with
!>  v(1:i-1) = 0 and v(i) = 1; v(i+1:n) is stored on exit in A(i,i+1:n),
!>  and tau in TAU(i).
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE SGELQ2(M,N,A,Lda,Tau,Work,Info)
      USE S_SLARF
      USE S_SLARFG
      USE S_XERBLA
      IMPLICIT NONE
!*--SGELQ2137
!
! PARAMETER definitions rewritten by SPAG
!
      REAL , PARAMETER  ::  ONE = 1.0E+0
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(*) :: Tau
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      REAL :: aii
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
         CALL XERBLA('SGELQ2',-Info)
         RETURN
      ENDIF
!
      k = MIN(M,N)
!
      DO i = 1 , k
!
!        Generate elementary reflector H(i) to annihilate A(i,i+1:n)
!
         CALL SLARFG(N-i+1,A(i,i),A(i,MIN(i+1,N)),Lda,Tau(i))
         IF ( i<M ) THEN
!
!           Apply H(i) to A(i+1:m,i:n) from the right
!
            aii = A(i,i)
            A(i,i) = ONE
            CALL SLARF('Right',M-i,N-i+1,A(i,i),Lda,Tau(i),A(i+1,i),Lda,&
     &                 Work)
            A(i,i) = aii
         ENDIF
      ENDDO
!
!     End of SGELQ2
!
      END SUBROUTINE SGELQ2
