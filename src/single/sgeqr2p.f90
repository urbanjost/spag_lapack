!*==sgeqr2p.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b SGEQR2P computes the QR factorization of a general rectangular matrix with non-negative diagonal elements using an unblocked algorithm.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SGEQR2P + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgeqr2p.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgeqr2p.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgeqr2p.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SGEQR2P( M, N, A, LDA, TAU, WORK, INFO )
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
!> SGEQR2P computes a QR factorization of a real m-by-n matrix A:
!>
!>    A = Q * ( R ),
!>            ( 0 )
!>
!> where:
!>
!>    Q is a m-by-m orthogonal matrix;
!>    R is an upper-triangular n-by-n matrix with nonnegative diagonal
!>    entries;
!>    0 is a (m-n)-by-n zero matrix, if m > n.
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
!>          On exit, the elements on and above the diagonal of the array
!>          contain the min(m,n) by n upper trapezoidal matrix R (R is
!>          upper triangular if m >= n). The diagonal entries of R
!>          are nonnegative; the elements below the diagonal,
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
!>          WORK is REAL array, dimension (N)
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
!>     Q = H(1) H(2) . . . H(k), where k = min(m,n).
!>
!>  Each H(i) has the form
!>
!>     H(i) = I - tau * v * v**T
!>
!>  where tau is a real scalar, and v is a real vector with
!>  v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i),
!>  and tau in TAU(i).
!>
!> See Lapack Working Note 203 for details
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE SGEQR2P(M,N,A,Lda,Tau,Work,Info)
      IMPLICIT NONE
!*--SGEQR2P138
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
      REAL A(Lda,*) , Tau(*) , Work(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ONE
      PARAMETER (ONE=1.0E+0)
!     ..
!     .. Local Scalars ..
      INTEGER i , k
      REAL aii
!     ..
!     .. External Subroutines ..
      EXTERNAL SLARF , SLARFGP , XERBLA
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
         CALL XERBLA('SGEQR2P',-Info)
         RETURN
      ENDIF
!
      k = MIN(M,N)
!
      DO i = 1 , k
!
!        Generate elementary reflector H(i) to annihilate A(i+1:m,i)
!
         CALL SLARFGP(M-i+1,A(i,i),A(MIN(i+1,M),i),1,Tau(i))
         IF ( i<N ) THEN
!
!           Apply H(i) to A(i:m,i+1:n) from the left
!
            aii = A(i,i)
            A(i,i) = ONE
            CALL SLARF('Left',M-i+1,N-i,A(i,i),1,Tau(i),A(i,i+1),Lda,   &
     &                 Work)
            A(i,i) = aii
         ENDIF
      ENDDO
!
!     End of SGEQR2P
!
      END SUBROUTINE SGEQR2P
