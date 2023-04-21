!*==dorg2l.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
 
!> \brief \b DORG2L generates all or part of the orthogonal matrix Q from a QL factorization determined by sgeqlf (unblocked algorithm).
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DORG2L + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dorg2l.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dorg2l.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dorg2l.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DORG2L( M, N, K, A, LDA, TAU, WORK, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, K, LDA, M, N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DORG2L generates an m by n real matrix Q with orthonormal columns,
!> which is defined as the last n columns of a product of k elementary
!> reflectors of order m
!>
!>       Q  =  H(k) . . . H(2) H(1)
!>
!> as returned by DGEQLF.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix Q. M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix Q. M >= N >= 0.
!> \endverbatim
!>
!> \param[in] K
!> \verbatim
!>          K is INTEGER
!>          The number of elementary reflectors whose product defines the
!>          matrix Q. N >= K >= 0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (LDA,N)
!>          On entry, the (n-k+i)-th column must contain the vector which
!>          defines the elementary reflector H(i), for i = 1,2,...,k, as
!>          returned by DGEQLF in the last k columns of its array
!>          argument A.
!>          On exit, the m by n matrix Q.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The first dimension of the array A. LDA >= max(1,M).
!> \endverbatim
!>
!> \param[in] TAU
!> \verbatim
!>          TAU is DOUBLE PRECISION array, dimension (K)
!>          TAU(i) must contain the scalar factor of the elementary
!>          reflector H(i), as returned by DGEQLF.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (N)
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0: successful exit
!>          < 0: if INFO = -i, the i-th argument has an illegal value
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
      SUBROUTINE DORG2L(M,N,K,A,Lda,Tau,Work,Info)
      IMPLICIT NONE
!*--DORG2L119
!
!  -- LAPACK computational routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Info , K , Lda , M , N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION A(Lda,*) , Tau(*) , Work(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ONE , ZERO
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
!     ..
!     .. Local Scalars ..
      INTEGER i , ii , j , l
!     ..
!     .. External Subroutines ..
      EXTERNAL DLARF , DSCAL , XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      Info = 0
      IF ( M<0 ) THEN
         Info = -1
      ELSEIF ( N<0 .OR. N>M ) THEN
         Info = -2
      ELSEIF ( K<0 .OR. K>N ) THEN
         Info = -3
      ELSEIF ( Lda<MAX(1,M) ) THEN
         Info = -5
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('DORG2L',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( N<=0 ) RETURN
!
!     Initialise columns 1:n-k to columns of the unit matrix
!
      DO j = 1 , N - K
         DO l = 1 , M
            A(l,j) = ZERO
         ENDDO
         A(M-N+j,j) = ONE
      ENDDO
!
      DO i = 1 , K
         ii = N - K + i
!
!        Apply H(i) to A(1:m-k+i,1:n-k+i) from the left
!
         A(M-N+ii,ii) = ONE
         CALL DLARF('Left',M-N+ii,ii-1,A(1,ii),1,Tau(i),A,Lda,Work)
         CALL DSCAL(M-N+ii-1,-Tau(i),A(1,ii),1)
         A(M-N+ii,ii) = ONE - Tau(i)
!
!        Set A(m-k+i+1:m,n-k+i) to zero
!
         DO l = M - N + ii + 1 , M
            A(l,ii) = ZERO
         ENDDO
      ENDDO
!
!     End of DORG2L
!
      END SUBROUTINE DORG2L
