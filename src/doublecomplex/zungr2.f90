!*==zungr2.f90  processed by SPAG 7.51RB at 20:09 on  3 Mar 2022
!> \brief \b ZUNGR2 generates all or part of the unitary matrix Q from an RQ factorization determined by cgerqf (unblocked algorithm).
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZUNGR2 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zungr2.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zungr2.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zungr2.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZUNGR2( M, N, K, A, LDA, TAU, WORK, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, K, LDA, M, N
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
!> ZUNGR2 generates an m by n complex matrix Q with orthonormal rows,
!> which is defined as the last m rows of a product of k elementary
!> reflectors of order n
!>
!>       Q  =  H(1)**H H(2)**H . . . H(k)**H
!>
!> as returned by ZGERQF.
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
!>          The number of columns of the matrix Q. N >= M.
!> \endverbatim
!>
!> \param[in] K
!> \verbatim
!>          K is INTEGER
!>          The number of elementary reflectors whose product defines the
!>          matrix Q. M >= K >= 0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension (LDA,N)
!>          On entry, the (m-k+i)-th row must contain the vector which
!>          defines the elementary reflector H(i), for i = 1,2,...,k, as
!>          returned by ZGERQF in the last k rows of its array argument
!>          A.
!>          On exit, the m-by-n matrix Q.
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
!>          TAU is COMPLEX*16 array, dimension (K)
!>          TAU(i) must contain the scalar factor of the elementary
!>          reflector H(i), as returned by ZGERQF.
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
!> \ingroup complex16OTHERcomputational
!
!  =====================================================================
      SUBROUTINE ZUNGR2(M,N,K,A,Lda,Tau,Work,Info)
      IMPLICIT NONE
!*--ZUNGR2118
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
      COMPLEX*16 A(Lda,*) , Tau(*) , Work(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      COMPLEX*16 ONE , ZERO
      PARAMETER (ONE=(1.0D+0,0.0D+0),ZERO=(0.0D+0,0.0D+0))
!     ..
!     .. Local Scalars ..
      INTEGER i , ii , j , l
!     ..
!     .. External Subroutines ..
      EXTERNAL XERBLA , ZLACGV , ZLARF , ZSCAL
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC DCONJG , MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      Info = 0
      IF ( M<0 ) THEN
         Info = -1
      ELSEIF ( N<M ) THEN
         Info = -2
      ELSEIF ( K<0 .OR. K>M ) THEN
         Info = -3
      ELSEIF ( Lda<MAX(1,M) ) THEN
         Info = -5
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('ZUNGR2',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( M<=0 ) RETURN
!
      IF ( K<M ) THEN
!
!        Initialise rows 1:m-k to rows of the unit matrix
!
         DO j = 1 , N
            DO l = 1 , M - K
               A(l,j) = ZERO
            ENDDO
            IF ( j>N-M .AND. j<=N-K ) A(M-N+j,j) = ONE
         ENDDO
      ENDIF
!
      DO i = 1 , K
         ii = M - K + i
!
!        Apply H(i)**H to A(1:m-k+i,1:n-k+i) from the right
!
         CALL ZLACGV(N-M+ii-1,A(ii,1),Lda)
         A(ii,N-M+ii) = ONE
         CALL ZLARF('Right',ii-1,N-M+ii,A(ii,1),Lda,DCONJG(Tau(i)),A,   &
     &              Lda,Work)
         CALL ZSCAL(N-M+ii-1,-Tau(i),A(ii,1),Lda)
         CALL ZLACGV(N-M+ii-1,A(ii,1),Lda)
         A(ii,N-M+ii) = ONE - DCONJG(Tau(i))
!
!        Set A(m-k+i,n-k+i+1:n) to zero
!
         DO l = N - M + ii + 1 , N
            A(ii,l) = ZERO
         ENDDO
      ENDDO
!
!     End of ZUNGR2
!
      END SUBROUTINE ZUNGR2
