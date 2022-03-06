!*==dlagsy.f90  processed by SPAG 7.51RB at 17:44 on  4 Mar 2022
!> \brief \b dlagsy
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE dlagsy( N, K, D, A, LDA, ISEED, WORK, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, K, LDA, N
!       ..
!       .. Array Arguments ..
!       INTEGER            ISEED( 4 )
!       DOUBLE PRECISION   A( LDA, * ), D( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> dlagsy generates a real symmetric matrix A, by pre- and post-
!> multiplying a real diagonal matrix D with a random orthogonal matrix:
!> A = U*D*U'. The semi-bandwidth may then be reduced to k by additional
!> orthogonal transformations.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in] K
!> \verbatim
!>          K is INTEGER
!>          The number of nonzero subdiagonals within the band of A.
!>          0 <= K <= N-1.
!> \endverbatim
!>
!> \param[in] D
!> \verbatim
!>          D is DOUBLE PRECISION array, dimension (N)
!>          The diagonal elements of the diagonal matrix D.
!> \endverbatim
!>
!> \param[out] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (LDA,N)
!>          The generated n by n symmetric matrix A (the full matrix is
!>          stored).
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= N.
!> \endverbatim
!>
!> \param[in,out] ISEED
!> \verbatim
!>          ISEED is INTEGER array, dimension (4)
!>          On entry, the seed of the random number generator; the array
!>          elements must be between 0 and 4095, and ISEED(4) must be
!>          odd.
!>          On exit, the seed is updated.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (2*N)
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
!> \ingroup double_matgen
!
!  =====================================================================
      SUBROUTINE dlagsy(N,K,D,A,Lda,Iseed,Work,Info)
      IMPLICIT NONE
!*--dlagsy105
!
!  -- LAPACK auxiliary routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Info , K , Lda , N
!     ..
!     .. Array Arguments ..
      INTEGER Iseed(4)
      DOUBLE PRECISION A(Lda,*) , D(*) , Work(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO , ONE , HALF
      PARAMETER (ZERO=0.0D+0,ONE=1.0D+0,HALF=0.5D+0)
!     ..
!     .. Local Scalars ..
      INTEGER i , j
      DOUBLE PRECISION alpha , tau , wa , wb , wn
!     ..
!     .. External Subroutines ..
      EXTERNAL DAXPY , DGEMV , DGER , DLARNV , DSCAL , DSYMV , DSYR2 ,  XERBLA
!     ..
!     .. External Functions ..
      DOUBLE PRECISION DDOT , DNRM2
      EXTERNAL DDOT , DNRM2
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX , SIGN
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      Info = 0
      IF ( N<0 ) THEN
         Info = -1
      ELSEIF ( K<0 .OR. K>N-1 ) THEN
         Info = -2
      ELSEIF ( Lda<MAX(1,N) ) THEN
         Info = -5
      ENDIF
      IF ( Info<0 ) THEN
         CALL XERBLA('dlagsy',-Info)
         RETURN
      ENDIF
!
!     initialize lower triangle of A to diagonal matrix
!
      DO j = 1 , N
         DO i = j + 1 , N
            A(i,j) = ZERO
         ENDDO
      ENDDO
      DO i = 1 , N
         A(i,i) = D(i)
      ENDDO
!
!     Generate lower triangle of symmetric matrix
!
      DO i = N - 1 , 1 , -1
!
!        generate random reflection
!
         CALL DLARNV(3,Iseed,N-i+1,Work)
         wn = DNRM2(N-i+1,Work,1)
         wa = SIGN(wn,Work(1))
         IF ( wn==ZERO ) THEN
            tau = ZERO
         ELSE
            wb = Work(1) + wa
            CALL DSCAL(N-i,ONE/wb,Work(2),1)
            Work(1) = ONE
            tau = wb/wa
         ENDIF
!
!        apply random reflection to A(i:n,i:n) from the left
!        and the right
!
!        compute  y := tau * A * u
!
         CALL DSYMV('Lower',N-i+1,tau,A(i,i),Lda,Work,1,ZERO,Work(N+1), &
     &              1)
!
!        compute  v := y - 1/2 * tau * ( y, u ) * u
!
         alpha = -HALF*tau*DDOT(N-i+1,Work(N+1),1,Work,1)
         CALL DAXPY(N-i+1,alpha,Work,1,Work(N+1),1)
!
!        apply the transformation as a rank-2 update to A(i:n,i:n)
!
         CALL DSYR2('Lower',N-i+1,-ONE,Work,1,Work(N+1),1,A(i,i),Lda)
      ENDDO
!
!     Reduce number of subdiagonals to K
!
      DO i = 1 , N - 1 - K
!
!        generate reflection to annihilate A(k+i+1:n,i)
!
         wn = DNRM2(N-K-i+1,A(K+i,i),1)
         wa = SIGN(wn,A(K+i,i))
         IF ( wn==ZERO ) THEN
            tau = ZERO
         ELSE
            wb = A(K+i,i) + wa
            CALL DSCAL(N-K-i,ONE/wb,A(K+i+1,i),1)
            A(K+i,i) = ONE
            tau = wb/wa
         ENDIF
!
!        apply reflection to A(k+i:n,i+1:k+i-1) from the left
!
         CALL DGEMV('Transpose',N-K-i+1,K-1,ONE,A(K+i,i+1),Lda,A(K+i,i),&
     &              1,ZERO,Work,1)
         CALL DGER(N-K-i+1,K-1,-tau,A(K+i,i),1,Work,1,A(K+i,i+1),Lda)
!
!        apply reflection to A(k+i:n,k+i:n) from the left and the right
!
!        compute  y := tau * A * u
!
         CALL DSYMV('Lower',N-K-i+1,tau,A(K+i,K+i),Lda,A(K+i,i),1,ZERO, &
     &              Work,1)
!
!        compute  v := y - 1/2 * tau * ( y, u ) * u
!
         alpha = -HALF*tau*DDOT(N-K-i+1,Work,1,A(K+i,i),1)
         CALL DAXPY(N-K-i+1,alpha,A(K+i,i),1,Work,1)
!
!        apply symmetric rank-2 update to A(k+i:n,k+i:n)
!
         CALL DSYR2('Lower',N-K-i+1,-ONE,A(K+i,i),1,Work,1,A(K+i,K+i),  &
     &              Lda)
!
         A(K+i,i) = -wa
         DO j = K + i + 1 , N
            A(j,i) = ZERO
         ENDDO
      ENDDO
!
!     Store full symmetric matrix
!
      DO j = 1 , N
         DO i = j + 1 , N
            A(j,i) = A(i,j)
         ENDDO
      ENDDO
!
!     End of dlagsy
!
      END SUBROUTINE dlagsy
