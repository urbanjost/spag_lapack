!*==claghe.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b CLAGHE
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE CLAGHE( N, K, D, A, LDA, ISEED, WORK, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, K, LDA, N
!       ..
!       .. Array Arguments ..
!       INTEGER            ISEED( 4 )
!       REAL               D( * )
!       COMPLEX            A( LDA, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CLAGHE generates a complex hermitian matrix A, by pre- and post-
!> multiplying a real diagonal matrix D with a random unitary matrix:
!> A = U*D*U'. The semi-bandwidth may then be reduced to k by additional
!> unitary transformations.
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
!>          D is REAL array, dimension (N)
!>          The diagonal elements of the diagonal matrix D.
!> \endverbatim
!>
!> \param[out] A
!> \verbatim
!>          A is COMPLEX array, dimension (LDA,N)
!>          The generated n by n hermitian matrix A (the full matrix is
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
!>          WORK is COMPLEX array, dimension (2*N)
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
!> \ingroup complex_matgen
!
!  =====================================================================
      SUBROUTINE CLAGHE(N,K,D,A,Lda,Iseed,Work,Info)
      IMPLICIT NONE
!*--CLAGHE106
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
      REAL D(*)
      COMPLEX A(Lda,*) , Work(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      COMPLEX ZERO , ONE , HALF
      PARAMETER (ZERO=(0.0E+0,0.0E+0),ONE=(1.0E+0,0.0E+0),              &
     &           HALF=(0.5E+0,0.0E+0))
!     ..
!     .. Local Scalars ..
      INTEGER i , j
      REAL wn
      COMPLEX alpha , tau , wa , wb
!     ..
!     .. External Subroutines ..
      EXTERNAL CAXPY , CGEMV , CGERC , CHEMV , CHER2 , CLARNV , CSCAL , &
     &         XERBLA
!     ..
!     .. External Functions ..
      REAL SCNRM2
      COMPLEX CDOTC
      EXTERNAL SCNRM2 , CDOTC
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , CONJG , MAX , REAL
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
         CALL XERBLA('CLAGHE',-Info)
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
!     Generate lower triangle of hermitian matrix
!
      DO i = N - 1 , 1 , -1
!
!        generate random reflection
!
         CALL CLARNV(3,Iseed,N-i+1,Work)
         wn = SCNRM2(N-i+1,Work,1)
         wa = (wn/ABS(Work(1)))*Work(1)
         IF ( wn==ZERO ) THEN
            tau = ZERO
         ELSE
            wb = Work(1) + wa
            CALL CSCAL(N-i,ONE/wb,Work(2),1)
            Work(1) = ONE
            tau = REAL(wb/wa)
         ENDIF
!
!        apply random reflection to A(i:n,i:n) from the left
!        and the right
!
!        compute  y := tau * A * u
!
         CALL CHEMV('Lower',N-i+1,tau,A(i,i),Lda,Work,1,ZERO,Work(N+1), &
     &              1)
!
!        compute  v := y - 1/2 * tau * ( y, u ) * u
!
         alpha = -HALF*tau*CDOTC(N-i+1,Work(N+1),1,Work,1)
         CALL CAXPY(N-i+1,alpha,Work,1,Work(N+1),1)
!
!        apply the transformation as a rank-2 update to A(i:n,i:n)
!
         CALL CHER2('Lower',N-i+1,-ONE,Work,1,Work(N+1),1,A(i,i),Lda)
      ENDDO
!
!     Reduce number of subdiagonals to K
!
      DO i = 1 , N - 1 - K
!
!        generate reflection to annihilate A(k+i+1:n,i)
!
         wn = SCNRM2(N-K-i+1,A(K+i,i),1)
         wa = (wn/ABS(A(K+i,i)))*A(K+i,i)
         IF ( wn==ZERO ) THEN
            tau = ZERO
         ELSE
            wb = A(K+i,i) + wa
            CALL CSCAL(N-K-i,ONE/wb,A(K+i+1,i),1)
            A(K+i,i) = ONE
            tau = REAL(wb/wa)
         ENDIF
!
!        apply reflection to A(k+i:n,i+1:k+i-1) from the left
!
         CALL CGEMV('Conjugate transpose',N-K-i+1,K-1,ONE,A(K+i,i+1),   &
     &              Lda,A(K+i,i),1,ZERO,Work,1)
         CALL CGERC(N-K-i+1,K-1,-tau,A(K+i,i),1,Work,1,A(K+i,i+1),Lda)
!
!        apply reflection to A(k+i:n,k+i:n) from the left and the right
!
!        compute  y := tau * A * u
!
         CALL CHEMV('Lower',N-K-i+1,tau,A(K+i,K+i),Lda,A(K+i,i),1,ZERO, &
     &              Work,1)
!
!        compute  v := y - 1/2 * tau * ( y, u ) * u
!
         alpha = -HALF*tau*CDOTC(N-K-i+1,Work,1,A(K+i,i),1)
         CALL CAXPY(N-K-i+1,alpha,A(K+i,i),1,Work,1)
!
!        apply hermitian rank-2 update to A(k+i:n,k+i:n)
!
         CALL CHER2('Lower',N-K-i+1,-ONE,A(K+i,i),1,Work,1,A(K+i,K+i),  &
     &              Lda)
!
         A(K+i,i) = -wa
         DO j = K + i + 1 , N
            A(j,i) = ZERO
         ENDDO
      ENDDO
!
!     Store full hermitian matrix
!
      DO j = 1 , N
         DO i = j + 1 , N
            A(j,i) = CONJG(A(i,j))
         ENDDO
      ENDDO
!
!     End of CLAGHE
!
      END SUBROUTINE CLAGHE
