!*==zlagge.f90  processed by SPAG 7.51RB at 17:44 on  4 Mar 2022
!> \brief \b zlagge
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE zlagge( M, N, KL, KU, D, A, LDA, ISEED, WORK, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, KL, KU, LDA, M, N
!       ..
!       .. Array Arguments ..
!       INTEGER            ISEED( 4 )
!       DOUBLE PRECISION   D( * )
!       COMPLEX*16         A( LDA, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> zlagge generates a complex general m by n matrix A, by pre- and post-
!> multiplying a real diagonal matrix D with random unitary matrices:
!> A = U*D*V. The lower and upper bandwidths may then be reduced to
!> kl and ku by additional unitary transformations.
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
!> \param[in] KL
!> \verbatim
!>          KL is INTEGER
!>          The number of nonzero subdiagonals within the band of A.
!>          0 <= KL <= M-1.
!> \endverbatim
!>
!> \param[in] KU
!> \verbatim
!>          KU is INTEGER
!>          The number of nonzero superdiagonals within the band of A.
!>          0 <= KU <= N-1.
!> \endverbatim
!>
!> \param[in] D
!> \verbatim
!>          D is DOUBLE PRECISION array, dimension (min(M,N))
!>          The diagonal elements of the diagonal matrix D.
!> \endverbatim
!>
!> \param[out] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension (LDA,N)
!>          The generated m by n matrix A.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= M.
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
!>          WORK is COMPLEX*16 array, dimension (M+N)
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
!> \ingroup complex16_matgen
!
!  =====================================================================
      SUBROUTINE zlagge(M,N,Kl,Ku,D,A,Lda,Iseed,Work,Info)
      IMPLICIT NONE
!*--zlagge118
!
!  -- LAPACK auxiliary routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Info , Kl , Ku , Lda , M , N
!     ..
!     .. Array Arguments ..
      INTEGER Iseed(4)
      DOUBLE PRECISION D(*)
      COMPLEX*16 A(Lda,*) , Work(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      COMPLEX*16 ZERO , ONE
      PARAMETER (ZERO=(0.0D+0,0.0D+0),ONE=(1.0D+0,0.0D+0))
!     ..
!     .. Local Scalars ..
      INTEGER i , j
      DOUBLE PRECISION wn
      COMPLEX*16 tau , wa , wb
!     ..
!     .. External Subroutines ..
      EXTERNAL XERBLA , ZGEMV , ZGERC , ZLACGV , ZLARNV , ZSCAL
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , DBLE , MAX , MIN
!     ..
!     .. External Functions ..
      DOUBLE PRECISION DZNRM2
      EXTERNAL DZNRM2
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
      ELSEIF ( Kl<0 .OR. Kl>M-1 ) THEN
         Info = -3
      ELSEIF ( Ku<0 .OR. Ku>N-1 ) THEN
         Info = -4
      ELSEIF ( Lda<MAX(1,M) ) THEN
         Info = -7
      ENDIF
      IF ( Info<0 ) THEN
         CALL XERBLA('zlagge',-Info)
         RETURN
      ENDIF
!
!     initialize A to diagonal matrix
!
      DO j = 1 , N
         DO i = 1 , M
            A(i,j) = ZERO
         ENDDO
      ENDDO
      DO i = 1 , MIN(M,N)
         A(i,i) = D(i)
      ENDDO
!
!     Quick exit if the user wants a diagonal matrix
!
      IF ( (Kl==0) .AND. (Ku==0) ) RETURN
!
!     pre- and post-multiply A by random unitary matrices
!
      DO i = MIN(M,N) , 1 , -1
         IF ( i<M ) THEN
!
!           generate random reflection
!
            CALL ZLARNV(3,Iseed,M-i+1,Work)
            wn = DZNRM2(M-i+1,Work,1)
            wa = (wn/ABS(Work(1)))*Work(1)
            IF ( wn==ZERO ) THEN
               tau = ZERO
            ELSE
               wb = Work(1) + wa
               CALL ZSCAL(M-i,ONE/wb,Work(2),1)
               Work(1) = ONE
               tau = DBLE(wb/wa)
            ENDIF
!
!           multiply A(i:m,i:n) by random reflection from the left
!
            CALL ZGEMV('Conjugate transpose',M-i+1,N-i+1,ONE,A(i,i),Lda,&
     &                 Work,1,ZERO,Work(M+1),1)
            CALL ZGERC(M-i+1,N-i+1,-tau,Work,1,Work(M+1),1,A(i,i),Lda)
         ENDIF
         IF ( i<N ) THEN
!
!           generate random reflection
!
            CALL ZLARNV(3,Iseed,N-i+1,Work)
            wn = DZNRM2(N-i+1,Work,1)
            wa = (wn/ABS(Work(1)))*Work(1)
            IF ( wn==ZERO ) THEN
               tau = ZERO
            ELSE
               wb = Work(1) + wa
               CALL ZSCAL(N-i,ONE/wb,Work(2),1)
               Work(1) = ONE
               tau = DBLE(wb/wa)
            ENDIF
!
!           multiply A(i:m,i:n) by random reflection from the right
!
            CALL ZGEMV('No transpose',M-i+1,N-i+1,ONE,A(i,i),Lda,Work,1,&
     &                 ZERO,Work(N+1),1)
            CALL ZGERC(M-i+1,N-i+1,-tau,Work(N+1),1,Work,1,A(i,i),Lda)
         ENDIF
      ENDDO
!
!     Reduce number of subdiagonals to KL and number of superdiagonals
!     to KU
!
      DO i = 1 , MAX(M-1-Kl,N-1-Ku)
         IF ( Kl<=Ku ) THEN
!
!           annihilate subdiagonal elements first (necessary if KL = 0)
!
            IF ( i<=MIN(M-1-Kl,N) ) THEN
!
!              generate reflection to annihilate A(kl+i+1:m,i)
!
               wn = DZNRM2(M-Kl-i+1,A(Kl+i,i),1)
               wa = (wn/ABS(A(Kl+i,i)))*A(Kl+i,i)
               IF ( wn==ZERO ) THEN
                  tau = ZERO
               ELSE
                  wb = A(Kl+i,i) + wa
                  CALL ZSCAL(M-Kl-i,ONE/wb,A(Kl+i+1,i),1)
                  A(Kl+i,i) = ONE
                  tau = DBLE(wb/wa)
               ENDIF
!
!              apply reflection to A(kl+i:m,i+1:n) from the left
!
               CALL ZGEMV('Conjugate transpose',M-Kl-i+1,N-i,ONE,       &
     &                    A(Kl+i,i+1),Lda,A(Kl+i,i),1,ZERO,Work,1)
               CALL ZGERC(M-Kl-i+1,N-i,-tau,A(Kl+i,i),1,Work,1,         &
     &                    A(Kl+i,i+1),Lda)
               A(Kl+i,i) = -wa
            ENDIF
!
            IF ( i<=MIN(N-1-Ku,M) ) THEN
!
!              generate reflection to annihilate A(i,ku+i+1:n)
!
               wn = DZNRM2(N-Ku-i+1,A(i,Ku+i),Lda)
               wa = (wn/ABS(A(i,Ku+i)))*A(i,Ku+i)
               IF ( wn==ZERO ) THEN
                  tau = ZERO
               ELSE
                  wb = A(i,Ku+i) + wa
                  CALL ZSCAL(N-Ku-i,ONE/wb,A(i,Ku+i+1),Lda)
                  A(i,Ku+i) = ONE
                  tau = DBLE(wb/wa)
               ENDIF
!
!              apply reflection to A(i+1:m,ku+i:n) from the right
!
               CALL ZLACGV(N-Ku-i+1,A(i,Ku+i),Lda)
               CALL ZGEMV('No transpose',M-i,N-Ku-i+1,ONE,A(i+1,Ku+i),  &
     &                    Lda,A(i,Ku+i),Lda,ZERO,Work,1)
               CALL ZGERC(M-i,N-Ku-i+1,-tau,Work,1,A(i,Ku+i),Lda,       &
     &                    A(i+1,Ku+i),Lda)
               A(i,Ku+i) = -wa
            ENDIF
         ELSE
!
!           annihilate superdiagonal elements first (necessary if
!           KU = 0)
!
            IF ( i<=MIN(N-1-Ku,M) ) THEN
!
!              generate reflection to annihilate A(i,ku+i+1:n)
!
               wn = DZNRM2(N-Ku-i+1,A(i,Ku+i),Lda)
               wa = (wn/ABS(A(i,Ku+i)))*A(i,Ku+i)
               IF ( wn==ZERO ) THEN
                  tau = ZERO
               ELSE
                  wb = A(i,Ku+i) + wa
                  CALL ZSCAL(N-Ku-i,ONE/wb,A(i,Ku+i+1),Lda)
                  A(i,Ku+i) = ONE
                  tau = DBLE(wb/wa)
               ENDIF
!
!              apply reflection to A(i+1:m,ku+i:n) from the right
!
               CALL ZLACGV(N-Ku-i+1,A(i,Ku+i),Lda)
               CALL ZGEMV('No transpose',M-i,N-Ku-i+1,ONE,A(i+1,Ku+i),  &
     &                    Lda,A(i,Ku+i),Lda,ZERO,Work,1)
               CALL ZGERC(M-i,N-Ku-i+1,-tau,Work,1,A(i,Ku+i),Lda,       &
     &                    A(i+1,Ku+i),Lda)
               A(i,Ku+i) = -wa
            ENDIF
!
            IF ( i<=MIN(M-1-Kl,N) ) THEN
!
!              generate reflection to annihilate A(kl+i+1:m,i)
!
               wn = DZNRM2(M-Kl-i+1,A(Kl+i,i),1)
               wa = (wn/ABS(A(Kl+i,i)))*A(Kl+i,i)
               IF ( wn==ZERO ) THEN
                  tau = ZERO
               ELSE
                  wb = A(Kl+i,i) + wa
                  CALL ZSCAL(M-Kl-i,ONE/wb,A(Kl+i+1,i),1)
                  A(Kl+i,i) = ONE
                  tau = DBLE(wb/wa)
               ENDIF
!
!              apply reflection to A(kl+i:m,i+1:n) from the left
!
               CALL ZGEMV('Conjugate transpose',M-Kl-i+1,N-i,ONE,       &
     &                    A(Kl+i,i+1),Lda,A(Kl+i,i),1,ZERO,Work,1)
               CALL ZGERC(M-Kl-i+1,N-i,-tau,A(Kl+i,i),1,Work,1,         &
     &                    A(Kl+i,i+1),Lda)
               A(Kl+i,i) = -wa
            ENDIF
         ENDIF
!
         IF ( i<=N ) THEN
            DO j = Kl + i + 1 , M
               A(j,i) = ZERO
            ENDDO
         ENDIF
!
         IF ( i<=M ) THEN
            DO j = Ku + i + 1 , N
               A(i,j) = ZERO
            ENDDO
         ENDIF
      ENDDO
!
!     End of zlagge
!
      END SUBROUTINE zlagge
