!*==slarge.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b SLARGE
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE SLARGE( N, A, LDA, ISEED, WORK, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, LDA, N
!       ..
!       .. Array Arguments ..
!       INTEGER            ISEED( 4 )
!       REAL               A( LDA, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SLARGE pre- and post-multiplies a real general n by n matrix A
!> with a random orthogonal matrix: A = U*D*U'.
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
!> \param[in,out] A
!> \verbatim
!>          A is REAL array, dimension (LDA,N)
!>          On entry, the original n by n matrix A.
!>          On exit, A is overwritten by U*A*U' for some random
!>          orthogonal matrix U.
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
!>          WORK is REAL array, dimension (2*N)
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
!> \ingroup real_matgen
!
!  =====================================================================
      SUBROUTINE SLARGE(N,A,Lda,Iseed,Work,Info)
      IMPLICIT NONE
!*--SLARGE91
!
!  -- LAPACK auxiliary routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Info , Lda , N
!     ..
!     .. Array Arguments ..
      INTEGER Iseed(4)
      REAL A(Lda,*) , Work(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ZERO , ONE
      PARAMETER (ZERO=0.0E+0,ONE=1.0E+0)
!     ..
!     .. Local Scalars ..
      INTEGER i
      REAL tau , wa , wb , wn
!     ..
!     .. External Subroutines ..
      EXTERNAL SGEMV , SGER , SLARNV , SSCAL , XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX , SIGN
!     ..
!     .. External Functions ..
      REAL SNRM2
      EXTERNAL SNRM2
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      Info = 0
      IF ( N<0 ) THEN
         Info = -1
      ELSEIF ( Lda<MAX(1,N) ) THEN
         Info = -3
      ENDIF
      IF ( Info<0 ) THEN
         CALL XERBLA('SLARGE',-Info)
         RETURN
      ENDIF
!
!     pre- and post-multiply A by random orthogonal matrix
!
      DO i = N , 1 , -1
!
!        generate random reflection
!
         CALL SLARNV(3,Iseed,N-i+1,Work)
         wn = SNRM2(N-i+1,Work,1)
         wa = SIGN(wn,Work(1))
         IF ( wn==ZERO ) THEN
            tau = ZERO
         ELSE
            wb = Work(1) + wa
            CALL SSCAL(N-i,ONE/wb,Work(2),1)
            Work(1) = ONE
            tau = wb/wa
         ENDIF
!
!        multiply A(i:n,1:n) by random reflection from the left
!
         CALL SGEMV('Transpose',N-i+1,N,ONE,A(i,1),Lda,Work,1,ZERO,     &
     &              Work(N+1),1)
         CALL SGER(N-i+1,N,-tau,Work,1,Work(N+1),1,A(i,1),Lda)
!
!        multiply A(1:n,i:n) by random reflection from the right
!
         CALL SGEMV('No transpose',N,N-i+1,ONE,A(1,i),Lda,Work,1,ZERO,  &
     &              Work(N+1),1)
         CALL SGER(N,N-i+1,-tau,Work(N+1),1,Work,1,A(1,i),Lda)
      ENDDO
!
!     End of SLARGE
!
      END SUBROUTINE SLARGE
