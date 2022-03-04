!*==cungr2.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b CUNGR2 generates all or part of the unitary matrix Q from an RQ factorization determined by cgerqf (unblocked algorithm).
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CUNGR2 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cungr2.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cungr2.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cungr2.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CUNGR2( M, N, K, A, LDA, TAU, WORK, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, K, LDA, M, N
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
!> CUNGR2 generates an m by n complex matrix Q with orthonormal rows,
!> which is defined as the last m rows of a product of k elementary
!> reflectors of order n
!>
!>       Q  =  H(1)**H H(2)**H . . . H(k)**H
!>
!> as returned by CGERQF.
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
!>          A is COMPLEX array, dimension (LDA,N)
!>          On entry, the (m-k+i)-th row must contain the vector which
!>          defines the elementary reflector H(i), for i = 1,2,...,k, as
!>          returned by CGERQF in the last k rows of its array argument
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
!>          TAU is COMPLEX array, dimension (K)
!>          TAU(i) must contain the scalar factor of the elementary
!>          reflector H(i), as returned by CGERQF.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX array, dimension (M)
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
!> \ingroup complexOTHERcomputational
!
!  =====================================================================
      SUBROUTINE CUNGR2(M,N,K,A,Lda,Tau,Work,Info)
      USE S_CLACGV
      USE S_CLARF
      USE S_CSCAL
      USE S_XERBLA
      IMPLICIT NONE
!*--CUNGR2122
!
! PARAMETER definitions rewritten by SPAG
!
      COMPLEX , PARAMETER  ::  ONE = (1.0E+0,0.0E+0) ,                  &
     &                         ZERO = (0.0E+0,0.0E+0)
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: K
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(*) :: Tau
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      INTEGER :: i , ii , j , l
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
      ELSEIF ( N<M ) THEN
         Info = -2
      ELSEIF ( K<0 .OR. K>M ) THEN
         Info = -3
      ELSEIF ( Lda<MAX(1,M) ) THEN
         Info = -5
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('CUNGR2',-Info)
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
         CALL CLACGV(N-M+ii-1,A(ii,1),Lda)
         A(ii,N-M+ii) = ONE
         CALL CLARF('Right',ii-1,N-M+ii,A(ii,1),Lda,CONJG(Tau(i)),A,Lda,&
     &              Work)
         CALL CSCAL(N-M+ii-1,-Tau(i),A(ii,1),Lda)
         CALL CLACGV(N-M+ii-1,A(ii,1),Lda)
         A(ii,N-M+ii) = ONE - CONJG(Tau(i))
!
!        Set A(m-k+i,n-k+i+1:n) to zero
!
         DO l = N - M + ii + 1 , N
            A(ii,l) = ZERO
         ENDDO
      ENDDO
!
!     End of CUNGR2
!
      END SUBROUTINE CUNGR2
