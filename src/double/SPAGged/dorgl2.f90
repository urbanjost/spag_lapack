!*==dorgl2.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b DORGL2
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DORGL2 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dorgl2.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dorgl2.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dorgl2.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DORGL2( M, N, K, A, LDA, TAU, WORK, INFO )
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
!> DORGL2 generates an m by n real matrix Q with orthonormal rows,
!> which is defined as the first m rows of a product of k elementary
!> reflectors of order n
!>
!>       Q  =  H(k) . . . H(2) H(1)
!>
!> as returned by DGELQF.
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
!>          A is DOUBLE PRECISION array, dimension (LDA,N)
!>          On entry, the i-th row must contain the vector which defines
!>          the elementary reflector H(i), for i = 1,2,...,k, as returned
!>          by DGELQF in the first k rows of its array argument A.
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
!>          TAU is DOUBLE PRECISION array, dimension (K)
!>          TAU(i) must contain the scalar factor of the elementary
!>          reflector H(i), as returned by DGELQF.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (M)
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
      SUBROUTINE DORGL2(M,N,K,A,Lda,Tau,Work,Info)
      USE F77KINDS                        
      USE S_DLARF
      USE S_DSCAL
      USE S_XERBLA
      IMPLICIT NONE
!*--DORGL2121
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: K
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(*) :: Tau
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      INTEGER :: i , j , l
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
         CALL XERBLA('DORGL2',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( M<=0 ) RETURN
!
      IF ( K<M ) THEN
!
!        Initialise rows k+1:m to rows of the unit matrix
!
         DO j = 1 , N
            DO l = K + 1 , M
               A(l,j) = ZERO
            ENDDO
            IF ( j>K .AND. j<=M ) A(j,j) = ONE
         ENDDO
      ENDIF
!
      DO i = K , 1 , -1
!
!        Apply H(i) to A(i:m,i:n) from the right
!
         IF ( i<N ) THEN
            IF ( i<M ) THEN
               A(i,i) = ONE
               CALL DLARF('Right',M-i,N-i+1,A(i,i),Lda,Tau(i),A(i+1,i), &
     &                    Lda,Work)
            ENDIF
            CALL DSCAL(N-i,-Tau(i),A(i,i+1),Lda)
         ENDIF
         A(i,i) = ONE - Tau(i)
!
!        Set A(i,1:i-1) to zero
!
         DO l = 1 , i - 1
            A(i,l) = ZERO
         ENDDO
      ENDDO
!
!     End of DORGL2
!
      END SUBROUTINE DORGL2
