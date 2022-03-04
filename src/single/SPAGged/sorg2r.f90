!*==sorg2r.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b SORG2R generates all or part of the orthogonal matrix Q from a QR factorization determined by sgeqrf (unblocked algorithm).
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SORG2R + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sorg2r.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sorg2r.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sorg2r.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SORG2R( M, N, K, A, LDA, TAU, WORK, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, K, LDA, M, N
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
!> SORG2R generates an m by n real matrix Q with orthonormal columns,
!> which is defined as the first n columns of a product of k elementary
!> reflectors of order m
!>
!>       Q  =  H(1) H(2) . . . H(k)
!>
!> as returned by SGEQRF.
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
!>          A is REAL array, dimension (LDA,N)
!>          On entry, the i-th column must contain the vector which
!>          defines the elementary reflector H(i), for i = 1,2,...,k, as
!>          returned by SGEQRF in the first k columns of its array
!>          argument A.
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
!>          TAU is REAL array, dimension (K)
!>          TAU(i) must contain the scalar factor of the elementary
!>          reflector H(i), as returned by SGEQRF.
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
!> \ingroup realOTHERcomputational
!
!  =====================================================================
      SUBROUTINE SORG2R(M,N,K,A,Lda,Tau,Work,Info)
      USE S_SLARF
      USE S_SSCAL
      USE S_XERBLA
      IMPLICIT NONE
!*--SORG2R121
!
! PARAMETER definitions rewritten by SPAG
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: K
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(*) :: Tau
      REAL , DIMENSION(*) :: Work
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
      ELSEIF ( N<0 .OR. N>M ) THEN
         Info = -2
      ELSEIF ( K<0 .OR. K>N ) THEN
         Info = -3
      ELSEIF ( Lda<MAX(1,M) ) THEN
         Info = -5
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('SORG2R',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( N<=0 ) RETURN
!
!     Initialise columns k+1:n to columns of the unit matrix
!
      DO j = K + 1 , N
         DO l = 1 , M
            A(l,j) = ZERO
         ENDDO
         A(j,j) = ONE
      ENDDO
!
      DO i = K , 1 , -1
!
!        Apply H(i) to A(i:m,i:n) from the left
!
         IF ( i<N ) THEN
            A(i,i) = ONE
            CALL SLARF('Left',M-i+1,N-i,A(i,i),1,Tau(i),A(i,i+1),Lda,   &
     &                 Work)
         ENDIF
         IF ( i<M ) CALL SSCAL(M-i,-Tau(i),A(i+1,i),1)
         A(i,i) = ONE - Tau(i)
!
!        Set A(1:i-1,i) to zero
!
         DO l = 1 , i - 1
            A(l,i) = ZERO
         ENDDO
      ENDDO
!
!     End of SORG2R
!
      END SUBROUTINE SORG2R
