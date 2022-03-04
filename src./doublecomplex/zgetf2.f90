!*==zgetf2.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b ZGETF2 computes the LU factorization of a general m-by-n matrix using partial pivoting with row interchanges (unblocked algorithm).
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZGETF2 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zgetf2.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zgetf2.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zgetf2.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZGETF2( M, N, A, LDA, IPIV, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, LDA, M, N
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * )
!       COMPLEX*16         A( LDA, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZGETF2 computes an LU factorization of a general m-by-n matrix A
!> using partial pivoting with row interchanges.
!>
!> The factorization has the form
!>    A = P * L * U
!> where P is a permutation matrix, L is lower triangular with unit
!> diagonal elements (lower trapezoidal if m > n), and U is upper
!> triangular (upper trapezoidal if m < n).
!>
!> This is the right-looking Level 2 BLAS version of the algorithm.
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
!>          A is COMPLEX*16 array, dimension (LDA,N)
!>          On entry, the m by n matrix to be factored.
!>          On exit, the factors L and U from the factorization
!>          A = P*L*U; the unit diagonal elements of L are not stored.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,M).
!> \endverbatim
!>
!> \param[out] IPIV
!> \verbatim
!>          IPIV is INTEGER array, dimension (min(M,N))
!>          The pivot indices; for 1 <= i <= min(M,N), row i of the
!>          matrix was interchanged with row IPIV(i).
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0: successful exit
!>          < 0: if INFO = -k, the k-th argument had an illegal value
!>          > 0: if INFO = k, U(k,k) is exactly zero. The factorization
!>               has been completed, but the factor U is exactly
!>               singular, and division by zero will occur if it is used
!>               to solve a system of equations.
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
!> \ingroup complex16GEcomputational
!
!  =====================================================================
      SUBROUTINE ZGETF2(M,N,A,Lda,Ipiv,Info)
      USE F77KINDS                        
      USE S_DLAMCH
      USE S_IZAMAX
      USE S_XERBLA
      USE S_ZGERU
      USE S_ZSCAL
      USE S_ZSWAP
      IMPLICIT NONE
!*--ZGETF2119
!
! PARAMETER definitions rewritten by SPAG
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ONE = (1.0D+0,0.0D+0) ,        &
     &                 ZERO = (0.0D+0,0.0D+0)
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER , INTENT(IN) :: M
      INTEGER :: N
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(OUT) , DIMENSION(*) :: Ipiv
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      INTEGER :: i , j , jp
      REAL(R8KIND) :: sfmin
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
!     .. External Functions ..
!     ..
!     .. External Subroutines ..
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
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
         CALL XERBLA('ZGETF2',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( M==0 .OR. N==0 ) RETURN
!
!     Compute machine safe minimum
!
      sfmin = DLAMCH('S')
!
      DO j = 1 , MIN(M,N)
!
!        Find pivot and test for singularity.
!
         jp = j - 1 + IZAMAX(M-j+1,A(j,j),1)
         Ipiv(j) = jp
         IF ( A(jp,j)/=ZERO ) THEN
!
!           Apply the interchange to columns 1:N.
!
            IF ( jp/=j ) CALL ZSWAP(N,A(j,1),Lda,A(jp,1),Lda)
!
!           Compute elements J+1:M of J-th column.
!
            IF ( j<M ) THEN
               IF ( ABS(A(j,j))>=sfmin ) THEN
                  CALL ZSCAL(M-j,ONE/A(j,j),A(j+1,j),1)
               ELSE
                  DO i = 1 , M - j
                     A(j+i,j) = A(j+i,j)/A(j,j)
                  ENDDO
               ENDIF
            ENDIF
!
         ELSEIF ( Info==0 ) THEN
!
            Info = j
         ENDIF
!
!
!           Update trailing submatrix.
!
         IF ( j<MIN(M,N) ) CALL ZGERU(M-j,N-j,-ONE,A(j+1,j),1,A(j,j+1), &
     &                                Lda,A(j+1,j+1),Lda)
      ENDDO
!
!     End of ZGETF2
!
      END SUBROUTINE ZGETF2
