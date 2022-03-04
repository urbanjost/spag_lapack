!*==sgetf2.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b SGETF2 computes the LU factorization of a general m-by-n matrix using partial pivoting with row interchanges (unblocked algorithm).
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SGETF2 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgetf2.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgetf2.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgetf2.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SGETF2( M, N, A, LDA, IPIV, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, LDA, M, N
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * )
!       REAL               A( LDA, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SGETF2 computes an LU factorization of a general m-by-n matrix A
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
!>          A is REAL array, dimension (LDA,N)
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
!> \ingroup realGEcomputational
!
!  =====================================================================
      SUBROUTINE SGETF2(M,N,A,Lda,Ipiv,Info)
      IMPLICIT NONE
!*--SGETF2112
!
!  -- LAPACK computational routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Info , Lda , M , N
!     ..
!     .. Array Arguments ..
      INTEGER Ipiv(*)
      REAL A(Lda,*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ONE , ZERO
      PARAMETER (ONE=1.0E+0,ZERO=0.0E+0)
!     ..
!     .. Local Scalars ..
      REAL sfmin
      INTEGER i , j , jp
!     ..
!     .. External Functions ..
      REAL SLAMCH
      INTEGER ISAMAX
      EXTERNAL SLAMCH , ISAMAX
!     ..
!     .. External Subroutines ..
      EXTERNAL SGER , SSCAL , SSWAP , XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX , MIN
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
         CALL XERBLA('SGETF2',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( M==0 .OR. N==0 ) RETURN
!
!     Compute machine safe minimum
!
      sfmin = SLAMCH('S')
!
      DO j = 1 , MIN(M,N)
!
!        Find pivot and test for singularity.
!
         jp = j - 1 + ISAMAX(M-j+1,A(j,j),1)
         Ipiv(j) = jp
         IF ( A(jp,j)/=ZERO ) THEN
!
!           Apply the interchange to columns 1:N.
!
            IF ( jp/=j ) CALL SSWAP(N,A(j,1),Lda,A(jp,1),Lda)
!
!           Compute elements J+1:M of J-th column.
!
            IF ( j<M ) THEN
               IF ( ABS(A(j,j))>=sfmin ) THEN
                  CALL SSCAL(M-j,ONE/A(j,j),A(j+1,j),1)
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
         IF ( j<MIN(M,N) ) CALL SGER(M-j,N-j,-ONE,A(j+1,j),1,A(j,j+1),  &
     &                               Lda,A(j+1,j+1),Lda)
      ENDDO
!
!     End of SGETF2
!
      END SUBROUTINE SGETF2
