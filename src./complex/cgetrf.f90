!*==cgetrf.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b CGETRF
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CGETRF + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgetrf.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgetrf.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgetrf.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CGETRF( M, N, A, LDA, IPIV, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, LDA, M, N
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * )
!       COMPLEX            A( LDA, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CGETRF computes an LU factorization of a general M-by-N matrix A
!> using partial pivoting with row interchanges.
!>
!> The factorization has the form
!>    A = P * L * U
!> where P is a permutation matrix, L is lower triangular with unit
!> diagonal elements (lower trapezoidal if m > n), and U is upper
!> triangular (upper trapezoidal if m < n).
!>
!> This is the right-looking Level 3 BLAS version of the algorithm.
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
!>          A is COMPLEX array, dimension (LDA,N)
!>          On entry, the M-by-N matrix to be factored.
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
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value
!>          > 0:  if INFO = i, U(i,i) is exactly zero. The factorization
!>                has been completed, but the factor U is exactly
!>                singular, and division by zero will occur if it is used
!>                to solve a system of equations.
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
!> \ingroup complexGEcomputational
!
!  =====================================================================
      SUBROUTINE CGETRF(M,N,A,Lda,Ipiv,Info)
      USE S_CGEMM
      USE S_CGETRF2
      USE S_CLASWP
      USE S_CTRSM
      USE S_ILAENV
      USE S_XERBLA
      IMPLICIT NONE
!*--CGETRF118
!
! PARAMETER definitions rewritten by SPAG
!
      COMPLEX , PARAMETER  ::  ONE = (1.0E+0,0.0E+0)
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER :: M
      INTEGER :: N
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Ipiv
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      INTEGER :: i , iinfo , j , jb , nb
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
!     .. External Functions ..
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
         CALL XERBLA('CGETRF',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( M==0 .OR. N==0 ) RETURN
!
!     Determine the block size for this environment.
!
      nb = ILAENV(1,'CGETRF',' ',M,N,-1,-1)
      IF ( nb<=1 .OR. nb>=MIN(M,N) ) THEN
!
!        Use unblocked code.
!
         CALL CGETRF2(M,N,A,Lda,Ipiv,Info)
      ELSE
!
!        Use blocked code.
!
         DO j = 1 , MIN(M,N) , nb
            jb = MIN(MIN(M,N)-j+1,nb)
!
!           Factor diagonal and subdiagonal blocks and test for exact
!           singularity.
!
            CALL CGETRF2(M-j+1,jb,A(j,j),Lda,Ipiv(j),iinfo)
!
!           Adjust INFO and the pivot indices.
!
            IF ( Info==0 .AND. iinfo>0 ) Info = iinfo + j - 1
            DO i = j , MIN(M,j+jb-1)
               Ipiv(i) = j - 1 + Ipiv(i)
            ENDDO
!
!           Apply interchanges to columns 1:J-1.
!
            CALL CLASWP(j-1,A,Lda,j,j+jb-1,Ipiv,1)
!
            IF ( j+jb<=N ) THEN
!
!              Apply interchanges to columns J+JB:N.
!
               CALL CLASWP(N-j-jb+1,A(1,j+jb),Lda,j,j+jb-1,Ipiv,1)
!
!              Compute block row of U.
!
               CALL CTRSM('Left','Lower','No transpose','Unit',jb,      &
     &                    N-j-jb+1,ONE,A(j,j),Lda,A(j,j+jb),Lda)
!
!                 Update trailing submatrix.
!
               IF ( j+jb<=M ) CALL CGEMM('No transpose','No transpose', &
     &              M-j-jb+1,N-j-jb+1,jb,-ONE,A(j+jb,j),Lda,A(j,j+jb),  &
     &              Lda,ONE,A(j+jb,j+jb),Lda)
            ENDIF
         ENDDO
      ENDIF
!
!     End of CGETRF
!
      END SUBROUTINE CGETRF
