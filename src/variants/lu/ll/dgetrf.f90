!*==dgetrf.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b DGETRF VARIANT: left-looking Level 3 BLAS version of the algorithm.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DGETRF ( M, N, A, LDA, IPIV, INFO)
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, LDA, M, N
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * )
!       DOUBLE PRECISION   A( LDA, * )
!       ..
!
!  Purpose
!  =======
!
!>\details \b Purpose:
!>\verbatim
!>
!> DGETRF computes an LU factorization of a general M-by-N matrix A
!> using partial pivoting with row interchanges.
!>
!> The factorization has the form
!>    A = P * L * U
!> where P is a permutation matrix, L is lower triangular with unit
!> diagonal elements (lower trapezoidal if m > n), and U is upper
!> triangular (upper trapezoidal if m < n).
!>
!> This is the left-looking Level 3 BLAS version of the algorithm.
!>
!>\endverbatim
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
!>          A is DOUBLE PRECISION array, dimension (LDA,N)
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
!>
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
!> \ingroup variantsGEcomputational
!
!  =====================================================================
      SUBROUTINE DGETRF(M,N,A,Lda,Ipiv,Info)
      IMPLICIT NONE
!*--DGETRF104
!
!  -- LAPACK computational routine (version 3.1) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Info , Lda , M , N
!     ..
!     .. Array Arguments ..
      INTEGER Ipiv(*)
      DOUBLE PRECISION A(Lda,*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ONE
      PARAMETER (ONE=1.0D+0)
!     ..
!     .. Local Scalars ..
      INTEGER i , iinfo , j , jb , k , nb
!     ..
!     .. External Subroutines ..
      EXTERNAL DGEMM , DGETF2 , DLASWP , DTRSM , XERBLA
!     ..
!     .. External Functions ..
      INTEGER ILAENV
      EXTERNAL ILAENV
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
         CALL XERBLA('DGETRF',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( M==0 .OR. N==0 ) RETURN
!
!     Determine the block size for this environment.
!
      nb = ILAENV(1,'DGETRF',' ',M,N,-1,-1)
      IF ( nb<=1 .OR. nb>=MIN(M,N) ) THEN
!
!        Use unblocked code.
!
         CALL DGETF2(M,N,A,Lda,Ipiv,Info)
 
      ELSE
!
!        Use blocked code.
!
         DO j = 1 , MIN(M,N) , nb
            jb = MIN(MIN(M,N)-j+1,nb)
!
!           Update before factoring the current panel
!
            DO k = 1 , j - nb , nb
!
!              Apply interchanges to rows K:K+NB-1.
!
               CALL DLASWP(jb,A(1,j),Lda,k,k+nb-1,Ipiv,1)
!
!              Compute block row of U.
!
               CALL DTRSM('Left','Lower','No transpose','Unit',nb,jb,   &
     &                    ONE,A(k,k),Lda,A(k,j),Lda)
!
!              Update trailing submatrix.
!
               CALL DGEMM('No transpose','No transpose',M-k-nb+1,jb,nb, &
     &                    -ONE,A(k+nb,k),Lda,A(k,j),Lda,ONE,A(k+nb,j),  &
     &                    Lda)
            ENDDO
!
!           Factor diagonal and subdiagonal blocks and test for exact
!           singularity.
!
            CALL DGETF2(M-j+1,jb,A(j,j),Lda,Ipiv(j),iinfo)
!
!           Adjust INFO and the pivot indices.
!
            IF ( Info==0 .AND. iinfo>0 ) Info = iinfo + j - 1
            DO i = j , MIN(M,j+jb-1)
               Ipiv(i) = j - 1 + Ipiv(i)
            ENDDO
!
         ENDDO
 
!
!        Apply interchanges to the left-overs
!
         DO k = 1 , MIN(M,N) , nb
            CALL DLASWP(k-1,A(1,1),Lda,k,MIN(k+nb-1,MIN(M,N)),Ipiv,1)
         ENDDO
!
!        Apply update to the M+1:N columns when N > M
!
         IF ( N>M ) THEN
 
            CALL DLASWP(N-M,A(1,M+1),Lda,1,M,Ipiv,1)
 
            DO k = 1 , M , nb
 
               jb = MIN(M-k+1,nb)
!
               CALL DTRSM('Left','Lower','No transpose','Unit',jb,N-M,  &
     &                    ONE,A(k,k),Lda,A(k,M+1),Lda)
 
!
               IF ( k+nb<=M ) CALL DGEMM('No transpose','No transpose', &
     &              M-k-nb+1,N-M,nb,-ONE,A(k+nb,k),Lda,A(k,M+1),Lda,ONE,&
     &              A(k+nb,M+1),Lda)
            ENDDO
         ENDIF
!
      ENDIF
!
!     End of DGETRF
!
      END SUBROUTINE DGETRF
