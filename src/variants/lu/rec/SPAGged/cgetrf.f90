!*==cgetrf.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b CGETRF VARIANT: iterative version of Sivan Toledo's recursive LU algorithm
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
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
!  Purpose
!  =======
!
!>\details \b Purpose:
!>\verbatim
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
!> This code implements an iterative version of Sivan Toledo's recursive
!> LU algorithm[1].  For square matrices, this iterative versions should
!> be within a factor of two of the optimum number of memory transfers.
!>
!> The pattern is as follows, with the large blocks of U being updated
!> in one call to DTRSM, and the dotted lines denoting sections that
!> have had all pending permutations applied:
!>
!>  1 2 3 4 5 6 7 8
!> +-+-+---+-------+------
!> | |1|   |       |
!> |.+-+ 2 |       |
!> | | |   |       |
!> |.|.+-+-+   4   |
!> | | | |1|       |
!> | | |.+-+       |
!> | | | | |       |
!> |.|.|.|.+-+-+---+  8
!> | | | | | |1|   |
!> | | | | |.+-+ 2 |
!> | | | | | | |   |
!> | | | | |.|.+-+-+
!> | | | | | | | |1|
!> | | | | | | |.+-+
!> | | | | | | | | |
!> |.|.|.|.|.|.|.|.+-----
!> | | | | | | | | |
!>
!> The 1-2-1-4-1-2-1-8-... pattern is the position of the last 1 bit in
!> the binary expansion of the current column.  Each Schur update is
!> applied as soon as the necessary portion of U is available.
!>
!> [1] Toledo, S. 1997. Locality of Reference in LU Decomposition with
!> Partial Pivoting. SIAM J. Matrix Anal. Appl. 18, 4 (Oct. 1997),
!> 1065-1081. http://dx.doi.org/10.1137/S0895479896297744
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
      SUBROUTINE CGETRF(M,N,A,Lda,Ipiv,Info)
      USE S_CGEMM
      USE S_CLASWP
      USE S_CSCAL
      USE S_CTRSM
      USE S_ICAMAX
      USE S_SISNAN
      USE S_SLAMCH
      USE S_XERBLA
      IMPLICIT NONE
!*--CGETRF146
!
! PARAMETER definitions rewritten by SPAG
!
      COMPLEX , PARAMETER  ::  ONE = (1.0E+0,0.0E+0) ,                  &
     &                         NEGONE = (-1.0E+0,0.0E+0)
      REAL , PARAMETER  ::  ZERO = 0.0E+0
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER :: M
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , DIMENSION(*) :: Ipiv
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      INTEGER :: i , ipivstart , j , jp , jpivstart , kahead , kcols ,  &
     &           kstart , npived , nstep , ntopiv
      REAL :: pivmag , sfmin
      COMPLEX :: tmp
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
         CALL XERBLA('CGETRF',-Info)
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
      nstep = MIN(M,N)
      DO j = 1 , nstep
         kahead = IAND(j,-j)
         kstart = j + 1 - kahead
         kcols = MIN(kahead,M-j)
!
!        Find pivot.
!
         jp = j - 1 + ICAMAX(M-j+1,A(j,j),1)
         Ipiv(j) = jp
 
!        Permute just this column.
         IF ( jp/=j ) THEN
            tmp = A(j,j)
            A(j,j) = A(jp,j)
            A(jp,j) = tmp
         ENDIF
 
!        Apply pending permutations to L
         ntopiv = 1
         ipivstart = j
         jpivstart = j - ntopiv
         DO WHILE ( ntopiv<kahead )
            CALL CLASWP(ntopiv,A(1,jpivstart),Lda,ipivstart,j,Ipiv,1)
            ipivstart = ipivstart - ntopiv
            ntopiv = ntopiv*2
            jpivstart = jpivstart - ntopiv
         ENDDO
 
!        Permute U block to match L
         CALL CLASWP(kcols,A(1,j+1),Lda,kstart,j,Ipiv,1)
 
!        Factor the current column
         pivmag = ABS(A(j,j))
         IF ( pivmag/=ZERO .AND. .NOT.SISNAN(pivmag) ) THEN
            IF ( pivmag>=sfmin ) THEN
               CALL CSCAL(M-j,ONE/A(j,j),A(j+1,j),1)
            ELSE
               DO i = 1 , M - j
                  A(j+i,j) = A(j+i,j)/A(j,j)
               ENDDO
            ENDIF
         ELSEIF ( pivmag==ZERO .AND. Info==0 ) THEN
            Info = j
         ENDIF
 
!        Solve for U block.
         CALL CTRSM('Left','Lower','No transpose','Unit',kahead,kcols,  &
     &              ONE,A(kstart,kstart),Lda,A(kstart,j+1),Lda)
!        Schur complement.
         CALL CGEMM('No transpose','No transpose',M-j,kcols,kahead,     &
     &              NEGONE,A(j+1,kstart),Lda,A(kstart,j+1),Lda,ONE,     &
     &              A(j+1,j+1),Lda)
      ENDDO
 
!     Handle pivot permutations on the way out of the recursion
      npived = IAND(nstep,-nstep)
      j = nstep - npived
      DO WHILE ( j>0 )
         ntopiv = IAND(j,-j)
         CALL CLASWP(ntopiv,A(1,j-ntopiv+1),Lda,j+1,nstep,Ipiv,1)
         j = j - ntopiv
      ENDDO
 
!     If short and wide, handle the rest of the columns.
      IF ( M<N ) THEN
         CALL CLASWP(N-M,A(1,M+kcols+1),Lda,1,M,Ipiv,1)
         CALL CTRSM('Left','Lower','No transpose','Unit',M,N-M,ONE,A,   &
     &              Lda,A(1,M+kcols+1),Lda)
      ENDIF
 
!
!     End of CGETRF
!
      END SUBROUTINE CGETRF
