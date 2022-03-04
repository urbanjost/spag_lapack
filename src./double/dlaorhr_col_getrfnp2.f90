!*==dlaorhr_col_getrfnp2.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b DLAORHR_COL_GETRFNP2
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DLAORHR_GETRF2NP + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaorhr_col_getrfnp2.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaorhr_col_getrfnp2.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaorhr_col_getrfnp2.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       RECURSIVE SUBROUTINE DLAORHR_COL_GETRFNP2( M, N, A, LDA, D, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, LDA, M, N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   A( LDA, * ), D( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLAORHR_COL_GETRFNP2 computes the modified LU factorization without
!> pivoting of a real general M-by-N matrix A. The factorization has
!> the form:
!>
!>     A - S = L * U,
!>
!> where:
!>    S is a m-by-n diagonal sign matrix with the diagonal D, so that
!>    D(i) = S(i,i), 1 <= i <= min(M,N). The diagonal D is constructed
!>    as D(i)=-SIGN(A(i,i)), where A(i,i) is the value after performing
!>    i-1 steps of Gaussian elimination. This means that the diagonal
!>    element at each step of "modified" Gaussian elimination is at
!>    least one in absolute value (so that division-by-zero not
!>    possible during the division by the diagonal element);
!>
!>    L is a M-by-N lower triangular matrix with unit diagonal elements
!>    (lower trapezoidal if M > N);
!>
!>    and U is a M-by-N upper triangular matrix
!>    (upper trapezoidal if M < N).
!>
!> This routine is an auxiliary routine used in the Householder
!> reconstruction routine DORHR_COL. In DORHR_COL, this routine is
!> applied to an M-by-N matrix A with orthonormal columns, where each
!> element is bounded by one in absolute value. With the choice of
!> the matrix S above, one can show that the diagonal element at each
!> step of Gaussian elimination is the largest (in absolute value) in
!> the column on or below the diagonal, so that no pivoting is required
!> for numerical stability [1].
!>
!> For more details on the Householder reconstruction algorithm,
!> including the modified LU factorization, see [1].
!>
!> This is the recursive version of the LU factorization algorithm.
!> Denote A - S by B. The algorithm divides the matrix B into four
!> submatrices:
!>
!>        [  B11 | B12  ]  where B11 is n1 by n1,
!>    B = [ -----|----- ]        B21 is (m-n1) by n1,
!>        [  B21 | B22  ]        B12 is n1 by n2,
!>                               B22 is (m-n1) by n2,
!>                               with n1 = min(m,n)/2, n2 = n-n1.
!>
!>
!> The subroutine calls itself to factor B11, solves for B21,
!> solves for B12, updates B22, then calls itself to factor B22.
!>
!> For more details on the recursive LU algorithm, see [2].
!>
!> DLAORHR_COL_GETRFNP2 is called to factorize a block by the blocked
!> routine DLAORHR_COL_GETRFNP, which uses blocked code calling
!> Level 3 BLAS to update the submatrix. However, DLAORHR_COL_GETRFNP2
!> is self-sufficient and can be used without DLAORHR_COL_GETRFNP.
!>
!> [1] "Reconstructing Householder vectors from tall-skinny QR",
!>     G. Ballard, J. Demmel, L. Grigori, M. Jacquelin, H.D. Nguyen,
!>     E. Solomonik, J. Parallel Distrib. Comput.,
!>     vol. 85, pp. 3-31, 2015.
!>
!> [2] "Recursion leads to automatic variable blocking for dense linear
!>     algebra algorithms", F. Gustavson, IBM J. of Res. and Dev.,
!>     vol. 41, no. 6, pp. 737-755, 1997.
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
!>          A is DOUBLE PRECISION array, dimension (LDA,N)
!>          On entry, the M-by-N matrix to be factored.
!>          On exit, the factors L and U from the factorization
!>          A-S=L*U; the unit diagonal elements of L are not stored.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,M).
!> \endverbatim
!>
!> \param[out] D
!> \verbatim
!>          D is DOUBLE PRECISION array, dimension min(M,N)
!>          The diagonal elements of the diagonal M-by-N sign matrix S,
!>          D(i) = S(i,i), where 1 <= i <= min(M,N). The elements can
!>          be only plus or minus one.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value
!> \endverbatim
!>
!  Authors:
!  ========
!
!> \author Univ. of Tennessee
!> \author Univ. of California Berkeley
!> \author Univ. of Colorado Denver
!> \author NAG Ltd.
!
!> \date November 2019
!
!> \ingroup doubleGEcomputational
!
!> \par Contributors:
!  ==================
!>
!> \verbatim
!>
!> November 2019, Igor Kozachenko,
!>                Computer Science Division,
!>                University of California, Berkeley
!>
!> \endverbatim
!
!  =====================================================================
      RECURSIVE SUBROUTINE DLAORHR_COL_GETRFNP2(M,N,A,Lda,D,Info)
      USE F77KINDS                        
      USE S_DGEMM
      USE S_DLAMCH
      USE S_DSCAL
      USE S_DTRSM
      USE S_XERBLA
      IMPLICIT NONE
!*--DLAORHR_COL_GETRFNP2177
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: D
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      INTEGER :: i , iinfo , n1 , n2
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
!     Test the input parameters
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
         CALL XERBLA('DLAORHR_COL_GETRFNP2',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( MIN(M,N)==0 ) RETURN
 
      IF ( M==1 ) THEN
!
!        One row case, (also recursion termination case),
!        use unblocked code
!
!        Transfer the sign
!
         D(1) = -DSIGN(ONE,A(1,1))
!
!        Construct the row of U
!
         A(1,1) = A(1,1) - D(1)
!
      ELSEIF ( N==1 ) THEN
!
!        One column case, (also recursion termination case),
!        use unblocked code
!
!        Transfer the sign
!
         D(1) = -DSIGN(ONE,A(1,1))
!
!        Construct the row of U
!
         A(1,1) = A(1,1) - D(1)
!
!        Scale the elements 2:M of the column
!
!        Determine machine safe minimum
!
         sfmin = DLAMCH('S')
!
!        Construct the subdiagonal elements of L
!
         IF ( ABS(A(1,1))>=sfmin ) THEN
            CALL DSCAL(M-1,ONE/A(1,1),A(2,1),1)
         ELSE
            DO i = 2 , M
               A(i,1) = A(i,1)/A(1,1)
            ENDDO
         ENDIF
!
      ELSE
!
!        Divide the matrix B into four submatrices
!
         n1 = MIN(M,N)/2
         n2 = N - n1
 
!
!        Factor B11, recursive call
!
         CALL DLAORHR_COL_GETRFNP2(n1,n1,A,Lda,D,iinfo)
!
!        Solve for B21
!
         CALL DTRSM('R','U','N','N',M-n1,n1,ONE,A,Lda,A(n1+1,1),Lda)
!
!        Solve for B12
!
         CALL DTRSM('L','L','N','U',n1,n2,ONE,A,Lda,A(1,n1+1),Lda)
!
!        Update B22, i.e. compute the Schur complement
!        B22 := B22 - B21*B12
!
         CALL DGEMM('N','N',M-n1,n2,n1,-ONE,A(n1+1,1),Lda,A(1,n1+1),Lda,&
     &              ONE,A(n1+1,n1+1),Lda)
!
!        Factor B22, recursive call
!
         CALL DLAORHR_COL_GETRFNP2(M-n1,n2,A(n1+1,n1+1),Lda,D(n1+1),    &
     &                             iinfo)
!
      ENDIF
!
!     End of DLAORHR_COL_GETRFNP2
!
      END SUBROUTINE DLAORHR_COL_GETRFNP2
