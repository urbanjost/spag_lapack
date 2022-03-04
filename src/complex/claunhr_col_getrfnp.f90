!*==claunhr_col_getrfnp.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b CLAUNHR_COL_GETRFNP
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CLAUNHR_COL_GETRFNP + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/claunhr_col_getrfnp.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/claunhr_col_getrfnp.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/claunhr_col_getrfnp.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CLAUNHR_COL_GETRFNP( M, N, A, LDA, D, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, LDA, M, N
!       ..
!       .. Array Arguments ..
!       COMPLEX            A( LDA, * ), D( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CLAUNHR_COL_GETRFNP computes the modified LU factorization without
!> pivoting of a complex general M-by-N matrix A. The factorization has
!> the form:
!>
!>     A - S = L * U,
!>
!> where:
!>    S is a m-by-n diagonal sign matrix with the diagonal D, so that
!>    D(i) = S(i,i), 1 <= i <= min(M,N). The diagonal D is constructed
!>    as D(i)=-SIGN(A(i,i)), where A(i,i) is the value after performing
!>    i-1 steps of Gaussian elimination. This means that the diagonal
!>    element at each step of "modified" Gaussian elimination is
!>    at least one in absolute value (so that division-by-zero not
!>    not possible during the division by the diagonal element);
!>
!>    L is a M-by-N lower triangular matrix with unit diagonal elements
!>    (lower trapezoidal if M > N);
!>
!>    and U is a M-by-N upper triangular matrix
!>    (upper trapezoidal if M < N).
!>
!> This routine is an auxiliary routine used in the Householder
!> reconstruction routine CUNHR_COL. In CUNHR_COL, this routine is
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
!> This is the blocked right-looking version of the algorithm,
!> calling Level 3 BLAS to update the submatrix. To factorize a block,
!> this routine calls the recursive routine CLAUNHR_COL_GETRFNP2.
!>
!> [1] "Reconstructing Householder vectors from tall-skinny QR",
!>     G. Ballard, J. Demmel, L. Grigori, M. Jacquelin, H.D. Nguyen,
!>     E. Solomonik, J. Parallel Distrib. Comput.,
!>     vol. 85, pp. 3-31, 2015.
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
!>          D is COMPLEX array, dimension min(M,N)
!>          The diagonal elements of the diagonal M-by-N sign matrix S,
!>          D(i) = S(i,i), where 1 <= i <= min(M,N). The elements can be
!>          only ( +1.0, 0.0 ) or (-1.0, 0.0 ).
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
!> \ingroup complexGEcomputational
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
      SUBROUTINE CLAUNHR_COL_GETRFNP(M,N,A,Lda,D,Info)
      IMPLICIT NONE
!*--CLAUNHR_COL_GETRFNP150
!
!  -- LAPACK computational routine (version 3.9.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2019
!
!     .. Scalar Arguments ..
      INTEGER Info , Lda , M , N
!     ..
!     .. Array Arguments ..
      COMPLEX A(Lda,*) , D(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      COMPLEX CONE
      PARAMETER (CONE=(1.0E+0,0.0E+0))
!     ..
!     .. Local Scalars ..
      INTEGER iinfo , j , jb , nb
!     ..
!     .. External Subroutines ..
      EXTERNAL CGEMM , CLAUNHR_COL_GETRFNP2 , CTRSM , XERBLA
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
         CALL XERBLA('CLAUNHR_COL_GETRFNP',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( MIN(M,N)==0 ) RETURN
!
!     Determine the block size for this environment.
!
 
      nb = ILAENV(1,'CLAUNHR_COL_GETRFNP',' ',M,N,-1,-1)
 
      IF ( nb<=1 .OR. nb>=MIN(M,N) ) THEN
!
!        Use unblocked code.
!
         CALL CLAUNHR_COL_GETRFNP2(M,N,A,Lda,D,Info)
      ELSE
!
!        Use blocked code.
!
         DO j = 1 , MIN(M,N) , nb
            jb = MIN(MIN(M,N)-j+1,nb)
!
!           Factor diagonal and subdiagonal blocks.
!
            CALL CLAUNHR_COL_GETRFNP2(M-j+1,jb,A(j,j),Lda,D(j),iinfo)
!
            IF ( j+jb<=N ) THEN
!
!              Compute block row of U.
!
               CALL CTRSM('Left','Lower','No transpose','Unit',jb,      &
     &                    N-j-jb+1,CONE,A(j,j),Lda,A(j,j+jb),Lda)
!
!                 Update trailing submatrix.
!
               IF ( j+jb<=M ) CALL CGEMM('No transpose','No transpose', &
     &              M-j-jb+1,N-j-jb+1,jb,-CONE,A(j+jb,j),Lda,A(j,j+jb), &
     &              Lda,CONE,A(j+jb,j+jb),Lda)
            ENDIF
         ENDDO
      ENDIF
!
!     End of CLAUNHR_COL_GETRFNP
!
      END SUBROUTINE CLAUNHR_COL_GETRFNP
