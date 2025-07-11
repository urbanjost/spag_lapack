!*==spotrf2.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b SPOTRF2
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       RECURSIVE SUBROUTINE SPOTRF2( UPLO, N, A, LDA, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            INFO, LDA, N
!       ..
!       .. Array Arguments ..
!       REAL               A( LDA, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SPOTRF2 computes the Cholesky factorization of a real symmetric
!> positive definite matrix A using the recursive algorithm.
!>
!> The factorization has the form
!>    A = U**T * U,  if UPLO = 'U', or
!>    A = L  * L**T,  if UPLO = 'L',
!> where U is an upper triangular matrix and L is lower triangular.
!>
!> This is the recursive version of the algorithm. It divides
!> the matrix into four submatrices:
!>
!>        [  A11 | A12  ]  where A11 is n1 by n1 and A22 is n2 by n2
!>    A = [ -----|----- ]  with n1 = n/2
!>        [  A21 | A22  ]       n2 = n-n1
!>
!> The subroutine calls itself to factor A11. Update and scale A21
!> or A12, update A22 then call itself to factor A22.
!>
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          = 'U':  Upper triangle of A is stored;
!>          = 'L':  Lower triangle of A is stored.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is REAL array, dimension (LDA,N)
!>          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
!>          N-by-N upper triangular part of A contains the upper
!>          triangular part of the matrix A, and the strictly lower
!>          triangular part of A is not referenced.  If UPLO = 'L', the
!>          leading N-by-N lower triangular part of A contains the lower
!>          triangular part of the matrix A, and the strictly upper
!>          triangular part of A is not referenced.
!>
!>          On exit, if INFO = 0, the factor U or L from the Cholesky
!>          factorization A = U**T*U or A = L*L**T.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,N).
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value
!>          > 0:  if INFO = i, the leading minor of order i is not
!>                positive definite, and the factorization could not be
!>                completed.
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
!> \date November 2017
!
!> \ingroup realPOcomputational
!
!  =====================================================================
      RECURSIVE SUBROUTINE SPOTRF2(Uplo,N,A,Lda,Info)
      IMPLICIT NONE
!*--SPOTRF2110
!
!  -- LAPACK computational routine (version 3.8.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2017
!
!     .. Scalar Arguments ..
      CHARACTER Uplo
      INTEGER Info , Lda , N
!     ..
!     .. Array Arguments ..
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
      LOGICAL upper
      INTEGER n1 , n2 , iinfo
!     ..
!     .. External Functions ..
      LOGICAL LSAME , SISNAN
      EXTERNAL LSAME , SISNAN
!     ..
!     .. External Subroutines ..
      EXTERNAL SSYRK , STRSM , XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX , SQRT
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters
!
      Info = 0
      upper = LSAME(Uplo,'U')
      IF ( .NOT.upper .AND. .NOT.LSAME(Uplo,'L') ) THEN
         Info = -1
      ELSEIF ( N<0 ) THEN
         Info = -2
      ELSEIF ( Lda<MAX(1,N) ) THEN
         Info = -4
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('SPOTRF2',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( N==0 ) RETURN
!
!     N=1 case
!
      IF ( N==1 ) THEN
!
!        Test for non-positive-definiteness
!
         IF ( A(1,1)<=ZERO .OR. SISNAN(A(1,1)) ) THEN
            Info = 1
            RETURN
         ENDIF
!
!        Factor
!
         A(1,1) = SQRT(A(1,1))
!
!     Use recursive code
!
      ELSE
         n1 = N/2
         n2 = N - n1
!
!        Factor A11
!
         CALL SPOTRF2(Uplo,n1,A(1,1),Lda,iinfo)
         IF ( iinfo/=0 ) THEN
            Info = iinfo
            RETURN
         ENDIF
!
!        Compute the Cholesky factorization A = U**T*U
!
         IF ( upper ) THEN
!
!           Update and scale A12
!
            CALL STRSM('L','U','T','N',n1,n2,ONE,A(1,1),Lda,A(1,n1+1),  &
     &                 Lda)
!
!           Update and factor A22
!
            CALL SSYRK(Uplo,'T',n2,n1,-ONE,A(1,n1+1),Lda,ONE,           &
     &                 A(n1+1,n1+1),Lda)
            CALL SPOTRF2(Uplo,n2,A(n1+1,n1+1),Lda,iinfo)
            IF ( iinfo/=0 ) THEN
               Info = iinfo + n1
               RETURN
            ENDIF
!
!        Compute the Cholesky factorization A = L*L**T
!
         ELSE
!
!           Update and scale A21
!
            CALL STRSM('R','L','T','N',n2,n1,ONE,A(1,1),Lda,A(n1+1,1),  &
     &                 Lda)
!
!           Update and factor A22
!
            CALL SSYRK(Uplo,'N',n2,n1,-ONE,A(n1+1,1),Lda,ONE,           &
     &                 A(n1+1,n1+1),Lda)
            CALL SPOTRF2(Uplo,n2,A(n1+1,n1+1),Lda,iinfo)
            IF ( iinfo/=0 ) THEN
               Info = iinfo + n1
               RETURN
            ENDIF
         ENDIF
      ENDIF
!
!     End of SPOTRF2
!
      END SUBROUTINE SPOTRF2
