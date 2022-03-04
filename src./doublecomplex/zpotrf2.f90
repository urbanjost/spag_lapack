!*==zpotrf2.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b ZPOTRF2
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       RECURSIVE SUBROUTINE ZPOTRF2( UPLO, N, A, LDA, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            INFO, LDA, N
!       ..
!       .. Array Arguments ..
!       COMPLEX*16         A( LDA, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZPOTRF2 computes the Cholesky factorization of a Hermitian
!> positive definite matrix A using the recursive algorithm.
!>
!> The factorization has the form
!>    A = U**H * U,  if UPLO = 'U', or
!>    A = L  * L**H,  if UPLO = 'L',
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
!>          A is COMPLEX*16 array, dimension (LDA,N)
!>          On entry, the Hermitian matrix A.  If UPLO = 'U', the leading
!>          N-by-N upper triangular part of A contains the upper
!>          triangular part of the matrix A, and the strictly lower
!>          triangular part of A is not referenced.  If UPLO = 'L', the
!>          leading N-by-N lower triangular part of A contains the lower
!>          triangular part of the matrix A, and the strictly upper
!>          triangular part of A is not referenced.
!>
!>          On exit, if INFO = 0, the factor U or L from the Cholesky
!>          factorization A = U**H*U or A = L*L**H.
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
!> \date December 2016
!
!> \ingroup complex16POcomputational
!
!  =====================================================================
      RECURSIVE SUBROUTINE ZPOTRF2(Uplo,N,A,Lda,Info)
      USE F77KINDS                        
      USE S_DISNAN
      USE S_LSAME
      USE S_XERBLA
      USE S_ZHERK
      USE S_ZTRSM
      IMPLICIT NONE
!*--ZPOTRF2116
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      COMPLEX(CX16KIND) , PARAMETER  ::  CONE = (1.0D+0,0.0D+0)
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      REAL(R8KIND) :: ajj
      INTEGER :: iinfo , n1 , n2
      LOGICAL :: upper
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
      upper = LSAME(Uplo,'U')
      IF ( .NOT.upper .AND. .NOT.LSAME(Uplo,'L') ) THEN
         Info = -1
      ELSEIF ( N<0 ) THEN
         Info = -2
      ELSEIF ( Lda<MAX(1,N) ) THEN
         Info = -4
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('ZPOTRF2',-Info)
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
         ajj = DBLE(A(1,1))
         IF ( ajj<=ZERO .OR. DISNAN(ajj) ) THEN
            Info = 1
            RETURN
         ENDIF
!
!        Factor
!
         A(1,1) = SQRT(ajj)
!
!     Use recursive code
!
      ELSE
         n1 = N/2
         n2 = N - n1
!
!        Factor A11
!
         CALL ZPOTRF2(Uplo,n1,A(1,1),Lda,iinfo)
         IF ( iinfo/=0 ) THEN
            Info = iinfo
            RETURN
         ENDIF
!
!        Compute the Cholesky factorization A = U**H*U
!
         IF ( upper ) THEN
!
!           Update and scale A12
!
            CALL ZTRSM('L','U','C','N',n1,n2,CONE,A(1,1),Lda,A(1,n1+1), &
     &                 Lda)
!
!           Update and factor A22
!
            CALL ZHERK(Uplo,'C',n2,n1,-ONE,A(1,n1+1),Lda,ONE,           &
     &                 A(n1+1,n1+1),Lda)
            CALL ZPOTRF2(Uplo,n2,A(n1+1,n1+1),Lda,iinfo)
            IF ( iinfo/=0 ) THEN
               Info = iinfo + n1
               RETURN
            ENDIF
!
!        Compute the Cholesky factorization A = L*L**H
!
         ELSE
!
!           Update and scale A21
!
            CALL ZTRSM('R','L','C','N',n2,n1,CONE,A(1,1),Lda,A(n1+1,1), &
     &                 Lda)
!
!           Update and factor A22
!
            CALL ZHERK(Uplo,'N',n2,n1,-ONE,A(n1+1,1),Lda,ONE,           &
     &                 A(n1+1,n1+1),Lda)
            CALL ZPOTRF2(Uplo,n2,A(n1+1,n1+1),Lda,iinfo)
            IF ( iinfo/=0 ) THEN
               Info = iinfo + n1
               RETURN
            ENDIF
         ENDIF
      ENDIF
!
!     End of ZPOTRF2
!
      END SUBROUTINE ZPOTRF2
