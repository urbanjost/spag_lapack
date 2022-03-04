!*==cpotrf.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b CPOTRF VARIANT: right looking block version of the algorithm, calling Level 3 BLAS.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE CPOTRF ( UPLO, N, A, LDA, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            INFO, LDA, N
!       ..
!       .. Array Arguments ..
!       COMPLEX            A( LDA, * )
!       ..
!
!  Purpose
!  =======
!
!>\details \b Purpose:
!>\verbatim
!>
!> CPOTRF computes the Cholesky factorization of a real Hermitian
!> positive definite matrix A.
!>
!> The factorization has the form
!>    A = U**H * U,  if UPLO = 'U', or
!>    A = L  * L**H,  if UPLO = 'L',
!> where U is an upper triangular matrix and L is lower triangular.
!>
!> This is the right looking block version of the algorithm, calling Level 3 BLAS.
!>
!>\endverbatim
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
!>          A is COMPLEX array, dimension (LDA,N)
!>          On entry, the Hermitian matrix A.  If UPLO = 'U', the leading
!>          N-by-N upper triangular part of A contains the upper
!>          triangular part of the matrix A, and the strictly lower
!>          triangular part of A is not referenced.  If UPLO = 'L', the
!>          leading N-by-N lower triangular part of A contains the lower
!>          triangular part of the matrix A, and the strictly upper
!>          triangular part of A is not referenced.
!> \endverbatim
!> \verbatim
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
!> \ingroup variantsPOcomputational
!
!  =====================================================================
      SUBROUTINE CPOTRF(Uplo,N,A,Lda,Info)
      USE S_CGEMM
      USE S_CHERK
      USE S_CPOTF2
      USE S_CTRSM
      USE S_ILAENV
      USE S_LSAME
      USE S_XERBLA
      IMPLICIT NONE
!*--CPOTRF111
!
! PARAMETER definitions rewritten by SPAG
!
      REAL , PARAMETER  ::  ONE = 1.0E+0
      COMPLEX , PARAMETER  ::  CONE = (1.0E+0,0.0E+0)
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      INTEGER :: j , jb , nb
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
!     Test the input parameters.
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
         CALL XERBLA('CPOTRF',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( N==0 ) RETURN
!
!     Determine the block size for this environment.
!
      nb = ILAENV(1,'CPOTRF',Uplo,N,-1,-1,-1)
      IF ( nb<=1 .OR. nb>=N ) THEN
!
!        Use unblocked code.
!
         CALL CPOTF2(Uplo,N,A,Lda,Info)
!
!        Use blocked code.
!
      ELSEIF ( upper ) THEN
!
!           Compute the Cholesky factorization A = U'*U.
!
         DO j = 1 , N , nb
!
!              Update and factorize the current diagonal block and test
!              for non-positive-definiteness.
!
            jb = MIN(nb,N-j+1)
 
            CALL CPOTF2('Upper',jb,A(j,j),Lda,Info)
 
            IF ( Info/=0 ) GOTO 100
 
            IF ( j+jb<=N ) THEN
!
!                 Updating the trailing submatrix.
!
               CALL CTRSM('Left','Upper','Conjugate Transpose',         &
     &                    'Non-unit',jb,N-j-jb+1,CONE,A(j,j),Lda,       &
     &                    A(j,j+jb),Lda)
               CALL CHERK('Upper','Conjugate transpose',N-j-jb+1,jb,    &
     &                    -ONE,A(j,j+jb),Lda,ONE,A(j+jb,j+jb),Lda)
            ENDIF
         ENDDO
!
      ELSE
!
!           Compute the Cholesky factorization A = L*L'.
!
         DO j = 1 , N , nb
!
!              Update and factorize the current diagonal block and test
!              for non-positive-definiteness.
!
            jb = MIN(nb,N-j+1)
 
            CALL CPOTF2('Lower',jb,A(j,j),Lda,Info)
 
            IF ( Info/=0 ) GOTO 100
 
            IF ( j+jb<=N ) THEN
!
!                Updating the trailing submatrix.
!
               CALL CTRSM('Right','Lower','Conjugate Transpose',        &
     &                    'Non-unit',N-j-jb+1,jb,CONE,A(j,j),Lda,       &
     &                    A(j+jb,j),Lda)
 
               CALL CHERK('Lower','No Transpose',N-j-jb+1,jb,-ONE,      &
     &                    A(j+jb,j),Lda,ONE,A(j+jb,j+jb),Lda)
            ENDIF
         ENDDO
      ENDIF
      GOTO 99999
!
 100  Info = Info + j - 1
!
!
!     End of CPOTRF
!
99999 END SUBROUTINE CPOTRF
