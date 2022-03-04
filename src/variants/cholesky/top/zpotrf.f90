!*==zpotrf.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b ZPOTRF VARIANT: top-looking block version of the algorithm, calling Level 3 BLAS.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZPOTRF ( UPLO, N, A, LDA, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            INFO, LDA, N
!       ..
!       .. Array Arguments ..
!       COMPLEX*16         A( LDA, * )
!       ..
!
!  Purpose
!  =======
!
!>\details \b Purpose:
!>\verbatim
!>
!> ZPOTRF computes the Cholesky factorization of a real symmetric
!> positive definite matrix A.
!>
!> The factorization has the form
!>    A = U**H * U,  if UPLO = 'U', or
!>    A = L  * L**H,  if UPLO = 'L',
!> where U is an upper triangular matrix and L is lower triangular.
!>
!> This is the top-looking block version of the algorithm, calling Level 3 BLAS.
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
!>          A is COMPLEX*16 array, dimension (LDA,N)
!>          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
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
      SUBROUTINE ZPOTRF(Uplo,N,A,Lda,Info)
      IMPLICIT NONE
!*--ZPOTRF104
!
!  -- LAPACK computational routine (version 3.1) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Uplo
      INTEGER Info , Lda , N
!     ..
!     .. Array Arguments ..
      COMPLEX*16 A(Lda,*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ONE
      COMPLEX*16 CONE
      PARAMETER (ONE=1.0D+0,CONE=(1.0D+0,0.0D+0))
!     ..
!     .. Local Scalars ..
      LOGICAL upper
      INTEGER j , jb , nb
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      INTEGER ILAENV
      EXTERNAL LSAME , ILAENV
!     ..
!     .. External Subroutines ..
      EXTERNAL ZGEMM , ZPOTF2 , ZHERK , ZTRSM , XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX , MIN
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
         CALL XERBLA('ZPOTRF',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( N==0 ) RETURN
!
!     Determine the block size for this environment.
!
      nb = ILAENV(1,'ZPOTRF',Uplo,N,-1,-1,-1)
      IF ( nb<=1 .OR. nb>=N ) THEN
!
!        Use unblocked code.
!
         CALL ZPOTF2(Uplo,N,A,Lda,Info)
!
!        Use blocked code.
!
      ELSEIF ( upper ) THEN
!
!           Compute the Cholesky factorization A = U'*U.
!
         DO j = 1 , N , nb
 
            jb = MIN(nb,N-j+1)
!
!              Compute the current block.
!
            CALL ZTRSM('Left','Upper','Conjugate Transpose','Non-unit', &
     &                 j-1,jb,CONE,A(1,1),Lda,A(1,j),Lda)
 
            CALL ZHERK('Upper','Conjugate Transpose',jb,j-1,-ONE,A(1,j),&
     &                 Lda,ONE,A(j,j),Lda)
!
!              Update and factorize the current diagonal block and test
!              for non-positive-definiteness.
!
            CALL ZPOTF2('Upper',jb,A(j,j),Lda,Info)
            IF ( Info/=0 ) GOTO 100
 
         ENDDO
!
      ELSE
!
!           Compute the Cholesky factorization A = L*L'.
!
         DO j = 1 , N , nb
 
            jb = MIN(nb,N-j+1)
!
!              Compute the current block.
!
            CALL ZTRSM('Right','Lower','Conjugate Transpose','Non-unit',&
     &                 jb,j-1,CONE,A(1,1),Lda,A(j,1),Lda)
 
            CALL ZHERK('Lower','No Transpose',jb,j-1,-ONE,A(j,1),Lda,   &
     &                 ONE,A(j,j),Lda)
!
!              Update and factorize the current diagonal block and test
!              for non-positive-definiteness.
!
            CALL ZPOTF2('Lower',jb,A(j,j),Lda,Info)
            IF ( Info/=0 ) GOTO 100
 
         ENDDO
      ENDIF
      GOTO 99999
!
 100  Info = Info + j - 1
!
!
!     End of ZPOTRF
!
99999 END SUBROUTINE ZPOTRF
