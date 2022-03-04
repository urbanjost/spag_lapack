!*==zpotf2.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b ZPOTF2 computes the Cholesky factorization of a symmetric/Hermitian positive definite matrix (unblocked algorithm).
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZPOTF2 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zpotf2.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zpotf2.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zpotf2.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZPOTF2( UPLO, N, A, LDA, INFO )
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
!> ZPOTF2 computes the Cholesky factorization of a complex Hermitian
!> positive definite matrix A.
!>
!> The factorization has the form
!>    A = U**H * U ,  if UPLO = 'U', or
!>    A = L  * L**H,  if UPLO = 'L',
!> where U is an upper triangular matrix and L is lower triangular.
!>
!> This is the unblocked version of the algorithm, calling Level 2 BLAS.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          Specifies whether the upper or lower triangular part of the
!>          Hermitian matrix A is stored.
!>          = 'U':  Upper triangular
!>          = 'L':  Lower triangular
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
!>          n by n upper triangular part of A contains the upper
!>          triangular part of the matrix A, and the strictly lower
!>          triangular part of A is not referenced.  If UPLO = 'L', the
!>          leading n by n lower triangular part of A contains the lower
!>          triangular part of the matrix A, and the strictly upper
!>          triangular part of A is not referenced.
!>
!>          On exit, if INFO = 0, the factor U or L from the Cholesky
!>          factorization A = U**H *U  or A = L*L**H.
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
!>          = 0: successful exit
!>          < 0: if INFO = -k, the k-th argument had an illegal value
!>          > 0: if INFO = k, the leading minor of order k is not
!>               positive definite, and the factorization could not be
!>               completed.
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
      SUBROUTINE ZPOTF2(Uplo,N,A,Lda,Info)
      USE F77KINDS                        
      USE S_DISNAN
      USE S_LSAME
      USE S_XERBLA
      USE S_ZDOTC
      USE S_ZDSCAL
      USE S_ZGEMV
      USE S_ZLACGV
      IMPLICIT NONE
!*--ZPOTF2121
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
      INTEGER :: j
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
         CALL XERBLA('ZPOTF2',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( N==0 ) RETURN
!
      IF ( upper ) THEN
!
!        Compute the Cholesky factorization A = U**H *U.
!
         DO j = 1 , N
!
!           Compute U(J,J) and test for non-positive-definiteness.
!
            ajj = DBLE(A(j,j)) - ZDOTC(j-1,A(1,j),1,A(1,j),1)
            IF ( ajj<=ZERO .OR. DISNAN(ajj) ) THEN
               A(j,j) = ajj
               GOTO 100
            ENDIF
            ajj = SQRT(ajj)
            A(j,j) = ajj
!
!           Compute elements J+1:N of row J.
!
            IF ( j<N ) THEN
               CALL ZLACGV(j-1,A(1,j),1)
               CALL ZGEMV('Transpose',j-1,N-j,-CONE,A(1,j+1),Lda,A(1,j),&
     &                    1,CONE,A(j,j+1),Lda)
               CALL ZLACGV(j-1,A(1,j),1)
               CALL ZDSCAL(N-j,ONE/ajj,A(j,j+1),Lda)
            ENDIF
         ENDDO
      ELSE
!
!        Compute the Cholesky factorization A = L*L**H.
!
         DO j = 1 , N
!
!           Compute L(J,J) and test for non-positive-definiteness.
!
            ajj = DBLE(A(j,j)) - ZDOTC(j-1,A(j,1),Lda,A(j,1),Lda)
            IF ( ajj<=ZERO .OR. DISNAN(ajj) ) THEN
               A(j,j) = ajj
               GOTO 100
            ENDIF
            ajj = SQRT(ajj)
            A(j,j) = ajj
!
!           Compute elements J+1:N of column J.
!
            IF ( j<N ) THEN
               CALL ZLACGV(j-1,A(j,1),Lda)
               CALL ZGEMV('No transpose',N-j,j-1,-CONE,A(j+1,1),Lda,    &
     &                    A(j,1),Lda,CONE,A(j+1,j),1)
               CALL ZLACGV(j-1,A(j,1),Lda)
               CALL ZDSCAL(N-j,ONE/ajj,A(j+1,j),1)
            ENDIF
         ENDDO
      ENDIF
      GOTO 99999
!
 100  Info = j
!
!
!     End of ZPOTF2
!
99999 END SUBROUTINE ZPOTF2
