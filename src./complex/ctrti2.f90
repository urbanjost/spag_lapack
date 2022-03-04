!*==ctrti2.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b CTRTI2 computes the inverse of a triangular matrix (unblocked algorithm).
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CTRTI2 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ctrti2.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ctrti2.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ctrti2.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CTRTI2( UPLO, DIAG, N, A, LDA, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          DIAG, UPLO
!       INTEGER            INFO, LDA, N
!       ..
!       .. Array Arguments ..
!       COMPLEX            A( LDA, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CTRTI2 computes the inverse of a complex upper or lower triangular
!> matrix.
!>
!> This is the Level 2 BLAS version of the algorithm.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          Specifies whether the matrix A is upper or lower triangular.
!>          = 'U':  Upper triangular
!>          = 'L':  Lower triangular
!> \endverbatim
!>
!> \param[in] DIAG
!> \verbatim
!>          DIAG is CHARACTER*1
!>          Specifies whether or not the matrix A is unit triangular.
!>          = 'N':  Non-unit triangular
!>          = 'U':  Unit triangular
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
!>          On entry, the triangular matrix A.  If UPLO = 'U', the
!>          leading n by n upper triangular part of the array A contains
!>          the upper triangular matrix, and the strictly lower
!>          triangular part of A is not referenced.  If UPLO = 'L', the
!>          leading n by n lower triangular part of the array A contains
!>          the lower triangular matrix, and the strictly upper
!>          triangular part of A is not referenced.  If DIAG = 'U', the
!>          diagonal elements of A are also not referenced and are
!>          assumed to be 1.
!>
!>          On exit, the (triangular) inverse of the original matrix, in
!>          the same storage format.
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
!> \ingroup complexOTHERcomputational
!
!  =====================================================================
      SUBROUTINE CTRTI2(Uplo,Diag,N,A,Lda,Info)
      USE S_CSCAL
      USE S_CTRMV
      USE S_LSAME
      USE S_XERBLA
      IMPLICIT NONE
!*--CTRTI2118
!
! PARAMETER definitions rewritten by SPAG
!
      COMPLEX , PARAMETER  ::  ONE = (1.0E+0,0.0E+0)
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Uplo
      CHARACTER :: Diag
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      COMPLEX :: ajj
      INTEGER :: j
      LOGICAL :: nounit , upper
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
      nounit = LSAME(Diag,'N')
      IF ( .NOT.upper .AND. .NOT.LSAME(Uplo,'L') ) THEN
         Info = -1
      ELSEIF ( .NOT.nounit .AND. .NOT.LSAME(Diag,'U') ) THEN
         Info = -2
      ELSEIF ( N<0 ) THEN
         Info = -3
      ELSEIF ( Lda<MAX(1,N) ) THEN
         Info = -5
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('CTRTI2',-Info)
         RETURN
      ENDIF
!
      IF ( upper ) THEN
!
!        Compute inverse of upper triangular matrix.
!
         DO j = 1 , N
            IF ( nounit ) THEN
               A(j,j) = ONE/A(j,j)
               ajj = -A(j,j)
            ELSE
               ajj = -ONE
            ENDIF
!
!           Compute elements 1:j-1 of j-th column.
!
            CALL CTRMV('Upper','No transpose',Diag,j-1,A,Lda,A(1,j),1)
            CALL CSCAL(j-1,ajj,A(1,j),1)
         ENDDO
      ELSE
!
!        Compute inverse of lower triangular matrix.
!
         DO j = N , 1 , -1
            IF ( nounit ) THEN
               A(j,j) = ONE/A(j,j)
               ajj = -A(j,j)
            ELSE
               ajj = -ONE
            ENDIF
            IF ( j<N ) THEN
!
!              Compute elements j+1:n of j-th column.
!
               CALL CTRMV('Lower','No transpose',Diag,N-j,A(j+1,j+1),   &
     &                    Lda,A(j+1,j),1)
               CALL CSCAL(N-j,ajj,A(j+1,j),1)
            ENDIF
         ENDDO
      ENDIF
!
!
!     End of CTRTI2
!
      END SUBROUTINE CTRTI2
