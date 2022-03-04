!*==claset.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b CLASET initializes the off-diagonal elements and the diagonal elements of a matrix to given values.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CLASET + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/claset.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/claset.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/claset.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CLASET( UPLO, M, N, ALPHA, BETA, A, LDA )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            LDA, M, N
!       COMPLEX            ALPHA, BETA
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
!> CLASET initializes a 2-D array A to BETA on the diagonal and
!> ALPHA on the offdiagonals.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          Specifies the part of the matrix A to be set.
!>          = 'U':      Upper triangular part is set. The lower triangle
!>                      is unchanged.
!>          = 'L':      Lower triangular part is set. The upper triangle
!>                      is unchanged.
!>          Otherwise:  All of the matrix A is set.
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          On entry, M specifies the number of rows of A.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          On entry, N specifies the number of columns of A.
!> \endverbatim
!>
!> \param[in] ALPHA
!> \verbatim
!>          ALPHA is COMPLEX
!>          All the offdiagonal array elements are set to ALPHA.
!> \endverbatim
!>
!> \param[in] BETA
!> \verbatim
!>          BETA is COMPLEX
!>          All the diagonal array elements are set to BETA.
!> \endverbatim
!>
!> \param[out] A
!> \verbatim
!>          A is COMPLEX array, dimension (LDA,N)
!>          On entry, the m by n matrix A.
!>          On exit, A(i,j) = ALPHA, 1 <= i <= m, 1 <= j <= n, i.ne.j;
!>                   A(i,i) = BETA , 1 <= i <= min(m,n)
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,M).
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
!> \ingroup complexOTHERauxiliary
!
!  =====================================================================
      SUBROUTINE CLASET(Uplo,M,N,Alpha,Beta,A,Lda)
      USE S_LSAME
      IMPLICIT NONE
!*--CLASET111
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(IN) :: Alpha
      COMPLEX , INTENT(IN) :: Beta
      COMPLEX , INTENT(OUT) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
!
! Local variable declarations rewritten by SPAG
!
      INTEGER :: i , j
!
! End of declarations rewritten by SPAG
!
!     ..
!     .. Array Arguments ..
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
!     ..
!     .. External Functions ..
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. Executable Statements ..
!
      IF ( LSAME(Uplo,'U') ) THEN
!
!        Set the diagonal to BETA and the strictly upper triangular
!        part of the array to ALPHA.
!
         DO j = 2 , N
            DO i = 1 , MIN(j-1,M)
               A(i,j) = Alpha
            ENDDO
         ENDDO
         DO i = 1 , MIN(N,M)
            A(i,i) = Beta
         ENDDO
!
      ELSEIF ( LSAME(Uplo,'L') ) THEN
!
!        Set the diagonal to BETA and the strictly lower triangular
!        part of the array to ALPHA.
!
         DO j = 1 , MIN(M,N)
            DO i = j + 1 , M
               A(i,j) = Alpha
            ENDDO
         ENDDO
         DO i = 1 , MIN(N,M)
            A(i,i) = Beta
         ENDDO
!
      ELSE
!
!        Set the array to BETA on the diagonal and ALPHA on the
!        offdiagonal.
!
         DO j = 1 , N
            DO i = 1 , M
               A(i,j) = Alpha
            ENDDO
         ENDDO
         DO i = 1 , MIN(M,N)
            A(i,i) = Beta
         ENDDO
      ENDIF
!
!
!     End of CLASET
!
      END SUBROUTINE CLASET
