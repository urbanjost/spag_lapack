!*==dlaset.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b DLASET initializes the off-diagonal elements and the diagonal elements of a matrix to given values.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DLASET + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaset.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaset.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaset.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLASET( UPLO, M, N, ALPHA, BETA, A, LDA )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            LDA, M, N
!       DOUBLE PRECISION   ALPHA, BETA
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   A( LDA, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLASET initializes an m-by-n matrix A to BETA on the diagonal and
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
!>          = 'U':      Upper triangular part is set; the strictly lower
!>                      triangular part of A is not changed.
!>          = 'L':      Lower triangular part is set; the strictly upper
!>                      triangular part of A is not changed.
!>          Otherwise:  All of the matrix A is set.
!> \endverbatim
!>
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
!> \param[in] ALPHA
!> \verbatim
!>          ALPHA is DOUBLE PRECISION
!>          The constant to which the offdiagonal elements are to be set.
!> \endverbatim
!>
!> \param[in] BETA
!> \verbatim
!>          BETA is DOUBLE PRECISION
!>          The constant to which the diagonal elements are to be set.
!> \endverbatim
!>
!> \param[out] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (LDA,N)
!>          On exit, the leading m-by-n submatrix of A is set as follows:
!>
!>          if UPLO = 'U', A(i,j) = ALPHA, 1<=i<=j-1, 1<=j<=n,
!>          if UPLO = 'L', A(i,j) = ALPHA, j+1<=i<=m, 1<=j<=n,
!>          otherwise,     A(i,j) = ALPHA, 1<=i<=m, 1<=j<=n, i.ne.j,
!>
!>          and, for all UPLO, A(i,i) = BETA, 1<=i<=min(m,n).
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
!> \ingroup OTHERauxiliary
!
!  =====================================================================
      SUBROUTINE DLASET(Uplo,M,N,Alpha,Beta,A,Lda)
      USE F77KINDS                        
      USE S_LSAME
      IMPLICIT NONE
!*--DLASET116
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(IN) :: Alpha
      REAL(R8KIND) , INTENT(IN) :: Beta
      REAL(R8KIND) , INTENT(OUT) , DIMENSION(Lda,*) :: A
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
! =====================================================================
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
!        Set the strictly upper triangular or trapezoidal part of the
!        array to ALPHA.
!
         DO j = 2 , N
            DO i = 1 , MIN(j-1,M)
               A(i,j) = Alpha
            ENDDO
         ENDDO
!
      ELSEIF ( LSAME(Uplo,'L') ) THEN
!
!        Set the strictly lower triangular or trapezoidal part of the
!        array to ALPHA.
!
         DO j = 1 , MIN(M,N)
            DO i = j + 1 , M
               A(i,j) = Alpha
            ENDDO
         ENDDO
!
      ELSE
!
!        Set the leading m-by-n submatrix to ALPHA.
!
         DO j = 1 , N
            DO i = 1 , M
               A(i,j) = Alpha
            ENDDO
         ENDDO
      ENDIF
!
!     Set the first min(M,N) diagonal elements to BETA.
!
      DO i = 1 , MIN(M,N)
         A(i,i) = Beta
      ENDDO
!
!
!     End of DLASET
!
      END SUBROUTINE DLASET
