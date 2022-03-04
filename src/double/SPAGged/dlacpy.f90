!*==dlacpy.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b DLACPY copies all or part of one two-dimensional array to another.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DLACPY + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlacpy.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlacpy.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlacpy.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLACPY( UPLO, M, N, A, LDA, B, LDB )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            LDA, LDB, M, N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLACPY copies all or part of a two-dimensional matrix A to another
!> matrix B.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          Specifies the part of the matrix A to be copied to B.
!>          = 'U':      Upper triangular part
!>          = 'L':      Lower triangular part
!>          Otherwise:  All of the matrix A
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
!> \param[in] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (LDA,N)
!>          The m by n matrix A.  If UPLO = 'U', only the upper triangle
!>          or trapezoid is accessed; if UPLO = 'L', only the lower
!>          triangle or trapezoid is accessed.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,M).
!> \endverbatim
!>
!> \param[out] B
!> \verbatim
!>          B is DOUBLE PRECISION array, dimension (LDB,N)
!>          On exit, B = A in the locations specified by UPLO.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of the array B.  LDB >= max(1,M).
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
      SUBROUTINE DLACPY(Uplo,M,N,A,Lda,B,Ldb)
      USE F77KINDS                        
      USE S_LSAME
      IMPLICIT NONE
!*--DLACPY109
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL(R8KIND) , INTENT(OUT) , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
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
         DO j = 1 , N
            DO i = 1 , MIN(j,M)
               B(i,j) = A(i,j)
            ENDDO
         ENDDO
      ELSEIF ( LSAME(Uplo,'L') ) THEN
         DO j = 1 , N
            DO i = j , M
               B(i,j) = A(i,j)
            ENDDO
         ENDDO
      ELSE
         DO j = 1 , N
            DO i = 1 , M
               B(i,j) = A(i,j)
            ENDDO
         ENDDO
      ENDIF
!
!     End of DLACPY
!
      END SUBROUTINE DLACPY
