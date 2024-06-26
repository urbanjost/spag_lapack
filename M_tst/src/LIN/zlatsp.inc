!*==zlatsp.f90  processed by SPAG 7.51RB at 17:35 on  4 Mar 2022
!> \brief \b zlatsp
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE zlatsp( UPLO, N, X, ISEED )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            N
!       ..
!       .. Array Arguments ..
!       INTEGER            ISEED( * )
!       COMPLEX*16         X( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> zlatsp generates a special test matrix for the complex symmetric
!> (indefinite) factorization for packed matrices.  The pivot blocks of
!> the generated matrix will be in the following order:
!>    2x2 pivot block, non diagonalizable
!>    1x1 pivot block
!>    2x2 pivot block, diagonalizable
!>    (cycle repeats)
!> A row interchange is required for each non-diagonalizable 2x2 block.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER
!>          Specifies whether the generated matrix is to be upper or
!>          lower triangular.
!>          = 'U':  Upper triangular
!>          = 'L':  Lower triangular
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The dimension of the matrix to be generated.
!> \endverbatim
!>
!> \param[out] X
!> \verbatim
!>          X is COMPLEX*16 array, dimension (N*(N+1)/2)
!>          The generated matrix in packed storage format.  The matrix
!>          consists of 3x3 and 2x2 diagonal blocks which result in the
!>          pivot sequence given above.  The matrix outside these
!>          diagonal blocks is zero.
!> \endverbatim
!>
!> \param[in,out] ISEED
!> \verbatim
!>          ISEED is INTEGER array, dimension (4)
!>          On entry, the seed for the random number generator.  The last
!>          of the four integers must be odd.  (modified on exit)
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
!> \ingroup complex16_lin
!
!  =====================================================================
      SUBROUTINE zlatsp(Uplo,N,X,Iseed)
      IMPLICIT NONE
!*--ZLATSP88
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Uplo
      INTEGER N
!     ..
!     .. Array Arguments ..
      INTEGER Iseed(*)
      COMPLEX*16 X(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      COMPLEX*16 EYE
      PARAMETER (EYE=(0.0D0,1.0D0))
!     ..
!     .. Local Scalars ..
      INTEGER j , jj , n5
      DOUBLE PRECISION alpha , alpha3 , beta
      COMPLEX*16 a , b , c , r
!     ..
!     .. External Functions ..
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , SQRT
!     ..
!     .. Executable Statements ..
!
!     Initialize constants
!
      alpha = (1.D0+SQRT(17.D0))/8.D0
      beta = alpha - 1.D0/1000.D0
      alpha3 = alpha*alpha*alpha
!
!     Fill the matrix with zeros.
!
      DO j = 1 , N*(N+1)/2
         X(j) = 0.0D0
      ENDDO
!
!     UPLO = 'U':  Upper triangular storage
!
      IF ( Uplo=='U' ) THEN
         n5 = N/5
         n5 = N - 5*n5 + 1
!
         jj = N*(N+1)/2
         DO j = N , n5 , -5
            a = alpha3*ZLARND(5,Iseed)
            b = ZLARND(5,Iseed)/alpha
            c = a - 2.D0*b*EYE
            r = c/beta
            X(jj) = a
            X(jj-2) = b
            jj = jj - j
            X(jj) = ZLARND(2,Iseed)
            X(jj-1) = r
            jj = jj - (j-1)
            X(jj) = c
            jj = jj - (j-2)
            X(jj) = ZLARND(2,Iseed)
            jj = jj - (j-3)
            X(jj) = ZLARND(2,Iseed)
            IF ( ABS(X(jj+(j-3)))>ABS(X(jj)) ) THEN
               X(jj+(j-4)) = 2.0D0*X(jj+(j-3))
            ELSE
               X(jj+(j-4)) = 2.0D0*X(jj)
            ENDIF
            jj = jj - (j-4)
         ENDDO
!
!        Clean-up for N not a multiple of 5.
!
         j = n5 - 1
         IF ( j>2 ) THEN
            a = alpha3*ZLARND(5,Iseed)
            b = ZLARND(5,Iseed)/alpha
            c = a - 2.D0*b*EYE
            r = c/beta
            X(jj) = a
            X(jj-2) = b
            jj = jj - j
            X(jj) = ZLARND(2,Iseed)
            X(jj-1) = r
            jj = jj - (j-1)
            X(jj) = c
            jj = jj - (j-2)
            j = j - 3
         ENDIF
         IF ( j>1 ) THEN
            X(jj) = ZLARND(2,Iseed)
            X(jj-j) = ZLARND(2,Iseed)
            IF ( ABS(X(jj))>ABS(X(jj-j)) ) THEN
               X(jj-1) = 2.0D0*X(jj)
            ELSE
               X(jj-1) = 2.0D0*X(jj-j)
            ENDIF
            jj = jj - j - (j-1)
            j = j - 2
         ELSEIF ( j==1 ) THEN
            X(jj) = ZLARND(2,Iseed)
            j = j - 1
         ENDIF
!
!     UPLO = 'L':  Lower triangular storage
!
      ELSE
         n5 = N/5
         n5 = n5*5
!
         jj = 1
         DO j = 1 , n5 , 5
            a = alpha3*ZLARND(5,Iseed)
            b = ZLARND(5,Iseed)/alpha
            c = a - 2.D0*b*EYE
            r = c/beta
            X(jj) = a
            X(jj+2) = b
            jj = jj + (N-j+1)
            X(jj) = ZLARND(2,Iseed)
            X(jj+1) = r
            jj = jj + (N-j)
            X(jj) = c
            jj = jj + (N-j-1)
            X(jj) = ZLARND(2,Iseed)
            jj = jj + (N-j-2)
            X(jj) = ZLARND(2,Iseed)
            IF ( ABS(X(jj-(N-j-2)))>ABS(X(jj)) ) THEN
               X(jj-(N-j-2)+1) = 2.0D0*X(jj-(N-j-2))
            ELSE
               X(jj-(N-j-2)+1) = 2.0D0*X(jj)
            ENDIF
            jj = jj + (N-j-3)
         ENDDO
!
!        Clean-up for N not a multiple of 5.
!
         j = n5 + 1
         IF ( j<N-1 ) THEN
            a = alpha3*ZLARND(5,Iseed)
            b = ZLARND(5,Iseed)/alpha
            c = a - 2.D0*b*EYE
            r = c/beta
            X(jj) = a
            X(jj+2) = b
            jj = jj + (N-j+1)
            X(jj) = ZLARND(2,Iseed)
            X(jj+1) = r
            jj = jj + (N-j)
            X(jj) = c
            jj = jj + (N-j-1)
            j = j + 3
         ENDIF
         IF ( j<N ) THEN
            X(jj) = ZLARND(2,Iseed)
            X(jj+(N-j+1)) = ZLARND(2,Iseed)
            IF ( ABS(X(jj))>ABS(X(jj+(N-j+1))) ) THEN
               X(jj+1) = 2.0D0*X(jj)
            ELSE
               X(jj+1) = 2.0D0*X(jj+(N-j+1))
            ENDIF
            jj = jj + (N-j+1) + (N-j)
            j = j + 2
         ELSEIF ( j==N ) THEN
            X(jj) = ZLARND(2,Iseed)
            jj = jj + (N-j+1)
            j = j + 1
         ENDIF
      ENDIF
!
!
!     End of zlatsp
!
      END SUBROUTINE zlatsp
