!*==clatsy.f90  processed by SPAG 7.51RB at 17:34 on  4 Mar 2022
!> \brief \b clatsy
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE clatsy( UPLO, N, X, LDX, ISEED )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            LDX, N
!       ..
!       .. Array Arguments ..
!       INTEGER            ISEED( * )
!       COMPLEX            X( LDX, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> clatsy generates a special test matrix for the complex symmetric
!> (indefinite) factorization.  The pivot blocks of the generated matrix
!> will be in the following order:
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
!>          X is COMPLEX array, dimension (LDX,N)
!>          The generated matrix, consisting of 3x3 and 2x2 diagonal
!>          blocks which result in the pivot sequence given above.
!>          The matrix outside of these diagonal blocks is zero.
!> \endverbatim
!>
!> \param[in] LDX
!> \verbatim
!>          LDX is INTEGER
!>          The leading dimension of the array X.
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
!> \ingroup complex_lin
!
!  =====================================================================
      SUBROUTINE clatsy(Uplo,N,X,Ldx,Iseed)
      IMPLICIT NONE
!*--CLATSY93
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Uplo
      INTEGER Ldx , N
!     ..
!     .. Array Arguments ..
      INTEGER Iseed(*)
      COMPLEX X(Ldx,*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      COMPLEX EYE
      PARAMETER (EYE=(0.0,1.0))
!     ..
!     .. Local Scalars ..
      INTEGER i , j , n5
      REAL alpha , alpha3 , beta
      COMPLEX a , b , c , r
!     ..
!     .. External Functions ..
      COMPLEX CLARND
EXTERNAL    CLARND
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , SQRT
!     ..
!     .. Executable Statements ..
!
!     Initialize constants
!
      alpha = (1.+SQRT(17.))/8.
      beta = alpha - 1./1000.
      alpha3 = alpha*alpha*alpha
!
!     UPLO = 'U':  Upper triangular storage
!
      IF ( Uplo=='U' ) THEN
!
!        Fill the upper triangle of the matrix with zeros.
!
         DO j = 1 , N
            DO i = 1 , j
               X(i,j) = 0.0
            ENDDO
         ENDDO
         n5 = N/5
         n5 = N - 5*n5 + 1
!
         DO i = N , n5 , -5
            a = alpha3*CLARND(5,Iseed)
            b = CLARND(5,Iseed)/alpha
            c = a - 2.*b*EYE
            r = c/beta
            X(i,i) = a
            X(i-2,i) = b
            X(i-2,i-1) = r
            X(i-2,i-2) = c
            X(i-1,i-1) = CLARND(2,Iseed)
            X(i-3,i-3) = CLARND(2,Iseed)
            X(i-4,i-4) = CLARND(2,Iseed)
            IF ( ABS(X(i-3,i-3))>ABS(X(i-4,i-4)) ) THEN
               X(i-4,i-3) = 2.0*X(i-3,i-3)
            ELSE
               X(i-4,i-3) = 2.0*X(i-4,i-4)
            ENDIF
         ENDDO
!
!        Clean-up for N not a multiple of 5.
!
         i = n5 - 1
         IF ( i>2 ) THEN
            a = alpha3*CLARND(5,Iseed)
            b = CLARND(5,Iseed)/alpha
            c = a - 2.*b*EYE
            r = c/beta
            X(i,i) = a
            X(i-2,i) = b
            X(i-2,i-1) = r
            X(i-2,i-2) = c
            X(i-1,i-1) = CLARND(2,Iseed)
            i = i - 3
         ENDIF
         IF ( i>1 ) THEN
            X(i,i) = CLARND(2,Iseed)
            X(i-1,i-1) = CLARND(2,Iseed)
            IF ( ABS(X(i,i))>ABS(X(i-1,i-1)) ) THEN
               X(i-1,i) = 2.0*X(i,i)
            ELSE
               X(i-1,i) = 2.0*X(i-1,i-1)
            ENDIF
            i = i - 2
         ELSEIF ( i==1 ) THEN
            X(i,i) = CLARND(2,Iseed)
            i = i - 1
         ENDIF
!
!     UPLO = 'L':  Lower triangular storage
!
      ELSE
!
!        Fill the lower triangle of the matrix with zeros.
!
         DO j = 1 , N
            DO i = j , N
               X(i,j) = 0.0
            ENDDO
         ENDDO
         n5 = N/5
         n5 = n5*5
!
         DO i = 1 , n5 , 5
            a = alpha3*CLARND(5,Iseed)
            b = CLARND(5,Iseed)/alpha
            c = a - 2.*b*EYE
            r = c/beta
            X(i,i) = a
            X(i+2,i) = b
            X(i+2,i+1) = r
            X(i+2,i+2) = c
            X(i+1,i+1) = CLARND(2,Iseed)
            X(i+3,i+3) = CLARND(2,Iseed)
            X(i+4,i+4) = CLARND(2,Iseed)
            IF ( ABS(X(i+3,i+3))>ABS(X(i+4,i+4)) ) THEN
               X(i+4,i+3) = 2.0*X(i+3,i+3)
            ELSE
               X(i+4,i+3) = 2.0*X(i+4,i+4)
            ENDIF
         ENDDO
!
!        Clean-up for N not a multiple of 5.
!
         i = n5 + 1
         IF ( i<N-1 ) THEN
            a = alpha3*CLARND(5,Iseed)
            b = CLARND(5,Iseed)/alpha
            c = a - 2.*b*EYE
            r = c/beta
            X(i,i) = a
            X(i+2,i) = b
            X(i+2,i+1) = r
            X(i+2,i+2) = c
            X(i+1,i+1) = CLARND(2,Iseed)
            i = i + 3
         ENDIF
         IF ( i<N ) THEN
            X(i,i) = CLARND(2,Iseed)
            X(i+1,i+1) = CLARND(2,Iseed)
            IF ( ABS(X(i,i))>ABS(X(i+1,i+1)) ) THEN
               X(i+1,i) = 2.0*X(i,i)
            ELSE
               X(i+1,i) = 2.0*X(i+1,i+1)
            ENDIF
            i = i + 2
         ELSEIF ( i==N ) THEN
            X(i,i) = CLARND(2,Iseed)
            i = i + 1
         ENDIF
      ENDIF
!
!
!     End of clatsy
!
      END SUBROUTINE clatsy
