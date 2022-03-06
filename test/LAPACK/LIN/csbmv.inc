!*==csbmv.f90  processed by SPAG 7.51RB at 17:34 on  4 Mar 2022
!> \brief \b csbmv
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE csbmv( UPLO, N, K, ALPHA, A, LDA, X, INCX, BETA, Y,
!                         INCY )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            INCX, INCY, K, LDA, N
!       COMPLEX            ALPHA, BETA
!       ..
!       .. Array Arguments ..
!       COMPLEX            A( LDA, * ), X( * ), Y( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> csbmv  performs the matrix-vector  operation
!>
!>    y := alpha*A*x + beta*y,
!>
!> where alpha and beta are scalars, x and y are n element vectors and
!> A is an n by n symmetric band matrix, with k super-diagonals.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \verbatim
!>  UPLO   - CHARACTER*1
!>           On entry, UPLO specifies whether the upper or lower
!>           triangular part of the band matrix A is being supplied as
!>           follows:
!>
!>              UPLO = 'U' or 'u'   The upper triangular part of A is
!>                                  being supplied.
!>
!>              UPLO = 'L' or 'l'   The lower triangular part of A is
!>                                  being supplied.
!>
!>           Unchanged on exit.
!>
!>  N      - INTEGER
!>           On entry, N specifies the order of the matrix A.
!>           N must be at least zero.
!>           Unchanged on exit.
!>
!>  K      - INTEGER
!>           On entry, K specifies the number of super-diagonals of the
!>           matrix A. K must satisfy  0 .le. K.
!>           Unchanged on exit.
!>
!>  ALPHA  - COMPLEX
!>           On entry, ALPHA specifies the scalar alpha.
!>           Unchanged on exit.
!>
!>  A      - COMPLEX array, dimension( LDA, N )
!>           Before entry with UPLO = 'U' or 'u', the leading ( k + 1 )
!>           by n part of the array A must contain the upper triangular
!>           band part of the symmetric matrix, supplied column by
!>           column, with the leading diagonal of the matrix in row
!>           ( k + 1 ) of the array, the first super-diagonal starting at
!>           position 2 in row k, and so on. The top left k by k triangle
!>           of the array A is not referenced.
!>           The following program segment will transfer the upper
!>           triangular part of a symmetric band matrix from conventional
!>           full matrix storage to band storage:
!>
!>                 DO 20, J = 1, N
!>                    M = K + 1 - J
!>                    DO 10, I = MAX( 1, J - K ), J
!>                       A( M + I, J ) = matrix( I, J )
!>              10    CONTINUE
!>              20 CONTINUE
!>
!>           Before entry with UPLO = 'L' or 'l', the leading ( k + 1 )
!>           by n part of the array A must contain the lower triangular
!>           band part of the symmetric matrix, supplied column by
!>           column, with the leading diagonal of the matrix in row 1 of
!>           the array, the first sub-diagonal starting at position 1 in
!>           row 2, and so on. The bottom right k by k triangle of the
!>           array A is not referenced.
!>           The following program segment will transfer the lower
!>           triangular part of a symmetric band matrix from conventional
!>           full matrix storage to band storage:
!>
!>                 DO 20, J = 1, N
!>                    M = 1 - J
!>                    DO 10, I = J, MIN( N, J + K )
!>                       A( M + I, J ) = matrix( I, J )
!>              10    CONTINUE
!>              20 CONTINUE
!>
!>           Unchanged on exit.
!>
!>  LDA    - INTEGER
!>           On entry, LDA specifies the first dimension of A as declared
!>           in the calling (sub) program. LDA must be at least
!>           ( k + 1 ).
!>           Unchanged on exit.
!>
!>  X      - COMPLEX array, dimension at least
!>           ( 1 + ( N - 1 )*abs( INCX ) ).
!>           Before entry, the incremented array X must contain the
!>           vector x.
!>           Unchanged on exit.
!>
!>  INCX   - INTEGER
!>           On entry, INCX specifies the increment for the elements of
!>           X. INCX must not be zero.
!>           Unchanged on exit.
!>
!>  BETA   - COMPLEX
!>           On entry, BETA specifies the scalar beta.
!>           Unchanged on exit.
!>
!>  Y      - COMPLEX array, dimension at least
!>           ( 1 + ( N - 1 )*abs( INCY ) ).
!>           Before entry, the incremented array Y must contain the
!>           vector y. On exit, Y is overwritten by the updated vector y.
!>
!>  INCY   - INTEGER
!>           On entry, INCY specifies the increment for the elements of
!>           Y. INCY must not be zero.
!>           Unchanged on exit.
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
      SUBROUTINE csbmv(Uplo,N,K,Alpha,A,Lda,X,Incx,Beta,Y,Incy)
      IMPLICIT NONE
!*--CSBMV155
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Uplo
      INTEGER Incx , Incy , K , Lda , N
      COMPLEX Alpha , Beta
!     ..
!     .. Array Arguments ..
      COMPLEX A(Lda,*) , X(*) , Y(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      COMPLEX ONE
      PARAMETER (ONE=(1.0E+0,0.0E+0))
      COMPLEX ZERO
      PARAMETER (ZERO=(0.0E+0,0.0E+0))
!     ..
!     .. Local Scalars ..
      INTEGER i , info , ix , iy , j , jx , jy , kplus1 , kx , ky , l
      COMPLEX temp1 , temp2
!     ..
!     .. External Functions ..
      LOGICAL LSAME
EXTERNAL    LSAME
!     ..
!     .. External Subroutines ..
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX , MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      info = 0
      IF ( .NOT.LSAME(Uplo,'U') .AND. .NOT.LSAME(Uplo,'L') ) THEN
         info = 1
      ELSEIF ( N<0 ) THEN
         info = 2
      ELSEIF ( K<0 ) THEN
         info = 3
      ELSEIF ( Lda<(K+1) ) THEN
         info = 6
      ELSEIF ( Incx==0 ) THEN
         info = 8
      ELSEIF ( Incy==0 ) THEN
         info = 11
      ENDIF
      IF ( info/=0 ) THEN
         CALL xerbla('csbmv ',info)
         RETURN
      ENDIF
!
!     Quick return if possible.
!
      IF ( (N==0) .OR. ((Alpha==ZERO) .AND. (Beta==ONE)) ) RETURN
!
!     Set up the start points in  X  and  Y.
!
      IF ( Incx>0 ) THEN
         kx = 1
      ELSE
         kx = 1 - (N-1)*Incx
      ENDIF
      IF ( Incy>0 ) THEN
         ky = 1
      ELSE
         ky = 1 - (N-1)*Incy
      ENDIF
!
!     Start the operations. In this version the elements of the array A
!     are accessed sequentially with one pass through A.
!
!     First form  y := beta*y.
!
      IF ( Beta/=ONE ) THEN
         IF ( Incy/=1 ) THEN
            iy = ky
            IF ( Beta==ZERO ) THEN
               DO i = 1 , N
                  Y(iy) = ZERO
                  iy = iy + Incy
               ENDDO
            ELSE
               DO i = 1 , N
                  Y(iy) = Beta*Y(iy)
                  iy = iy + Incy
               ENDDO
            ENDIF
         ELSEIF ( Beta==ZERO ) THEN
            DO i = 1 , N
               Y(i) = ZERO
            ENDDO
         ELSE
            DO i = 1 , N
               Y(i) = Beta*Y(i)
            ENDDO
         ENDIF
      ENDIF
      IF ( Alpha==ZERO ) RETURN
      IF ( LSAME(Uplo,'U') ) THEN
!
!        Form  y  when upper triangle of A is stored.
!
         kplus1 = K + 1
         IF ( (Incx==1) .AND. (Incy==1) ) THEN
            DO j = 1 , N
               temp1 = Alpha*X(j)
               temp2 = ZERO
               l = kplus1 - j
               DO i = MAX(1,j-K) , j - 1
                  Y(i) = Y(i) + temp1*A(l+i,j)
                  temp2 = temp2 + A(l+i,j)*X(i)
               ENDDO
               Y(j) = Y(j) + temp1*A(kplus1,j) + Alpha*temp2
            ENDDO
         ELSE
            jx = kx
            jy = ky
            DO j = 1 , N
               temp1 = Alpha*X(jx)
               temp2 = ZERO
               ix = kx
               iy = ky
               l = kplus1 - j
               DO i = MAX(1,j-K) , j - 1
                  Y(iy) = Y(iy) + temp1*A(l+i,j)
                  temp2 = temp2 + A(l+i,j)*X(ix)
                  ix = ix + Incx
                  iy = iy + Incy
               ENDDO
               Y(jy) = Y(jy) + temp1*A(kplus1,j) + Alpha*temp2
               jx = jx + Incx
               jy = jy + Incy
               IF ( j>K ) THEN
                  kx = kx + Incx
                  ky = ky + Incy
               ENDIF
            ENDDO
         ENDIF
!
!        Form  y  when lower triangle of A is stored.
!
      ELSEIF ( (Incx==1) .AND. (Incy==1) ) THEN
         DO j = 1 , N
            temp1 = Alpha*X(j)
            temp2 = ZERO
            Y(j) = Y(j) + temp1*A(1,j)
            l = 1 - j
            DO i = j + 1 , MIN(N,j+K)
               Y(i) = Y(i) + temp1*A(l+i,j)
               temp2 = temp2 + A(l+i,j)*X(i)
            ENDDO
            Y(j) = Y(j) + Alpha*temp2
         ENDDO
      ELSE
         jx = kx
         jy = ky
         DO j = 1 , N
            temp1 = Alpha*X(jx)
            temp2 = ZERO
            Y(jy) = Y(jy) + temp1*A(1,j)
            l = 1 - j
            ix = jx
            iy = jy
            DO i = j + 1 , MIN(N,j+K)
               ix = ix + Incx
               iy = iy + Incy
               Y(iy) = Y(iy) + temp1*A(l+i,j)
               temp2 = temp2 + A(l+i,j)*X(ix)
            ENDDO
            Y(jy) = Y(jy) + Alpha*temp2
            jx = jx + Incx
            jy = jy + Incy
         ENDDO
      ENDIF
!
!
!     End of csbmv
!
      END SUBROUTINE csbmv
