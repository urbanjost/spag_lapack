!*==chbmv.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b CHBMV
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE CHBMV(UPLO,N,K,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
!
!       .. Scalar Arguments ..
!       COMPLEX ALPHA,BETA
!       INTEGER INCX,INCY,K,LDA,N
!       CHARACTER UPLO
!       ..
!       .. Array Arguments ..
!       COMPLEX A(LDA,*),X(*),Y(*)
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CHBMV  performs the matrix-vector  operation
!>
!>    y := alpha*A*x + beta*y,
!>
!> where alpha and beta are scalars, x and y are n element vectors and
!> A is an n by n hermitian band matrix, with k super-diagonals.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>           On entry, UPLO specifies whether the upper or lower
!>           triangular part of the band matrix A is being supplied as
!>           follows:
!>
!>              UPLO = 'U' or 'u'   The upper triangular part of A is
!>                                  being supplied.
!>
!>              UPLO = 'L' or 'l'   The lower triangular part of A is
!>                                  being supplied.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>           On entry, N specifies the order of the matrix A.
!>           N must be at least zero.
!> \endverbatim
!>
!> \param[in] K
!> \verbatim
!>          K is INTEGER
!>           On entry, K specifies the number of super-diagonals of the
!>           matrix A. K must satisfy  0 .le. K.
!> \endverbatim
!>
!> \param[in] ALPHA
!> \verbatim
!>          ALPHA is COMPLEX
!>           On entry, ALPHA specifies the scalar alpha.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is COMPLEX array, dimension ( LDA, N )
!>           Before entry with UPLO = 'U' or 'u', the leading ( k + 1 )
!>           by n part of the array A must contain the upper triangular
!>           band part of the hermitian matrix, supplied column by
!>           column, with the leading diagonal of the matrix in row
!>           ( k + 1 ) of the array, the first super-diagonal starting at
!>           position 2 in row k, and so on. The top left k by k triangle
!>           of the array A is not referenced.
!>           The following program segment will transfer the upper
!>           triangular part of a hermitian band matrix from conventional
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
!>           band part of the hermitian matrix, supplied column by
!>           column, with the leading diagonal of the matrix in row 1 of
!>           the array, the first sub-diagonal starting at position 1 in
!>           row 2, and so on. The bottom right k by k triangle of the
!>           array A is not referenced.
!>           The following program segment will transfer the lower
!>           triangular part of a hermitian band matrix from conventional
!>           full matrix storage to band storage:
!>
!>                 DO 20, J = 1, N
!>                    M = 1 - J
!>                    DO 10, I = J, MIN( N, J + K )
!>                       A( M + I, J ) = matrix( I, J )
!>              10    CONTINUE
!>              20 CONTINUE
!>
!>           Note that the imaginary parts of the diagonal elements need
!>           not be set and are assumed to be zero.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>           On entry, LDA specifies the first dimension of A as declared
!>           in the calling (sub) program. LDA must be at least
!>           ( k + 1 ).
!> \endverbatim
!>
!> \param[in] X
!> \verbatim
!>          X is COMPLEX array, dimension at least
!>           ( 1 + ( n - 1 )*abs( INCX ) ).
!>           Before entry, the incremented array X must contain the
!>           vector x.
!> \endverbatim
!>
!> \param[in] INCX
!> \verbatim
!>          INCX is INTEGER
!>           On entry, INCX specifies the increment for the elements of
!>           X. INCX must not be zero.
!> \endverbatim
!>
!> \param[in] BETA
!> \verbatim
!>          BETA is COMPLEX
!>           On entry, BETA specifies the scalar beta.
!> \endverbatim
!>
!> \param[in,out] Y
!> \verbatim
!>          Y is COMPLEX array, dimension at least
!>           ( 1 + ( n - 1 )*abs( INCY ) ).
!>           Before entry, the incremented array Y must contain the
!>           vector y. On exit, Y is overwritten by the updated vector y.
!> \endverbatim
!>
!> \param[in] INCY
!> \verbatim
!>          INCY is INTEGER
!>           On entry, INCY specifies the increment for the elements of
!>           Y. INCY must not be zero.
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
!> \ingroup complex_blas_level2
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  Level 2 Blas routine.
!>  The vector and matrix arguments are not referenced when N = 0, or M = 0
!>
!>  -- Written on 22-October-1986.
!>     Jack Dongarra, Argonne National Lab.
!>     Jeremy Du Croz, Nag Central Office.
!>     Sven Hammarling, Nag Central Office.
!>     Richard Hanson, Sandia National Labs.
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE CHBMV(Uplo,N,K,Alpha,A,Lda,X,Incx,Beta,Y,Incy)
      USE S_LSAME
      USE S_XERBLA
      IMPLICIT NONE
!*--CHBMV193
!
! PARAMETER definitions rewritten by SPAG
!
      COMPLEX , PARAMETER  ::  ONE = (1.0E+0,0.0E+0) ,                  &
     &                         ZERO = (0.0E+0,0.0E+0)
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: K
      COMPLEX , INTENT(IN) :: Alpha
      COMPLEX , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      COMPLEX , INTENT(IN) , DIMENSION(*) :: X
      INTEGER , INTENT(IN) :: Incx
      COMPLEX , INTENT(IN) :: Beta
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Y
      INTEGER , INTENT(IN) :: Incy
!
! Local variable declarations rewritten by SPAG
!
      INTEGER :: i , info , ix , iy , j , jx , jy , kplus1 , kx , ky , l
      COMPLEX :: temp1 , temp2
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
         CALL XERBLA('CHBMV ',info)
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
                  temp2 = temp2 + CONJG(A(l+i,j))*X(i)
               ENDDO
               Y(j) = Y(j) + temp1*REAL(A(kplus1,j)) + Alpha*temp2
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
                  temp2 = temp2 + CONJG(A(l+i,j))*X(ix)
                  ix = ix + Incx
                  iy = iy + Incy
               ENDDO
               Y(jy) = Y(jy) + temp1*REAL(A(kplus1,j)) + Alpha*temp2
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
            Y(j) = Y(j) + temp1*REAL(A(1,j))
            l = 1 - j
            DO i = j + 1 , MIN(N,j+K)
               Y(i) = Y(i) + temp1*A(l+i,j)
               temp2 = temp2 + CONJG(A(l+i,j))*X(i)
            ENDDO
            Y(j) = Y(j) + Alpha*temp2
         ENDDO
      ELSE
         jx = kx
         jy = ky
         DO j = 1 , N
            temp1 = Alpha*X(jx)
            temp2 = ZERO
            Y(jy) = Y(jy) + temp1*REAL(A(1,j))
            l = 1 - j
            ix = jx
            iy = jy
            DO i = j + 1 , MIN(N,j+K)
               ix = ix + Incx
               iy = iy + Incy
               Y(iy) = Y(iy) + temp1*A(l+i,j)
               temp2 = temp2 + CONJG(A(l+i,j))*X(ix)
            ENDDO
            Y(jy) = Y(jy) + Alpha*temp2
            jx = jx + Incx
            jy = jy + Incy
         ENDDO
      ENDIF
!
!
!     End of CHBMV .
!
      END SUBROUTINE CHBMV
