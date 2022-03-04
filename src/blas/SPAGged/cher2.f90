!*==cher2.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b CHER2
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE CHER2(UPLO,N,ALPHA,X,INCX,Y,INCY,A,LDA)
!
!       .. Scalar Arguments ..
!       COMPLEX ALPHA
!       INTEGER INCX,INCY,LDA,N
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
!> CHER2  performs the hermitian rank 2 operation
!>
!>    A := alpha*x*y**H + conjg( alpha )*y*x**H + A,
!>
!> where alpha is a scalar, x and y are n element vectors and A is an n
!> by n hermitian matrix.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>           On entry, UPLO specifies whether the upper or lower
!>           triangular part of the array A is to be referenced as
!>           follows:
!>
!>              UPLO = 'U' or 'u'   Only the upper triangular part of A
!>                                  is to be referenced.
!>
!>              UPLO = 'L' or 'l'   Only the lower triangular part of A
!>                                  is to be referenced.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>           On entry, N specifies the order of the matrix A.
!>           N must be at least zero.
!> \endverbatim
!>
!> \param[in] ALPHA
!> \verbatim
!>          ALPHA is COMPLEX
!>           On entry, ALPHA specifies the scalar alpha.
!> \endverbatim
!>
!> \param[in] X
!> \verbatim
!>          X is COMPLEX array, dimension at least
!>           ( 1 + ( n - 1 )*abs( INCX ) ).
!>           Before entry, the incremented array X must contain the n
!>           element vector x.
!> \endverbatim
!>
!> \param[in] INCX
!> \verbatim
!>          INCX is INTEGER
!>           On entry, INCX specifies the increment for the elements of
!>           X. INCX must not be zero.
!> \endverbatim
!>
!> \param[in] Y
!> \verbatim
!>          Y is COMPLEX array, dimension at least
!>           ( 1 + ( n - 1 )*abs( INCY ) ).
!>           Before entry, the incremented array Y must contain the n
!>           element vector y.
!> \endverbatim
!>
!> \param[in] INCY
!> \verbatim
!>          INCY is INTEGER
!>           On entry, INCY specifies the increment for the elements of
!>           Y. INCY must not be zero.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is COMPLEX array, dimension ( LDA, N )
!>           Before entry with  UPLO = 'U' or 'u', the leading n by n
!>           upper triangular part of the array A must contain the upper
!>           triangular part of the hermitian matrix and the strictly
!>           lower triangular part of A is not referenced. On exit, the
!>           upper triangular part of the array A is overwritten by the
!>           upper triangular part of the updated matrix.
!>           Before entry with UPLO = 'L' or 'l', the leading n by n
!>           lower triangular part of the array A must contain the lower
!>           triangular part of the hermitian matrix and the strictly
!>           upper triangular part of A is not referenced. On exit, the
!>           lower triangular part of the array A is overwritten by the
!>           lower triangular part of the updated matrix.
!>           Note that the imaginary parts of the diagonal elements need
!>           not be set, they are assumed to be zero, and on exit they
!>           are set to zero.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>           On entry, LDA specifies the first dimension of A as declared
!>           in the calling (sub) program. LDA must be at least
!>           max( 1, n ).
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
!>
!>  -- Written on 22-October-1986.
!>     Jack Dongarra, Argonne National Lab.
!>     Jeremy Du Croz, Nag Central Office.
!>     Sven Hammarling, Nag Central Office.
!>     Richard Hanson, Sandia National Labs.
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE CHER2(Uplo,N,Alpha,X,Incx,Y,Incy,A,Lda)
      USE S_LSAME
      USE S_XERBLA
      IMPLICIT NONE
!*--CHER2156
!
! PARAMETER definitions rewritten by SPAG
!
      COMPLEX , PARAMETER  ::  ZERO = (0.0E+0,0.0E+0)
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(IN) :: Alpha
      COMPLEX , INTENT(IN) , DIMENSION(*) :: X
      INTEGER , INTENT(IN) :: Incx
      COMPLEX , INTENT(IN) , DIMENSION(*) :: Y
      INTEGER , INTENT(IN) :: Incy
      COMPLEX , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
!
! Local variable declarations rewritten by SPAG
!
      INTEGER :: i , info , ix , iy , j , jx , jy , kx , ky
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
      ELSEIF ( Incx==0 ) THEN
         info = 5
      ELSEIF ( Incy==0 ) THEN
         info = 7
      ELSEIF ( Lda<MAX(1,N) ) THEN
         info = 9
      ENDIF
      IF ( info/=0 ) THEN
         CALL XERBLA('CHER2 ',info)
         RETURN
      ENDIF
!
!     Quick return if possible.
!
      IF ( (N==0) .OR. (Alpha==ZERO) ) RETURN
!
!     Set up the start points in X and Y if the increments are not both
!     unity.
!
      IF ( (Incx/=1) .OR. (Incy/=1) ) THEN
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
         jx = kx
         jy = ky
      ENDIF
!
!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through the triangular part
!     of A.
!
      IF ( LSAME(Uplo,'U') ) THEN
!
!        Form  A  when A is stored in the upper triangle.
!
         IF ( (Incx==1) .AND. (Incy==1) ) THEN
            DO j = 1 , N
               IF ( (X(j)/=ZERO) .OR. (Y(j)/=ZERO) ) THEN
                  temp1 = Alpha*CONJG(Y(j))
                  temp2 = CONJG(Alpha*X(j))
                  DO i = 1 , j - 1
                     A(i,j) = A(i,j) + X(i)*temp1 + Y(i)*temp2
                  ENDDO
                  A(j,j) = REAL(A(j,j)) + REAL(X(j)*temp1+Y(j)*temp2)
               ELSE
                  A(j,j) = REAL(A(j,j))
               ENDIF
            ENDDO
         ELSE
            DO j = 1 , N
               IF ( (X(jx)/=ZERO) .OR. (Y(jy)/=ZERO) ) THEN
                  temp1 = Alpha*CONJG(Y(jy))
                  temp2 = CONJG(Alpha*X(jx))
                  ix = kx
                  iy = ky
                  DO i = 1 , j - 1
                     A(i,j) = A(i,j) + X(ix)*temp1 + Y(iy)*temp2
                     ix = ix + Incx
                     iy = iy + Incy
                  ENDDO
                  A(j,j) = REAL(A(j,j)) + REAL(X(jx)*temp1+Y(jy)*temp2)
               ELSE
                  A(j,j) = REAL(A(j,j))
               ENDIF
               jx = jx + Incx
               jy = jy + Incy
            ENDDO
         ENDIF
!
!        Form  A  when A is stored in the lower triangle.
!
      ELSEIF ( (Incx==1) .AND. (Incy==1) ) THEN
         DO j = 1 , N
            IF ( (X(j)/=ZERO) .OR. (Y(j)/=ZERO) ) THEN
               temp1 = Alpha*CONJG(Y(j))
               temp2 = CONJG(Alpha*X(j))
               A(j,j) = REAL(A(j,j)) + REAL(X(j)*temp1+Y(j)*temp2)
               DO i = j + 1 , N
                  A(i,j) = A(i,j) + X(i)*temp1 + Y(i)*temp2
               ENDDO
            ELSE
               A(j,j) = REAL(A(j,j))
            ENDIF
         ENDDO
      ELSE
         DO j = 1 , N
            IF ( (X(jx)/=ZERO) .OR. (Y(jy)/=ZERO) ) THEN
               temp1 = Alpha*CONJG(Y(jy))
               temp2 = CONJG(Alpha*X(jx))
               A(j,j) = REAL(A(j,j)) + REAL(X(jx)*temp1+Y(jy)*temp2)
               ix = jx
               iy = jy
               DO i = j + 1 , N
                  ix = ix + Incx
                  iy = iy + Incy
                  A(i,j) = A(i,j) + X(ix)*temp1 + Y(iy)*temp2
               ENDDO
            ELSE
               A(j,j) = REAL(A(j,j))
            ENDIF
            jx = jx + Incx
            jy = jy + Incy
         ENDDO
      ENDIF
!
!
!     End of CHER2 .
!
      END SUBROUTINE CHER2
