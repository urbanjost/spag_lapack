!*==cher.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b CHER
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE CHER(UPLO,N,ALPHA,X,INCX,A,LDA)
!
!       .. Scalar Arguments ..
!       REAL ALPHA
!       INTEGER INCX,LDA,N
!       CHARACTER UPLO
!       ..
!       .. Array Arguments ..
!       COMPLEX A(LDA,*),X(*)
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CHER   performs the hermitian rank 1 operation
!>
!>    A := alpha*x*x**H + A,
!>
!> where alpha is a real scalar, x is an n element vector and A is an
!> n by n hermitian matrix.
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
!>          ALPHA is REAL
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
      SUBROUTINE CHER(Uplo,N,Alpha,X,Incx,A,Lda)
      USE S_LSAME
      USE S_XERBLA
      IMPLICIT NONE
!*--CHER141
!
! PARAMETER definitions rewritten by SPAG
!
      COMPLEX , PARAMETER  ::  ZERO = (0.0E+0,0.0E+0)
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(IN) :: Alpha
      COMPLEX , INTENT(IN) , DIMENSION(*) :: X
      INTEGER , INTENT(IN) :: Incx
      COMPLEX , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
!
! Local variable declarations rewritten by SPAG
!
      INTEGER :: i , info , ix , j , jx , kx
      COMPLEX :: temp
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
      ELSEIF ( Lda<MAX(1,N) ) THEN
         info = 7
      ENDIF
      IF ( info/=0 ) THEN
         CALL XERBLA('CHER  ',info)
         RETURN
      ENDIF
!
!     Quick return if possible.
!
      IF ( (N==0) .OR. (Alpha==REAL(ZERO)) ) RETURN
!
!     Set the start point in X if the increment is not unity.
!
      IF ( Incx<=0 ) THEN
         kx = 1 - (N-1)*Incx
      ELSEIF ( Incx/=1 ) THEN
         kx = 1
      ENDIF
!
!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through the triangular part
!     of A.
!
      IF ( LSAME(Uplo,'U') ) THEN
!
!        Form  A  when A is stored in upper triangle.
!
         IF ( Incx==1 ) THEN
            DO j = 1 , N
               IF ( X(j)/=ZERO ) THEN
                  temp = Alpha*CONJG(X(j))
                  DO i = 1 , j - 1
                     A(i,j) = A(i,j) + X(i)*temp
                  ENDDO
                  A(j,j) = REAL(A(j,j)) + REAL(X(j)*temp)
               ELSE
                  A(j,j) = REAL(A(j,j))
               ENDIF
            ENDDO
         ELSE
            jx = kx
            DO j = 1 , N
               IF ( X(jx)/=ZERO ) THEN
                  temp = Alpha*CONJG(X(jx))
                  ix = kx
                  DO i = 1 , j - 1
                     A(i,j) = A(i,j) + X(ix)*temp
                     ix = ix + Incx
                  ENDDO
                  A(j,j) = REAL(A(j,j)) + REAL(X(jx)*temp)
               ELSE
                  A(j,j) = REAL(A(j,j))
               ENDIF
               jx = jx + Incx
            ENDDO
         ENDIF
!
!        Form  A  when A is stored in lower triangle.
!
      ELSEIF ( Incx==1 ) THEN
         DO j = 1 , N
            IF ( X(j)/=ZERO ) THEN
               temp = Alpha*CONJG(X(j))
               A(j,j) = REAL(A(j,j)) + REAL(temp*X(j))
               DO i = j + 1 , N
                  A(i,j) = A(i,j) + X(i)*temp
               ENDDO
            ELSE
               A(j,j) = REAL(A(j,j))
            ENDIF
         ENDDO
      ELSE
         jx = kx
         DO j = 1 , N
            IF ( X(jx)/=ZERO ) THEN
               temp = Alpha*CONJG(X(jx))
               A(j,j) = REAL(A(j,j)) + REAL(temp*X(jx))
               ix = jx
               DO i = j + 1 , N
                  ix = ix + Incx
                  A(i,j) = A(i,j) + X(ix)*temp
               ENDDO
            ELSE
               A(j,j) = REAL(A(j,j))
            ENDIF
            jx = jx + Incx
         ENDDO
      ENDIF
!
!
!     End of CHER  .
!
      END SUBROUTINE CHER
