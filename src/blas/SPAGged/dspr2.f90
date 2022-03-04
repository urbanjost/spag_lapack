!*==dspr2.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b DSPR2
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DSPR2(UPLO,N,ALPHA,X,INCX,Y,INCY,AP)
!
!       .. Scalar Arguments ..
!       DOUBLE PRECISION ALPHA
!       INTEGER INCX,INCY,N
!       CHARACTER UPLO
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION AP(*),X(*),Y(*)
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DSPR2  performs the symmetric rank 2 operation
!>
!>    A := alpha*x*y**T + alpha*y*x**T + A,
!>
!> where alpha is a scalar, x and y are n element vectors and A is an
!> n by n symmetric matrix, supplied in packed form.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>           On entry, UPLO specifies whether the upper or lower
!>           triangular part of the matrix A is supplied in the packed
!>           array AP as follows:
!>
!>              UPLO = 'U' or 'u'   The upper triangular part of A is
!>                                  supplied in AP.
!>
!>              UPLO = 'L' or 'l'   The lower triangular part of A is
!>                                  supplied in AP.
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
!>          ALPHA is DOUBLE PRECISION.
!>           On entry, ALPHA specifies the scalar alpha.
!> \endverbatim
!>
!> \param[in] X
!> \verbatim
!>          X is DOUBLE PRECISION array, dimension at least
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
!>          Y is DOUBLE PRECISION array, dimension at least
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
!> \param[in,out] AP
!> \verbatim
!>          AP is DOUBLE PRECISION array, dimension at least
!>           ( ( n*( n + 1 ) )/2 ).
!>           Before entry with  UPLO = 'U' or 'u', the array AP must
!>           contain the upper triangular part of the symmetric matrix
!>           packed sequentially, column by column, so that AP( 1 )
!>           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 1, 2 )
!>           and a( 2, 2 ) respectively, and so on. On exit, the array
!>           AP is overwritten by the upper triangular part of the
!>           updated matrix.
!>           Before entry with UPLO = 'L' or 'l', the array AP must
!>           contain the lower triangular part of the symmetric matrix
!>           packed sequentially, column by column, so that AP( 1 )
!>           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 2, 1 )
!>           and a( 3, 1 ) respectively, and so on. On exit, the array
!>           AP is overwritten by the lower triangular part of the
!>           updated matrix.
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
!> \ingroup double_blas_level2
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
      SUBROUTINE DSPR2(Uplo,N,Alpha,X,Incx,Y,Incy,Ap)
      USE F77KINDS                        
      USE S_LSAME
      USE S_XERBLA
      IMPLICIT NONE
!*--DSPR2149
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Uplo
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(IN) :: Alpha
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: X
      INTEGER , INTENT(IN) :: Incx
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: Y
      INTEGER , INTENT(IN) :: Incy
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Ap
!
! Local variable declarations rewritten by SPAG
!
      INTEGER :: i , info , ix , iy , j , jx , jy , k , kk , kx , ky
      REAL(R8KIND) :: temp1 , temp2
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
      ENDIF
      IF ( info/=0 ) THEN
         CALL XERBLA('DSPR2 ',info)
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
!     Start the operations. In this version the elements of the array AP
!     are accessed sequentially with one pass through AP.
!
      kk = 1
      IF ( LSAME(Uplo,'U') ) THEN
!
!        Form  A  when upper triangle is stored in AP.
!
         IF ( (Incx==1) .AND. (Incy==1) ) THEN
            DO j = 1 , N
               IF ( (X(j)/=ZERO) .OR. (Y(j)/=ZERO) ) THEN
                  temp1 = Alpha*Y(j)
                  temp2 = Alpha*X(j)
                  k = kk
                  DO i = 1 , j
                     Ap(k) = Ap(k) + X(i)*temp1 + Y(i)*temp2
                     k = k + 1
                  ENDDO
               ENDIF
               kk = kk + j
            ENDDO
         ELSE
            DO j = 1 , N
               IF ( (X(jx)/=ZERO) .OR. (Y(jy)/=ZERO) ) THEN
                  temp1 = Alpha*Y(jy)
                  temp2 = Alpha*X(jx)
                  ix = kx
                  iy = ky
                  DO k = kk , kk + j - 1
                     Ap(k) = Ap(k) + X(ix)*temp1 + Y(iy)*temp2
                     ix = ix + Incx
                     iy = iy + Incy
                  ENDDO
               ENDIF
               jx = jx + Incx
               jy = jy + Incy
               kk = kk + j
            ENDDO
         ENDIF
!
!        Form  A  when lower triangle is stored in AP.
!
      ELSEIF ( (Incx==1) .AND. (Incy==1) ) THEN
         DO j = 1 , N
            IF ( (X(j)/=ZERO) .OR. (Y(j)/=ZERO) ) THEN
               temp1 = Alpha*Y(j)
               temp2 = Alpha*X(j)
               k = kk
               DO i = j , N
                  Ap(k) = Ap(k) + X(i)*temp1 + Y(i)*temp2
                  k = k + 1
               ENDDO
            ENDIF
            kk = kk + N - j + 1
         ENDDO
      ELSE
         DO j = 1 , N
            IF ( (X(jx)/=ZERO) .OR. (Y(jy)/=ZERO) ) THEN
               temp1 = Alpha*Y(jy)
               temp2 = Alpha*X(jx)
               ix = jx
               iy = jy
               DO k = kk , kk + N - j
                  Ap(k) = Ap(k) + X(ix)*temp1 + Y(iy)*temp2
                  ix = ix + Incx
                  iy = iy + Incy
               ENDDO
            ENDIF
            jx = jx + Incx
            jy = jy + Incy
            kk = kk + N - j + 1
         ENDDO
      ENDIF
!
!
!     End of DSPR2 .
!
      END SUBROUTINE DSPR2
