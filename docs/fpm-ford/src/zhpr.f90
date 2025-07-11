!*==zhpr.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b ZHPR
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZHPR(UPLO,N,ALPHA,X,INCX,AP)
!
!       .. Scalar Arguments ..
!       DOUBLE PRECISION ALPHA
!       INTEGER INCX,N
!       CHARACTER UPLO
!       ..
!       .. Array Arguments ..
!       COMPLEX*16 AP(*),X(*)
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZHPR    performs the hermitian rank 1 operation
!>
!>    A := alpha*x*x**H + A,
!>
!> where alpha is a real scalar, x is an n element vector and A is an
!> n by n hermitian matrix, supplied in packed form.
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
!>          X is COMPLEX*16 array, dimension at least
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
!> \param[in,out] AP
!> \verbatim
!>          AP is COMPLEX*16 array, dimension at least
!>           ( ( n*( n + 1 ) )/2 ).
!>           Before entry with  UPLO = 'U' or 'u', the array AP must
!>           contain the upper triangular part of the hermitian matrix
!>           packed sequentially, column by column, so that AP( 1 )
!>           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 1, 2 )
!>           and a( 2, 2 ) respectively, and so on. On exit, the array
!>           AP is overwritten by the upper triangular part of the
!>           updated matrix.
!>           Before entry with UPLO = 'L' or 'l', the array AP must
!>           contain the lower triangular part of the hermitian matrix
!>           packed sequentially, column by column, so that AP( 1 )
!>           contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 2, 1 )
!>           and a( 3, 1 ) respectively, and so on. On exit, the array
!>           AP is overwritten by the lower triangular part of the
!>           updated matrix.
!>           Note that the imaginary parts of the diagonal elements need
!>           not be set, they are assumed to be zero, and on exit they
!>           are set to zero.
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
!> \ingroup complex16_blas_level2
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
      SUBROUTINE ZHPR(Uplo,N,Alpha,X,Incx,Ap)
      IMPLICIT NONE
!*--ZHPR134
!
!  -- Reference BLAS level2 routine (version 3.7.0) --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      DOUBLE PRECISION Alpha
      INTEGER Incx , N
      CHARACTER Uplo
!     ..
!     .. Array Arguments ..
      COMPLEX*16 Ap(*) , X(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      COMPLEX*16 ZERO
      PARAMETER (ZERO=(0.0D+0,0.0D+0))
!     ..
!     .. Local Scalars ..
      COMPLEX*16 temp
      INTEGER i , info , ix , j , jx , k , kk , kx
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC DBLE , DCONJG
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
      ENDIF
      IF ( info/=0 ) THEN
         CALL XERBLA('ZHPR  ',info)
         RETURN
      ENDIF
!
!     Quick return if possible.
!
      IF ( (N==0) .OR. (Alpha==DBLE(ZERO)) ) RETURN
!
!     Set the start point in X if the increment is not unity.
!
      IF ( Incx<=0 ) THEN
         kx = 1 - (N-1)*Incx
      ELSEIF ( Incx/=1 ) THEN
         kx = 1
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
         IF ( Incx==1 ) THEN
            DO j = 1 , N
               IF ( X(j)/=ZERO ) THEN
                  temp = Alpha*DCONJG(X(j))
                  k = kk
                  DO i = 1 , j - 1
                     Ap(k) = Ap(k) + X(i)*temp
                     k = k + 1
                  ENDDO
                  Ap(kk+j-1) = DBLE(Ap(kk+j-1)) + DBLE(X(j)*temp)
               ELSE
                  Ap(kk+j-1) = DBLE(Ap(kk+j-1))
               ENDIF
               kk = kk + j
            ENDDO
         ELSE
            jx = kx
            DO j = 1 , N
               IF ( X(jx)/=ZERO ) THEN
                  temp = Alpha*DCONJG(X(jx))
                  ix = kx
                  DO k = kk , kk + j - 2
                     Ap(k) = Ap(k) + X(ix)*temp
                     ix = ix + Incx
                  ENDDO
                  Ap(kk+j-1) = DBLE(Ap(kk+j-1)) + DBLE(X(jx)*temp)
               ELSE
                  Ap(kk+j-1) = DBLE(Ap(kk+j-1))
               ENDIF
               jx = jx + Incx
               kk = kk + j
            ENDDO
         ENDIF
!
!        Form  A  when lower triangle is stored in AP.
!
      ELSEIF ( Incx==1 ) THEN
         DO j = 1 , N
            IF ( X(j)/=ZERO ) THEN
               temp = Alpha*DCONJG(X(j))
               Ap(kk) = DBLE(Ap(kk)) + DBLE(temp*X(j))
               k = kk + 1
               DO i = j + 1 , N
                  Ap(k) = Ap(k) + X(i)*temp
                  k = k + 1
               ENDDO
            ELSE
               Ap(kk) = DBLE(Ap(kk))
            ENDIF
            kk = kk + N - j + 1
         ENDDO
      ELSE
         jx = kx
         DO j = 1 , N
            IF ( X(jx)/=ZERO ) THEN
               temp = Alpha*DCONJG(X(jx))
               Ap(kk) = DBLE(Ap(kk)) + DBLE(temp*X(jx))
               ix = jx
               DO k = kk + 1 , kk + N - j
                  ix = ix + Incx
                  Ap(k) = Ap(k) + X(ix)*temp
               ENDDO
            ELSE
               Ap(kk) = DBLE(Ap(kk))
            ENDIF
            jx = jx + Incx
            kk = kk + N - j + 1
         ENDDO
      ENDIF
!
!
!     End of ZHPR  .
!
      END SUBROUTINE ZHPR
