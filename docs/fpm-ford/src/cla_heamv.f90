!*==cla_heamv.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b CLA_HEAMV computes a matrix-vector product using a Hermitian indefinite matrix to calculate error bounds.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CLA_HEAMV + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cla_heamv.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cla_heamv.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cla_heamv.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CLA_HEAMV( UPLO, N, ALPHA, A, LDA, X, INCX, BETA, Y,
!                             INCY )
!
!       .. Scalar Arguments ..
!       REAL               ALPHA, BETA
!       INTEGER            INCX, INCY, LDA, N, UPLO
!       ..
!       .. Array Arguments ..
!       COMPLEX            A( LDA, * ), X( * )
!       REAL               Y( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CLA_SYAMV  performs the matrix-vector operation
!>
!>         y := alpha*abs(A)*abs(x) + beta*abs(y),
!>
!> where alpha and beta are scalars, x and y are vectors and A is an
!> n by n symmetric matrix.
!>
!> This function is primarily used in calculating error bounds.
!> To protect against underflow during evaluation, components in
!> the resulting vector are perturbed away from zero by (N+1)
!> times the underflow threshold.  To prevent unnecessarily large
!> errors for block-structure embedded in general matrices,
!> "symbolically" zero components are not perturbed.  A zero
!> entry is considered "symbolic" if all multiplications involved
!> in computing that entry have at least one zero multiplicand.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is INTEGER
!>           On entry, UPLO specifies whether the upper or lower
!>           triangular part of the array A is to be referenced as
!>           follows:
!>
!>              UPLO = BLAS_UPPER   Only the upper triangular part of A
!>                                  is to be referenced.
!>
!>              UPLO = BLAS_LOWER   Only the lower triangular part of A
!>                                  is to be referenced.
!>
!>           Unchanged on exit.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>           On entry, N specifies the number of columns of the matrix A.
!>           N must be at least zero.
!>           Unchanged on exit.
!> \endverbatim
!>
!> \param[in] ALPHA
!> \verbatim
!>          ALPHA is REAL .
!>           On entry, ALPHA specifies the scalar alpha.
!>           Unchanged on exit.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is COMPLEX array, dimension ( LDA, n ).
!>           Before entry, the leading m by n part of the array A must
!>           contain the matrix of coefficients.
!>           Unchanged on exit.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>           On entry, LDA specifies the first dimension of A as declared
!>           in the calling (sub) program. LDA must be at least
!>           max( 1, n ).
!>           Unchanged on exit.
!> \endverbatim
!>
!> \param[in] X
!> \verbatim
!>          X is COMPLEX array, dimension
!>           ( 1 + ( n - 1 )*abs( INCX ) )
!>           Before entry, the incremented array X must contain the
!>           vector x.
!>           Unchanged on exit.
!> \endverbatim
!>
!> \param[in] INCX
!> \verbatim
!>          INCX is INTEGER
!>           On entry, INCX specifies the increment for the elements of
!>           X. INCX must not be zero.
!>           Unchanged on exit.
!> \endverbatim
!>
!> \param[in] BETA
!> \verbatim
!>          BETA is REAL .
!>           On entry, BETA specifies the scalar beta. When BETA is
!>           supplied as zero then Y need not be set on input.
!>           Unchanged on exit.
!> \endverbatim
!>
!> \param[in,out] Y
!> \verbatim
!>          Y is REAL array, dimension
!>           ( 1 + ( n - 1 )*abs( INCY ) )
!>           Before entry with BETA non-zero, the incremented array Y
!>           must contain the vector y. On exit, Y is overwritten by the
!>           updated vector y.
!> \endverbatim
!>
!> \param[in] INCY
!> \verbatim
!>          INCY is INTEGER
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
!> \date June 2017
!
!> \ingroup complexHEcomputational
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
!>  -- Modified for the absolute-value product, April 2006
!>     Jason Riedy, UC Berkeley
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE CLA_HEAMV(Uplo,N,Alpha,A,Lda,X,Incx,Beta,Y,Incy)
      IMPLICIT NONE
!*--CLA_HEAMV181
!
!  -- LAPACK computational routine (version 3.7.1) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     June 2017
!
!     .. Scalar Arguments ..
      REAL Alpha , Beta
      INTEGER Incx , Incy , Lda , N , Uplo
!     ..
!     .. Array Arguments ..
      COMPLEX A(Lda,*) , X(*)
      REAL Y(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ONE , ZERO
      PARAMETER (ONE=1.0E+0,ZERO=0.0E+0)
!     ..
!     .. Local Scalars ..
      LOGICAL symb_zero
      REAL temp , safe1
      INTEGER i , info , iy , j , jx , kx , ky
      COMPLEX zdum
!     ..
!     .. External Subroutines ..
      EXTERNAL XERBLA , SLAMCH
      REAL SLAMCH
!     ..
!     .. External Functions ..
      EXTERNAL ILAUPLO
      INTEGER ILAUPLO
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX , ABS , SIGN , REAL , AIMAG
!     ..
!     .. Statement Functions ..
      REAL CABS1
!     ..
!     .. Statement Function Definitions ..
      CABS1(zdum) = ABS(REAL(zdum)) + ABS(AIMAG(zdum))
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      info = 0
      IF ( Uplo/=ILAUPLO('U') .AND. Uplo/=ILAUPLO('L') ) THEN
         info = 1
      ELSEIF ( N<0 ) THEN
         info = 2
      ELSEIF ( Lda<MAX(1,N) ) THEN
         info = 5
      ELSEIF ( Incx==0 ) THEN
         info = 7
      ELSEIF ( Incy==0 ) THEN
         info = 10
      ENDIF
      IF ( info/=0 ) THEN
         CALL XERBLA('CHEMV ',info)
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
!     Set SAFE1 essentially to be the underflow threshold times the
!     number of additions in each row.
!
      safe1 = SLAMCH('Safe minimum')
      safe1 = (N+1)*safe1
!
!     Form  y := alpha*abs(A)*abs(x) + beta*abs(y).
!
!     The O(N^2) SYMB_ZERO tests could be replaced by O(N) queries to
!     the inexact flag.  Still doesn't help change the iteration order
!     to per-column.
!
      iy = ky
      IF ( Incx==1 ) THEN
         IF ( Uplo==ILAUPLO('U') ) THEN
            DO i = 1 , N
               IF ( Beta==ZERO ) THEN
                  symb_zero = .TRUE.
                  Y(iy) = 0.0
               ELSEIF ( Y(iy)==ZERO ) THEN
                  symb_zero = .TRUE.
               ELSE
                  symb_zero = .FALSE.
                  Y(iy) = Beta*ABS(Y(iy))
               ENDIF
               IF ( Alpha/=ZERO ) THEN
                  DO j = 1 , i
                     temp = CABS1(A(j,i))
                     symb_zero = symb_zero .AND.                        &
     &                           (X(j)==ZERO .OR. temp==ZERO)
 
                     Y(iy) = Y(iy) + Alpha*CABS1(X(j))*temp
                  ENDDO
                  DO j = i + 1 , N
                     temp = CABS1(A(i,j))
                     symb_zero = symb_zero .AND.                        &
     &                           (X(j)==ZERO .OR. temp==ZERO)
 
                     Y(iy) = Y(iy) + Alpha*CABS1(X(j))*temp
                  ENDDO
               ENDIF
 
               IF ( .NOT.symb_zero ) Y(iy) = Y(iy) + SIGN(safe1,Y(iy))
 
               iy = iy + Incy
            ENDDO
         ELSE
            DO i = 1 , N
               IF ( Beta==ZERO ) THEN
                  symb_zero = .TRUE.
                  Y(iy) = 0.0
               ELSEIF ( Y(iy)==ZERO ) THEN
                  symb_zero = .TRUE.
               ELSE
                  symb_zero = .FALSE.
                  Y(iy) = Beta*ABS(Y(iy))
               ENDIF
               IF ( Alpha/=ZERO ) THEN
                  DO j = 1 , i
                     temp = CABS1(A(i,j))
                     symb_zero = symb_zero .AND.                        &
     &                           (X(j)==ZERO .OR. temp==ZERO)
 
                     Y(iy) = Y(iy) + Alpha*CABS1(X(j))*temp
                  ENDDO
                  DO j = i + 1 , N
                     temp = CABS1(A(j,i))
                     symb_zero = symb_zero .AND.                        &
     &                           (X(j)==ZERO .OR. temp==ZERO)
 
                     Y(iy) = Y(iy) + Alpha*CABS1(X(j))*temp
                  ENDDO
               ENDIF
 
               IF ( .NOT.symb_zero ) Y(iy) = Y(iy) + SIGN(safe1,Y(iy))
 
               iy = iy + Incy
            ENDDO
         ENDIF
      ELSEIF ( Uplo==ILAUPLO('U') ) THEN
         DO i = 1 , N
            IF ( Beta==ZERO ) THEN
               symb_zero = .TRUE.
               Y(iy) = 0.0
            ELSEIF ( Y(iy)==ZERO ) THEN
               symb_zero = .TRUE.
            ELSE
               symb_zero = .FALSE.
               Y(iy) = Beta*ABS(Y(iy))
            ENDIF
            jx = kx
            IF ( Alpha/=ZERO ) THEN
               DO j = 1 , i
                  temp = CABS1(A(j,i))
                  symb_zero = symb_zero .AND.                           &
     &                        (X(j)==ZERO .OR. temp==ZERO)
 
                  Y(iy) = Y(iy) + Alpha*CABS1(X(jx))*temp
                  jx = jx + Incx
               ENDDO
               DO j = i + 1 , N
                  temp = CABS1(A(i,j))
                  symb_zero = symb_zero .AND.                           &
     &                        (X(j)==ZERO .OR. temp==ZERO)
 
                  Y(iy) = Y(iy) + Alpha*CABS1(X(jx))*temp
                  jx = jx + Incx
               ENDDO
            ENDIF
 
            IF ( .NOT.symb_zero ) Y(iy) = Y(iy) + SIGN(safe1,Y(iy))
 
            iy = iy + Incy
         ENDDO
      ELSE
         DO i = 1 , N
            IF ( Beta==ZERO ) THEN
               symb_zero = .TRUE.
               Y(iy) = 0.0
            ELSEIF ( Y(iy)==ZERO ) THEN
               symb_zero = .TRUE.
            ELSE
               symb_zero = .FALSE.
               Y(iy) = Beta*ABS(Y(iy))
            ENDIF
            jx = kx
            IF ( Alpha/=ZERO ) THEN
               DO j = 1 , i
                  temp = CABS1(A(i,j))
                  symb_zero = symb_zero .AND.                           &
     &                        (X(j)==ZERO .OR. temp==ZERO)
 
                  Y(iy) = Y(iy) + Alpha*CABS1(X(jx))*temp
                  jx = jx + Incx
               ENDDO
               DO j = i + 1 , N
                  temp = CABS1(A(j,i))
                  symb_zero = symb_zero .AND.                           &
     &                        (X(j)==ZERO .OR. temp==ZERO)
 
                  Y(iy) = Y(iy) + Alpha*CABS1(X(jx))*temp
                  jx = jx + Incx
               ENDDO
            ENDIF
 
            IF ( .NOT.symb_zero ) Y(iy) = Y(iy) + SIGN(safe1,Y(iy))
 
            iy = iy + Incy
         ENDDO
 
      ENDIF
!
!
!     End of CLA_HEAMV
!
      END SUBROUTINE CLA_HEAMV
