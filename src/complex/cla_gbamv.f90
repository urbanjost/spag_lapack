!*==cla_gbamv.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b CLA_GBAMV performs a matrix-vector operation to calculate error bounds.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CLA_GBAMV + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cla_gbamv.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cla_gbamv.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cla_gbamv.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CLA_GBAMV( TRANS, M, N, KL, KU, ALPHA, AB, LDAB, X,
!                             INCX, BETA, Y, INCY )
!
!       .. Scalar Arguments ..
!       REAL               ALPHA, BETA
!       INTEGER            INCX, INCY, LDAB, M, N, KL, KU, TRANS
!       ..
!       .. Array Arguments ..
!       COMPLEX            AB( LDAB, * ), X( * )
!       REAL               Y( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CLA_GBAMV  performs one of the matrix-vector operations
!>
!>         y := alpha*abs(A)*abs(x) + beta*abs(y),
!>    or   y := alpha*abs(A)**T*abs(x) + beta*abs(y),
!>
!> where alpha and beta are scalars, x and y are vectors and A is an
!> m by n matrix.
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
!> \param[in] TRANS
!> \verbatim
!>          TRANS is INTEGER
!>           On entry, TRANS specifies the operation to be performed as
!>           follows:
!>
!>             BLAS_NO_TRANS      y := alpha*abs(A)*abs(x) + beta*abs(y)
!>             BLAS_TRANS         y := alpha*abs(A**T)*abs(x) + beta*abs(y)
!>             BLAS_CONJ_TRANS    y := alpha*abs(A**T)*abs(x) + beta*abs(y)
!>
!>           Unchanged on exit.
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>           On entry, M specifies the number of rows of the matrix A.
!>           M must be at least zero.
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
!> \param[in] KL
!> \verbatim
!>          KL is INTEGER
!>           The number of subdiagonals within the band of A.  KL >= 0.
!> \endverbatim
!>
!> \param[in] KU
!> \verbatim
!>          KU is INTEGER
!>           The number of superdiagonals within the band of A.  KU >= 0.
!> \endverbatim
!>
!> \param[in] ALPHA
!> \verbatim
!>          ALPHA is REAL
!>           On entry, ALPHA specifies the scalar alpha.
!>           Unchanged on exit.
!> \endverbatim
!>
!> \param[in] AB
!> \verbatim
!>          AB is COMPLEX array, dimension (LDAB,n)
!>           Before entry, the leading m by n part of the array AB must
!>           contain the matrix of coefficients.
!>           Unchanged on exit.
!> \endverbatim
!>
!> \param[in] LDAB
!> \verbatim
!>          LDAB is INTEGER
!>           On entry, LDAB specifies the first dimension of AB as declared
!>           in the calling (sub) program. LDAB must be at least
!>           max( 1, m ).
!>           Unchanged on exit.
!> \endverbatim
!>
!> \param[in] X
!> \verbatim
!>          X is COMPLEX array, dimension
!>           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
!>           and at least
!>           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
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
!>          BETA is REAL
!>           On entry, BETA specifies the scalar beta. When BETA is
!>           supplied as zero then Y need not be set on input.
!>           Unchanged on exit.
!> \endverbatim
!>
!> \param[in,out] Y
!> \verbatim
!>          Y is REAL array, dimension
!>           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
!>           and at least
!>           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
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
!>
!>  Level 2 Blas routine.
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
!> \date June 2016
!
!> \ingroup complexGBcomputational
!
!  =====================================================================
      SUBROUTINE CLA_GBAMV(Trans,M,N,Kl,Ku,Alpha,Ab,Ldab,X,Incx,Beta,Y, &
     &                     Incy)
      IMPLICIT NONE
!*--CLA_GBAMV190
!
!  -- LAPACK computational routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     June 2016
!
!     .. Scalar Arguments ..
      REAL Alpha , Beta
      INTEGER Incx , Incy , Ldab , M , N , Kl , Ku , Trans
!     ..
!     .. Array Arguments ..
      COMPLEX Ab(Ldab,*) , X(*)
      REAL Y(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      COMPLEX ONE , ZERO
      PARAMETER (ONE=1.0E+0,ZERO=0.0E+0)
!     ..
!     .. Local Scalars ..
      LOGICAL symb_zero
      REAL temp , safe1
      INTEGER i , info , iy , j , jx , kx , ky , lenx , leny , kd , ke
      COMPLEX cdum
!     ..
!     .. External Subroutines ..
      EXTERNAL XERBLA , SLAMCH
      REAL SLAMCH
!     ..
!     .. External Functions ..
      EXTERNAL ILATRANS
      INTEGER ILATRANS
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX , ABS , REAL , AIMAG , SIGN
!     ..
!     .. Statement Functions
      REAL CABS1
!     ..
!     .. Statement Function Definitions ..
      CABS1(cdum) = ABS(REAL(cdum)) + ABS(AIMAG(cdum))
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      info = 0
      IF ( (Trans/=ILATRANS('N')) .AND. (Trans/=ILATRANS('T')) .AND.    &
     &     (Trans/=ILATRANS('C')) ) THEN
         info = 1
      ELSEIF ( M<0 ) THEN
         info = 2
      ELSEIF ( N<0 ) THEN
         info = 3
      ELSEIF ( Kl<0 .OR. Kl>M-1 ) THEN
         info = 4
      ELSEIF ( Ku<0 .OR. Ku>N-1 ) THEN
         info = 5
      ELSEIF ( Ldab<Kl+Ku+1 ) THEN
         info = 6
      ELSEIF ( Incx==0 ) THEN
         info = 8
      ELSEIF ( Incy==0 ) THEN
         info = 11
      ENDIF
      IF ( info/=0 ) THEN
         CALL XERBLA('CLA_GBAMV ',info)
         RETURN
      ENDIF
!
!     Quick return if possible.
!
      IF ( (M==0) .OR. (N==0) .OR. ((Alpha==ZERO) .AND. (Beta==ONE)) )  &
     &     RETURN
!
!     Set  LENX  and  LENY, the lengths of the vectors x and y, and set
!     up the start points in  X  and  Y.
!
      IF ( Trans==ILATRANS('N') ) THEN
         lenx = N
         leny = M
      ELSE
         lenx = M
         leny = N
      ENDIF
      IF ( Incx>0 ) THEN
         kx = 1
      ELSE
         kx = 1 - (lenx-1)*Incx
      ENDIF
      IF ( Incy>0 ) THEN
         ky = 1
      ELSE
         ky = 1 - (leny-1)*Incy
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
!     The O(M*N) SYMB_ZERO tests could be replaced by O(N) queries to
!     the inexact flag.  Still doesn't help change the iteration order
!     to per-column.
!
      kd = Ku + 1
      ke = Kl + 1
      iy = ky
      IF ( Incx==1 ) THEN
         IF ( Trans==ILATRANS('N') ) THEN
            DO i = 1 , leny
               IF ( Beta==0.0 ) THEN
                  symb_zero = .TRUE.
                  Y(iy) = 0.0
               ELSEIF ( Y(iy)==0.0 ) THEN
                  symb_zero = .TRUE.
               ELSE
                  symb_zero = .FALSE.
                  Y(iy) = Beta*ABS(Y(iy))
               ENDIF
               IF ( Alpha/=0.0 ) THEN
                  DO j = MAX(i-Kl,1) , MIN(i+Ku,lenx)
                     temp = CABS1(Ab(kd+i-j,j))
                     symb_zero = symb_zero .AND.                        &
     &                           (X(j)==ZERO .OR. temp==ZERO)
 
                     Y(iy) = Y(iy) + Alpha*CABS1(X(j))*temp
                  ENDDO
               ENDIF
 
               IF ( .NOT.symb_zero ) Y(iy) = Y(iy) + SIGN(safe1,Y(iy))
 
               iy = iy + Incy
            ENDDO
         ELSE
            DO i = 1 , leny
               IF ( Beta==0.0 ) THEN
                  symb_zero = .TRUE.
                  Y(iy) = 0.0
               ELSEIF ( Y(iy)==0.0 ) THEN
                  symb_zero = .TRUE.
               ELSE
                  symb_zero = .FALSE.
                  Y(iy) = Beta*ABS(Y(iy))
               ENDIF
               IF ( Alpha/=0.0 ) THEN
                  DO j = MAX(i-Kl,1) , MIN(i+Ku,lenx)
                     temp = CABS1(Ab(ke-i+j,i))
                     symb_zero = symb_zero .AND.                        &
     &                           (X(j)==ZERO .OR. temp==ZERO)
 
                     Y(iy) = Y(iy) + Alpha*CABS1(X(j))*temp
                  ENDDO
               ENDIF
 
               IF ( .NOT.symb_zero ) Y(iy) = Y(iy) + SIGN(safe1,Y(iy))
 
               iy = iy + Incy
            ENDDO
         ENDIF
      ELSEIF ( Trans==ILATRANS('N') ) THEN
         DO i = 1 , leny
            IF ( Beta==0.0 ) THEN
               symb_zero = .TRUE.
               Y(iy) = 0.0
            ELSEIF ( Y(iy)==0.0 ) THEN
               symb_zero = .TRUE.
            ELSE
               symb_zero = .FALSE.
               Y(iy) = Beta*ABS(Y(iy))
            ENDIF
            IF ( Alpha/=0.0 ) THEN
               jx = kx
               DO j = MAX(i-Kl,1) , MIN(i+Ku,lenx)
                  temp = CABS1(Ab(kd+i-j,j))
                  symb_zero = symb_zero .AND.                           &
     &                        (X(jx)==ZERO .OR. temp==ZERO)
 
                  Y(iy) = Y(iy) + Alpha*CABS1(X(jx))*temp
                  jx = jx + Incx
               ENDDO
            ENDIF
 
            IF ( .NOT.symb_zero ) Y(iy) = Y(iy) + SIGN(safe1,Y(iy))
 
            iy = iy + Incy
         ENDDO
      ELSE
         DO i = 1 , leny
            IF ( Beta==0.0 ) THEN
               symb_zero = .TRUE.
               Y(iy) = 0.0
            ELSEIF ( Y(iy)==0.0 ) THEN
               symb_zero = .TRUE.
            ELSE
               symb_zero = .FALSE.
               Y(iy) = Beta*ABS(Y(iy))
            ENDIF
            IF ( Alpha/=0.0 ) THEN
               jx = kx
               DO j = MAX(i-Kl,1) , MIN(i+Ku,lenx)
                  temp = CABS1(Ab(ke-i+j,i))
                  symb_zero = symb_zero .AND.                           &
     &                        (X(jx)==ZERO .OR. temp==ZERO)
 
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
!     End of CLA_GBAMV
!
      END SUBROUTINE CLA_GBAMV
