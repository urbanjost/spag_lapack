!*==cgemv.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b CGEMV
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE CGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
!
!       .. Scalar Arguments ..
!       COMPLEX ALPHA,BETA
!       INTEGER INCX,INCY,LDA,M,N
!       CHARACTER TRANS
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
!> CGEMV performs one of the matrix-vector operations
!>
!>    y := alpha*A*x + beta*y,   or   y := alpha*A**T*x + beta*y,   or
!>
!>    y := alpha*A**H*x + beta*y,
!>
!> where alpha and beta are scalars, x and y are vectors and A is an
!> m by n matrix.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] TRANS
!> \verbatim
!>          TRANS is CHARACTER*1
!>           On entry, TRANS specifies the operation to be performed as
!>           follows:
!>
!>              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.
!>
!>              TRANS = 'T' or 't'   y := alpha*A**T*x + beta*y.
!>
!>              TRANS = 'C' or 'c'   y := alpha*A**H*x + beta*y.
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>           On entry, M specifies the number of rows of the matrix A.
!>           M must be at least zero.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>           On entry, N specifies the number of columns of the matrix A.
!>           N must be at least zero.
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
!>           Before entry, the leading m by n part of the array A must
!>           contain the matrix of coefficients.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>           On entry, LDA specifies the first dimension of A as declared
!>           in the calling (sub) program. LDA must be at least
!>           max( 1, m ).
!> \endverbatim
!>
!> \param[in] X
!> \verbatim
!>          X is COMPLEX array, dimension at least
!>           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
!>           and at least
!>           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
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
!>           On entry, BETA specifies the scalar beta. When BETA is
!>           supplied as zero then Y need not be set on input.
!> \endverbatim
!>
!> \param[in,out] Y
!> \verbatim
!>          Y is COMPLEX array, dimension at least
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
      SUBROUTINE CGEMV(Trans,M,N,Alpha,A,Lda,X,Incx,Beta,Y,Incy)
      IMPLICIT NONE
!*--CGEMV162
!
!  -- Reference BLAS level2 routine (version 3.7.0) --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      COMPLEX Alpha , Beta
      INTEGER Incx , Incy , Lda , M , N
      CHARACTER Trans
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
      COMPLEX temp
      INTEGER i , info , ix , iy , j , jx , jy , kx , ky , lenx , leny
      LOGICAL noconj
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC CONJG , MAX
!     ..
!
!     Test the input parameters.
!
      info = 0
      IF ( .NOT.LSAME(Trans,'N') .AND. .NOT.LSAME(Trans,'T') .AND.      &
     &     .NOT.LSAME(Trans,'C') ) THEN
         info = 1
      ELSEIF ( M<0 ) THEN
         info = 2
      ELSEIF ( N<0 ) THEN
         info = 3
      ELSEIF ( Lda<MAX(1,M) ) THEN
         info = 6
      ELSEIF ( Incx==0 ) THEN
         info = 8
      ELSEIF ( Incy==0 ) THEN
         info = 11
      ENDIF
      IF ( info/=0 ) THEN
         CALL XERBLA('CGEMV ',info)
         RETURN
      ENDIF
!
!     Quick return if possible.
!
      IF ( (M==0) .OR. (N==0) .OR. ((Alpha==ZERO) .AND. (Beta==ONE)) )  &
     &     RETURN
!
      noconj = LSAME(Trans,'T')
!
!     Set  LENX  and  LENY, the lengths of the vectors x and y, and set
!     up the start points in  X  and  Y.
!
      IF ( LSAME(Trans,'N') ) THEN
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
!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through A.
!
!     First form  y := beta*y.
!
      IF ( Beta/=ONE ) THEN
         IF ( Incy/=1 ) THEN
            iy = ky
            IF ( Beta==ZERO ) THEN
               DO i = 1 , leny
                  Y(iy) = ZERO
                  iy = iy + Incy
               ENDDO
            ELSE
               DO i = 1 , leny
                  Y(iy) = Beta*Y(iy)
                  iy = iy + Incy
               ENDDO
            ENDIF
         ELSEIF ( Beta==ZERO ) THEN
            DO i = 1 , leny
               Y(i) = ZERO
            ENDDO
         ELSE
            DO i = 1 , leny
               Y(i) = Beta*Y(i)
            ENDDO
         ENDIF
      ENDIF
      IF ( Alpha==ZERO ) RETURN
      IF ( LSAME(Trans,'N') ) THEN
!
!        Form  y := alpha*A*x + y.
!
         jx = kx
         IF ( Incy==1 ) THEN
            DO j = 1 , N
               temp = Alpha*X(jx)
               DO i = 1 , M
                  Y(i) = Y(i) + temp*A(i,j)
               ENDDO
               jx = jx + Incx
            ENDDO
         ELSE
            DO j = 1 , N
               temp = Alpha*X(jx)
               iy = ky
               DO i = 1 , M
                  Y(iy) = Y(iy) + temp*A(i,j)
                  iy = iy + Incy
               ENDDO
               jx = jx + Incx
            ENDDO
         ENDIF
      ELSE
!
!        Form  y := alpha*A**T*x + y  or  y := alpha*A**H*x + y.
!
         jy = ky
         IF ( Incx==1 ) THEN
            DO j = 1 , N
               temp = ZERO
               IF ( noconj ) THEN
                  DO i = 1 , M
                     temp = temp + A(i,j)*X(i)
                  ENDDO
               ELSE
                  DO i = 1 , M
                     temp = temp + CONJG(A(i,j))*X(i)
                  ENDDO
               ENDIF
               Y(jy) = Y(jy) + Alpha*temp
               jy = jy + Incy
            ENDDO
         ELSE
            DO j = 1 , N
               temp = ZERO
               ix = kx
               IF ( noconj ) THEN
                  DO i = 1 , M
                     temp = temp + A(i,j)*X(ix)
                     ix = ix + Incx
                  ENDDO
               ELSE
                  DO i = 1 , M
                     temp = temp + CONJG(A(i,j))*X(ix)
                     ix = ix + Incx
                  ENDDO
               ENDIF
               Y(jy) = Y(jy) + Alpha*temp
               jy = jy + Incy
            ENDDO
         ENDIF
      ENDIF
!
!
!     End of CGEMV .
!
      END SUBROUTINE CGEMV
