!*==dgbmv.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b DGBMV
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DGBMV(TRANS,M,N,KL,KU,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
!
!       .. Scalar Arguments ..
!       DOUBLE PRECISION ALPHA,BETA
!       INTEGER INCX,INCY,KL,KU,LDA,M,N
!       CHARACTER TRANS
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION A(LDA,*),X(*),Y(*)
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DGBMV  performs one of the matrix-vector operations
!>
!>    y := alpha*A*x + beta*y,   or   y := alpha*A**T*x + beta*y,
!>
!> where alpha and beta are scalars, x and y are vectors and A is an
!> m by n band matrix, with kl sub-diagonals and ku super-diagonals.
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
!>              TRANS = 'C' or 'c'   y := alpha*A**T*x + beta*y.
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
!> \param[in] KL
!> \verbatim
!>          KL is INTEGER
!>           On entry, KL specifies the number of sub-diagonals of the
!>           matrix A. KL must satisfy  0 .le. KL.
!> \endverbatim
!>
!> \param[in] KU
!> \verbatim
!>          KU is INTEGER
!>           On entry, KU specifies the number of super-diagonals of the
!>           matrix A. KU must satisfy  0 .le. KU.
!> \endverbatim
!>
!> \param[in] ALPHA
!> \verbatim
!>          ALPHA is DOUBLE PRECISION.
!>           On entry, ALPHA specifies the scalar alpha.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension ( LDA, N )
!>           Before entry, the leading ( kl + ku + 1 ) by n part of the
!>           array A must contain the matrix of coefficients, supplied
!>           column by column, with the leading diagonal of the matrix in
!>           row ( ku + 1 ) of the array, the first super-diagonal
!>           starting at position 2 in row ku, the first sub-diagonal
!>           starting at position 1 in row ( ku + 2 ), and so on.
!>           Elements in the array A that do not correspond to elements
!>           in the band matrix (such as the top left ku by ku triangle)
!>           are not referenced.
!>           The following program segment will transfer a band matrix
!>           from conventional full matrix storage to band storage:
!>
!>                 DO 20, J = 1, N
!>                    K = KU + 1 - J
!>                    DO 10, I = MAX( 1, J - KU ), MIN( M, J + KL )
!>                       A( K + I, J ) = matrix( I, J )
!>              10    CONTINUE
!>              20 CONTINUE
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>           On entry, LDA specifies the first dimension of A as declared
!>           in the calling (sub) program. LDA must be at least
!>           ( kl + ku + 1 ).
!> \endverbatim
!>
!> \param[in] X
!> \verbatim
!>          X is DOUBLE PRECISION array, dimension at least
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
!>          BETA is DOUBLE PRECISION.
!>           On entry, BETA specifies the scalar beta. When BETA is
!>           supplied as zero then Y need not be set on input.
!> \endverbatim
!>
!> \param[in,out] Y
!> \verbatim
!>          Y is DOUBLE PRECISION array, dimension at least
!>           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
!>           and at least
!>           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
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
!> \ingroup double_blas_level2
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
      SUBROUTINE DGBMV(Trans,M,N,Kl,Ku,Alpha,A,Lda,X,Incx,Beta,Y,Incy)
      IMPLICIT NONE
!*--DGBMV189
!
!  -- Reference BLAS level2 routine (version 3.7.0) --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      DOUBLE PRECISION Alpha , Beta
      INTEGER Incx , Incy , Kl , Ku , Lda , M , N
      CHARACTER Trans
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION A(Lda,*) , X(*) , Y(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ONE , ZERO
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
!     ..
!     .. Local Scalars ..
      DOUBLE PRECISION temp
      INTEGER i , info , ix , iy , j , jx , jy , k , kup1 , kx , ky ,   &
     &        lenx , leny
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX , MIN
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
      ELSEIF ( Kl<0 ) THEN
         info = 4
      ELSEIF ( Ku<0 ) THEN
         info = 5
      ELSEIF ( Lda<(Kl+Ku+1) ) THEN
         info = 8
      ELSEIF ( Incx==0 ) THEN
         info = 10
      ELSEIF ( Incy==0 ) THEN
         info = 13
      ENDIF
      IF ( info/=0 ) THEN
         CALL XERBLA('DGBMV ',info)
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
!     accessed sequentially with one pass through the band part of A.
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
      kup1 = Ku + 1
      IF ( LSAME(Trans,'N') ) THEN
!
!        Form  y := alpha*A*x + y.
!
         jx = kx
         IF ( Incy==1 ) THEN
            DO j = 1 , N
               temp = Alpha*X(jx)
               k = kup1 - j
               DO i = MAX(1,j-Ku) , MIN(M,j+Kl)
                  Y(i) = Y(i) + temp*A(k+i,j)
               ENDDO
               jx = jx + Incx
            ENDDO
         ELSE
            DO j = 1 , N
               temp = Alpha*X(jx)
               iy = ky
               k = kup1 - j
               DO i = MAX(1,j-Ku) , MIN(M,j+Kl)
                  Y(iy) = Y(iy) + temp*A(k+i,j)
                  iy = iy + Incy
               ENDDO
               jx = jx + Incx
               IF ( j>Ku ) ky = ky + Incy
            ENDDO
         ENDIF
      ELSE
!
!        Form  y := alpha*A**T*x + y.
!
         jy = ky
         IF ( Incx==1 ) THEN
            DO j = 1 , N
               temp = ZERO
               k = kup1 - j
               DO i = MAX(1,j-Ku) , MIN(M,j+Kl)
                  temp = temp + A(k+i,j)*X(i)
               ENDDO
               Y(jy) = Y(jy) + Alpha*temp
               jy = jy + Incy
            ENDDO
         ELSE
            DO j = 1 , N
               temp = ZERO
               ix = kx
               k = kup1 - j
               DO i = MAX(1,j-Ku) , MIN(M,j+Kl)
                  temp = temp + A(k+i,j)*X(ix)
                  ix = ix + Incx
               ENDDO
               Y(jy) = Y(jy) + Alpha*temp
               jy = jy + Incy
               IF ( j>Ku ) kx = kx + Incx
            ENDDO
         ENDIF
      ENDIF
!
!
!     End of DGBMV .
!
      END SUBROUTINE DGBMV
