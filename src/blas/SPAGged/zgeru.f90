!*==zgeru.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b ZGERU
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZGERU(M,N,ALPHA,X,INCX,Y,INCY,A,LDA)
!
!       .. Scalar Arguments ..
!       COMPLEX*16 ALPHA
!       INTEGER INCX,INCY,LDA,M,N
!       ..
!       .. Array Arguments ..
!       COMPLEX*16 A(LDA,*),X(*),Y(*)
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZGERU  performs the rank 1 operation
!>
!>    A := alpha*x*y**T + A,
!>
!> where alpha is a scalar, x is an m element vector, y is an n element
!> vector and A is an m by n matrix.
!> \endverbatim
!
!  Arguments:
!  ==========
!
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
!>          ALPHA is COMPLEX*16
!>           On entry, ALPHA specifies the scalar alpha.
!> \endverbatim
!>
!> \param[in] X
!> \verbatim
!>          X is COMPLEX*16 array, dimension at least
!>           ( 1 + ( m - 1 )*abs( INCX ) ).
!>           Before entry, the incremented array X must contain the m
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
!>          Y is COMPLEX*16 array, dimension at least
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
!>          A is COMPLEX*16 array, dimension ( LDA, N )
!>           Before entry, the leading m by n part of the array A must
!>           contain the matrix of coefficients. On exit, A is
!>           overwritten by the updated matrix.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>           On entry, LDA specifies the first dimension of A as declared
!>           in the calling (sub) program. LDA must be at least
!>           max( 1, m ).
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
      SUBROUTINE ZGERU(M,N,Alpha,X,Incx,Y,Incy,A,Lda)
      USE F77KINDS                        
      USE S_XERBLA
      IMPLICIT NONE
!*--ZGERU136
!
! PARAMETER definitions rewritten by SPAG
!
      COMPLEX(CX16KIND) , PARAMETER  ::  ZERO = (0.0D+0,0.0D+0)
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(IN) :: Alpha
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(*) :: X
      INTEGER , INTENT(IN) :: Incx
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(*) :: Y
      INTEGER , INTENT(IN) :: Incy
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
!
! Local variable declarations rewritten by SPAG
!
      INTEGER :: i , info , ix , j , jy , kx
      COMPLEX(CX16KIND) :: temp
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
!     .. External Subroutines ..
!     ..
!     .. Intrinsic Functions ..
!     ..
!
!     Test the input parameters.
!
      info = 0
      IF ( M<0 ) THEN
         info = 1
      ELSEIF ( N<0 ) THEN
         info = 2
      ELSEIF ( Incx==0 ) THEN
         info = 5
      ELSEIF ( Incy==0 ) THEN
         info = 7
      ELSEIF ( Lda<MAX(1,M) ) THEN
         info = 9
      ENDIF
      IF ( info/=0 ) THEN
         CALL XERBLA('ZGERU ',info)
         RETURN
      ENDIF
!
!     Quick return if possible.
!
      IF ( (M==0) .OR. (N==0) .OR. (Alpha==ZERO) ) RETURN
!
!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through A.
!
      IF ( Incy>0 ) THEN
         jy = 1
      ELSE
         jy = 1 - (N-1)*Incy
      ENDIF
      IF ( Incx==1 ) THEN
         DO j = 1 , N
            IF ( Y(jy)/=ZERO ) THEN
               temp = Alpha*Y(jy)
               DO i = 1 , M
                  A(i,j) = A(i,j) + X(i)*temp
               ENDDO
            ENDIF
            jy = jy + Incy
         ENDDO
      ELSE
         IF ( Incx>0 ) THEN
            kx = 1
         ELSE
            kx = 1 - (M-1)*Incx
         ENDIF
         DO j = 1 , N
            IF ( Y(jy)/=ZERO ) THEN
               temp = Alpha*Y(jy)
               ix = kx
               DO i = 1 , M
                  A(i,j) = A(i,j) + X(ix)*temp
                  ix = ix + Incx
               ENDDO
            ENDIF
            jy = jy + Incy
         ENDDO
      ENDIF
!
!
!     End of ZGERU .
!
      END SUBROUTINE ZGERU
