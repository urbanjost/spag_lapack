!*==strmv.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b STRMV
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE STRMV(UPLO,TRANS,DIAG,N,A,LDA,X,INCX)
!
!       .. Scalar Arguments ..
!       INTEGER INCX,LDA,N
!       CHARACTER DIAG,TRANS,UPLO
!       ..
!       .. Array Arguments ..
!       REAL A(LDA,*),X(*)
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> STRMV  performs one of the matrix-vector operations
!>
!>    x := A*x,   or   x := A**T*x,
!>
!> where x is an n element vector and  A is an n by n unit, or non-unit,
!> upper or lower triangular matrix.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>           On entry, UPLO specifies whether the matrix is an upper or
!>           lower triangular matrix as follows:
!>
!>              UPLO = 'U' or 'u'   A is an upper triangular matrix.
!>
!>              UPLO = 'L' or 'l'   A is a lower triangular matrix.
!> \endverbatim
!>
!> \param[in] TRANS
!> \verbatim
!>          TRANS is CHARACTER*1
!>           On entry, TRANS specifies the operation to be performed as
!>           follows:
!>
!>              TRANS = 'N' or 'n'   x := A*x.
!>
!>              TRANS = 'T' or 't'   x := A**T*x.
!>
!>              TRANS = 'C' or 'c'   x := A**T*x.
!> \endverbatim
!>
!> \param[in] DIAG
!> \verbatim
!>          DIAG is CHARACTER*1
!>           On entry, DIAG specifies whether or not A is unit
!>           triangular as follows:
!>
!>              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
!>
!>              DIAG = 'N' or 'n'   A is not assumed to be unit
!>                                  triangular.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>           On entry, N specifies the order of the matrix A.
!>           N must be at least zero.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is REAL array, dimension ( LDA, N )
!>           Before entry with  UPLO = 'U' or 'u', the leading n by n
!>           upper triangular part of the array A must contain the upper
!>           triangular matrix and the strictly lower triangular part of
!>           A is not referenced.
!>           Before entry with UPLO = 'L' or 'l', the leading n by n
!>           lower triangular part of the array A must contain the lower
!>           triangular matrix and the strictly upper triangular part of
!>           A is not referenced.
!>           Note that when  DIAG = 'U' or 'u', the diagonal elements of
!>           A are not referenced either, but are assumed to be unity.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>           On entry, LDA specifies the first dimension of A as declared
!>           in the calling (sub) program. LDA must be at least
!>           max( 1, n ).
!> \endverbatim
!>
!> \param[in,out] X
!> \verbatim
!>          X is REAL array, dimension at least
!>           ( 1 + ( n - 1 )*abs( INCX ) ).
!>           Before entry, the incremented array X must contain the n
!>           element vector x. On exit, X is overwritten with the
!>           transformed vector x.
!> \endverbatim
!>
!> \param[in] INCX
!> \verbatim
!>          INCX is INTEGER
!>           On entry, INCX specifies the increment for the elements of
!>           X. INCX must not be zero.
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
!> \ingroup single_blas_level2
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
      SUBROUTINE STRMV(Uplo,Trans,Diag,N,A,Lda,X,Incx)
      IMPLICIT NONE
!*--STRMV151
!
!  -- Reference BLAS level2 routine (version 3.7.0) --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Incx , Lda , N
      CHARACTER Diag , Trans , Uplo
!     ..
!     .. Array Arguments ..
      REAL A(Lda,*) , X(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ZERO
      PARAMETER (ZERO=0.0E+0)
!     ..
!     .. Local Scalars ..
      REAL temp
      INTEGER i , info , ix , j , jx , kx
      LOGICAL nounit
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX
!     ..
!
!     Test the input parameters.
!
      info = 0
      IF ( .NOT.LSAME(Uplo,'U') .AND. .NOT.LSAME(Uplo,'L') ) THEN
         info = 1
      ELSEIF ( .NOT.LSAME(Trans,'N') .AND. .NOT.LSAME(Trans,'T') .AND.  &
     &         .NOT.LSAME(Trans,'C') ) THEN
         info = 2
      ELSEIF ( .NOT.LSAME(Diag,'U') .AND. .NOT.LSAME(Diag,'N') ) THEN
         info = 3
      ELSEIF ( N<0 ) THEN
         info = 4
      ELSEIF ( Lda<MAX(1,N) ) THEN
         info = 6
      ELSEIF ( Incx==0 ) THEN
         info = 8
      ENDIF
      IF ( info/=0 ) THEN
         CALL XERBLA('STRMV ',info)
         RETURN
      ENDIF
!
!     Quick return if possible.
!
      IF ( N==0 ) RETURN
!
      nounit = LSAME(Diag,'N')
!
!     Set up the start point in X if the increment is not unity. This
!     will be  ( N - 1 )*INCX  too small for descending loops.
!
      IF ( Incx<=0 ) THEN
         kx = 1 - (N-1)*Incx
      ELSEIF ( Incx/=1 ) THEN
         kx = 1
      ENDIF
!
!     Start the operations. In this version the elements of A are
!     accessed sequentially with one pass through A.
!
      IF ( LSAME(Trans,'N') ) THEN
!
!        Form  x := A*x.
!
         IF ( LSAME(Uplo,'U') ) THEN
            IF ( Incx==1 ) THEN
               DO j = 1 , N
                  IF ( X(j)/=ZERO ) THEN
                     temp = X(j)
                     DO i = 1 , j - 1
                        X(i) = X(i) + temp*A(i,j)
                     ENDDO
                     IF ( nounit ) X(j) = X(j)*A(j,j)
                  ENDIF
               ENDDO
            ELSE
               jx = kx
               DO j = 1 , N
                  IF ( X(jx)/=ZERO ) THEN
                     temp = X(jx)
                     ix = kx
                     DO i = 1 , j - 1
                        X(ix) = X(ix) + temp*A(i,j)
                        ix = ix + Incx
                     ENDDO
                     IF ( nounit ) X(jx) = X(jx)*A(j,j)
                  ENDIF
                  jx = jx + Incx
               ENDDO
            ENDIF
         ELSEIF ( Incx==1 ) THEN
            DO j = N , 1 , -1
               IF ( X(j)/=ZERO ) THEN
                  temp = X(j)
                  DO i = N , j + 1 , -1
                     X(i) = X(i) + temp*A(i,j)
                  ENDDO
                  IF ( nounit ) X(j) = X(j)*A(j,j)
               ENDIF
            ENDDO
         ELSE
            kx = kx + (N-1)*Incx
            jx = kx
            DO j = N , 1 , -1
               IF ( X(jx)/=ZERO ) THEN
                  temp = X(jx)
                  ix = kx
                  DO i = N , j + 1 , -1
                     X(ix) = X(ix) + temp*A(i,j)
                     ix = ix - Incx
                  ENDDO
                  IF ( nounit ) X(jx) = X(jx)*A(j,j)
               ENDIF
               jx = jx - Incx
            ENDDO
         ENDIF
!
!        Form  x := A**T*x.
!
      ELSEIF ( LSAME(Uplo,'U') ) THEN
         IF ( Incx==1 ) THEN
            DO j = N , 1 , -1
               temp = X(j)
               IF ( nounit ) temp = temp*A(j,j)
               DO i = j - 1 , 1 , -1
                  temp = temp + A(i,j)*X(i)
               ENDDO
               X(j) = temp
            ENDDO
         ELSE
            jx = kx + (N-1)*Incx
            DO j = N , 1 , -1
               temp = X(jx)
               ix = jx
               IF ( nounit ) temp = temp*A(j,j)
               DO i = j - 1 , 1 , -1
                  ix = ix - Incx
                  temp = temp + A(i,j)*X(ix)
               ENDDO
               X(jx) = temp
               jx = jx - Incx
            ENDDO
         ENDIF
      ELSEIF ( Incx==1 ) THEN
         DO j = 1 , N
            temp = X(j)
            IF ( nounit ) temp = temp*A(j,j)
            DO i = j + 1 , N
               temp = temp + A(i,j)*X(i)
            ENDDO
            X(j) = temp
         ENDDO
      ELSE
         jx = kx
         DO j = 1 , N
            temp = X(jx)
            ix = jx
            IF ( nounit ) temp = temp*A(j,j)
            DO i = j + 1 , N
               ix = ix + Incx
               temp = temp + A(i,j)*X(ix)
            ENDDO
            X(jx) = temp
            jx = jx + Incx
         ENDDO
      ENDIF
!
!
!     End of STRMV .
!
      END SUBROUTINE STRMV
