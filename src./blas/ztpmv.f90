!*==ztpmv.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
!> \brief \b ZTPMV
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZTPMV(UPLO,TRANS,DIAG,N,AP,X,INCX)
!
!       .. Scalar Arguments ..
!       INTEGER INCX,N
!       CHARACTER DIAG,TRANS,UPLO
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
!> ZTPMV  performs one of the matrix-vector operations
!>
!>    x := A*x,   or   x := A**T*x,   or   x := A**H*x,
!>
!> where x is an n element vector and  A is an n by n unit, or non-unit,
!> upper or lower triangular matrix, supplied in packed form.
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
!>              TRANS = 'C' or 'c'   x := A**H*x.
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
!> \param[in] AP
!> \verbatim
!>          AP is COMPLEX*16 array, dimension at least
!>           ( ( n*( n + 1 ) )/2 ).
!>           Before entry with  UPLO = 'U' or 'u', the array AP must
!>           contain the upper triangular matrix packed sequentially,
!>           column by column, so that AP( 1 ) contains a( 1, 1 ),
!>           AP( 2 ) and AP( 3 ) contain a( 1, 2 ) and a( 2, 2 )
!>           respectively, and so on.
!>           Before entry with UPLO = 'L' or 'l', the array AP must
!>           contain the lower triangular matrix packed sequentially,
!>           column by column, so that AP( 1 ) contains a( 1, 1 ),
!>           AP( 2 ) and AP( 3 ) contain a( 2, 1 ) and a( 3, 1 )
!>           respectively, and so on.
!>           Note that when  DIAG = 'U' or 'u', the diagonal elements of
!>           A are not referenced, but are assumed to be unity.
!> \endverbatim
!>
!> \param[in,out] X
!> \verbatim
!>          X is COMPLEX*16 array, dimension at least
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
!> \ingroup complex16_blas_level2
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
      SUBROUTINE ZTPMV(Uplo,Trans,Diag,N,Ap,X,Incx)
      USE F77KINDS
      USE S_LSAME
      USE S_XERBLA
      IMPLICIT NONE
!*--ZTPMV149
!
! PARAMETER definitions rewritten by SPAG
!
      COMPLEX(cx16kind) , PARAMETER  ::  ZERO = (0.0D+0,0.0D+0)
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER , INTENT(IN) :: Uplo
      CHARACTER , INTENT(IN) :: Trans
      CHARACTER , INTENT(IN) :: Diag
      INTEGER , INTENT(IN) :: N
      COMPLEX(cx16kind) , INTENT(IN) , DIMENSION(*) :: Ap
      COMPLEX(cx16kind) , INTENT(INOUT) , DIMENSION(*) :: X
      INTEGER , INTENT(IN) :: Incx
!
! Local variable declarations rewritten by SPAG
!
      INTEGER :: i , info , ix , j , jx , k , kk , kx
      LOGICAL :: noconj , nounit
      COMPLEX(cx16kind) :: temp
!
! End of declarations rewritten by SPAG
!
!
! Dummy argument declarations rewritten by SPAG
!
!
! Local variable declarations rewritten by SPAG
!
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
      ELSEIF ( .NOT.LSAME(Trans,'N') .AND. .NOT.LSAME(Trans,'T') .AND.  &
     &         .NOT.LSAME(Trans,'C') ) THEN
         info = 2
      ELSEIF ( .NOT.LSAME(Diag,'U') .AND. .NOT.LSAME(Diag,'N') ) THEN
         info = 3
      ELSEIF ( N<0 ) THEN
         info = 4
      ELSEIF ( Incx==0 ) THEN
         info = 7
      ENDIF
      IF ( info/=0 ) THEN
         CALL XERBLA('ZTPMV ',info)
         RETURN
      ENDIF
!
!     Quick return if possible.
!
      IF ( N==0 ) RETURN
!
      noconj = LSAME(Trans,'T')
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
!     Start the operations. In this version the elements of AP are
!     accessed sequentially with one pass through AP.
!
      IF ( LSAME(Trans,'N') ) THEN
!
!        Form  x:= A*x.
!
         IF ( LSAME(Uplo,'U') ) THEN
            kk = 1
            IF ( Incx==1 ) THEN
               DO j = 1 , N
                  IF ( X(j)/=ZERO ) THEN
                     temp = X(j)
                     k = kk
                     DO i = 1 , j - 1
                        X(i) = X(i) + temp*Ap(k)
                        k = k + 1
                     ENDDO
                     IF ( nounit ) X(j) = X(j)*Ap(kk+j-1)
                  ENDIF
                  kk = kk + j
               ENDDO
            ELSE
               jx = kx
               DO j = 1 , N
                  IF ( X(jx)/=ZERO ) THEN
                     temp = X(jx)
                     ix = kx
                     DO k = kk , kk + j - 2
                        X(ix) = X(ix) + temp*Ap(k)
                        ix = ix + Incx
                     ENDDO
                     IF ( nounit ) X(jx) = X(jx)*Ap(kk+j-1)
                  ENDIF
                  jx = jx + Incx
                  kk = kk + j
               ENDDO
            ENDIF
         ELSE
            kk = (N*(N+1))/2
            IF ( Incx==1 ) THEN
               DO j = N , 1 , -1
                  IF ( X(j)/=ZERO ) THEN
                     temp = X(j)
                     k = kk
                     DO i = N , j + 1 , -1
                        X(i) = X(i) + temp*Ap(k)
                        k = k - 1
                     ENDDO
                     IF ( nounit ) X(j) = X(j)*Ap(kk-N+j)
                  ENDIF
                  kk = kk - (N-j+1)
               ENDDO
            ELSE
               kx = kx + (N-1)*Incx
               jx = kx
               DO j = N , 1 , -1
                  IF ( X(jx)/=ZERO ) THEN
                     temp = X(jx)
                     ix = kx
                     DO k = kk , kk - (N-(j+1)) , -1
                        X(ix) = X(ix) + temp*Ap(k)
                        ix = ix - Incx
                     ENDDO
                     IF ( nounit ) X(jx) = X(jx)*Ap(kk-N+j)
                  ENDIF
                  jx = jx - Incx
                  kk = kk - (N-j+1)
               ENDDO
            ENDIF
         ENDIF
!
!        Form  x := A**T*x  or  x := A**H*x.
!
      ELSEIF ( LSAME(Uplo,'U') ) THEN
         kk = (N*(N+1))/2
         IF ( Incx==1 ) THEN
            DO j = N , 1 , -1
               temp = X(j)
               k = kk - 1
               IF ( noconj ) THEN
                  IF ( nounit ) temp = temp*Ap(kk)
                  DO i = j - 1 , 1 , -1
                     temp = temp + Ap(k)*X(i)
                     k = k - 1
                  ENDDO
               ELSE
                  IF ( nounit ) temp = temp*DCONJG(Ap(kk))
                  DO i = j - 1 , 1 , -1
                     temp = temp + DCONJG(Ap(k))*X(i)
                     k = k - 1
                  ENDDO
               ENDIF
               X(j) = temp
               kk = kk - j
            ENDDO
         ELSE
            jx = kx + (N-1)*Incx
            DO j = N , 1 , -1
               temp = X(jx)
               ix = jx
               IF ( noconj ) THEN
                  IF ( nounit ) temp = temp*Ap(kk)
                  DO k = kk - 1 , kk - j + 1 , -1
                     ix = ix - Incx
                     temp = temp + Ap(k)*X(ix)
                  ENDDO
               ELSE
                  IF ( nounit ) temp = temp*DCONJG(Ap(kk))
                  DO k = kk - 1 , kk - j + 1 , -1
                     ix = ix - Incx
                     temp = temp + DCONJG(Ap(k))*X(ix)
                  ENDDO
               ENDIF
               X(jx) = temp
               jx = jx - Incx
               kk = kk - j
            ENDDO
         ENDIF
      ELSE
         kk = 1
         IF ( Incx==1 ) THEN
            DO j = 1 , N
               temp = X(j)
               k = kk + 1
               IF ( noconj ) THEN
                  IF ( nounit ) temp = temp*Ap(kk)
                  DO i = j + 1 , N
                     temp = temp + Ap(k)*X(i)
                     k = k + 1
                  ENDDO
               ELSE
                  IF ( nounit ) temp = temp*DCONJG(Ap(kk))
                  DO i = j + 1 , N
                     temp = temp + DCONJG(Ap(k))*X(i)
                     k = k + 1
                  ENDDO
               ENDIF
               X(j) = temp
               kk = kk + (N-j+1)
            ENDDO
         ELSE
            jx = kx
            DO j = 1 , N
               temp = X(jx)
               ix = jx
               IF ( noconj ) THEN
                  IF ( nounit ) temp = temp*Ap(kk)
                  DO k = kk + 1 , kk + N - j
                     ix = ix + Incx
                     temp = temp + Ap(k)*X(ix)
                  ENDDO
               ELSE
                  IF ( nounit ) temp = temp*DCONJG(Ap(kk))
                  DO k = kk + 1 , kk + N - j
                     ix = ix + Incx
                     temp = temp + DCONJG(Ap(k))*X(ix)
                  ENDDO
               ENDIF
               X(jx) = temp
               jx = jx + Incx
               kk = kk + (N-j+1)
            ENDDO
         ENDIF
      ENDIF
!
!
!     End of ZTPMV .
!
      END SUBROUTINE ZTPMV
