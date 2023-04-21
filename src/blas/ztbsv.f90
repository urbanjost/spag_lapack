!*==ztbsv.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b ZTBSV
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZTBSV(UPLO,TRANS,DIAG,N,K,A,LDA,X,INCX)
!
!       .. Scalar Arguments ..
!       INTEGER INCX,K,LDA,N
!       CHARACTER DIAG,TRANS,UPLO
!       ..
!       .. Array Arguments ..
!       COMPLEX*16 A(LDA,*),X(*)
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZTBSV  solves one of the systems of equations
!>
!>    A*x = b,   or   A**T*x = b,   or   A**H*x = b,
!>
!> where b and x are n element vectors and A is an n by n unit, or
!> non-unit, upper or lower triangular band matrix, with ( k + 1 )
!> diagonals.
!>
!> No test for singularity or near-singularity is included in this
!> routine. Such tests must be performed before calling this routine.
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
!>           On entry, TRANS specifies the equations to be solved as
!>           follows:
!>
!>              TRANS = 'N' or 'n'   A*x = b.
!>
!>              TRANS = 'T' or 't'   A**T*x = b.
!>
!>              TRANS = 'C' or 'c'   A**H*x = b.
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
!> \param[in] K
!> \verbatim
!>          K is INTEGER
!>           On entry with UPLO = 'U' or 'u', K specifies the number of
!>           super-diagonals of the matrix A.
!>           On entry with UPLO = 'L' or 'l', K specifies the number of
!>           sub-diagonals of the matrix A.
!>           K must satisfy  0 .le. K.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension ( LDA, N )
!>           Before entry with UPLO = 'U' or 'u', the leading ( k + 1 )
!>           by n part of the array A must contain the upper triangular
!>           band part of the matrix of coefficients, supplied column by
!>           column, with the leading diagonal of the matrix in row
!>           ( k + 1 ) of the array, the first super-diagonal starting at
!>           position 2 in row k, and so on. The top left k by k triangle
!>           of the array A is not referenced.
!>           The following program segment will transfer an upper
!>           triangular band matrix from conventional full matrix storage
!>           to band storage:
!>
!>                 DO 20, J = 1, N
!>                    M = K + 1 - J
!>                    DO 10, I = MAX( 1, J - K ), J
!>                       A( M + I, J ) = matrix( I, J )
!>              10    CONTINUE
!>              20 CONTINUE
!>
!>           Before entry with UPLO = 'L' or 'l', the leading ( k + 1 )
!>           by n part of the array A must contain the lower triangular
!>           band part of the matrix of coefficients, supplied column by
!>           column, with the leading diagonal of the matrix in row 1 of
!>           the array, the first sub-diagonal starting at position 1 in
!>           row 2, and so on. The bottom right k by k triangle of the
!>           array A is not referenced.
!>           The following program segment will transfer a lower
!>           triangular band matrix from conventional full matrix storage
!>           to band storage:
!>
!>                 DO 20, J = 1, N
!>                    M = 1 - J
!>                    DO 10, I = J, MIN( N, J + K )
!>                       A( M + I, J ) = matrix( I, J )
!>              10    CONTINUE
!>              20 CONTINUE
!>
!>           Note that when DIAG = 'U' or 'u' the elements of the array A
!>           corresponding to the diagonal elements of the matrix are not
!>           referenced, but are assumed to be unity.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>           On entry, LDA specifies the first dimension of A as declared
!>           in the calling (sub) program. LDA must be at least
!>           ( k + 1 ).
!> \endverbatim
!>
!> \param[in,out] X
!> \verbatim
!>          X is COMPLEX*16 array, dimension at least
!>           ( 1 + ( n - 1 )*abs( INCX ) ).
!>           Before entry, the incremented array X must contain the n
!>           element right-hand side vector b. On exit, X is overwritten
!>           with the solution vector x.
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
!>
!>  -- Written on 22-October-1986.
!>     Jack Dongarra, Argonne National Lab.
!>     Jeremy Du Croz, Nag Central Office.
!>     Sven Hammarling, Nag Central Office.
!>     Richard Hanson, Sandia National Labs.
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE ZTBSV(Uplo,Trans,Diag,N,K,A,Lda,X,Incx)
      IMPLICIT NONE
!*--ZTBSV193
!
!  -- Reference BLAS level2 routine (version 3.7.0) --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Incx , K , Lda , N
      CHARACTER Diag , Trans , Uplo
!     ..
!     .. Array Arguments ..
      COMPLEX*16 A(Lda,*) , X(*)
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
      INTEGER i , info , ix , j , jx , kplus1 , kx , l
      LOGICAL noconj , nounit
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC DCONJG , MAX , MIN
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
      ELSEIF ( K<0 ) THEN
         info = 5
      ELSEIF ( Lda<(K+1) ) THEN
         info = 7
      ELSEIF ( Incx==0 ) THEN
         info = 9
      ENDIF
      IF ( info/=0 ) THEN
         CALL XERBLA('ZTBSV ',info)
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
!     Start the operations. In this version the elements of A are
!     accessed by sequentially with one pass through A.
!
      IF ( LSAME(Trans,'N') ) THEN
!
!        Form  x := inv( A )*x.
!
         IF ( LSAME(Uplo,'U') ) THEN
            kplus1 = K + 1
            IF ( Incx==1 ) THEN
               DO j = N , 1 , -1
                  IF ( X(j)/=ZERO ) THEN
                     l = kplus1 - j
                     IF ( nounit ) X(j) = X(j)/A(kplus1,j)
                     temp = X(j)
                     DO i = j - 1 , MAX(1,j-K) , -1
                        X(i) = X(i) - temp*A(l+i,j)
                     ENDDO
                  ENDIF
               ENDDO
            ELSE
               kx = kx + (N-1)*Incx
               jx = kx
               DO j = N , 1 , -1
                  kx = kx - Incx
                  IF ( X(jx)/=ZERO ) THEN
                     ix = kx
                     l = kplus1 - j
                     IF ( nounit ) X(jx) = X(jx)/A(kplus1,j)
                     temp = X(jx)
                     DO i = j - 1 , MAX(1,j-K) , -1
                        X(ix) = X(ix) - temp*A(l+i,j)
                        ix = ix - Incx
                     ENDDO
                  ENDIF
                  jx = jx - Incx
               ENDDO
            ENDIF
         ELSEIF ( Incx==1 ) THEN
            DO j = 1 , N
               IF ( X(j)/=ZERO ) THEN
                  l = 1 - j
                  IF ( nounit ) X(j) = X(j)/A(1,j)
                  temp = X(j)
                  DO i = j + 1 , MIN(N,j+K)
                     X(i) = X(i) - temp*A(l+i,j)
                  ENDDO
               ENDIF
            ENDDO
         ELSE
            jx = kx
            DO j = 1 , N
               kx = kx + Incx
               IF ( X(jx)/=ZERO ) THEN
                  ix = kx
                  l = 1 - j
                  IF ( nounit ) X(jx) = X(jx)/A(1,j)
                  temp = X(jx)
                  DO i = j + 1 , MIN(N,j+K)
                     X(ix) = X(ix) - temp*A(l+i,j)
                     ix = ix + Incx
                  ENDDO
               ENDIF
               jx = jx + Incx
            ENDDO
         ENDIF
!
!        Form  x := inv( A**T )*x  or  x := inv( A**H )*x.
!
      ELSEIF ( LSAME(Uplo,'U') ) THEN
         kplus1 = K + 1
         IF ( Incx==1 ) THEN
            DO j = 1 , N
               temp = X(j)
               l = kplus1 - j
               IF ( noconj ) THEN
                  DO i = MAX(1,j-K) , j - 1
                     temp = temp - A(l+i,j)*X(i)
                  ENDDO
                  IF ( nounit ) temp = temp/A(kplus1,j)
               ELSE
                  DO i = MAX(1,j-K) , j - 1
                     temp = temp - DCONJG(A(l+i,j))*X(i)
                  ENDDO
                  IF ( nounit ) temp = temp/DCONJG(A(kplus1,j))
               ENDIF
               X(j) = temp
            ENDDO
         ELSE
            jx = kx
            DO j = 1 , N
               temp = X(jx)
               ix = kx
               l = kplus1 - j
               IF ( noconj ) THEN
                  DO i = MAX(1,j-K) , j - 1
                     temp = temp - A(l+i,j)*X(ix)
                     ix = ix + Incx
                  ENDDO
                  IF ( nounit ) temp = temp/A(kplus1,j)
               ELSE
                  DO i = MAX(1,j-K) , j - 1
                     temp = temp - DCONJG(A(l+i,j))*X(ix)
                     ix = ix + Incx
                  ENDDO
                  IF ( nounit ) temp = temp/DCONJG(A(kplus1,j))
               ENDIF
               X(jx) = temp
               jx = jx + Incx
               IF ( j>K ) kx = kx + Incx
            ENDDO
         ENDIF
      ELSEIF ( Incx==1 ) THEN
         DO j = N , 1 , -1
            temp = X(j)
            l = 1 - j
            IF ( noconj ) THEN
               DO i = MIN(N,j+K) , j + 1 , -1
                  temp = temp - A(l+i,j)*X(i)
               ENDDO
               IF ( nounit ) temp = temp/A(1,j)
            ELSE
               DO i = MIN(N,j+K) , j + 1 , -1
                  temp = temp - DCONJG(A(l+i,j))*X(i)
               ENDDO
               IF ( nounit ) temp = temp/DCONJG(A(1,j))
            ENDIF
            X(j) = temp
         ENDDO
      ELSE
         kx = kx + (N-1)*Incx
         jx = kx
         DO j = N , 1 , -1
            temp = X(jx)
            ix = kx
            l = 1 - j
            IF ( noconj ) THEN
               DO i = MIN(N,j+K) , j + 1 , -1
                  temp = temp - A(l+i,j)*X(ix)
                  ix = ix - Incx
               ENDDO
               IF ( nounit ) temp = temp/A(1,j)
            ELSE
               DO i = MIN(N,j+K) , j + 1 , -1
                  temp = temp - DCONJG(A(l+i,j))*X(ix)
                  ix = ix - Incx
               ENDDO
               IF ( nounit ) temp = temp/DCONJG(A(1,j))
            ENDIF
            X(jx) = temp
            jx = jx - Incx
            IF ( (N-j)>=K ) kx = kx - Incx
         ENDDO
      ENDIF
!
!
!     End of ZTBSV .
!
      END SUBROUTINE ZTBSV
