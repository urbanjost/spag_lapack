!*==ctrsm.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b CTRSM
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE CTRSM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
!
!       .. Scalar Arguments ..
!       COMPLEX ALPHA
!       INTEGER LDA,LDB,M,N
!       CHARACTER DIAG,SIDE,TRANSA,UPLO
!       ..
!       .. Array Arguments ..
!       COMPLEX A(LDA,*),B(LDB,*)
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CTRSM  solves one of the matrix equations
!>
!>    op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
!>
!> where alpha is a scalar, X and B are m by n matrices, A is a unit, or
!> non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
!>
!>    op( A ) = A   or   op( A ) = A**T   or   op( A ) = A**H.
!>
!> The matrix X is overwritten on B.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] SIDE
!> \verbatim
!>          SIDE is CHARACTER*1
!>           On entry, SIDE specifies whether op( A ) appears on the left
!>           or right of X as follows:
!>
!>              SIDE = 'L' or 'l'   op( A )*X = alpha*B.
!>
!>              SIDE = 'R' or 'r'   X*op( A ) = alpha*B.
!> \endverbatim
!>
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>           On entry, UPLO specifies whether the matrix A is an upper or
!>           lower triangular matrix as follows:
!>
!>              UPLO = 'U' or 'u'   A is an upper triangular matrix.
!>
!>              UPLO = 'L' or 'l'   A is a lower triangular matrix.
!> \endverbatim
!>
!> \param[in] TRANSA
!> \verbatim
!>          TRANSA is CHARACTER*1
!>           On entry, TRANSA specifies the form of op( A ) to be used in
!>           the matrix multiplication as follows:
!>
!>              TRANSA = 'N' or 'n'   op( A ) = A.
!>
!>              TRANSA = 'T' or 't'   op( A ) = A**T.
!>
!>              TRANSA = 'C' or 'c'   op( A ) = A**H.
!> \endverbatim
!>
!> \param[in] DIAG
!> \verbatim
!>          DIAG is CHARACTER*1
!>           On entry, DIAG specifies whether or not A is unit triangular
!>           as follows:
!>
!>              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
!>
!>              DIAG = 'N' or 'n'   A is not assumed to be unit
!>                                  triangular.
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>           On entry, M specifies the number of rows of B. M must be at
!>           least zero.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>           On entry, N specifies the number of columns of B.  N must be
!>           at least zero.
!> \endverbatim
!>
!> \param[in] ALPHA
!> \verbatim
!>          ALPHA is COMPLEX
!>           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
!>           zero then  A is not referenced and  B need not be set before
!>           entry.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is COMPLEX array, dimension ( LDA, k ),
!>           where k is m when SIDE = 'L' or 'l'
!>             and k is n when SIDE = 'R' or 'r'.
!>           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
!>           upper triangular part of the array  A must contain the upper
!>           triangular matrix  and the strictly lower triangular part of
!>           A is not referenced.
!>           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
!>           lower triangular part of the array  A must contain the lower
!>           triangular matrix  and the strictly upper triangular part of
!>           A is not referenced.
!>           Note that when  DIAG = 'U' or 'u',  the diagonal elements of
!>           A  are not referenced either,  but are assumed to be  unity.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>           On entry, LDA specifies the first dimension of A as declared
!>           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
!>           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
!>           then LDA must be at least max( 1, n ).
!> \endverbatim
!>
!> \param[in,out] B
!> \verbatim
!>          B is COMPLEX array, dimension ( LDB, N )
!>           Before entry,  the leading  m by n part of the array  B must
!>           contain  the  right-hand  side  matrix  B,  and  on exit  is
!>           overwritten by the solution matrix  X.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>           On entry, LDB specifies the first dimension of B as declared
!>           in  the  calling  (sub)  program.   LDB  must  be  at  least
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
!> \ingroup complex_blas_level3
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  Level 3 Blas routine.
!>
!>  -- Written on 8-February-1989.
!>     Jack Dongarra, Argonne National Laboratory.
!>     Iain Duff, AERE Harwell.
!>     Jeremy Du Croz, Numerical Algorithms Group Ltd.
!>     Sven Hammarling, Numerical Algorithms Group Ltd.
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE CTRSM(Side,Uplo,Transa,Diag,M,N,Alpha,A,Lda,B,Ldb)
      IMPLICIT NONE
!*--CTRSM184
!
!  -- Reference BLAS level3 routine (version 3.7.0) --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      COMPLEX Alpha
      INTEGER Lda , Ldb , M , N
      CHARACTER Diag , Side , Transa , Uplo
!     ..
!     .. Array Arguments ..
      COMPLEX A(Lda,*) , B(Ldb,*)
!     ..
!
!  =====================================================================
!
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
!     .. Local Scalars ..
      COMPLEX temp
      INTEGER i , info , j , k , nrowa
      LOGICAL lside , noconj , nounit , upper
!     ..
!     .. Parameters ..
      COMPLEX ONE
      PARAMETER (ONE=(1.0E+0,0.0E+0))
      COMPLEX ZERO
      PARAMETER (ZERO=(0.0E+0,0.0E+0))
!     ..
!
!     Test the input parameters.
!
      lside = LSAME(Side,'L')
      IF ( lside ) THEN
         nrowa = M
      ELSE
         nrowa = N
      ENDIF
      noconj = LSAME(Transa,'T')
      nounit = LSAME(Diag,'N')
      upper = LSAME(Uplo,'U')
!
      info = 0
      IF ( (.NOT.lside) .AND. (.NOT.LSAME(Side,'R')) ) THEN
         info = 1
      ELSEIF ( (.NOT.upper) .AND. (.NOT.LSAME(Uplo,'L')) ) THEN
         info = 2
      ELSEIF ( (.NOT.LSAME(Transa,'N')) .AND. (.NOT.LSAME(Transa,'T'))  &
     &         .AND. (.NOT.LSAME(Transa,'C')) ) THEN
         info = 3
      ELSEIF ( (.NOT.LSAME(Diag,'U')) .AND. (.NOT.LSAME(Diag,'N')) )    &
     &         THEN
         info = 4
      ELSEIF ( M<0 ) THEN
         info = 5
      ELSEIF ( N<0 ) THEN
         info = 6
      ELSEIF ( Lda<MAX(1,nrowa) ) THEN
         info = 9
      ELSEIF ( Ldb<MAX(1,M) ) THEN
         info = 11
      ENDIF
      IF ( info/=0 ) THEN
         CALL XERBLA('CTRSM ',info)
         RETURN
      ENDIF
!
!     Quick return if possible.
!
      IF ( M==0 .OR. N==0 ) RETURN
!
!     And when  alpha.eq.zero.
!
      IF ( Alpha==ZERO ) THEN
         DO j = 1 , N
            DO i = 1 , M
               B(i,j) = ZERO
            ENDDO
         ENDDO
         RETURN
      ENDIF
!
!     Start the operations.
!
      IF ( lside ) THEN
         IF ( LSAME(Transa,'N') ) THEN
!
!           Form  B := alpha*inv( A )*B.
!
            IF ( upper ) THEN
               DO j = 1 , N
                  IF ( Alpha/=ONE ) THEN
                     DO i = 1 , M
                        B(i,j) = Alpha*B(i,j)
                     ENDDO
                  ENDIF
                  DO k = M , 1 , -1
                     IF ( B(k,j)/=ZERO ) THEN
                        IF ( nounit ) B(k,j) = B(k,j)/A(k,k)
                        DO i = 1 , k - 1
                           B(i,j) = B(i,j) - B(k,j)*A(i,k)
                        ENDDO
                     ENDIF
                  ENDDO
               ENDDO
            ELSE
               DO j = 1 , N
                  IF ( Alpha/=ONE ) THEN
                     DO i = 1 , M
                        B(i,j) = Alpha*B(i,j)
                     ENDDO
                  ENDIF
                  DO k = 1 , M
                     IF ( B(k,j)/=ZERO ) THEN
                        IF ( nounit ) B(k,j) = B(k,j)/A(k,k)
                        DO i = k + 1 , M
                           B(i,j) = B(i,j) - B(k,j)*A(i,k)
                        ENDDO
                     ENDIF
                  ENDDO
               ENDDO
            ENDIF
!
!           Form  B := alpha*inv( A**T )*B
!           or    B := alpha*inv( A**H )*B.
!
         ELSEIF ( upper ) THEN
            DO j = 1 , N
               DO i = 1 , M
                  temp = Alpha*B(i,j)
                  IF ( noconj ) THEN
                     DO k = 1 , i - 1
                        temp = temp - A(k,i)*B(k,j)
                     ENDDO
                     IF ( nounit ) temp = temp/A(i,i)
                  ELSE
                     DO k = 1 , i - 1
                        temp = temp - CONJG(A(k,i))*B(k,j)
                     ENDDO
                     IF ( nounit ) temp = temp/CONJG(A(i,i))
                  ENDIF
                  B(i,j) = temp
               ENDDO
            ENDDO
         ELSE
            DO j = 1 , N
               DO i = M , 1 , -1
                  temp = Alpha*B(i,j)
                  IF ( noconj ) THEN
                     DO k = i + 1 , M
                        temp = temp - A(k,i)*B(k,j)
                     ENDDO
                     IF ( nounit ) temp = temp/A(i,i)
                  ELSE
                     DO k = i + 1 , M
                        temp = temp - CONJG(A(k,i))*B(k,j)
                     ENDDO
                     IF ( nounit ) temp = temp/CONJG(A(i,i))
                  ENDIF
                  B(i,j) = temp
               ENDDO
            ENDDO
         ENDIF
      ELSEIF ( LSAME(Transa,'N') ) THEN
!
!           Form  B := alpha*B*inv( A ).
!
         IF ( upper ) THEN
            DO j = 1 , N
               IF ( Alpha/=ONE ) THEN
                  DO i = 1 , M
                     B(i,j) = Alpha*B(i,j)
                  ENDDO
               ENDIF
               DO k = 1 , j - 1
                  IF ( A(k,j)/=ZERO ) THEN
                     DO i = 1 , M
                        B(i,j) = B(i,j) - A(k,j)*B(i,k)
                     ENDDO
                  ENDIF
               ENDDO
               IF ( nounit ) THEN
                  temp = ONE/A(j,j)
                  DO i = 1 , M
                     B(i,j) = temp*B(i,j)
                  ENDDO
               ENDIF
            ENDDO
         ELSE
            DO j = N , 1 , -1
               IF ( Alpha/=ONE ) THEN
                  DO i = 1 , M
                     B(i,j) = Alpha*B(i,j)
                  ENDDO
               ENDIF
               DO k = j + 1 , N
                  IF ( A(k,j)/=ZERO ) THEN
                     DO i = 1 , M
                        B(i,j) = B(i,j) - A(k,j)*B(i,k)
                     ENDDO
                  ENDIF
               ENDDO
               IF ( nounit ) THEN
                  temp = ONE/A(j,j)
                  DO i = 1 , M
                     B(i,j) = temp*B(i,j)
                  ENDDO
               ENDIF
            ENDDO
         ENDIF
!
!           Form  B := alpha*B*inv( A**T )
!           or    B := alpha*B*inv( A**H ).
!
      ELSEIF ( upper ) THEN
         DO k = N , 1 , -1
            IF ( nounit ) THEN
               IF ( noconj ) THEN
                  temp = ONE/A(k,k)
               ELSE
                  temp = ONE/CONJG(A(k,k))
               ENDIF
               DO i = 1 , M
                  B(i,k) = temp*B(i,k)
               ENDDO
            ENDIF
            DO j = 1 , k - 1
               IF ( A(j,k)/=ZERO ) THEN
                  IF ( noconj ) THEN
                     temp = A(j,k)
                  ELSE
                     temp = CONJG(A(j,k))
                  ENDIF
                  DO i = 1 , M
                     B(i,j) = B(i,j) - temp*B(i,k)
                  ENDDO
               ENDIF
            ENDDO
            IF ( Alpha/=ONE ) THEN
               DO i = 1 , M
                  B(i,k) = Alpha*B(i,k)
               ENDDO
            ENDIF
         ENDDO
      ELSE
         DO k = 1 , N
            IF ( nounit ) THEN
               IF ( noconj ) THEN
                  temp = ONE/A(k,k)
               ELSE
                  temp = ONE/CONJG(A(k,k))
               ENDIF
               DO i = 1 , M
                  B(i,k) = temp*B(i,k)
               ENDDO
            ENDIF
            DO j = k + 1 , N
               IF ( A(j,k)/=ZERO ) THEN
                  IF ( noconj ) THEN
                     temp = A(j,k)
                  ELSE
                     temp = CONJG(A(j,k))
                  ENDIF
                  DO i = 1 , M
                     B(i,j) = B(i,j) - temp*B(i,k)
                  ENDDO
               ENDIF
            ENDDO
            IF ( Alpha/=ONE ) THEN
               DO i = 1 , M
                  B(i,k) = Alpha*B(i,k)
               ENDDO
            ENDIF
         ENDDO
      ENDIF
!
!
!     End of CTRSM .
!
      END SUBROUTINE CTRSM
