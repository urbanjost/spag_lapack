!*==dtrmm.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b DTRMM
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DTRMM(SIDE,UPLO,TRANSA,DIAG,M,N,ALPHA,A,LDA,B,LDB)
!
!       .. Scalar Arguments ..
!       DOUBLE PRECISION ALPHA
!       INTEGER LDA,LDB,M,N
!       CHARACTER DIAG,SIDE,TRANSA,UPLO
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION A(LDA,*),B(LDB,*)
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DTRMM  performs one of the matrix-matrix operations
!>
!>    B := alpha*op( A )*B,   or   B := alpha*B*op( A ),
!>
!> where  alpha  is a scalar,  B  is an m by n matrix,  A  is a unit, or
!> non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
!>
!>    op( A ) = A   or   op( A ) = A**T.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] SIDE
!> \verbatim
!>          SIDE is CHARACTER*1
!>           On entry,  SIDE specifies whether  op( A ) multiplies B from
!>           the left or right as follows:
!>
!>              SIDE = 'L' or 'l'   B := alpha*op( A )*B.
!>
!>              SIDE = 'R' or 'r'   B := alpha*B*op( A ).
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
!>              TRANSA = 'C' or 'c'   op( A ) = A**T.
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
!>          ALPHA is DOUBLE PRECISION.
!>           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
!>           zero then  A is not referenced and  B need not be set before
!>           entry.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>           A is DOUBLE PRECISION array, dimension ( LDA, k ), where k is m
!>           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.
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
!>          B is DOUBLE PRECISION array, dimension ( LDB, N )
!>           Before entry,  the leading  m by n part of the array  B must
!>           contain the matrix  B,  and  on exit  is overwritten  by the
!>           transformed matrix.
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
!> \ingroup double_blas_level3
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
      SUBROUTINE DTRMM(Side,Uplo,Transa,Diag,M,N,Alpha,A,Lda,B,Ldb)
      IMPLICIT NONE
!*--DTRMM181
!
!  -- Reference BLAS level3 routine (version 3.7.0) --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      DOUBLE PRECISION Alpha
      INTEGER Lda , Ldb , M , N
      CHARACTER Diag , Side , Transa , Uplo
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION A(Lda,*) , B(Ldb,*)
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
      INTRINSIC MAX
!     ..
!     .. Local Scalars ..
      DOUBLE PRECISION temp
      INTEGER i , info , j , k , nrowa
      LOGICAL lside , nounit , upper
!     ..
!     .. Parameters ..
      DOUBLE PRECISION ONE , ZERO
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
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
         CALL XERBLA('DTRMM ',info)
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
!           Form  B := alpha*A*B.
!
            IF ( upper ) THEN
               DO j = 1 , N
                  DO k = 1 , M
                     IF ( B(k,j)/=ZERO ) THEN
                        temp = Alpha*B(k,j)
                        DO i = 1 , k - 1
                           B(i,j) = B(i,j) + temp*A(i,k)
                        ENDDO
                        IF ( nounit ) temp = temp*A(k,k)
                        B(k,j) = temp
                     ENDIF
                  ENDDO
               ENDDO
            ELSE
               DO j = 1 , N
                  DO k = M , 1 , -1
                     IF ( B(k,j)/=ZERO ) THEN
                        temp = Alpha*B(k,j)
                        B(k,j) = temp
                        IF ( nounit ) B(k,j) = B(k,j)*A(k,k)
                        DO i = k + 1 , M
                           B(i,j) = B(i,j) + temp*A(i,k)
                        ENDDO
                     ENDIF
                  ENDDO
               ENDDO
            ENDIF
!
!           Form  B := alpha*A**T*B.
!
         ELSEIF ( upper ) THEN
            DO j = 1 , N
               DO i = M , 1 , -1
                  temp = B(i,j)
                  IF ( nounit ) temp = temp*A(i,i)
                  DO k = 1 , i - 1
                     temp = temp + A(k,i)*B(k,j)
                  ENDDO
                  B(i,j) = Alpha*temp
               ENDDO
            ENDDO
         ELSE
            DO j = 1 , N
               DO i = 1 , M
                  temp = B(i,j)
                  IF ( nounit ) temp = temp*A(i,i)
                  DO k = i + 1 , M
                     temp = temp + A(k,i)*B(k,j)
                  ENDDO
                  B(i,j) = Alpha*temp
               ENDDO
            ENDDO
         ENDIF
      ELSEIF ( LSAME(Transa,'N') ) THEN
!
!           Form  B := alpha*B*A.
!
         IF ( upper ) THEN
            DO j = N , 1 , -1
               temp = Alpha
               IF ( nounit ) temp = temp*A(j,j)
               DO i = 1 , M
                  B(i,j) = temp*B(i,j)
               ENDDO
               DO k = 1 , j - 1
                  IF ( A(k,j)/=ZERO ) THEN
                     temp = Alpha*A(k,j)
                     DO i = 1 , M
                        B(i,j) = B(i,j) + temp*B(i,k)
                     ENDDO
                  ENDIF
               ENDDO
            ENDDO
         ELSE
            DO j = 1 , N
               temp = Alpha
               IF ( nounit ) temp = temp*A(j,j)
               DO i = 1 , M
                  B(i,j) = temp*B(i,j)
               ENDDO
               DO k = j + 1 , N
                  IF ( A(k,j)/=ZERO ) THEN
                     temp = Alpha*A(k,j)
                     DO i = 1 , M
                        B(i,j) = B(i,j) + temp*B(i,k)
                     ENDDO
                  ENDIF
               ENDDO
            ENDDO
         ENDIF
!
!           Form  B := alpha*B*A**T.
!
      ELSEIF ( upper ) THEN
         DO k = 1 , N
            DO j = 1 , k - 1
               IF ( A(j,k)/=ZERO ) THEN
                  temp = Alpha*A(j,k)
                  DO i = 1 , M
                     B(i,j) = B(i,j) + temp*B(i,k)
                  ENDDO
               ENDIF
            ENDDO
            temp = Alpha
            IF ( nounit ) temp = temp*A(k,k)
            IF ( temp/=ONE ) THEN
               DO i = 1 , M
                  B(i,k) = temp*B(i,k)
               ENDDO
            ENDIF
         ENDDO
      ELSE
         DO k = N , 1 , -1
            DO j = k + 1 , N
               IF ( A(j,k)/=ZERO ) THEN
                  temp = Alpha*A(j,k)
                  DO i = 1 , M
                     B(i,j) = B(i,j) + temp*B(i,k)
                  ENDDO
               ENDIF
            ENDDO
            temp = Alpha
            IF ( nounit ) temp = temp*A(k,k)
            IF ( temp/=ONE ) THEN
               DO i = 1 , M
                  B(i,k) = temp*B(i,k)
               ENDDO
            ENDIF
         ENDDO
      ENDIF
!
!
!     End of DTRMM .
!
      END SUBROUTINE DTRMM
