!*==chemm.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b CHEMM
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE CHEMM(SIDE,UPLO,M,N,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
!
!       .. Scalar Arguments ..
!       COMPLEX ALPHA,BETA
!       INTEGER LDA,LDB,LDC,M,N
!       CHARACTER SIDE,UPLO
!       ..
!       .. Array Arguments ..
!       COMPLEX A(LDA,*),B(LDB,*),C(LDC,*)
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CHEMM  performs one of the matrix-matrix operations
!>
!>    C := alpha*A*B + beta*C,
!>
!> or
!>
!>    C := alpha*B*A + beta*C,
!>
!> where alpha and beta are scalars, A is an hermitian matrix and  B and
!> C are m by n matrices.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] SIDE
!> \verbatim
!>          SIDE is CHARACTER*1
!>           On entry,  SIDE  specifies whether  the  hermitian matrix  A
!>           appears on the  left or right  in the  operation as follows:
!>
!>              SIDE = 'L' or 'l'   C := alpha*A*B + beta*C,
!>
!>              SIDE = 'R' or 'r'   C := alpha*B*A + beta*C,
!> \endverbatim
!>
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>           On  entry,   UPLO  specifies  whether  the  upper  or  lower
!>           triangular  part  of  the  hermitian  matrix   A  is  to  be
!>           referenced as follows:
!>
!>              UPLO = 'U' or 'u'   Only the upper triangular part of the
!>                                  hermitian matrix is to be referenced.
!>
!>              UPLO = 'L' or 'l'   Only the lower triangular part of the
!>                                  hermitian matrix is to be referenced.
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>           On entry,  M  specifies the number of rows of the matrix  C.
!>           M  must be at least zero.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>           On entry, N specifies the number of columns of the matrix C.
!>           N  must be at least zero.
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
!>          A is COMPLEX array, dimension ( LDA, ka ), where ka is
!>           m  when  SIDE = 'L' or 'l'  and is n  otherwise.
!>           Before entry  with  SIDE = 'L' or 'l',  the  m by m  part of
!>           the array  A  must contain the  hermitian matrix,  such that
!>           when  UPLO = 'U' or 'u', the leading m by m upper triangular
!>           part of the array  A  must contain the upper triangular part
!>           of the  hermitian matrix and the  strictly  lower triangular
!>           part of  A  is not referenced,  and when  UPLO = 'L' or 'l',
!>           the leading  m by m  lower triangular part  of the  array  A
!>           must  contain  the  lower triangular part  of the  hermitian
!>           matrix and the  strictly upper triangular part of  A  is not
!>           referenced.
!>           Before entry  with  SIDE = 'R' or 'r',  the  n by n  part of
!>           the array  A  must contain the  hermitian matrix,  such that
!>           when  UPLO = 'U' or 'u', the leading n by n upper triangular
!>           part of the array  A  must contain the upper triangular part
!>           of the  hermitian matrix and the  strictly  lower triangular
!>           part of  A  is not referenced,  and when  UPLO = 'L' or 'l',
!>           the leading  n by n  lower triangular part  of the  array  A
!>           must  contain  the  lower triangular part  of the  hermitian
!>           matrix and the  strictly upper triangular part of  A  is not
!>           referenced.
!>           Note that the imaginary parts  of the diagonal elements need
!>           not be set, they are assumed to be zero.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>           On entry, LDA specifies the first dimension of A as declared
!>           in the  calling (sub) program. When  SIDE = 'L' or 'l'  then
!>           LDA must be at least  max( 1, m ), otherwise  LDA must be at
!>           least max( 1, n ).
!> \endverbatim
!>
!> \param[in] B
!> \verbatim
!>          B is COMPLEX array, dimension ( LDB, N )
!>           Before entry, the leading  m by n part of the array  B  must
!>           contain the matrix B.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>           On entry, LDB specifies the first dimension of B as declared
!>           in  the  calling  (sub)  program.   LDB  must  be  at  least
!>           max( 1, m ).
!> \endverbatim
!>
!> \param[in] BETA
!> \verbatim
!>          BETA is COMPLEX
!>           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
!>           supplied as zero then C need not be set on input.
!> \endverbatim
!>
!> \param[in,out] C
!> \verbatim
!>          C is COMPLEX array, dimension ( LDC, N )
!>           Before entry, the leading  m by n  part of the array  C must
!>           contain the matrix  C,  except when  beta  is zero, in which
!>           case C need not be set on entry.
!>           On exit, the array  C  is overwritten by the  m by n updated
!>           matrix.
!> \endverbatim
!>
!> \param[in] LDC
!> \verbatim
!>          LDC is INTEGER
!>           On entry, LDC specifies the first dimension of C as declared
!>           in  the  calling  (sub)  program.   LDC  must  be  at  least
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
      SUBROUTINE CHEMM(Side,Uplo,M,N,Alpha,A,Lda,B,Ldb,Beta,C,Ldc)
      IMPLICIT NONE
!*--CHEMM195
!
!  -- Reference BLAS level3 routine (version 3.7.0) --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      COMPLEX Alpha , Beta
      INTEGER Lda , Ldb , Ldc , M , N
      CHARACTER Side , Uplo
!     ..
!     .. Array Arguments ..
      COMPLEX A(Lda,*) , B(Ldb,*) , C(Ldc,*)
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
      INTRINSIC CONJG , MAX , REAL
!     ..
!     .. Local Scalars ..
      COMPLEX temp1 , temp2
      INTEGER i , info , j , k , nrowa
      LOGICAL upper
!     ..
!     .. Parameters ..
      COMPLEX ONE
      PARAMETER (ONE=(1.0E+0,0.0E+0))
      COMPLEX ZERO
      PARAMETER (ZERO=(0.0E+0,0.0E+0))
!     ..
!
!     Set NROWA as the number of rows of A.
!
      IF ( LSAME(Side,'L') ) THEN
         nrowa = M
      ELSE
         nrowa = N
      ENDIF
      upper = LSAME(Uplo,'U')
!
!     Test the input parameters.
!
      info = 0
      IF ( (.NOT.LSAME(Side,'L')) .AND. (.NOT.LSAME(Side,'R')) ) THEN
         info = 1
      ELSEIF ( (.NOT.upper) .AND. (.NOT.LSAME(Uplo,'L')) ) THEN
         info = 2
      ELSEIF ( M<0 ) THEN
         info = 3
      ELSEIF ( N<0 ) THEN
         info = 4
      ELSEIF ( Lda<MAX(1,nrowa) ) THEN
         info = 7
      ELSEIF ( Ldb<MAX(1,M) ) THEN
         info = 9
      ELSEIF ( Ldc<MAX(1,M) ) THEN
         info = 12
      ENDIF
      IF ( info/=0 ) THEN
         CALL XERBLA('CHEMM ',info)
         RETURN
      ENDIF
!
!     Quick return if possible.
!
      IF ( (M==0) .OR. (N==0) .OR. ((Alpha==ZERO) .AND. (Beta==ONE)) )  &
     &     RETURN
!
!     And when  alpha.eq.zero.
!
      IF ( Alpha==ZERO ) THEN
         IF ( Beta==ZERO ) THEN
            DO j = 1 , N
               DO i = 1 , M
                  C(i,j) = ZERO
               ENDDO
            ENDDO
         ELSE
            DO j = 1 , N
               DO i = 1 , M
                  C(i,j) = Beta*C(i,j)
               ENDDO
            ENDDO
         ENDIF
         RETURN
      ENDIF
!
!     Start the operations.
!
      IF ( .NOT.(LSAME(Side,'L')) ) THEN
!
!        Form  C := alpha*B*A + beta*C.
!
         DO j = 1 , N
            temp1 = Alpha*REAL(A(j,j))
            IF ( Beta==ZERO ) THEN
               DO i = 1 , M
                  C(i,j) = temp1*B(i,j)
               ENDDO
            ELSE
               DO i = 1 , M
                  C(i,j) = Beta*C(i,j) + temp1*B(i,j)
               ENDDO
            ENDIF
            DO k = 1 , j - 1
               IF ( upper ) THEN
                  temp1 = Alpha*A(k,j)
               ELSE
                  temp1 = Alpha*CONJG(A(j,k))
               ENDIF
               DO i = 1 , M
                  C(i,j) = C(i,j) + temp1*B(i,k)
               ENDDO
            ENDDO
            DO k = j + 1 , N
               IF ( upper ) THEN
                  temp1 = Alpha*CONJG(A(j,k))
               ELSE
                  temp1 = Alpha*A(k,j)
               ENDIF
               DO i = 1 , M
                  C(i,j) = C(i,j) + temp1*B(i,k)
               ENDDO
            ENDDO
         ENDDO
!
!        Form  C := alpha*A*B + beta*C.
!
      ELSEIF ( upper ) THEN
         DO j = 1 , N
            DO i = 1 , M
               temp1 = Alpha*B(i,j)
               temp2 = ZERO
               DO k = 1 , i - 1
                  C(k,j) = C(k,j) + temp1*A(k,i)
                  temp2 = temp2 + B(k,j)*CONJG(A(k,i))
               ENDDO
               IF ( Beta==ZERO ) THEN
                  C(i,j) = temp1*REAL(A(i,i)) + Alpha*temp2
               ELSE
                  C(i,j) = Beta*C(i,j) + temp1*REAL(A(i,i))             &
     &                     + Alpha*temp2
               ENDIF
            ENDDO
         ENDDO
      ELSE
         DO j = 1 , N
            DO i = M , 1 , -1
               temp1 = Alpha*B(i,j)
               temp2 = ZERO
               DO k = i + 1 , M
                  C(k,j) = C(k,j) + temp1*A(k,i)
                  temp2 = temp2 + B(k,j)*CONJG(A(k,i))
               ENDDO
               IF ( Beta==ZERO ) THEN
                  C(i,j) = temp1*REAL(A(i,i)) + Alpha*temp2
               ELSE
                  C(i,j) = Beta*C(i,j) + temp1*REAL(A(i,i))             &
     &                     + Alpha*temp2
               ENDIF
            ENDDO
         ENDDO
      ENDIF
!
!
!     End of CHEMM .
!
      END SUBROUTINE CHEMM
