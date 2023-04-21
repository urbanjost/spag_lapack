!*==zher2k.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b ZHER2K
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZHER2K(UPLO,TRANS,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
!
!       .. Scalar Arguments ..
!       COMPLEX*16 ALPHA
!       DOUBLE PRECISION BETA
!       INTEGER K,LDA,LDB,LDC,N
!       CHARACTER TRANS,UPLO
!       ..
!       .. Array Arguments ..
!       COMPLEX*16 A(LDA,*),B(LDB,*),C(LDC,*)
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZHER2K  performs one of the hermitian rank 2k operations
!>
!>    C := alpha*A*B**H + conjg( alpha )*B*A**H + beta*C,
!>
!> or
!>
!>    C := alpha*A**H*B + conjg( alpha )*B**H*A + beta*C,
!>
!> where  alpha and beta  are scalars with  beta  real,  C is an  n by n
!> hermitian matrix and  A and B  are  n by k matrices in the first case
!> and  k by n  matrices in the second case.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>           On  entry,   UPLO  specifies  whether  the  upper  or  lower
!>           triangular  part  of the  array  C  is to be  referenced  as
!>           follows:
!>
!>              UPLO = 'U' or 'u'   Only the  upper triangular part of  C
!>                                  is to be referenced.
!>
!>              UPLO = 'L' or 'l'   Only the  lower triangular part of  C
!>                                  is to be referenced.
!> \endverbatim
!>
!> \param[in] TRANS
!> \verbatim
!>          TRANS is CHARACTER*1
!>           On entry,  TRANS  specifies the operation to be performed as
!>           follows:
!>
!>              TRANS = 'N' or 'n'    C := alpha*A*B**H          +
!>                                         conjg( alpha )*B*A**H +
!>                                         beta*C.
!>
!>              TRANS = 'C' or 'c'    C := alpha*A**H*B          +
!>                                         conjg( alpha )*B**H*A +
!>                                         beta*C.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>           On entry,  N specifies the order of the matrix C.  N must be
!>           at least zero.
!> \endverbatim
!>
!> \param[in] K
!> \verbatim
!>          K is INTEGER
!>           On entry with  TRANS = 'N' or 'n',  K  specifies  the number
!>           of  columns  of the  matrices  A and B,  and on  entry  with
!>           TRANS = 'C' or 'c',  K  specifies  the number of rows of the
!>           matrices  A and B.  K must be at least zero.
!> \endverbatim
!>
!> \param[in] ALPHA
!> \verbatim
!>          ALPHA is COMPLEX*16 .
!>           On entry, ALPHA specifies the scalar alpha.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension ( LDA, ka ), where ka is
!>           k  when  TRANS = 'N' or 'n',  and is  n  otherwise.
!>           Before entry with  TRANS = 'N' or 'n',  the  leading  n by k
!>           part of the array  A  must contain the matrix  A,  otherwise
!>           the leading  k by n  part of the array  A  must contain  the
!>           matrix A.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>           On entry, LDA specifies the first dimension of A as declared
!>           in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n'
!>           then  LDA must be at least  max( 1, n ), otherwise  LDA must
!>           be at least  max( 1, k ).
!> \endverbatim
!>
!> \param[in] B
!> \verbatim
!>          B is COMPLEX*16 array, dimension ( LDB, kb ), where kb is
!>           k  when  TRANS = 'N' or 'n',  and is  n  otherwise.
!>           Before entry with  TRANS = 'N' or 'n',  the  leading  n by k
!>           part of the array  B  must contain the matrix  B,  otherwise
!>           the leading  k by n  part of the array  B  must contain  the
!>           matrix B.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>           On entry, LDB specifies the first dimension of B as declared
!>           in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n'
!>           then  LDB must be at least  max( 1, n ), otherwise  LDB must
!>           be at least  max( 1, k ).
!>           Unchanged on exit.
!> \endverbatim
!>
!> \param[in] BETA
!> \verbatim
!>          BETA is DOUBLE PRECISION .
!>           On entry, BETA specifies the scalar beta.
!> \endverbatim
!>
!> \param[in,out] C
!> \verbatim
!>          C is COMPLEX*16 array, dimension ( LDC, N )
!>           Before entry  with  UPLO = 'U' or 'u',  the leading  n by n
!>           upper triangular part of the array C must contain the upper
!>           triangular part  of the  hermitian matrix  and the strictly
!>           lower triangular part of C is not referenced.  On exit, the
!>           upper triangular part of the array  C is overwritten by the
!>           upper triangular part of the updated matrix.
!>           Before entry  with  UPLO = 'L' or 'l',  the leading  n by n
!>           lower triangular part of the array C must contain the lower
!>           triangular part  of the  hermitian matrix  and the strictly
!>           upper triangular part of C is not referenced.  On exit, the
!>           lower triangular part of the array  C is overwritten by the
!>           lower triangular part of the updated matrix.
!>           Note that the imaginary parts of the diagonal elements need
!>           not be set,  they are assumed to be zero,  and on exit they
!>           are set to zero.
!> \endverbatim
!>
!> \param[in] LDC
!> \verbatim
!>          LDC is INTEGER
!>           On entry, LDC specifies the first dimension of C as declared
!>           in  the  calling  (sub)  program.   LDC  must  be  at  least
!>           max( 1, n ).
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
!> \ingroup complex16_blas_level3
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
!>
!>  -- Modified 8-Nov-93 to set C(J,J) to DBLE( C(J,J) ) when BETA = 1.
!>     Ed Anderson, Cray Research Inc.
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE ZHER2K(Uplo,Trans,N,K,Alpha,A,Lda,B,Ldb,Beta,C,Ldc)
      IMPLICIT NONE
!*--ZHER2K202
!
!  -- Reference BLAS level3 routine (version 3.7.0) --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      COMPLEX*16 Alpha
      DOUBLE PRECISION Beta
      INTEGER K , Lda , Ldb , Ldc , N
      CHARACTER Trans , Uplo
!     ..
!     .. Array Arguments ..
      COMPLEX*16 A(Lda,*) , B(Ldb,*) , C(Ldc,*)
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
      INTRINSIC DBLE , DCONJG , MAX
!     ..
!     .. Local Scalars ..
      COMPLEX*16 temp1 , temp2
      INTEGER i , info , j , l , nrowa
      LOGICAL upper
!     ..
!     .. Parameters ..
      DOUBLE PRECISION ONE
      PARAMETER (ONE=1.0D+0)
      COMPLEX*16 ZERO
      PARAMETER (ZERO=(0.0D+0,0.0D+0))
!     ..
!
!     Test the input parameters.
!
      IF ( LSAME(Trans,'N') ) THEN
         nrowa = N
      ELSE
         nrowa = K
      ENDIF
      upper = LSAME(Uplo,'U')
!
      info = 0
      IF ( (.NOT.upper) .AND. (.NOT.LSAME(Uplo,'L')) ) THEN
         info = 1
      ELSEIF ( (.NOT.LSAME(Trans,'N')) .AND. (.NOT.LSAME(Trans,'C')) )  &
     &         THEN
         info = 2
      ELSEIF ( N<0 ) THEN
         info = 3
      ELSEIF ( K<0 ) THEN
         info = 4
      ELSEIF ( Lda<MAX(1,nrowa) ) THEN
         info = 7
      ELSEIF ( Ldb<MAX(1,nrowa) ) THEN
         info = 9
      ELSEIF ( Ldc<MAX(1,N) ) THEN
         info = 12
      ENDIF
      IF ( info/=0 ) THEN
         CALL XERBLA('ZHER2K',info)
         RETURN
      ENDIF
!
!     Quick return if possible.
!
      IF ( (N==0) .OR. (((Alpha==ZERO) .OR. (K==0)) .AND. (Beta==ONE)) )&
     &     RETURN
!
!     And when  alpha.eq.zero.
!
      IF ( Alpha==ZERO ) THEN
         IF ( upper ) THEN
            IF ( Beta==DBLE(ZERO) ) THEN
               DO j = 1 , N
                  DO i = 1 , j
                     C(i,j) = ZERO
                  ENDDO
               ENDDO
            ELSE
               DO j = 1 , N
                  DO i = 1 , j - 1
                     C(i,j) = Beta*C(i,j)
                  ENDDO
                  C(j,j) = Beta*DBLE(C(j,j))
               ENDDO
            ENDIF
         ELSEIF ( Beta==DBLE(ZERO) ) THEN
            DO j = 1 , N
               DO i = j , N
                  C(i,j) = ZERO
               ENDDO
            ENDDO
         ELSE
            DO j = 1 , N
               C(j,j) = Beta*DBLE(C(j,j))
               DO i = j + 1 , N
                  C(i,j) = Beta*C(i,j)
               ENDDO
            ENDDO
         ENDIF
         RETURN
      ENDIF
!
!     Start the operations.
!
      IF ( LSAME(Trans,'N') ) THEN
!
!        Form  C := alpha*A*B**H + conjg( alpha )*B*A**H +
!                   C.
!
         IF ( upper ) THEN
            DO j = 1 , N
               IF ( Beta==DBLE(ZERO) ) THEN
                  DO i = 1 , j
                     C(i,j) = ZERO
                  ENDDO
               ELSEIF ( Beta/=ONE ) THEN
                  DO i = 1 , j - 1
                     C(i,j) = Beta*C(i,j)
                  ENDDO
                  C(j,j) = Beta*DBLE(C(j,j))
               ELSE
                  C(j,j) = DBLE(C(j,j))
               ENDIF
               DO l = 1 , K
                  IF ( (A(j,l)/=ZERO) .OR. (B(j,l)/=ZERO) ) THEN
                     temp1 = Alpha*DCONJG(B(j,l))
                     temp2 = DCONJG(Alpha*A(j,l))
                     DO i = 1 , j - 1
                        C(i,j) = C(i,j) + A(i,l)*temp1 + B(i,l)*temp2
                     ENDDO
                     C(j,j) = DBLE(C(j,j))                              &
     &                        + DBLE(A(j,l)*temp1+B(j,l)*temp2)
                  ENDIF
               ENDDO
            ENDDO
         ELSE
            DO j = 1 , N
               IF ( Beta==DBLE(ZERO) ) THEN
                  DO i = j , N
                     C(i,j) = ZERO
                  ENDDO
               ELSEIF ( Beta/=ONE ) THEN
                  DO i = j + 1 , N
                     C(i,j) = Beta*C(i,j)
                  ENDDO
                  C(j,j) = Beta*DBLE(C(j,j))
               ELSE
                  C(j,j) = DBLE(C(j,j))
               ENDIF
               DO l = 1 , K
                  IF ( (A(j,l)/=ZERO) .OR. (B(j,l)/=ZERO) ) THEN
                     temp1 = Alpha*DCONJG(B(j,l))
                     temp2 = DCONJG(Alpha*A(j,l))
                     DO i = j + 1 , N
                        C(i,j) = C(i,j) + A(i,l)*temp1 + B(i,l)*temp2
                     ENDDO
                     C(j,j) = DBLE(C(j,j))                              &
     &                        + DBLE(A(j,l)*temp1+B(j,l)*temp2)
                  ENDIF
               ENDDO
            ENDDO
         ENDIF
!
!        Form  C := alpha*A**H*B + conjg( alpha )*B**H*A +
!                   C.
!
      ELSEIF ( upper ) THEN
         DO j = 1 , N
            DO i = 1 , j
               temp1 = ZERO
               temp2 = ZERO
               DO l = 1 , K
                  temp1 = temp1 + DCONJG(A(l,i))*B(l,j)
                  temp2 = temp2 + DCONJG(B(l,i))*A(l,j)
               ENDDO
               IF ( i==j ) THEN
                  IF ( Beta==DBLE(ZERO) ) THEN
                     C(j,j) = DBLE(Alpha*temp1+DCONJG(Alpha)*temp2)
                  ELSE
                     C(j,j) = Beta*DBLE(C(j,j))                         &
     &                        + DBLE(Alpha*temp1+DCONJG(Alpha)*temp2)
                  ENDIF
               ELSEIF ( Beta==DBLE(ZERO) ) THEN
                  C(i,j) = Alpha*temp1 + DCONJG(Alpha)*temp2
               ELSE
                  C(i,j) = Beta*C(i,j) + Alpha*temp1 + DCONJG(Alpha)    &
     &                     *temp2
               ENDIF
            ENDDO
         ENDDO
      ELSE
         DO j = 1 , N
            DO i = j , N
               temp1 = ZERO
               temp2 = ZERO
               DO l = 1 , K
                  temp1 = temp1 + DCONJG(A(l,i))*B(l,j)
                  temp2 = temp2 + DCONJG(B(l,i))*A(l,j)
               ENDDO
               IF ( i==j ) THEN
                  IF ( Beta==DBLE(ZERO) ) THEN
                     C(j,j) = DBLE(Alpha*temp1+DCONJG(Alpha)*temp2)
                  ELSE
                     C(j,j) = Beta*DBLE(C(j,j))                         &
     &                        + DBLE(Alpha*temp1+DCONJG(Alpha)*temp2)
                  ENDIF
               ELSEIF ( Beta==DBLE(ZERO) ) THEN
                  C(i,j) = Alpha*temp1 + DCONJG(Alpha)*temp2
               ELSE
                  C(i,j) = Beta*C(i,j) + Alpha*temp1 + DCONJG(Alpha)    &
     &                     *temp2
               ENDIF
            ENDDO
         ENDDO
      ENDIF
!
!
!     End of ZHER2K.
!
      END SUBROUTINE ZHER2K
