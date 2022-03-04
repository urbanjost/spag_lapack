!*==zherk.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b ZHERK
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZHERK(UPLO,TRANS,N,K,ALPHA,A,LDA,BETA,C,LDC)
!
!       .. Scalar Arguments ..
!       DOUBLE PRECISION ALPHA,BETA
!       INTEGER K,LDA,LDC,N
!       CHARACTER TRANS,UPLO
!       ..
!       .. Array Arguments ..
!       COMPLEX*16 A(LDA,*),C(LDC,*)
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZHERK  performs one of the hermitian rank k operations
!>
!>    C := alpha*A*A**H + beta*C,
!>
!> or
!>
!>    C := alpha*A**H*A + beta*C,
!>
!> where  alpha and beta  are  real scalars,  C is an  n by n  hermitian
!> matrix and  A  is an  n by k  matrix in the  first case and a  k by n
!> matrix in the second case.
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
!>              TRANS = 'N' or 'n'   C := alpha*A*A**H + beta*C.
!>
!>              TRANS = 'C' or 'c'   C := alpha*A**H*A + beta*C.
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
!>           of  columns   of  the   matrix   A,   and  on   entry   with
!>           TRANS = 'C' or 'c',  K  specifies  the number of rows of the
!>           matrix A.  K must be at least zero.
!> \endverbatim
!>
!> \param[in] ALPHA
!> \verbatim
!>          ALPHA is DOUBLE PRECISION .
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
!> \param[in] BETA
!> \verbatim
!>          BETA is DOUBLE PRECISION.
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
      SUBROUTINE ZHERK(Uplo,Trans,N,K,Alpha,A,Lda,Beta,C,Ldc)
      USE F77KINDS                        
      USE S_LSAME
      USE S_XERBLA
      IMPLICIT NONE
!*--ZHERK180
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Uplo
      CHARACTER :: Trans
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: K
      REAL(R8KIND) , INTENT(IN) :: Alpha
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL(R8KIND) , INTENT(IN) :: Beta
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldc,*) :: C
      INTEGER , INTENT(IN) :: Ldc
!
! Local variable declarations rewritten by SPAG
!
      INTEGER :: i , info , j , l , nrowa
      REAL(R8KIND) :: rtemp
      COMPLEX(CX16KIND) :: temp
      LOGICAL :: upper
!
! End of declarations rewritten by SPAG
!
!     ..
!     .. Array Arguments ..
!     ..
!
!  =====================================================================
!
!     .. External Functions ..
!     ..
!     .. External Subroutines ..
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. Local Scalars ..
!     ..
!     .. Parameters ..
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
      ELSEIF ( Ldc<MAX(1,N) ) THEN
         info = 10
      ENDIF
      IF ( info/=0 ) THEN
         CALL XERBLA('ZHERK ',info)
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
            IF ( Beta==ZERO ) THEN
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
         ELSEIF ( Beta==ZERO ) THEN
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
!        Form  C := alpha*A*A**H + beta*C.
!
         IF ( upper ) THEN
            DO j = 1 , N
               IF ( Beta==ZERO ) THEN
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
                  IF ( A(j,l)/=DCMPLX(ZERO) ) THEN
                     temp = Alpha*DCONJG(A(j,l))
                     DO i = 1 , j - 1
                        C(i,j) = C(i,j) + temp*A(i,l)
                     ENDDO
                     C(j,j) = DBLE(C(j,j)) + DBLE(temp*A(i,l))
                  ENDIF
               ENDDO
            ENDDO
         ELSE
            DO j = 1 , N
               IF ( Beta==ZERO ) THEN
                  DO i = j , N
                     C(i,j) = ZERO
                  ENDDO
               ELSEIF ( Beta/=ONE ) THEN
                  C(j,j) = Beta*DBLE(C(j,j))
                  DO i = j + 1 , N
                     C(i,j) = Beta*C(i,j)
                  ENDDO
               ELSE
                  C(j,j) = DBLE(C(j,j))
               ENDIF
               DO l = 1 , K
                  IF ( A(j,l)/=DCMPLX(ZERO) ) THEN
                     temp = Alpha*DCONJG(A(j,l))
                     C(j,j) = DBLE(C(j,j)) + DBLE(temp*A(j,l))
                     DO i = j + 1 , N
                        C(i,j) = C(i,j) + temp*A(i,l)
                     ENDDO
                  ENDIF
               ENDDO
            ENDDO
         ENDIF
!
!        Form  C := alpha*A**H*A + beta*C.
!
      ELSEIF ( upper ) THEN
         DO j = 1 , N
            DO i = 1 , j - 1
               temp = ZERO
               DO l = 1 , K
                  temp = temp + DCONJG(A(l,i))*A(l,j)
               ENDDO
               IF ( Beta==ZERO ) THEN
                  C(i,j) = Alpha*temp
               ELSE
                  C(i,j) = Alpha*temp + Beta*C(i,j)
               ENDIF
            ENDDO
            rtemp = ZERO
            DO l = 1 , K
               rtemp = rtemp + DCONJG(A(l,j))*A(l,j)
            ENDDO
            IF ( Beta==ZERO ) THEN
               C(j,j) = Alpha*rtemp
            ELSE
               C(j,j) = Alpha*rtemp + Beta*DBLE(C(j,j))
            ENDIF
         ENDDO
      ELSE
         DO j = 1 , N
            rtemp = ZERO
            DO l = 1 , K
               rtemp = rtemp + DCONJG(A(l,j))*A(l,j)
            ENDDO
            IF ( Beta==ZERO ) THEN
               C(j,j) = Alpha*rtemp
            ELSE
               C(j,j) = Alpha*rtemp + Beta*DBLE(C(j,j))
            ENDIF
            DO i = j + 1 , N
               temp = ZERO
               DO l = 1 , K
                  temp = temp + DCONJG(A(l,i))*A(l,j)
               ENDDO
               IF ( Beta==ZERO ) THEN
                  C(i,j) = Alpha*temp
               ELSE
                  C(i,j) = Alpha*temp + Beta*C(i,j)
               ENDIF
            ENDDO
         ENDDO
      ENDIF
!
!
!     End of ZHERK .
!
      END SUBROUTINE ZHERK
