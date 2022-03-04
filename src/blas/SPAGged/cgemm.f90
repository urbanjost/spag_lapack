!*==cgemm.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b CGEMM
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE CGEMM(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
!
!       .. Scalar Arguments ..
!       COMPLEX ALPHA,BETA
!       INTEGER K,LDA,LDB,LDC,M,N
!       CHARACTER TRANSA,TRANSB
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
!> CGEMM  performs one of the matrix-matrix operations
!>
!>    C := alpha*op( A )*op( B ) + beta*C,
!>
!> where  op( X ) is one of
!>
!>    op( X ) = X   or   op( X ) = X**T   or   op( X ) = X**H,
!>
!> alpha and beta are scalars, and A, B and C are matrices, with op( A )
!> an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] TRANSA
!> \verbatim
!>          TRANSA is CHARACTER*1
!>           On entry, TRANSA specifies the form of op( A ) to be used in
!>           the matrix multiplication as follows:
!>
!>              TRANSA = 'N' or 'n',  op( A ) = A.
!>
!>              TRANSA = 'T' or 't',  op( A ) = A**T.
!>
!>              TRANSA = 'C' or 'c',  op( A ) = A**H.
!> \endverbatim
!>
!> \param[in] TRANSB
!> \verbatim
!>          TRANSB is CHARACTER*1
!>           On entry, TRANSB specifies the form of op( B ) to be used in
!>           the matrix multiplication as follows:
!>
!>              TRANSB = 'N' or 'n',  op( B ) = B.
!>
!>              TRANSB = 'T' or 't',  op( B ) = B**T.
!>
!>              TRANSB = 'C' or 'c',  op( B ) = B**H.
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>           On entry,  M  specifies  the number  of rows  of the  matrix
!>           op( A )  and of the  matrix  C.  M  must  be at least  zero.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>           On entry,  N  specifies the number  of columns of the matrix
!>           op( B ) and the number of columns of the matrix C. N must be
!>           at least zero.
!> \endverbatim
!>
!> \param[in] K
!> \verbatim
!>          K is INTEGER
!>           On entry,  K  specifies  the number of columns of the matrix
!>           op( A ) and the number of rows of the matrix op( B ). K must
!>           be at least  zero.
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
!>           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
!>           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
!>           part of the array  A  must contain the matrix  A,  otherwise
!>           the leading  k by m  part of the array  A  must contain  the
!>           matrix A.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>           On entry, LDA specifies the first dimension of A as declared
!>           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
!>           LDA must be at least  max( 1, m ), otherwise  LDA must be at
!>           least  max( 1, k ).
!> \endverbatim
!>
!> \param[in] B
!> \verbatim
!>          B is COMPLEX array, dimension ( LDB, kb ), where kb is
!>           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
!>           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
!>           part of the array  B  must contain the matrix  B,  otherwise
!>           the leading  n by k  part of the array  B  must contain  the
!>           matrix B.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>           On entry, LDB specifies the first dimension of B as declared
!>           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
!>           LDB must be at least  max( 1, k ), otherwise  LDB must be at
!>           least  max( 1, n ).
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
!>           On exit, the array  C  is overwritten by the  m by n  matrix
!>           ( alpha*op( A )*op( B ) + beta*C ).
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
      SUBROUTINE CGEMM(Transa,Transb,M,N,K,Alpha,A,Lda,B,Ldb,Beta,C,Ldc)
      USE S_LSAME
      USE S_XERBLA
      IMPLICIT NONE
!*--CGEMM193
!
! PARAMETER definitions rewritten by SPAG
!
      COMPLEX , PARAMETER  ::  ONE = (1.0E+0,0.0E+0) ,                  &
     &                         ZERO = (0.0E+0,0.0E+0)
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Transa
      CHARACTER :: Transb
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: K
      COMPLEX , INTENT(IN) :: Alpha
      COMPLEX , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      COMPLEX , INTENT(IN) , DIMENSION(Ldb,*) :: B
      INTEGER , INTENT(IN) :: Ldb
      COMPLEX , INTENT(IN) :: Beta
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldc,*) :: C
      INTEGER , INTENT(IN) :: Ldc
!
! Local variable declarations rewritten by SPAG
!
      LOGICAL :: conja , conjb , nota , notb
      INTEGER :: i , info , j , l , nrowa , nrowb
      COMPLEX :: temp
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
!     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
!     conjugated or transposed, set  CONJA and CONJB  as true if  A  and
!     B  respectively are to be  transposed but  not conjugated  and set
!     NROWA and  NROWB  as the number of rows of  A  and  B  respectively.
!
      nota = LSAME(Transa,'N')
      notb = LSAME(Transb,'N')
      conja = LSAME(Transa,'C')
      conjb = LSAME(Transb,'C')
      IF ( nota ) THEN
         nrowa = M
      ELSE
         nrowa = K
      ENDIF
      IF ( notb ) THEN
         nrowb = K
      ELSE
         nrowb = N
      ENDIF
!
!     Test the input parameters.
!
      info = 0
      IF ( (.NOT.nota) .AND. (.NOT.conja) .AND. (.NOT.LSAME(Transa,'T'))&
     &     ) THEN
         info = 1
      ELSEIF ( (.NOT.notb) .AND. (.NOT.conjb) .AND.                     &
     &         (.NOT.LSAME(Transb,'T')) ) THEN
         info = 2
      ELSEIF ( M<0 ) THEN
         info = 3
      ELSEIF ( N<0 ) THEN
         info = 4
      ELSEIF ( K<0 ) THEN
         info = 5
      ELSEIF ( Lda<MAX(1,nrowa) ) THEN
         info = 8
      ELSEIF ( Ldb<MAX(1,nrowb) ) THEN
         info = 10
      ELSEIF ( Ldc<MAX(1,M) ) THEN
         info = 13
      ENDIF
      IF ( info/=0 ) THEN
         CALL XERBLA('CGEMM ',info)
         RETURN
      ENDIF
!
!     Quick return if possible.
!
      IF ( (M==0) .OR. (N==0) .OR.                                      &
     &     (((Alpha==ZERO) .OR. (K==0)) .AND. (Beta==ONE)) ) RETURN
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
      IF ( notb ) THEN
         IF ( nota ) THEN
!
!           Form  C := alpha*A*B + beta*C.
!
            DO j = 1 , N
               IF ( Beta==ZERO ) THEN
                  DO i = 1 , M
                     C(i,j) = ZERO
                  ENDDO
               ELSEIF ( Beta/=ONE ) THEN
                  DO i = 1 , M
                     C(i,j) = Beta*C(i,j)
                  ENDDO
               ENDIF
               DO l = 1 , K
                  temp = Alpha*B(l,j)
                  DO i = 1 , M
                     C(i,j) = C(i,j) + temp*A(i,l)
                  ENDDO
               ENDDO
            ENDDO
         ELSEIF ( conja ) THEN
!
!           Form  C := alpha*A**H*B + beta*C.
!
            DO j = 1 , N
               DO i = 1 , M
                  temp = ZERO
                  DO l = 1 , K
                     temp = temp + CONJG(A(l,i))*B(l,j)
                  ENDDO
                  IF ( Beta==ZERO ) THEN
                     C(i,j) = Alpha*temp
                  ELSE
                     C(i,j) = Alpha*temp + Beta*C(i,j)
                  ENDIF
               ENDDO
            ENDDO
         ELSE
!
!           Form  C := alpha*A**T*B + beta*C
!
            DO j = 1 , N
               DO i = 1 , M
                  temp = ZERO
                  DO l = 1 , K
                     temp = temp + A(l,i)*B(l,j)
                  ENDDO
                  IF ( Beta==ZERO ) THEN
                     C(i,j) = Alpha*temp
                  ELSE
                     C(i,j) = Alpha*temp + Beta*C(i,j)
                  ENDIF
               ENDDO
            ENDDO
         ENDIF
      ELSEIF ( nota ) THEN
         IF ( conjb ) THEN
!
!           Form  C := alpha*A*B**H + beta*C.
!
            DO j = 1 , N
               IF ( Beta==ZERO ) THEN
                  DO i = 1 , M
                     C(i,j) = ZERO
                  ENDDO
               ELSEIF ( Beta/=ONE ) THEN
                  DO i = 1 , M
                     C(i,j) = Beta*C(i,j)
                  ENDDO
               ENDIF
               DO l = 1 , K
                  temp = Alpha*CONJG(B(j,l))
                  DO i = 1 , M
                     C(i,j) = C(i,j) + temp*A(i,l)
                  ENDDO
               ENDDO
            ENDDO
         ELSE
!
!           Form  C := alpha*A*B**T + beta*C
!
            DO j = 1 , N
               IF ( Beta==ZERO ) THEN
                  DO i = 1 , M
                     C(i,j) = ZERO
                  ENDDO
               ELSEIF ( Beta/=ONE ) THEN
                  DO i = 1 , M
                     C(i,j) = Beta*C(i,j)
                  ENDDO
               ENDIF
               DO l = 1 , K
                  temp = Alpha*B(j,l)
                  DO i = 1 , M
                     C(i,j) = C(i,j) + temp*A(i,l)
                  ENDDO
               ENDDO
            ENDDO
         ENDIF
      ELSEIF ( conja ) THEN
         IF ( conjb ) THEN
!
!           Form  C := alpha*A**H*B**H + beta*C.
!
            DO j = 1 , N
               DO i = 1 , M
                  temp = ZERO
                  DO l = 1 , K
                     temp = temp + CONJG(A(l,i))*CONJG(B(j,l))
                  ENDDO
                  IF ( Beta==ZERO ) THEN
                     C(i,j) = Alpha*temp
                  ELSE
                     C(i,j) = Alpha*temp + Beta*C(i,j)
                  ENDIF
               ENDDO
            ENDDO
         ELSE
!
!           Form  C := alpha*A**H*B**T + beta*C
!
            DO j = 1 , N
               DO i = 1 , M
                  temp = ZERO
                  DO l = 1 , K
                     temp = temp + CONJG(A(l,i))*B(j,l)
                  ENDDO
                  IF ( Beta==ZERO ) THEN
                     C(i,j) = Alpha*temp
                  ELSE
                     C(i,j) = Alpha*temp + Beta*C(i,j)
                  ENDIF
               ENDDO
            ENDDO
         ENDIF
      ELSEIF ( conjb ) THEN
!
!           Form  C := alpha*A**T*B**H + beta*C
!
         DO j = 1 , N
            DO i = 1 , M
               temp = ZERO
               DO l = 1 , K
                  temp = temp + A(l,i)*CONJG(B(j,l))
               ENDDO
               IF ( Beta==ZERO ) THEN
                  C(i,j) = Alpha*temp
               ELSE
                  C(i,j) = Alpha*temp + Beta*C(i,j)
               ENDIF
            ENDDO
         ENDDO
      ELSE
!
!           Form  C := alpha*A**T*B**T + beta*C
!
         DO j = 1 , N
            DO i = 1 , M
               temp = ZERO
               DO l = 1 , K
                  temp = temp + A(l,i)*B(j,l)
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
!     End of CGEMM .
!
      END SUBROUTINE CGEMM
