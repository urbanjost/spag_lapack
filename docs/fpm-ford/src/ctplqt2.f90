!*==ctplqt2.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b CTPLQT2
!
!  Definition:
!  ===========
!
!       SUBROUTINE CTPLQT2( M, N, L, A, LDA, B, LDB, T, LDT, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER   INFO, LDA, LDB, LDT, N, M, L
!       ..
!       .. Array Arguments ..
!       COMPLEX   A( LDA, * ), B( LDB, * ), T( LDT, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CTPLQT2 computes a LQ a factorization of a complex "triangular-pentagonal"
!> matrix C, which is composed of a triangular block A and pentagonal block B,
!> using the compact WY representation for Q.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The total number of rows of the matrix B.
!>          M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix B, and the order of
!>          the triangular matrix A.
!>          N >= 0.
!> \endverbatim
!>
!> \param[in] L
!> \verbatim
!>          L is INTEGER
!>          The number of rows of the lower trapezoidal part of B.
!>          MIN(M,N) >= L >= 0.  See Further Details.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is COMPLEX array, dimension (LDA,M)
!>          On entry, the lower triangular M-by-M matrix A.
!>          On exit, the elements on and below the diagonal of the array
!>          contain the lower triangular matrix L.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,M).
!> \endverbatim
!>
!> \param[in,out] B
!> \verbatim
!>          B is COMPLEX array, dimension (LDB,N)
!>          On entry, the pentagonal M-by-N matrix B.  The first N-L columns
!>          are rectangular, and the last L columns are lower trapezoidal.
!>          On exit, B contains the pentagonal matrix V.  See Further Details.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of the array B.  LDB >= max(1,M).
!> \endverbatim
!>
!> \param[out] T
!> \verbatim
!>          T is COMPLEX array, dimension (LDT,M)
!>          The N-by-N upper triangular factor T of the block reflector.
!>          See Further Details.
!> \endverbatim
!>
!> \param[in] LDT
!> \verbatim
!>          LDT is INTEGER
!>          The leading dimension of the array T.  LDT >= max(1,M)
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0: successful exit
!>          < 0: if INFO = -i, the i-th argument had an illegal value
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
!> \date June 2017
!
!> \ingroup doubleOTHERcomputational
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  The input matrix C is a M-by-(M+N) matrix
!>
!>               C = [ A ][ B ]
!>
!>
!>  where A is an lower triangular M-by-M matrix, and B is M-by-N pentagonal
!>  matrix consisting of a M-by-(N-L) rectangular matrix B1 left of a M-by-L
!>  upper trapezoidal matrix B2:
!>
!>               B = [ B1 ][ B2 ]
!>                   [ B1 ]  <-     M-by-(N-L) rectangular
!>                   [ B2 ]  <-     M-by-L lower trapezoidal.
!>
!>  The lower trapezoidal matrix B2 consists of the first L columns of a
!>  N-by-N lower triangular matrix, where 0 <= L <= MIN(M,N).  If L=0,
!>  B is rectangular M-by-N; if M=L=N, B is lower triangular.
!>
!>  The matrix W stores the elementary reflectors H(i) in the i-th row
!>  above the diagonal (of A) in the M-by-(M+N) input matrix C
!>
!>               C = [ A ][ B ]
!>                   [ A ]  <- lower triangular M-by-M
!>                   [ B ]  <- M-by-N pentagonal
!>
!>  so that W can be represented as
!>
!>               W = [ I ][ V ]
!>                   [ I ]  <- identity, M-by-M
!>                   [ V ]  <- M-by-N, same form as B.
!>
!>  Thus, all of information needed for W is contained on exit in B, which
!>  we call V above.  Note that V has the same form as B; that is,
!>
!>               W = [ V1 ][ V2 ]
!>                   [ V1 ] <-     M-by-(N-L) rectangular
!>                   [ V2 ] <-     M-by-L lower trapezoidal.
!>
!>  The rows of V represent the vectors which define the H(i)'s.
!>  The (M+N)-by-(M+N) block reflector H is then given by
!>
!>               H = I - W**T * T * W
!>
!>  where W^H is the conjugate transpose of W and T is the upper triangular
!>  factor of the block reflector.
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE CTPLQT2(M,N,L,A,Lda,B,Ldb,T,Ldt,Info)
      IMPLICIT NONE
!*--CTPLQT2166
!
!  -- LAPACK computational routine (version 3.7.1) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     June 2017
!
!     .. Scalar Arguments ..
      INTEGER Info , Lda , Ldb , Ldt , N , M , L
!     ..
!     .. Array Arguments ..
      COMPLEX A(Lda,*) , B(Ldb,*) , T(Ldt,*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      COMPLEX ONE , ZERO
      PARAMETER (ZERO=(0.0E+0,0.0E+0),ONE=(1.0E+0,0.0E+0))
!     ..
!     .. Local Scalars ..
      INTEGER i , j , p , mp , np
      COMPLEX alpha
!     ..
!     .. External Subroutines ..
      EXTERNAL CLARFG , CGEMV , CGERC , CTRMV , XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX , MIN
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      Info = 0
      IF ( M<0 ) THEN
         Info = -1
      ELSEIF ( N<0 ) THEN
         Info = -2
      ELSEIF ( L<0 .OR. L>MIN(M,N) ) THEN
         Info = -3
      ELSEIF ( Lda<MAX(1,M) ) THEN
         Info = -5
      ELSEIF ( Ldb<MAX(1,M) ) THEN
         Info = -7
      ELSEIF ( Ldt<MAX(1,M) ) THEN
         Info = -9
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('CTPLQT2',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( N==0 .OR. M==0 ) RETURN
!
      DO i = 1 , M
!
!        Generate elementary reflector H(I) to annihilate B(I,:)
!
         p = N - L + MIN(L,i)
         CALL CLARFG(p+1,A(i,i),B(i,1),Ldb,T(1,i))
         T(1,i) = CONJG(T(1,i))
         IF ( i<M ) THEN
            DO j = 1 , p
               B(i,j) = CONJG(B(i,j))
            ENDDO
!
!           W(M-I:1) := C(I+1:M,I:N) * C(I,I:N) [use W = T(M,:)]
!
            DO j = 1 , M - i
               T(M,j) = (A(i+j,i))
            ENDDO
            CALL CGEMV('N',M-i,p,ONE,B(i+1,1),Ldb,B(i,1),Ldb,ONE,T(M,1),&
     &                 Ldt)
!
!           C(I+1:M,I:N) = C(I+1:M,I:N) + alpha * C(I,I:N)*W(M-1:1)^H
!
            alpha = -(T(1,i))
            DO j = 1 , M - i
               A(i+j,i) = A(i+j,i) + alpha*(T(M,j))
            ENDDO
            CALL CGERC(M-i,p,(alpha),T(M,1),Ldt,B(i,1),Ldb,B(i+1,1),Ldb)
            DO j = 1 , p
               B(i,j) = CONJG(B(i,j))
            ENDDO
         ENDIF
      ENDDO
!
      DO i = 2 , M
!
!        T(I,1:I-1) := C(I:I-1,1:N)**H * (alpha * C(I,I:N))
!
         alpha = -(T(1,i))
         DO j = 1 , i - 1
            T(i,j) = ZERO
         ENDDO
         p = MIN(i-1,L)
         np = MIN(N-L+1,N)
         mp = MIN(p+1,M)
         DO j = 1 , N - L + p
            B(i,j) = CONJG(B(i,j))
         ENDDO
!
!        Triangular part of B2
!
         DO j = 1 , p
            T(i,j) = (alpha*B(i,N-L+j))
         ENDDO
         CALL CTRMV('L','N','N',p,B(1,np),Ldb,T(i,1),Ldt)
!
!        Rectangular part of B2
!
         CALL CGEMV('N',i-1-p,L,alpha,B(mp,np),Ldb,B(i,np),Ldb,ZERO,    &
     &              T(i,mp),Ldt)
!
!        B1
 
!
         CALL CGEMV('N',i-1,N-L,alpha,B,Ldb,B(i,1),Ldb,ONE,T(i,1),Ldt)
!
 
!
!        T(1:I-1,I) := T(1:I-1,1:I-1) * T(I,1:I-1)
!
         DO j = 1 , i - 1
            T(i,j) = CONJG(T(i,j))
         ENDDO
         CALL CTRMV('L','C','N',i-1,T,Ldt,T(i,1),Ldt)
         DO j = 1 , i - 1
            T(i,j) = CONJG(T(i,j))
         ENDDO
         DO j = 1 , N - L + p
            B(i,j) = CONJG(B(i,j))
         ENDDO
!
!        T(I,I) = tau(I)
!
         T(i,i) = T(1,i)
         T(1,i) = ZERO
      ENDDO
      DO i = 1 , M
         DO j = i + 1 , M
            T(i,j) = (T(j,i))
            T(j,i) = ZERO
         ENDDO
      ENDDO
 
!
!     End of CTPLQT2
!
      END SUBROUTINE CTPLQT2
