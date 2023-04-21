!*==stpqrt2.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b STPQRT2 computes a QR factorization of a real or complex "triangular-pentagonal" matrix, which is composed of a triangular block and a pentagonal block, using the compact WY representation for Q.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download STPQRT2 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/stpqrt2.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/stpqrt2.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/stpqrt2.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE STPQRT2( M, N, L, A, LDA, B, LDB, T, LDT, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER   INFO, LDA, LDB, LDT, N, M, L
!       ..
!       .. Array Arguments ..
!       REAL   A( LDA, * ), B( LDB, * ), T( LDT, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> STPQRT2 computes a QR factorization of a real "triangular-pentagonal"
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
!>          The number of rows of the upper trapezoidal part of B.
!>          MIN(M,N) >= L >= 0.  See Further Details.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is REAL array, dimension (LDA,N)
!>          On entry, the upper triangular N-by-N matrix A.
!>          On exit, the elements on and above the diagonal of the array
!>          contain the upper triangular matrix R.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,N).
!> \endverbatim
!>
!> \param[in,out] B
!> \verbatim
!>          B is REAL array, dimension (LDB,N)
!>          On entry, the pentagonal M-by-N matrix B.  The first M-L rows
!>          are rectangular, and the last L rows are upper trapezoidal.
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
!>          T is REAL array, dimension (LDT,N)
!>          The N-by-N upper triangular factor T of the block reflector.
!>          See Further Details.
!> \endverbatim
!>
!> \param[in] LDT
!> \verbatim
!>          LDT is INTEGER
!>          The leading dimension of the array T.  LDT >= max(1,N)
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
!> \date December 2016
!
!> \ingroup realOTHERcomputational
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  The input matrix C is a (N+M)-by-N matrix
!>
!>               C = [ A ]
!>                   [ B ]
!>
!>  where A is an upper triangular N-by-N matrix, and B is M-by-N pentagonal
!>  matrix consisting of a (M-L)-by-N rectangular matrix B1 on top of a L-by-N
!>  upper trapezoidal matrix B2:
!>
!>               B = [ B1 ]  <- (M-L)-by-N rectangular
!>                   [ B2 ]  <-     L-by-N upper trapezoidal.
!>
!>  The upper trapezoidal matrix B2 consists of the first L rows of a
!>  N-by-N upper triangular matrix, where 0 <= L <= MIN(M,N).  If L=0,
!>  B is rectangular M-by-N; if M=L=N, B is upper triangular.
!>
!>  The matrix W stores the elementary reflectors H(i) in the i-th column
!>  below the diagonal (of A) in the (N+M)-by-N input matrix C
!>
!>               C = [ A ]  <- upper triangular N-by-N
!>                   [ B ]  <- M-by-N pentagonal
!>
!>  so that W can be represented as
!>
!>               W = [ I ]  <- identity, N-by-N
!>                   [ V ]  <- M-by-N, same form as B.
!>
!>  Thus, all of information needed for W is contained on exit in B, which
!>  we call V above.  Note that V has the same form as B; that is,
!>
!>               V = [ V1 ] <- (M-L)-by-N rectangular
!>                   [ V2 ] <-     L-by-N upper trapezoidal.
!>
!>  The columns of V represent the vectors which define the H(i)'s.
!>  The (M+N)-by-(M+N) block reflector H is then given by
!>
!>               H = I - W * T * W^H
!>
!>  where W^H is the conjugate transpose of W and T is the upper triangular
!>  factor of the block reflector.
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE STPQRT2(M,N,L,A,Lda,B,Ldb,T,Ldt,Info)
      IMPLICIT NONE
!*--STPQRT2177
!
!  -- LAPACK computational routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Info , Lda , Ldb , Ldt , N , M , L
!     ..
!     .. Array Arguments ..
      REAL A(Lda,*) , B(Ldb,*) , T(Ldt,*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ONE , ZERO
      PARAMETER (ONE=1.0,ZERO=0.0)
!     ..
!     .. Local Scalars ..
      INTEGER i , j , p , mp , np
      REAL alpha
!     ..
!     .. External Subroutines ..
      EXTERNAL SLARFG , SGEMV , SGER , STRMV , XERBLA
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
      ELSEIF ( Lda<MAX(1,N) ) THEN
         Info = -5
      ELSEIF ( Ldb<MAX(1,M) ) THEN
         Info = -7
      ELSEIF ( Ldt<MAX(1,N) ) THEN
         Info = -9
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('STPQRT2',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( N==0 .OR. M==0 ) RETURN
!
      DO i = 1 , N
!
!        Generate elementary reflector H(I) to annihilate B(:,I)
!
         p = M - L + MIN(L,i)
         CALL SLARFG(p+1,A(i,i),B(1,i),1,T(i,1))
         IF ( i<N ) THEN
!
!           W(1:N-I) := C(I:M,I+1:N)^H * C(I:M,I) [use W = T(:,N)]
!
            DO j = 1 , N - i
               T(j,N) = (A(i,i+j))
            ENDDO
            CALL SGEMV('T',p,N-i,ONE,B(1,i+1),Ldb,B(1,i),1,ONE,T(1,N),1)
!
!           C(I:M,I+1:N) = C(I:m,I+1:N) + alpha*C(I:M,I)*W(1:N-1)^H
!
            alpha = -(T(i,1))
            DO j = 1 , N - i
               A(i,i+j) = A(i,i+j) + alpha*(T(j,N))
            ENDDO
            CALL SGER(p,N-i,alpha,B(1,i),1,T(1,N),1,B(1,i+1),Ldb)
         ENDIF
      ENDDO
!
      DO i = 2 , N
!
!        T(1:I-1,I) := C(I:M,1:I-1)^H * (alpha * C(I:M,I))
!
         alpha = -T(i,1)
 
         DO j = 1 , i - 1
            T(j,i) = ZERO
         ENDDO
         p = MIN(i-1,L)
         mp = MIN(M-L+1,M)
         np = MIN(p+1,N)
!
!        Triangular part of B2
!
         DO j = 1 , p
            T(j,i) = alpha*B(M-L+j,i)
         ENDDO
         CALL STRMV('U','T','N',p,B(mp,1),Ldb,T(1,i),1)
!
!        Rectangular part of B2
!
         CALL SGEMV('T',L,i-1-p,alpha,B(mp,np),Ldb,B(mp,i),1,ZERO,      &
     &              T(np,i),1)
!
!        B1
!
         CALL SGEMV('T',M-L,i-1,alpha,B,Ldb,B(1,i),1,ONE,T(1,i),1)
!
!        T(1:I-1,I) := T(1:I-1,1:I-1) * T(1:I-1,I)
!
         CALL STRMV('U','N','N',i-1,T,Ldt,T(1,i),1)
!
!        T(I,I) = tau(I)
!
         T(i,i) = T(i,1)
         T(i,1) = ZERO
      ENDDO
 
!
!     End of STPQRT2
!
      END SUBROUTINE STPQRT2
