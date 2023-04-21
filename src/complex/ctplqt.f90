!*==ctplqt.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b CTPLQT
!
!  Definition:
!  ===========
!
!       SUBROUTINE CTPLQT( M, N, L, MB, A, LDA, B, LDB, T, LDT, WORK,
!                          INFO )
!
!       .. Scalar Arguments ..
!       INTEGER         INFO, LDA, LDB, LDT, N, M, L, MB
!       ..
!       .. Array Arguments ..
!       COMPLEX      A( LDA, * ), B( LDB, * ), T( LDT, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CTPLQT computes a blocked LQ factorization of a complex
!> "triangular-pentagonal" matrix C, which is composed of a
!> triangular block A and pentagonal block B, using the compact
!> WY representation for Q.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix B, and the order of the
!>          triangular matrix A.
!>          M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix B.
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
!> \param[in] MB
!> \verbatim
!>          MB is INTEGER
!>          The block size to be used in the blocked QR.  M >= MB >= 1.
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
!>          T is COMPLEX array, dimension (LDT,N)
!>          The lower triangular block reflectors stored in compact form
!>          as a sequence of upper triangular blocks.  See Further Details.
!> \endverbatim
!>
!> \param[in] LDT
!> \verbatim
!>          LDT is INTEGER
!>          The leading dimension of the array T.  LDT >= MB.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX array, dimension (MB*M)
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value
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
!>               C = [ A ] [ B ]
!>
!>
!>  where A is an lower triangular M-by-M matrix, and B is M-by-N pentagonal
!>  matrix consisting of a M-by-(N-L) rectangular matrix B1 on left of a M-by-L
!>  upper trapezoidal matrix B2:
!>          [ B ] = [ B1 ] [ B2 ]
!>                   [ B1 ]  <- M-by-(N-L) rectangular
!>                   [ B2 ]  <-     M-by-L lower trapezoidal.
!>
!>  The lower trapezoidal matrix B2 consists of the first L columns of a
!>  M-by-M lower triangular matrix, where 0 <= L <= MIN(M,N).  If L=0,
!>  B is rectangular M-by-N; if M=L=N, B is lower triangular.
!>
!>  The matrix W stores the elementary reflectors H(i) in the i-th row
!>  above the diagonal (of A) in the M-by-(M+N) input matrix C
!>            [ C ] = [ A ] [ B ]
!>                   [ A ]  <- lower triangular M-by-M
!>                   [ B ]  <- M-by-N pentagonal
!>
!>  so that W can be represented as
!>            [ W ] = [ I ] [ V ]
!>                   [ I ]  <- identity, M-by-M
!>                   [ V ]  <- M-by-N, same form as B.
!>
!>  Thus, all of information needed for W is contained on exit in B, which
!>  we call V above.  Note that V has the same form as B; that is,
!>            [ V ] = [ V1 ] [ V2 ]
!>                   [ V1 ] <- M-by-(N-L) rectangular
!>                   [ V2 ] <-     M-by-L lower trapezoidal.
!>
!>  The rows of V represent the vectors which define the H(i)'s.
!>
!>  The number of blocks is B = ceiling(M/MB), where each
!>  block is of order MB except for the last block, which is of order
!>  IB = M - (M-1)*MB.  For each of the B blocks, a upper triangular block
!>  reflector factor is computed: T1, T2, ..., TB.  The MB-by-MB (and IB-by-IB
!>  for the last block) T's are stored in the MB-by-N matrix T as
!>
!>               T = [T1 T2 ... TB].
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE CTPLQT(M,N,L,Mb,A,Lda,B,Ldb,T,Ldt,Work,Info)
      IMPLICIT NONE
!*--CTPLQT177
!
!  -- LAPACK computational routine (version 3.7.1) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     June 2017
!
!     .. Scalar Arguments ..
      INTEGER Info , Lda , Ldb , Ldt , N , M , L , Mb
!     ..
!     .. Array Arguments ..
      COMPLEX A(Lda,*) , B(Ldb,*) , T(Ldt,*) , Work(*)
!     ..
!
! =====================================================================
!
!     ..
!     .. Local Scalars ..
      INTEGER i , ib , lb , nb , iinfo
!     ..
!     .. External Subroutines ..
      EXTERNAL CTPLQT2 , CTPRFB , XERBLA
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
      ELSEIF ( L<0 .OR. (L>MIN(M,N) .AND. MIN(M,N)>=0) ) THEN
         Info = -3
      ELSEIF ( Mb<1 .OR. (Mb>M .AND. M>0) ) THEN
         Info = -4
      ELSEIF ( Lda<MAX(1,M) ) THEN
         Info = -6
      ELSEIF ( Ldb<MAX(1,M) ) THEN
         Info = -8
      ELSEIF ( Ldt<Mb ) THEN
         Info = -10
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('CTPLQT',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( M==0 .OR. N==0 ) RETURN
!
      DO i = 1 , M , Mb
!
!     Compute the QR factorization of the current block
!
         ib = MIN(M-i+1,Mb)
         nb = MIN(N-L+i+ib-1,N)
         IF ( i>=L ) THEN
            lb = 0
         ELSE
            lb = nb - N + L - i + 1
         ENDIF
!
         CALL CTPLQT2(ib,nb,lb,A(i,i),Lda,B(i,1),Ldb,T(1,i),Ldt,iinfo)
!
!     Update by applying H**T to B(I+IB:M,:) from the right
!
         IF ( i+ib<=M ) CALL CTPRFB('R','N','F','R',M-i-ib+1,nb,ib,lb,  &
     &                              B(i,1),Ldb,T(1,i),Ldt,A(i+ib,i),Lda,&
     &                              B(i+ib,1),Ldb,Work,M-i-ib+1)
      ENDDO
!
!     End of CTPLQT
!
      END SUBROUTINE CTPLQT
