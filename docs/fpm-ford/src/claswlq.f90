!*==claswlq.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b CLASWLQ
!
!  Definition:
!  ===========
!
!       SUBROUTINE CLASWLQ( M, N, MB, NB, A, LDA, T, LDT, WORK,
!                            LWORK, INFO)
!
!       .. Scalar Arguments ..
!       INTEGER           INFO, LDA, M, N, MB, NB, LDT, LWORK
!       ..
!       .. Array Arguments ..
!       COMPLEX           A( LDA, * ), T( LDT, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CLASWLQ computes a blocked Tall-Skinny LQ factorization of
!> a complex M-by-N matrix A for M <= N:
!>
!>    A = ( L 0 ) *  Q,
!>
!> where:
!>
!>    Q is a n-by-N orthogonal matrix, stored on exit in an implicit
!>    form in the elements above the diagonal of the array A and in
!>    the elements of the array T;
!>    L is a lower-triangular M-by-M matrix stored on exit in
!>    the elements on and below the diagonal of the array A.
!>    0 is a M-by-(N-M) zero matrix, if M < N, and is not stored.
!>
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix A.  M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix A.  N >= M >= 0.
!> \endverbatim
!>
!> \param[in] MB
!> \verbatim
!>          MB is INTEGER
!>          The row block size to be used in the blocked QR.
!>          M >= MB >= 1
!> \endverbatim
!> \param[in] NB
!> \verbatim
!>          NB is INTEGER
!>          The column block size to be used in the blocked QR.
!>          NB > M.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is COMPLEX array, dimension (LDA,N)
!>          On entry, the M-by-N matrix A.
!>          On exit, the elements on and below the diagonal
!>          of the array contain the N-by-N lower triangular matrix L;
!>          the elements above the diagonal represent Q by the rows
!>          of blocked V (see Further Details).
!>
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,M).
!> \endverbatim
!>
!> \param[out] T
!> \verbatim
!>          T is COMPLEX array,
!>          dimension (LDT, N * Number_of_row_blocks)
!>          where Number_of_row_blocks = CEIL((N-M)/(NB-M))
!>          The blocked upper triangular block reflectors stored in compact form
!>          as a sequence of upper triangular blocks.
!>          See Further Details below.
!> \endverbatim
!>
!> \param[in] LDT
!> \verbatim
!>          LDT is INTEGER
!>          The leading dimension of the array T.  LDT >= MB.
!> \endverbatim
!>
!>
!> \param[out] WORK
!> \verbatim
!>         (workspace) COMPLEX array, dimension (MAX(1,LWORK))
!>
!> \endverbatim
!> \param[in] LWORK
!> \verbatim
!>          The dimension of the array WORK.  LWORK >= MB*M.
!>          If LWORK = -1, then a workspace query is assumed; the routine
!>          only calculates the optimal size of the WORK array, returns
!>          this value as the first entry of the WORK array, and no error
!>          message related to LWORK is issued by XERBLA.
!>
!> \endverbatim
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
!> \par Further Details:
!  =====================
!>
!> \verbatim
!> Short-Wide LQ (SWLQ) performs LQ by a sequence of orthogonal transformations,
!> representing Q as a product of other orthogonal matrices
!>   Q = Q(1) * Q(2) * . . . * Q(k)
!> where each Q(i) zeros out upper diagonal entries of a block of NB rows of A:
!>   Q(1) zeros out the upper diagonal entries of rows 1:NB of A
!>   Q(2) zeros out the bottom MB-N rows of rows [1:M,NB+1:2*NB-M] of A
!>   Q(3) zeros out the bottom MB-N rows of rows [1:M,2*NB-M+1:3*NB-2*M] of A
!>   . . .
!>
!> Q(1) is computed by GELQT, which represents Q(1) by Householder vectors
!> stored under the diagonal of rows 1:MB of A, and by upper triangular
!> block reflectors, stored in array T(1:LDT,1:N).
!> For more information see Further Details in GELQT.
!>
!> Q(i) for i>1 is computed by TPLQT, which represents Q(i) by Householder vectors
!> stored in columns [(i-1)*(NB-M)+M+1:i*(NB-M)+M] of A, and by upper triangular
!> block reflectors, stored in array T(1:LDT,(i-1)*M+1:i*M).
!> The last Q(k) may use fewer rows.
!> For more information see Further Details in TPQRT.
!>
!> For more details of the overall algorithm, see the description of
!> Sequential TSQR in Section 2.2 of [1].
!>
!> [1] “Communication-Optimal Parallel and Sequential QR and LU Factorizations,”
!>     J. Demmel, L. Grigori, M. Hoemmen, J. Langou,
!>     SIAM J. Sci. Comput, vol. 34, no. 1, 2012
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE CLASWLQ(M,N,Mb,Nb,A,Lda,T,Ldt,Work,Lwork,Info)
      IMPLICIT NONE
!*--CLASWLQ165
!
!  -- LAPACK computational routine (version 3.9.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd. --
!     June 2017
!
!     .. Scalar Arguments ..
      INTEGER Info , Lda , M , N , Mb , Nb , Lwork , Ldt
!     ..
!     .. Array Arguments ..
      COMPLEX A(Lda,*) , Work(*) , T(Ldt,*)
!     ..
!
!  =====================================================================
!
!     ..
!     .. Local Scalars ..
      LOGICAL lquery
      INTEGER i , ii , kk , ctr
!     ..
!     .. EXTERNAL FUNCTIONS ..
      LOGICAL LSAME
      EXTERNAL LSAME
!     .. EXTERNAL SUBROUTINES ..
      EXTERNAL CGELQT , CTPLQT , XERBLA
!     .. INTRINSIC FUNCTIONS ..
      INTRINSIC MAX , MIN , MOD
!     ..
!     .. EXTERNAL FUNCTIONS ..
      INTEGER ILAENV
      EXTERNAL ILAENV
!     ..
!     .. EXECUTABLE STATEMENTS ..
!
!     TEST THE INPUT ARGUMENTS
!
      Info = 0
!
      lquery = (Lwork==-1)
!
      IF ( M<0 ) THEN
         Info = -1
      ELSEIF ( N<0 .OR. N<M ) THEN
         Info = -2
      ELSEIF ( Mb<1 .OR. (Mb>M .AND. M>0) ) THEN
         Info = -3
      ELSEIF ( Nb<=M ) THEN
         Info = -4
      ELSEIF ( Lda<MAX(1,M) ) THEN
         Info = -5
      ELSEIF ( Ldt<Mb ) THEN
         Info = -8
      ELSEIF ( (Lwork<M*Mb) .AND. (.NOT.lquery) ) THEN
         Info = -10
      ENDIF
      IF ( Info==0 ) Work(1) = Mb*M
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('CLASWLQ',-Info)
         RETURN
      ELSEIF ( lquery ) THEN
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( MIN(M,N)==0 ) RETURN
!
!     The LQ Decomposition
!
      IF ( (M>=N) .OR. (Nb<=M) .OR. (Nb>=N) ) THEN
         CALL CGELQT(M,N,Mb,A,Lda,T,Ldt,Work,Info)
         RETURN
      ENDIF
!
      kk = MOD((N-M),(Nb-M))
      ii = N - kk + 1
!
!      Compute the LQ factorization of the first block A(1:M,1:NB)
!
      CALL CGELQT(M,Nb,Mb,A(1,1),Lda,T,Ldt,Work,Info)
      ctr = 1
!
      DO i = Nb + 1 , ii - Nb + M , (Nb-M)
!
!      Compute the QR factorization of the current block A(1:M,I:I+NB-M)
!
         CALL CTPLQT(M,Nb-M,0,Mb,A(1,1),Lda,A(1,i),Lda,T(1,ctr*M+1),Ldt,&
     &               Work,Info)
         ctr = ctr + 1
      ENDDO
!
!     Compute the QR factorization of the last block A(1:M,II:N)
!
      IF ( ii<=N ) CALL CTPLQT(M,kk,0,Mb,A(1,1),Lda,A(1,ii),Lda,        &
     &                         T(1,ctr*M+1),Ldt,Work,Info)
!
      Work(1) = M*Mb
!
!     End of CLASWLQ
!
      END SUBROUTINE CLASWLQ
