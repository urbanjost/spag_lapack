!*==slatsqr.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b SLATSQR
!
!  Definition:
!  ===========
!
!       SUBROUTINE SLATSQR( M, N, MB, NB, A, LDA, T, LDT, WORK,
!                           LWORK, INFO)
!
!       .. Scalar Arguments ..
!       INTEGER           INFO, LDA, M, N, MB, NB, LDT, LWORK
!       ..
!       .. Array Arguments ..
!       REAL              A( LDA, * ), T( LDT, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SLATSQR computes a blocked Tall-Skinny QR factorization of
!> a real M-by-N matrix A for M >= N:
!>
!>    A = Q * ( R ),
!>            ( 0 )
!>
!> where:
!>
!>    Q is a M-by-M orthogonal matrix, stored on exit in an implicit
!>    form in the elements below the diagonal of the array A and in
!>    the elements of the array T;
!>
!>    R is an upper-triangular N-by-N matrix, stored on exit in
!>    the elements on and above the diagonal of the array A.
!>
!>    0 is a (M-N)-by-N zero matrix, and is not stored.
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
!>          The number of columns of the matrix A. M >= N >= 0.
!> \endverbatim
!>
!> \param[in] MB
!> \verbatim
!>          MB is INTEGER
!>          The row block size to be used in the blocked QR.
!>          MB > N.
!> \endverbatim
!>
!> \param[in] NB
!> \verbatim
!>          NB is INTEGER
!>          The column block size to be used in the blocked QR.
!>          N >= NB >= 1.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is REAL array, dimension (LDA,N)
!>          On entry, the M-by-N matrix A.
!>          On exit, the elements on and above the diagonal
!>          of the array contain the N-by-N upper triangular matrix R;
!>          the elements below the diagonal represent Q by the columns
!>          of blocked V (see Further Details).
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
!>          T is REAL array,
!>          dimension (LDT, N * Number_of_row_blocks)
!>          where Number_of_row_blocks = CEIL((M-N)/(MB-N))
!>          The blocked upper triangular block reflectors stored in compact form
!>          as a sequence of upper triangular blocks.
!>          See Further Details below.
!> \endverbatim
!>
!> \param[in] LDT
!> \verbatim
!>          LDT is INTEGER
!>          The leading dimension of the array T.  LDT >= NB.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>         (workspace) REAL array, dimension (MAX(1,LWORK))
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          The dimension of the array WORK.  LWORK >= NB*N.
!>          If LWORK = -1, then a workspace query is assumed; the routine
!>          only calculates the optimal size of the WORK array, returns
!>          this value as the first entry of the WORK array, and no error
!>          message related to LWORK is issued by XERBLA.
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
!> \par Further Details:
!  =====================
!>
!> \verbatim
!> Tall-Skinny QR (TSQR) performs QR by a sequence of orthogonal transformations,
!> representing Q as a product of other orthogonal matrices
!>   Q = Q(1) * Q(2) * . . . * Q(k)
!> where each Q(i) zeros out subdiagonal entries of a block of MB rows of A:
!>   Q(1) zeros out the subdiagonal entries of rows 1:MB of A
!>   Q(2) zeros out the bottom MB-N rows of rows [1:N,MB+1:2*MB-N] of A
!>   Q(3) zeros out the bottom MB-N rows of rows [1:N,2*MB-N+1:3*MB-2*N] of A
!>   . . .
!>
!> Q(1) is computed by GEQRT, which represents Q(1) by Householder vectors
!> stored under the diagonal of rows 1:MB of A, and by upper triangular
!> block reflectors, stored in array T(1:LDT,1:N).
!> For more information see Further Details in GEQRT.
!>
!> Q(i) for i>1 is computed by TPQRT, which represents Q(i) by Householder vectors
!> stored in rows [(i-1)*(MB-N)+N+1:i*(MB-N)+N] of A, and by upper triangular
!> block reflectors, stored in array T(1:LDT,(i-1)*N+1:i*N).
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
      SUBROUTINE SLATSQR(M,N,Mb,Nb,A,Lda,T,Ldt,Work,Lwork,Info)
      USE S_LSAME
      USE S_SGEQRT
      USE S_STPQRT
      USE S_XERBLA
      IMPLICIT NONE
!*--SLATSQR171
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER :: M
      INTEGER :: N
      INTEGER :: Mb
      INTEGER :: Nb
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      REAL , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      INTEGER :: ctr , i , ii , kk
      LOGICAL :: lquery
!
! End of declarations rewritten by SPAG
!
!     ..
!     .. Array Arguments ..
!     ..
!
!  =====================================================================
!
!     ..
!     .. Local Scalars ..
!     ..
!     .. EXTERNAL FUNCTIONS ..
!     .. EXTERNAL SUBROUTINES ..
!     .. INTRINSIC FUNCTIONS ..
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
      ELSEIF ( N<0 .OR. M<N ) THEN
         Info = -2
      ELSEIF ( Mb<=N ) THEN
         Info = -3
      ELSEIF ( Nb<1 .OR. (Nb>N .AND. N>0) ) THEN
         Info = -4
      ELSEIF ( Lda<MAX(1,M) ) THEN
         Info = -5
      ELSEIF ( Ldt<Nb ) THEN
         Info = -8
      ELSEIF ( Lwork<(N*Nb) .AND. (.NOT.lquery) ) THEN
         Info = -10
      ENDIF
      IF ( Info==0 ) Work(1) = Nb*N
      IF ( Info/=0 ) THEN
         CALL XERBLA('SLATSQR',-Info)
         RETURN
      ELSEIF ( lquery ) THEN
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( MIN(M,N)==0 ) RETURN
!
!     The QR Decomposition
!
      IF ( (Mb<=N) .OR. (Mb>=M) ) THEN
         CALL SGEQRT(M,N,Nb,A,Lda,T,Ldt,Work,Info)
         RETURN
      ENDIF
      kk = MOD((M-N),(Mb-N))
      ii = M - kk + 1
!
!      Compute the QR factorization of the first block A(1:MB,1:N)
!
      CALL SGEQRT(Mb,N,Nb,A(1,1),Lda,T,Ldt,Work,Info)
!
      ctr = 1
      DO i = Mb + 1 , ii - Mb + N , (Mb-N)
!
!      Compute the QR factorization of the current block A(I:I+MB-N,1:N)
!
         CALL STPQRT(Mb-N,N,0,Nb,A(1,1),Lda,A(i,1),Lda,T(1,ctr*N+1),Ldt,&
     &               Work,Info)
         ctr = ctr + 1
      ENDDO
!
!      Compute the QR factorization of the last block A(II:M,1:N)
!
      IF ( ii<=M ) CALL STPQRT(kk,N,0,Nb,A(1,1),Lda,A(ii,1),Lda,        &
     &                         T(1,ctr*N+1),Ldt,Work,Info)
!
      Work(1) = N*Nb
!
!     End of SLATSQR
!
      END SUBROUTINE SLATSQR
