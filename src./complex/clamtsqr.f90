!*==clamtsqr.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b CLAMTSQR
!
!  Definition:
!  ===========
!
!      SUBROUTINE CLAMTSQR( SIDE, TRANS, M, N, K, MB, NB, A, LDA, T,
!     $                     LDT, C, LDC, WORK, LWORK, INFO )
!
!
!     .. Scalar Arguments ..
!      CHARACTER         SIDE, TRANS
!      INTEGER           INFO, LDA, M, N, K, MB, NB, LDT, LWORK, LDC
!     ..
!     .. Array Arguments ..
!      COMPLEX        A( LDA, * ), WORK( * ), C(LDC, * ),
!     $                  T( LDT, * )
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>      CLAMTSQR overwrites the general complex M-by-N matrix C with
!>
!>
!>                 SIDE = 'L'     SIDE = 'R'
!> TRANS = 'N':      Q * C          C * Q
!> TRANS = 'C':      Q**H * C       C * Q**H
!>      where Q is a real orthogonal matrix defined as the product
!>      of blocked elementary reflectors computed by tall skinny
!>      QR factorization (CLATSQR)
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] SIDE
!> \verbatim
!>          SIDE is CHARACTER*1
!>          = 'L': apply Q or Q**H from the Left;
!>          = 'R': apply Q or Q**H from the Right.
!> \endverbatim
!>
!> \param[in] TRANS
!> \verbatim
!>          TRANS is CHARACTER*1
!>          = 'N':  No transpose, apply Q;
!>          = 'C':  Conjugate Transpose, apply Q**H.
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix A.  M >=0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix C. M >= N >= 0.
!> \endverbatim
!>
!> \param[in] K
!> \verbatim
!>          K is INTEGER
!>          The number of elementary reflectors whose product defines
!>          the matrix Q.
!>          N >= K >= 0;
!>
!> \endverbatim
!>
!> \param[in] MB
!> \verbatim
!>          MB is INTEGER
!>          The block size to be used in the blocked QR.
!>          MB > N. (must be the same as DLATSQR)
!> \endverbatim
!>
!> \param[in] NB
!> \verbatim
!>          NB is INTEGER
!>          The column block size to be used in the blocked QR.
!>          N >= NB >= 1.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is COMPLEX array, dimension (LDA,K)
!>          The i-th column must contain the vector which defines the
!>          blockedelementary reflector H(i), for i = 1,2,...,k, as
!>          returned by DLATSQR in the first k columns of
!>          its array argument A.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.
!>          If SIDE = 'L', LDA >= max(1,M);
!>          if SIDE = 'R', LDA >= max(1,N).
!> \endverbatim
!>
!> \param[in] T
!> \verbatim
!>          T is COMPLEX array, dimension
!>          ( N * Number of blocks(CEIL(M-K/MB-K)),
!>          The blocked upper triangular block reflectors stored in compact form
!>          as a sequence of upper triangular blocks.  See below
!>          for further details.
!> \endverbatim
!>
!> \param[in] LDT
!> \verbatim
!>          LDT is INTEGER
!>          The leading dimension of the array T.  LDT >= NB.
!> \endverbatim
!>
!> \param[in,out] C
!> \verbatim
!>          C is COMPLEX array, dimension (LDC,N)
!>          On entry, the M-by-N matrix C.
!>          On exit, C is overwritten by Q*C or Q**H*C or C*Q**H or C*Q.
!> \endverbatim
!>
!> \param[in] LDC
!> \verbatim
!>          LDC is INTEGER
!>          The leading dimension of the array C. LDC >= max(1,M).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>         (workspace) COMPLEX array, dimension (MAX(1,LWORK))
!>
!> \endverbatim
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The dimension of the array WORK.
!>
!>          If SIDE = 'L', LWORK >= max(1,N)*NB;
!>          if SIDE = 'R', LWORK >= max(1,MB)*NB.
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
      SUBROUTINE CLAMTSQR(Side,Trans,M,N,K,Mb,Nb,A,Lda,T,Ldt,C,Ldc,Work,&
     &                    Lwork,Info)
      USE S_CGEMQRT
      USE S_CTPMQRT
      USE S_LSAME
      USE S_XERBLA
      IMPLICIT NONE
!*--CLAMTSQR204
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Side
      CHARACTER :: Trans
      INTEGER :: M
      INTEGER :: N
      INTEGER :: K
      INTEGER :: Mb
      INTEGER :: Nb
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      COMPLEX , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      INTEGER :: ctr , i , ii , kk , lw
      LOGICAL :: left , lquery , notran , right , tran
!
! End of declarations rewritten by SPAG
!
!     ..
!     .. Array Arguments ..
!     ..
!
! =====================================================================
!
!     ..
!     .. Local Scalars ..
!     ..
!     .. External Functions ..
!     .. External Subroutines ..
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      lquery = Lwork<0
      notran = LSAME(Trans,'N')
      tran = LSAME(Trans,'C')
      left = LSAME(Side,'L')
      right = LSAME(Side,'R')
      IF ( left ) THEN
         lw = N*Nb
      ELSE
         lw = M*Nb
      ENDIF
!
      Info = 0
      IF ( .NOT.left .AND. .NOT.right ) THEN
         Info = -1
      ELSEIF ( .NOT.tran .AND. .NOT.notran ) THEN
         Info = -2
      ELSEIF ( M<0 ) THEN
         Info = -3
      ELSEIF ( N<0 ) THEN
         Info = -4
      ELSEIF ( K<0 ) THEN
         Info = -5
      ELSEIF ( Lda<MAX(1,K) ) THEN
         Info = -9
      ELSEIF ( Ldt<MAX(1,Nb) ) THEN
         Info = -11
      ELSEIF ( Ldc<MAX(1,M) ) THEN
         Info = -13
      ELSEIF ( (Lwork<MAX(1,lw)) .AND. (.NOT.lquery) ) THEN
         Info = -15
      ENDIF
!
!     Determine the block size if it is tall skinny or short and wide
!
      IF ( Info==0 ) Work(1) = lw
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('CLAMTSQR',-Info)
         RETURN
      ELSEIF ( lquery ) THEN
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( MIN(M,N,K)==0 ) RETURN
!
      IF ( (Mb<=K) .OR. (Mb>=MAX(M,N,K)) ) THEN
         CALL CGEMQRT(Side,Trans,M,N,K,Nb,A,Lda,T,Ldt,C,Ldc,Work,Info)
         RETURN
      ENDIF
!
      IF ( left .AND. notran ) THEN
!
!         Multiply Q to the last block of C
!
         kk = MOD((M-K),(Mb-K))
         ctr = (M-K)/(Mb-K)
         IF ( kk>0 ) THEN
            ii = M - kk + 1
            CALL CTPMQRT('L','N',kk,N,K,0,Nb,A(ii,1),Lda,T(1,ctr*K+1),  &
     &                   Ldt,C(1,1),Ldc,C(ii,1),Ldc,Work,Info)
         ELSE
            ii = M + 1
         ENDIF
!
         DO i = ii - (Mb-K) , Mb + 1 , -(Mb-K)
!
!         Multiply Q to the current block of C (I:I+MB,1:N)
!
            ctr = ctr - 1
            CALL CTPMQRT('L','N',Mb-K,N,K,0,Nb,A(i,1),Lda,T(1,ctr*K+1), &
     &                   Ldt,C(1,1),Ldc,C(i,1),Ldc,Work,Info)
 
         ENDDO
!
!         Multiply Q to the first block of C (1:MB,1:N)
!
         CALL CGEMQRT('L','N',Mb,N,K,Nb,A(1,1),Lda,T,Ldt,C(1,1),Ldc,    &
     &                Work,Info)
!
      ELSEIF ( left .AND. tran ) THEN
!
!         Multiply Q to the first block of C
!
         kk = MOD((M-K),(Mb-K))
         ii = M - kk + 1
         ctr = 1
         CALL CGEMQRT('L','C',Mb,N,K,Nb,A(1,1),Lda,T,Ldt,C(1,1),Ldc,    &
     &                Work,Info)
!
         DO i = Mb + 1 , ii - Mb + K , (Mb-K)
!
!         Multiply Q to the current block of C (I:I+MB,1:N)
!
            CALL CTPMQRT('L','C',Mb-K,N,K,0,Nb,A(i,1),Lda,T(1,ctr*K+1), &
     &                   Ldt,C(1,1),Ldc,C(i,1),Ldc,Work,Info)
            ctr = ctr + 1
!
         ENDDO
!
!         Multiply Q to the last block of C
!
!
         IF ( ii<=M ) CALL CTPMQRT('L','C',kk,N,K,0,Nb,A(ii,1),Lda,     &
     &                             T(1,ctr*K+1),Ldt,C(1,1),Ldc,C(ii,1), &
     &                             Ldc,Work,Info)
!
      ELSEIF ( right .AND. tran ) THEN
!
!         Multiply Q to the last block of C
!
         kk = MOD((N-K),(Mb-K))
         ctr = (N-K)/(Mb-K)
         IF ( kk>0 ) THEN
            ii = N - kk + 1
            CALL CTPMQRT('R','C',M,kk,K,0,Nb,A(ii,1),Lda,T(1,ctr*K+1),  &
     &                   Ldt,C(1,1),Ldc,C(1,ii),Ldc,Work,Info)
         ELSE
            ii = N + 1
         ENDIF
!
         DO i = ii - (Mb-K) , Mb + 1 , -(Mb-K)
!
!         Multiply Q to the current block of C (1:M,I:I+MB)
!
            ctr = ctr - 1
            CALL CTPMQRT('R','C',M,Mb-K,K,0,Nb,A(i,1),Lda,T(1,ctr*K+1), &
     &                   Ldt,C(1,1),Ldc,C(1,i),Ldc,Work,Info)
         ENDDO
!
!         Multiply Q to the first block of C (1:M,1:MB)
!
         CALL CGEMQRT('R','C',M,Mb,K,Nb,A(1,1),Lda,T,Ldt,C(1,1),Ldc,    &
     &                Work,Info)
!
      ELSEIF ( right .AND. notran ) THEN
!
!         Multiply Q to the first block of C
!
         kk = MOD((N-K),(Mb-K))
         ii = N - kk + 1
         ctr = 1
         CALL CGEMQRT('R','N',M,Mb,K,Nb,A(1,1),Lda,T,Ldt,C(1,1),Ldc,    &
     &                Work,Info)
!
         DO i = Mb + 1 , ii - Mb + K , (Mb-K)
!
!         Multiply Q to the current block of C (1:M,I:I+MB)
!
            CALL CTPMQRT('R','N',M,Mb-K,K,0,Nb,A(i,1),Lda,T(1,ctr*K+1), &
     &                   Ldt,C(1,1),Ldc,C(1,i),Ldc,Work,Info)
            ctr = ctr + 1
!
         ENDDO
!
!         Multiply Q to the last block of C
!
!
         IF ( ii<=N ) CALL CTPMQRT('R','N',M,kk,K,0,Nb,A(ii,1),Lda,     &
     &                             T(1,ctr*K+1),Ldt,C(1,1),Ldc,C(1,ii), &
     &                             Ldc,Work,Info)
!
      ENDIF
!
      Work(1) = lw
!
!     End of CLAMTSQR
!
      END SUBROUTINE CLAMTSQR
