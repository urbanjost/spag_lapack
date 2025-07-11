!*==slamtsqr.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b SLAMTSQR
!
!  Definition:
!  ===========
!
!      SUBROUTINE SLAMTSQR( SIDE, TRANS, M, N, K, MB, NB, A, LDA, T,
!     $                     LDT, C, LDC, WORK, LWORK, INFO )
!
!
!     .. Scalar Arguments ..
!      CHARACTER         SIDE, TRANS
!      INTEGER           INFO, LDA, M, N, K, MB, NB, LDT, LWORK, LDC
!     ..
!     .. Array Arguments ..
!      DOUBLE        A( LDA, * ), WORK( * ), C(LDC, * ),
!     $                  T( LDT, * )
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>      SLAMTSQR overwrites the general real M-by-N matrix C with
!>
!>
!>                 SIDE = 'L'     SIDE = 'R'
!> TRANS = 'N':      Q * C          C * Q
!> TRANS = 'T':      Q**T * C       C * Q**T
!>      where Q is a real orthogonal matrix defined as the product
!>      of blocked elementary reflectors computed by tall skinny
!>      QR factorization (DLATSQR)
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] SIDE
!> \verbatim
!>          SIDE is CHARACTER*1
!>          = 'L': apply Q or Q**T from the Left;
!>          = 'R': apply Q or Q**T from the Right.
!> \endverbatim
!>
!> \param[in] TRANS
!> \verbatim
!>          TRANS is CHARACTER*1
!>          = 'N':  No transpose, apply Q;
!>          = 'T':  Transpose, apply Q**T.
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
!>          A is REAL array, dimension (LDA,K)
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
!>          T is REAL array, dimension
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
!>          C is REAL array, dimension (LDC,N)
!>          On entry, the M-by-N matrix C.
!>          On exit, C is overwritten by Q*C or Q**T*C or C*Q**T or C*Q.
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
!>         (workspace) REAL array, dimension (MAX(1,LWORK))
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
      SUBROUTINE SLAMTSQR(Side,Trans,M,N,K,Mb,Nb,A,Lda,T,Ldt,C,Ldc,Work,&
     &                    Lwork,Info)
      IMPLICIT NONE
!*--SLAMTSQR200
!
!  -- LAPACK computational routine (version 3.7.1) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     June 2017
!
!     .. Scalar Arguments ..
      CHARACTER Side , Trans
      INTEGER Info , Lda , M , N , K , Mb , Nb , Ldt , Lwork , Ldc
!     ..
!     .. Array Arguments ..
      REAL A(Lda,*) , Work(*) , C(Ldc,*) , T(Ldt,*)
!     ..
!
! =====================================================================
!
!     ..
!     .. Local Scalars ..
      LOGICAL left , right , tran , notran , lquery
      INTEGER i , ii , kk , lw , ctr
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
!     .. External Subroutines ..
      EXTERNAL SGEMQRT , STPMQRT , XERBLA
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      lquery = Lwork<0
      notran = LSAME(Trans,'N')
      tran = LSAME(Trans,'T')
      left = LSAME(Side,'L')
      right = LSAME(Side,'R')
      IF ( left ) THEN
         lw = N*Nb
      ELSE
         lw = Mb*Nb
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
         CALL XERBLA('SLAMTSQR',-Info)
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
         CALL SGEMQRT(Side,Trans,M,N,K,Nb,A,Lda,T,Ldt,C,Ldc,Work,Info)
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
            CALL STPMQRT('L','N',kk,N,K,0,Nb,A(ii,1),Lda,T(1,ctr*K+1),  &
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
            CALL STPMQRT('L','N',Mb-K,N,K,0,Nb,A(i,1),Lda,T(1,ctr*K+1), &
     &                   Ldt,C(1,1),Ldc,C(i,1),Ldc,Work,Info)
!
         ENDDO
!
!         Multiply Q to the first block of C (1:MB,1:N)
!
         CALL SGEMQRT('L','N',Mb,N,K,Nb,A(1,1),Lda,T,Ldt,C(1,1),Ldc,    &
     &                Work,Info)
!
      ELSEIF ( left .AND. tran ) THEN
!
!         Multiply Q to the first block of C
!
         kk = MOD((M-K),(Mb-K))
         ii = M - kk + 1
         ctr = 1
         CALL SGEMQRT('L','T',Mb,N,K,Nb,A(1,1),Lda,T,Ldt,C(1,1),Ldc,    &
     &                Work,Info)
!
         DO i = Mb + 1 , ii - Mb + K , (Mb-K)
!
!         Multiply Q to the current block of C (I:I+MB,1:N)
!
            CALL STPMQRT('L','T',Mb-K,N,K,0,Nb,A(i,1),Lda,T(1,ctr*K+1), &
     &                   Ldt,C(1,1),Ldc,C(i,1),Ldc,Work,Info)
            ctr = ctr + 1
!
         ENDDO
!
!         Multiply Q to the last block of C
!
!
         IF ( ii<=M ) CALL STPMQRT('L','T',kk,N,K,0,Nb,A(ii,1),Lda,     &
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
            CALL STPMQRT('R','T',M,kk,K,0,Nb,A(ii,1),Lda,T(1,ctr*K+1),  &
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
            CALL STPMQRT('R','T',M,Mb-K,K,0,Nb,A(i,1),Lda,T(1,ctr*K+1), &
     &                   Ldt,C(1,1),Ldc,C(1,i),Ldc,Work,Info)
!
         ENDDO
!
!         Multiply Q to the first block of C (1:M,1:MB)
!
         CALL SGEMQRT('R','T',M,Mb,K,Nb,A(1,1),Lda,T,Ldt,C(1,1),Ldc,    &
     &                Work,Info)
!
      ELSEIF ( right .AND. notran ) THEN
!
!         Multiply Q to the first block of C
!
         kk = MOD((N-K),(Mb-K))
         ii = N - kk + 1
         ctr = 1
         CALL SGEMQRT('R','N',M,Mb,K,Nb,A(1,1),Lda,T,Ldt,C(1,1),Ldc,    &
     &                Work,Info)
!
         DO i = Mb + 1 , ii - Mb + K , (Mb-K)
!
!         Multiply Q to the current block of C (1:M,I:I+MB)
!
            CALL STPMQRT('R','N',M,Mb-K,K,0,Nb,A(i,1),Lda,T(1,ctr*K+1), &
     &                   Ldt,C(1,1),Ldc,C(1,i),Ldc,Work,Info)
            ctr = ctr + 1
!
         ENDDO
!
!         Multiply Q to the last block of C
!
!
         IF ( ii<=N ) CALL STPMQRT('R','N',M,kk,K,0,Nb,A(ii,1),Lda,     &
     &                             T(1,ctr*K+1),Ldt,C(1,1),Ldc,C(1,ii), &
     &                             Ldc,Work,Info)
!
      ENDIF
!
      Work(1) = lw
!
!     End of SLAMTSQR
!
      END SUBROUTINE SLAMTSQR
