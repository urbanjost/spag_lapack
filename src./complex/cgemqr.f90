!*==cgemqr.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b CGEMQR
!
!  Definition:
!  ===========
!
!      SUBROUTINE CGEMQR( SIDE, TRANS, M, N, K, A, LDA, T,
!     $                   TSIZE, C, LDC, WORK, LWORK, INFO )
!
!
!     .. Scalar Arguments ..
!     CHARACTER         SIDE, TRANS
!     INTEGER           INFO, LDA, M, N, K, LDT, TSIZE, LWORK, LDC
!     ..
!     .. Array Arguments ..
!     COMPLEX           A( LDA, * ), T( * ), C( LDC, * ), WORK( * )
!     ..
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CGEMQR overwrites the general real M-by-N matrix C with
!>
!>                      SIDE = 'L'     SIDE = 'R'
!>      TRANS = 'N':      Q * C          C * Q
!>      TRANS = 'T':      Q**H * C       C * Q**H
!>
!> where Q is a complex unitary matrix defined as the product
!> of blocked elementary reflectors computed by tall skinny
!> QR factorization (CGEQR)
!>
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
!>          The number of columns of the matrix C. N >= 0.
!> \endverbatim
!>
!> \param[in] K
!> \verbatim
!>          K is INTEGER
!>          The number of elementary reflectors whose product defines
!>          the matrix Q.
!>          If SIDE = 'L', M >= K >= 0;
!>          if SIDE = 'R', N >= K >= 0.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is COMPLEX array, dimension (LDA,K)
!>          Part of the data structure to represent Q as returned by CGEQR.
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
!>          T is COMPLEX array, dimension (MAX(5,TSIZE)).
!>          Part of the data structure to represent Q as returned by CGEQR.
!> \endverbatim
!>
!> \param[in] TSIZE
!> \verbatim
!>          TSIZE is INTEGER
!>          The dimension of the array T. TSIZE >= 5.
!> \endverbatim
!>
!> \param[in,out] C
!> \verbatim
!>          C is COMPLEX array, dimension (LDC,N)
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
!>         (workspace) COMPLEX array, dimension (MAX(1,LWORK))
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The dimension of the array WORK.
!>          If LWORK = -1, then a workspace query is assumed. The routine
!>          only calculates the size of the WORK array, returns this
!>          value as WORK(1), and no error message related to WORK
!>          is issued by XERBLA.
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
!> \par Further Details
!  ====================
!>
!> \verbatim
!>
!> These details are particular for this LAPACK implementation. Users should not
!> take them for granted. These details may change in the future, and are not likely
!> true for another LAPACK implementation. These details are relevant if one wants
!> to try to understand the code. They are not part of the interface.
!>
!> In this version,
!>
!>          T(2): row block size (MB)
!>          T(3): column block size (NB)
!>          T(6:TSIZE): data structure needed for Q, computed by
!>                           CLATSQR or CGEQRT
!>
!>  Depending on the matrix dimensions M and N, and row and column
!>  block sizes MB and NB returned by ILAENV, CGEQR will use either
!>  CLATSQR (if the matrix is tall-and-skinny) or CGEQRT to compute
!>  the QR factorization.
!>  This version of CGEMQR will use either CLAMTSQR or CGEMQRT to
!>  multiply matrix Q by another matrix.
!>  Further Details in CLAMTSQR or CGEMQRT.
!>
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE CGEMQR(Side,Trans,M,N,K,A,Lda,T,Tsize,C,Ldc,Work,Lwork,&
     &                  Info)
      USE S_CGEMQRT
      USE S_CLAMTSQR
      USE S_LSAME
      USE S_XERBLA
      IMPLICIT NONE
!*--CGEMQR178
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Side
      CHARACTER :: Trans
      INTEGER :: M
      INTEGER :: N
      INTEGER :: K
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(*) :: T
      INTEGER , INTENT(IN) :: Tsize
      COMPLEX , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      COMPLEX , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      LOGICAL :: left , lquery , notran , right , tran
      INTEGER :: lw , mb , mn , nb , nblcks
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
!     ..
!     .. External Subroutines ..
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      lquery = Lwork== - 1
      notran = LSAME(Trans,'N')
      tran = LSAME(Trans,'C')
      left = LSAME(Side,'L')
      right = LSAME(Side,'R')
!
      mb = INT(T(2))
      nb = INT(T(3))
      IF ( left ) THEN
         lw = N*nb
         mn = M
      ELSE
         lw = mb*nb
         mn = N
      ENDIF
!
      IF ( (mb<=K) .OR. (mn<=K) ) THEN
         nblcks = 1
      ELSEIF ( MOD(mn-K,mb-K)==0 ) THEN
         nblcks = (mn-K)/(mb-K)
      ELSE
         nblcks = (mn-K)/(mb-K) + 1
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
      ELSEIF ( K<0 .OR. K>mn ) THEN
         Info = -5
      ELSEIF ( Lda<MAX(1,mn) ) THEN
         Info = -7
      ELSEIF ( Tsize<5 ) THEN
         Info = -9
      ELSEIF ( Ldc<MAX(1,M) ) THEN
         Info = -11
      ELSEIF ( (Lwork<MAX(1,lw)) .AND. (.NOT.lquery) ) THEN
         Info = -13
      ENDIF
!
      IF ( Info==0 ) Work(1) = lw
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('CGEMQR',-Info)
         RETURN
      ELSEIF ( lquery ) THEN
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( MIN(M,N,K)==0 ) RETURN
!
      IF ( (left .AND. M<=K) .OR. (right .AND. N<=K) .OR. (mb<=K) .OR.  &
     &     (mb>=MAX(M,N,K)) ) THEN
         CALL CGEMQRT(Side,Trans,M,N,K,nb,A,Lda,T(6),nb,C,Ldc,Work,Info)
      ELSE
         CALL CLAMTSQR(Side,Trans,M,N,K,mb,nb,A,Lda,T(6),nb,C,Ldc,Work, &
     &                 Lwork,Info)
      ENDIF
!
      Work(1) = lw
!
!
!     End of CGEMQR
!
      END SUBROUTINE CGEMQR
