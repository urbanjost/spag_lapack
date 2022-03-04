!*==zgemlq.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b ZGEMLQ
!
!  Definition:
!  ===========
!
!      SUBROUTINE ZGEMLQ( SIDE, TRANS, M, N, K, A, LDA, T,
!     $                   TSIZE, C, LDC, WORK, LWORK, INFO )
!
!
!     .. Scalar Arguments ..
!      CHARACTER          SIDE, TRANS
!      INTEGER            INFO, LDA, M, N, K, LDT, TSIZE, LWORK, LDC
!     ..
!     .. Array Arguments ..
!      COMPLEX*16         A( LDA, * ), T( * ), C(LDC, * ), WORK( * )
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>     ZGEMLQ overwrites the general real M-by-N matrix C with
!>
!>                      SIDE = 'L'     SIDE = 'R'
!>      TRANS = 'N':      Q * C          C * Q
!>      TRANS = 'C':      Q**H * C       C * Q**H
!>      where Q is a complex unitary matrix defined as the product
!>      of blocked elementary reflectors computed by short wide
!>      LQ factorization (ZGELQ)
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
!>
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension
!>                               (LDA,M) if SIDE = 'L',
!>                               (LDA,N) if SIDE = 'R'
!>          Part of the data structure to represent Q as returned by ZGELQ.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A. LDA >= max(1,K).
!> \endverbatim
!>
!> \param[in] T
!> \verbatim
!>          T is COMPLEX*16 array, dimension (MAX(5,TSIZE)).
!>          Part of the data structure to represent Q as returned by ZGELQ.
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
!>          C is COMPLEX*16 array, dimension (LDC,N)
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
!>         (workspace) COMPLEX*16 array, dimension (MAX(1,LWORK))
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
!>                           ZLASWLQ or ZGELQT
!>
!>  Depending on the matrix dimensions M and N, and row and column
!>  block sizes MB and NB returned by ILAENV, ZGELQ will use either
!>  ZLASWLQ (if the matrix is wide-and-short) or ZGELQT to compute
!>  the LQ factorization.
!>  This version of ZGEMLQ will use either ZLAMSWLQ or ZGEMLQT to
!>  multiply matrix Q by another matrix.
!>  Further Details in ZLAMSWLQ or ZGEMLQT.
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE ZGEMLQ(Side,Trans,M,N,K,A,Lda,T,Tsize,C,Ldc,Work,Lwork,&
     &                  Info)
      USE F77KINDS                        
      USE S_LSAME
      USE S_XERBLA
      USE S_ZGEMLQT
      USE S_ZLAMSWLQ
      IMPLICIT NONE
!*--ZGEMLQ176
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Side
      CHARACTER :: Trans
      INTEGER :: M
      INTEGER :: N
      INTEGER :: K
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(*) :: T
      INTEGER , INTENT(IN) :: Tsize
      COMPLEX(CX16KIND) , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
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
         lw = N*mb
         mn = M
      ELSE
         lw = M*mb
         mn = N
      ENDIF
!
      IF ( (nb<=K) .OR. (mn<=K) ) THEN
         nblcks = 1
      ELSEIF ( MOD(mn-K,nb-K)==0 ) THEN
         nblcks = (mn-K)/(nb-K)
      ELSE
         nblcks = (mn-K)/(nb-K) + 1
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
      ELSEIF ( Lda<MAX(1,K) ) THEN
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
         CALL XERBLA('ZGEMLQ',-Info)
         RETURN
      ELSEIF ( lquery ) THEN
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( MIN(M,N,K)==0 ) RETURN
!
      IF ( (left .AND. M<=K) .OR. (right .AND. N<=K) .OR. (nb<=K) .OR.  &
     &     (nb>=MAX(M,N,K)) ) THEN
         CALL ZGEMLQT(Side,Trans,M,N,K,mb,A,Lda,T(6),mb,C,Ldc,Work,Info)
      ELSE
         CALL ZLAMSWLQ(Side,Trans,M,N,K,mb,nb,A,Lda,T(6),mb,C,Ldc,Work, &
     &                 Lwork,Info)
      ENDIF
!
      Work(1) = lw
!
!
!     End of ZGEMLQ
!
      END SUBROUTINE ZGEMLQ
