!*==ztpmqrt.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b ZTPMQRT
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZTPMQRT + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ztpmqrt.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ztpmqrt.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ztpmqrt.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZTPMQRT( SIDE, TRANS, M, N, K, L, NB, V, LDV, T, LDT,
!                           A, LDA, B, LDB, WORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER SIDE, TRANS
!       INTEGER   INFO, K, LDV, LDA, LDB, M, N, L, NB, LDT
!       ..
!       .. Array Arguments ..
!       COMPLEX*16   V( LDV, * ), A( LDA, * ), B( LDB, * ), T( LDT, * ),
!      $          WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZTPMQRT applies a complex orthogonal matrix Q obtained from a
!> "triangular-pentagonal" complex block reflector H to a general
!> complex matrix C, which consists of two blocks A and B.
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
!>          = 'C':  Transpose, apply Q**H.
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix B. M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix B. N >= 0.
!> \endverbatim
!>
!> \param[in] K
!> \verbatim
!>          K is INTEGER
!>          The number of elementary reflectors whose product defines
!>          the matrix Q.
!> \endverbatim
!>
!> \param[in] L
!> \verbatim
!>          L is INTEGER
!>          The order of the trapezoidal part of V.
!>          K >= L >= 0.  See Further Details.
!> \endverbatim
!>
!> \param[in] NB
!> \verbatim
!>          NB is INTEGER
!>          The block size used for the storage of T.  K >= NB >= 1.
!>          This must be the same value of NB used to generate T
!>          in CTPQRT.
!> \endverbatim
!>
!> \param[in] V
!> \verbatim
!>          V is COMPLEX*16 array, dimension (LDV,K)
!>          The i-th column must contain the vector which defines the
!>          elementary reflector H(i), for i = 1,2,...,k, as returned by
!>          CTPQRT in B.  See Further Details.
!> \endverbatim
!>
!> \param[in] LDV
!> \verbatim
!>          LDV is INTEGER
!>          The leading dimension of the array V.
!>          If SIDE = 'L', LDV >= max(1,M);
!>          if SIDE = 'R', LDV >= max(1,N).
!> \endverbatim
!>
!> \param[in] T
!> \verbatim
!>          T is COMPLEX*16 array, dimension (LDT,K)
!>          The upper triangular factors of the block reflectors
!>          as returned by CTPQRT, stored as a NB-by-K matrix.
!> \endverbatim
!>
!> \param[in] LDT
!> \verbatim
!>          LDT is INTEGER
!>          The leading dimension of the array T.  LDT >= NB.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension
!>          (LDA,N) if SIDE = 'L' or
!>          (LDA,K) if SIDE = 'R'
!>          On entry, the K-by-N or M-by-K matrix A.
!>          On exit, A is overwritten by the corresponding block of
!>          Q*C or Q**H*C or C*Q or C*Q**H.  See Further Details.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.
!>          If SIDE = 'L', LDC >= max(1,K);
!>          If SIDE = 'R', LDC >= max(1,M).
!> \endverbatim
!>
!> \param[in,out] B
!> \verbatim
!>          B is COMPLEX*16 array, dimension (LDB,N)
!>          On entry, the M-by-N matrix B.
!>          On exit, B is overwritten by the corresponding block of
!>          Q*C or Q**H*C or C*Q or C*Q**H.  See Further Details.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of the array B.
!>          LDB >= max(1,M).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX*16 array. The dimension of WORK is
!>           N*NB if SIDE = 'L', or  M*NB if SIDE = 'R'.
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
!> \date November 2017
!
!> \ingroup complex16OTHERcomputational
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  The columns of the pentagonal matrix V contain the elementary reflectors
!>  H(1), H(2), ..., H(K); V is composed of a rectangular block V1 and a
!>  trapezoidal block V2:
!>
!>        V = [V1]
!>            [V2].
!>
!>  The size of the trapezoidal block V2 is determined by the parameter L,
!>  where 0 <= L <= K; V2 is upper trapezoidal, consisting of the first L
!>  rows of a K-by-K upper triangular matrix.  If L=K, V2 is upper triangular;
!>  if L=0, there is no trapezoidal block, hence V = V1 is rectangular.
!>
!>  If SIDE = 'L':  C = [A]  where A is K-by-N,  B is M-by-N and V is M-by-K.
!>                      [B]
!>
!>  If SIDE = 'R':  C = [A B]  where A is M-by-K, B is M-by-N and V is N-by-K.
!>
!>  The complex orthogonal matrix Q is formed from V and T.
!>
!>  If TRANS='N' and SIDE='L', C is on exit replaced with Q * C.
!>
!>  If TRANS='C' and SIDE='L', C is on exit replaced with Q**H * C.
!>
!>  If TRANS='N' and SIDE='R', C is on exit replaced with C * Q.
!>
!>  If TRANS='C' and SIDE='R', C is on exit replaced with C * Q**H.
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE ZTPMQRT(Side,Trans,M,N,K,L,Nb,V,Ldv,T,Ldt,A,Lda,B,Ldb, &
     &                   Work,Info)
      USE F77KINDS                        
      USE S_LSAME
      USE S_XERBLA
      USE S_ZTPRFB
      IMPLICIT NONE
!*--ZTPMQRT224
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Side
      CHARACTER :: Trans
      INTEGER :: M
      INTEGER :: N
      INTEGER , INTENT(IN) :: K
      INTEGER , INTENT(IN) :: L
      INTEGER , INTENT(IN) :: Nb
      COMPLEX(CX16KIND) , DIMENSION(Ldv,*) :: V
      INTEGER :: Ldv
      COMPLEX(CX16KIND) , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      INTEGER :: i , ib , kf , lb , ldaq , ldvq , mb
      LOGICAL :: left , notran , right , tran
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
!     .. External Functions ..
!     ..
!     .. External Subroutines ..
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. Executable Statements ..
!
!     .. Test the input arguments ..
!
      Info = 0
      left = LSAME(Side,'L')
      right = LSAME(Side,'R')
      tran = LSAME(Trans,'C')
      notran = LSAME(Trans,'N')
!
      IF ( left ) THEN
         ldvq = MAX(1,M)
         ldaq = MAX(1,K)
      ELSEIF ( right ) THEN
         ldvq = MAX(1,N)
         ldaq = MAX(1,M)
      ENDIF
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
      ELSEIF ( L<0 .OR. L>K ) THEN
         Info = -6
      ELSEIF ( Nb<1 .OR. (Nb>K .AND. K>0) ) THEN
         Info = -7
      ELSEIF ( Ldv<ldvq ) THEN
         Info = -9
      ELSEIF ( Ldt<Nb ) THEN
         Info = -11
      ELSEIF ( Lda<ldaq ) THEN
         Info = -13
      ELSEIF ( Ldb<MAX(1,M) ) THEN
         Info = -15
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('ZTPMQRT',-Info)
         RETURN
      ENDIF
!
!     .. Quick return if possible ..
!
      IF ( M==0 .OR. N==0 .OR. K==0 ) RETURN
!
      IF ( left .AND. tran ) THEN
!
         DO i = 1 , K , Nb
            ib = MIN(Nb,K-i+1)
            mb = MIN(M-L+i+ib-1,M)
            IF ( i>=L ) THEN
               lb = 0
            ELSE
               lb = mb - M + L - i + 1
            ENDIF
            CALL ZTPRFB('L','C','F','C',mb,N,ib,lb,V(1,i),Ldv,T(1,i),   &
     &                  Ldt,A(i,1),Lda,B,Ldb,Work,ib)
         ENDDO
!
      ELSEIF ( right .AND. notran ) THEN
!
         DO i = 1 , K , Nb
            ib = MIN(Nb,K-i+1)
            mb = MIN(N-L+i+ib-1,N)
            IF ( i>=L ) THEN
               lb = 0
            ELSE
               lb = mb - N + L - i + 1
            ENDIF
            CALL ZTPRFB('R','N','F','C',M,mb,ib,lb,V(1,i),Ldv,T(1,i),   &
     &                  Ldt,A(1,i),Lda,B,Ldb,Work,M)
         ENDDO
!
      ELSEIF ( left .AND. notran ) THEN
!
         kf = ((K-1)/Nb)*Nb + 1
         DO i = kf , 1 , -Nb
            ib = MIN(Nb,K-i+1)
            mb = MIN(M-L+i+ib-1,M)
            IF ( i>=L ) THEN
               lb = 0
            ELSE
               lb = mb - M + L - i + 1
            ENDIF
            CALL ZTPRFB('L','N','F','C',mb,N,ib,lb,V(1,i),Ldv,T(1,i),   &
     &                  Ldt,A(i,1),Lda,B,Ldb,Work,ib)
         ENDDO
!
      ELSEIF ( right .AND. tran ) THEN
!
         kf = ((K-1)/Nb)*Nb + 1
         DO i = kf , 1 , -Nb
            ib = MIN(Nb,K-i+1)
            mb = MIN(N-L+i+ib-1,N)
            IF ( i>=L ) THEN
               lb = 0
            ELSE
               lb = mb - N + L - i + 1
            ENDIF
            CALL ZTPRFB('R','C','F','C',M,mb,ib,lb,V(1,i),Ldv,T(1,i),   &
     &                  Ldt,A(1,i),Lda,B,Ldb,Work,M)
         ENDDO
!
      ENDIF
!
!
!     End of ZTPMQRT
!
      END SUBROUTINE ZTPMQRT
