!*==stpmqrt.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b STPMQRT
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download STPMQRT + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/stpmqrt.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/stpmqrt.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/stpmqrt.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE STPMQRT( SIDE, TRANS, M, N, K, L, NB, V, LDV, T, LDT,
!                           A, LDA, B, LDB, WORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER SIDE, TRANS
!       INTEGER   INFO, K, LDV, LDA, LDB, M, N, L, NB, LDT
!       ..
!       .. Array Arguments ..
!       REAL   V( LDV, * ), A( LDA, * ), B( LDB, * ), T( LDT, * ),
!      $          WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> STPMQRT applies a real orthogonal matrix Q obtained from a
!> "triangular-pentagonal" real block reflector H to a general
!> real matrix C, which consists of two blocks A and B.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] SIDE
!> \verbatim
!>          SIDE is CHARACTER*1
!>          = 'L': apply Q or Q^T from the Left;
!>          = 'R': apply Q or Q^T from the Right.
!> \endverbatim
!>
!> \param[in] TRANS
!> \verbatim
!>          TRANS is CHARACTER*1
!>          = 'N':  No transpose, apply Q;
!>          = 'T':  Transpose, apply Q^T.
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
!>          V is REAL array, dimension (LDV,K)
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
!>          T is REAL array, dimension (LDT,K)
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
!>          A is REAL array, dimension
!>          (LDA,N) if SIDE = 'L' or
!>          (LDA,K) if SIDE = 'R'
!>          On entry, the K-by-N or M-by-K matrix A.
!>          On exit, A is overwritten by the corresponding block of
!>          Q*C or Q^T*C or C*Q or C*Q^T.  See Further Details.
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
!>          B is REAL array, dimension (LDB,N)
!>          On entry, the M-by-N matrix B.
!>          On exit, B is overwritten by the corresponding block of
!>          Q*C or Q^T*C or C*Q or C*Q^T.  See Further Details.
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
!>          WORK is REAL array. The dimension of WORK is
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
!> \ingroup realOTHERcomputational
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
!>  The real orthogonal matrix Q is formed from V and T.
!>
!>  If TRANS='N' and SIDE='L', C is on exit replaced with Q * C.
!>
!>  If TRANS='T' and SIDE='L', C is on exit replaced with Q^T * C.
!>
!>  If TRANS='N' and SIDE='R', C is on exit replaced with C * Q.
!>
!>  If TRANS='T' and SIDE='R', C is on exit replaced with C * Q^T.
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE STPMQRT(Side,Trans,M,N,K,L,Nb,V,Ldv,T,Ldt,A,Lda,B,Ldb, &
     &                   Work,Info)
      IMPLICIT NONE
!*--STPMQRT220
!
!  -- LAPACK computational routine (version 3.8.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2017
!
!     .. Scalar Arguments ..
      CHARACTER Side , Trans
      INTEGER Info , K , Ldv , Lda , Ldb , M , N , L , Nb , Ldt
!     ..
!     .. Array Arguments ..
      REAL V(Ldv,*) , A(Lda,*) , B(Ldb,*) , T(Ldt,*) , Work(*)
!     ..
!
!  =====================================================================
!
!     ..
!     .. Local Scalars ..
      LOGICAL left , right , tran , notran
      INTEGER i , ib , mb , lb , kf , ldaq , ldvq
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL STPRFB , XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX , MIN
!     ..
!     .. Executable Statements ..
!
!     .. Test the input arguments ..
!
      Info = 0
      left = LSAME(Side,'L')
      right = LSAME(Side,'R')
      tran = LSAME(Trans,'T')
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
         CALL XERBLA('STPMQRT',-Info)
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
            CALL STPRFB('L','T','F','C',mb,N,ib,lb,V(1,i),Ldv,T(1,i),   &
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
            CALL STPRFB('R','N','F','C',M,mb,ib,lb,V(1,i),Ldv,T(1,i),   &
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
            CALL STPRFB('L','N','F','C',mb,N,ib,lb,V(1,i),Ldv,T(1,i),   &
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
            CALL STPRFB('R','T','F','C',M,mb,ib,lb,V(1,i),Ldv,T(1,i),   &
     &                  Ldt,A(1,i),Lda,B,Ldb,Work,M)
         ENDDO
!
      ENDIF
!
!
!     End of STPMQRT
!
      END SUBROUTINE STPMQRT
