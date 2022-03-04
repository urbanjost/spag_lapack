!*==zgemqrt.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b ZGEMQRT
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZGEMQRT + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zgemqrt.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zgemqrt.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zgemqrt.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZGEMQRT( SIDE, TRANS, M, N, K, NB, V, LDV, T, LDT,
!                          C, LDC, WORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER SIDE, TRANS
!       INTEGER   INFO, K, LDV, LDC, M, N, NB, LDT
!       ..
!       .. Array Arguments ..
!       COMPLEX*16   V( LDV, * ), C( LDC, * ), T( LDT, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZGEMQRT overwrites the general complex M-by-N matrix C with
!>
!>                 SIDE = 'L'     SIDE = 'R'
!> TRANS = 'N':      Q C            C Q
!> TRANS = 'C':    Q**H C            C Q**H
!>
!> where Q is a complex orthogonal matrix defined as the product of K
!> elementary reflectors:
!>
!>       Q = H(1) H(2) . . . H(K) = I - V T V**H
!>
!> generated using the compact WY representation as returned by ZGEQRT.
!>
!> Q is of order M if SIDE = 'L' and of order N  if SIDE = 'R'.
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
!>          The number of rows of the matrix C. M >= 0.
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
!> \param[in] NB
!> \verbatim
!>          NB is INTEGER
!>          The block size used for the storage of T.  K >= NB >= 1.
!>          This must be the same value of NB used to generate T
!>          in CGEQRT.
!> \endverbatim
!>
!> \param[in] V
!> \verbatim
!>          V is COMPLEX*16 array, dimension (LDV,K)
!>          The i-th column must contain the vector which defines the
!>          elementary reflector H(i), for i = 1,2,...,k, as returned by
!>          CGEQRT in the first K columns of its array argument A.
!> \endverbatim
!>
!> \param[in] LDV
!> \verbatim
!>          LDV is INTEGER
!>          The leading dimension of the array V.
!>          If SIDE = 'L', LDA >= max(1,M);
!>          if SIDE = 'R', LDA >= max(1,N).
!> \endverbatim
!>
!> \param[in] T
!> \verbatim
!>          T is COMPLEX*16 array, dimension (LDT,K)
!>          The upper triangular factors of the block reflectors
!>          as returned by CGEQRT, stored as a NB-by-N matrix.
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
!>          C is COMPLEX*16 array, dimension (LDC,N)
!>          On entry, the M-by-N matrix C.
!>          On exit, C is overwritten by Q C, Q**H C, C Q**H or C Q.
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
!> \date December 2016
!
!> \ingroup complex16GEcomputational
!
!  =====================================================================
      SUBROUTINE ZGEMQRT(Side,Trans,M,N,K,Nb,V,Ldv,T,Ldt,C,Ldc,Work,    &
     &                   Info)
      IMPLICIT NONE
!*--ZGEMQRT172
!
!  -- LAPACK computational routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Side , Trans
      INTEGER Info , K , Ldv , Ldc , M , N , Nb , Ldt
!     ..
!     .. Array Arguments ..
      COMPLEX*16 V(Ldv,*) , C(Ldc,*) , T(Ldt,*) , Work(*)
!     ..
!
!  =====================================================================
!
!     ..
!     .. Local Scalars ..
      LOGICAL left , right , tran , notran
      INTEGER i , ib , ldwork , kf , q
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL XERBLA , ZLARFB
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
      tran = LSAME(Trans,'C')
      notran = LSAME(Trans,'N')
!
      IF ( left ) THEN
         ldwork = MAX(1,N)
         q = M
      ELSEIF ( right ) THEN
         ldwork = MAX(1,M)
         q = N
      ENDIF
      IF ( .NOT.left .AND. .NOT.right ) THEN
         Info = -1
      ELSEIF ( .NOT.tran .AND. .NOT.notran ) THEN
         Info = -2
      ELSEIF ( M<0 ) THEN
         Info = -3
      ELSEIF ( N<0 ) THEN
         Info = -4
      ELSEIF ( K<0 .OR. K>q ) THEN
         Info = -5
      ELSEIF ( Nb<1 .OR. (Nb>K .AND. K>0) ) THEN
         Info = -6
      ELSEIF ( Ldv<MAX(1,q) ) THEN
         Info = -8
      ELSEIF ( Ldt<Nb ) THEN
         Info = -10
      ELSEIF ( Ldc<MAX(1,M) ) THEN
         Info = -12
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('ZGEMQRT',-Info)
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
            CALL ZLARFB('L','C','F','C',M-i+1,N,ib,V(i,i),Ldv,T(1,i),   &
     &                  Ldt,C(i,1),Ldc,Work,ldwork)
         ENDDO
!
      ELSEIF ( right .AND. notran ) THEN
!
         DO i = 1 , K , Nb
            ib = MIN(Nb,K-i+1)
            CALL ZLARFB('R','N','F','C',M,N-i+1,ib,V(i,i),Ldv,T(1,i),   &
     &                  Ldt,C(1,i),Ldc,Work,ldwork)
         ENDDO
!
      ELSEIF ( left .AND. notran ) THEN
!
         kf = ((K-1)/Nb)*Nb + 1
         DO i = kf , 1 , -Nb
            ib = MIN(Nb,K-i+1)
            CALL ZLARFB('L','N','F','C',M-i+1,N,ib,V(i,i),Ldv,T(1,i),   &
     &                  Ldt,C(i,1),Ldc,Work,ldwork)
         ENDDO
!
      ELSEIF ( right .AND. tran ) THEN
!
         kf = ((K-1)/Nb)*Nb + 1
         DO i = kf , 1 , -Nb
            ib = MIN(Nb,K-i+1)
            CALL ZLARFB('R','C','F','C',M,N-i+1,ib,V(i,i),Ldv,T(1,i),   &
     &                  Ldt,C(1,i),Ldc,Work,ldwork)
         ENDDO
!
      ENDIF
!
!
!     End of ZGEMQRT
!
      END SUBROUTINE ZGEMQRT
