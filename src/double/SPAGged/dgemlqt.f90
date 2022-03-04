!*==dgemlqt.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b DGEMLQT
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DGEMLQT + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgemlqt.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgemlqt.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgemlqt.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DGEMLQT( SIDE, TRANS, M, N, K, MB, V, LDV, T, LDT,
!                          C, LDC, WORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER SIDE, TRANS
!       INTEGER   INFO, K, LDV, LDC, M, N, MB, LDT
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION V( LDV, * ), C( LDC, * ), T( LDT, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DGEMLQT overwrites the general real M-by-N matrix C with
!>
!>                 SIDE = 'L'     SIDE = 'R'
!> TRANS = 'N':      Q C            C Q
!> TRANS = 'T':   Q**T C            C Q**T
!>
!> where Q is a real orthogonal matrix defined as the product of K
!> elementary reflectors:
!>
!>       Q = H(1) H(2) . . . H(K) = I - V T V**T
!>
!> generated using the compact WY representation as returned by DGELQT.
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
!>          = 'L': apply Q or Q**T from the Left;
!>          = 'R': apply Q or Q**T from the Right.
!> \endverbatim
!>
!> \param[in] TRANS
!> \verbatim
!>          TRANS is CHARACTER*1
!>          = 'N':  No transpose, apply Q;
!>          = 'C':  Transpose, apply Q**T.
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
!> \param[in] MB
!> \verbatim
!>          MB is INTEGER
!>          The block size used for the storage of T.  K >= MB >= 1.
!>          This must be the same value of MB used to generate T
!>          in DGELQT.
!> \endverbatim
!>
!> \param[in] V
!> \verbatim
!>          V is DOUBLE PRECISION array, dimension
!>                               (LDV,M) if SIDE = 'L',
!>                               (LDV,N) if SIDE = 'R'
!>          The i-th row must contain the vector which defines the
!>          elementary reflector H(i), for i = 1,2,...,k, as returned by
!>          DGELQT in the first K rows of its array argument A.
!> \endverbatim
!>
!> \param[in] LDV
!> \verbatim
!>          LDV is INTEGER
!>          The leading dimension of the array V. LDV >= max(1,K).
!> \endverbatim
!>
!> \param[in] T
!> \verbatim
!>          T is DOUBLE PRECISION array, dimension (LDT,K)
!>          The upper triangular factors of the block reflectors
!>          as returned by DGELQT, stored as a MB-by-K matrix.
!> \endverbatim
!>
!> \param[in] LDT
!> \verbatim
!>          LDT is INTEGER
!>          The leading dimension of the array T.  LDT >= MB.
!> \endverbatim
!>
!> \param[in,out] C
!> \verbatim
!>          C is DOUBLE PRECISION array, dimension (LDC,N)
!>          On entry, the M-by-N matrix C.
!>          On exit, C is overwritten by Q C, Q**T C, C Q**T or C Q.
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
!>          WORK is DOUBLE PRECISION array. The dimension of
!>          WORK is N*MB if SIDE = 'L', or  M*MB if SIDE = 'R'.
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
!> \ingroup doubleGEcomputational
!
!  =====================================================================
      SUBROUTINE DGEMLQT(Side,Trans,M,N,K,Mb,V,Ldv,T,Ldt,C,Ldc,Work,    &
     &                   Info)
      USE F77KINDS                        
      USE S_DLARFB
      USE S_LSAME
      USE S_XERBLA
      IMPLICIT NONE
!*--DGEMLQT176
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Side
      CHARACTER :: Trans
      INTEGER :: M
      INTEGER :: N
      INTEGER , INTENT(IN) :: K
      INTEGER , INTENT(IN) :: Mb
      REAL(R8KIND) , DIMENSION(Ldv,*) :: V
      INTEGER :: Ldv
      REAL(R8KIND) , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      REAL(R8KIND) , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      INTEGER :: i , ib , kf , ldwork
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
      tran = LSAME(Trans,'T')
      notran = LSAME(Trans,'N')
!
      IF ( left ) THEN
         ldwork = MAX(1,N)
      ELSEIF ( right ) THEN
         ldwork = MAX(1,M)
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
      ELSEIF ( Mb<1 .OR. (Mb>K .AND. K>0) ) THEN
         Info = -6
      ELSEIF ( Ldv<MAX(1,K) ) THEN
         Info = -8
      ELSEIF ( Ldt<Mb ) THEN
         Info = -10
      ELSEIF ( Ldc<MAX(1,M) ) THEN
         Info = -12
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('DGEMLQT',-Info)
         RETURN
      ENDIF
!
!     .. Quick return if possible ..
!
      IF ( M==0 .OR. N==0 .OR. K==0 ) RETURN
!
      IF ( left .AND. notran ) THEN
!
         DO i = 1 , K , Mb
            ib = MIN(Mb,K-i+1)
            CALL DLARFB('L','T','F','R',M-i+1,N,ib,V(i,i),Ldv,T(1,i),   &
     &                  Ldt,C(i,1),Ldc,Work,ldwork)
         ENDDO
!
      ELSEIF ( right .AND. tran ) THEN
!
         DO i = 1 , K , Mb
            ib = MIN(Mb,K-i+1)
            CALL DLARFB('R','N','F','R',M,N-i+1,ib,V(i,i),Ldv,T(1,i),   &
     &                  Ldt,C(1,i),Ldc,Work,ldwork)
         ENDDO
!
      ELSEIF ( left .AND. tran ) THEN
!
         kf = ((K-1)/Mb)*Mb + 1
         DO i = kf , 1 , -Mb
            ib = MIN(Mb,K-i+1)
            CALL DLARFB('L','N','F','R',M-i+1,N,ib,V(i,i),Ldv,T(1,i),   &
     &                  Ldt,C(i,1),Ldc,Work,ldwork)
         ENDDO
!
      ELSEIF ( right .AND. notran ) THEN
!
         kf = ((K-1)/Mb)*Mb + 1
         DO i = kf , 1 , -Mb
            ib = MIN(Mb,K-i+1)
            CALL DLARFB('R','T','F','R',M,N-i+1,ib,V(i,i),Ldv,T(1,i),   &
     &                  Ldt,C(1,i),Ldc,Work,ldwork)
         ENDDO
!
      ENDIF
!
!
!     End of DGEMLQT
!
      END SUBROUTINE DGEMLQT
