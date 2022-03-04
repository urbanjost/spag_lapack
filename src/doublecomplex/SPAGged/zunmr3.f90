!*==zunmr3.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b ZUNMR3 multiplies a general matrix by the unitary matrix from a RZ factorization determined by ctzrzf (unblocked algorithm).
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZUNMR3 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zunmr3.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zunmr3.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zunmr3.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZUNMR3( SIDE, TRANS, M, N, K, L, A, LDA, TAU, C, LDC,
!                          WORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          SIDE, TRANS
!       INTEGER            INFO, K, L, LDA, LDC, M, N
!       ..
!       .. Array Arguments ..
!       COMPLEX*16         A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZUNMR3 overwrites the general complex m by n matrix C with
!>
!>       Q * C  if SIDE = 'L' and TRANS = 'N', or
!>
!>       Q**H* C  if SIDE = 'L' and TRANS = 'C', or
!>
!>       C * Q  if SIDE = 'R' and TRANS = 'N', or
!>
!>       C * Q**H if SIDE = 'R' and TRANS = 'C',
!>
!> where Q is a complex unitary matrix defined as the product of k
!> elementary reflectors
!>
!>       Q = H(1) H(2) . . . H(k)
!>
!> as returned by ZTZRZF. Q is of order m if SIDE = 'L' and of order n
!> if SIDE = 'R'.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] SIDE
!> \verbatim
!>          SIDE is CHARACTER*1
!>          = 'L': apply Q or Q**H from the Left
!>          = 'R': apply Q or Q**H from the Right
!> \endverbatim
!>
!> \param[in] TRANS
!> \verbatim
!>          TRANS is CHARACTER*1
!>          = 'N': apply Q  (No transpose)
!>          = 'C': apply Q**H (Conjugate transpose)
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
!> \param[in] L
!> \verbatim
!>          L is INTEGER
!>          The number of columns of the matrix A containing
!>          the meaningful part of the Householder reflectors.
!>          If SIDE = 'L', M >= L >= 0, if SIDE = 'R', N >= L >= 0.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension
!>                               (LDA,M) if SIDE = 'L',
!>                               (LDA,N) if SIDE = 'R'
!>          The i-th row must contain the vector which defines the
!>          elementary reflector H(i), for i = 1,2,...,k, as returned by
!>          ZTZRZF in the last k rows of its array argument A.
!>          A is modified by the routine but restored on exit.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A. LDA >= max(1,K).
!> \endverbatim
!>
!> \param[in] TAU
!> \verbatim
!>          TAU is COMPLEX*16 array, dimension (K)
!>          TAU(i) must contain the scalar factor of the elementary
!>          reflector H(i), as returned by ZTZRZF.
!> \endverbatim
!>
!> \param[in,out] C
!> \verbatim
!>          C is COMPLEX*16 array, dimension (LDC,N)
!>          On entry, the m-by-n matrix C.
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
!>          WORK is COMPLEX*16 array, dimension
!>                                   (N) if SIDE = 'L',
!>                                   (M) if SIDE = 'R'
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0: successful exit
!>          < 0: if INFO = -i, the i-th argument had an illegal value
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
!> \ingroup complex16OTHERcomputational
!
!> \par Contributors:
!  ==================
!>
!>    A. Petitet, Computer Science Dept., Univ. of Tenn., Knoxville, USA
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE ZUNMR3(Side,Trans,M,N,K,L,A,Lda,Tau,C,Ldc,Work,Info)
      USE F77KINDS                        
      USE S_LSAME
      USE S_XERBLA
      USE S_ZLARZ
      IMPLICIT NONE
!*--ZUNMR3185
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Side
      CHARACTER :: Trans
      INTEGER , INTENT(IN) :: M
      INTEGER , INTENT(IN) :: N
      INTEGER , INTENT(IN) :: K
      INTEGER :: L
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(*) :: Tau
      COMPLEX(CX16KIND) , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      INTEGER :: i , i1 , i2 , i3 , ic , ja , jc , mi , ni , nq
      LOGICAL :: left , notran
      COMPLEX(CX16KIND) :: taui
!
! End of declarations rewritten by SPAG
!
!     ..
!     .. Array Arguments ..
!     ..
!
!  =====================================================================
!
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
      Info = 0
      left = LSAME(Side,'L')
      notran = LSAME(Trans,'N')
!
!     NQ is the order of Q
!
      IF ( left ) THEN
         nq = M
      ELSE
         nq = N
      ENDIF
      IF ( .NOT.left .AND. .NOT.LSAME(Side,'R') ) THEN
         Info = -1
      ELSEIF ( .NOT.notran .AND. .NOT.LSAME(Trans,'C') ) THEN
         Info = -2
      ELSEIF ( M<0 ) THEN
         Info = -3
      ELSEIF ( N<0 ) THEN
         Info = -4
      ELSEIF ( K<0 .OR. K>nq ) THEN
         Info = -5
      ELSEIF ( L<0 .OR. (left .AND. (L>M)) .OR. (.NOT.left .AND. (L>N)) &
     &         ) THEN
         Info = -6
      ELSEIF ( Lda<MAX(1,K) ) THEN
         Info = -8
      ELSEIF ( Ldc<MAX(1,M) ) THEN
         Info = -11
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('ZUNMR3',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( M==0 .OR. N==0 .OR. K==0 ) RETURN
!
      IF ( (left .AND. .NOT.notran .OR. .NOT.left .AND. notran) ) THEN
         i1 = 1
         i2 = K
         i3 = 1
      ELSE
         i1 = K
         i2 = 1
         i3 = -1
      ENDIF
!
      IF ( left ) THEN
         ni = N
         ja = M - L + 1
         jc = 1
      ELSE
         mi = M
         ja = N - L + 1
         ic = 1
      ENDIF
!
      DO i = i1 , i2 , i3
         IF ( left ) THEN
!
!           H(i) or H(i)**H is applied to C(i:m,1:n)
!
            mi = M - i + 1
            ic = i
         ELSE
!
!           H(i) or H(i)**H is applied to C(1:m,i:n)
!
            ni = N - i + 1
            jc = i
         ENDIF
!
!        Apply H(i) or H(i)**H
!
         IF ( notran ) THEN
            taui = Tau(i)
         ELSE
            taui = DCONJG(Tau(i))
         ENDIF
         CALL ZLARZ(Side,mi,ni,L,A(i,ja),Lda,taui,C(ic,jc),Ldc,Work)
!
      ENDDO
!
!
!     End of ZUNMR3
!
      END SUBROUTINE ZUNMR3
