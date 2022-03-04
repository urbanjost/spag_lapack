!*==cunmlq.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b CUNMLQ
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CUNMLQ + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cunmlq.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cunmlq.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cunmlq.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CUNMLQ( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC,
!                          WORK, LWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          SIDE, TRANS
!       INTEGER            INFO, K, LDA, LDC, LWORK, M, N
!       ..
!       .. Array Arguments ..
!       COMPLEX            A( LDA, * ), C( LDC, * ), TAU( * ),
!      $                   WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CUNMLQ overwrites the general complex M-by-N matrix C with
!>
!>                 SIDE = 'L'     SIDE = 'R'
!> TRANS = 'N':      Q * C          C * Q
!> TRANS = 'C':      Q**H * C       C * Q**H
!>
!> where Q is a complex unitary matrix defined as the product of k
!> elementary reflectors
!>
!>       Q = H(k)**H . . . H(2)**H H(1)**H
!>
!> as returned by CGELQF. Q is of order M if SIDE = 'L' and of order N
!> if SIDE = 'R'.
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
!>          = 'C':  Conjugate transpose, apply Q**H.
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
!> \param[in] A
!> \verbatim
!>          A is COMPLEX array, dimension
!>                               (LDA,M) if SIDE = 'L',
!>                               (LDA,N) if SIDE = 'R'
!>          The i-th row must contain the vector which defines the
!>          elementary reflector H(i), for i = 1,2,...,k, as returned by
!>          CGELQF in the first k rows of its array argument A.
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
!>          TAU is COMPLEX array, dimension (K)
!>          TAU(i) must contain the scalar factor of the elementary
!>          reflector H(i), as returned by CGELQF.
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
!>          WORK is COMPLEX array, dimension (MAX(1,LWORK))
!>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The dimension of the array WORK.
!>          If SIDE = 'L', LWORK >= max(1,N);
!>          if SIDE = 'R', LWORK >= max(1,M).
!>          For good performance, LWORK should generally be larger.
!>
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
!> \date December 2016
!
!> \ingroup complexOTHERcomputational
!
!  =====================================================================
      SUBROUTINE CUNMLQ(Side,Trans,M,N,K,A,Lda,Tau,C,Ldc,Work,Lwork,    &
     &                  Info)
      USE S_CLARFB
      USE S_CLARFT
      USE S_CUNML2
      USE S_ILAENV
      USE S_LSAME
      USE S_XERBLA
      IMPLICIT NONE
!*--CUNMLQ178
!
! PARAMETER definitions rewritten by SPAG
!
      INTEGER , PARAMETER  ::  NBMAX = 64 , LDT = NBMAX + 1 ,           &
     &                         TSIZE = LDT*NBMAX
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
      COMPLEX , DIMENSION(*) :: Tau
      COMPLEX , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      COMPLEX , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      INTEGER :: i , i1 , i2 , i3 , ib , ic , iinfo , iwt , jc ,        &
     &           ldwork , lwkopt , mi , nb , nbmin , ni , nq , nw
      LOGICAL :: left , lquery , notran
      CHARACTER :: transt
!
! End of declarations rewritten by SPAG
!
!     ..
!     .. Array Arguments ..
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
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
      Info = 0
      left = LSAME(Side,'L')
      notran = LSAME(Trans,'N')
      lquery = (Lwork==-1)
!
!     NQ is the order of Q and NW is the minimum dimension of WORK
!
      IF ( left ) THEN
         nq = M
         nw = N
      ELSE
         nq = N
         nw = M
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
      ELSEIF ( Lda<MAX(1,K) ) THEN
         Info = -7
      ELSEIF ( Ldc<MAX(1,M) ) THEN
         Info = -10
      ELSEIF ( Lwork<MAX(1,nw) .AND. .NOT.lquery ) THEN
         Info = -12
      ENDIF
!
      IF ( Info==0 ) THEN
!
!        Compute the workspace requirements
!
         IF ( M==0 .OR. N==0 .OR. K==0 ) THEN
            lwkopt = 1
         ELSE
            nb = MIN(NBMAX,ILAENV(1,'CUNMLQ',Side//Trans,M,N,K,-1))
            lwkopt = MAX(1,nw)*nb + TSIZE
         ENDIF
         Work(1) = lwkopt
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('CUNMLQ',-Info)
         RETURN
      ELSEIF ( lquery ) THEN
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( M==0 .OR. N==0 .OR. K==0 ) RETURN
!
!     Determine the block size
!
      nbmin = 2
      ldwork = nw
      IF ( nb>1 .AND. nb<K ) THEN
         IF ( Lwork<nw*nb+TSIZE ) THEN
            nb = (Lwork-TSIZE)/ldwork
            nbmin = MAX(2,ILAENV(2,'CUNMLQ',Side//Trans,M,N,K,-1))
         ENDIF
      ENDIF
!
      IF ( nb<nbmin .OR. nb>=K ) THEN
!
!        Use unblocked code
!
         CALL CUNML2(Side,Trans,M,N,K,A,Lda,Tau,C,Ldc,Work,iinfo)
      ELSE
!
!        Use blocked code
!
         iwt = 1 + nw*nb
         IF ( (left .AND. notran) .OR. (.NOT.left .AND. .NOT.notran) )  &
     &        THEN
            i1 = 1
            i2 = K
            i3 = nb
         ELSE
            i1 = ((K-1)/nb)*nb + 1
            i2 = 1
            i3 = -nb
         ENDIF
!
         IF ( left ) THEN
            ni = N
            jc = 1
         ELSE
            mi = M
            ic = 1
         ENDIF
!
         IF ( notran ) THEN
            transt = 'C'
         ELSE
            transt = 'N'
         ENDIF
!
         DO i = i1 , i2 , i3
            ib = MIN(nb,K-i+1)
!
!           Form the triangular factor of the block reflector
!           H = H(i) H(i+1) . . . H(i+ib-1)
!
            CALL CLARFT('Forward','Rowwise',nq-i+1,ib,A(i,i),Lda,Tau(i),&
     &                  Work(iwt),LDT)
            IF ( left ) THEN
!
!              H or H**H is applied to C(i:m,1:n)
!
               mi = M - i + 1
               ic = i
            ELSE
!
!              H or H**H is applied to C(1:m,i:n)
!
               ni = N - i + 1
               jc = i
            ENDIF
!
!           Apply H or H**H
!
            CALL CLARFB(Side,transt,'Forward','Rowwise',mi,ni,ib,A(i,i),&
     &                  Lda,Work(iwt),LDT,C(ic,jc),Ldc,Work,ldwork)
         ENDDO
      ENDIF
      Work(1) = lwkopt
!
!     End of CUNMLQ
!
      END SUBROUTINE CUNMLQ
