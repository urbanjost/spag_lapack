!*==zunmhr.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b ZUNMHR
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZUNMHR + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zunmhr.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zunmhr.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zunmhr.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZUNMHR( SIDE, TRANS, M, N, ILO, IHI, A, LDA, TAU, C,
!                          LDC, WORK, LWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          SIDE, TRANS
!       INTEGER            IHI, ILO, INFO, LDA, LDC, LWORK, M, N
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
!> ZUNMHR overwrites the general complex M-by-N matrix C with
!>
!>                 SIDE = 'L'     SIDE = 'R'
!> TRANS = 'N':      Q * C          C * Q
!> TRANS = 'C':      Q**H * C       C * Q**H
!>
!> where Q is a complex unitary matrix of order nq, with nq = m if
!> SIDE = 'L' and nq = n if SIDE = 'R'. Q is defined as the product of
!> IHI-ILO elementary reflectors, as returned by ZGEHRD:
!>
!> Q = H(ilo) H(ilo+1) . . . H(ihi-1).
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
!> \param[in] ILO
!> \verbatim
!>          ILO is INTEGER
!> \endverbatim
!>
!> \param[in] IHI
!> \verbatim
!>          IHI is INTEGER
!>
!>          ILO and IHI must have the same values as in the previous call
!>          of ZGEHRD. Q is equal to the unit matrix except in the
!>          submatrix Q(ilo+1:ihi,ilo+1:ihi).
!>          If SIDE = 'L', then 1 <= ILO <= IHI <= M, if M > 0, and
!>          ILO = 1 and IHI = 0, if M = 0;
!>          if SIDE = 'R', then 1 <= ILO <= IHI <= N, if N > 0, and
!>          ILO = 1 and IHI = 0, if N = 0.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension
!>                               (LDA,M) if SIDE = 'L'
!>                               (LDA,N) if SIDE = 'R'
!>          The vectors which define the elementary reflectors, as
!>          returned by ZGEHRD.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.
!>          LDA >= max(1,M) if SIDE = 'L'; LDA >= max(1,N) if SIDE = 'R'.
!> \endverbatim
!>
!> \param[in] TAU
!> \verbatim
!>          TAU is COMPLEX*16 array, dimension
!>                               (M-1) if SIDE = 'L'
!>                               (N-1) if SIDE = 'R'
!>          TAU(i) must contain the scalar factor of the elementary
!>          reflector H(i), as returned by ZGEHRD.
!> \endverbatim
!>
!> \param[in,out] C
!> \verbatim
!>          C is COMPLEX*16 array, dimension (LDC,N)
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
!>          WORK is COMPLEX*16 array, dimension (MAX(1,LWORK))
!>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The dimension of the array WORK.
!>          If SIDE = 'L', LWORK >= max(1,N);
!>          if SIDE = 'R', LWORK >= max(1,M).
!>          For optimum performance LWORK >= N*NB if SIDE = 'L', and
!>          LWORK >= M*NB if SIDE = 'R', where NB is the optimal
!>          blocksize.
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
!> \ingroup complex16OTHERcomputational
!
!  =====================================================================
      SUBROUTINE ZUNMHR(Side,Trans,M,N,Ilo,Ihi,A,Lda,Tau,C,Ldc,Work,    &
     &                  Lwork,Info)
      USE F77KINDS                        
      USE S_ILAENV
      USE S_LSAME
      USE S_XERBLA
      USE S_ZUNMQR
      IMPLICIT NONE
!*--ZUNMHR187
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Side
      CHARACTER :: Trans
      INTEGER :: M
      INTEGER :: N
      INTEGER , INTENT(IN) :: Ilo
      INTEGER , INTENT(IN) :: Ihi
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(*) :: Tau
      COMPLEX(CX16KIND) , DIMENSION(Ldc,*) :: C
      INTEGER :: Ldc
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      INTEGER :: i1 , i2 , iinfo , lwkopt , mi , nb , nh , ni , nq , nw
      LOGICAL :: left , lquery
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
      nh = Ihi - Ilo
      left = LSAME(Side,'L')
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
      ELSEIF ( .NOT.LSAME(Trans,'N') .AND. .NOT.LSAME(Trans,'C') ) THEN
         Info = -2
      ELSEIF ( M<0 ) THEN
         Info = -3
      ELSEIF ( N<0 ) THEN
         Info = -4
      ELSEIF ( Ilo<1 .OR. Ilo>MAX(1,nq) ) THEN
         Info = -5
      ELSEIF ( Ihi<MIN(Ilo,nq) .OR. Ihi>nq ) THEN
         Info = -6
      ELSEIF ( Lda<MAX(1,nq) ) THEN
         Info = -8
      ELSEIF ( Ldc<MAX(1,M) ) THEN
         Info = -11
      ELSEIF ( Lwork<MAX(1,nw) .AND. .NOT.lquery ) THEN
         Info = -13
      ENDIF
!
      IF ( Info==0 ) THEN
         IF ( left ) THEN
            nb = ILAENV(1,'ZUNMQR',Side//Trans,nh,N,nh,-1)
         ELSE
            nb = ILAENV(1,'ZUNMQR',Side//Trans,M,nh,nh,-1)
         ENDIF
         lwkopt = MAX(1,nw)*nb
         Work(1) = lwkopt
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('ZUNMHR',-Info)
         RETURN
      ELSEIF ( lquery ) THEN
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( M==0 .OR. N==0 .OR. nh==0 ) THEN
         Work(1) = 1
         RETURN
      ENDIF
!
      IF ( left ) THEN
         mi = nh
         ni = N
         i1 = Ilo + 1
         i2 = 1
      ELSE
         mi = M
         ni = nh
         i1 = 1
         i2 = Ilo + 1
      ENDIF
!
      CALL ZUNMQR(Side,Trans,mi,ni,nh,A(Ilo+1,Ilo),Lda,Tau(Ilo),C(i1,i2)&
     &            ,Ldc,Work,Lwork,iinfo)
!
      Work(1) = lwkopt
!
!     End of ZUNMHR
!
      END SUBROUTINE ZUNMHR
