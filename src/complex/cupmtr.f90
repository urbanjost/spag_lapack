!*==cupmtr.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b CUPMTR
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CUPMTR + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cupmtr.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cupmtr.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cupmtr.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CUPMTR( SIDE, UPLO, TRANS, M, N, AP, TAU, C, LDC, WORK,
!                          INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          SIDE, TRANS, UPLO
!       INTEGER            INFO, LDC, M, N
!       ..
!       .. Array Arguments ..
!       COMPLEX            AP( * ), C( LDC, * ), TAU( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CUPMTR overwrites the general complex M-by-N matrix C with
!>
!>                 SIDE = 'L'     SIDE = 'R'
!> TRANS = 'N':      Q * C          C * Q
!> TRANS = 'C':      Q**H * C       C * Q**H
!>
!> where Q is a complex unitary matrix of order nq, with nq = m if
!> SIDE = 'L' and nq = n if SIDE = 'R'. Q is defined as the product of
!> nq-1 elementary reflectors, as returned by CHPTRD using packed
!> storage:
!>
!> if UPLO = 'U', Q = H(nq-1) . . . H(2) H(1);
!>
!> if UPLO = 'L', Q = H(1) H(2) . . . H(nq-1).
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
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          = 'U': Upper triangular packed storage used in previous
!>                 call to CHPTRD;
!>          = 'L': Lower triangular packed storage used in previous
!>                 call to CHPTRD.
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
!> \param[in] AP
!> \verbatim
!>          AP is COMPLEX array, dimension
!>                               (M*(M+1)/2) if SIDE = 'L'
!>                               (N*(N+1)/2) if SIDE = 'R'
!>          The vectors which define the elementary reflectors, as
!>          returned by CHPTRD.  AP is modified by the routine but
!>          restored on exit.
!> \endverbatim
!>
!> \param[in] TAU
!> \verbatim
!>          TAU is COMPLEX array, dimension (M-1) if SIDE = 'L'
!>                                     or (N-1) if SIDE = 'R'
!>          TAU(i) must contain the scalar factor of the elementary
!>          reflector H(i), as returned by CHPTRD.
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
!>          WORK is COMPLEX array, dimension
!>                                   (N) if SIDE = 'L'
!>                                   (M) if SIDE = 'R'
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
      SUBROUTINE CUPMTR(Side,Uplo,Trans,M,N,Ap,Tau,C,Ldc,Work,Info)
      IMPLICIT NONE
!*--CUPMTR153
!
!  -- LAPACK computational routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Side , Trans , Uplo
      INTEGER Info , Ldc , M , N
!     ..
!     .. Array Arguments ..
      COMPLEX Ap(*) , C(Ldc,*) , Tau(*) , Work(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      COMPLEX ONE
      PARAMETER (ONE=(1.0E+0,0.0E+0))
!     ..
!     .. Local Scalars ..
      LOGICAL forwrd , left , notran , upper
      INTEGER i , i1 , i2 , i3 , ic , ii , jc , mi , ni , nq
      COMPLEX aii , taui
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL CLARF , XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC CONJG , MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      Info = 0
      left = LSAME(Side,'L')
      notran = LSAME(Trans,'N')
      upper = LSAME(Uplo,'U')
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
      ELSEIF ( .NOT.upper .AND. .NOT.LSAME(Uplo,'L') ) THEN
         Info = -2
      ELSEIF ( .NOT.notran .AND. .NOT.LSAME(Trans,'C') ) THEN
         Info = -3
      ELSEIF ( M<0 ) THEN
         Info = -4
      ELSEIF ( N<0 ) THEN
         Info = -5
      ELSEIF ( Ldc<MAX(1,M) ) THEN
         Info = -9
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('CUPMTR',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( M==0 .OR. N==0 ) RETURN
!
      IF ( upper ) THEN
!
!        Q was determined by a call to CHPTRD with UPLO = 'U'
!
         forwrd = (left .AND. notran) .OR. (.NOT.left .AND. .NOT.notran)
!
         IF ( forwrd ) THEN
            i1 = 1
            i2 = nq - 1
            i3 = 1
            ii = 2
         ELSE
            i1 = nq - 1
            i2 = 1
            i3 = -1
            ii = nq*(nq+1)/2 - 1
         ENDIF
!
         IF ( left ) THEN
            ni = N
         ELSE
            mi = M
         ENDIF
!
         DO i = i1 , i2 , i3
            IF ( left ) THEN
!
!              H(i) or H(i)**H is applied to C(1:i,1:n)
!
               mi = i
            ELSE
!
!              H(i) or H(i)**H is applied to C(1:m,1:i)
!
               ni = i
            ENDIF
!
!           Apply H(i) or H(i)**H
!
            IF ( notran ) THEN
               taui = Tau(i)
            ELSE
               taui = CONJG(Tau(i))
            ENDIF
            aii = Ap(ii)
            Ap(ii) = ONE
            CALL CLARF(Side,mi,ni,Ap(ii-i+1),1,taui,C,Ldc,Work)
            Ap(ii) = aii
!
            IF ( forwrd ) THEN
               ii = ii + i + 2
            ELSE
               ii = ii - i - 1
            ENDIF
         ENDDO
      ELSE
!
!        Q was determined by a call to CHPTRD with UPLO = 'L'.
!
         forwrd = (left .AND. .NOT.notran) .OR. (.NOT.left .AND. notran)
!
         IF ( forwrd ) THEN
            i1 = 1
            i2 = nq - 1
            i3 = 1
            ii = 2
         ELSE
            i1 = nq - 1
            i2 = 1
            i3 = -1
            ii = nq*(nq+1)/2 - 1
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
         DO i = i1 , i2 , i3
            aii = Ap(ii)
            Ap(ii) = ONE
            IF ( left ) THEN
!
!              H(i) or H(i)**H is applied to C(i+1:m,1:n)
!
               mi = M - i
               ic = i + 1
            ELSE
!
!              H(i) or H(i)**H is applied to C(1:m,i+1:n)
!
               ni = N - i
               jc = i + 1
            ENDIF
!
!           Apply H(i) or H(i)**H
!
            IF ( notran ) THEN
               taui = Tau(i)
            ELSE
               taui = CONJG(Tau(i))
            ENDIF
            CALL CLARF(Side,mi,ni,Ap(ii),1,taui,C(ic,jc),Ldc,Work)
            Ap(ii) = aii
!
            IF ( forwrd ) THEN
               ii = ii + nq - i + 1
            ELSE
               ii = ii - nq + i - 2
            ENDIF
         ENDDO
      ENDIF
!
!     End of CUPMTR
!
      END SUBROUTINE CUPMTR
