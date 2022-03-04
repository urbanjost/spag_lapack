!*==zunmtr.f90  processed by SPAG 7.51RB at 20:09 on  3 Mar 2022
!> \brief \b ZUNMTR
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZUNMTR + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zunmtr.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zunmtr.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zunmtr.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZUNMTR( SIDE, UPLO, TRANS, M, N, A, LDA, TAU, C, LDC,
!                          WORK, LWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          SIDE, TRANS, UPLO
!       INTEGER            INFO, LDA, LDC, LWORK, M, N
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
!> ZUNMTR overwrites the general complex M-by-N matrix C with
!>
!>                 SIDE = 'L'     SIDE = 'R'
!> TRANS = 'N':      Q * C          C * Q
!> TRANS = 'C':      Q**H * C       C * Q**H
!>
!> where Q is a complex unitary matrix of order nq, with nq = m if
!> SIDE = 'L' and nq = n if SIDE = 'R'. Q is defined as the product of
!> nq-1 elementary reflectors, as returned by ZHETRD:
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
!>          = 'U': Upper triangle of A contains elementary reflectors
!>                 from ZHETRD;
!>          = 'L': Lower triangle of A contains elementary reflectors
!>                 from ZHETRD.
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
!> \param[in] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension
!>                               (LDA,M) if SIDE = 'L'
!>                               (LDA,N) if SIDE = 'R'
!>          The vectors which define the elementary reflectors, as
!>          returned by ZHETRD.
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
!>          reflector H(i), as returned by ZHETRD.
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
!>          LWORK >=M*NB if SIDE = 'R', where NB is the optimal
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
      SUBROUTINE ZUNMTR(Side,Uplo,Trans,M,N,A,Lda,Tau,C,Ldc,Work,Lwork, &
     &                  Info)
      IMPLICIT NONE
!*--ZUNMTR175
!
!  -- LAPACK computational routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Side , Trans , Uplo
      INTEGER Info , Lda , Ldc , Lwork , M , N
!     ..
!     .. Array Arguments ..
      COMPLEX*16 A(Lda,*) , C(Ldc,*) , Tau(*) , Work(*)
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
      LOGICAL left , lquery , upper
      INTEGER i1 , i2 , iinfo , lwkopt , mi , nb , ni , nq , nw
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      INTEGER ILAENV
      EXTERNAL LSAME , ILAENV
!     ..
!     .. External Subroutines ..
      EXTERNAL XERBLA , ZUNMQL , ZUNMQR
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      Info = 0
      left = LSAME(Side,'L')
      upper = LSAME(Uplo,'U')
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
      ELSEIF ( .NOT.upper .AND. .NOT.LSAME(Uplo,'L') ) THEN
         Info = -2
      ELSEIF ( .NOT.LSAME(Trans,'N') .AND. .NOT.LSAME(Trans,'C') ) THEN
         Info = -3
      ELSEIF ( M<0 ) THEN
         Info = -4
      ELSEIF ( N<0 ) THEN
         Info = -5
      ELSEIF ( Lda<MAX(1,nq) ) THEN
         Info = -7
      ELSEIF ( Ldc<MAX(1,M) ) THEN
         Info = -10
      ELSEIF ( Lwork<MAX(1,nw) .AND. .NOT.lquery ) THEN
         Info = -12
      ENDIF
!
      IF ( Info==0 ) THEN
         IF ( upper ) THEN
            IF ( left ) THEN
               nb = ILAENV(1,'ZUNMQL',Side//Trans,M-1,N,M-1,-1)
            ELSE
               nb = ILAENV(1,'ZUNMQL',Side//Trans,M,N-1,N-1,-1)
            ENDIF
         ELSEIF ( left ) THEN
            nb = ILAENV(1,'ZUNMQR',Side//Trans,M-1,N,M-1,-1)
         ELSE
            nb = ILAENV(1,'ZUNMQR',Side//Trans,M,N-1,N-1,-1)
         ENDIF
         lwkopt = MAX(1,nw)*nb
         Work(1) = lwkopt
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('ZUNMTR',-Info)
         RETURN
      ELSEIF ( lquery ) THEN
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( M==0 .OR. N==0 .OR. nq==1 ) THEN
         Work(1) = 1
         RETURN
      ENDIF
!
      IF ( left ) THEN
         mi = M - 1
         ni = N
      ELSE
         mi = M
         ni = N - 1
      ENDIF
!
      IF ( upper ) THEN
!
!        Q was determined by a call to ZHETRD with UPLO = 'U'
!
         CALL ZUNMQL(Side,Trans,mi,ni,nq-1,A(1,2),Lda,Tau,C,Ldc,Work,   &
     &               Lwork,iinfo)
      ELSE
!
!        Q was determined by a call to ZHETRD with UPLO = 'L'
!
         IF ( left ) THEN
            i1 = 2
            i2 = 1
         ELSE
            i1 = 1
            i2 = 2
         ENDIF
         CALL ZUNMQR(Side,Trans,mi,ni,nq-1,A(2,1),Lda,Tau,C(i1,i2),Ldc, &
     &               Work,Lwork,iinfo)
      ENDIF
      Work(1) = lwkopt
!
!     End of ZUNMTR
!
      END SUBROUTINE ZUNMTR
