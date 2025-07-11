!*==dgesvdq.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief <b> DGESVDQ computes the singular value decomposition (SVD) with a QR-Preconditioned QR SVD Method for GE matrices</b>
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DGESVDQ + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgesvdq.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgesvdq.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgesvdq.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!      SUBROUTINE DGESVDQ( JOBA, JOBP, JOBR, JOBU, JOBV, M, N, A, LDA,
!                          S, U, LDU, V, LDV, NUMRANK, IWORK, LIWORK,
!                          WORK, LWORK, RWORK, LRWORK, INFO )
!
!     .. Scalar Arguments ..
!      IMPLICIT    NONE
!      CHARACTER   JOBA, JOBP, JOBR, JOBU, JOBV
!      INTEGER     M, N, LDA, LDU, LDV, NUMRANK, LIWORK, LWORK, LRWORK,
!                  INFO
!     ..
!     .. Array Arguments ..
!      DOUBLE PRECISION  A( LDA, * ), U( LDU, * ), V( LDV, * ), WORK( * )
!      DOUBLE PRECISION  S( * ), RWORK( * )
!      INTEGER     IWORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DGESVDQ computes the singular value decomposition (SVD) of a real
!> M-by-N matrix A, where M >= N. The SVD of A is written as
!>                                    [++]   [xx]   [x0]   [xx]
!>              A = U * SIGMA * V^*,  [++] = [xx] * [ox] * [xx]
!>                                    [++]   [xx]
!> where SIGMA is an N-by-N diagonal matrix, U is an M-by-N orthonormal
!> matrix, and V is an N-by-N orthogonal matrix. The diagonal elements
!> of SIGMA are the singular values of A. The columns of U and V are the
!> left and the right singular vectors of A, respectively.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] JOBA
!> \verbatim
!>  JOBA is CHARACTER*1
!>  Specifies the level of accuracy in the computed SVD
!>  = 'A' The requested accuracy corresponds to having the backward
!>        error bounded by || delta A ||_F <= f(m,n) * EPS * || A ||_F,
!>        where EPS = DLAMCH('Epsilon'). This authorises DGESVDQ to
!>        truncate the computed triangular factor in a rank revealing
!>        QR factorization whenever the truncated part is below the
!>        threshold of the order of EPS * ||A||_F. This is aggressive
!>        truncation level.
!>  = 'M' Similarly as with 'A', but the truncation is more gentle: it
!>        is allowed only when there is a drop on the diagonal of the
!>        triangular factor in the QR factorization. This is medium
!>        truncation level.
!>  = 'H' High accuracy requested. No numerical rank determination based
!>        on the rank revealing QR factorization is attempted.
!>  = 'E' Same as 'H', and in addition the condition number of column
!>        scaled A is estimated and returned in  RWORK(1).
!>        N^(-1/4)*RWORK(1) <= ||pinv(A_scaled)||_2 <= N^(1/4)*RWORK(1)
!> \endverbatim
!>
!> \param[in] JOBP
!> \verbatim
!>  JOBP is CHARACTER*1
!>  = 'P' The rows of A are ordered in decreasing order with respect to
!>        ||A(i,:)||_\infty. This enhances numerical accuracy at the cost
!>        of extra data movement. Recommended for numerical robustness.
!>  = 'N' No row pivoting.
!> \endverbatim
!>
!> \param[in] JOBR
!> \verbatim
!>          JOBR is CHARACTER*1
!>          = 'T' After the initial pivoted QR factorization, DGESVD is applied to
!>          the transposed R**T of the computed triangular factor R. This involves
!>          some extra data movement (matrix transpositions). Useful for
!>          experiments, research and development.
!>          = 'N' The triangular factor R is given as input to DGESVD. This may be
!>          preferred as it involves less data movement.
!> \endverbatim
!>
!> \param[in] JOBU
!> \verbatim
!>          JOBU is CHARACTER*1
!>          = 'A' All M left singular vectors are computed and returned in the
!>          matrix U. See the description of U.
!>          = 'S' or 'U' N = min(M,N) left singular vectors are computed and returned
!>          in the matrix U. See the description of U.
!>          = 'R' Numerical rank NUMRANK is determined and only NUMRANK left singular
!>          vectors are computed and returned in the matrix U.
!>          = 'F' The N left singular vectors are returned in factored form as the
!>          product of the Q factor from the initial QR factorization and the
!>          N left singular vectors of (R**T , 0)**T. If row pivoting is used,
!>          then the necessary information on the row pivoting is stored in
!>          IWORK(N+1:N+M-1).
!>          = 'N' The left singular vectors are not computed.
!> \endverbatim
!>
!> \param[in] JOBV
!> \verbatim
!>          JOBV is CHARACTER*1
!>          = 'A', 'V' All N right singular vectors are computed and returned in
!>          the matrix V.
!>          = 'R' Numerical rank NUMRANK is determined and only NUMRANK right singular
!>          vectors are computed and returned in the matrix V. This option is
!>          allowed only if JOBU = 'R' or JOBU = 'N'; otherwise it is illegal.
!>          = 'N' The right singular vectors are not computed.
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the input matrix A.  M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the input matrix A.  M >= N >= 0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is DOUBLE PRECISION array of dimensions LDA x N
!>          On entry, the input matrix A.
!>          On exit, if JOBU .NE. 'N' or JOBV .NE. 'N', the lower triangle of A contains
!>          the Householder vectors as stored by DGEQP3. If JOBU = 'F', these Householder
!>          vectors together with WORK(1:N) can be used to restore the Q factors from
!>          the initial pivoted QR factorization of A. See the description of U.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER.
!>          The leading dimension of the array A.  LDA >= max(1,M).
!> \endverbatim
!>
!> \param[out] S
!> \verbatim
!>          S is DOUBLE PRECISION array of dimension N.
!>          The singular values of A, ordered so that S(i) >= S(i+1).
!> \endverbatim
!>
!> \param[out] U
!> \verbatim
!>          U is DOUBLE PRECISION array, dimension
!>          LDU x M if JOBU = 'A'; see the description of LDU. In this case,
!>          on exit, U contains the M left singular vectors.
!>          LDU x N if JOBU = 'S', 'U', 'R' ; see the description of LDU. In this
!>          case, U contains the leading N or the leading NUMRANK left singular vectors.
!>          LDU x N if JOBU = 'F' ; see the description of LDU. In this case U
!>          contains N x N orthogonal matrix that can be used to form the left
!>          singular vectors.
!>          If JOBU = 'N', U is not referenced.
!> \endverbatim
!>
!> \param[in] LDU
!> \verbatim
!>          LDU is INTEGER.
!>          The leading dimension of the array U.
!>          If JOBU = 'A', 'S', 'U', 'R',  LDU >= max(1,M).
!>          If JOBU = 'F',                 LDU >= max(1,N).
!>          Otherwise,                     LDU >= 1.
!> \endverbatim
!>
!> \param[out] V
!> \verbatim
!>          V is DOUBLE PRECISION array, dimension
!>          LDV x N if JOBV = 'A', 'V', 'R' or if JOBA = 'E' .
!>          If JOBV = 'A', or 'V',  V contains the N-by-N orthogonal matrix  V**T;
!>          If JOBV = 'R', V contains the first NUMRANK rows of V**T (the right
!>          singular vectors, stored rowwise, of the NUMRANK largest singular values).
!>          If JOBV = 'N' and JOBA = 'E', V is used as a workspace.
!>          If JOBV = 'N', and JOBA.NE.'E', V is not referenced.
!> \endverbatim
!>
!> \param[in] LDV
!> \verbatim
!>          LDV is INTEGER
!>          The leading dimension of the array V.
!>          If JOBV = 'A', 'V', 'R',  or JOBA = 'E', LDV >= max(1,N).
!>          Otherwise,                               LDV >= 1.
!> \endverbatim
!>
!> \param[out] NUMRANK
!> \verbatim
!>          NUMRANK is INTEGER
!>          NUMRANK is the numerical rank first determined after the rank
!>          revealing QR factorization, following the strategy specified by the
!>          value of JOBA. If JOBV = 'R' and JOBU = 'R', only NUMRANK
!>          leading singular values and vectors are then requested in the call
!>          of DGESVD. The final value of NUMRANK might be further reduced if
!>          some singular values are computed as zeros.
!> \endverbatim
!>
!> \param[out] IWORK
!> \verbatim
!>          IWORK is INTEGER array, dimension (max(1, LIWORK)).
!>          On exit, IWORK(1:N) contains column pivoting permutation of the
!>          rank revealing QR factorization.
!>          If JOBP = 'P', IWORK(N+1:N+M-1) contains the indices of the sequence
!>          of row swaps used in row pivoting. These can be used to restore the
!>          left singular vectors in the case JOBU = 'F'.
!>
!>          If LIWORK, LWORK, or LRWORK = -1, then on exit, if INFO = 0,
!>          IWORK(1) returns the minimal LIWORK.
!> \endverbatim
!>
!> \param[in] LIWORK
!> \verbatim
!>          LIWORK is INTEGER
!>          The dimension of the array IWORK.
!>          LIWORK >= N + M - 1,     if JOBP = 'P' and JOBA .NE. 'E';
!>          LIWORK >= N              if JOBP = 'N' and JOBA .NE. 'E';
!>          LIWORK >= N + M - 1 + N, if JOBP = 'P' and JOBA = 'E';
!>          LIWORK >= N + N          if JOBP = 'N' and JOBA = 'E'.
!>
!>          If LIWORK = -1, then a workspace query is assumed; the routine
!>          only calculates and returns the optimal and minimal sizes
!>          for the WORK, IWORK, and RWORK arrays, and no error
!>          message related to LWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (max(2, LWORK)), used as a workspace.
!>          On exit, if, on entry, LWORK.NE.-1, WORK(1:N) contains parameters
!>          needed to recover the Q factor from the QR factorization computed by
!>          DGEQP3.
!>
!>          If LIWORK, LWORK, or LRWORK = -1, then on exit, if INFO = 0,
!>          WORK(1) returns the optimal LWORK, and
!>          WORK(2) returns the minimal LWORK.
!> \endverbatim
!>
!> \param[in,out] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The dimension of the array WORK. It is determined as follows:
!>          Let  LWQP3 = 3*N+1,  LWCON = 3*N, and let
!>          LWORQ = { MAX( N, 1 ),  if JOBU = 'R', 'S', or 'U'
!>                  { MAX( M, 1 ),  if JOBU = 'A'
!>          LWSVD = MAX( 5*N, 1 )
!>          LWLQF = MAX( N/2, 1 ), LWSVD2 = MAX( 5*(N/2), 1 ), LWORLQ = MAX( N, 1 ),
!>          LWQRF = MAX( N/2, 1 ), LWORQ2 = MAX( N, 1 )
!>          Then the minimal value of LWORK is:
!>          = MAX( N + LWQP3, LWSVD )        if only the singular values are needed;
!>          = MAX( N + LWQP3, LWCON, LWSVD ) if only the singular values are needed,
!>                                   and a scaled condition estimate requested;
!>
!>          = N + MAX( LWQP3, LWSVD, LWORQ ) if the singular values and the left
!>                                   singular vectors are requested;
!>          = N + MAX( LWQP3, LWCON, LWSVD, LWORQ ) if the singular values and the left
!>                                   singular vectors are requested, and also
!>                                   a scaled condition estimate requested;
!>
!>          = N + MAX( LWQP3, LWSVD )        if the singular values and the right
!>                                   singular vectors are requested;
!>          = N + MAX( LWQP3, LWCON, LWSVD ) if the singular values and the right
!>                                   singular vectors are requested, and also
!>                                   a scaled condition etimate requested;
!>
!>          = N + MAX( LWQP3, LWSVD, LWORQ ) if the full SVD is requested with JOBV = 'R';
!>                                   independent of JOBR;
!>          = N + MAX( LWQP3, LWCON, LWSVD, LWORQ ) if the full SVD is requested,
!>                                   JOBV = 'R' and, also a scaled condition
!>                                   estimate requested; independent of JOBR;
!>          = MAX( N + MAX( LWQP3, LWSVD, LWORQ ),
!>         N + MAX( LWQP3, N/2+LWLQF, N/2+LWSVD2, N/2+LWORLQ, LWORQ) ) if the
!>                         full SVD is requested with JOBV = 'A' or 'V', and
!>                         JOBR ='N'
!>          = MAX( N + MAX( LWQP3, LWCON, LWSVD, LWORQ ),
!>         N + MAX( LWQP3, LWCON, N/2+LWLQF, N/2+LWSVD2, N/2+LWORLQ, LWORQ ) )
!>                         if the full SVD is requested with JOBV = 'A' or 'V', and
!>                         JOBR ='N', and also a scaled condition number estimate
!>                         requested.
!>          = MAX( N + MAX( LWQP3, LWSVD, LWORQ ),
!>         N + MAX( LWQP3, N/2+LWQRF, N/2+LWSVD2, N/2+LWORQ2, LWORQ ) ) if the
!>                         full SVD is requested with JOBV = 'A', 'V', and JOBR ='T'
!>          = MAX( N + MAX( LWQP3, LWCON, LWSVD, LWORQ ),
!>         N + MAX( LWQP3, LWCON, N/2+LWQRF, N/2+LWSVD2, N/2+LWORQ2, LWORQ ) )
!>                         if the full SVD is requested with JOBV = 'A' or 'V', and
!>                         JOBR ='T', and also a scaled condition number estimate
!>                         requested.
!>          Finally, LWORK must be at least two: LWORK = MAX( 2, LWORK ).
!>
!>          If LWORK = -1, then a workspace query is assumed; the routine
!>          only calculates and returns the optimal and minimal sizes
!>          for the WORK, IWORK, and RWORK arrays, and no error
!>          message related to LWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is DOUBLE PRECISION array, dimension (max(1, LRWORK)).
!>          On exit,
!>          1. If JOBA = 'E', RWORK(1) contains an estimate of the condition
!>          number of column scaled A. If A = C * D where D is diagonal and C
!>          has unit columns in the Euclidean norm, then, assuming full column rank,
!>          N^(-1/4) * RWORK(1) <= ||pinv(C)||_2 <= N^(1/4) * RWORK(1).
!>          Otherwise, RWORK(1) = -1.
!>          2. RWORK(2) contains the number of singular values computed as
!>          exact zeros in DGESVD applied to the upper triangular or trapezoidal
!>          R (from the initial QR factorization). In case of early exit (no call to
!>          DGESVD, such as in the case of zero matrix) RWORK(2) = -1.
!>
!>          If LIWORK, LWORK, or LRWORK = -1, then on exit, if INFO = 0,
!>          RWORK(1) returns the minimal LRWORK.
!> \endverbatim
!>
!> \param[in] LRWORK
!> \verbatim
!>          LRWORK is INTEGER.
!>          The dimension of the array RWORK.
!>          If JOBP ='P', then LRWORK >= MAX(2, M).
!>          Otherwise, LRWORK >= 2
!>
!>          If LRWORK = -1, then a workspace query is assumed; the routine
!>          only calculates and returns the optimal and minimal sizes
!>          for the WORK, IWORK, and RWORK arrays, and no error
!>          message related to LWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit.
!>          < 0:  if INFO = -i, the i-th argument had an illegal value.
!>          > 0:  if DBDSQR did not converge, INFO specifies how many superdiagonals
!>          of an intermediate bidiagonal form B (computed in DGESVD) did not
!>          converge to zero.
!> \endverbatim
!
!> \par Further Details:
!  ========================
!>
!> \verbatim
!>
!>   1. The data movement (matrix transpose) is coded using simple nested
!>   DO-loops because BLAS and LAPACK do not provide corresponding subroutines.
!>   Those DO-loops are easily identified in this source code - by the CONTINUE
!>   statements labeled with 11**. In an optimized version of this code, the
!>   nested DO loops should be replaced with calls to an optimized subroutine.
!>   2. This code scales A by 1/SQRT(M) if the largest ABS(A(i,j)) could cause
!>   column norm overflow. This is the minial precaution and it is left to the
!>   SVD routine (CGESVD) to do its own preemptive scaling if potential over-
!>   or underflows are detected. To avoid repeated scanning of the array A,
!>   an optimal implementation would do all necessary scaling before calling
!>   CGESVD and the scaling in CGESVD can be switched off.
!>   3. Other comments related to code optimization are given in comments in the
!>   code, enlosed in [[double brackets]].
!> \endverbatim
!
!> \par Bugs, examples and comments
!  ===========================
!
!> \verbatim
!>  Please report all bugs and send interesting examples and/or comments to
!>  drmac@math.hr. Thank you.
!> \endverbatim
!
!> \par References
!  ===============
!
!> \verbatim
!>  [1] Zlatko Drmac, Algorithm 977: A QR-Preconditioned QR SVD Method for
!>      Computing the SVD with High Accuracy. ACM Trans. Math. Softw.
!>      44(1): 11:1-11:30 (2017)
!>
!>  SIGMA library, xGESVDQ section updated February 2016.
!>  Developed and coded by Zlatko Drmac, Department of Mathematics
!>  University of Zagreb, Croatia, drmac@math.hr
!> \endverbatim
!
!
!> \par Contributors:
!  ==================
!>
!> \verbatim
!> Developed and coded by Zlatko Drmac, Department of Mathematics
!>  University of Zagreb, Croatia, drmac@math.hr
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
!> \date November 2018
!
!> \ingroup doubleGEsing
!
!  =====================================================================
      SUBROUTINE DGESVDQ(Joba,Jobp,Jobr,Jobu,Jobv,M,N,A,Lda,S,U,Ldu,V,  &
     &                   Ldv,Numrank,Iwork,Liwork,Work,Lwork,Rwork,     &
     &                   Lrwork,Info)
!     .. Scalar Arguments ..
      IMPLICIT NONE
!*--DGESVDQ420
      CHARACTER Joba , Jobp , Jobr , Jobu , Jobv
      INTEGER M , N , Lda , Ldu , Ldv , Numrank , Liwork , Lwork ,      &
     &        Lrwork , Info
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION A(Lda,*) , U(Ldu,*) , V(Ldv,*) , Work(*)
      DOUBLE PRECISION S(*) , Rwork(*)
      INTEGER Iwork(*)
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO , ONE
      PARAMETER (ZERO=0.0D0,ONE=1.0D0)
!     .. Local Scalars ..
      INTEGER ierr , iwoff , nr , n1 , optratio , p , q
      INTEGER lwcon , lwqp3 , lwrk_dgelqf , lwrk_dgesvd , lwrk_dgesvd2 ,&
     &        lwrk_dgeqp3 , lwrk_dgeqrf , lwrk_dormlq , lwrk_dormqr ,   &
     &        lwrk_dormqr2 , lwlqf , lwqrf , lwsvd , lwsvd2 , lworq ,   &
     &        lworq2 , lworlq , minwrk , minwrk2 , optwrk , optwrk2 ,   &
     &        iminwrk , rminwrk
      LOGICAL accla , acclm , acclh , ascaled , conda , dntwu , dntwv , &
     &        lquery , lsvc0 , lsvec , rowprm , rsvec , rtrans , wntua ,&
     &        wntuf , wntur , wntus , wntva , wntvr
      DOUBLE PRECISION big , epsln , rtmp , sconda , sfmin
!     .. Local Arrays
      DOUBLE PRECISION rdummy(1)
!     ..
!     .. External Subroutines (BLAS, LAPACK)
      EXTERNAL DGELQF , DGEQP3 , DGEQRF , DGESVD , DLACPY , DLAPMT ,    &
     &         DLASCL , DLASET , DLASWP , DSCAL , DPOCON , DORMLQ ,     &
     &         DORMQR , XERBLA
!     ..
!     .. External Functions (BLAS, LAPACK)
      LOGICAL LSAME
      INTEGER IDAMAX
      DOUBLE PRECISION DLANGE , DNRM2 , DLAMCH
      EXTERNAL DLANGE , LSAME , IDAMAX , DNRM2 , DLAMCH
!     ..
!     .. Intrinsic Functions ..
!
      INTRINSIC ABS , MAX , MIN , DBLE , SQRT
!
!     Test the input arguments
!
      wntus = LSAME(Jobu,'S') .OR. LSAME(Jobu,'U')
      wntur = LSAME(Jobu,'R')
      wntua = LSAME(Jobu,'A')
      wntuf = LSAME(Jobu,'F')
      lsvc0 = wntus .OR. wntur .OR. wntua
      lsvec = lsvc0 .OR. wntuf
      dntwu = LSAME(Jobu,'N')
!
      wntvr = LSAME(Jobv,'R')
      wntva = LSAME(Jobv,'A') .OR. LSAME(Jobv,'V')
      rsvec = wntvr .OR. wntva
      dntwv = LSAME(Jobv,'N')
!
      accla = LSAME(Joba,'A')
      acclm = LSAME(Joba,'M')
      conda = LSAME(Joba,'E')
      acclh = LSAME(Joba,'H') .OR. conda
!
      rowprm = LSAME(Jobp,'P')
      rtrans = LSAME(Jobr,'T')
!
      IF ( rowprm ) THEN
         IF ( conda ) THEN
            iminwrk = MAX(1,N+M-1+N)
         ELSE
            iminwrk = MAX(1,N+M-1)
         ENDIF
         rminwrk = MAX(2,M)
      ELSE
         IF ( conda ) THEN
            iminwrk = MAX(1,N+N)
         ELSE
            iminwrk = MAX(1,N)
         ENDIF
         rminwrk = 2
      ENDIF
      lquery = (Liwork==-1 .OR. Lwork==-1 .OR. Lrwork==-1)
      Info = 0
      IF ( .NOT.(accla .OR. acclm .OR. acclh) ) THEN
         Info = -1
      ELSEIF ( .NOT.(rowprm .OR. LSAME(Jobp,'N')) ) THEN
         Info = -2
      ELSEIF ( .NOT.(rtrans .OR. LSAME(Jobr,'N')) ) THEN
         Info = -3
      ELSEIF ( .NOT.(lsvec .OR. dntwu) ) THEN
         Info = -4
      ELSEIF ( wntur .AND. wntva ) THEN
         Info = -5
      ELSEIF ( .NOT.(rsvec .OR. dntwv) ) THEN
         Info = -5
      ELSEIF ( M<0 ) THEN
         Info = -6
      ELSEIF ( (N<0) .OR. (N>M) ) THEN
         Info = -7
      ELSEIF ( Lda<MAX(1,M) ) THEN
         Info = -9
      ELSEIF ( Ldu<1 .OR. (lsvc0 .AND. Ldu<M) .OR. (wntuf .AND. Ldu<N) )&
     &         THEN
         Info = -12
      ELSEIF ( Ldv<1 .OR. (rsvec .AND. Ldv<N) .OR. (conda .AND. Ldv<N) )&
     &         THEN
         Info = -14
      ELSEIF ( Liwork<iminwrk .AND. .NOT.lquery ) THEN
         Info = -17
      ENDIF
!
!
      IF ( Info==0 ) THEN
!        .. compute the minimal and the optimal workspace lengths
!        [[The expressions for computing the minimal and the optimal
!        values of LWORK are written with a lot of redundancy and
!        can be simplified. However, this detailed form is easier for
!        maintenance and modifications of the code.]]
!
!        .. minimal workspace length for DGEQP3 of an M x N matrix
         lwqp3 = 3*N + 1
!        .. minimal workspace length for DORMQR to build left singular vectors
         IF ( wntus .OR. wntur ) THEN
            lworq = MAX(N,1)
         ELSEIF ( wntua ) THEN
            lworq = MAX(M,1)
         ENDIF
!        .. minimal workspace length for DPOCON of an N x N matrix
         lwcon = 3*N
!        .. DGESVD of an N x N matrix
         lwsvd = MAX(5*N,1)
         IF ( lquery ) THEN
            CALL DGEQP3(M,N,A,Lda,Iwork,rdummy,rdummy,-1,ierr)
            lwrk_dgeqp3 = INT(rdummy(1))
            IF ( wntus .OR. wntur ) THEN
               CALL DORMQR('L','N',M,N,N,A,Lda,rdummy,U,Ldu,rdummy,-1,  &
     &                     ierr)
               lwrk_dormqr = INT(rdummy(1))
            ELSEIF ( wntua ) THEN
               CALL DORMQR('L','N',M,M,N,A,Lda,rdummy,U,Ldu,rdummy,-1,  &
     &                     ierr)
               lwrk_dormqr = INT(rdummy(1))
            ELSE
               lwrk_dormqr = 0
            ENDIF
         ENDIF
         minwrk = 2
         optwrk = 2
         IF ( .NOT.(lsvec .OR. rsvec) ) THEN
!            .. minimal and optimal sizes of the workspace if
!            only the singular values are requested
            IF ( conda ) THEN
               minwrk = MAX(N+lwqp3,lwcon,lwsvd)
            ELSE
               minwrk = MAX(N+lwqp3,lwsvd)
            ENDIF
            IF ( lquery ) THEN
               CALL DGESVD('N','N',N,N,A,Lda,S,U,Ldu,V,Ldv,rdummy,-1,   &
     &                     ierr)
               lwrk_dgesvd = INT(rdummy(1))
               IF ( conda ) THEN
                  optwrk = MAX(N+lwrk_dgeqp3,N+lwcon,lwrk_dgesvd)
               ELSE
                  optwrk = MAX(N+lwrk_dgeqp3,lwrk_dgesvd)
               ENDIF
            ENDIF
         ELSEIF ( lsvec .AND. (.NOT.rsvec) ) THEN
!            .. minimal and optimal sizes of the workspace if the
!            singular values and the left singular vectors are requested
            IF ( conda ) THEN
               minwrk = N + MAX(lwqp3,lwcon,lwsvd,lworq)
            ELSE
               minwrk = N + MAX(lwqp3,lwsvd,lworq)
            ENDIF
            IF ( lquery ) THEN
               IF ( rtrans ) THEN
                  CALL DGESVD('N','O',N,N,A,Lda,S,U,Ldu,V,Ldv,rdummy,-1,&
     &                        ierr)
               ELSE
                  CALL DGESVD('O','N',N,N,A,Lda,S,U,Ldu,V,Ldv,rdummy,-1,&
     &                        ierr)
               ENDIF
               lwrk_dgesvd = INT(rdummy(1))
               IF ( conda ) THEN
                  optwrk = N + MAX(lwrk_dgeqp3,lwcon,lwrk_dgesvd,       &
     &                     lwrk_dormqr)
               ELSE
                  optwrk = N + MAX(lwrk_dgeqp3,lwrk_dgesvd,lwrk_dormqr)
               ENDIF
            ENDIF
         ELSEIF ( rsvec .AND. (.NOT.lsvec) ) THEN
!            .. minimal and optimal sizes of the workspace if the
!            singular values and the right singular vectors are requested
            IF ( conda ) THEN
               minwrk = N + MAX(lwqp3,lwcon,lwsvd)
            ELSE
               minwrk = N + MAX(lwqp3,lwsvd)
            ENDIF
            IF ( lquery ) THEN
               IF ( rtrans ) THEN
                  CALL DGESVD('O','N',N,N,A,Lda,S,U,Ldu,V,Ldv,rdummy,-1,&
     &                        ierr)
               ELSE
                  CALL DGESVD('N','O',N,N,A,Lda,S,U,Ldu,V,Ldv,rdummy,-1,&
     &                        ierr)
               ENDIF
               lwrk_dgesvd = INT(rdummy(1))
               IF ( conda ) THEN
                  optwrk = N + MAX(lwrk_dgeqp3,lwcon,lwrk_dgesvd)
               ELSE
                  optwrk = N + MAX(lwrk_dgeqp3,lwrk_dgesvd)
               ENDIF
            ENDIF
         ELSE
!            .. minimal and optimal sizes of the workspace if the
!            full SVD is requested
            IF ( rtrans ) THEN
               minwrk = MAX(lwqp3,lwsvd,lworq)
               IF ( conda ) minwrk = MAX(minwrk,lwcon)
               minwrk = minwrk + N
               IF ( wntva ) THEN
!                   .. minimal workspace length for N x N/2 DGEQRF
                  lwqrf = MAX(N/2,1)
!                   .. minimal workspace length for N/2 x N/2 DGESVD
                  lwsvd2 = MAX(5*(N/2),1)
                  lworq2 = MAX(N,1)
                  minwrk2 = MAX(lwqp3,N/2+lwqrf,N/2+lwsvd2,N/2+lworq2,  &
     &                      lworq)
                  IF ( conda ) minwrk2 = MAX(minwrk2,lwcon)
                  minwrk2 = N + minwrk2
                  minwrk = MAX(minwrk,minwrk2)
               ENDIF
            ELSE
               minwrk = MAX(lwqp3,lwsvd,lworq)
               IF ( conda ) minwrk = MAX(minwrk,lwcon)
               minwrk = minwrk + N
               IF ( wntva ) THEN
!                   .. minimal workspace length for N/2 x N DGELQF
                  lwlqf = MAX(N/2,1)
                  lwsvd2 = MAX(5*(N/2),1)
                  lworlq = MAX(N,1)
                  minwrk2 = MAX(lwqp3,N/2+lwlqf,N/2+lwsvd2,N/2+lworlq,  &
     &                      lworq)
                  IF ( conda ) minwrk2 = MAX(minwrk2,lwcon)
                  minwrk2 = N + minwrk2
                  minwrk = MAX(minwrk,minwrk2)
               ENDIF
            ENDIF
            IF ( lquery ) THEN
               IF ( rtrans ) THEN
                  CALL DGESVD('O','A',N,N,A,Lda,S,U,Ldu,V,Ldv,rdummy,-1,&
     &                        ierr)
                  lwrk_dgesvd = INT(rdummy(1))
                  optwrk = MAX(lwrk_dgeqp3,lwrk_dgesvd,lwrk_dormqr)
                  IF ( conda ) optwrk = MAX(optwrk,lwcon)
                  optwrk = N + optwrk
                  IF ( wntva ) THEN
                     CALL DGEQRF(N,N/2,U,Ldu,rdummy,rdummy,-1,ierr)
                     lwrk_dgeqrf = INT(rdummy(1))
                     CALL DGESVD('S','O',N/2,N/2,V,Ldv,S,U,Ldu,V,Ldv,   &
     &                           rdummy,-1,ierr)
                     lwrk_dgesvd2 = INT(rdummy(1))
                     CALL DORMQR('R','C',N,N,N/2,U,Ldu,rdummy,V,Ldv,    &
     &                           rdummy,-1,ierr)
                     lwrk_dormqr2 = INT(rdummy(1))
                     optwrk2 = MAX(lwrk_dgeqp3,N/2+lwrk_dgeqrf,         &
     &                         N/2+lwrk_dgesvd2,N/2+lwrk_dormqr2)
                     IF ( conda ) optwrk2 = MAX(optwrk2,lwcon)
                     optwrk2 = N + optwrk2
                     optwrk = MAX(optwrk,optwrk2)
                  ENDIF
               ELSE
                  CALL DGESVD('S','O',N,N,A,Lda,S,U,Ldu,V,Ldv,rdummy,-1,&
     &                        ierr)
                  lwrk_dgesvd = INT(rdummy(1))
                  optwrk = MAX(lwrk_dgeqp3,lwrk_dgesvd,lwrk_dormqr)
                  IF ( conda ) optwrk = MAX(optwrk,lwcon)
                  optwrk = N + optwrk
                  IF ( wntva ) THEN
                     CALL DGELQF(N/2,N,U,Ldu,rdummy,rdummy,-1,ierr)
                     lwrk_dgelqf = INT(rdummy(1))
                     CALL DGESVD('S','O',N/2,N/2,V,Ldv,S,U,Ldu,V,Ldv,   &
     &                           rdummy,-1,ierr)
                     lwrk_dgesvd2 = INT(rdummy(1))
                     CALL DORMLQ('R','N',N,N,N/2,U,Ldu,rdummy,V,Ldv,    &
     &                           rdummy,-1,ierr)
                     lwrk_dormlq = INT(rdummy(1))
                     optwrk2 = MAX(lwrk_dgeqp3,N/2+lwrk_dgelqf,         &
     &                         N/2+lwrk_dgesvd2,N/2+lwrk_dormlq)
                     IF ( conda ) optwrk2 = MAX(optwrk2,lwcon)
                     optwrk2 = N + optwrk2
                     optwrk = MAX(optwrk,optwrk2)
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
!
         minwrk = MAX(2,minwrk)
         optwrk = MAX(2,optwrk)
         IF ( Lwork<minwrk .AND. (.NOT.lquery) ) Info = -19
!
      ENDIF
!
      IF ( Info==0 .AND. Lrwork<rminwrk .AND. .NOT.lquery ) Info = -21
      IF ( Info/=0 ) THEN
         CALL XERBLA('DGESVDQ',-Info)
         RETURN
      ELSEIF ( lquery ) THEN
!
!     Return optimal workspace
!
         Iwork(1) = iminwrk
         Work(1) = optwrk
         Work(2) = minwrk
         Rwork(1) = rminwrk
         RETURN
      ENDIF
!
!     Quick return if the matrix is void.
!
!     .. all output is void.
      IF ( (M==0) .OR. (N==0) ) RETURN
!
      big = DLAMCH('O')
      ascaled = .FALSE.
      iwoff = 1
      IF ( rowprm ) THEN
         iwoff = M
!           .. reordering the rows in decreasing sequence in the
!           ell-infinity norm - this enhances numerical robustness in
!           the case of differently scaled rows.
         DO p = 1 , M
!               RWORK(p) = ABS( A(p,ICAMAX(N,A(p,1),LDA)) )
!               [[DLANGE will return NaN if an entry of the p-th row is Nan]]
            Rwork(p) = DLANGE('M',1,N,A(p,1),Lda,rdummy)
!               .. check for NaN's and Inf's
            IF ( (Rwork(p)/=Rwork(p)) .OR. ((Rwork(p)*ZERO)/=ZERO) )    &
     &           THEN
               Info = -8
               CALL XERBLA('DGESVDQ',-Info)
               RETURN
            ENDIF
         ENDDO
         DO p = 1 , M - 1
            q = IDAMAX(M-p+1,Rwork(p),1) + p - 1
            Iwork(N+p) = q
            IF ( p/=q ) THEN
               rtmp = Rwork(p)
               Rwork(p) = Rwork(q)
               Rwork(q) = rtmp
            ENDIF
         ENDDO
!
         IF ( Rwork(1)==ZERO ) THEN
!              Quick return: A is the M x N zero matrix.
            Numrank = 0
            CALL DLASET('G',N,1,ZERO,ZERO,S,N)
            IF ( wntus ) CALL DLASET('G',M,N,ZERO,ONE,U,Ldu)
            IF ( wntua ) CALL DLASET('G',M,M,ZERO,ONE,U,Ldu)
            IF ( wntva ) CALL DLASET('G',N,N,ZERO,ONE,V,Ldv)
            IF ( wntuf ) THEN
               CALL DLASET('G',N,1,ZERO,ZERO,Work,N)
               CALL DLASET('G',M,N,ZERO,ONE,U,Ldu)
            ENDIF
            DO p = 1 , N
               Iwork(p) = p
            ENDDO
            IF ( rowprm ) THEN
               DO p = N + 1 , N + M - 1
                  Iwork(p) = p - N
               ENDDO
            ENDIF
            IF ( conda ) Rwork(1) = -1
            Rwork(2) = -1
            RETURN
         ENDIF
!
         IF ( Rwork(1)>big/SQRT(DBLE(M)) ) THEN
!               .. to prevent overflow in the QR factorization, scale the
!               matrix by 1/sqrt(M) if too large entry detected
            CALL DLASCL('G',0,0,SQRT(DBLE(M)),ONE,M,N,A,Lda,ierr)
            ascaled = .TRUE.
         ENDIF
         CALL DLASWP(N,A,Lda,1,M-1,Iwork(N+1),1)
      ENDIF
!
!    .. At this stage, preemptive scaling is done only to avoid column
!    norms overflows during the QR factorization. The SVD procedure should
!    have its own scaling to save the singular values from overflows and
!    underflows. That depends on the SVD procedure.
!
      IF ( .NOT.rowprm ) THEN
         rtmp = DLANGE('M',M,N,A,Lda,rdummy)
         IF ( (rtmp/=rtmp) .OR. ((rtmp*ZERO)/=ZERO) ) THEN
            Info = -8
            CALL XERBLA('DGESVDQ',-Info)
            RETURN
         ENDIF
         IF ( rtmp>big/SQRT(DBLE(M)) ) THEN
!             .. to prevent overflow in the QR factorization, scale the
!             matrix by 1/sqrt(M) if too large entry detected
            CALL DLASCL('G',0,0,SQRT(DBLE(M)),ONE,M,N,A,Lda,ierr)
            ascaled = .TRUE.
         ENDIF
      ENDIF
!
!     .. QR factorization with column pivoting
!
!     A * P = Q * [ R ]
!                 [ 0 ]
!
      DO p = 1 , N
!        .. all columns are free columns
         Iwork(p) = 0
      ENDDO
      CALL DGEQP3(M,N,A,Lda,Iwork,Work,Work(N+1),Lwork-N,ierr)
!
!    If the user requested accuracy level allows truncation in the
!    computed upper triangular factor, the matrix R is examined and,
!    if possible, replaced with its leading upper trapezoidal part.
!
      epsln = DLAMCH('E')
      sfmin = DLAMCH('S')
!     SMALL = SFMIN / EPSLN
      nr = N
!
      IF ( accla ) THEN
!
!        Standard absolute error bound suffices. All sigma_i with
!        sigma_i < N*EPS*||A||_F are flushed to zero. This is an
!        aggressive enforcement of lower numerical rank by introducing a
!        backward error of the order of N*EPS*||A||_F.
         nr = 1
         rtmp = SQRT(DBLE(N))*epsln
         DO p = 2 , N
            IF ( ABS(A(p,p))<(rtmp*ABS(A(1,1))) ) EXIT
            nr = nr + 1
         ENDDO
!
      ELSEIF ( acclm ) THEN
!        .. similarly as above, only slightly more gentle (less aggressive).
!        Sudden drop on the diagonal of R is used as the criterion for being
!        close-to-rank-deficient. The threshold is set to EPSLN=DLAMCH('E').
!        [[This can be made more flexible by replacing this hard-coded value
!        with a user specified threshold.]] Also, the values that underflow
!        will be truncated.
         nr = 1
         DO p = 2 , N
            IF ( (ABS(A(p,p))<(epsln*ABS(A(p-1,p-1)))) .OR.             &
     &           (ABS(A(p,p))<sfmin) ) EXIT
            nr = nr + 1
         ENDDO
!
      ELSE
!        .. RRQR not authorized to determine numerical rank except in the
!        obvious case of zero pivots.
!        .. inspect R for exact zeros on the diagonal;
!        R(i,i)=0 => R(i:N,i:N)=0.
         nr = 1
         DO p = 2 , N
            IF ( ABS(A(p,p))==ZERO ) EXIT
            nr = nr + 1
         ENDDO
!
         IF ( conda ) THEN
!           Estimate the scaled condition number of A. Use the fact that it is
!           the same as the scaled condition number of R.
!              .. V is used as workspace
            CALL DLACPY('U',N,N,A,Lda,V,Ldv)
!              Only the leading NR x NR submatrix of the triangular factor
!              is considered. Only if NR=N will this give a reliable error
!              bound. However, even for NR < N, this can be used on an
!              expert level and obtain useful information in the sense of
!              perturbation theory.
            DO p = 1 , nr
               rtmp = DNRM2(p,V(1,p),1)
               CALL DSCAL(p,ONE/rtmp,V(1,p),1)
            ENDDO
            IF ( .NOT.(lsvec .OR. rsvec) ) THEN
               CALL DPOCON('U',nr,V,Ldv,ONE,rtmp,Work,Iwork(N+iwoff),   &
     &                     ierr)
            ELSE
               CALL DPOCON('U',nr,V,Ldv,ONE,rtmp,Work(N+1),             &
     &                     Iwork(N+iwoff),ierr)
            ENDIF
            sconda = ONE/SQRT(rtmp)
!           For NR=N, SCONDA is an estimate of SQRT(||(R^* * R)^(-1)||_1),
!           N^(-1/4) * SCONDA <= ||R^(-1)||_2 <= N^(1/4) * SCONDA
!           See the reference [1] for more details.
         ENDIF
!
      ENDIF
!
      IF ( wntur ) THEN
         n1 = nr
      ELSEIF ( wntus .OR. wntuf ) THEN
         n1 = N
      ELSEIF ( wntua ) THEN
         n1 = M
      ENDIF
!
      IF ( .NOT.(rsvec .OR. lsvec) ) THEN
!.......................................................................
!        .. only the singular values are requested
!.......................................................................
         IF ( rtrans ) THEN
!
!         .. compute the singular values of R**T = [A](1:NR,1:N)**T
!           .. set the lower triangle of [A] to [A](1:NR,1:N)**T and
!           the upper triangle of [A] to zero.
            DO p = 1 , MIN(N,nr)
               DO q = p + 1 , N
                  A(q,p) = A(p,q)
                  IF ( q<=nr ) A(p,q) = ZERO
               ENDDO
            ENDDO
!
            CALL DGESVD('N','N',N,nr,A,Lda,S,U,Ldu,V,Ldv,Work,Lwork,    &
     &                  Info)
!
         ELSE
!
!           .. compute the singular values of R = [A](1:NR,1:N)
!
            IF ( nr>1 ) CALL DLASET('L',nr-1,nr-1,ZERO,ZERO,A(2,1),Lda)
            CALL DGESVD('N','N',nr,N,A,Lda,S,U,Ldu,V,Ldv,Work,Lwork,    &
     &                  Info)
!
         ENDIF
!
      ELSEIF ( lsvec .AND. (.NOT.rsvec) ) THEN
!.......................................................................
!       .. the singular values and the left singular vectors requested
!.......................................................................""""""""
         IF ( rtrans ) THEN
!            .. apply DGESVD to R**T
!            .. copy R**T into [U] and overwrite [U] with the right singular
!            vectors of R
            DO p = 1 , nr
               DO q = p , N
                  U(q,p) = A(p,q)
               ENDDO
            ENDDO
            IF ( nr>1 ) CALL DLASET('U',nr-1,nr-1,ZERO,ZERO,U(1,2),Ldu)
!           .. the left singular vectors not computed, the NR right singular
!           vectors overwrite [U](1:NR,1:NR) as transposed. These
!           will be pre-multiplied by Q to build the left singular vectors of A.
            CALL DGESVD('N','O',N,nr,U,Ldu,S,U,Ldu,U,Ldu,Work(N+1),     &
     &                  Lwork-N,Info)
!
            DO p = 1 , nr
               DO q = p + 1 , nr
                  rtmp = U(q,p)
                  U(q,p) = U(p,q)
                  U(p,q) = rtmp
               ENDDO
            ENDDO
!
         ELSE
!            .. apply DGESVD to R
!            .. copy R into [U] and overwrite [U] with the left singular vectors
            CALL DLACPY('U',nr,N,A,Lda,U,Ldu)
            IF ( nr>1 ) CALL DLASET('L',nr-1,nr-1,ZERO,ZERO,U(2,1),Ldu)
!            .. the right singular vectors not computed, the NR left singular
!            vectors overwrite [U](1:NR,1:NR)
            CALL DGESVD('O','N',nr,N,U,Ldu,S,U,Ldu,V,Ldv,Work(N+1),     &
     &                  Lwork-N,Info)
!               .. now [U](1:NR,1:NR) contains the NR left singular vectors of
!               R. These will be pre-multiplied by Q to build the left singular
!               vectors of A.
         ENDIF
!
!           .. assemble the left singular vector matrix U of dimensions
!              (M x NR) or (M x N) or (M x M).
         IF ( (nr<M) .AND. (.NOT.wntuf) ) THEN
            CALL DLASET('A',M-nr,nr,ZERO,ZERO,U(nr+1,1),Ldu)
            IF ( nr<n1 ) THEN
               CALL DLASET('A',nr,n1-nr,ZERO,ZERO,U(1,nr+1),Ldu)
               CALL DLASET('A',M-nr,n1-nr,ZERO,ONE,U(nr+1,nr+1),Ldu)
            ENDIF
         ENDIF
!
!           The Q matrix from the first QRF is built into the left singular
!           vectors matrix U.
!
         IF ( .NOT.wntuf ) CALL DORMQR('L','N',M,n1,N,A,Lda,Work,U,Ldu, &
     &                                 Work(N+1),Lwork-N,ierr)
         IF ( rowprm .AND. .NOT.wntuf )                                 &
     &        CALL DLASWP(n1,U,Ldu,1,M-1,Iwork(N+1),-1)
!
      ELSEIF ( rsvec .AND. (.NOT.lsvec) ) THEN
!.......................................................................
!       .. the singular values and the right singular vectors requested
!.......................................................................
         IF ( rtrans ) THEN
!            .. apply DGESVD to R**T
!            .. copy R**T into V and overwrite V with the left singular vectors
            DO p = 1 , nr
               DO q = p , N
                  V(q,p) = (A(p,q))
               ENDDO
            ENDDO
            IF ( nr>1 ) CALL DLASET('U',nr-1,nr-1,ZERO,ZERO,V(1,2),Ldv)
!           .. the left singular vectors of R**T overwrite V, the right singular
!           vectors not computed
            IF ( wntvr .OR. (nr==N) ) THEN
               CALL DGESVD('O','N',N,nr,V,Ldv,S,U,Ldu,U,Ldu,Work(N+1),  &
     &                     Lwork-N,Info)
!
               DO p = 1 , nr
                  DO q = p + 1 , nr
                     rtmp = V(q,p)
                     V(q,p) = V(p,q)
                     V(p,q) = rtmp
                  ENDDO
               ENDDO
!
               IF ( nr<N ) THEN
                  DO p = 1 , nr
                     DO q = nr + 1 , N
                        V(p,q) = V(q,p)
                     ENDDO
                  ENDDO
               ENDIF
               CALL DLAPMT(.FALSE.,nr,N,V,Ldv,Iwork)
            ELSE
!               .. need all N right singular vectors and NR < N
!               [!] This is simple implementation that augments [V](1:N,1:NR)
!               by padding a zero block. In the case NR << N, a more efficient
!               way is to first use the QR factorization. For more details
!               how to implement this, see the " FULL SVD " branch.
               CALL DLASET('G',N,N-nr,ZERO,ZERO,V(1,nr+1),Ldv)
               CALL DGESVD('O','N',N,N,V,Ldv,S,U,Ldu,U,Ldu,Work(N+1),   &
     &                     Lwork-N,Info)
!
               DO p = 1 , N
                  DO q = p + 1 , N
                     rtmp = V(q,p)
                     V(q,p) = V(p,q)
                     V(p,q) = rtmp
                  ENDDO
               ENDDO
               CALL DLAPMT(.FALSE.,N,N,V,Ldv,Iwork)
            ENDIF
!
         ELSE
!            .. aply DGESVD to R
!            .. copy R into V and overwrite V with the right singular vectors
            CALL DLACPY('U',nr,N,A,Lda,V,Ldv)
            IF ( nr>1 ) CALL DLASET('L',nr-1,nr-1,ZERO,ZERO,V(2,1),Ldv)
!            .. the right singular vectors overwrite V, the NR left singular
!            vectors stored in U(1:NR,1:NR)
            IF ( wntvr .OR. (nr==N) ) THEN
               CALL DGESVD('N','O',nr,N,V,Ldv,S,U,Ldu,V,Ldv,Work(N+1),  &
     &                     Lwork-N,Info)
               CALL DLAPMT(.FALSE.,nr,N,V,Ldv,Iwork)
!               .. now [V](1:NR,1:N) contains V(1:N,1:NR)**T
            ELSE
!               .. need all N right singular vectors and NR < N
!               [!] This is simple implementation that augments [V](1:NR,1:N)
!               by padding a zero block. In the case NR << N, a more efficient
!               way is to first use the LQ factorization. For more details
!               how to implement this, see the " FULL SVD " branch.
               CALL DLASET('G',N-nr,N,ZERO,ZERO,V(nr+1,1),Ldv)
               CALL DGESVD('N','O',N,N,V,Ldv,S,U,Ldu,V,Ldv,Work(N+1),   &
     &                     Lwork-N,Info)
               CALL DLAPMT(.FALSE.,N,N,V,Ldv,Iwork)
            ENDIF
!            .. now [V] contains the transposed matrix of the right singular
!            vectors of A.
         ENDIF
!
      ELSE
!.......................................................................
!       .. FULL SVD requested
!.......................................................................
         IF ( rtrans ) THEN
!
!            .. apply DGESVD to R**T [[this option is left for R&D&T]]
!
            IF ( wntvr .OR. (nr==N) ) THEN
!            .. copy R**T into [V] and overwrite [V] with the left singular
!            vectors of R**T
               DO p = 1 , nr
                  DO q = p , N
                     V(q,p) = A(p,q)
                  ENDDO
               ENDDO
               IF ( nr>1 ) CALL DLASET('U',nr-1,nr-1,ZERO,ZERO,V(1,2),  &
     &                                 Ldv)
!
!           .. the left singular vectors of R**T overwrite [V], the NR right
!           singular vectors of R**T stored in [U](1:NR,1:NR) as transposed
               CALL DGESVD('O','A',N,nr,V,Ldv,S,V,Ldv,U,Ldu,Work(N+1),  &
     &                     Lwork-N,Info)
!              .. assemble V
               DO p = 1 , nr
                  DO q = p + 1 , nr
                     rtmp = V(q,p)
                     V(q,p) = V(p,q)
                     V(p,q) = rtmp
                  ENDDO
               ENDDO
               IF ( nr<N ) THEN
                  DO p = 1 , nr
                     DO q = nr + 1 , N
                        V(p,q) = V(q,p)
                     ENDDO
                  ENDDO
               ENDIF
               CALL DLAPMT(.FALSE.,nr,N,V,Ldv,Iwork)
!
               DO p = 1 , nr
                  DO q = p + 1 , nr
                     rtmp = U(q,p)
                     U(q,p) = U(p,q)
                     U(p,q) = rtmp
                  ENDDO
               ENDDO
!
               IF ( (nr<M) .AND. .NOT.(wntuf) ) THEN
                  CALL DLASET('A',M-nr,nr,ZERO,ZERO,U(nr+1,1),Ldu)
                  IF ( nr<n1 ) THEN
                     CALL DLASET('A',nr,n1-nr,ZERO,ZERO,U(1,nr+1),Ldu)
                     CALL DLASET('A',M-nr,n1-nr,ZERO,ONE,U(nr+1,nr+1),  &
     &                           Ldu)
                  ENDIF
               ENDIF
!
            ELSE
!               .. need all N right singular vectors and NR < N
!            .. copy R**T into [V] and overwrite [V] with the left singular
!            vectors of R**T
!               [[The optimal ratio N/NR for using QRF instead of padding
!                 with zeros. Here hard coded to 2; it must be at least
!                 two due to work space constraints.]]
!               OPTRATIO = ILAENV(6, 'DGESVD', 'S' // 'O', NR,N,0,0)
!               OPTRATIO = MAX( OPTRATIO, 2 )
               optratio = 2
               IF ( optratio*nr>N ) THEN
                  DO p = 1 , nr
                     DO q = p , N
                        V(q,p) = A(p,q)
                     ENDDO
                  ENDDO
                  IF ( nr>1 ) CALL DLASET('U',nr-1,nr-1,ZERO,ZERO,V(1,2)&
     &                 ,Ldv)
!
                  CALL DLASET('A',N,N-nr,ZERO,ZERO,V(1,nr+1),Ldv)
                  CALL DGESVD('O','A',N,N,V,Ldv,S,V,Ldv,U,Ldu,Work(N+1),&
     &                        Lwork-N,Info)
!
                  DO p = 1 , N
                     DO q = p + 1 , N
                        rtmp = V(q,p)
                        V(q,p) = V(p,q)
                        V(p,q) = rtmp
                     ENDDO
                  ENDDO
                  CALL DLAPMT(.FALSE.,N,N,V,Ldv,Iwork)
!              .. assemble the left singular vector matrix U of dimensions
!              (M x N1), i.e. (M x N) or (M x M).
!
                  DO p = 1 , N
                     DO q = p + 1 , N
                        rtmp = U(q,p)
                        U(q,p) = U(p,q)
                        U(p,q) = rtmp
                     ENDDO
                  ENDDO
!
                  IF ( (N<M) .AND. .NOT.(wntuf) ) THEN
                     CALL DLASET('A',M-N,N,ZERO,ZERO,U(N+1,1),Ldu)
                     IF ( N<n1 ) THEN
                        CALL DLASET('A',N,n1-N,ZERO,ZERO,U(1,N+1),Ldu)
                        CALL DLASET('A',M-N,n1-N,ZERO,ONE,U(N+1,N+1),   &
     &                              Ldu)
                     ENDIF
                  ENDIF
               ELSE
!                  .. copy R**T into [U] and overwrite [U] with the right
!                  singular vectors of R
                  DO p = 1 , nr
                     DO q = p , N
                        U(q,nr+p) = A(p,q)
                     ENDDO
                  ENDDO
                  IF ( nr>1 ) CALL DLASET('U',nr-1,nr-1,ZERO,ZERO,      &
     &                 U(1,nr+2),Ldu)
                  CALL DGEQRF(N,nr,U(1,nr+1),Ldu,Work(N+1),Work(N+nr+1),&
     &                        Lwork-N-nr,ierr)
                  DO p = 1 , nr
                     DO q = 1 , N
                        V(q,p) = U(p,nr+q)
                     ENDDO
                  ENDDO
                  CALL DLASET('U',nr-1,nr-1,ZERO,ZERO,V(1,2),Ldv)
                  CALL DGESVD('S','O',nr,nr,V,Ldv,S,U,Ldu,V,Ldv,        &
     &                        Work(N+nr+1),Lwork-N-nr,Info)
                  CALL DLASET('A',N-nr,nr,ZERO,ZERO,V(nr+1,1),Ldv)
                  CALL DLASET('A',nr,N-nr,ZERO,ZERO,V(1,nr+1),Ldv)
                  CALL DLASET('A',N-nr,N-nr,ZERO,ONE,V(nr+1,nr+1),Ldv)
                  CALL DORMQR('R','C',N,N,nr,U(1,nr+1),Ldu,Work(N+1),V, &
     &                        Ldv,Work(N+nr+1),Lwork-N-nr,ierr)
                  CALL DLAPMT(.FALSE.,N,N,V,Ldv,Iwork)
!                 .. assemble the left singular vector matrix U of dimensions
!                 (M x NR) or (M x N) or (M x M).
                  IF ( (nr<M) .AND. .NOT.(wntuf) ) THEN
                     CALL DLASET('A',M-nr,nr,ZERO,ZERO,U(nr+1,1),Ldu)
                     IF ( nr<n1 ) THEN
                        CALL DLASET('A',nr,n1-nr,ZERO,ZERO,U(1,nr+1),   &
     &                              Ldu)
                        CALL DLASET('A',M-nr,n1-nr,ZERO,ONE,U(nr+1,nr+1)&
     &                              ,Ldu)
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
!
!
!            .. apply DGESVD to R [[this is the recommended option]]
!
         ELSEIF ( wntvr .OR. (nr==N) ) THEN
!                .. copy R into [V] and overwrite V with the right singular vectors
            CALL DLACPY('U',nr,N,A,Lda,V,Ldv)
            IF ( nr>1 ) CALL DLASET('L',nr-1,nr-1,ZERO,ZERO,V(2,1),Ldv)
!               .. the right singular vectors of R overwrite [V], the NR left
!               singular vectors of R stored in [U](1:NR,1:NR)
            CALL DGESVD('S','O',nr,N,V,Ldv,S,U,Ldu,V,Ldv,Work(N+1),     &
     &                  Lwork-N,Info)
            CALL DLAPMT(.FALSE.,nr,N,V,Ldv,Iwork)
!               .. now [V](1:NR,1:N) contains V(1:N,1:NR)**T
!               .. assemble the left singular vector matrix U of dimensions
!              (M x NR) or (M x N) or (M x M).
            IF ( (nr<M) .AND. .NOT.(wntuf) ) THEN
               CALL DLASET('A',M-nr,nr,ZERO,ZERO,U(nr+1,1),Ldu)
               IF ( nr<n1 ) THEN
                  CALL DLASET('A',nr,n1-nr,ZERO,ZERO,U(1,nr+1),Ldu)
                  CALL DLASET('A',M-nr,n1-nr,ZERO,ONE,U(nr+1,nr+1),Ldu)
               ENDIF
            ENDIF
!
         ELSE
!              .. need all N right singular vectors and NR < N
!              .. the requested number of the left singular vectors
!               is then N1 (N or M)
!               [[The optimal ratio N/NR for using LQ instead of padding
!                 with zeros. Here hard coded to 2; it must be at least
!                 two due to work space constraints.]]
!               OPTRATIO = ILAENV(6, 'DGESVD', 'S' // 'O', NR,N,0,0)
!               OPTRATIO = MAX( OPTRATIO, 2 )
            optratio = 2
            IF ( optratio*nr>N ) THEN
               CALL DLACPY('U',nr,N,A,Lda,V,Ldv)
               IF ( nr>1 ) CALL DLASET('L',nr-1,nr-1,ZERO,ZERO,V(2,1),  &
     &                                 Ldv)
!              .. the right singular vectors of R overwrite [V], the NR left
!                 singular vectors of R stored in [U](1:NR,1:NR)
               CALL DLASET('A',N-nr,N,ZERO,ZERO,V(nr+1,1),Ldv)
               CALL DGESVD('S','O',N,N,V,Ldv,S,U,Ldu,V,Ldv,Work(N+1),   &
     &                     Lwork-N,Info)
               CALL DLAPMT(.FALSE.,N,N,V,Ldv,Iwork)
!                 .. now [V] contains the transposed matrix of the right
!                 singular vectors of A. The leading N left singular vectors
!                 are in [U](1:N,1:N)
!                 .. assemble the left singular vector matrix U of dimensions
!                 (M x N1), i.e. (M x N) or (M x M).
               IF ( (N<M) .AND. .NOT.(wntuf) ) THEN
                  CALL DLASET('A',M-N,N,ZERO,ZERO,U(N+1,1),Ldu)
                  IF ( N<n1 ) THEN
                     CALL DLASET('A',N,n1-N,ZERO,ZERO,U(1,N+1),Ldu)
                     CALL DLASET('A',M-N,n1-N,ZERO,ONE,U(N+1,N+1),Ldu)
                  ENDIF
               ENDIF
            ELSE
               CALL DLACPY('U',nr,N,A,Lda,U(nr+1,1),Ldu)
               IF ( nr>1 ) CALL DLASET('L',nr-1,nr-1,ZERO,ZERO,U(nr+2,1)&
     &                                 ,Ldu)
               CALL DGELQF(nr,N,U(nr+1,1),Ldu,Work(N+1),Work(N+nr+1),   &
     &                     Lwork-N-nr,ierr)
               CALL DLACPY('L',nr,nr,U(nr+1,1),Ldu,V,Ldv)
               IF ( nr>1 ) CALL DLASET('U',nr-1,nr-1,ZERO,ZERO,V(1,2),  &
     &                                 Ldv)
               CALL DGESVD('S','O',nr,nr,V,Ldv,S,U,Ldu,V,Ldv,           &
     &                     Work(N+nr+1),Lwork-N-nr,Info)
               CALL DLASET('A',N-nr,nr,ZERO,ZERO,V(nr+1,1),Ldv)
               CALL DLASET('A',nr,N-nr,ZERO,ZERO,V(1,nr+1),Ldv)
               CALL DLASET('A',N-nr,N-nr,ZERO,ONE,V(nr+1,nr+1),Ldv)
               CALL DORMLQ('R','N',N,N,nr,U(nr+1,1),Ldu,Work(N+1),V,Ldv,&
     &                     Work(N+nr+1),Lwork-N-nr,ierr)
               CALL DLAPMT(.FALSE.,N,N,V,Ldv,Iwork)
!               .. assemble the left singular vector matrix U of dimensions
!              (M x NR) or (M x N) or (M x M).
               IF ( (nr<M) .AND. .NOT.(wntuf) ) THEN
                  CALL DLASET('A',M-nr,nr,ZERO,ZERO,U(nr+1,1),Ldu)
                  IF ( nr<n1 ) THEN
                     CALL DLASET('A',nr,n1-nr,ZERO,ZERO,U(1,nr+1),Ldu)
                     CALL DLASET('A',M-nr,n1-nr,ZERO,ONE,U(nr+1,nr+1),  &
     &                           Ldu)
                  ENDIF
               ENDIF
            ENDIF
!        .. end of the "R**T or R" branch
         ENDIF
!
!           The Q matrix from the first QRF is built into the left singular
!           vectors matrix U.
!
         IF ( .NOT.wntuf ) CALL DORMQR('L','N',M,n1,N,A,Lda,Work,U,Ldu, &
     &                                 Work(N+1),Lwork-N,ierr)
         IF ( rowprm .AND. .NOT.wntuf )                                 &
     &        CALL DLASWP(n1,U,Ldu,1,M-1,Iwork(N+1),-1)
!
!     ... end of the "full SVD" branch
      ENDIF
!
!     Check whether some singular values are returned as zeros, e.g.
!     due to underflow, and update the numerical rank.
      p = nr
      DO q = p , 1 , -1
         IF ( S(q)>ZERO ) EXIT
         nr = nr - 1
      ENDDO
!
!     .. if numerical rank deficiency is detected, the truncated
!     singular values are set to zero.
      IF ( nr<N ) CALL DLASET('G',N-nr,1,ZERO,ZERO,S(nr+1),N)
!     .. undo scaling; this may cause overflow in the largest singular
!     values.
      IF ( ascaled ) CALL DLASCL('G',0,0,ONE,SQRT(DBLE(M)),nr,1,S,N,    &
     &                           ierr)
      IF ( conda ) Rwork(1) = sconda
      Rwork(2) = p - nr
!     .. p-NR is the number of singular values that are computed as
!     exact zeros in DGESVD() applied to the (possibly truncated)
!     full row rank triangular (trapezoidal) factor of A.
      Numrank = nr
!
!
!     End of DGESVDQ
!
      END SUBROUTINE DGESVDQ
