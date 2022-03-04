!*==cgesvdq.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief <b> CGESVDQ computes the singular value decomposition (SVD) with a QR-Preconditioned QR SVD Method for GE matrices</b>
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CGESVDQ + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgesvdq.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgesvdq.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgesvdq.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!      SUBROUTINE CGESVDQ( JOBA, JOBP, JOBR, JOBU, JOBV, M, N, A, LDA,
!                          S, U, LDU, V, LDV, NUMRANK, IWORK, LIWORK,
!                          CWORK, LCWORK, RWORK, LRWORK, INFO )
!
!     .. Scalar Arguments ..
!      IMPLICIT    NONE
!      CHARACTER   JOBA, JOBP, JOBR, JOBU, JOBV
!      INTEGER     M, N, LDA, LDU, LDV, NUMRANK, LIWORK, LCWORK, LRWORK,
!                  INFO
!     ..
!     .. Array Arguments ..
!      COMPLEX     A( LDA, * ), U( LDU, * ), V( LDV, * ), CWORK( * )
!      REAL        S( * ), RWORK( * )
!      INTEGER     IWORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CGESVDQ computes the singular value decomposition (SVD) of a complex
!> M-by-N matrix A, where M >= N. The SVD of A is written as
!>                                    [++]   [xx]   [x0]   [xx]
!>              A = U * SIGMA * V^*,  [++] = [xx] * [ox] * [xx]
!>                                    [++]   [xx]
!> where SIGMA is an N-by-N diagonal matrix, U is an M-by-N orthonormal
!> matrix, and V is an N-by-N unitary matrix. The diagonal elements
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
!>        where EPS = SLAMCH('Epsilon'). This authorises CGESVDQ to
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
!>          = 'T' After the initial pivoted QR factorization, CGESVD is applied to
!>          the adjoint R**H of the computed triangular factor R. This involves
!>          some extra data movement (matrix transpositions). Useful for
!>          experiments, research and development.
!>          = 'N' The triangular factor R is given as input to CGESVD. This may be
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
!>          N left singular vectors of (R**H , 0)**H. If row pivoting is used,
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
!>          A is COMPLEX array of dimensions LDA x N
!>          On entry, the input matrix A.
!>          On exit, if JOBU .NE. 'N' or JOBV .NE. 'N', the lower triangle of A contains
!>          the Householder vectors as stored by CGEQP3. If JOBU = 'F', these Householder
!>          vectors together with CWORK(1:N) can be used to restore the Q factors from
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
!>          S is REAL array of dimension N.
!>          The singular values of A, ordered so that S(i) >= S(i+1).
!> \endverbatim
!>
!> \param[out] U
!> \verbatim
!>          U is COMPLEX array, dimension
!>          LDU x M if JOBU = 'A'; see the description of LDU. In this case,
!>          on exit, U contains the M left singular vectors.
!>          LDU x N if JOBU = 'S', 'U', 'R' ; see the description of LDU. In this
!>          case, U contains the leading N or the leading NUMRANK left singular vectors.
!>          LDU x N if JOBU = 'F' ; see the description of LDU. In this case U
!>          contains N x N unitary matrix that can be used to form the left
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
!>          V is COMPLEX array, dimension
!>          LDV x N if JOBV = 'A', 'V', 'R' or if JOBA = 'E' .
!>          If JOBV = 'A', or 'V',  V contains the N-by-N unitary matrix  V**H;
!>          If JOBV = 'R', V contains the first NUMRANK rows of V**H (the right
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
!>          of CGESVD. The final value of NUMRANK might be further reduced if
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
!>          If LIWORK, LCWORK, or LRWORK = -1, then on exit, if INFO = 0,
!>          IWORK(1) returns the minimal LIWORK.
!> \endverbatim
!>
!> \param[in] LIWORK
!> \verbatim
!>          LIWORK is INTEGER
!>          The dimension of the array IWORK.
!>          LIWORK >= N + M - 1,  if JOBP = 'P';
!>          LIWORK >= N           if JOBP = 'N'.
!>
!>          If LIWORK = -1, then a workspace query is assumed; the routine
!>          only calculates and returns the optimal and minimal sizes
!>          for the CWORK, IWORK, and RWORK arrays, and no error
!>          message related to LCWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] CWORK
!> \verbatim
!>          CWORK is COMPLEX array, dimension (max(2, LCWORK)), used as a workspace.
!>          On exit, if, on entry, LCWORK.NE.-1, CWORK(1:N) contains parameters
!>          needed to recover the Q factor from the QR factorization computed by
!>          CGEQP3.
!>
!>          If LIWORK, LCWORK, or LRWORK = -1, then on exit, if INFO = 0,
!>          CWORK(1) returns the optimal LCWORK, and
!>          CWORK(2) returns the minimal LCWORK.
!> \endverbatim
!>
!> \param[in,out] LCWORK
!> \verbatim
!>          LCWORK is INTEGER
!>          The dimension of the array CWORK. It is determined as follows:
!>          Let  LWQP3 = N+1,  LWCON = 2*N, and let
!>          LWUNQ = { MAX( N, 1 ),  if JOBU = 'R', 'S', or 'U'
!>                  { MAX( M, 1 ),  if JOBU = 'A'
!>          LWSVD = MAX( 3*N, 1 )
!>          LWLQF = MAX( N/2, 1 ), LWSVD2 = MAX( 3*(N/2), 1 ), LWUNLQ = MAX( N, 1 ),
!>          LWQRF = MAX( N/2, 1 ), LWUNQ2 = MAX( N, 1 )
!>          Then the minimal value of LCWORK is:
!>          = MAX( N + LWQP3, LWSVD )        if only the singular values are needed;
!>          = MAX( N + LWQP3, LWCON, LWSVD ) if only the singular values are needed,
!>                                   and a scaled condition estimate requested;
!>
!>          = N + MAX( LWQP3, LWSVD, LWUNQ ) if the singular values and the left
!>                                   singular vectors are requested;
!>          = N + MAX( LWQP3, LWCON, LWSVD, LWUNQ ) if the singular values and the left
!>                                   singular vectors are requested, and also
!>                                   a scaled condition estimate requested;
!>
!>          = N + MAX( LWQP3, LWSVD )        if the singular values and the right
!>                                   singular vectors are requested;
!>          = N + MAX( LWQP3, LWCON, LWSVD ) if the singular values and the right
!>                                   singular vectors are requested, and also
!>                                   a scaled condition etimate requested;
!>
!>          = N + MAX( LWQP3, LWSVD, LWUNQ ) if the full SVD is requested with JOBV = 'R';
!>                                   independent of JOBR;
!>          = N + MAX( LWQP3, LWCON, LWSVD, LWUNQ ) if the full SVD is requested,
!>                                   JOBV = 'R' and, also a scaled condition
!>                                   estimate requested; independent of JOBR;
!>          = MAX( N + MAX( LWQP3, LWSVD, LWUNQ ),
!>         N + MAX( LWQP3, N/2+LWLQF, N/2+LWSVD2, N/2+LWUNLQ, LWUNQ) ) if the
!>                         full SVD is requested with JOBV = 'A' or 'V', and
!>                         JOBR ='N'
!>          = MAX( N + MAX( LWQP3, LWCON, LWSVD, LWUNQ ),
!>         N + MAX( LWQP3, LWCON, N/2+LWLQF, N/2+LWSVD2, N/2+LWUNLQ, LWUNQ ) )
!>                         if the full SVD is requested with JOBV = 'A' or 'V', and
!>                         JOBR ='N', and also a scaled condition number estimate
!>                         requested.
!>          = MAX( N + MAX( LWQP3, LWSVD, LWUNQ ),
!>         N + MAX( LWQP3, N/2+LWQRF, N/2+LWSVD2, N/2+LWUNQ2, LWUNQ ) ) if the
!>                         full SVD is requested with JOBV = 'A', 'V', and JOBR ='T'
!>          = MAX( N + MAX( LWQP3, LWCON, LWSVD, LWUNQ ),
!>         N + MAX( LWQP3, LWCON, N/2+LWQRF, N/2+LWSVD2, N/2+LWUNQ2, LWUNQ ) )
!>                         if the full SVD is requested with JOBV = 'A', 'V' and
!>                         JOBR ='T', and also a scaled condition number estimate
!>                         requested.
!>          Finally, LCWORK must be at least two: LCWORK = MAX( 2, LCWORK ).
!>
!>          If LCWORK = -1, then a workspace query is assumed; the routine
!>          only calculates and returns the optimal and minimal sizes
!>          for the CWORK, IWORK, and RWORK arrays, and no error
!>          message related to LCWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is REAL array, dimension (max(1, LRWORK)).
!>          On exit,
!>          1. If JOBA = 'E', RWORK(1) contains an estimate of the condition
!>          number of column scaled A. If A = C * D where D is diagonal and C
!>          has unit columns in the Euclidean norm, then, assuming full column rank,
!>          N^(-1/4) * RWORK(1) <= ||pinv(C)||_2 <= N^(1/4) * RWORK(1).
!>          Otherwise, RWORK(1) = -1.
!>          2. RWORK(2) contains the number of singular values computed as
!>          exact zeros in CGESVD applied to the upper triangular or trapezoidal
!>          R (from the initial QR factorization). In case of early exit (no call to
!>          CGESVD, such as in the case of zero matrix) RWORK(2) = -1.
!>
!>          If LIWORK, LCWORK, or LRWORK = -1, then on exit, if INFO = 0,
!>          RWORK(1) returns the minimal LRWORK.
!> \endverbatim
!>
!> \param[in] LRWORK
!> \verbatim
!>          LRWORK is INTEGER.
!>          The dimension of the array RWORK.
!>          If JOBP ='P', then LRWORK >= MAX(2, M, 5*N);
!>          Otherwise, LRWORK >= MAX(2, 5*N).
!>
!>          If LRWORK = -1, then a workspace query is assumed; the routine
!>          only calculates and returns the optimal and minimal sizes
!>          for the CWORK, IWORK, and RWORK arrays, and no error
!>          message related to LCWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit.
!>          < 0:  if INFO = -i, the i-th argument had an illegal value.
!>          > 0:  if CBDSQR did not converge, INFO specifies how many superdiagonals
!>          of an intermediate bidiagonal form B (computed in CGESVD) did not
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
!> \ingroup complexGEsing
!
!  =====================================================================
      SUBROUTINE CGESVDQ(Joba,Jobp,Jobr,Jobu,Jobv,M,N,A,Lda,S,U,Ldu,V,  &
     &                   Ldv,Numrank,Iwork,Liwork,Cwork,Lcwork,Rwork,   &
     &                   Lrwork,Info)
!     .. Scalar Arguments ..
      USE S_CGELQF
      USE S_CGEQP3
      USE S_CGEQRF
      USE S_CGESVD
      USE S_CLACPY
      USE S_CLANGE
      USE S_CLAPMT
      USE S_CLASCL
      USE S_CLASET
      USE S_CLASWP
      USE S_CPOCON
      USE S_CSSCAL
      USE S_CUNMLQ
      USE S_CUNMQR
      USE S_ISAMAX
      USE S_LSAME
      USE S_SCNRM2
      USE S_SLAMCH
      USE S_SLASCL
      USE S_SLASET
      USE S_XERBLA
      IMPLICIT NONE
!*--CGESVDQ439
!
! PARAMETER definitions rewritten by SPAG
!
      REAL , PARAMETER  ::  ZERO = 0.0E0 , ONE = 1.0E0
      COMPLEX , PARAMETER  ::  CZERO = (0.0E0,0.0E0) ,                  &
     &                         CONE = (1.0E0,0.0E0)
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Joba
      CHARACTER :: Jobp
      CHARACTER :: Jobr
      CHARACTER :: Jobu
      CHARACTER :: Jobv
      INTEGER :: M
      INTEGER :: N
      COMPLEX , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(*) :: S
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldu,*) :: U
      INTEGER :: Ldu
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldv,*) :: V
      INTEGER :: Ldv
      INTEGER , INTENT(OUT) :: Numrank
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(IN) :: Liwork
      COMPLEX , DIMENSION(*) :: Cwork
      INTEGER :: Lcwork
      REAL , INTENT(INOUT) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(IN) :: Lrwork
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      LOGICAL :: accla , acclh , acclm , ascaled , conda , dntwu ,      &
     &           dntwv , lquery , lsvc0 , lsvec , rowprm , rsvec ,      &
     &           rtrans , wntua , wntuf , wntur , wntus , wntva , wntvr
      REAL :: big , epsln , rtmp , sconda , sfmin
      COMPLEX , DIMENSION(1) :: cdummy
      COMPLEX :: ctmp
      INTEGER :: ierr , iminwrk , lwcon , lwlqf , lwqp3 , lwqrf ,       &
     &           lwrk_cgelqf , lwrk_cgeqp3 , lwrk_cgeqrf , lwrk_cgesvd ,&
     &           lwrk_cgesvd2 , lwrk_cunmlq , lwrk_cunmqr ,             &
     &           lwrk_cunmqr2 , lwsvd , lwsvd2 , lwunlq , lwunq ,       &
     &           lwunq2 , minwrk , minwrk2 , n1 , nr , optratio ,       &
     &           optwrk , optwrk2 , p , q , rminwrk
      REAL , DIMENSION(1) :: rdummy
!
! End of declarations rewritten by SPAG
!
!     ..
!     .. Array Arguments ..
!
!  =====================================================================
!
!     .. Parameters ..
!     ..
!     .. Local Scalars ..
!     ..
!     .. Local Arrays
!     ..
!     .. External Subroutines (BLAS, LAPACK)
!     ..
!     .. External Functions (BLAS, LAPACK)
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. Executable Statements ..
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
         iminwrk = MAX(1,N+M-1)
         rminwrk = MAX(2,M,5*N)
      ELSE
         iminwrk = MAX(1,N)
         rminwrk = MAX(2,5*N)
      ENDIF
      lquery = (Liwork==-1 .OR. Lcwork==-1 .OR. Lrwork==-1)
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
!
!     Compute workspace
!        .. compute the minimal and the optimal workspace lengths
!        [[The expressions for computing the minimal and the optimal
!        values of LCWORK are written with a lot of redundancy and
!        can be simplified. However, this detailed form is easier for
!        maintenance and modifications of the code.]]
!
!        .. minimal workspace length for CGEQP3 of an M x N matrix
         lwqp3 = N + 1
!        .. minimal workspace length for CUNMQR to build left singular vectors
         IF ( wntus .OR. wntur ) THEN
            lwunq = MAX(N,1)
         ELSEIF ( wntua ) THEN
            lwunq = MAX(M,1)
         ENDIF
!        .. minimal workspace length for CPOCON of an N x N matrix
         lwcon = 2*N
!        .. CGESVD of an N x N matrix
         lwsvd = MAX(3*N,1)
         IF ( lquery ) THEN
            CALL CGEQP3(M,N,A,Lda,Iwork,cdummy,cdummy,-1,rdummy,ierr)
            lwrk_cgeqp3 = INT(cdummy(1))
            IF ( wntus .OR. wntur ) THEN
               CALL CUNMQR('L','N',M,N,N,A,Lda,cdummy,U,Ldu,cdummy,-1,  &
     &                     ierr)
               lwrk_cunmqr = INT(cdummy(1))
            ELSEIF ( wntua ) THEN
               CALL CUNMQR('L','N',M,M,N,A,Lda,cdummy,U,Ldu,cdummy,-1,  &
     &                     ierr)
               lwrk_cunmqr = INT(cdummy(1))
            ELSE
               lwrk_cunmqr = 0
            ENDIF
         ENDIF
         minwrk = 2
         optwrk = 2
         IF ( .NOT.(lsvec .OR. rsvec) ) THEN
!            .. minimal and optimal sizes of the complex workspace if
!            only the singular values are requested
            IF ( conda ) THEN
               minwrk = MAX(N+lwqp3,lwcon,lwsvd)
            ELSE
               minwrk = MAX(N+lwqp3,lwsvd)
            ENDIF
            IF ( lquery ) THEN
               CALL CGESVD('N','N',N,N,A,Lda,S,U,Ldu,V,Ldv,cdummy,-1,   &
     &                     rdummy,ierr)
               lwrk_cgesvd = INT(cdummy(1))
               IF ( conda ) THEN
                  optwrk = MAX(N+lwrk_cgeqp3,N+lwcon,lwrk_cgesvd)
               ELSE
                  optwrk = MAX(N+lwrk_cgeqp3,lwrk_cgesvd)
               ENDIF
            ENDIF
         ELSEIF ( lsvec .AND. (.NOT.rsvec) ) THEN
!            .. minimal and optimal sizes of the complex workspace if the
!            singular values and the left singular vectors are requested
            IF ( conda ) THEN
               minwrk = N + MAX(lwqp3,lwcon,lwsvd,lwunq)
            ELSE
               minwrk = N + MAX(lwqp3,lwsvd,lwunq)
            ENDIF
            IF ( lquery ) THEN
               IF ( rtrans ) THEN
                  CALL CGESVD('N','O',N,N,A,Lda,S,U,Ldu,V,Ldv,cdummy,-1,&
     &                        rdummy,ierr)
               ELSE
                  CALL CGESVD('O','N',N,N,A,Lda,S,U,Ldu,V,Ldv,cdummy,-1,&
     &                        rdummy,ierr)
               ENDIF
               lwrk_cgesvd = INT(cdummy(1))
               IF ( conda ) THEN
                  optwrk = N + MAX(lwrk_cgeqp3,lwcon,lwrk_cgesvd,       &
     &                     lwrk_cunmqr)
               ELSE
                  optwrk = N + MAX(lwrk_cgeqp3,lwrk_cgesvd,lwrk_cunmqr)
               ENDIF
            ENDIF
         ELSEIF ( rsvec .AND. (.NOT.lsvec) ) THEN
!            .. minimal and optimal sizes of the complex workspace if the
!            singular values and the right singular vectors are requested
            IF ( conda ) THEN
               minwrk = N + MAX(lwqp3,lwcon,lwsvd)
            ELSE
               minwrk = N + MAX(lwqp3,lwsvd)
            ENDIF
            IF ( lquery ) THEN
               IF ( rtrans ) THEN
                  CALL CGESVD('O','N',N,N,A,Lda,S,U,Ldu,V,Ldv,cdummy,-1,&
     &                        rdummy,ierr)
               ELSE
                  CALL CGESVD('N','O',N,N,A,Lda,S,U,Ldu,V,Ldv,cdummy,-1,&
     &                        rdummy,ierr)
               ENDIF
               lwrk_cgesvd = INT(cdummy(1))
               IF ( conda ) THEN
                  optwrk = N + MAX(lwrk_cgeqp3,lwcon,lwrk_cgesvd)
               ELSE
                  optwrk = N + MAX(lwrk_cgeqp3,lwrk_cgesvd)
               ENDIF
            ENDIF
         ELSE
!            .. minimal and optimal sizes of the complex workspace if the
!            full SVD is requested
            IF ( rtrans ) THEN
               minwrk = MAX(lwqp3,lwsvd,lwunq)
               IF ( conda ) minwrk = MAX(minwrk,lwcon)
               minwrk = minwrk + N
               IF ( wntva ) THEN
!                   .. minimal workspace length for N x N/2 CGEQRF
                  lwqrf = MAX(N/2,1)
!                   .. minimal workspace length for N/2 x N/2 CGESVD
                  lwsvd2 = MAX(3*(N/2),1)
                  lwunq2 = MAX(N,1)
                  minwrk2 = MAX(lwqp3,N/2+lwqrf,N/2+lwsvd2,N/2+lwunq2,  &
     &                      lwunq)
                  IF ( conda ) minwrk2 = MAX(minwrk2,lwcon)
                  minwrk2 = N + minwrk2
                  minwrk = MAX(minwrk,minwrk2)
               ENDIF
            ELSE
               minwrk = MAX(lwqp3,lwsvd,lwunq)
               IF ( conda ) minwrk = MAX(minwrk,lwcon)
               minwrk = minwrk + N
               IF ( wntva ) THEN
!                   .. minimal workspace length for N/2 x N CGELQF
                  lwlqf = MAX(N/2,1)
                  lwsvd2 = MAX(3*(N/2),1)
                  lwunlq = MAX(N,1)
                  minwrk2 = MAX(lwqp3,N/2+lwlqf,N/2+lwsvd2,N/2+lwunlq,  &
     &                      lwunq)
                  IF ( conda ) minwrk2 = MAX(minwrk2,lwcon)
                  minwrk2 = N + minwrk2
                  minwrk = MAX(minwrk,minwrk2)
               ENDIF
            ENDIF
            IF ( lquery ) THEN
               IF ( rtrans ) THEN
                  CALL CGESVD('O','A',N,N,A,Lda,S,U,Ldu,V,Ldv,cdummy,-1,&
     &                        rdummy,ierr)
                  lwrk_cgesvd = INT(cdummy(1))
                  optwrk = MAX(lwrk_cgeqp3,lwrk_cgesvd,lwrk_cunmqr)
                  IF ( conda ) optwrk = MAX(optwrk,lwcon)
                  optwrk = N + optwrk
                  IF ( wntva ) THEN
                     CALL CGEQRF(N,N/2,U,Ldu,cdummy,cdummy,-1,ierr)
                     lwrk_cgeqrf = INT(cdummy(1))
                     CALL CGESVD('S','O',N/2,N/2,V,Ldv,S,U,Ldu,V,Ldv,   &
     &                           cdummy,-1,rdummy,ierr)
                     lwrk_cgesvd2 = INT(cdummy(1))
                     CALL CUNMQR('R','C',N,N,N/2,U,Ldu,cdummy,V,Ldv,    &
     &                           cdummy,-1,ierr)
                     lwrk_cunmqr2 = INT(cdummy(1))
                     optwrk2 = MAX(lwrk_cgeqp3,N/2+lwrk_cgeqrf,         &
     &                         N/2+lwrk_cgesvd2,N/2+lwrk_cunmqr2)
                     IF ( conda ) optwrk2 = MAX(optwrk2,lwcon)
                     optwrk2 = N + optwrk2
                     optwrk = MAX(optwrk,optwrk2)
                  ENDIF
               ELSE
                  CALL CGESVD('S','O',N,N,A,Lda,S,U,Ldu,V,Ldv,cdummy,-1,&
     &                        rdummy,ierr)
                  lwrk_cgesvd = INT(cdummy(1))
                  optwrk = MAX(lwrk_cgeqp3,lwrk_cgesvd,lwrk_cunmqr)
                  IF ( conda ) optwrk = MAX(optwrk,lwcon)
                  optwrk = N + optwrk
                  IF ( wntva ) THEN
                     CALL CGELQF(N/2,N,U,Ldu,cdummy,cdummy,-1,ierr)
                     lwrk_cgelqf = INT(cdummy(1))
                     CALL CGESVD('S','O',N/2,N/2,V,Ldv,S,U,Ldu,V,Ldv,   &
     &                           cdummy,-1,rdummy,ierr)
                     lwrk_cgesvd2 = INT(cdummy(1))
                     CALL CUNMLQ('R','N',N,N,N/2,U,Ldu,cdummy,V,Ldv,    &
     &                           cdummy,-1,ierr)
                     lwrk_cunmlq = INT(cdummy(1))
                     optwrk2 = MAX(lwrk_cgeqp3,N/2+lwrk_cgelqf,         &
     &                         N/2+lwrk_cgesvd2,N/2+lwrk_cunmlq)
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
         IF ( Lcwork<minwrk .AND. (.NOT.lquery) ) Info = -19
!
      ENDIF
!
      IF ( Info==0 .AND. Lrwork<rminwrk .AND. .NOT.lquery ) Info = -21
      IF ( Info/=0 ) THEN
         CALL XERBLA('CGESVDQ',-Info)
         RETURN
      ELSEIF ( lquery ) THEN
!
!     Return optimal workspace
!
         Iwork(1) = iminwrk
         Cwork(1) = optwrk
         Cwork(2) = minwrk
         Rwork(1) = rminwrk
         RETURN
      ENDIF
!
!     Quick return if the matrix is void.
!
!     .. all output is void.
      IF ( (M==0) .OR. (N==0) ) RETURN
!
      big = SLAMCH('O')
      ascaled = .FALSE.
      IF ( rowprm ) THEN
!           .. reordering the rows in decreasing sequence in the
!           ell-infinity norm - this enhances numerical robustness in
!           the case of differently scaled rows.
         DO p = 1 , M
!               RWORK(p) = ABS( A(p,ICAMAX(N,A(p,1),LDA)) )
!               [[CLANGE will return NaN if an entry of the p-th row is Nan]]
            Rwork(p) = CLANGE('M',1,N,A(p,1),Lda,rdummy)
!               .. check for NaN's and Inf's
            IF ( (Rwork(p)/=Rwork(p)) .OR. ((Rwork(p)*ZERO)/=ZERO) )    &
     &           THEN
               Info = -8
               CALL XERBLA('CGESVDQ',-Info)
               RETURN
            ENDIF
         ENDDO
         DO p = 1 , M - 1
            q = ISAMAX(M-p+1,Rwork(p),1) + p - 1
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
            CALL SLASET('G',N,1,ZERO,ZERO,S,N)
            IF ( wntus ) CALL CLASET('G',M,N,CZERO,CONE,U,Ldu)
            IF ( wntua ) CALL CLASET('G',M,M,CZERO,CONE,U,Ldu)
            IF ( wntva ) CALL CLASET('G',N,N,CZERO,CONE,V,Ldv)
            IF ( wntuf ) THEN
               CALL CLASET('G',N,1,CZERO,CZERO,Cwork,N)
               CALL CLASET('G',M,N,CZERO,CONE,U,Ldu)
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
         IF ( Rwork(1)>big/SQRT(REAL(M)) ) THEN
!               .. to prevent overflow in the QR factorization, scale the
!               matrix by 1/sqrt(M) if too large entry detected
            CALL CLASCL('G',0,0,SQRT(REAL(M)),ONE,M,N,A,Lda,ierr)
            ascaled = .TRUE.
         ENDIF
         CALL CLASWP(N,A,Lda,1,M-1,Iwork(N+1),1)
      ENDIF
!
!    .. At this stage, preemptive scaling is done only to avoid column
!    norms overflows during the QR factorization. The SVD procedure should
!    have its own scaling to save the singular values from overflows and
!    underflows. That depends on the SVD procedure.
!
      IF ( .NOT.rowprm ) THEN
         rtmp = CLANGE('M',M,N,A,Lda,Rwork)
         IF ( (rtmp/=rtmp) .OR. ((rtmp*ZERO)/=ZERO) ) THEN
            Info = -8
            CALL XERBLA('CGESVDQ',-Info)
            RETURN
         ENDIF
         IF ( rtmp>big/SQRT(REAL(M)) ) THEN
!             .. to prevent overflow in the QR factorization, scale the
!             matrix by 1/sqrt(M) if too large entry detected
            CALL CLASCL('G',0,0,SQRT(REAL(M)),ONE,M,N,A,Lda,ierr)
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
      CALL CGEQP3(M,N,A,Lda,Iwork,Cwork,Cwork(N+1),Lcwork-N,Rwork,ierr)
!
!    If the user requested accuracy level allows truncation in the
!    computed upper triangular factor, the matrix R is examined and,
!    if possible, replaced with its leading upper trapezoidal part.
!
      epsln = SLAMCH('E')
      sfmin = SLAMCH('S')
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
         rtmp = SQRT(REAL(N))*epsln
         DO p = 2 , N
            IF ( ABS(A(p,p))<(rtmp*ABS(A(1,1))) ) EXIT
            nr = nr + 1
         ENDDO
!
      ELSEIF ( acclm ) THEN
!        .. similarly as above, only slightly more gentle (less aggressive).
!        Sudden drop on the diagonal of R is used as the criterion for being
!        close-to-rank-deficient. The threshold is set to EPSLN=SLAMCH('E').
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
            CALL CLACPY('U',N,N,A,Lda,V,Ldv)
!              Only the leading NR x NR submatrix of the triangular factor
!              is considered. Only if NR=N will this give a reliable error
!              bound. However, even for NR < N, this can be used on an
!              expert level and obtain useful information in the sense of
!              perturbation theory.
            DO p = 1 , nr
               rtmp = SCNRM2(p,V(1,p),1)
               CALL CSSCAL(p,ONE/rtmp,V(1,p),1)
            ENDDO
            IF ( .NOT.(lsvec .OR. rsvec) ) THEN
               CALL CPOCON('U',nr,V,Ldv,ONE,rtmp,Cwork,Rwork,ierr)
            ELSE
               CALL CPOCON('U',nr,V,Ldv,ONE,rtmp,Cwork(N+1),Rwork,ierr)
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
!         .. compute the singular values of R**H = [A](1:NR,1:N)**H
!           .. set the lower triangle of [A] to [A](1:NR,1:N)**H and
!           the upper triangle of [A] to zero.
            DO p = 1 , MIN(N,nr)
               A(p,p) = CONJG(A(p,p))
               DO q = p + 1 , N
                  A(q,p) = CONJG(A(p,q))
                  IF ( q<=nr ) A(p,q) = CZERO
               ENDDO
            ENDDO
!
            CALL CGESVD('N','N',N,nr,A,Lda,S,U,Ldu,V,Ldv,Cwork,Lcwork,  &
     &                  Rwork,Info)
!
         ELSE
!
!           .. compute the singular values of R = [A](1:NR,1:N)
!
            IF ( nr>1 ) CALL CLASET('L',nr-1,nr-1,CZERO,CZERO,A(2,1),   &
     &                              Lda)
            CALL CGESVD('N','N',nr,N,A,Lda,S,U,Ldu,V,Ldv,Cwork,Lcwork,  &
     &                  Rwork,Info)
!
         ENDIF
!
      ELSEIF ( lsvec .AND. (.NOT.rsvec) ) THEN
!.......................................................................
!       .. the singular values and the left singular vectors requested
!.......................................................................""""""""
         IF ( rtrans ) THEN
!            .. apply CGESVD to R**H
!            .. copy R**H into [U] and overwrite [U] with the right singular
!            vectors of R
            DO p = 1 , nr
               DO q = p , N
                  U(q,p) = CONJG(A(p,q))
               ENDDO
            ENDDO
            IF ( nr>1 ) CALL CLASET('U',nr-1,nr-1,CZERO,CZERO,U(1,2),   &
     &                              Ldu)
!           .. the left singular vectors not computed, the NR right singular
!           vectors overwrite [U](1:NR,1:NR) as conjugate transposed. These
!           will be pre-multiplied by Q to build the left singular vectors of A.
            CALL CGESVD('N','O',N,nr,U,Ldu,S,U,Ldu,U,Ldu,Cwork(N+1),    &
     &                  Lcwork-N,Rwork,Info)
!
            DO p = 1 , nr
               U(p,p) = CONJG(U(p,p))
               DO q = p + 1 , nr
                  ctmp = CONJG(U(q,p))
                  U(q,p) = CONJG(U(p,q))
                  U(p,q) = ctmp
               ENDDO
            ENDDO
!
         ELSE
!            .. apply CGESVD to R
!            .. copy R into [U] and overwrite [U] with the left singular vectors
            CALL CLACPY('U',nr,N,A,Lda,U,Ldu)
            IF ( nr>1 ) CALL CLASET('L',nr-1,nr-1,CZERO,CZERO,U(2,1),   &
     &                              Ldu)
!            .. the right singular vectors not computed, the NR left singular
!            vectors overwrite [U](1:NR,1:NR)
            CALL CGESVD('O','N',nr,N,U,Ldu,S,U,Ldu,V,Ldv,Cwork(N+1),    &
     &                  Lcwork-N,Rwork,Info)
!               .. now [U](1:NR,1:NR) contains the NR left singular vectors of
!               R. These will be pre-multiplied by Q to build the left singular
!               vectors of A.
         ENDIF
!
!           .. assemble the left singular vector matrix U of dimensions
!              (M x NR) or (M x N) or (M x M).
         IF ( (nr<M) .AND. (.NOT.wntuf) ) THEN
            CALL CLASET('A',M-nr,nr,CZERO,CZERO,U(nr+1,1),Ldu)
            IF ( nr<n1 ) THEN
               CALL CLASET('A',nr,n1-nr,CZERO,CZERO,U(1,nr+1),Ldu)
               CALL CLASET('A',M-nr,n1-nr,CZERO,CONE,U(nr+1,nr+1),Ldu)
            ENDIF
         ENDIF
!
!           The Q matrix from the first QRF is built into the left singular
!           vectors matrix U.
!
         IF ( .NOT.wntuf ) CALL CUNMQR('L','N',M,n1,N,A,Lda,Cwork,U,Ldu,&
     &                                 Cwork(N+1),Lcwork-N,ierr)
         IF ( rowprm .AND. .NOT.wntuf )                                 &
     &        CALL CLASWP(n1,U,Ldu,1,M-1,Iwork(N+1),-1)
!
      ELSEIF ( rsvec .AND. (.NOT.lsvec) ) THEN
!.......................................................................
!       .. the singular values and the right singular vectors requested
!.......................................................................
         IF ( rtrans ) THEN
!            .. apply CGESVD to R**H
!            .. copy R**H into V and overwrite V with the left singular vectors
            DO p = 1 , nr
               DO q = p , N
                  V(q,p) = CONJG(A(p,q))
               ENDDO
            ENDDO
            IF ( nr>1 ) CALL CLASET('U',nr-1,nr-1,CZERO,CZERO,V(1,2),   &
     &                              Ldv)
!           .. the left singular vectors of R**H overwrite V, the right singular
!           vectors not computed
            IF ( wntvr .OR. (nr==N) ) THEN
               CALL CGESVD('O','N',N,nr,V,Ldv,S,U,Ldu,U,Ldu,Cwork(N+1), &
     &                     Lcwork-N,Rwork,Info)
!
               DO p = 1 , nr
                  V(p,p) = CONJG(V(p,p))
                  DO q = p + 1 , nr
                     ctmp = CONJG(V(q,p))
                     V(q,p) = CONJG(V(p,q))
                     V(p,q) = ctmp
                  ENDDO
               ENDDO
!
               IF ( nr<N ) THEN
                  DO p = 1 , nr
                     DO q = nr + 1 , N
                        V(p,q) = CONJG(V(q,p))
                     ENDDO
                  ENDDO
               ENDIF
               CALL CLAPMT(.FALSE.,nr,N,V,Ldv,Iwork)
            ELSE
!               .. need all N right singular vectors and NR < N
!               [!] This is simple implementation that augments [V](1:N,1:NR)
!               by padding a zero block. In the case NR << N, a more efficient
!               way is to first use the QR factorization. For more details
!               how to implement this, see the " FULL SVD " branch.
               CALL CLASET('G',N,N-nr,CZERO,CZERO,V(1,nr+1),Ldv)
               CALL CGESVD('O','N',N,N,V,Ldv,S,U,Ldu,U,Ldu,Cwork(N+1),  &
     &                     Lcwork-N,Rwork,Info)
!
               DO p = 1 , N
                  V(p,p) = CONJG(V(p,p))
                  DO q = p + 1 , N
                     ctmp = CONJG(V(q,p))
                     V(q,p) = CONJG(V(p,q))
                     V(p,q) = ctmp
                  ENDDO
               ENDDO
               CALL CLAPMT(.FALSE.,N,N,V,Ldv,Iwork)
            ENDIF
!
         ELSE
!            .. aply CGESVD to R
!            .. copy R into V and overwrite V with the right singular vectors
            CALL CLACPY('U',nr,N,A,Lda,V,Ldv)
            IF ( nr>1 ) CALL CLASET('L',nr-1,nr-1,CZERO,CZERO,V(2,1),   &
     &                              Ldv)
!            .. the right singular vectors overwrite V, the NR left singular
!            vectors stored in U(1:NR,1:NR)
            IF ( wntvr .OR. (nr==N) ) THEN
               CALL CGESVD('N','O',nr,N,V,Ldv,S,U,Ldu,V,Ldv,Cwork(N+1), &
     &                     Lcwork-N,Rwork,Info)
               CALL CLAPMT(.FALSE.,nr,N,V,Ldv,Iwork)
!               .. now [V](1:NR,1:N) contains V(1:N,1:NR)**H
            ELSE
!               .. need all N right singular vectors and NR < N
!               [!] This is simple implementation that augments [V](1:NR,1:N)
!               by padding a zero block. In the case NR << N, a more efficient
!               way is to first use the LQ factorization. For more details
!               how to implement this, see the " FULL SVD " branch.
               CALL CLASET('G',N-nr,N,CZERO,CZERO,V(nr+1,1),Ldv)
               CALL CGESVD('N','O',N,N,V,Ldv,S,U,Ldu,V,Ldv,Cwork(N+1),  &
     &                     Lcwork-N,Rwork,Info)
               CALL CLAPMT(.FALSE.,N,N,V,Ldv,Iwork)
            ENDIF
!            .. now [V] contains the adjoint of the matrix of the right singular
!            vectors of A.
         ENDIF
!
      ELSE
!.......................................................................
!       .. FULL SVD requested
!.......................................................................
         IF ( rtrans ) THEN
!
!            .. apply CGESVD to R**H [[this option is left for R&D&T]]
!
            IF ( wntvr .OR. (nr==N) ) THEN
!            .. copy R**H into [V] and overwrite [V] with the left singular
!            vectors of R**H
               DO p = 1 , nr
                  DO q = p , N
                     V(q,p) = CONJG(A(p,q))
                  ENDDO
               ENDDO
               IF ( nr>1 ) CALL CLASET('U',nr-1,nr-1,CZERO,CZERO,V(1,2),&
     &                                 Ldv)
!
!           .. the left singular vectors of R**H overwrite [V], the NR right
!           singular vectors of R**H stored in [U](1:NR,1:NR) as conjugate
!           transposed
               CALL CGESVD('O','A',N,nr,V,Ldv,S,V,Ldv,U,Ldu,Cwork(N+1), &
     &                     Lcwork-N,Rwork,Info)
!              .. assemble V
               DO p = 1 , nr
                  V(p,p) = CONJG(V(p,p))
                  DO q = p + 1 , nr
                     ctmp = CONJG(V(q,p))
                     V(q,p) = CONJG(V(p,q))
                     V(p,q) = ctmp
                  ENDDO
               ENDDO
               IF ( nr<N ) THEN
                  DO p = 1 , nr
                     DO q = nr + 1 , N
                        V(p,q) = CONJG(V(q,p))
                     ENDDO
                  ENDDO
               ENDIF
               CALL CLAPMT(.FALSE.,nr,N,V,Ldv,Iwork)
!
               DO p = 1 , nr
                  U(p,p) = CONJG(U(p,p))
                  DO q = p + 1 , nr
                     ctmp = CONJG(U(q,p))
                     U(q,p) = CONJG(U(p,q))
                     U(p,q) = ctmp
                  ENDDO
               ENDDO
!
               IF ( (nr<M) .AND. .NOT.(wntuf) ) THEN
                  CALL CLASET('A',M-nr,nr,CZERO,CZERO,U(nr+1,1),Ldu)
                  IF ( nr<n1 ) THEN
                     CALL CLASET('A',nr,n1-nr,CZERO,CZERO,U(1,nr+1),Ldu)
                     CALL CLASET('A',M-nr,n1-nr,CZERO,CONE,U(nr+1,nr+1),&
     &                           Ldu)
                  ENDIF
               ENDIF
!
            ELSE
!               .. need all N right singular vectors and NR < N
!            .. copy R**H into [V] and overwrite [V] with the left singular
!            vectors of R**H
!               [[The optimal ratio N/NR for using QRF instead of padding
!                 with zeros. Here hard coded to 2; it must be at least
!                 two due to work space constraints.]]
!               OPTRATIO = ILAENV(6, 'CGESVD', 'S' // 'O', NR,N,0,0)
!               OPTRATIO = MAX( OPTRATIO, 2 )
               optratio = 2
               IF ( optratio*nr>N ) THEN
                  DO p = 1 , nr
                     DO q = p , N
                        V(q,p) = CONJG(A(p,q))
                     ENDDO
                  ENDDO
                  IF ( nr>1 ) CALL CLASET('U',nr-1,nr-1,CZERO,CZERO,    &
     &                 V(1,2),Ldv)
!
                  CALL CLASET('A',N,N-nr,CZERO,CZERO,V(1,nr+1),Ldv)
                  CALL CGESVD('O','A',N,N,V,Ldv,S,V,Ldv,U,Ldu,Cwork(N+1)&
     &                        ,Lcwork-N,Rwork,Info)
!
                  DO p = 1 , N
                     V(p,p) = CONJG(V(p,p))
                     DO q = p + 1 , N
                        ctmp = CONJG(V(q,p))
                        V(q,p) = CONJG(V(p,q))
                        V(p,q) = ctmp
                     ENDDO
                  ENDDO
                  CALL CLAPMT(.FALSE.,N,N,V,Ldv,Iwork)
!              .. assemble the left singular vector matrix U of dimensions
!              (M x N1), i.e. (M x N) or (M x M).
!
                  DO p = 1 , N
                     U(p,p) = CONJG(U(p,p))
                     DO q = p + 1 , N
                        ctmp = CONJG(U(q,p))
                        U(q,p) = CONJG(U(p,q))
                        U(p,q) = ctmp
                     ENDDO
                  ENDDO
!
                  IF ( (N<M) .AND. .NOT.(wntuf) ) THEN
                     CALL CLASET('A',M-N,N,CZERO,CZERO,U(N+1,1),Ldu)
                     IF ( N<n1 ) THEN
                        CALL CLASET('A',N,n1-N,CZERO,CZERO,U(1,N+1),Ldu)
                        CALL CLASET('A',M-N,n1-N,CZERO,CONE,U(N+1,N+1), &
     &                              Ldu)
                     ENDIF
                  ENDIF
               ELSE
!                  .. copy R**H into [U] and overwrite [U] with the right
!                  singular vectors of R
                  DO p = 1 , nr
                     DO q = p , N
                        U(q,nr+p) = CONJG(A(p,q))
                     ENDDO
                  ENDDO
                  IF ( nr>1 ) CALL CLASET('U',nr-1,nr-1,CZERO,CZERO,    &
     &                 U(1,nr+2),Ldu)
                  CALL CGEQRF(N,nr,U(1,nr+1),Ldu,Cwork(N+1),            &
     &                        Cwork(N+nr+1),Lcwork-N-nr,ierr)
                  DO p = 1 , nr
                     DO q = 1 , N
                        V(q,p) = CONJG(U(p,nr+q))
                     ENDDO
                  ENDDO
                  CALL CLASET('U',nr-1,nr-1,CZERO,CZERO,V(1,2),Ldv)
                  CALL CGESVD('S','O',nr,nr,V,Ldv,S,U,Ldu,V,Ldv,        &
     &                        Cwork(N+nr+1),Lcwork-N-nr,Rwork,Info)
                  CALL CLASET('A',N-nr,nr,CZERO,CZERO,V(nr+1,1),Ldv)
                  CALL CLASET('A',nr,N-nr,CZERO,CZERO,V(1,nr+1),Ldv)
                  CALL CLASET('A',N-nr,N-nr,CZERO,CONE,V(nr+1,nr+1),Ldv)
                  CALL CUNMQR('R','C',N,N,nr,U(1,nr+1),Ldu,Cwork(N+1),V,&
     &                        Ldv,Cwork(N+nr+1),Lcwork-N-nr,ierr)
                  CALL CLAPMT(.FALSE.,N,N,V,Ldv,Iwork)
!                 .. assemble the left singular vector matrix U of dimensions
!                 (M x NR) or (M x N) or (M x M).
                  IF ( (nr<M) .AND. .NOT.(wntuf) ) THEN
                     CALL CLASET('A',M-nr,nr,CZERO,CZERO,U(nr+1,1),Ldu)
                     IF ( nr<n1 ) THEN
                        CALL CLASET('A',nr,n1-nr,CZERO,CZERO,U(1,nr+1), &
     &                              Ldu)
                        CALL CLASET('A',M-nr,n1-nr,CZERO,CONE,          &
     &                              U(nr+1,nr+1),Ldu)
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
!
!
!            .. apply CGESVD to R [[this is the recommended option]]
!
         ELSEIF ( wntvr .OR. (nr==N) ) THEN
!                .. copy R into [V] and overwrite V with the right singular vectors
            CALL CLACPY('U',nr,N,A,Lda,V,Ldv)
            IF ( nr>1 ) CALL CLASET('L',nr-1,nr-1,CZERO,CZERO,V(2,1),   &
     &                              Ldv)
!               .. the right singular vectors of R overwrite [V], the NR left
!               singular vectors of R stored in [U](1:NR,1:NR)
            CALL CGESVD('S','O',nr,N,V,Ldv,S,U,Ldu,V,Ldv,Cwork(N+1),    &
     &                  Lcwork-N,Rwork,Info)
            CALL CLAPMT(.FALSE.,nr,N,V,Ldv,Iwork)
!               .. now [V](1:NR,1:N) contains V(1:N,1:NR)**H
!               .. assemble the left singular vector matrix U of dimensions
!              (M x NR) or (M x N) or (M x M).
            IF ( (nr<M) .AND. .NOT.(wntuf) ) THEN
               CALL CLASET('A',M-nr,nr,CZERO,CZERO,U(nr+1,1),Ldu)
               IF ( nr<n1 ) THEN
                  CALL CLASET('A',nr,n1-nr,CZERO,CZERO,U(1,nr+1),Ldu)
                  CALL CLASET('A',M-nr,n1-nr,CZERO,CONE,U(nr+1,nr+1),   &
     &                        Ldu)
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
!               OPTRATIO = ILAENV(6, 'CGESVD', 'S' // 'O', NR,N,0,0)
!               OPTRATIO = MAX( OPTRATIO, 2 )
            optratio = 2
            IF ( optratio*nr>N ) THEN
               CALL CLACPY('U',nr,N,A,Lda,V,Ldv)
               IF ( nr>1 ) CALL CLASET('L',nr-1,nr-1,CZERO,CZERO,V(2,1),&
     &                                 Ldv)
!              .. the right singular vectors of R overwrite [V], the NR left
!                 singular vectors of R stored in [U](1:NR,1:NR)
               CALL CLASET('A',N-nr,N,CZERO,CZERO,V(nr+1,1),Ldv)
               CALL CGESVD('S','O',N,N,V,Ldv,S,U,Ldu,V,Ldv,Cwork(N+1),  &
     &                     Lcwork-N,Rwork,Info)
               CALL CLAPMT(.FALSE.,N,N,V,Ldv,Iwork)
!                 .. now [V] contains the adjoint of the matrix of the right
!                 singular vectors of A. The leading N left singular vectors
!                 are in [U](1:N,1:N)
!                 .. assemble the left singular vector matrix U of dimensions
!                 (M x N1), i.e. (M x N) or (M x M).
               IF ( (N<M) .AND. .NOT.(wntuf) ) THEN
                  CALL CLASET('A',M-N,N,CZERO,CZERO,U(N+1,1),Ldu)
                  IF ( N<n1 ) THEN
                     CALL CLASET('A',N,n1-N,CZERO,CZERO,U(1,N+1),Ldu)
                     CALL CLASET('A',M-N,n1-N,CZERO,CONE,U(N+1,N+1),Ldu)
                  ENDIF
               ENDIF
            ELSE
               CALL CLACPY('U',nr,N,A,Lda,U(nr+1,1),Ldu)
               IF ( nr>1 ) CALL CLASET('L',nr-1,nr-1,CZERO,CZERO,       &
     &                                 U(nr+2,1),Ldu)
               CALL CGELQF(nr,N,U(nr+1,1),Ldu,Cwork(N+1),Cwork(N+nr+1), &
     &                     Lcwork-N-nr,ierr)
               CALL CLACPY('L',nr,nr,U(nr+1,1),Ldu,V,Ldv)
               IF ( nr>1 ) CALL CLASET('U',nr-1,nr-1,CZERO,CZERO,V(1,2),&
     &                                 Ldv)
               CALL CGESVD('S','O',nr,nr,V,Ldv,S,U,Ldu,V,Ldv,           &
     &                     Cwork(N+nr+1),Lcwork-N-nr,Rwork,Info)
               CALL CLASET('A',N-nr,nr,CZERO,CZERO,V(nr+1,1),Ldv)
               CALL CLASET('A',nr,N-nr,CZERO,CZERO,V(1,nr+1),Ldv)
               CALL CLASET('A',N-nr,N-nr,CZERO,CONE,V(nr+1,nr+1),Ldv)
               CALL CUNMLQ('R','N',N,N,nr,U(nr+1,1),Ldu,Cwork(N+1),V,   &
     &                     Ldv,Cwork(N+nr+1),Lcwork-N-nr,ierr)
               CALL CLAPMT(.FALSE.,N,N,V,Ldv,Iwork)
!               .. assemble the left singular vector matrix U of dimensions
!              (M x NR) or (M x N) or (M x M).
               IF ( (nr<M) .AND. .NOT.(wntuf) ) THEN
                  CALL CLASET('A',M-nr,nr,CZERO,CZERO,U(nr+1,1),Ldu)
                  IF ( nr<n1 ) THEN
                     CALL CLASET('A',nr,n1-nr,CZERO,CZERO,U(1,nr+1),Ldu)
                     CALL CLASET('A',M-nr,n1-nr,CZERO,CONE,U(nr+1,nr+1),&
     &                           Ldu)
                  ENDIF
               ENDIF
            ENDIF
!        .. end of the "R**H or R" branch
         ENDIF
!
!           The Q matrix from the first QRF is built into the left singular
!           vectors matrix U.
!
         IF ( .NOT.wntuf ) CALL CUNMQR('L','N',M,n1,N,A,Lda,Cwork,U,Ldu,&
     &                                 Cwork(N+1),Lcwork-N,ierr)
         IF ( rowprm .AND. .NOT.wntuf )                                 &
     &        CALL CLASWP(n1,U,Ldu,1,M-1,Iwork(N+1),-1)
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
      IF ( nr<N ) CALL SLASET('G',N-nr,1,ZERO,ZERO,S(nr+1),N)
!     .. undo scaling; this may cause overflow in the largest singular
!     values.
      IF ( ascaled ) CALL SLASCL('G',0,0,ONE,SQRT(REAL(M)),nr,1,S,N,    &
     &                           ierr)
      IF ( conda ) Rwork(1) = sconda
      Rwork(2) = p - nr
!     .. p-NR is the number of singular values that are computed as
!     exact zeros in CGESVD() applied to the (possibly truncated)
!     full row rank triangular (trapezoidal) factor of A.
      Numrank = nr
!
!
!     End of CGESVDQ
!
      END SUBROUTINE CGESVDQ
