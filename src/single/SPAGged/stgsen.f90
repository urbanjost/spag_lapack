!*==stgsen.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b STGSEN
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download STGSEN + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/stgsen.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/stgsen.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/stgsen.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE STGSEN( IJOB, WANTQ, WANTZ, SELECT, N, A, LDA, B, LDB,
!                          ALPHAR, ALPHAI, BETA, Q, LDQ, Z, LDZ, M, PL,
!                          PR, DIF, WORK, LWORK, IWORK, LIWORK, INFO )
!
!       .. Scalar Arguments ..
!       LOGICAL            WANTQ, WANTZ
!       INTEGER            IJOB, INFO, LDA, LDB, LDQ, LDZ, LIWORK, LWORK,
!      $                   M, N
!       REAL               PL, PR
!       ..
!       .. Array Arguments ..
!       LOGICAL            SELECT( * )
!       INTEGER            IWORK( * )
!       REAL               A( LDA, * ), ALPHAI( * ), ALPHAR( * ),
!      $                   B( LDB, * ), BETA( * ), DIF( * ), Q( LDQ, * ),
!      $                   WORK( * ), Z( LDZ, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> STGSEN reorders the generalized real Schur decomposition of a real
!> matrix pair (A, B) (in terms of an orthonormal equivalence trans-
!> formation Q**T * (A, B) * Z), so that a selected cluster of eigenvalues
!> appears in the leading diagonal blocks of the upper quasi-triangular
!> matrix A and the upper triangular B. The leading columns of Q and
!> Z form orthonormal bases of the corresponding left and right eigen-
!> spaces (deflating subspaces). (A, B) must be in generalized real
!> Schur canonical form (as returned by SGGES), i.e. A is block upper
!> triangular with 1-by-1 and 2-by-2 diagonal blocks. B is upper
!> triangular.
!>
!> STGSEN also computes the generalized eigenvalues
!>
!>             w(j) = (ALPHAR(j) + i*ALPHAI(j))/BETA(j)
!>
!> of the reordered matrix pair (A, B).
!>
!> Optionally, STGSEN computes the estimates of reciprocal condition
!> numbers for eigenvalues and eigenspaces. These are Difu[(A11,B11),
!> (A22,B22)] and Difl[(A11,B11), (A22,B22)], i.e. the separation(s)
!> between the matrix pairs (A11, B11) and (A22,B22) that correspond to
!> the selected cluster and the eigenvalues outside the cluster, resp.,
!> and norms of "projections" onto left and right eigenspaces w.r.t.
!> the selected cluster in the (1,1)-block.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] IJOB
!> \verbatim
!>          IJOB is INTEGER
!>          Specifies whether condition numbers are required for the
!>          cluster of eigenvalues (PL and PR) or the deflating subspaces
!>          (Difu and Difl):
!>           =0: Only reorder w.r.t. SELECT. No extras.
!>           =1: Reciprocal of norms of "projections" onto left and right
!>               eigenspaces w.r.t. the selected cluster (PL and PR).
!>           =2: Upper bounds on Difu and Difl. F-norm-based estimate
!>               (DIF(1:2)).
!>           =3: Estimate of Difu and Difl. 1-norm-based estimate
!>               (DIF(1:2)).
!>               About 5 times as expensive as IJOB = 2.
!>           =4: Compute PL, PR and DIF (i.e. 0, 1 and 2 above): Economic
!>               version to get it all.
!>           =5: Compute PL, PR and DIF (i.e. 0, 1 and 3 above)
!> \endverbatim
!>
!> \param[in] WANTQ
!> \verbatim
!>          WANTQ is LOGICAL
!>          .TRUE. : update the left transformation matrix Q;
!>          .FALSE.: do not update Q.
!> \endverbatim
!>
!> \param[in] WANTZ
!> \verbatim
!>          WANTZ is LOGICAL
!>          .TRUE. : update the right transformation matrix Z;
!>          .FALSE.: do not update Z.
!> \endverbatim
!>
!> \param[in] SELECT
!> \verbatim
!>          SELECT is LOGICAL array, dimension (N)
!>          SELECT specifies the eigenvalues in the selected cluster.
!>          To select a real eigenvalue w(j), SELECT(j) must be set to
!>          .TRUE.. To select a complex conjugate pair of eigenvalues
!>          w(j) and w(j+1), corresponding to a 2-by-2 diagonal block,
!>          either SELECT(j) or SELECT(j+1) or both must be set to
!>          .TRUE.; a complex conjugate pair of eigenvalues must be
!>          either both included in the cluster or both excluded.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrices A and B. N >= 0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is REAL array, dimension(LDA,N)
!>          On entry, the upper quasi-triangular matrix A, with (A, B) in
!>          generalized real Schur canonical form.
!>          On exit, A is overwritten by the reordered matrix A.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A. LDA >= max(1,N).
!> \endverbatim
!>
!> \param[in,out] B
!> \verbatim
!>          B is REAL array, dimension(LDB,N)
!>          On entry, the upper triangular matrix B, with (A, B) in
!>          generalized real Schur canonical form.
!>          On exit, B is overwritten by the reordered matrix B.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of the array B. LDB >= max(1,N).
!> \endverbatim
!>
!> \param[out] ALPHAR
!> \verbatim
!>          ALPHAR is REAL array, dimension (N)
!> \endverbatim
!>
!> \param[out] ALPHAI
!> \verbatim
!>          ALPHAI is REAL array, dimension (N)
!> \endverbatim
!>
!> \param[out] BETA
!> \verbatim
!>          BETA is REAL array, dimension (N)
!>
!>          On exit, (ALPHAR(j) + ALPHAI(j)*i)/BETA(j), j=1,...,N, will
!>          be the generalized eigenvalues.  ALPHAR(j) + ALPHAI(j)*i
!>          and BETA(j),j=1,...,N  are the diagonals of the complex Schur
!>          form (S,T) that would result if the 2-by-2 diagonal blocks of
!>          the real generalized Schur form of (A,B) were further reduced
!>          to triangular form using complex unitary transformations.
!>          If ALPHAI(j) is zero, then the j-th eigenvalue is real; if
!>          positive, then the j-th and (j+1)-st eigenvalues are a
!>          complex conjugate pair, with ALPHAI(j+1) negative.
!> \endverbatim
!>
!> \param[in,out] Q
!> \verbatim
!>          Q is REAL array, dimension (LDQ,N)
!>          On entry, if WANTQ = .TRUE., Q is an N-by-N matrix.
!>          On exit, Q has been postmultiplied by the left orthogonal
!>          transformation matrix which reorder (A, B); The leading M
!>          columns of Q form orthonormal bases for the specified pair of
!>          left eigenspaces (deflating subspaces).
!>          If WANTQ = .FALSE., Q is not referenced.
!> \endverbatim
!>
!> \param[in] LDQ
!> \verbatim
!>          LDQ is INTEGER
!>          The leading dimension of the array Q.  LDQ >= 1;
!>          and if WANTQ = .TRUE., LDQ >= N.
!> \endverbatim
!>
!> \param[in,out] Z
!> \verbatim
!>          Z is REAL array, dimension (LDZ,N)
!>          On entry, if WANTZ = .TRUE., Z is an N-by-N matrix.
!>          On exit, Z has been postmultiplied by the left orthogonal
!>          transformation matrix which reorder (A, B); The leading M
!>          columns of Z form orthonormal bases for the specified pair of
!>          left eigenspaces (deflating subspaces).
!>          If WANTZ = .FALSE., Z is not referenced.
!> \endverbatim
!>
!> \param[in] LDZ
!> \verbatim
!>          LDZ is INTEGER
!>          The leading dimension of the array Z. LDZ >= 1;
!>          If WANTZ = .TRUE., LDZ >= N.
!> \endverbatim
!>
!> \param[out] M
!> \verbatim
!>          M is INTEGER
!>          The dimension of the specified pair of left and right eigen-
!>          spaces (deflating subspaces). 0 <= M <= N.
!> \endverbatim
!>
!> \param[out] PL
!> \verbatim
!>          PL is REAL
!> \endverbatim
!>
!> \param[out] PR
!> \verbatim
!>          PR is REAL
!>
!>          If IJOB = 1, 4 or 5, PL, PR are lower bounds on the
!>          reciprocal of the norm of "projections" onto left and right
!>          eigenspaces with respect to the selected cluster.
!>          0 < PL, PR <= 1.
!>          If M = 0 or M = N, PL = PR  = 1.
!>          If IJOB = 0, 2 or 3, PL and PR are not referenced.
!> \endverbatim
!>
!> \param[out] DIF
!> \verbatim
!>          DIF is REAL array, dimension (2).
!>          If IJOB >= 2, DIF(1:2) store the estimates of Difu and Difl.
!>          If IJOB = 2 or 4, DIF(1:2) are F-norm-based upper bounds on
!>          Difu and Difl. If IJOB = 3 or 5, DIF(1:2) are 1-norm-based
!>          estimates of Difu and Difl.
!>          If M = 0 or N, DIF(1:2) = F-norm([A, B]).
!>          If IJOB = 0 or 1, DIF is not referenced.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is REAL array, dimension (MAX(1,LWORK))
!>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The dimension of the array WORK. LWORK >=  4*N+16.
!>          If IJOB = 1, 2 or 4, LWORK >= MAX(4*N+16, 2*M*(N-M)).
!>          If IJOB = 3 or 5, LWORK >= MAX(4*N+16, 4*M*(N-M)).
!>
!>          If LWORK = -1, then a workspace query is assumed; the routine
!>          only calculates the optimal size of the WORK array, returns
!>          this value as the first entry of the WORK array, and no error
!>          message related to LWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] IWORK
!> \verbatim
!>          IWORK is INTEGER array, dimension (MAX(1,LIWORK))
!>          On exit, if INFO = 0, IWORK(1) returns the optimal LIWORK.
!> \endverbatim
!>
!> \param[in] LIWORK
!> \verbatim
!>          LIWORK is INTEGER
!>          The dimension of the array IWORK. LIWORK >= 1.
!>          If IJOB = 1, 2 or 4, LIWORK >=  N+6.
!>          If IJOB = 3 or 5, LIWORK >= MAX(2*M*(N-M), N+6).
!>
!>          If LIWORK = -1, then a workspace query is assumed; the
!>          routine only calculates the optimal size of the IWORK array,
!>          returns this value as the first entry of the IWORK array, and
!>          no error message related to LIWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>            =0: Successful exit.
!>            <0: If INFO = -i, the i-th argument had an illegal value.
!>            =1: Reordering of (A, B) failed because the transformed
!>                matrix pair (A, B) would be too far from generalized
!>                Schur form; the problem is very ill-conditioned.
!>                (A, B) may have been partially reordered.
!>                If requested, 0 is returned in DIF(*), PL and PR.
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
!> \date June 2016
!
!> \ingroup realOTHERcomputational
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  STGSEN first collects the selected eigenvalues by computing
!>  orthogonal U and W that move them to the top left corner of (A, B).
!>  In other words, the selected eigenvalues are the eigenvalues of
!>  (A11, B11) in:
!>
!>              U**T*(A, B)*W = (A11 A12) (B11 B12) n1
!>                              ( 0  A22),( 0  B22) n2
!>                                n1  n2    n1  n2
!>
!>  where N = n1+n2 and U**T means the transpose of U. The first n1 columns
!>  of U and W span the specified pair of left and right eigenspaces
!>  (deflating subspaces) of (A, B).
!>
!>  If (A, B) has been obtained from the generalized real Schur
!>  decomposition of a matrix pair (C, D) = Q*(A, B)*Z**T, then the
!>  reordered generalized real Schur form of (C, D) is given by
!>
!>           (C, D) = (Q*U)*(U**T*(A, B)*W)*(Z*W)**T,
!>
!>  and the first n1 columns of Q*U and Z*W span the corresponding
!>  deflating subspaces of (C, D) (Q and Z store Q*U and Z*W, resp.).
!>
!>  Note that if the selected eigenvalue is sufficiently ill-conditioned,
!>  then its value may differ significantly from its value before
!>  reordering.
!>
!>  The reciprocal condition numbers of the left and right eigenspaces
!>  spanned by the first n1 columns of U and W (or Q*U and Z*W) may
!>  be returned in DIF(1:2), corresponding to Difu and Difl, resp.
!>
!>  The Difu and Difl are defined as:
!>
!>       Difu[(A11, B11), (A22, B22)] = sigma-min( Zu )
!>  and
!>       Difl[(A11, B11), (A22, B22)] = Difu[(A22, B22), (A11, B11)],
!>
!>  where sigma-min(Zu) is the smallest singular value of the
!>  (2*n1*n2)-by-(2*n1*n2) matrix
!>
!>       Zu = [ kron(In2, A11)  -kron(A22**T, In1) ]
!>            [ kron(In2, B11)  -kron(B22**T, In1) ].
!>
!>  Here, Inx is the identity matrix of size nx and A22**T is the
!>  transpose of A22. kron(X, Y) is the Kronecker product between
!>  the matrices X and Y.
!>
!>  When DIF(2) is small, small changes in (A, B) can cause large changes
!>  in the deflating subspace. An approximate (asymptotic) bound on the
!>  maximum angular error in the computed deflating subspaces is
!>
!>       EPS * norm((A, B)) / DIF(2),
!>
!>  where EPS is the machine precision.
!>
!>  The reciprocal norm of the projectors on the left and right
!>  eigenspaces associated with (A11, B11) may be returned in PL and PR.
!>  They are computed as follows. First we compute L and R so that
!>  P*(A, B)*Q is block diagonal, where
!>
!>       P = ( I -L ) n1           Q = ( I R ) n1
!>           ( 0  I ) n2    and        ( 0 I ) n2
!>             n1 n2                    n1 n2
!>
!>  and (L, R) is the solution to the generalized Sylvester equation
!>
!>       A11*R - L*A22 = -A12
!>       B11*R - L*B22 = -B12
!>
!>  Then PL = (F-norm(L)**2+1)**(-1/2) and PR = (F-norm(R)**2+1)**(-1/2).
!>  An approximate (asymptotic) bound on the average absolute error of
!>  the selected eigenvalues is
!>
!>       EPS * norm((A, B)) / PL.
!>
!>  There are also global error bounds which valid for perturbations up
!>  to a certain restriction:  A lower bound (x) on the smallest
!>  F-norm(E,F) for which an eigenvalue of (A11, B11) may move and
!>  coalesce with an eigenvalue of (A22, B22) under perturbation (E,F),
!>  (i.e. (A + E, B + F), is
!>
!>   x = min(Difu,Difl)/((1/(PL*PL)+1/(PR*PR))**(1/2)+2*max(1/PL,1/PR)).
!>
!>  An approximate bound on x can be computed from DIF(1:2), PL and PR.
!>
!>  If y = ( F-norm(E,F) / x) <= 1, the angles between the perturbed
!>  (L', R') and unperturbed (L, R) left and right deflating subspaces
!>  associated with the selected cluster in the (1,1)-blocks can be
!>  bounded as
!>
!>   max-angle(L, L') <= arctan( y * PL / (1 - y * (1 - PL * PL)**(1/2))
!>   max-angle(R, R') <= arctan( y * PR / (1 - y * (1 - PR * PR)**(1/2))
!>
!>  See LAPACK User's Guide section 4.11 or the following references
!>  for more information.
!>
!>  Note that if the default method for computing the Frobenius-norm-
!>  based estimate DIF is not wanted (see SLATDF), then the parameter
!>  IDIFJB (see below) should be changed from 3 to 4 (routine SLATDF
!>  (IJOB = 2 will be used)). See STGSYL for more details.
!> \endverbatim
!
!> \par Contributors:
!  ==================
!>
!>     Bo Kagstrom and Peter Poromaa, Department of Computing Science,
!>     Umea University, S-901 87 Umea, Sweden.
!
!> \par References:
!  ================
!>
!> \verbatim
!>
!>  [1] B. Kagstrom; A Direct Method for Reordering Eigenvalues in the
!>      Generalized Real Schur Form of a Regular Matrix Pair (A, B), in
!>      M.S. Moonen et al (eds), Linear Algebra for Large Scale and
!>      Real-Time Applications, Kluwer Academic Publ. 1993, pp 195-218.
!>
!>  [2] B. Kagstrom and P. Poromaa; Computing Eigenspaces with Specified
!>      Eigenvalues of a Regular Matrix Pair (A, B) and Condition
!>      Estimation: Theory, Algorithms and Software,
!>      Report UMINF - 94.04, Department of Computing Science, Umea
!>      University, S-901 87 Umea, Sweden, 1994. Also as LAPACK Working
!>      Note 87. To appear in Numerical Algorithms, 1996.
!>
!>  [3] B. Kagstrom and P. Poromaa, LAPACK-Style Algorithms and Software
!>      for Solving the Generalized Sylvester Equation and Estimating the
!>      Separation between Regular Matrix Pairs, Report UMINF - 93.23,
!>      Department of Computing Science, Umea University, S-901 87 Umea,
!>      Sweden, December 1993, Revised April 1994, Also as LAPACK Working
!>      Note 75. To appear in ACM Trans. on Math. Software, Vol 22, No 1,
!>      1996.
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE STGSEN(Ijob,Wantq,Wantz,Select,N,A,Lda,B,Ldb,Alphar,   &
     &                  Alphai,Beta,Q,Ldq,Z,Ldz,M,Pl,Pr,Dif,Work,Lwork, &
     &                  Iwork,Liwork,Info)
      USE S_SLACN2
      USE S_SLACPY
      USE S_SLAG2
      USE S_SLAMCH
      USE S_SLASSQ
      USE S_STGEXC
      USE S_STGSYL
      USE S_XERBLA
      IMPLICIT NONE
!*--STGSEN463
!
! PARAMETER definitions rewritten by SPAG
!
      INTEGER , PARAMETER  ::  IDIFJB = 3
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER , INTENT(IN) :: Ijob
      LOGICAL :: Wantq
      LOGICAL :: Wantz
      LOGICAL , INTENT(IN) , DIMENSION(*) :: Select
      INTEGER :: N
      REAL , INTENT(INOUT) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL , DIMENSION(*) :: Alphar
      REAL , INTENT(INOUT) , DIMENSION(*) :: Alphai
      REAL , DIMENSION(*) :: Beta
      REAL , INTENT(INOUT) , DIMENSION(Ldq,*) :: Q
      INTEGER :: Ldq
      REAL , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      INTEGER , INTENT(INOUT) :: M
      REAL , INTENT(INOUT) :: Pl
      REAL , INTENT(INOUT) :: Pr
      REAL , INTENT(INOUT) , DIMENSION(*) :: Dif
      REAL , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(IN) :: Liwork
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      REAL :: dscale , dsum , eps , rdscal , smlnum
      INTEGER :: i , ierr , ijb , k , kase , kk , ks , liwmin , lwmin , &
     &           mn2 , n1 , n2
      INTEGER , DIMENSION(3) :: isave
      LOGICAL :: lquery , pair , swap , wantd , wantd1 , wantd2 , wantp
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
!     .. Local Arrays ..
!     ..
!     .. External Subroutines ..
!     ..
!     .. External Functions ..
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. Executable Statements ..
!
!     Decode and test the input parameters
!
      Info = 0
      lquery = (Lwork==-1 .OR. Liwork==-1)
!
      IF ( Ijob<0 .OR. Ijob>5 ) THEN
         Info = -1
      ELSEIF ( N<0 ) THEN
         Info = -5
      ELSEIF ( Lda<MAX(1,N) ) THEN
         Info = -7
      ELSEIF ( Ldb<MAX(1,N) ) THEN
         Info = -9
      ELSEIF ( Ldq<1 .OR. (Wantq .AND. Ldq<N) ) THEN
         Info = -14
      ELSEIF ( Ldz<1 .OR. (Wantz .AND. Ldz<N) ) THEN
         Info = -16
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('STGSEN',-Info)
         RETURN
      ENDIF
!
!     Get machine constants
!
      eps = SLAMCH('P')
      smlnum = SLAMCH('S')/eps
      ierr = 0
!
      wantp = Ijob==1 .OR. Ijob>=4
      wantd1 = Ijob==2 .OR. Ijob==4
      wantd2 = Ijob==3 .OR. Ijob==5
      wantd = wantd1 .OR. wantd2
!
!     Set M to the dimension of the specified pair of deflating
!     subspaces.
!
      M = 0
      pair = .FALSE.
      IF ( .NOT.lquery .OR. Ijob/=0 ) THEN
         DO k = 1 , N
            IF ( pair ) THEN
               pair = .FALSE.
            ELSEIF ( k<N ) THEN
               IF ( A(k+1,k)==ZERO ) THEN
                  IF ( Select(k) ) M = M + 1
               ELSE
                  pair = .TRUE.
                  IF ( Select(k) .OR. Select(k+1) ) M = M + 2
               ENDIF
            ELSE
               IF ( Select(N) ) M = M + 1
            ENDIF
         ENDDO
      ENDIF
!
      IF ( Ijob==1 .OR. Ijob==2 .OR. Ijob==4 ) THEN
         lwmin = MAX(1,4*N+16,2*M*(N-M))
         liwmin = MAX(1,N+6)
      ELSEIF ( Ijob==3 .OR. Ijob==5 ) THEN
         lwmin = MAX(1,4*N+16,4*M*(N-M))
         liwmin = MAX(1,2*M*(N-M),N+6)
      ELSE
         lwmin = MAX(1,4*N+16)
         liwmin = 1
      ENDIF
!
      Work(1) = lwmin
      Iwork(1) = liwmin
!
      IF ( Lwork<lwmin .AND. .NOT.lquery ) THEN
         Info = -22
      ELSEIF ( Liwork<liwmin .AND. .NOT.lquery ) THEN
         Info = -24
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('STGSEN',-Info)
         RETURN
      ELSEIF ( lquery ) THEN
         RETURN
      ENDIF
!
!     Quick return if possible.
!
      IF ( M==N .OR. M==0 ) THEN
         IF ( wantp ) THEN
            Pl = ONE
            Pr = ONE
         ENDIF
         IF ( wantd ) THEN
            dscale = ZERO
            dsum = ONE
            DO i = 1 , N
               CALL SLASSQ(N,A(1,i),1,dscale,dsum)
               CALL SLASSQ(N,B(1,i),1,dscale,dsum)
            ENDDO
            Dif(1) = dscale*SQRT(dsum)
            Dif(2) = Dif(1)
         ENDIF
         GOTO 100
      ENDIF
!
!     Collect the selected blocks at the top-left corner of (A, B).
!
      ks = 0
      pair = .FALSE.
      DO k = 1 , N
         IF ( pair ) THEN
            pair = .FALSE.
         ELSE
!
            swap = Select(k)
            IF ( k<N ) THEN
               IF ( A(k+1,k)/=ZERO ) THEN
                  pair = .TRUE.
                  swap = swap .OR. Select(k+1)
               ENDIF
            ENDIF
!
            IF ( swap ) THEN
               ks = ks + 1
!
!              Swap the K-th block to position KS.
!              Perform the reordering of diagonal blocks in (A, B)
!              by orthogonal transformation matrices and update
!              Q and Z accordingly (if requested):
!
               kk = k
               IF ( k/=ks ) CALL STGEXC(Wantq,Wantz,N,A,Lda,B,Ldb,Q,Ldq,&
     &                                  Z,Ldz,kk,ks,Work,Lwork,ierr)
!
               IF ( ierr>0 ) THEN
!
!                 Swap is rejected: exit.
!
                  Info = 1
                  IF ( wantp ) THEN
                     Pl = ZERO
                     Pr = ZERO
                  ENDIF
                  IF ( wantd ) THEN
                     Dif(1) = ZERO
                     Dif(2) = ZERO
                  ENDIF
                  GOTO 100
               ENDIF
!
               IF ( pair ) ks = ks + 1
            ENDIF
         ENDIF
      ENDDO
      IF ( wantp ) THEN
!
!        Solve generalized Sylvester equation for R and L
!        and compute PL and PR.
!
         n1 = M
         n2 = N - M
         i = n1 + 1
         ijb = 0
         CALL SLACPY('Full',n1,n2,A(1,i),Lda,Work,n1)
         CALL SLACPY('Full',n1,n2,B(1,i),Ldb,Work(n1*n2+1),n1)
         CALL STGSYL('N',ijb,n1,n2,A,Lda,A(i,i),Lda,Work,n1,B,Ldb,B(i,i)&
     &               ,Ldb,Work(n1*n2+1),n1,dscale,Dif(1),Work(n1*n2*2+1)&
     &               ,Lwork-2*n1*n2,Iwork,ierr)
!
!        Estimate the reciprocal of norms of "projections" onto left
!        and right eigenspaces.
!
         rdscal = ZERO
         dsum = ONE
         CALL SLASSQ(n1*n2,Work,1,rdscal,dsum)
         Pl = rdscal*SQRT(dsum)
         IF ( Pl==ZERO ) THEN
            Pl = ONE
         ELSE
            Pl = dscale/(SQRT(dscale*dscale/Pl+Pl)*SQRT(Pl))
         ENDIF
         rdscal = ZERO
         dsum = ONE
         CALL SLASSQ(n1*n2,Work(n1*n2+1),1,rdscal,dsum)
         Pr = rdscal*SQRT(dsum)
         IF ( Pr==ZERO ) THEN
            Pr = ONE
         ELSE
            Pr = dscale/(SQRT(dscale*dscale/Pr+Pr)*SQRT(Pr))
         ENDIF
      ENDIF
!
      IF ( wantd ) THEN
!
!        Compute estimates of Difu and Difl.
!
         IF ( wantd1 ) THEN
            n1 = M
            n2 = N - M
            i = n1 + 1
            ijb = IDIFJB
!
!           Frobenius norm-based Difu-estimate.
!
            CALL STGSYL('N',ijb,n1,n2,A,Lda,A(i,i),Lda,Work,n1,B,Ldb,   &
     &                  B(i,i),Ldb,Work(n1*n2+1),n1,dscale,Dif(1),      &
     &                  Work(2*n1*n2+1),Lwork-2*n1*n2,Iwork,ierr)
!
!           Frobenius norm-based Difl-estimate.
!
            CALL STGSYL('N',ijb,n2,n1,A(i,i),Lda,A,Lda,Work,n2,B(i,i),  &
     &                  Ldb,B,Ldb,Work(n1*n2+1),n2,dscale,Dif(2),       &
     &                  Work(2*n1*n2+1),Lwork-2*n1*n2,Iwork,ierr)
         ELSE
!
!
!           Compute 1-norm-based estimates of Difu and Difl using
!           reversed communication with SLACN2. In each step a
!           generalized Sylvester equation or a transposed variant
!           is solved.
!
            kase = 0
            n1 = M
            n2 = N - M
            i = n1 + 1
            ijb = 0
            mn2 = 2*n1*n2
            DO
!
!           1-norm-based estimate of Difu.
!
               CALL SLACN2(mn2,Work(mn2+1),Work,Iwork,Dif(1),kase,isave)
               IF ( kase/=0 ) THEN
                  IF ( kase==1 ) THEN
!
!                 Solve generalized Sylvester equation.
!
                     CALL STGSYL('N',ijb,n1,n2,A,Lda,A(i,i),Lda,Work,n1,&
     &                           B,Ldb,B(i,i),Ldb,Work(n1*n2+1),n1,     &
     &                           dscale,Dif(1),Work(2*n1*n2+1),         &
     &                           Lwork-2*n1*n2,Iwork,ierr)
                  ELSE
!
!                 Solve the transposed variant.
!
                     CALL STGSYL('T',ijb,n1,n2,A,Lda,A(i,i),Lda,Work,n1,&
     &                           B,Ldb,B(i,i),Ldb,Work(n1*n2+1),n1,     &
     &                           dscale,Dif(1),Work(2*n1*n2+1),         &
     &                           Lwork-2*n1*n2,Iwork,ierr)
                  ENDIF
                  CYCLE
               ENDIF
               Dif(1) = dscale/Dif(1)
               EXIT
            ENDDO
            DO
!
!           1-norm-based estimate of Difl.
!
               CALL SLACN2(mn2,Work(mn2+1),Work,Iwork,Dif(2),kase,isave)
               IF ( kase/=0 ) THEN
                  IF ( kase==1 ) THEN
!
!                 Solve generalized Sylvester equation.
!
                     CALL STGSYL('N',ijb,n2,n1,A(i,i),Lda,A,Lda,Work,n2,&
     &                           B(i,i),Ldb,B,Ldb,Work(n1*n2+1),n2,     &
     &                           dscale,Dif(2),Work(2*n1*n2+1),         &
     &                           Lwork-2*n1*n2,Iwork,ierr)
                  ELSE
!
!                 Solve the transposed variant.
!
                     CALL STGSYL('T',ijb,n2,n1,A(i,i),Lda,A,Lda,Work,n2,&
     &                           B(i,i),Ldb,B,Ldb,Work(n1*n2+1),n2,     &
     &                           dscale,Dif(2),Work(2*n1*n2+1),         &
     &                           Lwork-2*n1*n2,Iwork,ierr)
                  ENDIF
                  CYCLE
               ENDIF
               Dif(2) = dscale/Dif(2)
               EXIT
            ENDDO
!
         ENDIF
      ENDIF
!
!
!     Compute generalized eigenvalues of reordered pair (A, B) and
!     normalize the generalized Schur form.
!
 100  pair = .FALSE.
      DO k = 1 , N
         IF ( pair ) THEN
            pair = .FALSE.
         ELSE
!
            IF ( k<N ) THEN
               IF ( A(k+1,k)/=ZERO ) pair = .TRUE.
            ENDIF
!
            IF ( pair ) THEN
!
!             Compute the eigenvalue(s) at position K.
!
               Work(1) = A(k,k)
               Work(2) = A(k+1,k)
               Work(3) = A(k,k+1)
               Work(4) = A(k+1,k+1)
               Work(5) = B(k,k)
               Work(6) = B(k+1,k)
               Work(7) = B(k,k+1)
               Work(8) = B(k+1,k+1)
               CALL SLAG2(Work,2,Work(5),2,smlnum*eps,Beta(k),Beta(k+1),&
     &                    Alphar(k),Alphar(k+1),Alphai(k))
               Alphai(k+1) = -Alphai(k)
!
            ELSE
!
               IF ( SIGN(ONE,B(k,k))<ZERO ) THEN
!
!                 If B(K,K) is negative, make it positive
!
                  DO i = 1 , N
                     A(k,i) = -A(k,i)
                     B(k,i) = -B(k,i)
                     IF ( Wantq ) Q(i,k) = -Q(i,k)
                  ENDDO
               ENDIF
!
               Alphar(k) = A(k,k)
               Alphai(k) = ZERO
               Beta(k) = B(k,k)
!
            ENDIF
         ENDIF
      ENDDO
!
      Work(1) = lwmin
      Iwork(1) = liwmin
!
!
!     End of STGSEN
!
      END SUBROUTINE STGSEN
