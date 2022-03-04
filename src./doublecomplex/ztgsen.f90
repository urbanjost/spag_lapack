!*==ztgsen.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b ZTGSEN
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZTGSEN + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ztgsen.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ztgsen.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ztgsen.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZTGSEN( IJOB, WANTQ, WANTZ, SELECT, N, A, LDA, B, LDB,
!                          ALPHA, BETA, Q, LDQ, Z, LDZ, M, PL, PR, DIF,
!                          WORK, LWORK, IWORK, LIWORK, INFO )
!
!       .. Scalar Arguments ..
!       LOGICAL            WANTQ, WANTZ
!       INTEGER            IJOB, INFO, LDA, LDB, LDQ, LDZ, LIWORK, LWORK,
!      $                   M, N
!       DOUBLE PRECISION   PL, PR
!       ..
!       .. Array Arguments ..
!       LOGICAL            SELECT( * )
!       INTEGER            IWORK( * )
!       DOUBLE PRECISION   DIF( * )
!       COMPLEX*16         A( LDA, * ), ALPHA( * ), B( LDB, * ),
!      $                   BETA( * ), Q( LDQ, * ), WORK( * ), Z( LDZ, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZTGSEN reorders the generalized Schur decomposition of a complex
!> matrix pair (A, B) (in terms of an unitary equivalence trans-
!> formation Q**H * (A, B) * Z), so that a selected cluster of eigenvalues
!> appears in the leading diagonal blocks of the pair (A,B). The leading
!> columns of Q and Z form unitary bases of the corresponding left and
!> right eigenspaces (deflating subspaces). (A, B) must be in
!> generalized Schur canonical form, that is, A and B are both upper
!> triangular.
!>
!> ZTGSEN also computes the generalized eigenvalues
!>
!>          w(j)= ALPHA(j) / BETA(j)
!>
!> of the reordered matrix pair (A, B).
!>
!> Optionally, the routine computes estimates of reciprocal condition
!> numbers for eigenvalues and eigenspaces. These are Difu[(A11,B11),
!> (A22,B22)] and Difl[(A11,B11), (A22,B22)], i.e. the separation(s)
!> between the matrix pairs (A11, B11) and (A22,B22) that correspond to
!> the selected cluster and the eigenvalues outside the cluster, resp.,
!> and norms of "projections" onto left and right eigenspaces w.r.t.
!> the selected cluster in the (1,1)-block.
!>
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
!>          SELECT specifies the eigenvalues in the selected cluster. To
!>          select an eigenvalue w(j), SELECT(j) must be set to
!>          .TRUE..
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
!>          A is COMPLEX*16 array, dimension(LDA,N)
!>          On entry, the upper triangular matrix A, in generalized
!>          Schur canonical form.
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
!>          B is COMPLEX*16 array, dimension(LDB,N)
!>          On entry, the upper triangular matrix B, in generalized
!>          Schur canonical form.
!>          On exit, B is overwritten by the reordered matrix B.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of the array B. LDB >= max(1,N).
!> \endverbatim
!>
!> \param[out] ALPHA
!> \verbatim
!>          ALPHA is COMPLEX*16 array, dimension (N)
!> \endverbatim
!>
!> \param[out] BETA
!> \verbatim
!>          BETA is COMPLEX*16 array, dimension (N)
!>
!>          The diagonal elements of A and B, respectively,
!>          when the pair (A,B) has been reduced to generalized Schur
!>          form.  ALPHA(i)/BETA(i) i=1,...,N are the generalized
!>          eigenvalues.
!> \endverbatim
!>
!> \param[in,out] Q
!> \verbatim
!>          Q is COMPLEX*16 array, dimension (LDQ,N)
!>          On entry, if WANTQ = .TRUE., Q is an N-by-N matrix.
!>          On exit, Q has been postmultiplied by the left unitary
!>          transformation matrix which reorder (A, B); The leading M
!>          columns of Q form orthonormal bases for the specified pair of
!>          left eigenspaces (deflating subspaces).
!>          If WANTQ = .FALSE., Q is not referenced.
!> \endverbatim
!>
!> \param[in] LDQ
!> \verbatim
!>          LDQ is INTEGER
!>          The leading dimension of the array Q. LDQ >= 1.
!>          If WANTQ = .TRUE., LDQ >= N.
!> \endverbatim
!>
!> \param[in,out] Z
!> \verbatim
!>          Z is COMPLEX*16 array, dimension (LDZ,N)
!>          On entry, if WANTZ = .TRUE., Z is an N-by-N matrix.
!>          On exit, Z has been postmultiplied by the left unitary
!>          transformation matrix which reorder (A, B); The leading M
!>          columns of Z form orthonormal bases for the specified pair of
!>          left eigenspaces (deflating subspaces).
!>          If WANTZ = .FALSE., Z is not referenced.
!> \endverbatim
!>
!> \param[in] LDZ
!> \verbatim
!>          LDZ is INTEGER
!>          The leading dimension of the array Z. LDZ >= 1.
!>          If WANTZ = .TRUE., LDZ >= N.
!> \endverbatim
!>
!> \param[out] M
!> \verbatim
!>          M is INTEGER
!>          The dimension of the specified pair of left and right
!>          eigenspaces, (deflating subspaces) 0 <= M <= N.
!> \endverbatim
!>
!> \param[out] PL
!> \verbatim
!>          PL is DOUBLE PRECISION
!> \endverbatim
!>
!> \param[out] PR
!> \verbatim
!>          PR is DOUBLE PRECISION
!>
!>          If IJOB = 1, 4 or 5, PL, PR are lower bounds on the
!>          reciprocal  of the norm of "projections" onto left and right
!>          eigenspace with respect to the selected cluster.
!>          0 < PL, PR <= 1.
!>          If M = 0 or M = N, PL = PR  = 1.
!>          If IJOB = 0, 2 or 3 PL, PR are not referenced.
!> \endverbatim
!>
!> \param[out] DIF
!> \verbatim
!>          DIF is DOUBLE PRECISION array, dimension (2).
!>          If IJOB >= 2, DIF(1:2) store the estimates of Difu and Difl.
!>          If IJOB = 2 or 4, DIF(1:2) are F-norm-based upper bounds on
!>          Difu and Difl. If IJOB = 3 or 5, DIF(1:2) are 1-norm-based
!>          estimates of Difu and Difl, computed using reversed
!>          communication with ZLACN2.
!>          If M = 0 or N, DIF(1:2) = F-norm([A, B]).
!>          If IJOB = 0 or 1, DIF is not referenced.
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
!>          The dimension of the array WORK. LWORK >=  1
!>          If IJOB = 1, 2 or 4, LWORK >=  2*M*(N-M)
!>          If IJOB = 3 or 5, LWORK >=  4*M*(N-M)
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
!>          If IJOB = 1, 2 or 4, LIWORK >=  N+2;
!>          If IJOB = 3 or 5, LIWORK >= MAX(N+2, 2*M*(N-M));
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
!> \ingroup complex16OTHERcomputational
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  ZTGSEN first collects the selected eigenvalues by computing unitary
!>  U and W that move them to the top left corner of (A, B). In other
!>  words, the selected eigenvalues are the eigenvalues of (A11, B11) in
!>
!>              U**H*(A, B)*W = (A11 A12) (B11 B12) n1
!>                              ( 0  A22),( 0  B22) n2
!>                                n1  n2    n1  n2
!>
!>  where N = n1+n2 and U**H means the conjugate transpose of U. The first
!>  n1 columns of U and W span the specified pair of left and right
!>  eigenspaces (deflating subspaces) of (A, B).
!>
!>  If (A, B) has been obtained from the generalized real Schur
!>  decomposition of a matrix pair (C, D) = Q*(A, B)*Z**H, then the
!>  reordered generalized Schur form of (C, D) is given by
!>
!>           (C, D) = (Q*U)*(U**H *(A, B)*W)*(Z*W)**H,
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
!>       Zu = [ kron(In2, A11)  -kron(A22**H, In1) ]
!>            [ kron(In2, B11)  -kron(B22**H, In1) ].
!>
!>  Here, Inx is the identity matrix of size nx and A22**H is the
!>  conjugate transpose of A22. kron(X, Y) is the Kronecker product between
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
!>  based estimate DIF is not wanted (see ZLATDF), then the parameter
!>  IDIFJB (see below) should be changed from 3 to 4 (routine ZLATDF
!>  (IJOB = 2 will be used)). See ZTGSYL for more details.
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
!>  [1] B. Kagstrom; A Direct Method for Reordering Eigenvalues in the
!>      Generalized Real Schur Form of a Regular Matrix Pair (A, B), in
!>      M.S. Moonen et al (eds), Linear Algebra for Large Scale and
!>      Real-Time Applications, Kluwer Academic Publ. 1993, pp 195-218.
!> \n
!>  [2] B. Kagstrom and P. Poromaa; Computing Eigenspaces with Specified
!>      Eigenvalues of a Regular Matrix Pair (A, B) and Condition
!>      Estimation: Theory, Algorithms and Software, Report
!>      UMINF - 94.04, Department of Computing Science, Umea University,
!>      S-901 87 Umea, Sweden, 1994. Also as LAPACK Working Note 87.
!>      To appear in Numerical Algorithms, 1996.
!> \n
!>  [3] B. Kagstrom and P. Poromaa, LAPACK-Style Algorithms and Software
!>      for Solving the Generalized Sylvester Equation and Estimating the
!>      Separation between Regular Matrix Pairs, Report UMINF - 93.23,
!>      Department of Computing Science, Umea University, S-901 87 Umea,
!>      Sweden, December 1993, Revised April 1994, Also as LAPACK working
!>      Note 75. To appear in ACM Trans. on Math. Software, Vol 22, No 1,
!>      1996.
!>
!  =====================================================================
      SUBROUTINE ZTGSEN(Ijob,Wantq,Wantz,Select,N,A,Lda,B,Ldb,Alpha,    &
     &                  Beta,Q,Ldq,Z,Ldz,M,Pl,Pr,Dif,Work,Lwork,Iwork,  &
     &                  Liwork,Info)
      USE F77KINDS                        
      USE S_DLAMCH
      USE S_XERBLA
      USE S_ZLACN2
      USE S_ZLACPY
      USE S_ZLASSQ
      USE S_ZSCAL
      USE S_ZTGEXC
      USE S_ZTGSYL
      IMPLICIT NONE
!*--ZTGSEN446
!
! PARAMETER definitions rewritten by SPAG
!
      INTEGER , PARAMETER  ::  IDIFJB = 3
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER , INTENT(IN) :: Ijob
      LOGICAL :: Wantq
      LOGICAL :: Wantz
      LOGICAL , INTENT(IN) , DIMENSION(*) :: Select
      INTEGER :: N
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX(CX16KIND) , INTENT(OUT) , DIMENSION(*) :: Alpha
      COMPLEX(CX16KIND) , INTENT(OUT) , DIMENSION(*) :: Beta
      COMPLEX(CX16KIND) , DIMENSION(Ldq,*) :: Q
      INTEGER :: Ldq
      COMPLEX(CX16KIND) , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      INTEGER , INTENT(INOUT) :: M
      REAL(R8KIND) , INTENT(INOUT) :: Pl
      REAL(R8KIND) , INTENT(INOUT) :: Pr
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Dif
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(IN) :: Liwork
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      REAL(R8KIND) :: dscale , dsum , rdscal , safmin
      INTEGER :: i , ierr , ijb , k , kase , ks , liwmin , lwmin , mn2 ,&
     &           n1 , n2
      INTEGER , DIMENSION(3) :: isave
      LOGICAL :: lquery , swap , wantd , wantd1 , wantd2 , wantp
      COMPLEX(CX16KIND) :: temp1 , temp2
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
!     .. Intrinsic Functions ..
!     ..
!     .. External Functions ..
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
         Info = -13
      ELSEIF ( Ldz<1 .OR. (Wantz .AND. Ldz<N) ) THEN
         Info = -15
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('ZTGSEN',-Info)
         RETURN
      ENDIF
!
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
      IF ( .NOT.lquery .OR. Ijob/=0 ) THEN
         DO k = 1 , N
            Alpha(k) = A(k,k)
            Beta(k) = B(k,k)
            IF ( k<N ) THEN
               IF ( Select(k) ) M = M + 1
            ELSE
               IF ( Select(N) ) M = M + 1
            ENDIF
         ENDDO
      ENDIF
!
      IF ( Ijob==1 .OR. Ijob==2 .OR. Ijob==4 ) THEN
         lwmin = MAX(1,2*M*(N-M))
         liwmin = MAX(1,N+2)
      ELSEIF ( Ijob==3 .OR. Ijob==5 ) THEN
         lwmin = MAX(1,4*M*(N-M))
         liwmin = MAX(1,2*M*(N-M),N+2)
      ELSE
         lwmin = 1
         liwmin = 1
      ENDIF
!
      Work(1) = lwmin
      Iwork(1) = liwmin
!
      IF ( Lwork<lwmin .AND. .NOT.lquery ) THEN
         Info = -21
      ELSEIF ( Liwork<liwmin .AND. .NOT.lquery ) THEN
         Info = -23
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('ZTGSEN',-Info)
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
               CALL ZLASSQ(N,A(1,i),1,dscale,dsum)
               CALL ZLASSQ(N,B(1,i),1,dscale,dsum)
            ENDDO
            Dif(1) = dscale*SQRT(dsum)
            Dif(2) = Dif(1)
         ENDIF
         GOTO 100
      ENDIF
!
!     Get machine constant
!
      safmin = DLAMCH('S')
!
!     Collect the selected blocks at the top-left corner of (A, B).
!
      ks = 0
      DO k = 1 , N
         swap = Select(k)
         IF ( swap ) THEN
            ks = ks + 1
!
!           Swap the K-th block to position KS. Compute unitary Q
!           and Z that will swap adjacent diagonal blocks in (A, B).
!
            IF ( k/=ks ) CALL ZTGEXC(Wantq,Wantz,N,A,Lda,B,Ldb,Q,Ldq,Z, &
     &                               Ldz,k,ks,ierr)
!
            IF ( ierr>0 ) THEN
!
!              Swap is rejected: exit.
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
         ENDIF
      ENDDO
      IF ( wantp ) THEN
!
!        Solve generalized Sylvester equation for R and L:
!                   A11 * R - L * A22 = A12
!                   B11 * R - L * B22 = B12
!
         n1 = M
         n2 = N - M
         i = n1 + 1
         CALL ZLACPY('Full',n1,n2,A(1,i),Lda,Work,n1)
         CALL ZLACPY('Full',n1,n2,B(1,i),Ldb,Work(n1*n2+1),n1)
         ijb = 0
         CALL ZTGSYL('N',ijb,n1,n2,A,Lda,A(i,i),Lda,Work,n1,B,Ldb,B(i,i)&
     &               ,Ldb,Work(n1*n2+1),n1,dscale,Dif(1),Work(n1*n2*2+1)&
     &               ,Lwork-2*n1*n2,Iwork,ierr)
!
!        Estimate the reciprocal of norms of "projections" onto
!        left and right eigenspaces
!
         rdscal = ZERO
         dsum = ONE
         CALL ZLASSQ(n1*n2,Work,1,rdscal,dsum)
         Pl = rdscal*SQRT(dsum)
         IF ( Pl==ZERO ) THEN
            Pl = ONE
         ELSE
            Pl = dscale/(SQRT(dscale*dscale/Pl+Pl)*SQRT(Pl))
         ENDIF
         rdscal = ZERO
         dsum = ONE
         CALL ZLASSQ(n1*n2,Work(n1*n2+1),1,rdscal,dsum)
         Pr = rdscal*SQRT(dsum)
         IF ( Pr==ZERO ) THEN
            Pr = ONE
         ELSE
            Pr = dscale/(SQRT(dscale*dscale/Pr+Pr)*SQRT(Pr))
         ENDIF
      ENDIF
      IF ( wantd ) THEN
!
!        Compute estimates Difu and Difl.
!
         IF ( wantd1 ) THEN
            n1 = M
            n2 = N - M
            i = n1 + 1
            ijb = IDIFJB
!
!           Frobenius norm-based Difu estimate.
!
            CALL ZTGSYL('N',ijb,n1,n2,A,Lda,A(i,i),Lda,Work,n1,B,Ldb,   &
     &                  B(i,i),Ldb,Work(n1*n2+1),n1,dscale,Dif(1),      &
     &                  Work(n1*n2*2+1),Lwork-2*n1*n2,Iwork,ierr)
!
!           Frobenius norm-based Difl estimate.
!
            CALL ZTGSYL('N',ijb,n2,n1,A(i,i),Lda,A,Lda,Work,n2,B(i,i),  &
     &                  Ldb,B,Ldb,Work(n1*n2+1),n2,dscale,Dif(2),       &
     &                  Work(n1*n2*2+1),Lwork-2*n1*n2,Iwork,ierr)
         ELSE
!
!           Compute 1-norm-based estimates of Difu and Difl using
!           reversed communication with ZLACN2. In each step a
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
               CALL ZLACN2(mn2,Work(mn2+1),Work,Dif(1),kase,isave)
               IF ( kase/=0 ) THEN
                  IF ( kase==1 ) THEN
!
!                 Solve generalized Sylvester equation
!
                     CALL ZTGSYL('N',ijb,n1,n2,A,Lda,A(i,i),Lda,Work,n1,&
     &                           B,Ldb,B(i,i),Ldb,Work(n1*n2+1),n1,     &
     &                           dscale,Dif(1),Work(n1*n2*2+1),         &
     &                           Lwork-2*n1*n2,Iwork,ierr)
                  ELSE
!
!                 Solve the transposed variant.
!
                     CALL ZTGSYL('C',ijb,n1,n2,A,Lda,A(i,i),Lda,Work,n1,&
     &                           B,Ldb,B(i,i),Ldb,Work(n1*n2+1),n1,     &
     &                           dscale,Dif(1),Work(n1*n2*2+1),         &
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
               CALL ZLACN2(mn2,Work(mn2+1),Work,Dif(2),kase,isave)
               IF ( kase/=0 ) THEN
                  IF ( kase==1 ) THEN
!
!                 Solve generalized Sylvester equation
!
                     CALL ZTGSYL('N',ijb,n2,n1,A(i,i),Lda,A,Lda,Work,n2,&
     &                           B(i,i),Ldb,B,Ldb,Work(n1*n2+1),n2,     &
     &                           dscale,Dif(2),Work(n1*n2*2+1),         &
     &                           Lwork-2*n1*n2,Iwork,ierr)
                  ELSE
!
!                 Solve the transposed variant.
!
                     CALL ZTGSYL('C',ijb,n2,n1,A(i,i),Lda,A,Lda,Work,n2,&
     &                           B,Ldb,B(i,i),Ldb,Work(n1*n2+1),n2,     &
     &                           dscale,Dif(2),Work(n1*n2*2+1),         &
     &                           Lwork-2*n1*n2,Iwork,ierr)
                  ENDIF
                  CYCLE
               ENDIF
               Dif(2) = dscale/Dif(2)
               EXIT
            ENDDO
         ENDIF
      ENDIF
!
!     If B(K,K) is complex, make it real and positive (normalization
!     of the generalized Schur form) and Store the generalized
!     eigenvalues of reordered pair (A, B)
!
      DO k = 1 , N
         dscale = ABS(B(k,k))
         IF ( dscale>safmin ) THEN
            temp1 = DCONJG(B(k,k)/dscale)
            temp2 = B(k,k)/dscale
            B(k,k) = dscale
            CALL ZSCAL(N-k,temp1,B(k,k+1),Ldb)
            CALL ZSCAL(N-k+1,temp1,A(k,k),Lda)
            IF ( Wantq ) CALL ZSCAL(N,temp2,Q(1,k),1)
         ELSE
            B(k,k) = DCMPLX(ZERO,ZERO)
         ENDIF
!
         Alpha(k) = A(k,k)
         Beta(k) = B(k,k)
!
      ENDDO
!
!
 100  Work(1) = lwmin
      Iwork(1) = liwmin
!
!
!     End of ZTGSEN
!
      END SUBROUTINE ZTGSEN
