!*==ctgsen.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b CTGSEN
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CTGSEN + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ctgsen.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ctgsen.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ctgsen.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CTGSEN( IJOB, WANTQ, WANTZ, SELECT, N, A, LDA, B, LDB,
!                          ALPHA, BETA, Q, LDQ, Z, LDZ, M, PL, PR, DIF,
!                          WORK, LWORK, IWORK, LIWORK, INFO )
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
!       REAL               DIF( * )
!       COMPLEX            A( LDA, * ), ALPHA( * ), B( LDB, * ),
!      $                   BETA( * ), Q( LDQ, * ), WORK( * ), Z( LDZ, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CTGSEN reorders the generalized Schur decomposition of a complex
!> matrix pair (A, B) (in terms of an unitary equivalence trans-
!> formation Q**H * (A, B) * Z), so that a selected cluster of eigenvalues
!> appears in the leading diagonal blocks of the pair (A,B). The leading
!> columns of Q and Z form unitary bases of the corresponding left and
!> right eigenspaces (deflating subspaces). (A, B) must be in
!> generalized Schur canonical form, that is, A and B are both upper
!> triangular.
!>
!> CTGSEN also computes the generalized eigenvalues
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
!>          A is COMPLEX array, dimension(LDA,N)
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
!>          B is COMPLEX array, dimension(LDB,N)
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
!>          ALPHA is COMPLEX array, dimension (N)
!> \endverbatim
!>
!> \param[out] BETA
!> \verbatim
!>          BETA is COMPLEX array, dimension (N)
!>
!>          The diagonal elements of A and B, respectively,
!>          when the pair (A,B) has been reduced to generalized Schur
!>          form.  ALPHA(i)/BETA(i) i=1,...,N are the generalized
!>          eigenvalues.
!> \endverbatim
!>
!> \param[in,out] Q
!> \verbatim
!>          Q is COMPLEX array, dimension (LDQ,N)
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
!>          Z is COMPLEX array, dimension (LDZ,N)
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
!>          PL is REAL
!> \endverbatim
!>
!> \param[out] PR
!> \verbatim
!>          PR is REAL
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
!>          DIF is REAL array, dimension (2).
!>          If IJOB >= 2, DIF(1:2) store the estimates of Difu and Difl.
!>          If IJOB = 2 or 4, DIF(1:2) are F-norm-based upper bounds on
!>          Difu and Difl. If IJOB = 3 or 5, DIF(1:2) are 1-norm-based
!>          estimates of Difu and Difl, computed using reversed
!>          communication with CLACN2.
!>          If M = 0 or N, DIF(1:2) = F-norm([A, B]).
!>          If IJOB = 0 or 1, DIF is not referenced.
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
!> \ingroup complexOTHERcomputational
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  CTGSEN first collects the selected eigenvalues by computing unitary
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
!>  decomposition of a matrix pair (C, D) = Q*(A, B)*Z', then the
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
!>  conjuguate transpose of A22. kron(X, Y) is the Kronecker product between
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
!>  based estimate DIF is not wanted (see CLATDF), then the parameter
!>  IDIFJB (see below) should be changed from 3 to 4 (routine CLATDF
!>  (IJOB = 2 will be used)). See CTGSYL for more details.
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
      SUBROUTINE CTGSEN(Ijob,Wantq,Wantz,Select,N,A,Lda,B,Ldb,Alpha,    &
     &                  Beta,Q,Ldq,Z,Ldz,M,Pl,Pr,Dif,Work,Lwork,Iwork,  &
     &                  Liwork,Info)
      IMPLICIT NONE
!*--CTGSEN437
!
!  -- LAPACK computational routine (version 3.7.1) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     June 2016
!
!     .. Scalar Arguments ..
      LOGICAL Wantq , Wantz
      INTEGER Ijob , Info , Lda , Ldb , Ldq , Ldz , Liwork , Lwork , M ,&
     &        N
      REAL Pl , Pr
!     ..
!     .. Array Arguments ..
      LOGICAL Select(*)
      INTEGER Iwork(*)
      REAL Dif(*)
      COMPLEX A(Lda,*) , Alpha(*) , B(Ldb,*) , Beta(*) , Q(Ldq,*) ,     &
     &        Work(*) , Z(Ldz,*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      INTEGER IDIFJB
      PARAMETER (IDIFJB=3)
      REAL ZERO , ONE
      PARAMETER (ZERO=0.0E+0,ONE=1.0E+0)
!     ..
!     .. Local Scalars ..
      LOGICAL lquery , swap , wantd , wantd1 , wantd2 , wantp
      INTEGER i , ierr , ijb , k , kase , ks , liwmin , lwmin , mn2 ,   &
     &        n1 , n2
      REAL dscale , dsum , rdscal , safmin
      COMPLEX temp1 , temp2
!     ..
!     .. Local Arrays ..
      INTEGER isave(3)
!     ..
!     .. External Subroutines ..
      REAL SLAMCH
      EXTERNAL CLACN2 , CLACPY , CLASSQ , CSCAL , CTGEXC , CTGSYL ,     &
     &         SLAMCH , XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , CMPLX , CONJG , MAX , SQRT
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
         CALL XERBLA('CTGSEN',-Info)
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
         CALL XERBLA('CTGSEN',-Info)
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
               CALL CLASSQ(N,A(1,i),1,dscale,dsum)
               CALL CLASSQ(N,B(1,i),1,dscale,dsum)
            ENDDO
            Dif(1) = dscale*SQRT(dsum)
            Dif(2) = Dif(1)
         ENDIF
         GOTO 100
      ENDIF
!
!     Get machine constant
!
      safmin = SLAMCH('S')
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
            IF ( k/=ks ) CALL CTGEXC(Wantq,Wantz,N,A,Lda,B,Ldb,Q,Ldq,Z, &
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
         CALL CLACPY('Full',n1,n2,A(1,i),Lda,Work,n1)
         CALL CLACPY('Full',n1,n2,B(1,i),Ldb,Work(n1*n2+1),n1)
         ijb = 0
         CALL CTGSYL('N',ijb,n1,n2,A,Lda,A(i,i),Lda,Work,n1,B,Ldb,B(i,i)&
     &               ,Ldb,Work(n1*n2+1),n1,dscale,Dif(1),Work(n1*n2*2+1)&
     &               ,Lwork-2*n1*n2,Iwork,ierr)
!
!        Estimate the reciprocal of norms of "projections" onto
!        left and right eigenspaces
!
         rdscal = ZERO
         dsum = ONE
         CALL CLASSQ(n1*n2,Work,1,rdscal,dsum)
         Pl = rdscal*SQRT(dsum)
         IF ( Pl==ZERO ) THEN
            Pl = ONE
         ELSE
            Pl = dscale/(SQRT(dscale*dscale/Pl+Pl)*SQRT(Pl))
         ENDIF
         rdscal = ZERO
         dsum = ONE
         CALL CLASSQ(n1*n2,Work(n1*n2+1),1,rdscal,dsum)
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
            CALL CTGSYL('N',ijb,n1,n2,A,Lda,A(i,i),Lda,Work,n1,B,Ldb,   &
     &                  B(i,i),Ldb,Work(n1*n2+1),n1,dscale,Dif(1),      &
     &                  Work(n1*n2*2+1),Lwork-2*n1*n2,Iwork,ierr)
!
!           Frobenius norm-based Difl estimate.
!
            CALL CTGSYL('N',ijb,n2,n1,A(i,i),Lda,A,Lda,Work,n2,B(i,i),  &
     &                  Ldb,B,Ldb,Work(n1*n2+1),n2,dscale,Dif(2),       &
     &                  Work(n1*n2*2+1),Lwork-2*n1*n2,Iwork,ierr)
         ELSE
!
!           Compute 1-norm-based estimates of Difu and Difl using
!           reversed communication with CLACN2. In each step a
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
               CALL CLACN2(mn2,Work(mn2+1),Work,Dif(1),kase,isave)
               IF ( kase/=0 ) THEN
                  IF ( kase==1 ) THEN
!
!                 Solve generalized Sylvester equation
!
                     CALL CTGSYL('N',ijb,n1,n2,A,Lda,A(i,i),Lda,Work,n1,&
     &                           B,Ldb,B(i,i),Ldb,Work(n1*n2+1),n1,     &
     &                           dscale,Dif(1),Work(n1*n2*2+1),         &
     &                           Lwork-2*n1*n2,Iwork,ierr)
                  ELSE
!
!                 Solve the transposed variant.
!
                     CALL CTGSYL('C',ijb,n1,n2,A,Lda,A(i,i),Lda,Work,n1,&
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
               CALL CLACN2(mn2,Work(mn2+1),Work,Dif(2),kase,isave)
               IF ( kase/=0 ) THEN
                  IF ( kase==1 ) THEN
!
!                 Solve generalized Sylvester equation
!
                     CALL CTGSYL('N',ijb,n2,n1,A(i,i),Lda,A,Lda,Work,n2,&
     &                           B(i,i),Ldb,B,Ldb,Work(n1*n2+1),n2,     &
     &                           dscale,Dif(2),Work(n1*n2*2+1),         &
     &                           Lwork-2*n1*n2,Iwork,ierr)
                  ELSE
!
!                 Solve the transposed variant.
!
                     CALL CTGSYL('C',ijb,n2,n1,A(i,i),Lda,A,Lda,Work,n2,&
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
            temp1 = CONJG(B(k,k)/dscale)
            temp2 = B(k,k)/dscale
            B(k,k) = dscale
            CALL CSCAL(N-k,temp1,B(k,k+1),Ldb)
            CALL CSCAL(N-k+1,temp1,A(k,k),Lda)
            IF ( Wantq ) CALL CSCAL(N,temp2,Q(1,k),1)
         ELSE
            B(k,k) = CMPLX(ZERO,ZERO)
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
!     End of CTGSEN
!
      END SUBROUTINE CTGSEN
