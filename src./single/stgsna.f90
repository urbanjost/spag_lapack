!*==stgsna.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b STGSNA
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download STGSNA + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/stgsna.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/stgsna.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/stgsna.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE STGSNA( JOB, HOWMNY, SELECT, N, A, LDA, B, LDB, VL,
!                          LDVL, VR, LDVR, S, DIF, MM, M, WORK, LWORK,
!                          IWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          HOWMNY, JOB
!       INTEGER            INFO, LDA, LDB, LDVL, LDVR, LWORK, M, MM, N
!       ..
!       .. Array Arguments ..
!       LOGICAL            SELECT( * )
!       INTEGER            IWORK( * )
!       REAL               A( LDA, * ), B( LDB, * ), DIF( * ), S( * ),
!      $                   VL( LDVL, * ), VR( LDVR, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> STGSNA estimates reciprocal condition numbers for specified
!> eigenvalues and/or eigenvectors of a matrix pair (A, B) in
!> generalized real Schur canonical form (or of any matrix pair
!> (Q*A*Z**T, Q*B*Z**T) with orthogonal matrices Q and Z, where
!> Z**T denotes the transpose of Z.
!>
!> (A, B) must be in generalized real Schur form (as returned by SGGES),
!> i.e. A is block upper triangular with 1-by-1 and 2-by-2 diagonal
!> blocks. B is upper triangular.
!>
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] JOB
!> \verbatim
!>          JOB is CHARACTER*1
!>          Specifies whether condition numbers are required for
!>          eigenvalues (S) or eigenvectors (DIF):
!>          = 'E': for eigenvalues only (S);
!>          = 'V': for eigenvectors only (DIF);
!>          = 'B': for both eigenvalues and eigenvectors (S and DIF).
!> \endverbatim
!>
!> \param[in] HOWMNY
!> \verbatim
!>          HOWMNY is CHARACTER*1
!>          = 'A': compute condition numbers for all eigenpairs;
!>          = 'S': compute condition numbers for selected eigenpairs
!>                 specified by the array SELECT.
!> \endverbatim
!>
!> \param[in] SELECT
!> \verbatim
!>          SELECT is LOGICAL array, dimension (N)
!>          If HOWMNY = 'S', SELECT specifies the eigenpairs for which
!>          condition numbers are required. To select condition numbers
!>          for the eigenpair corresponding to a real eigenvalue w(j),
!>          SELECT(j) must be set to .TRUE.. To select condition numbers
!>          corresponding to a complex conjugate pair of eigenvalues w(j)
!>          and w(j+1), either SELECT(j) or SELECT(j+1) or both, must be
!>          set to .TRUE..
!>          If HOWMNY = 'A', SELECT is not referenced.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the square matrix pair (A, B). N >= 0.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is REAL array, dimension (LDA,N)
!>          The upper quasi-triangular matrix A in the pair (A,B).
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A. LDA >= max(1,N).
!> \endverbatim
!>
!> \param[in] B
!> \verbatim
!>          B is REAL array, dimension (LDB,N)
!>          The upper triangular matrix B in the pair (A,B).
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of the array B. LDB >= max(1,N).
!> \endverbatim
!>
!> \param[in] VL
!> \verbatim
!>          VL is REAL array, dimension (LDVL,M)
!>          If JOB = 'E' or 'B', VL must contain left eigenvectors of
!>          (A, B), corresponding to the eigenpairs specified by HOWMNY
!>          and SELECT. The eigenvectors must be stored in consecutive
!>          columns of VL, as returned by STGEVC.
!>          If JOB = 'V', VL is not referenced.
!> \endverbatim
!>
!> \param[in] LDVL
!> \verbatim
!>          LDVL is INTEGER
!>          The leading dimension of the array VL. LDVL >= 1.
!>          If JOB = 'E' or 'B', LDVL >= N.
!> \endverbatim
!>
!> \param[in] VR
!> \verbatim
!>          VR is REAL array, dimension (LDVR,M)
!>          If JOB = 'E' or 'B', VR must contain right eigenvectors of
!>          (A, B), corresponding to the eigenpairs specified by HOWMNY
!>          and SELECT. The eigenvectors must be stored in consecutive
!>          columns ov VR, as returned by STGEVC.
!>          If JOB = 'V', VR is not referenced.
!> \endverbatim
!>
!> \param[in] LDVR
!> \verbatim
!>          LDVR is INTEGER
!>          The leading dimension of the array VR. LDVR >= 1.
!>          If JOB = 'E' or 'B', LDVR >= N.
!> \endverbatim
!>
!> \param[out] S
!> \verbatim
!>          S is REAL array, dimension (MM)
!>          If JOB = 'E' or 'B', the reciprocal condition numbers of the
!>          selected eigenvalues, stored in consecutive elements of the
!>          array. For a complex conjugate pair of eigenvalues two
!>          consecutive elements of S are set to the same value. Thus
!>          S(j), DIF(j), and the j-th columns of VL and VR all
!>          correspond to the same eigenpair (but not in general the
!>          j-th eigenpair, unless all eigenpairs are selected).
!>          If JOB = 'V', S is not referenced.
!> \endverbatim
!>
!> \param[out] DIF
!> \verbatim
!>          DIF is REAL array, dimension (MM)
!>          If JOB = 'V' or 'B', the estimated reciprocal condition
!>          numbers of the selected eigenvectors, stored in consecutive
!>          elements of the array. For a complex eigenvector two
!>          consecutive elements of DIF are set to the same value. If
!>          the eigenvalues cannot be reordered to compute DIF(j), DIF(j)
!>          is set to 0; this can only occur when the true value would be
!>          very small anyway.
!>          If JOB = 'E', DIF is not referenced.
!> \endverbatim
!>
!> \param[in] MM
!> \verbatim
!>          MM is INTEGER
!>          The number of elements in the arrays S and DIF. MM >= M.
!> \endverbatim
!>
!> \param[out] M
!> \verbatim
!>          M is INTEGER
!>          The number of elements of the arrays S and DIF used to store
!>          the specified condition numbers; for each selected real
!>          eigenvalue one element is used, and for each selected complex
!>          conjugate pair of eigenvalues, two elements are used.
!>          If HOWMNY = 'A', M is set to N.
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
!>          The dimension of the array WORK. LWORK >= max(1,N).
!>          If JOB = 'V' or 'B' LWORK >= 2*N*(N+2)+16.
!>
!>          If LWORK = -1, then a workspace query is assumed; the routine
!>          only calculates the optimal size of the WORK array, returns
!>          this value as the first entry of the WORK array, and no error
!>          message related to LWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] IWORK
!> \verbatim
!>          IWORK is INTEGER array, dimension (N + 6)
!>          If JOB = 'E', IWORK is not referenced.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          =0: Successful exit
!>          <0: If INFO = -i, the i-th argument had an illegal value
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
!> \ingroup realOTHERcomputational
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  The reciprocal of the condition number of a generalized eigenvalue
!>  w = (a, b) is defined as
!>
!>       S(w) = (|u**TAv|**2 + |u**TBv|**2)**(1/2) / (norm(u)*norm(v))
!>
!>  where u and v are the left and right eigenvectors of (A, B)
!>  corresponding to w; |z| denotes the absolute value of the complex
!>  number, and norm(u) denotes the 2-norm of the vector u.
!>  The pair (a, b) corresponds to an eigenvalue w = a/b (= u**TAv/u**TBv)
!>  of the matrix pair (A, B). If both a and b equal zero, then (A B) is
!>  singular and S(I) = -1 is returned.
!>
!>  An approximate error bound on the chordal distance between the i-th
!>  computed generalized eigenvalue w and the corresponding exact
!>  eigenvalue lambda is
!>
!>       chord(w, lambda) <= EPS * norm(A, B) / S(I)
!>
!>  where EPS is the machine precision.
!>
!>  The reciprocal of the condition number DIF(i) of right eigenvector u
!>  and left eigenvector v corresponding to the generalized eigenvalue w
!>  is defined as follows:
!>
!>  a) If the i-th eigenvalue w = (a,b) is real
!>
!>     Suppose U and V are orthogonal transformations such that
!>
!>              U**T*(A, B)*V  = (S, T) = ( a   *  ) ( b  *  )  1
!>                                        ( 0  S22 ),( 0 T22 )  n-1
!>                                          1  n-1     1 n-1
!>
!>     Then the reciprocal condition number DIF(i) is
!>
!>                Difl((a, b), (S22, T22)) = sigma-min( Zl ),
!>
!>     where sigma-min(Zl) denotes the smallest singular value of the
!>     2(n-1)-by-2(n-1) matrix
!>
!>         Zl = [ kron(a, In-1)  -kron(1, S22) ]
!>              [ kron(b, In-1)  -kron(1, T22) ] .
!>
!>     Here In-1 is the identity matrix of size n-1. kron(X, Y) is the
!>     Kronecker product between the matrices X and Y.
!>
!>     Note that if the default method for computing DIF(i) is wanted
!>     (see SLATDF), then the parameter DIFDRI (see below) should be
!>     changed from 3 to 4 (routine SLATDF(IJOB = 2 will be used)).
!>     See STGSYL for more details.
!>
!>  b) If the i-th and (i+1)-th eigenvalues are complex conjugate pair,
!>
!>     Suppose U and V are orthogonal transformations such that
!>
!>              U**T*(A, B)*V = (S, T) = ( S11  *   ) ( T11  *  )  2
!>                                       ( 0    S22 ),( 0    T22) n-2
!>                                         2    n-2     2    n-2
!>
!>     and (S11, T11) corresponds to the complex conjugate eigenvalue
!>     pair (w, conjg(w)). There exist unitary matrices U1 and V1 such
!>     that
!>
!>       U1**T*S11*V1 = ( s11 s12 ) and U1**T*T11*V1 = ( t11 t12 )
!>                      (  0  s22 )                    (  0  t22 )
!>
!>     where the generalized eigenvalues w = s11/t11 and
!>     conjg(w) = s22/t22.
!>
!>     Then the reciprocal condition number DIF(i) is bounded by
!>
!>         min( d1, max( 1, |real(s11)/real(s22)| )*d2 )
!>
!>     where, d1 = Difl((s11, t11), (s22, t22)) = sigma-min(Z1), where
!>     Z1 is the complex 2-by-2 matrix
!>
!>              Z1 =  [ s11  -s22 ]
!>                    [ t11  -t22 ],
!>
!>     This is done by computing (using real arithmetic) the
!>     roots of the characteristical polynomial det(Z1**T * Z1 - lambda I),
!>     where Z1**T denotes the transpose of Z1 and det(X) denotes
!>     the determinant of X.
!>
!>     and d2 is an upper bound on Difl((S11, T11), (S22, T22)), i.e. an
!>     upper bound on sigma-min(Z2), where Z2 is (2n-2)-by-(2n-2)
!>
!>              Z2 = [ kron(S11**T, In-2)  -kron(I2, S22) ]
!>                   [ kron(T11**T, In-2)  -kron(I2, T22) ]
!>
!>     Note that if the default method for computing DIF is wanted (see
!>     SLATDF), then the parameter DIFDRI (see below) should be changed
!>     from 3 to 4 (routine SLATDF(IJOB = 2 will be used)). See STGSYL
!>     for more details.
!>
!>  For each eigenvalue/vector specified by SELECT, DIF stores a
!>  Frobenius norm-based estimate of Difl.
!>
!>  An approximate error bound for the i-th computed eigenvector VL(i) or
!>  VR(i) is given by
!>
!>             EPS * norm(A, B) / DIF(i).
!>
!>  See ref. [2-3] for more details and further references.
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
!>      Note 75.  To appear in ACM Trans. on Math. Software, Vol 22,
!>      No 1, 1996.
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE STGSNA(Job,Howmny,Select,N,A,Lda,B,Ldb,Vl,Ldvl,Vr,Ldvr,&
     &                  S,Dif,Mm,M,Work,Lwork,Iwork,Info)
      USE S_LSAME
      USE S_SDOT
      USE S_SGEMV
      USE S_SLACPY
      USE S_SLAG2
      USE S_SLAMCH
      USE S_SLAPY2
      USE S_SNRM2
      USE S_STGEXC
      USE S_STGSYL
      USE S_XERBLA
      IMPLICIT NONE
!*--STGSNA395
!
! PARAMETER definitions rewritten by SPAG
!
      INTEGER , PARAMETER  ::  DIFDRI = 3
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0 ,              &
     &                      TWO = 2.0E+0 , FOUR = 4.0E+0
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Job
      CHARACTER :: Howmny
      LOGICAL , INTENT(IN) , DIMENSION(*) :: Select
      INTEGER :: N
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL , DIMENSION(Ldvl,*) :: Vl
      INTEGER , INTENT(IN) :: Ldvl
      REAL , DIMENSION(Ldvr,*) :: Vr
      INTEGER , INTENT(IN) :: Ldvr
      REAL , INTENT(INOUT) , DIMENSION(*) :: S
      REAL , INTENT(INOUT) , DIMENSION(*) :: Dif
      INTEGER , INTENT(IN) :: Mm
      INTEGER , INTENT(INOUT) :: M
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      REAL :: alphai , alphar , alprqt , beta , c1 , c2 , cond , eps ,  &
     &        lnrm , rnrm , root1 , root2 , scale , smlnum , tmpii ,    &
     &        tmpir , tmpri , tmprr , uhav , uhavi , uhbv , uhbvi
      REAL , DIMENSION(1) :: dummy , dummy1
      INTEGER :: i , ierr , ifst , ilst , iz , k , ks , lwmin , n1 , n2
      LOGICAL :: lquery , pair , somcon , wantbh , wantdf , wants
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
!     .. External Functions ..
!     ..
!     .. External Subroutines ..
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. Executable Statements ..
!
!     Decode and test the input parameters
!
      wantbh = LSAME(Job,'B')
      wants = LSAME(Job,'E') .OR. wantbh
      wantdf = LSAME(Job,'V') .OR. wantbh
!
      somcon = LSAME(Howmny,'S')
!
      Info = 0
      lquery = (Lwork==-1)
!
      IF ( .NOT.wants .AND. .NOT.wantdf ) THEN
         Info = -1
      ELSEIF ( .NOT.LSAME(Howmny,'A') .AND. .NOT.somcon ) THEN
         Info = -2
      ELSEIF ( N<0 ) THEN
         Info = -4
      ELSEIF ( Lda<MAX(1,N) ) THEN
         Info = -6
      ELSEIF ( Ldb<MAX(1,N) ) THEN
         Info = -8
      ELSEIF ( wants .AND. Ldvl<N ) THEN
         Info = -10
      ELSEIF ( wants .AND. Ldvr<N ) THEN
         Info = -12
      ELSE
!
!        Set M to the number of eigenpairs for which condition numbers
!        are required, and test MM.
!
         IF ( somcon ) THEN
            M = 0
            pair = .FALSE.
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
         ELSE
            M = N
         ENDIF
!
         IF ( N==0 ) THEN
            lwmin = 1
         ELSEIF ( LSAME(Job,'V') .OR. LSAME(Job,'B') ) THEN
            lwmin = 2*N*(N+2) + 16
         ELSE
            lwmin = N
         ENDIF
         Work(1) = lwmin
!
         IF ( Mm<M ) THEN
            Info = -15
         ELSEIF ( Lwork<lwmin .AND. .NOT.lquery ) THEN
            Info = -18
         ENDIF
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('STGSNA',-Info)
         RETURN
      ELSEIF ( lquery ) THEN
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( N==0 ) RETURN
!
!     Get machine constants
!
      eps = SLAMCH('P')
      smlnum = SLAMCH('S')/eps
      ks = 0
      pair = .FALSE.
!
      DO k = 1 , N
!
!        Determine whether A(k,k) begins a 1-by-1 or 2-by-2 block.
!
         IF ( pair ) THEN
            pair = .FALSE.
            CYCLE
         ELSE
            IF ( k<N ) pair = A(k+1,k)/=ZERO
         ENDIF
!
!        Determine whether condition numbers are required for the k-th
!        eigenpair.
!
         IF ( somcon ) THEN
            IF ( pair ) THEN
               IF ( .NOT.Select(k) .AND. .NOT.Select(k+1) ) CYCLE
            ELSEIF ( .NOT.Select(k) ) THEN
               CYCLE
            ENDIF
         ENDIF
!
         ks = ks + 1
!
         IF ( wants ) THEN
!
!           Compute the reciprocal condition number of the k-th
!           eigenvalue.
!
            IF ( pair ) THEN
!
!              Complex eigenvalue pair.
!
               rnrm = SLAPY2(SNRM2(N,Vr(1,ks),1),SNRM2(N,Vr(1,ks+1),1))
               lnrm = SLAPY2(SNRM2(N,Vl(1,ks),1),SNRM2(N,Vl(1,ks+1),1))
               CALL SGEMV('N',N,N,ONE,A,Lda,Vr(1,ks),1,ZERO,Work,1)
               tmprr = SDOT(N,Work,1,Vl(1,ks),1)
               tmpri = SDOT(N,Work,1,Vl(1,ks+1),1)
               CALL SGEMV('N',N,N,ONE,A,Lda,Vr(1,ks+1),1,ZERO,Work,1)
               tmpii = SDOT(N,Work,1,Vl(1,ks+1),1)
               tmpir = SDOT(N,Work,1,Vl(1,ks),1)
               uhav = tmprr + tmpii
               uhavi = tmpir - tmpri
               CALL SGEMV('N',N,N,ONE,B,Ldb,Vr(1,ks),1,ZERO,Work,1)
               tmprr = SDOT(N,Work,1,Vl(1,ks),1)
               tmpri = SDOT(N,Work,1,Vl(1,ks+1),1)
               CALL SGEMV('N',N,N,ONE,B,Ldb,Vr(1,ks+1),1,ZERO,Work,1)
               tmpii = SDOT(N,Work,1,Vl(1,ks+1),1)
               tmpir = SDOT(N,Work,1,Vl(1,ks),1)
               uhbv = tmprr + tmpii
               uhbvi = tmpir - tmpri
               uhav = SLAPY2(uhav,uhavi)
               uhbv = SLAPY2(uhbv,uhbvi)
               cond = SLAPY2(uhav,uhbv)
               S(ks) = cond/(rnrm*lnrm)
               S(ks+1) = S(ks)
!
            ELSE
!
!              Real eigenvalue.
!
               rnrm = SNRM2(N,Vr(1,ks),1)
               lnrm = SNRM2(N,Vl(1,ks),1)
               CALL SGEMV('N',N,N,ONE,A,Lda,Vr(1,ks),1,ZERO,Work,1)
               uhav = SDOT(N,Work,1,Vl(1,ks),1)
               CALL SGEMV('N',N,N,ONE,B,Ldb,Vr(1,ks),1,ZERO,Work,1)
               uhbv = SDOT(N,Work,1,Vl(1,ks),1)
               cond = SLAPY2(uhav,uhbv)
               IF ( cond==ZERO ) THEN
                  S(ks) = -ONE
               ELSE
                  S(ks) = cond/(rnrm*lnrm)
               ENDIF
            ENDIF
         ENDIF
!
         IF ( wantdf ) THEN
            IF ( N==1 ) THEN
               Dif(ks) = SLAPY2(A(1,1),B(1,1))
               CYCLE
            ENDIF
!
!           Estimate the reciprocal condition number of the k-th
!           eigenvectors.
            IF ( pair ) THEN
!
!              Copy the  2-by 2 pencil beginning at (A(k,k), B(k, k)).
!              Compute the eigenvalue(s) at position K.
!
               Work(1) = A(k,k)
               Work(2) = A(k+1,k)
               Work(3) = A(k,k+1)
               Work(4) = A(k+1,k+1)
               Work(5) = B(k,k)
               Work(6) = B(k+1,k)
               Work(7) = B(k,k+1)
               Work(8) = B(k+1,k+1)
               CALL SLAG2(Work,2,Work(5),2,smlnum*eps,beta,dummy1(1),   &
     &                    alphar,dummy(1),alphai)
               alprqt = ONE
               c1 = TWO*(alphar*alphar+alphai*alphai+beta*beta)
               c2 = FOUR*beta*beta*alphai*alphai
               root1 = c1 + SQRT(c1*c1-4.0*c2)
               root2 = c2/root1
               root1 = root1/TWO
               cond = MIN(SQRT(root1),SQRT(root2))
            ENDIF
!
!           Copy the matrix (A, B) to the array WORK and swap the
!           diagonal block beginning at A(k,k) to the (1,1) position.
!
            CALL SLACPY('Full',N,N,A,Lda,Work,N)
            CALL SLACPY('Full',N,N,B,Ldb,Work(N*N+1),N)
            ifst = k
            ilst = 1
!
            CALL STGEXC(.FALSE.,.FALSE.,N,Work,N,Work(N*N+1),N,dummy,1, &
     &                  dummy1,1,ifst,ilst,Work(N*N*2+1),Lwork-2*N*N,   &
     &                  ierr)
!
            IF ( ierr>0 ) THEN
!
!              Ill-conditioned problem - swap rejected.
!
               Dif(ks) = ZERO
            ELSE
!
!              Reordering successful, solve generalized Sylvester
!              equation for R and L,
!                         A22 * R - L * A11 = A12
!                         B22 * R - L * B11 = B12,
!              and compute estimate of Difl((A11,B11), (A22, B22)).
!
               n1 = 1
               IF ( Work(2)/=ZERO ) n1 = 2
               n2 = N - n1
               IF ( n2==0 ) THEN
                  Dif(ks) = cond
               ELSE
                  i = N*N + 1
                  iz = 2*N*N + 1
                  CALL STGSYL('N',DIFDRI,n2,n1,Work(N*n1+n1+1),N,Work,N,&
     &                        Work(n1+1),N,Work(N*n1+n1+i),N,Work(i),N, &
     &                        Work(n1+i),N,scale,Dif(ks),Work(iz+1),    &
     &                        Lwork-2*N*N,Iwork,ierr)
!
                  IF ( pair ) Dif(ks) = MIN(MAX(ONE,alprqt)*Dif(ks),cond&
     &                                  )
               ENDIF
            ENDIF
            IF ( pair ) Dif(ks+1) = Dif(ks)
         ENDIF
         IF ( pair ) ks = ks + 1
!
      ENDDO
      Work(1) = lwmin
!
!     End of STGSNA
!
      END SUBROUTINE STGSNA