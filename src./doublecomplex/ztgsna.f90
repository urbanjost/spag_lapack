!*==ztgsna.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b ZTGSNA
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZTGSNA + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ztgsna.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ztgsna.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ztgsna.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZTGSNA( JOB, HOWMNY, SELECT, N, A, LDA, B, LDB, VL,
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
!       DOUBLE PRECISION   DIF( * ), S( * )
!       COMPLEX*16         A( LDA, * ), B( LDB, * ), VL( LDVL, * ),
!      $                   VR( LDVR, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZTGSNA estimates reciprocal condition numbers for specified
!> eigenvalues and/or eigenvectors of a matrix pair (A, B).
!>
!> (A, B) must be in generalized Schur canonical form, that is, A and
!> B are both upper triangular.
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
!>          for the corresponding j-th eigenvalue and/or eigenvector,
!>          SELECT(j) must be set to .TRUE..
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
!>          A is COMPLEX*16 array, dimension (LDA,N)
!>          The upper triangular matrix A in the pair (A,B).
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
!>          B is COMPLEX*16 array, dimension (LDB,N)
!>          The upper triangular matrix B in the pair (A, B).
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
!>          VL is COMPLEX*16 array, dimension (LDVL,M)
!>          IF JOB = 'E' or 'B', VL must contain left eigenvectors of
!>          (A, B), corresponding to the eigenpairs specified by HOWMNY
!>          and SELECT.  The eigenvectors must be stored in consecutive
!>          columns of VL, as returned by ZTGEVC.
!>          If JOB = 'V', VL is not referenced.
!> \endverbatim
!>
!> \param[in] LDVL
!> \verbatim
!>          LDVL is INTEGER
!>          The leading dimension of the array VL. LDVL >= 1; and
!>          If JOB = 'E' or 'B', LDVL >= N.
!> \endverbatim
!>
!> \param[in] VR
!> \verbatim
!>          VR is COMPLEX*16 array, dimension (LDVR,M)
!>          IF JOB = 'E' or 'B', VR must contain right eigenvectors of
!>          (A, B), corresponding to the eigenpairs specified by HOWMNY
!>          and SELECT.  The eigenvectors must be stored in consecutive
!>          columns of VR, as returned by ZTGEVC.
!>          If JOB = 'V', VR is not referenced.
!> \endverbatim
!>
!> \param[in] LDVR
!> \verbatim
!>          LDVR is INTEGER
!>          The leading dimension of the array VR. LDVR >= 1;
!>          If JOB = 'E' or 'B', LDVR >= N.
!> \endverbatim
!>
!> \param[out] S
!> \verbatim
!>          S is DOUBLE PRECISION array, dimension (MM)
!>          If JOB = 'E' or 'B', the reciprocal condition numbers of the
!>          selected eigenvalues, stored in consecutive elements of the
!>          array.
!>          If JOB = 'V', S is not referenced.
!> \endverbatim
!>
!> \param[out] DIF
!> \verbatim
!>          DIF is DOUBLE PRECISION array, dimension (MM)
!>          If JOB = 'V' or 'B', the estimated reciprocal condition
!>          numbers of the selected eigenvectors, stored in consecutive
!>          elements of the array.
!>          If the eigenvalues cannot be reordered to compute DIF(j),
!>          DIF(j) is set to 0; this can only occur when the true value
!>          would be very small anyway.
!>          For each eigenvalue/vector specified by SELECT, DIF stores
!>          a Frobenius norm-based estimate of Difl.
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
!>          the specified condition numbers; for each selected eigenvalue
!>          one element is used. If HOWMNY = 'A', M is set to N.
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
!>          The dimension of the array WORK. LWORK >= max(1,N).
!>          If JOB = 'V' or 'B', LWORK >= max(1,2*N*N).
!> \endverbatim
!>
!> \param[out] IWORK
!> \verbatim
!>          IWORK is INTEGER array, dimension (N+2)
!>          If JOB = 'E', IWORK is not referenced.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0: Successful exit
!>          < 0: If INFO = -i, the i-th argument had an illegal value
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
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  The reciprocal of the condition number of the i-th generalized
!>  eigenvalue w = (a, b) is defined as
!>
!>          S(I) = (|v**HAu|**2 + |v**HBu|**2)**(1/2) / (norm(u)*norm(v))
!>
!>  where u and v are the right and left eigenvectors of (A, B)
!>  corresponding to w; |z| denotes the absolute value of the complex
!>  number, and norm(u) denotes the 2-norm of the vector u. The pair
!>  (a, b) corresponds to an eigenvalue w = a/b (= v**HAu/v**HBu) of the
!>  matrix pair (A, B). If both a and b equal zero, then (A,B) is
!>  singular and S(I) = -1 is returned.
!>
!>  An approximate error bound on the chordal distance between the i-th
!>  computed generalized eigenvalue w and the corresponding exact
!>  eigenvalue lambda is
!>
!>          chord(w, lambda) <=   EPS * norm(A, B) / S(I),
!>
!>  where EPS is the machine precision.
!>
!>  The reciprocal of the condition number of the right eigenvector u
!>  and left eigenvector v corresponding to the generalized eigenvalue w
!>  is defined as follows. Suppose
!>
!>                   (A, B) = ( a   *  ) ( b  *  )  1
!>                            ( 0  A22 ),( 0 B22 )  n-1
!>                              1  n-1     1 n-1
!>
!>  Then the reciprocal condition number DIF(I) is
!>
!>          Difl[(a, b), (A22, B22)]  = sigma-min( Zl )
!>
!>  where sigma-min(Zl) denotes the smallest singular value of
!>
!>         Zl = [ kron(a, In-1) -kron(1, A22) ]
!>              [ kron(b, In-1) -kron(1, B22) ].
!>
!>  Here In-1 is the identity matrix of size n-1 and X**H is the conjugate
!>  transpose of X. kron(X, Y) is the Kronecker product between the
!>  matrices X and Y.
!>
!>  We approximate the smallest singular value of Zl with an upper
!>  bound. This is done by ZLATDF.
!>
!>  An approximate error bound for a computed eigenvector VL(i) or
!>  VR(i) is given by
!>
!>                      EPS * norm(A, B) / DIF(i).
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
!>      Estimation: Theory, Algorithms and Software, Report
!>      UMINF - 94.04, Department of Computing Science, Umea University,
!>      S-901 87 Umea, Sweden, 1994. Also as LAPACK Working Note 87.
!>      To appear in Numerical Algorithms, 1996.
!>
!>  [3] B. Kagstrom and P. Poromaa, LAPACK-Style Algorithms and Software
!>      for Solving the Generalized Sylvester Equation and Estimating the
!>      Separation between Regular Matrix Pairs, Report UMINF - 93.23,
!>      Department of Computing Science, Umea University, S-901 87 Umea,
!>      Sweden, December 1993, Revised April 1994, Also as LAPACK Working
!>      Note 75.
!>      To appear in ACM Trans. on Math. Software, Vol 22, No 1, 1996.
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE ZTGSNA(Job,Howmny,Select,N,A,Lda,B,Ldb,Vl,Ldvl,Vr,Ldvr,&
     &                  S,Dif,Mm,M,Work,Lwork,Iwork,Info)
      USE F77KINDS                        
      USE S_DLABAD
      USE S_DLAMCH
      USE S_DLAPY2
      USE S_DZNRM2
      USE S_LSAME
      USE S_XERBLA
      USE S_ZDOTC
      USE S_ZGEMV
      USE S_ZLACPY
      USE S_ZTGEXC
      USE S_ZTGSYL
      IMPLICIT NONE
!*--ZTGSNA326
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
      INTEGER , PARAMETER  ::  IDIFJB = 3
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Job
      CHARACTER :: Howmny
      LOGICAL , INTENT(IN) , DIMENSION(*) :: Select
      INTEGER :: N
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX(CX16KIND) , DIMENSION(Ldvl,*) :: Vl
      INTEGER , INTENT(IN) :: Ldvl
      COMPLEX(CX16KIND) , DIMENSION(Ldvr,*) :: Vr
      INTEGER , INTENT(IN) :: Ldvr
      REAL(R8KIND) , INTENT(OUT) , DIMENSION(*) :: S
      REAL(R8KIND) , DIMENSION(*) :: Dif
      INTEGER , INTENT(IN) :: Mm
      INTEGER , INTENT(INOUT) :: M
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      REAL(R8KIND) :: bignum , cond , eps , lnrm , rnrm , scale , smlnum
      COMPLEX(CX16KIND) , DIMENSION(1) :: dummy , dummy1
      INTEGER :: i , ierr , ifst , ilst , k , ks , lwmin , n1 , n2
      LOGICAL :: lquery , somcon , wantbh , wantdf , wants
      COMPLEX(CX16KIND) :: yhax , yhbx
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
            DO k = 1 , N
               IF ( Select(k) ) M = M + 1
            ENDDO
         ELSE
            M = N
         ENDIF
!
         IF ( N==0 ) THEN
            lwmin = 1
         ELSEIF ( LSAME(Job,'V') .OR. LSAME(Job,'B') ) THEN
            lwmin = 2*N*N
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
         CALL XERBLA('ZTGSNA',-Info)
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
      eps = DLAMCH('P')
      smlnum = DLAMCH('S')/eps
      bignum = ONE/smlnum
      CALL DLABAD(smlnum,bignum)
      ks = 0
      DO k = 1 , N
!
!        Determine whether condition numbers are required for the k-th
!        eigenpair.
!
         IF ( somcon ) THEN
            IF ( .NOT.Select(k) ) CYCLE
         ENDIF
!
         ks = ks + 1
!
         IF ( wants ) THEN
!
!           Compute the reciprocal condition number of the k-th
!           eigenvalue.
!
            rnrm = DZNRM2(N,Vr(1,ks),1)
            lnrm = DZNRM2(N,Vl(1,ks),1)
            CALL ZGEMV('N',N,N,DCMPLX(ONE,ZERO),A,Lda,Vr(1,ks),1,       &
     &                 DCMPLX(ZERO,ZERO),Work,1)
            yhax = ZDOTC(N,Work,1,Vl(1,ks),1)
            CALL ZGEMV('N',N,N,DCMPLX(ONE,ZERO),B,Ldb,Vr(1,ks),1,       &
     &                 DCMPLX(ZERO,ZERO),Work,1)
            yhbx = ZDOTC(N,Work,1,Vl(1,ks),1)
            cond = DLAPY2(ABS(yhax),ABS(yhbx))
            IF ( cond==ZERO ) THEN
               S(ks) = -ONE
            ELSE
               S(ks) = cond/(rnrm*lnrm)
            ENDIF
         ENDIF
!
         IF ( wantdf ) THEN
            IF ( N==1 ) THEN
               Dif(ks) = DLAPY2(ABS(A(1,1)),ABS(B(1,1)))
            ELSE
!
!              Estimate the reciprocal condition number of the k-th
!              eigenvectors.
!
!              Copy the matrix (A, B) to the array WORK and move the
!              (k,k)th pair to the (1,1) position.
!
               CALL ZLACPY('Full',N,N,A,Lda,Work,N)
               CALL ZLACPY('Full',N,N,B,Ldb,Work(N*N+1),N)
               ifst = k
               ilst = 1
!
               CALL ZTGEXC(.FALSE.,.FALSE.,N,Work,N,Work(N*N+1),N,dummy,&
     &                     1,dummy1,1,ifst,ilst,ierr)
!
               IF ( ierr>0 ) THEN
!
!                 Ill-conditioned problem - swap rejected.
!
                  Dif(ks) = ZERO
               ELSE
!
!                 Reordering successful, solve generalized Sylvester
!                 equation for R and L,
!                            A22 * R - L * A11 = A12
!                            B22 * R - L * B11 = B12,
!                 and compute estimate of Difl[(A11,B11), (A22, B22)].
!
                  n1 = 1
                  n2 = N - n1
                  i = N*N + 1
                  CALL ZTGSYL('N',IDIFJB,n2,n1,Work(N*n1+n1+1),N,Work,N,&
     &                        Work(n1+1),N,Work(N*n1+n1+i),N,Work(i),N, &
     &                        Work(n1+i),N,scale,Dif(ks),dummy,1,Iwork, &
     &                        ierr)
               ENDIF
            ENDIF
         ENDIF
!
      ENDDO
      Work(1) = lwmin
!
!     End of ZTGSNA
!
      END SUBROUTINE ZTGSNA
