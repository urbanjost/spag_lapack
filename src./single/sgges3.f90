!*==sgges3.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief <b> SGGES3 computes the eigenvalues, the Schur form, and, optionally, the matrix of Schur vectors for GE matrices (blocked algorithm)</b>
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SGGES3 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgges3.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgges3.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgges3.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SGGES3( JOBVSL, JOBVSR, SORT, SELCTG, N, A, LDA, B,
!      $                   LDB, SDIM, ALPHAR, ALPHAI, BETA, VSL, LDVSL,
!      $                   VSR, LDVSR, WORK, LWORK, BWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          JOBVSL, JOBVSR, SORT
!       INTEGER            INFO, LDA, LDB, LDVSL, LDVSR, LWORK, N, SDIM
!       ..
!       .. Array Arguments ..
!       LOGICAL            BWORK( * )
!       REAL               A( LDA, * ), ALPHAI( * ), ALPHAR( * ),
!      $                   B( LDB, * ), BETA( * ), VSL( LDVSL, * ),
!      $                   VSR( LDVSR, * ), WORK( * )
!       ..
!       .. Function Arguments ..
!       LOGICAL            SELCTG
!       EXTERNAL           SELCTG
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SGGES3 computes for a pair of N-by-N real nonsymmetric matrices (A,B),
!> the generalized eigenvalues, the generalized real Schur form (S,T),
!> optionally, the left and/or right matrices of Schur vectors (VSL and
!> VSR). This gives the generalized Schur factorization
!>
!>          (A,B) = ( (VSL)*S*(VSR)**T, (VSL)*T*(VSR)**T )
!>
!> Optionally, it also orders the eigenvalues so that a selected cluster
!> of eigenvalues appears in the leading diagonal blocks of the upper
!> quasi-triangular matrix S and the upper triangular matrix T.The
!> leading columns of VSL and VSR then form an orthonormal basis for the
!> corresponding left and right eigenspaces (deflating subspaces).
!>
!> (If only the generalized eigenvalues are needed, use the driver
!> SGGEV instead, which is faster.)
!>
!> A generalized eigenvalue for a pair of matrices (A,B) is a scalar w
!> or a ratio alpha/beta = w, such that  A - w*B is singular.  It is
!> usually represented as the pair (alpha,beta), as there is a
!> reasonable interpretation for beta=0 or both being zero.
!>
!> A pair of matrices (S,T) is in generalized real Schur form if T is
!> upper triangular with non-negative diagonal and S is block upper
!> triangular with 1-by-1 and 2-by-2 blocks.  1-by-1 blocks correspond
!> to real generalized eigenvalues, while 2-by-2 blocks of S will be
!> "standardized" by making the corresponding elements of T have the
!> form:
!>         [  a  0  ]
!>         [  0  b  ]
!>
!> and the pair of corresponding 2-by-2 blocks in S and T will have a
!> complex conjugate pair of generalized eigenvalues.
!>
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] JOBVSL
!> \verbatim
!>          JOBVSL is CHARACTER*1
!>          = 'N':  do not compute the left Schur vectors;
!>          = 'V':  compute the left Schur vectors.
!> \endverbatim
!>
!> \param[in] JOBVSR
!> \verbatim
!>          JOBVSR is CHARACTER*1
!>          = 'N':  do not compute the right Schur vectors;
!>          = 'V':  compute the right Schur vectors.
!> \endverbatim
!>
!> \param[in] SORT
!> \verbatim
!>          SORT is CHARACTER*1
!>          Specifies whether or not to order the eigenvalues on the
!>          diagonal of the generalized Schur form.
!>          = 'N':  Eigenvalues are not ordered;
!>          = 'S':  Eigenvalues are ordered (see SELCTG);
!> \endverbatim
!>
!> \param[in] SELCTG
!> \verbatim
!>          SELCTG is a LOGICAL FUNCTION of three REAL arguments
!>          SELCTG must be declared EXTERNAL in the calling subroutine.
!>          If SORT = 'N', SELCTG is not referenced.
!>          If SORT = 'S', SELCTG is used to select eigenvalues to sort
!>          to the top left of the Schur form.
!>          An eigenvalue (ALPHAR(j)+ALPHAI(j))/BETA(j) is selected if
!>          SELCTG(ALPHAR(j),ALPHAI(j),BETA(j)) is true; i.e. if either
!>          one of a complex conjugate pair of eigenvalues is selected,
!>          then both complex eigenvalues are selected.
!>
!>          Note that in the ill-conditioned case, a selected complex
!>          eigenvalue may no longer satisfy SELCTG(ALPHAR(j),ALPHAI(j),
!>          BETA(j)) = .TRUE. after ordering. INFO is to be set to N+2
!>          in this case.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrices A, B, VSL, and VSR.  N >= 0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is REAL array, dimension (LDA, N)
!>          On entry, the first of the pair of matrices.
!>          On exit, A has been overwritten by its generalized Schur
!>          form S.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of A.  LDA >= max(1,N).
!> \endverbatim
!>
!> \param[in,out] B
!> \verbatim
!>          B is REAL array, dimension (LDB, N)
!>          On entry, the second of the pair of matrices.
!>          On exit, B has been overwritten by its generalized Schur
!>          form T.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of B.  LDB >= max(1,N).
!> \endverbatim
!>
!> \param[out] SDIM
!> \verbatim
!>          SDIM is INTEGER
!>          If SORT = 'N', SDIM = 0.
!>          If SORT = 'S', SDIM = number of eigenvalues (after sorting)
!>          for which SELCTG is true.  (Complex conjugate pairs for which
!>          SELCTG is true for either eigenvalue count as 2.)
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
!>          On exit, (ALPHAR(j) + ALPHAI(j)*i)/BETA(j), j=1,...,N, will
!>          be the generalized eigenvalues.  ALPHAR(j) + ALPHAI(j)*i,
!>          and  BETA(j),j=1,...,N are the diagonals of the complex Schur
!>          form (S,T) that would result if the 2-by-2 diagonal blocks of
!>          the real Schur form of (A,B) were further reduced to
!>          triangular form using 2-by-2 complex unitary transformations.
!>          If ALPHAI(j) is zero, then the j-th eigenvalue is real; if
!>          positive, then the j-th and (j+1)-st eigenvalues are a
!>          complex conjugate pair, with ALPHAI(j+1) negative.
!>
!>          Note: the quotients ALPHAR(j)/BETA(j) and ALPHAI(j)/BETA(j)
!>          may easily over- or underflow, and BETA(j) may even be zero.
!>          Thus, the user should avoid naively computing the ratio.
!>          However, ALPHAR and ALPHAI will be always less than and
!>          usually comparable with norm(A) in magnitude, and BETA always
!>          less than and usually comparable with norm(B).
!> \endverbatim
!>
!> \param[out] VSL
!> \verbatim
!>          VSL is REAL array, dimension (LDVSL,N)
!>          If JOBVSL = 'V', VSL will contain the left Schur vectors.
!>          Not referenced if JOBVSL = 'N'.
!> \endverbatim
!>
!> \param[in] LDVSL
!> \verbatim
!>          LDVSL is INTEGER
!>          The leading dimension of the matrix VSL. LDVSL >=1, and
!>          if JOBVSL = 'V', LDVSL >= N.
!> \endverbatim
!>
!> \param[out] VSR
!> \verbatim
!>          VSR is REAL array, dimension (LDVSR,N)
!>          If JOBVSR = 'V', VSR will contain the right Schur vectors.
!>          Not referenced if JOBVSR = 'N'.
!> \endverbatim
!>
!> \param[in] LDVSR
!> \verbatim
!>          LDVSR is INTEGER
!>          The leading dimension of the matrix VSR. LDVSR >= 1, and
!>          if JOBVSR = 'V', LDVSR >= N.
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
!>          The dimension of the array WORK.
!>
!>          If LWORK = -1, then a workspace query is assumed; the routine
!>          only calculates the optimal size of the WORK array, returns
!>          this value as the first entry of the WORK array, and no error
!>          message related to LWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] BWORK
!> \verbatim
!>          BWORK is LOGICAL array, dimension (N)
!>          Not referenced if SORT = 'N'.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value.
!>          = 1,...,N:
!>                The QZ iteration failed.  (A,B) are not in Schur
!>                form, but ALPHAR(j), ALPHAI(j), and BETA(j) should
!>                be correct for j=INFO+1,...,N.
!>          > N:  =N+1: other than QZ iteration failed in SHGEQZ.
!>                =N+2: after reordering, roundoff changed values of
!>                      some complex eigenvalues so that leading
!>                      eigenvalues in the Generalized Schur form no
!>                      longer satisfy SELCTG=.TRUE.  This could also
!>                      be caused due to scaling.
!>                =N+3: reordering failed in STGSEN.
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
!> \date January 2015
!
!> \ingroup realGEeigen
!
!  =====================================================================
      SUBROUTINE SGGES3(Jobvsl,Jobvsr,Sort,SELCTG,N,A,Lda,B,Ldb,Sdim,   &
     &                  Alphar,Alphai,Beta,Vsl,Ldvsl,Vsr,Ldvsr,Work,    &
     &                  Lwork,Bwork,Info)
      USE S_LSAME
      USE S_SGEQRF
      USE S_SGGBAK
      USE S_SGGBAL
      USE S_SGGHD3
      USE S_SHGEQZ
      USE S_SLABAD
      USE S_SLACPY
      USE S_SLAMCH
      USE S_SLANGE
      USE S_SLASCL
      USE S_SLASET
      USE S_SORGQR
      USE S_SORMQR
      USE S_STGSEN
      USE S_XERBLA
      IMPLICIT NONE
!*--SGGES3302
!
! PARAMETER definitions rewritten by SPAG
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Jobvsl
      CHARACTER :: Jobvsr
      CHARACTER :: Sort
      LOGICAL , EXTERNAL :: SELCTG
      INTEGER :: N
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Sdim
      REAL , INTENT(INOUT) , DIMENSION(*) :: Alphar
      REAL , INTENT(INOUT) , DIMENSION(*) :: Alphai
      REAL , INTENT(INOUT) , DIMENSION(*) :: Beta
      REAL , DIMENSION(Ldvsl,*) :: Vsl
      INTEGER :: Ldvsl
      REAL , DIMENSION(Ldvsr,*) :: Vsr
      INTEGER :: Ldvsr
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      LOGICAL , DIMENSION(*) :: Bwork
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      REAL :: anrm , anrmto , bignum , bnrm , bnrmto , eps , pvsl ,     &
     &        pvsr , safmax , safmin , smlnum
      LOGICAL :: cursl , ilascl , ilbscl , ilvsl , ilvsr , lastsl ,     &
     &           lquery , lst2sl , wantst
      REAL , DIMENSION(2) :: dif
      INTEGER :: i , icols , ierr , ihi , ijobvl , ijobvr , ileft ,     &
     &           ilo , ip , iright , irows , itau , iwrk , lwkopt
      INTEGER , DIMENSION(1) :: idum
!
! End of declarations rewritten by SPAG
!
!     ..
!     .. Array Arguments ..
!     ..
!     .. Function Arguments ..
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
!     Decode the input arguments
!
      IF ( LSAME(Jobvsl,'N') ) THEN
         ijobvl = 1
         ilvsl = .FALSE.
      ELSEIF ( LSAME(Jobvsl,'V') ) THEN
         ijobvl = 2
         ilvsl = .TRUE.
      ELSE
         ijobvl = -1
         ilvsl = .FALSE.
      ENDIF
!
      IF ( LSAME(Jobvsr,'N') ) THEN
         ijobvr = 1
         ilvsr = .FALSE.
      ELSEIF ( LSAME(Jobvsr,'V') ) THEN
         ijobvr = 2
         ilvsr = .TRUE.
      ELSE
         ijobvr = -1
         ilvsr = .FALSE.
      ENDIF
!
      wantst = LSAME(Sort,'S')
!
!     Test the input arguments
!
      Info = 0
      lquery = (Lwork==-1)
      IF ( ijobvl<=0 ) THEN
         Info = -1
      ELSEIF ( ijobvr<=0 ) THEN
         Info = -2
      ELSEIF ( (.NOT.wantst) .AND. (.NOT.LSAME(Sort,'N')) ) THEN
         Info = -3
      ELSEIF ( N<0 ) THEN
         Info = -5
      ELSEIF ( Lda<MAX(1,N) ) THEN
         Info = -7
      ELSEIF ( Ldb<MAX(1,N) ) THEN
         Info = -9
      ELSEIF ( Ldvsl<1 .OR. (ilvsl .AND. Ldvsl<N) ) THEN
         Info = -15
      ELSEIF ( Ldvsr<1 .OR. (ilvsr .AND. Ldvsr<N) ) THEN
         Info = -17
      ELSEIF ( Lwork<6*N+16 .AND. .NOT.lquery ) THEN
         Info = -19
      ENDIF
!
!     Compute workspace
!
      IF ( Info==0 ) THEN
         CALL SGEQRF(N,N,B,Ldb,Work,Work,-1,ierr)
         lwkopt = MAX(6*N+16,3*N+INT(Work(1)))
         CALL SORMQR('L','T',N,N,N,B,Ldb,Work,A,Lda,Work,-1,ierr)
         lwkopt = MAX(lwkopt,3*N+INT(Work(1)))
         IF ( ilvsl ) THEN
            CALL SORGQR(N,N,N,Vsl,Ldvsl,Work,Work,-1,ierr)
            lwkopt = MAX(lwkopt,3*N+INT(Work(1)))
         ENDIF
         CALL SGGHD3(Jobvsl,Jobvsr,N,1,N,A,Lda,B,Ldb,Vsl,Ldvsl,Vsr,     &
     &               Ldvsr,Work,-1,ierr)
         lwkopt = MAX(lwkopt,3*N+INT(Work(1)))
         CALL SHGEQZ('S',Jobvsl,Jobvsr,N,1,N,A,Lda,B,Ldb,Alphar,Alphai, &
     &               Beta,Vsl,Ldvsl,Vsr,Ldvsr,Work,-1,ierr)
         lwkopt = MAX(lwkopt,2*N+INT(Work(1)))
         IF ( wantst ) THEN
            CALL STGSEN(0,ilvsl,ilvsr,Bwork,N,A,Lda,B,Ldb,Alphar,Alphai,&
     &                  Beta,Vsl,Ldvsl,Vsr,Ldvsr,Sdim,pvsl,pvsr,dif,    &
     &                  Work,-1,idum,1,ierr)
            lwkopt = MAX(lwkopt,2*N+INT(Work(1)))
         ENDIF
         Work(1) = lwkopt
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('SGGES3 ',-Info)
         RETURN
      ELSEIF ( lquery ) THEN
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( N==0 ) THEN
         Sdim = 0
         RETURN
      ENDIF
!
!     Get machine constants
!
      eps = SLAMCH('P')
      safmin = SLAMCH('S')
      safmax = ONE/safmin
      CALL SLABAD(safmin,safmax)
      smlnum = SQRT(safmin)/eps
      bignum = ONE/smlnum
!
!     Scale A if max element outside range [SMLNUM,BIGNUM]
!
      anrm = SLANGE('M',N,N,A,Lda,Work)
      ilascl = .FALSE.
      IF ( anrm>ZERO .AND. anrm<smlnum ) THEN
         anrmto = smlnum
         ilascl = .TRUE.
      ELSEIF ( anrm>bignum ) THEN
         anrmto = bignum
         ilascl = .TRUE.
      ENDIF
      IF ( ilascl ) CALL SLASCL('G',0,0,anrm,anrmto,N,N,A,Lda,ierr)
!
!     Scale B if max element outside range [SMLNUM,BIGNUM]
!
      bnrm = SLANGE('M',N,N,B,Ldb,Work)
      ilbscl = .FALSE.
      IF ( bnrm>ZERO .AND. bnrm<smlnum ) THEN
         bnrmto = smlnum
         ilbscl = .TRUE.
      ELSEIF ( bnrm>bignum ) THEN
         bnrmto = bignum
         ilbscl = .TRUE.
      ENDIF
      IF ( ilbscl ) CALL SLASCL('G',0,0,bnrm,bnrmto,N,N,B,Ldb,ierr)
!
!     Permute the matrix to make it more nearly triangular
!
      ileft = 1
      iright = N + 1
      iwrk = iright + N
      CALL SGGBAL('P',N,A,Lda,B,Ldb,ilo,ihi,Work(ileft),Work(iright),   &
     &            Work(iwrk),ierr)
!
!     Reduce B to triangular form (QR decomposition of B)
!
      irows = ihi + 1 - ilo
      icols = N + 1 - ilo
      itau = iwrk
      iwrk = itau + irows
      CALL SGEQRF(irows,icols,B(ilo,ilo),Ldb,Work(itau),Work(iwrk),     &
     &            Lwork+1-iwrk,ierr)
!
!     Apply the orthogonal transformation to matrix A
!
      CALL SORMQR('L','T',irows,icols,irows,B(ilo,ilo),Ldb,Work(itau),  &
     &            A(ilo,ilo),Lda,Work(iwrk),Lwork+1-iwrk,ierr)
!
!     Initialize VSL
!
      IF ( ilvsl ) THEN
         CALL SLASET('Full',N,N,ZERO,ONE,Vsl,Ldvsl)
         IF ( irows>1 ) CALL SLACPY('L',irows-1,irows-1,B(ilo+1,ilo),   &
     &                              Ldb,Vsl(ilo+1,ilo),Ldvsl)
         CALL SORGQR(irows,irows,irows,Vsl(ilo,ilo),Ldvsl,Work(itau),   &
     &               Work(iwrk),Lwork+1-iwrk,ierr)
      ENDIF
!
!     Initialize VSR
!
      IF ( ilvsr ) CALL SLASET('Full',N,N,ZERO,ONE,Vsr,Ldvsr)
!
!     Reduce to generalized Hessenberg form
!
      CALL SGGHD3(Jobvsl,Jobvsr,N,ilo,ihi,A,Lda,B,Ldb,Vsl,Ldvsl,Vsr,    &
     &            Ldvsr,Work(iwrk),Lwork+1-iwrk,ierr)
!
!     Perform QZ algorithm, computing Schur vectors if desired
!
      iwrk = itau
      CALL SHGEQZ('S',Jobvsl,Jobvsr,N,ilo,ihi,A,Lda,B,Ldb,Alphar,Alphai,&
     &            Beta,Vsl,Ldvsl,Vsr,Ldvsr,Work(iwrk),Lwork+1-iwrk,ierr)
      IF ( ierr/=0 ) THEN
         IF ( ierr>0 .AND. ierr<=N ) THEN
            Info = ierr
         ELSEIF ( ierr>N .AND. ierr<=2*N ) THEN
            Info = ierr - N
         ELSE
            Info = N + 1
         ENDIF
         GOTO 100
      ENDIF
!
!     Sort eigenvalues ALPHA/BETA if desired
!
      Sdim = 0
      IF ( wantst ) THEN
!
!        Undo scaling on eigenvalues before SELCTGing
!
         IF ( ilascl ) THEN
            CALL SLASCL('G',0,0,anrmto,anrm,N,1,Alphar,N,ierr)
            CALL SLASCL('G',0,0,anrmto,anrm,N,1,Alphai,N,ierr)
         ENDIF
         IF ( ilbscl ) CALL SLASCL('G',0,0,bnrmto,bnrm,N,1,Beta,N,ierr)
!
!        Select eigenvalues
!
         DO i = 1 , N
            Bwork(i) = SELCTG(Alphar(i),Alphai(i),Beta(i))
         ENDDO
!
         CALL STGSEN(0,ilvsl,ilvsr,Bwork,N,A,Lda,B,Ldb,Alphar,Alphai,   &
     &               Beta,Vsl,Ldvsl,Vsr,Ldvsr,Sdim,pvsl,pvsr,dif,       &
     &               Work(iwrk),Lwork-iwrk+1,idum,1,ierr)
         IF ( ierr==1 ) Info = N + 3
!
      ENDIF
!
!     Apply back-permutation to VSL and VSR
!
      IF ( ilvsl ) CALL SGGBAK('P','L',N,ilo,ihi,Work(ileft),           &
     &                         Work(iright),N,Vsl,Ldvsl,ierr)
!
      IF ( ilvsr ) CALL SGGBAK('P','R',N,ilo,ihi,Work(ileft),           &
     &                         Work(iright),N,Vsr,Ldvsr,ierr)
!
!     Check if unscaling would cause over/underflow, if so, rescale
!     (ALPHAR(I),ALPHAI(I),BETA(I)) so BETA(I) is on the order of
!     B(I,I) and ALPHAR(I) and ALPHAI(I) are on the order of A(I,I)
!
      IF ( ilascl ) THEN
         DO i = 1 , N
            IF ( Alphai(i)/=ZERO ) THEN
               IF ( (Alphar(i)/safmax)>(anrmto/anrm) .OR.               &
     &              (safmin/Alphar(i))>(anrm/anrmto) ) THEN
                  Work(1) = ABS(A(i,i)/Alphar(i))
                  Beta(i) = Beta(i)*Work(1)
                  Alphar(i) = Alphar(i)*Work(1)
                  Alphai(i) = Alphai(i)*Work(1)
               ELSEIF ( (Alphai(i)/safmax)>(anrmto/anrm) .OR.           &
     &                  (safmin/Alphai(i))>(anrm/anrmto) ) THEN
                  Work(1) = ABS(A(i,i+1)/Alphai(i))
                  Beta(i) = Beta(i)*Work(1)
                  Alphar(i) = Alphar(i)*Work(1)
                  Alphai(i) = Alphai(i)*Work(1)
               ENDIF
            ENDIF
         ENDDO
      ENDIF
!
      IF ( ilbscl ) THEN
         DO i = 1 , N
            IF ( Alphai(i)/=ZERO ) THEN
               IF ( (Beta(i)/safmax)>(bnrmto/bnrm) .OR. (safmin/Beta(i))&
     &              >(bnrm/bnrmto) ) THEN
                  Work(1) = ABS(B(i,i)/Beta(i))
                  Beta(i) = Beta(i)*Work(1)
                  Alphar(i) = Alphar(i)*Work(1)
                  Alphai(i) = Alphai(i)*Work(1)
               ENDIF
            ENDIF
         ENDDO
      ENDIF
!
!     Undo scaling
!
      IF ( ilascl ) THEN
         CALL SLASCL('H',0,0,anrmto,anrm,N,N,A,Lda,ierr)
         CALL SLASCL('G',0,0,anrmto,anrm,N,1,Alphar,N,ierr)
         CALL SLASCL('G',0,0,anrmto,anrm,N,1,Alphai,N,ierr)
      ENDIF
!
      IF ( ilbscl ) THEN
         CALL SLASCL('U',0,0,bnrmto,bnrm,N,N,B,Ldb,ierr)
         CALL SLASCL('G',0,0,bnrmto,bnrm,N,1,Beta,N,ierr)
      ENDIF
!
      IF ( wantst ) THEN
!
!        Check if reordering is correct
!
         lastsl = .TRUE.
         lst2sl = .TRUE.
         Sdim = 0
         ip = 0
         DO i = 1 , N
            cursl = SELCTG(Alphar(i),Alphai(i),Beta(i))
            IF ( Alphai(i)==ZERO ) THEN
               IF ( cursl ) Sdim = Sdim + 1
               ip = 0
               IF ( cursl .AND. .NOT.lastsl ) Info = N + 2
            ELSEIF ( ip==1 ) THEN
!
!                 Last eigenvalue of conjugate pair
!
               cursl = cursl .OR. lastsl
               lastsl = cursl
               IF ( cursl ) Sdim = Sdim + 2
               ip = -1
               IF ( cursl .AND. .NOT.lst2sl ) Info = N + 2
            ELSE
!
!                 First eigenvalue of conjugate pair
!
               ip = 1
            ENDIF
            lst2sl = lastsl
            lastsl = cursl
         ENDDO
!
      ENDIF
!
!
 100  Work(1) = lwkopt
!
!
!     End of SGGES3
!
      END SUBROUTINE SGGES3
