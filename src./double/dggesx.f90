!*==dggesx.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief <b> DGGESX computes the eigenvalues, the Schur form, and, optionally, the matrix of Schur vectors for GE matrices</b>
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DGGESX + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dggesx.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dggesx.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dggesx.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DGGESX( JOBVSL, JOBVSR, SORT, SELCTG, SENSE, N, A, LDA,
!                          B, LDB, SDIM, ALPHAR, ALPHAI, BETA, VSL, LDVSL,
!                          VSR, LDVSR, RCONDE, RCONDV, WORK, LWORK, IWORK,
!                          LIWORK, BWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          JOBVSL, JOBVSR, SENSE, SORT
!       INTEGER            INFO, LDA, LDB, LDVSL, LDVSR, LIWORK, LWORK, N,
!      $                   SDIM
!       ..
!       .. Array Arguments ..
!       LOGICAL            BWORK( * )
!       INTEGER            IWORK( * )
!       DOUBLE PRECISION   A( LDA, * ), ALPHAI( * ), ALPHAR( * ),
!      $                   B( LDB, * ), BETA( * ), RCONDE( 2 ),
!      $                   RCONDV( 2 ), VSL( LDVSL, * ), VSR( LDVSR, * ),
!      $                   WORK( * )
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
!> DGGESX computes for a pair of N-by-N real nonsymmetric matrices
!> (A,B), the generalized eigenvalues, the real Schur form (S,T), and,
!> optionally, the left and/or right matrices of Schur vectors (VSL and
!> VSR).  This gives the generalized Schur factorization
!>
!>      (A,B) = ( (VSL) S (VSR)**T, (VSL) T (VSR)**T )
!>
!> Optionally, it also orders the eigenvalues so that a selected cluster
!> of eigenvalues appears in the leading diagonal blocks of the upper
!> quasi-triangular matrix S and the upper triangular matrix T; computes
!> a reciprocal condition number for the average of the selected
!> eigenvalues (RCONDE); and computes a reciprocal condition number for
!> the right and left deflating subspaces corresponding to the selected
!> eigenvalues (RCONDV). The leading columns of VSL and VSR then form
!> an orthonormal basis for the corresponding left and right eigenspaces
!> (deflating subspaces).
!>
!> A generalized eigenvalue for a pair of matrices (A,B) is a scalar w
!> or a ratio alpha/beta = w, such that  A - w*B is singular.  It is
!> usually represented as the pair (alpha,beta), as there is a
!> reasonable interpretation for beta=0 or for both being zero.
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
!>          = 'S':  Eigenvalues are ordered (see SELCTG).
!> \endverbatim
!>
!> \param[in] SELCTG
!> \verbatim
!>          SELCTG is a LOGICAL FUNCTION of three DOUBLE PRECISION arguments
!>          SELCTG must be declared EXTERNAL in the calling subroutine.
!>          If SORT = 'N', SELCTG is not referenced.
!>          If SORT = 'S', SELCTG is used to select eigenvalues to sort
!>          to the top left of the Schur form.
!>          An eigenvalue (ALPHAR(j)+ALPHAI(j))/BETA(j) is selected if
!>          SELCTG(ALPHAR(j),ALPHAI(j),BETA(j)) is true; i.e. if either
!>          one of a complex conjugate pair of eigenvalues is selected,
!>          then both complex eigenvalues are selected.
!>          Note that a selected complex eigenvalue may no longer satisfy
!>          SELCTG(ALPHAR(j),ALPHAI(j),BETA(j)) = .TRUE. after ordering,
!>          since ordering may change the value of complex eigenvalues
!>          (especially if the eigenvalue is ill-conditioned), in this
!>          case INFO is set to N+3.
!> \endverbatim
!>
!> \param[in] SENSE
!> \verbatim
!>          SENSE is CHARACTER*1
!>          Determines which reciprocal condition numbers are computed.
!>          = 'N':  None are computed;
!>          = 'E':  Computed for average of selected eigenvalues only;
!>          = 'V':  Computed for selected deflating subspaces only;
!>          = 'B':  Computed for both.
!>          If SENSE = 'E', 'V', or 'B', SORT must equal 'S'.
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
!>          A is DOUBLE PRECISION array, dimension (LDA, N)
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
!>          B is DOUBLE PRECISION array, dimension (LDB, N)
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
!>          ALPHAR is DOUBLE PRECISION array, dimension (N)
!> \endverbatim
!>
!> \param[out] ALPHAI
!> \verbatim
!>          ALPHAI is DOUBLE PRECISION array, dimension (N)
!> \endverbatim
!>
!> \param[out] BETA
!> \verbatim
!>          BETA is DOUBLE PRECISION array, dimension (N)
!>          On exit, (ALPHAR(j) + ALPHAI(j)*i)/BETA(j), j=1,...,N, will
!>          be the generalized eigenvalues.  ALPHAR(j) + ALPHAI(j)*i
!>          and BETA(j),j=1,...,N  are the diagonals of the complex Schur
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
!>          VSL is DOUBLE PRECISION array, dimension (LDVSL,N)
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
!>          VSR is DOUBLE PRECISION array, dimension (LDVSR,N)
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
!> \param[out] RCONDE
!> \verbatim
!>          RCONDE is DOUBLE PRECISION array, dimension ( 2 )
!>          If SENSE = 'E' or 'B', RCONDE(1) and RCONDE(2) contain the
!>          reciprocal condition numbers for the average of the selected
!>          eigenvalues.
!>          Not referenced if SENSE = 'N' or 'V'.
!> \endverbatim
!>
!> \param[out] RCONDV
!> \verbatim
!>          RCONDV is DOUBLE PRECISION array, dimension ( 2 )
!>          If SENSE = 'V' or 'B', RCONDV(1) and RCONDV(2) contain the
!>          reciprocal condition numbers for the selected deflating
!>          subspaces.
!>          Not referenced if SENSE = 'N' or 'E'.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))
!>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The dimension of the array WORK.
!>          If N = 0, LWORK >= 1, else if SENSE = 'E', 'V', or 'B',
!>          LWORK >= max( 8*N, 6*N+16, 2*SDIM*(N-SDIM) ), else
!>          LWORK >= max( 8*N, 6*N+16 ).
!>          Note that 2*SDIM*(N-SDIM) <= N*N/2.
!>          Note also that an error is only returned if
!>          LWORK < max( 8*N, 6*N+16), but if SENSE = 'E' or 'V' or 'B'
!>          this may not be large enough.
!>
!>          If LWORK = -1, then a workspace query is assumed; the routine
!>          only calculates the bound on the optimal size of the WORK
!>          array and the minimum size of the IWORK array, returns these
!>          values as the first entries of the WORK and IWORK arrays, and
!>          no error message related to LWORK or LIWORK is issued by
!>          XERBLA.
!> \endverbatim
!>
!> \param[out] IWORK
!> \verbatim
!>          IWORK is INTEGER array, dimension (MAX(1,LIWORK))
!>          On exit, if INFO = 0, IWORK(1) returns the minimum LIWORK.
!> \endverbatim
!>
!> \param[in] LIWORK
!> \verbatim
!>          LIWORK is INTEGER
!>          The dimension of the array IWORK.
!>          If SENSE = 'N' or N = 0, LIWORK >= 1, otherwise
!>          LIWORK >= N+6.
!>
!>          If LIWORK = -1, then a workspace query is assumed; the
!>          routine only calculates the bound on the optimal size of the
!>          WORK array and the minimum size of the IWORK array, returns
!>          these values as the first entries of the WORK and IWORK
!>          arrays, and no error message related to LWORK or LIWORK is
!>          issued by XERBLA.
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
!>          > N:  =N+1: other than QZ iteration failed in DHGEQZ
!>                =N+2: after reordering, roundoff changed values of
!>                      some complex eigenvalues so that leading
!>                      eigenvalues in the Generalized Schur form no
!>                      longer satisfy SELCTG=.TRUE.  This could also
!>                      be caused due to scaling.
!>                =N+3: reordering failed in DTGSEN.
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
!> \date June 2017
!
!> \ingroup doubleGEeigen
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  An approximate (asymptotic) bound on the average absolute error of
!>  the selected eigenvalues is
!>
!>       EPS * norm((A, B)) / RCONDE( 1 ).
!>
!>  An approximate (asymptotic) bound on the maximum angular error in
!>  the computed deflating subspaces is
!>
!>       EPS * norm((A, B)) / RCONDV( 2 ).
!>
!>  See LAPACK User's Guide, section 4.11 for more information.
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE DGGESX(Jobvsl,Jobvsr,Sort,SELCTG,Sense,N,A,Lda,B,Ldb,  &
     &                  Sdim,Alphar,Alphai,Beta,Vsl,Ldvsl,Vsr,Ldvsr,    &
     &                  Rconde,Rcondv,Work,Lwork,Iwork,Liwork,Bwork,    &
     &                  Info)
      USE F77KINDS                        
      USE S_DGEQRF
      USE S_DGGBAK
      USE S_DGGBAL
      USE S_DGGHRD
      USE S_DHGEQZ
      USE S_DLABAD
      USE S_DLACPY
      USE S_DLAMCH
      USE S_DLANGE
      USE S_DLASCL
      USE S_DLASET
      USE S_DORGQR
      USE S_DORMQR
      USE S_DTGSEN
      USE S_ILAENV
      USE S_LSAME
      USE S_XERBLA
      IMPLICIT NONE
!*--DGGESX387
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Jobvsl
      CHARACTER :: Jobvsr
      CHARACTER :: Sort
      LOGICAL , EXTERNAL :: SELCTG
      CHARACTER :: Sense
      INTEGER :: N
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Sdim
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Alphar
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Alphai
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Beta
      REAL(R8KIND) , DIMENSION(Ldvsl,*) :: Vsl
      INTEGER :: Ldvsl
      REAL(R8KIND) , DIMENSION(Ldvsr,*) :: Vsr
      INTEGER :: Ldvsr
      REAL(R8KIND) , INTENT(OUT) , DIMENSION(2) :: Rconde
      REAL(R8KIND) , INTENT(OUT) , DIMENSION(2) :: Rcondv
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER :: Liwork
      LOGICAL , DIMENSION(*) :: Bwork
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      REAL(R8KIND) :: anrm , anrmto , bignum , bnrm , bnrmto , eps ,    &
     &                pl , pr , safmax , safmin , smlnum
      LOGICAL :: cursl , ilascl , ilbscl , ilvsl , ilvsr , lastsl ,     &
     &           lquery , lst2sl , wantsb , wantse , wantsn , wantst ,  &
     &           wantsv
      REAL(R8KIND) , DIMENSION(2) :: dif
      INTEGER :: i , icols , ierr , ihi , ijob , ijobvl , ijobvr ,      &
     &           ileft , ilo , ip , iright , irows , itau , iwrk ,      &
     &           liwmin , lwrk , maxwrk , minwrk
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
      wantsn = LSAME(Sense,'N')
      wantse = LSAME(Sense,'E')
      wantsv = LSAME(Sense,'V')
      wantsb = LSAME(Sense,'B')
      lquery = (Lwork==-1 .OR. Liwork==-1)
      IF ( wantsn ) THEN
         ijob = 0
      ELSEIF ( wantse ) THEN
         ijob = 1
      ELSEIF ( wantsv ) THEN
         ijob = 2
      ELSEIF ( wantsb ) THEN
         ijob = 4
      ENDIF
!
!     Test the input arguments
!
      Info = 0
      IF ( ijobvl<=0 ) THEN
         Info = -1
      ELSEIF ( ijobvr<=0 ) THEN
         Info = -2
      ELSEIF ( (.NOT.wantst) .AND. (.NOT.LSAME(Sort,'N')) ) THEN
         Info = -3
      ELSEIF ( .NOT.(wantsn .OR. wantse .OR. wantsv .OR. wantsb) .OR.   &
     &         (.NOT.wantst .AND. .NOT.wantsn) ) THEN
         Info = -5
      ELSEIF ( N<0 ) THEN
         Info = -6
      ELSEIF ( Lda<MAX(1,N) ) THEN
         Info = -8
      ELSEIF ( Ldb<MAX(1,N) ) THEN
         Info = -10
      ELSEIF ( Ldvsl<1 .OR. (ilvsl .AND. Ldvsl<N) ) THEN
         Info = -16
      ELSEIF ( Ldvsr<1 .OR. (ilvsr .AND. Ldvsr<N) ) THEN
         Info = -18
      ENDIF
!
!     Compute workspace
!      (Note: Comments in the code beginning "Workspace:" describe the
!       minimal amount of workspace needed at that point in the code,
!       as well as the preferred amount for good performance.
!       NB refers to the optimal block size for the immediately
!       following subroutine, as returned by ILAENV.)
!
      IF ( Info==0 ) THEN
         IF ( N>0 ) THEN
            minwrk = MAX(8*N,6*N+16)
            maxwrk = minwrk - N + N*ILAENV(1,'DGEQRF',' ',N,1,N,0)
            maxwrk = MAX(maxwrk,                                        &
     &               minwrk-N+N*ILAENV(1,'DORMQR',' ',N,1,N,-1))
            IF ( ilvsl ) maxwrk = MAX(maxwrk,minwrk-N+N*ILAENV(1,       &
     &                            'DORGQR',' ',N,1,N,-1))
            lwrk = maxwrk
            IF ( ijob>=1 ) lwrk = MAX(lwrk,N*N/2)
         ELSE
            minwrk = 1
            maxwrk = 1
            lwrk = 1
         ENDIF
         Work(1) = lwrk
         IF ( wantsn .OR. N==0 ) THEN
            liwmin = 1
         ELSE
            liwmin = N + 6
         ENDIF
         Iwork(1) = liwmin
!
         IF ( Lwork<minwrk .AND. .NOT.lquery ) THEN
            Info = -22
         ELSEIF ( Liwork<liwmin .AND. .NOT.lquery ) THEN
            Info = -24
         ENDIF
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('DGGESX',-Info)
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
      eps = DLAMCH('P')
      safmin = DLAMCH('S')
      safmax = ONE/safmin
      CALL DLABAD(safmin,safmax)
      smlnum = SQRT(safmin)/eps
      bignum = ONE/smlnum
!
!     Scale A if max element outside range [SMLNUM,BIGNUM]
!
      anrm = DLANGE('M',N,N,A,Lda,Work)
      ilascl = .FALSE.
      IF ( anrm>ZERO .AND. anrm<smlnum ) THEN
         anrmto = smlnum
         ilascl = .TRUE.
      ELSEIF ( anrm>bignum ) THEN
         anrmto = bignum
         ilascl = .TRUE.
      ENDIF
      IF ( ilascl ) CALL DLASCL('G',0,0,anrm,anrmto,N,N,A,Lda,ierr)
!
!     Scale B if max element outside range [SMLNUM,BIGNUM]
!
      bnrm = DLANGE('M',N,N,B,Ldb,Work)
      ilbscl = .FALSE.
      IF ( bnrm>ZERO .AND. bnrm<smlnum ) THEN
         bnrmto = smlnum
         ilbscl = .TRUE.
      ELSEIF ( bnrm>bignum ) THEN
         bnrmto = bignum
         ilbscl = .TRUE.
      ENDIF
      IF ( ilbscl ) CALL DLASCL('G',0,0,bnrm,bnrmto,N,N,B,Ldb,ierr)
!
!     Permute the matrix to make it more nearly triangular
!     (Workspace: need 6*N + 2*N for permutation parameters)
!
      ileft = 1
      iright = N + 1
      iwrk = iright + N
      CALL DGGBAL('P',N,A,Lda,B,Ldb,ilo,ihi,Work(ileft),Work(iright),   &
     &            Work(iwrk),ierr)
!
!     Reduce B to triangular form (QR decomposition of B)
!     (Workspace: need N, prefer N*NB)
!
      irows = ihi + 1 - ilo
      icols = N + 1 - ilo
      itau = iwrk
      iwrk = itau + irows
      CALL DGEQRF(irows,icols,B(ilo,ilo),Ldb,Work(itau),Work(iwrk),     &
     &            Lwork+1-iwrk,ierr)
!
!     Apply the orthogonal transformation to matrix A
!     (Workspace: need N, prefer N*NB)
!
      CALL DORMQR('L','T',irows,icols,irows,B(ilo,ilo),Ldb,Work(itau),  &
     &            A(ilo,ilo),Lda,Work(iwrk),Lwork+1-iwrk,ierr)
!
!     Initialize VSL
!     (Workspace: need N, prefer N*NB)
!
      IF ( ilvsl ) THEN
         CALL DLASET('Full',N,N,ZERO,ONE,Vsl,Ldvsl)
         IF ( irows>1 ) CALL DLACPY('L',irows-1,irows-1,B(ilo+1,ilo),   &
     &                              Ldb,Vsl(ilo+1,ilo),Ldvsl)
         CALL DORGQR(irows,irows,irows,Vsl(ilo,ilo),Ldvsl,Work(itau),   &
     &               Work(iwrk),Lwork+1-iwrk,ierr)
      ENDIF
!
!     Initialize VSR
!
      IF ( ilvsr ) CALL DLASET('Full',N,N,ZERO,ONE,Vsr,Ldvsr)
!
!     Reduce to generalized Hessenberg form
!     (Workspace: none needed)
!
      CALL DGGHRD(Jobvsl,Jobvsr,N,ilo,ihi,A,Lda,B,Ldb,Vsl,Ldvsl,Vsr,    &
     &            Ldvsr,ierr)
!
      Sdim = 0
!
!     Perform QZ algorithm, computing Schur vectors if desired
!     (Workspace: need N)
!
      iwrk = itau
      CALL DHGEQZ('S',Jobvsl,Jobvsr,N,ilo,ihi,A,Lda,B,Ldb,Alphar,Alphai,&
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
!     Sort eigenvalues ALPHA/BETA and compute the reciprocal of
!     condition number(s)
!     (Workspace: If IJOB >= 1, need MAX( 8*(N+1), 2*SDIM*(N-SDIM) )
!                 otherwise, need 8*(N+1) )
!
      IF ( wantst ) THEN
!
!        Undo scaling on eigenvalues before SELCTGing
!
         IF ( ilascl ) THEN
            CALL DLASCL('G',0,0,anrmto,anrm,N,1,Alphar,N,ierr)
            CALL DLASCL('G',0,0,anrmto,anrm,N,1,Alphai,N,ierr)
         ENDIF
         IF ( ilbscl ) CALL DLASCL('G',0,0,bnrmto,bnrm,N,1,Beta,N,ierr)
!
!        Select eigenvalues
!
         DO i = 1 , N
            Bwork(i) = SELCTG(Alphar(i),Alphai(i),Beta(i))
         ENDDO
!
!        Reorder eigenvalues, transform Generalized Schur vectors, and
!        compute reciprocal condition numbers
!
         CALL DTGSEN(ijob,ilvsl,ilvsr,Bwork,N,A,Lda,B,Ldb,Alphar,Alphai,&
     &               Beta,Vsl,Ldvsl,Vsr,Ldvsr,Sdim,pl,pr,dif,Work(iwrk),&
     &               Lwork-iwrk+1,Iwork,Liwork,ierr)
!
         IF ( ijob>=1 ) maxwrk = MAX(maxwrk,2*Sdim*(N-Sdim))
         IF ( ierr==-22 ) THEN
!
!            not enough real workspace
!
            Info = -22
         ELSE
            IF ( ijob==1 .OR. ijob==4 ) THEN
               Rconde(1) = pl
               Rconde(2) = pr
            ENDIF
            IF ( ijob==2 .OR. ijob==4 ) THEN
               Rcondv(1) = dif(1)
               Rcondv(2) = dif(2)
            ENDIF
            IF ( ierr==1 ) Info = N + 3
         ENDIF
!
      ENDIF
!
!     Apply permutation to VSL and VSR
!     (Workspace: none needed)
!
      IF ( ilvsl ) CALL DGGBAK('P','L',N,ilo,ihi,Work(ileft),           &
     &                         Work(iright),N,Vsl,Ldvsl,ierr)
!
      IF ( ilvsr ) CALL DGGBAK('P','R',N,ilo,ihi,Work(ileft),           &
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
         CALL DLASCL('H',0,0,anrmto,anrm,N,N,A,Lda,ierr)
         CALL DLASCL('G',0,0,anrmto,anrm,N,1,Alphar,N,ierr)
         CALL DLASCL('G',0,0,anrmto,anrm,N,1,Alphai,N,ierr)
      ENDIF
!
      IF ( ilbscl ) THEN
         CALL DLASCL('U',0,0,bnrmto,bnrm,N,N,B,Ldb,ierr)
         CALL DLASCL('G',0,0,bnrmto,bnrm,N,1,Beta,N,ierr)
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
 100  Work(1) = maxwrk
      Iwork(1) = liwmin
!
!
!     End of DGGESX
!
      END SUBROUTINE DGGESX
