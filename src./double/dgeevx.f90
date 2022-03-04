!*==dgeevx.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief <b> DGEEVX computes the eigenvalues and, optionally, the left and/or right eigenvectors for GE matrices</b>
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DGEEVX + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgeevx.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgeevx.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgeevx.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DGEEVX( BALANC, JOBVL, JOBVR, SENSE, N, A, LDA, WR, WI,
!                          VL, LDVL, VR, LDVR, ILO, IHI, SCALE, ABNRM,
!                          RCONDE, RCONDV, WORK, LWORK, IWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          BALANC, JOBVL, JOBVR, SENSE
!       INTEGER            IHI, ILO, INFO, LDA, LDVL, LDVR, LWORK, N
!       DOUBLE PRECISION   ABNRM
!       ..
!       .. Array Arguments ..
!       INTEGER            IWORK( * )
!       DOUBLE PRECISION   A( LDA, * ), RCONDE( * ), RCONDV( * ),
!      $                   SCALE( * ), VL( LDVL, * ), VR( LDVR, * ),
!      $                   WI( * ), WORK( * ), WR( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DGEEVX computes for an N-by-N real nonsymmetric matrix A, the
!> eigenvalues and, optionally, the left and/or right eigenvectors.
!>
!> Optionally also, it computes a balancing transformation to improve
!> the conditioning of the eigenvalues and eigenvectors (ILO, IHI,
!> SCALE, and ABNRM), reciprocal condition numbers for the eigenvalues
!> (RCONDE), and reciprocal condition numbers for the right
!> eigenvectors (RCONDV).
!>
!> The right eigenvector v(j) of A satisfies
!>                  A * v(j) = lambda(j) * v(j)
!> where lambda(j) is its eigenvalue.
!> The left eigenvector u(j) of A satisfies
!>               u(j)**H * A = lambda(j) * u(j)**H
!> where u(j)**H denotes the conjugate-transpose of u(j).
!>
!> The computed eigenvectors are normalized to have Euclidean norm
!> equal to 1 and largest component real.
!>
!> Balancing a matrix means permuting the rows and columns to make it
!> more nearly upper triangular, and applying a diagonal similarity
!> transformation D * A * D**(-1), where D is a diagonal matrix, to
!> make its rows and columns closer in norm and the condition numbers
!> of its eigenvalues and eigenvectors smaller.  The computed
!> reciprocal condition numbers correspond to the balanced matrix.
!> Permuting rows and columns will not change the condition numbers
!> (in exact arithmetic) but diagonal scaling will.  For further
!> explanation of balancing, see section 4.10.2 of the LAPACK
!> Users' Guide.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] BALANC
!> \verbatim
!>          BALANC is CHARACTER*1
!>          Indicates how the input matrix should be diagonally scaled
!>          and/or permuted to improve the conditioning of its
!>          eigenvalues.
!>          = 'N': Do not diagonally scale or permute;
!>          = 'P': Perform permutations to make the matrix more nearly
!>                 upper triangular. Do not diagonally scale;
!>          = 'S': Diagonally scale the matrix, i.e. replace A by
!>                 D*A*D**(-1), where D is a diagonal matrix chosen
!>                 to make the rows and columns of A more equal in
!>                 norm. Do not permute;
!>          = 'B': Both diagonally scale and permute A.
!>
!>          Computed reciprocal condition numbers will be for the matrix
!>          after balancing and/or permuting. Permuting does not change
!>          condition numbers (in exact arithmetic), but balancing does.
!> \endverbatim
!>
!> \param[in] JOBVL
!> \verbatim
!>          JOBVL is CHARACTER*1
!>          = 'N': left eigenvectors of A are not computed;
!>          = 'V': left eigenvectors of A are computed.
!>          If SENSE = 'E' or 'B', JOBVL must = 'V'.
!> \endverbatim
!>
!> \param[in] JOBVR
!> \verbatim
!>          JOBVR is CHARACTER*1
!>          = 'N': right eigenvectors of A are not computed;
!>          = 'V': right eigenvectors of A are computed.
!>          If SENSE = 'E' or 'B', JOBVR must = 'V'.
!> \endverbatim
!>
!> \param[in] SENSE
!> \verbatim
!>          SENSE is CHARACTER*1
!>          Determines which reciprocal condition numbers are computed.
!>          = 'N': None are computed;
!>          = 'E': Computed for eigenvalues only;
!>          = 'V': Computed for right eigenvectors only;
!>          = 'B': Computed for eigenvalues and right eigenvectors.
!>
!>          If SENSE = 'E' or 'B', both left and right eigenvectors
!>          must also be computed (JOBVL = 'V' and JOBVR = 'V').
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix A. N >= 0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (LDA,N)
!>          On entry, the N-by-N matrix A.
!>          On exit, A has been overwritten.  If JOBVL = 'V' or
!>          JOBVR = 'V', A contains the real Schur form of the balanced
!>          version of the input matrix A.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,N).
!> \endverbatim
!>
!> \param[out] WR
!> \verbatim
!>          WR is DOUBLE PRECISION array, dimension (N)
!> \endverbatim
!>
!> \param[out] WI
!> \verbatim
!>          WI is DOUBLE PRECISION array, dimension (N)
!>          WR and WI contain the real and imaginary parts,
!>          respectively, of the computed eigenvalues.  Complex
!>          conjugate pairs of eigenvalues will appear consecutively
!>          with the eigenvalue having the positive imaginary part
!>          first.
!> \endverbatim
!>
!> \param[out] VL
!> \verbatim
!>          VL is DOUBLE PRECISION array, dimension (LDVL,N)
!>          If JOBVL = 'V', the left eigenvectors u(j) are stored one
!>          after another in the columns of VL, in the same order
!>          as their eigenvalues.
!>          If JOBVL = 'N', VL is not referenced.
!>          If the j-th eigenvalue is real, then u(j) = VL(:,j),
!>          the j-th column of VL.
!>          If the j-th and (j+1)-st eigenvalues form a complex
!>          conjugate pair, then u(j) = VL(:,j) + i*VL(:,j+1) and
!>          u(j+1) = VL(:,j) - i*VL(:,j+1).
!> \endverbatim
!>
!> \param[in] LDVL
!> \verbatim
!>          LDVL is INTEGER
!>          The leading dimension of the array VL.  LDVL >= 1; if
!>          JOBVL = 'V', LDVL >= N.
!> \endverbatim
!>
!> \param[out] VR
!> \verbatim
!>          VR is DOUBLE PRECISION array, dimension (LDVR,N)
!>          If JOBVR = 'V', the right eigenvectors v(j) are stored one
!>          after another in the columns of VR, in the same order
!>          as their eigenvalues.
!>          If JOBVR = 'N', VR is not referenced.
!>          If the j-th eigenvalue is real, then v(j) = VR(:,j),
!>          the j-th column of VR.
!>          If the j-th and (j+1)-st eigenvalues form a complex
!>          conjugate pair, then v(j) = VR(:,j) + i*VR(:,j+1) and
!>          v(j+1) = VR(:,j) - i*VR(:,j+1).
!> \endverbatim
!>
!> \param[in] LDVR
!> \verbatim
!>          LDVR is INTEGER
!>          The leading dimension of the array VR.  LDVR >= 1, and if
!>          JOBVR = 'V', LDVR >= N.
!> \endverbatim
!>
!> \param[out] ILO
!> \verbatim
!>          ILO is INTEGER
!> \endverbatim
!>
!> \param[out] IHI
!> \verbatim
!>          IHI is INTEGER
!>          ILO and IHI are integer values determined when A was
!>          balanced.  The balanced A(i,j) = 0 if I > J and
!>          J = 1,...,ILO-1 or I = IHI+1,...,N.
!> \endverbatim
!>
!> \param[out] SCALE
!> \verbatim
!>          SCALE is DOUBLE PRECISION array, dimension (N)
!>          Details of the permutations and scaling factors applied
!>          when balancing A.  If P(j) is the index of the row and column
!>          interchanged with row and column j, and D(j) is the scaling
!>          factor applied to row and column j, then
!>          SCALE(J) = P(J),    for J = 1,...,ILO-1
!>                   = D(J),    for J = ILO,...,IHI
!>                   = P(J)     for J = IHI+1,...,N.
!>          The order in which the interchanges are made is N to IHI+1,
!>          then 1 to ILO-1.
!> \endverbatim
!>
!> \param[out] ABNRM
!> \verbatim
!>          ABNRM is DOUBLE PRECISION
!>          The one-norm of the balanced matrix (the maximum
!>          of the sum of absolute values of elements of any column).
!> \endverbatim
!>
!> \param[out] RCONDE
!> \verbatim
!>          RCONDE is DOUBLE PRECISION array, dimension (N)
!>          RCONDE(j) is the reciprocal condition number of the j-th
!>          eigenvalue.
!> \endverbatim
!>
!> \param[out] RCONDV
!> \verbatim
!>          RCONDV is DOUBLE PRECISION array, dimension (N)
!>          RCONDV(j) is the reciprocal condition number of the j-th
!>          right eigenvector.
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
!>          The dimension of the array WORK.   If SENSE = 'N' or 'E',
!>          LWORK >= max(1,2*N), and if JOBVL = 'V' or JOBVR = 'V',
!>          LWORK >= 3*N.  If SENSE = 'V' or 'B', LWORK >= N*(N+6).
!>          For good performance, LWORK must generally be larger.
!>
!>          If LWORK = -1, then a workspace query is assumed; the routine
!>          only calculates the optimal size of the WORK array, returns
!>          this value as the first entry of the WORK array, and no error
!>          message related to LWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] IWORK
!> \verbatim
!>          IWORK is INTEGER array, dimension (2*N-2)
!>          If SENSE = 'N' or 'E', not referenced.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value.
!>          > 0:  if INFO = i, the QR algorithm failed to compute all the
!>                eigenvalues, and no eigenvectors or condition numbers
!>                have been computed; elements 1:ILO-1 and i+1:N of WR
!>                and WI contain eigenvalues which have converged.
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
!  @precisions fortran d -> s
!
!> \ingroup doubleGEeigen
!
!  =====================================================================
      SUBROUTINE DGEEVX(Balanc,Jobvl,Jobvr,Sense,N,A,Lda,Wr,Wi,Vl,Ldvl, &
     &                  Vr,Ldvr,Ilo,Ihi,Scale,Abnrm,Rconde,Rcondv,Work, &
     &                  Lwork,Iwork,Info)
      USE F77KINDS                        
      USE S_DGEBAK
      USE S_DGEBAL
      USE S_DGEHRD
      USE S_DHSEQR
      USE S_DLABAD
      USE S_DLACPY
      USE S_DLAMCH
      USE S_DLANGE
      USE S_DLAPY2
      USE S_DLARTG
      USE S_DLASCL
      USE S_DNRM2
      USE S_DORGHR
      USE S_DROT
      USE S_DSCAL
      USE S_DTREVC3
      USE S_DTRSNA
      USE S_IDAMAX
      USE S_ILAENV
      USE S_LSAME
      USE S_XERBLA
      IMPLICIT NONE
!*--DGEEVX331
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D0 , ONE = 1.0D0
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Balanc
      CHARACTER :: Jobvl
      CHARACTER :: Jobvr
      CHARACTER :: Sense
      INTEGER :: N
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(*) :: Wr
      REAL(R8KIND) , DIMENSION(*) :: Wi
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldvl,*) :: Vl
      INTEGER :: Ldvl
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldvr,*) :: Vr
      INTEGER :: Ldvr
      INTEGER :: Ilo
      INTEGER :: Ihi
      REAL(R8KIND) , DIMENSION(*) :: Scale
      REAL(R8KIND) , INTENT(INOUT) :: Abnrm
      REAL(R8KIND) , DIMENSION(*) :: Rconde
      REAL(R8KIND) , DIMENSION(*) :: Rcondv
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      REAL(R8KIND) :: anrm , bignum , cs , cscale , eps , r , scl ,     &
     &                smlnum , sn
      REAL(R8KIND) , DIMENSION(1) :: dum
      INTEGER :: hswork , i , icond , ierr , itau , iwrk , k ,          &
     &           lwork_trevc , maxwrk , minwrk , nout
      CHARACTER :: job , side
      LOGICAL :: lquery , scalea , wantvl , wantvr , wntsnb , wntsne ,  &
     &           wntsnn , wntsnv
      LOGICAL , DIMENSION(1) :: select
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
!     Test the input arguments
!
      Info = 0
      lquery = (Lwork==-1)
      wantvl = LSAME(Jobvl,'V')
      wantvr = LSAME(Jobvr,'V')
      wntsnn = LSAME(Sense,'N')
      wntsne = LSAME(Sense,'E')
      wntsnv = LSAME(Sense,'V')
      wntsnb = LSAME(Sense,'B')
      IF ( .NOT.(LSAME(Balanc,'N') .OR. LSAME(Balanc,'S') .OR.          &
     &     LSAME(Balanc,'P') .OR. LSAME(Balanc,'B')) ) THEN
         Info = -1
      ELSEIF ( (.NOT.wantvl) .AND. (.NOT.LSAME(Jobvl,'N')) ) THEN
         Info = -2
      ELSEIF ( (.NOT.wantvr) .AND. (.NOT.LSAME(Jobvr,'N')) ) THEN
         Info = -3
      ELSEIF ( .NOT.(wntsnn .OR. wntsne .OR. wntsnb .OR. wntsnv) .OR.   &
     &         ((wntsne .OR. wntsnb) .AND. .NOT.(wantvl .AND. wantvr)) )&
     &         THEN
         Info = -4
      ELSEIF ( N<0 ) THEN
         Info = -5
      ELSEIF ( Lda<MAX(1,N) ) THEN
         Info = -7
      ELSEIF ( Ldvl<1 .OR. (wantvl .AND. Ldvl<N) ) THEN
         Info = -11
      ELSEIF ( Ldvr<1 .OR. (wantvr .AND. Ldvr<N) ) THEN
         Info = -13
      ENDIF
!
!     Compute workspace
!      (Note: Comments in the code beginning "Workspace:" describe the
!       minimal amount of workspace needed at that point in the code,
!       as well as the preferred amount for good performance.
!       NB refers to the optimal block size for the immediately
!       following subroutine, as returned by ILAENV.
!       HSWORK refers to the workspace preferred by DHSEQR, as
!       calculated below. HSWORK is computed assuming ILO=1 and IHI=N,
!       the worst case.)
!
      IF ( Info==0 ) THEN
         IF ( N==0 ) THEN
            minwrk = 1
            maxwrk = 1
         ELSE
            maxwrk = N + N*ILAENV(1,'DGEHRD',' ',N,1,N,0)
!
            IF ( wantvl ) THEN
               CALL DTREVC3('L','B',select,N,A,Lda,Vl,Ldvl,Vr,Ldvr,N,   &
     &                      nout,Work,-1,ierr)
               lwork_trevc = INT(Work(1))
               maxwrk = MAX(maxwrk,N+lwork_trevc)
               CALL DHSEQR('S','V',N,1,N,A,Lda,Wr,Wi,Vl,Ldvl,Work,-1,   &
     &                     Info)
            ELSEIF ( wantvr ) THEN
               CALL DTREVC3('R','B',select,N,A,Lda,Vl,Ldvl,Vr,Ldvr,N,   &
     &                      nout,Work,-1,ierr)
               lwork_trevc = INT(Work(1))
               maxwrk = MAX(maxwrk,N+lwork_trevc)
               CALL DHSEQR('S','V',N,1,N,A,Lda,Wr,Wi,Vr,Ldvr,Work,-1,   &
     &                     Info)
            ELSEIF ( wntsnn ) THEN
               CALL DHSEQR('E','N',N,1,N,A,Lda,Wr,Wi,Vr,Ldvr,Work,-1,   &
     &                     Info)
            ELSE
               CALL DHSEQR('S','N',N,1,N,A,Lda,Wr,Wi,Vr,Ldvr,Work,-1,   &
     &                     Info)
            ENDIF
            hswork = INT(Work(1))
!
            IF ( (.NOT.wantvl) .AND. (.NOT.wantvr) ) THEN
               minwrk = 2*N
               IF ( .NOT.wntsnn ) minwrk = MAX(minwrk,N*N+6*N)
               maxwrk = MAX(maxwrk,hswork)
               IF ( .NOT.wntsnn ) maxwrk = MAX(maxwrk,N*N+6*N)
            ELSE
               minwrk = 3*N
               IF ( (.NOT.wntsnn) .AND. (.NOT.wntsne) )                 &
     &              minwrk = MAX(minwrk,N*N+6*N)
               maxwrk = MAX(maxwrk,hswork)
               maxwrk = MAX(maxwrk,N+(N-1)                              &
     &                  *ILAENV(1,'DORGHR',' ',N,1,N,-1))
               IF ( (.NOT.wntsnn) .AND. (.NOT.wntsne) )                 &
     &              maxwrk = MAX(maxwrk,N*N+6*N)
               maxwrk = MAX(maxwrk,3*N)
            ENDIF
            maxwrk = MAX(maxwrk,minwrk)
         ENDIF
         Work(1) = maxwrk
!
         IF ( Lwork<minwrk .AND. .NOT.lquery ) Info = -21
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('DGEEVX',-Info)
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
      smlnum = DLAMCH('S')
      bignum = ONE/smlnum
      CALL DLABAD(smlnum,bignum)
      smlnum = SQRT(smlnum)/eps
      bignum = ONE/smlnum
!
!     Scale A if max element outside range [SMLNUM,BIGNUM]
!
      icond = 0
      anrm = DLANGE('M',N,N,A,Lda,dum)
      scalea = .FALSE.
      IF ( anrm>ZERO .AND. anrm<smlnum ) THEN
         scalea = .TRUE.
         cscale = smlnum
      ELSEIF ( anrm>bignum ) THEN
         scalea = .TRUE.
         cscale = bignum
      ENDIF
      IF ( scalea ) CALL DLASCL('G',0,0,anrm,cscale,N,N,A,Lda,ierr)
!
!     Balance the matrix and compute ABNRM
!
      CALL DGEBAL(Balanc,N,A,Lda,Ilo,Ihi,Scale,ierr)
      Abnrm = DLANGE('1',N,N,A,Lda,dum)
      IF ( scalea ) THEN
         dum(1) = Abnrm
         CALL DLASCL('G',0,0,cscale,anrm,1,1,dum,1,ierr)
         Abnrm = dum(1)
      ENDIF
!
!     Reduce to upper Hessenberg form
!     (Workspace: need 2*N, prefer N+N*NB)
!
      itau = 1
      iwrk = itau + N
      CALL DGEHRD(N,Ilo,Ihi,A,Lda,Work(itau),Work(iwrk),Lwork-iwrk+1,   &
     &            ierr)
!
      IF ( wantvl ) THEN
!
!        Want left eigenvectors
!        Copy Householder vectors to VL
!
         side = 'L'
         CALL DLACPY('L',N,N,A,Lda,Vl,Ldvl)
!
!        Generate orthogonal matrix in VL
!        (Workspace: need 2*N-1, prefer N+(N-1)*NB)
!
         CALL DORGHR(N,Ilo,Ihi,Vl,Ldvl,Work(itau),Work(iwrk),           &
     &               Lwork-iwrk+1,ierr)
!
!        Perform QR iteration, accumulating Schur vectors in VL
!        (Workspace: need 1, prefer HSWORK (see comments) )
!
         iwrk = itau
         CALL DHSEQR('S','V',N,Ilo,Ihi,A,Lda,Wr,Wi,Vl,Ldvl,Work(iwrk),  &
     &               Lwork-iwrk+1,Info)
!
         IF ( wantvr ) THEN
!
!           Want left and right eigenvectors
!           Copy Schur vectors to VR
!
            side = 'B'
            CALL DLACPY('F',N,N,Vl,Ldvl,Vr,Ldvr)
         ENDIF
!
      ELSEIF ( wantvr ) THEN
!
!        Want right eigenvectors
!        Copy Householder vectors to VR
!
         side = 'R'
         CALL DLACPY('L',N,N,A,Lda,Vr,Ldvr)
!
!        Generate orthogonal matrix in VR
!        (Workspace: need 2*N-1, prefer N+(N-1)*NB)
!
         CALL DORGHR(N,Ilo,Ihi,Vr,Ldvr,Work(itau),Work(iwrk),           &
     &               Lwork-iwrk+1,ierr)
!
!        Perform QR iteration, accumulating Schur vectors in VR
!        (Workspace: need 1, prefer HSWORK (see comments) )
!
         iwrk = itau
         CALL DHSEQR('S','V',N,Ilo,Ihi,A,Lda,Wr,Wi,Vr,Ldvr,Work(iwrk),  &
     &               Lwork-iwrk+1,Info)
!
      ELSE
!
!        Compute eigenvalues only
!        If condition numbers desired, compute Schur form
!
         IF ( wntsnn ) THEN
            job = 'E'
         ELSE
            job = 'S'
         ENDIF
!
!        (Workspace: need 1, prefer HSWORK (see comments) )
!
         iwrk = itau
         CALL DHSEQR(job,'N',N,Ilo,Ihi,A,Lda,Wr,Wi,Vr,Ldvr,Work(iwrk),  &
     &               Lwork-iwrk+1,Info)
      ENDIF
!
!     If INFO .NE. 0 from DHSEQR, then quit
!
      IF ( Info==0 ) THEN
!
!
!        Compute left and/or right eigenvectors
!        (Workspace: need 3*N, prefer N + 2*N*NB)
!
         IF ( wantvl .OR. wantvr ) CALL DTREVC3(side,'B',select,N,A,Lda,&
     &        Vl,Ldvl,Vr,Ldvr,N,nout,Work(iwrk),Lwork-iwrk+1,ierr)
!
!     Compute condition numbers if desired
!     (Workspace: need N*N+6*N unless SENSE = 'E')
!
         IF ( .NOT.wntsnn ) CALL DTRSNA(Sense,'A',select,N,A,Lda,Vl,    &
     &                                  Ldvl,Vr,Ldvr,Rconde,Rcondv,N,   &
     &                                  nout,Work(iwrk),N,Iwork,icond)
!
         IF ( wantvl ) THEN
!
!        Undo balancing of left eigenvectors
!
            CALL DGEBAK(Balanc,'L',N,Ilo,Ihi,Scale,N,Vl,Ldvl,ierr)
!
!        Normalize left eigenvectors and make largest component real
!
            DO i = 1 , N
               IF ( Wi(i)==ZERO ) THEN
                  scl = ONE/DNRM2(N,Vl(1,i),1)
                  CALL DSCAL(N,scl,Vl(1,i),1)
               ELSEIF ( Wi(i)>ZERO ) THEN
                  scl = ONE/DLAPY2(DNRM2(N,Vl(1,i),1),                  &
     &                  DNRM2(N,Vl(1,i+1),1))
                  CALL DSCAL(N,scl,Vl(1,i),1)
                  CALL DSCAL(N,scl,Vl(1,i+1),1)
                  DO k = 1 , N
                     Work(k) = Vl(k,i)**2 + Vl(k,i+1)**2
                  ENDDO
                  k = IDAMAX(N,Work,1)
                  CALL DLARTG(Vl(k,i),Vl(k,i+1),cs,sn,r)
                  CALL DROT(N,Vl(1,i),1,Vl(1,i+1),1,cs,sn)
                  Vl(k,i+1) = ZERO
               ENDIF
            ENDDO
         ENDIF
!
         IF ( wantvr ) THEN
!
!        Undo balancing of right eigenvectors
!
            CALL DGEBAK(Balanc,'R',N,Ilo,Ihi,Scale,N,Vr,Ldvr,ierr)
!
!        Normalize right eigenvectors and make largest component real
!
            DO i = 1 , N
               IF ( Wi(i)==ZERO ) THEN
                  scl = ONE/DNRM2(N,Vr(1,i),1)
                  CALL DSCAL(N,scl,Vr(1,i),1)
               ELSEIF ( Wi(i)>ZERO ) THEN
                  scl = ONE/DLAPY2(DNRM2(N,Vr(1,i),1),                  &
     &                  DNRM2(N,Vr(1,i+1),1))
                  CALL DSCAL(N,scl,Vr(1,i),1)
                  CALL DSCAL(N,scl,Vr(1,i+1),1)
                  DO k = 1 , N
                     Work(k) = Vr(k,i)**2 + Vr(k,i+1)**2
                  ENDDO
                  k = IDAMAX(N,Work,1)
                  CALL DLARTG(Vr(k,i),Vr(k,i+1),cs,sn,r)
                  CALL DROT(N,Vr(1,i),1,Vr(1,i+1),1,cs,sn)
                  Vr(k,i+1) = ZERO
               ENDIF
            ENDDO
         ENDIF
      ENDIF
!
!     Undo scaling if necessary
!
      IF ( scalea ) THEN
         CALL DLASCL('G',0,0,cscale,anrm,N-Info,1,Wr(Info+1),           &
     &               MAX(N-Info,1),ierr)
         CALL DLASCL('G',0,0,cscale,anrm,N-Info,1,Wi(Info+1),           &
     &               MAX(N-Info,1),ierr)
         IF ( Info==0 ) THEN
            IF ( (wntsnv .OR. wntsnb) .AND. icond==0 )                  &
     &           CALL DLASCL('G',0,0,cscale,anrm,N,1,Rcondv,N,ierr)
         ELSE
            CALL DLASCL('G',0,0,cscale,anrm,Ilo-1,1,Wr,N,ierr)
            CALL DLASCL('G',0,0,cscale,anrm,Ilo-1,1,Wi,N,ierr)
         ENDIF
      ENDIF
!
      Work(1) = maxwrk
!
!     End of DGEEVX
!
      END SUBROUTINE DGEEVX
