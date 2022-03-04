!*==zggevx.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief <b> ZGGEVX computes the eigenvalues and, optionally, the left and/or right eigenvectors for GE matrices</b>
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZGGEVX + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zggevx.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zggevx.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zggevx.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZGGEVX( BALANC, JOBVL, JOBVR, SENSE, N, A, LDA, B, LDB,
!                          ALPHA, BETA, VL, LDVL, VR, LDVR, ILO, IHI,
!                          LSCALE, RSCALE, ABNRM, BBNRM, RCONDE, RCONDV,
!                          WORK, LWORK, RWORK, IWORK, BWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          BALANC, JOBVL, JOBVR, SENSE
!       INTEGER            IHI, ILO, INFO, LDA, LDB, LDVL, LDVR, LWORK, N
!       DOUBLE PRECISION   ABNRM, BBNRM
!       ..
!       .. Array Arguments ..
!       LOGICAL            BWORK( * )
!       INTEGER            IWORK( * )
!       DOUBLE PRECISION   LSCALE( * ), RCONDE( * ), RCONDV( * ),
!      $                   RSCALE( * ), RWORK( * )
!       COMPLEX*16         A( LDA, * ), ALPHA( * ), B( LDB, * ),
!      $                   BETA( * ), VL( LDVL, * ), VR( LDVR, * ),
!      $                   WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZGGEVX computes for a pair of N-by-N complex nonsymmetric matrices
!> (A,B) the generalized eigenvalues, and optionally, the left and/or
!> right generalized eigenvectors.
!>
!> Optionally, it also computes a balancing transformation to improve
!> the conditioning of the eigenvalues and eigenvectors (ILO, IHI,
!> LSCALE, RSCALE, ABNRM, and BBNRM), reciprocal condition numbers for
!> the eigenvalues (RCONDE), and reciprocal condition numbers for the
!> right eigenvectors (RCONDV).
!>
!> A generalized eigenvalue for a pair of matrices (A,B) is a scalar
!> lambda or a ratio alpha/beta = lambda, such that A - lambda*B is
!> singular. It is usually represented as the pair (alpha,beta), as
!> there is a reasonable interpretation for beta=0, and even for both
!> being zero.
!>
!> The right eigenvector v(j) corresponding to the eigenvalue lambda(j)
!> of (A,B) satisfies
!>                  A * v(j) = lambda(j) * B * v(j) .
!> The left eigenvector u(j) corresponding to the eigenvalue lambda(j)
!> of (A,B) satisfies
!>                  u(j)**H * A  = lambda(j) * u(j)**H * B.
!> where u(j)**H is the conjugate-transpose of u(j).
!>
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] BALANC
!> \verbatim
!>          BALANC is CHARACTER*1
!>          Specifies the balance option to be performed:
!>          = 'N':  do not diagonally scale or permute;
!>          = 'P':  permute only;
!>          = 'S':  scale only;
!>          = 'B':  both permute and scale.
!>          Computed reciprocal condition numbers will be for the
!>          matrices after permuting and/or balancing. Permuting does
!>          not change condition numbers (in exact arithmetic), but
!>          balancing does.
!> \endverbatim
!>
!> \param[in] JOBVL
!> \verbatim
!>          JOBVL is CHARACTER*1
!>          = 'N':  do not compute the left generalized eigenvectors;
!>          = 'V':  compute the left generalized eigenvectors.
!> \endverbatim
!>
!> \param[in] JOBVR
!> \verbatim
!>          JOBVR is CHARACTER*1
!>          = 'N':  do not compute the right generalized eigenvectors;
!>          = 'V':  compute the right generalized eigenvectors.
!> \endverbatim
!>
!> \param[in] SENSE
!> \verbatim
!>          SENSE is CHARACTER*1
!>          Determines which reciprocal condition numbers are computed.
!>          = 'N': none are computed;
!>          = 'E': computed for eigenvalues only;
!>          = 'V': computed for eigenvectors only;
!>          = 'B': computed for eigenvalues and eigenvectors.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrices A, B, VL, and VR.  N >= 0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension (LDA, N)
!>          On entry, the matrix A in the pair (A,B).
!>          On exit, A has been overwritten. If JOBVL='V' or JOBVR='V'
!>          or both, then A contains the first part of the complex Schur
!>          form of the "balanced" versions of the input A and B.
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
!>          B is COMPLEX*16 array, dimension (LDB, N)
!>          On entry, the matrix B in the pair (A,B).
!>          On exit, B has been overwritten. If JOBVL='V' or JOBVR='V'
!>          or both, then B contains the second part of the complex
!>          Schur form of the "balanced" versions of the input A and B.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of B.  LDB >= max(1,N).
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
!>          On exit, ALPHA(j)/BETA(j), j=1,...,N, will be the generalized
!>          eigenvalues.
!>
!>          Note: the quotient ALPHA(j)/BETA(j) ) may easily over- or
!>          underflow, and BETA(j) may even be zero.  Thus, the user
!>          should avoid naively computing the ratio ALPHA/BETA.
!>          However, ALPHA will be always less than and usually
!>          comparable with norm(A) in magnitude, and BETA always less
!>          than and usually comparable with norm(B).
!> \endverbatim
!>
!> \param[out] VL
!> \verbatim
!>          VL is COMPLEX*16 array, dimension (LDVL,N)
!>          If JOBVL = 'V', the left generalized eigenvectors u(j) are
!>          stored one after another in the columns of VL, in the same
!>          order as their eigenvalues.
!>          Each eigenvector will be scaled so the largest component
!>          will have abs(real part) + abs(imag. part) = 1.
!>          Not referenced if JOBVL = 'N'.
!> \endverbatim
!>
!> \param[in] LDVL
!> \verbatim
!>          LDVL is INTEGER
!>          The leading dimension of the matrix VL. LDVL >= 1, and
!>          if JOBVL = 'V', LDVL >= N.
!> \endverbatim
!>
!> \param[out] VR
!> \verbatim
!>          VR is COMPLEX*16 array, dimension (LDVR,N)
!>          If JOBVR = 'V', the right generalized eigenvectors v(j) are
!>          stored one after another in the columns of VR, in the same
!>          order as their eigenvalues.
!>          Each eigenvector will be scaled so the largest component
!>          will have abs(real part) + abs(imag. part) = 1.
!>          Not referenced if JOBVR = 'N'.
!> \endverbatim
!>
!> \param[in] LDVR
!> \verbatim
!>          LDVR is INTEGER
!>          The leading dimension of the matrix VR. LDVR >= 1, and
!>          if JOBVR = 'V', LDVR >= N.
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
!>          ILO and IHI are integer values such that on exit
!>          A(i,j) = 0 and B(i,j) = 0 if i > j and
!>          j = 1,...,ILO-1 or i = IHI+1,...,N.
!>          If BALANC = 'N' or 'S', ILO = 1 and IHI = N.
!> \endverbatim
!>
!> \param[out] LSCALE
!> \verbatim
!>          LSCALE is DOUBLE PRECISION array, dimension (N)
!>          Details of the permutations and scaling factors applied
!>          to the left side of A and B.  If PL(j) is the index of the
!>          row interchanged with row j, and DL(j) is the scaling
!>          factor applied to row j, then
!>            LSCALE(j) = PL(j)  for j = 1,...,ILO-1
!>                      = DL(j)  for j = ILO,...,IHI
!>                      = PL(j)  for j = IHI+1,...,N.
!>          The order in which the interchanges are made is N to IHI+1,
!>          then 1 to ILO-1.
!> \endverbatim
!>
!> \param[out] RSCALE
!> \verbatim
!>          RSCALE is DOUBLE PRECISION array, dimension (N)
!>          Details of the permutations and scaling factors applied
!>          to the right side of A and B.  If PR(j) is the index of the
!>          column interchanged with column j, and DR(j) is the scaling
!>          factor applied to column j, then
!>            RSCALE(j) = PR(j)  for j = 1,...,ILO-1
!>                      = DR(j)  for j = ILO,...,IHI
!>                      = PR(j)  for j = IHI+1,...,N
!>          The order in which the interchanges are made is N to IHI+1,
!>          then 1 to ILO-1.
!> \endverbatim
!>
!> \param[out] ABNRM
!> \verbatim
!>          ABNRM is DOUBLE PRECISION
!>          The one-norm of the balanced matrix A.
!> \endverbatim
!>
!> \param[out] BBNRM
!> \verbatim
!>          BBNRM is DOUBLE PRECISION
!>          The one-norm of the balanced matrix B.
!> \endverbatim
!>
!> \param[out] RCONDE
!> \verbatim
!>          RCONDE is DOUBLE PRECISION array, dimension (N)
!>          If SENSE = 'E' or 'B', the reciprocal condition numbers of
!>          the eigenvalues, stored in consecutive elements of the array.
!>          If SENSE = 'N' or 'V', RCONDE is not referenced.
!> \endverbatim
!>
!> \param[out] RCONDV
!> \verbatim
!>          RCONDV is DOUBLE PRECISION array, dimension (N)
!>          If JOB = 'V' or 'B', the estimated reciprocal condition
!>          numbers of the eigenvectors, stored in consecutive elements
!>          of the array. If the eigenvalues cannot be reordered to
!>          compute RCONDV(j), RCONDV(j) is set to 0; this can only occur
!>          when the true value would be very small anyway.
!>          If SENSE = 'N' or 'E', RCONDV is not referenced.
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
!>          The dimension of the array WORK. LWORK >= max(1,2*N).
!>          If SENSE = 'E', LWORK >= max(1,4*N).
!>          If SENSE = 'V' or 'B', LWORK >= max(1,2*N*N+2*N).
!>
!>          If LWORK = -1, then a workspace query is assumed; the routine
!>          only calculates the optimal size of the WORK array, returns
!>          this value as the first entry of the WORK array, and no error
!>          message related to LWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is DOUBLE PRECISION array, dimension (lrwork)
!>          lrwork must be at least max(1,6*N) if BALANC = 'S' or 'B',
!>          and at least max(1,2*N) otherwise.
!>          Real workspace.
!> \endverbatim
!>
!> \param[out] IWORK
!> \verbatim
!>          IWORK is INTEGER array, dimension (N+2)
!>          If SENSE = 'E', IWORK is not referenced.
!> \endverbatim
!>
!> \param[out] BWORK
!> \verbatim
!>          BWORK is LOGICAL array, dimension (N)
!>          If SENSE = 'N', BWORK is not referenced.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value.
!>          = 1,...,N:
!>                The QZ iteration failed.  No eigenvectors have been
!>                calculated, but ALPHA(j) and BETA(j) should be correct
!>                for j=INFO+1,...,N.
!>          > N:  =N+1: other than QZ iteration failed in ZHGEQZ.
!>                =N+2: error return from ZTGEVC.
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
!> \date April 2012
!
!> \ingroup complex16GEeigen
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  Balancing a matrix pair (A,B) includes, first, permuting rows and
!>  columns to isolate eigenvalues, second, applying diagonal similarity
!>  transformation to the rows and columns to make the rows and columns
!>  as close in norm as possible. The computed reciprocal condition
!>  numbers correspond to the balanced matrix. Permuting rows and columns
!>  will not change the condition numbers (in exact arithmetic) but
!>  diagonal scaling will.  For further explanation of balancing, see
!>  section 4.11.1.2 of LAPACK Users' Guide.
!>
!>  An approximate error bound on the chordal distance between the i-th
!>  computed generalized eigenvalue w and the corresponding exact
!>  eigenvalue lambda is
!>
!>       chord(w, lambda) <= EPS * norm(ABNRM, BBNRM) / RCONDE(I)
!>
!>  An approximate error bound for the angle between the i-th computed
!>  eigenvector VL(i) or VR(i) is given by
!>
!>       EPS * norm(ABNRM, BBNRM) / DIF(i).
!>
!>  For further explanation of the reciprocal condition numbers RCONDE
!>  and RCONDV, see section 4.11 of LAPACK User's Guide.
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE ZGGEVX(Balanc,Jobvl,Jobvr,Sense,N,A,Lda,B,Ldb,Alpha,   &
     &                  Beta,Vl,Ldvl,Vr,Ldvr,Ilo,Ihi,Lscale,Rscale,     &
     &                  Abnrm,Bbnrm,Rconde,Rcondv,Work,Lwork,Rwork,     &
     &                  Iwork,Bwork,Info)
      IMPLICIT NONE
!*--ZGGEVX378
!
!  -- LAPACK driver routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     April 2012
!
!     .. Scalar Arguments ..
      CHARACTER Balanc , Jobvl , Jobvr , Sense
      INTEGER Ihi , Ilo , Info , Lda , Ldb , Ldvl , Ldvr , Lwork , N
      DOUBLE PRECISION Abnrm , Bbnrm
!     ..
!     .. Array Arguments ..
      LOGICAL Bwork(*)
      INTEGER Iwork(*)
      DOUBLE PRECISION Lscale(*) , Rconde(*) , Rcondv(*) , Rscale(*) ,  &
     &                 Rwork(*)
      COMPLEX*16 A(Lda,*) , Alpha(*) , B(Ldb,*) , Beta(*) , Vl(Ldvl,*) ,&
     &           Vr(Ldvr,*) , Work(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO , ONE
      PARAMETER (ZERO=0.0D+0,ONE=1.0D+0)
      COMPLEX*16 CZERO , CONE
      PARAMETER (CZERO=(0.0D+0,0.0D+0),CONE=(1.0D+0,0.0D+0))
!     ..
!     .. Local Scalars ..
      LOGICAL ilascl , ilbscl , ilv , ilvl , ilvr , lquery , noscl ,    &
     &        wantsb , wantse , wantsn , wantsv
      CHARACTER chtemp
      INTEGER i , icols , ierr , ijobvl , ijobvr , in , irows , itau ,  &
     &        iwrk , iwrk1 , j , jc , jr , m , maxwrk , minwrk
      DOUBLE PRECISION anrm , anrmto , bignum , bnrm , bnrmto , eps ,   &
     &                 smlnum , temp
      COMPLEX*16 x
!     ..
!     .. Local Arrays ..
      LOGICAL ldumma(1)
!     ..
!     .. External Subroutines ..
      EXTERNAL DLABAD , DLASCL , XERBLA , ZGEQRF , ZGGBAK , ZGGBAL ,    &
     &         ZGGHRD , ZHGEQZ , ZLACPY , ZLASCL , ZLASET , ZTGEVC ,    &
     &         ZTGSNA , ZUNGQR , ZUNMQR
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      INTEGER ILAENV
      DOUBLE PRECISION DLAMCH , ZLANGE
      EXTERNAL LSAME , ILAENV , DLAMCH , ZLANGE
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , DBLE , DIMAG , MAX , SQRT
!     ..
!     .. Statement Functions ..
      DOUBLE PRECISION ABS1
!     ..
!     .. Statement Function definitions ..
      ABS1(x) = ABS(DBLE(x)) + ABS(DIMAG(x))
!     ..
!     .. Executable Statements ..
!
!     Decode the input arguments
!
      IF ( LSAME(Jobvl,'N') ) THEN
         ijobvl = 1
         ilvl = .FALSE.
      ELSEIF ( LSAME(Jobvl,'V') ) THEN
         ijobvl = 2
         ilvl = .TRUE.
      ELSE
         ijobvl = -1
         ilvl = .FALSE.
      ENDIF
!
      IF ( LSAME(Jobvr,'N') ) THEN
         ijobvr = 1
         ilvr = .FALSE.
      ELSEIF ( LSAME(Jobvr,'V') ) THEN
         ijobvr = 2
         ilvr = .TRUE.
      ELSE
         ijobvr = -1
         ilvr = .FALSE.
      ENDIF
      ilv = ilvl .OR. ilvr
!
      noscl = LSAME(Balanc,'N') .OR. LSAME(Balanc,'P')
      wantsn = LSAME(Sense,'N')
      wantse = LSAME(Sense,'E')
      wantsv = LSAME(Sense,'V')
      wantsb = LSAME(Sense,'B')
!
!     Test the input arguments
!
      Info = 0
      lquery = (Lwork==-1)
      IF ( .NOT.(noscl .OR. LSAME(Balanc,'S') .OR. LSAME(Balanc,'B')) ) &
     &     THEN
         Info = -1
      ELSEIF ( ijobvl<=0 ) THEN
         Info = -2
      ELSEIF ( ijobvr<=0 ) THEN
         Info = -3
      ELSEIF ( .NOT.(wantsn .OR. wantse .OR. wantsb .OR. wantsv) ) THEN
         Info = -4
      ELSEIF ( N<0 ) THEN
         Info = -5
      ELSEIF ( Lda<MAX(1,N) ) THEN
         Info = -7
      ELSEIF ( Ldb<MAX(1,N) ) THEN
         Info = -9
      ELSEIF ( Ldvl<1 .OR. (ilvl .AND. Ldvl<N) ) THEN
         Info = -13
      ELSEIF ( Ldvr<1 .OR. (ilvr .AND. Ldvr<N) ) THEN
         Info = -15
      ENDIF
!
!     Compute workspace
!      (Note: Comments in the code beginning "Workspace:" describe the
!       minimal amount of workspace needed at that point in the code,
!       as well as the preferred amount for good performance.
!       NB refers to the optimal block size for the immediately
!       following subroutine, as returned by ILAENV. The workspace is
!       computed assuming ILO = 1 and IHI = N, the worst case.)
!
      IF ( Info==0 ) THEN
         IF ( N==0 ) THEN
            minwrk = 1
            maxwrk = 1
         ELSE
            minwrk = 2*N
            IF ( wantse ) THEN
               minwrk = 4*N
            ELSEIF ( wantsv .OR. wantsb ) THEN
               minwrk = 2*N*(N+1)
            ENDIF
            maxwrk = minwrk
            maxwrk = MAX(maxwrk,N+N*ILAENV(1,'ZGEQRF',' ',N,1,N,0))
            maxwrk = MAX(maxwrk,N+N*ILAENV(1,'ZUNMQR',' ',N,1,N,0))
            IF ( ilvl ) maxwrk = MAX(maxwrk,N+N*ILAENV(1,'ZUNGQR',' ',N,&
     &                           1,N,0))
         ENDIF
         Work(1) = maxwrk
!
         IF ( Lwork<minwrk .AND. .NOT.lquery ) Info = -25
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('ZGGEVX',-Info)
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
      anrm = ZLANGE('M',N,N,A,Lda,Rwork)
      ilascl = .FALSE.
      IF ( anrm>ZERO .AND. anrm<smlnum ) THEN
         anrmto = smlnum
         ilascl = .TRUE.
      ELSEIF ( anrm>bignum ) THEN
         anrmto = bignum
         ilascl = .TRUE.
      ENDIF
      IF ( ilascl ) CALL ZLASCL('G',0,0,anrm,anrmto,N,N,A,Lda,ierr)
!
!     Scale B if max element outside range [SMLNUM,BIGNUM]
!
      bnrm = ZLANGE('M',N,N,B,Ldb,Rwork)
      ilbscl = .FALSE.
      IF ( bnrm>ZERO .AND. bnrm<smlnum ) THEN
         bnrmto = smlnum
         ilbscl = .TRUE.
      ELSEIF ( bnrm>bignum ) THEN
         bnrmto = bignum
         ilbscl = .TRUE.
      ENDIF
      IF ( ilbscl ) CALL ZLASCL('G',0,0,bnrm,bnrmto,N,N,B,Ldb,ierr)
!
!     Permute and/or balance the matrix pair (A,B)
!     (Real Workspace: need 6*N if BALANC = 'S' or 'B', 1 otherwise)
!
      CALL ZGGBAL(Balanc,N,A,Lda,B,Ldb,Ilo,Ihi,Lscale,Rscale,Rwork,ierr)
!
!     Compute ABNRM and BBNRM
!
      Abnrm = ZLANGE('1',N,N,A,Lda,Rwork(1))
      IF ( ilascl ) THEN
         Rwork(1) = Abnrm
         CALL DLASCL('G',0,0,anrmto,anrm,1,1,Rwork(1),1,ierr)
         Abnrm = Rwork(1)
      ENDIF
!
      Bbnrm = ZLANGE('1',N,N,B,Ldb,Rwork(1))
      IF ( ilbscl ) THEN
         Rwork(1) = Bbnrm
         CALL DLASCL('G',0,0,bnrmto,bnrm,1,1,Rwork(1),1,ierr)
         Bbnrm = Rwork(1)
      ENDIF
!
!     Reduce B to triangular form (QR decomposition of B)
!     (Complex Workspace: need N, prefer N*NB )
!
      irows = Ihi + 1 - Ilo
      IF ( ilv .OR. .NOT.wantsn ) THEN
         icols = N + 1 - Ilo
      ELSE
         icols = irows
      ENDIF
      itau = 1
      iwrk = itau + irows
      CALL ZGEQRF(irows,icols,B(Ilo,Ilo),Ldb,Work(itau),Work(iwrk),     &
     &            Lwork+1-iwrk,ierr)
!
!     Apply the unitary transformation to A
!     (Complex Workspace: need N, prefer N*NB)
!
      CALL ZUNMQR('L','C',irows,icols,irows,B(Ilo,Ilo),Ldb,Work(itau),  &
     &            A(Ilo,Ilo),Lda,Work(iwrk),Lwork+1-iwrk,ierr)
!
!     Initialize VL and/or VR
!     (Workspace: need N, prefer N*NB)
!
      IF ( ilvl ) THEN
         CALL ZLASET('Full',N,N,CZERO,CONE,Vl,Ldvl)
         IF ( irows>1 ) CALL ZLACPY('L',irows-1,irows-1,B(Ilo+1,Ilo),   &
     &                              Ldb,Vl(Ilo+1,Ilo),Ldvl)
         CALL ZUNGQR(irows,irows,irows,Vl(Ilo,Ilo),Ldvl,Work(itau),     &
     &               Work(iwrk),Lwork+1-iwrk,ierr)
      ENDIF
!
      IF ( ilvr ) CALL ZLASET('Full',N,N,CZERO,CONE,Vr,Ldvr)
!
!     Reduce to generalized Hessenberg form
!     (Workspace: none needed)
!
      IF ( ilv .OR. .NOT.wantsn ) THEN
!
!        Eigenvectors requested -- work on whole matrix.
!
         CALL ZGGHRD(Jobvl,Jobvr,N,Ilo,Ihi,A,Lda,B,Ldb,Vl,Ldvl,Vr,Ldvr, &
     &               ierr)
      ELSE
         CALL ZGGHRD('N','N',irows,1,irows,A(Ilo,Ilo),Lda,B(Ilo,Ilo),   &
     &               Ldb,Vl,Ldvl,Vr,Ldvr,ierr)
      ENDIF
!
!     Perform QZ algorithm (Compute eigenvalues, and optionally, the
!     Schur forms and Schur vectors)
!     (Complex Workspace: need N)
!     (Real Workspace: need N)
!
      iwrk = itau
      IF ( ilv .OR. .NOT.wantsn ) THEN
         chtemp = 'S'
      ELSE
         chtemp = 'E'
      ENDIF
!
      CALL ZHGEQZ(chtemp,Jobvl,Jobvr,N,Ilo,Ihi,A,Lda,B,Ldb,Alpha,Beta,  &
     &            Vl,Ldvl,Vr,Ldvr,Work(iwrk),Lwork+1-iwrk,Rwork,ierr)
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
!     Compute Eigenvectors and estimate condition numbers if desired
!     ZTGEVC: (Complex Workspace: need 2*N )
!             (Real Workspace:    need 2*N )
!     ZTGSNA: (Complex Workspace: need 2*N*N if SENSE='V' or 'B')
!             (Integer Workspace: need N+2 )
!
      IF ( ilv .OR. .NOT.wantsn ) THEN
         IF ( ilv ) THEN
            IF ( .NOT.(ilvl) ) THEN
               chtemp = 'R'
            ELSEIF ( ilvr ) THEN
               chtemp = 'B'
            ELSE
               chtemp = 'L'
            ENDIF
!
            CALL ZTGEVC(chtemp,'B',ldumma,N,A,Lda,B,Ldb,Vl,Ldvl,Vr,Ldvr,&
     &                  N,in,Work(iwrk),Rwork,ierr)
            IF ( ierr/=0 ) THEN
               Info = N + 2
               GOTO 100
            ENDIF
         ENDIF
!
         IF ( .NOT.wantsn ) THEN
!
!           compute eigenvectors (DTGEVC) and estimate condition
!           numbers (DTGSNA). Note that the definition of the condition
!           number is not invariant under transformation (u,v) to
!           (Q*u, Z*v), where (u,v) are eigenvectors of the generalized
!           Schur form (S,T), Q and Z are orthogonal matrices. In order
!           to avoid using extra 2*N*N workspace, we have to
!           re-calculate eigenvectors and estimate the condition numbers
!           one at a time.
!
            DO i = 1 , N
!
               DO j = 1 , N
                  Bwork(j) = .FALSE.
               ENDDO
               Bwork(i) = .TRUE.
!
               iwrk = N + 1
               iwrk1 = iwrk + N
!
               IF ( wantse .OR. wantsb ) THEN
                  CALL ZTGEVC('B','S',Bwork,N,A,Lda,B,Ldb,Work(1),N,    &
     &                        Work(iwrk),N,1,m,Work(iwrk1),Rwork,ierr)
                  IF ( ierr/=0 ) THEN
                     Info = N + 2
                     GOTO 100
                  ENDIF
               ENDIF
!
               CALL ZTGSNA(Sense,'S',Bwork,N,A,Lda,B,Ldb,Work(1),N,     &
     &                     Work(iwrk),N,Rconde(i),Rcondv(i),1,m,        &
     &                     Work(iwrk1),Lwork-iwrk1+1,Iwork,ierr)
!
            ENDDO
         ENDIF
      ENDIF
!
!     Undo balancing on VL and VR and normalization
!     (Workspace: none needed)
!
      IF ( ilvl ) THEN
         CALL ZGGBAK(Balanc,'L',N,Ilo,Ihi,Lscale,Rscale,N,Vl,Ldvl,ierr)
!
         DO jc = 1 , N
            temp = ZERO
            DO jr = 1 , N
               temp = MAX(temp,ABS1(Vl(jr,jc)))
            ENDDO
            IF ( temp>=smlnum ) THEN
               temp = ONE/temp
               DO jr = 1 , N
                  Vl(jr,jc) = Vl(jr,jc)*temp
               ENDDO
            ENDIF
         ENDDO
      ENDIF
!
      IF ( ilvr ) THEN
         CALL ZGGBAK(Balanc,'R',N,Ilo,Ihi,Lscale,Rscale,N,Vr,Ldvr,ierr)
         DO jc = 1 , N
            temp = ZERO
            DO jr = 1 , N
               temp = MAX(temp,ABS1(Vr(jr,jc)))
            ENDDO
            IF ( temp>=smlnum ) THEN
               temp = ONE/temp
               DO jr = 1 , N
                  Vr(jr,jc) = Vr(jr,jc)*temp
               ENDDO
            ENDIF
         ENDDO
      ENDIF
!
!     Undo scaling if necessary
!
!
 100  IF ( ilascl ) CALL ZLASCL('G',0,0,anrmto,anrm,N,1,Alpha,N,ierr)
!
      IF ( ilbscl ) CALL ZLASCL('G',0,0,bnrmto,bnrm,N,1,Beta,N,ierr)
!
      Work(1) = maxwrk
!
!     End of ZGGEVX
!
      END SUBROUTINE ZGGEVX
