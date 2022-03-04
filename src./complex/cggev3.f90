!*==cggev3.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief <b> CGGEV3 computes the eigenvalues and, optionally, the left and/or right eigenvectors for GE matrices (blocked algorithm)</b>
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CGGEV3 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cggev3.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cggev3.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cggev3.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CGGEV3( JOBVL, JOBVR, N, A, LDA, B, LDB, ALPHA, BETA,
!      $                   VL, LDVL, VR, LDVR, WORK, LWORK, RWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          JOBVL, JOBVR
!       INTEGER            INFO, LDA, LDB, LDVL, LDVR, LWORK, N
!       ..
!       .. Array Arguments ..
!       REAL               RWORK( * )
!       COMPLEX            A( LDA, * ), ALPHA( * ), B( LDB, * ),
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
!> CGGEV3 computes for a pair of N-by-N complex nonsymmetric matrices
!> (A,B), the generalized eigenvalues, and optionally, the left and/or
!> right generalized eigenvectors.
!>
!> A generalized eigenvalue for a pair of matrices (A,B) is a scalar
!> lambda or a ratio alpha/beta = lambda, such that A - lambda*B is
!> singular. It is usually represented as the pair (alpha,beta), as
!> there is a reasonable interpretation for beta=0, and even for both
!> being zero.
!>
!> The right generalized eigenvector v(j) corresponding to the
!> generalized eigenvalue lambda(j) of (A,B) satisfies
!>
!>              A * v(j) = lambda(j) * B * v(j).
!>
!> The left generalized eigenvector u(j) corresponding to the
!> generalized eigenvalues lambda(j) of (A,B) satisfies
!>
!>              u(j)**H * A = lambda(j) * u(j)**H * B
!>
!> where u(j)**H is the conjugate-transpose of u(j).
!> \endverbatim
!
!  Arguments:
!  ==========
!
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
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrices A, B, VL, and VR.  N >= 0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is COMPLEX array, dimension (LDA, N)
!>          On entry, the matrix A in the pair (A,B).
!>          On exit, A has been overwritten.
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
!>          B is COMPLEX array, dimension (LDB, N)
!>          On entry, the matrix B in the pair (A,B).
!>          On exit, B has been overwritten.
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
!>          ALPHA is COMPLEX array, dimension (N)
!> \endverbatim
!>
!> \param[out] BETA
!> \verbatim
!>          BETA is COMPLEX array, dimension (N)
!>          On exit, ALPHA(j)/BETA(j), j=1,...,N, will be the
!>          generalized eigenvalues.
!>
!>          Note: the quotients ALPHA(j)/BETA(j) may easily over- or
!>          underflow, and BETA(j) may even be zero.  Thus, the user
!>          should avoid naively computing the ratio alpha/beta.
!>          However, ALPHA will be always less than and usually
!>          comparable with norm(A) in magnitude, and BETA always less
!>          than and usually comparable with norm(B).
!> \endverbatim
!>
!> \param[out] VL
!> \verbatim
!>          VL is COMPLEX array, dimension (LDVL,N)
!>          If JOBVL = 'V', the left generalized eigenvectors u(j) are
!>          stored one after another in the columns of VL, in the same
!>          order as their eigenvalues.
!>          Each eigenvector is scaled so the largest component has
!>          abs(real part) + abs(imag. part) = 1.
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
!>          VR is COMPLEX array, dimension (LDVR,N)
!>          If JOBVR = 'V', the right generalized eigenvectors v(j) are
!>          stored one after another in the columns of VR, in the same
!>          order as their eigenvalues.
!>          Each eigenvector is scaled so the largest component has
!>          abs(real part) + abs(imag. part) = 1.
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
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX array, dimension (MAX(1,LWORK))
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
!> \param[out] RWORK
!> \verbatim
!>          RWORK is REAL array, dimension (8*N)
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value.
!>          =1,...,N:
!>                The QZ iteration failed.  No eigenvectors have been
!>                calculated, but ALPHA(j) and BETA(j) should be
!>                correct for j=INFO+1,...,N.
!>          > N:  =N+1: other then QZ iteration failed in SHGEQZ,
!>                =N+2: error return from STGEVC.
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
!> \ingroup complexGEeigen
!
!  =====================================================================
      SUBROUTINE CGGEV3(Jobvl,Jobvr,N,A,Lda,B,Ldb,Alpha,Beta,Vl,Ldvl,Vr,&
     &                  Ldvr,Work,Lwork,Rwork,Info)
      USE S_CGEQRF
      USE S_CGGBAK
      USE S_CGGBAL
      USE S_CGGHD3
      USE S_CHGEQZ
      USE S_CLACPY
      USE S_CLANGE
      USE S_CLASCL
      USE S_CLASET
      USE S_CTGEVC
      USE S_CUNGQR
      USE S_CUNMQR
      USE S_LSAME
      USE S_SLABAD
      USE S_SLAMCH
      USE S_XERBLA
      IMPLICIT NONE
!*--CGGEV3236
!
! PARAMETER definitions rewritten by SPAG
!
      REAL , PARAMETER  ::  ZERO = 0.0E0 , ONE = 1.0E0
      COMPLEX , PARAMETER  ::  CZERO = (0.0E0,0.0E0) ,                  &
     &                         CONE = (1.0E0,0.0E0)
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Jobvl
      CHARACTER :: Jobvr
      INTEGER :: N
      COMPLEX , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX , DIMENSION(*) :: Alpha
      COMPLEX , DIMENSION(*) :: Beta
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldvl,*) :: Vl
      INTEGER :: Ldvl
      COMPLEX , INTENT(INOUT) , DIMENSION(Ldvr,*) :: Vr
      INTEGER :: Ldvr
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      REAL , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      REAL :: ABS1
      REAL :: anrm , anrmto , bignum , bnrm , bnrmto , eps , smlnum ,   &
     &        temp
      CHARACTER :: chtemp
      INTEGER :: icols , ierr , ihi , ijobvl , ijobvr , ileft , ilo ,   &
     &           in , iright , irows , irwrk , itau , iwrk , jc , jr ,  &
     &           lwkopt
      LOGICAL :: ilascl , ilbscl , ilv , ilvl , ilvr , lquery
      LOGICAL , DIMENSION(1) :: ldumma
      COMPLEX :: x
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
!     .. Statement Functions ..
!     ..
!     .. Statement Function definitions ..
      ABS1(x) = ABS(REAL(x)) + ABS(AIMAG(x))
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
!     Test the input arguments
!
      Info = 0
      lquery = (Lwork==-1)
      IF ( ijobvl<=0 ) THEN
         Info = -1
      ELSEIF ( ijobvr<=0 ) THEN
         Info = -2
      ELSEIF ( N<0 ) THEN
         Info = -3
      ELSEIF ( Lda<MAX(1,N) ) THEN
         Info = -5
      ELSEIF ( Ldb<MAX(1,N) ) THEN
         Info = -7
      ELSEIF ( Ldvl<1 .OR. (ilvl .AND. Ldvl<N) ) THEN
         Info = -11
      ELSEIF ( Ldvr<1 .OR. (ilvr .AND. Ldvr<N) ) THEN
         Info = -13
      ELSEIF ( Lwork<MAX(1,2*N) .AND. .NOT.lquery ) THEN
         Info = -15
      ENDIF
!
!     Compute workspace
!
      IF ( Info==0 ) THEN
         CALL CGEQRF(N,N,B,Ldb,Work,Work,-1,ierr)
         lwkopt = MAX(N,N+INT(Work(1)))
         CALL CUNMQR('L','C',N,N,N,B,Ldb,Work,A,Lda,Work,-1,ierr)
         lwkopt = MAX(lwkopt,N+INT(Work(1)))
         IF ( ilvl ) THEN
            CALL CUNGQR(N,N,N,Vl,Ldvl,Work,Work,-1,ierr)
            lwkopt = MAX(lwkopt,N+INT(Work(1)))
         ENDIF
         IF ( ilv ) THEN
            CALL CGGHD3(Jobvl,Jobvr,N,1,N,A,Lda,B,Ldb,Vl,Ldvl,Vr,Ldvr,  &
     &                  Work,-1,ierr)
            lwkopt = MAX(lwkopt,N+INT(Work(1)))
            CALL CHGEQZ('S',Jobvl,Jobvr,N,1,N,A,Lda,B,Ldb,Alpha,Beta,Vl,&
     &                  Ldvl,Vr,Ldvr,Work,-1,Rwork,ierr)
            lwkopt = MAX(lwkopt,N+INT(Work(1)))
         ELSE
            CALL CGGHD3('N','N',N,1,N,A,Lda,B,Ldb,Vl,Ldvl,Vr,Ldvr,Work, &
     &                  -1,ierr)
            lwkopt = MAX(lwkopt,N+INT(Work(1)))
            CALL CHGEQZ('E',Jobvl,Jobvr,N,1,N,A,Lda,B,Ldb,Alpha,Beta,Vl,&
     &                  Ldvl,Vr,Ldvr,Work,-1,Rwork,ierr)
            lwkopt = MAX(lwkopt,N+INT(Work(1)))
         ENDIF
         Work(1) = CMPLX(lwkopt)
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('CGGEV3 ',-Info)
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
      eps = SLAMCH('E')*SLAMCH('B')
      smlnum = SLAMCH('S')
      bignum = ONE/smlnum
      CALL SLABAD(smlnum,bignum)
      smlnum = SQRT(smlnum)/eps
      bignum = ONE/smlnum
!
!     Scale A if max element outside range [SMLNUM,BIGNUM]
!
      anrm = CLANGE('M',N,N,A,Lda,Rwork)
      ilascl = .FALSE.
      IF ( anrm>ZERO .AND. anrm<smlnum ) THEN
         anrmto = smlnum
         ilascl = .TRUE.
      ELSEIF ( anrm>bignum ) THEN
         anrmto = bignum
         ilascl = .TRUE.
      ENDIF
      IF ( ilascl ) CALL CLASCL('G',0,0,anrm,anrmto,N,N,A,Lda,ierr)
!
!     Scale B if max element outside range [SMLNUM,BIGNUM]
!
      bnrm = CLANGE('M',N,N,B,Ldb,Rwork)
      ilbscl = .FALSE.
      IF ( bnrm>ZERO .AND. bnrm<smlnum ) THEN
         bnrmto = smlnum
         ilbscl = .TRUE.
      ELSEIF ( bnrm>bignum ) THEN
         bnrmto = bignum
         ilbscl = .TRUE.
      ENDIF
      IF ( ilbscl ) CALL CLASCL('G',0,0,bnrm,bnrmto,N,N,B,Ldb,ierr)
!
!     Permute the matrices A, B to isolate eigenvalues if possible
!
      ileft = 1
      iright = N + 1
      irwrk = iright + N
      CALL CGGBAL('P',N,A,Lda,B,Ldb,ilo,ihi,Rwork(ileft),Rwork(iright), &
     &            Rwork(irwrk),ierr)
!
!     Reduce B to triangular form (QR decomposition of B)
!
      irows = ihi + 1 - ilo
      IF ( ilv ) THEN
         icols = N + 1 - ilo
      ELSE
         icols = irows
      ENDIF
      itau = 1
      iwrk = itau + irows
      CALL CGEQRF(irows,icols,B(ilo,ilo),Ldb,Work(itau),Work(iwrk),     &
     &            Lwork+1-iwrk,ierr)
!
!     Apply the orthogonal transformation to matrix A
!
      CALL CUNMQR('L','C',irows,icols,irows,B(ilo,ilo),Ldb,Work(itau),  &
     &            A(ilo,ilo),Lda,Work(iwrk),Lwork+1-iwrk,ierr)
!
!     Initialize VL
!
      IF ( ilvl ) THEN
         CALL CLASET('Full',N,N,CZERO,CONE,Vl,Ldvl)
         IF ( irows>1 ) CALL CLACPY('L',irows-1,irows-1,B(ilo+1,ilo),   &
     &                              Ldb,Vl(ilo+1,ilo),Ldvl)
         CALL CUNGQR(irows,irows,irows,Vl(ilo,ilo),Ldvl,Work(itau),     &
     &               Work(iwrk),Lwork+1-iwrk,ierr)
      ENDIF
!
!     Initialize VR
!
      IF ( ilvr ) CALL CLASET('Full',N,N,CZERO,CONE,Vr,Ldvr)
!
!     Reduce to generalized Hessenberg form
!
      IF ( ilv ) THEN
!
!        Eigenvectors requested -- work on whole matrix.
!
         CALL CGGHD3(Jobvl,Jobvr,N,ilo,ihi,A,Lda,B,Ldb,Vl,Ldvl,Vr,Ldvr, &
     &               Work(iwrk),Lwork+1-iwrk,ierr)
      ELSE
         CALL CGGHD3('N','N',irows,1,irows,A(ilo,ilo),Lda,B(ilo,ilo),   &
     &               Ldb,Vl,Ldvl,Vr,Ldvr,Work(iwrk),Lwork+1-iwrk,ierr)
      ENDIF
!
!     Perform QZ algorithm (Compute eigenvalues, and optionally, the
!     Schur form and Schur vectors)
!
      iwrk = itau
      IF ( ilv ) THEN
         chtemp = 'S'
      ELSE
         chtemp = 'E'
      ENDIF
      CALL CHGEQZ(chtemp,Jobvl,Jobvr,N,ilo,ihi,A,Lda,B,Ldb,Alpha,Beta,  &
     &            Vl,Ldvl,Vr,Ldvr,Work(iwrk),Lwork+1-iwrk,Rwork(irwrk), &
     &            ierr)
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
!     Compute Eigenvectors
!
      IF ( ilv ) THEN
         IF ( .NOT.(ilvl) ) THEN
            chtemp = 'R'
         ELSEIF ( ilvr ) THEN
            chtemp = 'B'
         ELSE
            chtemp = 'L'
         ENDIF
!
         CALL CTGEVC(chtemp,'B',ldumma,N,A,Lda,B,Ldb,Vl,Ldvl,Vr,Ldvr,N, &
     &               in,Work(iwrk),Rwork(irwrk),ierr)
         IF ( ierr/=0 ) THEN
            Info = N + 2
            GOTO 100
         ENDIF
!
!        Undo balancing on VL and VR and normalization
!
         IF ( ilvl ) THEN
            CALL CGGBAK('P','L',N,ilo,ihi,Rwork(ileft),Rwork(iright),N, &
     &                  Vl,Ldvl,ierr)
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
         IF ( ilvr ) THEN
            CALL CGGBAK('P','R',N,ilo,ihi,Rwork(ileft),Rwork(iright),N, &
     &                  Vr,Ldvr,ierr)
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
      ENDIF
!
!     Undo scaling if necessary
!
!
 100  IF ( ilascl ) CALL CLASCL('G',0,0,anrmto,anrm,N,1,Alpha,N,ierr)
!
      IF ( ilbscl ) CALL CLASCL('G',0,0,bnrmto,bnrm,N,1,Beta,N,ierr)
!
      Work(1) = CMPLX(lwkopt)
!
!     End of CGGEV3
!
      END SUBROUTINE CGGEV3
