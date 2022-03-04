!*==sggev3.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief <b> SGGEV3 computes the eigenvalues and, optionally, the left and/or right eigenvectors for GE matrices (blocked algorithm)</b>
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SGGEV3 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sggev3.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sggev3.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sggev3.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SGGEV3( JOBVL, JOBVR, N, A, LDA, B, LDB, ALPHAR,
!      $                   ALPHAI, BETA, VL, LDVL, VR, LDVR, WORK, LWORK,
!      $                   INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          JOBVL, JOBVR
!       INTEGER            INFO, LDA, LDB, LDVL, LDVR, LWORK, N
!       ..
!       .. Array Arguments ..
!       REAL               A( LDA, * ), ALPHAI( * ), ALPHAR( * ),
!      $                   B( LDB, * ), BETA( * ), VL( LDVL, * ),
!      $                   VR( LDVR, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SGGEV3 computes for a pair of N-by-N real nonsymmetric matrices (A,B)
!> the generalized eigenvalues, and optionally, the left and/or right
!> generalized eigenvectors.
!>
!> A generalized eigenvalue for a pair of matrices (A,B) is a scalar
!> lambda or a ratio alpha/beta = lambda, such that A - lambda*B is
!> singular. It is usually represented as the pair (alpha,beta), as
!> there is a reasonable interpretation for beta=0, and even for both
!> being zero.
!>
!> The right eigenvector v(j) corresponding to the eigenvalue lambda(j)
!> of (A,B) satisfies
!>
!>                  A * v(j) = lambda(j) * B * v(j).
!>
!> The left eigenvector u(j) corresponding to the eigenvalue lambda(j)
!> of (A,B) satisfies
!>
!>                  u(j)**H * A  = lambda(j) * u(j)**H * B .
!>
!> where u(j)**H is the conjugate-transpose of u(j).
!>
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
!>          A is REAL array, dimension (LDA, N)
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
!>          B is REAL array, dimension (LDB, N)
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
!>          be the generalized eigenvalues.  If ALPHAI(j) is zero, then
!>          the j-th eigenvalue is real; if positive, then the j-th and
!>          (j+1)-st eigenvalues are a complex conjugate pair, with
!>          ALPHAI(j+1) negative.
!>
!>          Note: the quotients ALPHAR(j)/BETA(j) and ALPHAI(j)/BETA(j)
!>          may easily over- or underflow, and BETA(j) may even be zero.
!>          Thus, the user should avoid naively computing the ratio
!>          alpha/beta.  However, ALPHAR and ALPHAI will be always less
!>          than and usually comparable with norm(A) in magnitude, and
!>          BETA always less than and usually comparable with norm(B).
!> \endverbatim
!>
!> \param[out] VL
!> \verbatim
!>          VL is REAL array, dimension (LDVL,N)
!>          If JOBVL = 'V', the left eigenvectors u(j) are stored one
!>          after another in the columns of VL, in the same order as
!>          their eigenvalues. If the j-th eigenvalue is real, then
!>          u(j) = VL(:,j), the j-th column of VL. If the j-th and
!>          (j+1)-th eigenvalues form a complex conjugate pair, then
!>          u(j) = VL(:,j)+i*VL(:,j+1) and u(j+1) = VL(:,j)-i*VL(:,j+1).
!>          Each eigenvector is scaled so the largest component has
!>          abs(real part)+abs(imag. part)=1.
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
!>          VR is REAL array, dimension (LDVR,N)
!>          If JOBVR = 'V', the right eigenvectors v(j) are stored one
!>          after another in the columns of VR, in the same order as
!>          their eigenvalues. If the j-th eigenvalue is real, then
!>          v(j) = VR(:,j), the j-th column of VR. If the j-th and
!>          (j+1)-th eigenvalues form a complex conjugate pair, then
!>          v(j) = VR(:,j)+i*VR(:,j+1) and v(j+1) = VR(:,j)-i*VR(:,j+1).
!>          Each eigenvector is scaled so the largest component has
!>          abs(real part)+abs(imag. part)=1.
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
!>          WORK is REAL array, dimension (MAX(1,LWORK))
!>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>
!>          If LWORK = -1, then a workspace query is assumed; the routine
!>          only calculates the optimal size of the WORK array, returns
!>          this value as the first entry of the WORK array, and no error
!>          message related to LWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value.
!>          = 1,...,N:
!>                The QZ iteration failed.  No eigenvectors have been
!>                calculated, but ALPHAR(j), ALPHAI(j), and BETA(j)
!>                should be correct for j=INFO+1,...,N.
!>          > N:  =N+1: other than QZ iteration failed in SHGEQZ.
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
!> \ingroup realGEeigen
!
!  =====================================================================
      SUBROUTINE SGGEV3(Jobvl,Jobvr,N,A,Lda,B,Ldb,Alphar,Alphai,Beta,Vl,&
     &                  Ldvl,Vr,Ldvr,Work,Lwork,Info)
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
      USE S_STGEVC
      USE S_XERBLA
      IMPLICIT NONE
!*--SGGEV3245
!
! PARAMETER definitions rewritten by SPAG
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Jobvl
      CHARACTER :: Jobvr
      INTEGER :: N
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL , DIMENSION(*) :: Alphar
      REAL , DIMENSION(*) :: Alphai
      REAL , DIMENSION(*) :: Beta
      REAL , INTENT(INOUT) , DIMENSION(Ldvl,*) :: Vl
      INTEGER :: Ldvl
      REAL , INTENT(INOUT) , DIMENSION(Ldvr,*) :: Vr
      INTEGER :: Ldvr
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      REAL :: anrm , anrmto , bignum , bnrm , bnrmto , eps , smlnum ,   &
     &        temp
      CHARACTER :: chtemp
      INTEGER :: icols , ierr , ihi , ijobvl , ijobvr , ileft , ilo ,   &
     &           in , iright , irows , itau , iwrk , jc , jr , lwkopt
      LOGICAL :: ilascl , ilbscl , ilv , ilvl , ilvr , lquery
      LOGICAL , DIMENSION(1) :: ldumma
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
         Info = -12
      ELSEIF ( Ldvr<1 .OR. (ilvr .AND. Ldvr<N) ) THEN
         Info = -14
      ELSEIF ( Lwork<MAX(1,8*N) .AND. .NOT.lquery ) THEN
         Info = -16
      ENDIF
!
!     Compute workspace
!
      IF ( Info==0 ) THEN
         CALL SGEQRF(N,N,B,Ldb,Work,Work,-1,ierr)
         lwkopt = MAX(1,8*N,3*N+INT(Work(1)))
         CALL SORMQR('L','T',N,N,N,B,Ldb,Work,A,Lda,Work,-1,ierr)
         lwkopt = MAX(lwkopt,3*N+INT(Work(1)))
         CALL SGGHD3(Jobvl,Jobvr,N,1,N,A,Lda,B,Ldb,Vl,Ldvl,Vr,Ldvr,Work,&
     &               -1,ierr)
         lwkopt = MAX(lwkopt,3*N+INT(Work(1)))
         IF ( ilvl ) THEN
            CALL SORGQR(N,N,N,Vl,Ldvl,Work,Work,-1,ierr)
            lwkopt = MAX(lwkopt,3*N+INT(Work(1)))
            CALL SHGEQZ('S',Jobvl,Jobvr,N,1,N,A,Lda,B,Ldb,Alphar,Alphai,&
     &                  Beta,Vl,Ldvl,Vr,Ldvr,Work,-1,ierr)
            lwkopt = MAX(lwkopt,2*N+INT(Work(1)))
         ELSE
            CALL SHGEQZ('E',Jobvl,Jobvr,N,1,N,A,Lda,B,Ldb,Alphar,Alphai,&
     &                  Beta,Vl,Ldvl,Vr,Ldvr,Work,-1,ierr)
            lwkopt = MAX(lwkopt,2*N+INT(Work(1)))
         ENDIF
         Work(1) = REAL(lwkopt)
!
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('SGGEV3 ',-Info)
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
      smlnum = SLAMCH('S')
      bignum = ONE/smlnum
      CALL SLABAD(smlnum,bignum)
      smlnum = SQRT(smlnum)/eps
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
!     Permute the matrices A, B to isolate eigenvalues if possible
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
      IF ( ilv ) THEN
         icols = N + 1 - ilo
      ELSE
         icols = irows
      ENDIF
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
!     Initialize VL
!
      IF ( ilvl ) THEN
         CALL SLASET('Full',N,N,ZERO,ONE,Vl,Ldvl)
         IF ( irows>1 ) CALL SLACPY('L',irows-1,irows-1,B(ilo+1,ilo),   &
     &                              Ldb,Vl(ilo+1,ilo),Ldvl)
         CALL SORGQR(irows,irows,irows,Vl(ilo,ilo),Ldvl,Work(itau),     &
     &               Work(iwrk),Lwork+1-iwrk,ierr)
      ENDIF
!
!     Initialize VR
!
      IF ( ilvr ) CALL SLASET('Full',N,N,ZERO,ONE,Vr,Ldvr)
!
!     Reduce to generalized Hessenberg form
!
      IF ( ilv ) THEN
!
!        Eigenvectors requested -- work on whole matrix.
!
         CALL SGGHD3(Jobvl,Jobvr,N,ilo,ihi,A,Lda,B,Ldb,Vl,Ldvl,Vr,Ldvr, &
     &               Work(iwrk),Lwork+1-iwrk,ierr)
      ELSE
         CALL SGGHD3('N','N',irows,1,irows,A(ilo,ilo),Lda,B(ilo,ilo),   &
     &               Ldb,Vl,Ldvl,Vr,Ldvr,Work(iwrk),Lwork+1-iwrk,ierr)
      ENDIF
!
!     Perform QZ algorithm (Compute eigenvalues, and optionally, the
!     Schur forms and Schur vectors)
!
      iwrk = itau
      IF ( ilv ) THEN
         chtemp = 'S'
      ELSE
         chtemp = 'E'
      ENDIF
      CALL SHGEQZ(chtemp,Jobvl,Jobvr,N,ilo,ihi,A,Lda,B,Ldb,Alphar,      &
     &            Alphai,Beta,Vl,Ldvl,Vr,Ldvr,Work(iwrk),Lwork+1-iwrk,  &
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
         CALL STGEVC(chtemp,'B',ldumma,N,A,Lda,B,Ldb,Vl,Ldvl,Vr,Ldvr,N, &
     &               in,Work(iwrk),ierr)
         IF ( ierr/=0 ) THEN
            Info = N + 2
            GOTO 100
         ENDIF
!
!        Undo balancing on VL and VR and normalization
!
         IF ( ilvl ) THEN
            CALL SGGBAK('P','L',N,ilo,ihi,Work(ileft),Work(iright),N,Vl,&
     &                  Ldvl,ierr)
            DO jc = 1 , N
               IF ( Alphai(jc)>=ZERO ) THEN
                  temp = ZERO
                  IF ( Alphai(jc)==ZERO ) THEN
                     DO jr = 1 , N
                        temp = MAX(temp,ABS(Vl(jr,jc)))
                     ENDDO
                  ELSE
                     DO jr = 1 , N
                        temp = MAX(temp,ABS(Vl(jr,jc))+ABS(Vl(jr,jc+1)))
                     ENDDO
                  ENDIF
                  IF ( temp>=smlnum ) THEN
                     temp = ONE/temp
                     IF ( Alphai(jc)==ZERO ) THEN
                        DO jr = 1 , N
                           Vl(jr,jc) = Vl(jr,jc)*temp
                        ENDDO
                     ELSE
                        DO jr = 1 , N
                           Vl(jr,jc) = Vl(jr,jc)*temp
                           Vl(jr,jc+1) = Vl(jr,jc+1)*temp
                        ENDDO
                     ENDIF
                  ENDIF
               ENDIF
            ENDDO
         ENDIF
         IF ( ilvr ) THEN
            CALL SGGBAK('P','R',N,ilo,ihi,Work(ileft),Work(iright),N,Vr,&
     &                  Ldvr,ierr)
            DO jc = 1 , N
               IF ( Alphai(jc)>=ZERO ) THEN
                  temp = ZERO
                  IF ( Alphai(jc)==ZERO ) THEN
                     DO jr = 1 , N
                        temp = MAX(temp,ABS(Vr(jr,jc)))
                     ENDDO
                  ELSE
                     DO jr = 1 , N
                        temp = MAX(temp,ABS(Vr(jr,jc))+ABS(Vr(jr,jc+1)))
                     ENDDO
                  ENDIF
                  IF ( temp>=smlnum ) THEN
                     temp = ONE/temp
                     IF ( Alphai(jc)==ZERO ) THEN
                        DO jr = 1 , N
                           Vr(jr,jc) = Vr(jr,jc)*temp
                        ENDDO
                     ELSE
                        DO jr = 1 , N
                           Vr(jr,jc) = Vr(jr,jc)*temp
                           Vr(jr,jc+1) = Vr(jr,jc+1)*temp
                        ENDDO
                     ENDIF
                  ENDIF
               ENDIF
            ENDDO
         ENDIF
!
!        End of eigenvector calculation
!
      ENDIF
!
!     Undo scaling if necessary
!
!
 100  IF ( ilascl ) THEN
         CALL SLASCL('G',0,0,anrmto,anrm,N,1,Alphar,N,ierr)
         CALL SLASCL('G',0,0,anrmto,anrm,N,1,Alphai,N,ierr)
      ENDIF
!
      IF ( ilbscl ) CALL SLASCL('G',0,0,bnrmto,bnrm,N,1,Beta,N,ierr)
!
      Work(1) = REAL(lwkopt)
!
!     End of SGGEV3
!
      END SUBROUTINE SGGEV3