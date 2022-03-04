!*==dgegv.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief <b> DGEEVX computes the eigenvalues and, optionally, the left and/or right eigenvectors for GE matrices</b>
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DGEGV + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgegv.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgegv.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgegv.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DGEGV( JOBVL, JOBVR, N, A, LDA, B, LDB, ALPHAR, ALPHAI,
!                         BETA, VL, LDVL, VR, LDVR, WORK, LWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          JOBVL, JOBVR
!       INTEGER            INFO, LDA, LDB, LDVL, LDVR, LWORK, N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   A( LDA, * ), ALPHAI( * ), ALPHAR( * ),
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
!> This routine is deprecated and has been replaced by routine DGGEV.
!>
!> DGEGV computes the eigenvalues and, optionally, the left and/or right
!> eigenvectors of a real matrix pair (A,B).
!> Given two square matrices A and B,
!> the generalized nonsymmetric eigenvalue problem (GNEP) is to find the
!> eigenvalues lambda and corresponding (non-zero) eigenvectors x such
!> that
!>
!>    A*x = lambda*B*x.
!>
!> An alternate form is to find the eigenvalues mu and corresponding
!> eigenvectors y such that
!>
!>    mu*A*y = B*y.
!>
!> These two forms are equivalent with mu = 1/lambda and x = y if
!> neither lambda nor mu is zero.  In order to deal with the case that
!> lambda or mu is zero or small, two values alpha and beta are returned
!> for each eigenvalue, such that lambda = alpha/beta and
!> mu = beta/alpha.
!>
!> The vectors x and y in the above equations are right eigenvectors of
!> the matrix pair (A,B).  Vectors u and v satisfying
!>
!>    u**H*A = lambda*u**H*B  or  mu*v**H*A = v**H*B
!>
!> are left eigenvectors of (A,B).
!>
!> Note: this routine performs "full balancing" on A and B
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] JOBVL
!> \verbatim
!>          JOBVL is CHARACTER*1
!>          = 'N':  do not compute the left generalized eigenvectors;
!>          = 'V':  compute the left generalized eigenvectors (returned
!>                  in VL).
!> \endverbatim
!>
!> \param[in] JOBVR
!> \verbatim
!>          JOBVR is CHARACTER*1
!>          = 'N':  do not compute the right generalized eigenvectors;
!>          = 'V':  compute the right generalized eigenvectors (returned
!>                  in VR).
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
!>          A is DOUBLE PRECISION array, dimension (LDA, N)
!>          On entry, the matrix A.
!>          If JOBVL = 'V' or JOBVR = 'V', then on exit A
!>          contains the real Schur form of A from the generalized Schur
!>          factorization of the pair (A,B) after balancing.
!>          If no eigenvectors were computed, then only the diagonal
!>          blocks from the Schur form will be correct.  See DGGHRD and
!>          DHGEQZ for details.
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
!>          On entry, the matrix B.
!>          If JOBVL = 'V' or JOBVR = 'V', then on exit B contains the
!>          upper triangular matrix obtained from B in the generalized
!>          Schur factorization of the pair (A,B) after balancing.
!>          If no eigenvectors were computed, then only those elements of
!>          B corresponding to the diagonal blocks from the Schur form of
!>          A will be correct.  See DGGHRD and DHGEQZ for details.
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
!>          ALPHAR is DOUBLE PRECISION array, dimension (N)
!>          The real parts of each scalar alpha defining an eigenvalue of
!>          GNEP.
!> \endverbatim
!>
!> \param[out] ALPHAI
!> \verbatim
!>          ALPHAI is DOUBLE PRECISION array, dimension (N)
!>          The imaginary parts of each scalar alpha defining an
!>          eigenvalue of GNEP.  If ALPHAI(j) is zero, then the j-th
!>          eigenvalue is real; if positive, then the j-th and
!>          (j+1)-st eigenvalues are a complex conjugate pair, with
!>          ALPHAI(j+1) = -ALPHAI(j).
!> \endverbatim
!>
!> \param[out] BETA
!> \verbatim
!>          BETA is DOUBLE PRECISION array, dimension (N)
!>          The scalars beta that define the eigenvalues of GNEP.
!>
!>          Together, the quantities alpha = (ALPHAR(j),ALPHAI(j)) and
!>          beta = BETA(j) represent the j-th eigenvalue of the matrix
!>          pair (A,B), in one of the forms lambda = alpha/beta or
!>          mu = beta/alpha.  Since either lambda or mu may overflow,
!>          they should not, in general, be computed.
!> \endverbatim
!>
!> \param[out] VL
!> \verbatim
!>          VL is DOUBLE PRECISION array, dimension (LDVL,N)
!>          If JOBVL = 'V', the left eigenvectors u(j) are stored
!>          in the columns of VL, in the same order as their eigenvalues.
!>          If the j-th eigenvalue is real, then u(j) = VL(:,j).
!>          If the j-th and (j+1)-st eigenvalues form a complex conjugate
!>          pair, then
!>             u(j) = VL(:,j) + i*VL(:,j+1)
!>          and
!>            u(j+1) = VL(:,j) - i*VL(:,j+1).
!>
!>          Each eigenvector is scaled so that its largest component has
!>          abs(real part) + abs(imag. part) = 1, except for eigenvectors
!>          corresponding to an eigenvalue with alpha = beta = 0, which
!>          are set to zero.
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
!>          VR is DOUBLE PRECISION array, dimension (LDVR,N)
!>          If JOBVR = 'V', the right eigenvectors x(j) are stored
!>          in the columns of VR, in the same order as their eigenvalues.
!>          If the j-th eigenvalue is real, then x(j) = VR(:,j).
!>          If the j-th and (j+1)-st eigenvalues form a complex conjugate
!>          pair, then
!>            x(j) = VR(:,j) + i*VR(:,j+1)
!>          and
!>            x(j+1) = VR(:,j) - i*VR(:,j+1).
!>
!>          Each eigenvector is scaled so that its largest component has
!>          abs(real part) + abs(imag. part) = 1, except for eigenvalues
!>          corresponding to an eigenvalue with alpha = beta = 0, which
!>          are set to zero.
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
!>          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))
!>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The dimension of the array WORK.  LWORK >= max(1,8*N).
!>          For good performance, LWORK must generally be larger.
!>          To compute the optimal value of LWORK, call ILAENV to get
!>          blocksizes (for DGEQRF, DORMQR, and DORGQR.)  Then compute:
!>          NB  -- MAX of the blocksizes for DGEQRF, DORMQR, and DORGQR;
!>          The optimal LWORK is:
!>              2*N + MAX( 6*N, N*(NB+1) ).
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
!>          > N:  errors that usually indicate LAPACK problems:
!>                =N+1: error return from DGGBAL
!>                =N+2: error return from DGEQRF
!>                =N+3: error return from DORMQR
!>                =N+4: error return from DORGQR
!>                =N+5: error return from DGGHRD
!>                =N+6: error return from DHGEQZ (other than failed
!>                                                iteration)
!>                =N+7: error return from DTGEVC
!>                =N+8: error return from DGGBAK (computing VL)
!>                =N+9: error return from DGGBAK (computing VR)
!>                =N+10: error return from DLASCL (various calls)
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
!> \ingroup doubleGEeigen
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  Balancing
!>  ---------
!>
!>  This driver calls DGGBAL to both permute and scale rows and columns
!>  of A and B.  The permutations PL and PR are chosen so that PL*A*PR
!>  and PL*B*R will be upper triangular except for the diagonal blocks
!>  A(i:j,i:j) and B(i:j,i:j), with i and j as close together as
!>  possible.  The diagonal scaling matrices DL and DR are chosen so
!>  that the pair  DL*PL*A*PR*DR, DL*PL*B*PR*DR have elements close to
!>  one (except for the elements that start out zero.)
!>
!>  After the eigenvalues and eigenvectors of the balanced matrices
!>  have been computed, DGGBAK transforms the eigenvectors back to what
!>  they would have been (in perfect arithmetic) if they had not been
!>  balanced.
!>
!>  Contents of A and B on Exit
!>  -------- -- - --- - -- ----
!>
!>  If any eigenvectors are computed (either JOBVL='V' or JOBVR='V' or
!>  both), then on exit the arrays A and B will contain the real Schur
!>  form[*] of the "balanced" versions of A and B.  If no eigenvectors
!>  are computed, then only the diagonal blocks will be correct.
!>
!>  [*] See DHGEQZ, DGEGS, or read the book "Matrix Computations",
!>      by Golub & van Loan, pub. by Johns Hopkins U. Press.
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE DGEGV(Jobvl,Jobvr,N,A,Lda,B,Ldb,Alphar,Alphai,Beta,Vl, &
     &                 Ldvl,Vr,Ldvr,Work,Lwork,Info)
      USE F77KINDS                        
      USE S_DGEQRF
      USE S_DGGBAK
      USE S_DGGBAL
      USE S_DGGHRD
      USE S_DHGEQZ
      USE S_DLACPY
      USE S_DLAMCH
      USE S_DLANGE
      USE S_DLASCL
      USE S_DLASET
      USE S_DORGQR
      USE S_DORMQR
      USE S_DTGEVC
      USE S_ILAENV
      USE S_LSAME
      USE S_XERBLA
      IMPLICIT NONE
!*--DGEGV327
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D0 , ONE = 1.0D0
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Jobvl
      CHARACTER :: Jobvr
      INTEGER :: N
      REAL(R8KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Alphar
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Alphai
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Beta
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldvl,*) :: Vl
      INTEGER :: Ldvl
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldvr,*) :: Vr
      INTEGER :: Ldvr
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      REAL(R8KIND) :: absai , absar , absb , anrm , anrm1 , anrm2 ,     &
     &                bnrm , bnrm1 , bnrm2 , eps , onepls , safmax ,    &
     &                safmin , salfai , salfar , sbeta , scale , temp
      CHARACTER :: chtemp
      INTEGER :: icols , ihi , iinfo , ijobvl , ijobvr , ileft , ilo ,  &
     &           in , iright , irows , itau , iwork , jc , jr , lopt ,  &
     &           lwkmin , lwkopt , nb , nb1 , nb2 , nb3
      LOGICAL :: ilimit , ilv , ilvl , ilvr , lquery
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
      lwkmin = MAX(8*N,1)
      lwkopt = lwkmin
      Work(1) = lwkopt
      lquery = (Lwork==-1)
      Info = 0
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
      ELSEIF ( Lwork<lwkmin .AND. .NOT.lquery ) THEN
         Info = -16
      ENDIF
!
      IF ( Info==0 ) THEN
         nb1 = ILAENV(1,'DGEQRF',' ',N,N,-1,-1)
         nb2 = ILAENV(1,'DORMQR',' ',N,N,N,-1)
         nb3 = ILAENV(1,'DORGQR',' ',N,N,N,-1)
         nb = MAX(nb1,nb2,nb3)
         lopt = 2*N + MAX(6*N,N*(nb+1))
         Work(1) = lopt
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('DGEGV ',-Info)
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
      eps = DLAMCH('E')*DLAMCH('B')
      safmin = DLAMCH('S')
      safmin = safmin + safmin
      safmax = ONE/safmin
      onepls = ONE + (4*eps)
!
!     Scale A
!
      anrm = DLANGE('M',N,N,A,Lda,Work)
      anrm1 = anrm
      anrm2 = ONE
      IF ( anrm<ONE ) THEN
         IF ( safmax*anrm<ONE ) THEN
            anrm1 = safmin
            anrm2 = safmax*anrm
         ENDIF
      ENDIF
!
      IF ( anrm>ZERO ) THEN
         CALL DLASCL('G',-1,-1,anrm,ONE,N,N,A,Lda,iinfo)
         IF ( iinfo/=0 ) THEN
            Info = N + 10
            RETURN
         ENDIF
      ENDIF
!
!     Scale B
!
      bnrm = DLANGE('M',N,N,B,Ldb,Work)
      bnrm1 = bnrm
      bnrm2 = ONE
      IF ( bnrm<ONE ) THEN
         IF ( safmax*bnrm<ONE ) THEN
            bnrm1 = safmin
            bnrm2 = safmax*bnrm
         ENDIF
      ENDIF
!
      IF ( bnrm>ZERO ) THEN
         CALL DLASCL('G',-1,-1,bnrm,ONE,N,N,B,Ldb,iinfo)
         IF ( iinfo/=0 ) THEN
            Info = N + 10
            RETURN
         ENDIF
      ENDIF
!
!     Permute the matrix to make it more nearly triangular
!     Workspace layout:  (8*N words -- "work" requires 6*N words)
!        left_permutation, right_permutation, work...
!
      ileft = 1
      iright = N + 1
      iwork = iright + N
      CALL DGGBAL('P',N,A,Lda,B,Ldb,ilo,ihi,Work(ileft),Work(iright),   &
     &            Work(iwork),iinfo)
      IF ( iinfo/=0 ) THEN
         Info = N + 1
         GOTO 100
      ENDIF
!
!     Reduce B to triangular form, and initialize VL and/or VR
!     Workspace layout:  ("work..." must have at least N words)
!        left_permutation, right_permutation, tau, work...
!
      irows = ihi + 1 - ilo
      IF ( ilv ) THEN
         icols = N + 1 - ilo
      ELSE
         icols = irows
      ENDIF
      itau = iwork
      iwork = itau + irows
      CALL DGEQRF(irows,icols,B(ilo,ilo),Ldb,Work(itau),Work(iwork),    &
     &            Lwork+1-iwork,iinfo)
      IF ( iinfo>=0 ) lwkopt = MAX(lwkopt,INT(Work(iwork))+iwork-1)
      IF ( iinfo/=0 ) THEN
         Info = N + 2
         GOTO 100
      ENDIF
!
      CALL DORMQR('L','T',irows,icols,irows,B(ilo,ilo),Ldb,Work(itau),  &
     &            A(ilo,ilo),Lda,Work(iwork),Lwork+1-iwork,iinfo)
      IF ( iinfo>=0 ) lwkopt = MAX(lwkopt,INT(Work(iwork))+iwork-1)
      IF ( iinfo/=0 ) THEN
         Info = N + 3
         GOTO 100
      ENDIF
!
      IF ( ilvl ) THEN
         CALL DLASET('Full',N,N,ZERO,ONE,Vl,Ldvl)
         CALL DLACPY('L',irows-1,irows-1,B(ilo+1,ilo),Ldb,Vl(ilo+1,ilo),&
     &               Ldvl)
         CALL DORGQR(irows,irows,irows,Vl(ilo,ilo),Ldvl,Work(itau),     &
     &               Work(iwork),Lwork+1-iwork,iinfo)
         IF ( iinfo>=0 ) lwkopt = MAX(lwkopt,INT(Work(iwork))+iwork-1)
         IF ( iinfo/=0 ) THEN
            Info = N + 4
            GOTO 100
         ENDIF
      ENDIF
!
      IF ( ilvr ) CALL DLASET('Full',N,N,ZERO,ONE,Vr,Ldvr)
!
!     Reduce to generalized Hessenberg form
!
      IF ( ilv ) THEN
!
!        Eigenvectors requested -- work on whole matrix.
!
         CALL DGGHRD(Jobvl,Jobvr,N,ilo,ihi,A,Lda,B,Ldb,Vl,Ldvl,Vr,Ldvr, &
     &               iinfo)
      ELSE
         CALL DGGHRD('N','N',irows,1,irows,A(ilo,ilo),Lda,B(ilo,ilo),   &
     &               Ldb,Vl,Ldvl,Vr,Ldvr,iinfo)
      ENDIF
      IF ( iinfo/=0 ) THEN
         Info = N + 5
         GOTO 100
      ENDIF
!
!     Perform QZ algorithm
!     Workspace layout:  ("work..." must have at least 1 word)
!        left_permutation, right_permutation, work...
!
      iwork = itau
      IF ( ilv ) THEN
         chtemp = 'S'
      ELSE
         chtemp = 'E'
      ENDIF
      CALL DHGEQZ(chtemp,Jobvl,Jobvr,N,ilo,ihi,A,Lda,B,Ldb,Alphar,      &
     &            Alphai,Beta,Vl,Ldvl,Vr,Ldvr,Work(iwork),Lwork+1-iwork,&
     &            iinfo)
      IF ( iinfo>=0 ) lwkopt = MAX(lwkopt,INT(Work(iwork))+iwork-1)
      IF ( iinfo/=0 ) THEN
         IF ( iinfo>0 .AND. iinfo<=N ) THEN
            Info = iinfo
         ELSEIF ( iinfo>N .AND. iinfo<=2*N ) THEN
            Info = iinfo - N
         ELSE
            Info = N + 6
         ENDIF
         GOTO 100
      ENDIF
!
      IF ( ilv ) THEN
!
!        Compute Eigenvectors  (DTGEVC requires 6*N words of workspace)
!
         IF ( .NOT.(ilvl) ) THEN
            chtemp = 'R'
         ELSEIF ( ilvr ) THEN
            chtemp = 'B'
         ELSE
            chtemp = 'L'
         ENDIF
!
         CALL DTGEVC(chtemp,'B',ldumma,N,A,Lda,B,Ldb,Vl,Ldvl,Vr,Ldvr,N, &
     &               in,Work(iwork),iinfo)
         IF ( iinfo/=0 ) THEN
            Info = N + 7
            GOTO 100
         ENDIF
!
!        Undo balancing on VL and VR, rescale
!
         IF ( ilvl ) THEN
            CALL DGGBAK('P','L',N,ilo,ihi,Work(ileft),Work(iright),N,Vl,&
     &                  Ldvl,iinfo)
            IF ( iinfo/=0 ) THEN
               Info = N + 8
               GOTO 100
            ENDIF
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
                  IF ( temp>=safmin ) THEN
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
            CALL DGGBAK('P','R',N,ilo,ihi,Work(ileft),Work(iright),N,Vr,&
     &                  Ldvr,iinfo)
            IF ( iinfo/=0 ) THEN
               Info = N + 9
               GOTO 100
            ENDIF
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
                  IF ( temp>=safmin ) THEN
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
!     Undo scaling in alpha, beta
!
!     Note: this does not give the alpha and beta for the unscaled
!     problem.
!
!     Un-scaling is limited to avoid underflow in alpha and beta
!     if they are significant.
!
      DO jc = 1 , N
         absar = ABS(Alphar(jc))
         absai = ABS(Alphai(jc))
         absb = ABS(Beta(jc))
         salfar = anrm*Alphar(jc)
         salfai = anrm*Alphai(jc)
         sbeta = bnrm*Beta(jc)
         ilimit = .FALSE.
         scale = ONE
!
!        Check for significant underflow in ALPHAI
!
         IF ( ABS(salfai)<safmin .AND.                                  &
     &        absai>=MAX(safmin,eps*absar,eps*absb) ) THEN
            ilimit = .TRUE.
            scale = (onepls*safmin/anrm1)/MAX(onepls*safmin,anrm2*absai)
!
         ELSEIF ( salfai==ZERO ) THEN
!
!           If insignificant underflow in ALPHAI, then make the
!           conjugate eigenvalue real.
!
            IF ( Alphai(jc)<ZERO .AND. jc>1 ) THEN
               Alphai(jc-1) = ZERO
            ELSEIF ( Alphai(jc)>ZERO .AND. jc<N ) THEN
               Alphai(jc+1) = ZERO
            ENDIF
         ENDIF
!
!        Check for significant underflow in ALPHAR
!
         IF ( ABS(salfar)<safmin .AND.                                  &
     &        absar>=MAX(safmin,eps*absai,eps*absb) ) THEN
            ilimit = .TRUE.
            scale = MAX(scale,(onepls*safmin/anrm1)                     &
     &              /MAX(onepls*safmin,anrm2*absar))
         ENDIF
!
!        Check for significant underflow in BETA
!
         IF ( ABS(sbeta)<safmin .AND.                                   &
     &        absb>=MAX(safmin,eps*absar,eps*absai) ) THEN
            ilimit = .TRUE.
            scale = MAX(scale,(onepls*safmin/bnrm1)                     &
     &              /MAX(onepls*safmin,bnrm2*absb))
         ENDIF
!
!        Check for possible overflow when limiting scaling
!
         IF ( ilimit ) THEN
            temp = (scale*safmin)                                       &
     &             *MAX(ABS(salfar),ABS(salfai),ABS(sbeta))
            IF ( temp>ONE ) scale = scale/temp
            IF ( scale<ONE ) ilimit = .FALSE.
         ENDIF
!
!        Recompute un-scaled ALPHAR, ALPHAI, BETA if necessary.
!
         IF ( ilimit ) THEN
            salfar = (scale*Alphar(jc))*anrm
            salfai = (scale*Alphai(jc))*anrm
            sbeta = (scale*Beta(jc))*bnrm
         ENDIF
         Alphar(jc) = salfar
         Alphai(jc) = salfai
         Beta(jc) = sbeta
      ENDDO
!
 100  Work(1) = lwkopt
!
!
!     End of DGEGV
!
      END SUBROUTINE DGEGV
