!*==cgges3.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief <b> CGGES3 computes the eigenvalues, the Schur form, and, optionally, the matrix of Schur vectors for GE matrices (blocked algorithm)</b>
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CGGES3 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgges3.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgges3.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgges3.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CGGES3( JOBVSL, JOBVSR, SORT, SELCTG, N, A, LDA, B,
!      $                   LDB, SDIM, ALPHA, BETA, VSL, LDVSL, VSR, LDVSR,
!      $                   WORK, LWORK, RWORK, BWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          JOBVSL, JOBVSR, SORT
!       INTEGER            INFO, LDA, LDB, LDVSL, LDVSR, LWORK, N, SDIM
!       ..
!       .. Array Arguments ..
!       LOGICAL            BWORK( * )
!       REAL               RWORK( * )
!       COMPLEX            A( LDA, * ), ALPHA( * ), B( LDB, * ),
!      $                   BETA( * ), VSL( LDVSL, * ), VSR( LDVSR, * ),
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
!> CGGES3 computes for a pair of N-by-N complex nonsymmetric matrices
!> (A,B), the generalized eigenvalues, the generalized complex Schur
!> form (S, T), and optionally left and/or right Schur vectors (VSL
!> and VSR). This gives the generalized Schur factorization
!>
!>         (A,B) = ( (VSL)*S*(VSR)**H, (VSL)*T*(VSR)**H )
!>
!> where (VSR)**H is the conjugate-transpose of VSR.
!>
!> Optionally, it also orders the eigenvalues so that a selected cluster
!> of eigenvalues appears in the leading diagonal blocks of the upper
!> triangular matrix S and the upper triangular matrix T. The leading
!> columns of VSL and VSR then form an unitary basis for the
!> corresponding left and right eigenspaces (deflating subspaces).
!>
!> (If only the generalized eigenvalues are needed, use the driver
!> CGGEV instead, which is faster.)
!>
!> A generalized eigenvalue for a pair of matrices (A,B) is a scalar w
!> or a ratio alpha/beta = w, such that  A - w*B is singular.  It is
!> usually represented as the pair (alpha,beta), as there is a
!> reasonable interpretation for beta=0, and even for both being zero.
!>
!> A pair of matrices (S,T) is in generalized complex Schur form if S
!> and T are upper triangular and, in addition, the diagonal elements
!> of T are non-negative real numbers.
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
!>          SELCTG is a LOGICAL FUNCTION of two COMPLEX arguments
!>          SELCTG must be declared EXTERNAL in the calling subroutine.
!>          If SORT = 'N', SELCTG is not referenced.
!>          If SORT = 'S', SELCTG is used to select eigenvalues to sort
!>          to the top left of the Schur form.
!>          An eigenvalue ALPHA(j)/BETA(j) is selected if
!>          SELCTG(ALPHA(j),BETA(j)) is true.
!>
!>          Note that a selected complex eigenvalue may no longer satisfy
!>          SELCTG(ALPHA(j),BETA(j)) = .TRUE. after ordering, since
!>          ordering may change the value of complex eigenvalues
!>          (especially if the eigenvalue is ill-conditioned), in this
!>          case INFO is set to N+2 (See INFO below).
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
!>          A is COMPLEX array, dimension (LDA, N)
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
!>          B is COMPLEX array, dimension (LDB, N)
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
!>          for which SELCTG is true.
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
!>          On exit,  ALPHA(j)/BETA(j), j=1,...,N, will be the
!>          generalized eigenvalues.  ALPHA(j), j=1,...,N  and  BETA(j),
!>          j=1,...,N  are the diagonals of the complex Schur form (A,B)
!>          output by CGGES3. The  BETA(j) will be non-negative real.
!>
!>          Note: the quotients ALPHA(j)/BETA(j) may easily over- or
!>          underflow, and BETA(j) may even be zero.  Thus, the user
!>          should avoid naively computing the ratio alpha/beta.
!>          However, ALPHA will be always less than and usually
!>          comparable with norm(A) in magnitude, and BETA always less
!>          than and usually comparable with norm(B).
!> \endverbatim
!>
!> \param[out] VSL
!> \verbatim
!>          VSL is COMPLEX array, dimension (LDVSL,N)
!>          If JOBVSL = 'V', VSL will contain the left Schur vectors.
!>          Not referenced if JOBVSL = 'N'.
!> \endverbatim
!>
!> \param[in] LDVSL
!> \verbatim
!>          LDVSL is INTEGER
!>          The leading dimension of the matrix VSL. LDVSL >= 1, and
!>          if JOBVSL = 'V', LDVSL >= N.
!> \endverbatim
!>
!> \param[out] VSR
!> \verbatim
!>          VSR is COMPLEX array, dimension (LDVSR,N)
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
!>          =1,...,N:
!>                The QZ iteration failed.  (A,B) are not in Schur
!>                form, but ALPHA(j) and BETA(j) should be correct for
!>                j=INFO+1,...,N.
!>          > N:  =N+1: other than QZ iteration failed in CHGEQZ
!>                =N+2: after reordering, roundoff changed values of
!>                      some complex eigenvalues so that leading
!>                      eigenvalues in the Generalized Schur form no
!>                      longer satisfy SELCTG=.TRUE.  This could also
!>                      be caused due to scaling.
!>                =N+3: reordering failed in CTGSEN.
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
      SUBROUTINE CGGES3(Jobvsl,Jobvsr,Sort,SELCTG,N,A,Lda,B,Ldb,Sdim,   &
     &                  Alpha,Beta,Vsl,Ldvsl,Vsr,Ldvsr,Work,Lwork,Rwork,&
     &                  Bwork,Info)
      IMPLICIT NONE
!*--CGGES3273
!
!  -- LAPACK driver routine (version 3.6.1) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     January 2015
!
!     .. Scalar Arguments ..
      CHARACTER Jobvsl , Jobvsr , Sort
      INTEGER Info , Lda , Ldb , Ldvsl , Ldvsr , Lwork , N , Sdim
!     ..
!     .. Array Arguments ..
      LOGICAL Bwork(*)
      REAL Rwork(*)
      COMPLEX A(Lda,*) , Alpha(*) , B(Ldb,*) , Beta(*) , Vsl(Ldvsl,*) , &
     &        Vsr(Ldvsr,*) , Work(*)
!     ..
!     .. Function Arguments ..
      LOGICAL SELCTG
      EXTERNAL SELCTG
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ZERO , ONE
      PARAMETER (ZERO=0.0E0,ONE=1.0E0)
      COMPLEX CZERO , CONE
      PARAMETER (CZERO=(0.0E0,0.0E0),CONE=(1.0E0,0.0E0))
!     ..
!     .. Local Scalars ..
      LOGICAL cursl , ilascl , ilbscl , ilvsl , ilvsr , lastsl ,        &
     &        lquery , wantst
      INTEGER i , icols , ierr , ihi , ijobvl , ijobvr , ileft , ilo ,  &
     &        iright , irows , irwrk , itau , iwrk , lwkopt
      REAL anrm , anrmto , bignum , bnrm , bnrmto , eps , pvsl , pvsr , &
     &     smlnum
!     ..
!     .. Local Arrays ..
      INTEGER idum(1)
      REAL dif(2)
!     ..
!     .. External Subroutines ..
      EXTERNAL CGEQRF , CGGBAK , CGGBAL , CGGHD3 , CHGEQZ , CLACPY ,    &
     &         CLASCL , CLASET , CTGSEN , CUNGQR , CUNMQR , SLABAD ,    &
     &         XERBLA
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      REAL CLANGE , SLAMCH
      EXTERNAL LSAME , CLANGE , SLAMCH
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX , SQRT
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
         Info = -14
      ELSEIF ( Ldvsr<1 .OR. (ilvsr .AND. Ldvsr<N) ) THEN
         Info = -16
      ELSEIF ( Lwork<MAX(1,2*N) .AND. .NOT.lquery ) THEN
         Info = -18
      ENDIF
!
!     Compute workspace
!
      IF ( Info==0 ) THEN
         CALL CGEQRF(N,N,B,Ldb,Work,Work,-1,ierr)
         lwkopt = MAX(1,N+INT(Work(1)))
         CALL CUNMQR('L','C',N,N,N,B,Ldb,Work,A,Lda,Work,-1,ierr)
         lwkopt = MAX(lwkopt,N+INT(Work(1)))
         IF ( ilvsl ) THEN
            CALL CUNGQR(N,N,N,Vsl,Ldvsl,Work,Work,-1,ierr)
            lwkopt = MAX(lwkopt,N+INT(Work(1)))
         ENDIF
         CALL CGGHD3(Jobvsl,Jobvsr,N,1,N,A,Lda,B,Ldb,Vsl,Ldvsl,Vsr,     &
     &               Ldvsr,Work,-1,ierr)
         lwkopt = MAX(lwkopt,N+INT(Work(1)))
         CALL CHGEQZ('S',Jobvsl,Jobvsr,N,1,N,A,Lda,B,Ldb,Alpha,Beta,Vsl,&
     &               Ldvsl,Vsr,Ldvsr,Work,-1,Rwork,ierr)
         lwkopt = MAX(lwkopt,INT(Work(1)))
         IF ( wantst ) THEN
            CALL CTGSEN(0,ilvsl,ilvsr,Bwork,N,A,Lda,B,Ldb,Alpha,Beta,   &
     &                  Vsl,Ldvsl,Vsr,Ldvsr,Sdim,pvsl,pvsr,dif,Work,-1, &
     &                  idum,1,ierr)
            lwkopt = MAX(lwkopt,INT(Work(1)))
         ENDIF
         Work(1) = CMPLX(lwkopt)
      ENDIF
 
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('CGGES3 ',-Info)
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
!
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
!
      IF ( ilbscl ) CALL CLASCL('G',0,0,bnrm,bnrmto,N,N,B,Ldb,ierr)
!
!     Permute the matrix to make it more nearly triangular
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
      icols = N + 1 - ilo
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
!     Initialize VSL
!
      IF ( ilvsl ) THEN
         CALL CLASET('Full',N,N,CZERO,CONE,Vsl,Ldvsl)
         IF ( irows>1 ) CALL CLACPY('L',irows-1,irows-1,B(ilo+1,ilo),   &
     &                              Ldb,Vsl(ilo+1,ilo),Ldvsl)
         CALL CUNGQR(irows,irows,irows,Vsl(ilo,ilo),Ldvsl,Work(itau),   &
     &               Work(iwrk),Lwork+1-iwrk,ierr)
      ENDIF
!
!     Initialize VSR
!
      IF ( ilvsr ) CALL CLASET('Full',N,N,CZERO,CONE,Vsr,Ldvsr)
!
!     Reduce to generalized Hessenberg form
!
      CALL CGGHD3(Jobvsl,Jobvsr,N,ilo,ihi,A,Lda,B,Ldb,Vsl,Ldvsl,Vsr,    &
     &            Ldvsr,Work(iwrk),Lwork+1-iwrk,ierr)
!
      Sdim = 0
!
!     Perform QZ algorithm, computing Schur vectors if desired
!
      iwrk = itau
      CALL CHGEQZ('S',Jobvsl,Jobvsr,N,ilo,ihi,A,Lda,B,Ldb,Alpha,Beta,   &
     &            Vsl,Ldvsl,Vsr,Ldvsr,Work(iwrk),Lwork+1-iwrk,          &
     &            Rwork(irwrk),ierr)
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
      IF ( wantst ) THEN
!
!        Undo scaling on eigenvalues before selecting
!
         IF ( ilascl ) CALL CLASCL('G',0,0,anrm,anrmto,N,1,Alpha,N,ierr)
         IF ( ilbscl ) CALL CLASCL('G',0,0,bnrm,bnrmto,N,1,Beta,N,ierr)
!
!        Select eigenvalues
!
         DO i = 1 , N
            Bwork(i) = SELCTG(Alpha(i),Beta(i))
         ENDDO
!
         CALL CTGSEN(0,ilvsl,ilvsr,Bwork,N,A,Lda,B,Ldb,Alpha,Beta,Vsl,  &
     &               Ldvsl,Vsr,Ldvsr,Sdim,pvsl,pvsr,dif,Work(iwrk),     &
     &               Lwork-iwrk+1,idum,1,ierr)
         IF ( ierr==1 ) Info = N + 3
!
      ENDIF
!
!     Apply back-permutation to VSL and VSR
!
      IF ( ilvsl ) CALL CGGBAK('P','L',N,ilo,ihi,Rwork(ileft),          &
     &                         Rwork(iright),N,Vsl,Ldvsl,ierr)
      IF ( ilvsr ) CALL CGGBAK('P','R',N,ilo,ihi,Rwork(ileft),          &
     &                         Rwork(iright),N,Vsr,Ldvsr,ierr)
!
!     Undo scaling
!
      IF ( ilascl ) THEN
         CALL CLASCL('U',0,0,anrmto,anrm,N,N,A,Lda,ierr)
         CALL CLASCL('G',0,0,anrmto,anrm,N,1,Alpha,N,ierr)
      ENDIF
!
      IF ( ilbscl ) THEN
         CALL CLASCL('U',0,0,bnrmto,bnrm,N,N,B,Ldb,ierr)
         CALL CLASCL('G',0,0,bnrmto,bnrm,N,1,Beta,N,ierr)
      ENDIF
!
      IF ( wantst ) THEN
!
!        Check if reordering is correct
!
         lastsl = .TRUE.
         Sdim = 0
         DO i = 1 , N
            cursl = SELCTG(Alpha(i),Beta(i))
            IF ( cursl ) Sdim = Sdim + 1
            IF ( cursl .AND. .NOT.lastsl ) Info = N + 2
            lastsl = cursl
         ENDDO
!
      ENDIF
!
!
 100  Work(1) = CMPLX(lwkopt)
!
!
!     End of CGGES3
!
      END SUBROUTINE CGGES3
