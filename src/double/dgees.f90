!*==dgees.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief <b> DGEES computes the eigenvalues, the Schur form, and, optionally, the matrix of Schur vectors for GE matrices</b>
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DGEES + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgees.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgees.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgees.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DGEES( JOBVS, SORT, SELECT, N, A, LDA, SDIM, WR, WI,
!                         VS, LDVS, WORK, LWORK, BWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          JOBVS, SORT
!       INTEGER            INFO, LDA, LDVS, LWORK, N, SDIM
!       ..
!       .. Array Arguments ..
!       LOGICAL            BWORK( * )
!       DOUBLE PRECISION   A( LDA, * ), VS( LDVS, * ), WI( * ), WORK( * ),
!      $                   WR( * )
!       ..
!       .. Function Arguments ..
!       LOGICAL            SELECT
!       EXTERNAL           SELECT
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DGEES computes for an N-by-N real nonsymmetric matrix A, the
!> eigenvalues, the real Schur form T, and, optionally, the matrix of
!> Schur vectors Z.  This gives the Schur factorization A = Z*T*(Z**T).
!>
!> Optionally, it also orders the eigenvalues on the diagonal of the
!> real Schur form so that selected eigenvalues are at the top left.
!> The leading columns of Z then form an orthonormal basis for the
!> invariant subspace corresponding to the selected eigenvalues.
!>
!> A matrix is in real Schur form if it is upper quasi-triangular with
!> 1-by-1 and 2-by-2 blocks. 2-by-2 blocks will be standardized in the
!> form
!>         [  a  b  ]
!>         [  c  a  ]
!>
!> where b*c < 0. The eigenvalues of such a block are a +- sqrt(bc).
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] JOBVS
!> \verbatim
!>          JOBVS is CHARACTER*1
!>          = 'N': Schur vectors are not computed;
!>          = 'V': Schur vectors are computed.
!> \endverbatim
!>
!> \param[in] SORT
!> \verbatim
!>          SORT is CHARACTER*1
!>          Specifies whether or not to order the eigenvalues on the
!>          diagonal of the Schur form.
!>          = 'N': Eigenvalues are not ordered;
!>          = 'S': Eigenvalues are ordered (see SELECT).
!> \endverbatim
!>
!> \param[in] SELECT
!> \verbatim
!>          SELECT is a LOGICAL FUNCTION of two DOUBLE PRECISION arguments
!>          SELECT must be declared EXTERNAL in the calling subroutine.
!>          If SORT = 'S', SELECT is used to select eigenvalues to sort
!>          to the top left of the Schur form.
!>          If SORT = 'N', SELECT is not referenced.
!>          An eigenvalue WR(j)+sqrt(-1)*WI(j) is selected if
!>          SELECT(WR(j),WI(j)) is true; i.e., if either one of a complex
!>          conjugate pair of eigenvalues is selected, then both complex
!>          eigenvalues are selected.
!>          Note that a selected complex eigenvalue may no longer
!>          satisfy SELECT(WR(j),WI(j)) = .TRUE. after ordering, since
!>          ordering may change the value of complex eigenvalues
!>          (especially if the eigenvalue is ill-conditioned); in this
!>          case INFO is set to N+2 (see INFO below).
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
!>          On exit, A has been overwritten by its real Schur form T.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,N).
!> \endverbatim
!>
!> \param[out] SDIM
!> \verbatim
!>          SDIM is INTEGER
!>          If SORT = 'N', SDIM = 0.
!>          If SORT = 'S', SDIM = number of eigenvalues (after sorting)
!>                         for which SELECT is true. (Complex conjugate
!>                         pairs for which SELECT is true for either
!>                         eigenvalue count as 2.)
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
!>          respectively, of the computed eigenvalues in the same order
!>          that they appear on the diagonal of the output Schur form T.
!>          Complex conjugate pairs of eigenvalues will appear
!>          consecutively with the eigenvalue having the positive
!>          imaginary part first.
!> \endverbatim
!>
!> \param[out] VS
!> \verbatim
!>          VS is DOUBLE PRECISION array, dimension (LDVS,N)
!>          If JOBVS = 'V', VS contains the orthogonal matrix Z of Schur
!>          vectors.
!>          If JOBVS = 'N', VS is not referenced.
!> \endverbatim
!>
!> \param[in] LDVS
!> \verbatim
!>          LDVS is INTEGER
!>          The leading dimension of the array VS.  LDVS >= 1; if
!>          JOBVS = 'V', LDVS >= N.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))
!>          On exit, if INFO = 0, WORK(1) contains the optimal LWORK.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The dimension of the array WORK.  LWORK >= max(1,3*N).
!>          For good performance, LWORK must generally be larger.
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
!>          = 0: successful exit
!>          < 0: if INFO = -i, the i-th argument had an illegal value.
!>          > 0: if INFO = i, and i is
!>             <= N: the QR algorithm failed to compute all the
!>                   eigenvalues; elements 1:ILO-1 and i+1:N of WR and WI
!>                   contain those eigenvalues which have converged; if
!>                   JOBVS = 'V', VS contains the matrix which reduces A
!>                   to its partially converged Schur form.
!>             = N+1: the eigenvalues could not be reordered because some
!>                   eigenvalues were too close to separate (the problem
!>                   is very ill-conditioned);
!>             = N+2: after reordering, roundoff changed values of some
!>                   complex eigenvalues so that leading eigenvalues in
!>                   the Schur form no longer satisfy SELECT=.TRUE.  This
!>                   could also be caused by underflow due to scaling.
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
!  =====================================================================
      SUBROUTINE DGEES(Jobvs,Sort,SELECT,N,A,Lda,Sdim,Wr,Wi,Vs,Ldvs,    &
     &                 Work,Lwork,Bwork,Info)
      IMPLICIT NONE
!*--DGEES220
!
!  -- LAPACK driver routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Jobvs , Sort
      INTEGER Info , Lda , Ldvs , Lwork , N , Sdim
!     ..
!     .. Array Arguments ..
      LOGICAL Bwork(*)
      DOUBLE PRECISION A(Lda,*) , Vs(Ldvs,*) , Wi(*) , Work(*) , Wr(*)
!     ..
!     .. Function Arguments ..
      LOGICAL SELECT
      EXTERNAL SELECT
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO , ONE
      PARAMETER (ZERO=0.0D0,ONE=1.0D0)
!     ..
!     .. Local Scalars ..
      LOGICAL cursl , lastsl , lquery , lst2sl , scalea , wantst ,      &
     &        wantvs
      INTEGER hswork , i , i1 , i2 , ibal , icond , ierr , ieval , ihi ,&
     &        ilo , inxt , ip , itau , iwrk , maxwrk , minwrk
      DOUBLE PRECISION anrm , bignum , cscale , eps , s , sep , smlnum
!     ..
!     .. Local Arrays ..
      INTEGER idum(1)
      DOUBLE PRECISION dum(1)
!     ..
!     .. External Subroutines ..
      EXTERNAL DCOPY , DGEBAK , DGEBAL , DGEHRD , DHSEQR , DLACPY ,     &
     &         DLABAD , DLASCL , DORGHR , DSWAP , DTRSEN , XERBLA
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      INTEGER ILAENV
      DOUBLE PRECISION DLAMCH , DLANGE
      EXTERNAL LSAME , ILAENV , DLAMCH , DLANGE
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX , SQRT
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      Info = 0
      lquery = (Lwork==-1)
      wantvs = LSAME(Jobvs,'V')
      wantst = LSAME(Sort,'S')
      IF ( (.NOT.wantvs) .AND. (.NOT.LSAME(Jobvs,'N')) ) THEN
         Info = -1
      ELSEIF ( (.NOT.wantst) .AND. (.NOT.LSAME(Sort,'N')) ) THEN
         Info = -2
      ELSEIF ( N<0 ) THEN
         Info = -4
      ELSEIF ( Lda<MAX(1,N) ) THEN
         Info = -6
      ELSEIF ( Ldvs<1 .OR. (wantvs .AND. Ldvs<N) ) THEN
         Info = -11
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
            maxwrk = 2*N + N*ILAENV(1,'DGEHRD',' ',N,1,N,0)
            minwrk = 3*N
!
            CALL DHSEQR('S',Jobvs,N,1,N,A,Lda,Wr,Wi,Vs,Ldvs,Work,-1,    &
     &                  ieval)
            hswork = Work(1)
!
            IF ( .NOT.wantvs ) THEN
               maxwrk = MAX(maxwrk,N+hswork)
            ELSE
               maxwrk = MAX(maxwrk,2*N+(N-1)                            &
     &                  *ILAENV(1,'DORGHR',' ',N,1,N,-1))
               maxwrk = MAX(maxwrk,N+hswork)
            ENDIF
         ENDIF
         Work(1) = maxwrk
!
         IF ( Lwork<minwrk .AND. .NOT.lquery ) Info = -13
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('DGEES ',-Info)
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
      smlnum = DLAMCH('S')
      bignum = ONE/smlnum
      CALL DLABAD(smlnum,bignum)
      smlnum = SQRT(smlnum)/eps
      bignum = ONE/smlnum
!
!     Scale A if max element outside range [SMLNUM,BIGNUM]
!
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
!     Permute the matrix to make it more nearly triangular
!     (Workspace: need N)
!
      ibal = 1
      CALL DGEBAL('P',N,A,Lda,ilo,ihi,Work(ibal),ierr)
!
!     Reduce to upper Hessenberg form
!     (Workspace: need 3*N, prefer 2*N+N*NB)
!
      itau = N + ibal
      iwrk = N + itau
      CALL DGEHRD(N,ilo,ihi,A,Lda,Work(itau),Work(iwrk),Lwork-iwrk+1,   &
     &            ierr)
!
      IF ( wantvs ) THEN
!
!        Copy Householder vectors to VS
!
         CALL DLACPY('L',N,N,A,Lda,Vs,Ldvs)
!
!        Generate orthogonal matrix in VS
!        (Workspace: need 3*N-1, prefer 2*N+(N-1)*NB)
!
         CALL DORGHR(N,ilo,ihi,Vs,Ldvs,Work(itau),Work(iwrk),           &
     &               Lwork-iwrk+1,ierr)
      ENDIF
!
      Sdim = 0
!
!     Perform QR iteration, accumulating Schur vectors in VS if desired
!     (Workspace: need N+1, prefer N+HSWORK (see comments) )
!
      iwrk = itau
      CALL DHSEQR('S',Jobvs,N,ilo,ihi,A,Lda,Wr,Wi,Vs,Ldvs,Work(iwrk),   &
     &            Lwork-iwrk+1,ieval)
      IF ( ieval>0 ) Info = ieval
!
!     Sort eigenvalues if desired
!
      IF ( wantst .AND. Info==0 ) THEN
         IF ( scalea ) THEN
            CALL DLASCL('G',0,0,cscale,anrm,N,1,Wr,N,ierr)
            CALL DLASCL('G',0,0,cscale,anrm,N,1,Wi,N,ierr)
         ENDIF
         DO i = 1 , N
            Bwork(i) = SELECT(Wr(i),Wi(i))
         ENDDO
!
!        Reorder eigenvalues and transform Schur vectors
!        (Workspace: none needed)
!
         CALL DTRSEN('N',Jobvs,Bwork,N,A,Lda,Vs,Ldvs,Wr,Wi,Sdim,s,sep,  &
     &               Work(iwrk),Lwork-iwrk+1,idum,1,icond)
         IF ( icond>0 ) Info = N + icond
      ENDIF
!
!
!        Undo balancing
!        (Workspace: need N)
!
      IF ( wantvs ) CALL DGEBAK('P','R',N,ilo,ihi,Work(ibal),N,Vs,Ldvs, &
     &                          ierr)
!
      IF ( scalea ) THEN
!
!        Undo scaling for the Schur form of A
!
         CALL DLASCL('H',0,0,cscale,anrm,N,N,A,Lda,ierr)
         CALL DCOPY(N,A,Lda+1,Wr,1)
         IF ( cscale==smlnum ) THEN
!
!           If scaling back towards underflow, adjust WI if an
!           offdiagonal element of a 2-by-2 block in the Schur form
!           underflows.
!
            IF ( ieval>0 ) THEN
               i1 = ieval + 1
               i2 = ihi - 1
               CALL DLASCL('G',0,0,cscale,anrm,ilo-1,1,Wi,MAX(ilo-1,1), &
     &                     ierr)
            ELSEIF ( wantst ) THEN
               i1 = 1
               i2 = N - 1
            ELSE
               i1 = ilo
               i2 = ihi - 1
            ENDIF
            inxt = i1 - 1
            DO i = i1 , i2
               IF ( i>=inxt ) THEN
                  IF ( Wi(i)==ZERO ) THEN
                     inxt = i + 1
                  ELSE
                     IF ( A(i+1,i)==ZERO ) THEN
                        Wi(i) = ZERO
                        Wi(i+1) = ZERO
                     ELSEIF ( A(i+1,i)/=ZERO .AND. A(i,i+1)==ZERO ) THEN
                        Wi(i) = ZERO
                        Wi(i+1) = ZERO
                        IF ( i>1 ) CALL DSWAP(i-1,A(1,i),1,A(1,i+1),1)
                        IF ( N>i+1 )                                    &
     &                       CALL DSWAP(N-i-1,A(i,i+2),Lda,A(i+1,i+2),  &
     &                       Lda)
                        IF ( wantvs )                                   &
     &                       CALL DSWAP(N,Vs(1,i),1,Vs(1,i+1),1)
                        A(i,i+1) = A(i+1,i)
                        A(i+1,i) = ZERO
                     ENDIF
                     inxt = i + 2
                  ENDIF
               ENDIF
            ENDDO
         ENDIF
!
!        Undo scaling for the imaginary part of the eigenvalues
!
         CALL DLASCL('G',0,0,cscale,anrm,N-ieval,1,Wi(ieval+1),         &
     &               MAX(N-ieval,1),ierr)
      ENDIF
!
      IF ( wantst .AND. Info==0 ) THEN
!
!        Check if reordering successful
!
         lastsl = .TRUE.
         lst2sl = .TRUE.
         Sdim = 0
         ip = 0
         DO i = 1 , N
            cursl = SELECT(Wr(i),Wi(i))
            IF ( Wi(i)==ZERO ) THEN
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
      ENDIF
!
      Work(1) = maxwrk
!
!     End of DGEES
!
      END SUBROUTINE DGEES
