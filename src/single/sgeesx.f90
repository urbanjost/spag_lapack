!*==sgeesx.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief <b> SGEESX computes the eigenvalues, the Schur form, and, optionally, the matrix of Schur vectors for GE matrices</b>
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SGEESX + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgeesx.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgeesx.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgeesx.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SGEESX( JOBVS, SORT, SELECT, SENSE, N, A, LDA, SDIM,
!                          WR, WI, VS, LDVS, RCONDE, RCONDV, WORK, LWORK,
!                          IWORK, LIWORK, BWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          JOBVS, SENSE, SORT
!       INTEGER            INFO, LDA, LDVS, LIWORK, LWORK, N, SDIM
!       REAL               RCONDE, RCONDV
!       ..
!       .. Array Arguments ..
!       LOGICAL            BWORK( * )
!       INTEGER            IWORK( * )
!       REAL               A( LDA, * ), VS( LDVS, * ), WI( * ), WORK( * ),
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
!> SGEESX computes for an N-by-N real nonsymmetric matrix A, the
!> eigenvalues, the real Schur form T, and, optionally, the matrix of
!> Schur vectors Z.  This gives the Schur factorization A = Z*T*(Z**T).
!>
!> Optionally, it also orders the eigenvalues on the diagonal of the
!> real Schur form so that selected eigenvalues are at the top left;
!> computes a reciprocal condition number for the average of the
!> selected eigenvalues (RCONDE); and computes a reciprocal condition
!> number for the right invariant subspace corresponding to the
!> selected eigenvalues (RCONDV).  The leading columns of Z form an
!> orthonormal basis for this invariant subspace.
!>
!> For further explanation of the reciprocal condition numbers RCONDE
!> and RCONDV, see Section 4.10 of the LAPACK Users' Guide (where
!> these quantities are called s and sep respectively).
!>
!> A real matrix is in real Schur form if it is upper quasi-triangular
!> with 1-by-1 and 2-by-2 blocks. 2-by-2 blocks will be standardized in
!> the form
!>           [  a  b  ]
!>           [  c  a  ]
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
!>          SELECT is a LOGICAL FUNCTION of two REAL arguments
!>          SELECT must be declared EXTERNAL in the calling subroutine.
!>          If SORT = 'S', SELECT is used to select eigenvalues to sort
!>          to the top left of the Schur form.
!>          If SORT = 'N', SELECT is not referenced.
!>          An eigenvalue WR(j)+sqrt(-1)*WI(j) is selected if
!>          SELECT(WR(j),WI(j)) is true; i.e., if either one of a
!>          complex conjugate pair of eigenvalues is selected, then both
!>          are.  Note that a selected complex eigenvalue may no longer
!>          satisfy SELECT(WR(j),WI(j)) = .TRUE. after ordering, since
!>          ordering may change the value of complex eigenvalues
!>          (especially if the eigenvalue is ill-conditioned); in this
!>          case INFO may be set to N+3 (see INFO below).
!> \endverbatim
!>
!> \param[in] SENSE
!> \verbatim
!>          SENSE is CHARACTER*1
!>          Determines which reciprocal condition numbers are computed.
!>          = 'N': None are computed;
!>          = 'E': Computed for average of selected eigenvalues only;
!>          = 'V': Computed for selected right invariant subspace only;
!>          = 'B': Computed for both.
!>          If SENSE = 'E', 'V' or 'B', SORT must equal 'S'.
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
!>          A is REAL array, dimension (LDA, N)
!>          On entry, the N-by-N matrix A.
!>          On exit, A is overwritten by its real Schur form T.
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
!>          WR is REAL array, dimension (N)
!> \endverbatim
!>
!> \param[out] WI
!> \verbatim
!>          WI is REAL array, dimension (N)
!>          WR and WI contain the real and imaginary parts, respectively,
!>          of the computed eigenvalues, in the same order that they
!>          appear on the diagonal of the output Schur form T.  Complex
!>          conjugate pairs of eigenvalues appear consecutively with the
!>          eigenvalue having the positive imaginary part first.
!> \endverbatim
!>
!> \param[out] VS
!> \verbatim
!>          VS is REAL array, dimension (LDVS,N)
!>          If JOBVS = 'V', VS contains the orthogonal matrix Z of Schur
!>          vectors.
!>          If JOBVS = 'N', VS is not referenced.
!> \endverbatim
!>
!> \param[in] LDVS
!> \verbatim
!>          LDVS is INTEGER
!>          The leading dimension of the array VS.  LDVS >= 1, and if
!>          JOBVS = 'V', LDVS >= N.
!> \endverbatim
!>
!> \param[out] RCONDE
!> \verbatim
!>          RCONDE is REAL
!>          If SENSE = 'E' or 'B', RCONDE contains the reciprocal
!>          condition number for the average of the selected eigenvalues.
!>          Not referenced if SENSE = 'N' or 'V'.
!> \endverbatim
!>
!> \param[out] RCONDV
!> \verbatim
!>          RCONDV is REAL
!>          If SENSE = 'V' or 'B', RCONDV contains the reciprocal
!>          condition number for the selected right invariant subspace.
!>          Not referenced if SENSE = 'N' or 'E'.
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
!>          The dimension of the array WORK.  LWORK >= max(1,3*N).
!>          Also, if SENSE = 'E' or 'V' or 'B',
!>          LWORK >= N+2*SDIM*(N-SDIM), where SDIM is the number of
!>          selected eigenvalues computed by this routine.  Note that
!>          N+2*SDIM*(N-SDIM) <= N+N*N/2. Note also that an error is only
!>          returned if LWORK < max(1,3*N), but if SENSE = 'E' or 'V' or
!>          'B' this may not be large enough.
!>          For good performance, LWORK must generally be larger.
!>
!>          If LWORK = -1, then a workspace query is assumed; the routine
!>          only calculates upper bounds on the optimal sizes of the
!>          arrays WORK and IWORK, returns these values as the first
!>          entries of the WORK and IWORK arrays, and no error messages
!>          related to LWORK or LIWORK are issued by XERBLA.
!> \endverbatim
!>
!> \param[out] IWORK
!> \verbatim
!>          IWORK is INTEGER array, dimension (MAX(1,LIWORK))
!>          On exit, if INFO = 0, IWORK(1) returns the optimal LIWORK.
!> \endverbatim
!>
!> \param[in] LIWORK
!> \verbatim
!>          LIWORK is INTEGER
!>          The dimension of the array IWORK.
!>          LIWORK >= 1; if SENSE = 'V' or 'B', LIWORK >= SDIM*(N-SDIM).
!>          Note that SDIM*(N-SDIM) <= N*N/4. Note also that an error is
!>          only returned if LIWORK < 1, but if SENSE = 'V' or 'B' this
!>          may not be large enough.
!>
!>          If LIWORK = -1, then a workspace query is assumed; the
!>          routine only calculates upper bounds on the optimal sizes of
!>          the arrays WORK and IWORK, returns these values as the first
!>          entries of the WORK and IWORK arrays, and no error messages
!>          related to LWORK or LIWORK are issued by XERBLA.
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
!>                   JOBVS = 'V', VS contains the transformation which
!>                   reduces A to its partially converged Schur form.
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
!> \date June 2016
!
!> \ingroup realGEeigen
!
!  =====================================================================
      SUBROUTINE SGEESX(Jobvs,Sort,SELECT,Sense,N,A,Lda,Sdim,Wr,Wi,Vs,  &
     &                  Ldvs,Rconde,Rcondv,Work,Lwork,Iwork,Liwork,     &
     &                  Bwork,Info)
      IMPLICIT NONE
!*--SGEESX285
!
!  -- LAPACK driver routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     June 2016
!
!     .. Scalar Arguments ..
      CHARACTER Jobvs , Sense , Sort
      INTEGER Info , Lda , Ldvs , Liwork , Lwork , N , Sdim
      REAL Rconde , Rcondv
!     ..
!     .. Array Arguments ..
      LOGICAL Bwork(*)
      INTEGER Iwork(*)
      REAL A(Lda,*) , Vs(Ldvs,*) , Wi(*) , Work(*) , Wr(*)
!     ..
!     .. Function Arguments ..
      LOGICAL SELECT
      EXTERNAL SELECT
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ZERO , ONE
      PARAMETER (ZERO=0.0E0,ONE=1.0E0)
!     ..
!     .. Local Scalars ..
      LOGICAL cursl , lastsl , lquery , lst2sl , scalea , wantsb ,      &
     &        wantse , wantsn , wantst , wantsv , wantvs
      INTEGER hswork , i , i1 , i2 , ibal , icond , ierr , ieval , ihi ,&
     &        ilo , inxt , ip , itau , iwrk , lwrk , liwrk , maxwrk ,   &
     &        minwrk
      REAL anrm , bignum , cscale , eps , smlnum
!     ..
!     .. Local Arrays ..
      REAL dum(1)
!     ..
!     .. External Subroutines ..
      EXTERNAL SCOPY , SGEBAK , SGEBAL , SGEHRD , SHSEQR , SLABAD ,     &
     &         SLACPY , SLASCL , SORGHR , SSWAP , STRSEN , XERBLA
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      INTEGER ILAENV
      REAL SLAMCH , SLANGE
      EXTERNAL LSAME , ILAENV , SLAMCH , SLANGE
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX , SQRT
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      Info = 0
      wantvs = LSAME(Jobvs,'V')
      wantst = LSAME(Sort,'S')
      wantsn = LSAME(Sense,'N')
      wantse = LSAME(Sense,'E')
      wantsv = LSAME(Sense,'V')
      wantsb = LSAME(Sense,'B')
      lquery = (Lwork==-1 .OR. Liwork==-1)
!
      IF ( (.NOT.wantvs) .AND. (.NOT.LSAME(Jobvs,'N')) ) THEN
         Info = -1
      ELSEIF ( (.NOT.wantst) .AND. (.NOT.LSAME(Sort,'N')) ) THEN
         Info = -2
      ELSEIF ( .NOT.(wantsn .OR. wantse .OR. wantsv .OR. wantsb) .OR.   &
     &         (.NOT.wantst .AND. .NOT.wantsn) ) THEN
         Info = -4
      ELSEIF ( N<0 ) THEN
         Info = -5
      ELSEIF ( Lda<MAX(1,N) ) THEN
         Info = -7
      ELSEIF ( Ldvs<1 .OR. (wantvs .AND. Ldvs<N) ) THEN
         Info = -12
      ENDIF
!
!     Compute workspace
!      (Note: Comments in the code beginning "RWorkspace:" describe the
!       minimal amount of real workspace needed at that point in the
!       code, as well as the preferred amount for good performance.
!       IWorkspace refers to integer workspace.
!       NB refers to the optimal block size for the immediately
!       following subroutine, as returned by ILAENV.
!       HSWORK refers to the workspace preferred by SHSEQR, as
!       calculated below. HSWORK is computed assuming ILO=1 and IHI=N,
!       the worst case.
!       If SENSE = 'E', 'V' or 'B', then the amount of workspace needed
!       depends on SDIM, which is computed by the routine STRSEN later
!       in the code.)
!
      IF ( Info==0 ) THEN
         liwrk = 1
         IF ( N==0 ) THEN
            minwrk = 1
            lwrk = 1
         ELSE
            maxwrk = 2*N + N*ILAENV(1,'SGEHRD',' ',N,1,N,0)
            minwrk = 3*N
!
            CALL SHSEQR('S',Jobvs,N,1,N,A,Lda,Wr,Wi,Vs,Ldvs,Work,-1,    &
     &                  ieval)
            hswork = Work(1)
!
            IF ( .NOT.wantvs ) THEN
               maxwrk = MAX(maxwrk,N+hswork)
            ELSE
               maxwrk = MAX(maxwrk,2*N+(N-1)                            &
     &                  *ILAENV(1,'SORGHR',' ',N,1,N,-1))
               maxwrk = MAX(maxwrk,N+hswork)
            ENDIF
            lwrk = maxwrk
            IF ( .NOT.wantsn ) lwrk = MAX(lwrk,N+(N*N)/2)
            IF ( wantsv .OR. wantsb ) liwrk = (N*N)/4
         ENDIF
         Iwork(1) = liwrk
         Work(1) = lwrk
!
         IF ( Lwork<minwrk .AND. .NOT.lquery ) THEN
            Info = -16
         ELSEIF ( Liwork<1 .AND. .NOT.lquery ) THEN
            Info = -18
         ENDIF
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('SGEESX',-Info)
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
      anrm = SLANGE('M',N,N,A,Lda,dum)
      scalea = .FALSE.
      IF ( anrm>ZERO .AND. anrm<smlnum ) THEN
         scalea = .TRUE.
         cscale = smlnum
      ELSEIF ( anrm>bignum ) THEN
         scalea = .TRUE.
         cscale = bignum
      ENDIF
      IF ( scalea ) CALL SLASCL('G',0,0,anrm,cscale,N,N,A,Lda,ierr)
!
!     Permute the matrix to make it more nearly triangular
!     (RWorkspace: need N)
!
      ibal = 1
      CALL SGEBAL('P',N,A,Lda,ilo,ihi,Work(ibal),ierr)
!
!     Reduce to upper Hessenberg form
!     (RWorkspace: need 3*N, prefer 2*N+N*NB)
!
      itau = N + ibal
      iwrk = N + itau
      CALL SGEHRD(N,ilo,ihi,A,Lda,Work(itau),Work(iwrk),Lwork-iwrk+1,   &
     &            ierr)
!
      IF ( wantvs ) THEN
!
!        Copy Householder vectors to VS
!
         CALL SLACPY('L',N,N,A,Lda,Vs,Ldvs)
!
!        Generate orthogonal matrix in VS
!        (RWorkspace: need 3*N-1, prefer 2*N+(N-1)*NB)
!
         CALL SORGHR(N,ilo,ihi,Vs,Ldvs,Work(itau),Work(iwrk),           &
     &               Lwork-iwrk+1,ierr)
      ENDIF
!
      Sdim = 0
!
!     Perform QR iteration, accumulating Schur vectors in VS if desired
!     (RWorkspace: need N+1, prefer N+HSWORK (see comments) )
!
      iwrk = itau
      CALL SHSEQR('S',Jobvs,N,ilo,ihi,A,Lda,Wr,Wi,Vs,Ldvs,Work(iwrk),   &
     &            Lwork-iwrk+1,ieval)
      IF ( ieval>0 ) Info = ieval
!
!     Sort eigenvalues if desired
!
      IF ( wantst .AND. Info==0 ) THEN
         IF ( scalea ) THEN
            CALL SLASCL('G',0,0,cscale,anrm,N,1,Wr,N,ierr)
            CALL SLASCL('G',0,0,cscale,anrm,N,1,Wi,N,ierr)
         ENDIF
         DO i = 1 , N
            Bwork(i) = SELECT(Wr(i),Wi(i))
         ENDDO
!
!        Reorder eigenvalues, transform Schur vectors, and compute
!        reciprocal condition numbers
!        (RWorkspace: if SENSE is not 'N', need N+2*SDIM*(N-SDIM)
!                     otherwise, need N )
!        (IWorkspace: if SENSE is 'V' or 'B', need SDIM*(N-SDIM)
!                     otherwise, need 0 )
!
         CALL STRSEN(Sense,Jobvs,Bwork,N,A,Lda,Vs,Ldvs,Wr,Wi,Sdim,      &
     &               Rconde,Rcondv,Work(iwrk),Lwork-iwrk+1,Iwork,Liwork,&
     &               icond)
         IF ( .NOT.wantsn ) maxwrk = MAX(maxwrk,N+2*Sdim*(N-Sdim))
         IF ( icond==-15 ) THEN
!
!           Not enough real workspace
!
            Info = -16
         ELSEIF ( icond==-17 ) THEN
!
!           Not enough integer workspace
!
            Info = -18
         ELSEIF ( icond>0 ) THEN
!
!           STRSEN failed to reorder or to restore standard Schur form
!
            Info = icond + N
         ENDIF
      ENDIF
!
!
!        Undo balancing
!        (RWorkspace: need N)
!
      IF ( wantvs ) CALL SGEBAK('P','R',N,ilo,ihi,Work(ibal),N,Vs,Ldvs, &
     &                          ierr)
!
      IF ( scalea ) THEN
!
!        Undo scaling for the Schur form of A
!
         CALL SLASCL('H',0,0,cscale,anrm,N,N,A,Lda,ierr)
         CALL SCOPY(N,A,Lda+1,Wr,1)
         IF ( (wantsv .OR. wantsb) .AND. Info==0 ) THEN
            dum(1) = Rcondv
            CALL SLASCL('G',0,0,cscale,anrm,1,1,dum,1,ierr)
            Rcondv = dum(1)
         ENDIF
         IF ( cscale==smlnum ) THEN
!
!           If scaling back towards underflow, adjust WI if an
!           offdiagonal element of a 2-by-2 block in the Schur form
!           underflows.
!
            IF ( ieval>0 ) THEN
               i1 = ieval + 1
               i2 = ihi - 1
               CALL SLASCL('G',0,0,cscale,anrm,ilo-1,1,Wi,N,ierr)
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
                        IF ( i>1 ) CALL SSWAP(i-1,A(1,i),1,A(1,i+1),1)
                        IF ( N>i+1 )                                    &
     &                       CALL SSWAP(N-i-1,A(i,i+2),Lda,A(i+1,i+2),  &
     &                       Lda)
                        IF ( wantvs )                                   &
     &                       CALL SSWAP(N,Vs(1,i),1,Vs(1,i+1),1)
                        A(i,i+1) = A(i+1,i)
                        A(i+1,i) = ZERO
                     ENDIF
                     inxt = i + 2
                  ENDIF
               ENDIF
            ENDDO
         ENDIF
         CALL SLASCL('G',0,0,cscale,anrm,N-ieval,1,Wi(ieval+1),         &
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
      IF ( wantsv .OR. wantsb ) THEN
         Iwork(1) = Sdim*(N-Sdim)
      ELSE
         Iwork(1) = 1
      ENDIF
!
!
!     End of SGEESX
!
      END SUBROUTINE SGEESX
