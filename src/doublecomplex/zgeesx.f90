!*==zgeesx.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief <b> ZGEESX computes the eigenvalues, the Schur form, and, optionally, the matrix of Schur vectors for GE matrices</b>
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZGEESX + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zgeesx.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zgeesx.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zgeesx.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZGEESX( JOBVS, SORT, SELECT, SENSE, N, A, LDA, SDIM, W,
!                          VS, LDVS, RCONDE, RCONDV, WORK, LWORK, RWORK,
!                          BWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          JOBVS, SENSE, SORT
!       INTEGER            INFO, LDA, LDVS, LWORK, N, SDIM
!       DOUBLE PRECISION   RCONDE, RCONDV
!       ..
!       .. Array Arguments ..
!       LOGICAL            BWORK( * )
!       DOUBLE PRECISION   RWORK( * )
!       COMPLEX*16         A( LDA, * ), VS( LDVS, * ), W( * ), WORK( * )
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
!> ZGEESX computes for an N-by-N complex nonsymmetric matrix A, the
!> eigenvalues, the Schur form T, and, optionally, the matrix of Schur
!> vectors Z.  This gives the Schur factorization A = Z*T*(Z**H).
!>
!> Optionally, it also orders the eigenvalues on the diagonal of the
!> Schur form so that selected eigenvalues are at the top left;
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
!> A complex matrix is in Schur form if it is upper triangular.
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
!>          SELECT is a LOGICAL FUNCTION of one COMPLEX*16 argument
!>          SELECT must be declared EXTERNAL in the calling subroutine.
!>          If SORT = 'S', SELECT is used to select eigenvalues to order
!>          to the top left of the Schur form.
!>          If SORT = 'N', SELECT is not referenced.
!>          An eigenvalue W(j) is selected if SELECT(W(j)) is true.
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
!>          A is COMPLEX*16 array, dimension (LDA, N)
!>          On entry, the N-by-N matrix A.
!>          On exit, A is overwritten by its Schur form T.
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
!>          If SORT = 'S', SDIM = number of eigenvalues for which
!>                         SELECT is true.
!> \endverbatim
!>
!> \param[out] W
!> \verbatim
!>          W is COMPLEX*16 array, dimension (N)
!>          W contains the computed eigenvalues, in the same order
!>          that they appear on the diagonal of the output Schur form T.
!> \endverbatim
!>
!> \param[out] VS
!> \verbatim
!>          VS is COMPLEX*16 array, dimension (LDVS,N)
!>          If JOBVS = 'V', VS contains the unitary matrix Z of Schur
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
!>          RCONDE is DOUBLE PRECISION
!>          If SENSE = 'E' or 'B', RCONDE contains the reciprocal
!>          condition number for the average of the selected eigenvalues.
!>          Not referenced if SENSE = 'N' or 'V'.
!> \endverbatim
!>
!> \param[out] RCONDV
!> \verbatim
!>          RCONDV is DOUBLE PRECISION
!>          If SENSE = 'V' or 'B', RCONDV contains the reciprocal
!>          condition number for the selected right invariant subspace.
!>          Not referenced if SENSE = 'N' or 'E'.
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
!>          The dimension of the array WORK.  LWORK >= max(1,2*N).
!>          Also, if SENSE = 'E' or 'V' or 'B', LWORK >= 2*SDIM*(N-SDIM),
!>          where SDIM is the number of selected eigenvalues computed by
!>          this routine.  Note that 2*SDIM*(N-SDIM) <= N*N/2. Note also
!>          that an error is only returned if LWORK < max(1,2*N), but if
!>          SENSE = 'E' or 'V' or 'B' this may not be large enough.
!>          For good performance, LWORK must generally be larger.
!>
!>          If LWORK = -1, then a workspace query is assumed; the routine
!>          only calculates upper bound on the optimal size of the
!>          array WORK, returns this value as the first entry of the WORK
!>          array, and no error message related to LWORK is issued by
!>          XERBLA.
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is DOUBLE PRECISION array, dimension (N)
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
!>                   eigenvalues; elements 1:ILO-1 and i+1:N of W
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
!> \ingroup complex16GEeigen
!
!  =====================================================================
      SUBROUTINE ZGEESX(Jobvs,Sort,SELECT,Sense,N,A,Lda,Sdim,W,Vs,Ldvs, &
     &                  Rconde,Rcondv,Work,Lwork,Rwork,Bwork,Info)
      IMPLICIT NONE
!*--ZGEESX242
!
!  -- LAPACK driver routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     June 2016
!
!     .. Scalar Arguments ..
      CHARACTER Jobvs , Sense , Sort
      INTEGER Info , Lda , Ldvs , Lwork , N , Sdim
      DOUBLE PRECISION Rconde , Rcondv
!     ..
!     .. Array Arguments ..
      LOGICAL Bwork(*)
      DOUBLE PRECISION Rwork(*)
      COMPLEX*16 A(Lda,*) , Vs(Ldvs,*) , W(*) , Work(*)
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
      LOGICAL lquery , scalea , wantsb , wantse , wantsn , wantst ,     &
     &        wantsv , wantvs
      INTEGER hswork , i , ibal , icond , ierr , ieval , ihi , ilo ,    &
     &        itau , iwrk , lwrk , maxwrk , minwrk
      DOUBLE PRECISION anrm , bignum , cscale , eps , smlnum
!     ..
!     .. Local Arrays ..
      DOUBLE PRECISION dum(1)
!     ..
!     .. External Subroutines ..
      EXTERNAL DLABAD , DLASCL , XERBLA , ZCOPY , ZGEBAK , ZGEBAL ,     &
     &         ZGEHRD , ZHSEQR , ZLACPY , ZLASCL , ZTRSEN , ZUNGHR
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      INTEGER ILAENV
      DOUBLE PRECISION DLAMCH , ZLANGE
      EXTERNAL LSAME , ILAENV , DLAMCH , ZLANGE
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
      lquery = (Lwork==-1)
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
         Info = -11
      ENDIF
!
!     Compute workspace
!      (Note: Comments in the code beginning "Workspace:" describe the
!       minimal amount of real workspace needed at that point in the
!       code, as well as the preferred amount for good performance.
!       CWorkspace refers to complex workspace, and RWorkspace to real
!       workspace. NB refers to the optimal block size for the
!       immediately following subroutine, as returned by ILAENV.
!       HSWORK refers to the workspace preferred by ZHSEQR, as
!       calculated below. HSWORK is computed assuming ILO=1 and IHI=N,
!       the worst case.
!       If SENSE = 'E', 'V' or 'B', then the amount of workspace needed
!       depends on SDIM, which is computed by the routine ZTRSEN later
!       in the code.)
!
      IF ( Info==0 ) THEN
         IF ( N==0 ) THEN
            minwrk = 1
            lwrk = 1
         ELSE
            maxwrk = N + N*ILAENV(1,'ZGEHRD',' ',N,1,N,0)
            minwrk = 2*N
!
            CALL ZHSEQR('S',Jobvs,N,1,N,A,Lda,W,Vs,Ldvs,Work,-1,ieval)
            hswork = Work(1)
!
            IF ( .NOT.wantvs ) THEN
               maxwrk = MAX(maxwrk,hswork)
            ELSE
               maxwrk = MAX(maxwrk,N+(N-1)                              &
     &                  *ILAENV(1,'ZUNGHR',' ',N,1,N,-1))
               maxwrk = MAX(maxwrk,hswork)
            ENDIF
            lwrk = maxwrk
            IF ( .NOT.wantsn ) lwrk = MAX(lwrk,(N*N)/2)
         ENDIF
         Work(1) = lwrk
!
         IF ( Lwork<minwrk .AND. .NOT.lquery ) Info = -15
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('ZGEESX',-Info)
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
      anrm = ZLANGE('M',N,N,A,Lda,dum)
      scalea = .FALSE.
      IF ( anrm>ZERO .AND. anrm<smlnum ) THEN
         scalea = .TRUE.
         cscale = smlnum
      ELSEIF ( anrm>bignum ) THEN
         scalea = .TRUE.
         cscale = bignum
      ENDIF
      IF ( scalea ) CALL ZLASCL('G',0,0,anrm,cscale,N,N,A,Lda,ierr)
!
!
!     Permute the matrix to make it more nearly triangular
!     (CWorkspace: none)
!     (RWorkspace: need N)
!
      ibal = 1
      CALL ZGEBAL('P',N,A,Lda,ilo,ihi,Rwork(ibal),ierr)
!
!     Reduce to upper Hessenberg form
!     (CWorkspace: need 2*N, prefer N+N*NB)
!     (RWorkspace: none)
!
      itau = 1
      iwrk = N + itau
      CALL ZGEHRD(N,ilo,ihi,A,Lda,Work(itau),Work(iwrk),Lwork-iwrk+1,   &
     &            ierr)
!
      IF ( wantvs ) THEN
!
!        Copy Householder vectors to VS
!
         CALL ZLACPY('L',N,N,A,Lda,Vs,Ldvs)
!
!        Generate unitary matrix in VS
!        (CWorkspace: need 2*N-1, prefer N+(N-1)*NB)
!        (RWorkspace: none)
!
         CALL ZUNGHR(N,ilo,ihi,Vs,Ldvs,Work(itau),Work(iwrk),           &
     &               Lwork-iwrk+1,ierr)
      ENDIF
!
      Sdim = 0
!
!     Perform QR iteration, accumulating Schur vectors in VS if desired
!     (CWorkspace: need 1, prefer HSWORK (see comments) )
!     (RWorkspace: none)
!
      iwrk = itau
      CALL ZHSEQR('S',Jobvs,N,ilo,ihi,A,Lda,W,Vs,Ldvs,Work(iwrk),       &
     &            Lwork-iwrk+1,ieval)
      IF ( ieval>0 ) Info = ieval
!
!     Sort eigenvalues if desired
!
      IF ( wantst .AND. Info==0 ) THEN
         IF ( scalea ) CALL ZLASCL('G',0,0,cscale,anrm,N,1,W,N,ierr)
         DO i = 1 , N
            Bwork(i) = SELECT(W(i))
         ENDDO
!
!        Reorder eigenvalues, transform Schur vectors, and compute
!        reciprocal condition numbers
!        (CWorkspace: if SENSE is not 'N', need 2*SDIM*(N-SDIM)
!                     otherwise, need none )
!        (RWorkspace: none)
!
         CALL ZTRSEN(Sense,Jobvs,Bwork,N,A,Lda,Vs,Ldvs,W,Sdim,Rconde,   &
     &               Rcondv,Work(iwrk),Lwork-iwrk+1,icond)
         IF ( .NOT.wantsn ) maxwrk = MAX(maxwrk,2*Sdim*(N-Sdim))
!
!           Not enough complex workspace
!
         IF ( icond==-14 ) Info = -15
      ENDIF
!
!
!        Undo balancing
!        (CWorkspace: none)
!        (RWorkspace: need N)
!
      IF ( wantvs ) CALL ZGEBAK('P','R',N,ilo,ihi,Rwork(ibal),N,Vs,Ldvs,&
     &                          ierr)
!
      IF ( scalea ) THEN
!
!        Undo scaling for the Schur form of A
!
         CALL ZLASCL('U',0,0,cscale,anrm,N,N,A,Lda,ierr)
         CALL ZCOPY(N,A,Lda+1,W,1)
         IF ( (wantsv .OR. wantsb) .AND. Info==0 ) THEN
            dum(1) = Rcondv
            CALL DLASCL('G',0,0,cscale,anrm,1,1,dum,1,ierr)
            Rcondv = dum(1)
         ENDIF
      ENDIF
!
      Work(1) = maxwrk
!
!     End of ZGEESX
!
      END SUBROUTINE ZGEESX
