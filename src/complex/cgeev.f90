!*==cgeev.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief <b> CGEEV computes the eigenvalues and, optionally, the left and/or right eigenvectors for GE matrices</b>
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CGEEV + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgeev.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgeev.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgeev.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CGEEV( JOBVL, JOBVR, N, A, LDA, W, VL, LDVL, VR, LDVR,
!                         WORK, LWORK, RWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          JOBVL, JOBVR
!       INTEGER            INFO, LDA, LDVL, LDVR, LWORK, N
!       ..
!       .. Array Arguments ..
!       REAL   RWORK( * )
!       COMPLEX         A( LDA, * ), VL( LDVL, * ), VR( LDVR, * ),
!      $                   W( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CGEEV computes for an N-by-N complex nonsymmetric matrix A, the
!> eigenvalues and, optionally, the left and/or right eigenvectors.
!>
!> The right eigenvector v(j) of A satisfies
!>                  A * v(j) = lambda(j) * v(j)
!> where lambda(j) is its eigenvalue.
!> The left eigenvector u(j) of A satisfies
!>               u(j)**H * A = lambda(j) * u(j)**H
!> where u(j)**H denotes the conjugate transpose of u(j).
!>
!> The computed eigenvectors are normalized to have Euclidean norm
!> equal to 1 and largest component real.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] JOBVL
!> \verbatim
!>          JOBVL is CHARACTER*1
!>          = 'N': left eigenvectors of A are not computed;
!>          = 'V': left eigenvectors of are computed.
!> \endverbatim
!>
!> \param[in] JOBVR
!> \verbatim
!>          JOBVR is CHARACTER*1
!>          = 'N': right eigenvectors of A are not computed;
!>          = 'V': right eigenvectors of A are computed.
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
!>          A is COMPLEX array, dimension (LDA,N)
!>          On entry, the N-by-N matrix A.
!>          On exit, A has been overwritten.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,N).
!> \endverbatim
!>
!> \param[out] W
!> \verbatim
!>          W is COMPLEX array, dimension (N)
!>          W contains the computed eigenvalues.
!> \endverbatim
!>
!> \param[out] VL
!> \verbatim
!>          VL is COMPLEX array, dimension (LDVL,N)
!>          If JOBVL = 'V', the left eigenvectors u(j) are stored one
!>          after another in the columns of VL, in the same order
!>          as their eigenvalues.
!>          If JOBVL = 'N', VL is not referenced.
!>          u(j) = VL(:,j), the j-th column of VL.
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
!>          VR is COMPLEX array, dimension (LDVR,N)
!>          If JOBVR = 'V', the right eigenvectors v(j) are stored one
!>          after another in the columns of VR, in the same order
!>          as their eigenvalues.
!>          If JOBVR = 'N', VR is not referenced.
!>          v(j) = VR(:,j), the j-th column of VR.
!> \endverbatim
!>
!> \param[in] LDVR
!> \verbatim
!>          LDVR is INTEGER
!>          The leading dimension of the array VR.  LDVR >= 1; if
!>          JOBVR = 'V', LDVR >= N.
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
!>          The dimension of the array WORK.  LWORK >= max(1,2*N).
!>          For good performance, LWORK must generally be larger.
!>
!>          If LWORK = -1, then a workspace query is assumed; the routine
!>          only calculates the optimal size of the WORK array, returns
!>          this value as the first entry of the WORK array, and no error
!>          message related to LWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is REAL array, dimension (2*N)
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value.
!>          > 0:  if INFO = i, the QR algorithm failed to compute all the
!>                eigenvalues, and no eigenvectors have been computed;
!>                elements i+1:N of W contain eigenvalues which have
!>                converged.
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
!  @generated from zgeev.f, fortran z -> c, Tue Apr 19 01:47:44 2016
!
!> \ingroup complexGEeigen
!
!  =====================================================================
      SUBROUTINE CGEEV(Jobvl,Jobvr,N,A,Lda,W,Vl,Ldvl,Vr,Ldvr,Work,Lwork,&
     &                 Rwork,Info)
      IMPLICIT NONE
!*--CGEEV183
!
!  -- LAPACK driver routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     June 2016
!
!     .. Scalar Arguments ..
      CHARACTER Jobvl , Jobvr
      INTEGER Info , Lda , Ldvl , Ldvr , Lwork , N
!     ..
!     .. Array Arguments ..
      REAL Rwork(*)
      COMPLEX A(Lda,*) , Vl(Ldvl,*) , Vr(Ldvr,*) , W(*) , Work(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ZERO , ONE
      PARAMETER (ZERO=0.0E0,ONE=1.0E0)
!     ..
!     .. Local Scalars ..
      LOGICAL lquery , scalea , wantvl , wantvr
      CHARACTER side
      INTEGER hswork , i , ibal , ierr , ihi , ilo , irwork , itau ,    &
     &        iwrk , k , lwork_trevc , maxwrk , minwrk , nout
      REAL anrm , bignum , cscale , eps , scl , smlnum
      COMPLEX tmp
!     ..
!     .. Local Arrays ..
      LOGICAL select(1)
      REAL dum(1)
!     ..
!     .. External Subroutines ..
      EXTERNAL SLABAD , XERBLA , CSSCAL , CGEBAK , CGEBAL , CGEHRD ,    &
     &         CHSEQR , CLACPY , CLASCL , CSCAL , CTREVC3 , CUNGHR
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      INTEGER ISAMAX , ILAENV
      REAL SLAMCH , SCNRM2 , CLANGE
      EXTERNAL LSAME , ISAMAX , ILAENV , SLAMCH , SCNRM2 , CLANGE
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC REAL , CMPLX , CONJG , AIMAG , MAX , SQRT
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments
!
      Info = 0
      lquery = (Lwork==-1)
      wantvl = LSAME(Jobvl,'V')
      wantvr = LSAME(Jobvr,'V')
      IF ( (.NOT.wantvl) .AND. (.NOT.LSAME(Jobvl,'N')) ) THEN
         Info = -1
      ELSEIF ( (.NOT.wantvr) .AND. (.NOT.LSAME(Jobvr,'N')) ) THEN
         Info = -2
      ELSEIF ( N<0 ) THEN
         Info = -3
      ELSEIF ( Lda<MAX(1,N) ) THEN
         Info = -5
      ELSEIF ( Ldvl<1 .OR. (wantvl .AND. Ldvl<N) ) THEN
         Info = -8
      ELSEIF ( Ldvr<1 .OR. (wantvr .AND. Ldvr<N) ) THEN
         Info = -10
      ENDIF
!
!     Compute workspace
!      (Note: Comments in the code beginning "Workspace:" describe the
!       minimal amount of workspace needed at that point in the code,
!       as well as the preferred amount for good performance.
!       CWorkspace refers to complex workspace, and RWorkspace to real
!       workspace. NB refers to the optimal block size for the
!       immediately following subroutine, as returned by ILAENV.
!       HSWORK refers to the workspace preferred by CHSEQR, as
!       calculated below. HSWORK is computed assuming ILO=1 and IHI=N,
!       the worst case.)
!
      IF ( Info==0 ) THEN
         IF ( N==0 ) THEN
            minwrk = 1
            maxwrk = 1
         ELSE
            maxwrk = N + N*ILAENV(1,'CGEHRD',' ',N,1,N,0)
            minwrk = 2*N
            IF ( wantvl ) THEN
               maxwrk = MAX(maxwrk,N+(N-1)                              &
     &                  *ILAENV(1,'CUNGHR',' ',N,1,N,-1))
               CALL CTREVC3('L','B',select,N,A,Lda,Vl,Ldvl,Vr,Ldvr,N,   &
     &                      nout,Work,-1,Rwork,-1,ierr)
               lwork_trevc = INT(Work(1))
               maxwrk = MAX(maxwrk,N+lwork_trevc)
               CALL CHSEQR('S','V',N,1,N,A,Lda,W,Vl,Ldvl,Work,-1,Info)
            ELSEIF ( wantvr ) THEN
               maxwrk = MAX(maxwrk,N+(N-1)                              &
     &                  *ILAENV(1,'CUNGHR',' ',N,1,N,-1))
               CALL CTREVC3('R','B',select,N,A,Lda,Vl,Ldvl,Vr,Ldvr,N,   &
     &                      nout,Work,-1,Rwork,-1,ierr)
               lwork_trevc = INT(Work(1))
               maxwrk = MAX(maxwrk,N+lwork_trevc)
               CALL CHSEQR('S','V',N,1,N,A,Lda,W,Vr,Ldvr,Work,-1,Info)
            ELSE
               CALL CHSEQR('E','N',N,1,N,A,Lda,W,Vr,Ldvr,Work,-1,Info)
            ENDIF
            hswork = INT(Work(1))
            maxwrk = MAX(maxwrk,hswork,minwrk)
         ENDIF
         Work(1) = maxwrk
!
         IF ( Lwork<minwrk .AND. .NOT.lquery ) Info = -12
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('CGEEV ',-Info)
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
      anrm = CLANGE('M',N,N,A,Lda,dum)
      scalea = .FALSE.
      IF ( anrm>ZERO .AND. anrm<smlnum ) THEN
         scalea = .TRUE.
         cscale = smlnum
      ELSEIF ( anrm>bignum ) THEN
         scalea = .TRUE.
         cscale = bignum
      ENDIF
      IF ( scalea ) CALL CLASCL('G',0,0,anrm,cscale,N,N,A,Lda,ierr)
!
!     Balance the matrix
!     (CWorkspace: none)
!     (RWorkspace: need N)
!
      ibal = 1
      CALL CGEBAL('B',N,A,Lda,ilo,ihi,Rwork(ibal),ierr)
!
!     Reduce to upper Hessenberg form
!     (CWorkspace: need 2*N, prefer N+N*NB)
!     (RWorkspace: none)
!
      itau = 1
      iwrk = itau + N
      CALL CGEHRD(N,ilo,ihi,A,Lda,Work(itau),Work(iwrk),Lwork-iwrk+1,   &
     &            ierr)
!
      IF ( wantvl ) THEN
!
!        Want left eigenvectors
!        Copy Householder vectors to VL
!
         side = 'L'
         CALL CLACPY('L',N,N,A,Lda,Vl,Ldvl)
!
!        Generate unitary matrix in VL
!        (CWorkspace: need 2*N-1, prefer N+(N-1)*NB)
!        (RWorkspace: none)
!
         CALL CUNGHR(N,ilo,ihi,Vl,Ldvl,Work(itau),Work(iwrk),           &
     &               Lwork-iwrk+1,ierr)
!
!        Perform QR iteration, accumulating Schur vectors in VL
!        (CWorkspace: need 1, prefer HSWORK (see comments) )
!        (RWorkspace: none)
!
         iwrk = itau
         CALL CHSEQR('S','V',N,ilo,ihi,A,Lda,W,Vl,Ldvl,Work(iwrk),      &
     &               Lwork-iwrk+1,Info)
!
         IF ( wantvr ) THEN
!
!           Want left and right eigenvectors
!           Copy Schur vectors to VR
!
            side = 'B'
            CALL CLACPY('F',N,N,Vl,Ldvl,Vr,Ldvr)
         ENDIF
!
      ELSEIF ( wantvr ) THEN
!
!        Want right eigenvectors
!        Copy Householder vectors to VR
!
         side = 'R'
         CALL CLACPY('L',N,N,A,Lda,Vr,Ldvr)
!
!        Generate unitary matrix in VR
!        (CWorkspace: need 2*N-1, prefer N+(N-1)*NB)
!        (RWorkspace: none)
!
         CALL CUNGHR(N,ilo,ihi,Vr,Ldvr,Work(itau),Work(iwrk),           &
     &               Lwork-iwrk+1,ierr)
!
!        Perform QR iteration, accumulating Schur vectors in VR
!        (CWorkspace: need 1, prefer HSWORK (see comments) )
!        (RWorkspace: none)
!
         iwrk = itau
         CALL CHSEQR('S','V',N,ilo,ihi,A,Lda,W,Vr,Ldvr,Work(iwrk),      &
     &               Lwork-iwrk+1,Info)
!
      ELSE
!
!        Compute eigenvalues only
!        (CWorkspace: need 1, prefer HSWORK (see comments) )
!        (RWorkspace: none)
!
         iwrk = itau
         CALL CHSEQR('E','N',N,ilo,ihi,A,Lda,W,Vr,Ldvr,Work(iwrk),      &
     &               Lwork-iwrk+1,Info)
      ENDIF
!
!     If INFO .NE. 0 from CHSEQR, then quit
!
      IF ( Info==0 ) THEN
!
         IF ( wantvl .OR. wantvr ) THEN
!
!        Compute left and/or right eigenvectors
!        (CWorkspace: need 2*N, prefer N + 2*N*NB)
!        (RWorkspace: need 2*N)
!
            irwork = ibal + N
            CALL CTREVC3(side,'B',select,N,A,Lda,Vl,Ldvl,Vr,Ldvr,N,nout,&
     &                   Work(iwrk),Lwork-iwrk+1,Rwork(irwork),N,ierr)
         ENDIF
!
         IF ( wantvl ) THEN
!
!        Undo balancing of left eigenvectors
!        (CWorkspace: none)
!        (RWorkspace: need N)
!
            CALL CGEBAK('B','L',N,ilo,ihi,Rwork(ibal),N,Vl,Ldvl,ierr)
!
!        Normalize left eigenvectors and make largest component real
!
            DO i = 1 , N
               scl = ONE/SCNRM2(N,Vl(1,i),1)
               CALL CSSCAL(N,scl,Vl(1,i),1)
               DO k = 1 , N
                  Rwork(irwork+k-1) = REAL(Vl(k,i))**2 + AIMAG(Vl(k,i)) &
     &                                **2
               ENDDO
               k = ISAMAX(N,Rwork(irwork),1)
               tmp = CONJG(Vl(k,i))/SQRT(Rwork(irwork+k-1))
               CALL CSCAL(N,tmp,Vl(1,i),1)
               Vl(k,i) = CMPLX(REAL(Vl(k,i)),ZERO)
            ENDDO
         ENDIF
!
         IF ( wantvr ) THEN
!
!        Undo balancing of right eigenvectors
!        (CWorkspace: none)
!        (RWorkspace: need N)
!
            CALL CGEBAK('B','R',N,ilo,ihi,Rwork(ibal),N,Vr,Ldvr,ierr)
!
!        Normalize right eigenvectors and make largest component real
!
            DO i = 1 , N
               scl = ONE/SCNRM2(N,Vr(1,i),1)
               CALL CSSCAL(N,scl,Vr(1,i),1)
               DO k = 1 , N
                  Rwork(irwork+k-1) = REAL(Vr(k,i))**2 + AIMAG(Vr(k,i)) &
     &                                **2
               ENDDO
               k = ISAMAX(N,Rwork(irwork),1)
               tmp = CONJG(Vr(k,i))/SQRT(Rwork(irwork+k-1))
               CALL CSCAL(N,tmp,Vr(1,i),1)
               Vr(k,i) = CMPLX(REAL(Vr(k,i)),ZERO)
            ENDDO
         ENDIF
      ENDIF
!
!     Undo scaling if necessary
!
      IF ( scalea ) THEN
         CALL CLASCL('G',0,0,cscale,anrm,N-Info,1,W(Info+1),            &
     &               MAX(N-Info,1),ierr)
         IF ( Info>0 ) CALL CLASCL('G',0,0,cscale,anrm,ilo-1,1,W,N,ierr)
      ENDIF
!
      Work(1) = maxwrk
!
!     End of CGEEV
!
      END SUBROUTINE CGEEV
