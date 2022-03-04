!*==sgeev.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief <b> SGEEV computes the eigenvalues and, optionally, the left and/or right eigenvectors for GE matrices</b>
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SGEEV + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgeev.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgeev.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgeev.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SGEEV( JOBVL, JOBVR, N, A, LDA, WR, WI, VL, LDVL, VR,
!                         LDVR, WORK, LWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          JOBVL, JOBVR
!       INTEGER            INFO, LDA, LDVL, LDVR, LWORK, N
!       ..
!       .. Array Arguments ..
!       REAL   A( LDA, * ), VL( LDVL, * ), VR( LDVR, * ),
!      $                   WI( * ), WORK( * ), WR( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SGEEV computes for an N-by-N real nonsymmetric matrix A, the
!> eigenvalues and, optionally, the left and/or right eigenvectors.
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
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] JOBVL
!> \verbatim
!>          JOBVL is CHARACTER*1
!>          = 'N': left eigenvectors of A are not computed;
!>          = 'V': left eigenvectors of A are computed.
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
!>          A is REAL array, dimension (LDA,N)
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
!> \param[out] WR
!> \verbatim
!>          WR is REAL array, dimension (N)
!> \endverbatim
!>
!> \param[out] WI
!> \verbatim
!>          WI is REAL array, dimension (N)
!>          WR and WI contain the real and imaginary parts,
!>          respectively, of the computed eigenvalues.  Complex
!>          conjugate pairs of eigenvalues appear consecutively
!>          with the eigenvalue having the positive imaginary part
!>          first.
!> \endverbatim
!>
!> \param[out] VL
!> \verbatim
!>          VL is REAL array, dimension (LDVL,N)
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
!>          VR is REAL array, dimension (LDVR,N)
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
!>          The leading dimension of the array VR.  LDVR >= 1; if
!>          JOBVR = 'V', LDVR >= N.
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
!>          The dimension of the array WORK.  LWORK >= max(1,3*N), and
!>          if JOBVL = 'V' or JOBVR = 'V', LWORK >= 4*N.  For good
!>          performance, LWORK must generally be larger.
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
!>          > 0:  if INFO = i, the QR algorithm failed to compute all the
!>                eigenvalues, and no eigenvectors have been computed;
!>                elements i+1:N of WR and WI contain eigenvalues which
!>                have converged.
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
!  @generated from dgeev.f, fortran d -> s, Tue Apr 19 01:47:44 2016
!
!> \ingroup realGEeigen
!
!  =====================================================================
      SUBROUTINE SGEEV(Jobvl,Jobvr,N,A,Lda,Wr,Wi,Vl,Ldvl,Vr,Ldvr,Work,  &
     &                 Lwork,Info)
      USE S_ILAENV
      USE S_ISAMAX
      USE S_LSAME
      USE S_SGEBAK
      USE S_SGEBAL
      USE S_SGEHRD
      USE S_SHSEQR
      USE S_SLABAD
      USE S_SLACPY
      USE S_SLAMCH
      USE S_SLANGE
      USE S_SLAPY2
      USE S_SLARTG
      USE S_SLASCL
      USE S_SNRM2
      USE S_SORGHR
      USE S_SROT
      USE S_SSCAL
      USE S_STREVC3
      USE S_XERBLA
      IMPLICIT NONE
!*--SGEEV215
!
! PARAMETER definitions rewritten by SPAG
!
      REAL , PARAMETER  ::  ZERO = 0.0E0 , ONE = 1.0E0
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Jobvl
      CHARACTER :: Jobvr
      INTEGER :: N
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(*) :: Wr
      REAL , DIMENSION(*) :: Wi
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
      REAL :: anrm , bignum , cs , cscale , eps , r , scl , smlnum , sn
      REAL , DIMENSION(1) :: dum
      INTEGER :: hswork , i , ibal , ierr , ihi , ilo , itau , iwrk ,   &
     &           k , lwork_trevc , maxwrk , minwrk , nout
      LOGICAL :: lquery , scalea , wantvl , wantvr
      LOGICAL , DIMENSION(1) :: select
      CHARACTER :: side
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
      IF ( (.NOT.wantvl) .AND. (.NOT.LSAME(Jobvl,'N')) ) THEN
         Info = -1
      ELSEIF ( (.NOT.wantvr) .AND. (.NOT.LSAME(Jobvr,'N')) ) THEN
         Info = -2
      ELSEIF ( N<0 ) THEN
         Info = -3
      ELSEIF ( Lda<MAX(1,N) ) THEN
         Info = -5
      ELSEIF ( Ldvl<1 .OR. (wantvl .AND. Ldvl<N) ) THEN
         Info = -9
      ELSEIF ( Ldvr<1 .OR. (wantvr .AND. Ldvr<N) ) THEN
         Info = -11
      ENDIF
!
!     Compute workspace
!      (Note: Comments in the code beginning "Workspace:" describe the
!       minimal amount of workspace needed at that point in the code,
!       as well as the preferred amount for good performance.
!       NB refers to the optimal block size for the immediately
!       following subroutine, as returned by ILAENV.
!       HSWORK refers to the workspace preferred by SHSEQR, as
!       calculated below. HSWORK is computed assuming ILO=1 and IHI=N,
!       the worst case.)
!
      IF ( Info==0 ) THEN
         IF ( N==0 ) THEN
            minwrk = 1
            maxwrk = 1
         ELSE
            maxwrk = 2*N + N*ILAENV(1,'SGEHRD',' ',N,1,N,0)
            IF ( wantvl ) THEN
               minwrk = 4*N
               maxwrk = MAX(maxwrk,2*N+(N-1)                            &
     &                  *ILAENV(1,'SORGHR',' ',N,1,N,-1))
               CALL SHSEQR('S','V',N,1,N,A,Lda,Wr,Wi,Vl,Ldvl,Work,-1,   &
     &                     Info)
               hswork = INT(Work(1))
               maxwrk = MAX(maxwrk,N+1,N+hswork)
               CALL STREVC3('L','B',select,N,A,Lda,Vl,Ldvl,Vr,Ldvr,N,   &
     &                      nout,Work,-1,ierr)
               lwork_trevc = INT(Work(1))
               maxwrk = MAX(maxwrk,N+lwork_trevc)
               maxwrk = MAX(maxwrk,4*N)
            ELSEIF ( wantvr ) THEN
               minwrk = 4*N
               maxwrk = MAX(maxwrk,2*N+(N-1)                            &
     &                  *ILAENV(1,'SORGHR',' ',N,1,N,-1))
               CALL SHSEQR('S','V',N,1,N,A,Lda,Wr,Wi,Vr,Ldvr,Work,-1,   &
     &                     Info)
               hswork = INT(Work(1))
               maxwrk = MAX(maxwrk,N+1,N+hswork)
               CALL STREVC3('R','B',select,N,A,Lda,Vl,Ldvl,Vr,Ldvr,N,   &
     &                      nout,Work,-1,ierr)
               lwork_trevc = INT(Work(1))
               maxwrk = MAX(maxwrk,N+lwork_trevc)
               maxwrk = MAX(maxwrk,4*N)
            ELSE
               minwrk = 3*N
               CALL SHSEQR('E','N',N,1,N,A,Lda,Wr,Wi,Vr,Ldvr,Work,-1,   &
     &                     Info)
               hswork = INT(Work(1))
               maxwrk = MAX(maxwrk,N+1,N+hswork)
            ENDIF
            maxwrk = MAX(maxwrk,minwrk)
         ENDIF
         Work(1) = maxwrk
!
         IF ( Lwork<minwrk .AND. .NOT.lquery ) Info = -13
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('SGEEV ',-Info)
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
!     Balance the matrix
!     (Workspace: need N)
!
      ibal = 1
      CALL SGEBAL('B',N,A,Lda,ilo,ihi,Work(ibal),ierr)
!
!     Reduce to upper Hessenberg form
!     (Workspace: need 3*N, prefer 2*N+N*NB)
!
      itau = ibal + N
      iwrk = itau + N
      CALL SGEHRD(N,ilo,ihi,A,Lda,Work(itau),Work(iwrk),Lwork-iwrk+1,   &
     &            ierr)
!
      IF ( wantvl ) THEN
!
!        Want left eigenvectors
!        Copy Householder vectors to VL
!
         side = 'L'
         CALL SLACPY('L',N,N,A,Lda,Vl,Ldvl)
!
!        Generate orthogonal matrix in VL
!        (Workspace: need 3*N-1, prefer 2*N+(N-1)*NB)
!
         CALL SORGHR(N,ilo,ihi,Vl,Ldvl,Work(itau),Work(iwrk),           &
     &               Lwork-iwrk+1,ierr)
!
!        Perform QR iteration, accumulating Schur vectors in VL
!        (Workspace: need N+1, prefer N+HSWORK (see comments) )
!
         iwrk = itau
         CALL SHSEQR('S','V',N,ilo,ihi,A,Lda,Wr,Wi,Vl,Ldvl,Work(iwrk),  &
     &               Lwork-iwrk+1,Info)
!
         IF ( wantvr ) THEN
!
!           Want left and right eigenvectors
!           Copy Schur vectors to VR
!
            side = 'B'
            CALL SLACPY('F',N,N,Vl,Ldvl,Vr,Ldvr)
         ENDIF
!
      ELSEIF ( wantvr ) THEN
!
!        Want right eigenvectors
!        Copy Householder vectors to VR
!
         side = 'R'
         CALL SLACPY('L',N,N,A,Lda,Vr,Ldvr)
!
!        Generate orthogonal matrix in VR
!        (Workspace: need 3*N-1, prefer 2*N+(N-1)*NB)
!
         CALL SORGHR(N,ilo,ihi,Vr,Ldvr,Work(itau),Work(iwrk),           &
     &               Lwork-iwrk+1,ierr)
!
!        Perform QR iteration, accumulating Schur vectors in VR
!        (Workspace: need N+1, prefer N+HSWORK (see comments) )
!
         iwrk = itau
         CALL SHSEQR('S','V',N,ilo,ihi,A,Lda,Wr,Wi,Vr,Ldvr,Work(iwrk),  &
     &               Lwork-iwrk+1,Info)
!
      ELSE
!
!        Compute eigenvalues only
!        (Workspace: need N+1, prefer N+HSWORK (see comments) )
!
         iwrk = itau
         CALL SHSEQR('E','N',N,ilo,ihi,A,Lda,Wr,Wi,Vr,Ldvr,Work(iwrk),  &
     &               Lwork-iwrk+1,Info)
      ENDIF
!
!     If INFO .NE. 0 from SHSEQR, then quit
!
      IF ( Info==0 ) THEN
!
!
!        Compute left and/or right eigenvectors
!        (Workspace: need 4*N, prefer N + N + 2*N*NB)
!
         IF ( wantvl .OR. wantvr ) CALL STREVC3(side,'B',select,N,A,Lda,&
     &        Vl,Ldvl,Vr,Ldvr,N,nout,Work(iwrk),Lwork-iwrk+1,ierr)
!
         IF ( wantvl ) THEN
!
!        Undo balancing of left eigenvectors
!        (Workspace: need N)
!
            CALL SGEBAK('B','L',N,ilo,ihi,Work(ibal),N,Vl,Ldvl,ierr)
!
!        Normalize left eigenvectors and make largest component real
!
            DO i = 1 , N
               IF ( Wi(i)==ZERO ) THEN
                  scl = ONE/SNRM2(N,Vl(1,i),1)
                  CALL SSCAL(N,scl,Vl(1,i),1)
               ELSEIF ( Wi(i)>ZERO ) THEN
                  scl = ONE/SLAPY2(SNRM2(N,Vl(1,i),1),                  &
     &                  SNRM2(N,Vl(1,i+1),1))
                  CALL SSCAL(N,scl,Vl(1,i),1)
                  CALL SSCAL(N,scl,Vl(1,i+1),1)
                  DO k = 1 , N
                     Work(iwrk+k-1) = Vl(k,i)**2 + Vl(k,i+1)**2
                  ENDDO
                  k = ISAMAX(N,Work(iwrk),1)
                  CALL SLARTG(Vl(k,i),Vl(k,i+1),cs,sn,r)
                  CALL SROT(N,Vl(1,i),1,Vl(1,i+1),1,cs,sn)
                  Vl(k,i+1) = ZERO
               ENDIF
            ENDDO
         ENDIF
!
         IF ( wantvr ) THEN
!
!        Undo balancing of right eigenvectors
!        (Workspace: need N)
!
            CALL SGEBAK('B','R',N,ilo,ihi,Work(ibal),N,Vr,Ldvr,ierr)
!
!        Normalize right eigenvectors and make largest component real
!
            DO i = 1 , N
               IF ( Wi(i)==ZERO ) THEN
                  scl = ONE/SNRM2(N,Vr(1,i),1)
                  CALL SSCAL(N,scl,Vr(1,i),1)
               ELSEIF ( Wi(i)>ZERO ) THEN
                  scl = ONE/SLAPY2(SNRM2(N,Vr(1,i),1),                  &
     &                  SNRM2(N,Vr(1,i+1),1))
                  CALL SSCAL(N,scl,Vr(1,i),1)
                  CALL SSCAL(N,scl,Vr(1,i+1),1)
                  DO k = 1 , N
                     Work(iwrk+k-1) = Vr(k,i)**2 + Vr(k,i+1)**2
                  ENDDO
                  k = ISAMAX(N,Work(iwrk),1)
                  CALL SLARTG(Vr(k,i),Vr(k,i+1),cs,sn,r)
                  CALL SROT(N,Vr(1,i),1,Vr(1,i+1),1,cs,sn)
                  Vr(k,i+1) = ZERO
               ENDIF
            ENDDO
         ENDIF
      ENDIF
!
!     Undo scaling if necessary
!
      IF ( scalea ) THEN
         CALL SLASCL('G',0,0,cscale,anrm,N-Info,1,Wr(Info+1),           &
     &               MAX(N-Info,1),ierr)
         CALL SLASCL('G',0,0,cscale,anrm,N-Info,1,Wi(Info+1),           &
     &               MAX(N-Info,1),ierr)
         IF ( Info>0 ) THEN
            CALL SLASCL('G',0,0,cscale,anrm,ilo-1,1,Wr,N,ierr)
            CALL SLASCL('G',0,0,cscale,anrm,ilo-1,1,Wi,N,ierr)
         ENDIF
      ENDIF
!
      Work(1) = maxwrk
!
!     End of SGEEV
!
      END SUBROUTINE SGEEV
