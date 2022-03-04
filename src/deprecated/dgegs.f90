!*==dgegs.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief <b> DGEEVX computes the eigenvalues and, optionally, the left and/or right eigenvectors for GE matrices</b>
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DGEGS + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgegs.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgegs.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgegs.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DGEGS( JOBVSL, JOBVSR, N, A, LDA, B, LDB, ALPHAR,
!                         ALPHAI, BETA, VSL, LDVSL, VSR, LDVSR, WORK,
!                         LWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          JOBVSL, JOBVSR
!       INTEGER            INFO, LDA, LDB, LDVSL, LDVSR, LWORK, N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   A( LDA, * ), ALPHAI( * ), ALPHAR( * ),
!      $                   B( LDB, * ), BETA( * ), VSL( LDVSL, * ),
!      $                   VSR( LDVSR, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> This routine is deprecated and has been replaced by routine DGGES.
!>
!> DGEGS computes the eigenvalues, real Schur form, and, optionally,
!> left and or/right Schur vectors of a real matrix pair (A,B).
!> Given two square matrices A and B, the generalized real Schur
!> factorization has the form
!>
!>   A = Q*S*Z**T,  B = Q*T*Z**T
!>
!> where Q and Z are orthogonal matrices, T is upper triangular, and S
!> is an upper quasi-triangular matrix with 1-by-1 and 2-by-2 diagonal
!> blocks, the 2-by-2 blocks corresponding to complex conjugate pairs
!> of eigenvalues of (A,B).  The columns of Q are the left Schur vectors
!> and the columns of Z are the right Schur vectors.
!>
!> If only the eigenvalues of (A,B) are needed, the driver routine
!> DGEGV should be used instead.  See DGEGV for a description of the
!> eigenvalues of the generalized nonsymmetric eigenvalue problem
!> (GNEP).
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] JOBVSL
!> \verbatim
!>          JOBVSL is CHARACTER*1
!>          = 'N':  do not compute the left Schur vectors;
!>          = 'V':  compute the left Schur vectors (returned in VSL).
!> \endverbatim
!>
!> \param[in] JOBVSR
!> \verbatim
!>          JOBVSR is CHARACTER*1
!>          = 'N':  do not compute the right Schur vectors;
!>          = 'V':  compute the right Schur vectors (returned in VSR).
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
!>          A is DOUBLE PRECISION array, dimension (LDA, N)
!>          On entry, the matrix A.
!>          On exit, the upper quasi-triangular matrix S from the
!>          generalized real Schur factorization.
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
!>          On exit, the upper triangular matrix T from the generalized
!>          real Schur factorization.
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
!>          The real parts of each scalar alpha defining an eigenvalue
!>          of GNEP.
!> \endverbatim
!>
!> \param[out] ALPHAI
!> \verbatim
!>          ALPHAI is DOUBLE PRECISION array, dimension (N)
!>          The imaginary parts of each scalar alpha defining an
!>          eigenvalue of GNEP.  If ALPHAI(j) is zero, then the j-th
!>          eigenvalue is real; if positive, then the j-th and (j+1)-st
!>          eigenvalues are a complex conjugate pair, with
!>          ALPHAI(j+1) = -ALPHAI(j).
!> \endverbatim
!>
!> \param[out] BETA
!> \verbatim
!>          BETA is DOUBLE PRECISION array, dimension (N)
!>          The scalars beta that define the eigenvalues of GNEP.
!>          Together, the quantities alpha = (ALPHAR(j),ALPHAI(j)) and
!>          beta = BETA(j) represent the j-th eigenvalue of the matrix
!>          pair (A,B), in one of the forms lambda = alpha/beta or
!>          mu = beta/alpha.  Since either lambda or mu may overflow,
!>          they should not, in general, be computed.
!> \endverbatim
!>
!> \param[out] VSL
!> \verbatim
!>          VSL is DOUBLE PRECISION array, dimension (LDVSL,N)
!>          If JOBVSL = 'V', the matrix of left Schur vectors Q.
!>          Not referenced if JOBVSL = 'N'.
!> \endverbatim
!>
!> \param[in] LDVSL
!> \verbatim
!>          LDVSL is INTEGER
!>          The leading dimension of the matrix VSL. LDVSL >=1, and
!>          if JOBVSL = 'V', LDVSL >= N.
!> \endverbatim
!>
!> \param[out] VSR
!> \verbatim
!>          VSR is DOUBLE PRECISION array, dimension (LDVSR,N)
!>          If JOBVSR = 'V', the matrix of right Schur vectors Z.
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
!>          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))
!>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The dimension of the array WORK.  LWORK >= max(1,4*N).
!>          For good performance, LWORK must generally be larger.
!>          To compute the optimal value of LWORK, call ILAENV to get
!>          blocksizes (for DGEQRF, DORMQR, and DORGQR.)  Then compute:
!>          NB  -- MAX of the blocksizes for DGEQRF, DORMQR, and DORGQR
!>          The optimal LWORK is  2*N + N*(NB+1).
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
!>                The QZ iteration failed.  (A,B) are not in Schur
!>                form, but ALPHAR(j), ALPHAI(j), and BETA(j) should
!>                be correct for j=INFO+1,...,N.
!>          > N:  errors that usually indicate LAPACK problems:
!>                =N+1: error return from DGGBAL
!>                =N+2: error return from DGEQRF
!>                =N+3: error return from DORMQR
!>                =N+4: error return from DORGQR
!>                =N+5: error return from DGGHRD
!>                =N+6: error return from DHGEQZ (other than failed
!>                                                iteration)
!>                =N+7: error return from DGGBAK (computing VSL)
!>                =N+8: error return from DGGBAK (computing VSR)
!>                =N+9: error return from DLASCL (various places)
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
      SUBROUTINE DGEGS(Jobvsl,Jobvsr,N,A,Lda,B,Ldb,Alphar,Alphai,Beta,  &
     &                 Vsl,Ldvsl,Vsr,Ldvsr,Work,Lwork,Info)
      IMPLICIT NONE
!*--DGEGS230
!
!  -- LAPACK driver routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Jobvsl , Jobvsr
      INTEGER Info , Lda , Ldb , Ldvsl , Ldvsr , Lwork , N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION A(Lda,*) , Alphai(*) , Alphar(*) , B(Ldb,*) ,    &
     &                 Beta(*) , Vsl(Ldvsl,*) , Vsr(Ldvsr,*) , Work(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO , ONE
      PARAMETER (ZERO=0.0D0,ONE=1.0D0)
!     ..
!     .. Local Scalars ..
      LOGICAL ilascl , ilbscl , ilvsl , ilvsr , lquery
      INTEGER icols , ihi , iinfo , ijobvl , ijobvr , ileft , ilo ,     &
     &        iright , irows , itau , iwork , lopt , lwkmin , lwkopt ,  &
     &        nb , nb1 , nb2 , nb3
      DOUBLE PRECISION anrm , anrmto , bignum , bnrm , bnrmto , eps ,   &
     &                 safmin , smlnum
!     ..
!     .. External Subroutines ..
      EXTERNAL DGEQRF , DGGBAK , DGGBAL , DGGHRD , DHGEQZ , DLACPY ,    &
     &         DLASCL , DLASET , DORGQR , DORMQR , XERBLA
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      INTEGER ILAENV
      DOUBLE PRECISION DLAMCH , DLANGE
      EXTERNAL LSAME , ILAENV , DLAMCH , DLANGE
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC INT , MAX
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
!     Test the input arguments
!
      lwkmin = MAX(4*N,1)
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
      ELSEIF ( Ldvsl<1 .OR. (ilvsl .AND. Ldvsl<N) ) THEN
         Info = -12
      ELSEIF ( Ldvsr<1 .OR. (ilvsr .AND. Ldvsr<N) ) THEN
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
         lopt = 2*N + N*(nb+1)
         Work(1) = lopt
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('DGEGS ',-Info)
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
      smlnum = N*safmin/eps
      bignum = ONE/smlnum
!
!     Scale A if max element outside range [SMLNUM,BIGNUM]
!
      anrm = DLANGE('M',N,N,A,Lda,Work)
      ilascl = .FALSE.
      IF ( anrm>ZERO .AND. anrm<smlnum ) THEN
         anrmto = smlnum
         ilascl = .TRUE.
      ELSEIF ( anrm>bignum ) THEN
         anrmto = bignum
         ilascl = .TRUE.
      ENDIF
!
      IF ( ilascl ) THEN
         CALL DLASCL('G',-1,-1,anrm,anrmto,N,N,A,Lda,iinfo)
         IF ( iinfo/=0 ) THEN
            Info = N + 9
            RETURN
         ENDIF
      ENDIF
!
!     Scale B if max element outside range [SMLNUM,BIGNUM]
!
      bnrm = DLANGE('M',N,N,B,Ldb,Work)
      ilbscl = .FALSE.
      IF ( bnrm>ZERO .AND. bnrm<smlnum ) THEN
         bnrmto = smlnum
         ilbscl = .TRUE.
      ELSEIF ( bnrm>bignum ) THEN
         bnrmto = bignum
         ilbscl = .TRUE.
      ENDIF
!
      IF ( ilbscl ) THEN
         CALL DLASCL('G',-1,-1,bnrm,bnrmto,N,N,B,Ldb,iinfo)
         IF ( iinfo/=0 ) THEN
            Info = N + 9
            RETURN
         ENDIF
      ENDIF
!
!     Permute the matrix to make it more nearly triangular
!     Workspace layout:  (2*N words -- "work..." not actually used)
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
!     Reduce B to triangular form, and initialize VSL and/or VSR
!     Workspace layout:  ("work..." must have at least N words)
!        left_permutation, right_permutation, tau, work...
!
      irows = ihi + 1 - ilo
      icols = N + 1 - ilo
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
      IF ( ilvsl ) THEN
         CALL DLASET('Full',N,N,ZERO,ONE,Vsl,Ldvsl)
         CALL DLACPY('L',irows-1,irows-1,B(ilo+1,ilo),Ldb,Vsl(ilo+1,ilo)&
     &               ,Ldvsl)
         CALL DORGQR(irows,irows,irows,Vsl(ilo,ilo),Ldvsl,Work(itau),   &
     &               Work(iwork),Lwork+1-iwork,iinfo)
         IF ( iinfo>=0 ) lwkopt = MAX(lwkopt,INT(Work(iwork))+iwork-1)
         IF ( iinfo/=0 ) THEN
            Info = N + 4
            GOTO 100
         ENDIF
      ENDIF
!
      IF ( ilvsr ) CALL DLASET('Full',N,N,ZERO,ONE,Vsr,Ldvsr)
!
!     Reduce to generalized Hessenberg form
!
      CALL DGGHRD(Jobvsl,Jobvsr,N,ilo,ihi,A,Lda,B,Ldb,Vsl,Ldvsl,Vsr,    &
     &            Ldvsr,iinfo)
      IF ( iinfo/=0 ) THEN
         Info = N + 5
         GOTO 100
      ENDIF
!
!     Perform QZ algorithm, computing Schur vectors if desired
!     Workspace layout:  ("work..." must have at least 1 word)
!        left_permutation, right_permutation, work...
!
      iwork = itau
      CALL DHGEQZ('S',Jobvsl,Jobvsr,N,ilo,ihi,A,Lda,B,Ldb,Alphar,Alphai,&
     &            Beta,Vsl,Ldvsl,Vsr,Ldvsr,Work(iwork),Lwork+1-iwork,   &
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
!     Apply permutation to VSL and VSR
!
      IF ( ilvsl ) THEN
         CALL DGGBAK('P','L',N,ilo,ihi,Work(ileft),Work(iright),N,Vsl,  &
     &               Ldvsl,iinfo)
         IF ( iinfo/=0 ) THEN
            Info = N + 7
            GOTO 100
         ENDIF
      ENDIF
      IF ( ilvsr ) THEN
         CALL DGGBAK('P','R',N,ilo,ihi,Work(ileft),Work(iright),N,Vsr,  &
     &               Ldvsr,iinfo)
         IF ( iinfo/=0 ) THEN
            Info = N + 8
            GOTO 100
         ENDIF
      ENDIF
!
!     Undo scaling
!
      IF ( ilascl ) THEN
         CALL DLASCL('H',-1,-1,anrmto,anrm,N,N,A,Lda,iinfo)
         IF ( iinfo/=0 ) THEN
            Info = N + 9
            RETURN
         ENDIF
         CALL DLASCL('G',-1,-1,anrmto,anrm,N,1,Alphar,N,iinfo)
         IF ( iinfo/=0 ) THEN
            Info = N + 9
            RETURN
         ENDIF
         CALL DLASCL('G',-1,-1,anrmto,anrm,N,1,Alphai,N,iinfo)
         IF ( iinfo/=0 ) THEN
            Info = N + 9
            RETURN
         ENDIF
      ENDIF
!
      IF ( ilbscl ) THEN
         CALL DLASCL('U',-1,-1,bnrmto,bnrm,N,N,B,Ldb,iinfo)
         IF ( iinfo/=0 ) THEN
            Info = N + 9
            RETURN
         ENDIF
         CALL DLASCL('G',-1,-1,bnrmto,bnrm,N,1,Beta,N,iinfo)
         IF ( iinfo/=0 ) THEN
            Info = N + 9
            RETURN
         ENDIF
      ENDIF
!
 100  Work(1) = lwkopt
!
!
!     End of DGEGS
!
      END SUBROUTINE DGEGS
