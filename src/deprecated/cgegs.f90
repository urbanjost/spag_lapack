!*==cgegs.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief <b> CGEEVX computes the eigenvalues and, optionally, the left and/or right eigenvectors for GE matrices</b>
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CGEGS + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgegs.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgegs.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgegs.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CGEGS( JOBVSL, JOBVSR, N, A, LDA, B, LDB, ALPHA, BETA,
!                         VSL, LDVSL, VSR, LDVSR, WORK, LWORK, RWORK,
!                         INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          JOBVSL, JOBVSR
!       INTEGER            INFO, LDA, LDB, LDVSL, LDVSR, LWORK, N
!       ..
!       .. Array Arguments ..
!       REAL               RWORK( * )
!       COMPLEX            A( LDA, * ), ALPHA( * ), B( LDB, * ),
!      $                   BETA( * ), VSL( LDVSL, * ), VSR( LDVSR, * ),
!      $                   WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> This routine is deprecated and has been replaced by routine CGGES.
!>
!> CGEGS computes the eigenvalues, Schur form, and, optionally, the
!> left and or/right Schur vectors of a complex matrix pair (A,B).
!> Given two square matrices A and B, the generalized Schur
!> factorization has the form
!>
!>    A = Q*S*Z**H,  B = Q*T*Z**H
!>
!> where Q and Z are unitary matrices and S and T are upper triangular.
!> The columns of Q are the left Schur vectors
!> and the columns of Z are the right Schur vectors.
!>
!> If only the eigenvalues of (A,B) are needed, the driver routine
!> CGEGV should be used instead.  See CGEGV for a description of the
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
!>          A is COMPLEX array, dimension (LDA, N)
!>          On entry, the matrix A.
!>          On exit, the upper triangular matrix S from the generalized
!>          Schur factorization.
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
!>          On entry, the matrix B.
!>          On exit, the upper triangular matrix T from the generalized
!>          Schur factorization.
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
!>          The complex scalars alpha that define the eigenvalues of
!>          GNEP.  ALPHA(j) = S(j,j), the diagonal element of the Schur
!>          form of A.
!> \endverbatim
!>
!> \param[out] BETA
!> \verbatim
!>          BETA is COMPLEX array, dimension (N)
!>          The non-negative real scalars beta that define the
!>          eigenvalues of GNEP.  BETA(j) = T(j,j), the diagonal element
!>          of the triangular factor T.
!>
!>          Together, the quantities alpha = ALPHA(j) and beta = BETA(j)
!>          represent the j-th eigenvalue of the matrix pair (A,B), in
!>          one of the forms lambda = alpha/beta or mu = beta/alpha.
!>          Since either lambda or mu may overflow, they should not,
!>          in general, be computed.
!> \endverbatim
!>
!> \param[out] VSL
!> \verbatim
!>          VSL is COMPLEX array, dimension (LDVSL,N)
!>          If JOBVSL = 'V', the matrix of left Schur vectors Q.
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
!>          WORK is COMPLEX array, dimension (MAX(1,LWORK))
!>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The dimension of the array WORK.  LWORK >= max(1,2*N).
!>          For good performance, LWORK must generally be larger.
!>          To compute the optimal value of LWORK, call ILAENV to get
!>          blocksizes (for CGEQRF, CUNMQR, and CUNGQR.)  Then compute:
!>          NB  -- MAX of the blocksizes for CGEQRF, CUNMQR, and CUNGQR;
!>          the optimal LWORK is N*(NB+1).
!>
!>          If LWORK = -1, then a workspace query is assumed; the routine
!>          only calculates the optimal size of the WORK array, returns
!>          this value as the first entry of the WORK array, and no error
!>          message related to LWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is REAL array, dimension (3*N)
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
!>          > N:  errors that usually indicate LAPACK problems:
!>                =N+1: error return from CGGBAL
!>                =N+2: error return from CGEQRF
!>                =N+3: error return from CUNMQR
!>                =N+4: error return from CUNGQR
!>                =N+5: error return from CGGHRD
!>                =N+6: error return from CHGEQZ (other than failed
!>                                               iteration)
!>                =N+7: error return from CGGBAK (computing VSL)
!>                =N+8: error return from CGGBAK (computing VSR)
!>                =N+9: error return from CLASCL (various places)
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
!> \ingroup complexGEeigen
!
!  =====================================================================
      SUBROUTINE CGEGS(Jobvsl,Jobvsr,N,A,Lda,B,Ldb,Alpha,Beta,Vsl,Ldvsl,&
     &                 Vsr,Ldvsr,Work,Lwork,Rwork,Info)
      IMPLICIT NONE
!*--CGEGS228
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
      REAL Rwork(*)
      COMPLEX A(Lda,*) , Alpha(*) , B(Ldb,*) , Beta(*) , Vsl(Ldvsl,*) , &
     &        Vsr(Ldvsr,*) , Work(*)
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
      LOGICAL ilascl , ilbscl , ilvsl , ilvsr , lquery
      INTEGER icols , ihi , iinfo , ijobvl , ijobvr , ileft , ilo ,     &
     &        iright , irows , irwork , itau , iwork , lopt , lwkmin ,  &
     &        lwkopt , nb , nb1 , nb2 , nb3
      REAL anrm , anrmto , bignum , bnrm , bnrmto , eps , safmin ,      &
     &     smlnum
!     ..
!     .. External Subroutines ..
      EXTERNAL CGEQRF , CGGBAK , CGGBAL , CGGHRD , CHGEQZ , CLACPY ,    &
     &         CLASCL , CLASET , CUNGQR , CUNMQR , XERBLA
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      INTEGER ILAENV
      REAL CLANGE , SLAMCH
      EXTERNAL ILAENV , LSAME , CLANGE , SLAMCH
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
      lwkmin = MAX(2*N,1)
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
         Info = -11
      ELSEIF ( Ldvsr<1 .OR. (ilvsr .AND. Ldvsr<N) ) THEN
         Info = -13
      ELSEIF ( Lwork<lwkmin .AND. .NOT.lquery ) THEN
         Info = -15
      ENDIF
!
      IF ( Info==0 ) THEN
         nb1 = ILAENV(1,'CGEQRF',' ',N,N,-1,-1)
         nb2 = ILAENV(1,'CUNMQR',' ',N,N,N,-1)
         nb3 = ILAENV(1,'CUNGQR',' ',N,N,N,-1)
         nb = MAX(nb1,nb2,nb3)
         lopt = N*(nb+1)
         Work(1) = lopt
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('CGEGS ',-Info)
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
      safmin = SLAMCH('S')
      smlnum = N*safmin/eps
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
      IF ( ilascl ) THEN
         CALL CLASCL('G',-1,-1,anrm,anrmto,N,N,A,Lda,iinfo)
         IF ( iinfo/=0 ) THEN
            Info = N + 9
            RETURN
         ENDIF
      ENDIF
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
      IF ( ilbscl ) THEN
         CALL CLASCL('G',-1,-1,bnrm,bnrmto,N,N,B,Ldb,iinfo)
         IF ( iinfo/=0 ) THEN
            Info = N + 9
            RETURN
         ENDIF
      ENDIF
!
!     Permute the matrix to make it more nearly triangular
!
      ileft = 1
      iright = N + 1
      irwork = iright + N
      iwork = 1
      CALL CGGBAL('P',N,A,Lda,B,Ldb,ilo,ihi,Rwork(ileft),Rwork(iright), &
     &            Rwork(irwork),iinfo)
      IF ( iinfo/=0 ) THEN
         Info = N + 1
         GOTO 100
      ENDIF
!
!     Reduce B to triangular form, and initialize VSL and/or VSR
!
      irows = ihi + 1 - ilo
      icols = N + 1 - ilo
      itau = iwork
      iwork = itau + irows
      CALL CGEQRF(irows,icols,B(ilo,ilo),Ldb,Work(itau),Work(iwork),    &
     &            Lwork+1-iwork,iinfo)
      IF ( iinfo>=0 ) lwkopt = MAX(lwkopt,INT(Work(iwork))+iwork-1)
      IF ( iinfo/=0 ) THEN
         Info = N + 2
         GOTO 100
      ENDIF
!
      CALL CUNMQR('L','C',irows,icols,irows,B(ilo,ilo),Ldb,Work(itau),  &
     &            A(ilo,ilo),Lda,Work(iwork),Lwork+1-iwork,iinfo)
      IF ( iinfo>=0 ) lwkopt = MAX(lwkopt,INT(Work(iwork))+iwork-1)
      IF ( iinfo/=0 ) THEN
         Info = N + 3
         GOTO 100
      ENDIF
!
      IF ( ilvsl ) THEN
         CALL CLASET('Full',N,N,CZERO,CONE,Vsl,Ldvsl)
         CALL CLACPY('L',irows-1,irows-1,B(ilo+1,ilo),Ldb,Vsl(ilo+1,ilo)&
     &               ,Ldvsl)
         CALL CUNGQR(irows,irows,irows,Vsl(ilo,ilo),Ldvsl,Work(itau),   &
     &               Work(iwork),Lwork+1-iwork,iinfo)
         IF ( iinfo>=0 ) lwkopt = MAX(lwkopt,INT(Work(iwork))+iwork-1)
         IF ( iinfo/=0 ) THEN
            Info = N + 4
            GOTO 100
         ENDIF
      ENDIF
!
      IF ( ilvsr ) CALL CLASET('Full',N,N,CZERO,CONE,Vsr,Ldvsr)
!
!     Reduce to generalized Hessenberg form
!
      CALL CGGHRD(Jobvsl,Jobvsr,N,ilo,ihi,A,Lda,B,Ldb,Vsl,Ldvsl,Vsr,    &
     &            Ldvsr,iinfo)
      IF ( iinfo/=0 ) THEN
         Info = N + 5
         GOTO 100
      ENDIF
!
!     Perform QZ algorithm, computing Schur vectors if desired
!
      iwork = itau
      CALL CHGEQZ('S',Jobvsl,Jobvsr,N,ilo,ihi,A,Lda,B,Ldb,Alpha,Beta,   &
     &            Vsl,Ldvsl,Vsr,Ldvsr,Work(iwork),Lwork+1-iwork,        &
     &            Rwork(irwork),iinfo)
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
         CALL CGGBAK('P','L',N,ilo,ihi,Rwork(ileft),Rwork(iright),N,Vsl,&
     &               Ldvsl,iinfo)
         IF ( iinfo/=0 ) THEN
            Info = N + 7
            GOTO 100
         ENDIF
      ENDIF
      IF ( ilvsr ) THEN
         CALL CGGBAK('P','R',N,ilo,ihi,Rwork(ileft),Rwork(iright),N,Vsr,&
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
         CALL CLASCL('U',-1,-1,anrmto,anrm,N,N,A,Lda,iinfo)
         IF ( iinfo/=0 ) THEN
            Info = N + 9
            RETURN
         ENDIF
         CALL CLASCL('G',-1,-1,anrmto,anrm,N,1,Alpha,N,iinfo)
         IF ( iinfo/=0 ) THEN
            Info = N + 9
            RETURN
         ENDIF
      ENDIF
!
      IF ( ilbscl ) THEN
         CALL CLASCL('U',-1,-1,bnrmto,bnrm,N,N,B,Ldb,iinfo)
         IF ( iinfo/=0 ) THEN
            Info = N + 9
            RETURN
         ENDIF
         CALL CLASCL('G',-1,-1,bnrmto,bnrm,N,1,Beta,N,iinfo)
         IF ( iinfo/=0 ) THEN
            Info = N + 9
            RETURN
         ENDIF
      ENDIF
!
 100  Work(1) = lwkopt
!
!
!     End of CGEGS
!
      END SUBROUTINE CGEGS
