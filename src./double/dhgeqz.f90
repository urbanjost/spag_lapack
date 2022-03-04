!*==dhgeqz.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b DHGEQZ
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DHGEQZ + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dhgeqz.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dhgeqz.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dhgeqz.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DHGEQZ( JOB, COMPQ, COMPZ, N, ILO, IHI, H, LDH, T, LDT,
!                          ALPHAR, ALPHAI, BETA, Q, LDQ, Z, LDZ, WORK,
!                          LWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          COMPQ, COMPZ, JOB
!       INTEGER            IHI, ILO, INFO, LDH, LDQ, LDT, LDZ, LWORK, N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   ALPHAI( * ), ALPHAR( * ), BETA( * ),
!      $                   H( LDH, * ), Q( LDQ, * ), T( LDT, * ),
!      $                   WORK( * ), Z( LDZ, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DHGEQZ computes the eigenvalues of a real matrix pair (H,T),
!> where H is an upper Hessenberg matrix and T is upper triangular,
!> using the double-shift QZ method.
!> Matrix pairs of this type are produced by the reduction to
!> generalized upper Hessenberg form of a real matrix pair (A,B):
!>
!>    A = Q1*H*Z1**T,  B = Q1*T*Z1**T,
!>
!> as computed by DGGHRD.
!>
!> If JOB='S', then the Hessenberg-triangular pair (H,T) is
!> also reduced to generalized Schur form,
!>
!>    H = Q*S*Z**T,  T = Q*P*Z**T,
!>
!> where Q and Z are orthogonal matrices, P is an upper triangular
!> matrix, and S is a quasi-triangular matrix with 1-by-1 and 2-by-2
!> diagonal blocks.
!>
!> The 1-by-1 blocks correspond to real eigenvalues of the matrix pair
!> (H,T) and the 2-by-2 blocks correspond to complex conjugate pairs of
!> eigenvalues.
!>
!> Additionally, the 2-by-2 upper triangular diagonal blocks of P
!> corresponding to 2-by-2 blocks of S are reduced to positive diagonal
!> form, i.e., if S(j+1,j) is non-zero, then P(j+1,j) = P(j,j+1) = 0,
!> P(j,j) > 0, and P(j+1,j+1) > 0.
!>
!> Optionally, the orthogonal matrix Q from the generalized Schur
!> factorization may be postmultiplied into an input matrix Q1, and the
!> orthogonal matrix Z may be postmultiplied into an input matrix Z1.
!> If Q1 and Z1 are the orthogonal matrices from DGGHRD that reduced
!> the matrix pair (A,B) to generalized upper Hessenberg form, then the
!> output matrices Q1*Q and Z1*Z are the orthogonal factors from the
!> generalized Schur factorization of (A,B):
!>
!>    A = (Q1*Q)*S*(Z1*Z)**T,  B = (Q1*Q)*P*(Z1*Z)**T.
!>
!> To avoid overflow, eigenvalues of the matrix pair (H,T) (equivalently,
!> of (A,B)) are computed as a pair of values (alpha,beta), where alpha is
!> complex and beta real.
!> If beta is nonzero, lambda = alpha / beta is an eigenvalue of the
!> generalized nonsymmetric eigenvalue problem (GNEP)
!>    A*x = lambda*B*x
!> and if alpha is nonzero, mu = beta / alpha is an eigenvalue of the
!> alternate form of the GNEP
!>    mu*A*y = B*y.
!> Real eigenvalues can be read directly from the generalized Schur
!> form:
!>   alpha = S(i,i), beta = P(i,i).
!>
!> Ref: C.B. Moler & G.W. Stewart, "An Algorithm for Generalized Matrix
!>      Eigenvalue Problems", SIAM J. Numer. Anal., 10(1973),
!>      pp. 241--256.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] JOB
!> \verbatim
!>          JOB is CHARACTER*1
!>          = 'E': Compute eigenvalues only;
!>          = 'S': Compute eigenvalues and the Schur form.
!> \endverbatim
!>
!> \param[in] COMPQ
!> \verbatim
!>          COMPQ is CHARACTER*1
!>          = 'N': Left Schur vectors (Q) are not computed;
!>          = 'I': Q is initialized to the unit matrix and the matrix Q
!>                 of left Schur vectors of (H,T) is returned;
!>          = 'V': Q must contain an orthogonal matrix Q1 on entry and
!>                 the product Q1*Q is returned.
!> \endverbatim
!>
!> \param[in] COMPZ
!> \verbatim
!>          COMPZ is CHARACTER*1
!>          = 'N': Right Schur vectors (Z) are not computed;
!>          = 'I': Z is initialized to the unit matrix and the matrix Z
!>                 of right Schur vectors of (H,T) is returned;
!>          = 'V': Z must contain an orthogonal matrix Z1 on entry and
!>                 the product Z1*Z is returned.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrices H, T, Q, and Z.  N >= 0.
!> \endverbatim
!>
!> \param[in] ILO
!> \verbatim
!>          ILO is INTEGER
!> \endverbatim
!>
!> \param[in] IHI
!> \verbatim
!>          IHI is INTEGER
!>          ILO and IHI mark the rows and columns of H which are in
!>          Hessenberg form.  It is assumed that A is already upper
!>          triangular in rows and columns 1:ILO-1 and IHI+1:N.
!>          If N > 0, 1 <= ILO <= IHI <= N; if N = 0, ILO=1 and IHI=0.
!> \endverbatim
!>
!> \param[in,out] H
!> \verbatim
!>          H is DOUBLE PRECISION array, dimension (LDH, N)
!>          On entry, the N-by-N upper Hessenberg matrix H.
!>          On exit, if JOB = 'S', H contains the upper quasi-triangular
!>          matrix S from the generalized Schur factorization.
!>          If JOB = 'E', the diagonal blocks of H match those of S, but
!>          the rest of H is unspecified.
!> \endverbatim
!>
!> \param[in] LDH
!> \verbatim
!>          LDH is INTEGER
!>          The leading dimension of the array H.  LDH >= max( 1, N ).
!> \endverbatim
!>
!> \param[in,out] T
!> \verbatim
!>          T is DOUBLE PRECISION array, dimension (LDT, N)
!>          On entry, the N-by-N upper triangular matrix T.
!>          On exit, if JOB = 'S', T contains the upper triangular
!>          matrix P from the generalized Schur factorization;
!>          2-by-2 diagonal blocks of P corresponding to 2-by-2 blocks of S
!>          are reduced to positive diagonal form, i.e., if H(j+1,j) is
!>          non-zero, then T(j+1,j) = T(j,j+1) = 0, T(j,j) > 0, and
!>          T(j+1,j+1) > 0.
!>          If JOB = 'E', the diagonal blocks of T match those of P, but
!>          the rest of T is unspecified.
!> \endverbatim
!>
!> \param[in] LDT
!> \verbatim
!>          LDT is INTEGER
!>          The leading dimension of the array T.  LDT >= max( 1, N ).
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
!>          eigenvalue of GNEP.
!>          If ALPHAI(j) is zero, then the j-th eigenvalue is real; if
!>          positive, then the j-th and (j+1)-st eigenvalues are a
!>          complex conjugate pair, with ALPHAI(j+1) = -ALPHAI(j).
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
!> \param[in,out] Q
!> \verbatim
!>          Q is DOUBLE PRECISION array, dimension (LDQ, N)
!>          On entry, if COMPQ = 'V', the orthogonal matrix Q1 used in
!>          the reduction of (A,B) to generalized Hessenberg form.
!>          On exit, if COMPQ = 'I', the orthogonal matrix of left Schur
!>          vectors of (H,T), and if COMPQ = 'V', the orthogonal matrix
!>          of left Schur vectors of (A,B).
!>          Not referenced if COMPQ = 'N'.
!> \endverbatim
!>
!> \param[in] LDQ
!> \verbatim
!>          LDQ is INTEGER
!>          The leading dimension of the array Q.  LDQ >= 1.
!>          If COMPQ='V' or 'I', then LDQ >= N.
!> \endverbatim
!>
!> \param[in,out] Z
!> \verbatim
!>          Z is DOUBLE PRECISION array, dimension (LDZ, N)
!>          On entry, if COMPZ = 'V', the orthogonal matrix Z1 used in
!>          the reduction of (A,B) to generalized Hessenberg form.
!>          On exit, if COMPZ = 'I', the orthogonal matrix of
!>          right Schur vectors of (H,T), and if COMPZ = 'V', the
!>          orthogonal matrix of right Schur vectors of (A,B).
!>          Not referenced if COMPZ = 'N'.
!> \endverbatim
!>
!> \param[in] LDZ
!> \verbatim
!>          LDZ is INTEGER
!>          The leading dimension of the array Z.  LDZ >= 1.
!>          If COMPZ='V' or 'I', then LDZ >= N.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))
!>          On exit, if INFO >= 0, WORK(1) returns the optimal LWORK.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The dimension of the array WORK.  LWORK >= max(1,N).
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
!>          = 0: successful exit
!>          < 0: if INFO = -i, the i-th argument had an illegal value
!>          = 1,...,N: the QZ iteration did not converge.  (H,T) is not
!>                     in Schur form, but ALPHAR(i), ALPHAI(i), and
!>                     BETA(i), i=INFO+1,...,N should be correct.
!>          = N+1,...,2*N: the shift calculation failed.  (H,T) is not
!>                     in Schur form, but ALPHAR(i), ALPHAI(i), and
!>                     BETA(i), i=INFO-N+1,...,N should be correct.
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
!> \ingroup doubleGEcomputational
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  Iteration counters:
!>
!>  JITER  -- counts iterations.
!>  IITER  -- counts iterations run since ILAST was last
!>            changed.  This is therefore reset only when a 1-by-1 or
!>            2-by-2 block deflates off the bottom.
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE DHGEQZ(Job,Compq,Compz,N,Ilo,Ihi,H,Ldh,T,Ldt,Alphar,   &
     &                  Alphai,Beta,Q,Ldq,Z,Ldz,Work,Lwork,Info)
      USE F77KINDS                        
      USE S_DLAG2
      USE S_DLAMCH
      USE S_DLANHS
      USE S_DLAPY2
      USE S_DLAPY3
      USE S_DLARFG
      USE S_DLARTG
      USE S_DLASET
      USE S_DLASV2
      USE S_DROT
      USE S_LSAME
      USE S_XERBLA
      IMPLICIT NONE
!*--DHGEQZ320
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  HALF = 0.5D+0 , ZERO = 0.0D+0 ,     &
     &                              ONE = 1.0D+0 , SAFETY = 1.0D+2
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Job
      CHARACTER :: Compq
      CHARACTER :: Compz
      INTEGER :: N
      INTEGER , INTENT(IN) :: Ilo
      INTEGER , INTENT(IN) :: Ihi
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldh,*) :: H
      INTEGER :: Ldh
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      REAL(R8KIND) , INTENT(OUT) , DIMENSION(*) :: Alphar
      REAL(R8KIND) , INTENT(OUT) , DIMENSION(*) :: Alphai
      REAL(R8KIND) , INTENT(OUT) , DIMENSION(*) :: Beta
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldq,*) :: Q
      INTEGER :: Ldq
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      REAL(R8KIND) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      REAL(R8KIND) :: a11 , a12 , a1i , a1r , a21 , a22 , a2i , a2r ,   &
     &                ad11 , ad11l , ad12 , ad12l , ad21 , ad21l ,      &
     &                ad22 , ad22l , ad32l , an , anorm , ascale ,      &
     &                atol , b11 , b1a , b1i , b1r , b22 , b2a , b2i ,  &
     &                b2r , bn , bnorm , bscale , btol , c , c11i ,     &
     &                c11r , c12 , c21 , c22i , c22r , cl , cq , cr ,   &
     &                cz , eshift , s , s1 , s1inv , s2 , safmax ,      &
     &                safmin , scale , sl , sqi , sqr , sr , szi , szr ,&
     &                t1 , tau , temp , temp2 , tempi , tempr , u1 ,    &
     &                u12 , u12l , u2 , ulp , vs , w11 , w12 , w21 ,    &
     &                w22 , wabs , wi , wr , wr2
      INTEGER :: icompq , icompz , ifirst , ifrstm , iiter , ilast ,    &
     &           ilastm , in , ischur , istart , j , jc , jch , jiter , &
     &           jr , maxit
      LOGICAL :: ilazr2 , ilazro , ilpivt , ilq , ilschr , ilz , lquery
      REAL(R8KIND) , DIMENSION(3) :: v
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
!    $                     SAFETY = 1.0E+0 )
!     ..
!     .. Local Scalars ..
!     ..
!     .. Local Arrays ..
!     ..
!     .. External Functions ..
!     ..
!     .. External Subroutines ..
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. Executable Statements ..
!
!     Decode JOB, COMPQ, COMPZ
!
      IF ( LSAME(Job,'E') ) THEN
         ilschr = .FALSE.
         ischur = 1
      ELSEIF ( LSAME(Job,'S') ) THEN
         ilschr = .TRUE.
         ischur = 2
      ELSE
         ischur = 0
      ENDIF
!
      IF ( LSAME(Compq,'N') ) THEN
         ilq = .FALSE.
         icompq = 1
      ELSEIF ( LSAME(Compq,'V') ) THEN
         ilq = .TRUE.
         icompq = 2
      ELSEIF ( LSAME(Compq,'I') ) THEN
         ilq = .TRUE.
         icompq = 3
      ELSE
         icompq = 0
      ENDIF
!
      IF ( LSAME(Compz,'N') ) THEN
         ilz = .FALSE.
         icompz = 1
      ELSEIF ( LSAME(Compz,'V') ) THEN
         ilz = .TRUE.
         icompz = 2
      ELSEIF ( LSAME(Compz,'I') ) THEN
         ilz = .TRUE.
         icompz = 3
      ELSE
         icompz = 0
      ENDIF
!
!     Check Argument Values
!
      Info = 0
      Work(1) = MAX(1,N)
      lquery = (Lwork==-1)
      IF ( ischur==0 ) THEN
         Info = -1
      ELSEIF ( icompq==0 ) THEN
         Info = -2
      ELSEIF ( icompz==0 ) THEN
         Info = -3
      ELSEIF ( N<0 ) THEN
         Info = -4
      ELSEIF ( Ilo<1 ) THEN
         Info = -5
      ELSEIF ( Ihi>N .OR. Ihi<Ilo-1 ) THEN
         Info = -6
      ELSEIF ( Ldh<N ) THEN
         Info = -8
      ELSEIF ( Ldt<N ) THEN
         Info = -10
      ELSEIF ( Ldq<1 .OR. (ilq .AND. Ldq<N) ) THEN
         Info = -15
      ELSEIF ( Ldz<1 .OR. (ilz .AND. Ldz<N) ) THEN
         Info = -17
      ELSEIF ( Lwork<MAX(1,N) .AND. .NOT.lquery ) THEN
         Info = -19
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('DHGEQZ',-Info)
         RETURN
      ELSEIF ( lquery ) THEN
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( N<=0 ) THEN
         Work(1) = DBLE(1)
         RETURN
      ENDIF
!
!     Initialize Q and Z
!
      IF ( icompq==3 ) CALL DLASET('Full',N,N,ZERO,ONE,Q,Ldq)
      IF ( icompz==3 ) CALL DLASET('Full',N,N,ZERO,ONE,Z,Ldz)
!
!     Machine Constants
!
      in = Ihi + 1 - Ilo
      safmin = DLAMCH('S')
      safmax = ONE/safmin
      ulp = DLAMCH('E')*DLAMCH('B')
      anorm = DLANHS('F',in,H(Ilo,Ilo),Ldh,Work)
      bnorm = DLANHS('F',in,T(Ilo,Ilo),Ldt,Work)
      atol = MAX(safmin,ulp*anorm)
      btol = MAX(safmin,ulp*bnorm)
      ascale = ONE/MAX(safmin,anorm)
      bscale = ONE/MAX(safmin,bnorm)
!
!     Set Eigenvalues IHI+1:N
!
      DO j = Ihi + 1 , N
         IF ( T(j,j)<ZERO ) THEN
            IF ( ilschr ) THEN
               DO jr = 1 , j
                  H(jr,j) = -H(jr,j)
                  T(jr,j) = -T(jr,j)
               ENDDO
            ELSE
               H(j,j) = -H(j,j)
               T(j,j) = -T(j,j)
            ENDIF
            IF ( ilz ) THEN
               DO jr = 1 , N
                  Z(jr,j) = -Z(jr,j)
               ENDDO
            ENDIF
         ENDIF
         Alphar(j) = H(j,j)
         Alphai(j) = ZERO
         Beta(j) = T(j,j)
      ENDDO
!
!     If IHI < ILO, skip QZ steps
!
      IF ( Ihi>=Ilo ) THEN
!
!     MAIN QZ ITERATION LOOP
!
!     Initialize dynamic indices
!
!     Eigenvalues ILAST+1:N have been found.
!        Column operations modify rows IFRSTM:whatever.
!        Row operations modify columns whatever:ILASTM.
!
!     If only eigenvalues are being computed, then
!        IFRSTM is the row of the last splitting row above row ILAST;
!        this is always at least ILO.
!     IITER counts iterations since the last eigenvalue was found,
!        to tell when to use an extraordinary shift.
!     MAXIT is the maximum number of QZ sweeps allowed.
!
         ilast = Ihi
         IF ( ilschr ) THEN
            ifrstm = 1
            ilastm = N
         ELSE
            ifrstm = Ilo
            ilastm = Ihi
         ENDIF
         iiter = 0
         eshift = ZERO
         maxit = 30*(Ihi-Ilo+1)
!
         DO jiter = 1 , maxit
!
!        Split the matrix if possible.
!
!        Two tests:
!           1: H(j,j-1)=0  or  j=ILO
!           2: T(j,j)=0
!
!
!           Special case: j=ILAST
!
            IF ( ilast==Ilo ) GOTO 40
            IF ( ABS(H(ilast,ilast-1))                                  &
     &           <=MAX(safmin,ulp*(ABS(H(ilast,ilast))                  &
     &           +ABS(H(ilast-1,ilast-1)))) ) THEN
               H(ilast,ilast-1) = ZERO
               GOTO 40
            ENDIF
!
            IF ( ABS(T(ilast,ilast))                                    &
     &           <=MAX(safmin,ulp*(ABS(T(ilast-1,ilast))                &
     &           +ABS(T(ilast-1,ilast-1)))) ) THEN
               T(ilast,ilast) = ZERO
               GOTO 20
            ENDIF
!
!        General case: j<ILAST
!
            DO j = ilast - 1 , Ilo , -1
!
!           Test 1: for H(j,j-1)=0 or j=ILO
!
               IF ( j==Ilo ) THEN
                  ilazro = .TRUE.
               ELSEIF ( ABS(H(j,j-1))                                   &
     &                  <=MAX(safmin,ulp*(ABS(H(j,j))+ABS(H(j-1,j-1)))) &
     &                  ) THEN
                  H(j,j-1) = ZERO
                  ilazro = .TRUE.
               ELSE
                  ilazro = .FALSE.
               ENDIF
!
!           Test 2: for T(j,j)=0
!
               temp = ABS(T(j,j+1))
               IF ( j>Ilo ) temp = temp + ABS(T(j-1,j))
               IF ( ABS(T(j,j))<MAX(safmin,ulp*temp) ) THEN
                  T(j,j) = ZERO
!
!              Test 1a: Check for 2 consecutive small subdiagonals in A
!
                  ilazr2 = .FALSE.
                  IF ( .NOT.ilazro ) THEN
                     temp = ABS(H(j,j-1))
                     temp2 = ABS(H(j,j))
                     tempr = MAX(temp,temp2)
                     IF ( tempr<ONE .AND. tempr/=ZERO ) THEN
                        temp = temp/tempr
                        temp2 = temp2/tempr
                     ENDIF
                     IF ( temp*(ascale*ABS(H(j+1,j)))                   &
     &                    <=temp2*(ascale*atol) ) ilazr2 = .TRUE.
                  ENDIF
!
!              If both tests pass (1 & 2), i.e., the leading diagonal
!              element of B in the block is zero, split a 1x1 block off
!              at the top. (I.e., at the J-th row/column) The leading
!              diagonal element of the remainder can also be zero, so
!              this may have to be done repeatedly.
!
                  IF ( ilazro .OR. ilazr2 ) THEN
                     DO jch = j , ilast - 1
                        temp = H(jch,jch)
                        CALL DLARTG(temp,H(jch+1,jch),c,s,H(jch,jch))
                        H(jch+1,jch) = ZERO
                        CALL DROT(ilastm-jch,H(jch,jch+1),Ldh,          &
     &                            H(jch+1,jch+1),Ldh,c,s)
                        CALL DROT(ilastm-jch,T(jch,jch+1),Ldt,          &
     &                            T(jch+1,jch+1),Ldt,c,s)
                        IF ( ilq ) CALL DROT(N,Q(1,jch),1,Q(1,jch+1),1, &
     &                       c,s)
                        IF ( ilazr2 ) H(jch,jch-1) = H(jch,jch-1)*c
                        ilazr2 = .FALSE.
                        IF ( ABS(T(jch+1,jch+1))>=btol ) THEN
                           IF ( jch+1>=ilast ) GOTO 40
                           ifirst = jch + 1
                           GOTO 60
                        ENDIF
                        T(jch+1,jch+1) = ZERO
                     ENDDO
                  ELSE
!
!                 Only test 2 passed -- chase the zero to T(ILAST,ILAST)
!                 Then process as in the case T(ILAST,ILAST)=0
!
                     DO jch = j , ilast - 1
                        temp = T(jch,jch+1)
                        CALL DLARTG(temp,T(jch+1,jch+1),c,s,T(jch,jch+1)&
     &                              )
                        T(jch+1,jch+1) = ZERO
                        IF ( jch<ilastm-1 )                             &
     &                       CALL DROT(ilastm-jch-1,T(jch,jch+2),Ldt,   &
     &                       T(jch+1,jch+2),Ldt,c,s)
                        CALL DROT(ilastm-jch+2,H(jch,jch-1),Ldh,        &
     &                            H(jch+1,jch-1),Ldh,c,s)
                        IF ( ilq ) CALL DROT(N,Q(1,jch),1,Q(1,jch+1),1, &
     &                       c,s)
                        temp = H(jch+1,jch)
                        CALL DLARTG(temp,H(jch+1,jch-1),c,s,H(jch+1,jch)&
     &                              )
                        H(jch+1,jch-1) = ZERO
                        CALL DROT(jch+1-ifrstm,H(ifrstm,jch),1,         &
     &                            H(ifrstm,jch-1),1,c,s)
                        CALL DROT(jch-ifrstm,T(ifrstm,jch),1,           &
     &                            T(ifrstm,jch-1),1,c,s)
                        IF ( ilz ) CALL DROT(N,Z(1,jch),1,Z(1,jch-1),1, &
     &                       c,s)
                     ENDDO
                  ENDIF
                  GOTO 20
               ELSEIF ( ilazro ) THEN
!
!              Only test 1 passed -- work on J:ILAST
!
                  ifirst = j
                  GOTO 60
               ENDIF
!
!           Neither test passed -- try next J
!
            ENDDO
!
!        (Drop-through is "impossible")
!
            Info = N + 1
            GOTO 200
!
!        T(ILAST,ILAST)=0 -- clear H(ILAST,ILAST-1) to split off a
!        1x1 block.
!
 20         temp = H(ilast,ilast)
            CALL DLARTG(temp,H(ilast,ilast-1),c,s,H(ilast,ilast))
            H(ilast,ilast-1) = ZERO
            CALL DROT(ilast-ifrstm,H(ifrstm,ilast),1,H(ifrstm,ilast-1), &
     &                1,c,s)
            CALL DROT(ilast-ifrstm,T(ifrstm,ilast),1,T(ifrstm,ilast-1), &
     &                1,c,s)
            IF ( ilz ) CALL DROT(N,Z(1,ilast),1,Z(1,ilast-1),1,c,s)
!
!        H(ILAST,ILAST-1)=0 -- Standardize B, set ALPHAR, ALPHAI,
!                              and BETA
!
 40         IF ( T(ilast,ilast)<ZERO ) THEN
               IF ( ilschr ) THEN
                  DO j = ifrstm , ilast
                     H(j,ilast) = -H(j,ilast)
                     T(j,ilast) = -T(j,ilast)
                  ENDDO
               ELSE
                  H(ilast,ilast) = -H(ilast,ilast)
                  T(ilast,ilast) = -T(ilast,ilast)
               ENDIF
               IF ( ilz ) THEN
                  DO j = 1 , N
                     Z(j,ilast) = -Z(j,ilast)
                  ENDDO
               ENDIF
            ENDIF
            Alphar(ilast) = H(ilast,ilast)
            Alphai(ilast) = ZERO
            Beta(ilast) = T(ilast,ilast)
!
!        Go to next block -- exit if finished.
!
            ilast = ilast - 1
            IF ( ilast<Ilo ) GOTO 100
!
!        Reset counters
!
            iiter = 0
            eshift = ZERO
            IF ( .NOT.ilschr ) THEN
               ilastm = ilast
               IF ( ifrstm>ilast ) ifrstm = Ilo
            ENDIF
            CYCLE
!
!        QZ step
!
!        This iteration only involves rows/columns IFIRST:ILAST. We
!        assume IFIRST < ILAST, and that the diagonal of B is non-zero.
!
 60         iiter = iiter + 1
            IF ( .NOT.ilschr ) ifrstm = ifirst
!
!        Compute single shifts.
!
!        At this point, IFIRST < ILAST, and the diagonal elements of
!        T(IFIRST:ILAST,IFIRST,ILAST) are larger than BTOL (in
!        magnitude)
!
            IF ( (iiter/10)*10==iiter ) THEN
!
!           Exceptional shift.  Chosen for no particularly good reason.
!           (Single shift only.)
!
               IF ( (DBLE(maxit)*safmin)*ABS(H(ilast,ilast-1))          &
     &              <ABS(T(ilast-1,ilast-1)) ) THEN
                  eshift = H(ilast,ilast-1)/T(ilast-1,ilast-1)
               ELSE
                  eshift = eshift + ONE/(safmin*DBLE(maxit))
               ENDIF
               s1 = ONE
               wr = eshift
!
            ELSE
!
!           Shifts based on the generalized eigenvalues of the
!           bottom-right 2x2 block of A and B. The first eigenvalue
!           returned by DLAG2 is the Wilkinson shift (AEP p.512),
!
               CALL DLAG2(H(ilast-1,ilast-1),Ldh,T(ilast-1,ilast-1),Ldt,&
     &                    safmin*SAFETY,s1,s2,wr,wr2,wi)
!
               IF ( ABS((wr/s1)*T(ilast,ilast)-H(ilast,ilast))          &
     &              >ABS((wr2/s2)*T(ilast,ilast)-H(ilast,ilast)) ) THEN
                  temp = wr
                  wr = wr2
                  wr2 = temp
                  temp = s1
                  s1 = s2
                  s2 = temp
               ENDIF
               temp = MAX(s1,safmin*MAX(ONE,ABS(wr),ABS(wi)))
               IF ( wi/=ZERO ) THEN
!
!        Use Francis double-shift
!
!        Note: the Francis double-shift should work with real shifts,
!              but only if the block is at least 3x3.
!              This code may break if this point is reached with
!              a 2x2 block with real eigenvalues.
!
                  IF ( ifirst+1==ilast ) THEN
!
!           Special case -- 2x2 block with complex eigenvectors
!
!           Step 1: Standardize, that is, rotate so that
!
!                       ( B11  0  )
!                   B = (         )  with B11 non-negative.
!                       (  0  B22 )
!
                     CALL DLASV2(T(ilast-1,ilast-1),T(ilast-1,ilast),   &
     &                           T(ilast,ilast),b22,b11,sr,cr,sl,cl)
!
                     IF ( b11<ZERO ) THEN
                        cr = -cr
                        sr = -sr
                        b11 = -b11
                        b22 = -b22
                     ENDIF
!
                     CALL DROT(ilastm+1-ifirst,H(ilast-1,ilast-1),Ldh,  &
     &                         H(ilast,ilast-1),Ldh,cl,sl)
                     CALL DROT(ilast+1-ifrstm,H(ifrstm,ilast-1),1,      &
     &                         H(ifrstm,ilast),1,cr,sr)
!
                     IF ( ilast<ilastm )                                &
     &                    CALL DROT(ilastm-ilast,T(ilast-1,ilast+1),Ldt,&
     &                    T(ilast,ilast+1),Ldt,cl,sl)
                     IF ( ifrstm<ilast-1 )                              &
     &                    CALL DROT(ifirst-ifrstm,T(ifrstm,ilast-1),1,  &
     &                    T(ifrstm,ilast),1,cr,sr)
!
                     IF ( ilq ) CALL DROT(N,Q(1,ilast-1),1,Q(1,ilast),1,&
     &                    cl,sl)
                     IF ( ilz ) CALL DROT(N,Z(1,ilast-1),1,Z(1,ilast),1,&
     &                    cr,sr)
!
                     T(ilast-1,ilast-1) = b11
                     T(ilast-1,ilast) = ZERO
                     T(ilast,ilast-1) = ZERO
                     T(ilast,ilast) = b22
!
!           If B22 is negative, negate column ILAST
!
                     IF ( b22<ZERO ) THEN
                        DO j = ifrstm , ilast
                           H(j,ilast) = -H(j,ilast)
                           T(j,ilast) = -T(j,ilast)
                        ENDDO
!
                        IF ( ilz ) THEN
                           DO j = 1 , N
                              Z(j,ilast) = -Z(j,ilast)
                           ENDDO
                        ENDIF
                        b22 = -b22
                     ENDIF
!
!           Step 2: Compute ALPHAR, ALPHAI, and BETA (see refs.)
!
!           Recompute shift
!
                     CALL DLAG2(H(ilast-1,ilast-1),Ldh,                 &
     &                          T(ilast-1,ilast-1),Ldt,safmin*SAFETY,s1,&
     &                          temp,wr,temp2,wi)
!
!           If standardization has perturbed the shift onto real line,
!           do another (real single-shift) QR step.
!
                     IF ( wi/=ZERO ) THEN
                        s1inv = ONE/s1
!
!           Do EISPACK (QZVAL) computation of alpha and beta
!
                        a11 = H(ilast-1,ilast-1)
                        a21 = H(ilast,ilast-1)
                        a12 = H(ilast-1,ilast)
                        a22 = H(ilast,ilast)
!
!           Compute complex Givens rotation on right
!           (Assume some element of C = (sA - wB) > unfl )
!                            __
!           (sA - wB) ( CZ   -SZ )
!                     ( SZ    CZ )
!
                        c11r = s1*a11 - wr*b11
                        c11i = -wi*b11
                        c12 = s1*a12
                        c21 = s1*a21
                        c22r = s1*a22 - wr*b22
                        c22i = -wi*b22
!
                        IF ( ABS(c11r)+ABS(c11i)+ABS(c12)>ABS(c21)      &
     &                       +ABS(c22r)+ABS(c22i) ) THEN
                           t1 = DLAPY3(c12,c11r,c11i)
                           cz = c12/t1
                           szr = -c11r/t1
                           szi = -c11i/t1
                        ELSE
                           cz = DLAPY2(c22r,c22i)
                           IF ( cz<=safmin ) THEN
                              cz = ZERO
                              szr = ONE
                              szi = ZERO
                           ELSE
                              tempr = c22r/cz
                              tempi = c22i/cz
                              t1 = DLAPY2(cz,c21)
                              cz = cz/t1
                              szr = -c21*tempr/t1
                              szi = c21*tempi/t1
                           ENDIF
                        ENDIF
!
!           Compute Givens rotation on left
!
!           (  CQ   SQ )
!           (  __      )  A or B
!           ( -SQ   CQ )
!
                        an = ABS(a11) + ABS(a12) + ABS(a21) + ABS(a22)
                        bn = ABS(b11) + ABS(b22)
                        wabs = ABS(wr) + ABS(wi)
                        IF ( s1*an>wabs*bn ) THEN
                           cq = cz*b11
                           sqr = szr*b22
                           sqi = -szi*b22
                        ELSE
                           a1r = cz*a11 + szr*a12
                           a1i = szi*a12
                           a2r = cz*a21 + szr*a22
                           a2i = szi*a22
                           cq = DLAPY2(a1r,a1i)
                           IF ( cq<=safmin ) THEN
                              cq = ZERO
                              sqr = ONE
                              sqi = ZERO
                           ELSE
                              tempr = a1r/cq
                              tempi = a1i/cq
                              sqr = tempr*a2r + tempi*a2i
                              sqi = tempi*a2r - tempr*a2i
                           ENDIF
                        ENDIF
                        t1 = DLAPY3(cq,sqr,sqi)
                        cq = cq/t1
                        sqr = sqr/t1
                        sqi = sqi/t1
!
!           Compute diagonal elements of QBZ
!
                        tempr = sqr*szr - sqi*szi
                        tempi = sqr*szi + sqi*szr
                        b1r = cq*cz*b11 + tempr*b22
                        b1i = tempi*b22
                        b1a = DLAPY2(b1r,b1i)
                        b2r = cq*cz*b22 + tempr*b11
                        b2i = -tempi*b11
                        b2a = DLAPY2(b2r,b2i)
!
!           Normalize so beta > 0, and Im( alpha1 ) > 0
!
                        Beta(ilast-1) = b1a
                        Beta(ilast) = b2a
                        Alphar(ilast-1) = (wr*b1a)*s1inv
                        Alphai(ilast-1) = (wi*b1a)*s1inv
                        Alphar(ilast) = (wr*b2a)*s1inv
                        Alphai(ilast) = -(wi*b2a)*s1inv
!
!           Step 3: Go to next block -- exit if finished.
!
                        ilast = ifirst - 1
                        IF ( ilast<Ilo ) GOTO 100
!
!           Reset counters
!
                        iiter = 0
                        eshift = ZERO
                        IF ( .NOT.ilschr ) THEN
                           ilastm = ilast
                           IF ( ifrstm>ilast ) ifrstm = Ilo
                        ENDIF
                     ENDIF
                  ELSE
!
!           Usual case: 3x3 or larger block, using Francis implicit
!                       double-shift
!
!                                    2
!           Eigenvalue equation is  w  - c w + d = 0,
!
!                                         -1 2        -1
!           so compute 1st column of  (A B  )  - c A B   + d
!           using the formula in QZIT (from EISPACK)
!
!           We assume that the block is at least 3x3
!
                     ad11 = (ascale*H(ilast-1,ilast-1))                 &
     &                      /(bscale*T(ilast-1,ilast-1))
                     ad21 = (ascale*H(ilast,ilast-1))                   &
     &                      /(bscale*T(ilast-1,ilast-1))
                     ad12 = (ascale*H(ilast-1,ilast))                   &
     &                      /(bscale*T(ilast,ilast))
                     ad22 = (ascale*H(ilast,ilast))                     &
     &                      /(bscale*T(ilast,ilast))
                     u12 = T(ilast-1,ilast)/T(ilast,ilast)
                     ad11l = (ascale*H(ifirst,ifirst))                  &
     &                       /(bscale*T(ifirst,ifirst))
                     ad21l = (ascale*H(ifirst+1,ifirst))                &
     &                       /(bscale*T(ifirst,ifirst))
                     ad12l = (ascale*H(ifirst,ifirst+1))                &
     &                       /(bscale*T(ifirst+1,ifirst+1))
                     ad22l = (ascale*H(ifirst+1,ifirst+1))              &
     &                       /(bscale*T(ifirst+1,ifirst+1))
                     ad32l = (ascale*H(ifirst+2,ifirst+1))              &
     &                       /(bscale*T(ifirst+1,ifirst+1))
                     u12l = T(ifirst,ifirst+1)/T(ifirst+1,ifirst+1)
!
                     v(1) = (ad11-ad11l)*(ad22-ad11l) - ad12*ad21 +     &
     &                      ad21*u12*ad11l + (ad12l-ad11l*u12l)*ad21l
                     v(2) = ((ad22l-ad11l)-ad21l*u12l-(ad11-ad11l)      &
     &                      -(ad22-ad11l)+ad21*u12)*ad21l
                     v(3) = ad32l*ad21l
!
                     istart = ifirst
!
                     CALL DLARFG(3,v(1),v(2),1,tau)
                     v(1) = ONE
!
!           Sweep
!
                     DO j = istart , ilast - 2
!
!              All but last elements: use 3x3 Householder transforms.
!
!              Zero (j-1)st column of A
!
                        IF ( j>istart ) THEN
                           v(1) = H(j,j-1)
                           v(2) = H(j+1,j-1)
                           v(3) = H(j+2,j-1)
!
                           CALL DLARFG(3,H(j,j-1),v(2),1,tau)
                           v(1) = ONE
                           H(j+1,j-1) = ZERO
                           H(j+2,j-1) = ZERO
                        ENDIF
!
                        DO jc = j , ilastm
                           temp = tau*(H(j,jc)+v(2)*H(j+1,jc)+v(3)      &
     &                            *H(j+2,jc))
                           H(j,jc) = H(j,jc) - temp
                           H(j+1,jc) = H(j+1,jc) - temp*v(2)
                           H(j+2,jc) = H(j+2,jc) - temp*v(3)
                           temp2 = tau*(T(j,jc)+v(2)*T(j+1,jc)+v(3)     &
     &                             *T(j+2,jc))
                           T(j,jc) = T(j,jc) - temp2
                           T(j+1,jc) = T(j+1,jc) - temp2*v(2)
                           T(j+2,jc) = T(j+2,jc) - temp2*v(3)
                        ENDDO
                        IF ( ilq ) THEN
                           DO jr = 1 , N
                              temp = tau*(Q(jr,j)+v(2)*Q(jr,j+1)+v(3)   &
     &                               *Q(jr,j+2))
                              Q(jr,j) = Q(jr,j) - temp
                              Q(jr,j+1) = Q(jr,j+1) - temp*v(2)
                              Q(jr,j+2) = Q(jr,j+2) - temp*v(3)
                           ENDDO
                        ENDIF
!
!              Zero j-th column of B (see DLAGBC for details)
!
!              Swap rows to pivot
!
                        ilpivt = .FALSE.
                        temp = MAX(ABS(T(j+1,j+1)),ABS(T(j+1,j+2)))
                        temp2 = MAX(ABS(T(j+2,j+1)),ABS(T(j+2,j+2)))
                        IF ( MAX(temp,temp2)<safmin ) THEN
                           scale = ZERO
                           u1 = ONE
                           u2 = ZERO
                           GOTO 62
                        ELSEIF ( temp>=temp2 ) THEN
                           w11 = T(j+1,j+1)
                           w21 = T(j+2,j+1)
                           w12 = T(j+1,j+2)
                           w22 = T(j+2,j+2)
                           u1 = T(j+1,j)
                           u2 = T(j+2,j)
                        ELSE
                           w21 = T(j+1,j+1)
                           w11 = T(j+2,j+1)
                           w22 = T(j+1,j+2)
                           w12 = T(j+2,j+2)
                           u2 = T(j+1,j)
                           u1 = T(j+2,j)
                        ENDIF
!
!              Swap columns if nec.
!
                        IF ( ABS(w12)>ABS(w11) ) THEN
                           ilpivt = .TRUE.
                           temp = w12
                           temp2 = w22
                           w12 = w11
                           w22 = w21
                           w11 = temp
                           w21 = temp2
                        ENDIF
!
!              LU-factor
!
                        temp = w21/w11
                        u2 = u2 - temp*u1
                        w22 = w22 - temp*w12
                        w21 = ZERO
!
!              Compute SCALE
!
                        scale = ONE
                        IF ( ABS(w22)<safmin ) THEN
                           scale = ZERO
                           u2 = ONE
                           u1 = -w12/w11
                           GOTO 62
                        ENDIF
                        IF ( ABS(w22)<ABS(u2) ) scale = ABS(w22/u2)
                        IF ( ABS(w11)<ABS(u1) )                         &
     &                       scale = MIN(scale,ABS(w11/u1))
!
!              Solve
!
                        u2 = (scale*u2)/w22
                        u1 = (scale*u1-w12*u2)/w11
!
 62                     IF ( ilpivt ) THEN
                           temp = u2
                           u2 = u1
                           u1 = temp
                        ENDIF
!
!              Compute Householder Vector
!
                        t1 = SQRT(scale**2+u1**2+u2**2)
                        tau = ONE + scale/t1
                        vs = -ONE/(scale+t1)
                        v(1) = ONE
                        v(2) = vs*u1
                        v(3) = vs*u2
!
!              Apply transformations from the right.
!
                        DO jr = ifrstm , MIN(j+3,ilast)
                           temp = tau*(H(jr,j)+v(2)*H(jr,j+1)+v(3)      &
     &                            *H(jr,j+2))
                           H(jr,j) = H(jr,j) - temp
                           H(jr,j+1) = H(jr,j+1) - temp*v(2)
                           H(jr,j+2) = H(jr,j+2) - temp*v(3)
                        ENDDO
                        DO jr = ifrstm , j + 2
                           temp = tau*(T(jr,j)+v(2)*T(jr,j+1)+v(3)      &
     &                            *T(jr,j+2))
                           T(jr,j) = T(jr,j) - temp
                           T(jr,j+1) = T(jr,j+1) - temp*v(2)
                           T(jr,j+2) = T(jr,j+2) - temp*v(3)
                        ENDDO
                        IF ( ilz ) THEN
                           DO jr = 1 , N
                              temp = tau*(Z(jr,j)+v(2)*Z(jr,j+1)+v(3)   &
     &                               *Z(jr,j+2))
                              Z(jr,j) = Z(jr,j) - temp
                              Z(jr,j+1) = Z(jr,j+1) - temp*v(2)
                              Z(jr,j+2) = Z(jr,j+2) - temp*v(3)
                           ENDDO
                        ENDIF
                        T(j+1,j) = ZERO
                        T(j+2,j) = ZERO
                     ENDDO
!
!           Last elements: Use Givens rotations
!
!           Rotations from the left
!
                     j = ilast - 1
                     temp = H(j,j-1)
                     CALL DLARTG(temp,H(j+1,j-1),c,s,H(j,j-1))
                     H(j+1,j-1) = ZERO
!
                     DO jc = j , ilastm
                        temp = c*H(j,jc) + s*H(j+1,jc)
                        H(j+1,jc) = -s*H(j,jc) + c*H(j+1,jc)
                        H(j,jc) = temp
                        temp2 = c*T(j,jc) + s*T(j+1,jc)
                        T(j+1,jc) = -s*T(j,jc) + c*T(j+1,jc)
                        T(j,jc) = temp2
                     ENDDO
                     IF ( ilq ) THEN
                        DO jr = 1 , N
                           temp = c*Q(jr,j) + s*Q(jr,j+1)
                           Q(jr,j+1) = -s*Q(jr,j) + c*Q(jr,j+1)
                           Q(jr,j) = temp
                        ENDDO
                     ENDIF
!
!           Rotations from the right.
!
                     temp = T(j+1,j+1)
                     CALL DLARTG(temp,T(j+1,j),c,s,T(j+1,j+1))
                     T(j+1,j) = ZERO
!
                     DO jr = ifrstm , ilast
                        temp = c*H(jr,j+1) + s*H(jr,j)
                        H(jr,j) = -s*H(jr,j+1) + c*H(jr,j)
                        H(jr,j+1) = temp
                     ENDDO
                     DO jr = ifrstm , ilast - 1
                        temp = c*T(jr,j+1) + s*T(jr,j)
                        T(jr,j) = -s*T(jr,j+1) + c*T(jr,j)
                        T(jr,j+1) = temp
                     ENDDO
                     IF ( ilz ) THEN
                        DO jr = 1 , N
                           temp = c*Z(jr,j+1) + s*Z(jr,j)
                           Z(jr,j) = -s*Z(jr,j+1) + c*Z(jr,j)
                           Z(jr,j+1) = temp
                        ENDDO
                     ENDIF
!
!           End of Double-Shift code
!
                  ENDIF
!
                  CYCLE
               ENDIF
            ENDIF
!
!        Fiddle with shift to avoid overflow
!
            temp = MIN(ascale,ONE)*(HALF*safmax)
            IF ( s1>temp ) THEN
               scale = temp/s1
            ELSE
               scale = ONE
            ENDIF
!
            temp = MIN(bscale,ONE)*(HALF*safmax)
            IF ( ABS(wr)>temp ) scale = MIN(scale,temp/ABS(wr))
            s1 = scale*s1
            wr = scale*wr
!
!        Now check for two consecutive small subdiagonals.
!
            DO j = ilast - 1 , ifirst + 1 , -1
               istart = j
               temp = ABS(s1*H(j,j-1))
               temp2 = ABS(s1*H(j,j)-wr*T(j,j))
               tempr = MAX(temp,temp2)
               IF ( tempr<ONE .AND. tempr/=ZERO ) THEN
                  temp = temp/tempr
                  temp2 = temp2/tempr
               ENDIF
               IF ( ABS((ascale*H(j+1,j))*temp)<=(ascale*atol)*temp2 )  &
     &              GOTO 80
            ENDDO
!
            istart = ifirst
!
!        Do an implicit single-shift QZ sweep.
!
!        Initial Q
!
 80         temp = s1*H(istart,istart) - wr*T(istart,istart)
            temp2 = s1*H(istart+1,istart)
            CALL DLARTG(temp,temp2,c,s,tempr)
!
!        Sweep
!
            DO j = istart , ilast - 1
               IF ( j>istart ) THEN
                  temp = H(j,j-1)
                  CALL DLARTG(temp,H(j+1,j-1),c,s,H(j,j-1))
                  H(j+1,j-1) = ZERO
               ENDIF
!
               DO jc = j , ilastm
                  temp = c*H(j,jc) + s*H(j+1,jc)
                  H(j+1,jc) = -s*H(j,jc) + c*H(j+1,jc)
                  H(j,jc) = temp
                  temp2 = c*T(j,jc) + s*T(j+1,jc)
                  T(j+1,jc) = -s*T(j,jc) + c*T(j+1,jc)
                  T(j,jc) = temp2
               ENDDO
               IF ( ilq ) THEN
                  DO jr = 1 , N
                     temp = c*Q(jr,j) + s*Q(jr,j+1)
                     Q(jr,j+1) = -s*Q(jr,j) + c*Q(jr,j+1)
                     Q(jr,j) = temp
                  ENDDO
               ENDIF
!
               temp = T(j+1,j+1)
               CALL DLARTG(temp,T(j+1,j),c,s,T(j+1,j+1))
               T(j+1,j) = ZERO
!
               DO jr = ifrstm , MIN(j+2,ilast)
                  temp = c*H(jr,j+1) + s*H(jr,j)
                  H(jr,j) = -s*H(jr,j+1) + c*H(jr,j)
                  H(jr,j+1) = temp
               ENDDO
               DO jr = ifrstm , j
                  temp = c*T(jr,j+1) + s*T(jr,j)
                  T(jr,j) = -s*T(jr,j+1) + c*T(jr,j)
                  T(jr,j+1) = temp
               ENDDO
               IF ( ilz ) THEN
                  DO jr = 1 , N
                     temp = c*Z(jr,j+1) + s*Z(jr,j)
                     Z(jr,j) = -s*Z(jr,j+1) + c*Z(jr,j)
                     Z(jr,j+1) = temp
                  ENDDO
               ENDIF
!
            ENDDO
!
!        End of iteration loop
!
         ENDDO
!
!     Drop-through = non-convergence
!
         Info = ilast
         GOTO 200
      ENDIF
!
!     Successful completion of all QZ steps
!
!
!     Set Eigenvalues 1:ILO-1
!
 100  DO j = 1 , Ilo - 1
         IF ( T(j,j)<ZERO ) THEN
            IF ( ilschr ) THEN
               DO jr = 1 , j
                  H(jr,j) = -H(jr,j)
                  T(jr,j) = -T(jr,j)
               ENDDO
            ELSE
               H(j,j) = -H(j,j)
               T(j,j) = -T(j,j)
            ENDIF
            IF ( ilz ) THEN
               DO jr = 1 , N
                  Z(jr,j) = -Z(jr,j)
               ENDDO
            ENDIF
         ENDIF
         Alphar(j) = H(j,j)
         Alphai(j) = ZERO
         Beta(j) = T(j,j)
      ENDDO
!
!     Normal Termination
!
      Info = 0
!
!     Exit (other than argument error) -- return optimal workspace size
!
 200  Work(1) = DBLE(N)
!
!     End of DHGEQZ
!
      END SUBROUTINE DHGEQZ
