!*==zhgeqz.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b ZHGEQZ
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZHGEQZ + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhgeqz.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhgeqz.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhgeqz.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZHGEQZ( JOB, COMPQ, COMPZ, N, ILO, IHI, H, LDH, T, LDT,
!                          ALPHA, BETA, Q, LDQ, Z, LDZ, WORK, LWORK,
!                          RWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          COMPQ, COMPZ, JOB
!       INTEGER            IHI, ILO, INFO, LDH, LDQ, LDT, LDZ, LWORK, N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   RWORK( * )
!       COMPLEX*16         ALPHA( * ), BETA( * ), H( LDH, * ),
!      $                   Q( LDQ, * ), T( LDT, * ), WORK( * ),
!      $                   Z( LDZ, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZHGEQZ computes the eigenvalues of a complex matrix pair (H,T),
!> where H is an upper Hessenberg matrix and T is upper triangular,
!> using the single-shift QZ method.
!> Matrix pairs of this type are produced by the reduction to
!> generalized upper Hessenberg form of a complex matrix pair (A,B):
!>
!>    A = Q1*H*Z1**H,  B = Q1*T*Z1**H,
!>
!> as computed by ZGGHRD.
!>
!> If JOB='S', then the Hessenberg-triangular pair (H,T) is
!> also reduced to generalized Schur form,
!>
!>    H = Q*S*Z**H,  T = Q*P*Z**H,
!>
!> where Q and Z are unitary matrices and S and P are upper triangular.
!>
!> Optionally, the unitary matrix Q from the generalized Schur
!> factorization may be postmultiplied into an input matrix Q1, and the
!> unitary matrix Z may be postmultiplied into an input matrix Z1.
!> If Q1 and Z1 are the unitary matrices from ZGGHRD that reduced
!> the matrix pair (A,B) to generalized Hessenberg form, then the output
!> matrices Q1*Q and Z1*Z are the unitary factors from the generalized
!> Schur factorization of (A,B):
!>
!>    A = (Q1*Q)*S*(Z1*Z)**H,  B = (Q1*Q)*P*(Z1*Z)**H.
!>
!> To avoid overflow, eigenvalues of the matrix pair (H,T)
!> (equivalently, of (A,B)) are computed as a pair of complex values
!> (alpha,beta).  If beta is nonzero, lambda = alpha / beta is an
!> eigenvalue of the generalized nonsymmetric eigenvalue problem (GNEP)
!>    A*x = lambda*B*x
!> and if alpha is nonzero, mu = beta / alpha is an eigenvalue of the
!> alternate form of the GNEP
!>    mu*A*y = B*y.
!> The values of alpha and beta for the i-th eigenvalue can be read
!> directly from the generalized Schur form:  alpha = S(i,i),
!> beta = P(i,i).
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
!>          = 'S': Computer eigenvalues and the Schur form.
!> \endverbatim
!>
!> \param[in] COMPQ
!> \verbatim
!>          COMPQ is CHARACTER*1
!>          = 'N': Left Schur vectors (Q) are not computed;
!>          = 'I': Q is initialized to the unit matrix and the matrix Q
!>                 of left Schur vectors of (H,T) is returned;
!>          = 'V': Q must contain a unitary matrix Q1 on entry and
!>                 the product Q1*Q is returned.
!> \endverbatim
!>
!> \param[in] COMPZ
!> \verbatim
!>          COMPZ is CHARACTER*1
!>          = 'N': Right Schur vectors (Z) are not computed;
!>          = 'I': Q is initialized to the unit matrix and the matrix Z
!>                 of right Schur vectors of (H,T) is returned;
!>          = 'V': Z must contain a unitary matrix Z1 on entry and
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
!>          H is COMPLEX*16 array, dimension (LDH, N)
!>          On entry, the N-by-N upper Hessenberg matrix H.
!>          On exit, if JOB = 'S', H contains the upper triangular
!>          matrix S from the generalized Schur factorization.
!>          If JOB = 'E', the diagonal of H matches that of S, but
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
!>          T is COMPLEX*16 array, dimension (LDT, N)
!>          On entry, the N-by-N upper triangular matrix T.
!>          On exit, if JOB = 'S', T contains the upper triangular
!>          matrix P from the generalized Schur factorization.
!>          If JOB = 'E', the diagonal of T matches that of P, but
!>          the rest of T is unspecified.
!> \endverbatim
!>
!> \param[in] LDT
!> \verbatim
!>          LDT is INTEGER
!>          The leading dimension of the array T.  LDT >= max( 1, N ).
!> \endverbatim
!>
!> \param[out] ALPHA
!> \verbatim
!>          ALPHA is COMPLEX*16 array, dimension (N)
!>          The complex scalars alpha that define the eigenvalues of
!>          GNEP.  ALPHA(i) = S(i,i) in the generalized Schur
!>          factorization.
!> \endverbatim
!>
!> \param[out] BETA
!> \verbatim
!>          BETA is COMPLEX*16 array, dimension (N)
!>          The real non-negative scalars beta that define the
!>          eigenvalues of GNEP.  BETA(i) = P(i,i) in the generalized
!>          Schur factorization.
!>
!>          Together, the quantities alpha = ALPHA(j) and beta = BETA(j)
!>          represent the j-th eigenvalue of the matrix pair (A,B), in
!>          one of the forms lambda = alpha/beta or mu = beta/alpha.
!>          Since either lambda or mu may overflow, they should not,
!>          in general, be computed.
!> \endverbatim
!>
!> \param[in,out] Q
!> \verbatim
!>          Q is COMPLEX*16 array, dimension (LDQ, N)
!>          On entry, if COMPQ = 'V', the unitary matrix Q1 used in the
!>          reduction of (A,B) to generalized Hessenberg form.
!>          On exit, if COMPQ = 'I', the unitary matrix of left Schur
!>          vectors of (H,T), and if COMPQ = 'V', the unitary matrix of
!>          left Schur vectors of (A,B).
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
!>          Z is COMPLEX*16 array, dimension (LDZ, N)
!>          On entry, if COMPZ = 'V', the unitary matrix Z1 used in the
!>          reduction of (A,B) to generalized Hessenberg form.
!>          On exit, if COMPZ = 'I', the unitary matrix of right Schur
!>          vectors of (H,T), and if COMPZ = 'V', the unitary matrix of
!>          right Schur vectors of (A,B).
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
!>          WORK is COMPLEX*16 array, dimension (MAX(1,LWORK))
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
!> \param[out] RWORK
!> \verbatim
!>          RWORK is DOUBLE PRECISION array, dimension (N)
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0: successful exit
!>          < 0: if INFO = -i, the i-th argument had an illegal value
!>          = 1,...,N: the QZ iteration did not converge.  (H,T) is not
!>                     in Schur form, but ALPHA(i) and BETA(i),
!>                     i=INFO+1,...,N should be correct.
!>          = N+1,...,2*N: the shift calculation failed.  (H,T) is not
!>                     in Schur form, but ALPHA(i) and BETA(i),
!>                     i=INFO-N+1,...,N should be correct.
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
!> \date April 2012
!
!> \ingroup complex16GEcomputational
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  We assume that complex ABS works as long as its value is less than
!>  overflow.
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE ZHGEQZ(Job,Compq,Compz,N,Ilo,Ihi,H,Ldh,T,Ldt,Alpha,    &
     &                  Beta,Q,Ldq,Z,Ldz,Work,Lwork,Rwork,Info)
      USE F77KINDS                        
      USE S_DLAMCH
      USE S_LSAME
      USE S_XERBLA
      USE S_ZLADIV
      USE S_ZLANHS
      USE S_ZLARTG
      USE S_ZLASET
      USE S_ZROT
      USE S_ZSCAL
      IMPLICIT NONE
!*--ZHGEQZ297
!
! PARAMETER definitions rewritten by SPAG
!
      COMPLEX(CX16KIND) , PARAMETER  ::  CZERO = (0.0D+0,0.0D+0) ,      &
     &                 CONE = (1.0D+0,0.0D+0)
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0 ,      &
     &                              HALF = 0.5D+0
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Job
      CHARACTER :: Compq
      CHARACTER :: Compz
      INTEGER :: N
      INTEGER , INTENT(IN) :: Ilo
      INTEGER , INTENT(IN) :: Ihi
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldh,*) :: H
      INTEGER :: Ldh
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      COMPLEX(CX16KIND) , INTENT(OUT) , DIMENSION(*) :: Alpha
      COMPLEX(CX16KIND) , INTENT(OUT) , DIMENSION(*) :: Beta
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldq,*) :: Q
      INTEGER :: Ldq
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      COMPLEX(CX16KIND) , INTENT(OUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      REAL(R8KIND) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      COMPLEX(CX16KIND) :: abi12 , abi22 , ad11 , ad12 , ad21 , ad22 ,  &
     &                     ctemp , ctemp2 , ctemp3 , eshift , s ,       &
     &                     shift , signbc , u12 , x , y
      REAL(R8KIND) :: ABS1
      REAL(R8KIND) :: absb , anorm , ascale , atol , bnorm , bscale ,   &
     &                btol , c , safmin , temp , temp2 , tempr , ulp
      INTEGER :: icompq , icompz , ifirst , ifrstm , iiter , ilast ,    &
     &           ilastm , in , ischur , istart , j , jc , jch , jiter , &
     &           jr , maxit
      LOGICAL :: ilazr2 , ilazro , ilq , ilschr , ilz , lquery
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
!     .. External Functions ..
!     ..
!     .. External Subroutines ..
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. Statement Functions ..
!     ..
!     .. Statement Function definitions ..
      ABS1(x) = ABS(DBLE(x)) + ABS(DIMAG(x))
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
         ilschr = .TRUE.
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
         ilq = .TRUE.
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
         ilz = .TRUE.
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
         Info = -14
      ELSEIF ( Ldz<1 .OR. (ilz .AND. Ldz<N) ) THEN
         Info = -16
      ELSEIF ( Lwork<MAX(1,N) .AND. .NOT.lquery ) THEN
         Info = -18
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('ZHGEQZ',-Info)
         RETURN
      ELSEIF ( lquery ) THEN
         RETURN
      ENDIF
!
!     Quick return if possible
!
!     WORK( 1 ) = CMPLX( 1 )
      IF ( N<=0 ) THEN
         Work(1) = DCMPLX(1)
         RETURN
      ENDIF
!
!     Initialize Q and Z
!
      IF ( icompq==3 ) CALL ZLASET('Full',N,N,CZERO,CONE,Q,Ldq)
      IF ( icompz==3 ) CALL ZLASET('Full',N,N,CZERO,CONE,Z,Ldz)
!
!     Machine Constants
!
      in = Ihi + 1 - Ilo
      safmin = DLAMCH('S')
      ulp = DLAMCH('E')*DLAMCH('B')
      anorm = ZLANHS('F',in,H(Ilo,Ilo),Ldh,Rwork)
      bnorm = ZLANHS('F',in,T(Ilo,Ilo),Ldt,Rwork)
      atol = MAX(safmin,ulp*anorm)
      btol = MAX(safmin,ulp*bnorm)
      ascale = ONE/MAX(safmin,anorm)
      bscale = ONE/MAX(safmin,bnorm)
!
!
!     Set Eigenvalues IHI+1:N
!
      DO j = Ihi + 1 , N
         absb = ABS(T(j,j))
         IF ( absb>safmin ) THEN
            signbc = DCONJG(T(j,j)/absb)
            T(j,j) = absb
            IF ( ilschr ) THEN
               CALL ZSCAL(j-1,signbc,T(1,j),1)
               CALL ZSCAL(j,signbc,H(1,j),1)
            ELSE
               CALL ZSCAL(1,signbc,H(j,j),1)
            ENDIF
            IF ( ilz ) CALL ZSCAL(N,signbc,Z(1,j),1)
         ELSE
            T(j,j) = CZERO
         ENDIF
         Alpha(j) = H(j,j)
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
!        Column operations modify rows IFRSTM:whatever
!        Row operations modify columns whatever:ILASTM
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
         eshift = CZERO
         maxit = 30*(Ihi-Ilo+1)
!
         DO jiter = 1 , maxit
!
!        Check for too many iterations.
!
            IF ( jiter>maxit ) EXIT
!
!        Split the matrix if possible.
!
!        Two tests:
!           1: H(j,j-1)=0  or  j=ILO
!           2: T(j,j)=0
!
!        Special case: j=ILAST
!
            IF ( ilast==Ilo ) GOTO 40
            IF ( ABS1(H(ilast,ilast-1))                                 &
     &           <=MAX(safmin,ulp*(ABS1(H(ilast,ilast))                 &
     &           +ABS1(H(ilast-1,ilast-1)))) ) THEN
               H(ilast,ilast-1) = CZERO
               GOTO 40
            ENDIF
!
            IF ( ABS(T(ilast,ilast))                                    &
     &           <=MAX(safmin,ulp*(ABS(T(ilast-1,ilast))                &
     &           +ABS(T(ilast-1,ilast-1)))) ) THEN
               T(ilast,ilast) = CZERO
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
               ELSEIF ( ABS1(H(j,j-1))                                  &
     &                  <=MAX(safmin,ulp*(ABS1(H(j,j))+ABS1(H(j-1,j-1)))&
     &                  ) ) THEN
                  H(j,j-1) = CZERO
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
                  T(j,j) = CZERO
!
!              Test 1a: Check for 2 consecutive small subdiagonals in A
!
                  ilazr2 = .FALSE.
                  IF ( .NOT.ilazro ) THEN
                     IF ( ABS1(H(j,j-1))*(ascale*ABS1(H(j+1,j)))        &
     &                    <=ABS1(H(j,j))*(ascale*atol) ) ilazr2 = .TRUE.
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
                        ctemp = H(jch,jch)
                        CALL ZLARTG(ctemp,H(jch+1,jch),c,s,H(jch,jch))
                        H(jch+1,jch) = CZERO
                        CALL ZROT(ilastm-jch,H(jch,jch+1),Ldh,          &
     &                            H(jch+1,jch+1),Ldh,c,s)
                        CALL ZROT(ilastm-jch,T(jch,jch+1),Ldt,          &
     &                            T(jch+1,jch+1),Ldt,c,s)
                        IF ( ilq ) CALL ZROT(N,Q(1,jch),1,Q(1,jch+1),1, &
     &                       c,DCONJG(s))
                        IF ( ilazr2 ) H(jch,jch-1) = H(jch,jch-1)*c
                        ilazr2 = .FALSE.
                        IF ( ABS1(T(jch+1,jch+1))>=btol ) THEN
                           IF ( jch+1>=ilast ) GOTO 40
                           ifirst = jch + 1
                           GOTO 60
                        ENDIF
                        T(jch+1,jch+1) = CZERO
                     ENDDO
                  ELSE
!
!                 Only test 2 passed -- chase the zero to T(ILAST,ILAST)
!                 Then process as in the case T(ILAST,ILAST)=0
!
                     DO jch = j , ilast - 1
                        ctemp = T(jch,jch+1)
                        CALL ZLARTG(ctemp,T(jch+1,jch+1),c,s,           &
     &                              T(jch,jch+1))
                        T(jch+1,jch+1) = CZERO
                        IF ( jch<ilastm-1 )                             &
     &                       CALL ZROT(ilastm-jch-1,T(jch,jch+2),Ldt,   &
     &                       T(jch+1,jch+2),Ldt,c,s)
                        CALL ZROT(ilastm-jch+2,H(jch,jch-1),Ldh,        &
     &                            H(jch+1,jch-1),Ldh,c,s)
                        IF ( ilq ) CALL ZROT(N,Q(1,jch),1,Q(1,jch+1),1, &
     &                       c,DCONJG(s))
                        ctemp = H(jch+1,jch)
                        CALL ZLARTG(ctemp,H(jch+1,jch-1),c,s,           &
     &                              H(jch+1,jch))
                        H(jch+1,jch-1) = CZERO
                        CALL ZROT(jch+1-ifrstm,H(ifrstm,jch),1,         &
     &                            H(ifrstm,jch-1),1,c,s)
                        CALL ZROT(jch-ifrstm,T(ifrstm,jch),1,           &
     &                            T(ifrstm,jch-1),1,c,s)
                        IF ( ilz ) CALL ZROT(N,Z(1,jch),1,Z(1,jch-1),1, &
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
            Info = 2*N + 1
            GOTO 200
!
!        T(ILAST,ILAST)=0 -- clear H(ILAST,ILAST-1) to split off a
!        1x1 block.
!
 20         ctemp = H(ilast,ilast)
            CALL ZLARTG(ctemp,H(ilast,ilast-1),c,s,H(ilast,ilast))
            H(ilast,ilast-1) = CZERO
            CALL ZROT(ilast-ifrstm,H(ifrstm,ilast),1,H(ifrstm,ilast-1), &
     &                1,c,s)
            CALL ZROT(ilast-ifrstm,T(ifrstm,ilast),1,T(ifrstm,ilast-1), &
     &                1,c,s)
            IF ( ilz ) CALL ZROT(N,Z(1,ilast),1,Z(1,ilast-1),1,c,s)
!
!        H(ILAST,ILAST-1)=0 -- Standardize B, set ALPHA and BETA
!
 40         absb = ABS(T(ilast,ilast))
            IF ( absb>safmin ) THEN
               signbc = DCONJG(T(ilast,ilast)/absb)
               T(ilast,ilast) = absb
               IF ( ilschr ) THEN
                  CALL ZSCAL(ilast-ifrstm,signbc,T(ifrstm,ilast),1)
                  CALL ZSCAL(ilast+1-ifrstm,signbc,H(ifrstm,ilast),1)
               ELSE
                  CALL ZSCAL(1,signbc,H(ilast,ilast),1)
               ENDIF
               IF ( ilz ) CALL ZSCAL(N,signbc,Z(1,ilast),1)
            ELSE
               T(ilast,ilast) = CZERO
            ENDIF
            Alpha(ilast) = H(ilast,ilast)
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
            eshift = CZERO
            IF ( .NOT.ilschr ) THEN
               ilastm = ilast
               IF ( ifrstm>ilast ) ifrstm = Ilo
            ENDIF
            CYCLE
!
!        QZ step
!
!        This iteration only involves rows/columns IFIRST:ILAST.  We
!        assume IFIRST < ILAST, and that the diagonal of B is non-zero.
!
 60         iiter = iiter + 1
            IF ( .NOT.ilschr ) ifrstm = ifirst
!
!        Compute the Shift.
!
!        At this point, IFIRST < ILAST, and the diagonal elements of
!        T(IFIRST:ILAST,IFIRST,ILAST) are larger than BTOL (in
!        magnitude)
!
            IF ( (iiter/10)*10/=iiter ) THEN
!
!           The Wilkinson shift (AEP p.512), i.e., the eigenvalue of
!           the bottom-right 2x2 block of A inv(B) which is nearest to
!           the bottom-right element.
!
!           We factor B as U*D, where U has unit diagonals, and
!           compute (A*inv(D))*inv(U).
!
               u12 = (bscale*T(ilast-1,ilast))/(bscale*T(ilast,ilast))
               ad11 = (ascale*H(ilast-1,ilast-1))                       &
     &                /(bscale*T(ilast-1,ilast-1))
               ad21 = (ascale*H(ilast,ilast-1))                         &
     &                /(bscale*T(ilast-1,ilast-1))
               ad12 = (ascale*H(ilast-1,ilast))/(bscale*T(ilast,ilast))
               ad22 = (ascale*H(ilast,ilast))/(bscale*T(ilast,ilast))
               abi22 = ad22 - u12*ad21
               abi12 = ad12 - u12*ad11
!
               shift = abi22
               ctemp = SQRT(abi12)*SQRT(ad21)
               temp = ABS1(ctemp)
               IF ( ctemp/=ZERO ) THEN
                  x = HALF*(ad11-shift)
                  temp2 = ABS1(x)
                  temp = MAX(temp,ABS1(x))
                  y = temp*SQRT((x/temp)**2+(ctemp/temp)**2)
                  IF ( temp2>ZERO ) THEN
                     IF ( DBLE(x/temp2)*DBLE(y)+DIMAG(x/temp2)*DIMAG(y) &
     &                    <ZERO ) y = -y
                  ENDIF
                  shift = shift - ctemp*ZLADIV(ctemp,(x+y))
               ENDIF
            ELSE
!
!           Exceptional shift.  Chosen for no particularly good reason.
!
               IF ( (iiter/20)*20==iiter .AND.                          &
     &              bscale*ABS1(T(ilast,ilast))>safmin ) THEN
                  eshift = eshift + (ascale*H(ilast,ilast))             &
     &                     /(bscale*T(ilast,ilast))
               ELSE
                  eshift = eshift + (ascale*H(ilast,ilast-1))           &
     &                     /(bscale*T(ilast-1,ilast-1))
               ENDIF
               shift = eshift
            ENDIF
!
!        Now check for two consecutive small subdiagonals.
!
            DO j = ilast - 1 , ifirst + 1 , -1
               istart = j
               ctemp = ascale*H(j,j) - shift*(bscale*T(j,j))
               temp = ABS1(ctemp)
               temp2 = ascale*ABS1(H(j+1,j))
               tempr = MAX(temp,temp2)
               IF ( tempr<ONE .AND. tempr/=ZERO ) THEN
                  temp = temp/tempr
                  temp2 = temp2/tempr
               ENDIF
               IF ( ABS1(H(j,j-1))*temp2<=temp*atol ) GOTO 80
            ENDDO
!
            istart = ifirst
            ctemp = ascale*H(ifirst,ifirst)                             &
     &              - shift*(bscale*T(ifirst,ifirst))
!
!        Do an implicit-shift QZ sweep.
!
!        Initial Q
!
 80         ctemp2 = ascale*H(istart+1,istart)
            CALL ZLARTG(ctemp,ctemp2,c,s,ctemp3)
!
!        Sweep
!
            DO j = istart , ilast - 1
               IF ( j>istart ) THEN
                  ctemp = H(j,j-1)
                  CALL ZLARTG(ctemp,H(j+1,j-1),c,s,H(j,j-1))
                  H(j+1,j-1) = CZERO
               ENDIF
!
               DO jc = j , ilastm
                  ctemp = c*H(j,jc) + s*H(j+1,jc)
                  H(j+1,jc) = -DCONJG(s)*H(j,jc) + c*H(j+1,jc)
                  H(j,jc) = ctemp
                  ctemp2 = c*T(j,jc) + s*T(j+1,jc)
                  T(j+1,jc) = -DCONJG(s)*T(j,jc) + c*T(j+1,jc)
                  T(j,jc) = ctemp2
               ENDDO
               IF ( ilq ) THEN
                  DO jr = 1 , N
                     ctemp = c*Q(jr,j) + DCONJG(s)*Q(jr,j+1)
                     Q(jr,j+1) = -s*Q(jr,j) + c*Q(jr,j+1)
                     Q(jr,j) = ctemp
                  ENDDO
               ENDIF
!
               ctemp = T(j+1,j+1)
               CALL ZLARTG(ctemp,T(j+1,j),c,s,T(j+1,j+1))
               T(j+1,j) = CZERO
!
               DO jr = ifrstm , MIN(j+2,ilast)
                  ctemp = c*H(jr,j+1) + s*H(jr,j)
                  H(jr,j) = -DCONJG(s)*H(jr,j+1) + c*H(jr,j)
                  H(jr,j+1) = ctemp
               ENDDO
               DO jr = ifrstm , j
                  ctemp = c*T(jr,j+1) + s*T(jr,j)
                  T(jr,j) = -DCONJG(s)*T(jr,j+1) + c*T(jr,j)
                  T(jr,j+1) = ctemp
               ENDDO
               IF ( ilz ) THEN
                  DO jr = 1 , N
                     ctemp = c*Z(jr,j+1) + s*Z(jr,j)
                     Z(jr,j) = -DCONJG(s)*Z(jr,j+1) + c*Z(jr,j)
                     Z(jr,j+1) = ctemp
                  ENDDO
               ENDIF
            ENDDO
!
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
         absb = ABS(T(j,j))
         IF ( absb>safmin ) THEN
            signbc = DCONJG(T(j,j)/absb)
            T(j,j) = absb
            IF ( ilschr ) THEN
               CALL ZSCAL(j-1,signbc,T(1,j),1)
               CALL ZSCAL(j,signbc,H(1,j),1)
            ELSE
               CALL ZSCAL(1,signbc,H(j,j),1)
            ENDIF
            IF ( ilz ) CALL ZSCAL(N,signbc,Z(1,j),1)
         ELSE
            T(j,j) = CZERO
         ENDIF
         Alpha(j) = H(j,j)
         Beta(j) = T(j,j)
      ENDDO
!
!     Normal Termination
!
      Info = 0
!
!     Exit (other than argument error) -- return optimal workspace size
!
 200  Work(1) = DCMPLX(N)
!
!     End of ZHGEQZ
!
      END SUBROUTINE ZHGEQZ
