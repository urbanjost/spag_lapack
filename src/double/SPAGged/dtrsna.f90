!*==dtrsna.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b DTRSNA
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DTRSNA + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dtrsna.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dtrsna.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dtrsna.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DTRSNA( JOB, HOWMNY, SELECT, N, T, LDT, VL, LDVL, VR,
!                          LDVR, S, SEP, MM, M, WORK, LDWORK, IWORK,
!                          INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          HOWMNY, JOB
!       INTEGER            INFO, LDT, LDVL, LDVR, LDWORK, M, MM, N
!       ..
!       .. Array Arguments ..
!       LOGICAL            SELECT( * )
!       INTEGER            IWORK( * )
!       DOUBLE PRECISION   S( * ), SEP( * ), T( LDT, * ), VL( LDVL, * ),
!      $                   VR( LDVR, * ), WORK( LDWORK, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DTRSNA estimates reciprocal condition numbers for specified
!> eigenvalues and/or right eigenvectors of a real upper
!> quasi-triangular matrix T (or of any matrix Q*T*Q**T with Q
!> orthogonal).
!>
!> T must be in Schur canonical form (as returned by DHSEQR), that is,
!> block upper triangular with 1-by-1 and 2-by-2 diagonal blocks; each
!> 2-by-2 diagonal block has its diagonal elements equal and its
!> off-diagonal elements of opposite sign.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] JOB
!> \verbatim
!>          JOB is CHARACTER*1
!>          Specifies whether condition numbers are required for
!>          eigenvalues (S) or eigenvectors (SEP):
!>          = 'E': for eigenvalues only (S);
!>          = 'V': for eigenvectors only (SEP);
!>          = 'B': for both eigenvalues and eigenvectors (S and SEP).
!> \endverbatim
!>
!> \param[in] HOWMNY
!> \verbatim
!>          HOWMNY is CHARACTER*1
!>          = 'A': compute condition numbers for all eigenpairs;
!>          = 'S': compute condition numbers for selected eigenpairs
!>                 specified by the array SELECT.
!> \endverbatim
!>
!> \param[in] SELECT
!> \verbatim
!>          SELECT is LOGICAL array, dimension (N)
!>          If HOWMNY = 'S', SELECT specifies the eigenpairs for which
!>          condition numbers are required. To select condition numbers
!>          for the eigenpair corresponding to a real eigenvalue w(j),
!>          SELECT(j) must be set to .TRUE.. To select condition numbers
!>          corresponding to a complex conjugate pair of eigenvalues w(j)
!>          and w(j+1), either SELECT(j) or SELECT(j+1) or both, must be
!>          set to .TRUE..
!>          If HOWMNY = 'A', SELECT is not referenced.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix T. N >= 0.
!> \endverbatim
!>
!> \param[in] T
!> \verbatim
!>          T is DOUBLE PRECISION array, dimension (LDT,N)
!>          The upper quasi-triangular matrix T, in Schur canonical form.
!> \endverbatim
!>
!> \param[in] LDT
!> \verbatim
!>          LDT is INTEGER
!>          The leading dimension of the array T. LDT >= max(1,N).
!> \endverbatim
!>
!> \param[in] VL
!> \verbatim
!>          VL is DOUBLE PRECISION array, dimension (LDVL,M)
!>          If JOB = 'E' or 'B', VL must contain left eigenvectors of T
!>          (or of any Q*T*Q**T with Q orthogonal), corresponding to the
!>          eigenpairs specified by HOWMNY and SELECT. The eigenvectors
!>          must be stored in consecutive columns of VL, as returned by
!>          DHSEIN or DTREVC.
!>          If JOB = 'V', VL is not referenced.
!> \endverbatim
!>
!> \param[in] LDVL
!> \verbatim
!>          LDVL is INTEGER
!>          The leading dimension of the array VL.
!>          LDVL >= 1; and if JOB = 'E' or 'B', LDVL >= N.
!> \endverbatim
!>
!> \param[in] VR
!> \verbatim
!>          VR is DOUBLE PRECISION array, dimension (LDVR,M)
!>          If JOB = 'E' or 'B', VR must contain right eigenvectors of T
!>          (or of any Q*T*Q**T with Q orthogonal), corresponding to the
!>          eigenpairs specified by HOWMNY and SELECT. The eigenvectors
!>          must be stored in consecutive columns of VR, as returned by
!>          DHSEIN or DTREVC.
!>          If JOB = 'V', VR is not referenced.
!> \endverbatim
!>
!> \param[in] LDVR
!> \verbatim
!>          LDVR is INTEGER
!>          The leading dimension of the array VR.
!>          LDVR >= 1; and if JOB = 'E' or 'B', LDVR >= N.
!> \endverbatim
!>
!> \param[out] S
!> \verbatim
!>          S is DOUBLE PRECISION array, dimension (MM)
!>          If JOB = 'E' or 'B', the reciprocal condition numbers of the
!>          selected eigenvalues, stored in consecutive elements of the
!>          array. For a complex conjugate pair of eigenvalues two
!>          consecutive elements of S are set to the same value. Thus
!>          S(j), SEP(j), and the j-th columns of VL and VR all
!>          correspond to the same eigenpair (but not in general the
!>          j-th eigenpair, unless all eigenpairs are selected).
!>          If JOB = 'V', S is not referenced.
!> \endverbatim
!>
!> \param[out] SEP
!> \verbatim
!>          SEP is DOUBLE PRECISION array, dimension (MM)
!>          If JOB = 'V' or 'B', the estimated reciprocal condition
!>          numbers of the selected eigenvectors, stored in consecutive
!>          elements of the array. For a complex eigenvector two
!>          consecutive elements of SEP are set to the same value. If
!>          the eigenvalues cannot be reordered to compute SEP(j), SEP(j)
!>          is set to 0; this can only occur when the true value would be
!>          very small anyway.
!>          If JOB = 'E', SEP is not referenced.
!> \endverbatim
!>
!> \param[in] MM
!> \verbatim
!>          MM is INTEGER
!>          The number of elements in the arrays S (if JOB = 'E' or 'B')
!>           and/or SEP (if JOB = 'V' or 'B'). MM >= M.
!> \endverbatim
!>
!> \param[out] M
!> \verbatim
!>          M is INTEGER
!>          The number of elements of the arrays S and/or SEP actually
!>          used to store the estimated condition numbers.
!>          If HOWMNY = 'A', M is set to N.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (LDWORK,N+6)
!>          If JOB = 'E', WORK is not referenced.
!> \endverbatim
!>
!> \param[in] LDWORK
!> \verbatim
!>          LDWORK is INTEGER
!>          The leading dimension of the array WORK.
!>          LDWORK >= 1; and if JOB = 'V' or 'B', LDWORK >= N.
!> \endverbatim
!>
!> \param[out] IWORK
!> \verbatim
!>          IWORK is INTEGER array, dimension (2*(N-1))
!>          If JOB = 'E', IWORK is not referenced.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0: successful exit
!>          < 0: if INFO = -i, the i-th argument had an illegal value
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
!> \date November 2017
!
!> \ingroup doubleOTHERcomputational
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  The reciprocal of the condition number of an eigenvalue lambda is
!>  defined as
!>
!>          S(lambda) = |v**T*u| / (norm(u)*norm(v))
!>
!>  where u and v are the right and left eigenvectors of T corresponding
!>  to lambda; v**T denotes the transpose of v, and norm(u)
!>  denotes the Euclidean norm. These reciprocal condition numbers always
!>  lie between zero (very badly conditioned) and one (very well
!>  conditioned). If n = 1, S(lambda) is defined to be 1.
!>
!>  An approximate error bound for a computed eigenvalue W(i) is given by
!>
!>                      EPS * norm(T) / S(i)
!>
!>  where EPS is the machine precision.
!>
!>  The reciprocal of the condition number of the right eigenvector u
!>  corresponding to lambda is defined as follows. Suppose
!>
!>              T = ( lambda  c  )
!>                  (   0    T22 )
!>
!>  Then the reciprocal condition number is
!>
!>          SEP( lambda, T22 ) = sigma-min( T22 - lambda*I )
!>
!>  where sigma-min denotes the smallest singular value. We approximate
!>  the smallest singular value by the reciprocal of an estimate of the
!>  one-norm of the inverse of T22 - lambda*I. If n = 1, SEP(1) is
!>  defined to be abs(T(1,1)).
!>
!>  An approximate error bound for a computed right eigenvector VR(i)
!>  is given by
!>
!>                      EPS * norm(T) / SEP(i)
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE DTRSNA(Job,Howmny,Select,N,T,Ldt,Vl,Ldvl,Vr,Ldvr,S,Sep,&
     &                  Mm,M,Work,Ldwork,Iwork,Info)
      USE F77KINDS                        
      USE S_DDOT
      USE S_DLABAD
      USE S_DLACN2
      USE S_DLACPY
      USE S_DLAMCH
      USE S_DLAPY2
      USE S_DLAQTR
      USE S_DNRM2
      USE S_DTREXC
      USE S_LSAME
      USE S_XERBLA
      IMPLICIT NONE
!*--DTRSNA280
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0 ,      &
     &                              TWO = 2.0D+0
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Job
      CHARACTER :: Howmny
      LOGICAL , INTENT(IN) , DIMENSION(*) :: Select
      INTEGER :: N
      REAL(R8KIND) , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      REAL(R8KIND) , DIMENSION(Ldvl,*) :: Vl
      INTEGER , INTENT(IN) :: Ldvl
      REAL(R8KIND) , DIMENSION(Ldvr,*) :: Vr
      INTEGER , INTENT(IN) :: Ldvr
      REAL(R8KIND) , INTENT(OUT) , DIMENSION(*) :: S
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Sep
      INTEGER , INTENT(IN) :: Mm
      INTEGER , INTENT(INOUT) :: M
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(Ldwork,*) :: Work
      INTEGER :: Ldwork
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      REAL(R8KIND) :: bignum , cond , cs , delta , dumm , eps , est ,   &
     &                lnrm , mu , prod , prod1 , prod2 , rnrm , scale , &
     &                smlnum , sn
      REAL(R8KIND) , DIMENSION(1) :: dummy
      INTEGER :: i , ierr , ifst , ilst , j , k , kase , ks , n2 , nn
      INTEGER , DIMENSION(3) :: isave
      LOGICAL :: pair , somcon , wantbh , wants , wantsp
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
!     .. External Functions ..
!     ..
!     .. External Subroutines ..
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. Executable Statements ..
!
!     Decode and test the input parameters
!
      wantbh = LSAME(Job,'B')
      wants = LSAME(Job,'E') .OR. wantbh
      wantsp = LSAME(Job,'V') .OR. wantbh
!
      somcon = LSAME(Howmny,'S')
!
      Info = 0
      IF ( .NOT.wants .AND. .NOT.wantsp ) THEN
         Info = -1
      ELSEIF ( .NOT.LSAME(Howmny,'A') .AND. .NOT.somcon ) THEN
         Info = -2
      ELSEIF ( N<0 ) THEN
         Info = -4
      ELSEIF ( Ldt<MAX(1,N) ) THEN
         Info = -6
      ELSEIF ( Ldvl<1 .OR. (wants .AND. Ldvl<N) ) THEN
         Info = -8
      ELSEIF ( Ldvr<1 .OR. (wants .AND. Ldvr<N) ) THEN
         Info = -10
      ELSE
!
!        Set M to the number of eigenpairs for which condition numbers
!        are required, and test MM.
!
         IF ( somcon ) THEN
            M = 0
            pair = .FALSE.
            DO k = 1 , N
               IF ( pair ) THEN
                  pair = .FALSE.
               ELSEIF ( k<N ) THEN
                  IF ( T(k+1,k)==ZERO ) THEN
                     IF ( Select(k) ) M = M + 1
                  ELSE
                     pair = .TRUE.
                     IF ( Select(k) .OR. Select(k+1) ) M = M + 2
                  ENDIF
               ELSE
                  IF ( Select(N) ) M = M + 1
               ENDIF
            ENDDO
         ELSE
            M = N
         ENDIF
!
         IF ( Mm<M ) THEN
            Info = -13
         ELSEIF ( Ldwork<1 .OR. (wantsp .AND. Ldwork<N) ) THEN
            Info = -16
         ENDIF
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('DTRSNA',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( N==0 ) RETURN
!
      IF ( N==1 ) THEN
         IF ( somcon ) THEN
            IF ( .NOT.Select(1) ) RETURN
         ENDIF
         IF ( wants ) S(1) = ONE
         IF ( wantsp ) Sep(1) = ABS(T(1,1))
         RETURN
      ENDIF
!
!     Get machine constants
!
      eps = DLAMCH('P')
      smlnum = DLAMCH('S')/eps
      bignum = ONE/smlnum
      CALL DLABAD(smlnum,bignum)
!
      ks = 0
      pair = .FALSE.
      DO k = 1 , N
!
!        Determine whether T(k,k) begins a 1-by-1 or 2-by-2 block.
!
         IF ( pair ) THEN
            pair = .FALSE.
            CYCLE
         ELSE
            IF ( k<N ) pair = T(k+1,k)/=ZERO
         ENDIF
!
!        Determine whether condition numbers are required for the k-th
!        eigenpair.
!
         IF ( somcon ) THEN
            IF ( pair ) THEN
               IF ( .NOT.Select(k) .AND. .NOT.Select(k+1) ) CYCLE
            ELSEIF ( .NOT.Select(k) ) THEN
               CYCLE
            ENDIF
         ENDIF
!
         ks = ks + 1
!
         IF ( wants ) THEN
!
!           Compute the reciprocal condition number of the k-th
!           eigenvalue.
!
            IF ( .NOT.pair ) THEN
!
!              Real eigenvalue.
!
               prod = DDOT(N,Vr(1,ks),1,Vl(1,ks),1)
               rnrm = DNRM2(N,Vr(1,ks),1)
               lnrm = DNRM2(N,Vl(1,ks),1)
               S(ks) = ABS(prod)/(rnrm*lnrm)
            ELSE
!
!              Complex eigenvalue.
!
               prod1 = DDOT(N,Vr(1,ks),1,Vl(1,ks),1)
               prod1 = prod1 + DDOT(N,Vr(1,ks+1),1,Vl(1,ks+1),1)
               prod2 = DDOT(N,Vl(1,ks),1,Vr(1,ks+1),1)
               prod2 = prod2 - DDOT(N,Vl(1,ks+1),1,Vr(1,ks),1)
               rnrm = DLAPY2(DNRM2(N,Vr(1,ks),1),DNRM2(N,Vr(1,ks+1),1))
               lnrm = DLAPY2(DNRM2(N,Vl(1,ks),1),DNRM2(N,Vl(1,ks+1),1))
               cond = DLAPY2(prod1,prod2)/(rnrm*lnrm)
               S(ks) = cond
               S(ks+1) = cond
            ENDIF
         ENDIF
!
         IF ( wantsp ) THEN
!
!           Estimate the reciprocal condition number of the k-th
!           eigenvector.
!
!           Copy the matrix T to the array WORK and swap the diagonal
!           block beginning at T(k,k) to the (1,1) position.
!
            CALL DLACPY('Full',N,N,T,Ldt,Work,Ldwork)
            ifst = k
            ilst = 1
            CALL DTREXC('No Q',N,Work,Ldwork,dummy,1,ifst,ilst,         &
     &                  Work(1,N+1),ierr)
!
            IF ( ierr==1 .OR. ierr==2 ) THEN
!
!              Could not swap because blocks not well separated
!
               scale = ONE
               est = bignum
            ELSE
!
!              Reordering successful
!
               IF ( Work(2,1)==ZERO ) THEN
!
!                 Form C = T22 - lambda*I in WORK(2:N,2:N).
!
                  DO i = 2 , N
                     Work(i,i) = Work(i,i) - Work(1,1)
                  ENDDO
                  n2 = 1
                  nn = N - 1
               ELSE
!
!                 Triangularize the 2 by 2 block by unitary
!                 transformation U = [  cs   i*ss ]
!                                    [ i*ss   cs  ].
!                 such that the (1,1) position of WORK is complex
!                 eigenvalue lambda with positive imaginary part. (2,2)
!                 position of WORK is the complex eigenvalue lambda
!                 with negative imaginary  part.
!
                  mu = SQRT(ABS(Work(1,2)))*SQRT(ABS(Work(2,1)))
                  delta = DLAPY2(mu,Work(2,1))
                  cs = mu/delta
                  sn = -Work(2,1)/delta
!
!                 Form
!
!                 C**T = WORK(2:N,2:N) + i*[rwork(1) ..... rwork(n-1) ]
!                                          [   mu                     ]
!                                          [         ..               ]
!                                          [             ..           ]
!                                          [                  mu      ]
!                 where C**T is transpose of matrix C,
!                 and RWORK is stored starting in the N+1-st column of
!                 WORK.
!
                  DO j = 3 , N
                     Work(2,j) = cs*Work(2,j)
                     Work(j,j) = Work(j,j) - Work(1,1)
                  ENDDO
                  Work(2,2) = ZERO
!
                  Work(1,N+1) = TWO*mu
                  DO i = 2 , N - 1
                     Work(i,N+1) = sn*Work(1,i+1)
                  ENDDO
                  n2 = 2
                  nn = 2*(N-1)
               ENDIF
!
!              Estimate norm(inv(C**T))
!
               est = ZERO
               kase = 0
               DO
                  CALL DLACN2(nn,Work(1,N+2),Work(1,N+4),Iwork,est,kase,&
     &                        isave)
                  IF ( kase/=0 ) THEN
                     IF ( kase==1 ) THEN
                        IF ( n2==1 ) THEN
!
!                       Real eigenvalue: solve C**T*x = scale*c.
!
                           CALL DLAQTR(.TRUE.,.TRUE.,N-1,Work(2,2),     &
     &                                 Ldwork,dummy,dumm,scale,         &
     &                                 Work(1,N+4),Work(1,N+6),ierr)
                        ELSE
!
!                       Complex eigenvalue: solve
!                       C**T*(p+iq) = scale*(c+id) in real arithmetic.
!
                           CALL DLAQTR(.TRUE.,.FALSE.,N-1,Work(2,2),    &
     &                                 Ldwork,Work(1,N+1),mu,scale,     &
     &                                 Work(1,N+4),Work(1,N+6),ierr)
                        ENDIF
                     ELSEIF ( n2==1 ) THEN
!
!                       Real eigenvalue: solve C*x = scale*c.
!
                        CALL DLAQTR(.FALSE.,.TRUE.,N-1,Work(2,2),Ldwork,&
     &                              dummy,dumm,scale,Work(1,N+4),       &
     &                              Work(1,N+6),ierr)
                     ELSE
!
!                       Complex eigenvalue: solve
!                       C*(p+iq) = scale*(c+id) in real arithmetic.
!
                        CALL DLAQTR(.FALSE.,.FALSE.,N-1,Work(2,2),      &
     &                              Ldwork,Work(1,N+1),mu,scale,        &
     &                              Work(1,N+4),Work(1,N+6),ierr)
!
                     ENDIF
!
                     CYCLE
                  ENDIF
                  EXIT
               ENDDO
            ENDIF
!
            Sep(ks) = scale/MAX(est,smlnum)
            IF ( pair ) Sep(ks+1) = Sep(ks)
         ENDIF
!
         IF ( pair ) ks = ks + 1
!
      ENDDO
!
!     End of DTRSNA
!
      END SUBROUTINE DTRSNA
