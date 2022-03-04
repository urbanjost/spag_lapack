!*==ztrsna.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b ZTRSNA
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZTRSNA + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ztrsna.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ztrsna.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ztrsna.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZTRSNA( JOB, HOWMNY, SELECT, N, T, LDT, VL, LDVL, VR,
!                          LDVR, S, SEP, MM, M, WORK, LDWORK, RWORK,
!                          INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          HOWMNY, JOB
!       INTEGER            INFO, LDT, LDVL, LDVR, LDWORK, M, MM, N
!       ..
!       .. Array Arguments ..
!       LOGICAL            SELECT( * )
!       DOUBLE PRECISION   RWORK( * ), S( * ), SEP( * )
!       COMPLEX*16         T( LDT, * ), VL( LDVL, * ), VR( LDVR, * ),
!      $                   WORK( LDWORK, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZTRSNA estimates reciprocal condition numbers for specified
!> eigenvalues and/or right eigenvectors of a complex upper triangular
!> matrix T (or of any matrix Q*T*Q**H with Q unitary).
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
!>          for the j-th eigenpair, SELECT(j) must be set to .TRUE..
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
!>          T is COMPLEX*16 array, dimension (LDT,N)
!>          The upper triangular matrix T.
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
!>          VL is COMPLEX*16 array, dimension (LDVL,M)
!>          If JOB = 'E' or 'B', VL must contain left eigenvectors of T
!>          (or of any Q*T*Q**H with Q unitary), corresponding to the
!>          eigenpairs specified by HOWMNY and SELECT. The eigenvectors
!>          must be stored in consecutive columns of VL, as returned by
!>          ZHSEIN or ZTREVC.
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
!>          VR is COMPLEX*16 array, dimension (LDVR,M)
!>          If JOB = 'E' or 'B', VR must contain right eigenvectors of T
!>          (or of any Q*T*Q**H with Q unitary), corresponding to the
!>          eigenpairs specified by HOWMNY and SELECT. The eigenvectors
!>          must be stored in consecutive columns of VR, as returned by
!>          ZHSEIN or ZTREVC.
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
!>          array. Thus S(j), SEP(j), and the j-th columns of VL and VR
!>          all correspond to the same eigenpair (but not in general the
!>          j-th eigenpair, unless all eigenpairs are selected).
!>          If JOB = 'V', S is not referenced.
!> \endverbatim
!>
!> \param[out] SEP
!> \verbatim
!>          SEP is DOUBLE PRECISION array, dimension (MM)
!>          If JOB = 'V' or 'B', the estimated reciprocal condition
!>          numbers of the selected eigenvectors, stored in consecutive
!>          elements of the array.
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
!>          WORK is COMPLEX*16 array, dimension (LDWORK,N+6)
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
!> \param[out] RWORK
!> \verbatim
!>          RWORK is DOUBLE PRECISION array, dimension (N)
!>          If JOB = 'E', RWORK is not referenced.
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
!> \ingroup complex16OTHERcomputational
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  The reciprocal of the condition number of an eigenvalue lambda is
!>  defined as
!>
!>          S(lambda) = |v**H*u| / (norm(u)*norm(v))
!>
!>  where u and v are the right and left eigenvectors of T corresponding
!>  to lambda; v**H denotes the conjugate transpose of v, and norm(u)
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
      SUBROUTINE ZTRSNA(Job,Howmny,Select,N,T,Ldt,Vl,Ldvl,Vr,Ldvr,S,Sep,&
     &                  Mm,M,Work,Ldwork,Rwork,Info)
      USE F77KINDS                        
      USE S_DLABAD
      USE S_DLAMCH
      USE S_DZNRM2
      USE S_IZAMAX
      USE S_LSAME
      USE S_XERBLA
      USE S_ZDOTC
      USE S_ZDRSCL
      USE S_ZLACN2
      USE S_ZLACPY
      USE S_ZLATRS
      USE S_ZTREXC
      IMPLICIT NONE
!*--ZTRSNA265
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D0 + 0
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Job
      CHARACTER :: Howmny
      LOGICAL , INTENT(IN) , DIMENSION(*) :: Select
      INTEGER :: N
      COMPLEX(CX16KIND) , DIMENSION(Ldt,*) :: T
      INTEGER :: Ldt
      COMPLEX(CX16KIND) , DIMENSION(Ldvl,*) :: Vl
      INTEGER , INTENT(IN) :: Ldvl
      COMPLEX(CX16KIND) , DIMENSION(Ldvr,*) :: Vr
      INTEGER , INTENT(IN) :: Ldvr
      REAL(R8KIND) , INTENT(OUT) , DIMENSION(*) :: S
      REAL(R8KIND) , INTENT(OUT) , DIMENSION(*) :: Sep
      INTEGER , INTENT(IN) :: Mm
      INTEGER , INTENT(INOUT) :: M
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(Ldwork,*) :: Work
      INTEGER :: Ldwork
      REAL(R8KIND) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      REAL(R8KIND) :: bignum , eps , est , lnrm , rnrm , scale ,        &
     &                smlnum , xnorm
      REAL(R8KIND) :: CABS1
      COMPLEX(CX16KIND) :: cdum , prod
      COMPLEX(CX16KIND) , DIMENSION(1) :: dummy
      INTEGER :: i , ierr , ix , j , k , kase , ks
      INTEGER , DIMENSION(3) :: isave
      CHARACTER :: normin
      LOGICAL :: somcon , wantbh , wants , wantsp
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
!     .. Statement Functions ..
!     ..
!     .. Statement Function definitions ..
      CABS1(cdum) = ABS(DBLE(cdum)) + ABS(DIMAG(cdum))
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
!     Set M to the number of eigenpairs for which condition numbers are
!     to be computed.
!
      IF ( somcon ) THEN
         M = 0
         DO j = 1 , N
            IF ( Select(j) ) M = M + 1
         ENDDO
      ELSE
         M = N
      ENDIF
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
      ELSEIF ( Mm<M ) THEN
         Info = -13
      ELSEIF ( Ldwork<1 .OR. (wantsp .AND. Ldwork<N) ) THEN
         Info = -16
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('ZTRSNA',-Info)
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
      ks = 1
      DO k = 1 , N
!
         IF ( somcon ) THEN
            IF ( .NOT.Select(k) ) CYCLE
         ENDIF
!
         IF ( wants ) THEN
!
!           Compute the reciprocal condition number of the k-th
!           eigenvalue.
!
            prod = ZDOTC(N,Vr(1,ks),1,Vl(1,ks),1)
            rnrm = DZNRM2(N,Vr(1,ks),1)
            lnrm = DZNRM2(N,Vl(1,ks),1)
            S(ks) = ABS(prod)/(rnrm*lnrm)
!
         ENDIF
!
         IF ( wantsp ) THEN
!
!           Estimate the reciprocal condition number of the k-th
!           eigenvector.
!
!           Copy the matrix T to the array WORK and swap the k-th
!           diagonal element to the (1,1) position.
!
            CALL ZLACPY('Full',N,N,T,Ldt,Work,Ldwork)
            CALL ZTREXC('No Q',N,Work,Ldwork,dummy,1,k,1,ierr)
!
!           Form  C = T22 - lambda*I in WORK(2:N,2:N).
!
            DO i = 2 , N
               Work(i,i) = Work(i,i) - Work(1,1)
            ENDDO
!
!           Estimate a lower bound for the 1-norm of inv(C**H). The 1st
!           and (N+1)th columns of WORK are used to store work vectors.
!
            Sep(ks) = ZERO
            est = ZERO
            kase = 0
            normin = 'N'
            DO
               CALL ZLACN2(N-1,Work(1,N+1),Work,est,kase,isave)
!
               IF ( kase/=0 ) THEN
                  IF ( kase==1 ) THEN
!
!                 Solve C**H*x = scale*b
!
                     CALL ZLATRS('Upper','Conjugate transpose',         &
     &                           'Nonunit',normin,N-1,Work(2,2),Ldwork, &
     &                           Work,scale,Rwork,ierr)
                  ELSE
!
!                 Solve C*x = scale*b
!
                     CALL ZLATRS('Upper','No transpose','Nonunit',      &
     &                           normin,N-1,Work(2,2),Ldwork,Work,scale,&
     &                           Rwork,ierr)
                  ENDIF
                  normin = 'Y'
                  IF ( scale/=ONE ) THEN
!
!                 Multiply by 1/SCALE if doing so will not cause
!                 overflow.
!
                     ix = IZAMAX(N-1,Work,1)
                     xnorm = CABS1(Work(ix,1))
                     IF ( scale<xnorm*smlnum .OR. scale==ZERO ) EXIT
                     CALL ZDRSCL(N,scale,Work,1)
                  ENDIF
                  CYCLE
               ENDIF
!
               Sep(ks) = ONE/MAX(est,smlnum)
               EXIT
            ENDDO
         ENDIF
!
         ks = ks + 1
      ENDDO
!
!     End of ZTRSNA
!
      END SUBROUTINE ZTRSNA
