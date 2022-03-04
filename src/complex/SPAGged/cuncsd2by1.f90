!*==cuncsd2by1.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
 
!> \brief \b CUNCSD2BY1
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CUNCSD2BY1 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cuncsd2by1.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cuncsd2by1.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cuncsd2by1.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CUNCSD2BY1( JOBU1, JOBU2, JOBV1T, M, P, Q, X11, LDX11,
!                              X21, LDX21, THETA, U1, LDU1, U2, LDU2, V1T,
!                              LDV1T, WORK, LWORK, RWORK, LRWORK, IWORK,
!                              INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          JOBU1, JOBU2, JOBV1T
!       INTEGER            INFO, LDU1, LDU2, LDV1T, LWORK, LDX11, LDX21,
!      $                   M, P, Q
!       INTEGER            LRWORK, LRWORKMIN, LRWORKOPT
!       ..
!       .. Array Arguments ..
!       REAL               RWORK(*)
!       REAL               THETA(*)
!       COMPLEX            U1(LDU1,*), U2(LDU2,*), V1T(LDV1T,*), WORK(*),
!      $                   X11(LDX11,*), X21(LDX21,*)
!       INTEGER            IWORK(*)
!       ..
!
!
!> \par Purpose:
!  =============
!>
!>\verbatim
!>
!> CUNCSD2BY1 computes the CS decomposition of an M-by-Q matrix X with
!> orthonormal columns that has been partitioned into a 2-by-1 block
!> structure:
!>
!>                                [  I1 0  0 ]
!>                                [  0  C  0 ]
!>          [ X11 ]   [ U1 |    ] [  0  0  0 ]
!>      X = [-----] = [---------] [----------] V1**T .
!>          [ X21 ]   [    | U2 ] [  0  0  0 ]
!>                                [  0  S  0 ]
!>                                [  0  0  I2]
!>
!> X11 is P-by-Q. The unitary matrices U1, U2, and V1 are P-by-P,
!> (M-P)-by-(M-P), and Q-by-Q, respectively. C and S are R-by-R
!> nonnegative diagonal matrices satisfying C^2 + S^2 = I, in which
!> R = MIN(P,M-P,Q,M-Q). I1 is a K1-by-K1 identity matrix and I2 is a
!> K2-by-K2 identity matrix, where K1 = MAX(Q+P-M,0), K2 = MAX(Q-P,0).
!>
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] JOBU1
!> \verbatim
!>          JOBU1 is CHARACTER
!>          = 'Y':      U1 is computed;
!>          otherwise:  U1 is not computed.
!> \endverbatim
!>
!> \param[in] JOBU2
!> \verbatim
!>          JOBU2 is CHARACTER
!>          = 'Y':      U2 is computed;
!>          otherwise:  U2 is not computed.
!> \endverbatim
!>
!> \param[in] JOBV1T
!> \verbatim
!>          JOBV1T is CHARACTER
!>          = 'Y':      V1T is computed;
!>          otherwise:  V1T is not computed.
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows in X.
!> \endverbatim
!>
!> \param[in] P
!> \verbatim
!>          P is INTEGER
!>          The number of rows in X11. 0 <= P <= M.
!> \endverbatim
!>
!> \param[in] Q
!> \verbatim
!>          Q is INTEGER
!>          The number of columns in X11 and X21. 0 <= Q <= M.
!> \endverbatim
!>
!> \param[in,out] X11
!> \verbatim
!>          X11 is COMPLEX array, dimension (LDX11,Q)
!>          On entry, part of the unitary matrix whose CSD is desired.
!> \endverbatim
!>
!> \param[in] LDX11
!> \verbatim
!>          LDX11 is INTEGER
!>          The leading dimension of X11. LDX11 >= MAX(1,P).
!> \endverbatim
!>
!> \param[in,out] X21
!> \verbatim
!>          X21 is COMPLEX array, dimension (LDX21,Q)
!>          On entry, part of the unitary matrix whose CSD is desired.
!> \endverbatim
!>
!> \param[in] LDX21
!> \verbatim
!>          LDX21 is INTEGER
!>          The leading dimension of X21. LDX21 >= MAX(1,M-P).
!> \endverbatim
!>
!> \param[out] THETA
!> \verbatim
!>          THETA is REAL array, dimension (R), in which R =
!>          MIN(P,M-P,Q,M-Q).
!>          C = DIAG( COS(THETA(1)), ... , COS(THETA(R)) ) and
!>          S = DIAG( SIN(THETA(1)), ... , SIN(THETA(R)) ).
!> \endverbatim
!>
!> \param[out] U1
!> \verbatim
!>          U1 is COMPLEX array, dimension (P)
!>          If JOBU1 = 'Y', U1 contains the P-by-P unitary matrix U1.
!> \endverbatim
!>
!> \param[in] LDU1
!> \verbatim
!>          LDU1 is INTEGER
!>          The leading dimension of U1. If JOBU1 = 'Y', LDU1 >=
!>          MAX(1,P).
!> \endverbatim
!>
!> \param[out] U2
!> \verbatim
!>          U2 is COMPLEX array, dimension (M-P)
!>          If JOBU2 = 'Y', U2 contains the (M-P)-by-(M-P) unitary
!>          matrix U2.
!> \endverbatim
!>
!> \param[in] LDU2
!> \verbatim
!>          LDU2 is INTEGER
!>          The leading dimension of U2. If JOBU2 = 'Y', LDU2 >=
!>          MAX(1,M-P).
!> \endverbatim
!>
!> \param[out] V1T
!> \verbatim
!>          V1T is COMPLEX array, dimension (Q)
!>          If JOBV1T = 'Y', V1T contains the Q-by-Q matrix unitary
!>          matrix V1**T.
!> \endverbatim
!>
!> \param[in] LDV1T
!> \verbatim
!>          LDV1T is INTEGER
!>          The leading dimension of V1T. If JOBV1T = 'Y', LDV1T >=
!>          MAX(1,Q).
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
!>          The dimension of the array WORK.
!>
!>          If LWORK = -1, then a workspace query is assumed; the routine
!>          only calculates the optimal size of the WORK and RWORK
!>          arrays, returns this value as the first entry of the WORK
!>          and RWORK array, respectively, and no error message related
!>          to LWORK or LRWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is REAL array, dimension (MAX(1,LRWORK))
!>          On exit, if INFO = 0, RWORK(1) returns the optimal LRWORK.
!>          If INFO > 0 on exit, RWORK(2:R) contains the values PHI(1),
!>          ..., PHI(R-1) that, together with THETA(1), ..., THETA(R),
!>          define the matrix in intermediate bidiagonal-block form
!>          remaining after nonconvergence. INFO specifies the number
!>          of nonzero PHI's.
!> \endverbatim
!>
!> \param[in] LRWORK
!> \verbatim
!>          LRWORK is INTEGER
!>          The dimension of the array RWORK.
!>
!>          If LRWORK=-1, then a workspace query is assumed; the routine
!>          only calculates the optimal size of the WORK and RWORK
!>          arrays, returns this value as the first entry of the WORK
!>          and RWORK array, respectively, and no error message related
!>          to LWORK or LRWORK is issued by XERBLA.
!> \endverbatim
!
!> \param[out] IWORK
!> \verbatim
!>          IWORK is INTEGER array, dimension (M-MIN(P,M-P,Q,M-Q))
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit.
!>          < 0:  if INFO = -i, the i-th argument had an illegal value.
!>          > 0:  CBBCSD did not converge. See the description of WORK
!>                above for details.
!> \endverbatim
!
!> \par References:
!  ================
!>
!>  [1] Brian D. Sutton. Computing the complete CS decomposition. Numer.
!>      Algorithms, 50(1):33-65, 2009.
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
!> \ingroup complexOTHERcomputational
!
!  =====================================================================
      SUBROUTINE CUNCSD2BY1(Jobu1,Jobu2,Jobv1t,M,P,Q,X11,Ldx11,X21,     &
     &                      Ldx21,Theta,U1,Ldu1,U2,Ldu2,V1t,Ldv1t,Work, &
     &                      Lwork,Rwork,Lrwork,Iwork,Info)
      USE S_CBBCSD
      USE S_CCOPY
      USE S_CLACPY
      USE S_CLAPMR
      USE S_CLAPMT
      USE S_CUNBDB1
      USE S_CUNBDB2
      USE S_CUNBDB3
      USE S_CUNBDB4
      USE S_CUNGLQ
      USE S_CUNGQR
      USE S_LSAME
      USE S_XERBLA
      IMPLICIT NONE
!*--CUNCSD2BY1274
!
! PARAMETER definitions rewritten by SPAG
!
      COMPLEX , PARAMETER  ::  ONE = (1.0E0,0.0E0) ,                    &
     &                         ZERO = (0.0E0,0.0E0)
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Jobu1
      CHARACTER :: Jobu2
      CHARACTER :: Jobv1t
      INTEGER :: M
      INTEGER :: P
      INTEGER :: Q
      COMPLEX , DIMENSION(Ldx11,*) :: X11
      INTEGER :: Ldx11
      COMPLEX , DIMENSION(Ldx21,*) :: X21
      INTEGER :: Ldx21
      REAL , DIMENSION(*) :: Theta
      COMPLEX , DIMENSION(Ldu1,*) :: U1
      INTEGER :: Ldu1
      COMPLEX , DIMENSION(Ldu2,*) :: U2
      INTEGER :: Ldu2
      COMPLEX , DIMENSION(Ldv1t,*) :: V1t
      INTEGER :: Ldv1t
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , INTENT(IN) :: Lwork
      REAL , INTENT(INOUT) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(IN) :: Lrwork
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      COMPLEX , DIMENSION(1,1) :: cdum
      INTEGER :: childinfo , i , ib11d , ib11e , ib12d , ib12e , ib21d ,&
     &           ib21e , ib22d , ib22e , ibbcsd , iorbdb , iorglq ,     &
     &           iorgqr , iphi , itaup1 , itaup2 , itauq1 , j , lbbcsd ,&
     &           lorbdb , lorglq , lorglqmin , lorglqopt , lorgqr ,     &
     &           lorgqrmin , lorgqropt , lrworkmin , lrworkopt ,        &
     &           lworkmin , lworkopt , r
      REAL , DIMENSION(1) :: dum
      LOGICAL :: lquery , wantu1 , wantu2 , wantv1t
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
!     .. Intrinsic Function ..
!     ..
!     .. Executable Statements ..
!
!     Test input arguments
!
      Info = 0
      wantu1 = LSAME(Jobu1,'Y')
      wantu2 = LSAME(Jobu2,'Y')
      wantv1t = LSAME(Jobv1t,'Y')
      lquery = (Lwork==-1) .OR. (Lrwork==-1)
!
      IF ( M<0 ) THEN
         Info = -4
      ELSEIF ( P<0 .OR. P>M ) THEN
         Info = -5
      ELSEIF ( Q<0 .OR. Q>M ) THEN
         Info = -6
      ELSEIF ( Ldx11<MAX(1,P) ) THEN
         Info = -8
      ELSEIF ( Ldx21<MAX(1,M-P) ) THEN
         Info = -10
      ELSEIF ( wantu1 .AND. Ldu1<MAX(1,P) ) THEN
         Info = -13
      ELSEIF ( wantu2 .AND. Ldu2<MAX(1,M-P) ) THEN
         Info = -15
      ELSEIF ( wantv1t .AND. Ldv1t<MAX(1,Q) ) THEN
         Info = -17
      ENDIF
!
      r = MIN(P,M-P,Q,M-Q)
!
!     Compute workspace
!
!       WORK layout:
!     |-----------------------------------------|
!     | LWORKOPT (1)                            |
!     |-----------------------------------------|
!     | TAUP1 (MAX(1,P))                        |
!     | TAUP2 (MAX(1,M-P))                      |
!     | TAUQ1 (MAX(1,Q))                        |
!     |-----------------------------------------|
!     | CUNBDB WORK | CUNGQR WORK | CUNGLQ WORK |
!     |             |             |             |
!     |             |             |             |
!     |             |             |             |
!     |             |             |             |
!     |-----------------------------------------|
!       RWORK layout:
!     |------------------|
!     | LRWORKOPT (1)    |
!     |------------------|
!     | PHI (MAX(1,R-1)) |
!     |------------------|
!     | B11D (R)         |
!     | B11E (R-1)       |
!     | B12D (R)         |
!     | B12E (R-1)       |
!     | B21D (R)         |
!     | B21E (R-1)       |
!     | B22D (R)         |
!     | B22E (R-1)       |
!     | CBBCSD RWORK     |
!     |------------------|
!
      IF ( Info==0 ) THEN
         iphi = 2
         ib11d = iphi + MAX(1,r-1)
         ib11e = ib11d + MAX(1,r)
         ib12d = ib11e + MAX(1,r-1)
         ib12e = ib12d + MAX(1,r)
         ib21d = ib12e + MAX(1,r-1)
         ib21e = ib21d + MAX(1,r)
         ib22d = ib21e + MAX(1,r-1)
         ib22e = ib22d + MAX(1,r)
         ibbcsd = ib22e + MAX(1,r-1)
         itaup1 = 2
         itaup2 = itaup1 + MAX(1,P)
         itauq1 = itaup2 + MAX(1,M-P)
         iorbdb = itauq1 + MAX(1,Q)
         iorgqr = itauq1 + MAX(1,Q)
         iorglq = itauq1 + MAX(1,Q)
         lorgqrmin = 1
         lorgqropt = 1
         lorglqmin = 1
         lorglqopt = 1
         IF ( r==Q ) THEN
            CALL CUNBDB1(M,P,Q,X11,Ldx11,X21,Ldx21,Theta,dum,cdum,cdum, &
     &                   cdum,Work,-1,childinfo)
            lorbdb = INT(Work(1))
            IF ( wantu1 .AND. P>0 ) THEN
               CALL CUNGQR(P,P,Q,U1,Ldu1,cdum,Work(1),-1,childinfo)
               lorgqrmin = MAX(lorgqrmin,P)
               lorgqropt = MAX(lorgqropt,INT(Work(1)))
            ENDIF
            IF ( wantu2 .AND. M>P ) THEN
               CALL CUNGQR(M-P,M-P,Q,U2,Ldu2,cdum,Work(1),-1,childinfo)
               lorgqrmin = MAX(lorgqrmin,M-P)
               lorgqropt = MAX(lorgqropt,INT(Work(1)))
            ENDIF
            IF ( wantv1t .AND. Q>0 ) THEN
               CALL CUNGLQ(Q-1,Q-1,Q-1,V1t,Ldv1t,cdum,Work(1),-1,       &
     &                     childinfo)
               lorglqmin = MAX(lorglqmin,Q-1)
               lorglqopt = MAX(lorglqopt,INT(Work(1)))
            ENDIF
            CALL CBBCSD(Jobu1,Jobu2,Jobv1t,'N','N',M,P,Q,Theta,dum(1),  &
     &                  U1,Ldu1,U2,Ldu2,V1t,Ldv1t,cdum,1,dum,dum,dum,   &
     &                  dum,dum,dum,dum,dum,Rwork(1),-1,childinfo)
            lbbcsd = INT(Rwork(1))
         ELSEIF ( r==P ) THEN
            CALL CUNBDB2(M,P,Q,X11,Ldx11,X21,Ldx21,Theta,dum,cdum,cdum, &
     &                   cdum,Work(1),-1,childinfo)
            lorbdb = INT(Work(1))
            IF ( wantu1 .AND. P>0 ) THEN
               CALL CUNGQR(P-1,P-1,P-1,U1(2,2),Ldu1,cdum,Work(1),-1,    &
     &                     childinfo)
               lorgqrmin = MAX(lorgqrmin,P-1)
               lorgqropt = MAX(lorgqropt,INT(Work(1)))
            ENDIF
            IF ( wantu2 .AND. M>P ) THEN
               CALL CUNGQR(M-P,M-P,Q,U2,Ldu2,cdum,Work(1),-1,childinfo)
               lorgqrmin = MAX(lorgqrmin,M-P)
               lorgqropt = MAX(lorgqropt,INT(Work(1)))
            ENDIF
            IF ( wantv1t .AND. Q>0 ) THEN
               CALL CUNGLQ(Q,Q,r,V1t,Ldv1t,cdum,Work(1),-1,childinfo)
               lorglqmin = MAX(lorglqmin,Q)
               lorglqopt = MAX(lorglqopt,INT(Work(1)))
            ENDIF
            CALL CBBCSD(Jobv1t,'N',Jobu1,Jobu2,'T',M,Q,P,Theta,dum,V1t, &
     &                  Ldv1t,cdum,1,U1,Ldu1,U2,Ldu2,dum,dum,dum,dum,   &
     &                  dum,dum,dum,dum,Rwork(1),-1,childinfo)
            lbbcsd = INT(Rwork(1))
         ELSEIF ( r==M-P ) THEN
            CALL CUNBDB3(M,P,Q,X11,Ldx11,X21,Ldx21,Theta,dum,cdum,cdum, &
     &                   cdum,Work(1),-1,childinfo)
            lorbdb = INT(Work(1))
            IF ( wantu1 .AND. P>0 ) THEN
               CALL CUNGQR(P,P,Q,U1,Ldu1,cdum,Work(1),-1,childinfo)
               lorgqrmin = MAX(lorgqrmin,P)
               lorgqropt = MAX(lorgqropt,INT(Work(1)))
            ENDIF
            IF ( wantu2 .AND. M>P ) THEN
               CALL CUNGQR(M-P-1,M-P-1,M-P-1,U2(2,2),Ldu2,cdum,Work(1), &
     &                     -1,childinfo)
               lorgqrmin = MAX(lorgqrmin,M-P-1)
               lorgqropt = MAX(lorgqropt,INT(Work(1)))
            ENDIF
            IF ( wantv1t .AND. Q>0 ) THEN
               CALL CUNGLQ(Q,Q,r,V1t,Ldv1t,cdum,Work(1),-1,childinfo)
               lorglqmin = MAX(lorglqmin,Q)
               lorglqopt = MAX(lorglqopt,INT(Work(1)))
            ENDIF
            CALL CBBCSD('N',Jobv1t,Jobu2,Jobu1,'T',M,M-Q,M-P,Theta,dum, &
     &                  cdum,1,V1t,Ldv1t,U2,Ldu2,U1,Ldu1,dum,dum,dum,   &
     &                  dum,dum,dum,dum,dum,Rwork(1),-1,childinfo)
            lbbcsd = INT(Rwork(1))
         ELSE
            CALL CUNBDB4(M,P,Q,X11,Ldx11,X21,Ldx21,Theta,dum,cdum,cdum, &
     &                   cdum,cdum,Work(1),-1,childinfo)
            lorbdb = M + INT(Work(1))
            IF ( wantu1 .AND. P>0 ) THEN
               CALL CUNGQR(P,P,M-Q,U1,Ldu1,cdum,Work(1),-1,childinfo)
               lorgqrmin = MAX(lorgqrmin,P)
               lorgqropt = MAX(lorgqropt,INT(Work(1)))
            ENDIF
            IF ( wantu2 .AND. M>P ) THEN
               CALL CUNGQR(M-P,M-P,M-Q,U2,Ldu2,cdum,Work(1),-1,         &
     &                     childinfo)
               lorgqrmin = MAX(lorgqrmin,M-P)
               lorgqropt = MAX(lorgqropt,INT(Work(1)))
            ENDIF
            IF ( wantv1t .AND. Q>0 ) THEN
               CALL CUNGLQ(Q,Q,Q,V1t,Ldv1t,cdum,Work(1),-1,childinfo)
               lorglqmin = MAX(lorglqmin,Q)
               lorglqopt = MAX(lorglqopt,INT(Work(1)))
            ENDIF
            CALL CBBCSD(Jobu2,Jobu1,'N',Jobv1t,'N',M,M-P,M-Q,Theta,dum, &
     &                  U2,Ldu2,U1,Ldu1,cdum,1,V1t,Ldv1t,dum,dum,dum,   &
     &                  dum,dum,dum,dum,dum,Rwork(1),-1,childinfo)
            lbbcsd = INT(Rwork(1))
         ENDIF
         lrworkmin = ibbcsd + lbbcsd - 1
         lrworkopt = lrworkmin
         Rwork(1) = lrworkopt
         lworkmin = MAX(iorbdb+lorbdb-1,iorgqr+lorgqrmin-1,             &
     &              iorglq+lorglqmin-1)
         lworkopt = MAX(iorbdb+lorbdb-1,iorgqr+lorgqropt-1,             &
     &              iorglq+lorglqopt-1)
         Work(1) = lworkopt
         IF ( Lwork<lworkmin .AND. .NOT.lquery ) Info = -19
         IF ( Lrwork<lrworkmin .AND. .NOT.lquery ) Info = -21
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('CUNCSD2BY1',-Info)
         RETURN
      ELSEIF ( lquery ) THEN
         RETURN
      ENDIF
      lorgqr = Lwork - iorgqr + 1
      lorglq = Lwork - iorglq + 1
!
!     Handle four cases separately: R = Q, R = P, R = M-P, and R = M-Q,
!     in which R = MIN(P,M-P,Q,M-Q)
!
      IF ( r==Q ) THEN
!
!        Case 1: R = Q
!
!        Simultaneously bidiagonalize X11 and X21
!
         CALL CUNBDB1(M,P,Q,X11,Ldx11,X21,Ldx21,Theta,Rwork(iphi),      &
     &                Work(itaup1),Work(itaup2),Work(itauq1),           &
     &                Work(iorbdb),lorbdb,childinfo)
!
!        Accumulate Householder reflectors
!
         IF ( wantu1 .AND. P>0 ) THEN
            CALL CLACPY('L',P,Q,X11,Ldx11,U1,Ldu1)
            CALL CUNGQR(P,P,Q,U1,Ldu1,Work(itaup1),Work(iorgqr),lorgqr, &
     &                  childinfo)
         ENDIF
         IF ( wantu2 .AND. M>P ) THEN
            CALL CLACPY('L',M-P,Q,X21,Ldx21,U2,Ldu2)
            CALL CUNGQR(M-P,M-P,Q,U2,Ldu2,Work(itaup2),Work(iorgqr),    &
     &                  lorgqr,childinfo)
         ENDIF
         IF ( wantv1t .AND. Q>0 ) THEN
            V1t(1,1) = ONE
            DO j = 2 , Q
               V1t(1,j) = ZERO
               V1t(j,1) = ZERO
            ENDDO
            CALL CLACPY('U',Q-1,Q-1,X21(1,2),Ldx21,V1t(2,2),Ldv1t)
            CALL CUNGLQ(Q-1,Q-1,Q-1,V1t(2,2),Ldv1t,Work(itauq1),        &
     &                  Work(iorglq),lorglq,childinfo)
         ENDIF
!
!        Simultaneously diagonalize X11 and X21.
!
         CALL CBBCSD(Jobu1,Jobu2,Jobv1t,'N','N',M,P,Q,Theta,Rwork(iphi),&
     &               U1,Ldu1,U2,Ldu2,V1t,Ldv1t,cdum,1,Rwork(ib11d),     &
     &               Rwork(ib11e),Rwork(ib12d),Rwork(ib12e),Rwork(ib21d)&
     &               ,Rwork(ib21e),Rwork(ib22d),Rwork(ib22e),           &
     &               Rwork(ibbcsd),Lrwork-ibbcsd+1,childinfo)
!
!        Permute rows and columns to place zero submatrices in
!        preferred positions
!
         IF ( Q>0 .AND. wantu2 ) THEN
            DO i = 1 , Q
               Iwork(i) = M - P - Q + i
            ENDDO
            DO i = Q + 1 , M - P
               Iwork(i) = i - Q
            ENDDO
            CALL CLAPMT(.FALSE.,M-P,M-P,U2,Ldu2,Iwork)
         ENDIF
      ELSEIF ( r==P ) THEN
!
!        Case 2: R = P
!
!        Simultaneously bidiagonalize X11 and X21
!
         CALL CUNBDB2(M,P,Q,X11,Ldx11,X21,Ldx21,Theta,Rwork(iphi),      &
     &                Work(itaup1),Work(itaup2),Work(itauq1),           &
     &                Work(iorbdb),lorbdb,childinfo)
!
!        Accumulate Householder reflectors
!
         IF ( wantu1 .AND. P>0 ) THEN
            U1(1,1) = ONE
            DO j = 2 , P
               U1(1,j) = ZERO
               U1(j,1) = ZERO
            ENDDO
            CALL CLACPY('L',P-1,P-1,X11(2,1),Ldx11,U1(2,2),Ldu1)
            CALL CUNGQR(P-1,P-1,P-1,U1(2,2),Ldu1,Work(itaup1),          &
     &                  Work(iorgqr),lorgqr,childinfo)
         ENDIF
         IF ( wantu2 .AND. M>P ) THEN
            CALL CLACPY('L',M-P,Q,X21,Ldx21,U2,Ldu2)
            CALL CUNGQR(M-P,M-P,Q,U2,Ldu2,Work(itaup2),Work(iorgqr),    &
     &                  lorgqr,childinfo)
         ENDIF
         IF ( wantv1t .AND. Q>0 ) THEN
            CALL CLACPY('U',P,Q,X11,Ldx11,V1t,Ldv1t)
            CALL CUNGLQ(Q,Q,r,V1t,Ldv1t,Work(itauq1),Work(iorglq),      &
     &                  lorglq,childinfo)
         ENDIF
!
!        Simultaneously diagonalize X11 and X21.
!
         CALL CBBCSD(Jobv1t,'N',Jobu1,Jobu2,'T',M,Q,P,Theta,Rwork(iphi),&
     &               V1t,Ldv1t,cdum,1,U1,Ldu1,U2,Ldu2,Rwork(ib11d),     &
     &               Rwork(ib11e),Rwork(ib12d),Rwork(ib12e),Rwork(ib21d)&
     &               ,Rwork(ib21e),Rwork(ib22d),Rwork(ib22e),           &
     &               Rwork(ibbcsd),lbbcsd,childinfo)
!
!        Permute rows and columns to place identity submatrices in
!        preferred positions
!
         IF ( Q>0 .AND. wantu2 ) THEN
            DO i = 1 , Q
               Iwork(i) = M - P - Q + i
            ENDDO
            DO i = Q + 1 , M - P
               Iwork(i) = i - Q
            ENDDO
            CALL CLAPMT(.FALSE.,M-P,M-P,U2,Ldu2,Iwork)
         ENDIF
      ELSEIF ( r==M-P ) THEN
!
!        Case 3: R = M-P
!
!        Simultaneously bidiagonalize X11 and X21
!
         CALL CUNBDB3(M,P,Q,X11,Ldx11,X21,Ldx21,Theta,Rwork(iphi),      &
     &                Work(itaup1),Work(itaup2),Work(itauq1),           &
     &                Work(iorbdb),lorbdb,childinfo)
!
!        Accumulate Householder reflectors
!
         IF ( wantu1 .AND. P>0 ) THEN
            CALL CLACPY('L',P,Q,X11,Ldx11,U1,Ldu1)
            CALL CUNGQR(P,P,Q,U1,Ldu1,Work(itaup1),Work(iorgqr),lorgqr, &
     &                  childinfo)
         ENDIF
         IF ( wantu2 .AND. M>P ) THEN
            U2(1,1) = ONE
            DO j = 2 , M - P
               U2(1,j) = ZERO
               U2(j,1) = ZERO
            ENDDO
            CALL CLACPY('L',M-P-1,M-P-1,X21(2,1),Ldx21,U2(2,2),Ldu2)
            CALL CUNGQR(M-P-1,M-P-1,M-P-1,U2(2,2),Ldu2,Work(itaup2),    &
     &                  Work(iorgqr),lorgqr,childinfo)
         ENDIF
         IF ( wantv1t .AND. Q>0 ) THEN
            CALL CLACPY('U',M-P,Q,X21,Ldx21,V1t,Ldv1t)
            CALL CUNGLQ(Q,Q,r,V1t,Ldv1t,Work(itauq1),Work(iorglq),      &
     &                  lorglq,childinfo)
         ENDIF
!
!        Simultaneously diagonalize X11 and X21.
!
         CALL CBBCSD('N',Jobv1t,Jobu2,Jobu1,'T',M,M-Q,M-P,Theta,        &
     &               Rwork(iphi),cdum,1,V1t,Ldv1t,U2,Ldu2,U1,Ldu1,      &
     &               Rwork(ib11d),Rwork(ib11e),Rwork(ib12d),Rwork(ib12e)&
     &               ,Rwork(ib21d),Rwork(ib21e),Rwork(ib22d),           &
     &               Rwork(ib22e),Rwork(ibbcsd),lbbcsd,childinfo)
!
!        Permute rows and columns to place identity submatrices in
!        preferred positions
!
         IF ( Q>r ) THEN
            DO i = 1 , r
               Iwork(i) = Q - r + i
            ENDDO
            DO i = r + 1 , Q
               Iwork(i) = i - r
            ENDDO
            IF ( wantu1 ) CALL CLAPMT(.FALSE.,P,Q,U1,Ldu1,Iwork)
            IF ( wantv1t ) CALL CLAPMR(.FALSE.,Q,Q,V1t,Ldv1t,Iwork)
         ENDIF
      ELSE
!
!        Case 4: R = M-Q
!
!        Simultaneously bidiagonalize X11 and X21
!
         CALL CUNBDB4(M,P,Q,X11,Ldx11,X21,Ldx21,Theta,Rwork(iphi),      &
     &                Work(itaup1),Work(itaup2),Work(itauq1),           &
     &                Work(iorbdb),Work(iorbdb+M),lorbdb-M,childinfo)
!
!        Accumulate Householder reflectors
!
 
         IF ( wantu2 .AND. M>P ) CALL CCOPY(M-P,Work(iorbdb+P),1,U2,1)
         IF ( wantu1 .AND. P>0 ) THEN
            CALL CCOPY(P,Work(iorbdb),1,U1,1)
            DO j = 2 , P
               U1(1,j) = ZERO
            ENDDO
            CALL CLACPY('L',P-1,M-Q-1,X11(2,1),Ldx11,U1(2,2),Ldu1)
            CALL CUNGQR(P,P,M-Q,U1,Ldu1,Work(itaup1),Work(iorgqr),      &
     &                  lorgqr,childinfo)
         ENDIF
         IF ( wantu2 .AND. M>P ) THEN
            DO j = 2 , M - P
               U2(1,j) = ZERO
            ENDDO
            CALL CLACPY('L',M-P-1,M-Q-1,X21(2,1),Ldx21,U2(2,2),Ldu2)
            CALL CUNGQR(M-P,M-P,M-Q,U2,Ldu2,Work(itaup2),Work(iorgqr),  &
     &                  lorgqr,childinfo)
         ENDIF
         IF ( wantv1t .AND. Q>0 ) THEN
            CALL CLACPY('U',M-Q,Q,X21,Ldx21,V1t,Ldv1t)
            CALL CLACPY('U',P-(M-Q),Q-(M-Q),X11(M-Q+1,M-Q+1),Ldx11,     &
     &                  V1t(M-Q+1,M-Q+1),Ldv1t)
            CALL CLACPY('U',-P+Q,Q-P,X21(M-Q+1,P+1),Ldx21,V1t(P+1,P+1), &
     &                  Ldv1t)
            CALL CUNGLQ(Q,Q,Q,V1t,Ldv1t,Work(itauq1),Work(iorglq),      &
     &                  lorglq,childinfo)
         ENDIF
!
!        Simultaneously diagonalize X11 and X21.
!
         CALL CBBCSD(Jobu2,Jobu1,'N',Jobv1t,'N',M,M-P,M-Q,Theta,        &
     &               Rwork(iphi),U2,Ldu2,U1,Ldu1,cdum,1,V1t,Ldv1t,      &
     &               Rwork(ib11d),Rwork(ib11e),Rwork(ib12d),Rwork(ib12e)&
     &               ,Rwork(ib21d),Rwork(ib21e),Rwork(ib22d),           &
     &               Rwork(ib22e),Rwork(ibbcsd),lbbcsd,childinfo)
!
!        Permute rows and columns to place identity submatrices in
!        preferred positions
!
         IF ( P>r ) THEN
            DO i = 1 , r
               Iwork(i) = P - r + i
            ENDDO
            DO i = r + 1 , P
               Iwork(i) = i - r
            ENDDO
            IF ( wantu1 ) CALL CLAPMT(.FALSE.,P,P,U1,Ldu1,Iwork)
            IF ( wantv1t ) CALL CLAPMR(.FALSE.,P,Q,V1t,Ldv1t,Iwork)
         ENDIF
      ENDIF
!
!
!     End of CUNCSD2BY1
!
      END SUBROUTINE CUNCSD2BY1
