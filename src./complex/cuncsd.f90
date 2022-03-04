!*==cuncsd.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
 
!> \brief \b CUNCSD
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CUNCSD + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cuncsd.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cuncsd.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cuncsd.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       RECURSIVE SUBROUTINE CUNCSD( JOBU1, JOBU2, JOBV1T, JOBV2T, TRANS,
!                                    SIGNS, M, P, Q, X11, LDX11, X12,
!                                    LDX12, X21, LDX21, X22, LDX22, THETA,
!                                    U1, LDU1, U2, LDU2, V1T, LDV1T, V2T,
!                                    LDV2T, WORK, LWORK, RWORK, LRWORK,
!                                    IWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          JOBU1, JOBU2, JOBV1T, JOBV2T, SIGNS, TRANS
!       INTEGER            INFO, LDU1, LDU2, LDV1T, LDV2T, LDX11, LDX12,
!      $                   LDX21, LDX22, LRWORK, LWORK, M, P, Q
!       ..
!       .. Array Arguments ..
!       INTEGER            IWORK( * )
!       REAL               THETA( * )
!       REAL               RWORK( * )
!       COMPLEX            U1( LDU1, * ), U2( LDU2, * ), V1T( LDV1T, * ),
!      $                   V2T( LDV2T, * ), WORK( * ), X11( LDX11, * ),
!      $                   X12( LDX12, * ), X21( LDX21, * ), X22( LDX22,
!      $                   * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CUNCSD computes the CS decomposition of an M-by-M partitioned
!> unitary matrix X:
!>
!>                                 [  I  0  0 |  0  0  0 ]
!>                                 [  0  C  0 |  0 -S  0 ]
!>     [ X11 | X12 ]   [ U1 |    ] [  0  0  0 |  0  0 -I ] [ V1 |    ]**H
!> X = [-----------] = [---------] [---------------------] [---------]   .
!>     [ X21 | X22 ]   [    | U2 ] [  0  0  0 |  I  0  0 ] [    | V2 ]
!>                                 [  0  S  0 |  0  C  0 ]
!>                                 [  0  0  I |  0  0  0 ]
!>
!> X11 is P-by-Q. The unitary matrices U1, U2, V1, and V2 are P-by-P,
!> (M-P)-by-(M-P), Q-by-Q, and (M-Q)-by-(M-Q), respectively. C and S are
!> R-by-R nonnegative diagonal matrices satisfying C^2 + S^2 = I, in
!> which R = MIN(P,M-P,Q,M-Q).
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
!> \param[in] JOBV2T
!> \verbatim
!>          JOBV2T is CHARACTER
!>          = 'Y':      V2T is computed;
!>          otherwise:  V2T is not computed.
!> \endverbatim
!>
!> \param[in] TRANS
!> \verbatim
!>          TRANS is CHARACTER
!>          = 'T':      X, U1, U2, V1T, and V2T are stored in row-major
!>                      order;
!>          otherwise:  X, U1, U2, V1T, and V2T are stored in column-
!>                      major order.
!> \endverbatim
!>
!> \param[in] SIGNS
!> \verbatim
!>          SIGNS is CHARACTER
!>          = 'O':      The lower-left block is made nonpositive (the
!>                      "other" convention);
!>          otherwise:  The upper-right block is made nonpositive (the
!>                      "default" convention).
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows and columns in X.
!> \endverbatim
!>
!> \param[in] P
!> \verbatim
!>          P is INTEGER
!>          The number of rows in X11 and X12. 0 <= P <= M.
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
!> \param[in,out] X12
!> \verbatim
!>          X12 is COMPLEX array, dimension (LDX12,M-Q)
!>          On entry, part of the unitary matrix whose CSD is desired.
!> \endverbatim
!>
!> \param[in] LDX12
!> \verbatim
!>          LDX12 is INTEGER
!>          The leading dimension of X12. LDX12 >= MAX(1,P).
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
!>          The leading dimension of X11. LDX21 >= MAX(1,M-P).
!> \endverbatim
!>
!> \param[in,out] X22
!> \verbatim
!>          X22 is COMPLEX array, dimension (LDX22,M-Q)
!>          On entry, part of the unitary matrix whose CSD is desired.
!> \endverbatim
!>
!> \param[in] LDX22
!> \verbatim
!>          LDX22 is INTEGER
!>          The leading dimension of X11. LDX22 >= MAX(1,M-P).
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
!>          U1 is COMPLEX array, dimension (LDU1,P)
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
!>          U2 is COMPLEX array, dimension (LDU2,M-P)
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
!>          V1T is COMPLEX array, dimension (LDV1T,Q)
!>          If JOBV1T = 'Y', V1T contains the Q-by-Q matrix unitary
!>          matrix V1**H.
!> \endverbatim
!>
!> \param[in] LDV1T
!> \verbatim
!>          LDV1T is INTEGER
!>          The leading dimension of V1T. If JOBV1T = 'Y', LDV1T >=
!>          MAX(1,Q).
!> \endverbatim
!>
!> \param[out] V2T
!> \verbatim
!>          V2T is COMPLEX array, dimension (LDV2T,M-Q)
!>          If JOBV2T = 'Y', V2T contains the (M-Q)-by-(M-Q) unitary
!>          matrix V2**H.
!> \endverbatim
!>
!> \param[in] LDV2T
!> \verbatim
!>          LDV2T is INTEGER
!>          The leading dimension of V2T. If JOBV2T = 'Y', LDV2T >=
!>          MAX(1,M-Q).
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
!>          only calculates the optimal size of the WORK array, returns
!>          this value as the first entry of the work array, and no error
!>          message related to LWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is REAL array, dimension MAX(1,LRWORK)
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
!>          If LRWORK = -1, then a workspace query is assumed; the routine
!>          only calculates the optimal size of the RWORK array, returns
!>          this value as the first entry of the work array, and no error
!>          message related to LRWORK is issued by XERBLA.
!> \endverbatim
!>
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
!>          > 0:  CBBCSD did not converge. See the description of RWORK
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
      RECURSIVE SUBROUTINE CUNCSD(Jobu1,Jobu2,Jobv1t,Jobv2t,Trans,Signs,&
     &                            M,P,Q,X11,Ldx11,X12,Ldx12,X21,Ldx21,  &
     &                            X22,Ldx22,Theta,U1,Ldu1,U2,Ldu2,V1t,  &
     &                            Ldv1t,V2t,Ldv2t,Work,Lwork,Rwork,     &
     &                            Lrwork,Iwork,Info)
      USE S_CBBCSD
      USE S_CLACPY
      USE S_CLAPMR
      USE S_CLAPMT
      USE S_CUNBDB
      USE S_CUNGLQ
      USE S_CUNGQR
      USE S_LSAME
      USE S_XERBLA
      IMPLICIT NONE
!*--CUNCSD333
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
      CHARACTER :: Jobv2t
      CHARACTER :: Trans
      CHARACTER :: Signs
      INTEGER :: M
      INTEGER :: P
      INTEGER :: Q
      COMPLEX , DIMENSION(Ldx11,*) :: X11
      INTEGER :: Ldx11
      COMPLEX , DIMENSION(Ldx12,*) :: X12
      INTEGER :: Ldx12
      COMPLEX , DIMENSION(Ldx21,*) :: X21
      INTEGER :: Ldx21
      COMPLEX , DIMENSION(Ldx22,*) :: X22
      INTEGER :: Ldx22
      REAL , DIMENSION(*) :: Theta
      COMPLEX , DIMENSION(Ldu1,*) :: U1
      INTEGER :: Ldu1
      COMPLEX , DIMENSION(Ldu2,*) :: U2
      INTEGER :: Ldu2
      COMPLEX , DIMENSION(Ldv1t,*) :: V1t
      INTEGER :: Ldv1t
      COMPLEX , DIMENSION(Ldv2t,*) :: V2t
      INTEGER :: Ldv2t
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      REAL , INTENT(INOUT) , DIMENSION(*) :: Rwork
      INTEGER :: Lrwork
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      INTEGER :: childinfo , i , ib11d , ib11e , ib12d , ib12e , ib21d ,&
     &           ib21e , ib22d , ib22e , ibbcsd , iorbdb , iorglq ,     &
     &           iorgqr , iphi , itaup1 , itaup2 , itauq1 , itauq2 , j ,&
     &           lbbcsdwork , lbbcsdworkmin , lbbcsdworkopt ,           &
     &           lorbdbwork , lorbdbworkmin , lorbdbworkopt ,           &
     &           lorglqwork , lorglqworkmin , lorglqworkopt ,           &
     &           lorgqrwork , lorgqrworkmin , lorgqrworkopt ,           &
     &           lrworkmin , lrworkopt , lworkmin , lworkopt , p1 , q1
      LOGICAL :: colmajor , defaultsigns , lquery , lrquery , wantu1 ,  &
     &           wantu2 , wantv1t , wantv2t
      CHARACTER :: signst , transt
!
! End of declarations rewritten by SPAG
!
!     ..
!     .. Array Arguments ..
!     ..
!
!  ===================================================================
!
!     .. Parameters ..
!     ..
!     .. Local Scalars ..
!     ..
!     .. External Subroutines ..
!     ..
!     .. External Functions ..
!     ..
!     .. Intrinsic Functions
!     ..
!     .. Executable Statements ..
!
!     Test input arguments
!
      Info = 0
      wantu1 = LSAME(Jobu1,'Y')
      wantu2 = LSAME(Jobu2,'Y')
      wantv1t = LSAME(Jobv1t,'Y')
      wantv2t = LSAME(Jobv2t,'Y')
      colmajor = .NOT.LSAME(Trans,'T')
      defaultsigns = .NOT.LSAME(Signs,'O')
      lquery = Lwork== - 1
      lrquery = Lrwork== - 1
      IF ( M<0 ) THEN
         Info = -7
      ELSEIF ( P<0 .OR. P>M ) THEN
         Info = -8
      ELSEIF ( Q<0 .OR. Q>M ) THEN
         Info = -9
      ELSEIF ( colmajor .AND. Ldx11<MAX(1,P) ) THEN
         Info = -11
      ELSEIF ( .NOT.colmajor .AND. Ldx11<MAX(1,Q) ) THEN
         Info = -11
      ELSEIF ( colmajor .AND. Ldx12<MAX(1,P) ) THEN
         Info = -13
      ELSEIF ( .NOT.colmajor .AND. Ldx12<MAX(1,M-Q) ) THEN
         Info = -13
      ELSEIF ( colmajor .AND. Ldx21<MAX(1,M-P) ) THEN
         Info = -15
      ELSEIF ( .NOT.colmajor .AND. Ldx21<MAX(1,Q) ) THEN
         Info = -15
      ELSEIF ( colmajor .AND. Ldx22<MAX(1,M-P) ) THEN
         Info = -17
      ELSEIF ( .NOT.colmajor .AND. Ldx22<MAX(1,M-Q) ) THEN
         Info = -17
      ELSEIF ( wantu1 .AND. Ldu1<P ) THEN
         Info = -20
      ELSEIF ( wantu2 .AND. Ldu2<M-P ) THEN
         Info = -22
      ELSEIF ( wantv1t .AND. Ldv1t<Q ) THEN
         Info = -24
      ELSEIF ( wantv2t .AND. Ldv2t<M-Q ) THEN
         Info = -26
      ENDIF
!
!     Work with transpose if convenient
!
      IF ( Info==0 .AND. MIN(P,M-P)<MIN(Q,M-Q) ) THEN
         IF ( colmajor ) THEN
            transt = 'T'
         ELSE
            transt = 'N'
         ENDIF
         IF ( defaultsigns ) THEN
            signst = 'O'
         ELSE
            signst = 'D'
         ENDIF
         CALL CUNCSD(Jobv1t,Jobv2t,Jobu1,Jobu2,transt,signst,M,Q,P,X11, &
     &               Ldx11,X21,Ldx21,X12,Ldx12,X22,Ldx22,Theta,V1t,     &
     &               Ldv1t,V2t,Ldv2t,U1,Ldu1,U2,Ldu2,Work,Lwork,Rwork,  &
     &               Lrwork,Iwork,Info)
         RETURN
      ENDIF
!
!     Work with permutation [ 0 I; I 0 ] * X * [ 0 I; I 0 ] if
!     convenient
!
      IF ( Info==0 .AND. M-Q<Q ) THEN
         IF ( defaultsigns ) THEN
            signst = 'O'
         ELSE
            signst = 'D'
         ENDIF
         CALL CUNCSD(Jobu2,Jobu1,Jobv2t,Jobv1t,Trans,signst,M,M-P,M-Q,  &
     &               X22,Ldx22,X21,Ldx21,X12,Ldx12,X11,Ldx11,Theta,U2,  &
     &               Ldu2,U1,Ldu1,V2t,Ldv2t,V1t,Ldv1t,Work,Lwork,Rwork, &
     &               Lrwork,Iwork,Info)
         RETURN
      ENDIF
!
!     Compute workspace
!
      IF ( Info==0 ) THEN
!
!        Real workspace
!
         iphi = 2
         ib11d = iphi + MAX(1,Q-1)
         ib11e = ib11d + MAX(1,Q)
         ib12d = ib11e + MAX(1,Q-1)
         ib12e = ib12d + MAX(1,Q)
         ib21d = ib12e + MAX(1,Q-1)
         ib21e = ib21d + MAX(1,Q)
         ib22d = ib21e + MAX(1,Q-1)
         ib22e = ib22d + MAX(1,Q)
         ibbcsd = ib22e + MAX(1,Q-1)
         CALL CBBCSD(Jobu1,Jobu2,Jobv1t,Jobv2t,Trans,M,P,Q,Theta,Theta, &
     &               U1,Ldu1,U2,Ldu2,V1t,Ldv1t,V2t,Ldv2t,Theta,Theta,   &
     &               Theta,Theta,Theta,Theta,Theta,Theta,Rwork,-1,      &
     &               childinfo)
         lbbcsdworkopt = INT(Rwork(1))
         lbbcsdworkmin = lbbcsdworkopt
         lrworkopt = ibbcsd + lbbcsdworkopt - 1
         lrworkmin = ibbcsd + lbbcsdworkmin - 1
         Rwork(1) = lrworkopt
!
!        Complex workspace
!
         itaup1 = 2
         itaup2 = itaup1 + MAX(1,P)
         itauq1 = itaup2 + MAX(1,M-P)
         itauq2 = itauq1 + MAX(1,Q)
         iorgqr = itauq2 + MAX(1,M-Q)
         CALL CUNGQR(M-Q,M-Q,M-Q,U1,MAX(1,M-Q),U1,Work,-1,childinfo)
         lorgqrworkopt = INT(Work(1))
         lorgqrworkmin = MAX(1,M-Q)
         iorglq = itauq2 + MAX(1,M-Q)
         CALL CUNGLQ(M-Q,M-Q,M-Q,U1,MAX(1,M-Q),U1,Work,-1,childinfo)
         lorglqworkopt = INT(Work(1))
         lorglqworkmin = MAX(1,M-Q)
         iorbdb = itauq2 + MAX(1,M-Q)
         CALL CUNBDB(Trans,Signs,M,P,Q,X11,Ldx11,X12,Ldx12,X21,Ldx21,   &
     &               X22,Ldx22,Theta,Theta,U1,U2,V1t,V2t,Work,-1,       &
     &               childinfo)
         lorbdbworkopt = INT(Work(1))
         lorbdbworkmin = lorbdbworkopt
         lworkopt = MAX(iorgqr+lorgqrworkopt,iorglq+lorglqworkopt,      &
     &              iorbdb+lorbdbworkopt) - 1
         lworkmin = MAX(iorgqr+lorgqrworkmin,iorglq+lorglqworkmin,      &
     &              iorbdb+lorbdbworkmin) - 1
         Work(1) = MAX(lworkopt,lworkmin)
!
         IF ( Lwork<lworkmin .AND. .NOT.(lquery .OR. lrquery) ) THEN
            Info = -22
         ELSEIF ( Lrwork<lrworkmin .AND. .NOT.(lquery .OR. lrquery) )   &
     &            THEN
            Info = -24
         ELSE
            lorgqrwork = Lwork - iorgqr + 1
            lorglqwork = Lwork - iorglq + 1
            lorbdbwork = Lwork - iorbdb + 1
            lbbcsdwork = Lrwork - ibbcsd + 1
         ENDIF
      ENDIF
!
!     Abort if any illegal arguments
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('CUNCSD',-Info)
         RETURN
      ELSEIF ( lquery .OR. lrquery ) THEN
         RETURN
      ENDIF
!
!     Transform to bidiagonal block form
!
      CALL CUNBDB(Trans,Signs,M,P,Q,X11,Ldx11,X12,Ldx12,X21,Ldx21,X22,  &
     &            Ldx22,Theta,Rwork(iphi),Work(itaup1),Work(itaup2),    &
     &            Work(itauq1),Work(itauq2),Work(iorbdb),lorbdbwork,    &
     &            childinfo)
!
!     Accumulate Householder reflectors
!
      IF ( colmajor ) THEN
         IF ( wantu1 .AND. P>0 ) THEN
            CALL CLACPY('L',P,Q,X11,Ldx11,U1,Ldu1)
            CALL CUNGQR(P,P,Q,U1,Ldu1,Work(itaup1),Work(iorgqr),        &
     &                  lorgqrwork,Info)
         ENDIF
         IF ( wantu2 .AND. M>P ) THEN
            CALL CLACPY('L',M-P,Q,X21,Ldx21,U2,Ldu2)
            CALL CUNGQR(M-P,M-P,Q,U2,Ldu2,Work(itaup2),Work(iorgqr),    &
     &                  lorgqrwork,Info)
         ENDIF
         IF ( wantv1t .AND. Q>0 ) THEN
            CALL CLACPY('U',Q-1,Q-1,X11(1,2),Ldx11,V1t(2,2),Ldv1t)
            V1t(1,1) = ONE
            DO j = 2 , Q
               V1t(1,j) = ZERO
               V1t(j,1) = ZERO
            ENDDO
            CALL CUNGLQ(Q-1,Q-1,Q-1,V1t(2,2),Ldv1t,Work(itauq1),        &
     &                  Work(iorglq),lorglqwork,Info)
         ENDIF
         IF ( wantv2t .AND. M>Q ) THEN
            CALL CLACPY('U',P,M-Q,X12,Ldx12,V2t,Ldv2t)
            IF ( M-P>Q ) CALL CLACPY('U',M-P-Q,M-P-Q,X22(Q+1,P+1),Ldx22,&
     &                               V2t(P+1,P+1),Ldv2t)
            IF ( M>Q ) CALL CUNGLQ(M-Q,M-Q,M-Q,V2t,Ldv2t,Work(itauq2),  &
     &                             Work(iorglq),lorglqwork,Info)
         ENDIF
      ELSE
         IF ( wantu1 .AND. P>0 ) THEN
            CALL CLACPY('U',Q,P,X11,Ldx11,U1,Ldu1)
            CALL CUNGLQ(P,P,Q,U1,Ldu1,Work(itaup1),Work(iorglq),        &
     &                  lorglqwork,Info)
         ENDIF
         IF ( wantu2 .AND. M>P ) THEN
            CALL CLACPY('U',Q,M-P,X21,Ldx21,U2,Ldu2)
            CALL CUNGLQ(M-P,M-P,Q,U2,Ldu2,Work(itaup2),Work(iorglq),    &
     &                  lorglqwork,Info)
         ENDIF
         IF ( wantv1t .AND. Q>0 ) THEN
            CALL CLACPY('L',Q-1,Q-1,X11(2,1),Ldx11,V1t(2,2),Ldv1t)
            V1t(1,1) = ONE
            DO j = 2 , Q
               V1t(1,j) = ZERO
               V1t(j,1) = ZERO
            ENDDO
            CALL CUNGQR(Q-1,Q-1,Q-1,V1t(2,2),Ldv1t,Work(itauq1),        &
     &                  Work(iorgqr),lorgqrwork,Info)
         ENDIF
         IF ( wantv2t .AND. M>Q ) THEN
            p1 = MIN(P+1,M)
            q1 = MIN(Q+1,M)
            CALL CLACPY('L',M-Q,P,X12,Ldx12,V2t,Ldv2t)
            IF ( M>P+Q ) CALL CLACPY('L',M-P-Q,M-P-Q,X22(p1,q1),Ldx22,  &
     &                               V2t(P+1,P+1),Ldv2t)
            CALL CUNGQR(M-Q,M-Q,M-Q,V2t,Ldv2t,Work(itauq2),Work(iorgqr),&
     &                  lorgqrwork,Info)
         ENDIF
      ENDIF
!
!     Compute the CSD of the matrix in bidiagonal-block form
!
      CALL CBBCSD(Jobu1,Jobu2,Jobv1t,Jobv2t,Trans,M,P,Q,Theta,          &
     &            Rwork(iphi),U1,Ldu1,U2,Ldu2,V1t,Ldv1t,V2t,Ldv2t,      &
     &            Rwork(ib11d),Rwork(ib11e),Rwork(ib12d),Rwork(ib12e),  &
     &            Rwork(ib21d),Rwork(ib21e),Rwork(ib22d),Rwork(ib22e),  &
     &            Rwork(ibbcsd),lbbcsdwork,Info)
!
!     Permute rows and columns to place identity submatrices in top-
!     left corner of (1,1)-block and/or bottom-right corner of (1,2)-
!     block and/or bottom-right corner of (2,1)-block and/or top-left
!     corner of (2,2)-block
!
      IF ( Q>0 .AND. wantu2 ) THEN
         DO i = 1 , Q
            Iwork(i) = M - P - Q + i
         ENDDO
         DO i = Q + 1 , M - P
            Iwork(i) = i - Q
         ENDDO
         IF ( colmajor ) THEN
            CALL CLAPMT(.FALSE.,M-P,M-P,U2,Ldu2,Iwork)
         ELSE
            CALL CLAPMR(.FALSE.,M-P,M-P,U2,Ldu2,Iwork)
         ENDIF
      ENDIF
      IF ( M>0 .AND. wantv2t ) THEN
         DO i = 1 , P
            Iwork(i) = M - P - Q + i
         ENDDO
         DO i = P + 1 , M - Q
            Iwork(i) = i - P
         ENDDO
         IF ( .NOT.colmajor ) THEN
            CALL CLAPMT(.FALSE.,M-Q,M-Q,V2t,Ldv2t,Iwork)
         ELSE
            CALL CLAPMR(.FALSE.,M-Q,M-Q,V2t,Ldv2t,Iwork)
         ENDIF
      ENDIF
!
!
!     End CUNCSD
!
      END SUBROUTINE CUNCSD
