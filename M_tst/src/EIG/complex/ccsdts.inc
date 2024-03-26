!*==ccsdts.f90  processed by SPAG 7.51RB at 17:37 on  4 Mar 2022
!> \brief \b ccsdts
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ccsdts( M, P, Q, X, XF, LDX, U1, LDU1, U2, LDU2, V1T,
!                          LDV1T, V2T, LDV2T, THETA, IWORK, WORK, LWORK,
!                          RWORK, RESULT )
!
!       .. Scalar Arguments ..
!       INTEGER            LDX, LDU1, LDU2, LDV1T, LDV2T, LWORK, M, P, Q
!       ..
!       .. Array Arguments ..
!       INTEGER            IWORK( * )
!       REAL               RESULT( 15 ), RWORK( * ), THETA( * )
!       COMPLEX            U1( LDU1, * ), U2( LDU2, * ), V1T( LDV1T, * ),
!      $                   V2T( LDV2T, * ), WORK( LWORK ), X( LDX, * ),
!      $                   XF( LDX, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ccsdts tests CUNCSD, which, given an M-by-M partitioned unitary
!> matrix X,
!>              Q  M-Q
!>       X = [ X11 X12 ] P   ,
!>           [ X21 X22 ] M-P
!>
!> computes the CSD
!>
!>       [ U1    ]**T * [ X11 X12 ] * [ V1    ]
!>       [    U2 ]      [ X21 X22 ]   [    V2 ]
!>
!>                             [  I  0  0 |  0  0  0 ]
!>                             [  0  C  0 |  0 -S  0 ]
!>                             [  0  0  0 |  0  0 -I ]
!>                           = [---------------------] = [ D11 D12 ] .
!>                             [  0  0  0 |  I  0  0 ]   [ D21 D22 ]
!>                             [  0  S  0 |  0  C  0 ]
!>                             [  0  0  I |  0  0  0 ]
!>
!> and also SORCSD2BY1, which, given
!>          Q
!>       [ X11 ] P   ,
!>       [ X21 ] M-P
!>
!> computes the 2-by-1 CSD
!>
!>                                     [  I  0  0 ]
!>                                     [  0  C  0 ]
!>                                     [  0  0  0 ]
!>       [ U1    ]**T * [ X11 ] * V1 = [----------] = [ D11 ] ,
!>       [    U2 ]      [ X21 ]        [  0  0  0 ]   [ D21 ]
!>                                     [  0  S  0 ]
!>                                     [  0  0  I ]
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix X.  M >= 0.
!> \endverbatim
!>
!> \param[in] P
!> \verbatim
!>          P is INTEGER
!>          The number of rows of the matrix X11.  P >= 0.
!> \endverbatim
!>
!> \param[in] Q
!> \verbatim
!>          Q is INTEGER
!>          The number of columns of the matrix X11.  Q >= 0.
!> \endverbatim
!>
!> \param[in] X
!> \verbatim
!>          X is COMPLEX array, dimension (LDX,M)
!>          The M-by-M matrix X.
!> \endverbatim
!>
!> \param[out] XF
!> \verbatim
!>          XF is COMPLEX array, dimension (LDX,M)
!>          Details of the CSD of X, as returned by CUNCSD;
!>          see CUNCSD for further details.
!> \endverbatim
!>
!> \param[in] LDX
!> \verbatim
!>          LDX is INTEGER
!>          The leading dimension of the arrays X and XF.
!>          LDX >= max( 1,M ).
!> \endverbatim
!>
!> \param[out] U1
!> \verbatim
!>          U1 is COMPLEX array, dimension(LDU1,P)
!>          The P-by-P unitary matrix U1.
!> \endverbatim
!>
!> \param[in] LDU1
!> \verbatim
!>          LDU1 is INTEGER
!>          The leading dimension of the array U1. LDU >= max(1,P).
!> \endverbatim
!>
!> \param[out] U2
!> \verbatim
!>          U2 is COMPLEX array, dimension(LDU2,M-P)
!>          The (M-P)-by-(M-P) unitary matrix U2.
!> \endverbatim
!>
!> \param[in] LDU2
!> \verbatim
!>          LDU2 is INTEGER
!>          The leading dimension of the array U2. LDU >= max(1,M-P).
!> \endverbatim
!>
!> \param[out] V1T
!> \verbatim
!>          V1T is COMPLEX array, dimension(LDV1T,Q)
!>          The Q-by-Q unitary matrix V1T.
!> \endverbatim
!>
!> \param[in] LDV1T
!> \verbatim
!>          LDV1T is INTEGER
!>          The leading dimension of the array V1T. LDV1T >=
!>          max(1,Q).
!> \endverbatim
!>
!> \param[out] V2T
!> \verbatim
!>          V2T is COMPLEX array, dimension(LDV2T,M-Q)
!>          The (M-Q)-by-(M-Q) unitary matrix V2T.
!> \endverbatim
!>
!> \param[in] LDV2T
!> \verbatim
!>          LDV2T is INTEGER
!>          The leading dimension of the array V2T. LDV2T >=
!>          max(1,M-Q).
!> \endverbatim
!>
!> \param[out] THETA
!> \verbatim
!>          THETA is REAL array, dimension MIN(P,M-P,Q,M-Q)
!>          The CS values of X; the essentially diagonal matrices C and
!>          S are constructed from THETA; see subroutine CUNCSD for
!>          details.
!> \endverbatim
!>
!> \param[out] IWORK
!> \verbatim
!>          IWORK is INTEGER array, dimension (M)
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX array, dimension (LWORK)
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The dimension of the array WORK
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is REAL array
!> \endverbatim
!>
!> \param[out] RESULT
!> \verbatim
!>          RESULT is REAL array, dimension (15)
!>          The test ratios:
!>          First, the 2-by-2 CSD:
!>          RESULT(1) = norm( U1'*X11*V1 - D11 ) / ( MAX(1,P,Q)*EPS2 )
!>          RESULT(2) = norm( U1'*X12*V2 - D12 ) / ( MAX(1,P,M-Q)*EPS2 )
!>          RESULT(3) = norm( U2'*X21*V1 - D21 ) / ( MAX(1,M-P,Q)*EPS2 )
!>          RESULT(4) = norm( U2'*X22*V2 - D22 ) / ( MAX(1,M-P,M-Q)*EPS2 )
!>          RESULT(5) = norm( I - U1'*U1 ) / ( MAX(1,P)*ULP )
!>          RESULT(6) = norm( I - U2'*U2 ) / ( MAX(1,M-P)*ULP )
!>          RESULT(7) = norm( I - V1T'*V1T ) / ( MAX(1,Q)*ULP )
!>          RESULT(8) = norm( I - V2T'*V2T ) / ( MAX(1,M-Q)*ULP )
!>          RESULT(9) = 0        if THETA is in increasing order and
!>                               all angles are in [0,pi/2];
!>                    = ULPINV   otherwise.
!>          Then, the 2-by-1 CSD:
!>          RESULT(10) = norm( U1'*X11*V1 - D11 ) / ( MAX(1,P,Q)*EPS2 )
!>          RESULT(11) = norm( U2'*X21*V1 - D21 ) / ( MAX(1,M-P,Q)*EPS2 )
!>          RESULT(12) = norm( I - U1'*U1 ) / ( MAX(1,P)*ULP )
!>          RESULT(13) = norm( I - U2'*U2 ) / ( MAX(1,M-P)*ULP )
!>          RESULT(14) = norm( I - V1T'*V1T ) / ( MAX(1,Q)*ULP )
!>          RESULT(15) = 0        if THETA is in increasing order and
!>                                all angles are in [0,pi/2];
!>                     = ULPINV   otherwise.
!>          ( EPS2 = MAX( norm( I - X'*X ) / M, ULP ). )
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
!> \ingroup complex_eig
!
!  =====================================================================
      SUBROUTINE ccsdts(M,P,Q,X,Xf,Ldx,U1,Ldu1,U2,Ldu2,V1t,Ldv1t,V2t,   &
     &                  Ldv2t,Theta,Iwork,Work,Lwork,Rwork,Result)
      IMPLICIT NONE
!*--CCSDTS232
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Ldx , Ldu1 , Ldu2 , Ldv1t , Ldv2t , Lwork , M , P , Q
!     ..
!     .. Array Arguments ..
      INTEGER Iwork(*)
      REAL Result(15) , Rwork(*) , Theta(*)
      COMPLEX U1(Ldu1,*) , U2(Ldu2,*) , V1t(Ldv1t,*) , V2t(Ldv2t,*) ,   &
     &        Work(Lwork) , X(Ldx,*) , Xf(Ldx,*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL REALONE , REALZERO
      PARAMETER (REALONE=1.0E0,REALZERO=0.0E0)
      COMPLEX ZERO , ONE
      PARAMETER (ZERO=(0.0E0,0.0E0),ONE=(1.0E0,0.0E0))
      REAL PIOVER2
      PARAMETER (PIOVER2=1.57079632679489661923132169163975144210E0)
!     ..
!     .. Local Scalars ..
      INTEGER i , info , r
      REAL eps2 , resid , ulp , ulpinv
!     ..
!     .. External Functions ..
      REAL SLAMCH , CLANGE , CLANHE
      EXTERNAL SLAMCH , CLANGE , CLANHE
!     ..
!     .. External Subroutines ..
      EXTERNAL CGEMM , CHERK , CLACPY , CLASET , CUNCSD , CUNCSD2BY1
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC CMPLX , COS , MAX , MIN , REAL , SIN
!     ..
!     .. Executable Statements ..
!
      ulp = SLAMCH('Precision')
      ulpinv = REALONE/ulp
!
!     The first half of the routine checks the 2-by-2 CSD
!
      CALL CLASET('Full',M,M,ZERO,ONE,Work,Ldx)
      CALL CHERK('Upper','Conjugate transpose',M,M,-REALONE,X,Ldx,      &
     &           REALONE,Work,Ldx)
      IF ( M>0 ) THEN
         eps2 = MAX(ulp,CLANGE('1',M,M,Work,Ldx,Rwork)/REAL(M))
      ELSE
         eps2 = ulp
      ENDIF
      r = MIN(P,M-P,Q,M-Q)
!
!     Copy the matrix X to the array XF.
!
      CALL CLACPY('Full',M,M,X,Ldx,Xf,Ldx)
!
!     Compute the CSD
!
      CALL CUNCSD('Y','Y','Y','Y','N','D',M,P,Q,Xf(1,1),Ldx,Xf(1,Q+1),  &
     &            Ldx,Xf(P+1,1),Ldx,Xf(P+1,Q+1),Ldx,Theta,U1,Ldu1,U2,   &
     &            Ldu2,V1t,Ldv1t,V2t,Ldv2t,Work,Lwork,Rwork,17*(r+2),   &
     &            Iwork,info)
!
!     Compute XF := diag(U1,U2)'*X*diag(V1,V2) - [D11 D12; D21 D22]
!
      CALL CLACPY('Full',M,M,X,Ldx,Xf,Ldx)
!
      CALL CGEMM('No transpose','Conjugate transpose',P,Q,Q,ONE,Xf,Ldx, &
     &           V1t,Ldv1t,ZERO,Work,Ldx)
!
      CALL CGEMM('Conjugate transpose','No transpose',P,Q,P,ONE,U1,Ldu1,&
     &           Work,Ldx,ZERO,Xf,Ldx)
!
      DO i = 1 , MIN(P,Q) - r
         Xf(i,i) = Xf(i,i) - ONE
      ENDDO
      DO i = 1 , r
         Xf(MIN(P,Q)-r+i,MIN(P,Q)-r+i) = Xf(MIN(P,Q)-r+i,MIN(P,Q)-r+i)  &
     &      - CMPLX(COS(Theta(i)),0.0E0)
      ENDDO
!
      CALL CGEMM('No transpose','Conjugate transpose',P,M-Q,M-Q,ONE,    &
     &           Xf(1,Q+1),Ldx,V2t,Ldv2t,ZERO,Work,Ldx)
!
      CALL CGEMM('Conjugate transpose','No transpose',P,M-Q,P,ONE,U1,   &
     &           Ldu1,Work,Ldx,ZERO,Xf(1,Q+1),Ldx)
!
      DO i = 1 , MIN(P,M-Q) - r
         Xf(P-i+1,M-i+1) = Xf(P-i+1,M-i+1) + ONE
      ENDDO
      DO i = 1 , r
         Xf(P-(MIN(P,M-Q)-r)+1-i,M-(MIN(P,M-Q)-r)+1-i)                  &
     &      = Xf(P-(MIN(P,M-Q)-r)+1-i,M-(MIN(P,M-Q)-r)+1-i)             &
     &      + CMPLX(SIN(Theta(r-i+1)),0.0E0)
      ENDDO
!
      CALL CGEMM('No transpose','Conjugate transpose',M-P,Q,Q,ONE,      &
     &           Xf(P+1,1),Ldx,V1t,Ldv1t,ZERO,Work,Ldx)
!
      CALL CGEMM('Conjugate transpose','No transpose',M-P,Q,M-P,ONE,U2, &
     &           Ldu2,Work,Ldx,ZERO,Xf(P+1,1),Ldx)
!
      DO i = 1 , MIN(M-P,Q) - r
         Xf(M-i+1,Q-i+1) = Xf(M-i+1,Q-i+1) - ONE
      ENDDO
      DO i = 1 , r
         Xf(M-(MIN(M-P,Q)-r)+1-i,Q-(MIN(M-P,Q)-r)+1-i)                  &
     &      = Xf(M-(MIN(M-P,Q)-r)+1-i,Q-(MIN(M-P,Q)-r)+1-i)             &
     &      - CMPLX(SIN(Theta(r-i+1)),0.0E0)
      ENDDO
!
      CALL CGEMM('No transpose','Conjugate transpose',M-P,M-Q,M-Q,ONE,  &
     &           Xf(P+1,Q+1),Ldx,V2t,Ldv2t,ZERO,Work,Ldx)
!
      CALL CGEMM('Conjugate transpose','No transpose',M-P,M-Q,M-P,ONE,  &
     &           U2,Ldu2,Work,Ldx,ZERO,Xf(P+1,Q+1),Ldx)
!
      DO i = 1 , MIN(M-P,M-Q) - r
         Xf(P+i,Q+i) = Xf(P+i,Q+i) - ONE
      ENDDO
      DO i = 1 , r
         Xf(P+(MIN(M-P,M-Q)-r)+i,Q+(MIN(M-P,M-Q)-r)+i)                  &
     &      = Xf(P+(MIN(M-P,M-Q)-r)+i,Q+(MIN(M-P,M-Q)-r)+i)             &
     &      - CMPLX(COS(Theta(i)),0.0E0)
      ENDDO
!
!     Compute norm( U1'*X11*V1 - D11 ) / ( MAX(1,P,Q)*EPS2 ) .
!
      resid = CLANGE('1',P,Q,Xf,Ldx,Rwork)
      Result(1) = (resid/REAL(MAX(1,P,Q)))/eps2
!
!     Compute norm( U1'*X12*V2 - D12 ) / ( MAX(1,P,M-Q)*EPS2 ) .
!
      resid = CLANGE('1',P,M-Q,Xf(1,Q+1),Ldx,Rwork)
      Result(2) = (resid/REAL(MAX(1,P,M-Q)))/eps2
!
!     Compute norm( U2'*X21*V1 - D21 ) / ( MAX(1,M-P,Q)*EPS2 ) .
!
      resid = CLANGE('1',M-P,Q,Xf(P+1,1),Ldx,Rwork)
      Result(3) = (resid/REAL(MAX(1,M-P,Q)))/eps2
!
!     Compute norm( U2'*X22*V2 - D22 ) / ( MAX(1,M-P,M-Q)*EPS2 ) .
!
      resid = CLANGE('1',M-P,M-Q,Xf(P+1,Q+1),Ldx,Rwork)
      Result(4) = (resid/REAL(MAX(1,M-P,M-Q)))/eps2
!
!     Compute I - U1'*U1
!
      CALL CLASET('Full',P,P,ZERO,ONE,Work,Ldu1)
      CALL CHERK('Upper','Conjugate transpose',P,P,-REALONE,U1,Ldu1,    &
     &           REALONE,Work,Ldu1)
!
!     Compute norm( I - U'*U ) / ( MAX(1,P) * ULP ) .
!
      resid = CLANHE('1','Upper',P,Work,Ldu1,Rwork)
      Result(5) = (resid/REAL(MAX(1,P)))/ulp
!
!     Compute I - U2'*U2
!
      CALL CLASET('Full',M-P,M-P,ZERO,ONE,Work,Ldu2)
      CALL CHERK('Upper','Conjugate transpose',M-P,M-P,-REALONE,U2,Ldu2,&
     &           REALONE,Work,Ldu2)
!
!     Compute norm( I - U2'*U2 ) / ( MAX(1,M-P) * ULP ) .
!
      resid = CLANHE('1','Upper',M-P,Work,Ldu2,Rwork)
      Result(6) = (resid/REAL(MAX(1,M-P)))/ulp
!
!     Compute I - V1T*V1T'
!
      CALL CLASET('Full',Q,Q,ZERO,ONE,Work,Ldv1t)
      CALL CHERK('Upper','No transpose',Q,Q,-REALONE,V1t,Ldv1t,REALONE, &
     &           Work,Ldv1t)
!
!     Compute norm( I - V1T*V1T' ) / ( MAX(1,Q) * ULP ) .
!
      resid = CLANHE('1','Upper',Q,Work,Ldv1t,Rwork)
      Result(7) = (resid/REAL(MAX(1,Q)))/ulp
!
!     Compute I - V2T*V2T'
!
      CALL CLASET('Full',M-Q,M-Q,ZERO,ONE,Work,Ldv2t)
      CALL CHERK('Upper','No transpose',M-Q,M-Q,-REALONE,V2t,Ldv2t,     &
     &           REALONE,Work,Ldv2t)
!
!     Compute norm( I - V2T*V2T' ) / ( MAX(1,M-Q) * ULP ) .
!
      resid = CLANHE('1','Upper',M-Q,Work,Ldv2t,Rwork)
      Result(8) = (resid/REAL(MAX(1,M-Q)))/ulp
!
!     Check sorting
!
      Result(9) = REALZERO
      DO i = 1 , r
         IF ( Theta(i)<REALZERO .OR. Theta(i)>PIOVER2 ) Result(9)       &
     &        = ulpinv
         IF ( i>1 ) THEN
            IF ( Theta(i)<Theta(i-1) ) Result(9) = ulpinv
         ENDIF
      ENDDO
!
!     The second half of the routine checks the 2-by-1 CSD
!
      CALL CLASET('Full',Q,Q,ZERO,ONE,Work,Ldx)
      CALL CHERK('Upper','Conjugate transpose',Q,M,-REALONE,X,Ldx,      &
     &           REALONE,Work,Ldx)
      IF ( M>0 ) THEN
         eps2 = MAX(ulp,CLANGE('1',Q,Q,Work,Ldx,Rwork)/REAL(M))
      ELSE
         eps2 = ulp
      ENDIF
      r = MIN(P,M-P,Q,M-Q)
!
!     Copy the matrix X to the array XF.
!
      CALL CLACPY('Full',M,Q,X,Ldx,Xf,Ldx)
!
!     Compute the CSD
!
      CALL CUNCSD2BY1('Y','Y','Y',M,P,Q,Xf(1,1),Ldx,Xf(P+1,1),Ldx,Theta,&
     &                U1,Ldu1,U2,Ldu2,V1t,Ldv1t,Work,Lwork,Rwork,       &
     &                17*(r+2),Iwork,info)
!
!     Compute [X11;X21] := diag(U1,U2)'*[X11;X21]*V1 - [D11;D21]
!
      CALL CGEMM('No transpose','Conjugate transpose',P,Q,Q,ONE,X,Ldx,  &
     &           V1t,Ldv1t,ZERO,Work,Ldx)
!
      CALL CGEMM('Conjugate transpose','No transpose',P,Q,P,ONE,U1,Ldu1,&
     &           Work,Ldx,ZERO,X,Ldx)
!
      DO i = 1 , MIN(P,Q) - r
         X(i,i) = X(i,i) - ONE
      ENDDO
      DO i = 1 , r
         X(MIN(P,Q)-r+i,MIN(P,Q)-r+i) = X(MIN(P,Q)-r+i,MIN(P,Q)-r+i)    &
     &                                  - CMPLX(COS(Theta(i)),0.0E0)
      ENDDO
!
      CALL CGEMM('No transpose','Conjugate transpose',M-P,Q,Q,ONE,      &
     &           X(P+1,1),Ldx,V1t,Ldv1t,ZERO,Work,Ldx)
!
      CALL CGEMM('Conjugate transpose','No transpose',M-P,Q,M-P,ONE,U2, &
     &           Ldu2,Work,Ldx,ZERO,X(P+1,1),Ldx)
!
      DO i = 1 , MIN(M-P,Q) - r
         X(M-i+1,Q-i+1) = X(M-i+1,Q-i+1) - ONE
      ENDDO
      DO i = 1 , r
         X(M-(MIN(M-P,Q)-r)+1-i,Q-(MIN(M-P,Q)-r)+1-i)                   &
     &      = X(M-(MIN(M-P,Q)-r)+1-i,Q-(MIN(M-P,Q)-r)+1-i)              &
     &      - CMPLX(SIN(Theta(r-i+1)),0.0E0)
      ENDDO
!
!     Compute norm( U1'*X11*V1 - D11 ) / ( MAX(1,P,Q)*EPS2 ) .
!
      resid = CLANGE('1',P,Q,X,Ldx,Rwork)
      Result(10) = (resid/REAL(MAX(1,P,Q)))/eps2
!
!     Compute norm( U2'*X21*V1 - D21 ) / ( MAX(1,M-P,Q)*EPS2 ) .
!
      resid = CLANGE('1',M-P,Q,X(P+1,1),Ldx,Rwork)
      Result(11) = (resid/REAL(MAX(1,M-P,Q)))/eps2
!
!     Compute I - U1'*U1
!
      CALL CLASET('Full',P,P,ZERO,ONE,Work,Ldu1)
      CALL CHERK('Upper','Conjugate transpose',P,P,-REALONE,U1,Ldu1,    &
     &           REALONE,Work,Ldu1)
!
!     Compute norm( I - U1'*U1 ) / ( MAX(1,P) * ULP ) .
!
      resid = CLANHE('1','Upper',P,Work,Ldu1,Rwork)
      Result(12) = (resid/REAL(MAX(1,P)))/ulp
!
!     Compute I - U2'*U2
!
      CALL CLASET('Full',M-P,M-P,ZERO,ONE,Work,Ldu2)
      CALL CHERK('Upper','Conjugate transpose',M-P,M-P,-REALONE,U2,Ldu2,&
     &           REALONE,Work,Ldu2)
!
!     Compute norm( I - U2'*U2 ) / ( MAX(1,M-P) * ULP ) .
!
      resid = CLANHE('1','Upper',M-P,Work,Ldu2,Rwork)
      Result(13) = (resid/REAL(MAX(1,M-P)))/ulp
!
!     Compute I - V1T*V1T'
!
      CALL CLASET('Full',Q,Q,ZERO,ONE,Work,Ldv1t)
      CALL CHERK('Upper','No transpose',Q,Q,-REALONE,V1t,Ldv1t,REALONE, &
     &           Work,Ldv1t)
!
!     Compute norm( I - V1T*V1T' ) / ( MAX(1,Q) * ULP ) .
!
      resid = CLANHE('1','Upper',Q,Work,Ldv1t,Rwork)
      Result(14) = (resid/REAL(MAX(1,Q)))/ulp
!
!     Check sorting
!
      Result(15) = REALZERO
      DO i = 1 , r
         IF ( Theta(i)<REALZERO .OR. Theta(i)>PIOVER2 ) Result(15)      &
     &        = ulpinv
         IF ( i>1 ) THEN
            IF ( Theta(i)<Theta(i-1) ) Result(15) = ulpinv
         ENDIF
      ENDDO
!
!
!     End of ccsdts
!
      END SUBROUTINE ccsdts
