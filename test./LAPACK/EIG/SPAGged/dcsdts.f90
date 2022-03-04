!*==dcsdts.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b DCSDTS
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DCSDTS( M, P, Q, X, XF, LDX, U1, LDU1, U2, LDU2, V1T,
!                          LDV1T, V2T, LDV2T, THETA, IWORK, WORK, LWORK,
!                          RWORK, RESULT )
!
!       .. Scalar Arguments ..
!       INTEGER            LDX, LDU1, LDU2, LDV1T, LDV2T, LWORK, M, P, Q
!       ..
!       .. Array Arguments ..
!       INTEGER            IWORK( * )
!       DOUBLE PRECISION   RESULT( 15 ), RWORK( * ), THETA( * )
!       DOUBLE PRECISION   U1( LDU1, * ), U2( LDU2, * ), V1T( LDV1T, * ),
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
!> DCSDTS tests DORCSD, which, given an M-by-M partitioned orthogonal
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
!>                           = [---------------------] = [ D11 D12 ] ,
!>                             [  0  0  0 |  I  0  0 ]   [ D21 D22 ]
!>                             [  0  S  0 |  0  C  0 ]
!>                             [  0  0  I |  0  0  0 ]
!>
!> and also DORCSD2BY1, which, given
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
!>          X is DOUBLE PRECISION array, dimension (LDX,M)
!>          The M-by-M matrix X.
!> \endverbatim
!>
!> \param[out] XF
!> \verbatim
!>          XF is DOUBLE PRECISION array, dimension (LDX,M)
!>          Details of the CSD of X, as returned by DORCSD;
!>          see DORCSD for further details.
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
!>          U1 is DOUBLE PRECISION array, dimension(LDU1,P)
!>          The P-by-P orthogonal matrix U1.
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
!>          U2 is DOUBLE PRECISION array, dimension(LDU2,M-P)
!>          The (M-P)-by-(M-P) orthogonal matrix U2.
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
!>          V1T is DOUBLE PRECISION array, dimension(LDV1T,Q)
!>          The Q-by-Q orthogonal matrix V1T.
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
!>          V2T is DOUBLE PRECISION array, dimension(LDV2T,M-Q)
!>          The (M-Q)-by-(M-Q) orthogonal matrix V2T.
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
!>          THETA is DOUBLE PRECISION array, dimension MIN(P,M-P,Q,M-Q)
!>          The CS values of X; the essentially diagonal matrices C and
!>          S are constructed from THETA; see subroutine DORCSD for
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
!>          WORK is DOUBLE PRECISION array, dimension (LWORK)
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
!>          RWORK is DOUBLE PRECISION array
!> \endverbatim
!>
!> \param[out] RESULT
!> \verbatim
!>          RESULT is DOUBLE PRECISION array, dimension (15)
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
!> \ingroup double_eig
!
!  =====================================================================
      SUBROUTINE DCSDTS(M,P,Q,X,Xf,Ldx,U1,Ldu1,U2,Ldu2,V1t,Ldv1t,V2t,   &
     &                  Ldv2t,Theta,Iwork,Work,Lwork,Rwork,Result)
      IMPLICIT NONE
!*--DCSDTS232
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
      DOUBLE PRECISION Result(15) , Rwork(*) , Theta(*)
      DOUBLE PRECISION U1(Ldu1,*) , U2(Ldu2,*) , V1t(Ldv1t,*) ,         &
     &                 V2t(Ldv2t,*) , Work(Lwork) , X(Ldx,*) , Xf(Ldx,*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION REALONE , REALZERO
      PARAMETER (REALONE=1.0D0,REALZERO=0.0D0)
      DOUBLE PRECISION ZERO , ONE
      PARAMETER (ZERO=0.0D0,ONE=1.0D0)
      DOUBLE PRECISION PIOVER2
      PARAMETER (PIOVER2=1.57079632679489661923132169163975144210D0)
!     ..
!     .. Local Scalars ..
      INTEGER i , info , r
      DOUBLE PRECISION eps2 , resid , ulp , ulpinv
!     ..
!     .. External Functions ..
      DOUBLE PRECISION DLAMCH , DLANGE , DLANSY
      EXTERNAL DLAMCH , DLANGE , DLANSY
!     ..
!     .. External Subroutines ..
      EXTERNAL DGEMM , DLACPY , DLASET , DORCSD , DORCSD2BY1 , DSYRK
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC COS , DBLE , MAX , MIN , SIN
!     ..
!     .. Executable Statements ..
!
      ulp = DLAMCH('Precision')
      ulpinv = REALONE/ulp
!
!     The first half of the routine checks the 2-by-2 CSD
!
      CALL DLASET('Full',M,M,ZERO,ONE,Work,Ldx)
      CALL DSYRK('Upper','Conjugate transpose',M,M,-ONE,X,Ldx,ONE,Work, &
     &           Ldx)
      IF ( M>0 ) THEN
         eps2 = MAX(ulp,DLANGE('1',M,M,Work,Ldx,Rwork)/DBLE(M))
      ELSE
         eps2 = ulp
      ENDIF
      r = MIN(P,M-P,Q,M-Q)
!
!     Copy the matrix X to the array XF.
!
      CALL DLACPY('Full',M,M,X,Ldx,Xf,Ldx)
!
!     Compute the CSD
!
      CALL DORCSD('Y','Y','Y','Y','N','D',M,P,Q,Xf(1,1),Ldx,Xf(1,Q+1),  &
     &            Ldx,Xf(P+1,1),Ldx,Xf(P+1,Q+1),Ldx,Theta,U1,Ldu1,U2,   &
     &            Ldu2,V1t,Ldv1t,V2t,Ldv2t,Work,Lwork,Iwork,info)
!
!     Compute XF := diag(U1,U2)'*X*diag(V1,V2) - [D11 D12; D21 D22]
!
      CALL DLACPY('Full',M,M,X,Ldx,Xf,Ldx)
!
      CALL DGEMM('No transpose','Conjugate transpose',P,Q,Q,ONE,Xf,Ldx, &
     &           V1t,Ldv1t,ZERO,Work,Ldx)
!
      CALL DGEMM('Conjugate transpose','No transpose',P,Q,P,ONE,U1,Ldu1,&
     &           Work,Ldx,ZERO,Xf,Ldx)
!
      DO i = 1 , MIN(P,Q) - r
         Xf(i,i) = Xf(i,i) - ONE
      ENDDO
      DO i = 1 , r
         Xf(MIN(P,Q)-r+i,MIN(P,Q)-r+i) = Xf(MIN(P,Q)-r+i,MIN(P,Q)-r+i)  &
     &      - COS(Theta(i))
      ENDDO
!
      CALL DGEMM('No transpose','Conjugate transpose',P,M-Q,M-Q,ONE,    &
     &           Xf(1,Q+1),Ldx,V2t,Ldv2t,ZERO,Work,Ldx)
!
      CALL DGEMM('Conjugate transpose','No transpose',P,M-Q,P,ONE,U1,   &
     &           Ldu1,Work,Ldx,ZERO,Xf(1,Q+1),Ldx)
!
      DO i = 1 , MIN(P,M-Q) - r
         Xf(P-i+1,M-i+1) = Xf(P-i+1,M-i+1) + ONE
      ENDDO
      DO i = 1 , r
         Xf(P-(MIN(P,M-Q)-r)+1-i,M-(MIN(P,M-Q)-r)+1-i)                  &
     &      = Xf(P-(MIN(P,M-Q)-r)+1-i,M-(MIN(P,M-Q)-r)+1-i)             &
     &      + SIN(Theta(r-i+1))
      ENDDO
!
      CALL DGEMM('No transpose','Conjugate transpose',M-P,Q,Q,ONE,      &
     &           Xf(P+1,1),Ldx,V1t,Ldv1t,ZERO,Work,Ldx)
!
      CALL DGEMM('Conjugate transpose','No transpose',M-P,Q,M-P,ONE,U2, &
     &           Ldu2,Work,Ldx,ZERO,Xf(P+1,1),Ldx)
!
      DO i = 1 , MIN(M-P,Q) - r
         Xf(M-i+1,Q-i+1) = Xf(M-i+1,Q-i+1) - ONE
      ENDDO
      DO i = 1 , r
         Xf(M-(MIN(M-P,Q)-r)+1-i,Q-(MIN(M-P,Q)-r)+1-i)                  &
     &      = Xf(M-(MIN(M-P,Q)-r)+1-i,Q-(MIN(M-P,Q)-r)+1-i)             &
     &      - SIN(Theta(r-i+1))
      ENDDO
!
      CALL DGEMM('No transpose','Conjugate transpose',M-P,M-Q,M-Q,ONE,  &
     &           Xf(P+1,Q+1),Ldx,V2t,Ldv2t,ZERO,Work,Ldx)
!
      CALL DGEMM('Conjugate transpose','No transpose',M-P,M-Q,M-P,ONE,  &
     &           U2,Ldu2,Work,Ldx,ZERO,Xf(P+1,Q+1),Ldx)
!
      DO i = 1 , MIN(M-P,M-Q) - r
         Xf(P+i,Q+i) = Xf(P+i,Q+i) - ONE
      ENDDO
      DO i = 1 , r
         Xf(P+(MIN(M-P,M-Q)-r)+i,Q+(MIN(M-P,M-Q)-r)+i)                  &
     &      = Xf(P+(MIN(M-P,M-Q)-r)+i,Q+(MIN(M-P,M-Q)-r)+i)             &
     &      - COS(Theta(i))
      ENDDO
!
!     Compute norm( U1'*X11*V1 - D11 ) / ( MAX(1,P,Q)*EPS2 ) .
!
      resid = DLANGE('1',P,Q,Xf,Ldx,Rwork)
      Result(1) = (resid/DBLE(MAX(1,P,Q)))/eps2
!
!     Compute norm( U1'*X12*V2 - D12 ) / ( MAX(1,P,M-Q)*EPS2 ) .
!
      resid = DLANGE('1',P,M-Q,Xf(1,Q+1),Ldx,Rwork)
      Result(2) = (resid/DBLE(MAX(1,P,M-Q)))/eps2
!
!     Compute norm( U2'*X21*V1 - D21 ) / ( MAX(1,M-P,Q)*EPS2 ) .
!
      resid = DLANGE('1',M-P,Q,Xf(P+1,1),Ldx,Rwork)
      Result(3) = (resid/DBLE(MAX(1,M-P,Q)))/eps2
!
!     Compute norm( U2'*X22*V2 - D22 ) / ( MAX(1,M-P,M-Q)*EPS2 ) .
!
      resid = DLANGE('1',M-P,M-Q,Xf(P+1,Q+1),Ldx,Rwork)
      Result(4) = (resid/DBLE(MAX(1,M-P,M-Q)))/eps2
!
!     Compute I - U1'*U1
!
      CALL DLASET('Full',P,P,ZERO,ONE,Work,Ldu1)
      CALL DSYRK('Upper','Conjugate transpose',P,P,-ONE,U1,Ldu1,ONE,    &
     &           Work,Ldu1)
!
!     Compute norm( I - U'*U ) / ( MAX(1,P) * ULP ) .
!
      resid = DLANSY('1','Upper',P,Work,Ldu1,Rwork)
      Result(5) = (resid/DBLE(MAX(1,P)))/ulp
!
!     Compute I - U2'*U2
!
      CALL DLASET('Full',M-P,M-P,ZERO,ONE,Work,Ldu2)
      CALL DSYRK('Upper','Conjugate transpose',M-P,M-P,-ONE,U2,Ldu2,ONE,&
     &           Work,Ldu2)
!
!     Compute norm( I - U2'*U2 ) / ( MAX(1,M-P) * ULP ) .
!
      resid = DLANSY('1','Upper',M-P,Work,Ldu2,Rwork)
      Result(6) = (resid/DBLE(MAX(1,M-P)))/ulp
!
!     Compute I - V1T*V1T'
!
      CALL DLASET('Full',Q,Q,ZERO,ONE,Work,Ldv1t)
      CALL DSYRK('Upper','No transpose',Q,Q,-ONE,V1t,Ldv1t,ONE,Work,    &
     &           Ldv1t)
!
!     Compute norm( I - V1T*V1T' ) / ( MAX(1,Q) * ULP ) .
!
      resid = DLANSY('1','Upper',Q,Work,Ldv1t,Rwork)
      Result(7) = (resid/DBLE(MAX(1,Q)))/ulp
!
!     Compute I - V2T*V2T'
!
      CALL DLASET('Full',M-Q,M-Q,ZERO,ONE,Work,Ldv2t)
      CALL DSYRK('Upper','No transpose',M-Q,M-Q,-ONE,V2t,Ldv2t,ONE,Work,&
     &           Ldv2t)
!
!     Compute norm( I - V2T*V2T' ) / ( MAX(1,M-Q) * ULP ) .
!
      resid = DLANSY('1','Upper',M-Q,Work,Ldv2t,Rwork)
      Result(8) = (resid/DBLE(MAX(1,M-Q)))/ulp
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
      CALL DLASET('Full',Q,Q,ZERO,ONE,Work,Ldx)
      CALL DSYRK('Upper','Conjugate transpose',Q,M,-ONE,X,Ldx,ONE,Work, &
     &           Ldx)
      IF ( M>0 ) THEN
         eps2 = MAX(ulp,DLANGE('1',Q,Q,Work,Ldx,Rwork)/DBLE(M))
      ELSE
         eps2 = ulp
      ENDIF
      r = MIN(P,M-P,Q,M-Q)
!
!     Copy the matrix [ X11; X21 ] to the array XF.
!
      CALL DLACPY('Full',M,Q,X,Ldx,Xf,Ldx)
!
!     Compute the CSD
!
      CALL DORCSD2BY1('Y','Y','Y',M,P,Q,Xf(1,1),Ldx,Xf(P+1,1),Ldx,Theta,&
     &                U1,Ldu1,U2,Ldu2,V1t,Ldv1t,Work,Lwork,Iwork,info)
!
!     Compute [X11;X21] := diag(U1,U2)'*[X11;X21]*V1 - [D11;D21]
!
      CALL DGEMM('No transpose','Conjugate transpose',P,Q,Q,ONE,X,Ldx,  &
     &           V1t,Ldv1t,ZERO,Work,Ldx)
!
      CALL DGEMM('Conjugate transpose','No transpose',P,Q,P,ONE,U1,Ldu1,&
     &           Work,Ldx,ZERO,X,Ldx)
!
      DO i = 1 , MIN(P,Q) - r
         X(i,i) = X(i,i) - ONE
      ENDDO
      DO i = 1 , r
         X(MIN(P,Q)-r+i,MIN(P,Q)-r+i) = X(MIN(P,Q)-r+i,MIN(P,Q)-r+i)    &
     &                                  - COS(Theta(i))
      ENDDO
!
      CALL DGEMM('No transpose','Conjugate transpose',M-P,Q,Q,ONE,      &
     &           X(P+1,1),Ldx,V1t,Ldv1t,ZERO,Work,Ldx)
!
      CALL DGEMM('Conjugate transpose','No transpose',M-P,Q,M-P,ONE,U2, &
     &           Ldu2,Work,Ldx,ZERO,X(P+1,1),Ldx)
!
      DO i = 1 , MIN(M-P,Q) - r
         X(M-i+1,Q-i+1) = X(M-i+1,Q-i+1) - ONE
      ENDDO
      DO i = 1 , r
         X(M-(MIN(M-P,Q)-r)+1-i,Q-(MIN(M-P,Q)-r)+1-i)                   &
     &      = X(M-(MIN(M-P,Q)-r)+1-i,Q-(MIN(M-P,Q)-r)+1-i)              &
     &      - SIN(Theta(r-i+1))
      ENDDO
!
!     Compute norm( U1'*X11*V1 - D11 ) / ( MAX(1,P,Q)*EPS2 ) .
!
      resid = DLANGE('1',P,Q,X,Ldx,Rwork)
      Result(10) = (resid/DBLE(MAX(1,P,Q)))/eps2
!
!     Compute norm( U2'*X21*V1 - D21 ) / ( MAX(1,M-P,Q)*EPS2 ) .
!
      resid = DLANGE('1',M-P,Q,X(P+1,1),Ldx,Rwork)
      Result(11) = (resid/DBLE(MAX(1,M-P,Q)))/eps2
!
!     Compute I - U1'*U1
!
      CALL DLASET('Full',P,P,ZERO,ONE,Work,Ldu1)
      CALL DSYRK('Upper','Conjugate transpose',P,P,-ONE,U1,Ldu1,ONE,    &
     &           Work,Ldu1)
!
!     Compute norm( I - U1'*U1 ) / ( MAX(1,P) * ULP ) .
!
      resid = DLANSY('1','Upper',P,Work,Ldu1,Rwork)
      Result(12) = (resid/DBLE(MAX(1,P)))/ulp
!
!     Compute I - U2'*U2
!
      CALL DLASET('Full',M-P,M-P,ZERO,ONE,Work,Ldu2)
      CALL DSYRK('Upper','Conjugate transpose',M-P,M-P,-ONE,U2,Ldu2,ONE,&
     &           Work,Ldu2)
!
!     Compute norm( I - U2'*U2 ) / ( MAX(1,M-P) * ULP ) .
!
      resid = DLANSY('1','Upper',M-P,Work,Ldu2,Rwork)
      Result(13) = (resid/DBLE(MAX(1,M-P)))/ulp
!
!     Compute I - V1T*V1T'
!
      CALL DLASET('Full',Q,Q,ZERO,ONE,Work,Ldv1t)
      CALL DSYRK('Upper','No transpose',Q,Q,-ONE,V1t,Ldv1t,ONE,Work,    &
     &           Ldv1t)
!
!     Compute norm( I - V1T*V1T' ) / ( MAX(1,Q) * ULP ) .
!
      resid = DLANSY('1','Upper',Q,Work,Ldv1t,Rwork)
      Result(14) = (resid/DBLE(MAX(1,Q)))/ulp
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
!     End of DCSDTS
!
      END SUBROUTINE DCSDTS
