!*==cbbcsd.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b CBBCSD
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CBBCSD + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cbbcsd.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cbbcsd.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cbbcsd.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CBBCSD( JOBU1, JOBU2, JOBV1T, JOBV2T, TRANS, M, P, Q,
!                          THETA, PHI, U1, LDU1, U2, LDU2, V1T, LDV1T,
!                          V2T, LDV2T, B11D, B11E, B12D, B12E, B21D, B21E,
!                          B22D, B22E, RWORK, LRWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          JOBU1, JOBU2, JOBV1T, JOBV2T, TRANS
!       INTEGER            INFO, LDU1, LDU2, LDV1T, LDV2T, LRWORK, M, P, Q
!       ..
!       .. Array Arguments ..
!       REAL               B11D( * ), B11E( * ), B12D( * ), B12E( * ),
!      $                   B21D( * ), B21E( * ), B22D( * ), B22E( * ),
!      $                   PHI( * ), THETA( * ), RWORK( * )
!       COMPLEX            U1( LDU1, * ), U2( LDU2, * ), V1T( LDV1T, * ),
!      $                   V2T( LDV2T, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CBBCSD computes the CS decomposition of a unitary matrix in
!> bidiagonal-block form,
!>
!>
!>     [ B11 | B12 0  0 ]
!>     [  0  |  0 -I  0 ]
!> X = [----------------]
!>     [ B21 | B22 0  0 ]
!>     [  0  |  0  0  I ]
!>
!>                               [  C | -S  0  0 ]
!>                   [ U1 |    ] [  0 |  0 -I  0 ] [ V1 |    ]**H
!>                 = [---------] [---------------] [---------]   .
!>                   [    | U2 ] [  S |  C  0  0 ] [    | V2 ]
!>                               [  0 |  0  0  I ]
!>
!> X is M-by-M, its top-left block is P-by-Q, and Q must be no larger
!> than P, M-P, or M-Q. (If Q is not the smallest index, then X must be
!> transposed and/or permuted. This can be done in constant time using
!> the TRANS and SIGNS options. See CUNCSD for details.)
!>
!> The bidiagonal matrices B11, B12, B21, and B22 are represented
!> implicitly by angles THETA(1:Q) and PHI(1:Q-1).
!>
!> The unitary matrices U1, U2, V1T, and V2T are input/output.
!> The input matrices are pre- or post-multiplied by the appropriate
!> singular vector matrices.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] JOBU1
!> \verbatim
!>          JOBU1 is CHARACTER
!>          = 'Y':      U1 is updated;
!>          otherwise:  U1 is not updated.
!> \endverbatim
!>
!> \param[in] JOBU2
!> \verbatim
!>          JOBU2 is CHARACTER
!>          = 'Y':      U2 is updated;
!>          otherwise:  U2 is not updated.
!> \endverbatim
!>
!> \param[in] JOBV1T
!> \verbatim
!>          JOBV1T is CHARACTER
!>          = 'Y':      V1T is updated;
!>          otherwise:  V1T is not updated.
!> \endverbatim
!>
!> \param[in] JOBV2T
!> \verbatim
!>          JOBV2T is CHARACTER
!>          = 'Y':      V2T is updated;
!>          otherwise:  V2T is not updated.
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
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows and columns in X, the unitary matrix in
!>          bidiagonal-block form.
!> \endverbatim
!>
!> \param[in] P
!> \verbatim
!>          P is INTEGER
!>          The number of rows in the top-left block of X. 0 <= P <= M.
!> \endverbatim
!>
!> \param[in] Q
!> \verbatim
!>          Q is INTEGER
!>          The number of columns in the top-left block of X.
!>          0 <= Q <= MIN(P,M-P,M-Q).
!> \endverbatim
!>
!> \param[in,out] THETA
!> \verbatim
!>          THETA is REAL array, dimension (Q)
!>          On entry, the angles THETA(1),...,THETA(Q) that, along with
!>          PHI(1), ...,PHI(Q-1), define the matrix in bidiagonal-block
!>          form. On exit, the angles whose cosines and sines define the
!>          diagonal blocks in the CS decomposition.
!> \endverbatim
!>
!> \param[in,out] PHI
!> \verbatim
!>          PHI is REAL array, dimension (Q-1)
!>          The angles PHI(1),...,PHI(Q-1) that, along with THETA(1),...,
!>          THETA(Q), define the matrix in bidiagonal-block form.
!> \endverbatim
!>
!> \param[in,out] U1
!> \verbatim
!>          U1 is COMPLEX array, dimension (LDU1,P)
!>          On entry, a P-by-P matrix. On exit, U1 is postmultiplied
!>          by the left singular vector matrix common to [ B11 ; 0 ] and
!>          [ B12 0 0 ; 0 -I 0 0 ].
!> \endverbatim
!>
!> \param[in] LDU1
!> \verbatim
!>          LDU1 is INTEGER
!>          The leading dimension of the array U1, LDU1 >= MAX(1,P).
!> \endverbatim
!>
!> \param[in,out] U2
!> \verbatim
!>          U2 is COMPLEX array, dimension (LDU2,M-P)
!>          On entry, an (M-P)-by-(M-P) matrix. On exit, U2 is
!>          postmultiplied by the left singular vector matrix common to
!>          [ B21 ; 0 ] and [ B22 0 0 ; 0 0 I ].
!> \endverbatim
!>
!> \param[in] LDU2
!> \verbatim
!>          LDU2 is INTEGER
!>          The leading dimension of the array U2, LDU2 >= MAX(1,M-P).
!> \endverbatim
!>
!> \param[in,out] V1T
!> \verbatim
!>          V1T is COMPLEX array, dimension (LDV1T,Q)
!>          On entry, a Q-by-Q matrix. On exit, V1T is premultiplied
!>          by the conjugate transpose of the right singular vector
!>          matrix common to [ B11 ; 0 ] and [ B21 ; 0 ].
!> \endverbatim
!>
!> \param[in] LDV1T
!> \verbatim
!>          LDV1T is INTEGER
!>          The leading dimension of the array V1T, LDV1T >= MAX(1,Q).
!> \endverbatim
!>
!> \param[in,out] V2T
!> \verbatim
!>          V2T is COMPLEX array, dimension (LDV2T,M-Q)
!>          On entry, an (M-Q)-by-(M-Q) matrix. On exit, V2T is
!>          premultiplied by the conjugate transpose of the right
!>          singular vector matrix common to [ B12 0 0 ; 0 -I 0 ] and
!>          [ B22 0 0 ; 0 0 I ].
!> \endverbatim
!>
!> \param[in] LDV2T
!> \verbatim
!>          LDV2T is INTEGER
!>          The leading dimension of the array V2T, LDV2T >= MAX(1,M-Q).
!> \endverbatim
!>
!> \param[out] B11D
!> \verbatim
!>          B11D is REAL array, dimension (Q)
!>          When CBBCSD converges, B11D contains the cosines of THETA(1),
!>          ..., THETA(Q). If CBBCSD fails to converge, then B11D
!>          contains the diagonal of the partially reduced top-left
!>          block.
!> \endverbatim
!>
!> \param[out] B11E
!> \verbatim
!>          B11E is REAL array, dimension (Q-1)
!>          When CBBCSD converges, B11E contains zeros. If CBBCSD fails
!>          to converge, then B11E contains the superdiagonal of the
!>          partially reduced top-left block.
!> \endverbatim
!>
!> \param[out] B12D
!> \verbatim
!>          B12D is REAL array, dimension (Q)
!>          When CBBCSD converges, B12D contains the negative sines of
!>          THETA(1), ..., THETA(Q). If CBBCSD fails to converge, then
!>          B12D contains the diagonal of the partially reduced top-right
!>          block.
!> \endverbatim
!>
!> \param[out] B12E
!> \verbatim
!>          B12E is REAL array, dimension (Q-1)
!>          When CBBCSD converges, B12E contains zeros. If CBBCSD fails
!>          to converge, then B12E contains the subdiagonal of the
!>          partially reduced top-right block.
!> \endverbatim
!>
!> \param[out] B21D
!> \verbatim
!>          B21D is REAL array, dimension (Q)
!>          When CBBCSD converges, B21D contains the negative sines of
!>          THETA(1), ..., THETA(Q). If CBBCSD fails to converge, then
!>          B21D contains the diagonal of the partially reduced bottom-left
!>          block.
!> \endverbatim
!>
!> \param[out] B21E
!> \verbatim
!>          B21E is REAL array, dimension (Q-1)
!>          When CBBCSD converges, B21E contains zeros. If CBBCSD fails
!>          to converge, then B21E contains the subdiagonal of the
!>          partially reduced bottom-left block.
!> \endverbatim
!>
!> \param[out] B22D
!> \verbatim
!>          B22D is REAL array, dimension (Q)
!>          When CBBCSD converges, B22D contains the negative sines of
!>          THETA(1), ..., THETA(Q). If CBBCSD fails to converge, then
!>          B22D contains the diagonal of the partially reduced bottom-right
!>          block.
!> \endverbatim
!>
!> \param[out] B22E
!> \verbatim
!>          B22E is REAL array, dimension (Q-1)
!>          When CBBCSD converges, B22E contains zeros. If CBBCSD fails
!>          to converge, then B22E contains the subdiagonal of the
!>          partially reduced bottom-right block.
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is REAL array, dimension (MAX(1,LRWORK))
!>          On exit, if INFO = 0, RWORK(1) returns the optimal LRWORK.
!> \endverbatim
!>
!> \param[in] LRWORK
!> \verbatim
!>          LRWORK is INTEGER
!>          The dimension of the array RWORK. LRWORK >= MAX(1,8*Q).
!>
!>          If LRWORK = -1, then a workspace query is assumed; the
!>          routine only calculates the optimal size of the RWORK array,
!>          returns this value as the first entry of the work array, and
!>          no error message related to LRWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit.
!>          < 0:  if INFO = -i, the i-th argument had an illegal value.
!>          > 0:  if CBBCSD did not converge, INFO specifies the number
!>                of nonzero entries in PHI, and B11D, B11E, etc.,
!>                contain the partially reduced matrix.
!> \endverbatim
!
!> \par Internal Parameters:
!  =========================
!>
!> \verbatim
!>  TOLMUL  REAL, default = MAX(10,MIN(100,EPS**(-1/8)))
!>          TOLMUL controls the convergence criterion of the QR loop.
!>          Angles THETA(i), PHI(i) are rounded to 0 or PI/2 when they
!>          are within TOLMUL*EPS of either bound.
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
      SUBROUTINE CBBCSD(Jobu1,Jobu2,Jobv1t,Jobv2t,Trans,M,P,Q,Theta,Phi,&
     &                  U1,Ldu1,U2,Ldu2,V1t,Ldv1t,V2t,Ldv2t,B11d,B11e,  &
     &                  B12d,B12e,B21d,B21e,B22d,B22e,Rwork,Lrwork,Info)
      IMPLICIT NONE
!*--CBBCSD335
!
!  -- LAPACK computational routine (version 3.7.1) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     June 2016
!
!     .. Scalar Arguments ..
      CHARACTER Jobu1 , Jobu2 , Jobv1t , Jobv2t , Trans
      INTEGER Info , Ldu1 , Ldu2 , Ldv1t , Ldv2t , Lrwork , M , P , Q
!     ..
!     .. Array Arguments ..
      REAL B11d(*) , B11e(*) , B12d(*) , B12e(*) , B21d(*) , B21e(*) ,  &
     &     B22d(*) , B22e(*) , Phi(*) , Theta(*) , Rwork(*)
      COMPLEX U1(Ldu1,*) , U2(Ldu2,*) , V1t(Ldv1t,*) , V2t(Ldv2t,*)
!     ..
!
!  ===================================================================
!
!     .. Parameters ..
      INTEGER MAXITR
      PARAMETER (MAXITR=6)
      REAL HUNDRED , MEIGHTH , ONE , TEN , ZERO
      PARAMETER (HUNDRED=100.0E0,MEIGHTH=-0.125E0,ONE=1.0E0,TEN=10.0E0, &
     &           ZERO=0.0E0)
      COMPLEX NEGONECOMPLEX
      PARAMETER (NEGONECOMPLEX=(-1.0E0,0.0E0))
      REAL PIOVER2
      PARAMETER (PIOVER2=1.57079632679489661923132169163975144210E0)
!     ..
!     .. Local Scalars ..
      LOGICAL colmajor , lquery , restart11 , restart12 , restart21 ,   &
     &        restart22 , wantu1 , wantu2 , wantv1t , wantv2t
      INTEGER i , imin , imax , iter , iu1cs , iu1sn , iu2cs , iu2sn ,  &
     &        iv1tcs , iv1tsn , iv2tcs , iv2tsn , j , lrworkmin ,       &
     &        lrworkopt , maxit , mini
      REAL b11bulge , b12bulge , b21bulge , b22bulge , dummy , eps ,    &
     &     mu , nu , r , sigma11 , sigma21 , temp , thetamax ,          &
     &     thetamin , thresh , tol , tolmul , unfl , x1 , x2 , y1 , y2
!
!     .. External Subroutines ..
      EXTERNAL CLASR , CSCAL , CSWAP , SLARTGP , SLARTGS , SLAS2 ,      &
     &         XERBLA
!     ..
!     .. External Functions ..
      REAL SLAMCH
      LOGICAL LSAME
      EXTERNAL LSAME , SLAMCH
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , ATAN2 , COS , MAX , MIN , SIN , SQRT
!     ..
!     .. Executable Statements ..
!
!     Test input arguments
!
      Info = 0
      lquery = Lrwork== - 1
      wantu1 = LSAME(Jobu1,'Y')
      wantu2 = LSAME(Jobu2,'Y')
      wantv1t = LSAME(Jobv1t,'Y')
      wantv2t = LSAME(Jobv2t,'Y')
      colmajor = .NOT.LSAME(Trans,'T')
!
      IF ( M<0 ) THEN
         Info = -6
      ELSEIF ( P<0 .OR. P>M ) THEN
         Info = -7
      ELSEIF ( Q<0 .OR. Q>M ) THEN
         Info = -8
      ELSEIF ( Q>P .OR. Q>M-P .OR. Q>M-Q ) THEN
         Info = -8
      ELSEIF ( wantu1 .AND. Ldu1<P ) THEN
         Info = -12
      ELSEIF ( wantu2 .AND. Ldu2<M-P ) THEN
         Info = -14
      ELSEIF ( wantv1t .AND. Ldv1t<Q ) THEN
         Info = -16
      ELSEIF ( wantv2t .AND. Ldv2t<M-Q ) THEN
         Info = -18
      ENDIF
!
!     Quick return if Q = 0
!
      IF ( Info==0 .AND. Q==0 ) THEN
         lrworkmin = 1
         Rwork(1) = lrworkmin
         RETURN
      ENDIF
!
!     Compute workspace
!
      IF ( Info==0 ) THEN
         iu1cs = 1
         iu1sn = iu1cs + Q
         iu2cs = iu1sn + Q
         iu2sn = iu2cs + Q
         iv1tcs = iu2sn + Q
         iv1tsn = iv1tcs + Q
         iv2tcs = iv1tsn + Q
         iv2tsn = iv2tcs + Q
         lrworkopt = iv2tsn + Q - 1
         lrworkmin = lrworkopt
         Rwork(1) = lrworkopt
         IF ( Lrwork<lrworkmin .AND. .NOT.lquery ) Info = -28
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('CBBCSD',-Info)
         RETURN
      ELSEIF ( lquery ) THEN
         RETURN
      ENDIF
!
!     Get machine constants
!
      eps = SLAMCH('Epsilon')
      unfl = SLAMCH('Safe minimum')
      tolmul = MAX(TEN,MIN(HUNDRED,eps**MEIGHTH))
      tol = tolmul*eps
      thresh = MAX(tol,MAXITR*Q*Q*unfl)
!
!     Test for negligible sines or cosines
!
      DO i = 1 , Q
         IF ( Theta(i)<thresh ) THEN
            Theta(i) = ZERO
         ELSEIF ( Theta(i)>PIOVER2-thresh ) THEN
            Theta(i) = PIOVER2
         ENDIF
      ENDDO
      DO i = 1 , Q - 1
         IF ( Phi(i)<thresh ) THEN
            Phi(i) = ZERO
         ELSEIF ( Phi(i)>PIOVER2-thresh ) THEN
            Phi(i) = PIOVER2
         ENDIF
      ENDDO
!
!     Initial deflation
!
      imax = Q
      DO WHILE ( imax>1 )
         IF ( Phi(imax-1)/=ZERO ) EXIT
         imax = imax - 1
      ENDDO
      imin = imax - 1
      IF ( imin>1 ) THEN
         DO WHILE ( Phi(imin-1)/=ZERO )
            imin = imin - 1
            IF ( imin<=1 ) EXIT
         ENDDO
      ENDIF
!
!     Initialize iteration counter
!
      maxit = MAXITR*Q*Q
      iter = 0
!
!     Begin main iteration loop
!
      DO WHILE ( imax>1 )
!
!        Compute the matrix entries
!
         B11d(imin) = COS(Theta(imin))
         B21d(imin) = -SIN(Theta(imin))
         DO i = imin , imax - 1
            B11e(i) = -SIN(Theta(i))*SIN(Phi(i))
            B11d(i+1) = COS(Theta(i+1))*COS(Phi(i))
            B12d(i) = SIN(Theta(i))*COS(Phi(i))
            B12e(i) = COS(Theta(i+1))*SIN(Phi(i))
            B21e(i) = -COS(Theta(i))*SIN(Phi(i))
            B21d(i+1) = -SIN(Theta(i+1))*COS(Phi(i))
            B22d(i) = COS(Theta(i))*COS(Phi(i))
            B22e(i) = -SIN(Theta(i+1))*SIN(Phi(i))
         ENDDO
         B12d(imax) = SIN(Theta(imax))
         B22d(imax) = COS(Theta(imax))
!
!        Abort if not converging; otherwise, increment ITER
!
         IF ( iter>maxit ) THEN
            Info = 0
            DO i = 1 , Q
               IF ( Phi(i)/=ZERO ) Info = Info + 1
            ENDDO
            RETURN
         ENDIF
!
         iter = iter + imax - imin
!
!        Compute shifts
!
         thetamax = Theta(imin)
         thetamin = Theta(imin)
         DO i = imin + 1 , imax
            IF ( Theta(i)>thetamax ) thetamax = Theta(i)
            IF ( Theta(i)<thetamin ) thetamin = Theta(i)
         ENDDO
!
         IF ( thetamax>PIOVER2-thresh ) THEN
!
!           Zero on diagonals of B11 and B22; induce deflation with a
!           zero shift
!
            mu = ZERO
            nu = ONE
!
         ELSEIF ( thetamin<thresh ) THEN
!
!           Zero on diagonals of B12 and B22; induce deflation with a
!           zero shift
!
            mu = ONE
            nu = ZERO
!
         ELSE
!
!           Compute shifts for B11 and B21 and use the lesser
!
            CALL SLAS2(B11d(imax-1),B11e(imax-1),B11d(imax),sigma11,    &
     &                 dummy)
            CALL SLAS2(B21d(imax-1),B21e(imax-1),B21d(imax),sigma21,    &
     &                 dummy)
!
            IF ( sigma11<=sigma21 ) THEN
               mu = sigma11
               nu = SQRT(ONE-mu**2)
               IF ( mu<thresh ) THEN
                  mu = ZERO
                  nu = ONE
               ENDIF
            ELSE
               nu = sigma21
               mu = SQRT(1.0-nu**2)
               IF ( nu<thresh ) THEN
                  mu = ONE
                  nu = ZERO
               ENDIF
            ENDIF
         ENDIF
!
!        Rotate to produce bulges in B11 and B21
!
         IF ( mu<=nu ) THEN
            CALL SLARTGS(B11d(imin),B11e(imin),mu,Rwork(iv1tcs+imin-1), &
     &                   Rwork(iv1tsn+imin-1))
         ELSE
            CALL SLARTGS(B21d(imin),B21e(imin),nu,Rwork(iv1tcs+imin-1), &
     &                   Rwork(iv1tsn+imin-1))
         ENDIF
!
         temp = Rwork(iv1tcs+imin-1)*B11d(imin) + Rwork(iv1tsn+imin-1)  &
     &          *B11e(imin)
         B11e(imin) = Rwork(iv1tcs+imin-1)*B11e(imin)                   &
     &                - Rwork(iv1tsn+imin-1)*B11d(imin)
         B11d(imin) = temp
         b11bulge = Rwork(iv1tsn+imin-1)*B11d(imin+1)
         B11d(imin+1) = Rwork(iv1tcs+imin-1)*B11d(imin+1)
         temp = Rwork(iv1tcs+imin-1)*B21d(imin) + Rwork(iv1tsn+imin-1)  &
     &          *B21e(imin)
         B21e(imin) = Rwork(iv1tcs+imin-1)*B21e(imin)                   &
     &                - Rwork(iv1tsn+imin-1)*B21d(imin)
         B21d(imin) = temp
         b21bulge = Rwork(iv1tsn+imin-1)*B21d(imin+1)
         B21d(imin+1) = Rwork(iv1tcs+imin-1)*B21d(imin+1)
!
!        Compute THETA(IMIN)
!
         Theta(imin) = ATAN2(SQRT(B21d(imin)**2+b21bulge**2),           &
     &                 SQRT(B11d(imin)**2+b11bulge**2))
!
!        Chase the bulges in B11(IMIN+1,IMIN) and B21(IMIN+1,IMIN)
!
         IF ( B11d(imin)**2+b11bulge**2>thresh**2 ) THEN
            CALL SLARTGP(b11bulge,B11d(imin),Rwork(iu1sn+imin-1),       &
     &                   Rwork(iu1cs+imin-1),r)
         ELSEIF ( mu<=nu ) THEN
            CALL SLARTGS(B11e(imin),B11d(imin+1),mu,Rwork(iu1cs+imin-1),&
     &                   Rwork(iu1sn+imin-1))
         ELSE
            CALL SLARTGS(B12d(imin),B12e(imin),nu,Rwork(iu1cs+imin-1),  &
     &                   Rwork(iu1sn+imin-1))
         ENDIF
         IF ( B21d(imin)**2+b21bulge**2>thresh**2 ) THEN
            CALL SLARTGP(b21bulge,B21d(imin),Rwork(iu2sn+imin-1),       &
     &                   Rwork(iu2cs+imin-1),r)
         ELSEIF ( nu<mu ) THEN
            CALL SLARTGS(B21e(imin),B21d(imin+1),nu,Rwork(iu2cs+imin-1),&
     &                   Rwork(iu2sn+imin-1))
         ELSE
            CALL SLARTGS(B22d(imin),B22e(imin),mu,Rwork(iu2cs+imin-1),  &
     &                   Rwork(iu2sn+imin-1))
         ENDIF
         Rwork(iu2cs+imin-1) = -Rwork(iu2cs+imin-1)
         Rwork(iu2sn+imin-1) = -Rwork(iu2sn+imin-1)
!
         temp = Rwork(iu1cs+imin-1)*B11e(imin) + Rwork(iu1sn+imin-1)    &
     &          *B11d(imin+1)
         B11d(imin+1) = Rwork(iu1cs+imin-1)*B11d(imin+1)                &
     &                  - Rwork(iu1sn+imin-1)*B11e(imin)
         B11e(imin) = temp
         IF ( imax>imin+1 ) THEN
            b11bulge = Rwork(iu1sn+imin-1)*B11e(imin+1)
            B11e(imin+1) = Rwork(iu1cs+imin-1)*B11e(imin+1)
         ENDIF
         temp = Rwork(iu1cs+imin-1)*B12d(imin) + Rwork(iu1sn+imin-1)    &
     &          *B12e(imin)
         B12e(imin) = Rwork(iu1cs+imin-1)*B12e(imin)                    &
     &                - Rwork(iu1sn+imin-1)*B12d(imin)
         B12d(imin) = temp
         b12bulge = Rwork(iu1sn+imin-1)*B12d(imin+1)
         B12d(imin+1) = Rwork(iu1cs+imin-1)*B12d(imin+1)
         temp = Rwork(iu2cs+imin-1)*B21e(imin) + Rwork(iu2sn+imin-1)    &
     &          *B21d(imin+1)
         B21d(imin+1) = Rwork(iu2cs+imin-1)*B21d(imin+1)                &
     &                  - Rwork(iu2sn+imin-1)*B21e(imin)
         B21e(imin) = temp
         IF ( imax>imin+1 ) THEN
            b21bulge = Rwork(iu2sn+imin-1)*B21e(imin+1)
            B21e(imin+1) = Rwork(iu2cs+imin-1)*B21e(imin+1)
         ENDIF
         temp = Rwork(iu2cs+imin-1)*B22d(imin) + Rwork(iu2sn+imin-1)    &
     &          *B22e(imin)
         B22e(imin) = Rwork(iu2cs+imin-1)*B22e(imin)                    &
     &                - Rwork(iu2sn+imin-1)*B22d(imin)
         B22d(imin) = temp
         b22bulge = Rwork(iu2sn+imin-1)*B22d(imin+1)
         B22d(imin+1) = Rwork(iu2cs+imin-1)*B22d(imin+1)
!
!        Inner loop: chase bulges from B11(IMIN,IMIN+2),
!        B12(IMIN,IMIN+1), B21(IMIN,IMIN+2), and B22(IMIN,IMIN+1) to
!        bottom-right
!
         DO i = imin + 1 , imax - 1
!
!           Compute PHI(I-1)
!
            x1 = SIN(Theta(i-1))*B11e(i-1) + COS(Theta(i-1))*B21e(i-1)
            x2 = SIN(Theta(i-1))*b11bulge + COS(Theta(i-1))*b21bulge
            y1 = SIN(Theta(i-1))*B12d(i-1) + COS(Theta(i-1))*B22d(i-1)
            y2 = SIN(Theta(i-1))*b12bulge + COS(Theta(i-1))*b22bulge
!
            Phi(i-1) = ATAN2(SQRT(x1**2+x2**2),SQRT(y1**2+y2**2))
!
!           Determine if there are bulges to chase or if a new direct
!           summand has been reached
!
            restart11 = B11e(i-1)**2 + b11bulge**2<=thresh**2
            restart21 = B21e(i-1)**2 + b21bulge**2<=thresh**2
            restart12 = B12d(i-1)**2 + b12bulge**2<=thresh**2
            restart22 = B22d(i-1)**2 + b22bulge**2<=thresh**2
!
!           If possible, chase bulges from B11(I-1,I+1), B12(I-1,I),
!           B21(I-1,I+1), and B22(I-1,I). If necessary, restart bulge-
!           chasing by applying the original shift again.
!
            IF ( .NOT.restart11 .AND. .NOT.restart21 ) THEN
               CALL SLARTGP(x2,x1,Rwork(iv1tsn+i-1),Rwork(iv1tcs+i-1),r)
            ELSEIF ( .NOT.restart11 .AND. restart21 ) THEN
               CALL SLARTGP(b11bulge,B11e(i-1),Rwork(iv1tsn+i-1),       &
     &                      Rwork(iv1tcs+i-1),r)
            ELSEIF ( restart11 .AND. .NOT.restart21 ) THEN
               CALL SLARTGP(b21bulge,B21e(i-1),Rwork(iv1tsn+i-1),       &
     &                      Rwork(iv1tcs+i-1),r)
            ELSEIF ( mu<=nu ) THEN
               CALL SLARTGS(B11d(i),B11e(i),mu,Rwork(iv1tcs+i-1),       &
     &                      Rwork(iv1tsn+i-1))
            ELSE
               CALL SLARTGS(B21d(i),B21e(i),nu,Rwork(iv1tcs+i-1),       &
     &                      Rwork(iv1tsn+i-1))
            ENDIF
            Rwork(iv1tcs+i-1) = -Rwork(iv1tcs+i-1)
            Rwork(iv1tsn+i-1) = -Rwork(iv1tsn+i-1)
            IF ( .NOT.restart12 .AND. .NOT.restart22 ) THEN
               CALL SLARTGP(y2,y1,Rwork(iv2tsn+i-1-1),                  &
     &                      Rwork(iv2tcs+i-1-1),r)
            ELSEIF ( .NOT.restart12 .AND. restart22 ) THEN
               CALL SLARTGP(b12bulge,B12d(i-1),Rwork(iv2tsn+i-1-1),     &
     &                      Rwork(iv2tcs+i-1-1),r)
            ELSEIF ( restart12 .AND. .NOT.restart22 ) THEN
               CALL SLARTGP(b22bulge,B22d(i-1),Rwork(iv2tsn+i-1-1),     &
     &                      Rwork(iv2tcs+i-1-1),r)
            ELSEIF ( nu<mu ) THEN
               CALL SLARTGS(B12e(i-1),B12d(i),nu,Rwork(iv2tcs+i-1-1),   &
     &                      Rwork(iv2tsn+i-1-1))
            ELSE
               CALL SLARTGS(B22e(i-1),B22d(i),mu,Rwork(iv2tcs+i-1-1),   &
     &                      Rwork(iv2tsn+i-1-1))
            ENDIF
!
            temp = Rwork(iv1tcs+i-1)*B11d(i) + Rwork(iv1tsn+i-1)*B11e(i)
            B11e(i) = Rwork(iv1tcs+i-1)*B11e(i) - Rwork(iv1tsn+i-1)     &
     &                *B11d(i)
            B11d(i) = temp
            b11bulge = Rwork(iv1tsn+i-1)*B11d(i+1)
            B11d(i+1) = Rwork(iv1tcs+i-1)*B11d(i+1)
            temp = Rwork(iv1tcs+i-1)*B21d(i) + Rwork(iv1tsn+i-1)*B21e(i)
            B21e(i) = Rwork(iv1tcs+i-1)*B21e(i) - Rwork(iv1tsn+i-1)     &
     &                *B21d(i)
            B21d(i) = temp
            b21bulge = Rwork(iv1tsn+i-1)*B21d(i+1)
            B21d(i+1) = Rwork(iv1tcs+i-1)*B21d(i+1)
            temp = Rwork(iv2tcs+i-1-1)*B12e(i-1) + Rwork(iv2tsn+i-1-1)  &
     &             *B12d(i)
            B12d(i) = Rwork(iv2tcs+i-1-1)*B12d(i) - Rwork(iv2tsn+i-1-1) &
     &                *B12e(i-1)
            B12e(i-1) = temp
            b12bulge = Rwork(iv2tsn+i-1-1)*B12e(i)
            B12e(i) = Rwork(iv2tcs+i-1-1)*B12e(i)
            temp = Rwork(iv2tcs+i-1-1)*B22e(i-1) + Rwork(iv2tsn+i-1-1)  &
     &             *B22d(i)
            B22d(i) = Rwork(iv2tcs+i-1-1)*B22d(i) - Rwork(iv2tsn+i-1-1) &
     &                *B22e(i-1)
            B22e(i-1) = temp
            b22bulge = Rwork(iv2tsn+i-1-1)*B22e(i)
            B22e(i) = Rwork(iv2tcs+i-1-1)*B22e(i)
!
!           Compute THETA(I)
!
            x1 = COS(Phi(i-1))*B11d(i) + SIN(Phi(i-1))*B12e(i-1)
            x2 = COS(Phi(i-1))*b11bulge + SIN(Phi(i-1))*b12bulge
            y1 = COS(Phi(i-1))*B21d(i) + SIN(Phi(i-1))*B22e(i-1)
            y2 = COS(Phi(i-1))*b21bulge + SIN(Phi(i-1))*b22bulge
!
            Theta(i) = ATAN2(SQRT(y1**2+y2**2),SQRT(x1**2+x2**2))
!
!           Determine if there are bulges to chase or if a new direct
!           summand has been reached
!
            restart11 = B11d(i)**2 + b11bulge**2<=thresh**2
            restart12 = B12e(i-1)**2 + b12bulge**2<=thresh**2
            restart21 = B21d(i)**2 + b21bulge**2<=thresh**2
            restart22 = B22e(i-1)**2 + b22bulge**2<=thresh**2
!
!           If possible, chase bulges from B11(I+1,I), B12(I+1,I-1),
!           B21(I+1,I), and B22(I+1,I-1). If necessary, restart bulge-
!           chasing by applying the original shift again.
!
            IF ( .NOT.restart11 .AND. .NOT.restart12 ) THEN
               CALL SLARTGP(x2,x1,Rwork(iu1sn+i-1),Rwork(iu1cs+i-1),r)
            ELSEIF ( .NOT.restart11 .AND. restart12 ) THEN
               CALL SLARTGP(b11bulge,B11d(i),Rwork(iu1sn+i-1),          &
     &                      Rwork(iu1cs+i-1),r)
            ELSEIF ( restart11 .AND. .NOT.restart12 ) THEN
               CALL SLARTGP(b12bulge,B12e(i-1),Rwork(iu1sn+i-1),        &
     &                      Rwork(iu1cs+i-1),r)
            ELSEIF ( mu<=nu ) THEN
               CALL SLARTGS(B11e(i),B11d(i+1),mu,Rwork(iu1cs+i-1),      &
     &                      Rwork(iu1sn+i-1))
            ELSE
               CALL SLARTGS(B12d(i),B12e(i),nu,Rwork(iu1cs+i-1),        &
     &                      Rwork(iu1sn+i-1))
            ENDIF
            IF ( .NOT.restart21 .AND. .NOT.restart22 ) THEN
               CALL SLARTGP(y2,y1,Rwork(iu2sn+i-1),Rwork(iu2cs+i-1),r)
            ELSEIF ( .NOT.restart21 .AND. restart22 ) THEN
               CALL SLARTGP(b21bulge,B21d(i),Rwork(iu2sn+i-1),          &
     &                      Rwork(iu2cs+i-1),r)
            ELSEIF ( restart21 .AND. .NOT.restart22 ) THEN
               CALL SLARTGP(b22bulge,B22e(i-1),Rwork(iu2sn+i-1),        &
     &                      Rwork(iu2cs+i-1),r)
            ELSEIF ( nu<mu ) THEN
               CALL SLARTGS(B21e(i),B21e(i+1),nu,Rwork(iu2cs+i-1),      &
     &                      Rwork(iu2sn+i-1))
            ELSE
               CALL SLARTGS(B22d(i),B22e(i),mu,Rwork(iu2cs+i-1),        &
     &                      Rwork(iu2sn+i-1))
            ENDIF
            Rwork(iu2cs+i-1) = -Rwork(iu2cs+i-1)
            Rwork(iu2sn+i-1) = -Rwork(iu2sn+i-1)
!
            temp = Rwork(iu1cs+i-1)*B11e(i) + Rwork(iu1sn+i-1)*B11d(i+1)
            B11d(i+1) = Rwork(iu1cs+i-1)*B11d(i+1) - Rwork(iu1sn+i-1)   &
     &                  *B11e(i)
            B11e(i) = temp
            IF ( i<imax-1 ) THEN
               b11bulge = Rwork(iu1sn+i-1)*B11e(i+1)
               B11e(i+1) = Rwork(iu1cs+i-1)*B11e(i+1)
            ENDIF
            temp = Rwork(iu2cs+i-1)*B21e(i) + Rwork(iu2sn+i-1)*B21d(i+1)
            B21d(i+1) = Rwork(iu2cs+i-1)*B21d(i+1) - Rwork(iu2sn+i-1)   &
     &                  *B21e(i)
            B21e(i) = temp
            IF ( i<imax-1 ) THEN
               b21bulge = Rwork(iu2sn+i-1)*B21e(i+1)
               B21e(i+1) = Rwork(iu2cs+i-1)*B21e(i+1)
            ENDIF
            temp = Rwork(iu1cs+i-1)*B12d(i) + Rwork(iu1sn+i-1)*B12e(i)
            B12e(i) = Rwork(iu1cs+i-1)*B12e(i) - Rwork(iu1sn+i-1)       &
     &                *B12d(i)
            B12d(i) = temp
            b12bulge = Rwork(iu1sn+i-1)*B12d(i+1)
            B12d(i+1) = Rwork(iu1cs+i-1)*B12d(i+1)
            temp = Rwork(iu2cs+i-1)*B22d(i) + Rwork(iu2sn+i-1)*B22e(i)
            B22e(i) = Rwork(iu2cs+i-1)*B22e(i) - Rwork(iu2sn+i-1)       &
     &                *B22d(i)
            B22d(i) = temp
            b22bulge = Rwork(iu2sn+i-1)*B22d(i+1)
            B22d(i+1) = Rwork(iu2cs+i-1)*B22d(i+1)
!
         ENDDO
!
!        Compute PHI(IMAX-1)
!
         x1 = SIN(Theta(imax-1))*B11e(imax-1) + COS(Theta(imax-1))      &
     &        *B21e(imax-1)
         y1 = SIN(Theta(imax-1))*B12d(imax-1) + COS(Theta(imax-1))      &
     &        *B22d(imax-1)
         y2 = SIN(Theta(imax-1))*b12bulge + COS(Theta(imax-1))*b22bulge
!
         Phi(imax-1) = ATAN2(ABS(x1),SQRT(y1**2+y2**2))
!
!        Chase bulges from B12(IMAX-1,IMAX) and B22(IMAX-1,IMAX)
!
         restart12 = B12d(imax-1)**2 + b12bulge**2<=thresh**2
         restart22 = B22d(imax-1)**2 + b22bulge**2<=thresh**2
!
         IF ( .NOT.restart12 .AND. .NOT.restart22 ) THEN
            CALL SLARTGP(y2,y1,Rwork(iv2tsn+imax-1-1),                  &
     &                   Rwork(iv2tcs+imax-1-1),r)
         ELSEIF ( .NOT.restart12 .AND. restart22 ) THEN
            CALL SLARTGP(b12bulge,B12d(imax-1),Rwork(iv2tsn+imax-1-1),  &
     &                   Rwork(iv2tcs+imax-1-1),r)
         ELSEIF ( restart12 .AND. .NOT.restart22 ) THEN
            CALL SLARTGP(b22bulge,B22d(imax-1),Rwork(iv2tsn+imax-1-1),  &
     &                   Rwork(iv2tcs+imax-1-1),r)
         ELSEIF ( nu<mu ) THEN
            CALL SLARTGS(B12e(imax-1),B12d(imax),nu,                    &
     &                   Rwork(iv2tcs+imax-1-1),Rwork(iv2tsn+imax-1-1))
         ELSE
            CALL SLARTGS(B22e(imax-1),B22d(imax),mu,                    &
     &                   Rwork(iv2tcs+imax-1-1),Rwork(iv2tsn+imax-1-1))
         ENDIF
!
         temp = Rwork(iv2tcs+imax-1-1)*B12e(imax-1)                     &
     &          + Rwork(iv2tsn+imax-1-1)*B12d(imax)
         B12d(imax) = Rwork(iv2tcs+imax-1-1)*B12d(imax)                 &
     &                - Rwork(iv2tsn+imax-1-1)*B12e(imax-1)
         B12e(imax-1) = temp
         temp = Rwork(iv2tcs+imax-1-1)*B22e(imax-1)                     &
     &          + Rwork(iv2tsn+imax-1-1)*B22d(imax)
         B22d(imax) = Rwork(iv2tcs+imax-1-1)*B22d(imax)                 &
     &                - Rwork(iv2tsn+imax-1-1)*B22e(imax-1)
         B22e(imax-1) = temp
!
!        Update singular vectors
!
         IF ( wantu1 ) THEN
            IF ( colmajor ) THEN
               CALL CLASR('R','V','F',P,imax-imin+1,Rwork(iu1cs+imin-1),&
     &                    Rwork(iu1sn+imin-1),U1(1,imin),Ldu1)
            ELSE
               CALL CLASR('L','V','F',imax-imin+1,P,Rwork(iu1cs+imin-1),&
     &                    Rwork(iu1sn+imin-1),U1(imin,1),Ldu1)
            ENDIF
         ENDIF
         IF ( wantu2 ) THEN
            IF ( colmajor ) THEN
               CALL CLASR('R','V','F',M-P,imax-imin+1,                  &
     &                    Rwork(iu2cs+imin-1),Rwork(iu2sn+imin-1),      &
     &                    U2(1,imin),Ldu2)
            ELSE
               CALL CLASR('L','V','F',imax-imin+1,M-P,                  &
     &                    Rwork(iu2cs+imin-1),Rwork(iu2sn+imin-1),      &
     &                    U2(imin,1),Ldu2)
            ENDIF
         ENDIF
         IF ( wantv1t ) THEN
            IF ( colmajor ) THEN
               CALL CLASR('L','V','F',imax-imin+1,Q,Rwork(iv1tcs+imin-1)&
     &                    ,Rwork(iv1tsn+imin-1),V1t(imin,1),Ldv1t)
            ELSE
               CALL CLASR('R','V','F',Q,imax-imin+1,Rwork(iv1tcs+imin-1)&
     &                    ,Rwork(iv1tsn+imin-1),V1t(1,imin),Ldv1t)
            ENDIF
         ENDIF
         IF ( wantv2t ) THEN
            IF ( colmajor ) THEN
               CALL CLASR('L','V','F',imax-imin+1,M-Q,                  &
     &                    Rwork(iv2tcs+imin-1),Rwork(iv2tsn+imin-1),    &
     &                    V2t(imin,1),Ldv2t)
            ELSE
               CALL CLASR('R','V','F',M-Q,imax-imin+1,                  &
     &                    Rwork(iv2tcs+imin-1),Rwork(iv2tsn+imin-1),    &
     &                    V2t(1,imin),Ldv2t)
            ENDIF
         ENDIF
!
!        Fix signs on B11(IMAX-1,IMAX) and B21(IMAX-1,IMAX)
!
         IF ( B11e(imax-1)+B21e(imax-1)>0 ) THEN
            B11d(imax) = -B11d(imax)
            B21d(imax) = -B21d(imax)
            IF ( wantv1t ) THEN
               IF ( colmajor ) THEN
                  CALL CSCAL(Q,NEGONECOMPLEX,V1t(imax,1),Ldv1t)
               ELSE
                  CALL CSCAL(Q,NEGONECOMPLEX,V1t(1,imax),1)
               ENDIF
            ENDIF
         ENDIF
!
!        Compute THETA(IMAX)
!
         x1 = COS(Phi(imax-1))*B11d(imax) + SIN(Phi(imax-1))            &
     &        *B12e(imax-1)
         y1 = COS(Phi(imax-1))*B21d(imax) + SIN(Phi(imax-1))            &
     &        *B22e(imax-1)
!
         Theta(imax) = ATAN2(ABS(y1),ABS(x1))
!
!        Fix signs on B11(IMAX,IMAX), B12(IMAX,IMAX-1), B21(IMAX,IMAX),
!        and B22(IMAX,IMAX-1)
!
         IF ( B11d(imax)+B12e(imax-1)<0 ) THEN
            B12d(imax) = -B12d(imax)
            IF ( wantu1 ) THEN
               IF ( colmajor ) THEN
                  CALL CSCAL(P,NEGONECOMPLEX,U1(1,imax),1)
               ELSE
                  CALL CSCAL(P,NEGONECOMPLEX,U1(imax,1),Ldu1)
               ENDIF
            ENDIF
         ENDIF
         IF ( B21d(imax)+B22e(imax-1)>0 ) THEN
            B22d(imax) = -B22d(imax)
            IF ( wantu2 ) THEN
               IF ( colmajor ) THEN
                  CALL CSCAL(M-P,NEGONECOMPLEX,U2(1,imax),1)
               ELSE
                  CALL CSCAL(M-P,NEGONECOMPLEX,U2(imax,1),Ldu2)
               ENDIF
            ENDIF
         ENDIF
!
!        Fix signs on B12(IMAX,IMAX) and B22(IMAX,IMAX)
!
         IF ( B12d(imax)+B22d(imax)<0 ) THEN
            IF ( wantv2t ) THEN
               IF ( colmajor ) THEN
                  CALL CSCAL(M-Q,NEGONECOMPLEX,V2t(imax,1),Ldv2t)
               ELSE
                  CALL CSCAL(M-Q,NEGONECOMPLEX,V2t(1,imax),1)
               ENDIF
            ENDIF
         ENDIF
!
!        Test for negligible sines or cosines
!
         DO i = imin , imax
            IF ( Theta(i)<thresh ) THEN
               Theta(i) = ZERO
            ELSEIF ( Theta(i)>PIOVER2-thresh ) THEN
               Theta(i) = PIOVER2
            ENDIF
         ENDDO
         DO i = imin , imax - 1
            IF ( Phi(i)<thresh ) THEN
               Phi(i) = ZERO
            ELSEIF ( Phi(i)>PIOVER2-thresh ) THEN
               Phi(i) = PIOVER2
            ENDIF
         ENDDO
!
!        Deflate
!
         IF ( imax>1 ) THEN
            DO WHILE ( Phi(imax-1)==ZERO )
               imax = imax - 1
               IF ( imax<=1 ) EXIT
            ENDDO
         ENDIF
         IF ( imin>imax-1 ) imin = imax - 1
         IF ( imin>1 ) THEN
            DO WHILE ( Phi(imin-1)/=ZERO )
               imin = imin - 1
               IF ( imin<=1 ) EXIT
            ENDDO
         ENDIF
!
!        Repeat main iteration loop
!
      ENDDO
!
!     Postprocessing: order THETA from least to greatest
!
      DO i = 1 , Q
!
         mini = i
         thetamin = Theta(i)
         DO j = i + 1 , Q
            IF ( Theta(j)<thetamin ) THEN
               mini = j
               thetamin = Theta(j)
            ENDIF
         ENDDO
!
         IF ( mini/=i ) THEN
            Theta(mini) = Theta(i)
            Theta(i) = thetamin
            IF ( colmajor ) THEN
               IF ( wantu1 ) CALL CSWAP(P,U1(1,i),1,U1(1,mini),1)
               IF ( wantu2 ) CALL CSWAP(M-P,U2(1,i),1,U2(1,mini),1)
               IF ( wantv1t ) CALL CSWAP(Q,V1t(i,1),Ldv1t,V1t(mini,1),  &
     &              Ldv1t)
               IF ( wantv2t ) CALL CSWAP(M-Q,V2t(i,1),Ldv2t,V2t(mini,1),&
     &              Ldv2t)
            ELSE
               IF ( wantu1 ) CALL CSWAP(P,U1(i,1),Ldu1,U1(mini,1),Ldu1)
               IF ( wantu2 ) CALL CSWAP(M-P,U2(i,1),Ldu2,U2(mini,1),    &
     &                                  Ldu2)
               IF ( wantv1t ) CALL CSWAP(Q,V1t(1,i),1,V1t(1,mini),1)
               IF ( wantv2t ) CALL CSWAP(M-Q,V2t(1,i),1,V2t(1,mini),1)
            ENDIF
         ENDIF
!
      ENDDO
!
!
!     End of CBBCSD
!
      END SUBROUTINE CBBCSD
