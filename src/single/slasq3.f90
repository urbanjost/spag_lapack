!*==slasq3.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b SLASQ3 checks for deflation, computes a shift and calls dqds. Used by sbdsqr.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SLASQ3 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slasq3.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slasq3.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slasq3.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SLASQ3( I0, N0, Z, PP, DMIN, SIGMA, DESIG, QMAX, NFAIL,
!                          ITER, NDIV, IEEE, TTYPE, DMIN1, DMIN2, DN, DN1,
!                          DN2, G, TAU )
!
!       .. Scalar Arguments ..
!       LOGICAL            IEEE
!       INTEGER            I0, ITER, N0, NDIV, NFAIL, PP
!       REAL               DESIG, DMIN, DMIN1, DMIN2, DN, DN1, DN2, G,
!      $                   QMAX, SIGMA, TAU
!       ..
!       .. Array Arguments ..
!       REAL               Z( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SLASQ3 checks for deflation, computes a shift (TAU) and calls dqds.
!> In case of failure it changes shifts, and tries again until output
!> is positive.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] I0
!> \verbatim
!>          I0 is INTEGER
!>         First index.
!> \endverbatim
!>
!> \param[in,out] N0
!> \verbatim
!>          N0 is INTEGER
!>         Last index.
!> \endverbatim
!>
!> \param[in,out] Z
!> \verbatim
!>          Z is REAL array, dimension ( 4*N0 )
!>         Z holds the qd array.
!> \endverbatim
!>
!> \param[in,out] PP
!> \verbatim
!>          PP is INTEGER
!>         PP=0 for ping, PP=1 for pong.
!>         PP=2 indicates that flipping was applied to the Z array
!>         and that the initial tests for deflation should not be
!>         performed.
!> \endverbatim
!>
!> \param[out] DMIN
!> \verbatim
!>          DMIN is REAL
!>         Minimum value of d.
!> \endverbatim
!>
!> \param[out] SIGMA
!> \verbatim
!>          SIGMA is REAL
!>         Sum of shifts used in current segment.
!> \endverbatim
!>
!> \param[in,out] DESIG
!> \verbatim
!>          DESIG is REAL
!>         Lower order part of SIGMA
!> \endverbatim
!>
!> \param[in] QMAX
!> \verbatim
!>          QMAX is REAL
!>         Maximum value of q.
!> \endverbatim
!>
!> \param[in,out] NFAIL
!> \verbatim
!>          NFAIL is INTEGER
!>         Increment NFAIL by 1 each time the shift was too big.
!> \endverbatim
!>
!> \param[in,out] ITER
!> \verbatim
!>          ITER is INTEGER
!>         Increment ITER by 1 for each iteration.
!> \endverbatim
!>
!> \param[in,out] NDIV
!> \verbatim
!>          NDIV is INTEGER
!>         Increment NDIV by 1 for each division.
!> \endverbatim
!>
!> \param[in] IEEE
!> \verbatim
!>          IEEE is LOGICAL
!>         Flag for IEEE or non IEEE arithmetic (passed to SLASQ5).
!> \endverbatim
!>
!> \param[in,out] TTYPE
!> \verbatim
!>          TTYPE is INTEGER
!>         Shift type.
!> \endverbatim
!>
!> \param[in,out] DMIN1
!> \verbatim
!>          DMIN1 is REAL
!> \endverbatim
!>
!> \param[in,out] DMIN2
!> \verbatim
!>          DMIN2 is REAL
!> \endverbatim
!>
!> \param[in,out] DN
!> \verbatim
!>          DN is REAL
!> \endverbatim
!>
!> \param[in,out] DN1
!> \verbatim
!>          DN1 is REAL
!> \endverbatim
!>
!> \param[in,out] DN2
!> \verbatim
!>          DN2 is REAL
!> \endverbatim
!>
!> \param[in,out] G
!> \verbatim
!>          G is REAL
!> \endverbatim
!>
!> \param[in,out] TAU
!> \verbatim
!>          TAU is REAL
!>
!>         These are passed as arguments in order to save their values
!>         between calls to SLASQ3.
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
!> \ingroup auxOTHERcomputational
!
!  =====================================================================
      SUBROUTINE SLASQ3(I0,N0,Z,Pp,Dmin,Sigma,Desig,Qmax,Nfail,Iter,    &
     &                  Ndiv,Ieee,Ttype,Dmin1,Dmin2,Dn,Dn1,Dn2,G,Tau)
      IMPLICIT NONE
!*--SLASQ3185
!
!  -- LAPACK computational routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     June 2016
!
!     .. Scalar Arguments ..
      LOGICAL Ieee
      INTEGER I0 , Iter , N0 , Ndiv , Nfail , Pp
      REAL Desig , Dmin , Dmin1 , Dmin2 , Dn , Dn1 , Dn2 , G , Qmax ,   &
     &     Sigma , Tau
!     ..
!     .. Array Arguments ..
      REAL Z(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL CBIAS
      PARAMETER (CBIAS=1.50E0)
      REAL ZERO , QURTR , HALF , ONE , TWO , HUNDRD
      PARAMETER (ZERO=0.0E0,QURTR=0.250E0,HALF=0.5E0,ONE=1.0E0,         &
     &           TWO=2.0E0,HUNDRD=100.0E0)
!     ..
!     .. Local Scalars ..
      INTEGER ipn4 , j4 , n0in , nn , Ttype
      REAL eps , s , t , temp , tol , tol2
!     ..
!     .. External Subroutines ..
      EXTERNAL SLASQ4 , SLASQ5 , SLASQ6
!     ..
!     .. External Function ..
      REAL SLAMCH
      LOGICAL SISNAN
      EXTERNAL SISNAN , SLAMCH
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , MAX , MIN , SQRT
!     ..
!     .. Executable Statements ..
!
      n0in = N0
      eps = SLAMCH('Precision')
      tol = eps*HUNDRD
      tol2 = tol**2
!
!     Check for deflation.
!
!
 100  IF ( N0<I0 ) RETURN
      IF ( N0/=I0 ) THEN
         nn = 4*N0 + Pp
         IF ( N0==(I0+1) ) GOTO 200
!
!     Check whether E(N0-1) is negligible, 1 eigenvalue.
!
         IF ( Z(nn-5)>tol2*(Sigma+Z(nn-3)) .AND. Z(nn-2*Pp-4)           &
     &        >tol2*Z(nn-7) ) THEN
!
!     Check  whether E(N0-2) is negligible, 2 eigenvalues.
!
!
            IF ( Z(nn-9)<=tol2*Sigma .OR. Z(nn-2*Pp-8)<=tol2*Z(nn-11) ) &
     &           GOTO 200
!
            IF ( Pp==2 ) Pp = 0
!
!     Reverse the qd-array, if warranted.
!
            IF ( Dmin<=ZERO .OR. N0<n0in ) THEN
               IF ( CBIAS*Z(4*I0+Pp-3)<Z(4*N0+Pp-3) ) THEN
                  ipn4 = 4*(I0+N0)
                  DO j4 = 4*I0 , 2*(I0+N0-1) , 4
                     temp = Z(j4-3)
                     Z(j4-3) = Z(ipn4-j4-3)
                     Z(ipn4-j4-3) = temp
                     temp = Z(j4-2)
                     Z(j4-2) = Z(ipn4-j4-2)
                     Z(ipn4-j4-2) = temp
                     temp = Z(j4-1)
                     Z(j4-1) = Z(ipn4-j4-5)
                     Z(ipn4-j4-5) = temp
                     temp = Z(j4)
                     Z(j4) = Z(ipn4-j4-4)
                     Z(ipn4-j4-4) = temp
                  ENDDO
                  IF ( N0-I0<=4 ) THEN
                     Z(4*N0+Pp-1) = Z(4*I0+Pp-1)
                     Z(4*N0-Pp) = Z(4*I0-Pp)
                  ENDIF
                  Dmin2 = MIN(Dmin2,Z(4*N0+Pp-1))
                  Z(4*N0+Pp-1) = MIN(Z(4*N0+Pp-1),Z(4*I0+Pp-1),Z(4*I0+Pp&
     &                           +3))
                  Z(4*N0-Pp) = MIN(Z(4*N0-Pp),Z(4*I0-Pp),Z(4*I0-Pp+4))
                  Qmax = MAX(Qmax,Z(4*I0+Pp-3),Z(4*I0+Pp+1))
                  Dmin = -ZERO
               ENDIF
            ENDIF
!
!     Choose a shift.
!
            CALL SLASQ4(I0,N0,Z,Pp,n0in,Dmin,Dmin1,Dmin2,Dn,Dn1,Dn2,Tau,&
     &                  Ttype,G)
            DO
!
!     Call dqds until DMIN > 0.
!
!
               CALL SLASQ5(I0,N0,Z,Pp,Tau,Sigma,Dmin,Dmin1,Dmin2,Dn,Dn1,&
     &                     Dn2,Ieee,eps)
!
               Ndiv = Ndiv + (N0-I0+2)
               Iter = Iter + 1
!
!     Check status.
!
!
!        Success.
!
               IF ( Dmin>=ZERO .AND. Dmin1>=ZERO ) GOTO 400
!
               IF ( Dmin<ZERO .AND. Dmin1>ZERO .AND. Z(4*(N0-1)-Pp)     &
     &              <tol*(Sigma+Dn1) .AND. ABS(Dn)<tol*Sigma ) THEN
!
!        Convergence hidden by negative DN.
!
                  Z(4*(N0-1)-Pp+2) = ZERO
                  Dmin = ZERO
                  GOTO 400
               ELSEIF ( Dmin<ZERO ) THEN
!
!        TAU too big. Select new TAU and try again.
!
                  Nfail = Nfail + 1
                  IF ( Ttype<-22 ) THEN
!
!           Failed twice. Play it safe.
!
                     Tau = ZERO
                  ELSEIF ( Dmin1>ZERO ) THEN
!
!           Late failure. Gives excellent shift.
!
                     Tau = (Tau+Dmin)*(ONE-TWO*eps)
                     Ttype = Ttype - 11
                  ELSE
!
!           Early failure. Divide by 4.
!
                     Tau = QURTR*Tau
                     Ttype = Ttype - 12
                  ENDIF
               ELSEIF ( SISNAN(Dmin) ) THEN
!
!        NaN.
!
                  IF ( Tau==ZERO ) GOTO 300
                  Tau = ZERO
               ELSE
!
!        Possible underflow. Play it safe.
!
                  GOTO 300
               ENDIF
            ENDDO
         ENDIF
      ENDIF
!
!
      Z(4*N0-3) = Z(4*N0+Pp-3) + Sigma
      N0 = N0 - 1
      GOTO 100
!
!
 200  IF ( Z(nn-3)>Z(nn-7) ) THEN
         s = Z(nn-3)
         Z(nn-3) = Z(nn-7)
         Z(nn-7) = s
      ENDIF
      t = HALF*((Z(nn-7)-Z(nn-3))+Z(nn-5))
      IF ( Z(nn-5)>Z(nn-3)*tol2 .AND. t/=ZERO ) THEN
         s = Z(nn-3)*(Z(nn-5)/t)
         IF ( s<=t ) THEN
            s = Z(nn-3)*(Z(nn-5)/(t*(ONE+SQRT(ONE+s/t))))
         ELSE
            s = Z(nn-3)*(Z(nn-5)/(t+SQRT(t)*SQRT(t+s)))
         ENDIF
         t = Z(nn-7) + (s+Z(nn-5))
         Z(nn-3) = Z(nn-3)*(Z(nn-7)/t)
         Z(nn-7) = t
      ENDIF
      Z(4*N0-7) = Z(nn-7) + Sigma
      Z(4*N0-3) = Z(nn-3) + Sigma
      N0 = N0 - 2
      GOTO 100
!
!     Risk of underflow.
!
 300  CALL SLASQ6(I0,N0,Z,Pp,Dmin,Dmin1,Dmin2,Dn,Dn1,Dn2)
      Ndiv = Ndiv + (N0-I0+2)
      Iter = Iter + 1
      Tau = ZERO
!
 400  IF ( Tau<Sigma ) THEN
         Desig = Desig + Tau
         t = Sigma + Desig
         Desig = Desig - (t-Sigma)
      ELSE
         t = Sigma + Tau
         Desig = Sigma - (t-Tau) + Desig
      ENDIF
      Sigma = t
!
!
!     End of SLASQ3
!
      END SUBROUTINE SLASQ3
