!*==dlaic1.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b DLAIC1 applies one step of incremental condition estimation.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DLAIC1 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlaic1.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlaic1.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlaic1.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLAIC1( JOB, J, X, SEST, W, GAMMA, SESTPR, S, C )
!
!       .. Scalar Arguments ..
!       INTEGER            J, JOB
!       DOUBLE PRECISION   C, GAMMA, S, SEST, SESTPR
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   W( J ), X( J )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLAIC1 applies one step of incremental condition estimation in
!> its simplest version:
!>
!> Let x, twonorm(x) = 1, be an approximate singular vector of an j-by-j
!> lower triangular matrix L, such that
!>          twonorm(L*x) = sest
!> Then DLAIC1 computes sestpr, s, c such that
!> the vector
!>                 [ s*x ]
!>          xhat = [  c  ]
!> is an approximate singular vector of
!>                 [ L       0  ]
!>          Lhat = [ w**T gamma ]
!> in the sense that
!>          twonorm(Lhat*xhat) = sestpr.
!>
!> Depending on JOB, an estimate for the largest or smallest singular
!> value is computed.
!>
!> Note that [s c]**T and sestpr**2 is an eigenpair of the system
!>
!>     diag(sest*sest, 0) + [alpha  gamma] * [ alpha ]
!>                                           [ gamma ]
!>
!> where  alpha =  x**T*w.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] JOB
!> \verbatim
!>          JOB is INTEGER
!>          = 1: an estimate for the largest singular value is computed.
!>          = 2: an estimate for the smallest singular value is computed.
!> \endverbatim
!>
!> \param[in] J
!> \verbatim
!>          J is INTEGER
!>          Length of X and W
!> \endverbatim
!>
!> \param[in] X
!> \verbatim
!>          X is DOUBLE PRECISION array, dimension (J)
!>          The j-vector x.
!> \endverbatim
!>
!> \param[in] SEST
!> \verbatim
!>          SEST is DOUBLE PRECISION
!>          Estimated singular value of j by j matrix L
!> \endverbatim
!>
!> \param[in] W
!> \verbatim
!>          W is DOUBLE PRECISION array, dimension (J)
!>          The j-vector w.
!> \endverbatim
!>
!> \param[in] GAMMA
!> \verbatim
!>          GAMMA is DOUBLE PRECISION
!>          The diagonal element gamma.
!> \endverbatim
!>
!> \param[out] SESTPR
!> \verbatim
!>          SESTPR is DOUBLE PRECISION
!>          Estimated singular value of (j+1) by (j+1) matrix Lhat.
!> \endverbatim
!>
!> \param[out] S
!> \verbatim
!>          S is DOUBLE PRECISION
!>          Sine needed in forming xhat.
!> \endverbatim
!>
!> \param[out] C
!> \verbatim
!>          C is DOUBLE PRECISION
!>          Cosine needed in forming xhat.
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
!> \ingroup doubleOTHERauxiliary
!
!  =====================================================================
      SUBROUTINE DLAIC1(Job,J,X,Sest,W,Gamma,Sestpr,S,C)
      IMPLICIT NONE
!*--DLAIC1138
!
!  -- LAPACK auxiliary routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER J , Job
      DOUBLE PRECISION C , Gamma , S , Sest , Sestpr
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION W(J) , X(J)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO , ONE , TWO
      PARAMETER (ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0)
      DOUBLE PRECISION HALF , FOUR
      PARAMETER (HALF=0.5D0,FOUR=4.0D0)
!     ..
!     .. Local Scalars ..
      DOUBLE PRECISION absalp , absest , absgam , alpha , b , cosine ,  &
     &                 eps , norma , s1 , s2 , sine , t , test , tmp ,  &
     &                 zeta1 , zeta2
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , MAX , SIGN , SQRT
!     ..
!     .. External Functions ..
      DOUBLE PRECISION DDOT , DLAMCH
      EXTERNAL DDOT , DLAMCH
!     ..
!     .. Executable Statements ..
!
      eps = DLAMCH('Epsilon')
      alpha = DDOT(J,X,1,W,1)
!
      absalp = ABS(alpha)
      absgam = ABS(Gamma)
      absest = ABS(Sest)
!
      IF ( Job==1 ) THEN
!
!        Estimating largest singular value
!
!        special cases
!
         IF ( Sest==ZERO ) THEN
            s1 = MAX(absgam,absalp)
            IF ( s1==ZERO ) THEN
               S = ZERO
               C = ONE
               Sestpr = ZERO
            ELSE
               S = alpha/s1
               C = Gamma/s1
               tmp = SQRT(S*S+C*C)
               S = S/tmp
               C = C/tmp
               Sestpr = s1*tmp
            ENDIF
            RETURN
         ELSEIF ( absgam<=eps*absest ) THEN
            S = ONE
            C = ZERO
            tmp = MAX(absest,absalp)
            s1 = absest/tmp
            s2 = absalp/tmp
            Sestpr = tmp*SQRT(s1*s1+s2*s2)
            RETURN
         ELSEIF ( absalp<=eps*absest ) THEN
            s1 = absgam
            s2 = absest
            IF ( s1<=s2 ) THEN
               S = ONE
               C = ZERO
               Sestpr = s2
            ELSE
               S = ZERO
               C = ONE
               Sestpr = s1
            ENDIF
            RETURN
         ELSEIF ( absest<=eps*absalp .OR. absest<=eps*absgam ) THEN
            s1 = absgam
            s2 = absalp
            IF ( s1<=s2 ) THEN
               tmp = s1/s2
               S = SQRT(ONE+tmp*tmp)
               Sestpr = s2*S
               C = (Gamma/s2)/S
               S = SIGN(ONE,alpha)/S
            ELSE
               tmp = s2/s1
               C = SQRT(ONE+tmp*tmp)
               Sestpr = s1*C
               S = (alpha/s1)/C
               C = SIGN(ONE,Gamma)/C
            ENDIF
            RETURN
         ELSE
!
!           normal case
!
            zeta1 = alpha/absest
            zeta2 = Gamma/absest
!
            b = (ONE-zeta1*zeta1-zeta2*zeta2)*HALF
            C = zeta1*zeta1
            IF ( b>ZERO ) THEN
               t = C/(b+SQRT(b*b+C))
            ELSE
               t = SQRT(b*b+C) - b
            ENDIF
!
            sine = -zeta1/t
            cosine = -zeta2/(ONE+t)
            tmp = SQRT(sine*sine+cosine*cosine)
            S = sine/tmp
            C = cosine/tmp
            Sestpr = SQRT(t+ONE)*absest
            RETURN
         ENDIF
!
      ELSEIF ( Job==2 ) THEN
!
!        Estimating smallest singular value
!
!        special cases
!
         IF ( Sest==ZERO ) THEN
            Sestpr = ZERO
            IF ( MAX(absgam,absalp)==ZERO ) THEN
               sine = ONE
               cosine = ZERO
            ELSE
               sine = -Gamma
               cosine = alpha
            ENDIF
            s1 = MAX(ABS(sine),ABS(cosine))
            S = sine/s1
            C = cosine/s1
            tmp = SQRT(S*S+C*C)
            S = S/tmp
            C = C/tmp
            RETURN
         ELSEIF ( absgam<=eps*absest ) THEN
            S = ZERO
            C = ONE
            Sestpr = absgam
            RETURN
         ELSEIF ( absalp<=eps*absest ) THEN
            s1 = absgam
            s2 = absest
            IF ( s1<=s2 ) THEN
               S = ZERO
               C = ONE
               Sestpr = s1
            ELSE
               S = ONE
               C = ZERO
               Sestpr = s2
            ENDIF
            RETURN
         ELSEIF ( absest<=eps*absalp .OR. absest<=eps*absgam ) THEN
            s1 = absgam
            s2 = absalp
            IF ( s1<=s2 ) THEN
               tmp = s1/s2
               C = SQRT(ONE+tmp*tmp)
               Sestpr = absest*(tmp/C)
               S = -(Gamma/s2)/C
               C = SIGN(ONE,alpha)/C
            ELSE
               tmp = s2/s1
               S = SQRT(ONE+tmp*tmp)
               Sestpr = absest/S
               C = (alpha/s1)/S
               S = -SIGN(ONE,Gamma)/S
            ENDIF
            RETURN
         ELSE
!
!           normal case
!
            zeta1 = alpha/absest
            zeta2 = Gamma/absest
!
            norma = MAX(ONE+zeta1*zeta1+ABS(zeta1*zeta2),               &
     &              ABS(zeta1*zeta2)+zeta2*zeta2)
!
!           See if root is closer to zero or to ONE
!
            test = ONE + TWO*(zeta1-zeta2)*(zeta1+zeta2)
            IF ( test>=ZERO ) THEN
!
!              root is close to zero, compute directly
!
               b = (zeta1*zeta1+zeta2*zeta2+ONE)*HALF
               C = zeta2*zeta2
               t = C/(b+SQRT(ABS(b*b-C)))
               sine = zeta1/(ONE-t)
               cosine = -zeta2/t
               Sestpr = SQRT(t+FOUR*eps*eps*norma)*absest
            ELSE
!
!              root is closer to ONE, shift by that amount
!
               b = (zeta2*zeta2+zeta1*zeta1-ONE)*HALF
               C = zeta1*zeta1
               IF ( b>=ZERO ) THEN
                  t = -C/(b+SQRT(b*b+C))
               ELSE
                  t = b - SQRT(b*b+C)
               ENDIF
               sine = -zeta1/t
               cosine = -zeta2/(ONE+t)
               Sestpr = SQRT(ONE+t+FOUR*eps*eps*norma)*absest
            ENDIF
            tmp = SQRT(sine*sine+cosine*cosine)
            S = sine/tmp
            C = cosine/tmp
            RETURN
!
         ENDIF
      ENDIF
!
!     End of DLAIC1
!
      END SUBROUTINE DLAIC1
