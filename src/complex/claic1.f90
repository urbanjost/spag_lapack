!*==claic1.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b CLAIC1 applies one step of incremental condition estimation.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CLAIC1 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/claic1.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/claic1.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/claic1.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CLAIC1( JOB, J, X, SEST, W, GAMMA, SESTPR, S, C )
!
!       .. Scalar Arguments ..
!       INTEGER            J, JOB
!       REAL               SEST, SESTPR
!       COMPLEX            C, GAMMA, S
!       ..
!       .. Array Arguments ..
!       COMPLEX            W( J ), X( J )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CLAIC1 applies one step of incremental condition estimation in
!> its simplest version:
!>
!> Let x, twonorm(x) = 1, be an approximate singular vector of an j-by-j
!> lower triangular matrix L, such that
!>          twonorm(L*x) = sest
!> Then CLAIC1 computes sestpr, s, c such that
!> the vector
!>                 [ s*x ]
!>          xhat = [  c  ]
!> is an approximate singular vector of
!>                 [ L      0  ]
!>          Lhat = [ w**H gamma ]
!> in the sense that
!>          twonorm(Lhat*xhat) = sestpr.
!>
!> Depending on JOB, an estimate for the largest or smallest singular
!> value is computed.
!>
!> Note that [s c]**H and sestpr**2 is an eigenpair of the system
!>
!>     diag(sest*sest, 0) + [alpha  gamma] * [ conjg(alpha) ]
!>                                           [ conjg(gamma) ]
!>
!> where  alpha =  x**H*w.
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
!>          X is COMPLEX array, dimension (J)
!>          The j-vector x.
!> \endverbatim
!>
!> \param[in] SEST
!> \verbatim
!>          SEST is REAL
!>          Estimated singular value of j by j matrix L
!> \endverbatim
!>
!> \param[in] W
!> \verbatim
!>          W is COMPLEX array, dimension (J)
!>          The j-vector w.
!> \endverbatim
!>
!> \param[in] GAMMA
!> \verbatim
!>          GAMMA is COMPLEX
!>          The diagonal element gamma.
!> \endverbatim
!>
!> \param[out] SESTPR
!> \verbatim
!>          SESTPR is REAL
!>          Estimated singular value of (j+1) by (j+1) matrix Lhat.
!> \endverbatim
!>
!> \param[out] S
!> \verbatim
!>          S is COMPLEX
!>          Sine needed in forming xhat.
!> \endverbatim
!>
!> \param[out] C
!> \verbatim
!>          C is COMPLEX
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
!> \ingroup complexOTHERauxiliary
!
!  =====================================================================
      SUBROUTINE CLAIC1(Job,J,X,Sest,W,Gamma,Sestpr,S,C)
      IMPLICIT NONE
!*--CLAIC1139
!
!  -- LAPACK auxiliary routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER J , Job
      REAL Sest , Sestpr
      COMPLEX C , Gamma , S
!     ..
!     .. Array Arguments ..
      COMPLEX W(J) , X(J)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ZERO , ONE , TWO
      PARAMETER (ZERO=0.0E0,ONE=1.0E0,TWO=2.0E0)
      REAL HALF , FOUR
      PARAMETER (HALF=0.5E0,FOUR=4.0E0)
!     ..
!     .. Local Scalars ..
      REAL absalp , absest , absgam , b , eps , norma , s1 , s2 , scl , &
     &     t , test , tmp , zeta1 , zeta2
      COMPLEX alpha , cosine , sine
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , CONJG , MAX , SQRT
!     ..
!     .. External Functions ..
      REAL SLAMCH
      COMPLEX CDOTC
      EXTERNAL SLAMCH , CDOTC
!     ..
!     .. Executable Statements ..
!
      eps = SLAMCH('Epsilon')
      alpha = CDOTC(J,X,1,W,1)
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
               tmp = SQRT(S*CONJG(S)+C*CONJG(C))
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
               scl = SQRT(ONE+tmp*tmp)
               Sestpr = s2*scl
               S = (alpha/s2)/scl
               C = (Gamma/s2)/scl
            ELSE
               tmp = s2/s1
               scl = SQRT(ONE+tmp*tmp)
               Sestpr = s1*scl
               S = (alpha/s1)/scl
               C = (Gamma/s1)/scl
            ENDIF
            RETURN
         ELSE
!
!           normal case
!
            zeta1 = absalp/absest
            zeta2 = absgam/absest
!
            b = (ONE-zeta1*zeta1-zeta2*zeta2)*HALF
            C = zeta1*zeta1
            IF ( b>ZERO ) THEN
               t = C/(b+SQRT(b*b+C))
            ELSE
               t = SQRT(b*b+C) - b
            ENDIF
!
            sine = -(alpha/absest)/t
            cosine = -(Gamma/absest)/(ONE+t)
            tmp = SQRT(sine*CONJG(sine)+cosine*CONJG(cosine))
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
               sine = -CONJG(Gamma)
               cosine = CONJG(alpha)
            ENDIF
            s1 = MAX(ABS(sine),ABS(cosine))
            S = sine/s1
            C = cosine/s1
            tmp = SQRT(S*CONJG(S)+C*CONJG(C))
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
               scl = SQRT(ONE+tmp*tmp)
               Sestpr = absest*(tmp/scl)
               S = -(CONJG(Gamma)/s2)/scl
               C = (CONJG(alpha)/s2)/scl
            ELSE
               tmp = s2/s1
               scl = SQRT(ONE+tmp*tmp)
               Sestpr = absest/scl
               S = -(CONJG(Gamma)/s1)/scl
               C = (CONJG(alpha)/s1)/scl
            ENDIF
            RETURN
         ELSE
!
!           normal case
!
            zeta1 = absalp/absest
            zeta2 = absgam/absest
!
            norma = MAX(ONE+zeta1*zeta1+zeta1*zeta2,                    &
     &              zeta1*zeta2+zeta2*zeta2)
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
               sine = (alpha/absest)/(ONE-t)
               cosine = -(Gamma/absest)/t
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
               sine = -(alpha/absest)/t
               cosine = -(Gamma/absest)/(ONE+t)
               Sestpr = SQRT(ONE+t+FOUR*eps*eps*norma)*absest
            ENDIF
            tmp = SQRT(sine*CONJG(sine)+cosine*CONJG(cosine))
            S = sine/tmp
            C = cosine/tmp
            RETURN
!
         ENDIF
      ENDIF
!
!     End of CLAIC1
!
      END SUBROUTINE CLAIC1
