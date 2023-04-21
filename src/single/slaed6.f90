!*==slaed6.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b SLAED6 used by sstedc. Computes one Newton step in solution of the secular equation.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SLAED6 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slaed6.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slaed6.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slaed6.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SLAED6( KNITER, ORGATI, RHO, D, Z, FINIT, TAU, INFO )
!
!       .. Scalar Arguments ..
!       LOGICAL            ORGATI
!       INTEGER            INFO, KNITER
!       REAL               FINIT, RHO, TAU
!       ..
!       .. Array Arguments ..
!       REAL               D( 3 ), Z( 3 )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SLAED6 computes the positive or negative root (closest to the origin)
!> of
!>                  z(1)        z(2)        z(3)
!> f(x) =   rho + --------- + ---------- + ---------
!>                 d(1)-x      d(2)-x      d(3)-x
!>
!> It is assumed that
!>
!>       if ORGATI = .true. the root is between d(2) and d(3);
!>       otherwise it is between d(1) and d(2)
!>
!> This routine will be called by SLAED4 when necessary. In most cases,
!> the root sought is the smallest in magnitude, though it might not be
!> in some extremely rare situations.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] KNITER
!> \verbatim
!>          KNITER is INTEGER
!>               Refer to SLAED4 for its significance.
!> \endverbatim
!>
!> \param[in] ORGATI
!> \verbatim
!>          ORGATI is LOGICAL
!>               If ORGATI is true, the needed root is between d(2) and
!>               d(3); otherwise it is between d(1) and d(2).  See
!>               SLAED4 for further details.
!> \endverbatim
!>
!> \param[in] RHO
!> \verbatim
!>          RHO is REAL
!>               Refer to the equation f(x) above.
!> \endverbatim
!>
!> \param[in] D
!> \verbatim
!>          D is REAL array, dimension (3)
!>               D satisfies d(1) < d(2) < d(3).
!> \endverbatim
!>
!> \param[in] Z
!> \verbatim
!>          Z is REAL array, dimension (3)
!>               Each of the elements in z must be positive.
!> \endverbatim
!>
!> \param[in] FINIT
!> \verbatim
!>          FINIT is REAL
!>               The value of f at 0. It is more accurate than the one
!>               evaluated inside this routine (if someone wants to do
!>               so).
!> \endverbatim
!>
!> \param[out] TAU
!> \verbatim
!>          TAU is REAL
!>               The root of the equation f(x).
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>               = 0: successful exit
!>               > 0: if INFO = 1, failure to converge
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
!> \ingroup auxOTHERcomputational
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  10/02/03: This version has a few statements commented out for thread
!>  safety (machine parameters are computed on each entry). SJH.
!>
!>  05/10/06: Modified from a new version of Ren-Cang Li, use
!>     Gragg-Thornton-Warner cubic convergent scheme for better stability.
!> \endverbatim
!
!> \par Contributors:
!  ==================
!>
!>     Ren-Cang Li, Computer Science Division, University of California
!>     at Berkeley, USA
!>
!  =====================================================================
      SUBROUTINE SLAED6(Kniter,Orgati,Rho,D,Z,Finit,Tau,Info)
      IMPLICIT NONE
!*--SLAED6144
!
!  -- LAPACK computational routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      LOGICAL Orgati
      INTEGER Info , Kniter
      REAL Finit , Rho , Tau
!     ..
!     .. Array Arguments ..
      REAL D(3) , Z(3)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      INTEGER MAXIT
      PARAMETER (MAXIT=40)
      REAL ZERO , ONE , TWO , THREE , FOUR , EIGHT
      PARAMETER (ZERO=0.0E0,ONE=1.0E0,TWO=2.0E0,THREE=3.0E0,FOUR=4.0E0, &
     &           EIGHT=8.0E0)
!     ..
!     .. External Functions ..
      REAL SLAMCH
      EXTERNAL SLAMCH
!     ..
!     .. Local Arrays ..
      REAL dscale(3) , zscale(3)
!     ..
!     .. Local Scalars ..
      LOGICAL scale
      INTEGER i , iter , niter
      REAL a , b , base , c , ddf , df , eps , erretm , eta , f , fc ,  &
     &     sclfac , sclinv , small1 , small2 , sminv1 , sminv2 , temp , &
     &     temp1 , temp2 , temp3 , temp4 , lbd , ubd
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , INT , LOG , MAX , MIN , SQRT
!     ..
!     .. Executable Statements ..
!
      Info = 0
!
      IF ( Orgati ) THEN
         lbd = D(2)
         ubd = D(3)
      ELSE
         lbd = D(1)
         ubd = D(2)
      ENDIF
      IF ( Finit<ZERO ) THEN
         lbd = ZERO
      ELSE
         ubd = ZERO
      ENDIF
!
      niter = 1
      Tau = ZERO
      IF ( Kniter==2 ) THEN
         IF ( Orgati ) THEN
            temp = (D(3)-D(2))/TWO
            c = Rho + Z(1)/((D(1)-D(2))-temp)
            a = c*(D(2)+D(3)) + Z(2) + Z(3)
            b = c*D(2)*D(3) + Z(2)*D(3) + Z(3)*D(2)
         ELSE
            temp = (D(1)-D(2))/TWO
            c = Rho + Z(3)/((D(3)-D(2))-temp)
            a = c*(D(1)+D(2)) + Z(1) + Z(2)
            b = c*D(1)*D(2) + Z(1)*D(2) + Z(2)*D(1)
         ENDIF
         temp = MAX(ABS(a),ABS(b),ABS(c))
         a = a/temp
         b = b/temp
         c = c/temp
         IF ( c==ZERO ) THEN
            Tau = b/a
         ELSEIF ( a<=ZERO ) THEN
            Tau = (a-SQRT(ABS(a*a-FOUR*b*c)))/(TWO*c)
         ELSE
            Tau = TWO*b/(a+SQRT(ABS(a*a-FOUR*b*c)))
         ENDIF
         IF ( Tau<lbd .OR. Tau>ubd ) Tau = (lbd+ubd)/TWO
         IF ( D(1)==Tau .OR. D(2)==Tau .OR. D(3)==Tau ) THEN
            Tau = ZERO
         ELSE
            temp = Finit + Tau*Z(1)/(D(1)*(D(1)-Tau)) + Tau*Z(2)        &
     &             /(D(2)*(D(2)-Tau)) + Tau*Z(3)/(D(3)*(D(3)-Tau))
            IF ( temp<=ZERO ) THEN
               lbd = Tau
            ELSE
               ubd = Tau
            ENDIF
            IF ( ABS(Finit)<=ABS(temp) ) Tau = ZERO
         ENDIF
      ENDIF
!
!     get machine parameters for possible scaling to avoid overflow
!
!     modified by Sven: parameters SMALL1, SMINV1, SMALL2,
!     SMINV2, EPS are not SAVEd anymore between one call to the
!     others but recomputed at each call
!
      eps = SLAMCH('Epsilon')
      base = SLAMCH('Base')
      small1 = base**(INT(LOG(SLAMCH('SafMin'))/LOG(base)/THREE))
      sminv1 = ONE/small1
      small2 = small1*small1
      sminv2 = sminv1*sminv1
!
!     Determine if scaling of inputs necessary to avoid overflow
!     when computing 1/TEMP**3
!
      IF ( Orgati ) THEN
         temp = MIN(ABS(D(2)-Tau),ABS(D(3)-Tau))
      ELSE
         temp = MIN(ABS(D(1)-Tau),ABS(D(2)-Tau))
      ENDIF
      scale = .FALSE.
      IF ( temp<=small1 ) THEN
         scale = .TRUE.
         IF ( temp<=small2 ) THEN
!
!        Scale up by power of radix nearest 1/SAFMIN**(2/3)
!
            sclfac = sminv2
            sclinv = small2
         ELSE
!
!        Scale up by power of radix nearest 1/SAFMIN**(1/3)
!
            sclfac = sminv1
            sclinv = small1
         ENDIF
!
!        Scaling up safe because D, Z, TAU scaled elsewhere to be O(1)
!
         DO i = 1 , 3
            dscale(i) = D(i)*sclfac
            zscale(i) = Z(i)*sclfac
         ENDDO
         Tau = Tau*sclfac
         lbd = lbd*sclfac
         ubd = ubd*sclfac
      ELSE
!
!        Copy D and Z to DSCALE and ZSCALE
!
         DO i = 1 , 3
            dscale(i) = D(i)
            zscale(i) = Z(i)
         ENDDO
      ENDIF
!
      fc = ZERO
      df = ZERO
      ddf = ZERO
      DO i = 1 , 3
         temp = ONE/(dscale(i)-Tau)
         temp1 = zscale(i)*temp
         temp2 = temp1*temp
         temp3 = temp2*temp
         fc = fc + temp1/dscale(i)
         df = df + temp2
         ddf = ddf + temp3
      ENDDO
      f = Finit + Tau*fc
!
      IF ( ABS(f)>ZERO ) THEN
         IF ( f<=ZERO ) THEN
            lbd = Tau
         ELSE
            ubd = Tau
         ENDIF
!
!        Iteration begins -- Use Gragg-Thornton-Warner cubic convergent
!                            scheme
!
!     It is not hard to see that
!
!           1) Iterations will go up monotonically
!              if FINIT < 0;
!
!           2) Iterations will go down monotonically
!              if FINIT > 0.
!
         iter = niter + 1
!
         DO niter = iter , MAXIT
!
            IF ( Orgati ) THEN
               temp1 = dscale(2) - Tau
               temp2 = dscale(3) - Tau
            ELSE
               temp1 = dscale(1) - Tau
               temp2 = dscale(2) - Tau
            ENDIF
            a = (temp1+temp2)*f - temp1*temp2*df
            b = temp1*temp2*f
            c = f - (temp1+temp2)*df + temp1*temp2*ddf
            temp = MAX(ABS(a),ABS(b),ABS(c))
            a = a/temp
            b = b/temp
            c = c/temp
            IF ( c==ZERO ) THEN
               eta = b/a
            ELSEIF ( a<=ZERO ) THEN
               eta = (a-SQRT(ABS(a*a-FOUR*b*c)))/(TWO*c)
            ELSE
               eta = TWO*b/(a+SQRT(ABS(a*a-FOUR*b*c)))
            ENDIF
            IF ( f*eta>=ZERO ) eta = -f/df
!
            Tau = Tau + eta
            IF ( Tau<lbd .OR. Tau>ubd ) Tau = (lbd+ubd)/TWO
!
            fc = ZERO
            erretm = ZERO
            df = ZERO
            ddf = ZERO
            DO i = 1 , 3
               IF ( (dscale(i)-Tau)==ZERO ) GOTO 100
               temp = ONE/(dscale(i)-Tau)
               temp1 = zscale(i)*temp
               temp2 = temp1*temp
               temp3 = temp2*temp
               temp4 = temp1/dscale(i)
               fc = fc + temp4
               erretm = erretm + ABS(temp4)
               df = df + temp2
               ddf = ddf + temp3
            ENDDO
            f = Finit + Tau*fc
            erretm = EIGHT*(ABS(Finit)+ABS(Tau)*erretm) + ABS(Tau)*df
            IF ( (ABS(f)<=FOUR*eps*erretm) .OR.                         &
     &           ((ubd-lbd)<=FOUR*eps*ABS(Tau)) ) GOTO 100
            IF ( f<=ZERO ) THEN
               lbd = Tau
            ELSE
               ubd = Tau
            ENDIF
         ENDDO
         Info = 1
      ENDIF
!
!     Undo scaling
!
 100  IF ( scale ) Tau = Tau*sclinv
!
!     End of SLAED6
!
      END SUBROUTINE SLAED6
