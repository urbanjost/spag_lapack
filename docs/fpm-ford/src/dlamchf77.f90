!*==dlamch.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!
!***********************************************************************
!> \brief \b DLAMCHF77 deprecated
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!      DOUBLE PRECISION FUNCTION DLAMCH( CMACH )
!
!     .. Scalar Arguments ..
!     CHARACTER          CMACH
!     ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLAMCHF77 determines double precision machine parameters.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] CMACH
!> \verbatim
!>          Specifies the value to be returned by DLAMCH:
!>          = 'E' or 'e',   DLAMCH := eps
!>          = 'S' or 's ,   DLAMCH := sfmin
!>          = 'B' or 'b',   DLAMCH := base
!>          = 'P' or 'p',   DLAMCH := eps*base
!>          = 'N' or 'n',   DLAMCH := t
!>          = 'R' or 'r',   DLAMCH := rnd
!>          = 'M' or 'm',   DLAMCH := emin
!>          = 'U' or 'u',   DLAMCH := rmin
!>          = 'L' or 'l',   DLAMCH := emax
!>          = 'O' or 'o',   DLAMCH := rmax
!>          where
!>          eps   = relative machine precision
!>          sfmin = safe minimum, such that 1/sfmin does not overflow
!>          base  = base of the machine
!>          prec  = eps*base
!>          t     = number of (base) digits in the mantissa
!>          rnd   = 1.0 when rounding occurs in addition, 0.0 otherwise
!>          emin  = minimum exponent before (gradual) underflow
!>          rmin  = underflow threshold - base**(emin-1)
!>          emax  = largest exponent before overflow
!>          rmax  = overflow threshold  - (base**emax)*(1-eps)
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
!> \date April 2012
!
!> \ingroup auxOTHERauxiliary
!
!  =====================================================================
      DOUBLE PRECISION FUNCTION DLAMCH(Cmach)
      IMPLICIT NONE
!*--DLAMCH73
!
!  -- LAPACK auxiliary routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     April 2012
!
!     .. Scalar Arguments ..
      CHARACTER Cmach
!     ..
!     .. Parameters ..
      DOUBLE PRECISION ONE , ZERO
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
!     ..
!     .. Local Scalars ..
      LOGICAL first , lrnd
      INTEGER beta , imax , imin , it
      DOUBLE PRECISION base , emax , emin , eps , prec , rmach , rmax , &
     &                 rmin , rnd , sfmin , small , t
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL DLAMC2
!     ..
!     .. Save statement ..
      SAVE first , eps , sfmin , base , t , rnd , emin , rmin , emax ,  &
     &   rmax , prec
!     ..
!     .. Data statements ..
      DATA first/.TRUE./
!     ..
!     .. Executable Statements ..
!
      IF ( first ) THEN
         CALL DLAMC2(beta,it,lrnd,eps,imin,rmin,imax,rmax)
         base = beta
         t = it
         IF ( lrnd ) THEN
            rnd = ONE
            eps = (base**(1-it))/2
         ELSE
            rnd = ZERO
            eps = base**(1-it)
         ENDIF
         prec = eps*base
         emin = imin
         emax = imax
         sfmin = rmin
         small = ONE/rmax
!
!           Use SMALL plus a bit, to avoid the possibility of rounding
!           causing overflow when computing  1/sfmin.
!
         IF ( small>=sfmin ) sfmin = small*(ONE+eps)
      ENDIF
!
      IF ( LSAME(Cmach,'E') ) THEN
         rmach = eps
      ELSEIF ( LSAME(Cmach,'S') ) THEN
         rmach = sfmin
      ELSEIF ( LSAME(Cmach,'B') ) THEN
         rmach = base
      ELSEIF ( LSAME(Cmach,'P') ) THEN
         rmach = prec
      ELSEIF ( LSAME(Cmach,'N') ) THEN
         rmach = t
      ELSEIF ( LSAME(Cmach,'R') ) THEN
         rmach = rnd
      ELSEIF ( LSAME(Cmach,'M') ) THEN
         rmach = emin
      ELSEIF ( LSAME(Cmach,'U') ) THEN
         rmach = rmin
      ELSEIF ( LSAME(Cmach,'L') ) THEN
         rmach = emax
      ELSEIF ( LSAME(Cmach,'O') ) THEN
         rmach = rmax
      ENDIF
!
      DLAMCH = rmach
      first = .FALSE.
!
!     End of DLAMCH
!
      END FUNCTION DLAMCH
!*==dlamc1.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!
!***********************************************************************
!
!> \brief \b DLAMC1
!> \details
!> \b Purpose:
!> \verbatim
!> DLAMC1 determines the machine parameters given by BETA, T, RND, and
!> IEEE1.
!> \endverbatim
!>
!> \param[out] BETA
!> \verbatim
!>          The base of the machine.
!> \endverbatim
!>
!> \param[out] T
!> \verbatim
!>          The number of ( BETA ) digits in the mantissa.
!> \endverbatim
!>
!> \param[out] RND
!> \verbatim
!>          Specifies whether proper rounding  ( RND = .TRUE. )  or
!>          chopping  ( RND = .FALSE. )  occurs in addition. This may not
!>          be a reliable guide to the way in which the machine performs
!>          its arithmetic.
!> \endverbatim
!>
!> \param[out] IEEE1
!> \verbatim
!>          Specifies whether rounding appears to be done in the IEEE
!>          'round to nearest' style.
!> \endverbatim
!> \author LAPACK is a software package provided by Univ. of Tennessee, Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..
!> \date April 2012
!> \ingroup auxOTHERauxiliary
!>
!> \details \b Further \b Details
!> \verbatim
!>
!>  The routine is based on the routine  ENVRON  by Malcolm and
!>  incorporates suggestions by Gentleman and Marovich. See
!>
!>     Malcolm M. A. (1972) Algorithms to reveal properties of
!>        floating-point arithmetic. Comms. of the ACM, 15, 949-951.
!>
!>     Gentleman W. M. and Marovich S. B. (1974) More on algorithms
!>        that reveal properties of floating point arithmetic units.
!>        Comms. of the ACM, 17, 276-277.
!> \endverbatim
!>
      SUBROUTINE DLAMC1(Beta,T,Rnd,Ieee1)
      IMPLICIT NONE
!*--DLAMC1215
!
!  -- LAPACK auxiliary routine (version 3.7.0) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2010
!
!     .. Scalar Arguments ..
      LOGICAL Ieee1 , Rnd
      INTEGER Beta , T
!     ..
! =====================================================================
!
!     .. Local Scalars ..
      LOGICAL first , lieee1 , lrnd
      INTEGER lbeta , lt
      DOUBLE PRECISION a , b , c , f , one , qtr , savec , t1 , t2
!     ..
!     .. External Functions ..
      DOUBLE PRECISION DLAMC3
      EXTERNAL DLAMC3
!     ..
!     .. Save statement ..
      SAVE first , lieee1 , lbeta , lrnd , lt
!     ..
!     .. Data statements ..
      DATA first/.TRUE./
!     ..
!     .. Executable Statements ..
!
      IF ( first ) THEN
         one = 1
!
!        LBETA,  LIEEE1,  LT and  LRND  are the  local values  of  BETA,
!        IEEE1, T and RND.
!
!        Throughout this routine  we use the function  DLAMC3  to ensure
!        that relevant values are  stored and not held in registers,  or
!        are not affected by optimizers.
!
!        Compute  a = 2.0**m  with the  smallest positive integer m such
!        that
!
!           fl( a + 1.0 ) = a.
!
         a = 1
         c = 1
         DO
!
!+       WHILE( C.EQ.ONE )LOOP
            IF ( c==one ) THEN
               a = 2*a
               c = DLAMC3(a,one)
               c = DLAMC3(c,-a)
               CYCLE
            ENDIF
!+       END WHILE
!
!        Now compute  b = 2.0**m  with the smallest positive integer m
!        such that
!
!           fl( a + b ) .gt. a.
!
            b = 1
            c = DLAMC3(a,b)
            EXIT
         ENDDO
         DO
!
!+       WHILE( C.EQ.A )LOOP
            IF ( c==a ) THEN
               b = 2*b
               c = DLAMC3(a,b)
               CYCLE
            ENDIF
!+       END WHILE
!
!        Now compute the base.  a and c  are neighbouring floating point
!        numbers  in the  interval  ( beta**t, beta**( t + 1 ) )  and so
!        their difference is beta. Adding 0.25 to c is to ensure that it
!        is truncated to beta and not ( beta - 1 ).
!
            qtr = one/4
            savec = c
            c = DLAMC3(c,-a)
            lbeta = c + qtr
!
!        Now determine whether rounding or chopping occurs,  by adding a
!        bit  less  than  beta/2  and a  bit  more  than  beta/2  to  a.
!
            b = lbeta
            f = DLAMC3(b/2,-b/100)
            c = DLAMC3(f,a)
            IF ( c==a ) THEN
               lrnd = .TRUE.
            ELSE
               lrnd = .FALSE.
            ENDIF
            f = DLAMC3(b/2,b/100)
            c = DLAMC3(f,a)
            IF ( (lrnd) .AND. (c==a) ) lrnd = .FALSE.
!
!        Try and decide whether rounding is done in the  IEEE  'round to
!        nearest' style. B/2 is half a unit in the last place of the two
!        numbers A and SAVEC. Furthermore, A is even, i.e. has last  bit
!        zero, and SAVEC is odd. Thus adding B/2 to A should not  change
!        A, but adding B/2 to SAVEC should change SAVEC.
!
            t1 = DLAMC3(b/2,a)
            t2 = DLAMC3(b/2,savec)
            lieee1 = (t1==a) .AND. (t2>savec) .AND. lrnd
!
!        Now find  the  mantissa, t.  It should  be the  integer part of
!        log to the base beta of a,  however it is safer to determine  t
!        by powering.  So we find t as the smallest positive integer for
!        which
!
!           fl( beta**t + 1.0 ) = 1.0.
!
            lt = 0
            a = 1
            c = 1
            EXIT
         ENDDO
!
!+       WHILE( C.EQ.ONE )LOOP
         DO WHILE ( c==one )
            lt = lt + 1
            a = a*lbeta
            c = DLAMC3(a,one)
            c = DLAMC3(c,-a)
         ENDDO
!+       END WHILE
!
      ENDIF
!
      Beta = lbeta
      T = lt
      Rnd = lrnd
      Ieee1 = lieee1
      first = .FALSE.
!
!     End of DLAMC1
!
      END SUBROUTINE DLAMC1
!*==dlamc2.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!
!***********************************************************************
!
!> \brief \b DLAMC2
!> \details
!> \b Purpose:
!> \verbatim
!> DLAMC2 determines the machine parameters specified in its argument
!> list.
!> \endverbatim
!> \author LAPACK is a software package provided by Univ. of Tennessee, Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..
!> \date April 2012
!> \ingroup auxOTHERauxiliary
!>
!> \param[out] BETA
!> \verbatim
!>          The base of the machine.
!> \endverbatim
!>
!> \param[out] T
!> \verbatim
!>          The number of ( BETA ) digits in the mantissa.
!> \endverbatim
!>
!> \param[out] RND
!> \verbatim
!>          Specifies whether proper rounding  ( RND = .TRUE. )  or
!>          chopping  ( RND = .FALSE. )  occurs in addition. This may not
!>          be a reliable guide to the way in which the machine performs
!>          its arithmetic.
!> \endverbatim
!>
!> \param[out] EPS
!> \verbatim
!>          The smallest positive number such that
!>             fl( 1.0 - EPS ) .LT. 1.0,
!>          where fl denotes the computed value.
!> \endverbatim
!>
!> \param[out] EMIN
!> \verbatim
!>          The minimum exponent before (gradual) underflow occurs.
!> \endverbatim
!>
!> \param[out] RMIN
!> \verbatim
!>          The smallest normalized number for the machine, given by
!>          BASE**( EMIN - 1 ), where  BASE  is the floating point value
!>          of BETA.
!> \endverbatim
!>
!> \param[out] EMAX
!> \verbatim
!>          The maximum exponent before overflow occurs.
!> \endverbatim
!>
!> \param[out] RMAX
!> \verbatim
!>          The largest positive number for the machine, given by
!>          BASE**EMAX * ( 1 - EPS ), where  BASE  is the floating point
!>          value of BETA.
!> \endverbatim
!>
!> \details \b Further \b Details
!> \verbatim
!>
!>  The computation of  EPS  is based on a routine PARANOIA by
!>  W. Kahan of the University of California at Berkeley.
!> \endverbatim
      SUBROUTINE DLAMC2(Beta,T,Rnd,Eps,Emin,Rmin,Emax,Rmax)
      IMPLICIT NONE
!*--DLAMC2431
!
!  -- LAPACK auxiliary routine (version 3.7.0) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2010
!
!     .. Scalar Arguments ..
      LOGICAL Rnd
      INTEGER Beta , Emax , Emin , T
      DOUBLE PRECISION Eps , Rmax , Rmin
!     ..
! =====================================================================
!
!     .. Local Scalars ..
      LOGICAL first , ieee , iwarn , lieee1 , lrnd
      INTEGER gnmin , gpmin , i , lbeta , lemax , lemin , lt , ngnmin , &
     &        ngpmin
      DOUBLE PRECISION a , b , c , half , leps , lrmax , lrmin , one ,  &
     &                 rbase , sixth , small , third , two , zero
!     ..
!     .. External Functions ..
      DOUBLE PRECISION DLAMC3
      EXTERNAL DLAMC3
!     ..
!     .. External Subroutines ..
      EXTERNAL DLAMC1 , DLAMC4 , DLAMC5
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , MAX , MIN
!     ..
!     .. Save statement ..
      SAVE first , iwarn , lbeta , lemax , lemin , leps , lrmax ,       &
     &   lrmin , lt
!     ..
!     .. Data statements ..
      DATA first/.TRUE./ , iwarn/.FALSE./
!     ..
!     .. Executable Statements ..
!
      IF ( first ) THEN
         zero = 0
         one = 1
         two = 2
!
!        LBETA, LT, LRND, LEPS, LEMIN and LRMIN  are the local values of
!        BETA, T, RND, EPS, EMIN and RMIN.
!
!        Throughout this routine  we use the function  DLAMC3  to ensure
!        that relevant values are stored  and not held in registers,  or
!        are not affected by optimizers.
!
!        DLAMC1 returns the parameters  LBETA, LT, LRND and LIEEE1.
!
         CALL DLAMC1(lbeta,lt,lrnd,lieee1)
!
!        Start to find EPS.
!
         b = lbeta
         a = b**(-lt)
         leps = a
!
!        Try some tricks to see whether or not this is the correct  EPS.
!
         b = two/3
         half = one/2
         sixth = DLAMC3(b,-half)
         third = DLAMC3(sixth,sixth)
         b = DLAMC3(third,-half)
         b = DLAMC3(b,sixth)
         b = ABS(b)
         IF ( b<leps ) b = leps
!
         leps = 1
         DO
!
!+       WHILE( ( LEPS.GT.B ).AND.( B.GT.ZERO ) )LOOP
            IF ( (leps>b) .AND. (b>zero) ) THEN
               leps = b
               c = DLAMC3(half*leps,(two**5)*(leps**2))
               c = DLAMC3(half,-c)
               b = DLAMC3(half,c)
               c = DLAMC3(half,-b)
               b = DLAMC3(half,c)
               CYCLE
            ENDIF
!+       END WHILE
!
            IF ( a<leps ) leps = a
!
!        Computation of EPS complete.
!
!        Now find  EMIN.  Let A = + or - 1, and + or - (1 + BASE**(-3)).
!        Keep dividing  A by BETA until (gradual) underflow occurs. This
!        is detected when we cannot recover the previous A.
!
            rbase = one/lbeta
            small = one
            DO i = 1 , 3
               small = DLAMC3(small*rbase,zero)
            ENDDO
            a = DLAMC3(one,small)
            CALL DLAMC4(ngpmin,one,lbeta)
            CALL DLAMC4(ngnmin,-one,lbeta)
            CALL DLAMC4(gpmin,a,lbeta)
            CALL DLAMC4(gnmin,-a,lbeta)
            ieee = .FALSE.
!
            IF ( (ngpmin==ngnmin) .AND. (gpmin==gnmin) ) THEN
               IF ( ngpmin==gpmin ) THEN
                  lemin = ngpmin
!            ( Non twos-complement machines, no gradual underflow;
!              e.g.,  VAX )
               ELSEIF ( (gpmin-ngpmin)==3 ) THEN
                  lemin = ngpmin - 1 + lt
                  ieee = .TRUE.
!            ( Non twos-complement machines, with gradual underflow;
!              e.g., IEEE standard followers )
               ELSE
                  lemin = MIN(ngpmin,gpmin)
!            ( A guess; no known machine )
                  iwarn = .TRUE.
               ENDIF
!
            ELSEIF ( (ngpmin==gpmin) .AND. (ngnmin==gnmin) ) THEN
               IF ( ABS(ngpmin-ngnmin)==1 ) THEN
                  lemin = MAX(ngpmin,ngnmin)
!            ( Twos-complement machines, no gradual underflow;
!              e.g., CYBER 205 )
               ELSE
                  lemin = MIN(ngpmin,ngnmin)
!            ( A guess; no known machine )
                  iwarn = .TRUE.
               ENDIF
!
            ELSEIF ( (ABS(ngpmin-ngnmin)==1) .AND. (gpmin==gnmin) ) THEN
               IF ( (gpmin-MIN(ngpmin,ngnmin))==3 ) THEN
                  lemin = MAX(ngpmin,ngnmin) - 1 + lt
!            ( Twos-complement machines with gradual underflow;
!              no known machine )
               ELSE
                  lemin = MIN(ngpmin,ngnmin)
!            ( A guess; no known machine )
                  iwarn = .TRUE.
               ENDIF
!
            ELSE
               lemin = MIN(ngpmin,ngnmin,gpmin,gnmin)
!         ( A guess; no known machine )
               iwarn = .TRUE.
            ENDIF
            first = .FALSE.
!**
! Comment out this if block if EMIN is ok
            IF ( iwarn ) THEN
               first = .TRUE.
               WRITE (6,FMT=99001) lemin
            ENDIF
!**
!
!        Assume IEEE arithmetic if we found denormalised  numbers above,
!        or if arithmetic seems to round in the  IEEE style,  determined
!        in routine DLAMC1. A true IEEE machine should have both  things
!        true; however, faulty machines may have one or the other.
!
            ieee = ieee .OR. lieee1
!
!        Compute  RMIN by successive division by  BETA. We could compute
!        RMIN as BASE**( EMIN - 1 ),  but some machines underflow during
!        this computation.
!
            lrmin = 1
            DO i = 1 , 1 - lemin
               lrmin = DLAMC3(lrmin*rbase,zero)
            ENDDO
!
!        Finally, call DLAMC5 to compute EMAX and RMAX.
!
            CALL DLAMC5(lbeta,lt,lemin,ieee,lemax,lrmax)
            EXIT
         ENDDO
      ENDIF
!
      Beta = lbeta
      T = lt
      Rnd = lrnd
      Eps = leps
      Emin = lemin
      Rmin = lrmin
      Emax = lemax
      Rmax = lrmax
!
      RETURN
!
99001 FORMAT (//' WARNING. The value EMIN may be incorrect:-',          &
     &        '  EMIN = ',I8,                                           &
     &        /' If, after inspection, the value EMIN looks',           &
     &        ' acceptable please comment out ',                        &
     &        /' the IF block as marked within the code of routine',    &
     &        ' DLAMC2,',/' otherwise supply EMIN explicitly.',/)
!
!     End of DLAMC2
!
      END SUBROUTINE DLAMC2
!*==dlamc3.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!
!***********************************************************************
!
!> \brief \b DLAMC3
!> \details
!> \b Purpose:
!> \verbatim
!> DLAMC3  is intended to force  A  and  B  to be stored prior to doing
!> the addition of  A  and  B ,  for use in situations where optimizers
!> might hold one of these in a register.
!> \endverbatim
!>
!> \param[in] A
!>
!> \param[in] B
!> \verbatim
!>          The values A and B.
!> \endverbatim
 
      DOUBLE PRECISION FUNCTION DLAMC3(A,B)
      IMPLICIT NONE
!*--DLAMC3656
!
!  -- LAPACK auxiliary routine (version 3.7.0) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2010
!
!     .. Scalar Arguments ..
      DOUBLE PRECISION A , B
!     ..
! =====================================================================
!
!     .. Executable Statements ..
!
      DLAMC3 = A + B
!
!
!     End of DLAMC3
!
      END FUNCTION DLAMC3
!*==dlamc4.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!
!***********************************************************************
!
!> \brief \b DLAMC4
!> \details
!> \b Purpose:
!> \verbatim
!> DLAMC4 is a service routine for DLAMC2.
!> \endverbatim
!>
!> \param[out] EMIN
!> \verbatim
!>          The minimum exponent before (gradual) underflow, computed by
!>          setting A = START and dividing by BASE until the previous A
!>          can not be recovered.
!> \endverbatim
!>
!> \param[in] START
!> \verbatim
!>          The starting point for determining EMIN.
!> \endverbatim
!>
!> \param[in] BASE
!> \verbatim
!>          The base of the machine.
!> \endverbatim
!>
      SUBROUTINE DLAMC4(Emin,Start,Base)
      IMPLICIT NONE
!*--DLAMC4705
!
!  -- LAPACK auxiliary routine (version 3.7.0) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2010
!
!     .. Scalar Arguments ..
      INTEGER Base , Emin
      DOUBLE PRECISION Start
!     ..
! =====================================================================
!
!     .. Local Scalars ..
      INTEGER i
      DOUBLE PRECISION a , b1 , b2 , c1 , c2 , d1 , d2 , one , rbase ,  &
     &                 zero
!     ..
!     .. External Functions ..
      DOUBLE PRECISION DLAMC3
      EXTERNAL DLAMC3
!     ..
!     .. Executable Statements ..
!
      a = Start
      one = 1
      rbase = one/Base
      zero = 0
      Emin = 1
      b1 = DLAMC3(a*rbase,zero)
      c1 = a
      c2 = a
      d1 = a
      d2 = a
!+    WHILE( ( C1.EQ.A ).AND.( C2.EQ.A ).AND.
!    $       ( D1.EQ.A ).AND.( D2.EQ.A )      )LOOP
      DO WHILE ( (c1==a) .AND. (c2==a) .AND. (d1==a) .AND. (d2==a) )
         Emin = Emin - 1
         a = b1
         b1 = DLAMC3(a/Base,zero)
         c1 = DLAMC3(b1*Base,zero)
         d1 = zero
         DO i = 1 , Base
            d1 = d1 + b1
         ENDDO
         b2 = DLAMC3(a*rbase,zero)
         c2 = DLAMC3(b2/rbase,zero)
         d2 = zero
         DO i = 1 , Base
            d2 = d2 + b2
         ENDDO
      ENDDO
!+    END WHILE
!
!
!     End of DLAMC4
!
      END SUBROUTINE DLAMC4
!*==dlamc5.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!
!***********************************************************************
!
!> \brief \b DLAMC5
!> \details
!> \b Purpose:
!> \verbatim
!> DLAMC5 attempts to compute RMAX, the largest machine floating-point
!> number, without overflow.  It assumes that EMAX + abs(EMIN) sum
!> approximately to a power of 2.  It will fail on machines where this
!> assumption does not hold, for example, the Cyber 205 (EMIN = -28625,
!> EMAX = 28718).  It will also fail if the value supplied for EMIN is
!> too large (i.e. too close to zero), probably with overflow.
!> \endverbatim
!>
!> \param[in] BETA
!> \verbatim
!>          The base of floating-point arithmetic.
!> \endverbatim
!>
!> \param[in] P
!> \verbatim
!>          The number of base BETA digits in the mantissa of a
!>          floating-point value.
!> \endverbatim
!>
!> \param[in] EMIN
!> \verbatim
!>          The minimum exponent before (gradual) underflow.
!> \endverbatim
!>
!> \param[in] IEEE
!> \verbatim
!>          A logical flag specifying whether or not the arithmetic
!>          system is thought to comply with the IEEE standard.
!> \endverbatim
!>
!> \param[out] EMAX
!> \verbatim
!>          The largest exponent before overflow
!> \endverbatim
!>
!> \param[out] RMAX
!> \verbatim
!>          The largest machine floating-point number.
!> \endverbatim
!>
      SUBROUTINE DLAMC5(Beta,P,Emin,Ieee,Emax,Rmax)
      IMPLICIT NONE
!*--DLAMC5812
!
!  -- LAPACK auxiliary routine (version 3.7.0) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2010
!
!     .. Scalar Arguments ..
      LOGICAL Ieee
      INTEGER Beta , Emax , Emin , P
      DOUBLE PRECISION Rmax
!     ..
! =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO , ONE
      PARAMETER (ZERO=0.0D0,ONE=1.0D0)
!     ..
!     .. Local Scalars ..
      INTEGER exbits , expsum , i , lexp , nbits , try , uexp
      DOUBLE PRECISION oldy , recbas , y , z
!     ..
!     .. External Functions ..
      DOUBLE PRECISION DLAMC3
      EXTERNAL DLAMC3
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MOD
!     ..
!     .. Executable Statements ..
!
!     First compute LEXP and UEXP, two powers of 2 that bound
!     abs(EMIN). We then assume that EMAX + abs(EMIN) will sum
!     approximately to the bound that is closest to abs(EMIN).
!     (EMAX is the exponent of the required number RMAX).
!
      lexp = 1
      exbits = 1
      DO
         try = lexp*2
         IF ( try<=(-Emin) ) THEN
            lexp = try
            exbits = exbits + 1
            CYCLE
         ENDIF
         IF ( lexp==-Emin ) THEN
            uexp = lexp
         ELSE
            uexp = try
            exbits = exbits + 1
         ENDIF
!
!     Now -LEXP is less than or equal to EMIN, and -UEXP is greater
!     than or equal to EMIN. EXBITS is the number of bits needed to
!     store the exponent.
!
         IF ( (uexp+Emin)>(-lexp-Emin) ) THEN
            expsum = 2*lexp
         ELSE
            expsum = 2*uexp
         ENDIF
!
!     EXPSUM is the exponent range, approximately equal to
!     EMAX - EMIN + 1 .
!
         Emax = expsum + Emin - 1
         nbits = 1 + exbits + P
!
!     NBITS is the total number of bits needed to store a
!     floating-point number.
!
!
!        Either there are an odd number of bits used to store a
!        floating-point number, which is unlikely, or some bits are
!        not used in the representation of numbers, which is possible,
!        (e.g. Cray machines) or the mantissa has an implicit bit,
!        (e.g. IEEE machines, Dec Vax machines), which is perhaps the
!        most likely. We have to assume the last alternative.
!        If this is true, then we need to reduce EMAX by one because
!        there must be some way of representing zero in an implicit-bit
!        system. On machines like Cray, we are reducing EMAX by one
!        unnecessarily.
!
         IF ( (MOD(nbits,2)==1) .AND. (Beta==2) ) Emax = Emax - 1
!
!
!        Assume we are on an IEEE machine which reserves one exponent
!        for infinity and NaN.
!
         IF ( Ieee ) Emax = Emax - 1
!
!     Now create RMAX, the largest machine number, which should
!     be equal to (1.0 - BETA**(-P)) * BETA**EMAX .
!
!     First compute 1.0 - BETA**(-P), being careful that the
!     result is less than 1.0 .
!
         recbas = ONE/Beta
         z = Beta - ONE
         y = ZERO
         DO i = 1 , P
            z = z*recbas
            IF ( y<ONE ) oldy = y
            y = DLAMC3(y,z)
         ENDDO
         IF ( y>=ONE ) y = oldy
!
!     Now multiply by BETA**EMAX to get RMAX.
!
         DO i = 1 , Emax
            y = DLAMC3(y*Beta,ZERO)
         ENDDO
!
         Rmax = y
         EXIT
      ENDDO
!
!     End of DLAMC5
!
      END SUBROUTINE DLAMC5
