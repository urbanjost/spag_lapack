!*==dlamch.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b DLAMCH
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
!> DLAMCH determines double precision machine parameters.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] CMACH
!> \verbatim
!>          CMACH is CHARACTER*1
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
!> \date December 2016
!
!> \ingroup auxOTHERauxiliary
!
!  =====================================================================
      DOUBLE PRECISION FUNCTION DLAMCH(Cmach)
      IMPLICIT NONE
!*--DLAMCH72
!
!  -- LAPACK auxiliary routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Cmach
!     ..
!
! =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ONE , ZERO
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
!     ..
!     .. Local Scalars ..
      DOUBLE PRECISION rnd , eps , sfmin , small , rmach
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC DIGITS , EPSILON , HUGE , MAXEXPONENT , MINEXPONENT ,   &
     &          RADIX , TINY
!     ..
!     .. Executable Statements ..
!
!
!     Assume rounding, not chopping. Always.
!
      rnd = ONE
!
      IF ( ONE==rnd ) THEN
         eps = EPSILON(ZERO)*0.5
      ELSE
         eps = EPSILON(ZERO)
      ENDIF
!
      IF ( LSAME(Cmach,'E') ) THEN
         rmach = eps
      ELSEIF ( LSAME(Cmach,'S') ) THEN
         sfmin = TINY(ZERO)
         small = ONE/HUGE(ZERO)
!
!           Use SMALL plus a bit, to avoid the possibility of rounding
!           causing overflow when computing  1/sfmin.
!
         IF ( small>=sfmin ) sfmin = small*(ONE+eps)
         rmach = sfmin
      ELSEIF ( LSAME(Cmach,'B') ) THEN
         rmach = RADIX(ZERO)
      ELSEIF ( LSAME(Cmach,'P') ) THEN
         rmach = eps*RADIX(ZERO)
      ELSEIF ( LSAME(Cmach,'N') ) THEN
         rmach = DIGITS(ZERO)
      ELSEIF ( LSAME(Cmach,'R') ) THEN
         rmach = rnd
      ELSEIF ( LSAME(Cmach,'M') ) THEN
         rmach = MINEXPONENT(ZERO)
      ELSEIF ( LSAME(Cmach,'U') ) THEN
         rmach = TINY(ZERO)
      ELSEIF ( LSAME(Cmach,'L') ) THEN
         rmach = MAXEXPONENT(ZERO)
      ELSEIF ( LSAME(Cmach,'O') ) THEN
         rmach = HUGE(ZERO)
      ELSE
         rmach = ZERO
      ENDIF
!
      DLAMCH = rmach
!
!     End of DLAMCH
!
      END FUNCTION DLAMCH
!*==dlamc3.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!***********************************************************************
!> \brief \b DLAMC3
!> \details
!> \b Purpose:
!> \verbatim
!> DLAMC3  is intended to force  A  and  B  to be stored prior to doing
!> the addition of  A  and  B ,  for use in situations where optimizers
!> might hold one of these in a register.
!> \endverbatim
!> \author LAPACK is a software package provided by Univ. of Tennessee, Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..
!> \date December 2016
!> \ingroup auxOTHERauxiliary
!>
!> \param[in] A
!> \verbatim
!>          A is a DOUBLE PRECISION
!> \endverbatim
!>
!> \param[in] B
!> \verbatim
!>          B is a DOUBLE PRECISION
!>          The values A and B.
!> \endverbatim
!>
      DOUBLE PRECISION FUNCTION DLAMC3(A,B)
      IMPLICIT NONE
!*--DLAMC3176
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
