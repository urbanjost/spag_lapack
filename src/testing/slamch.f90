!*==slamch.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b SLAMCH
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!      REAL             FUNCTION SLAMCH( CMACH )
!
!     .. Scalar Arguments ..
!      CHARACTER          CMACH
!     ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SLAMCH determines single precision machine parameters.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] CMACH
!> \verbatim
!>          CMACH is CHARACTER*1
!>          Specifies the value to be returned by SLAMCH:
!>          = 'E' or 'e',   SLAMCH := eps
!>          = 'S' or 's ,   SLAMCH := sfmin
!>          = 'B' or 'b',   SLAMCH := base
!>          = 'P' or 'p',   SLAMCH := eps*base
!>          = 'N' or 'n',   SLAMCH := t
!>          = 'R' or 'r',   SLAMCH := rnd
!>          = 'M' or 'm',   SLAMCH := emin
!>          = 'U' or 'u',   SLAMCH := rmin
!>          = 'L' or 'l',   SLAMCH := emax
!>          = 'O' or 'o',   SLAMCH := rmax
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
      REAL FUNCTION SLAMCH(Cmach)
      IMPLICIT NONE
!*--SLAMCH72
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
      REAL ONE , ZERO
      PARAMETER (ONE=1.0E+0,ZERO=0.0E+0)
!     ..
!     .. Local Scalars ..
      REAL rnd , eps , sfmin , small , rmach
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
      SLAMCH = rmach
!
!     End of SLAMCH
!
      END FUNCTION SLAMCH
!*==slamc3.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!***********************************************************************
!> \brief \b SLAMC3
!> \details
!> \b Purpose:
!> \verbatim
!> SLAMC3  is intended to force  A  and  B  to be stored prior to doing
!> the addition of  A  and  B ,  for use in situations where optimizers
!> might hold one of these in a register.
!> \endverbatim
!> \author LAPACK is a software package provided by Univ. of Tennessee, Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..
!> \date December 2016
!> \ingroup auxOTHERauxiliary
!>
!> \param[in] A
!> \verbatim
!> \endverbatim
!>
!> \param[in] B
!> \verbatim
!>          The values A and B.
!> \endverbatim
!>
!
      REAL FUNCTION SLAMC3(A,B)
      IMPLICIT NONE
!*--SLAMC3175
!
!  -- LAPACK auxiliary routine (version 3.7.0) --
!     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
!     November 2010
!
!     .. Scalar Arguments ..
      REAL A , B
!     ..
! =====================================================================
!
!     .. Executable Statements ..
!
      SLAMC3 = A + B
!
!
!     End of SLAMC3
!
      END FUNCTION SLAMC3
