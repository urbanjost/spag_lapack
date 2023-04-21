PROGRAM DLAMCHTST
!*==aa0001.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
      IMPLICIT NONE
!*--AA00013
!> \brief \b DLAMCHTST
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!      PROGRAM DLAMCHTST
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
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
! =====================================================================
!
!     .. Local Scalars ..
      DOUBLE PRECISION base , emax , emin , eps , prec , rmax , rmin ,  &
     &                 rnd , sfmin , t
!     ..
!     .. External Functions ..
      DOUBLE PRECISION DLAMCH
      EXTERNAL DLAMCH
!     ..
!     .. Executable Statements ..
!
      eps = DLAMCH('Epsilon')
      sfmin = DLAMCH('Safe minimum')
      base = DLAMCH('Base')
      prec = DLAMCH('Precision')
      t = DLAMCH('Number of digits in mantissa')
      rnd = DLAMCH('Rounding mode')
      emin = DLAMCH('Minimum exponent')
      rmin = DLAMCH('Underflow threshold')
      emax = DLAMCH('Largest exponent')
      rmax = DLAMCH('Overflow threshold')
!
      WRITE (6,*) ' Epsilon                      = ' , eps
      WRITE (6,*) ' Safe minimum                 = ' , sfmin
      WRITE (6,*) ' Base                         = ' , base
      WRITE (6,*) ' Precision                    = ' , prec
      WRITE (6,*) ' Number of digits in mantissa = ' , t
      WRITE (6,*) ' Rounding mode                = ' , rnd
      WRITE (6,*) ' Minimum exponent             = ' , emin
      WRITE (6,*) ' Underflow threshold          = ' , rmin
      WRITE (6,*) ' Largest exponent             = ' , emax
      WRITE (6,*) ' Overflow threshold           = ' , rmax
      WRITE (6,*) ' Reciprocal of safe minimum   = ' , 1/sfmin
!
      END PROGRAM DLAMCHTST
