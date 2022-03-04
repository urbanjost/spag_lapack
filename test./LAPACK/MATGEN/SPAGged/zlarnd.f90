!*==zlarnd.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b ZLARND
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       COMPLEX*16   FUNCTION ZLARND( IDIST, ISEED )
!
!       .. Scalar Arguments ..
!       INTEGER            IDIST
!       ..
!       .. Array Arguments ..
!       INTEGER            ISEED( 4 )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZLARND returns a random complex number from a uniform or normal
!> distribution.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] IDIST
!> \verbatim
!>          IDIST is INTEGER
!>          Specifies the distribution of the random numbers:
!>          = 1:  real and imaginary parts each uniform (0,1)
!>          = 2:  real and imaginary parts each uniform (-1,1)
!>          = 3:  real and imaginary parts each normal (0,1)
!>          = 4:  uniformly distributed on the disc abs(z) <= 1
!>          = 5:  uniformly distributed on the circle abs(z) = 1
!> \endverbatim
!>
!> \param[in,out] ISEED
!> \verbatim
!>          ISEED is INTEGER array, dimension (4)
!>          On entry, the seed of the random number generator; the array
!>          elements must be between 0 and 4095, and ISEED(4) must be
!>          odd.
!>          On exit, the seed is updated.
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
!> \ingroup complex16_matgen
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  This routine calls the auxiliary routine DLARAN to generate a random
!>  real number from a uniform (0,1) distribution. The Box-Muller method
!>  is used to transform numbers from a uniform to a normal distribution.
!> \endverbatim
!>
!  =====================================================================
      COMPLEX*16 FUNCTION ZLARND(Idist,Iseed)
      IMPLICIT NONE
!*--ZLARND79
!
!  -- LAPACK auxiliary routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Idist
!     ..
!     .. Array Arguments ..
      INTEGER Iseed(4)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO , ONE , TWO
      PARAMETER (ZERO=0.0D+0,ONE=1.0D+0,TWO=2.0D+0)
      DOUBLE PRECISION TWOPI
      PARAMETER (TWOPI=6.28318530717958647692528676655900576839D+0)
!     ..
!     .. Local Scalars ..
      DOUBLE PRECISION t1 , t2
!     ..
!     .. External Functions ..
      DOUBLE PRECISION DLARAN
      EXTERNAL DLARAN
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC DCMPLX , EXP , LOG , SQRT
!     ..
!     .. Executable Statements ..
!
!     Generate a pair of real random numbers from a uniform (0,1)
!     distribution
!
      t1 = DLARAN(Iseed)
      t2 = DLARAN(Iseed)
!
      IF ( Idist==1 ) THEN
!
!        real and imaginary parts each uniform (0,1)
!
         ZLARND = DCMPLX(t1,t2)
      ELSEIF ( Idist==2 ) THEN
!
!        real and imaginary parts each uniform (-1,1)
!
         ZLARND = DCMPLX(TWO*t1-ONE,TWO*t2-ONE)
      ELSEIF ( Idist==3 ) THEN
!
!        real and imaginary parts each normal (0,1)
!
         ZLARND = SQRT(-TWO*LOG(t1))*EXP(DCMPLX(ZERO,TWOPI*t2))
      ELSEIF ( Idist==4 ) THEN
!
!        uniform distribution on the unit disc abs(z) <= 1
!
         ZLARND = SQRT(t1)*EXP(DCMPLX(ZERO,TWOPI*t2))
      ELSEIF ( Idist==5 ) THEN
!
!        uniform distribution on the unit circle abs(z) = 1
!
         ZLARND = EXP(DCMPLX(ZERO,TWOPI*t2))
      ENDIF
!
!     End of ZLARND
!
      END FUNCTION ZLARND
