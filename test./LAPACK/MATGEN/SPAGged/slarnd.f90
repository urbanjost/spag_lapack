!*==slarnd.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b SLARND
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       REAL             FUNCTION SLARND( IDIST, ISEED )
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
!> SLARND returns a random real number from a uniform or normal
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
!>          = 1:  uniform (0,1)
!>          = 2:  uniform (-1,1)
!>          = 3:  normal (0,1)
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
!> \ingroup real_matgen
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  This routine calls the auxiliary routine SLARAN to generate a random
!>  real number from a uniform (0,1) distribution. The Box-Muller method
!>  is used to transform numbers from a uniform to a normal distribution.
!> \endverbatim
!>
!  =====================================================================
      REAL FUNCTION SLARND(Idist,Iseed)
      IMPLICIT NONE
!*--SLARND77
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
      REAL ONE , TWO
      PARAMETER (ONE=1.0E+0,TWO=2.0E+0)
      REAL TWOPI
      PARAMETER (TWOPI=6.28318530717958647692528676655900576839E+0)
!     ..
!     .. Local Scalars ..
      REAL t1 , t2
!     ..
!     .. External Functions ..
      REAL SLARAN
      EXTERNAL SLARAN
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC COS , LOG , SQRT
!     ..
!     .. Executable Statements ..
!
!     Generate a real random number from a uniform (0,1) distribution
!
      t1 = SLARAN(Iseed)
!
      IF ( Idist==1 ) THEN
!
!        uniform (0,1)
!
         SLARND = t1
      ELSEIF ( Idist==2 ) THEN
!
!        uniform (-1,1)
!
         SLARND = TWO*t1 - ONE
      ELSEIF ( Idist==3 ) THEN
!
!        normal (0,1)
!
         t2 = SLARAN(Iseed)
         SLARND = SQRT(-TWO*LOG(t1))*COS(TWOPI*t2)
      ENDIF
!
!     End of SLARND
!
      END FUNCTION SLARND