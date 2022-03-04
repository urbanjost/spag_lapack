!*==slaran.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b SLARAN
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       REAL FUNCTION SLARAN( ISEED )
!
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
!> SLARAN returns a random real number from a uniform (0,1)
!> distribution.
!> \endverbatim
!
!  Arguments:
!  ==========
!
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
!>  This routine uses a multiplicative congruential method with modulus
!>  2**48 and multiplier 33952834046453 (see G.S.Fishman,
!>  'Multiplicative congruential random number generators with modulus
!>  2**b: an exhaustive analysis for b = 32 and a partial analysis for
!>  b = 48', Math. Comp. 189, pp 331-344, 1990).
!>
!>  48-bit integers are stored in 4 integer array elements with 12 bits
!>  per element. Hence the routine is portable across machines with
!>  integers of 32 bits or more.
!> \endverbatim
!>
!  =====================================================================
      REAL FUNCTION SLARAN(Iseed)
      IMPLICIT NONE
!*--SLARAN71
!
!  -- LAPACK auxiliary routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Array Arguments ..
      INTEGER Iseed(4)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      INTEGER M1 , M2 , M3 , M4
      PARAMETER (M1=494,M2=322,M3=2508,M4=2549)
      REAL ONE
      PARAMETER (ONE=1.0E+0)
      INTEGER IPW2
      REAL R
      PARAMETER (IPW2=4096,R=ONE/IPW2)
!     ..
!     .. Local Scalars ..
      INTEGER it1 , it2 , it3 , it4
      REAL rndout
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MOD , REAL
      DO
!     ..
!     .. Executable Statements ..
!
!     multiply the seed by the multiplier modulo 2**48
!
         it4 = Iseed(4)*M4
         it3 = it4/IPW2
         it4 = it4 - IPW2*it3
         it3 = it3 + Iseed(3)*M4 + Iseed(4)*M3
         it2 = it3/IPW2
         it3 = it3 - IPW2*it2
         it2 = it2 + Iseed(2)*M4 + Iseed(3)*M3 + Iseed(4)*M2
         it1 = it2/IPW2
         it2 = it2 - IPW2*it1
         it1 = it1 + Iseed(1)*M4 + Iseed(2)*M3 + Iseed(3)*M2 + Iseed(4) &
     &         *M1
         it1 = MOD(it1,IPW2)
!
!     return updated seed
!
         Iseed(1) = it1
         Iseed(2) = it2
         Iseed(3) = it3
         Iseed(4) = it4
!
!     convert 48-bit integer to a real number in the interval (0,1)
!
         rndout = R*(REAL(it1)+R*(REAL(it2)+R*(REAL(it3)+R*(REAL(it4))))&
     &            )
!
!        If a real number has n bits of precision, and the first
!        n bits of the 48-bit integer above happen to be all 1 (which
!        will occur about once every 2**n calls), then SLARAN will
!        be rounded to exactly 1.0. In IEEE single precision arithmetic,
!        this will happen relatively often since n = 24.
!        Since SLARAN is not supposed to return exactly 0.0 or 1.0
!        (and some callers of SLARAN, such as CLARND, depend on that),
!        the statistically correct thing to do in this situation is
!        simply to iterate again.
!        N.B. the case SLARAN = 0.0 should not be possible.
!
         IF ( rndout==1.0 ) CYCLE
!
         SLARAN = rndout
         EXIT
      ENDDO
!
!     End of SLARAN
!
      END FUNCTION SLARAN
