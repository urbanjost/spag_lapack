!*==slarnv.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b SLARNV returns a vector of random numbers from a uniform or normal distribution.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SLARNV + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slarnv.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slarnv.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slarnv.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SLARNV( IDIST, ISEED, N, X )
!
!       .. Scalar Arguments ..
!       INTEGER            IDIST, N
!       ..
!       .. Array Arguments ..
!       INTEGER            ISEED( 4 )
!       REAL               X( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SLARNV returns a vector of n random real numbers from a uniform or
!> normal distribution.
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
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of random numbers to be generated.
!> \endverbatim
!>
!> \param[out] X
!> \verbatim
!>          X is REAL array, dimension (N)
!>          The generated random numbers.
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
!> \ingroup OTHERauxiliary
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  This routine calls the auxiliary routine SLARUV to generate random
!>  real numbers from a uniform (0,1) distribution, in batches of up to
!>  128 using vectorisable code. The Box-Muller method is used to
!>  transform numbers from a uniform to a normal distribution.
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE SLARNV(Idist,Iseed,N,X)
      IMPLICIT NONE
!*--SLARNV101
!
!  -- LAPACK auxiliary routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Idist , N
!     ..
!     .. Array Arguments ..
      INTEGER Iseed(4)
      REAL X(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ONE , TWO
      PARAMETER (ONE=1.0E+0,TWO=2.0E+0)
      INTEGER LV
      PARAMETER (LV=128)
      REAL TWOPI
      PARAMETER (TWOPI=6.28318530717958647692528676655900576839E+0)
!     ..
!     .. Local Scalars ..
      INTEGER i , il , il2 , iv
!     ..
!     .. Local Arrays ..
      REAL u(LV)
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC COS , LOG , MIN , SQRT
!     ..
!     .. External Subroutines ..
      EXTERNAL SLARUV
!     ..
!     .. Executable Statements ..
!
      DO iv = 1 , N , LV/2
         il = MIN(LV/2,N-iv+1)
         IF ( Idist==3 ) THEN
            il2 = 2*il
         ELSE
            il2 = il
         ENDIF
!
!        Call SLARUV to generate IL2 numbers from a uniform (0,1)
!        distribution (IL2 <= LV)
!
         CALL SLARUV(Iseed,il2,u)
!
         IF ( Idist==1 ) THEN
!
!           Copy generated numbers
!
            DO i = 1 , il
               X(iv+i-1) = u(i)
            ENDDO
         ELSEIF ( Idist==2 ) THEN
!
!           Convert generated numbers to uniform (-1,1) distribution
!
            DO i = 1 , il
               X(iv+i-1) = TWO*u(i) - ONE
            ENDDO
         ELSEIF ( Idist==3 ) THEN
!
!           Convert generated numbers to normal (0,1) distribution
!
            DO i = 1 , il
               X(iv+i-1) = SQRT(-TWO*LOG(u(2*i-1)))*COS(TWOPI*u(2*i))
            ENDDO
         ENDIF
      ENDDO
!
!     End of SLARNV
!
      END SUBROUTINE SLARNV
