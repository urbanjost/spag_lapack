!*==clarnv.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b CLARNV returns a vector of random numbers from a uniform or normal distribution.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CLARNV + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clarnv.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clarnv.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clarnv.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CLARNV( IDIST, ISEED, N, X )
!
!       .. Scalar Arguments ..
!       INTEGER            IDIST, N
!       ..
!       .. Array Arguments ..
!       INTEGER            ISEED( 4 )
!       COMPLEX            X( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CLARNV returns a vector of n random complex numbers from a uniform or
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
!>          = 1:  real and imaginary parts each uniform (0,1)
!>          = 2:  real and imaginary parts each uniform (-1,1)
!>          = 3:  real and imaginary parts each normal (0,1)
!>          = 4:  uniformly distributed on the disc abs(z) < 1
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
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of random numbers to be generated.
!> \endverbatim
!>
!> \param[out] X
!> \verbatim
!>          X is COMPLEX array, dimension (N)
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
!> \ingroup complexOTHERauxiliary
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
      SUBROUTINE CLARNV(Idist,Iseed,N,X)
      USE S_SLARUV
      IMPLICIT NONE
!*--CLARNV104
!
! PARAMETER definitions rewritten by SPAG
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0 , TWO = 2.0E+0
      INTEGER , PARAMETER  ::  LV = 128
      REAL , PARAMETER  ::  TWOPI =                                     &
     &                      6.28318530717958647692528676655900576839E+0
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER , INTENT(IN) :: Idist
      INTEGER , DIMENSION(4) :: Iseed
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(OUT) , DIMENSION(*) :: X
!
! Local variable declarations rewritten by SPAG
!
      INTEGER :: i , il , iv
      REAL , DIMENSION(LV) :: u
!
! End of declarations rewritten by SPAG
!
!     ..
!     .. Array Arguments ..
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
!     ..
!     .. Local Scalars ..
!     ..
!     .. Local Arrays ..
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. External Subroutines ..
!     ..
!     .. Executable Statements ..
!
      DO iv = 1 , N , LV/2
         il = MIN(LV/2,N-iv+1)
!
!        Call SLARUV to generate 2*IL real numbers from a uniform (0,1)
!        distribution (2*IL <= LV)
!
         CALL SLARUV(Iseed,2*il,u)
!
         IF ( Idist==1 ) THEN
!
!           Copy generated numbers
!
            DO i = 1 , il
               X(iv+i-1) = CMPLX(u(2*i-1),u(2*i))
            ENDDO
         ELSEIF ( Idist==2 ) THEN
!
!           Convert generated numbers to uniform (-1,1) distribution
!
            DO i = 1 , il
               X(iv+i-1) = CMPLX(TWO*u(2*i-1)-ONE,TWO*u(2*i)-ONE)
            ENDDO
         ELSEIF ( Idist==3 ) THEN
!
!           Convert generated numbers to normal (0,1) distribution
!
            DO i = 1 , il
               X(iv+i-1) = SQRT(-TWO*LOG(u(2*i-1)))                     &
     &                     *EXP(CMPLX(ZERO,TWOPI*u(2*i)))
            ENDDO
         ELSEIF ( Idist==4 ) THEN
!
!           Convert generated numbers to complex numbers uniformly
!           distributed on the unit disk
!
            DO i = 1 , il
               X(iv+i-1) = SQRT(u(2*i-1))*EXP(CMPLX(ZERO,TWOPI*u(2*i)))
            ENDDO
         ELSEIF ( Idist==5 ) THEN
!
!           Convert generated numbers to complex numbers uniformly
!           distributed on the unit circle
!
            DO i = 1 , il
               X(iv+i-1) = EXP(CMPLX(ZERO,TWOPI*u(2*i)))
            ENDDO
         ENDIF
      ENDDO
!
!     End of CLARNV
!
      END SUBROUTINE CLARNV
