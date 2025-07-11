!*==srscl.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b SRSCL multiplies a vector by the reciprocal of a real scalar.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SRSCL + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/srscl.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/srscl.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/srscl.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SRSCL( N, SA, SX, INCX )
!
!       .. Scalar Arguments ..
!       INTEGER            INCX, N
!       REAL               SA
!       ..
!       .. Array Arguments ..
!       REAL               SX( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SRSCL multiplies an n-element real vector x by the real scalar 1/a.
!> This is done without overflow or underflow as long as
!> the final result x/a does not overflow or underflow.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of components of the vector x.
!> \endverbatim
!>
!> \param[in] SA
!> \verbatim
!>          SA is REAL
!>          The scalar a which is used to divide each component of x.
!>          SA must be >= 0, or the subroutine will divide by zero.
!> \endverbatim
!>
!> \param[in,out] SX
!> \verbatim
!>          SX is REAL array, dimension
!>                         (1+(N-1)*abs(INCX))
!>          The n-element vector x.
!> \endverbatim
!>
!> \param[in] INCX
!> \verbatim
!>          INCX is INTEGER
!>          The increment between successive values of the vector SX.
!>          > 0:  SX(1) = X(1) and SX(1+(i-1)*INCX) = x(i),     1< i<= n
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
!> \ingroup realOTHERauxiliary
!
!  =====================================================================
      SUBROUTINE SRSCL(N,Sa,Sx,Incx)
      IMPLICIT NONE
!*--SRSCL88
!
!  -- LAPACK auxiliary routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Incx , N
      REAL Sa
!     ..
!     .. Array Arguments ..
      REAL Sx(*)
!     ..
!
! =====================================================================
!
!     .. Parameters ..
      REAL ONE , ZERO
      PARAMETER (ONE=1.0E+0,ZERO=0.0E+0)
!     ..
!     .. Local Scalars ..
      LOGICAL done
      REAL bignum , cden , cden1 , cnum , cnum1 , mul , smlnum
!     ..
!     .. External Functions ..
      REAL SLAMCH
      EXTERNAL SLAMCH
!     ..
!     .. External Subroutines ..
      EXTERNAL SLABAD , SSCAL
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS
!     ..
!     .. Executable Statements ..
!
!     Quick return if possible
!
      IF ( N<=0 ) RETURN
!
!     Get machine parameters
!
      smlnum = SLAMCH('S')
      bignum = ONE/smlnum
      CALL SLABAD(smlnum,bignum)
!
!     Initialize the denominator to SA and the numerator to 1.
!
      cden = Sa
      cnum = ONE
      DO
!
         cden1 = cden*smlnum
         cnum1 = cnum/bignum
         IF ( ABS(cden1)>ABS(cnum) .AND. cnum/=ZERO ) THEN
!
!        Pre-multiply X by SMLNUM if CDEN is large compared to CNUM.
!
            mul = smlnum
            done = .FALSE.
            cden = cden1
         ELSEIF ( ABS(cnum1)>ABS(cden) ) THEN
!
!        Pre-multiply X by BIGNUM if CDEN is small compared to CNUM.
!
            mul = bignum
            done = .FALSE.
            cnum = cnum1
         ELSE
!
!        Multiply X by CNUM / CDEN and return.
!
            mul = cnum/cden
            done = .TRUE.
         ENDIF
!
!     Scale the vector X by MUL
!
         CALL SSCAL(N,mul,Sx,Incx)
!
         IF ( done ) EXIT
      ENDDO
!
!
!     End of SRSCL
!
      END SUBROUTINE SRSCL
