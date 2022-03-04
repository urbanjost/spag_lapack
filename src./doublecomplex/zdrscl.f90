!*==zdrscl.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b ZDRSCL multiplies a vector by the reciprocal of a real scalar.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZDRSCL + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zdrscl.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zdrscl.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zdrscl.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZDRSCL( N, SA, SX, INCX )
!
!       .. Scalar Arguments ..
!       INTEGER            INCX, N
!       DOUBLE PRECISION   SA
!       ..
!       .. Array Arguments ..
!       COMPLEX*16         SX( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZDRSCL multiplies an n-element complex vector x by the real scalar
!> 1/a.  This is done without overflow or underflow as long as
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
!>          SA is DOUBLE PRECISION
!>          The scalar a which is used to divide each component of x.
!>          SA must be >= 0, or the subroutine will divide by zero.
!> \endverbatim
!>
!> \param[in,out] SX
!> \verbatim
!>          SX is COMPLEX*16 array, dimension
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
!> \ingroup complex16OTHERauxiliary
!
!  =====================================================================
      SUBROUTINE ZDRSCL(N,Sa,Sx,Incx)
      USE F77KINDS                        
      USE S_DLABAD
      USE S_DLAMCH
      USE S_ZDSCAL
      IMPLICIT NONE
!*--ZDRSCL92
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER :: N
      REAL(R8KIND) , INTENT(IN) :: Sa
      COMPLEX(CX16KIND) , DIMENSION(*) :: Sx
      INTEGER :: Incx
!
! Local variable declarations rewritten by SPAG
!
      REAL(R8KIND) :: bignum , cden , cden1 , cnum , cnum1 , mul ,      &
     &                smlnum
      LOGICAL :: done
!
! End of declarations rewritten by SPAG
!
!     ..
!     .. Array Arguments ..
!     ..
!
! =====================================================================
!
!     .. Parameters ..
!     ..
!     .. Local Scalars ..
!     ..
!     .. External Functions ..
!     ..
!     .. External Subroutines ..
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. Executable Statements ..
!
!     Quick return if possible
!
      IF ( N<=0 ) RETURN
!
!     Get machine parameters
!
      smlnum = DLAMCH('S')
      bignum = ONE/smlnum
      CALL DLABAD(smlnum,bignum)
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
         CALL ZDSCAL(N,mul,Sx,Incx)
!
         IF ( done ) EXIT
      ENDDO
!
!
!     End of ZDRSCL
!
      END SUBROUTINE ZDRSCL
