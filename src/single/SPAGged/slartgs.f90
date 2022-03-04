!*==slartgs.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b SLARTGS generates a plane rotation designed to introduce a bulge in implicit QR iteration for the bidiagonal SVD problem.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SLARTGS + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slartgs.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slartgs.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slartgs.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SLARTGS( X, Y, SIGMA, CS, SN )
!
!       .. Scalar Arguments ..
!       REAL                    CS, SIGMA, SN, X, Y
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SLARTGS generates a plane rotation designed to introduce a bulge in
!> Golub-Reinsch-style implicit QR iteration for the bidiagonal SVD
!> problem. X and Y are the top-row entries, and SIGMA is the shift.
!> The computed CS and SN define a plane rotation satisfying
!>
!>    [  CS  SN  ]  .  [ X^2 - SIGMA ]  =  [ R ],
!>    [ -SN  CS  ]     [    X * Y    ]     [ 0 ]
!>
!> with R nonnegative.  If X^2 - SIGMA and X * Y are 0, then the
!> rotation is by PI/2.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] X
!> \verbatim
!>          X is REAL
!>          The (1,1) entry of an upper bidiagonal matrix.
!> \endverbatim
!>
!> \param[in] Y
!> \verbatim
!>          Y is REAL
!>          The (1,2) entry of an upper bidiagonal matrix.
!> \endverbatim
!>
!> \param[in] SIGMA
!> \verbatim
!>          SIGMA is REAL
!>          The shift.
!> \endverbatim
!>
!> \param[out] CS
!> \verbatim
!>          CS is REAL
!>          The cosine of the rotation.
!> \endverbatim
!>
!> \param[out] SN
!> \verbatim
!>          SN is REAL
!>          The sine of the rotation.
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
!> \date November 2017
!
!> \ingroup auxOTHERcomputational
!
!  =====================================================================
      SUBROUTINE SLARTGS(X,Y,Sigma,Cs,Sn)
      USE S_SLAMCH
      USE S_SLARTGP
      IMPLICIT NONE
!*--SLARTGS96
!
! PARAMETER definitions rewritten by SPAG
!
      REAL , PARAMETER  ::  NEGONE = -1.0E0 , ONE = 1.0E0 , ZERO = 0.0E0
!
! Dummy argument declarations rewritten by SPAG
!
      REAL , INTENT(IN) :: X
      REAL , INTENT(IN) :: Y
      REAL , INTENT(IN) :: Sigma
      REAL :: Cs
      REAL :: Sn
!
! Local variable declarations rewritten by SPAG
!
      REAL :: r , s , thresh , w , z
!
! End of declarations rewritten by SPAG
!
!     ..
!
!  ===================================================================
!
!     .. Parameters ..
!     ..
!     .. Local Scalars ..
!     ..
!     .. External Subroutines ..
!     ..
!     .. External Functions ..
!     .. Executable Statements ..
!
      thresh = SLAMCH('E')
!
!     Compute the first column of B**T*B - SIGMA^2*I, up to a scale
!     factor.
!
      IF ( (Sigma==ZERO .AND. ABS(X)<thresh) .OR.                       &
     &     (ABS(X)==Sigma .AND. Y==ZERO) ) THEN
         z = ZERO
         w = ZERO
      ELSEIF ( Sigma==ZERO ) THEN
         IF ( X>=ZERO ) THEN
            z = X
            w = Y
         ELSE
            z = -X
            w = -Y
         ENDIF
      ELSEIF ( ABS(X)<thresh ) THEN
         z = -Sigma*Sigma
         w = ZERO
      ELSE
         IF ( X>=ZERO ) THEN
            s = ONE
         ELSE
            s = NEGONE
         ENDIF
         z = s*(ABS(X)-Sigma)*(s+Sigma/X)
         w = s*Y
      ENDIF
!
!     Generate the rotation.
!     CALL SLARTGP( Z, W, CS, SN, R ) might seem more natural;
!     reordering the arguments ensures that if Z = 0 then the rotation
!     is by PI/2.
!
      CALL SLARTGP(w,z,Sn,Cs,r)
!
!
!     End SLARTGS
!
      END SUBROUTINE SLARTGS
