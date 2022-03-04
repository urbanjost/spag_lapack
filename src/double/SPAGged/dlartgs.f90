!*==dlartgs.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b DLARTGS generates a plane rotation designed to introduce a bulge in implicit QR iteration for the bidiagonal SVD problem.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DLARTGS + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlartgs.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlartgs.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlartgs.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLARTGS( X, Y, SIGMA, CS, SN )
!
!       .. Scalar Arguments ..
!       DOUBLE PRECISION        CS, SIGMA, SN, X, Y
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLARTGS generates a plane rotation designed to introduce a bulge in
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
!>          X is DOUBLE PRECISION
!>          The (1,1) entry of an upper bidiagonal matrix.
!> \endverbatim
!>
!> \param[in] Y
!> \verbatim
!>          Y is DOUBLE PRECISION
!>          The (1,2) entry of an upper bidiagonal matrix.
!> \endverbatim
!>
!> \param[in] SIGMA
!> \verbatim
!>          SIGMA is DOUBLE PRECISION
!>          The shift.
!> \endverbatim
!>
!> \param[out] CS
!> \verbatim
!>          CS is DOUBLE PRECISION
!>          The cosine of the rotation.
!> \endverbatim
!>
!> \param[out] SN
!> \verbatim
!>          SN is DOUBLE PRECISION
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
      SUBROUTINE DLARTGS(X,Y,Sigma,Cs,Sn)
      USE F77KINDS                        
      USE S_DLAMCH
      USE S_DLARTGP
      IMPLICIT NONE
!*--DLARTGS97
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  NEGONE = -1.0D0 , ONE = 1.0D0 ,     &
     &                              ZERO = 0.0D0
!
! Dummy argument declarations rewritten by SPAG
!
      REAL(R8KIND) , INTENT(IN) :: X
      REAL(R8KIND) , INTENT(IN) :: Y
      REAL(R8KIND) , INTENT(IN) :: Sigma
      REAL(R8KIND) :: Cs
      REAL(R8KIND) :: Sn
!
! Local variable declarations rewritten by SPAG
!
      REAL(R8KIND) :: r , s , thresh , w , z
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
      thresh = DLAMCH('E')
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
!     CALL DLARTGP( Z, W, CS, SN, R ) might seem more natural;
!     reordering the arguments ensures that if Z = 0 then the rotation
!     is by PI/2.
!
      CALL DLARTGP(w,z,Sn,Cs,r)
!
!
!     End DLARTGS
!
      END SUBROUTINE DLARTGS
