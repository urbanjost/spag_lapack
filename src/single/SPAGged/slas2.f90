!*==slas2.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b SLAS2 computes singular values of a 2-by-2 triangular matrix.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SLAS2 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slas2.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slas2.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slas2.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SLAS2( F, G, H, SSMIN, SSMAX )
!
!       .. Scalar Arguments ..
!       REAL               F, G, H, SSMAX, SSMIN
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SLAS2  computes the singular values of the 2-by-2 matrix
!>    [  F   G  ]
!>    [  0   H  ].
!> On return, SSMIN is the smaller singular value and SSMAX is the
!> larger singular value.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] F
!> \verbatim
!>          F is REAL
!>          The (1,1) element of the 2-by-2 matrix.
!> \endverbatim
!>
!> \param[in] G
!> \verbatim
!>          G is REAL
!>          The (1,2) element of the 2-by-2 matrix.
!> \endverbatim
!>
!> \param[in] H
!> \verbatim
!>          H is REAL
!>          The (2,2) element of the 2-by-2 matrix.
!> \endverbatim
!>
!> \param[out] SSMIN
!> \verbatim
!>          SSMIN is REAL
!>          The smaller singular value.
!> \endverbatim
!>
!> \param[out] SSMAX
!> \verbatim
!>          SSMAX is REAL
!>          The larger singular value.
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
!>  Barring over/underflow, all output quantities are correct to within
!>  a few units in the last place (ulps), even in the absence of a guard
!>  digit in addition/subtraction.
!>
!>  In IEEE arithmetic, the code works correctly if one matrix element is
!>  infinite.
!>
!>  Overflow will not occur unless the largest singular value itself
!>  overflows, or is within a few ulps of overflow. (On machines with
!>  partial overflow, like the Cray, overflow may occur if the largest
!>  singular value is within a factor of 2 of overflow.)
!>
!>  Underflow is harmless if underflow is gradual. Otherwise, results
!>  may correspond to a matrix modified by perturbations of size near
!>  the underflow threshold.
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE SLAS2(F,G,H,Ssmin,Ssmax)
      IMPLICIT NONE
!*--SLAS2111
!
! PARAMETER definitions rewritten by SPAG
!
      REAL , PARAMETER  ::  ZERO = 0.0E0 , ONE = 1.0E0 , TWO = 2.0E0
!
! Dummy argument declarations rewritten by SPAG
!
      REAL , INTENT(IN) :: F
      REAL , INTENT(IN) :: G
      REAL , INTENT(IN) :: H
      REAL , INTENT(INOUT) :: Ssmin
      REAL , INTENT(OUT) :: Ssmax
!
! Local variable declarations rewritten by SPAG
!
      REAL :: as , at , au , c , fa , fhmn , fhmx , ga , ha
!
! End of declarations rewritten by SPAG
!
!     ..
!
!  ====================================================================
!
!     .. Parameters ..
!     ..
!     .. Local Scalars ..
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. Executable Statements ..
!
      fa = ABS(F)
      ga = ABS(G)
      ha = ABS(H)
      fhmn = MIN(fa,ha)
      fhmx = MAX(fa,ha)
      IF ( fhmn==ZERO ) THEN
         Ssmin = ZERO
         IF ( fhmx==ZERO ) THEN
            Ssmax = ga
         ELSE
            Ssmax = MAX(fhmx,ga)                                        &
     &              *SQRT(ONE+(MIN(fhmx,ga)/MAX(fhmx,ga))**2)
         ENDIF
      ELSEIF ( ga<fhmx ) THEN
         as = ONE + fhmn/fhmx
         at = (fhmx-fhmn)/fhmx
         au = (ga/fhmx)**2
         c = TWO/(SQRT(as*as+au)+SQRT(at*at+au))
         Ssmin = fhmn*c
         Ssmax = fhmx/c
      ELSE
         au = fhmx/ga
         IF ( au==ZERO ) THEN
!
!              Avoid possible harmful underflow if exponent range
!              asymmetric (true SSMIN may not underflow even if
!              AU underflows)
!
            Ssmin = (fhmn*fhmx)/ga
            Ssmax = ga
         ELSE
            as = ONE + fhmn/fhmx
            at = (fhmx-fhmn)/fhmx
            c = ONE/(SQRT(ONE+(as*au)**2)+SQRT(ONE+(at*au)**2))
            Ssmin = (fhmn*c)*au
            Ssmin = Ssmin + Ssmin
            Ssmax = ga/(c+c)
         ENDIF
      ENDIF
!
!     End of SLAS2
!
      END SUBROUTINE SLAS2
