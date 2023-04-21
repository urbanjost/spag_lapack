!*==dlae2.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b DLAE2 computes the eigenvalues of a 2-by-2 symmetric matrix.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DLAE2 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlae2.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlae2.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlae2.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLAE2( A, B, C, RT1, RT2 )
!
!       .. Scalar Arguments ..
!       DOUBLE PRECISION   A, B, C, RT1, RT2
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLAE2  computes the eigenvalues of a 2-by-2 symmetric matrix
!>    [  A   B  ]
!>    [  B   C  ].
!> On return, RT1 is the eigenvalue of larger absolute value, and RT2
!> is the eigenvalue of smaller absolute value.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] A
!> \verbatim
!>          A is DOUBLE PRECISION
!>          The (1,1) element of the 2-by-2 matrix.
!> \endverbatim
!>
!> \param[in] B
!> \verbatim
!>          B is DOUBLE PRECISION
!>          The (1,2) and (2,1) elements of the 2-by-2 matrix.
!> \endverbatim
!>
!> \param[in] C
!> \verbatim
!>          C is DOUBLE PRECISION
!>          The (2,2) element of the 2-by-2 matrix.
!> \endverbatim
!>
!> \param[out] RT1
!> \verbatim
!>          RT1 is DOUBLE PRECISION
!>          The eigenvalue of larger absolute value.
!> \endverbatim
!>
!> \param[out] RT2
!> \verbatim
!>          RT2 is DOUBLE PRECISION
!>          The eigenvalue of smaller absolute value.
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
!>  RT1 is accurate to a few ulps barring over/underflow.
!>
!>  RT2 may be inaccurate if there is massive cancellation in the
!>  determinant A*C-B*B; higher precision or correctly rounded or
!>  correctly truncated arithmetic would be needed to compute RT2
!>  accurately in all cases.
!>
!>  Overflow is possible only if RT1 is within a factor of 5 of overflow.
!>  Underflow is harmless if the input data is 0 or exceeds
!>     underflow_threshold / macheps.
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE DLAE2(A,B,C,Rt1,Rt2)
      IMPLICIT NONE
!*--DLAE2106
!
!  -- LAPACK auxiliary routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      DOUBLE PRECISION A , B , C , Rt1 , Rt2
!     ..
!
! =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ONE
      PARAMETER (ONE=1.0D0)
      DOUBLE PRECISION TWO
      PARAMETER (TWO=2.0D0)
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
      DOUBLE PRECISION HALF
      PARAMETER (HALF=0.5D0)
!     ..
!     .. Local Scalars ..
      DOUBLE PRECISION ab , acmn , acmx , adf , df , rt , sm , tb
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , SQRT
!     ..
!     .. Executable Statements ..
!
!     Compute the eigenvalues
!
      sm = A + C
      df = A - C
      adf = ABS(df)
      tb = B + B
      ab = ABS(tb)
      IF ( ABS(A)>ABS(C) ) THEN
         acmx = A
         acmn = C
      ELSE
         acmx = C
         acmn = A
      ENDIF
      IF ( adf>ab ) THEN
         rt = adf*SQRT(ONE+(ab/adf)**2)
      ELSEIF ( adf<ab ) THEN
         rt = ab*SQRT(ONE+(adf/ab)**2)
      ELSE
!
!        Includes case AB=ADF=0
!
         rt = ab*SQRT(TWO)
      ENDIF
      IF ( sm<ZERO ) THEN
         Rt1 = HALF*(sm-rt)
!
!        Order of execution important.
!        To get fully accurate smaller eigenvalue,
!        next line needs to be executed in higher precision.
!
         Rt2 = (acmx/Rt1)*acmn - (B/Rt1)*B
      ELSEIF ( sm>ZERO ) THEN
         Rt1 = HALF*(sm+rt)
!
!        Order of execution important.
!        To get fully accurate smaller eigenvalue,
!        next line needs to be executed in higher precision.
!
         Rt2 = (acmx/Rt1)*acmn - (B/Rt1)*B
      ELSE
!
!        Includes case RT1 = RT2 = 0
!
         Rt1 = HALF*rt
         Rt2 = -HALF*rt
      ENDIF
!
!     End of DLAE2
!
      END SUBROUTINE DLAE2
