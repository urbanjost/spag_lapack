!*==slae2.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b SLAE2 computes the eigenvalues of a 2-by-2 symmetric matrix.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SLAE2 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slae2.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slae2.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slae2.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SLAE2( A, B, C, RT1, RT2 )
!
!       .. Scalar Arguments ..
!       REAL               A, B, C, RT1, RT2
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SLAE2  computes the eigenvalues of a 2-by-2 symmetric matrix
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
!>          A is REAL
!>          The (1,1) element of the 2-by-2 matrix.
!> \endverbatim
!>
!> \param[in] B
!> \verbatim
!>          B is REAL
!>          The (1,2) and (2,1) elements of the 2-by-2 matrix.
!> \endverbatim
!>
!> \param[in] C
!> \verbatim
!>          C is REAL
!>          The (2,2) element of the 2-by-2 matrix.
!> \endverbatim
!>
!> \param[out] RT1
!> \verbatim
!>          RT1 is REAL
!>          The eigenvalue of larger absolute value.
!> \endverbatim
!>
!> \param[out] RT2
!> \verbatim
!>          RT2 is REAL
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
      SUBROUTINE SLAE2(A,B,C,Rt1,Rt2)
      IMPLICIT NONE
!*--SLAE2106
!
! PARAMETER definitions rewritten by SPAG
!
      REAL , PARAMETER  ::  ONE = 1.0E0 , TWO = 2.0E0 , ZERO = 0.0E0 ,  &
     &                      HALF = 0.5E0
!
! Dummy argument declarations rewritten by SPAG
!
      REAL , INTENT(IN) :: A
      REAL , INTENT(IN) :: B
      REAL , INTENT(IN) :: C
      REAL , INTENT(INOUT) :: Rt1
      REAL , INTENT(OUT) :: Rt2
!
! Local variable declarations rewritten by SPAG
!
      REAL :: ab , acmn , acmx , adf , df , rt , sm , tb
!
! End of declarations rewritten by SPAG
!
!     ..
!
! =====================================================================
!
!     .. Parameters ..
!     ..
!     .. Local Scalars ..
!     ..
!     .. Intrinsic Functions ..
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
!     End of SLAE2
!
      END SUBROUTINE SLAE2
