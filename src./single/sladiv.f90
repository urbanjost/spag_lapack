!*==sladiv.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b SLADIV performs complex division in real arithmetic, avoiding unnecessary overflow.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SLADIV + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sladiv.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sladiv.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sladiv.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SLADIV( A, B, C, D, P, Q )
!
!       .. Scalar Arguments ..
!       REAL               A, B, C, D, P, Q
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SLADIV performs complex division in  real arithmetic
!>
!>                       a + i*b
!>            p + i*q = ---------
!>                       c + i*d
!>
!> The algorithm is due to Michael Baudin and Robert L. Smith
!> and can be found in the paper
!> "A Robust Complex Division in Scilab"
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] A
!> \verbatim
!>          A is REAL
!> \endverbatim
!>
!> \param[in] B
!> \verbatim
!>          B is REAL
!> \endverbatim
!>
!> \param[in] C
!> \verbatim
!>          C is REAL
!> \endverbatim
!>
!> \param[in] D
!> \verbatim
!>          D is REAL
!>          The scalars a, b, c, and d in the above expression.
!> \endverbatim
!>
!> \param[out] P
!> \verbatim
!>          P is REAL
!> \endverbatim
!>
!> \param[out] Q
!> \verbatim
!>          Q is REAL
!>          The scalars p and q in the above expression.
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
!> \date January 2013
!
!> \ingroup realOTHERauxiliary
!
!  =====================================================================
      SUBROUTINE SLADIV(A,B,C,D,P,Q)
      USE S_SLADIV1
      USE S_SLAMCH
      IMPLICIT NONE
!*--SLADIV97
!
! PARAMETER definitions rewritten by SPAG
!
      REAL , PARAMETER  ::  BS = 2.0E0 , HALF = 0.5E0 , TWO = 2.0E0
!
! Dummy argument declarations rewritten by SPAG
!
      REAL , INTENT(IN) :: A
      REAL , INTENT(IN) :: B
      REAL , INTENT(IN) :: C
      REAL , INTENT(IN) :: D
      REAL , INTENT(INOUT) :: P
      REAL , INTENT(INOUT) :: Q
!
! Local variable declarations rewritten by SPAG
!
      REAL :: aa , ab , bb , be , cc , cd , dd , eps , ov , s , un
!
! End of declarations rewritten by SPAG
!
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
!
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
      aa = A
      bb = B
      cc = C
      dd = D
      ab = MAX(ABS(A),ABS(B))
      cd = MAX(ABS(C),ABS(D))
      s = 1.0E0
 
      ov = SLAMCH('Overflow threshold')
      un = SLAMCH('Safe minimum')
      eps = SLAMCH('Epsilon')
      be = BS/(eps*eps)
 
      IF ( ab>=HALF*ov ) THEN
         aa = HALF*aa
         bb = HALF*bb
         s = TWO*s
      ENDIF
      IF ( cd>=HALF*ov ) THEN
         cc = HALF*cc
         dd = HALF*dd
         s = HALF*s
      ENDIF
      IF ( ab<=un*BS/eps ) THEN
         aa = aa*be
         bb = bb*be
         s = s/be
      ENDIF
      IF ( cd<=un*BS/eps ) THEN
         cc = cc*be
         dd = dd*be
         s = s*be
      ENDIF
      IF ( ABS(D)<=ABS(C) ) THEN
         CALL SLADIV1(aa,bb,cc,dd,P,Q)
      ELSE
         CALL SLADIV1(bb,aa,dd,cc,P,Q)
         Q = -Q
      ENDIF
      P = P*s
      Q = Q*s
!
!
!     End of SLADIV
!
      END SUBROUTINE SLADIV
!*==sladiv1.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
 
!> \ingroup realOTHERauxiliary
 
 
      SUBROUTINE SLADIV1(A,B,C,D,P,Q)
      USE S_SLADIV2
      IMPLICIT NONE
!*--SLADIV1188
!
! PARAMETER definitions rewritten by SPAG
!
      REAL , PARAMETER  ::  ONE = 1.0E0
!
! Dummy argument declarations rewritten by SPAG
!
      REAL , INTENT(INOUT) :: A
      REAL :: B
      REAL :: C
      REAL :: D
      REAL , INTENT(OUT) :: P
      REAL , INTENT(OUT) :: Q
!
! Local variable declarations rewritten by SPAG
!
      REAL :: r , t
!
! End of declarations rewritten by SPAG
!
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
!
!     .. Local Scalars ..
!     ..
!     .. External Functions ..
!     ..
!     .. Executable Statements ..
!
      r = D/C
      t = ONE/(C+D*r)
      P = SLADIV2(A,B,C,D,r,t)
      A = -A
      Q = SLADIV2(B,A,C,D,r,t)
!
!
!     End of SLADIV1
!
      END SUBROUTINE SLADIV1
!*==sladiv2.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
 
!> \ingroup realOTHERauxiliary
 
      FUNCTION SLADIV2(A,B,C,D,R,T)
      IMPLICIT NONE
!*--SLADIV2237
!
! PARAMETER definitions rewritten by SPAG
!
      REAL , PARAMETER  ::  ZERO = 0.0E0
      REAL :: SLADIV2
!
! Dummy argument declarations rewritten by SPAG
!
      REAL , INTENT(IN) :: A
      REAL , INTENT(IN) :: B
      REAL , INTENT(IN) :: C
      REAL , INTENT(IN) :: D
      REAL , INTENT(IN) :: R
      REAL , INTENT(IN) :: T
!
! Local variable declarations rewritten by SPAG
!
      REAL :: br
!
! End of declarations rewritten by SPAG
!
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
!
!     .. Local Scalars ..
!     ..
!     .. Executable Statements ..
!
      IF ( R/=ZERO ) THEN
         br = B*R
         IF ( br/=ZERO ) THEN
            SLADIV2 = (A+br)*T
         ELSE
            SLADIV2 = A*T + (B*T)*R
         ENDIF
      ELSE
         SLADIV2 = (A+D*(B/C))*T
      ENDIF
!
!
!     End of SLADIV
!
      END FUNCTION SLADIV2
