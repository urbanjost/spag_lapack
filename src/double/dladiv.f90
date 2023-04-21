!*==dladiv.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b DLADIV performs complex division in real arithmetic, avoiding unnecessary overflow.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DLADIV + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dladiv.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dladiv.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dladiv.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLADIV( A, B, C, D, P, Q )
!
!       .. Scalar Arguments ..
!       DOUBLE PRECISION   A, B, C, D, P, Q
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLADIV performs complex division in  real arithmetic
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
!>          A is DOUBLE PRECISION
!> \endverbatim
!>
!> \param[in] B
!> \verbatim
!>          B is DOUBLE PRECISION
!> \endverbatim
!>
!> \param[in] C
!> \verbatim
!>          C is DOUBLE PRECISION
!> \endverbatim
!>
!> \param[in] D
!> \verbatim
!>          D is DOUBLE PRECISION
!>          The scalars a, b, c, and d in the above expression.
!> \endverbatim
!>
!> \param[out] P
!> \verbatim
!>          P is DOUBLE PRECISION
!> \endverbatim
!>
!> \param[out] Q
!> \verbatim
!>          Q is DOUBLE PRECISION
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
!> \ingroup doubleOTHERauxiliary
!
!  =====================================================================
      SUBROUTINE DLADIV(A,B,C,D,P,Q)
      IMPLICIT NONE
!*--DLADIV95
!
!  -- LAPACK auxiliary routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     January 2013
!
!     .. Scalar Arguments ..
      DOUBLE PRECISION A , B , C , D , P , Q
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION BS
      PARAMETER (BS=2.0D0)
      DOUBLE PRECISION HALF
      PARAMETER (HALF=0.5D0)
      DOUBLE PRECISION TWO
      PARAMETER (TWO=2.0D0)
!
!     .. Local Scalars ..
      DOUBLE PRECISION aa , bb , cc , dd , ab , cd , s , ov , un , be , &
     &                 eps
!     ..
!     .. External Functions ..
      DOUBLE PRECISION DLAMCH
      EXTERNAL DLAMCH
!     ..
!     .. External Subroutines ..
      EXTERNAL DLADIV1
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , MAX
!     ..
!     .. Executable Statements ..
!
      aa = A
      bb = B
      cc = C
      dd = D
      ab = MAX(ABS(A),ABS(B))
      cd = MAX(ABS(C),ABS(D))
      s = 1.0D0
 
      ov = DLAMCH('Overflow threshold')
      un = DLAMCH('Safe minimum')
      eps = DLAMCH('Epsilon')
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
         CALL DLADIV1(aa,bb,cc,dd,P,Q)
      ELSE
         CALL DLADIV1(bb,aa,dd,cc,P,Q)
         Q = -Q
      ENDIF
      P = P*s
      Q = Q*s
!
!
!     End of DLADIV
!
      END SUBROUTINE DLADIV
!*==dladiv1.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
 
!> \ingroup doubleOTHERauxiliary
 
 
      SUBROUTINE DLADIV1(A,B,C,D,P,Q)
      IMPLICIT NONE
!*--DLADIV1185
!
!  -- LAPACK auxiliary routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     January 2013
!
!     .. Scalar Arguments ..
      DOUBLE PRECISION A , B , C , D , P , Q
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ONE
      PARAMETER (ONE=1.0D0)
!
!     .. Local Scalars ..
      DOUBLE PRECISION r , t
!     ..
!     .. External Functions ..
      DOUBLE PRECISION DLADIV2
      EXTERNAL DLADIV2
!     ..
!     .. Executable Statements ..
!
      r = D/C
      t = ONE/(C+D*r)
      P = DLADIV2(A,B,C,D,r,t)
      A = -A
      Q = DLADIV2(B,A,C,D,r,t)
!
!
!     End of DLADIV1
!
      END SUBROUTINE DLADIV1
!*==dladiv2.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
 
!> \ingroup doubleOTHERauxiliary
 
      DOUBLE PRECISION FUNCTION DLADIV2(A,B,C,D,R,T)
      IMPLICIT NONE
!*--DLADIV2227
!
!  -- LAPACK auxiliary routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     January 2013
!
!     .. Scalar Arguments ..
      DOUBLE PRECISION A , B , C , D , R , T
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
!
!     .. Local Scalars ..
      DOUBLE PRECISION br
!     ..
!     .. Executable Statements ..
!
      IF ( R/=ZERO ) THEN
         br = B*R
         IF ( br/=ZERO ) THEN
            DLADIV2 = (A+br)*T
         ELSE
            DLADIV2 = A*T + (B*T)*R
         ENDIF
      ELSE
         DLADIV2 = (A+D*(B/C))*T
      ENDIF
!
!
!     End of DLADIV12
!
      END FUNCTION DLADIV2
