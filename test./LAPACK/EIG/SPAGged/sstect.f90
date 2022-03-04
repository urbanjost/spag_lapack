!*==sstect.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b SSTECT
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE SSTECT( N, A, B, SHIFT, NUM )
!
!       .. Scalar Arguments ..
!       INTEGER            N, NUM
!       REAL               SHIFT
!       ..
!       .. Array Arguments ..
!       REAL               A( * ), B( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    SSTECT counts the number NUM of eigenvalues of a tridiagonal
!>    matrix T which are less than or equal to SHIFT. T has
!>    diagonal entries A(1), ... , A(N), and offdiagonal entries
!>    B(1), ..., B(N-1).
!>    See W. Kahan "Accurate Eigenvalues of a Symmetric Tridiagonal
!>    Matrix", Report CS41, Computer Science Dept., Stanford
!>    University, July 21, 1966
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The dimension of the tridiagonal matrix T.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is REAL array, dimension (N)
!>          The diagonal entries of the tridiagonal matrix T.
!> \endverbatim
!>
!> \param[in] B
!> \verbatim
!>          B is REAL array, dimension (N-1)
!>          The offdiagonal entries of the tridiagonal matrix T.
!> \endverbatim
!>
!> \param[in] SHIFT
!> \verbatim
!>          SHIFT is REAL
!>          The shift, used as described under Purpose.
!> \endverbatim
!>
!> \param[out] NUM
!> \verbatim
!>          NUM is INTEGER
!>          The number of eigenvalues of T less than or equal
!>          to SHIFT.
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
!> \ingroup single_eig
!
!  =====================================================================
      SUBROUTINE SSTECT(N,A,B,Shift,Num)
      IMPLICIT NONE
!*--SSTECT86
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER N , Num
      REAL Shift
!     ..
!     .. Array Arguments ..
      REAL A(*) , B(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ZERO , ONE , THREE
      PARAMETER (ZERO=0.0E0,ONE=1.0E0,THREE=3.0E0)
!     ..
!     .. Local Scalars ..
      INTEGER i
      REAL m1 , m2 , mx , ovfl , sov , sshift , ssun , sun , tmp , tom ,&
     &     u , unfl
!     ..
!     .. External Functions ..
      REAL SLAMCH
      EXTERNAL SLAMCH
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , MAX , SQRT
!     ..
!     .. Executable Statements ..
!
!     Get machine constants
!
      unfl = SLAMCH('Safe minimum')
      ovfl = SLAMCH('Overflow')
!
!     Find largest entry
!
      mx = ABS(A(1))
      DO i = 1 , N - 1
         mx = MAX(mx,ABS(A(i+1)),ABS(B(i)))
      ENDDO
!
!     Handle easy cases, including zero matrix
!
      IF ( Shift>=THREE*mx ) THEN
         Num = N
         RETURN
      ENDIF
      IF ( Shift<-THREE*mx ) THEN
         Num = 0
         RETURN
      ENDIF
!
!     Compute scale factors as in Kahan's report
!     At this point, MX .NE. 0 so we can divide by it
!
      sun = SQRT(unfl)
      ssun = SQRT(sun)
      sov = SQRT(ovfl)
      tom = ssun*sov
      IF ( mx<=ONE ) THEN
         m1 = ONE/mx
         m2 = tom
      ELSE
         m1 = ONE
         m2 = tom/mx
      ENDIF
!
!     Begin counting
!
      Num = 0
      sshift = (Shift*m1)*m2
      u = (A(1)*m1)*m2 - sshift
      IF ( u<=sun ) THEN
         IF ( u<=ZERO ) THEN
            Num = Num + 1
            IF ( u>-sun ) u = -sun
         ELSE
            u = sun
         ENDIF
      ENDIF
      DO i = 2 , N
         tmp = (B(i-1)*m1)*m2
         u = ((A(i)*m1)*m2-tmp*(tmp/u)) - sshift
         IF ( u<=sun ) THEN
            IF ( u<=ZERO ) THEN
               Num = Num + 1
               IF ( u>-sun ) u = -sun
            ELSE
               u = sun
            ENDIF
         ENDIF
      ENDDO
!
!     End of SSTECT
!
      END SUBROUTINE SSTECT