!*==dsvdct.f90  processed by SPAG 7.51RB at 17:37 on  4 Mar 2022
!> \brief \b dsvdct
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE dsvdct( N, S, E, SHIFT, NUM )
!
!       .. Scalar Arguments ..
!       INTEGER            N, NUM
!       DOUBLE PRECISION   SHIFT
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   E( * ), S( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> dsvdct counts the number NUM of eigenvalues of a 2*N by 2*N
!> tridiagonal matrix T which are less than or equal to SHIFT.  T is
!> formed by putting zeros on the diagonal and making the off-diagonals
!> equal to S(1), E(1), S(2), E(2), ... , E(N-1), S(N).  If SHIFT is
!> positive, NUM is equal to N plus the number of singular values of a
!> bidiagonal matrix B less than or equal to SHIFT.  Here B has diagonal
!> entries S(1), ..., S(N) and superdiagonal entries E(1), ... E(N-1).
!> If SHIFT is negative, NUM is equal to the number of singular values
!> of B greater than or equal to -SHIFT.
!>
!> See W. Kahan "Accurate Eigenvalues of a Symmetric Tridiagonal
!> Matrix", Report CS41, Computer Science Dept., Stanford University,
!> July 21, 1966
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The dimension of the bidiagonal matrix B.
!> \endverbatim
!>
!> \param[in] S
!> \verbatim
!>          S is DOUBLE PRECISION array, dimension (N)
!>          The diagonal entries of the bidiagonal matrix B.
!> \endverbatim
!>
!> \param[in] E
!> \verbatim
!>          E is DOUBLE PRECISION array of dimension (N-1)
!>          The superdiagonal entries of the bidiagonal matrix B.
!> \endverbatim
!>
!> \param[in] SHIFT
!> \verbatim
!>          SHIFT is DOUBLE PRECISION
!>          The shift, used as described under Purpose.
!> \endverbatim
!>
!> \param[out] NUM
!> \verbatim
!>          NUM is INTEGER
!>          The number of eigenvalues of T less than or equal to SHIFT.
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
!> \ingroup double_eig
!
!  =====================================================================
      SUBROUTINE dsvdct(N,S,E,Shift,Num)
      IMPLICIT NONE
!*--DSVDCT91
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER N , Num
      DOUBLE PRECISION Shift
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION E(*) , S(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ONE
      PARAMETER (ONE=1.0D0)
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
!     ..
!     .. Local Scalars ..
      INTEGER i
      DOUBLE PRECISION m1 , m2 , mx , ovfl , sov , sshift , ssun , sun ,&
     &                 tmp , tom , u , unfl
!     ..
!     .. External Functions ..
      DOUBLE PRECISION DLAMCH
      EXTERNAL DLAMCH
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , MAX , SQRT
!     ..
!     .. Executable Statements ..
!
!     Get machine constants
!
      unfl = 2*DLAMCH('Safe minimum')
      ovfl = ONE/unfl
!
!     Find largest entry
!
      mx = ABS(S(1))
      DO i = 1 , N - 1
         mx = MAX(mx,ABS(S(i+1)),ABS(E(i)))
      ENDDO
!
      IF ( mx==ZERO ) THEN
         IF ( Shift<ZERO ) THEN
            Num = 0
         ELSE
            Num = 2*N
         ENDIF
         RETURN
      ENDIF
!
!     Compute scale factors as in Kahan's report
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
      u = ONE
      Num = 0
      sshift = (Shift*m1)*m2
      u = -sshift
      IF ( u<=sun ) THEN
         IF ( u<=ZERO ) THEN
            Num = Num + 1
            IF ( u>-sun ) u = -sun
         ELSE
            u = sun
         ENDIF
      ENDIF
      tmp = (S(1)*m1)*m2
      u = -tmp*(tmp/u) - sshift
      IF ( u<=sun ) THEN
         IF ( u<=ZERO ) THEN
            Num = Num + 1
            IF ( u>-sun ) u = -sun
         ELSE
            u = sun
         ENDIF
      ENDIF
      DO i = 1 , N - 1
         tmp = (E(i)*m1)*m2
         u = -tmp*(tmp/u) - sshift
         IF ( u<=sun ) THEN
            IF ( u<=ZERO ) THEN
               Num = Num + 1
               IF ( u>-sun ) u = -sun
            ELSE
               u = sun
            ENDIF
         ENDIF
         tmp = (S(i+1)*m1)*m2
         u = -tmp*(tmp/u) - sshift
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
!     End of dsvdct
!
      END SUBROUTINE dsvdct
