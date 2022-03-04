!*==dsxt1.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b DSXT1
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       DOUBLE PRECISION FUNCTION DSXT1( IJOB, D1, N1, D2, N2, ABSTOL,
!                        ULP, UNFL )
!
!       .. Scalar Arguments ..
!       INTEGER            IJOB, N1, N2
!       DOUBLE PRECISION   ABSTOL, ULP, UNFL
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   D1( * ), D2( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DSXT1  computes the difference between a set of eigenvalues.
!> The result is returned as the function value.
!>
!> IJOB = 1:   Computes   max { min | D1(i)-D2(j) | }
!>                         i     j
!>
!> IJOB = 2:   Computes   max { min | D1(i)-D2(j) | /
!>                         i     j
!>                              ( ABSTOL + |D1(i)|*ULP ) }
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] IJOB
!> \verbatim
!>          IJOB is INTEGER
!>          Specifies the type of tests to be performed.  (See above.)
!> \endverbatim
!>
!> \param[in] D1
!> \verbatim
!>          D1 is DOUBLE PRECISION array, dimension (N1)
!>          The first array.  D1 should be in increasing order, i.e.,
!>          D1(j) <= D1(j+1).
!> \endverbatim
!>
!> \param[in] N1
!> \verbatim
!>          N1 is INTEGER
!>          The length of D1.
!> \endverbatim
!>
!> \param[in] D2
!> \verbatim
!>          D2 is DOUBLE PRECISION array, dimension (N2)
!>          The second array.  D2 should be in increasing order, i.e.,
!>          D2(j) <= D2(j+1).
!> \endverbatim
!>
!> \param[in] N2
!> \verbatim
!>          N2 is INTEGER
!>          The length of D2.
!> \endverbatim
!>
!> \param[in] ABSTOL
!> \verbatim
!>          ABSTOL is DOUBLE PRECISION
!>          The absolute tolerance, used as a measure of the error.
!> \endverbatim
!>
!> \param[in] ULP
!> \verbatim
!>          ULP is DOUBLE PRECISION
!>          Machine precision.
!> \endverbatim
!>
!> \param[in] UNFL
!> \verbatim
!>          UNFL is DOUBLE PRECISION
!>          The smallest positive number whose reciprocal does not
!>          overflow.
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
      DOUBLE PRECISION FUNCTION DSXT1(Ijob,D1,N1,D2,N2,Abstol,Ulp,Unfl)
      IMPLICIT NONE
!*--DSXT1109
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Ijob , N1 , N2
      DOUBLE PRECISION Abstol , Ulp , Unfl
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION D1(*) , D2(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
!     ..
!     .. Local Scalars ..
      INTEGER i , j
      DOUBLE PRECISION temp1 , temp2
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , MAX , MIN
!     ..
!     .. Executable Statements ..
!
      temp1 = ZERO
!
      j = 1
      DO i = 1 , N1
         DO
            IF ( D2(j)<D1(i) .AND. j<N2 ) THEN
               j = j + 1
               CYCLE
            ENDIF
            IF ( j==1 ) THEN
               temp2 = ABS(D2(j)-D1(i))
               IF ( Ijob==2 ) temp2 = temp2/MAX(Unfl,Abstol+Ulp*ABS(D1(i&
     &                                )))
            ELSE
               temp2 = MIN(ABS(D2(j)-D1(i)),ABS(D1(i)-D2(j-1)))
               IF ( Ijob==2 ) temp2 = temp2/MAX(Unfl,Abstol+Ulp*ABS(D1(i&
     &                                )))
            ENDIF
            temp1 = MAX(temp1,temp2)
            EXIT
         ENDDO
      ENDDO
!
      DSXT1 = temp1
!
!     End of DSXT1
!
      END FUNCTION DSXT1
