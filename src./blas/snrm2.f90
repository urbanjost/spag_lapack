!*==snrm2.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
!> \brief \b SNRM2
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       REAL FUNCTION SNRM2(N,X,INCX)
!
!       .. Scalar Arguments ..
!       INTEGER INCX,N
!       ..
!       .. Array Arguments ..
!       REAL X(*)
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SNRM2 returns the euclidean norm of a vector via the function
!> name, so that
!>
!>    SNRM2 := sqrt( x'*x ).
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>         number of elements in input vector(s)
!> \endverbatim
!>
!> \param[in] X
!> \verbatim
!>          X is REAL array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
!> \endverbatim
!>
!> \param[in] INCX
!> \verbatim
!>          INCX is INTEGER
!>         storage spacing between elements of SX
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
!> \ingroup single_blas_level1
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  -- This version written on 25-October-1982.
!>     Modified on 14-October-1993 to inline the call to SLASSQ.
!>     Sven Hammarling, Nag Ltd.
!> \endverbatim
!>
!  =====================================================================
      FUNCTION SNRM2(N,X,Incx)
      IMPLICIT NONE
!*--SNRM278
!
! PARAMETER definitions rewritten by SPAG
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
      REAL :: SNRM2
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(IN) , DIMENSION(*) :: X
      INTEGER , INTENT(IN) :: Incx
!
! Local variable declarations rewritten by SPAG
!
      REAL :: absxi , norm , scale , ssq
      INTEGER :: ix
!
! End of declarations rewritten by SPAG
!
!
! Dummy argument declarations rewritten by SPAG
!
!
! Local variable declarations rewritten by SPAG
!
!
! End of declarations rewritten by SPAG
!
!     ..
!     .. Array Arguments ..
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
!     ..
!     .. Local Scalars ..
!     ..
!     .. Intrinsic Functions ..
!     ..
      IF ( N<1 .OR. Incx<1 ) THEN
         norm = ZERO
      ELSEIF ( N==1 ) THEN
         norm = ABS(X(1))
      ELSE
         scale = ZERO
         ssq = ONE
!        The following loop is equivalent to this call to the LAPACK
!        auxiliary routine:
!        CALL SLASSQ( N, X, INCX, SCALE, SSQ )
!
         DO ix = 1 , 1 + (N-1)*Incx , Incx
            IF ( X(ix)/=ZERO ) THEN
               absxi = ABS(X(ix))
               IF ( scale<absxi ) THEN
                  ssq = ONE + ssq*(scale/absxi)**2
                  scale = absxi
               ELSE
                  ssq = ssq + (absxi/scale)**2
               ENDIF
            ENDIF
         ENDDO
         norm = scale*SQRT(ssq)
      ENDIF
!
      SNRM2 = norm
!
!     End of SNRM2.
!
      END FUNCTION SNRM2
