!*==dznrm2.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b DZNRM2
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       DOUBLE PRECISION FUNCTION DZNRM2(N,X,INCX)
!
!       .. Scalar Arguments ..
!       INTEGER INCX,N
!       ..
!       .. Array Arguments ..
!       COMPLEX*16 X(*)
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DZNRM2 returns the euclidean norm of a vector via the function
!> name, so that
!>
!>    DZNRM2 := sqrt( x**H*x )
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
!>          X is COMPLEX*16 array, dimension (N)
!>         complex vector with N elements
!> \endverbatim
!>
!> \param[in] INCX
!> \verbatim
!>          INCX is INTEGER
!>         storage spacing between elements of X
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
!> \ingroup double_blas_level1
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  -- This version written on 25-October-1982.
!>     Modified on 14-October-1993 to inline the call to ZLASSQ.
!>     Sven Hammarling, Nag Ltd.
!> \endverbatim
!>
!  =====================================================================
      FUNCTION DZNRM2(N,X,Incx)
      USE F77KINDS                        
      IMPLICIT NONE
!*--DZNRM280
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
      REAL(R8KIND) :: DZNRM2
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(*) :: X
      INTEGER , INTENT(IN) :: Incx
!
! Local variable declarations rewritten by SPAG
!
      INTEGER :: ix
      REAL(R8KIND) :: norm , scale , ssq , temp
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
      ELSE
         scale = ZERO
         ssq = ONE
!        The following loop is equivalent to this call to the LAPACK
!        auxiliary routine:
!        CALL ZLASSQ( N, X, INCX, SCALE, SSQ )
!
         DO ix = 1 , 1 + (N-1)*Incx , Incx
            IF ( DBLE(X(ix))/=ZERO ) THEN
               temp = ABS(DBLE(X(ix)))
               IF ( scale<temp ) THEN
                  ssq = ONE + ssq*(scale/temp)**2
                  scale = temp
               ELSE
                  ssq = ssq + (temp/scale)**2
               ENDIF
            ENDIF
            IF ( DIMAG(X(ix))/=ZERO ) THEN
               temp = ABS(DIMAG(X(ix)))
               IF ( scale<temp ) THEN
                  ssq = ONE + ssq*(scale/temp)**2
                  scale = temp
               ELSE
                  ssq = ssq + (temp/scale)**2
               ENDIF
            ENDIF
         ENDDO
         norm = scale*SQRT(ssq)
      ENDIF
!
      DZNRM2 = norm
!
!     End of DZNRM2.
!
      END FUNCTION DZNRM2
