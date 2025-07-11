!*==scnrm2.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b SCNRM2
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       REAL FUNCTION SCNRM2(N,X,INCX)
!
!       .. Scalar Arguments ..
!       INTEGER INCX,N
!       ..
!       .. Array Arguments ..
!       COMPLEX X(*)
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SCNRM2 returns the euclidean norm of a vector via the function
!> name, so that
!>
!>    SCNRM2 := sqrt( x**H*x )
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
!>          X is COMPLEX array, dimension (N)
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
!> \ingroup single_blas_level1
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  -- This version written on 25-October-1982.
!>     Modified on 14-October-1993 to inline the call to CLASSQ.
!>     Sven Hammarling, Nag Ltd.
!> \endverbatim
!>
!  =====================================================================
      REAL FUNCTION SCNRM2(N,X,Incx)
      IMPLICIT NONE
!*--SCNRM279
!
!  -- Reference BLAS level1 routine (version 3.8.0) --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2017
!
!     .. Scalar Arguments ..
      INTEGER Incx , N
!     ..
!     .. Array Arguments ..
      COMPLEX X(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ONE , ZERO
      PARAMETER (ONE=1.0E+0,ZERO=0.0E+0)
!     ..
!     .. Local Scalars ..
      REAL norm , scale , ssq , temp
      INTEGER ix
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , AIMAG , REAL , SQRT
!     ..
      IF ( N<1 .OR. Incx<1 ) THEN
         norm = ZERO
      ELSE
         scale = ZERO
         ssq = ONE
!        The following loop is equivalent to this call to the LAPACK
!        auxiliary routine:
!        CALL CLASSQ( N, X, INCX, SCALE, SSQ )
!
         DO ix = 1 , 1 + (N-1)*Incx , Incx
            IF ( REAL(X(ix))/=ZERO ) THEN
               temp = ABS(REAL(X(ix)))
               IF ( scale<temp ) THEN
                  ssq = ONE + ssq*(scale/temp)**2
                  scale = temp
               ELSE
                  ssq = ssq + (temp/scale)**2
               ENDIF
            ENDIF
            IF ( AIMAG(X(ix))/=ZERO ) THEN
               temp = ABS(AIMAG(X(ix)))
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
      SCNRM2 = norm
!
!     End of SCNRM2.
!
      END FUNCTION SCNRM2
