!*==zlassq.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b ZLASSQ updates a sum of squares represented in scaled form.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZLASSQ + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlassq.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlassq.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlassq.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZLASSQ( N, X, INCX, SCALE, SUMSQ )
!
!       .. Scalar Arguments ..
!       INTEGER            INCX, N
!       DOUBLE PRECISION   SCALE, SUMSQ
!       ..
!       .. Array Arguments ..
!       COMPLEX*16         X( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZLASSQ returns the values scl and ssq such that
!>
!>    ( scl**2 )*ssq = x( 1 )**2 +...+ x( n )**2 + ( scale**2 )*sumsq,
!>
!> where x( i ) = abs( X( 1 + ( i - 1 )*INCX ) ). The value of sumsq is
!> assumed to be at least unity and the value of ssq will then satisfy
!>
!>    1.0 <= ssq <= ( sumsq + 2*n ).
!>
!> scale is assumed to be non-negative and scl returns the value
!>
!>    scl = max( scale, abs( real( x( i ) ) ), abs( aimag( x( i ) ) ) ),
!>           i
!>
!> scale and sumsq must be supplied in SCALE and SUMSQ respectively.
!> SCALE and SUMSQ are overwritten by scl and ssq respectively.
!>
!> The routine makes only one pass through the vector X.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of elements to be used from the vector X.
!> \endverbatim
!>
!> \param[in] X
!> \verbatim
!>          X is COMPLEX*16 array, dimension (1+(N-1)*INCX)
!>          The vector x as described above.
!>             x( i )  = X( 1 + ( i - 1 )*INCX ), 1 <= i <= n.
!> \endverbatim
!>
!> \param[in] INCX
!> \verbatim
!>          INCX is INTEGER
!>          The increment between successive values of the vector X.
!>          INCX > 0.
!> \endverbatim
!>
!> \param[in,out] SCALE
!> \verbatim
!>          SCALE is DOUBLE PRECISION
!>          On entry, the value  scale  in the equation above.
!>          On exit, SCALE is overwritten with the value  scl .
!> \endverbatim
!>
!> \param[in,out] SUMSQ
!> \verbatim
!>          SUMSQ is DOUBLE PRECISION
!>          On entry, the value  sumsq  in the equation above.
!>          On exit, SUMSQ is overwritten with the value  ssq .
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
!> \ingroup complex16OTHERauxiliary
!
!  =====================================================================
      SUBROUTINE ZLASSQ(N,X,Incx,Scale,Sumsq)
      USE F77KINDS                        
      USE S_DISNAN
      IMPLICIT NONE
!*--ZLASSQ112
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(*) :: X
      INTEGER , INTENT(IN) :: Incx
      REAL(R8KIND) , INTENT(INOUT) :: Scale
      REAL(R8KIND) , INTENT(INOUT) :: Sumsq
!
! Local variable declarations rewritten by SPAG
!
      INTEGER :: ix
      REAL(R8KIND) :: temp1
!
! End of declarations rewritten by SPAG
!
!     ..
!     .. Array Arguments ..
!     ..
!
! =====================================================================
!
!     .. Parameters ..
!     ..
!     .. Local Scalars ..
!     ..
!     .. External Functions ..
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. Executable Statements ..
!
      IF ( N>0 ) THEN
         DO ix = 1 , 1 + (N-1)*Incx , Incx
            temp1 = ABS(DBLE(X(ix)))
            IF ( temp1>ZERO .OR. DISNAN(temp1) ) THEN
               IF ( Scale<temp1 ) THEN
                  Sumsq = 1 + Sumsq*(Scale/temp1)**2
                  Scale = temp1
               ELSE
                  Sumsq = Sumsq + (temp1/Scale)**2
               ENDIF
            ENDIF
            temp1 = ABS(DIMAG(X(ix)))
            IF ( temp1>ZERO .OR. DISNAN(temp1) ) THEN
               IF ( Scale<temp1 ) THEN
                  Sumsq = 1 + Sumsq*(Scale/temp1)**2
                  Scale = temp1
               ELSE
                  Sumsq = Sumsq + (temp1/Scale)**2
               ENDIF
            ENDIF
         ENDDO
      ENDIF
!
!
!     End of ZLASSQ
!
      END SUBROUTINE ZLASSQ
