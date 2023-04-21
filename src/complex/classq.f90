!*==classq.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b CLASSQ updates a sum of squares represented in scaled form.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CLASSQ + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/classq.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/classq.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/classq.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CLASSQ( N, X, INCX, SCALE, SUMSQ )
!
!       .. Scalar Arguments ..
!       INTEGER            INCX, N
!       REAL               SCALE, SUMSQ
!       ..
!       .. Array Arguments ..
!       COMPLEX            X( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CLASSQ returns the values scl and ssq such that
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
!>          X is COMPLEX array, dimension (1+(N-1)*INCX)
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
!>          SCALE is REAL
!>          On entry, the value  scale  in the equation above.
!>          On exit, SCALE is overwritten with the value  scl .
!> \endverbatim
!>
!> \param[in,out] SUMSQ
!> \verbatim
!>          SUMSQ is REAL
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
!> \ingroup complexOTHERauxiliary
!
!  =====================================================================
      SUBROUTINE CLASSQ(N,X,Incx,Scale,Sumsq)
      IMPLICIT NONE
!*--CLASSQ110
!
!  -- LAPACK auxiliary routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Incx , N
      REAL Scale , Sumsq
!     ..
!     .. Array Arguments ..
      COMPLEX X(*)
!     ..
!
! =====================================================================
!
!     .. Parameters ..
      REAL ZERO
      PARAMETER (ZERO=0.0E+0)
!     ..
!     .. Local Scalars ..
      INTEGER ix
      REAL temp1
!     ..
!     .. External Functions ..
      LOGICAL SISNAN
      EXTERNAL SISNAN
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , AIMAG , REAL
!     ..
!     .. Executable Statements ..
!
      IF ( N>0 ) THEN
         DO ix = 1 , 1 + (N-1)*Incx , Incx
            temp1 = ABS(REAL(X(ix)))
            IF ( temp1>ZERO .OR. SISNAN(temp1) ) THEN
               IF ( Scale<temp1 ) THEN
                  Sumsq = 1 + Sumsq*(Scale/temp1)**2
                  Scale = temp1
               ELSE
                  Sumsq = Sumsq + (temp1/Scale)**2
               ENDIF
            ENDIF
            temp1 = ABS(AIMAG(X(ix)))
            IF ( temp1>ZERO .OR. SISNAN(temp1) ) THEN
               IF ( Scale<temp1 .OR. SISNAN(temp1) ) THEN
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
!     End of CLASSQ
!
      END SUBROUTINE CLASSQ
