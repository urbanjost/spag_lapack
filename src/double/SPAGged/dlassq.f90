!*==dlassq.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b DLASSQ updates a sum of squares represented in scaled form.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DLASSQ + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlassq.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlassq.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlassq.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLASSQ( N, X, INCX, SCALE, SUMSQ )
!
!       .. Scalar Arguments ..
!       INTEGER            INCX, N
!       DOUBLE PRECISION   SCALE, SUMSQ
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   X( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLASSQ  returns the values  scl  and  smsq  such that
!>
!>    ( scl**2 )*smsq = x( 1 )**2 +...+ x( n )**2 + ( scale**2 )*sumsq,
!>
!> where  x( i ) = X( 1 + ( i - 1 )*INCX ). The value of  sumsq  is
!> assumed to be non-negative and  scl  returns the value
!>
!>    scl = max( scale, abs( x( i ) ) ).
!>
!> scale and sumsq must be supplied in SCALE and SUMSQ and
!> scl and smsq are overwritten on SCALE and SUMSQ respectively.
!>
!> The routine makes only one pass through the vector x.
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
!>          X is DOUBLE PRECISION array, dimension (1+(N-1)*INCX)
!>          The vector for which a scaled sum of squares is computed.
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
!>          On exit, SCALE is overwritten with  scl , the scaling factor
!>          for the sum of squares.
!> \endverbatim
!>
!> \param[in,out] SUMSQ
!> \verbatim
!>          SUMSQ is DOUBLE PRECISION
!>          On entry, the value  sumsq  in the equation above.
!>          On exit, SUMSQ is overwritten with  smsq , the basic sum of
!>          squares from which  scl  has been factored out.
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
!> \ingroup OTHERauxiliary
!
!  =====================================================================
      SUBROUTINE DLASSQ(N,X,Incx,Scale,Sumsq)
      USE F77KINDS                        
      USE S_DISNAN
      IMPLICIT NONE
!*--DLASSQ109
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: X
      INTEGER , INTENT(IN) :: Incx
      REAL(R8KIND) , INTENT(INOUT) :: Scale
      REAL(R8KIND) , INTENT(INOUT) :: Sumsq
!
! Local variable declarations rewritten by SPAG
!
      REAL(R8KIND) :: absxi
      INTEGER :: ix
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
            absxi = ABS(X(ix))
            IF ( absxi>ZERO .OR. DISNAN(absxi) ) THEN
               IF ( Scale<absxi ) THEN
                  Sumsq = 1 + Sumsq*(Scale/absxi)**2
                  Scale = absxi
               ELSE
                  Sumsq = Sumsq + (absxi/Scale)**2
               ENDIF
            ENDIF
         ENDDO
      ENDIF
!
!     End of DLASSQ
!
      END SUBROUTINE DLASSQ
