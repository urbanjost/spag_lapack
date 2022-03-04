!*==dlarfg.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b DLARFG generates an elementary reflector (Householder matrix).
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DLARFG + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlarfg.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlarfg.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlarfg.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLARFG( N, ALPHA, X, INCX, TAU )
!
!       .. Scalar Arguments ..
!       INTEGER            INCX, N
!       DOUBLE PRECISION   ALPHA, TAU
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
!> DLARFG generates a real elementary reflector H of order n, such
!> that
!>
!>       H * ( alpha ) = ( beta ),   H**T * H = I.
!>           (   x   )   (   0  )
!>
!> where alpha and beta are scalars, and x is an (n-1)-element real
!> vector. H is represented in the form
!>
!>       H = I - tau * ( 1 ) * ( 1 v**T ) ,
!>                     ( v )
!>
!> where tau is a real scalar and v is a real (n-1)-element
!> vector.
!>
!> If the elements of x are all zero, then tau = 0 and H is taken to be
!> the unit matrix.
!>
!> Otherwise  1 <= tau <= 2.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the elementary reflector.
!> \endverbatim
!>
!> \param[in,out] ALPHA
!> \verbatim
!>          ALPHA is DOUBLE PRECISION
!>          On entry, the value alpha.
!>          On exit, it is overwritten with the value beta.
!> \endverbatim
!>
!> \param[in,out] X
!> \verbatim
!>          X is DOUBLE PRECISION array, dimension
!>                         (1+(N-2)*abs(INCX))
!>          On entry, the vector x.
!>          On exit, it is overwritten with the vector v.
!> \endverbatim
!>
!> \param[in] INCX
!> \verbatim
!>          INCX is INTEGER
!>          The increment between elements of X. INCX > 0.
!> \endverbatim
!>
!> \param[out] TAU
!> \verbatim
!>          TAU is DOUBLE PRECISION
!>          The value tau.
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
!> \ingroup doubleOTHERauxiliary
!
!  =====================================================================
      SUBROUTINE DLARFG(N,Alpha,X,Incx,Tau)
      USE F77KINDS                        
      USE S_DLAMCH
      USE S_DLAPY2
      USE S_DNRM2
      USE S_DSCAL
      IMPLICIT NONE
!*--DLARFG115
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(INOUT) :: Alpha
      REAL(R8KIND) , DIMENSION(*) :: X
      INTEGER :: Incx
      REAL(R8KIND) , INTENT(OUT) :: Tau
!
! Local variable declarations rewritten by SPAG
!
      REAL(R8KIND) :: beta , rsafmn , safmin , xnorm
      INTEGER :: j , knt
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
!     .. External Functions ..
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. External Subroutines ..
!     ..
!     .. Executable Statements ..
!
      IF ( N<=1 ) THEN
         Tau = ZERO
         RETURN
      ENDIF
!
      xnorm = DNRM2(N-1,X,Incx)
!
      IF ( xnorm==ZERO ) THEN
!
!        H  =  I
!
         Tau = ZERO
      ELSE
!
!        general case
!
         beta = -SIGN(DLAPY2(Alpha,xnorm),Alpha)
         safmin = DLAMCH('S')/DLAMCH('E')
         knt = 0
         IF ( ABS(beta)<safmin ) THEN
!
!           XNORM, BETA may be inaccurate; scale X and recompute them
!
            rsafmn = ONE/safmin
            DO
               knt = knt + 1
               CALL DSCAL(N-1,rsafmn,X,Incx)
               beta = beta*rsafmn
               Alpha = Alpha*rsafmn
               IF ( (ABS(beta)>=safmin) .OR. (knt>=20) ) THEN
!
!           New BETA is at most 1, at least SAFMIN
!
                  xnorm = DNRM2(N-1,X,Incx)
                  beta = -SIGN(DLAPY2(Alpha,xnorm),Alpha)
                  EXIT
               ENDIF
            ENDDO
         ENDIF
         Tau = (beta-Alpha)/beta
         CALL DSCAL(N-1,ONE/(Alpha-beta),X,Incx)
!
!        If ALPHA is subnormal, it may lose relative accuracy
!
         DO j = 1 , knt
            beta = beta*safmin
         ENDDO
         Alpha = beta
      ENDIF
!
!
!     End of DLARFG
!
      END SUBROUTINE DLARFG
