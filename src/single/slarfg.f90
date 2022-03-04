!*==slarfg.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b SLARFG generates an elementary reflector (Householder matrix).
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SLARFG + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slarfg.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slarfg.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slarfg.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SLARFG( N, ALPHA, X, INCX, TAU )
!
!       .. Scalar Arguments ..
!       INTEGER            INCX, N
!       REAL               ALPHA, TAU
!       ..
!       .. Array Arguments ..
!       REAL               X( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SLARFG generates a real elementary reflector H of order n, such
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
!>          ALPHA is REAL
!>          On entry, the value alpha.
!>          On exit, it is overwritten with the value beta.
!> \endverbatim
!>
!> \param[in,out] X
!> \verbatim
!>          X is REAL array, dimension
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
!>          TAU is REAL
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
!> \ingroup realOTHERauxiliary
!
!  =====================================================================
      SUBROUTINE SLARFG(N,Alpha,X,Incx,Tau)
      IMPLICIT NONE
!*--SLARFG110
!
!  -- LAPACK auxiliary routine (version 3.8.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2017
!
!     .. Scalar Arguments ..
      INTEGER Incx , N
      REAL Alpha , Tau
!     ..
!     .. Array Arguments ..
      REAL X(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ONE , ZERO
      PARAMETER (ONE=1.0E+0,ZERO=0.0E+0)
!     ..
!     .. Local Scalars ..
      INTEGER j , knt
      REAL beta , rsafmn , safmin , xnorm
!     ..
!     .. External Functions ..
      REAL SLAMCH , SLAPY2 , SNRM2
      EXTERNAL SLAMCH , SLAPY2 , SNRM2
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , SIGN
!     ..
!     .. External Subroutines ..
      EXTERNAL SSCAL
!     ..
!     .. Executable Statements ..
!
      IF ( N<=1 ) THEN
         Tau = ZERO
         RETURN
      ENDIF
!
      xnorm = SNRM2(N-1,X,Incx)
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
         beta = -SIGN(SLAPY2(Alpha,xnorm),Alpha)
         safmin = SLAMCH('S')/SLAMCH('E')
         knt = 0
         IF ( ABS(beta)<safmin ) THEN
!
!           XNORM, BETA may be inaccurate; scale X and recompute them
!
            rsafmn = ONE/safmin
            DO
               knt = knt + 1
               CALL SSCAL(N-1,rsafmn,X,Incx)
               beta = beta*rsafmn
               Alpha = Alpha*rsafmn
               IF ( (ABS(beta)>=safmin) .OR. (knt>=20) ) THEN
!
!           New BETA is at most 1, at least SAFMIN
!
                  xnorm = SNRM2(N-1,X,Incx)
                  beta = -SIGN(SLAPY2(Alpha,xnorm),Alpha)
                  EXIT
               ENDIF
            ENDDO
         ENDIF
         Tau = (beta-Alpha)/beta
         CALL SSCAL(N-1,ONE/(Alpha-beta),X,Incx)
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
!     End of SLARFG
!
      END SUBROUTINE SLARFG
