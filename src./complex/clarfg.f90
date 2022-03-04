!*==clarfg.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b CLARFG generates an elementary reflector (Householder matrix).
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CLARFG + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clarfg.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clarfg.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clarfg.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CLARFG( N, ALPHA, X, INCX, TAU )
!
!       .. Scalar Arguments ..
!       INTEGER            INCX, N
!       COMPLEX            ALPHA, TAU
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
!> CLARFG generates a complex elementary reflector H of order n, such
!> that
!>
!>       H**H * ( alpha ) = ( beta ),   H**H * H = I.
!>              (   x   )   (   0  )
!>
!> where alpha and beta are scalars, with beta real, and x is an
!> (n-1)-element complex vector. H is represented in the form
!>
!>       H = I - tau * ( 1 ) * ( 1 v**H ) ,
!>                     ( v )
!>
!> where tau is a complex scalar and v is a complex (n-1)-element
!> vector. Note that H is not hermitian.
!>
!> If the elements of x are all zero and alpha is real, then tau = 0
!> and H is taken to be the unit matrix.
!>
!> Otherwise  1 <= real(tau) <= 2  and  abs(tau-1) <= 1 .
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
!>          ALPHA is COMPLEX
!>          On entry, the value alpha.
!>          On exit, it is overwritten with the value beta.
!> \endverbatim
!>
!> \param[in,out] X
!> \verbatim
!>          X is COMPLEX array, dimension
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
!>          TAU is COMPLEX
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
!> \ingroup complexOTHERauxiliary
!
!  =====================================================================
      SUBROUTINE CLARFG(N,Alpha,X,Incx,Tau)
      USE S_CLADIV
      USE S_CSCAL
      USE S_CSSCAL
      USE S_SCNRM2
      USE S_SLAMCH
      USE S_SLAPY3
      IMPLICIT NONE
!*--CLARFG116
!
! PARAMETER definitions rewritten by SPAG
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(INOUT) :: Alpha
      COMPLEX , DIMENSION(*) :: X
      INTEGER :: Incx
      COMPLEX , INTENT(OUT) :: Tau
!
! Local variable declarations rewritten by SPAG
!
      REAL :: alphi , alphr , beta , rsafmn , safmin , xnorm
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
      IF ( N<=0 ) THEN
         Tau = ZERO
         RETURN
      ENDIF
!
      xnorm = SCNRM2(N-1,X,Incx)
      alphr = REAL(Alpha)
      alphi = AIMAG(Alpha)
!
      IF ( xnorm==ZERO .AND. alphi==ZERO ) THEN
!
!        H  =  I
!
         Tau = ZERO
      ELSE
!
!        general case
!
         beta = -SIGN(SLAPY3(alphr,alphi,xnorm),alphr)
         safmin = SLAMCH('S')/SLAMCH('E')
         rsafmn = ONE/safmin
!
         knt = 0
         IF ( ABS(beta)<safmin ) THEN
            DO
!
!           XNORM, BETA may be inaccurate; scale X and recompute them
!
               knt = knt + 1
               CALL CSSCAL(N-1,rsafmn,X,Incx)
               beta = beta*rsafmn
               alphi = alphi*rsafmn
               alphr = alphr*rsafmn
               IF ( (ABS(beta)>=safmin) .OR. (knt>=20) ) THEN
!
!           New BETA is at most 1, at least SAFMIN
!
                  xnorm = SCNRM2(N-1,X,Incx)
                  Alpha = CMPLX(alphr,alphi)
                  beta = -SIGN(SLAPY3(alphr,alphi,xnorm),alphr)
                  EXIT
               ENDIF
            ENDDO
         ENDIF
         Tau = CMPLX((beta-alphr)/beta,-alphi/beta)
         Alpha = CLADIV(CMPLX(ONE),Alpha-beta)
         CALL CSCAL(N-1,Alpha,X,Incx)
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
!     End of CLARFG
!
      END SUBROUTINE CLARFG
