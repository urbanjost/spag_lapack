!*==clarfgp.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b CLARFGP generates an elementary reflector (Householder matrix) with non-negative beta.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CLARFGP + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clarfgp.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clarfgp.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clarfgp.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CLARFGP( N, ALPHA, X, INCX, TAU )
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
!> CLARFGP generates a complex elementary reflector H of order n, such
!> that
!>
!>       H**H * ( alpha ) = ( beta ),   H**H * H = I.
!>              (   x   )   (   0  )
!>
!> where alpha and beta are scalars, beta is real and non-negative, and
!> x is an (n-1)-element complex vector.  H is represented in the form
!>
!>       H = I - tau * ( 1 ) * ( 1 v**H ) ,
!>                     ( v )
!>
!> where tau is a complex scalar and v is a complex (n-1)-element
!> vector. Note that H is not hermitian.
!>
!> If the elements of x are all zero and alpha is real, then tau = 0
!> and H is taken to be the unit matrix.
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
      SUBROUTINE CLARFGP(N,Alpha,X,Incx,Tau)
      IMPLICIT NONE
!*--CLARFGP108
!
!  -- LAPACK auxiliary routine (version 3.8.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2017
!
!     .. Scalar Arguments ..
      INTEGER Incx , N
      COMPLEX Alpha , Tau
!     ..
!     .. Array Arguments ..
      COMPLEX X(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL TWO , ONE , ZERO
      PARAMETER (TWO=2.0E+0,ONE=1.0E+0,ZERO=0.0E+0)
!     ..
!     .. Local Scalars ..
      INTEGER j , knt
      REAL alphi , alphr , beta , bignum , smlnum , xnorm
      COMPLEX savealpha
!     ..
!     .. External Functions ..
      REAL SCNRM2 , SLAMCH , SLAPY3 , SLAPY2
      COMPLEX CLADIV
      EXTERNAL SCNRM2 , SLAMCH , SLAPY3 , SLAPY2 , CLADIV
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , AIMAG , CMPLX , REAL , SIGN
!     ..
!     .. External Subroutines ..
      EXTERNAL CSCAL , CSSCAL
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
      IF ( xnorm/=ZERO ) THEN
!
!        general case
!
         beta = SIGN(SLAPY3(alphr,alphi,xnorm),alphr)
         smlnum = SLAMCH('S')/SLAMCH('E')
         bignum = ONE/smlnum
!
         knt = 0
         IF ( ABS(beta)<smlnum ) THEN
            DO
!
!           XNORM, BETA may be inaccurate; scale X and recompute them
!
               knt = knt + 1
               CALL CSSCAL(N-1,bignum,X,Incx)
               beta = beta*bignum
               alphi = alphi*bignum
               alphr = alphr*bignum
               IF ( (ABS(beta)>=smlnum) .OR. (knt>=20) ) THEN
!
!           New BETA is at most 1, at least SMLNUM
!
                  xnorm = SCNRM2(N-1,X,Incx)
                  Alpha = CMPLX(alphr,alphi)
                  beta = SIGN(SLAPY3(alphr,alphi,xnorm),alphr)
                  EXIT
               ENDIF
            ENDDO
         ENDIF
         savealpha = Alpha
         Alpha = Alpha + beta
         IF ( beta<ZERO ) THEN
            beta = -beta
            Tau = -Alpha/beta
         ELSE
            alphr = alphi*(alphi/REAL(Alpha))
            alphr = alphr + xnorm*(xnorm/REAL(Alpha))
            Tau = CMPLX(alphr/beta,-alphi/beta)
            Alpha = CMPLX(-alphr,alphi)
         ENDIF
         Alpha = CLADIV(CMPLX(ONE),Alpha)
!
         IF ( ABS(Tau)<=smlnum ) THEN
!
!           In the case where the computed TAU ends up being a denormalized number,
!           it loses relative accuracy. This is a BIG problem. Solution: flush TAU
!           to ZERO (or TWO or whatever makes a nonnegative real number for BETA).
!
!           (Bug report provided by Pat Quillen from MathWorks on Jul 29, 2009.)
!           (Thanks Pat. Thanks MathWorks.)
!
            alphr = REAL(savealpha)
            alphi = AIMAG(savealpha)
            IF ( alphi/=ZERO ) THEN
               xnorm = SLAPY2(alphr,alphi)
               Tau = CMPLX(ONE-alphr/xnorm,-alphi/xnorm)
               DO j = 1 , N - 1
                  X(1+(j-1)*Incx) = ZERO
               ENDDO
               beta = xnorm
            ELSEIF ( alphr>=ZERO ) THEN
               Tau = ZERO
            ELSE
               Tau = TWO
               DO j = 1 , N - 1
                  X(1+(j-1)*Incx) = ZERO
               ENDDO
               beta = -savealpha
            ENDIF
!
         ELSE
!
!           This is the general case.
!
            CALL CSCAL(N-1,Alpha,X,Incx)
!
         ENDIF
!
!        If BETA is subnormal, it may lose relative accuracy
!
         DO j = 1 , knt
            beta = beta*smlnum
         ENDDO
         Alpha = beta
!
!        H  =  [1-alpha/abs(alpha) 0; 0 I], sign chosen so ALPHA >= 0.
!
      ELSEIF ( alphi/=ZERO ) THEN
!           Only "reflecting" the diagonal entry to be real and non-negative.
         xnorm = SLAPY2(alphr,alphi)
         Tau = CMPLX(ONE-alphr/xnorm,-alphi/xnorm)
         DO j = 1 , N - 1
            X(1+(j-1)*Incx) = ZERO
         ENDDO
         Alpha = xnorm
      ELSEIF ( alphr>=ZERO ) THEN
!              When TAU.eq.ZERO, the vector is special-cased to be
!              all zeros in the application routines.  We do not need
!              to clear it.
         Tau = ZERO
      ELSE
!              However, the application routines rely on explicit
!              zero checks when TAU.ne.ZERO, and we must clear X.
         Tau = TWO
         DO j = 1 , N - 1
            X(1+(j-1)*Incx) = ZERO
         ENDDO
         Alpha = -Alpha
      ENDIF
!
!
!     End of CLARFGP
!
      END SUBROUTINE CLARFGP
