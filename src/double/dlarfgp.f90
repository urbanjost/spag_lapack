!*==dlarfgp.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b DLARFGP generates an elementary reflector (Householder matrix) with non-negative beta.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DLARFGP + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlarfgp.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlarfgp.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlarfgp.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLARFGP( N, ALPHA, X, INCX, TAU )
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
!> DLARFGP generates a real elementary reflector H of order n, such
!> that
!>
!>       H * ( alpha ) = ( beta ),   H**T * H = I.
!>           (   x   )   (   0  )
!>
!> where alpha and beta are scalars, beta is non-negative, and x is
!> an (n-1)-element real vector.  H is represented in the form
!>
!>       H = I - tau * ( 1 ) * ( 1 v**T ) ,
!>                     ( v )
!>
!> where tau is a real scalar and v is a real (n-1)-element
!> vector.
!>
!> If the elements of x are all zero, then tau = 0 and H is taken to be
!> the unit matrix.
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
      SUBROUTINE DLARFGP(N,Alpha,X,Incx,Tau)
      IMPLICIT NONE
!*--DLARFGP108
!
!  -- LAPACK auxiliary routine (version 3.8.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2017
!
!     .. Scalar Arguments ..
      INTEGER Incx , N
      DOUBLE PRECISION Alpha , Tau
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION X(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION TWO , ONE , ZERO
      PARAMETER (TWO=2.0D+0,ONE=1.0D+0,ZERO=0.0D+0)
!     ..
!     .. Local Scalars ..
      INTEGER j , knt
      DOUBLE PRECISION beta , bignum , savealpha , smlnum , xnorm
!     ..
!     .. External Functions ..
      DOUBLE PRECISION DLAMCH , DLAPY2 , DNRM2
      EXTERNAL DLAMCH , DLAPY2 , DNRM2
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , SIGN
!     ..
!     .. External Subroutines ..
      EXTERNAL DSCAL
!     ..
!     .. Executable Statements ..
!
      IF ( N<=0 ) THEN
         Tau = ZERO
         RETURN
      ENDIF
!
      xnorm = DNRM2(N-1,X,Incx)
!
      IF ( xnorm/=ZERO ) THEN
!
!        general case
!
         beta = SIGN(DLAPY2(Alpha,xnorm),Alpha)
         smlnum = DLAMCH('S')/DLAMCH('E')
         knt = 0
         IF ( ABS(beta)<smlnum ) THEN
!
!           XNORM, BETA may be inaccurate; scale X and recompute them
!
            bignum = ONE/smlnum
            DO
               knt = knt + 1
               CALL DSCAL(N-1,bignum,X,Incx)
               beta = beta*bignum
               Alpha = Alpha*bignum
               IF ( (ABS(beta)>=smlnum) .OR. (knt>=20) ) THEN
!
!           New BETA is at most 1, at least SMLNUM
!
                  xnorm = DNRM2(N-1,X,Incx)
                  beta = SIGN(DLAPY2(Alpha,xnorm),Alpha)
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
            Alpha = xnorm*(xnorm/Alpha)
            Tau = Alpha/beta
            Alpha = -Alpha
         ENDIF
!
         IF ( ABS(Tau)>smlnum ) THEN
!
!           This is the general case.
!
            CALL DSCAL(N-1,ONE/Alpha,X,Incx)
!
!           In the case where the computed TAU ends up being a denormalized number,
!           it loses relative accuracy. This is a BIG problem. Solution: flush TAU
!           to ZERO. This explains the next IF statement.
!
!           (Bug report provided by Pat Quillen from MathWorks on Jul 29, 2009.)
!           (Thanks Pat. Thanks MathWorks.)
!
         ELSEIF ( savealpha>=ZERO ) THEN
            Tau = ZERO
         ELSE
            Tau = TWO
            DO j = 1 , N - 1
               X(1+(j-1)*Incx) = 0
            ENDDO
            beta = -savealpha
!
         ENDIF
!
!        If BETA is subnormal, it may lose relative accuracy
!
         DO j = 1 , knt
!
            beta = beta*smlnum
         ENDDO
         Alpha = beta
!
!        H  =  [+/-1, 0; I], sign chosen so ALPHA >= 0
!
      ELSEIF ( Alpha>=ZERO ) THEN
!           When TAU.eq.ZERO, the vector is special-cased to be
!           all zeros in the application routines.  We do not need
!           to clear it.
         Tau = ZERO
      ELSE
!           However, the application routines rely on explicit
!           zero checks when TAU.ne.ZERO, and we must clear X.
         Tau = TWO
         DO j = 1 , N - 1
            X(1+(j-1)*Incx) = 0
         ENDDO
         Alpha = -Alpha
      ENDIF
!
!
!     End of DLARFGP
!
      END SUBROUTINE DLARFGP
