!*==zlarfg.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b ZLARFG generates an elementary reflector (Householder matrix).
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZLARFG + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlarfg.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlarfg.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlarfg.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZLARFG( N, ALPHA, X, INCX, TAU )
!
!       .. Scalar Arguments ..
!       INTEGER            INCX, N
!       COMPLEX*16         ALPHA, TAU
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
!> ZLARFG generates a complex elementary reflector H of order n, such
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
!>          ALPHA is COMPLEX*16
!>          On entry, the value alpha.
!>          On exit, it is overwritten with the value beta.
!> \endverbatim
!>
!> \param[in,out] X
!> \verbatim
!>          X is COMPLEX*16 array, dimension
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
!>          TAU is COMPLEX*16
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
!> \ingroup complex16OTHERauxiliary
!
!  =====================================================================
      SUBROUTINE ZLARFG(N,Alpha,X,Incx,Tau)
      USE F77KINDS                        
      USE S_DLAMCH
      USE S_DLAPY3
      USE S_DZNRM2
      USE S_ZDSCAL
      USE S_ZLADIV
      USE S_ZSCAL
      IMPLICIT NONE
!*--ZLARFG117
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(INOUT) :: Alpha
      COMPLEX(CX16KIND) , DIMENSION(*) :: X
      INTEGER :: Incx
      COMPLEX(CX16KIND) , INTENT(OUT) :: Tau
!
! Local variable declarations rewritten by SPAG
!
      REAL(R8KIND) :: alphi , alphr , beta , rsafmn , safmin , xnorm
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
      xnorm = DZNRM2(N-1,X,Incx)
      alphr = DBLE(Alpha)
      alphi = DIMAG(Alpha)
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
         beta = -SIGN(DLAPY3(alphr,alphi,xnorm),alphr)
         safmin = DLAMCH('S')/DLAMCH('E')
         rsafmn = ONE/safmin
!
         knt = 0
         IF ( ABS(beta)<safmin ) THEN
            DO
!
!           XNORM, BETA may be inaccurate; scale X and recompute them
!
               knt = knt + 1
               CALL ZDSCAL(N-1,rsafmn,X,Incx)
               beta = beta*rsafmn
               alphi = alphi*rsafmn
               alphr = alphr*rsafmn
               IF ( (ABS(beta)>=safmin) .OR. (knt>=20) ) THEN
!
!           New BETA is at most 1, at least SAFMIN
!
                  xnorm = DZNRM2(N-1,X,Incx)
                  Alpha = DCMPLX(alphr,alphi)
                  beta = -SIGN(DLAPY3(alphr,alphi,xnorm),alphr)
                  EXIT
               ENDIF
            ENDDO
         ENDIF
         Tau = DCMPLX((beta-alphr)/beta,-alphi/beta)
         Alpha = ZLADIV(DCMPLX(ONE),Alpha-beta)
         CALL ZSCAL(N-1,Alpha,X,Incx)
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
!     End of ZLARFG
!
      END SUBROUTINE ZLARFG
