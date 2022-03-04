!*==slasd5.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b SLASD5 computes the square root of the i-th eigenvalue of a positive symmetric rank-one modification of a 2-by-2 diagonal matrix. Used by sbdsdc.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SLASD5 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slasd5.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slasd5.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slasd5.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SLASD5( I, D, Z, DELTA, RHO, DSIGMA, WORK )
!
!       .. Scalar Arguments ..
!       INTEGER            I
!       REAL               DSIGMA, RHO
!       ..
!       .. Array Arguments ..
!       REAL               D( 2 ), DELTA( 2 ), WORK( 2 ), Z( 2 )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> This subroutine computes the square root of the I-th eigenvalue
!> of a positive symmetric rank-one modification of a 2-by-2 diagonal
!> matrix
!>
!>            diag( D ) * diag( D ) +  RHO * Z * transpose(Z) .
!>
!> The diagonal entries in the array D are assumed to satisfy
!>
!>            0 <= D(i) < D(j)  for  i < j .
!>
!> We also assume RHO > 0 and that the Euclidean norm of the vector
!> Z is one.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] I
!> \verbatim
!>          I is INTEGER
!>         The index of the eigenvalue to be computed.  I = 1 or I = 2.
!> \endverbatim
!>
!> \param[in] D
!> \verbatim
!>          D is REAL array, dimension (2)
!>         The original eigenvalues.  We assume 0 <= D(1) < D(2).
!> \endverbatim
!>
!> \param[in] Z
!> \verbatim
!>          Z is REAL array, dimension (2)
!>         The components of the updating vector.
!> \endverbatim
!>
!> \param[out] DELTA
!> \verbatim
!>          DELTA is REAL array, dimension (2)
!>         Contains (D(j) - sigma_I) in its  j-th component.
!>         The vector DELTA contains the information necessary
!>         to construct the eigenvectors.
!> \endverbatim
!>
!> \param[in] RHO
!> \verbatim
!>          RHO is REAL
!>         The scalar in the symmetric updating formula.
!> \endverbatim
!>
!> \param[out] DSIGMA
!> \verbatim
!>          DSIGMA is REAL
!>         The computed sigma_I, the I-th updated eigenvalue.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is REAL array, dimension (2)
!>         WORK contains (D(j) + sigma_I) in its  j-th component.
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
!> \par Contributors:
!  ==================
!>
!>     Ren-Cang Li, Computer Science Division, University of California
!>     at Berkeley, USA
!>
!  =====================================================================
      SUBROUTINE SLASD5(I,D,Z,Delta,Rho,Dsigma,Work)
      IMPLICIT NONE
!*--SLASD5120
!
! PARAMETER definitions rewritten by SPAG
!
      REAL , PARAMETER  ::  ZERO = 0.0E+0 , ONE = 1.0E+0 ,              &
     &                      TWO = 2.0E+0 , THREE = 3.0E+0 ,             &
     &                      FOUR = 4.0E+0
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER , INTENT(IN) :: I
      REAL , INTENT(IN) , DIMENSION(2) :: D
      REAL , INTENT(IN) , DIMENSION(2) :: Z
      REAL , INTENT(OUT) , DIMENSION(2) :: Delta
      REAL , INTENT(IN) :: Rho
      REAL , INTENT(OUT) :: Dsigma
      REAL , INTENT(OUT) , DIMENSION(2) :: Work
!
! Local variable declarations rewritten by SPAG
!
      REAL :: b , c , del , delsq , tau , w
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
!     .. Intrinsic Functions ..
!     ..
!     .. Executable Statements ..
!
      del = D(2) - D(1)
      delsq = del*(D(2)+D(1))
      IF ( I==1 ) THEN
         w = ONE + FOUR*Rho*(Z(2)*Z(2)/(D(1)+THREE*D(2))-Z(1)*Z(1)      &
     &       /(THREE*D(1)+D(2)))/del
         IF ( w>ZERO ) THEN
            b = delsq + Rho*(Z(1)*Z(1)+Z(2)*Z(2))
            c = Rho*Z(1)*Z(1)*delsq
!
!           B > ZERO, always
!
!           The following TAU is DSIGMA * DSIGMA - D( 1 ) * D( 1 )
!
            tau = TWO*c/(b+SQRT(ABS(b*b-FOUR*c)))
!
!           The following TAU is DSIGMA - D( 1 )
!
            tau = tau/(D(1)+SQRT(D(1)*D(1)+tau))
            Dsigma = D(1) + tau
            Delta(1) = -tau
            Delta(2) = del - tau
            Work(1) = TWO*D(1) + tau
            Work(2) = (D(1)+tau) + D(2)
!           DELTA( 1 ) = -Z( 1 ) / TAU
!           DELTA( 2 ) = Z( 2 ) / ( DEL-TAU )
         ELSE
            b = -delsq + Rho*(Z(1)*Z(1)+Z(2)*Z(2))
            c = Rho*Z(2)*Z(2)*delsq
!
!           The following TAU is DSIGMA * DSIGMA - D( 2 ) * D( 2 )
!
            IF ( b>ZERO ) THEN
               tau = -TWO*c/(b+SQRT(b*b+FOUR*c))
            ELSE
               tau = (b-SQRT(b*b+FOUR*c))/TWO
            ENDIF
!
!           The following TAU is DSIGMA - D( 2 )
!
            tau = tau/(D(2)+SQRT(ABS(D(2)*D(2)+tau)))
            Dsigma = D(2) + tau
            Delta(1) = -(del+tau)
            Delta(2) = -tau
            Work(1) = D(1) + tau + D(2)
            Work(2) = TWO*D(2) + tau
!           DELTA( 1 ) = -Z( 1 ) / ( DEL+TAU )
!           DELTA( 2 ) = -Z( 2 ) / TAU
         ENDIF
!        TEMP = SQRT( DELTA( 1 )*DELTA( 1 )+DELTA( 2 )*DELTA( 2 ) )
!        DELTA( 1 ) = DELTA( 1 ) / TEMP
!        DELTA( 2 ) = DELTA( 2 ) / TEMP
      ELSE
!
!        Now I=2
!
         b = -delsq + Rho*(Z(1)*Z(1)+Z(2)*Z(2))
         c = Rho*Z(2)*Z(2)*delsq
!
!        The following TAU is DSIGMA * DSIGMA - D( 2 ) * D( 2 )
!
         IF ( b>ZERO ) THEN
            tau = (b+SQRT(b*b+FOUR*c))/TWO
         ELSE
            tau = TWO*c/(-b+SQRT(b*b+FOUR*c))
         ENDIF
!
!        The following TAU is DSIGMA - D( 2 )
!
         tau = tau/(D(2)+SQRT(D(2)*D(2)+tau))
         Dsigma = D(2) + tau
         Delta(1) = -(del+tau)
         Delta(2) = -tau
         Work(1) = D(1) + tau + D(2)
         Work(2) = TWO*D(2) + tau
!        DELTA( 1 ) = -Z( 1 ) / ( DEL+TAU )
!        DELTA( 2 ) = -Z( 2 ) / TAU
!        TEMP = SQRT( DELTA( 1 )*DELTA( 1 )+DELTA( 2 )*DELTA( 2 ) )
!        DELTA( 1 ) = DELTA( 1 ) / TEMP
!        DELTA( 2 ) = DELTA( 2 ) / TEMP
      ENDIF
!
!     End of SLASD5
!
      END SUBROUTINE SLASD5
