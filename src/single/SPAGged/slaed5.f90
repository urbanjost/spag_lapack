!*==slaed5.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b SLAED5 used by sstedc. Solves the 2-by-2 secular equation.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SLAED5 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slaed5.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slaed5.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slaed5.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SLAED5( I, D, Z, DELTA, RHO, DLAM )
!
!       .. Scalar Arguments ..
!       INTEGER            I
!       REAL               DLAM, RHO
!       ..
!       .. Array Arguments ..
!       REAL               D( 2 ), DELTA( 2 ), Z( 2 )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> This subroutine computes the I-th eigenvalue of a symmetric rank-one
!> modification of a 2-by-2 diagonal matrix
!>
!>            diag( D )  +  RHO * Z * transpose(Z) .
!>
!> The diagonal elements in the array D are assumed to satisfy
!>
!>            D(i) < D(j)  for  i < j .
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
!>         The original eigenvalues.  We assume D(1) < D(2).
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
!> \param[out] DLAM
!> \verbatim
!>          DLAM is REAL
!>         The computed lambda_I, the I-th updated eigenvalue.
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
!> \ingroup auxOTHERcomputational
!
!> \par Contributors:
!  ==================
!>
!>     Ren-Cang Li, Computer Science Division, University of California
!>     at Berkeley, USA
!>
!  =====================================================================
      SUBROUTINE SLAED5(I,D,Z,Delta,Rho,Dlam)
      IMPLICIT NONE
!*--SLAED5112
!
! PARAMETER definitions rewritten by SPAG
!
      REAL , PARAMETER  ::  ZERO = 0.0E0 , ONE = 1.0E0 , TWO = 2.0E0 ,  &
     &                      FOUR = 4.0E0
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER , INTENT(IN) :: I
      REAL , INTENT(IN) , DIMENSION(2) :: D
      REAL , INTENT(IN) , DIMENSION(2) :: Z
      REAL , INTENT(INOUT) , DIMENSION(2) :: Delta
      REAL , INTENT(IN) :: Rho
      REAL , INTENT(OUT) :: Dlam
!
! Local variable declarations rewritten by SPAG
!
      REAL :: b , c , del , tau , temp , w
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
      IF ( I==1 ) THEN
         w = ONE + TWO*Rho*(Z(2)*Z(2)-Z(1)*Z(1))/del
         IF ( w>ZERO ) THEN
            b = del + Rho*(Z(1)*Z(1)+Z(2)*Z(2))
            c = Rho*Z(1)*Z(1)*del
!
!           B > ZERO, always
!
            tau = TWO*c/(b+SQRT(ABS(b*b-FOUR*c)))
            Dlam = D(1) + tau
            Delta(1) = -Z(1)/tau
            Delta(2) = Z(2)/(del-tau)
         ELSE
            b = -del + Rho*(Z(1)*Z(1)+Z(2)*Z(2))
            c = Rho*Z(2)*Z(2)*del
            IF ( b>ZERO ) THEN
               tau = -TWO*c/(b+SQRT(b*b+FOUR*c))
            ELSE
               tau = (b-SQRT(b*b+FOUR*c))/TWO
            ENDIF
            Dlam = D(2) + tau
            Delta(1) = -Z(1)/(del+tau)
            Delta(2) = -Z(2)/tau
         ENDIF
         temp = SQRT(Delta(1)*Delta(1)+Delta(2)*Delta(2))
         Delta(1) = Delta(1)/temp
         Delta(2) = Delta(2)/temp
      ELSE
!
!     Now I=2
!
         b = -del + Rho*(Z(1)*Z(1)+Z(2)*Z(2))
         c = Rho*Z(2)*Z(2)*del
         IF ( b>ZERO ) THEN
            tau = (b+SQRT(b*b+FOUR*c))/TWO
         ELSE
            tau = TWO*c/(-b+SQRT(b*b+FOUR*c))
         ENDIF
         Dlam = D(2) + tau
         Delta(1) = -Z(1)/(del+tau)
         Delta(2) = -Z(2)/tau
         temp = SQRT(Delta(1)*Delta(1)+Delta(2)*Delta(2))
         Delta(1) = Delta(1)/temp
         Delta(2) = Delta(2)/temp
      ENDIF
!
!     End OF SLAED5
!
      END SUBROUTINE SLAED5
