!*==dget04.f90  processed by SPAG 7.51RB at 17:34 on  4 Mar 2022
!> \brief \b dget04
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE dget04( N, NRHS, X, LDX, XACT, LDXACT, RCOND, RESID )
!
!       .. Scalar Arguments ..
!       INTEGER            LDX, LDXACT, N, NRHS
!       DOUBLE PRECISION   RCOND, RESID
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   X( LDX, * ), XACT( LDXACT, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> dget04 computes the difference between a computed solution and the
!> true solution to a system of linear equations.
!>
!> RESID =  ( norm(X-XACT) * RCOND ) / ( norm(XACT) * EPS ),
!> where RCOND is the reciprocal of the condition number and EPS is the
!> machine epsilon.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of rows of the matrices X and XACT.  N >= 0.
!> \endverbatim
!>
!> \param[in] NRHS
!> \verbatim
!>          NRHS is INTEGER
!>          The number of columns of the matrices X and XACT.  NRHS >= 0.
!> \endverbatim
!>
!> \param[in] X
!> \verbatim
!>          X is DOUBLE PRECISION array, dimension (LDX,NRHS)
!>          The computed solution vectors.  Each vector is stored as a
!>          column of the matrix X.
!> \endverbatim
!>
!> \param[in] LDX
!> \verbatim
!>          LDX is INTEGER
!>          The leading dimension of the array X.  LDX >= max(1,N).
!> \endverbatim
!>
!> \param[in] XACT
!> \verbatim
!>          XACT is DOUBLE PRECISION array, dimension( LDX, NRHS )
!>          The exact solution vectors.  Each vector is stored as a
!>          column of the matrix XACT.
!> \endverbatim
!>
!> \param[in] LDXACT
!> \verbatim
!>          LDXACT is INTEGER
!>          The leading dimension of the array XACT.  LDXACT >= max(1,N).
!> \endverbatim
!>
!> \param[in] RCOND
!> \verbatim
!>          RCOND is DOUBLE PRECISION
!>          The reciprocal of the condition number of the coefficient
!>          matrix in the system of equations.
!> \endverbatim
!>
!> \param[out] RESID
!> \verbatim
!>          RESID is DOUBLE PRECISION
!>          The maximum over the NRHS solution vectors of
!>          ( norm(X-XACT) * RCOND ) / ( norm(XACT) * EPS )
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
!> \ingroup double_lin
!
!  =====================================================================
      SUBROUTINE dget04(N,Nrhs,X,Ldx,Xact,Ldxact,Rcond,Resid)
      IMPLICIT NONE
!*--DGET04106
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Ldx , Ldxact , N , Nrhs
      DOUBLE PRECISION Rcond , Resid
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION X(Ldx,*) , Xact(Ldxact,*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D+0)
!     ..
!     .. Local Scalars ..
      INTEGER i , ix , j
      DOUBLE PRECISION diffnm , eps , xnorm
!     ..
!     .. External Functions ..
      INTEGER IDAMAX
      DOUBLE PRECISION DLAMCH
EXTERNAL    IDAMAX,   DLAMCH
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , MAX
!     ..
!     .. Executable Statements ..
!
!     Quick exit if N = 0 or NRHS = 0.
!
      IF ( N<=0 .OR. Nrhs<=0 ) THEN
         Resid = ZERO
         RETURN
      ENDIF
!
!     Exit with RESID = 1/EPS if RCOND is invalid.
!
      eps = DLAMCH('Epsilon')
      IF ( Rcond<ZERO ) THEN
         Resid = 1.0D0/eps
         RETURN
      ENDIF
!
!     Compute the maximum of
!        norm(X - XACT) / ( norm(XACT) * EPS )
!     over all the vectors X and XACT .
!
      Resid = ZERO
      DO j = 1 , Nrhs
         ix = IDAMAX(N,Xact(1,j),1)
         xnorm = ABS(Xact(ix,j))
         diffnm = ZERO
         DO i = 1 , N
            diffnm = MAX(diffnm,ABS(X(i,j)-Xact(i,j)))
         ENDDO
         IF ( xnorm<=ZERO ) THEN
            IF ( diffnm>ZERO ) Resid = 1.0D0/eps
         ELSE
            Resid = MAX(Resid,(diffnm/xnorm)*Rcond)
         ENDIF
      ENDDO
      IF ( Resid*eps<1.0D0 ) Resid = Resid/eps
!
!
!     End of dget04
!
      END SUBROUTINE dget04
