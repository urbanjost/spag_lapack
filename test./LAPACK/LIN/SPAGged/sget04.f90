!*==sget04.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b SGET04
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE SGET04( N, NRHS, X, LDX, XACT, LDXACT, RCOND, RESID )
!
!       .. Scalar Arguments ..
!       INTEGER            LDX, LDXACT, N, NRHS
!       REAL               RCOND, RESID
!       ..
!       .. Array Arguments ..
!       REAL               X( LDX, * ), XACT( LDXACT, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SGET04 computes the difference between a computed solution and the
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
!>          X is REAL array, dimension (LDX,NRHS)
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
!>          XACT is REAL array, dimension( LDX, NRHS )
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
!>          RCOND is REAL
!>          The reciprocal of the condition number of the coefficient
!>          matrix in the system of equations.
!> \endverbatim
!>
!> \param[out] RESID
!> \verbatim
!>          RESID is REAL
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
!> \ingroup single_lin
!
!  =====================================================================
      SUBROUTINE SGET04(N,Nrhs,X,Ldx,Xact,Ldxact,Rcond,Resid)
      IMPLICIT NONE
!*--SGET04106
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Ldx , Ldxact , N , Nrhs
      REAL Rcond , Resid
!     ..
!     .. Array Arguments ..
      REAL X(Ldx,*) , Xact(Ldxact,*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ZERO
      PARAMETER (ZERO=0.0E+0)
!     ..
!     .. Local Scalars ..
      INTEGER i , ix , j
      REAL diffnm , eps , xnorm
!     ..
!     .. External Functions ..
      INTEGER ISAMAX
      REAL SLAMCH
      EXTERNAL ISAMAX , SLAMCH
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
      eps = SLAMCH('Epsilon')
      IF ( Rcond<ZERO ) THEN
         Resid = 1.0/eps
         RETURN
      ENDIF
!
!     Compute the maximum of
!        norm(X - XACT) / ( norm(XACT) * EPS )
!     over all the vectors X and XACT .
!
      Resid = ZERO
      DO j = 1 , Nrhs
         ix = ISAMAX(N,Xact(1,j),1)
         xnorm = ABS(Xact(ix,j))
         diffnm = ZERO
         DO i = 1 , N
            diffnm = MAX(diffnm,ABS(X(i,j)-Xact(i,j)))
         ENDDO
         IF ( xnorm<=ZERO ) THEN
            IF ( diffnm>ZERO ) Resid = 1.0/eps
         ELSE
            Resid = MAX(Resid,(diffnm/xnorm)*Rcond)
         ENDIF
      ENDDO
      IF ( Resid*eps<1.0 ) Resid = Resid/eps
!
!
!     End of SGET04
!
      END SUBROUTINE SGET04