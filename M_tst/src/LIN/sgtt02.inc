!*==sgtt02.f90  processed by SPAG 7.51RB at 17:34 on  4 Mar 2022
!> \brief \b sgtt02
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE sgtt02( TRANS, N, NRHS, DL, D, DU, X, LDX, B, LDB,
!                          RESID )
!
!       .. Scalar Arguments ..
!       CHARACTER          TRANS
!       INTEGER            LDB, LDX, N, NRHS
!       REAL               RESID
!       ..
!       .. Array Arguments ..
!       REAL               B( LDB, * ), D( * ), DL( * ), DU( * ),
!      $                   X( LDX, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> sgtt02 computes the residual for the solution to a tridiagonal
!> system of equations:
!>    RESID = norm(B - op(A)*X) / (norm(A) * norm(X) * EPS),
!> where EPS is the machine epsilon.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] TRANS
!> \verbatim
!>          TRANS is CHARACTER
!>          Specifies the form of the residual.
!>          = 'N':  B - A * X  (No transpose)
!>          = 'T':  B - A'* X  (Transpose)
!>          = 'C':  B - A'* X  (Conjugate transpose = Transpose)
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGTER
!>          The order of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in] NRHS
!> \verbatim
!>          NRHS is INTEGER
!>          The number of right hand sides, i.e., the number of columns
!>          of the matrices B and X.  NRHS >= 0.
!> \endverbatim
!>
!> \param[in] DL
!> \verbatim
!>          DL is REAL array, dimension (N-1)
!>          The (n-1) sub-diagonal elements of A.
!> \endverbatim
!>
!> \param[in] D
!> \verbatim
!>          D is REAL array, dimension (N)
!>          The diagonal elements of A.
!> \endverbatim
!>
!> \param[in] DU
!> \verbatim
!>          DU is REAL array, dimension (N-1)
!>          The (n-1) super-diagonal elements of A.
!> \endverbatim
!>
!> \param[in] X
!> \verbatim
!>          X is REAL array, dimension (LDX,NRHS)
!>          The computed solution vectors X.
!> \endverbatim
!>
!> \param[in] LDX
!> \verbatim
!>          LDX is INTEGER
!>          The leading dimension of the array X.  LDX >= max(1,N).
!> \endverbatim
!>
!> \param[in,out] B
!> \verbatim
!>          B is REAL array, dimension (LDB,NRHS)
!>          On entry, the right hand side vectors for the system of
!>          linear equations.
!>          On exit, B is overwritten with the difference B - op(A)*X.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of the array B.  LDB >= max(1,N).
!> \endverbatim
!>
!> \param[out] RESID
!> \verbatim
!>          RESID is REAL
!>          norm(B - op(A)*X) / (norm(A) * norm(X) * EPS)
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
      SUBROUTINE sgtt02(Trans,N,Nrhs,Dl,D,Du,X,Ldx,B,Ldb,Resid)
      IMPLICIT NONE
!*--SGTT02127
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Trans
      INTEGER Ldb , Ldx , N , Nrhs
      REAL Resid
!     ..
!     .. Array Arguments ..
      REAL B(Ldb,*) , D(*) , Dl(*) , Du(*) , X(Ldx,*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ONE , ZERO
      PARAMETER (ONE=1.0E+0,ZERO=0.0E+0)
!     ..
!     .. Local Scalars ..
      INTEGER j
      REAL anorm , bnorm , eps , xnorm
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      REAL SASUM , SLAMCH , SLANGT
EXTERNAL    LSAME,   SASUM,   SLAMCH,   SLANGT
!     ..
!     .. External Subroutines ..
EXTERNAL    SLAGTM
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX
!     ..
!     .. Executable Statements ..
!
!     Quick exit if N = 0 or NRHS = 0
!
      Resid = ZERO
      IF ( N<=0 .OR. Nrhs==0 ) RETURN
!
!     Compute the maximum over the number of right hand sides of
!        norm(B - op(A)*X) / ( norm(A) * norm(X) * EPS ).
!
      IF ( LSAME(Trans,'N') ) THEN
         anorm = SLANGT('1',N,Dl,D,Du)
      ELSE
         anorm = SLANGT('I',N,Dl,D,Du)
      ENDIF
!
!     Exit with RESID = 1/EPS if ANORM = 0.
!
      eps = SLAMCH('Epsilon')
      IF ( anorm<=ZERO ) THEN
         Resid = ONE/eps
         RETURN
      ENDIF
!
!     Compute B - op(A)*X.
!
      CALL SLAGTM(Trans,N,Nrhs,-ONE,Dl,D,Du,X,Ldx,ONE,B,Ldb)
!
      DO j = 1 , Nrhs
         bnorm = SASUM(N,B(1,j),1)
         xnorm = SASUM(N,X(1,j),1)
         IF ( xnorm<=ZERO ) THEN
            Resid = ONE/eps
         ELSE
            Resid = MAX(Resid,((bnorm/anorm)/xnorm)/eps)
         ENDIF
      ENDDO
!
!
!     End of sgtt02
!
      END SUBROUTINE sgtt02
