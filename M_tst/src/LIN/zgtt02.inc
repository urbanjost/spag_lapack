!*==zgtt02.f90  processed by SPAG 7.51RB at 17:35 on  4 Mar 2022
!> \brief \b zgtt02
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE zgtt02( TRANS, N, NRHS, DL, D, DU, X, LDX, B, LDB,
!                          RESID )
!
!       .. Scalar Arguments ..
!       CHARACTER          TRANS
!       INTEGER            LDB, LDX, N, NRHS
!       DOUBLE PRECISION   RESID
!       ..
!       .. Array Arguments ..
!       COMPLEX*16         B( LDB, * ), D( * ), DL( * ), DU( * ),
!      $                   X( LDX, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> zgtt02 computes the residual for the solution to a tridiagonal
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
!>          = 'N':  B - A * X     (No transpose)
!>          = 'T':  B - A**T * X  (Transpose)
!>          = 'C':  B - A**H * X  (Conjugate transpose)
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
!>          DL is COMPLEX*16 array, dimension (N-1)
!>          The (n-1) sub-diagonal elements of A.
!> \endverbatim
!>
!> \param[in] D
!> \verbatim
!>          D is COMPLEX*16 array, dimension (N)
!>          The diagonal elements of A.
!> \endverbatim
!>
!> \param[in] DU
!> \verbatim
!>          DU is COMPLEX*16 array, dimension (N-1)
!>          The (n-1) super-diagonal elements of A.
!> \endverbatim
!>
!> \param[in] X
!> \verbatim
!>          X is COMPLEX*16 array, dimension (LDX,NRHS)
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
!>          B is COMPLEX*16 array, dimension (LDB,NRHS)
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
!>          RESID is DOUBLE PRECISION
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
!> \ingroup complex16_lin
!
!  =====================================================================
      SUBROUTINE zgtt02(Trans,N,Nrhs,Dl,D,Du,X,Ldx,B,Ldb,Resid)
      IMPLICIT NONE
!*--ZGTT02127
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Trans
      INTEGER Ldb , Ldx , N , Nrhs
      DOUBLE PRECISION Resid
!     ..
!     .. Array Arguments ..
      COMPLEX*16 B(Ldb,*) , D(*) , Dl(*) , Du(*) , X(Ldx,*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ONE , ZERO
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
!     ..
!     .. Local Scalars ..
      INTEGER j
      DOUBLE PRECISION anorm , bnorm , eps , xnorm
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      DOUBLE PRECISION DLAMCH , DZASUM , ZLANGT
EXTERNAL    LSAME,   DLAMCH,   DZASUM,   ZLANGT
!     ..
!     .. External Subroutines ..
EXTERNAL    ZLAGTM
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
         anorm = ZLANGT('1',N,Dl,D,Du)
      ELSE
         anorm = ZLANGT('I',N,Dl,D,Du)
      ENDIF
!
!     Exit with RESID = 1/EPS if ANORM = 0.
!
      eps = DLAMCH('Epsilon')
      IF ( anorm<=ZERO ) THEN
         Resid = ONE/eps
         RETURN
      ENDIF
!
!     Compute B - op(A)*X.
!
      CALL ZLAGTM(Trans,N,Nrhs,-ONE,Dl,D,Du,X,Ldx,ONE,B,Ldb)
!
      DO j = 1 , Nrhs
         bnorm = DZASUM(N,B(1,j),1)
         xnorm = DZASUM(N,X(1,j),1)
         IF ( xnorm<=ZERO ) THEN
            Resid = ONE/eps
         ELSE
            Resid = MAX(Resid,((bnorm/anorm)/xnorm)/eps)
         ENDIF
      ENDDO
!
!
!     End of zgtt02
!
      END SUBROUTINE zgtt02
