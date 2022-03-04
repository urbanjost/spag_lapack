!*==zspt02.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b ZSPT02
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZSPT02( UPLO, N, NRHS, A, X, LDX, B, LDB, RWORK,
!                          RESID )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            LDB, LDX, N, NRHS
!       DOUBLE PRECISION   RESID
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   RWORK( * )
!       COMPLEX*16         A( * ), B( LDB, * ), X( LDX, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZSPT02 computes the residual in the solution of a complex symmetric
!> system of linear equations  A*x = b  when packed storage is used for
!> the coefficient matrix.  The ratio computed is
!>
!>    RESID = norm( B - A*X ) / ( norm(A) * norm(X) * EPS).
!>
!> where EPS is the machine precision.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          Specifies whether the upper or lower triangular part of the
!>          complex symmetric matrix A is stored:
!>          = 'U':  Upper triangular
!>          = 'L':  Lower triangular
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of rows and columns of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in] NRHS
!> \verbatim
!>          NRHS is INTEGER
!>          The number of columns of B, the matrix of right hand sides.
!>          NRHS >= 0.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension (N*(N+1)/2)
!>          The original complex symmetric matrix A, stored as a packed
!>          triangular matrix.
!> \endverbatim
!>
!> \param[in] X
!> \verbatim
!>          X is COMPLEX*16 array, dimension (LDX,NRHS)
!>          The computed solution vectors for the system of linear
!>          equations.
!> \endverbatim
!>
!> \param[in] LDX
!> \verbatim
!>          LDX is INTEGER
!>          The leading dimension of the array X.   LDX >= max(1,N).
!> \endverbatim
!>
!> \param[in,out] B
!> \verbatim
!>          B is COMPLEX*16 array, dimension (LDB,NRHS)
!>          On entry, the right hand side vectors for the system of
!>          linear equations.
!>          On exit, B is overwritten with the difference B - A*X.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of the array B.  LDB >= max(1,N).
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is DOUBLE PRECISION array, dimension (N)
!> \endverbatim
!>
!> \param[out] RESID
!> \verbatim
!>          RESID is DOUBLE PRECISION
!>          The maximum over the number of right hand sides of
!>          norm(B - A*X) / ( norm(A) * norm(X) * EPS ).
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
      SUBROUTINE ZSPT02(Uplo,N,Nrhs,A,X,Ldx,B,Ldb,Rwork,Resid)
      IMPLICIT NONE
!*--ZSPT02126
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Uplo
      INTEGER Ldb , Ldx , N , Nrhs
      DOUBLE PRECISION Resid
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION Rwork(*)
      COMPLEX*16 A(*) , B(Ldb,*) , X(Ldx,*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO , ONE
      PARAMETER (ZERO=0.0D+0,ONE=1.0D+0)
      COMPLEX*16 CONE
      PARAMETER (CONE=(1.0D+0,0.0D+0))
!     ..
!     .. Local Scalars ..
      INTEGER j
      DOUBLE PRECISION anorm , bnorm , eps , xnorm
!     ..
!     .. External Functions ..
      DOUBLE PRECISION DLAMCH , DZASUM , ZLANSP
      EXTERNAL DLAMCH , DZASUM , ZLANSP
!     ..
!     .. External Subroutines ..
      EXTERNAL ZSPMV
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX
!     ..
!     .. Executable Statements ..
!
!     Quick exit if N = 0 or NRHS = 0
!
      IF ( N<=0 .OR. Nrhs<=0 ) THEN
         Resid = ZERO
         RETURN
      ENDIF
!
!     Exit with RESID = 1/EPS if ANORM = 0.
!
      eps = DLAMCH('Epsilon')
      anorm = ZLANSP('1',Uplo,N,A,Rwork)
      IF ( anorm<=ZERO ) THEN
         Resid = ONE/eps
         RETURN
      ENDIF
!
!     Compute  B - A*X  for the matrix of right hand sides B.
!
      DO j = 1 , Nrhs
         CALL ZSPMV(Uplo,N,-CONE,A,X(1,j),1,CONE,B(1,j),1)
      ENDDO
!
!     Compute the maximum over the number of right hand sides of
!        norm( B - A*X ) / ( norm(A) * norm(X) * EPS ) .
!
      Resid = ZERO
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
!     End of ZSPT02
!
      END SUBROUTINE ZSPT02
