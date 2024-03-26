!*==dtpt03.f90  processed by SPAG 7.51RB at 17:34 on  4 Mar 2022
!> \brief \b dtpt03
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE dtpt03( UPLO, TRANS, DIAG, N, NRHS, AP, SCALE, CNORM,
!                          TSCAL, X, LDX, B, LDB, WORK, RESID )
!
!       .. Scalar Arguments ..
!       CHARACTER          DIAG, TRANS, UPLO
!       INTEGER            LDB, LDX, N, NRHS
!       DOUBLE PRECISION   RESID, SCALE, TSCAL
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   AP( * ), B( LDB, * ), CNORM( * ), WORK( * ),
!      $                   X( LDX, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> dtpt03 computes the residual for the solution to a scaled triangular
!> system of equations A*x = s*b  or  A'*x = s*b  when the triangular
!> matrix A is stored in packed format.  Here A' is the transpose of A,
!> s is a scalar, and x and b are N by NRHS matrices.  The test ratio is
!> the maximum over the number of right hand sides of
!>    norm(s*b - op(A)*x) / ( norm(op(A)) * norm(x) * EPS ),
!> where op(A) denotes A or A' and EPS is the machine epsilon.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          Specifies whether the matrix A is upper or lower triangular.
!>          = 'U':  Upper triangular
!>          = 'L':  Lower triangular
!> \endverbatim
!>
!> \param[in] TRANS
!> \verbatim
!>          TRANS is CHARACTER*1
!>          Specifies the operation applied to A.
!>          = 'N':  A *x = s*b  (No transpose)
!>          = 'T':  A'*x = s*b  (Transpose)
!>          = 'C':  A'*x = s*b  (Conjugate transpose = Transpose)
!> \endverbatim
!>
!> \param[in] DIAG
!> \verbatim
!>          DIAG is CHARACTER*1
!>          Specifies whether or not the matrix A is unit triangular.
!>          = 'N':  Non-unit triangular
!>          = 'U':  Unit triangular
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in] NRHS
!> \verbatim
!>          NRHS is INTEGER
!>          The number of right hand sides, i.e., the number of columns
!>          of the matrices X and B.  NRHS >= 0.
!> \endverbatim
!>
!> \param[in] AP
!> \verbatim
!>          AP is DOUBLE PRECISION array, dimension (N*(N+1)/2)
!>          The upper or lower triangular matrix A, packed columnwise in
!>          a linear array.  The j-th column of A is stored in the array
!>          AP as follows:
!>          if UPLO = 'U', AP((j-1)*j/2 + i) = A(i,j) for 1<=i<=j;
!>          if UPLO = 'L',
!>             AP((j-1)*(n-j) + j*(j+1)/2 + i-j) = A(i,j) for j<=i<=n.
!> \endverbatim
!>
!> \param[in] SCALE
!> \verbatim
!>          SCALE is DOUBLE PRECISION
!>          The scaling factor s used in solving the triangular system.
!> \endverbatim
!>
!> \param[in] CNORM
!> \verbatim
!>          CNORM is DOUBLE PRECISION array, dimension (N)
!>          The 1-norms of the columns of A, not counting the diagonal.
!> \endverbatim
!>
!> \param[in] TSCAL
!> \verbatim
!>          TSCAL is DOUBLE PRECISION
!>          The scaling factor used in computing the 1-norms in CNORM.
!>          CNORM actually contains the column norms of TSCAL*A.
!> \endverbatim
!>
!> \param[in] X
!> \verbatim
!>          X is DOUBLE PRECISION array, dimension (LDX,NRHS)
!>          The computed solution vectors for the system of linear
!>          equations.
!> \endverbatim
!>
!> \param[in] LDX
!> \verbatim
!>          LDX is INTEGER
!>          The leading dimension of the array X.  LDX >= max(1,N).
!> \endverbatim
!>
!> \param[in] B
!> \verbatim
!>          B is DOUBLE PRECISION array, dimension (LDB,NRHS)
!>          The right hand side vectors for the system of linear
!>          equations.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of the array B.  LDB >= max(1,N).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (N)
!> \endverbatim
!>
!> \param[out] RESID
!> \verbatim
!>          RESID is DOUBLE PRECISION
!>          The maximum over the number of right hand sides of
!>          norm(op(A)*x - s*b) / ( norm(op(A)) * norm(x) * EPS ).
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
      SUBROUTINE dtpt03(Uplo,Trans,Diag,N,Nrhs,Ap,Scale,Cnorm,Tscal,X,  &
     &                  Ldx,B,Ldb,Work,Resid)
      IMPLICIT NONE
!*--DTPT03165
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Diag , Trans , Uplo
      INTEGER Ldb , Ldx , N , Nrhs
      DOUBLE PRECISION Resid , Scale , Tscal
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION Ap(*) , B(Ldb,*) , Cnorm(*) , Work(*) , X(Ldx,*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ONE , ZERO
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
!     ..
!     .. Local Scalars ..
      INTEGER ix , j , jj
      DOUBLE PRECISION bignum , eps , err , smlnum , tnorm , xnorm ,    &
     &                 xscal
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      INTEGER IDAMAX
      DOUBLE PRECISION DLAMCH
EXTERNAL    LSAME,   IDAMAX,   DLAMCH
!     ..
!     .. External Subroutines ..
EXTERNAL    DAXPY,   DCOPY,   DLABAD,   DSCAL,   DTPMV
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , DBLE , MAX
!     ..
!     .. Executable Statements ..
!
!     Quick exit if N = 0.
!
      IF ( N<=0 .OR. Nrhs<=0 ) THEN
         Resid = ZERO
         RETURN
      ENDIF
      eps = DLAMCH('Epsilon')
      smlnum = DLAMCH('Safe minimum')
      bignum = ONE/smlnum
      CALL DLABAD(smlnum,bignum)
!
!     Compute the norm of the triangular matrix A using the column
!     norms already computed by DLATPS.
!
      tnorm = ZERO
      IF ( .NOT.(LSAME(Diag,'N')) ) THEN
         DO j = 1 , N
            tnorm = MAX(tnorm,Tscal+Cnorm(j))
         ENDDO
      ELSEIF ( LSAME(Uplo,'U') ) THEN
         jj = 1
         DO j = 1 , N
            tnorm = MAX(tnorm,Tscal*ABS(Ap(jj))+Cnorm(j))
            jj = jj + j + 1
         ENDDO
      ELSE
         jj = 1
         DO j = 1 , N
            tnorm = MAX(tnorm,Tscal*ABS(Ap(jj))+Cnorm(j))
            jj = jj + N - j + 1
         ENDDO
      ENDIF
!
!     Compute the maximum over the number of right hand sides of
!        norm(op(A)*x - s*b) / ( norm(op(A)) * norm(x) * EPS ).
!
      Resid = ZERO
      DO j = 1 , Nrhs
         CALL DCOPY(N,X(1,j),1,Work,1)
         ix = IDAMAX(N,Work,1)
         xnorm = MAX(ONE,ABS(X(ix,j)))
         xscal = (ONE/xnorm)/DBLE(N)
         CALL DSCAL(N,xscal,Work,1)
         CALL DTPMV(Uplo,Trans,Diag,N,Ap,Work,1)
         CALL DAXPY(N,-Scale*xscal,B(1,j),1,Work,1)
         ix = IDAMAX(N,Work,1)
         err = Tscal*ABS(Work(ix))
         ix = IDAMAX(N,X(1,j),1)
         xnorm = ABS(X(ix,j))
         IF ( err*smlnum<=xnorm ) THEN
            IF ( xnorm>ZERO ) err = err/xnorm
         ELSE
            IF ( err>ZERO ) err = ONE/eps
         ENDIF
         IF ( err*smlnum<=tnorm ) THEN
            IF ( tnorm>ZERO ) err = err/tnorm
         ELSE
            IF ( err>ZERO ) err = ONE/eps
         ENDIF
         Resid = MAX(Resid,err)
      ENDDO
!
!
!     End of dtpt03
!
      END SUBROUTINE dtpt03
