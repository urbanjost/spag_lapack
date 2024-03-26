!*==ztrt03.f90  processed by SPAG 7.51RB at 17:35 on  4 Mar 2022
!> \brief \b ztrt03
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ztrt03( UPLO, TRANS, DIAG, N, NRHS, A, LDA, SCALE,
!                          CNORM, TSCAL, X, LDX, B, LDB, WORK, RESID )
!
!       .. Scalar Arguments ..
!       CHARACTER          DIAG, TRANS, UPLO
!       INTEGER            LDA, LDB, LDX, N, NRHS
!       DOUBLE PRECISION   RESID, SCALE, TSCAL
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   CNORM( * )
!       COMPLEX*16         A( LDA, * ), B( LDB, * ), WORK( * ),
!      $                   X( LDX, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ztrt03 computes the residual for the solution to a scaled triangular
!> system of equations A*x = s*b,  A**T *x = s*b,  or  A**H *x = s*b.
!> Here A is a triangular matrix, A**T denotes the transpose of A, A**H
!> denotes the conjugate transpose of A, s is a scalar, and x and b are
!> N by NRHS matrices.  The test ratio is the maximum over the number of
!> right hand sides of
!>    norm(s*b - op(A)*x) / ( norm(op(A)) * norm(x) * EPS ),
!> where op(A) denotes A, A**T, or A**H, and EPS is the machine epsilon.
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
!>          = 'N':  A *x = s*b     (No transpose)
!>          = 'T':  A**T *x = s*b  (Transpose)
!>          = 'C':  A**H *x = s*b  (Conjugate transpose)
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
!> \param[in] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension (LDA,N)
!>          The triangular matrix A.  If UPLO = 'U', the leading n by n
!>          upper triangular part of the array A contains the upper
!>          triangular matrix, and the strictly lower triangular part of
!>          A is not referenced.  If UPLO = 'L', the leading n by n lower
!>          triangular part of the array A contains the lower triangular
!>          matrix, and the strictly upper triangular part of A is not
!>          referenced.  If DIAG = 'U', the diagonal elements of A are
!>          also not referenced and are assumed to be 1.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,N).
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
!>          X is COMPLEX*16 array, dimension (LDX,NRHS)
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
!>          B is COMPLEX*16 array, dimension (LDB,NRHS)
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
!>          WORK is COMPLEX*16 array, dimension (N)
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
!> \ingroup complex16_lin
!
!  =====================================================================
      SUBROUTINE ztrt03(Uplo,Trans,Diag,N,Nrhs,A,Lda,Scale,Cnorm,Tscal, &
     &                  X,Ldx,B,Ldb,Work,Resid)
      IMPLICIT NONE
!*--ZTRT03175
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Diag , Trans , Uplo
      INTEGER Lda , Ldb , Ldx , N , Nrhs
      DOUBLE PRECISION Resid , Scale , Tscal
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION Cnorm(*)
      COMPLEX*16 A(Lda,*) , B(Ldb,*) , Work(*) , X(Ldx,*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ONE , ZERO
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
!     ..
!     .. Local Scalars ..
      INTEGER ix , j
      DOUBLE PRECISION eps , err , smlnum , tnorm , xnorm , xscal
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      INTEGER IZAMAX
      DOUBLE PRECISION DLAMCH
EXTERNAL    LSAME,   IZAMAX,   DLAMCH
!     ..
!     .. External Subroutines ..
EXTERNAL    ZAXPY,   ZCOPY,   ZDSCAL,   ZTRMV
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , DBLE , DCMPLX , MAX
!     ..
!     .. Executable Statements ..
!
!     Quick exit if N = 0
!
      IF ( N<=0 .OR. Nrhs<=0 ) THEN
         Resid = ZERO
         RETURN
      ENDIF
      eps = DLAMCH('Epsilon')
      smlnum = DLAMCH('Safe minimum')
!
!     Compute the norm of the triangular matrix A using the column
!     norms already computed by ZLATRS.
!
      tnorm = ZERO
      IF ( LSAME(Diag,'N') ) THEN
         DO j = 1 , N
            tnorm = MAX(tnorm,Tscal*ABS(A(j,j))+Cnorm(j))
         ENDDO
      ELSE
         DO j = 1 , N
            tnorm = MAX(tnorm,Tscal+Cnorm(j))
         ENDDO
      ENDIF
!
!     Compute the maximum over the number of right hand sides of
!        norm(op(A)*x - s*b) / ( norm(op(A)) * norm(x) * EPS ).
!
      Resid = ZERO
      DO j = 1 , Nrhs
         CALL ZCOPY(N,X(1,j),1,Work,1)
         ix = IZAMAX(N,Work,1)
         xnorm = MAX(ONE,ABS(X(ix,j)))
         xscal = (ONE/xnorm)/DBLE(N)
         CALL ZDSCAL(N,xscal,Work,1)
         CALL ZTRMV(Uplo,Trans,Diag,N,A,Lda,Work,1)
         CALL ZAXPY(N,DCMPLX(-Scale*xscal),B(1,j),1,Work,1)
         ix = IZAMAX(N,Work,1)
         err = Tscal*ABS(Work(ix))
         ix = IZAMAX(N,X(1,j),1)
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
!     End of ztrt03
!
      END SUBROUTINE ztrt03
