!*==ztrt06.f90  processed by SPAG 7.51RB at 17:35 on  4 Mar 2022
!> \brief \b ztrt06
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ztrt06( RCOND, RCONDC, UPLO, DIAG, N, A, LDA, RWORK,
!                          RAT )
!
!       .. Scalar Arguments ..
!       CHARACTER          DIAG, UPLO
!       INTEGER            LDA, N
!       DOUBLE PRECISION   RAT, RCOND, RCONDC
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   RWORK( * )
!       COMPLEX*16         A( LDA, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ztrt06 computes a test ratio comparing RCOND (the reciprocal
!> condition number of a triangular matrix A) and RCONDC, the estimate
!> computed by ZTRCON.  Information about the triangular matrix A is
!> used if one estimate is zero and the other is non-zero to decide if
!> underflow in the estimate is justified.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] RCOND
!> \verbatim
!>          RCOND is DOUBLE PRECISION
!>          The estimate of the reciprocal condition number obtained by
!>          forming the explicit inverse of the matrix A and computing
!>          RCOND = 1/( norm(A) * norm(inv(A)) ).
!> \endverbatim
!>
!> \param[in] RCONDC
!> \verbatim
!>          RCONDC is DOUBLE PRECISION
!>          The estimate of the reciprocal condition number computed by
!>          ZTRCON.
!> \endverbatim
!>
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER
!>          Specifies whether the matrix A is upper or lower triangular.
!>          = 'U':  Upper triangular
!>          = 'L':  Lower triangular
!> \endverbatim
!>
!> \param[in] DIAG
!> \verbatim
!>          DIAG is CHARACTER
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
!> \param[out] RWORK
!> \verbatim
!>          RWORK is DOUBLE PRECISION array, dimension (N)
!> \endverbatim
!>
!> \param[out] RAT
!> \verbatim
!>          RAT is DOUBLE PRECISION
!>          The test ratio.  If both RCOND and RCONDC are nonzero,
!>             RAT = MAX( RCOND, RCONDC )/MIN( RCOND, RCONDC ) - 1.
!>          If RAT = 0, the two estimates are exactly the same.
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
      SUBROUTINE ztrt06(Rcond,Rcondc,Uplo,Diag,N,A,Lda,Rwork,Rat)
      IMPLICIT NONE
!*--ZTRT06125
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Diag , Uplo
      INTEGER Lda , N
      DOUBLE PRECISION Rat , Rcond , Rcondc
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION Rwork(*)
      COMPLEX*16 A(Lda,*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO , ONE
      PARAMETER (ZERO=0.0D+0,ONE=1.0D+0)
!     ..
!     .. Local Scalars ..
      DOUBLE PRECISION anorm , bignum , eps , rmax , rmin
!     ..
!     .. External Functions ..
      DOUBLE PRECISION DLAMCH , ZLANTR
EXTERNAL    DLAMCH,   ZLANTR
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX , MIN
!     ..
!     .. Executable Statements ..
!
      eps = DLAMCH('Epsilon')
      rmax = MAX(Rcond,Rcondc)
      rmin = MIN(Rcond,Rcondc)
!
!     Do the easy cases first.
!
      IF ( rmin<ZERO ) THEN
!
!        Invalid value for RCOND or RCONDC, return 1/EPS.
!
         Rat = ONE/eps
!
      ELSEIF ( rmin>ZERO ) THEN
!
!        Both estimates are positive, return RMAX/RMIN - 1.
!
         Rat = rmax/rmin - ONE
!
      ELSEIF ( rmax==ZERO ) THEN
!
!        Both estimates zero.
!
         Rat = ZERO
!
      ELSE
!
!        One estimate is zero, the other is non-zero.  If the matrix is
!        ill-conditioned, return the nonzero estimate multiplied by
!        1/EPS; if the matrix is badly scaled, return the nonzero
!        estimate multiplied by BIGNUM/TMAX, where TMAX is the maximum
!        element in absolute value in A.
!
         bignum = ONE/DLAMCH('Safe minimum')
         anorm = ZLANTR('M',Uplo,Diag,N,N,A,Lda,Rwork)
!
         Rat = rmax*(MIN(bignum/MAX(ONE,anorm),ONE/eps))
      ENDIF
!
!
!     End of ztrt06
!
      END SUBROUTINE ztrt06
