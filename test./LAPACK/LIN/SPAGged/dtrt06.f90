!*==dtrt06.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b DTRT06
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DTRT06( RCOND, RCONDC, UPLO, DIAG, N, A, LDA, WORK,
!                          RAT )
!
!       .. Scalar Arguments ..
!       CHARACTER          DIAG, UPLO
!       INTEGER            LDA, N
!       DOUBLE PRECISION   RAT, RCOND, RCONDC
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   A( LDA, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DTRT06 computes a test ratio comparing RCOND (the reciprocal
!> condition number of a triangular matrix A) and RCONDC, the estimate
!> computed by DTRCON.  Information about the triangular matrix A is
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
!>          DTRCON.
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
!>          A is DOUBLE PRECISION array, dimension (LDA,N)
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
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (N)
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
!> \ingroup double_lin
!
!  =====================================================================
      SUBROUTINE DTRT06(Rcond,Rcondc,Uplo,Diag,N,A,Lda,Work,Rat)
      IMPLICIT NONE
!*--DTRT06124
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
      DOUBLE PRECISION A(Lda,*) , Work(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO , ONE
      PARAMETER (ZERO=0.0D+0,ONE=1.0D+0)
!     ..
!     .. Local Scalars ..
      DOUBLE PRECISION anorm , bignum , eps , rmax , rmin , smlnum
!     ..
!     .. External Functions ..
      DOUBLE PRECISION DLAMCH , DLANTR
      EXTERNAL DLAMCH , DLANTR
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX , MIN
!     ..
!     .. External Subroutines ..
      EXTERNAL DLABAD
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
         smlnum = DLAMCH('Safe minimum')
         bignum = ONE/smlnum
         CALL DLABAD(smlnum,bignum)
         anorm = DLANTR('M',Uplo,Diag,N,N,A,Lda,Work)
!
         Rat = rmax*(MIN(bignum/MAX(ONE,anorm),ONE/eps))
      ENDIF
!
!
!     End of DTRT06
!
      END SUBROUTINE DTRT06