!*==dtpt06.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b DTPT06
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DTPT06( RCOND, RCONDC, UPLO, DIAG, N, AP, WORK, RAT )
!
!       .. Scalar Arguments ..
!       CHARACTER          DIAG, UPLO
!       INTEGER            N
!       DOUBLE PRECISION   RAT, RCOND, RCONDC
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   AP( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DTPT06 computes a test ratio comparing RCOND (the reciprocal
!> condition number of a triangular matrix A) and RCONDC, the estimate
!> computed by DTPCON.  Information about the triangular matrix A is
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
!>          DTPCON.
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
      SUBROUTINE DTPT06(Rcond,Rcondc,Uplo,Diag,N,Ap,Work,Rat)
      IMPLICIT NONE
!*--DTPT06115
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Diag , Uplo
      INTEGER N
      DOUBLE PRECISION Rat , Rcond , Rcondc
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION Ap(*) , Work(*)
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
      DOUBLE PRECISION DLAMCH , DLANTP
      EXTERNAL DLAMCH , DLANTP
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
         anorm = DLANTP('M',Uplo,Diag,N,Ap,Work)
!
         Rat = rmax*(MIN(bignum/MAX(ONE,anorm),ONE/eps))
      ENDIF
!
!
!     End of DTPT06
!
      END SUBROUTINE DTPT06
