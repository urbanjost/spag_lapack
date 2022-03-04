!*==stbt06.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b STBT06
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE STBT06( RCOND, RCONDC, UPLO, DIAG, N, KD, AB, LDAB,
!                          WORK, RAT )
!
!       .. Scalar Arguments ..
!       CHARACTER          DIAG, UPLO
!       INTEGER            KD, LDAB, N
!       REAL               RAT, RCOND, RCONDC
!       ..
!       .. Array Arguments ..
!       REAL               AB( LDAB, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> STBT06 computes a test ratio comparing RCOND (the reciprocal
!> condition number of a triangular matrix A) and RCONDC, the estimate
!> computed by STBCON.  Information about the triangular matrix A is
!> used if one estimate is zero and the other is non-zero to decide if
!> underflow in the estimate is justified.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] RCOND
!> \verbatim
!>          RCOND is REAL
!>          The estimate of the reciprocal condition number obtained by
!>          forming the explicit inverse of the matrix A and computing
!>          RCOND = 1/( norm(A) * norm(inv(A)) ).
!> \endverbatim
!>
!> \param[in] RCONDC
!> \verbatim
!>          RCONDC is REAL
!>          The estimate of the reciprocal condition number computed by
!>          STBCON.
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
!> \param[in] KD
!> \verbatim
!>          KD is INTEGER
!>          The number of superdiagonals or subdiagonals of the
!>          triangular band matrix A.  KD >= 0.
!> \endverbatim
!>
!> \param[in] AB
!> \verbatim
!>          AB is REAL array, dimension (LDAB,N)
!>          The upper or lower triangular band matrix A, stored in the
!>          first kd+1 rows of the array. The j-th column of A is stored
!>          in the j-th column of the array AB as follows:
!>          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j;
!>          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd).
!> \endverbatim
!>
!> \param[in] LDAB
!> \verbatim
!>          LDAB is INTEGER
!>          The leading dimension of the array AB.  LDAB >= KD+1.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is REAL array, dimension (N)
!> \endverbatim
!>
!> \param[out] RAT
!> \verbatim
!>          RAT is REAL
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
!> \ingroup single_lin
!
!  =====================================================================
      SUBROUTINE STBT06(Rcond,Rcondc,Uplo,Diag,N,Kd,Ab,Ldab,Work,Rat)
      IMPLICIT NONE
!*--STBT06128
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Diag , Uplo
      INTEGER Kd , Ldab , N
      REAL Rat , Rcond , Rcondc
!     ..
!     .. Array Arguments ..
      REAL Ab(Ldab,*) , Work(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ZERO , ONE
      PARAMETER (ZERO=0.0E+0,ONE=1.0E+0)
!     ..
!     .. Local Scalars ..
      REAL anorm , bignum , eps , rmax , rmin , smlnum
!     ..
!     .. External Functions ..
      REAL SLAMCH , SLANTB
      EXTERNAL SLAMCH , SLANTB
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX , MIN
!     ..
!     .. External Subroutines ..
      EXTERNAL SLABAD
!     ..
!     .. Executable Statements ..
!
      eps = SLAMCH('Epsilon')
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
         smlnum = SLAMCH('Safe minimum')
         bignum = ONE/smlnum
         CALL SLABAD(smlnum,bignum)
         anorm = SLANTB('M',Uplo,Diag,N,Kd,Ab,Ldab,Work)
!
         Rat = rmax*(MIN(bignum/MAX(ONE,anorm),ONE/eps))
      ENDIF
!
!
!     End of STBT06
!
      END SUBROUTINE STBT06
