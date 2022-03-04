!*==cspt01.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b CSPT01
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE CSPT01( UPLO, N, A, AFAC, IPIV, C, LDC, RWORK, RESID )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            LDC, N
!       REAL               RESID
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * )
!       REAL               RWORK( * )
!       COMPLEX            A( * ), AFAC( * ), C( LDC, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CSPT01 reconstructs a symmetric indefinite packed matrix A from its
!> diagonal pivoting factorization A = U*D*U' or A = L*D*L' and computes
!> the residual
!>    norm( C - A ) / ( N * norm(A) * EPS ),
!> where C is the reconstructed matrix and EPS is the machine epsilon.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          Specifies whether the upper or lower triangular part of the
!>          Hermitian matrix A is stored:
!>          = 'U':  Upper triangular
!>          = 'L':  Lower triangular
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
!>          A is COMPLEX array, dimension (N*(N+1)/2)
!>          The original symmetric matrix A, stored as a packed
!>          triangular matrix.
!> \endverbatim
!>
!> \param[in] AFAC
!> \verbatim
!>          AFAC is COMPLEX array, dimension (N*(N+1)/2)
!>          The factored form of the matrix A, stored as a packed
!>          triangular matrix.  AFAC contains the block diagonal matrix D
!>          and the multipliers used to obtain the factor L or U from the
!>          L*D*L' or U*D*U' factorization as computed by CSPTRF.
!> \endverbatim
!>
!> \param[in] IPIV
!> \verbatim
!>          IPIV is INTEGER array, dimension (N)
!>          The pivot indices from CSPTRF.
!> \endverbatim
!>
!> \param[out] C
!> \verbatim
!>          C is COMPLEX array, dimension (LDC,N)
!> \endverbatim
!>
!> \param[in] LDC
!> \verbatim
!>          LDC is INTEGER
!>          The leading dimension of the array C.  LDC >= max(1,N).
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is REAL array, dimension (N)
!> \endverbatim
!>
!> \param[out] RESID
!> \verbatim
!>          RESID is REAL
!>          If UPLO = 'L', norm(L*D*L' - A) / ( N * norm(A) * EPS )
!>          If UPLO = 'U', norm(U*D*U' - A) / ( N * norm(A) * EPS )
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
!> \ingroup complex_lin
!
!  =====================================================================
      SUBROUTINE CSPT01(Uplo,N,A,Afac,Ipiv,C,Ldc,Rwork,Resid)
      IMPLICIT NONE
!*--CSPT01116
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Uplo
      INTEGER Ldc , N
      REAL Resid
!     ..
!     .. Array Arguments ..
      INTEGER Ipiv(*)
      REAL Rwork(*)
      COMPLEX A(*) , Afac(*) , C(Ldc,*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ZERO , ONE
      PARAMETER (ZERO=0.0E+0,ONE=1.0E+0)
      COMPLEX CZERO , CONE
      PARAMETER (CZERO=(0.0E+0,0.0E+0),CONE=(1.0E+0,0.0E+0))
!     ..
!     .. Local Scalars ..
      INTEGER i , info , j , jc
      REAL anorm , eps
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      REAL CLANSP , CLANSY , SLAMCH
      EXTERNAL LSAME , CLANSP , CLANSY , SLAMCH
!     ..
!     .. External Subroutines ..
      EXTERNAL CLAVSP , CLASET
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC REAL
!     ..
!     .. Executable Statements ..
!
!     Quick exit if N = 0.
!
      IF ( N<=0 ) THEN
         Resid = ZERO
         RETURN
      ENDIF
!
!     Determine EPS and the norm of A.
!
      eps = SLAMCH('Epsilon')
      anorm = CLANSP('1',Uplo,N,A,Rwork)
!
!     Initialize C to the identity matrix.
!
      CALL CLASET('Full',N,N,CZERO,CONE,C,Ldc)
!
!     Call CLAVSP to form the product D * U' (or D * L' ).
!
      CALL CLAVSP(Uplo,'Transpose','Non-unit',N,N,Afac,Ipiv,C,Ldc,info)
!
!     Call CLAVSP again to multiply by U ( or L ).
!
      CALL CLAVSP(Uplo,'No transpose','Unit',N,N,Afac,Ipiv,C,Ldc,info)
!
!     Compute the difference  C - A .
!
      IF ( LSAME(Uplo,'U') ) THEN
         jc = 0
         DO j = 1 , N
            DO i = 1 , j
               C(i,j) = C(i,j) - A(jc+i)
            ENDDO
            jc = jc + j
         ENDDO
      ELSE
         jc = 1
         DO j = 1 , N
            DO i = j , N
               C(i,j) = C(i,j) - A(jc+i-j)
            ENDDO
            jc = jc + N - j + 1
         ENDDO
      ENDIF
!
!     Compute norm( C - A ) / ( N * norm(A) * EPS )
!
      Resid = CLANSY('1',Uplo,N,C,Ldc,Rwork)
!
      IF ( anorm<=ZERO ) THEN
         IF ( Resid/=ZERO ) Resid = ONE/eps
      ELSE
         Resid = ((Resid/REAL(N))/anorm)/eps
      ENDIF
!
!
!     End of CSPT01
!
      END SUBROUTINE CSPT01
