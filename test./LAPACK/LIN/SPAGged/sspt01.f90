!*==sspt01.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b SSPT01
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE SSPT01( UPLO, N, A, AFAC, IPIV, C, LDC, RWORK, RESID )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            LDC, N
!       REAL               RESID
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * )
!       REAL               A( * ), AFAC( * ), C( LDC, * ), RWORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SSPT01 reconstructs a symmetric indefinite packed matrix A from its
!> block L*D*L' or U*D*U' factorization and computes the residual
!>      norm( C - A ) / ( N * norm(A) * EPS ),
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
!>          symmetric matrix A is stored:
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
!> \param[in] A
!> \verbatim
!>          A is REAL array, dimension (N*(N+1)/2)
!>          The original symmetric matrix A, stored as a packed
!>          triangular matrix.
!> \endverbatim
!>
!> \param[in] AFAC
!> \verbatim
!>          AFAC is REAL array, dimension (N*(N+1)/2)
!>          The factored form of the matrix A, stored as a packed
!>          triangular matrix.  AFAC contains the block diagonal matrix D
!>          and the multipliers used to obtain the factor L or U from the
!>          block L*D*L' or U*D*U' factorization as computed by SSPTRF.
!> \endverbatim
!>
!> \param[in] IPIV
!> \verbatim
!>          IPIV is INTEGER array, dimension (N)
!>          The pivot indices from SSPTRF.
!> \endverbatim
!>
!> \param[out] C
!> \verbatim
!>          C is REAL array, dimension (LDC,N)
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
!> \ingroup single_lin
!
!  =====================================================================
      SUBROUTINE SSPT01(Uplo,N,A,Afac,Ipiv,C,Ldc,Rwork,Resid)
      IMPLICIT NONE
!*--SSPT01114
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
      REAL A(*) , Afac(*) , C(Ldc,*) , Rwork(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ZERO , ONE
      PARAMETER (ZERO=0.0E+0,ONE=1.0E+0)
!     ..
!     .. Local Scalars ..
      INTEGER i , info , j , jc
      REAL anorm , eps
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      REAL SLAMCH , SLANSP , SLANSY
      EXTERNAL LSAME , SLAMCH , SLANSP , SLANSY
!     ..
!     .. External Subroutines ..
      EXTERNAL SLAVSP , SLASET
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
      anorm = SLANSP('1',Uplo,N,A,Rwork)
!
!     Initialize C to the identity matrix.
!
      CALL SLASET('Full',N,N,ZERO,ONE,C,Ldc)
!
!     Call SLAVSP to form the product D * U' (or D * L' ).
!
      CALL SLAVSP(Uplo,'Transpose','Non-unit',N,N,Afac,Ipiv,C,Ldc,info)
!
!     Call SLAVSP again to multiply by U ( or L ).
!
      CALL SLAVSP(Uplo,'No transpose','Unit',N,N,Afac,Ipiv,C,Ldc,info)
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
      Resid = SLANSY('1',Uplo,N,C,Ldc,Rwork)
!
      IF ( anorm<=ZERO ) THEN
         IF ( Resid/=ZERO ) Resid = ONE/eps
      ELSE
         Resid = ((Resid/REAL(N))/anorm)/eps
      ENDIF
!
!
!     End of SSPT01
!
      END SUBROUTINE SSPT01
