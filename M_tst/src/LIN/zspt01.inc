!*==zspt01.f90  processed by SPAG 7.51RB at 17:35 on  4 Mar 2022
!> \brief \b zspt01
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE zspt01( UPLO, N, A, AFAC, IPIV, C, LDC, RWORK, RESID )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            LDC, N
!       DOUBLE PRECISION   RESID
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * )
!       DOUBLE PRECISION   RWORK( * )
!       COMPLEX*16         A( * ), AFAC( * ), C( LDC, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> zspt01 reconstructs a symmetric indefinite packed matrix A from its
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
!>          A is COMPLEX*16 array, dimension (N*(N+1)/2)
!>          The original symmetric matrix A, stored as a packed
!>          triangular matrix.
!> \endverbatim
!>
!> \param[in] AFAC
!> \verbatim
!>          AFAC is COMPLEX*16 array, dimension (N*(N+1)/2)
!>          The factored form of the matrix A, stored as a packed
!>          triangular matrix.  AFAC contains the block diagonal matrix D
!>          and the multipliers used to obtain the factor L or U from the
!>          L*D*L' or U*D*U' factorization as computed by ZSPTRF.
!> \endverbatim
!>
!> \param[in] IPIV
!> \verbatim
!>          IPIV is INTEGER array, dimension (N)
!>          The pivot indices from ZSPTRF.
!> \endverbatim
!>
!> \param[out] C
!> \verbatim
!>          C is COMPLEX*16 array, dimension (LDC,N)
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
!>          RWORK is DOUBLE PRECISION array, dimension (N)
!> \endverbatim
!>
!> \param[out] RESID
!> \verbatim
!>          RESID is DOUBLE PRECISION
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
!> \ingroup complex16_lin
!
!  =====================================================================
      SUBROUTINE zspt01(Uplo,N,A,Afac,Ipiv,C,Ldc,Rwork,Resid)
      IMPLICIT NONE
!*--ZSPT01116
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Uplo
      INTEGER Ldc , N
      DOUBLE PRECISION Resid
!     ..
!     .. Array Arguments ..
      INTEGER Ipiv(*)
      DOUBLE PRECISION Rwork(*)
      COMPLEX*16 A(*) , Afac(*) , C(Ldc,*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO , ONE
      PARAMETER (ZERO=0.0D+0,ONE=1.0D+0)
      COMPLEX*16 CZERO , CONE
      PARAMETER (CZERO=(0.0D+0,0.0D+0),CONE=(1.0D+0,0.0D+0))
!     ..
!     .. Local Scalars ..
      INTEGER i , info , j , jc
      DOUBLE PRECISION anorm , eps
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      DOUBLE PRECISION DLAMCH , ZLANSP , ZLANSY
EXTERNAL    LSAME,   DLAMCH,   ZLANSP,   ZLANSY
!     ..
!     .. External Subroutines ..
EXTERNAL    ZLASET
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC DBLE
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
      eps = DLAMCH('Epsilon')
      anorm = ZLANSP('1',Uplo,N,A,Rwork)
!
!     Initialize C to the identity matrix.
!
      CALL ZLASET('Full',N,N,CZERO,CONE,C,Ldc)
!
!     Call zlavsp to form the product D * U' (or D * L' ).
!
      CALL zlavsp(Uplo,'Transpose','Non-unit',N,N,Afac,Ipiv,C,Ldc,info)
!
!     Call zlavsp again to multiply by U ( or L ).
!
      CALL zlavsp(Uplo,'No transpose','Unit',N,N,Afac,Ipiv,C,Ldc,info)
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
      Resid = ZLANSY('1',Uplo,N,C,Ldc,Rwork)
!
      IF ( anorm<=ZERO ) THEN
         IF ( Resid/=ZERO ) Resid = ONE/eps
      ELSE
         Resid = ((Resid/DBLE(N))/anorm)/eps
      ENDIF
!
!
!     End of zspt01
!
      END SUBROUTINE zspt01
