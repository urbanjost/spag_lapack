!*==csyt01.f90  processed by SPAG 7.51RB at 17:34 on  4 Mar 2022
!> \brief \b csyt01
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE csyt01( UPLO, N, A, LDA, AFAC, LDAFAC, IPIV, C, LDC,
!                          RWORK, RESID )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            LDA, LDAFAC, LDC, N
!       REAL               RESID
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * )
!       REAL               RWORK( * )
!       COMPLEX            A( LDA, * ), AFAC( LDAFAC, * ), C( LDC, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> csyt01 reconstructs a complex symmetric indefinite matrix A from its
!> block L*D*L' or U*D*U' factorization and computes the residual
!>    norm( C - A ) / ( N * norm(A) * EPS ),
!> where C is the reconstructed matrix, EPS is the machine epsilon,
!> L' is the transpose of L, and U' is the transpose of U.
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
!> \param[in] A
!> \verbatim
!>          A is COMPLEX array, dimension (LDA,N)
!>          The original complex symmetric matrix A.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,N)
!> \endverbatim
!>
!> \param[in] AFAC
!> \verbatim
!>          AFAC is COMPLEX array, dimension (LDAFAC,N)
!>          The factored form of the matrix A.  AFAC contains the block
!>          diagonal matrix D and the multipliers used to obtain the
!>          factor L or U from the block L*D*L' or U*D*U' factorization
!>          as computed by CSYTRF.
!> \endverbatim
!>
!> \param[in] LDAFAC
!> \verbatim
!>          LDAFAC is INTEGER
!>          The leading dimension of the array AFAC.  LDAFAC >= max(1,N).
!> \endverbatim
!>
!> \param[in] IPIV
!> \verbatim
!>          IPIV is INTEGER array, dimension (N)
!>          The pivot indices from CSYTRF.
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
!> \date November 2013
!
!> \ingroup complex_lin
!
!  =====================================================================
      SUBROUTINE csyt01(Uplo,N,A,Lda,Afac,Ldafac,Ipiv,C,Ldc,Rwork,Resid)
      IMPLICIT NONE
!*--CSYT01128
!
!  -- LAPACK test routine (version 3.5.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2013
!
!     .. Scalar Arguments ..
      CHARACTER Uplo
      INTEGER Lda , Ldafac , Ldc , N
      REAL Resid
!     ..
!     .. Array Arguments ..
      INTEGER Ipiv(*)
      REAL Rwork(*)
      COMPLEX A(Lda,*) , Afac(Ldafac,*) , C(Ldc,*)
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
      INTEGER i , info , j
      REAL anorm , eps
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      REAL CLANSY , SLAMCH
EXTERNAL    LSAME,   CLANSY,   SLAMCH
!     ..
!     .. External Subroutines ..
EXTERNAL    CLASET
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
      anorm = CLANSY('1',Uplo,N,A,Lda,Rwork)
!
!     Initialize C to the identity matrix.
!
      CALL CLASET('Full',N,N,CZERO,CONE,C,Ldc)
!
!     Call clavsy to form the product D * U' (or D * L' ).
!
      CALL clavsy(Uplo,'Transpose','Non-unit',N,N,Afac,Ldafac,Ipiv,C,   &
     &            Ldc,info)
!
!     Call clavsy again to multiply by U (or L ).
!
      CALL clavsy(Uplo,'No transpose','Unit',N,N,Afac,Ldafac,Ipiv,C,Ldc,&
     &            info)
!
!     Compute the difference  C - A .
!
      IF ( LSAME(Uplo,'U') ) THEN
         DO j = 1 , N
            DO i = 1 , j
               C(i,j) = C(i,j) - A(i,j)
            ENDDO
         ENDDO
      ELSE
         DO j = 1 , N
            DO i = j , N
               C(i,j) = C(i,j) - A(i,j)
            ENDDO
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
!     End of csyt01
!
      END SUBROUTINE csyt01