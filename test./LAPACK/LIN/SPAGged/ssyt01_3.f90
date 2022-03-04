!*==ssyt01_3.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b SSYT01_3
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE SSYT01_3( UPLO, N, A, LDA, AFAC, LDAFAC, E, IPIV, C,
!                            LDC, RWORK, RESID )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            LDA, LDAFAC, LDC, N
!       DOUBLE PRECISION   RESID
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * )
!       DOUBLE PRECISION   A( LDA, * ), AFAC( LDAFAC, * ), C( LDC, * ),
!      $                   E( * ), RWORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SSYT01_3 reconstructs a symmetric indefinite matrix A from its
!> block L*D*L' or U*D*U' factorization computed by SSYTRF_RK
!> (or SSYTRF_BK) and computes the residual
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
!>          A is DOUBLE PRECISION array, dimension (LDA,N)
!>          The original symmetric matrix A.
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
!>          AFAC is DOUBLE PRECISION array, dimension (LDAFAC,N)
!>          Diagonal of the block diagonal matrix D and factors U or L
!>          as computed by SSYTRF_RK and SSYTRF_BK:
!>            a) ONLY diagonal elements of the symmetric block diagonal
!>               matrix D on the diagonal of A, i.e. D(k,k) = A(k,k);
!>               (superdiagonal (or subdiagonal) elements of D
!>                should be provided on entry in array E), and
!>            b) If UPLO = 'U': factor U in the superdiagonal part of A.
!>               If UPLO = 'L': factor L in the subdiagonal part of A.
!> \endverbatim
!>
!> \param[in] LDAFAC
!> \verbatim
!>          LDAFAC is INTEGER
!>          The leading dimension of the array AFAC.
!>          LDAFAC >= max(1,N).
!> \endverbatim
!>
!> \param[in] E
!> \verbatim
!>          E is DOUBLE PRECISION array, dimension (N)
!>          On entry, contains the superdiagonal (or subdiagonal)
!>          elements of the symmetric block diagonal matrix D
!>          with 1-by-1 or 2-by-2 diagonal blocks, where
!>          If UPLO = 'U': E(i) = D(i-1,i),i=2:N, E(1) not referenced;
!>          If UPLO = 'L': E(i) = D(i+1,i),i=1:N-1, E(N) not referenced.
!> \endverbatim
!>
!> \param[in] IPIV
!> \verbatim
!>          IPIV is INTEGER array, dimension (N)
!>          The pivot indices from SSYTRF_RK (or SSYTRF_BK).
!> \endverbatim
!>
!> \param[out] C
!> \verbatim
!>          C is DOUBLE PRECISION array, dimension (LDC,N)
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
!> \date June 2017
!
!> \ingroup single_lin
!
!  =====================================================================
      SUBROUTINE SSYT01_3(Uplo,N,A,Lda,Afac,Ldafac,E,Ipiv,C,Ldc,Rwork,  &
     &                    Resid)
      IMPLICIT NONE
!*--SSYT01_3144
!
!  -- LAPACK test routine (version 3.7.1) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     June 2017
!
!     .. Scalar Arguments ..
      CHARACTER Uplo
      INTEGER Lda , Ldafac , Ldc , N
      REAL Resid
!     ..
!     .. Array Arguments ..
      INTEGER Ipiv(*)
      REAL A(Lda,*) , Afac(Ldafac,*) , C(Ldc,*) , E(*) , Rwork(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ZERO , ONE
      PARAMETER (ZERO=0.0E+0,ONE=1.0E+0)
!     ..
!     .. Local Scalars ..
      INTEGER i , info , j
      REAL anorm , eps
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      REAL SLAMCH , SLANSY
      EXTERNAL LSAME , SLAMCH , SLANSY
!     ..
!     .. External Subroutines ..
      EXTERNAL SLASET , SLAVSY_ROOK , SSYCONVF_ROOK
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
!     a) Revert to multiplyers of L
!
      CALL SSYCONVF_ROOK(Uplo,'R',N,Afac,Ldafac,E,Ipiv,info)
!
!     1) Determine EPS and the norm of A.
!
      eps = SLAMCH('Epsilon')
      anorm = SLANSY('1',Uplo,N,A,Lda,Rwork)
!
!     2) Initialize C to the identity matrix.
!
      CALL SLASET('Full',N,N,ZERO,ONE,C,Ldc)
!
!     3) Call SLAVSY_ROOK to form the product D * U' (or D * L' ).
!
      CALL SLAVSY_ROOK(Uplo,'Transpose','Non-unit',N,N,Afac,Ldafac,Ipiv,&
     &                 C,Ldc,info)
!
!     4) Call SLAVSY_ROOK again to multiply by U (or L ).
!
      CALL SLAVSY_ROOK(Uplo,'No transpose','Unit',N,N,Afac,Ldafac,Ipiv, &
     &                 C,Ldc,info)
!
!     5) Compute the difference  C - A.
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
!     6) Compute norm( C - A ) / ( N * norm(A) * EPS )
!
      Resid = SLANSY('1',Uplo,N,C,Ldc,Rwork)
!
      IF ( anorm<=ZERO ) THEN
         IF ( Resid/=ZERO ) Resid = ONE/eps
      ELSE
         Resid = ((Resid/REAL(N))/anorm)/eps
      ENDIF
 
!
!     b) Convert to factor of L (or U)
!
      CALL SSYCONVF_ROOK(Uplo,'C',N,Afac,Ldafac,E,Ipiv,info)
!
!
!     End of SSYT01_3
!
      END SUBROUTINE SSYT01_3
