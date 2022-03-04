!*==zhet01_3.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b ZHET01_3
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZHET01_3( UPLO, N, A, LDA, AFAC, LDAFAC, E, IPIV, C,
!                            LDC, RWORK, RESID )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            LDA, LDAFAC, LDC, N
!       DOUBLE PRECISION   RESID
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * )
!       DOUBLE PRECISION   RWORK( * )
!       COMPLEX*16         A( LDA, * ), AFAC( LDAFAC, * ), C( LDC, * ),
!                          E( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZHET01_3 reconstructs a Hermitian indefinite matrix A from its
!> block L*D*L' or U*D*U' factorization computed by ZHETRF_RK
!> (or ZHETRF_BK) and computes the residual
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
!>          The number of rows and columns of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension (LDA,N)
!>          The original Hermitian matrix A.
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
!>          AFAC is COMPLEX*16 array, dimension (LDAFAC,N)
!>          Diagonal of the block diagonal matrix D and factors U or L
!>          as computed by ZHETRF_RK and ZHETRF_BK:
!>            a) ONLY diagonal elements of the Hermitian block diagonal
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
!>          E is COMPLEX*16 array, dimension (N)
!>          On entry, contains the superdiagonal (or subdiagonal)
!>          elements of the Hermitian block diagonal matrix D
!>          with 1-by-1 or 2-by-2 diagonal blocks, where
!>          If UPLO = 'U': E(i) = D(i-1,i),i=2:N, E(1) not referenced;
!>          If UPLO = 'L': E(i) = D(i+1,i),i=1:N-1, E(N) not referenced.
!> \endverbatim
!>
!> \param[in] IPIV
!> \verbatim
!>          IPIV is INTEGER array, dimension (N)
!>          The pivot indices from ZHETRF_RK (or ZHETRF_BK).
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
!> \date June 2017
!
!> \ingroup complex16_lin
!
!  =====================================================================
      SUBROUTINE ZHET01_3(Uplo,N,A,Lda,Afac,Ldafac,E,Ipiv,C,Ldc,Rwork,  &
     &                    Resid)
      IMPLICIT NONE
!*--ZHET01_3145
!
!  -- LAPACK test routine (version 3.7.1) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     June 2017
!
!     .. Scalar Arguments ..
      CHARACTER Uplo
      INTEGER Lda , Ldafac , Ldc , N
      DOUBLE PRECISION Resid
!     ..
!     .. Array Arguments ..
      INTEGER Ipiv(*)
      DOUBLE PRECISION Rwork(*)
      COMPLEX*16 A(Lda,*) , Afac(Ldafac,*) , C(Ldc,*) , E(*)
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
      INTEGER i , info , j
      DOUBLE PRECISION anorm , eps
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      DOUBLE PRECISION ZLANHE , DLAMCH
      EXTERNAL LSAME , ZLANHE , DLAMCH
!     ..
!     .. External Subroutines ..
      EXTERNAL ZLASET , ZLAVHE_ROOK , ZSYCONVF_ROOK
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC DIMAG , DBLE
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
      CALL ZSYCONVF_ROOK(Uplo,'R',N,Afac,Ldafac,E,Ipiv,info)
!
!     1) Determine EPS and the norm of A.
!
      eps = DLAMCH('Epsilon')
      anorm = ZLANHE('1',Uplo,N,A,Lda,Rwork)
!
!     Check the imaginary parts of the diagonal elements and return with
!     an error code if any are nonzero.
!
      DO j = 1 , N
         IF ( DIMAG(Afac(j,j))/=ZERO ) THEN
            Resid = ONE/eps
            RETURN
         ENDIF
      ENDDO
!
!     2) Initialize C to the identity matrix.
!
      CALL ZLASET('Full',N,N,CZERO,CONE,C,Ldc)
!
!     3) Call ZLAVHE_ROOK to form the product D * U' (or D * L' ).
!
      CALL ZLAVHE_ROOK(Uplo,'Conjugate','Non-unit',N,N,Afac,Ldafac,Ipiv,&
     &                 C,Ldc,info)
!
!     4) Call ZLAVHE_RK again to multiply by U (or L ).
!
      CALL ZLAVHE_ROOK(Uplo,'No transpose','Unit',N,N,Afac,Ldafac,Ipiv, &
     &                 C,Ldc,info)
!
!     5) Compute the difference  C - A .
!
      IF ( LSAME(Uplo,'U') ) THEN
         DO j = 1 , N
            DO i = 1 , j - 1
               C(i,j) = C(i,j) - A(i,j)
            ENDDO
            C(j,j) = C(j,j) - DBLE(A(j,j))
         ENDDO
      ELSE
         DO j = 1 , N
            C(j,j) = C(j,j) - DBLE(A(j,j))
            DO i = j + 1 , N
               C(i,j) = C(i,j) - A(i,j)
            ENDDO
         ENDDO
      ENDIF
!
!     6) Compute norm( C - A ) / ( N * norm(A) * EPS )
!
      Resid = ZLANHE('1',Uplo,N,C,Ldc,Rwork)
!
      IF ( anorm<=ZERO ) THEN
         IF ( Resid/=ZERO ) Resid = ONE/eps
      ELSE
         Resid = ((Resid/DBLE(N))/anorm)/eps
      ENDIF
!
!     b) Convert to factor of L (or U)
!
      CALL ZSYCONVF_ROOK(Uplo,'C',N,Afac,Ldafac,E,Ipiv,info)
!
!
!     End of ZHET01_3
!
      END SUBROUTINE ZHET01_3
