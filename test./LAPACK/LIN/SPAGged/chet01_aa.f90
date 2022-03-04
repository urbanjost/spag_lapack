!*==chet01_aa.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b CHET01_AA
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE CHET01_AA( UPLO, N, A, LDA, AFAC, LDAFAC, IPIV,
!                             C, LDC, RWORK, RESID )
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
!> CHET01_AA reconstructs a hermitian indefinite matrix A from its
!> block L*D*L' or U*D*U' factorization and computes the residual
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
!>          hermitian matrix A is stored:
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
!>          The original hermitian matrix A.
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
!>          as computed by CHETRF.
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
!>          The pivot indices from CHETRF.
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
!>          RWORK is COMPLEX array, dimension (N)
!> \endverbatim
!>
!> \param[out] RESID
!> \verbatim
!>          RESID is COMPLEX
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
!
!> \ingroup complex_lin
!
!  =====================================================================
      SUBROUTINE CHET01_AA(Uplo,N,A,Lda,Afac,Ldafac,Ipiv,C,Ldc,Rwork,   &
     &                     Resid)
      IMPLICIT NONE
!*--CHET01_AA129
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
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
      COMPLEX CZERO , CONE
      PARAMETER (CZERO=(0.0E+0,0.0E+0),CONE=(1.0E+0,0.0E+0))
      REAL ZERO , ONE
      PARAMETER (ZERO=0.0E+0,ONE=1.0E+0)
!     ..
!     .. Local Scalars ..
      INTEGER i , j
      REAL anorm , eps
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      REAL SLAMCH , CLANHE
      EXTERNAL LSAME , SLAMCH , CLANHE
!     ..
!     .. External Subroutines ..
      EXTERNAL CLASET , CLAVHE
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
      eps = SLAMCH('Epsilon')
      anorm = CLANHE('1',Uplo,N,A,Lda,Rwork)
!
!     Initialize C to the tridiagonal matrix T.
!
      CALL CLASET('Full',N,N,CZERO,CZERO,C,Ldc)
      CALL CLACPY('F',1,N,Afac(1,1),Ldafac+1,C(1,1),Ldc+1)
      IF ( N>1 ) THEN
         IF ( LSAME(Uplo,'U') ) THEN
            CALL CLACPY('F',1,N-1,Afac(1,2),Ldafac+1,C(1,2),Ldc+1)
            CALL CLACPY('F',1,N-1,Afac(1,2),Ldafac+1,C(2,1),Ldc+1)
            CALL CLACGV(N-1,C(2,1),Ldc+1)
         ELSE
            CALL CLACPY('F',1,N-1,Afac(2,1),Ldafac+1,C(1,2),Ldc+1)
            CALL CLACPY('F',1,N-1,Afac(2,1),Ldafac+1,C(2,1),Ldc+1)
            CALL CLACGV(N-1,C(1,2),Ldc+1)
         ENDIF
!
!        Call CTRMM to form the product U' * D (or L * D ).
!
         IF ( LSAME(Uplo,'U') ) THEN
            CALL CTRMM('Left',Uplo,'Conjugate transpose','Unit',N-1,N,  &
     &                 CONE,Afac(1,2),Ldafac,C(2,1),Ldc)
         ELSE
            CALL CTRMM('Left',Uplo,'No transpose','Unit',N-1,N,CONE,    &
     &                 Afac(2,1),Ldafac,C(2,1),Ldc)
         ENDIF
!
!        Call CTRMM again to multiply by U (or L ).
!
         IF ( LSAME(Uplo,'U') ) THEN
            CALL CTRMM('Right',Uplo,'No transpose','Unit',N,N-1,CONE,   &
     &                 Afac(1,2),Ldafac,C(1,2),Ldc)
         ELSE
            CALL CTRMM('Right',Uplo,'Conjugate transpose','Unit',N,N-1, &
     &                 CONE,Afac(2,1),Ldafac,C(1,2),Ldc)
         ENDIF
      ENDIF
!
!     Apply hermitian pivots
!
      DO j = N , 1 , -1
         i = Ipiv(j)
         IF ( i/=j ) CALL CSWAP(N,C(j,1),Ldc,C(i,1),Ldc)
      ENDDO
      DO j = N , 1 , -1
         i = Ipiv(j)
         IF ( i/=j ) CALL CSWAP(N,C(1,j),1,C(1,i),1)
      ENDDO
!
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
      Resid = CLANHE('1',Uplo,N,C,Ldc,Rwork)
!
      IF ( anorm<=ZERO ) THEN
         IF ( Resid/=ZERO ) Resid = ONE/eps
      ELSE
         Resid = ((Resid/DBLE(N))/anorm)/eps
      ENDIF
!
!
!     End of CHET01_AA
!
      END SUBROUTINE CHET01_AA
