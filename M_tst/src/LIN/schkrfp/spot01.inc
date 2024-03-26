!*==spot01.f90  processed by SPAG 7.51RB at 17:34 on  4 Mar 2022
!> \brief \b spot01
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE spot01( UPLO, N, A, LDA, AFAC, LDAFAC, RWORK, RESID )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            LDA, LDAFAC, N
!       REAL               RESID
!       ..
!       .. Array Arguments ..
!       REAL               A( LDA, * ), AFAC( LDAFAC, * ), RWORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> spot01 reconstructs a symmetric positive definite matrix  A  from
!> its L*L' or U'*U factorization and computes the residual
!>    norm( L*L' - A ) / ( N * norm(A) * EPS ) or
!>    norm( U'*U - A ) / ( N * norm(A) * EPS ),
!> where EPS is the machine epsilon.
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
!>          A is REAL array, dimension (LDA,N)
!>          The original symmetric matrix A.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,N)
!> \endverbatim
!>
!> \param[in,out] AFAC
!> \verbatim
!>          AFAC is REAL array, dimension (LDAFAC,N)
!>          On entry, the factor L or U from the L*L' or U'*U
!>          factorization of A.
!>          Overwritten with the reconstructed matrix, and then with the
!>          difference L*L' - A (or U'*U - A).
!> \endverbatim
!>
!> \param[in] LDAFAC
!> \verbatim
!>          LDAFAC is INTEGER
!>          The leading dimension of the array AFAC.  LDAFAC >= max(1,N).
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
!>          If UPLO = 'L', norm(L*L' - A) / ( N * norm(A) * EPS )
!>          If UPLO = 'U', norm(U'*U - A) / ( N * norm(A) * EPS )
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
      SUBROUTINE spot01(Uplo,N,A,Lda,Afac,Ldafac,Rwork,Resid)
      IMPLICIT NONE
!*--SPOT01108
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Uplo
      INTEGER Lda , Ldafac , N
      REAL Resid
!     ..
!     .. Array Arguments ..
      REAL A(Lda,*) , Afac(Ldafac,*) , Rwork(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ZERO , ONE
      PARAMETER (ZERO=0.0E+0,ONE=1.0E+0)
!     ..
!     .. Local Scalars ..
      INTEGER i , j , k
      REAL anorm , eps , t
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      REAL SDOT , SLAMCH , SLANSY
      EXTERNAL LSAME , SDOT , SLAMCH , SLANSY
!     ..
!     .. External Subroutines ..
      EXTERNAL SSCAL , SSYR , STRMV
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
!     Exit with RESID = 1/EPS if ANORM = 0.
!
      eps = SLAMCH('Epsilon')
      anorm = SLANSY('1',Uplo,N,A,Lda,Rwork)
      IF ( anorm<=ZERO ) THEN
         Resid = ONE/eps
         RETURN
      ENDIF
!
!     Compute the product U'*U, overwriting U.
!
      IF ( LSAME(Uplo,'U') ) THEN
         DO k = N , 1 , -1
!
!           Compute the (K,K) element of the result.
!
            t = SDOT(k,Afac(1,k),1,Afac(1,k),1)
            Afac(k,k) = t
!
!           Compute the rest of column K.
!
            CALL STRMV('Upper','Transpose','Non-unit',k-1,Afac,Ldafac,  &
     &                 Afac(1,k),1)
!
         ENDDO
!
!     Compute the product L*L', overwriting L.
!
      ELSE
         DO k = N , 1 , -1
!
!           Add a multiple of column K of the factor L to each of
!           columns K+1 through N.
!
            IF ( k+1<=N ) CALL SSYR('Lower',N-k,ONE,Afac(k+1,k),1,      &
     &                              Afac(k+1,k+1),Ldafac)
!
!           Scale column K by the diagonal element.
!
            t = Afac(k,k)
            CALL SSCAL(N-k+1,t,Afac(k,k),1)
!
         ENDDO
      ENDIF
!
!     Compute the difference  L*L' - A (or U'*U - A).
!
      IF ( LSAME(Uplo,'U') ) THEN
         DO j = 1 , N
            DO i = 1 , j
               Afac(i,j) = Afac(i,j) - A(i,j)
            ENDDO
         ENDDO
      ELSE
         DO j = 1 , N
            DO i = j , N
               Afac(i,j) = Afac(i,j) - A(i,j)
            ENDDO
         ENDDO
      ENDIF
!
!     Compute norm( L*U - A ) / ( N * norm(A) * EPS )
!
      Resid = SLANSY('1',Uplo,N,Afac,Ldafac,Rwork)
!
      Resid = ((Resid/REAL(N))/anorm)/eps
!
!
!     End of spot01
!
      END SUBROUTINE spot01
