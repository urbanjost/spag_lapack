!*==cppt01.f90  processed by SPAG 7.51RB at 17:34 on  4 Mar 2022
!> \brief \b cppt01
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE cppt01( UPLO, N, A, AFAC, RWORK, RESID )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            N
!       REAL               RESID
!       ..
!       .. Array Arguments ..
!       REAL               RWORK( * )
!       COMPLEX            A( * ), AFAC( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> cppt01 reconstructs a Hermitian positive definite packed matrix A
!> from its L*L' or U'*U factorization and computes the residual
!>    norm( L*L' - A ) / ( N * norm(A) * EPS ) or
!>    norm( U'*U - A ) / ( N * norm(A) * EPS ),
!> where EPS is the machine epsilon, L' is the conjugate transpose of
!> L, and U' is the conjugate transpose of U.
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
!>          A is COMPLEX array, dimension (N*(N+1)/2)
!>          The original Hermitian matrix A, stored as a packed
!>          triangular matrix.
!> \endverbatim
!>
!> \param[in,out] AFAC
!> \verbatim
!>          AFAC is COMPLEX array, dimension (N*(N+1)/2)
!>          On entry, the factor L or U from the L*L' or U'*U
!>          factorization of A, stored as a packed triangular matrix.
!>          Overwritten with the reconstructed matrix, and then with the
!>          difference L*L' - A (or U'*U - A).
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
!> \ingroup complex_lin
!
!  =====================================================================
      SUBROUTINE cppt01(Uplo,N,A,Afac,Rwork,Resid)
      IMPLICIT NONE
!*--CPPT0199
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Uplo
      INTEGER N
      REAL Resid
!     ..
!     .. Array Arguments ..
      REAL Rwork(*)
      COMPLEX A(*) , Afac(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ZERO , ONE
      PARAMETER (ZERO=0.0E+0,ONE=1.0E+0)
!     ..
!     .. Local Scalars ..
      INTEGER i , k , kc
      REAL anorm , eps , tr
      COMPLEX tc
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      REAL CLANHP , SLAMCH
      COMPLEX CDOTC
EXTERNAL    LSAME,   CLANHP,   SLAMCH,   CDOTC
!     ..
!     .. External Subroutines ..
EXTERNAL    CHPR,   CSCAL,   CTPMV
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC AIMAG , REAL
!     ..
!     .. Executable Statements ..
!
!     Quick exit if N = 0
!
      IF ( N<=0 ) THEN
         Resid = ZERO
         RETURN
      ENDIF
!
!     Exit with RESID = 1/EPS if ANORM = 0.
!
      eps = SLAMCH('Epsilon')
      anorm = CLANHP('1',Uplo,N,A,Rwork)
      IF ( anorm<=ZERO ) THEN
         Resid = ONE/eps
         RETURN
      ENDIF
!
!     Check the imaginary parts of the diagonal elements and return with
!     an error code if any are nonzero.
!
      kc = 1
      IF ( LSAME(Uplo,'U') ) THEN
         DO k = 1 , N
            IF ( AIMAG(Afac(kc))/=ZERO ) THEN
               Resid = ONE/eps
               RETURN
            ENDIF
            kc = kc + k + 1
         ENDDO
      ELSE
         DO k = 1 , N
            IF ( AIMAG(Afac(kc))/=ZERO ) THEN
               Resid = ONE/eps
               RETURN
            ENDIF
            kc = kc + N - k + 1
         ENDDO
      ENDIF
!
!     Compute the product U'*U, overwriting U.
!
      IF ( LSAME(Uplo,'U') ) THEN
         kc = (N*(N-1))/2 + 1
         DO k = N , 1 , -1
!
!           Compute the (K,K) element of the result.
!
            tr = CDOTC(k,Afac(kc),1,Afac(kc),1)
            Afac(kc+k-1) = tr
!
!           Compute the rest of column K.
!
            IF ( k>1 ) THEN
               CALL CTPMV('Upper','Conjugate','Non-unit',k-1,Afac,      &
     &                    Afac(kc),1)
               kc = kc - (k-1)
            ENDIF
         ENDDO
!
!        Compute the difference  L*L' - A
!
         kc = 1
         DO k = 1 , N
            DO i = 1 , k - 1
               Afac(kc+i-1) = Afac(kc+i-1) - A(kc+i-1)
            ENDDO
            Afac(kc+k-1) = Afac(kc+k-1) - REAL(A(kc+k-1))
            kc = kc + k
         ENDDO
!
!     Compute the product L*L', overwriting L.
!
      ELSE
         kc = (N*(N+1))/2
         DO k = N , 1 , -1
!
!           Add a multiple of column K of the factor L to each of
!           columns K+1 through N.
!
            IF ( k<N ) CALL CHPR('Lower',N-k,ONE,Afac(kc+1),1,          &
     &                           Afac(kc+N-k+1))
!
!           Scale column K by the diagonal element.
!
            tc = Afac(kc)
            CALL CSCAL(N-k+1,tc,Afac(kc),1)
!
            kc = kc - (N-k+2)
         ENDDO
!
!        Compute the difference  U'*U - A
!
         kc = 1
         DO k = 1 , N
            Afac(kc) = Afac(kc) - REAL(A(kc))
            DO i = k + 1 , N
               Afac(kc+i-k) = Afac(kc+i-k) - A(kc+i-k)
            ENDDO
            kc = kc + N - k + 1
         ENDDO
      ENDIF
!
!     Compute norm( L*U - A ) / ( N * norm(A) * EPS )
!
      Resid = CLANHP('1',Uplo,N,Afac,Rwork)
!
      Resid = ((Resid/REAL(N))/anorm)/eps
!
!
!     End of cppt01
!
      END SUBROUTINE cppt01
