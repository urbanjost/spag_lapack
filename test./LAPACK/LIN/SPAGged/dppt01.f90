!*==dppt01.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b DPPT01
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DPPT01( UPLO, N, A, AFAC, RWORK, RESID )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            N
!       DOUBLE PRECISION   RESID
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   A( * ), AFAC( * ), RWORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DPPT01 reconstructs a symmetric positive definite packed matrix A
!> from its L*L' or U'*U factorization and computes the residual
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
!>          A is DOUBLE PRECISION array, dimension (N*(N+1)/2)
!>          The original symmetric matrix A, stored as a packed
!>          triangular matrix.
!> \endverbatim
!>
!> \param[in,out] AFAC
!> \verbatim
!>          AFAC is DOUBLE PRECISION array, dimension (N*(N+1)/2)
!>          On entry, the factor L or U from the L*L' or U'*U
!>          factorization of A, stored as a packed triangular matrix.
!>          Overwritten with the reconstructed matrix, and then with the
!>          difference L*L' - A (or U'*U - A).
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
!> \ingroup double_lin
!
!  =====================================================================
      SUBROUTINE DPPT01(Uplo,N,A,Afac,Rwork,Resid)
      IMPLICIT NONE
!*--DPPT0197
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Uplo
      INTEGER N
      DOUBLE PRECISION Resid
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION A(*) , Afac(*) , Rwork(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO , ONE
      PARAMETER (ZERO=0.0D+0,ONE=1.0D+0)
!     ..
!     .. Local Scalars ..
      INTEGER i , k , kc , npp
      DOUBLE PRECISION anorm , eps , t
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      DOUBLE PRECISION DDOT , DLAMCH , DLANSP
      EXTERNAL LSAME , DDOT , DLAMCH , DLANSP
!     ..
!     .. External Subroutines ..
      EXTERNAL DSCAL , DSPR , DTPMV
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC DBLE
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
      eps = DLAMCH('Epsilon')
      anorm = DLANSP('1',Uplo,N,A,Rwork)
      IF ( anorm<=ZERO ) THEN
         Resid = ONE/eps
         RETURN
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
            t = DDOT(k,Afac(kc),1,Afac(kc),1)
            Afac(kc+k-1) = t
!
!           Compute the rest of column K.
!
            IF ( k>1 ) THEN
               CALL DTPMV('Upper','Transpose','Non-unit',k-1,Afac,      &
     &                    Afac(kc),1)
               kc = kc - (k-1)
            ENDIF
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
            IF ( k<N ) CALL DSPR('Lower',N-k,ONE,Afac(kc+1),1,          &
     &                           Afac(kc+N-k+1))
!
!           Scale column K by the diagonal element.
!
            t = Afac(kc)
            CALL DSCAL(N-k+1,t,Afac(kc),1)
!
            kc = kc - (N-k+2)
         ENDDO
      ENDIF
!
!     Compute the difference  L*L' - A (or U'*U - A).
!
      npp = N*(N+1)/2
      DO i = 1 , npp
         Afac(i) = Afac(i) - A(i)
      ENDDO
!
!     Compute norm( L*U - A ) / ( N * norm(A) * EPS )
!
      Resid = DLANSP('1',Uplo,N,Afac,Rwork)
!
      Resid = ((Resid/DBLE(N))/anorm)/eps
!
!
!     End of DPPT01
!
      END SUBROUTINE DPPT01
