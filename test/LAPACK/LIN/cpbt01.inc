!*==cpbt01.f90  processed by SPAG 7.51RB at 17:34 on  4 Mar 2022
!> \brief \b cpbt01
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE cpbt01( UPLO, N, KD, A, LDA, AFAC, LDAFAC, RWORK,
!                          RESID )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            KD, LDA, LDAFAC, N
!       REAL               RESID
!       ..
!       .. Array Arguments ..
!       REAL               RWORK( * )
!       COMPLEX            A( LDA, * ), AFAC( LDAFAC, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> cpbt01 reconstructs a Hermitian positive definite band matrix A from
!> its L*L' or U'*U factorization and computes the residual
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
!> \param[in] KD
!> \verbatim
!>          KD is INTEGER
!>          The number of super-diagonals of the matrix A if UPLO = 'U',
!>          or the number of sub-diagonals if UPLO = 'L'.  KD >= 0.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is COMPLEX array, dimension (LDA,N)
!>          The original Hermitian band matrix A.  If UPLO = 'U', the
!>          upper triangular part of A is stored as a band matrix; if
!>          UPLO = 'L', the lower triangular part of A is stored.  The
!>          columns of the appropriate triangle are stored in the columns
!>          of A and the diagonals of the triangle are stored in the rows
!>          of A.  See CPBTRF for further details.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER.
!>          The leading dimension of the array A.  LDA >= max(1,KD+1).
!> \endverbatim
!>
!> \param[in] AFAC
!> \verbatim
!>          AFAC is COMPLEX array, dimension (LDAFAC,N)
!>          The factored form of the matrix A.  AFAC contains the factor
!>          L or U from the L*L' or U'*U factorization in band storage
!>          format, as computed by CPBTRF.
!> \endverbatim
!>
!> \param[in] LDAFAC
!> \verbatim
!>          LDAFAC is INTEGER
!>          The leading dimension of the array AFAC.
!>          LDAFAC >= max(1,KD+1).
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
      SUBROUTINE cpbt01(Uplo,N,Kd,A,Lda,Afac,Ldafac,Rwork,Resid)
      IMPLICIT NONE
!*--CPBT01123
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Uplo
      INTEGER Kd , Lda , Ldafac , N
      REAL Resid
!     ..
!     .. Array Arguments ..
      REAL Rwork(*)
      COMPLEX A(Lda,*) , Afac(Ldafac,*)
!     ..
!
!  =====================================================================
!
!
!     .. Parameters ..
      REAL ZERO , ONE
      PARAMETER (ZERO=0.0E+0,ONE=1.0E+0)
!     ..
!     .. Local Scalars ..
      INTEGER i , j , k , kc , klen , ml , mu
      REAL akk , anorm , eps
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      REAL CLANHB , SLAMCH
      COMPLEX CDOTC
EXTERNAL    LSAME,   CLANHB,   SLAMCH,   CDOTC
!     ..
!     .. External Subroutines ..
EXTERNAL    CHER,   CSSCAL,   CTRMV
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC AIMAG , MAX , MIN , REAL
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
      anorm = CLANHB('1',Uplo,N,Kd,A,Lda,Rwork)
      IF ( anorm<=ZERO ) THEN
         Resid = ONE/eps
         RETURN
      ENDIF
!
!     Check the imaginary parts of the diagonal elements and return with
!     an error code if any are nonzero.
!
      IF ( LSAME(Uplo,'U') ) THEN
         DO j = 1 , N
            IF ( AIMAG(Afac(Kd+1,j))/=ZERO ) THEN
               Resid = ONE/eps
               RETURN
            ENDIF
         ENDDO
      ELSE
         DO j = 1 , N
            IF ( AIMAG(Afac(1,j))/=ZERO ) THEN
               Resid = ONE/eps
               RETURN
            ENDIF
         ENDDO
      ENDIF
!
!     Compute the product U'*U, overwriting U.
!
      IF ( LSAME(Uplo,'U') ) THEN
         DO k = N , 1 , -1
            kc = MAX(1,Kd+2-k)
            klen = Kd + 1 - kc
!
!           Compute the (K,K) element of the result.
!
            akk = CDOTC(klen+1,Afac(kc,k),1,Afac(kc,k),1)
            Afac(Kd+1,k) = akk
!
!           Compute the rest of column K.
!
            IF ( klen>0 ) CALL CTRMV('Upper','Conjugate','Non-unit',    &
     &                               klen,Afac(Kd+1,k-klen),Ldafac-1,   &
     &                               Afac(kc,k),1)
!
         ENDDO
!
!     UPLO = 'L':  Compute the product L*L', overwriting L.
!
      ELSE
         DO k = N , 1 , -1
            klen = MIN(Kd,N-k)
!
!           Add a multiple of column K of the factor L to each of
!           columns K+1 through N.
!
            IF ( klen>0 ) CALL CHER('Lower',klen,ONE,Afac(2,k),1,       &
     &                              Afac(1,k+1),Ldafac-1)
!
!           Scale column K by the diagonal element.
!
            akk = Afac(1,k)
            CALL CSSCAL(klen+1,akk,Afac(1,k),1)
!
         ENDDO
      ENDIF
!
!     Compute the difference  L*L' - A  or  U'*U - A.
!
      IF ( LSAME(Uplo,'U') ) THEN
         DO j = 1 , N
            mu = MAX(1,Kd+2-j)
            DO i = mu , Kd + 1
               Afac(i,j) = Afac(i,j) - A(i,j)
            ENDDO
         ENDDO
      ELSE
         DO j = 1 , N
            ml = MIN(Kd+1,N-j+1)
            DO i = 1 , ml
               Afac(i,j) = Afac(i,j) - A(i,j)
            ENDDO
         ENDDO
      ENDIF
!
!     Compute norm( L*L' - A ) / ( N * norm(A) * EPS )
!
      Resid = CLANHB('1',Uplo,N,Kd,Afac,Ldafac,Rwork)
!
      Resid = ((Resid/REAL(N))/anorm)/eps
!
!
!     End of cpbt01
!
      END SUBROUTINE cpbt01
