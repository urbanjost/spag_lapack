!*==cpst01.f90  processed by SPAG 7.51RB at 17:34 on  4 Mar 2022
!> \brief \b cpst01
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE cpst01( UPLO, N, A, LDA, AFAC, LDAFAC, PERM, LDPERM,
!                          PIV, RWORK, RESID, RANK )
!
!       .. Scalar Arguments ..
!       REAL               RESID
!       INTEGER            LDA, LDAFAC, LDPERM, N, RANK
!       CHARACTER          UPLO
!       ..
!       .. Array Arguments ..
!       COMPLEX            A( LDA, * ), AFAC( LDAFAC, * ),
!      $                   PERM( LDPERM, * )
!       REAL               RWORK( * )
!       INTEGER            PIV( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> cpst01 reconstructs an Hermitian positive semidefinite matrix A
!> from its L or U factors and the permutation matrix P and computes
!> the residual
!>    norm( P*L*L'*P' - A ) / ( N * norm(A) * EPS ) or
!>    norm( P*U'*U*P' - A ) / ( N * norm(A) * EPS ),
!> where EPS is the machine epsilon, L' is the conjugate transpose of L,
!> and U' is the conjugate transpose of U.
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
!>          A is COMPLEX array, dimension (LDA,N)
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
!>          AFAC is COMPLEX array, dimension (LDAFAC,N)
!>          The factor L or U from the L*L' or U'*U
!>          factorization of A.
!> \endverbatim
!>
!> \param[in] LDAFAC
!> \verbatim
!>          LDAFAC is INTEGER
!>          The leading dimension of the array AFAC.  LDAFAC >= max(1,N).
!> \endverbatim
!>
!> \param[out] PERM
!> \verbatim
!>          PERM is COMPLEX array, dimension (LDPERM,N)
!>          Overwritten with the reconstructed matrix, and then with the
!>          difference P*L*L'*P' - A (or P*U'*U*P' - A)
!> \endverbatim
!>
!> \param[in] LDPERM
!> \verbatim
!>          LDPERM is INTEGER
!>          The leading dimension of the array PERM.
!>          LDAPERM >= max(1,N).
!> \endverbatim
!>
!> \param[in] PIV
!> \verbatim
!>          PIV is INTEGER array, dimension (N)
!>          PIV is such that the nonzero entries are
!>          P( PIV( K ), K ) = 1.
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
!>
!> \param[in] RANK
!> \verbatim
!>          RANK is INTEGER
!>          number of nonzero singular values of A.
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
      SUBROUTINE cpst01(Uplo,N,A,Lda,Afac,Ldafac,Perm,Ldperm,Piv,Rwork, &
     &                  Resid,Rank)
      IMPLICIT NONE
!*--CPST01140
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      REAL Resid
      INTEGER Lda , Ldafac , Ldperm , N , Rank
      CHARACTER Uplo
!     ..
!     .. Array Arguments ..
      COMPLEX A(Lda,*) , Afac(Ldafac,*) , Perm(Ldperm,*)
      REAL Rwork(*)
      INTEGER Piv(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ZERO , ONE
      PARAMETER (ZERO=0.0E+0,ONE=1.0E+0)
      COMPLEX CZERO
      PARAMETER (CZERO=(0.0E+0,0.0E+0))
!     ..
!     .. Local Scalars ..
      COMPLEX tc
      REAL anorm , eps , tr
      INTEGER i , j , k
!     ..
!     .. External Functions ..
      COMPLEX CDOTC
      REAL CLANHE , SLAMCH
      LOGICAL LSAME
EXTERNAL    CDOTC,   CLANHE,   SLAMCH,   LSAME
!     ..
!     .. External Subroutines ..
EXTERNAL    CHER,   CSCAL,   CTRMV
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC AIMAG , CONJG , REAL
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
      anorm = CLANHE('1',Uplo,N,A,Lda,Rwork)
      IF ( anorm<=ZERO ) THEN
         Resid = ONE/eps
         RETURN
      ENDIF
!
!     Check the imaginary parts of the diagonal elements and return with
!     an error code if any are nonzero.
!
      DO j = 1 , N
         IF ( AIMAG(Afac(j,j))/=ZERO ) THEN
            Resid = ONE/eps
            RETURN
         ENDIF
      ENDDO
!
!     Compute the product U'*U, overwriting U.
!
      IF ( LSAME(Uplo,'U') ) THEN
!
         IF ( Rank<N ) THEN
            DO j = Rank + 1 , N
               DO i = Rank + 1 , j
                  Afac(i,j) = CZERO
               ENDDO
            ENDDO
         ENDIF
!
         DO k = N , 1 , -1
!
!           Compute the (K,K) element of the result.
!
            tr = CDOTC(k,Afac(1,k),1,Afac(1,k),1)
            Afac(k,k) = tr
!
!           Compute the rest of column K.
!
            CALL CTRMV('Upper','Conjugate','Non-unit',k-1,Afac,Ldafac,  &
     &                 Afac(1,k),1)
!
         ENDDO
!
!     Compute the product L*L', overwriting L.
!
      ELSE
!
         IF ( Rank<N ) THEN
            DO j = Rank + 1 , N
               DO i = j , N
                  Afac(i,j) = CZERO
               ENDDO
            ENDDO
         ENDIF
!
         DO k = N , 1 , -1
!           Add a multiple of column K of the factor L to each of
!           columns K+1 through N.
!
            IF ( k+1<=N ) CALL CHER('Lower',N-k,ONE,Afac(k+1,k),1,      &
     &                              Afac(k+1,k+1),Ldafac)
!
!           Scale column K by the diagonal element.
!
            tc = Afac(k,k)
            CALL CSCAL(N-k+1,tc,Afac(k,k),1)
         ENDDO
!
      ENDIF
!
!        Form P*L*L'*P' or P*U'*U*P'
!
      IF ( LSAME(Uplo,'U') ) THEN
!
         DO j = 1 , N
            DO i = 1 , N
               IF ( Piv(i)<=Piv(j) ) THEN
                  IF ( i<=j ) THEN
                     Perm(Piv(i),Piv(j)) = Afac(i,j)
                  ELSE
                     Perm(Piv(i),Piv(j)) = CONJG(Afac(j,i))
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
!
!
      ELSE
!
         DO j = 1 , N
            DO i = 1 , N
               IF ( Piv(i)>=Piv(j) ) THEN
                  IF ( i>=j ) THEN
                     Perm(Piv(i),Piv(j)) = Afac(i,j)
                  ELSE
                     Perm(Piv(i),Piv(j)) = CONJG(Afac(j,i))
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
!
      ENDIF
!
!     Compute the difference  P*L*L'*P' - A (or P*U'*U*P' - A).
!
      IF ( LSAME(Uplo,'U') ) THEN
         DO j = 1 , N
            DO i = 1 , j - 1
               Perm(i,j) = Perm(i,j) - A(i,j)
            ENDDO
            Perm(j,j) = Perm(j,j) - REAL(A(j,j))
         ENDDO
      ELSE
         DO j = 1 , N
            Perm(j,j) = Perm(j,j) - REAL(A(j,j))
            DO i = j + 1 , N
               Perm(i,j) = Perm(i,j) - A(i,j)
            ENDDO
         ENDDO
      ENDIF
!
!     Compute norm( P*L*L'P - A ) / ( N * norm(A) * EPS ), or
!     ( P*U'*U*P' - A )/ ( N * norm(A) * EPS ).
!
      Resid = CLANHE('1',Uplo,N,Perm,Ldafac,Rwork)
!
      Resid = ((Resid/REAL(N))/anorm)/eps
!
!
!     End of cpst01
!
      END SUBROUTINE cpst01
