!*==spst01.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b SPST01
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE SPST01( UPLO, N, A, LDA, AFAC, LDAFAC, PERM, LDPERM,
!                          PIV, RWORK, RESID, RANK )
!
!       .. Scalar Arguments ..
!       REAL               RESID
!       INTEGER            LDA, LDAFAC, LDPERM, N, RANK
!       CHARACTER          UPLO
!       ..
!       .. Array Arguments ..
!       REAL               A( LDA, * ), AFAC( LDAFAC, * ),
!      $                   PERM( LDPERM, * ), RWORK( * )
!       INTEGER            PIV( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SPST01 reconstructs a symmetric positive semidefinite matrix A
!> from its L or U factors and the permutation matrix P and computes
!> the residual
!>    norm( P*L*L'*P' - A ) / ( N * norm(A) * EPS ) or
!>    norm( P*U'*U*P' - A ) / ( N * norm(A) * EPS ),
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
!> \param[in] AFAC
!> \verbatim
!>          AFAC is REAL array, dimension (LDAFAC,N)
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
!>          PERM is REAL array, dimension (LDPERM,N)
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
!> \ingroup single_lin
!
!  =====================================================================
      SUBROUTINE SPST01(Uplo,N,A,Lda,Afac,Ldafac,Perm,Ldperm,Piv,Rwork, &
     &                  Resid,Rank)
      IMPLICIT NONE
!*--SPST01138
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
      REAL A(Lda,*) , Afac(Ldafac,*) , Perm(Ldperm,*) , Rwork(*)
      INTEGER Piv(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ZERO , ONE
      PARAMETER (ZERO=0.0E+0,ONE=1.0E+0)
!     ..
!     .. Local Scalars ..
      REAL anorm , eps , t
      INTEGER i , j , k
!     ..
!     .. External Functions ..
      REAL SDOT , SLAMCH , SLANSY
      LOGICAL LSAME
      EXTERNAL SDOT , SLAMCH , SLANSY , LSAME
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
!
         IF ( Rank<N ) THEN
            DO j = Rank + 1 , N
               DO i = Rank + 1 , j
                  Afac(i,j) = ZERO
               ENDDO
            ENDDO
         ENDIF
!
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
!
         IF ( Rank<N ) THEN
            DO j = Rank + 1 , N
               DO i = j , N
                  Afac(i,j) = ZERO
               ENDDO
            ENDDO
         ENDIF
!
         DO k = N , 1 , -1
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
                     Perm(Piv(i),Piv(j)) = Afac(j,i)
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
                     Perm(Piv(i),Piv(j)) = Afac(j,i)
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
            DO i = 1 , j
               Perm(i,j) = Perm(i,j) - A(i,j)
            ENDDO
         ENDDO
      ELSE
         DO j = 1 , N
            DO i = j , N
               Perm(i,j) = Perm(i,j) - A(i,j)
            ENDDO
         ENDDO
      ENDIF
!
!     Compute norm( P*L*L'P - A ) / ( N * norm(A) * EPS ), or
!     ( P*U'*U*P' - A )/ ( N * norm(A) * EPS ).
!
      Resid = SLANSY('1',Uplo,N,Perm,Ldafac,Rwork)
!
      Resid = ((Resid/REAL(N))/anorm)/eps
!
!
!     End of SPST01
!
      END SUBROUTINE SPST01
