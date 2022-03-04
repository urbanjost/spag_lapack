!*==zpst01.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b ZPST01
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZPST01( UPLO, N, A, LDA, AFAC, LDAFAC, PERM, LDPERM,
!                          PIV, RWORK, RESID, RANK )
!
!       .. Scalar Arguments ..
!       DOUBLE PRECISION   RESID
!       INTEGER            LDA, LDAFAC, LDPERM, N, RANK
!       CHARACTER          UPLO
!       ..
!       .. Array Arguments ..
!       COMPLEX*16         A( LDA, * ), AFAC( LDAFAC, * ),
!      $                   PERM( LDPERM, * )
!       DOUBLE PRECISION   RWORK( * )
!       INTEGER            PIV( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZPST01 reconstructs an Hermitian positive semidefinite matrix A
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
!>          PERM is COMPLEX*16 array, dimension (LDPERM,N)
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
!>          RWORK is DOUBLE PRECISION array, dimension (N)
!> \endverbatim
!>
!> \param[out] RESID
!> \verbatim
!>          RESID is DOUBLE PRECISION
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
!> \ingroup complex16_lin
!
!  =====================================================================
      SUBROUTINE ZPST01(Uplo,N,A,Lda,Afac,Ldafac,Perm,Ldperm,Piv,Rwork, &
     &                  Resid,Rank)
      IMPLICIT NONE
!*--ZPST01140
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      DOUBLE PRECISION Resid
      INTEGER Lda , Ldafac , Ldperm , N , Rank
      CHARACTER Uplo
!     ..
!     .. Array Arguments ..
      COMPLEX*16 A(Lda,*) , Afac(Ldafac,*) , Perm(Ldperm,*)
      DOUBLE PRECISION Rwork(*)
      INTEGER Piv(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO , ONE
      PARAMETER (ZERO=0.0D+0,ONE=1.0D+0)
      COMPLEX*16 CZERO
      PARAMETER (CZERO=(0.0D+0,0.0D+0))
!     ..
!     .. Local Scalars ..
      COMPLEX*16 tc
      DOUBLE PRECISION anorm , eps , tr
      INTEGER i , j , k
!     ..
!     .. External Functions ..
      COMPLEX*16 ZDOTC
      DOUBLE PRECISION DLAMCH , ZLANHE
      LOGICAL LSAME
      EXTERNAL ZDOTC , DLAMCH , ZLANHE , LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL ZHER , ZSCAL , ZTRMV
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC DBLE , DCONJG , DIMAG
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
      eps = DLAMCH('Epsilon')
      anorm = ZLANHE('1',Uplo,N,A,Lda,Rwork)
      IF ( anorm<=ZERO ) THEN
         Resid = ONE/eps
         RETURN
      ENDIF
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
            tr = ZDOTC(k,Afac(1,k),1,Afac(1,k),1)
            Afac(k,k) = tr
!
!           Compute the rest of column K.
!
            CALL ZTRMV('Upper','Conjugate','Non-unit',k-1,Afac,Ldafac,  &
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
            IF ( k+1<=N ) CALL ZHER('Lower',N-k,ONE,Afac(k+1,k),1,      &
     &                              Afac(k+1,k+1),Ldafac)
!
!           Scale column K by the diagonal element.
!
            tc = Afac(k,k)
            CALL ZSCAL(N-k+1,tc,Afac(k,k),1)
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
                     Perm(Piv(i),Piv(j)) = DCONJG(Afac(j,i))
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
                     Perm(Piv(i),Piv(j)) = DCONJG(Afac(j,i))
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
            Perm(j,j) = Perm(j,j) - DBLE(A(j,j))
         ENDDO
      ELSE
         DO j = 1 , N
            Perm(j,j) = Perm(j,j) - DBLE(A(j,j))
            DO i = j + 1 , N
               Perm(i,j) = Perm(i,j) - A(i,j)
            ENDDO
         ENDDO
      ENDIF
!
!     Compute norm( P*L*L'P - A ) / ( N * norm(A) * EPS ), or
!     ( P*U'*U*P' - A )/ ( N * norm(A) * EPS ).
!
      Resid = ZLANHE('1',Uplo,N,Perm,Ldafac,Rwork)
!
      Resid = ((Resid/DBLE(N))/anorm)/eps
!
!
!     End of ZPST01
!
      END SUBROUTINE ZPST01
