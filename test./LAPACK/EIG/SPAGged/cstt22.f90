!*==cstt22.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b CSTT22
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE CSTT22( N, M, KBAND, AD, AE, SD, SE, U, LDU, WORK,
!                          LDWORK, RWORK, RESULT )
!
!       .. Scalar Arguments ..
!       INTEGER            KBAND, LDU, LDWORK, M, N
!       ..
!       .. Array Arguments ..
!       REAL               AD( * ), AE( * ), RESULT( 2 ), RWORK( * ),
!      $                   SD( * ), SE( * )
!       COMPLEX            U( LDU, * ), WORK( LDWORK, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CSTT22  checks a set of M eigenvalues and eigenvectors,
!>
!>     A U = U S
!>
!> where A is Hermitian tridiagonal, the columns of U are unitary,
!> and S is diagonal (if KBAND=0) or Hermitian tridiagonal (if KBAND=1).
!> Two tests are performed:
!>
!>    RESULT(1) = | U* A U - S | / ( |A| m ulp )
!>
!>    RESULT(2) = | I - U*U | / ( m ulp )
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The size of the matrix.  If it is zero, CSTT22 does nothing.
!>          It must be at least zero.
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of eigenpairs to check.  If it is zero, CSTT22
!>          does nothing.  It must be at least zero.
!> \endverbatim
!>
!> \param[in] KBAND
!> \verbatim
!>          KBAND is INTEGER
!>          The bandwidth of the matrix S.  It may only be zero or one.
!>          If zero, then S is diagonal, and SE is not referenced.  If
!>          one, then S is Hermitian tri-diagonal.
!> \endverbatim
!>
!> \param[in] AD
!> \verbatim
!>          AD is REAL array, dimension (N)
!>          The diagonal of the original (unfactored) matrix A.  A is
!>          assumed to be Hermitian tridiagonal.
!> \endverbatim
!>
!> \param[in] AE
!> \verbatim
!>          AE is REAL array, dimension (N)
!>          The off-diagonal of the original (unfactored) matrix A.  A
!>          is assumed to be Hermitian tridiagonal.  AE(1) is ignored,
!>          AE(2) is the (1,2) and (2,1) element, etc.
!> \endverbatim
!>
!> \param[in] SD
!> \verbatim
!>          SD is REAL array, dimension (N)
!>          The diagonal of the (Hermitian tri-) diagonal matrix S.
!> \endverbatim
!>
!> \param[in] SE
!> \verbatim
!>          SE is REAL array, dimension (N)
!>          The off-diagonal of the (Hermitian tri-) diagonal matrix S.
!>          Not referenced if KBSND=0.  If KBAND=1, then AE(1) is
!>          ignored, SE(2) is the (1,2) and (2,1) element, etc.
!> \endverbatim
!>
!> \param[in] U
!> \verbatim
!>          U is REAL array, dimension (LDU, N)
!>          The unitary matrix in the decomposition.
!> \endverbatim
!>
!> \param[in] LDU
!> \verbatim
!>          LDU is INTEGER
!>          The leading dimension of U.  LDU must be at least N.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX array, dimension (LDWORK, M+1)
!> \endverbatim
!>
!> \param[in] LDWORK
!> \verbatim
!>          LDWORK is INTEGER
!>          The leading dimension of WORK.  LDWORK must be at least
!>          max(1,M).
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is REAL array, dimension (N)
!> \endverbatim
!>
!> \param[out] RESULT
!> \verbatim
!>          RESULT is REAL array, dimension (2)
!>          The values computed by the two tests described above.  The
!>          values are currently limited to 1/ulp, to avoid overflow.
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
!> \ingroup complex_eig
!
!  =====================================================================
      SUBROUTINE CSTT22(N,M,Kband,Ad,Ae,Sd,Se,U,Ldu,Work,Ldwork,Rwork,  &
     &                  Result)
      IMPLICIT NONE
!*--CSTT22149
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Kband , Ldu , Ldwork , M , N
!     ..
!     .. Array Arguments ..
      REAL Ad(*) , Ae(*) , Result(2) , Rwork(*) , Sd(*) , Se(*)
      COMPLEX U(Ldu,*) , Work(Ldwork,*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ZERO , ONE
      PARAMETER (ZERO=0.0E0,ONE=1.0E0)
      COMPLEX CZERO , CONE
      PARAMETER (CZERO=(0.0E+0,0.0E+0),CONE=(1.0E+0,0.0E+0))
!     ..
!     .. Local Scalars ..
      INTEGER i , j , k
      REAL anorm , ulp , unfl , wnorm
      COMPLEX aukj
!     ..
!     .. External Functions ..
      REAL CLANGE , CLANSY , SLAMCH
      EXTERNAL CLANGE , CLANSY , SLAMCH
!     ..
!     .. External Subroutines ..
      EXTERNAL CGEMM
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , MAX , MIN , REAL
!     ..
!     .. Executable Statements ..
!
      Result(1) = ZERO
      Result(2) = ZERO
      IF ( N<=0 .OR. M<=0 ) RETURN
!
      unfl = SLAMCH('Safe minimum')
      ulp = SLAMCH('Epsilon')
!
!     Do Test 1
!
!     Compute the 1-norm of A.
!
      IF ( N>1 ) THEN
         anorm = ABS(Ad(1)) + ABS(Ae(1))
         DO j = 2 , N - 1
            anorm = MAX(anorm,ABS(Ad(j))+ABS(Ae(j))+ABS(Ae(j-1)))
         ENDDO
         anorm = MAX(anorm,ABS(Ad(N))+ABS(Ae(N-1)))
      ELSE
         anorm = ABS(Ad(1))
      ENDIF
      anorm = MAX(anorm,unfl)
!
!     Norm of U*AU - S
!
      DO i = 1 , M
         DO j = 1 , M
            Work(i,j) = CZERO
            DO k = 1 , N
               aukj = Ad(k)*U(k,j)
               IF ( k/=N ) aukj = aukj + Ae(k)*U(k+1,j)
               IF ( k/=1 ) aukj = aukj + Ae(k-1)*U(k-1,j)
               Work(i,j) = Work(i,j) + U(k,i)*aukj
            ENDDO
         ENDDO
         Work(i,i) = Work(i,i) - Sd(i)
         IF ( Kband==1 ) THEN
            IF ( i/=1 ) Work(i,i-1) = Work(i,i-1) - Se(i-1)
            IF ( i/=N ) Work(i,i+1) = Work(i,i+1) - Se(i)
         ENDIF
      ENDDO
!
      wnorm = CLANSY('1','L',M,Work,M,Rwork)
!
      IF ( anorm>wnorm ) THEN
         Result(1) = (wnorm/anorm)/(M*ulp)
      ELSEIF ( anorm<ONE ) THEN
         Result(1) = (MIN(wnorm,M*anorm)/anorm)/(M*ulp)
      ELSE
         Result(1) = MIN(wnorm/anorm,REAL(M))/(M*ulp)
      ENDIF
!
!     Do Test 2
!
!     Compute  U*U - I
!
      CALL CGEMM('T','N',M,M,N,CONE,U,Ldu,U,Ldu,CZERO,Work,M)
!
      DO j = 1 , M
         Work(j,j) = Work(j,j) - ONE
      ENDDO
!
      Result(2) = MIN(REAL(M),CLANGE('1',M,M,Work,M,Rwork))/(M*ulp)
!
!
!     End of CSTT22
!
      END SUBROUTINE CSTT22
