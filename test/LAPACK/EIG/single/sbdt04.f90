!*==sbdt04.f90  processed by SPAG 7.51RB at 17:37 on  4 Mar 2022
!> \brief \b SBDT04
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE SBDT04( UPLO, N, D, E, S, NS, U, LDU, VT, LDVT,
!                          WORK, RESID )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            LDU, LDVT, N, NS
!       REAL               RESID
!       ..
!       .. Array Arguments ..
!       REAL               D( * ), E( * ), S( * ), U( LDU, * ),
!      $                   VT( LDVT, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SBDT04 reconstructs a bidiagonal matrix B from its (partial) SVD:
!>    S = U' * B * V
!> where U and V are orthogonal matrices and S is diagonal.
!>
!> The test ratio to test the singular value decomposition is
!>    RESID = norm( S - U' * B * V ) / ( n * norm(B) * EPS )
!> where VT = V' and EPS is the machine precision.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          Specifies whether the matrix B is upper or lower bidiagonal.
!>          = 'U':  Upper bidiagonal
!>          = 'L':  Lower bidiagonal
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix B.
!> \endverbatim
!>
!> \param[in] D
!> \verbatim
!>          D is REAL array, dimension (N)
!>          The n diagonal elements of the bidiagonal matrix B.
!> \endverbatim
!>
!> \param[in] E
!> \verbatim
!>          E is REAL array, dimension (N-1)
!>          The (n-1) superdiagonal elements of the bidiagonal matrix B
!>          if UPLO = 'U', or the (n-1) subdiagonal elements of B if
!>          UPLO = 'L'.
!> \endverbatim
!>
!> \param[in] S
!> \verbatim
!>          S is REAL array, dimension (NS)
!>          The singular values from the (partial) SVD of B, sorted in
!>          decreasing order.
!> \endverbatim
!>
!> \param[in] NS
!> \verbatim
!>          NS is INTEGER
!>          The number of singular values/vectors from the (partial)
!>          SVD of B.
!> \endverbatim
!>
!> \param[in] U
!> \verbatim
!>          U is REAL array, dimension (LDU,NS)
!>          The n by ns orthogonal matrix U in S = U' * B * V.
!> \endverbatim
!>
!> \param[in] LDU
!> \verbatim
!>          LDU is INTEGER
!>          The leading dimension of the array U.  LDU >= max(1,N)
!> \endverbatim
!>
!> \param[in] VT
!> \verbatim
!>          VT is REAL array, dimension (LDVT,N)
!>          The n by ns orthogonal matrix V in S = U' * B * V.
!> \endverbatim
!>
!> \param[in] LDVT
!> \verbatim
!>          LDVT is INTEGER
!>          The leading dimension of the array VT.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is REAL array, dimension (2*N)
!> \endverbatim
!>
!> \param[out] RESID
!> \verbatim
!>          RESID is REAL
!>          The test ratio:  norm(S - U' * B * V) / ( n * norm(B) * EPS )
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
!> \ingroup double_eig
!
!  =====================================================================
      SUBROUTINE SBDT04(Uplo,N,D,E,S,Ns,U,Ldu,Vt,Ldvt,Work,Resid)
      IMPLICIT NONE
!*--SBDT04134
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Uplo
      INTEGER Ldu , Ldvt , N , Ns
      REAL Resid
!     ..
!     .. Array Arguments ..
      REAL D(*) , E(*) , S(*) , U(Ldu,*) , Vt(Ldvt,*) , Work(*)
!     ..
!
! ======================================================================
!
!     .. Parameters ..
      REAL ZERO , ONE
      PARAMETER (ZERO=0.0E+0,ONE=1.0E+0)
!     ..
!     .. Local Scalars ..
      INTEGER i , j , k
      REAL bnorm , eps
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      INTEGER ISAMAX
      REAL SASUM , SLAMCH
      EXTERNAL LSAME , ISAMAX , SASUM , SLAMCH
!     ..
!     .. External Subroutines ..
      EXTERNAL SGEMM
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , REAL , MAX , MIN
!     ..
!     .. Executable Statements ..
!
!     Quick return if possible.
!
      Resid = ZERO
      IF ( N<=0 .OR. Ns<=0 ) RETURN
!
      eps = SLAMCH('Precision')
!
!     Compute S - U' * B * V.
!
      bnorm = ZERO
!
      IF ( LSAME(Uplo,'U') ) THEN
!
!        B is upper bidiagonal.
!
         k = 0
         DO i = 1 , Ns
            DO j = 1 , N - 1
               k = k + 1
               Work(k) = D(j)*Vt(i,j) + E(j)*Vt(i,j+1)
            ENDDO
            k = k + 1
            Work(k) = D(N)*Vt(i,N)
         ENDDO
         bnorm = ABS(D(1))
         DO i = 2 , N
            bnorm = MAX(bnorm,ABS(D(i))+ABS(E(i-1)))
         ENDDO
      ELSE
!
!        B is lower bidiagonal.
!
         k = 0
         DO i = 1 , Ns
            k = k + 1
            Work(k) = D(1)*Vt(i,1)
            DO j = 1 , N - 1
               k = k + 1
               Work(k) = E(j)*Vt(i,j) + D(j+1)*Vt(i,j+1)
            ENDDO
         ENDDO
         bnorm = ABS(D(N))
         DO i = 1 , N - 1
            bnorm = MAX(bnorm,ABS(D(i))+ABS(E(i)))
         ENDDO
      ENDIF
!
      CALL SGEMM('T','N',Ns,Ns,N,-ONE,U,Ldu,Work(1),N,ZERO,Work(1+N*Ns),&
     &           Ns)
!
!     norm(S - U' * B * V)
!
      k = N*Ns
      DO i = 1 , Ns
         Work(k+i) = Work(k+i) + S(i)
         Resid = MAX(Resid,SASUM(Ns,Work(k+1),1))
         k = k + Ns
      ENDDO
!
      IF ( bnorm<=ZERO ) THEN
         IF ( Resid/=ZERO ) Resid = ONE/eps
      ELSEIF ( bnorm>=Resid ) THEN
         Resid = (Resid/bnorm)/(REAL(N)*eps)
      ELSEIF ( bnorm<ONE ) THEN
         Resid = (MIN(Resid,REAL(N)*bnorm)/bnorm)/(REAL(N)*eps)
      ELSE
         Resid = MIN(Resid/bnorm,REAL(N))/(REAL(N)*eps)
      ENDIF
!
!
!     End of SBDT04
!
      END SUBROUTINE SBDT04
