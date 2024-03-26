!*==sbdt03.f90  processed by SPAG 7.51RB at 17:37 on  4 Mar 2022
!> \brief \b sbdt03
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE sbdt03( UPLO, N, KD, D, E, U, LDU, S, VT, LDVT, WORK,
!                          RESID )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            KD, LDU, LDVT, N
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
!> sbdt03 reconstructs a bidiagonal matrix B from its SVD:
!>    S = U' * B * V
!> where U and V are orthogonal matrices and S is diagonal.
!>
!> The test ratio to test the singular value decomposition is
!>    RESID = norm( B - U * S * VT ) / ( n * norm(B) * EPS )
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
!> \param[in] KD
!> \verbatim
!>          KD is INTEGER
!>          The bandwidth of the bidiagonal matrix B.  If KD = 1, the
!>          matrix B is bidiagonal, and if KD = 0, B is diagonal and E is
!>          not referenced.  If KD is greater than 1, it is assumed to be
!>          1, and if KD is less than 0, it is assumed to be 0.
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
!> \param[in] U
!> \verbatim
!>          U is REAL array, dimension (LDU,N)
!>          The n by n orthogonal matrix U in the reduction B = U'*A*P.
!> \endverbatim
!>
!> \param[in] LDU
!> \verbatim
!>          LDU is INTEGER
!>          The leading dimension of the array U.  LDU >= max(1,N)
!> \endverbatim
!>
!> \param[in] S
!> \verbatim
!>          S is REAL array, dimension (N)
!>          The singular values from the SVD of B, sorted in decreasing
!>          order.
!> \endverbatim
!>
!> \param[in] VT
!> \verbatim
!>          VT is REAL array, dimension (LDVT,N)
!>          The n by n orthogonal matrix V' in the reduction
!>          B = U * S * V'.
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
!>          The test ratio:  norm(B - U * S * V') / ( n * norm(A) * EPS )
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
!> \ingroup single_eig
!
!  =====================================================================
      SUBROUTINE sbdt03(Uplo,N,Kd,D,E,U,Ldu,S,Vt,Ldvt,Work,Resid)
      IMPLICIT NONE
!*--SBDT03138
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Uplo
      INTEGER Kd , Ldu , Ldvt , N
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
      INTEGER i , j
      REAL bnorm , eps
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      INTEGER ISAMAX
      REAL SASUM , SLAMCH
      EXTERNAL LSAME , ISAMAX , SASUM , SLAMCH
!     ..
!     .. External Subroutines ..
      EXTERNAL SGEMV
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , MAX , MIN , REAL
!     ..
!     .. Executable Statements ..
!
!     Quick return if possible
!
      Resid = ZERO
      IF ( N<=0 ) RETURN
!
!     Compute B - U * S * V' one column at a time.
!
      bnorm = ZERO
      IF ( Kd<1 ) THEN
!
!        B is diagonal.
!
         DO j = 1 , N
            DO i = 1 , N
               Work(N+i) = S(i)*Vt(i,j)
            ENDDO
            CALL SGEMV('No transpose',N,N,-ONE,U,Ldu,Work(N+1),1,ZERO,  &
     &                 Work,1)
            Work(j) = Work(j) + D(j)
            Resid = MAX(Resid,SASUM(N,Work,1))
         ENDDO
         j = ISAMAX(N,D,1)
         bnorm = ABS(D(j))
!
!        B is bidiagonal.
!
      ELSEIF ( LSAME(Uplo,'U') ) THEN
!
!           B is upper bidiagonal.
!
         DO j = 1 , N
            DO i = 1 , N
               Work(N+i) = S(i)*Vt(i,j)
            ENDDO
            CALL SGEMV('No transpose',N,N,-ONE,U,Ldu,Work(N+1),1,ZERO,  &
     &                 Work,1)
            Work(j) = Work(j) + D(j)
            IF ( j>1 ) THEN
               Work(j-1) = Work(j-1) + E(j-1)
               bnorm = MAX(bnorm,ABS(D(j))+ABS(E(j-1)))
            ELSE
               bnorm = MAX(bnorm,ABS(D(j)))
            ENDIF
            Resid = MAX(Resid,SASUM(N,Work,1))
         ENDDO
      ELSE
!
!           B is lower bidiagonal.
!
         DO j = 1 , N
            DO i = 1 , N
               Work(N+i) = S(i)*Vt(i,j)
            ENDDO
            CALL SGEMV('No transpose',N,N,-ONE,U,Ldu,Work(N+1),1,ZERO,  &
     &                 Work,1)
            Work(j) = Work(j) + D(j)
            IF ( j<N ) THEN
               Work(j+1) = Work(j+1) + E(j)
               bnorm = MAX(bnorm,ABS(D(j))+ABS(E(j)))
            ELSE
               bnorm = MAX(bnorm,ABS(D(j)))
            ENDIF
            Resid = MAX(Resid,SASUM(N,Work,1))
         ENDDO
      ENDIF
!
!     Compute norm(B - U * S * V') / ( n * norm(B) * EPS )
!
      eps = SLAMCH('Precision')
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
!     End of sbdt03
!
      END SUBROUTINE sbdt03
