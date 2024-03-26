!*==zbdt03.f90  processed by SPAG 7.51RB at 17:37 on  4 Mar 2022
!> \brief \b zbdt03
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE zbdt03( UPLO, N, KD, D, E, U, LDU, S, VT, LDVT, WORK,
!                          RESID )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            KD, LDU, LDVT, N
!       DOUBLE PRECISION   RESID
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   D( * ), E( * ), S( * )
!       COMPLEX*16         U( LDU, * ), VT( LDVT, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> zbdt03 reconstructs a bidiagonal matrix B from its SVD:
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
!>          D is DOUBLE PRECISION array, dimension (N)
!>          The n diagonal elements of the bidiagonal matrix B.
!> \endverbatim
!>
!> \param[in] E
!> \verbatim
!>          E is DOUBLE PRECISION array, dimension (N-1)
!>          The (n-1) superdiagonal elements of the bidiagonal matrix B
!>          if UPLO = 'U', or the (n-1) subdiagonal elements of B if
!>          UPLO = 'L'.
!> \endverbatim
!>
!> \param[in] U
!> \verbatim
!>          U is COMPLEX*16 array, dimension (LDU,N)
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
!>          S is DOUBLE PRECISION array, dimension (N)
!>          The singular values from the SVD of B, sorted in decreasing
!>          order.
!> \endverbatim
!>
!> \param[in] VT
!> \verbatim
!>          VT is COMPLEX*16 array, dimension (LDVT,N)
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
!>          WORK is COMPLEX*16 array, dimension (2*N)
!> \endverbatim
!>
!> \param[out] RESID
!> \verbatim
!>          RESID is DOUBLE PRECISION
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
!> \ingroup complex16_eig
!
!  =====================================================================
      SUBROUTINE zbdt03(Uplo,N,Kd,D,E,U,Ldu,S,Vt,Ldvt,Work,Resid)
      IMPLICIT NONE
!*--ZBDT03138
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Uplo
      INTEGER Kd , Ldu , Ldvt , N
      DOUBLE PRECISION Resid
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION D(*) , E(*) , S(*)
      COMPLEX*16 U(Ldu,*) , Vt(Ldvt,*) , Work(*)
!     ..
!
! ======================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO , ONE
      PARAMETER (ZERO=0.0D+0,ONE=1.0D+0)
!     ..
!     .. Local Scalars ..
      INTEGER i , j
      DOUBLE PRECISION bnorm , eps
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      INTEGER IDAMAX
      DOUBLE PRECISION DLAMCH , DZASUM
      EXTERNAL LSAME , IDAMAX , DLAMCH , DZASUM
!     ..
!     .. External Subroutines ..
      EXTERNAL ZGEMV
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , DBLE , DCMPLX , MAX , MIN
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
            CALL ZGEMV('No transpose',N,N,-DCMPLX(ONE),U,Ldu,Work(N+1), &
     &                 1,DCMPLX(ZERO),Work,1)
            Work(j) = Work(j) + D(j)
            Resid = MAX(Resid,DZASUM(N,Work,1))
         ENDDO
         j = IDAMAX(N,D,1)
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
            CALL ZGEMV('No transpose',N,N,-DCMPLX(ONE),U,Ldu,Work(N+1), &
     &                 1,DCMPLX(ZERO),Work,1)
            Work(j) = Work(j) + D(j)
            IF ( j>1 ) THEN
               Work(j-1) = Work(j-1) + E(j-1)
               bnorm = MAX(bnorm,ABS(D(j))+ABS(E(j-1)))
            ELSE
               bnorm = MAX(bnorm,ABS(D(j)))
            ENDIF
            Resid = MAX(Resid,DZASUM(N,Work,1))
         ENDDO
      ELSE
!
!           B is lower bidiagonal.
!
         DO j = 1 , N
            DO i = 1 , N
               Work(N+i) = S(i)*Vt(i,j)
            ENDDO
            CALL ZGEMV('No transpose',N,N,-DCMPLX(ONE),U,Ldu,Work(N+1), &
     &                 1,DCMPLX(ZERO),Work,1)
            Work(j) = Work(j) + D(j)
            IF ( j<N ) THEN
               Work(j+1) = Work(j+1) + E(j)
               bnorm = MAX(bnorm,ABS(D(j))+ABS(E(j)))
            ELSE
               bnorm = MAX(bnorm,ABS(D(j)))
            ENDIF
            Resid = MAX(Resid,DZASUM(N,Work,1))
         ENDDO
      ENDIF
!
!     Compute norm(B - U * S * V') / ( n * norm(B) * EPS )
!
      eps = DLAMCH('Precision')
!
      IF ( bnorm<=ZERO ) THEN
         IF ( Resid/=ZERO ) Resid = ONE/eps
      ELSEIF ( bnorm>=Resid ) THEN
         Resid = (Resid/bnorm)/(DBLE(N)*eps)
      ELSEIF ( bnorm<ONE ) THEN
         Resid = (MIN(Resid,DBLE(N)*bnorm)/bnorm)/(DBLE(N)*eps)
      ELSE
         Resid = MIN(Resid/bnorm,DBLE(N))/(DBLE(N)*eps)
      ENDIF
!
!
!     End of zbdt03
!
      END SUBROUTINE zbdt03
