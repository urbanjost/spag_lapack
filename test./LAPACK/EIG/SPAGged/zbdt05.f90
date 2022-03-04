!*==zbdt05.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b ZBDT05
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZBDT05( M, N, A, LDA, S, NS, U, LDU,
!                          VT, LDVT, WORK, RESID )
!
!       .. Scalar Arguments ..
!       INTEGER            LDA, LDU, LDVT, N, NS
!       DOUBLE PRECISION   RESID
!       ..
!       .. Array Arguments ..
!      DOUBLE PRECISION   S( * )
!      COMPLEX*16         A( LDA, * ), U( * ), VT( LDVT, * ), WORK( * )
!       ..
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZBDT05 reconstructs a bidiagonal matrix B from its (partial) SVD:
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
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrices A and U.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrices A and VT.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension (LDA,N)
!>          The m by n matrix A.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,M).
!> \endverbatim
!>
!> \param[in] S
!> \verbatim
!>          S is DOUBLE PRECISION array, dimension (NS)
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
!>          U is COMPLEX*16 array, dimension (LDU,NS)
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
!>          VT is COMPLEX*16 array, dimension (LDVT,N)
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
!>          WORK is COMPLEX*16 array, dimension (M,N)
!> \endverbatim
!>
!> \param[out] RESID
!> \verbatim
!>          RESID is DOUBLE PRECISION
!>          The test ratio:  norm(S - U' * A * V) / ( n * norm(A) * EPS )
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
      SUBROUTINE ZBDT05(M,N,A,Lda,S,Ns,U,Ldu,Vt,Ldvt,Work,Resid)
      IMPLICIT NONE
!*--ZBDT05128
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Lda , Ldu , Ldvt , M , N , Ns
      DOUBLE PRECISION Resid
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION S(*)
      COMPLEX*16 A(Lda,*) , U(*) , Vt(Ldvt,*) , Work(*)
!     ..
!
! ======================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO , ONE
      PARAMETER (ZERO=0.0D+0,ONE=1.0D+0)
      COMPLEX*16 CZERO , CONE
      PARAMETER (CZERO=(0.0D+0,0.0D+0),CONE=(1.0D+0,0.0D+0))
!     ..
!     .. Local Scalars ..
      INTEGER i , j
      DOUBLE PRECISION anorm , eps
!     ..
!     .. Local Arrays ..
      DOUBLE PRECISION dum(1)
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      INTEGER IDAMAX
      DOUBLE PRECISION DASUM , DLAMCH , ZLANGE
      EXTERNAL LSAME , IDAMAX , DASUM , DLAMCH , ZLANGE
      DOUBLE PRECISION DZASUM
!     ..
!     .. External Subroutines ..
      EXTERNAL ZGEMM
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , DBLE , MAX , MIN
!     ..
!     .. Executable Statements ..
!
!     Quick return if possible.
!
      Resid = ZERO
      IF ( MIN(M,N)<=0 .OR. Ns<=0 ) RETURN
!
      eps = DLAMCH('Precision')
      anorm = ZLANGE('M',M,N,A,Lda,dum)
!
!     Compute U' * A * V.
!
      CALL ZGEMM('N','C',M,Ns,N,CONE,A,Lda,Vt,Ldvt,CZERO,Work(1+Ns*Ns), &
     &           M)
      CALL ZGEMM('C','N',Ns,Ns,M,-CONE,U,Ldu,Work(1+Ns*Ns),M,CZERO,Work,&
     &           Ns)
!
!     norm(S - U' * B * V)
!
      j = 0
      DO i = 1 , Ns
         Work(j+i) = Work(j+i) + DCMPLX(S(i),ZERO)
         Resid = MAX(Resid,DZASUM(Ns,Work(j+1),1))
         j = j + Ns
      ENDDO
!
      IF ( anorm<=ZERO ) THEN
         IF ( Resid/=ZERO ) Resid = ONE/eps
      ELSEIF ( anorm>=Resid ) THEN
         Resid = (Resid/anorm)/(DBLE(N)*eps)
      ELSEIF ( anorm<ONE ) THEN
         Resid = (MIN(Resid,DBLE(N)*anorm)/anorm)/(DBLE(N)*eps)
      ELSE
         Resid = MIN(Resid/anorm,DBLE(N))/(DBLE(N)*eps)
      ENDIF
!
!
!     End of ZBDT05
!
      END SUBROUTINE ZBDT05
