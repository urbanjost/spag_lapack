!*==sbdt05.f90  processed by SPAG 7.51RB at 17:37 on  4 Mar 2022
!> \brief \b SBDT05
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE SBDT05( M, N, A, LDA, S, NS, U, LDU,
!                          VT, LDVT, WORK, RESID )
!
!       .. Scalar Arguments ..
!       INTEGER            LDA, LDU, LDVT, N, NS
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
!> SBDT05 reconstructs a bidiagonal matrix B from its (partial) SVD:
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
!>          A is REAL array, dimension (LDA,N)
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
!>          WORK is REAL array, dimension (M,N)
!> \endverbatim
!>
!> \param[out] RESID
!> \verbatim
!>          RESID is REAL
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
      SUBROUTINE SBDT05(M,N,A,Lda,S,Ns,U,Ldu,Vt,Ldvt,Work,Resid)
      IMPLICIT NONE
!*--SBDT05130
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Lda , Ldu , Ldvt , M , N , Ns
      REAL Resid
!     ..
!     .. Array Arguments ..
      REAL A(Lda,*) , S(*) , U(Ldu,*) , Vt(Ldvt,*) , Work(*)
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
      REAL anorm , eps
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      INTEGER ISAMAX
      REAL SASUM , SLAMCH , SLANGE
      EXTERNAL LSAME , ISAMAX , SASUM , SLAMCH , SLANGE
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
      IF ( MIN(M,N)<=0 .OR. Ns<=0 ) RETURN
!
      eps = SLAMCH('Precision')
      anorm = SLANGE('M',M,N,A,Lda,Work)
!
!     Compute U' * A * V.
!
      CALL SGEMM('N','T',M,Ns,N,ONE,A,Lda,Vt,Ldvt,ZERO,Work(1+Ns*Ns),M)
      CALL SGEMM('T','N',Ns,Ns,M,-ONE,U,Ldu,Work(1+Ns*Ns),M,ZERO,Work,  &
     &           Ns)
!
!     norm(S - U' * B * V)
!
      j = 0
      DO i = 1 , Ns
         Work(j+i) = Work(j+i) + S(i)
         Resid = MAX(Resid,SASUM(Ns,Work(j+1),1))
         j = j + Ns
      ENDDO
!
      IF ( anorm<=ZERO ) THEN
         IF ( Resid/=ZERO ) Resid = ONE/eps
      ELSEIF ( anorm>=Resid ) THEN
         Resid = (Resid/anorm)/(REAL(N)*eps)
      ELSEIF ( anorm<ONE ) THEN
         Resid = (MIN(Resid,REAL(N)*anorm)/anorm)/(REAL(N)*eps)
      ELSE
         Resid = MIN(Resid/anorm,REAL(N))/(REAL(N)*eps)
      ENDIF
!
!
!     End of SBDT05
!
      END SUBROUTINE SBDT05
