!*==sbdt02.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b SBDT02
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE SBDT02( M, N, B, LDB, C, LDC, U, LDU, WORK, RESID )
!
!       .. Scalar Arguments ..
!       INTEGER            LDB, LDC, LDU, M, N
!       REAL               RESID
!       ..
!       .. Array Arguments ..
!       REAL               B( LDB, * ), C( LDC, * ), U( LDU, * ),
!      $                   WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SBDT02 tests the change of basis C = U' * B by computing the residual
!>
!>    RESID = norm( B - U * C ) / ( max(m,n) * norm(B) * EPS ),
!>
!> where B and C are M by N matrices, U is an M by M orthogonal matrix,
!> and EPS is the machine precision.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrices B and C and the order of
!>          the matrix Q.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrices B and C.
!> \endverbatim
!>
!> \param[in] B
!> \verbatim
!>          B is REAL array, dimension (LDB,N)
!>          The m by n matrix B.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of the array B.  LDB >= max(1,M).
!> \endverbatim
!>
!> \param[in] C
!> \verbatim
!>          C is REAL array, dimension (LDC,N)
!>          The m by n matrix C, assumed to contain U' * B.
!> \endverbatim
!>
!> \param[in] LDC
!> \verbatim
!>          LDC is INTEGER
!>          The leading dimension of the array C.  LDC >= max(1,M).
!> \endverbatim
!>
!> \param[in] U
!> \verbatim
!>          U is REAL array, dimension (LDU,M)
!>          The m by m orthogonal matrix U.
!> \endverbatim
!>
!> \param[in] LDU
!> \verbatim
!>          LDU is INTEGER
!>          The leading dimension of the array U.  LDU >= max(1,M).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is REAL array, dimension (M)
!> \endverbatim
!>
!> \param[out] RESID
!> \verbatim
!>          RESID is REAL
!>          RESID = norm( B - U * C ) / ( max(m,n) * norm(B) * EPS ),
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
      SUBROUTINE SBDT02(M,N,B,Ldb,C,Ldc,U,Ldu,Work,Resid)
      IMPLICIT NONE
!*--SBDT02115
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Ldb , Ldc , Ldu , M , N
      REAL Resid
!     ..
!     .. Array Arguments ..
      REAL B(Ldb,*) , C(Ldc,*) , U(Ldu,*) , Work(*)
!     ..
!
! ======================================================================
!
!     .. Parameters ..
      REAL ZERO , ONE
      PARAMETER (ZERO=0.0E+0,ONE=1.0E+0)
!     ..
!     .. Local Scalars ..
      INTEGER j
      REAL bnorm , eps , realmn
!     ..
!     .. External Functions ..
      REAL SASUM , SLAMCH , SLANGE
      EXTERNAL SASUM , SLAMCH , SLANGE
!     ..
!     .. External Subroutines ..
      EXTERNAL SCOPY , SGEMV
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX , MIN , REAL
!     ..
!     .. Executable Statements ..
!
!     Quick return if possible
!
      Resid = ZERO
      IF ( M<=0 .OR. N<=0 ) RETURN
      realmn = REAL(MAX(M,N))
      eps = SLAMCH('Precision')
!
!     Compute norm( B - U * C )
!
      DO j = 1 , N
         CALL SCOPY(M,B(1,j),1,Work,1)
         CALL SGEMV('No transpose',M,M,-ONE,U,Ldu,C(1,j),1,ONE,Work,1)
         Resid = MAX(Resid,SASUM(M,Work,1))
      ENDDO
!
!     Compute norm of B.
!
      bnorm = SLANGE('1',M,N,B,Ldb,Work)
!
      IF ( bnorm<=ZERO ) THEN
         IF ( Resid/=ZERO ) Resid = ONE/eps
      ELSEIF ( bnorm>=Resid ) THEN
         Resid = (Resid/bnorm)/(realmn*eps)
      ELSEIF ( bnorm<ONE ) THEN
         Resid = (MIN(Resid,realmn*bnorm)/bnorm)/(realmn*eps)
      ELSE
         Resid = MIN(Resid/bnorm,realmn)/(realmn*eps)
      ENDIF
!
!     End of SBDT02
!
      END SUBROUTINE SBDT02
