!*==zbdt02.f90  processed by SPAG 7.51RB at 17:37 on  4 Mar 2022
!> \brief \b ZBDT02
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZBDT02( M, N, B, LDB, C, LDC, U, LDU, WORK, RWORK,
!                          RESID )
!
!       .. Scalar Arguments ..
!       INTEGER            LDB, LDC, LDU, M, N
!       DOUBLE PRECISION   RESID
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   RWORK( * )
!       COMPLEX*16         B( LDB, * ), C( LDC, * ), U( LDU, * ),
!      $                   WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZBDT02 tests the change of basis C = U' * B by computing the residual
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
!>          B is COMPLEX*16 array, dimension (LDB,N)
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
!>          C is COMPLEX*16 array, dimension (LDC,N)
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
!>          U is COMPLEX*16 array, dimension (LDU,M)
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
!>          WORK is COMPLEX*16 array, dimension (M)
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is DOUBLE PRECISION array, dimension (M)
!> \endverbatim
!>
!> \param[out] RESID
!> \verbatim
!>          RESID is DOUBLE PRECISION
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
!> \ingroup complex16_eig
!
!  =====================================================================
      SUBROUTINE ZBDT02(M,N,B,Ldb,C,Ldc,U,Ldu,Work,Rwork,Resid)
      IMPLICIT NONE
!*--ZBDT02122
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Ldb , Ldc , Ldu , M , N
      DOUBLE PRECISION Resid
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION Rwork(*)
      COMPLEX*16 B(Ldb,*) , C(Ldc,*) , U(Ldu,*) , Work(*)
!     ..
!
! ======================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO , ONE
      PARAMETER (ZERO=0.0D+0,ONE=1.0D+0)
!     ..
!     .. Local Scalars ..
      INTEGER j
      DOUBLE PRECISION bnorm , eps , realmn
!     ..
!     .. External Functions ..
      DOUBLE PRECISION DLAMCH , DZASUM , ZLANGE
      EXTERNAL DLAMCH , DZASUM , ZLANGE
!     ..
!     .. External Subroutines ..
      EXTERNAL ZCOPY , ZGEMV
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC DBLE , DCMPLX , MAX , MIN
!     ..
!     .. Executable Statements ..
!
!     Quick return if possible
!
      Resid = ZERO
      IF ( M<=0 .OR. N<=0 ) RETURN
      realmn = DBLE(MAX(M,N))
      eps = DLAMCH('Precision')
!
!     Compute norm( B - U * C )
!
      DO j = 1 , N
         CALL ZCOPY(M,B(1,j),1,Work,1)
         CALL ZGEMV('No transpose',M,M,-DCMPLX(ONE),U,Ldu,C(1,j),1,     &
     &              DCMPLX(ONE),Work,1)
         Resid = MAX(Resid,DZASUM(M,Work,1))
      ENDDO
!
!     Compute norm of B.
!
      bnorm = ZLANGE('1',M,N,B,Ldb,Rwork)
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
!     End of ZBDT02
!
      END SUBROUTINE ZBDT02
