!*==cbdt02.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b CBDT02
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE CBDT02( M, N, B, LDB, C, LDC, U, LDU, WORK, RWORK,
!                          RESID )
!
!       .. Scalar Arguments ..
!       INTEGER            LDB, LDC, LDU, M, N
!       REAL               RESID
!       ..
!       .. Array Arguments ..
!       REAL               RWORK( * )
!       COMPLEX            B( LDB, * ), C( LDC, * ), U( LDU, * ),
!      $                   WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CBDT02 tests the change of basis C = U' * B by computing the residual
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
!>          B is COMPLEX array, dimension (LDB,N)
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
!>          C is COMPLEX array, dimension (LDC,N)
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
!>          U is COMPLEX array, dimension (LDU,M)
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
!>          WORK is COMPLEX array, dimension (M)
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is REAL array, dimension (M)
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
!> \ingroup complex_eig
!
!  =====================================================================
      SUBROUTINE CBDT02(M,N,B,Ldb,C,Ldc,U,Ldu,Work,Rwork,Resid)
      IMPLICIT NONE
!*--CBDT02122
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
      REAL Rwork(*)
      COMPLEX B(Ldb,*) , C(Ldc,*) , U(Ldu,*) , Work(*)
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
      REAL CLANGE , SCASUM , SLAMCH
      EXTERNAL CLANGE , SCASUM , SLAMCH
!     ..
!     .. External Subroutines ..
      EXTERNAL CCOPY , CGEMV
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC CMPLX , MAX , MIN , REAL
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
         CALL CCOPY(M,B(1,j),1,Work,1)
         CALL CGEMV('No transpose',M,M,-CMPLX(ONE),U,Ldu,C(1,j),1,      &
     &              CMPLX(ONE),Work,1)
         Resid = MAX(Resid,SCASUM(M,Work,1))
      ENDDO
!
!     Compute norm of B.
!
      bnorm = CLANGE('1',M,N,B,Ldb,Rwork)
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
!     End of CBDT02
!
      END SUBROUTINE CBDT02
