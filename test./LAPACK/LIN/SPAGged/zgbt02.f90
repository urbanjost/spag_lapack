!*==zgbt02.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b ZGBT02
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZGBT02( TRANS, M, N, KL, KU, NRHS, A, LDA, X, LDX, B,
!                          LDB, RESID )
!
!       .. Scalar Arguments ..
!       CHARACTER          TRANS
!       INTEGER            KL, KU, LDA, LDB, LDX, M, N, NRHS
!       DOUBLE PRECISION   RESID
!       ..
!       .. Array Arguments ..
!       COMPLEX*16         A( LDA, * ), B( LDB, * ), X( LDX, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZGBT02 computes the residual for a solution of a banded system of
!> equations  A*x = b  or  A'*x = b:
!>    RESID = norm( B - A*X ) / ( norm(A) * norm(X) * EPS).
!> where EPS is the machine precision.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] TRANS
!> \verbatim
!>          TRANS is CHARACTER*1
!>          Specifies the form of the system of equations:
!>          = 'N':  A *x = b
!>          = 'T':  A'*x = b, where A' is the transpose of A
!>          = 'C':  A'*x = b, where A' is the transpose of A
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix A.  M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in] KL
!> \verbatim
!>          KL is INTEGER
!>          The number of subdiagonals within the band of A.  KL >= 0.
!> \endverbatim
!>
!> \param[in] KU
!> \verbatim
!>          KU is INTEGER
!>          The number of superdiagonals within the band of A.  KU >= 0.
!> \endverbatim
!>
!> \param[in] NRHS
!> \verbatim
!>          NRHS is INTEGER
!>          The number of columns of B.  NRHS >= 0.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension (LDA,N)
!>          The original matrix A in band storage, stored in rows 1 to
!>          KL+KU+1.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,KL+KU+1).
!> \endverbatim
!>
!> \param[in] X
!> \verbatim
!>          X is COMPLEX*16 array, dimension (LDX,NRHS)
!>          The computed solution vectors for the system of linear
!>          equations.
!> \endverbatim
!>
!> \param[in] LDX
!> \verbatim
!>          LDX is INTEGER
!>          The leading dimension of the array X.  If TRANS = 'N',
!>          LDX >= max(1,N); if TRANS = 'T' or 'C', LDX >= max(1,M).
!> \endverbatim
!>
!> \param[in,out] B
!> \verbatim
!>          B is COMPLEX*16 array, dimension (LDB,NRHS)
!>          On entry, the right hand side vectors for the system of
!>          linear equations.
!>          On exit, B is overwritten with the difference B - A*X.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of the array B.  IF TRANS = 'N',
!>          LDB >= max(1,M); if TRANS = 'T' or 'C', LDB >= max(1,N).
!> \endverbatim
!>
!> \param[out] RESID
!> \verbatim
!>          RESID is DOUBLE PRECISION
!>          The maximum over the number of right hand sides of
!>          norm(B - A*X) / ( norm(A) * norm(X) * EPS ).
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
      SUBROUTINE ZGBT02(Trans,M,N,Kl,Ku,Nrhs,A,Lda,X,Ldx,B,Ldb,Resid)
      IMPLICIT NONE
!*--ZGBT02142
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Trans
      INTEGER Kl , Ku , Lda , Ldb , Ldx , M , N , Nrhs
      DOUBLE PRECISION Resid
!     ..
!     .. Array Arguments ..
      COMPLEX*16 A(Lda,*) , B(Ldb,*) , X(Ldx,*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO , ONE
      PARAMETER (ZERO=0.0D+0,ONE=1.0D+0)
      COMPLEX*16 CONE
      PARAMETER (CONE=(1.0D+0,0.0D+0))
!     ..
!     .. Local Scalars ..
      INTEGER i1 , i2 , j , kd , n1
      DOUBLE PRECISION anorm , bnorm , eps , xnorm
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      DOUBLE PRECISION DLAMCH , DZASUM
      EXTERNAL LSAME , DLAMCH , DZASUM
!     ..
!     .. External Subroutines ..
      EXTERNAL ZGBMV
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX , MIN
!     ..
!     .. Executable Statements ..
!
!     Quick return if N = 0 pr NRHS = 0
!
      IF ( M<=0 .OR. N<=0 .OR. Nrhs<=0 ) THEN
         Resid = ZERO
         RETURN
      ENDIF
!
!     Exit with RESID = 1/EPS if ANORM = 0.
!
      eps = DLAMCH('Epsilon')
      kd = Ku + 1
      anorm = ZERO
      DO j = 1 , N
         i1 = MAX(kd+1-j,1)
         i2 = MIN(kd+M-j,Kl+kd)
         anorm = MAX(anorm,DZASUM(i2-i1+1,A(i1,j),1))
      ENDDO
      IF ( anorm<=ZERO ) THEN
         Resid = ONE/eps
         RETURN
      ENDIF
!
      IF ( LSAME(Trans,'T') .OR. LSAME(Trans,'C') ) THEN
         n1 = N
      ELSE
         n1 = M
      ENDIF
!
!     Compute  B - A*X (or  B - A'*X )
!
      DO j = 1 , Nrhs
         CALL ZGBMV(Trans,M,N,Kl,Ku,-CONE,A,Lda,X(1,j),1,CONE,B(1,j),1)
      ENDDO
!
!     Compute the maximum over the number of right hand sides of
!        norm(B - A*X) / ( norm(A) * norm(X) * EPS ).
!
      Resid = ZERO
      DO j = 1 , Nrhs
         bnorm = DZASUM(n1,B(1,j),1)
         xnorm = DZASUM(n1,X(1,j),1)
         IF ( xnorm<=ZERO ) THEN
            Resid = ONE/eps
         ELSE
            Resid = MAX(Resid,((bnorm/anorm)/xnorm)/eps)
         ENDIF
      ENDDO
!
!
!     End of ZGBT02
!
      END SUBROUTINE ZGBT02
