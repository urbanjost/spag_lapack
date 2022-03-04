!*==dqrt16.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b DQRT16
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DQRT16( TRANS, M, N, NRHS, A, LDA, X, LDX, B, LDB,
!                          RWORK, RESID )
!
!       .. Scalar Arguments ..
!       CHARACTER          TRANS
!       INTEGER            LDA, LDB, LDX, M, N, NRHS
!       DOUBLE PRECISION   RESID
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), RWORK( * ),
!      $                   X( LDX, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DQRT16 computes the residual for a solution of a system of linear
!> equations  A*x = b  or  A'*x = b:
!>    RESID = norm(B - A*X) / ( max(m,n) * norm(A) * norm(X) * EPS ),
!> where EPS is the machine epsilon.
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
!> \param[in] NRHS
!> \verbatim
!>          NRHS is INTEGER
!>          The number of columns of B, the matrix of right hand sides.
!>          NRHS >= 0.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (LDA,N)
!>          The original M x N matrix A.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,M).
!> \endverbatim
!>
!> \param[in] X
!> \verbatim
!>          X is DOUBLE PRECISION array, dimension (LDX,NRHS)
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
!>          B is DOUBLE PRECISION array, dimension (LDB,NRHS)
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
!> \param[out] RWORK
!> \verbatim
!>          RWORK is DOUBLE PRECISION array, dimension (M)
!> \endverbatim
!>
!> \param[out] RESID
!> \verbatim
!>          RESID is DOUBLE PRECISION
!>          The maximum over the number of right hand sides of
!>          norm(B - A*X) / ( max(m,n) * norm(A) * norm(X) * EPS ).
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
!> \ingroup double_lin
!
!  =====================================================================
      SUBROUTINE DQRT16(Trans,M,N,Nrhs,A,Lda,X,Ldx,B,Ldb,Rwork,Resid)
      IMPLICIT NONE
!*--DQRT16136
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Trans
      INTEGER Lda , Ldb , Ldx , M , N , Nrhs
      DOUBLE PRECISION Resid
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION A(Lda,*) , B(Ldb,*) , Rwork(*) , X(Ldx,*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO , ONE
      PARAMETER (ZERO=0.0D+0,ONE=1.0D+0)
!     ..
!     .. Local Scalars ..
      INTEGER j , n1 , n2
      DOUBLE PRECISION anorm , bnorm , eps , xnorm
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      DOUBLE PRECISION DASUM , DLAMCH , DLANGE
      EXTERNAL LSAME , DASUM , DLAMCH , DLANGE
!     ..
!     .. External Subroutines ..
      EXTERNAL DGEMM
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX
!     ..
!     .. Executable Statements ..
!
!     Quick exit if M = 0 or N = 0 or NRHS = 0
!
      IF ( M<=0 .OR. N<=0 .OR. Nrhs==0 ) THEN
         Resid = ZERO
         RETURN
      ENDIF
!
      IF ( LSAME(Trans,'T') .OR. LSAME(Trans,'C') ) THEN
         anorm = DLANGE('I',M,N,A,Lda,Rwork)
         n1 = N
         n2 = M
      ELSE
         anorm = DLANGE('1',M,N,A,Lda,Rwork)
         n1 = M
         n2 = N
      ENDIF
!
      eps = DLAMCH('Epsilon')
!
!     Compute  B - A*X  (or  B - A'*X ) and store in B.
!
      CALL DGEMM(Trans,'No transpose',n1,Nrhs,n2,-ONE,A,Lda,X,Ldx,ONE,B,&
     &           Ldb)
!
!     Compute the maximum over the number of right hand sides of
!        norm(B - A*X) / ( max(m,n) * norm(A) * norm(X) * EPS ) .
!
      Resid = ZERO
      DO j = 1 , Nrhs
         bnorm = DASUM(n1,B(1,j),1)
         xnorm = DASUM(n2,X(1,j),1)
         IF ( anorm==ZERO .AND. bnorm==ZERO ) THEN
            Resid = ZERO
         ELSEIF ( anorm<=ZERO .OR. xnorm<=ZERO ) THEN
            Resid = ONE/eps
         ELSE
            Resid = MAX(Resid,((bnorm/anorm)/xnorm)/(MAX(M,N)*eps))
         ENDIF
      ENDDO
!
!
!     End of DQRT16
!
      END SUBROUTINE DQRT16
