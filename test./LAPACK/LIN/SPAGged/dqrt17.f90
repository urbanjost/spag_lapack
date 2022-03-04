!*==dqrt17.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b DQRT17
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       DOUBLE PRECISION FUNCTION DQRT17( TRANS, IRESID, M, N, NRHS, A,
!                        LDA, X, LDX, B, LDB, C, WORK, LWORK )
!
!       .. Scalar Arguments ..
!       CHARACTER          TRANS
!       INTEGER            IRESID, LDA, LDB, LDX, LWORK, M, N, NRHS
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), C( LDB, * ),
!      $                   WORK( LWORK ), X( LDX, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DQRT17 computes the ratio
!>
!>    || R'*op(A) ||/(||A||*alpha*max(M,N,NRHS)*eps)
!>
!> where R = op(A)*X - B, op(A) is A or A', and
!>
!>    alpha = ||B|| if IRESID = 1 (zero-residual problem)
!>    alpha = ||R|| if IRESID = 2 (otherwise).
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] TRANS
!> \verbatim
!>          TRANS is CHARACTER*1
!>          Specifies whether or not the transpose of A is used.
!>          = 'N':  No transpose, op(A) = A.
!>          = 'T':  Transpose, op(A) = A'.
!> \endverbatim
!>
!> \param[in] IRESID
!> \verbatim
!>          IRESID is INTEGER
!>          IRESID = 1 indicates zero-residual problem.
!>          IRESID = 2 indicates non-zero residual.
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix A.
!>          If TRANS = 'N', the number of rows of the matrix B.
!>          If TRANS = 'T', the number of rows of the matrix X.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix  A.
!>          If TRANS = 'N', the number of rows of the matrix X.
!>          If TRANS = 'T', the number of rows of the matrix B.
!> \endverbatim
!>
!> \param[in] NRHS
!> \verbatim
!>          NRHS is INTEGER
!>          The number of columns of the matrices X and B.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (LDA,N)
!>          The m-by-n matrix A.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A. LDA >= M.
!> \endverbatim
!>
!> \param[in] X
!> \verbatim
!>          X is DOUBLE PRECISION array, dimension (LDX,NRHS)
!>          If TRANS = 'N', the n-by-nrhs matrix X.
!>          If TRANS = 'T', the m-by-nrhs matrix X.
!> \endverbatim
!>
!> \param[in] LDX
!> \verbatim
!>          LDX is INTEGER
!>          The leading dimension of the array X.
!>          If TRANS = 'N', LDX >= N.
!>          If TRANS = 'T', LDX >= M.
!> \endverbatim
!>
!> \param[in] B
!> \verbatim
!>          B is DOUBLE PRECISION array, dimension (LDB,NRHS)
!>          If TRANS = 'N', the m-by-nrhs matrix B.
!>          If TRANS = 'T', the n-by-nrhs matrix B.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of the array B.
!>          If TRANS = 'N', LDB >= M.
!>          If TRANS = 'T', LDB >= N.
!> \endverbatim
!>
!> \param[out] C
!> \verbatim
!>          C is DOUBLE PRECISION array, dimension (LDB,NRHS)
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (LWORK)
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The length of the array WORK.  LWORK >= NRHS*(M+N).
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
      DOUBLE PRECISION FUNCTION DQRT17(Trans,Iresid,M,N,Nrhs,A,Lda,X,   &
     &                                 Ldx,B,Ldb,C,Work,Lwork)
      IMPLICIT NONE
!*--DQRT17154
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Trans
      INTEGER Iresid , Lda , Ldb , Ldx , Lwork , M , N , Nrhs
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION A(Lda,*) , B(Ldb,*) , C(Ldb,*) , Work(Lwork) ,   &
     &                 X(Ldx,*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO , ONE
      PARAMETER (ZERO=0.0D0,ONE=1.0D0)
!     ..
!     .. Local Scalars ..
      INTEGER info , iscl , ncols , nrows
      DOUBLE PRECISION bignum , err , norma , normb , normrs , smlnum
!     ..
!     .. Local Arrays ..
      DOUBLE PRECISION rwork(1)
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      DOUBLE PRECISION DLAMCH , DLANGE
      EXTERNAL LSAME , DLAMCH , DLANGE
!     ..
!     .. External Subroutines ..
      EXTERNAL DGEMM , DLACPY , DLASCL , XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC DBLE , MAX
!     ..
!     .. Executable Statements ..
!
      DQRT17 = ZERO
!
      IF ( LSAME(Trans,'N') ) THEN
         nrows = M
         ncols = N
      ELSEIF ( LSAME(Trans,'T') ) THEN
         nrows = N
         ncols = M
      ELSE
         CALL XERBLA('DQRT17',1)
         RETURN
      ENDIF
!
      IF ( Lwork<ncols*Nrhs ) THEN
         CALL XERBLA('DQRT17',13)
         RETURN
      ENDIF
!
      IF ( M<=0 .OR. N<=0 .OR. Nrhs<=0 ) RETURN
!
      norma = DLANGE('One-norm',M,N,A,Lda,rwork)
      smlnum = DLAMCH('Safe minimum')/DLAMCH('Precision')
      bignum = ONE/smlnum
      iscl = 0
!
!     compute residual and scale it
!
      CALL DLACPY('All',nrows,Nrhs,B,Ldb,C,Ldb)
      CALL DGEMM(Trans,'No transpose',nrows,Nrhs,ncols,-ONE,A,Lda,X,Ldx,&
     &           ONE,C,Ldb)
      normrs = DLANGE('Max',nrows,Nrhs,C,Ldb,rwork)
      IF ( normrs>smlnum ) THEN
         iscl = 1
         CALL DLASCL('General',0,0,normrs,ONE,nrows,Nrhs,C,Ldb,info)
      ENDIF
!
!     compute R'*A
!
      CALL DGEMM('Transpose',Trans,Nrhs,ncols,nrows,ONE,C,Ldb,A,Lda,    &
     &           ZERO,Work,Nrhs)
!
!     compute and properly scale error
!
      err = DLANGE('One-norm',Nrhs,ncols,Work,Nrhs,rwork)
      IF ( norma/=ZERO ) err = err/norma
!
      IF ( iscl==1 ) err = err*normrs
!
      IF ( Iresid==1 ) THEN
         normb = DLANGE('One-norm',nrows,Nrhs,B,Ldb,rwork)
         IF ( normb/=ZERO ) err = err/normb
      ELSE
         IF ( normrs/=ZERO ) err = err/normrs
      ENDIF
!
      DQRT17 = err/(DLAMCH('Epsilon')*DBLE(MAX(M,N,Nrhs)))
!
!     End of DQRT17
!
      END FUNCTION DQRT17
