!*==zqrt17.f90  processed by SPAG 7.51RB at 17:35 on  4 Mar 2022
!> \brief \b zqrt17
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       DOUBLE PRECISION FUNCTION zqrt17( TRANS, IRESID, M, N, NRHS, A,
!                        LDA, X, LDX, B, LDB, C, WORK, LWORK )
!
!       .. Scalar Arguments ..
!       CHARACTER          TRANS
!       INTEGER            IRESID, LDA, LDB, LDX, LWORK, M, N, NRHS
!       ..
!       .. Array Arguments ..
!       COMPLEX*16         A( LDA, * ), B( LDB, * ), C( LDB, * ),
!      $                   WORK( LWORK ), X( LDX, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> zqrt17 computes the ratio
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
!>          = 'C':  Conjugate transpose, op(A) = A'.
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
!>          If TRANS = 'C', the number of rows of the matrix X.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix  A.
!>          If TRANS = 'N', the number of rows of the matrix X.
!>          If TRANS = 'C', the number of rows of the matrix B.
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
!>          A is COMPLEX*16 array, dimension (LDA,N)
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
!>          X is COMPLEX*16 array, dimension (LDX,NRHS)
!>          If TRANS = 'N', the n-by-nrhs matrix X.
!>          If TRANS = 'C', the m-by-nrhs matrix X.
!> \endverbatim
!>
!> \param[in] LDX
!> \verbatim
!>          LDX is INTEGER
!>          The leading dimension of the array X.
!>          If TRANS = 'N', LDX >= N.
!>          If TRANS = 'C', LDX >= M.
!> \endverbatim
!>
!> \param[in] B
!> \verbatim
!>          B is COMPLEX*16 array, dimension (LDB,NRHS)
!>          If TRANS = 'N', the m-by-nrhs matrix B.
!>          If TRANS = 'C', the n-by-nrhs matrix B.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of the array B.
!>          If TRANS = 'N', LDB >= M.
!>          If TRANS = 'C', LDB >= N.
!> \endverbatim
!>
!> \param[out] C
!> \verbatim
!>          C is COMPLEX*16 array, dimension (LDB,NRHS)
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX*16 array, dimension (LWORK)
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
!> \ingroup complex16_lin
!
!  =====================================================================
      DOUBLE PRECISION FUNCTION zqrt17(Trans,Iresid,M,N,Nrhs,A,Lda,X,   &
     &                                 Ldx,B,Ldb,C,Work,Lwork)
      IMPLICIT NONE
!*--ZQRT17154
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
      COMPLEX*16 A(Lda,*) , B(Ldb,*) , C(Ldb,*) , Work(Lwork) , X(Ldx,*)
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
      DOUBLE PRECISION DLAMCH , ZLANGE
EXTERNAL    LSAME,   DLAMCH,   ZLANGE
!     ..
!     .. External Subroutines ..
EXTERNAL    ZGEMM,   ZLACPY,   ZLASCL
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC DBLE , DCMPLX , MAX
!     ..
!     .. Executable Statements ..
!
      zqrt17 = ZERO
!
      IF ( LSAME(Trans,'N') ) THEN
         nrows = M
         ncols = N
      ELSEIF ( LSAME(Trans,'C') ) THEN
         nrows = N
         ncols = M
      ELSE
         CALL xerbla('zqrt17',1)
         RETURN
      ENDIF
!
      IF ( Lwork<ncols*Nrhs ) THEN
         CALL xerbla('zqrt17',13)
         RETURN
      ENDIF
!
      IF ( M<=0 .OR. N<=0 .OR. Nrhs<=0 ) RETURN
!
      norma = ZLANGE('One-norm',M,N,A,Lda,rwork)
      smlnum = DLAMCH('Safe minimum')/DLAMCH('Precision')
      bignum = ONE/smlnum
      iscl = 0
!
!     compute residual and scale it
!
      CALL ZLACPY('All',nrows,Nrhs,B,Ldb,C,Ldb)
      CALL ZGEMM(Trans,'No transpose',nrows,Nrhs,ncols,DCMPLX(-ONE),A,  &
     &           Lda,X,Ldx,DCMPLX(ONE),C,Ldb)
      normrs = ZLANGE('Max',nrows,Nrhs,C,Ldb,rwork)
      IF ( normrs>smlnum ) THEN
         iscl = 1
         CALL ZLASCL('General',0,0,normrs,ONE,nrows,Nrhs,C,Ldb,info)
      ENDIF
!
!     compute R'*A
!
      CALL ZGEMM('Conjugate transpose',Trans,Nrhs,ncols,nrows,          &
     &           DCMPLX(ONE),C,Ldb,A,Lda,DCMPLX(ZERO),Work,Nrhs)
!
!     compute and properly scale error
!
      err = ZLANGE('One-norm',Nrhs,ncols,Work,Nrhs,rwork)
      IF ( norma/=ZERO ) err = err/norma
!
      IF ( iscl==1 ) err = err*normrs
!
      IF ( Iresid==1 ) THEN
         normb = ZLANGE('One-norm',nrows,Nrhs,B,Ldb,rwork)
         IF ( normb/=ZERO ) err = err/normb
      ELSE
         IF ( normrs/=ZERO ) err = err/normrs
      ENDIF
!
      zqrt17 = err/(DLAMCH('Epsilon')*DBLE(MAX(M,N,Nrhs)))
!
!     End of zqrt17
!
      END FUNCTION zqrt17
