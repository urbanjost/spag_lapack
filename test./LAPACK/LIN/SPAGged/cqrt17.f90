!*==cqrt17.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b CQRT17
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       REAL             FUNCTION CQRT17( TRANS, IRESID, M, N, NRHS, A,
!                        LDA, X, LDX, B, LDB, C, WORK, LWORK )
!
!       .. Scalar Arguments ..
!       CHARACTER          TRANS
!       INTEGER            IRESID, LDA, LDB, LDX, LWORK, M, N, NRHS
!       ..
!       .. Array Arguments ..
!       COMPLEX            A( LDA, * ), B( LDB, * ), C( LDB, * ),
!      $                   WORK( LWORK ), X( LDX, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CQRT17 computes the ratio
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
!>          A is COMPLEX array, dimension (LDA,N)
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
!>          X is COMPLEX array, dimension (LDX,NRHS)
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
!>          B is COMPLEX array, dimension (LDB,NRHS)
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
!>          C is COMPLEX array, dimension (LDB,NRHS)
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX array, dimension (LWORK)
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
!> \ingroup complex_lin
!
!  =====================================================================
      REAL FUNCTION CQRT17(Trans,Iresid,M,N,Nrhs,A,Lda,X,Ldx,B,Ldb,C,   &
     &                     Work,Lwork)
      IMPLICIT NONE
!*--CQRT17154
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
      COMPLEX A(Lda,*) , B(Ldb,*) , C(Ldb,*) , Work(Lwork) , X(Ldx,*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ZERO , ONE
      PARAMETER (ZERO=0.0E0,ONE=1.0E0)
!     ..
!     .. Local Scalars ..
      INTEGER info , iscl , ncols , nrows
      REAL bignum , err , norma , normb , normrs , smlnum
!     ..
!     .. Local Arrays ..
      REAL rwork(1)
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      REAL CLANGE , SLAMCH
      EXTERNAL LSAME , CLANGE , SLAMCH
!     ..
!     .. External Subroutines ..
      EXTERNAL CGEMM , CLACPY , CLASCL , XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC CMPLX , MAX , REAL
!     ..
!     .. Executable Statements ..
!
      CQRT17 = ZERO
!
      IF ( LSAME(Trans,'N') ) THEN
         nrows = M
         ncols = N
      ELSEIF ( LSAME(Trans,'C') ) THEN
         nrows = N
         ncols = M
      ELSE
         CALL XERBLA('CQRT17',1)
         RETURN
      ENDIF
!
      IF ( Lwork<ncols*Nrhs ) THEN
         CALL XERBLA('CQRT17',13)
         RETURN
      ENDIF
!
      IF ( M<=0 .OR. N<=0 .OR. Nrhs<=0 ) RETURN
!
      norma = CLANGE('One-norm',M,N,A,Lda,rwork)
      smlnum = SLAMCH('Safe minimum')/SLAMCH('Precision')
      bignum = ONE/smlnum
      iscl = 0
!
!     compute residual and scale it
!
      CALL CLACPY('All',nrows,Nrhs,B,Ldb,C,Ldb)
      CALL CGEMM(Trans,'No transpose',nrows,Nrhs,ncols,CMPLX(-ONE),A,   &
     &           Lda,X,Ldx,CMPLX(ONE),C,Ldb)
      normrs = CLANGE('Max',nrows,Nrhs,C,Ldb,rwork)
      IF ( normrs>smlnum ) THEN
         iscl = 1
         CALL CLASCL('General',0,0,normrs,ONE,nrows,Nrhs,C,Ldb,info)
      ENDIF
!
!     compute R'*A
!
      CALL CGEMM('Conjugate transpose',Trans,Nrhs,ncols,nrows,CMPLX(ONE)&
     &           ,C,Ldb,A,Lda,CMPLX(ZERO),Work,Nrhs)
!
!     compute and properly scale error
!
      err = CLANGE('One-norm',Nrhs,ncols,Work,Nrhs,rwork)
      IF ( norma/=ZERO ) err = err/norma
!
      IF ( iscl==1 ) err = err*normrs
!
      IF ( Iresid==1 ) THEN
         normb = CLANGE('One-norm',nrows,Nrhs,B,Ldb,rwork)
         IF ( normb/=ZERO ) err = err/normb
      ELSE
         IF ( normrs/=ZERO ) err = err/normrs
      ENDIF
!
      CQRT17 = err/(SLAMCH('Epsilon')*REAL(MAX(M,N,Nrhs)))
!
!     End of CQRT17
!
      END FUNCTION CQRT17
