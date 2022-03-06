!*==zlarhs.f90  processed by SPAG 7.51RB at 17:35 on  4 Mar 2022
!> \brief \b zlarhs
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE zlarhs( PATH, XTYPE, UPLO, TRANS, M, N, KL, KU, NRHS,
!                          A, LDA, X, LDX, B, LDB, ISEED, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          TRANS, UPLO, XTYPE
!       CHARACTER*3        PATH
!       INTEGER            INFO, KL, KU, LDA, LDB, LDX, M, N, NRHS
!       ..
!       .. Array Arguments ..
!       INTEGER            ISEED( 4 )
!       COMPLEX*16         A( LDA, * ), B( LDB, * ), X( LDX, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> zlarhs chooses a set of NRHS random solution vectors and sets
!> up the right hand sides for the linear system
!>    op( A ) * X = B,
!> where op( A ) may be A, A**T (transpose of A), or A**H (conjugate
!> transpose of A).
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] PATH
!> \verbatim
!>          PATH is CHARACTER*3
!>          The type of the complex matrix A.  PATH may be given in any
!>          combination of upper and lower case.  Valid paths include
!>             xGE:  General m x n matrix
!>             xGB:  General banded matrix
!>             xPO:  Hermitian positive definite, 2-D storage
!>             xPP:  Hermitian positive definite packed
!>             xPB:  Hermitian positive definite banded
!>             xHE:  Hermitian indefinite, 2-D storage
!>             xHP:  Hermitian indefinite packed
!>             xHB:  Hermitian indefinite banded
!>             xSY:  Symmetric indefinite, 2-D storage
!>             xSP:  Symmetric indefinite packed
!>             xSB:  Symmetric indefinite banded
!>             xTR:  Triangular
!>             xTP:  Triangular packed
!>             xTB:  Triangular banded
!>             xQR:  General m x n matrix
!>             xLQ:  General m x n matrix
!>             xQL:  General m x n matrix
!>             xRQ:  General m x n matrix
!>          where the leading character indicates the precision.
!> \endverbatim
!>
!> \param[in] XTYPE
!> \verbatim
!>          XTYPE is CHARACTER*1
!>          Specifies how the exact solution X will be determined:
!>          = 'N':  New solution; generate a random X.
!>          = 'C':  Computed; use value of X on entry.
!> \endverbatim
!>
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          Used only if A is symmetric or triangular; specifies whether
!>          the upper or lower triangular part of the matrix A is stored.
!>          = 'U':  Upper triangular
!>          = 'L':  Lower triangular
!> \endverbatim
!>
!> \param[in] TRANS
!> \verbatim
!>          TRANS is CHARACTER*1
!>          Used only if A is nonsymmetric; specifies the operation
!>          applied to the matrix A.
!>          = 'N':  B := A    * X
!>          = 'T':  B := A**T * X
!>          = 'C':  B := A**H * X
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
!>          Used only if A is a band matrix; specifies the number of
!>          subdiagonals of A if A is a general band matrix or if A is
!>          symmetric or triangular and UPLO = 'L'; specifies the number
!>          of superdiagonals of A if A is symmetric or triangular and
!>          UPLO = 'U'.  0 <= KL <= M-1.
!> \endverbatim
!>
!> \param[in] KU
!> \verbatim
!>          KU is INTEGER
!>          Used only if A is a general band matrix or if A is
!>          triangular.
!>
!>          If PATH = xGB, specifies the number of superdiagonals of A,
!>          and 0 <= KU <= N-1.
!>
!>          If PATH = xTR, xTP, or xTB, specifies whether or not the
!>          matrix has unit diagonal:
!>          = 1:  matrix has non-unit diagonal (default)
!>          = 2:  matrix has unit diagonal
!> \endverbatim
!>
!> \param[in] NRHS
!> \verbatim
!>          NRHS is INTEGER
!>          The number of right hand side vectors in the system A*X = B.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension (LDA,N)
!>          The test matrix whose type is given by PATH.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.
!>          If PATH = xGB, LDA >= KL+KU+1.
!>          If PATH = xPB, xSB, xHB, or xTB, LDA >= KL+1.
!>          Otherwise, LDA >= max(1,M).
!> \endverbatim
!>
!> \param[in,out] X
!> \verbatim
!>          X is or output) COMPLEX*16 array, dimension (LDX,NRHS)
!>          On entry, if XTYPE = 'C' (for 'Computed'), then X contains
!>          the exact solution to the system of linear equations.
!>          On exit, if XTYPE = 'N' (for 'New'), then X is initialized
!>          with random values.
!> \endverbatim
!>
!> \param[in] LDX
!> \verbatim
!>          LDX is INTEGER
!>          The leading dimension of the array X.  If TRANS = 'N',
!>          LDX >= max(1,N); if TRANS = 'T', LDX >= max(1,M).
!> \endverbatim
!>
!> \param[out] B
!> \verbatim
!>          B is COMPLEX*16 array, dimension (LDB,NRHS)
!>          The right hand side vector(s) for the system of equations,
!>          computed from B = op(A) * X, where op(A) is determined by
!>          TRANS.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of the array B.  If TRANS = 'N',
!>          LDB >= max(1,M); if TRANS = 'T', LDB >= max(1,N).
!> \endverbatim
!>
!> \param[in,out] ISEED
!> \verbatim
!>          ISEED is INTEGER array, dimension (4)
!>          The seed vector for the random number generator (used in
!>          ZLATMS).  Modified on exit.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0: successful exit
!>          < 0: if INFO = -i, the i-th argument had an illegal value
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
!> \date June 2017
!
!> \ingroup complex16_lin
!
!  =====================================================================
      SUBROUTINE zlarhs(Path,Xtype,Uplo,Trans,M,N,Kl,Ku,Nrhs,A,Lda,X,   &
     &                  Ldx,B,Ldb,Iseed,Info)
      IMPLICIT NONE
!*--ZLARHS213
!
!  -- LAPACK test routine (version 3.7.1) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     June 2017
!
!     .. Scalar Arguments ..
      CHARACTER Trans , Uplo , Xtype
      CHARACTER*3 Path
      INTEGER Info , Kl , Ku , Lda , Ldb , Ldx , M , N , Nrhs
!     ..
!     .. Array Arguments ..
      INTEGER Iseed(4)
      COMPLEX*16 A(Lda,*) , B(Ldb,*) , X(Ldx,*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      COMPLEX*16 ONE , ZERO
      PARAMETER (ONE=(1.0D+0,0.0D+0),ZERO=(0.0D+0,0.0D+0))
!     ..
!     .. Local Scalars ..
      LOGICAL band , gen , notran , qrs , sym , tran , tri
      CHARACTER c1 , diag
      CHARACTER*2 c2
      INTEGER j , mb , nx
!     ..
!     .. External Functions ..
      LOGICAL LSAME , LSAMEN
EXTERNAL    LSAME,   LSAMEN
!     ..
!     .. External Subroutines ..
EXTERNAL    ZGBMV,   ZGEMM,   ZHBMV,   ZHEMM,   ZHPMV,   ZLACPY
EXTERNAL    ZLARNV,   ZSPMV,   ZSYMM,   ZTBMV,   ZTPMV,   ZTRMM
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      Info = 0
      c1 = Path(1:1)
      c2 = Path(2:3)
      tran = LSAME(Trans,'T') .OR. LSAME(Trans,'C')
      notran = .NOT.tran
      gen = LSAME(Path(2:2),'G')
      qrs = LSAME(Path(2:2),'Q') .OR. LSAME(Path(3:3),'Q')
      sym = LSAME(Path(2:2),'P') .OR. LSAME(Path(2:2),'S') .OR.         &
     &      LSAME(Path(2:2),'H')
      tri = LSAME(Path(2:2),'T')
      band = LSAME(Path(3:3),'B')
      IF ( .NOT.LSAME(c1,'Zomplex precision') ) THEN
         Info = -1
      ELSEIF ( .NOT.(LSAME(Xtype,'N') .OR. LSAME(Xtype,'C')) ) THEN
         Info = -2
      ELSEIF ( (sym .OR. tri) .AND.                                     &
     &         .NOT.(LSAME(Uplo,'U') .OR. LSAME(Uplo,'L')) ) THEN
         Info = -3
      ELSEIF ( (gen .OR. qrs) .AND. .NOT.(tran .OR. LSAME(Trans,'N')) ) &
     &         THEN
         Info = -4
      ELSEIF ( M<0 ) THEN
         Info = -5
      ELSEIF ( N<0 ) THEN
         Info = -6
      ELSEIF ( band .AND. Kl<0 ) THEN
         Info = -7
      ELSEIF ( band .AND. Ku<0 ) THEN
         Info = -8
      ELSEIF ( Nrhs<0 ) THEN
         Info = -9
      ELSEIF ( (.NOT.band .AND. Lda<MAX(1,M)) .OR.                      &
     &         (band .AND. (sym .OR. tri) .AND. Lda<Kl+1) .OR.          &
     &         (band .AND. gen .AND. Lda<Kl+Ku+1) ) THEN
         Info = -11
      ELSEIF ( (notran .AND. Ldx<MAX(1,N)) .OR.                         &
     &         (tran .AND. Ldx<MAX(1,M)) ) THEN
         Info = -13
      ELSEIF ( (notran .AND. Ldb<MAX(1,M)) .OR.                         &
     &         (tran .AND. Ldb<MAX(1,N)) ) THEN
         Info = -15
      ENDIF
      IF ( Info/=0 ) THEN
         CALL xerbla('zlarhs',-Info)
         RETURN
      ENDIF
!
!     Initialize X to NRHS random vectors unless XTYPE = 'C'.
!
      IF ( tran ) THEN
         nx = M
         mb = N
      ELSE
         nx = N
         mb = M
      ENDIF
      IF ( .NOT.LSAME(Xtype,'C') ) THEN
         DO j = 1 , Nrhs
            CALL ZLARNV(2,Iseed,N,X(1,j))
         ENDDO
      ENDIF
!
!     Multiply X by op( A ) using an appropriate
!     matrix multiply routine.
!
      IF ( LSAMEN(2,c2,'GE') .OR. LSAMEN(2,c2,'QR') .OR.                &
     &     LSAMEN(2,c2,'LQ') .OR. LSAMEN(2,c2,'QL') .OR.                &
     &     LSAMEN(2,c2,'RQ') ) THEN
!
!        General matrix
!
         CALL ZGEMM(Trans,'N',mb,Nrhs,nx,ONE,A,Lda,X,Ldx,ZERO,B,Ldb)
!
      ELSEIF ( LSAMEN(2,c2,'PO') .OR. LSAMEN(2,c2,'HE') ) THEN
!
!        Hermitian matrix, 2-D storage
!
         CALL ZHEMM('Left',Uplo,N,Nrhs,ONE,A,Lda,X,Ldx,ZERO,B,Ldb)
!
      ELSEIF ( LSAMEN(2,c2,'SY') ) THEN
!
!        Symmetric matrix, 2-D storage
!
         CALL ZSYMM('Left',Uplo,N,Nrhs,ONE,A,Lda,X,Ldx,ZERO,B,Ldb)
!
      ELSEIF ( LSAMEN(2,c2,'GB') ) THEN
!
!        General matrix, band storage
!
         DO j = 1 , Nrhs
            CALL ZGBMV(Trans,M,N,Kl,Ku,ONE,A,Lda,X(1,j),1,ZERO,B(1,j),1)
         ENDDO
!
      ELSEIF ( LSAMEN(2,c2,'PB') .OR. LSAMEN(2,c2,'HB') ) THEN
!
!        Hermitian matrix, band storage
!
         DO j = 1 , Nrhs
            CALL ZHBMV(Uplo,N,Kl,ONE,A,Lda,X(1,j),1,ZERO,B(1,j),1)
         ENDDO
!
      ELSEIF ( LSAMEN(2,c2,'SB') ) THEN
!
!        Symmetric matrix, band storage
!
         DO j = 1 , Nrhs
            CALL zsbmv(Uplo,N,Kl,ONE,A,Lda,X(1,j),1,ZERO,B(1,j),1)
         ENDDO
!
      ELSEIF ( LSAMEN(2,c2,'PP') .OR. LSAMEN(2,c2,'HP') ) THEN
!
!        Hermitian matrix, packed storage
!
         DO j = 1 , Nrhs
            CALL ZHPMV(Uplo,N,ONE,A,X(1,j),1,ZERO,B(1,j),1)
         ENDDO
!
      ELSEIF ( LSAMEN(2,c2,'SP') ) THEN
!
!        Symmetric matrix, packed storage
!
         DO j = 1 , Nrhs
            CALL ZSPMV(Uplo,N,ONE,A,X(1,j),1,ZERO,B(1,j),1)
         ENDDO
!
      ELSEIF ( LSAMEN(2,c2,'TR') ) THEN
!
!        Triangular matrix.  Note that for triangular matrices,
!           KU = 1 => non-unit triangular
!           KU = 2 => unit triangular
!
         CALL ZLACPY('Full',N,Nrhs,X,Ldx,B,Ldb)
         IF ( Ku==2 ) THEN
            diag = 'U'
         ELSE
            diag = 'N'
         ENDIF
         CALL ZTRMM('Left',Uplo,Trans,diag,N,Nrhs,ONE,A,Lda,B,Ldb)
!
      ELSEIF ( LSAMEN(2,c2,'TP') ) THEN
!
!        Triangular matrix, packed storage
!
         CALL ZLACPY('Full',N,Nrhs,X,Ldx,B,Ldb)
         IF ( Ku==2 ) THEN
            diag = 'U'
         ELSE
            diag = 'N'
         ENDIF
         DO j = 1 , Nrhs
            CALL ZTPMV(Uplo,Trans,diag,N,A,B(1,j),1)
         ENDDO
!
      ELSEIF ( LSAMEN(2,c2,'TB') ) THEN
!
!        Triangular matrix, banded storage
!
         CALL ZLACPY('Full',N,Nrhs,X,Ldx,B,Ldb)
         IF ( Ku==2 ) THEN
            diag = 'U'
         ELSE
            diag = 'N'
         ENDIF
         DO j = 1 , Nrhs
            CALL ZTBMV(Uplo,Trans,diag,N,Kl,A,Lda,B(1,j),1)
         ENDDO
!
      ELSE
!
!        If none of the above, set INFO = -1 and return
!
         Info = -1
         CALL xerbla('zlarhs',-Info)
      ENDIF
!
!
!     End of zlarhs
!
      END SUBROUTINE zlarhs
