!*==dlarhs.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b DLARHS
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLARHS( PATH, XTYPE, UPLO, TRANS, M, N, KL, KU, NRHS,
!                          A, LDA, X, LDX, B, LDB, ISEED, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          TRANS, UPLO, XTYPE
!       CHARACTER*3        PATH
!       INTEGER            INFO, KL, KU, LDA, LDB, LDX, M, N, NRHS
!       ..
!       .. Array Arguments ..
!       INTEGER            ISEED( 4 )
!       DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), X( LDX, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLARHS chooses a set of NRHS random solution vectors and sets
!> up the right hand sides for the linear system
!>    op( A ) * X = B,
!> where op( A ) may be A or A' (transpose of A).
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] PATH
!> \verbatim
!>          PATH is CHARACTER*3
!>          The type of the real matrix A.  PATH may be given in any
!>          combination of upper and lower case.  Valid types include
!>             xGE:  General m x n matrix
!>             xGB:  General banded matrix
!>             xPO:  Symmetric positive definite, 2-D storage
!>             xPP:  Symmetric positive definite packed
!>             xPB:  Symmetric positive definite banded
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
!>          Specifies whether the upper or lower triangular part of the
!>          matrix A is stored, if A is symmetric.
!>          = 'U':  Upper triangular
!>          = 'L':  Lower triangular
!> \endverbatim
!>
!> \param[in] TRANS
!> \verbatim
!>          TRANS is CHARACTER*1
!>          Specifies the operation applied to the matrix A.
!>          = 'N':  System is  A * x = b
!>          = 'T':  System is  A'* x = b
!>          = 'C':  System is  A'* x = b
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number or rows of the matrix A.  M >= 0.
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
!>          A is DOUBLE PRECISION array, dimension (LDA,N)
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
!>          X is or output) DOUBLE PRECISION array, dimension(LDX,NRHS)
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
!>          B is DOUBLE PRECISION array, dimension (LDB,NRHS)
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
!>          DLATMS).  Modified on exit.
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
!> \date December 2016
!
!> \ingroup double_lin
!
!  =====================================================================
      SUBROUTINE DLARHS(Path,Xtype,Uplo,Trans,M,N,Kl,Ku,Nrhs,A,Lda,X,   &
     &                  Ldx,B,Ldb,Iseed,Info)
      IMPLICIT NONE
!*--DLARHS208
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Trans , Uplo , Xtype
      CHARACTER*3 Path
      INTEGER Info , Kl , Ku , Lda , Ldb , Ldx , M , N , Nrhs
!     ..
!     .. Array Arguments ..
      INTEGER Iseed(4)
      DOUBLE PRECISION A(Lda,*) , B(Ldb,*) , X(Ldx,*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ONE , ZERO
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
!     ..
!     .. Local Scalars ..
      LOGICAL band , gen , notran , qrs , sym , tran , tri
      CHARACTER c1 , diag
      CHARACTER*2 c2
      INTEGER j , mb , nx
!     ..
!     .. External Functions ..
      LOGICAL LSAME , LSAMEN
      EXTERNAL LSAME , LSAMEN
!     ..
!     .. External Subroutines ..
      EXTERNAL DGBMV , DGEMM , DLACPY , DLARNV , DSBMV , DSPMV , DSYMM ,&
     &         DTBMV , DTPMV , DTRMM , XERBLA
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
      sym = LSAME(Path(2:2),'P') .OR. LSAME(Path(2:2),'S')
      tri = LSAME(Path(2:2),'T')
      band = LSAME(Path(3:3),'B')
      IF ( .NOT.LSAME(c1,'Double precision') ) THEN
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
         CALL XERBLA('DLARHS',-Info)
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
            CALL DLARNV(2,Iseed,N,X(1,j))
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
         CALL DGEMM(Trans,'N',mb,Nrhs,nx,ONE,A,Lda,X,Ldx,ZERO,B,Ldb)
!
      ELSEIF ( LSAMEN(2,c2,'PO') .OR. LSAMEN(2,c2,'SY') ) THEN
!
!        Symmetric matrix, 2-D storage
!
         CALL DSYMM('Left',Uplo,N,Nrhs,ONE,A,Lda,X,Ldx,ZERO,B,Ldb)
!
      ELSEIF ( LSAMEN(2,c2,'GB') ) THEN
!
!        General matrix, band storage
!
         DO j = 1 , Nrhs
            CALL DGBMV(Trans,mb,nx,Kl,Ku,ONE,A,Lda,X(1,j),1,ZERO,B(1,j),&
     &                 1)
         ENDDO
!
      ELSEIF ( LSAMEN(2,c2,'PB') ) THEN
!
!        Symmetric matrix, band storage
!
         DO j = 1 , Nrhs
            CALL DSBMV(Uplo,N,Kl,ONE,A,Lda,X(1,j),1,ZERO,B(1,j),1)
         ENDDO
!
      ELSEIF ( LSAMEN(2,c2,'PP') .OR. LSAMEN(2,c2,'SP') ) THEN
!
!        Symmetric matrix, packed storage
!
         DO j = 1 , Nrhs
            CALL DSPMV(Uplo,N,ONE,A,X(1,j),1,ZERO,B(1,j),1)
         ENDDO
!
      ELSEIF ( LSAMEN(2,c2,'TR') ) THEN
!
!        Triangular matrix.  Note that for triangular matrices,
!           KU = 1 => non-unit triangular
!           KU = 2 => unit triangular
!
         CALL DLACPY('Full',N,Nrhs,X,Ldx,B,Ldb)
         IF ( Ku==2 ) THEN
            diag = 'U'
         ELSE
            diag = 'N'
         ENDIF
         CALL DTRMM('Left',Uplo,Trans,diag,N,Nrhs,ONE,A,Lda,B,Ldb)
!
      ELSEIF ( LSAMEN(2,c2,'TP') ) THEN
!
!        Triangular matrix, packed storage
!
         CALL DLACPY('Full',N,Nrhs,X,Ldx,B,Ldb)
         IF ( Ku==2 ) THEN
            diag = 'U'
         ELSE
            diag = 'N'
         ENDIF
         DO j = 1 , Nrhs
            CALL DTPMV(Uplo,Trans,diag,N,A,B(1,j),1)
         ENDDO
!
      ELSEIF ( LSAMEN(2,c2,'TB') ) THEN
!
!        Triangular matrix, banded storage
!
         CALL DLACPY('Full',N,Nrhs,X,Ldx,B,Ldb)
         IF ( Ku==2 ) THEN
            diag = 'U'
         ELSE
            diag = 'N'
         ENDIF
         DO j = 1 , Nrhs
            CALL DTBMV(Uplo,Trans,diag,N,Kl,A,Lda,B(1,j),1)
         ENDDO
!
      ELSE
!
!        If PATH is none of the above, return with an error code.
!
         Info = -1
         CALL XERBLA('DLARHS',-Info)
      ENDIF
!
!
!     End of DLARHS
!
      END SUBROUTINE DLARHS
