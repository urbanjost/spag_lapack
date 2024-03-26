!*==dtrt01.f90  processed by SPAG 7.51RB at 17:34 on  4 Mar 2022
!> \brief \b dtrt01
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE dtrt01( UPLO, DIAG, N, A, LDA, AINV, LDAINV, RCOND,
!                          WORK, RESID )
!
!       .. Scalar Arguments ..
!       CHARACTER          DIAG, UPLO
!       INTEGER            LDA, LDAINV, N
!       DOUBLE PRECISION   RCOND, RESID
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   A( LDA, * ), AINV( LDAINV, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> dtrt01 computes the residual for a triangular matrix A times its
!> inverse:
!>    RESID = norm( A*AINV - I ) / ( N * norm(A) * norm(AINV) * EPS ),
!> where EPS is the machine epsilon.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          Specifies whether the matrix A is upper or lower triangular.
!>          = 'U':  Upper triangular
!>          = 'L':  Lower triangular
!> \endverbatim
!>
!> \param[in] DIAG
!> \verbatim
!>          DIAG is CHARACTER*1
!>          Specifies whether or not the matrix A is unit triangular.
!>          = 'N':  Non-unit triangular
!>          = 'U':  Unit triangular
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (LDA,N)
!>          The triangular matrix A.  If UPLO = 'U', the leading n by n
!>          upper triangular part of the array A contains the upper
!>          triangular matrix, and the strictly lower triangular part of
!>          A is not referenced.  If UPLO = 'L', the leading n by n lower
!>          triangular part of the array A contains the lower triangular
!>          matrix, and the strictly upper triangular part of A is not
!>          referenced.  If DIAG = 'U', the diagonal elements of A are
!>          also not referenced and are assumed to be 1.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,N).
!> \endverbatim
!>
!> \param[in,out] AINV
!> \verbatim
!>          AINV is DOUBLE PRECISION array, dimension (LDAINV,N)
!>          On entry, the (triangular) inverse of the matrix A, in the
!>          same storage format as A.
!>          On exit, the contents of AINV are destroyed.
!> \endverbatim
!>
!> \param[in] LDAINV
!> \verbatim
!>          LDAINV is INTEGER
!>          The leading dimension of the array AINV.  LDAINV >= max(1,N).
!> \endverbatim
!>
!> \param[out] RCOND
!> \verbatim
!>          RCOND is DOUBLE PRECISION
!>          The reciprocal condition number of A, computed as
!>          1/(norm(A) * norm(AINV)).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (N)
!> \endverbatim
!>
!> \param[out] RESID
!> \verbatim
!>          RESID is DOUBLE PRECISION
!>          norm(A*AINV - I) / ( N * norm(A) * norm(AINV) * EPS )
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
      SUBROUTINE dtrt01(Uplo,Diag,N,A,Lda,Ainv,Ldainv,Rcond,Work,Resid)
      IMPLICIT NONE
!*--DTRT01127
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Diag , Uplo
      INTEGER Lda , Ldainv , N
      DOUBLE PRECISION Rcond , Resid
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION A(Lda,*) , Ainv(Ldainv,*) , Work(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO , ONE
      PARAMETER (ZERO=0.0D+0,ONE=1.0D+0)
!     ..
!     .. Local Scalars ..
      INTEGER j
      DOUBLE PRECISION ainvnm , anorm , eps
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      DOUBLE PRECISION DLAMCH , DLANTR
EXTERNAL    LSAME,   DLAMCH,   DLANTR
!     ..
!     .. External Subroutines ..
EXTERNAL    DTRMV
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC DBLE
!     ..
!     .. Executable Statements ..
!
!     Quick exit if N = 0
!
      IF ( N<=0 ) THEN
         Rcond = ONE
         Resid = ZERO
         RETURN
      ENDIF
!
!     Exit with RESID = 1/EPS if ANORM = 0 or AINVNM = 0.
!
      eps = DLAMCH('Epsilon')
      anorm = DLANTR('1',Uplo,Diag,N,N,A,Lda,Work)
      ainvnm = DLANTR('1',Uplo,Diag,N,N,Ainv,Ldainv,Work)
      IF ( anorm<=ZERO .OR. ainvnm<=ZERO ) THEN
         Rcond = ZERO
         Resid = ONE/eps
         RETURN
      ENDIF
      Rcond = (ONE/anorm)/ainvnm
!
!     Set the diagonal of AINV to 1 if AINV has unit diagonal.
!
      IF ( LSAME(Diag,'U') ) THEN
         DO j = 1 , N
            Ainv(j,j) = ONE
         ENDDO
      ENDIF
!
!     Compute A * AINV, overwriting AINV.
!
      IF ( LSAME(Uplo,'U') ) THEN
         DO j = 1 , N
            CALL DTRMV('Upper','No transpose',Diag,j,A,Lda,Ainv(1,j),1)
         ENDDO
      ELSE
         DO j = 1 , N
            CALL DTRMV('Lower','No transpose',Diag,N-j+1,A(j,j),Lda,    &
     &                 Ainv(j,j),1)
         ENDDO
      ENDIF
!
!     Subtract 1 from each diagonal element to form A*AINV - I.
!
      DO j = 1 , N
         Ainv(j,j) = Ainv(j,j) - ONE
      ENDDO
!
!     Compute norm(A*AINV - I) / (N * norm(A) * norm(AINV) * EPS)
!
      Resid = DLANTR('1',Uplo,'Non-unit',N,N,Ainv,Ldainv,Work)
!
      Resid = ((Resid*Rcond)/DBLE(N))/eps
!
!
!     End of dtrt01
!
      END SUBROUTINE dtrt01
