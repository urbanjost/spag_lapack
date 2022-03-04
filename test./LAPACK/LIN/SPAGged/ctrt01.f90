!*==ctrt01.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b CTRT01
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE CTRT01( UPLO, DIAG, N, A, LDA, AINV, LDAINV, RCOND,
!                          RWORK, RESID )
!
!       .. Scalar Arguments ..
!       CHARACTER          DIAG, UPLO
!       INTEGER            LDA, LDAINV, N
!       REAL               RCOND, RESID
!       ..
!       .. Array Arguments ..
!       REAL               RWORK( * )
!       COMPLEX            A( LDA, * ), AINV( LDAINV, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CTRT01 computes the residual for a triangular matrix A times its
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
!>          A is COMPLEX array, dimension (LDA,N)
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
!> \param[in] AINV
!> \verbatim
!>          AINV is COMPLEX array, dimension (LDAINV,N)
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
!>          RCOND is REAL
!>          The reciprocal condition number of A, computed as
!>          1/(norm(A) * norm(AINV)).
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is REAL array, dimension (N)
!> \endverbatim
!>
!> \param[out] RESID
!> \verbatim
!>          RESID is REAL
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
!> \ingroup complex_lin
!
!  =====================================================================
      SUBROUTINE CTRT01(Uplo,Diag,N,A,Lda,Ainv,Ldainv,Rcond,Rwork,Resid)
      IMPLICIT NONE
!*--CTRT01128
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Diag , Uplo
      INTEGER Lda , Ldainv , N
      REAL Rcond , Resid
!     ..
!     .. Array Arguments ..
      REAL Rwork(*)
      COMPLEX A(Lda,*) , Ainv(Ldainv,*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ZERO , ONE
      PARAMETER (ZERO=0.0E+0,ONE=1.0E+0)
!     ..
!     .. Local Scalars ..
      INTEGER j
      REAL ainvnm , anorm , eps
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      REAL CLANTR , SLAMCH
      EXTERNAL LSAME , CLANTR , SLAMCH
!     ..
!     .. External Subroutines ..
      EXTERNAL CTRMV
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC REAL
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
      eps = SLAMCH('Epsilon')
      anorm = CLANTR('1',Uplo,Diag,N,N,A,Lda,Rwork)
      ainvnm = CLANTR('1',Uplo,Diag,N,N,Ainv,Ldainv,Rwork)
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
            CALL CTRMV('Upper','No transpose',Diag,j,A,Lda,Ainv(1,j),1)
         ENDDO
      ELSE
         DO j = 1 , N
            CALL CTRMV('Lower','No transpose',Diag,N-j+1,A(j,j),Lda,    &
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
      Resid = CLANTR('1',Uplo,'Non-unit',N,N,Ainv,Ldainv,Rwork)
!
      Resid = ((Resid*Rcond)/REAL(N))/eps
!
!
!     End of CTRT01
!
      END SUBROUTINE CTRT01
