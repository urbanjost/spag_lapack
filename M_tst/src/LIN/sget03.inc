!*==sget03.f90  processed by SPAG 7.51RB at 17:34 on  4 Mar 2022
!> \brief \b sget03
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE sget03( N, A, LDA, AINV, LDAINV, WORK, LDWORK, RWORK,
!                          RCOND, RESID )
!
!       .. Scalar Arguments ..
!       INTEGER            LDA, LDAINV, LDWORK, N
!       REAL               RCOND, RESID
!       ..
!       .. Array Arguments ..
!       REAL               A( LDA, * ), AINV( LDAINV, * ), RWORK( * ),
!      $                   WORK( LDWORK, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> sget03 computes the residual for a general matrix times its inverse:
!>    norm( I - AINV*A ) / ( N * norm(A) * norm(AINV) * EPS ),
!> where EPS is the machine epsilon.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of rows and columns of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is REAL array, dimension (LDA,N)
!>          The original N x N matrix A.
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
!>          AINV is REAL array, dimension (LDAINV,N)
!>          The inverse of the matrix A.
!> \endverbatim
!>
!> \param[in] LDAINV
!> \verbatim
!>          LDAINV is INTEGER
!>          The leading dimension of the array AINV.  LDAINV >= max(1,N).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is REAL array, dimension (LDWORK,N)
!> \endverbatim
!>
!> \param[in] LDWORK
!> \verbatim
!>          LDWORK is INTEGER
!>          The leading dimension of the array WORK.  LDWORK >= max(1,N).
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is REAL array, dimension (N)
!> \endverbatim
!>
!> \param[out] RCOND
!> \verbatim
!>          RCOND is REAL
!>          The reciprocal of the condition number of A, computed as
!>          ( 1/norm(A) ) / norm(AINV).
!> \endverbatim
!>
!> \param[out] RESID
!> \verbatim
!>          RESID is REAL
!>          norm(I - AINV*A) / ( N * norm(A) * norm(AINV) * EPS )
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
!> \ingroup single_lin
!
!  =====================================================================
      SUBROUTINE sget03(N,A,Lda,Ainv,Ldainv,Work,Ldwork,Rwork,Rcond,    &
     &                  Resid)
      IMPLICIT NONE
!*--SGET03113
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Lda , Ldainv , Ldwork , N
      REAL Rcond , Resid
!     ..
!     .. Array Arguments ..
      REAL A(Lda,*) , Ainv(Ldainv,*) , Rwork(*) , Work(Ldwork,*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ZERO , ONE
      PARAMETER (ZERO=0.0E+0,ONE=1.0E+0)
!     ..
!     .. Local Scalars ..
      INTEGER i
      REAL ainvnm , anorm , eps
!     ..
!     .. External Functions ..
      REAL SLAMCH , SLANGE
EXTERNAL    SLAMCH,   SLANGE
!     ..
!     .. External Subroutines ..
EXTERNAL    SGEMM
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC REAL
!     ..
!     .. Executable Statements ..
!
!     Quick exit if N = 0.
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
      anorm = SLANGE('1',N,N,A,Lda,Rwork)
      ainvnm = SLANGE('1',N,N,Ainv,Ldainv,Rwork)
      IF ( anorm<=ZERO .OR. ainvnm<=ZERO ) THEN
         Rcond = ZERO
         Resid = ONE/eps
         RETURN
      ENDIF
      Rcond = (ONE/anorm)/ainvnm
!
!     Compute I - A * AINV
!
      CALL SGEMM('No transpose','No transpose',N,N,N,-ONE,Ainv,Ldainv,A,&
     &           Lda,ZERO,Work,Ldwork)
      DO i = 1 , N
         Work(i,i) = ONE + Work(i,i)
      ENDDO
!
!     Compute norm(I - AINV*A) / (N * norm(A) * norm(AINV) * EPS)
!
      Resid = SLANGE('1',N,N,Work,Ldwork,Rwork)
!
      Resid = ((Resid*Rcond)/eps)/REAL(N)
!
!
!     End of sget03
!
      END SUBROUTINE sget03
