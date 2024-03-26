!*==zget03.f90  processed by SPAG 7.51RB at 17:35 on  4 Mar 2022
!> \brief \b zget03
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE zget03( N, A, LDA, AINV, LDAINV, WORK, LDWORK, RWORK,
!                          RCOND, RESID )
!
!       .. Scalar Arguments ..
!       INTEGER            LDA, LDAINV, LDWORK, N
!       DOUBLE PRECISION   RCOND, RESID
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   RWORK( * )
!       COMPLEX*16         A( LDA, * ), AINV( LDAINV, * ),
!      $                   WORK( LDWORK, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> zget03 computes the residual for a general matrix times its inverse:
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
!>          A is COMPLEX*16 array, dimension (LDA,N)
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
!>          AINV is COMPLEX*16 array, dimension (LDAINV,N)
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
!>          WORK is COMPLEX*16 array, dimension (LDWORK,N)
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
!>          RWORK is DOUBLE PRECISION array, dimension (N)
!> \endverbatim
!>
!> \param[out] RCOND
!> \verbatim
!>          RCOND is DOUBLE PRECISION
!>          The reciprocal of the condition number of A, computed as
!>          ( 1/norm(A) ) / norm(AINV).
!> \endverbatim
!>
!> \param[out] RESID
!> \verbatim
!>          RESID is DOUBLE PRECISION
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
!> \ingroup complex16_lin
!
!  =====================================================================
      SUBROUTINE zget03(N,A,Lda,Ainv,Ldainv,Work,Ldwork,Rwork,Rcond,    &
     &                  Resid)
      IMPLICIT NONE
!*--ZGET03114
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Lda , Ldainv , Ldwork , N
      DOUBLE PRECISION Rcond , Resid
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION Rwork(*)
      COMPLEX*16 A(Lda,*) , Ainv(Ldainv,*) , Work(Ldwork,*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO , ONE
      PARAMETER (ZERO=0.0D+0,ONE=1.0D+0)
      COMPLEX*16 CZERO , CONE
      PARAMETER (CZERO=(0.0D+0,0.0D+0),CONE=(1.0D+0,0.0D+0))
!     ..
!     .. Local Scalars ..
      INTEGER i
      DOUBLE PRECISION ainvnm , anorm , eps
!     ..
!     .. External Functions ..
      DOUBLE PRECISION DLAMCH , ZLANGE
EXTERNAL    DLAMCH,   ZLANGE
!     ..
!     .. External Subroutines ..
EXTERNAL    ZGEMM
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC DBLE
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
      eps = DLAMCH('Epsilon')
      anorm = ZLANGE('1',N,N,A,Lda,Rwork)
      ainvnm = ZLANGE('1',N,N,Ainv,Ldainv,Rwork)
      IF ( anorm<=ZERO .OR. ainvnm<=ZERO ) THEN
         Rcond = ZERO
         Resid = ONE/eps
         RETURN
      ENDIF
      Rcond = (ONE/anorm)/ainvnm
!
!     Compute I - A * AINV
!
      CALL ZGEMM('No transpose','No transpose',N,N,N,-CONE,Ainv,Ldainv, &
     &           A,Lda,CZERO,Work,Ldwork)
      DO i = 1 , N
         Work(i,i) = CONE + Work(i,i)
      ENDDO
!
!     Compute norm(I - AINV*A) / (N * norm(A) * norm(AINV) * EPS)
!
      Resid = ZLANGE('1',N,N,Work,Ldwork,Rwork)
!
      Resid = ((Resid*Rcond)/eps)/DBLE(N)
!
!
!     End of zget03
!
      END SUBROUTINE zget03
