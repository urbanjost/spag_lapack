!*==zpot03.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b ZPOT03
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZPOT03( UPLO, N, A, LDA, AINV, LDAINV, WORK, LDWORK,
!                          RWORK, RCOND, RESID )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
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
!> ZPOT03 computes the residual for a Hermitian matrix times its
!> inverse:
!>    norm( I - A*AINV ) / ( N * norm(A) * norm(AINV) * EPS ),
!> where EPS is the machine epsilon.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          Specifies whether the upper or lower triangular part of the
!>          Hermitian matrix A is stored:
!>          = 'U':  Upper triangular
!>          = 'L':  Lower triangular
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of rows and columns of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension (LDA,N)
!>          The original Hermitian matrix A.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,N)
!> \endverbatim
!>
!> \param[in,out] AINV
!> \verbatim
!>          AINV is COMPLEX*16 array, dimension (LDAINV,N)
!>          On entry, the inverse of the matrix A, stored as a Hermitian
!>          matrix in the same format as A.
!>          In this version, AINV is expanded into a full matrix and
!>          multiplied by A, so the opposing triangle of AINV will be
!>          changed; i.e., if the upper triangular part of AINV is
!>          stored, the lower triangular part will be used as work space.
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
!>          norm(I - A*AINV) / ( N * norm(A) * norm(AINV) * EPS )
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
      SUBROUTINE ZPOT03(Uplo,N,A,Lda,Ainv,Ldainv,Work,Ldwork,Rwork,     &
     &                  Rcond,Resid)
      IMPLICIT NONE
!*--ZPOT03130
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Uplo
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
      INTEGER i , j
      DOUBLE PRECISION ainvnm , anorm , eps
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      DOUBLE PRECISION DLAMCH , ZLANGE , ZLANHE
      EXTERNAL LSAME , DLAMCH , ZLANGE , ZLANHE
!     ..
!     .. External Subroutines ..
      EXTERNAL ZHEMM
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC DBLE , DCONJG
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
      anorm = ZLANHE('1',Uplo,N,A,Lda,Rwork)
      ainvnm = ZLANHE('1',Uplo,N,Ainv,Ldainv,Rwork)
      IF ( anorm<=ZERO .OR. ainvnm<=ZERO ) THEN
         Rcond = ZERO
         Resid = ONE/eps
         RETURN
      ENDIF
      Rcond = (ONE/anorm)/ainvnm
!
!     Expand AINV into a full matrix and call ZHEMM to multiply
!     AINV on the left by A.
!
      IF ( LSAME(Uplo,'U') ) THEN
         DO j = 1 , N
            DO i = 1 , j - 1
               Ainv(j,i) = DCONJG(Ainv(i,j))
            ENDDO
         ENDDO
      ELSE
         DO j = 1 , N
            DO i = j + 1 , N
               Ainv(j,i) = DCONJG(Ainv(i,j))
            ENDDO
         ENDDO
      ENDIF
      CALL ZHEMM('Left',Uplo,N,N,-CONE,A,Lda,Ainv,Ldainv,CZERO,Work,    &
     &           Ldwork)
!
!     Add the identity matrix to WORK .
!
      DO i = 1 , N
         Work(i,i) = Work(i,i) + CONE
      ENDDO
!
!     Compute norm(I - A*AINV) / (N * norm(A) * norm(AINV) * EPS)
!
      Resid = ZLANGE('1',N,N,Work,Ldwork,Rwork)
!
      Resid = ((Resid*Rcond)/eps)/DBLE(N)
!
!
!     End of ZPOT03
!
      END SUBROUTINE ZPOT03
