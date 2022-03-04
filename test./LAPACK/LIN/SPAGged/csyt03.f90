!*==csyt03.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b CSYT03
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE CSYT03( UPLO, N, A, LDA, AINV, LDAINV, WORK, LDWORK,
!                          RWORK, RCOND, RESID )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            LDA, LDAINV, LDWORK, N
!       REAL               RCOND, RESID
!       ..
!       .. Array Arguments ..
!       REAL               RWORK( * )
!       COMPLEX            A( LDA, * ), AINV( LDAINV, * ),
!      $                   WORK( LDWORK, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CSYT03 computes the residual for a complex symmetric matrix times
!> its inverse:
!>    norm( I - A*AINV ) / ( N * norm(A) * norm(AINV) * EPS )
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
!>          complex symmetric matrix A is stored:
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
!>          A is COMPLEX array, dimension (LDA,N)
!>          The original complex symmetric matrix A.
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
!>          AINV is COMPLEX array, dimension (LDAINV,N)
!>          On entry, the inverse of the matrix A, stored as a symmetric
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
!>          WORK is COMPLEX array, dimension (LDWORK,N)
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
!>          RCOND = 1/ (norm(A) * norm(AINV)).
!> \endverbatim
!>
!> \param[out] RESID
!> \verbatim
!>          RESID is REAL
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
!> \ingroup complex_lin
!
!  =====================================================================
      SUBROUTINE CSYT03(Uplo,N,A,Lda,Ainv,Ldainv,Work,Ldwork,Rwork,     &
     &                  Rcond,Resid)
      IMPLICIT NONE
!*--CSYT03130
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Uplo
      INTEGER Lda , Ldainv , Ldwork , N
      REAL Rcond , Resid
!     ..
!     .. Array Arguments ..
      REAL Rwork(*)
      COMPLEX A(Lda,*) , Ainv(Ldainv,*) , Work(Ldwork,*)
!     ..
!
!  =====================================================================
!
!
!     .. Parameters ..
      REAL ZERO , ONE
      PARAMETER (ZERO=0.0E+0,ONE=1.0E+0)
      COMPLEX CZERO , CONE
      PARAMETER (CZERO=(0.0E+0,0.0E+0),CONE=(1.0E+0,0.0E+0))
!     ..
!     .. Local Scalars ..
      INTEGER i , j
      REAL ainvnm , anorm , eps
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      REAL CLANGE , CLANSY , SLAMCH
      EXTERNAL LSAME , CLANGE , CLANSY , SLAMCH
!     ..
!     .. External Subroutines ..
      EXTERNAL CSYMM
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
      anorm = CLANSY('1',Uplo,N,A,Lda,Rwork)
      ainvnm = CLANSY('1',Uplo,N,Ainv,Ldainv,Rwork)
      IF ( anorm<=ZERO .OR. ainvnm<=ZERO ) THEN
         Rcond = ZERO
         Resid = ONE/eps
         RETURN
      ENDIF
      Rcond = (ONE/anorm)/ainvnm
!
!     Expand AINV into a full matrix and call CSYMM to multiply
!     AINV on the left by A (store the result in WORK).
!
      IF ( LSAME(Uplo,'U') ) THEN
         DO j = 1 , N
            DO i = 1 , j - 1
               Ainv(j,i) = Ainv(i,j)
            ENDDO
         ENDDO
      ELSE
         DO j = 1 , N
            DO i = j + 1 , N
               Ainv(j,i) = Ainv(i,j)
            ENDDO
         ENDDO
      ENDIF
      CALL CSYMM('Left',Uplo,N,N,-CONE,A,Lda,Ainv,Ldainv,CZERO,Work,    &
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
      Resid = CLANGE('1',N,N,Work,Ldwork,Rwork)
!
      Resid = ((Resid*Rcond)/eps)/REAL(N)
!
!
!     End of CSYT03
!
      END SUBROUTINE CSYT03
