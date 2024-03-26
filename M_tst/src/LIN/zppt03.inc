!*==zppt03.f90  processed by SPAG 7.51RB at 17:35 on  4 Mar 2022
!> \brief \b zppt03
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE zppt03( UPLO, N, A, AINV, WORK, LDWORK, RWORK, RCOND,
!                          RESID )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            LDWORK, N
!       DOUBLE PRECISION   RCOND, RESID
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   RWORK( * )
!       COMPLEX*16         A( * ), AINV( * ), WORK( LDWORK, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> zppt03 computes the residual for a Hermitian packed matrix times its
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
!>          A is COMPLEX*16 array, dimension (N*(N+1)/2)
!>          The original Hermitian matrix A, stored as a packed
!>          triangular matrix.
!> \endverbatim
!>
!> \param[in] AINV
!> \verbatim
!>          AINV is COMPLEX*16 array, dimension (N*(N+1)/2)
!>          The (Hermitian) inverse of the matrix A, stored as a packed
!>          triangular matrix.
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
      SUBROUTINE zppt03(Uplo,N,A,Ainv,Work,Ldwork,Rwork,Rcond,Resid)
      IMPLICIT NONE
!*--ZPPT03113
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Uplo
      INTEGER Ldwork , N
      DOUBLE PRECISION Rcond , Resid
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION Rwork(*)
      COMPLEX*16 A(*) , Ainv(*) , Work(Ldwork,*)
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
      INTEGER i , j , jj
      DOUBLE PRECISION ainvnm , anorm , eps
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      DOUBLE PRECISION DLAMCH , ZLANGE , ZLANHP
EXTERNAL    LSAME,   DLAMCH,   ZLANGE,   ZLANHP
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC DBLE , DCONJG
!     ..
!     .. External Subroutines ..
EXTERNAL    ZCOPY,   ZHPMV
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
      anorm = ZLANHP('1',Uplo,N,A,Rwork)
      ainvnm = ZLANHP('1',Uplo,N,Ainv,Rwork)
      IF ( anorm<=ZERO .OR. ainvnm<=ZERO ) THEN
         Rcond = ZERO
         Resid = ONE/eps
         RETURN
      ENDIF
      Rcond = (ONE/anorm)/ainvnm
!
!     UPLO = 'U':
!     Copy the leading N-1 x N-1 submatrix of AINV to WORK(1:N,2:N) and
!     expand it to a full matrix, then multiply by A one column at a
!     time, moving the result one column to the left.
!
      IF ( LSAME(Uplo,'U') ) THEN
!
!        Copy AINV
!
         jj = 1
         DO j = 1 , N - 1
            CALL ZCOPY(j,Ainv(jj),1,Work(1,j+1),1)
            DO i = 1 , j - 1
               Work(j,i+1) = DCONJG(Ainv(jj+i-1))
            ENDDO
            jj = jj + j
         ENDDO
         jj = ((N-1)*N)/2 + 1
         DO i = 1 , N - 1
            Work(N,i+1) = DCONJG(Ainv(jj+i-1))
         ENDDO
!
!        Multiply by A
!
         DO j = 1 , N - 1
            CALL ZHPMV('Upper',N,-CONE,A,Work(1,j+1),1,CZERO,Work(1,j), &
     &                 1)
         ENDDO
         CALL ZHPMV('Upper',N,-CONE,A,Ainv(jj),1,CZERO,Work(1,N),1)
!
!     UPLO = 'L':
!     Copy the trailing N-1 x N-1 submatrix of AINV to WORK(1:N,1:N-1)
!     and multiply by A, moving each column to the right.
!
      ELSE
!
!        Copy AINV
!
         DO i = 1 , N - 1
            Work(1,i) = DCONJG(Ainv(i+1))
         ENDDO
         jj = N + 1
         DO j = 2 , N
            CALL ZCOPY(N-j+1,Ainv(jj),1,Work(j,j-1),1)
            DO i = 1 , N - j
               Work(j,j+i-1) = DCONJG(Ainv(jj+i))
            ENDDO
            jj = jj + N - j + 1
         ENDDO
!
!        Multiply by A
!
         DO j = N , 2 , -1
            CALL ZHPMV('Lower',N,-CONE,A,Work(1,j-1),1,CZERO,Work(1,j), &
     &                 1)
         ENDDO
         CALL ZHPMV('Lower',N,-CONE,A,Ainv(1),1,CZERO,Work(1,1),1)
!
      ENDIF
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
!     End of zppt03
!
      END SUBROUTINE zppt03
