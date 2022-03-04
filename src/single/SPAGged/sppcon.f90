!*==sppcon.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b SPPCON
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SPPCON + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sppcon.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sppcon.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sppcon.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SPPCON( UPLO, N, AP, ANORM, RCOND, WORK, IWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            INFO, N
!       REAL               ANORM, RCOND
!       ..
!       .. Array Arguments ..
!       INTEGER            IWORK( * )
!       REAL               AP( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SPPCON estimates the reciprocal of the condition number (in the
!> 1-norm) of a real symmetric positive definite packed matrix using
!> the Cholesky factorization A = U**T*U or A = L*L**T computed by
!> SPPTRF.
!>
!> An estimate is obtained for norm(inv(A)), and the reciprocal of the
!> condition number is computed as RCOND = 1 / (ANORM * norm(inv(A))).
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          = 'U':  Upper triangle of A is stored;
!>          = 'L':  Lower triangle of A is stored.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in] AP
!> \verbatim
!>          AP is REAL array, dimension (N*(N+1)/2)
!>          The triangular factor U or L from the Cholesky factorization
!>          A = U**T*U or A = L*L**T, packed columnwise in a linear
!>          array.  The j-th column of U or L is stored in the array AP
!>          as follows:
!>          if UPLO = 'U', AP(i + (j-1)*j/2) = U(i,j) for 1<=i<=j;
!>          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = L(i,j) for j<=i<=n.
!> \endverbatim
!>
!> \param[in] ANORM
!> \verbatim
!>          ANORM is REAL
!>          The 1-norm (or infinity-norm) of the symmetric matrix A.
!> \endverbatim
!>
!> \param[out] RCOND
!> \verbatim
!>          RCOND is REAL
!>          The reciprocal of the condition number of the matrix A,
!>          computed as RCOND = 1/(ANORM * AINVNM), where AINVNM is an
!>          estimate of the 1-norm of inv(A) computed in this routine.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is REAL array, dimension (3*N)
!> \endverbatim
!>
!> \param[out] IWORK
!> \verbatim
!>          IWORK is INTEGER array, dimension (N)
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value
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
!> \ingroup realOTHERcomputational
!
!  =====================================================================
      SUBROUTINE SPPCON(Uplo,N,Ap,Anorm,Rcond,Work,Iwork,Info)
      USE S_ISAMAX
      USE S_LSAME
      USE S_SLACN2
      USE S_SLAMCH
      USE S_SLATPS
      USE S_SRSCL
      USE S_XERBLA
      IMPLICIT NONE
!*--SPPCON129
!
! PARAMETER definitions rewritten by SPAG
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Uplo
      INTEGER :: N
      REAL , DIMENSION(*) :: Ap
      REAL , INTENT(IN) :: Anorm
      REAL , INTENT(OUT) :: Rcond
      REAL , DIMENSION(*) :: Work
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      REAL :: ainvnm , scale , scalel , scaleu , smlnum
      INTEGER , DIMENSION(3) :: isave
      INTEGER :: ix , kase
      CHARACTER :: normin
      LOGICAL :: upper
!
! End of declarations rewritten by SPAG
!
!     ..
!     .. Array Arguments ..
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
!     ..
!     .. Local Scalars ..
!     ..
!     .. Local Arrays ..
!     ..
!     .. External Functions ..
!     ..
!     .. External Subroutines ..
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      Info = 0
      upper = LSAME(Uplo,'U')
      IF ( .NOT.upper .AND. .NOT.LSAME(Uplo,'L') ) THEN
         Info = -1
      ELSEIF ( N<0 ) THEN
         Info = -2
      ELSEIF ( Anorm<ZERO ) THEN
         Info = -4
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('SPPCON',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      Rcond = ZERO
      IF ( N==0 ) THEN
         Rcond = ONE
         RETURN
      ELSEIF ( Anorm==ZERO ) THEN
         RETURN
      ENDIF
!
      smlnum = SLAMCH('Safe minimum')
!
!     Estimate the 1-norm of the inverse.
!
      kase = 0
      normin = 'N'
      DO
         CALL SLACN2(N,Work(N+1),Work,Iwork,ainvnm,kase,isave)
         IF ( kase/=0 ) THEN
            IF ( upper ) THEN
!
!           Multiply by inv(U**T).
!
               CALL SLATPS('Upper','Transpose','Non-unit',normin,N,Ap,  &
     &                     Work,scalel,Work(2*N+1),Info)
               normin = 'Y'
!
!           Multiply by inv(U).
!
               CALL SLATPS('Upper','No transpose','Non-unit',normin,N,  &
     &                     Ap,Work,scaleu,Work(2*N+1),Info)
            ELSE
!
!           Multiply by inv(L).
!
               CALL SLATPS('Lower','No transpose','Non-unit',normin,N,  &
     &                     Ap,Work,scalel,Work(2*N+1),Info)
               normin = 'Y'
!
!           Multiply by inv(L**T).
!
               CALL SLATPS('Lower','Transpose','Non-unit',normin,N,Ap,  &
     &                     Work,scaleu,Work(2*N+1),Info)
            ENDIF
!
!        Multiply by 1/SCALE if doing so will not cause overflow.
!
            scale = scalel*scaleu
            IF ( scale/=ONE ) THEN
               ix = ISAMAX(N,Work,1)
               IF ( scale<ABS(Work(ix))*smlnum .OR. scale==ZERO ) EXIT
               CALL SRSCL(N,scale,Work,1)
            ENDIF
            CYCLE
         ENDIF
!
!     Compute the estimate of the reciprocal condition number.
!
         IF ( ainvnm/=ZERO ) Rcond = (ONE/ainvnm)/Anorm
         EXIT
      ENDDO
!
!
!     End of SPPCON
!
      END SUBROUTINE SPPCON
