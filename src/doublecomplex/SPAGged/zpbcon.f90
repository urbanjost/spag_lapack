!*==zpbcon.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b ZPBCON
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZPBCON + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zpbcon.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zpbcon.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zpbcon.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZPBCON( UPLO, N, KD, AB, LDAB, ANORM, RCOND, WORK,
!                          RWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            INFO, KD, LDAB, N
!       DOUBLE PRECISION   ANORM, RCOND
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   RWORK( * )
!       COMPLEX*16         AB( LDAB, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZPBCON estimates the reciprocal of the condition number (in the
!> 1-norm) of a complex Hermitian positive definite band matrix using
!> the Cholesky factorization A = U**H*U or A = L*L**H computed by
!> ZPBTRF.
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
!>          = 'U':  Upper triangular factor stored in AB;
!>          = 'L':  Lower triangular factor stored in AB.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in] KD
!> \verbatim
!>          KD is INTEGER
!>          The number of superdiagonals of the matrix A if UPLO = 'U',
!>          or the number of sub-diagonals if UPLO = 'L'.  KD >= 0.
!> \endverbatim
!>
!> \param[in] AB
!> \verbatim
!>          AB is COMPLEX*16 array, dimension (LDAB,N)
!>          The triangular factor U or L from the Cholesky factorization
!>          A = U**H*U or A = L*L**H of the band matrix A, stored in the
!>          first KD+1 rows of the array.  The j-th column of U or L is
!>          stored in the j-th column of the array AB as follows:
!>          if UPLO ='U', AB(kd+1+i-j,j) = U(i,j) for max(1,j-kd)<=i<=j;
!>          if UPLO ='L', AB(1+i-j,j)    = L(i,j) for j<=i<=min(n,j+kd).
!> \endverbatim
!>
!> \param[in] LDAB
!> \verbatim
!>          LDAB is INTEGER
!>          The leading dimension of the array AB.  LDAB >= KD+1.
!> \endverbatim
!>
!> \param[in] ANORM
!> \verbatim
!>          ANORM is DOUBLE PRECISION
!>          The 1-norm (or infinity-norm) of the Hermitian band matrix A.
!> \endverbatim
!>
!> \param[out] RCOND
!> \verbatim
!>          RCOND is DOUBLE PRECISION
!>          The reciprocal of the condition number of the matrix A,
!>          computed as RCOND = 1/(ANORM * AINVNM), where AINVNM is an
!>          estimate of the 1-norm of inv(A) computed in this routine.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX*16 array, dimension (2*N)
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is DOUBLE PRECISION array, dimension (N)
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
!> \ingroup complex16OTHERcomputational
!
!  =====================================================================
      SUBROUTINE ZPBCON(Uplo,N,Kd,Ab,Ldab,Anorm,Rcond,Work,Rwork,Info)
      USE F77KINDS                        
      USE S_DLAMCH
      USE S_IZAMAX
      USE S_LSAME
      USE S_XERBLA
      USE S_ZDRSCL
      USE S_ZLACN2
      USE S_ZLATBS
      IMPLICIT NONE
!*--ZPBCON144
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Kd
      COMPLEX(CX16KIND) , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      REAL(R8KIND) , INTENT(IN) :: Anorm
      REAL(R8KIND) , INTENT(OUT) :: Rcond
      COMPLEX(CX16KIND) , DIMENSION(*) :: Work
      REAL(R8KIND) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      REAL(R8KIND) :: ainvnm , scale , scalel , scaleu , smlnum
      REAL(R8KIND) :: CABS1
      INTEGER , DIMENSION(3) :: isave
      INTEGER :: ix , kase
      CHARACTER :: normin
      LOGICAL :: upper
      COMPLEX(CX16KIND) :: zdum
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
!     .. Statement Functions ..
!     ..
!     .. Statement Function definitions ..
      CABS1(zdum) = ABS(DBLE(zdum)) + ABS(DIMAG(zdum))
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
      ELSEIF ( Kd<0 ) THEN
         Info = -3
      ELSEIF ( Ldab<Kd+1 ) THEN
         Info = -5
      ELSEIF ( Anorm<ZERO ) THEN
         Info = -6
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('ZPBCON',-Info)
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
      smlnum = DLAMCH('Safe minimum')
!
!     Estimate the 1-norm of the inverse.
!
      kase = 0
      normin = 'N'
      DO
         CALL ZLACN2(N,Work(N+1),Work,ainvnm,kase,isave)
         IF ( kase/=0 ) THEN
            IF ( upper ) THEN
!
!           Multiply by inv(U**H).
!
               CALL ZLATBS('Upper','Conjugate transpose','Non-unit',    &
     &                     normin,N,Kd,Ab,Ldab,Work,scalel,Rwork,Info)
               normin = 'Y'
!
!           Multiply by inv(U).
!
               CALL ZLATBS('Upper','No transpose','Non-unit',normin,N,  &
     &                     Kd,Ab,Ldab,Work,scaleu,Rwork,Info)
            ELSE
!
!           Multiply by inv(L).
!
               CALL ZLATBS('Lower','No transpose','Non-unit',normin,N,  &
     &                     Kd,Ab,Ldab,Work,scalel,Rwork,Info)
               normin = 'Y'
!
!           Multiply by inv(L**H).
!
               CALL ZLATBS('Lower','Conjugate transpose','Non-unit',    &
     &                     normin,N,Kd,Ab,Ldab,Work,scaleu,Rwork,Info)
            ENDIF
!
!        Multiply by 1/SCALE if doing so will not cause overflow.
!
            scale = scalel*scaleu
            IF ( scale/=ONE ) THEN
               ix = IZAMAX(N,Work,1)
               IF ( scale<CABS1(Work(ix))*smlnum .OR. scale==ZERO ) EXIT
               CALL ZDRSCL(N,scale,Work,1)
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
!
!     End of ZPBCON
!
      END SUBROUTINE ZPBCON
