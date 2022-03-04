!*==dpbcon.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b DPBCON
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DPBCON + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dpbcon.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dpbcon.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dpbcon.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DPBCON( UPLO, N, KD, AB, LDAB, ANORM, RCOND, WORK,
!                          IWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            INFO, KD, LDAB, N
!       DOUBLE PRECISION   ANORM, RCOND
!       ..
!       .. Array Arguments ..
!       INTEGER            IWORK( * )
!       DOUBLE PRECISION   AB( LDAB, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DPBCON estimates the reciprocal of the condition number (in the
!> 1-norm) of a real symmetric positive definite band matrix using the
!> Cholesky factorization A = U**T*U or A = L*L**T computed by DPBTRF.
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
!>          or the number of subdiagonals if UPLO = 'L'.  KD >= 0.
!> \endverbatim
!>
!> \param[in] AB
!> \verbatim
!>          AB is DOUBLE PRECISION array, dimension (LDAB,N)
!>          The triangular factor U or L from the Cholesky factorization
!>          A = U**T*U or A = L*L**T of the band matrix A, stored in the
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
!>          The 1-norm (or infinity-norm) of the symmetric band matrix A.
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
!>          WORK is DOUBLE PRECISION array, dimension (3*N)
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
!> \ingroup doubleOTHERcomputational
!
!  =====================================================================
      SUBROUTINE DPBCON(Uplo,N,Kd,Ab,Ldab,Anorm,Rcond,Work,Iwork,Info)
      IMPLICIT NONE
!*--DPBCON135
!
!  -- LAPACK computational routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Uplo
      INTEGER Info , Kd , Ldab , N
      DOUBLE PRECISION Anorm , Rcond
!     ..
!     .. Array Arguments ..
      INTEGER Iwork(*)
      DOUBLE PRECISION Ab(Ldab,*) , Work(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ONE , ZERO
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
!     ..
!     .. Local Scalars ..
      LOGICAL upper
      CHARACTER normin
      INTEGER ix , kase
      DOUBLE PRECISION ainvnm , scale , scalel , scaleu , smlnum
!     ..
!     .. Local Arrays ..
      INTEGER isave(3)
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      INTEGER IDAMAX
      DOUBLE PRECISION DLAMCH
      EXTERNAL LSAME , IDAMAX , DLAMCH
!     ..
!     .. External Subroutines ..
      EXTERNAL DLACN2 , DLATBS , DRSCL , XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS
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
         CALL XERBLA('DPBCON',-Info)
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
         CALL DLACN2(N,Work(N+1),Work,Iwork,ainvnm,kase,isave)
         IF ( kase/=0 ) THEN
            IF ( upper ) THEN
!
!           Multiply by inv(U**T).
!
               CALL DLATBS('Upper','Transpose','Non-unit',normin,N,Kd,  &
     &                     Ab,Ldab,Work,scalel,Work(2*N+1),Info)
               normin = 'Y'
!
!           Multiply by inv(U).
!
               CALL DLATBS('Upper','No transpose','Non-unit',normin,N,  &
     &                     Kd,Ab,Ldab,Work,scaleu,Work(2*N+1),Info)
            ELSE
!
!           Multiply by inv(L).
!
               CALL DLATBS('Lower','No transpose','Non-unit',normin,N,  &
     &                     Kd,Ab,Ldab,Work,scalel,Work(2*N+1),Info)
               normin = 'Y'
!
!           Multiply by inv(L**T).
!
               CALL DLATBS('Lower','Transpose','Non-unit',normin,N,Kd,  &
     &                     Ab,Ldab,Work,scaleu,Work(2*N+1),Info)
            ENDIF
!
!        Multiply by 1/SCALE if doing so will not cause overflow.
!
            scale = scalel*scaleu
            IF ( scale/=ONE ) THEN
               ix = IDAMAX(N,Work,1)
               IF ( scale<ABS(Work(ix))*smlnum .OR. scale==ZERO ) EXIT
               CALL DRSCL(N,scale,Work,1)
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
!     End of DPBCON
!
      END SUBROUTINE DPBCON