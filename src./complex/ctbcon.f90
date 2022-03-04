!*==ctbcon.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b CTBCON
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CTBCON + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ctbcon.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ctbcon.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ctbcon.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CTBCON( NORM, UPLO, DIAG, N, KD, AB, LDAB, RCOND, WORK,
!                          RWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          DIAG, NORM, UPLO
!       INTEGER            INFO, KD, LDAB, N
!       REAL               RCOND
!       ..
!       .. Array Arguments ..
!       REAL               RWORK( * )
!       COMPLEX            AB( LDAB, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CTBCON estimates the reciprocal of the condition number of a
!> triangular band matrix A, in either the 1-norm or the infinity-norm.
!>
!> The norm of A is computed and an estimate is obtained for
!> norm(inv(A)), then the reciprocal of the condition number is
!> computed as
!>    RCOND = 1 / ( norm(A) * norm(inv(A)) ).
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] NORM
!> \verbatim
!>          NORM is CHARACTER*1
!>          Specifies whether the 1-norm condition number or the
!>          infinity-norm condition number is required:
!>          = '1' or 'O':  1-norm;
!>          = 'I':         Infinity-norm.
!> \endverbatim
!>
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          = 'U':  A is upper triangular;
!>          = 'L':  A is lower triangular.
!> \endverbatim
!>
!> \param[in] DIAG
!> \verbatim
!>          DIAG is CHARACTER*1
!>          = 'N':  A is non-unit triangular;
!>          = 'U':  A is unit triangular.
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
!>          The number of superdiagonals or subdiagonals of the
!>          triangular band matrix A.  KD >= 0.
!> \endverbatim
!>
!> \param[in] AB
!> \verbatim
!>          AB is COMPLEX array, dimension (LDAB,N)
!>          The upper or lower triangular band matrix A, stored in the
!>          first kd+1 rows of the array. The j-th column of A is stored
!>          in the j-th column of the array AB as follows:
!>          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j;
!>          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd).
!>          If DIAG = 'U', the diagonal elements of A are not referenced
!>          and are assumed to be 1.
!> \endverbatim
!>
!> \param[in] LDAB
!> \verbatim
!>          LDAB is INTEGER
!>          The leading dimension of the array AB.  LDAB >= KD+1.
!> \endverbatim
!>
!> \param[out] RCOND
!> \verbatim
!>          RCOND is REAL
!>          The reciprocal of the condition number of the matrix A,
!>          computed as RCOND = 1/(norm(A) * norm(inv(A))).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX array, dimension (2*N)
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is REAL array, dimension (N)
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
!> \ingroup complexOTHERcomputational
!
!  =====================================================================
      SUBROUTINE CTBCON(Norm,Uplo,Diag,N,Kd,Ab,Ldab,Rcond,Work,Rwork,   &
     &                  Info)
      USE S_CLACN2
      USE S_CLANTB
      USE S_CLATBS
      USE S_CSRSCL
      USE S_ICAMAX
      USE S_LSAME
      USE S_SLAMCH
      USE S_XERBLA
      IMPLICIT NONE
!*--CTBCON155
!
! PARAMETER definitions rewritten by SPAG
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , ZERO = 0.0E+0
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Norm
      CHARACTER :: Uplo
      CHARACTER :: Diag
      INTEGER :: N
      INTEGER :: Kd
      COMPLEX , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      REAL , INTENT(OUT) :: Rcond
      COMPLEX , DIMENSION(*) :: Work
      REAL , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      REAL :: ainvnm , anorm , scale , smlnum , xnorm
      REAL :: CABS1
      INTEGER , DIMENSION(3) :: isave
      INTEGER :: ix , kase , kase1
      CHARACTER :: normin
      LOGICAL :: nounit , onenrm , upper
      COMPLEX :: zdum
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
      CABS1(zdum) = ABS(REAL(zdum)) + ABS(AIMAG(zdum))
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      Info = 0
      upper = LSAME(Uplo,'U')
      onenrm = Norm=='1' .OR. LSAME(Norm,'O')
      nounit = LSAME(Diag,'N')
!
      IF ( .NOT.onenrm .AND. .NOT.LSAME(Norm,'I') ) THEN
         Info = -1
      ELSEIF ( .NOT.upper .AND. .NOT.LSAME(Uplo,'L') ) THEN
         Info = -2
      ELSEIF ( .NOT.nounit .AND. .NOT.LSAME(Diag,'U') ) THEN
         Info = -3
      ELSEIF ( N<0 ) THEN
         Info = -4
      ELSEIF ( Kd<0 ) THEN
         Info = -5
      ELSEIF ( Ldab<Kd+1 ) THEN
         Info = -7
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('CTBCON',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( N==0 ) THEN
         Rcond = ONE
         RETURN
      ENDIF
!
      Rcond = ZERO
      smlnum = SLAMCH('Safe minimum')*REAL(MAX(N,1))
!
!     Compute the 1-norm of the triangular matrix A or A**H.
!
      anorm = CLANTB(Norm,Uplo,Diag,N,Kd,Ab,Ldab,Rwork)
!
!     Continue only if ANORM > 0.
!
      IF ( anorm>ZERO ) THEN
!
!        Estimate the 1-norm of the inverse of A.
!
         ainvnm = ZERO
         normin = 'N'
         IF ( onenrm ) THEN
            kase1 = 1
         ELSE
            kase1 = 2
         ENDIF
         kase = 0
         DO
            CALL CLACN2(N,Work(N+1),Work,ainvnm,kase,isave)
            IF ( kase/=0 ) THEN
               IF ( kase==kase1 ) THEN
!
!              Multiply by inv(A).
!
                  CALL CLATBS(Uplo,'No transpose',Diag,normin,N,Kd,Ab,  &
     &                        Ldab,Work,scale,Rwork,Info)
               ELSE
!
!              Multiply by inv(A**H).
!
                  CALL CLATBS(Uplo,'Conjugate transpose',Diag,normin,N, &
     &                        Kd,Ab,Ldab,Work,scale,Rwork,Info)
               ENDIF
               normin = 'Y'
!
!           Multiply by 1/SCALE if doing so will not cause overflow.
!
               IF ( scale/=ONE ) THEN
                  ix = ICAMAX(N,Work,1)
                  xnorm = CABS1(Work(ix))
                  IF ( scale<xnorm*smlnum .OR. scale==ZERO ) EXIT
                  CALL CSRSCL(N,scale,Work,1)
               ENDIF
               CYCLE
            ENDIF
!
!        Compute the estimate of the reciprocal condition number.
!
            IF ( ainvnm/=ZERO ) Rcond = (ONE/anorm)/ainvnm
            EXIT
         ENDDO
      ENDIF
!
!
!     End of CTBCON
!
      END SUBROUTINE CTBCON
