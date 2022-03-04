!*==dgbcon.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b DGBCON
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DGBCON + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgbcon.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgbcon.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgbcon.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DGBCON( NORM, N, KL, KU, AB, LDAB, IPIV, ANORM, RCOND,
!                          WORK, IWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          NORM
!       INTEGER            INFO, KL, KU, LDAB, N
!       DOUBLE PRECISION   ANORM, RCOND
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * ), IWORK( * )
!       DOUBLE PRECISION   AB( LDAB, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DGBCON estimates the reciprocal of the condition number of a real
!> general band matrix A, in either the 1-norm or the infinity-norm,
!> using the LU factorization computed by DGBTRF.
!>
!> An estimate is obtained for norm(inv(A)), and the reciprocal of the
!> condition number is computed as
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
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in] KL
!> \verbatim
!>          KL is INTEGER
!>          The number of subdiagonals within the band of A.  KL >= 0.
!> \endverbatim
!>
!> \param[in] KU
!> \verbatim
!>          KU is INTEGER
!>          The number of superdiagonals within the band of A.  KU >= 0.
!> \endverbatim
!>
!> \param[in] AB
!> \verbatim
!>          AB is DOUBLE PRECISION array, dimension (LDAB,N)
!>          Details of the LU factorization of the band matrix A, as
!>          computed by DGBTRF.  U is stored as an upper triangular band
!>          matrix with KL+KU superdiagonals in rows 1 to KL+KU+1, and
!>          the multipliers used during the factorization are stored in
!>          rows KL+KU+2 to 2*KL+KU+1.
!> \endverbatim
!>
!> \param[in] LDAB
!> \verbatim
!>          LDAB is INTEGER
!>          The leading dimension of the array AB.  LDAB >= 2*KL+KU+1.
!> \endverbatim
!>
!> \param[in] IPIV
!> \verbatim
!>          IPIV is INTEGER array, dimension (N)
!>          The pivot indices; for 1 <= i <= N, row i of the matrix was
!>          interchanged with row IPIV(i).
!> \endverbatim
!>
!> \param[in] ANORM
!> \verbatim
!>          ANORM is DOUBLE PRECISION
!>          If NORM = '1' or 'O', the 1-norm of the original matrix A.
!>          If NORM = 'I', the infinity-norm of the original matrix A.
!> \endverbatim
!>
!> \param[out] RCOND
!> \verbatim
!>          RCOND is DOUBLE PRECISION
!>          The reciprocal of the condition number of the matrix A,
!>          computed as RCOND = 1/(norm(A) * norm(inv(A))).
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
!>          < 0: if INFO = -i, the i-th argument had an illegal value
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
!> \ingroup doubleGBcomputational
!
!  =====================================================================
      SUBROUTINE DGBCON(Norm,N,Kl,Ku,Ab,Ldab,Ipiv,Anorm,Rcond,Work,     &
     &                  Iwork,Info)
      USE F77KINDS                        
      USE S_DAXPY
      USE S_DDOT
      USE S_DLACN2
      USE S_DLAMCH
      USE S_DLATBS
      USE S_DRSCL
      USE S_IDAMAX
      USE S_LSAME
      USE S_XERBLA
      IMPLICIT NONE
!*--DGBCON160
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Norm
      INTEGER :: N
      INTEGER , INTENT(IN) :: Kl
      INTEGER , INTENT(IN) :: Ku
      REAL(R8KIND) , DIMENSION(Ldab,*) :: Ab
      INTEGER :: Ldab
      INTEGER , INTENT(IN) , DIMENSION(*) :: Ipiv
      REAL(R8KIND) , INTENT(IN) :: Anorm
      REAL(R8KIND) , INTENT(OUT) :: Rcond
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      REAL(R8KIND) :: ainvnm , scale , smlnum , t
      INTEGER , DIMENSION(3) :: isave
      INTEGER :: ix , j , jp , kase , kase1 , kd , lm
      LOGICAL :: lnoti , onenrm
      CHARACTER :: normin
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
      onenrm = Norm=='1' .OR. LSAME(Norm,'O')
      IF ( .NOT.onenrm .AND. .NOT.LSAME(Norm,'I') ) THEN
         Info = -1
      ELSEIF ( N<0 ) THEN
         Info = -2
      ELSEIF ( Kl<0 ) THEN
         Info = -3
      ELSEIF ( Ku<0 ) THEN
         Info = -4
      ELSEIF ( Ldab<2*Kl+Ku+1 ) THEN
         Info = -6
      ELSEIF ( Anorm<ZERO ) THEN
         Info = -8
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('DGBCON',-Info)
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
!     Estimate the norm of inv(A).
!
      ainvnm = ZERO
      normin = 'N'
      IF ( onenrm ) THEN
         kase1 = 1
      ELSE
         kase1 = 2
      ENDIF
      kd = Kl + Ku + 1
      lnoti = Kl>0
      kase = 0
      DO
         CALL DLACN2(N,Work(N+1),Work,Iwork,ainvnm,kase,isave)
         IF ( kase/=0 ) THEN
            IF ( kase==kase1 ) THEN
!
!           Multiply by inv(L).
!
               IF ( lnoti ) THEN
                  DO j = 1 , N - 1
                     lm = MIN(Kl,N-j)
                     jp = Ipiv(j)
                     t = Work(jp)
                     IF ( jp/=j ) THEN
                        Work(jp) = Work(j)
                        Work(j) = t
                     ENDIF
                     CALL DAXPY(lm,-t,Ab(kd+1,j),1,Work(j+1),1)
                  ENDDO
               ENDIF
!
!           Multiply by inv(U).
!
               CALL DLATBS('Upper','No transpose','Non-unit',normin,N,  &
     &                     Kl+Ku,Ab,Ldab,Work,scale,Work(2*N+1),Info)
            ELSE
!
!           Multiply by inv(U**T).
!
               CALL DLATBS('Upper','Transpose','Non-unit',normin,N,     &
     &                     Kl+Ku,Ab,Ldab,Work,scale,Work(2*N+1),Info)
!
!           Multiply by inv(L**T).
!
               IF ( lnoti ) THEN
                  DO j = N - 1 , 1 , -1
                     lm = MIN(Kl,N-j)
                     Work(j) = Work(j)                                  &
     &                         - DDOT(lm,Ab(kd+1,j),1,Work(j+1),1)
                     jp = Ipiv(j)
                     IF ( jp/=j ) THEN
                        t = Work(jp)
                        Work(jp) = Work(j)
                        Work(j) = t
                     ENDIF
                  ENDDO
               ENDIF
            ENDIF
!
!        Divide X by 1/SCALE if doing so will not cause overflow.
!
            normin = 'Y'
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
!     End of DGBCON
!
      END SUBROUTINE DGBCON
