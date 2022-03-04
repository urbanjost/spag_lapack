!*==sla_porcond.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b SLA_PORCOND estimates the Skeel condition number for a symmetric positive-definite matrix.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SLA_PORCOND + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sla_porcond.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sla_porcond.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sla_porcond.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       REAL FUNCTION SLA_PORCOND( UPLO, N, A, LDA, AF, LDAF, CMODE, C,
!                                  INFO, WORK, IWORK )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            N, LDA, LDAF, INFO, CMODE
!       REAL               A( LDA, * ), AF( LDAF, * ), WORK( * ),
!      $                   C( * )
!       ..
!       .. Array Arguments ..
!       INTEGER            IWORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    SLA_PORCOND Estimates the Skeel condition number of  op(A) * op2(C)
!>    where op2 is determined by CMODE as follows
!>    CMODE =  1    op2(C) = C
!>    CMODE =  0    op2(C) = I
!>    CMODE = -1    op2(C) = inv(C)
!>    The Skeel condition number  cond(A) = norminf( |inv(A)||A| )
!>    is computed by computing scaling factors R such that
!>    diag(R)*A*op2(C) is row equilibrated and computing the standard
!>    infinity-norm condition number.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>       = 'U':  Upper triangle of A is stored;
!>       = 'L':  Lower triangle of A is stored.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>     The number of linear equations, i.e., the order of the
!>     matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is REAL array, dimension (LDA,N)
!>     On entry, the N-by-N matrix A.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>     The leading dimension of the array A.  LDA >= max(1,N).
!> \endverbatim
!>
!> \param[in] AF
!> \verbatim
!>          AF is REAL array, dimension (LDAF,N)
!>     The triangular factor U or L from the Cholesky factorization
!>     A = U**T*U or A = L*L**T, as computed by SPOTRF.
!> \endverbatim
!>
!> \param[in] LDAF
!> \verbatim
!>          LDAF is INTEGER
!>     The leading dimension of the array AF.  LDAF >= max(1,N).
!> \endverbatim
!>
!> \param[in] CMODE
!> \verbatim
!>          CMODE is INTEGER
!>     Determines op2(C) in the formula op(A) * op2(C) as follows:
!>     CMODE =  1    op2(C) = C
!>     CMODE =  0    op2(C) = I
!>     CMODE = -1    op2(C) = inv(C)
!> \endverbatim
!>
!> \param[in] C
!> \verbatim
!>          C is REAL array, dimension (N)
!>     The vector C in the formula op(A) * op2(C).
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>       = 0:  Successful exit.
!>     i > 0:  The ith argument is invalid.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is REAL array, dimension (3*N).
!>     Workspace.
!> \endverbatim
!>
!> \param[out] IWORK
!> \verbatim
!>          IWORK is INTEGER array, dimension (N).
!>     Workspace.
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
!> \ingroup realPOcomputational
!
!  =====================================================================
      FUNCTION SLA_PORCOND(Uplo,N,A,Lda,Af,Ldaf,Cmode,C,Info,Work,Iwork)
      USE S_LSAME
      USE S_SLACN2
      USE S_SPOTRS
      USE S_XERBLA
      IMPLICIT NONE
!*--SLA_PORCOND147
      REAL :: SLA_PORCOND
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Uplo
      INTEGER :: N
      REAL , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL , DIMENSION(Ldaf,*) :: Af
      INTEGER :: Ldaf
      INTEGER , INTENT(IN) :: Cmode
      REAL , INTENT(IN) , DIMENSION(*) :: C
      INTEGER , INTENT(INOUT) :: Info
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER , DIMENSION(*) :: Iwork
!
! Local variable declarations rewritten by SPAG
!
      REAL :: ainvnm , tmp
      INTEGER :: i , j , kase
      INTEGER , DIMENSION(3) :: isave
      LOGICAL :: up
!
! End of declarations rewritten by SPAG
!
!     ..
!     .. Array Arguments ..
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
!     ..
!     .. Array Arguments ..
!     ..
!     .. External Functions ..
!     ..
!     .. External Subroutines ..
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. Executable Statements ..
!
      SLA_PORCOND = 0.0
!
      Info = 0
      IF ( N<0 ) Info = -2
      IF ( Info/=0 ) THEN
         CALL XERBLA('SLA_PORCOND',-Info)
         RETURN
      ENDIF
 
      IF ( N==0 ) THEN
         SLA_PORCOND = 1.0
         RETURN
      ENDIF
      up = .FALSE.
      IF ( LSAME(Uplo,'U') ) up = .TRUE.
!
!     Compute the equilibration matrix R such that
!     inv(R)*A*C has unit 1-norm.
!
      IF ( up ) THEN
         DO i = 1 , N
            tmp = 0.0
            IF ( Cmode==1 ) THEN
               DO j = 1 , i
                  tmp = tmp + ABS(A(j,i)*C(j))
               ENDDO
               DO j = i + 1 , N
                  tmp = tmp + ABS(A(i,j)*C(j))
               ENDDO
            ELSEIF ( Cmode==0 ) THEN
               DO j = 1 , i
                  tmp = tmp + ABS(A(j,i))
               ENDDO
               DO j = i + 1 , N
                  tmp = tmp + ABS(A(i,j))
               ENDDO
            ELSE
               DO j = 1 , i
                  tmp = tmp + ABS(A(j,i)/C(j))
               ENDDO
               DO j = i + 1 , N
                  tmp = tmp + ABS(A(i,j)/C(j))
               ENDDO
            ENDIF
            Work(2*N+i) = tmp
         ENDDO
      ELSE
         DO i = 1 , N
            tmp = 0.0
            IF ( Cmode==1 ) THEN
               DO j = 1 , i
                  tmp = tmp + ABS(A(i,j)*C(j))
               ENDDO
               DO j = i + 1 , N
                  tmp = tmp + ABS(A(j,i)*C(j))
               ENDDO
            ELSEIF ( Cmode==0 ) THEN
               DO j = 1 , i
                  tmp = tmp + ABS(A(i,j))
               ENDDO
               DO j = i + 1 , N
                  tmp = tmp + ABS(A(j,i))
               ENDDO
            ELSE
               DO j = 1 , i
                  tmp = tmp + ABS(A(i,j)/C(j))
               ENDDO
               DO j = i + 1 , N
                  tmp = tmp + ABS(A(j,i)/C(j))
               ENDDO
            ENDIF
            Work(2*N+i) = tmp
         ENDDO
      ENDIF
!
!     Estimate the norm of inv(op(A)).
!
      ainvnm = 0.0
 
      kase = 0
      DO
         CALL SLACN2(N,Work(N+1),Work,Iwork,ainvnm,kase,isave)
         IF ( kase/=0 ) THEN
            IF ( kase==2 ) THEN
!
!           Multiply by R.
!
               DO i = 1 , N
                  Work(i) = Work(i)*Work(2*N+i)
               ENDDO
 
               IF ( up ) THEN
                  CALL SPOTRS('Upper',N,1,Af,Ldaf,Work,N,Info)
               ELSE
                  CALL SPOTRS('Lower',N,1,Af,Ldaf,Work,N,Info)
               ENDIF
!
!           Multiply by inv(C).
!
               IF ( Cmode==1 ) THEN
                  DO i = 1 , N
                     Work(i) = Work(i)/C(i)
                  ENDDO
               ELSEIF ( Cmode==-1 ) THEN
                  DO i = 1 , N
                     Work(i) = Work(i)*C(i)
                  ENDDO
               ENDIF
            ELSE
!
!           Multiply by inv(C**T).
!
               IF ( Cmode==1 ) THEN
                  DO i = 1 , N
                     Work(i) = Work(i)/C(i)
                  ENDDO
               ELSEIF ( Cmode==-1 ) THEN
                  DO i = 1 , N
                     Work(i) = Work(i)*C(i)
                  ENDDO
               ENDIF
 
               IF ( up ) THEN
                  CALL SPOTRS('Upper',N,1,Af,Ldaf,Work,N,Info)
               ELSE
                  CALL SPOTRS('Lower',N,1,Af,Ldaf,Work,N,Info)
               ENDIF
!
!           Multiply by R.
!
               DO i = 1 , N
                  Work(i) = Work(i)*Work(2*N+i)
               ENDDO
            ENDIF
            CYCLE
         ENDIF
!
!     Compute the estimate of the reciprocal condition number.
!
         IF ( ainvnm/=0.0 ) SLA_PORCOND = (1.0/ainvnm)
         EXIT
      ENDDO
!
!
      END FUNCTION SLA_PORCOND
