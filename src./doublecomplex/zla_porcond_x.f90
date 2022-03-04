!*==zla_porcond_x.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b ZLA_PORCOND_X computes the infinity norm condition number of op(A)*diag(x) for Hermitian positive-definite matrices.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZLA_PORCOND_X + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zla_porcond_x.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zla_porcond_x.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zla_porcond_x.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       DOUBLE PRECISION FUNCTION ZLA_PORCOND_X( UPLO, N, A, LDA, AF,
!                                                LDAF, X, INFO, WORK,
!                                                RWORK )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            N, LDA, LDAF, INFO
!       ..
!       .. Array Arguments ..
!       COMPLEX*16         A( LDA, * ), AF( LDAF, * ), WORK( * ), X( * )
!       DOUBLE PRECISION   RWORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    ZLA_PORCOND_X Computes the infinity norm condition number of
!>    op(A) * diag(X) where X is a COMPLEX*16 vector.
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
!>          A is COMPLEX*16 array, dimension (LDA,N)
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
!>          AF is COMPLEX*16 array, dimension (LDAF,N)
!>     The triangular factor U or L from the Cholesky factorization
!>     A = U**H*U or A = L*L**H, as computed by ZPOTRF.
!> \endverbatim
!>
!> \param[in] LDAF
!> \verbatim
!>          LDAF is INTEGER
!>     The leading dimension of the array AF.  LDAF >= max(1,N).
!> \endverbatim
!>
!> \param[in] X
!> \verbatim
!>          X is COMPLEX*16 array, dimension (N)
!>     The vector X in the formula op(A) * diag(X).
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
!>          WORK is COMPLEX*16 array, dimension (2*N).
!>     Workspace.
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is DOUBLE PRECISION array, dimension (N).
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
!> \ingroup complex16POcomputational
!
!  =====================================================================
      FUNCTION ZLA_PORCOND_X(Uplo,N,A,Lda,Af,Ldaf,X,Info,Work,Rwork)
      USE F77KINDS                        
      USE S_LSAME
      USE S_XERBLA
      USE S_ZLACN2
      USE S_ZPOTRS
      IMPLICIT NONE
!*--ZLA_PORCOND_X132
      REAL(R8KIND) :: ZLA_PORCOND_X
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Uplo
      INTEGER :: N
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      COMPLEX(CX16KIND) , DIMENSION(Ldaf,*) :: Af
      INTEGER :: Ldaf
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(*) :: X
      INTEGER , INTENT(INOUT) :: Info
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Rwork
!
! Local variable declarations rewritten by SPAG
!
      REAL(R8KIND) :: ainvnm , anorm , tmp
      REAL(R8KIND) :: CABS1
      INTEGER :: i , j , kase
      INTEGER , DIMENSION(3) :: isave
      LOGICAL :: up , upper
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
!     .. Statement Function Definitions ..
      CABS1(zdum) = ABS(DBLE(zdum)) + ABS(DIMAG(zdum))
!     ..
!     .. Executable Statements ..
!
      ZLA_PORCOND_X = 0.0D+0
!
      Info = 0
      upper = LSAME(Uplo,'U')
      IF ( .NOT.upper .AND. .NOT.LSAME(Uplo,'L') ) THEN
         Info = -1
      ELSEIF ( N<0 ) THEN
         Info = -2
      ELSEIF ( Lda<MAX(1,N) ) THEN
         Info = -4
      ELSEIF ( Ldaf<MAX(1,N) ) THEN
         Info = -6
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('ZLA_PORCOND_X',-Info)
         RETURN
      ENDIF
      up = .FALSE.
      IF ( LSAME(Uplo,'U') ) up = .TRUE.
!
!     Compute norm of op(A)*op2(C).
!
      anorm = 0.0D+0
      IF ( up ) THEN
         DO i = 1 , N
            tmp = 0.0D+0
            DO j = 1 , i
               tmp = tmp + CABS1(A(j,i)*X(j))
            ENDDO
            DO j = i + 1 , N
               tmp = tmp + CABS1(A(i,j)*X(j))
            ENDDO
            Rwork(i) = tmp
            anorm = MAX(anorm,tmp)
         ENDDO
      ELSE
         DO i = 1 , N
            tmp = 0.0D+0
            DO j = 1 , i
               tmp = tmp + CABS1(A(i,j)*X(j))
            ENDDO
            DO j = i + 1 , N
               tmp = tmp + CABS1(A(j,i)*X(j))
            ENDDO
            Rwork(i) = tmp
            anorm = MAX(anorm,tmp)
         ENDDO
      ENDIF
!
!     Quick return if possible.
!
      IF ( N==0 ) THEN
         ZLA_PORCOND_X = 1.0D+0
         RETURN
      ELSEIF ( anorm==0.0D+0 ) THEN
         RETURN
      ENDIF
!
!     Estimate the norm of inv(op(A)).
!
      ainvnm = 0.0D+0
!
      kase = 0
      DO
         CALL ZLACN2(N,Work(N+1),Work,ainvnm,kase,isave)
         IF ( kase/=0 ) THEN
            IF ( kase==2 ) THEN
!
!           Multiply by R.
!
               DO i = 1 , N
                  Work(i) = Work(i)*Rwork(i)
               ENDDO
!
               IF ( up ) THEN
                  CALL ZPOTRS('U',N,1,Af,Ldaf,Work,N,Info)
               ELSE
                  CALL ZPOTRS('L',N,1,Af,Ldaf,Work,N,Info)
               ENDIF
!
!           Multiply by inv(X).
!
               DO i = 1 , N
                  Work(i) = Work(i)/X(i)
               ENDDO
            ELSE
!
!           Multiply by inv(X**H).
!
               DO i = 1 , N
                  Work(i) = Work(i)/X(i)
               ENDDO
!
               IF ( up ) THEN
                  CALL ZPOTRS('U',N,1,Af,Ldaf,Work,N,Info)
               ELSE
                  CALL ZPOTRS('L',N,1,Af,Ldaf,Work,N,Info)
               ENDIF
!
!           Multiply by R.
!
               DO i = 1 , N
                  Work(i) = Work(i)*Rwork(i)
               ENDDO
            ENDIF
            CYCLE
         ENDIF
!
!     Compute the estimate of the reciprocal condition number.
!
         IF ( ainvnm/=0.0D+0 ) ZLA_PORCOND_X = 1.0D+0/ainvnm
         EXIT
      ENDDO
!
!
      END FUNCTION ZLA_PORCOND_X
