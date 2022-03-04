!*==cla_gercond_c.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b CLA_GERCOND_C computes the infinity norm condition number of op(A)*inv(diag(c)) for general matrices.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CLA_GERCOND_C + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cla_gercond_c.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cla_gercond_c.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cla_gercond_c.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       REAL FUNCTION CLA_GERCOND_C( TRANS, N, A, LDA, AF, LDAF, IPIV, C,
!                                    CAPPLY, INFO, WORK, RWORK )
!
!       .. Scalar Arguments ..
!       CHARACTER          TRANS
!       LOGICAL            CAPPLY
!       INTEGER            N, LDA, LDAF, INFO
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * )
!       COMPLEX            A( LDA, * ), AF( LDAF, * ), WORK( * )
!       REAL               C( * ), RWORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>
!>    CLA_GERCOND_C computes the infinity norm condition number of
!>    op(A) * inv(diag(C)) where C is a REAL vector.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] TRANS
!> \verbatim
!>          TRANS is CHARACTER*1
!>     Specifies the form of the system of equations:
!>       = 'N':  A * X = B     (No transpose)
!>       = 'T':  A**T * X = B  (Transpose)
!>       = 'C':  A**H * X = B  (Conjugate Transpose = Transpose)
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
!>          A is COMPLEX array, dimension (LDA,N)
!>     On entry, the N-by-N matrix A
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
!>          AF is COMPLEX array, dimension (LDAF,N)
!>     The factors L and U from the factorization
!>     A = P*L*U as computed by CGETRF.
!> \endverbatim
!>
!> \param[in] LDAF
!> \verbatim
!>          LDAF is INTEGER
!>     The leading dimension of the array AF.  LDAF >= max(1,N).
!> \endverbatim
!>
!> \param[in] IPIV
!> \verbatim
!>          IPIV is INTEGER array, dimension (N)
!>     The pivot indices from the factorization A = P*L*U
!>     as computed by CGETRF; row i of the matrix was interchanged
!>     with row IPIV(i).
!> \endverbatim
!>
!> \param[in] C
!> \verbatim
!>          C is REAL array, dimension (N)
!>     The vector C in the formula op(A) * inv(diag(C)).
!> \endverbatim
!>
!> \param[in] CAPPLY
!> \verbatim
!>          CAPPLY is LOGICAL
!>     If .TRUE. then access the vector C in the formula above.
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
!>          WORK is COMPLEX array, dimension (2*N).
!>     Workspace.
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is REAL array, dimension (N).
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
!> \ingroup complexGEcomputational
!
!  =====================================================================
      FUNCTION CLA_GERCOND_C(Trans,N,A,Lda,Af,Ldaf,Ipiv,C,Capply,Info,  &
     &                       Work,Rwork)
      USE S_CGETRS
      USE S_CLACN2
      USE S_LSAME
      USE S_XERBLA
      IMPLICIT NONE
!*--CLA_GERCOND_C150
      REAL :: CLA_GERCOND_C
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Trans
      INTEGER :: N
      COMPLEX , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      COMPLEX , DIMENSION(Ldaf,*) :: Af
      INTEGER :: Ldaf
      INTEGER , DIMENSION(*) :: Ipiv
      REAL , INTENT(IN) , DIMENSION(*) :: C
      LOGICAL , INTENT(IN) :: Capply
      INTEGER , INTENT(INOUT) :: Info
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Work
      REAL , INTENT(INOUT) , DIMENSION(*) :: Rwork
!
! Local variable declarations rewritten by SPAG
!
      REAL :: ainvnm , anorm , tmp
      REAL :: CABS1
      INTEGER :: i , j , kase
      INTEGER , DIMENSION(3) :: isave
      LOGICAL :: notrans
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
      CABS1(zdum) = ABS(REAL(zdum)) + ABS(AIMAG(zdum))
!     ..
!     .. Executable Statements ..
      CLA_GERCOND_C = 0.0E+0
!
      Info = 0
      notrans = LSAME(Trans,'N')
      IF ( .NOT.notrans .AND. .NOT.LSAME(Trans,'T') .AND.               &
     &     .NOT.LSAME(Trans,'C') ) THEN
         Info = -1
      ELSEIF ( N<0 ) THEN
         Info = -2
      ELSEIF ( Lda<MAX(1,N) ) THEN
         Info = -4
      ELSEIF ( Ldaf<MAX(1,N) ) THEN
         Info = -6
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('CLA_GERCOND_C',-Info)
         RETURN
      ENDIF
!
!     Compute norm of op(A)*op2(C).
!
      anorm = 0.0E+0
      IF ( notrans ) THEN
         DO i = 1 , N
            tmp = 0.0E+0
            IF ( Capply ) THEN
               DO j = 1 , N
                  tmp = tmp + CABS1(A(i,j))/C(j)
               ENDDO
            ELSE
               DO j = 1 , N
                  tmp = tmp + CABS1(A(i,j))
               ENDDO
            ENDIF
            Rwork(i) = tmp
            anorm = MAX(anorm,tmp)
         ENDDO
      ELSE
         DO i = 1 , N
            tmp = 0.0E+0
            IF ( Capply ) THEN
               DO j = 1 , N
                  tmp = tmp + CABS1(A(j,i))/C(j)
               ENDDO
            ELSE
               DO j = 1 , N
                  tmp = tmp + CABS1(A(j,i))
               ENDDO
            ENDIF
            Rwork(i) = tmp
            anorm = MAX(anorm,tmp)
         ENDDO
      ENDIF
!
!     Quick return if possible.
!
      IF ( N==0 ) THEN
         CLA_GERCOND_C = 1.0E+0
         RETURN
      ELSEIF ( anorm==0.0E+0 ) THEN
         RETURN
      ENDIF
!
!     Estimate the norm of inv(op(A)).
!
      ainvnm = 0.0E+0
!
      kase = 0
      DO
         CALL CLACN2(N,Work(N+1),Work,ainvnm,kase,isave)
         IF ( kase/=0 ) THEN
            IF ( kase==2 ) THEN
!
!           Multiply by R.
!
               DO i = 1 , N
                  Work(i) = Work(i)*Rwork(i)
               ENDDO
!
               IF ( notrans ) THEN
                  CALL CGETRS('No transpose',N,1,Af,Ldaf,Ipiv,Work,N,   &
     &                        Info)
               ELSE
                  CALL CGETRS('Conjugate transpose',N,1,Af,Ldaf,Ipiv,   &
     &                        Work,N,Info)
               ENDIF
!
!           Multiply by inv(C).
!
               IF ( Capply ) THEN
                  DO i = 1 , N
                     Work(i) = Work(i)*C(i)
                  ENDDO
               ENDIF
            ELSE
!
!           Multiply by inv(C**H).
!
               IF ( Capply ) THEN
                  DO i = 1 , N
                     Work(i) = Work(i)*C(i)
                  ENDDO
               ENDIF
!
               IF ( notrans ) THEN
                  CALL CGETRS('Conjugate transpose',N,1,Af,Ldaf,Ipiv,   &
     &                        Work,N,Info)
               ELSE
                  CALL CGETRS('No transpose',N,1,Af,Ldaf,Ipiv,Work,N,   &
     &                        Info)
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
         IF ( ainvnm/=0.0E+0 ) CLA_GERCOND_C = 1.0E+0/ainvnm
         EXIT
      ENDDO
!
!
      END FUNCTION CLA_GERCOND_C