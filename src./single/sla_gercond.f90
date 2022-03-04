!*==sla_gercond.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b SLA_GERCOND estimates the Skeel condition number for a general matrix.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SLA_GERCOND + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sla_gercond.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sla_gercond.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sla_gercond.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       REAL FUNCTION SLA_GERCOND ( TRANS, N, A, LDA, AF, LDAF, IPIV,
!                                   CMODE, C, INFO, WORK, IWORK )
!
!       .. Scalar Arguments ..
!       CHARACTER          TRANS
!       INTEGER            N, LDA, LDAF, INFO, CMODE
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * ), IWORK( * )
!       REAL               A( LDA, * ), AF( LDAF, * ), WORK( * ),
!      $                   C( * )
!      ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    SLA_GERCOND estimates the Skeel condition number of op(A) * op2(C)
!>    where op2 is determined by CMODE as follows
!>    CMODE =  1    op2(C) = C
!>    CMODE =  0    op2(C) = I
!>    CMODE = -1    op2(C) = inv(C)
!>    The Skeel condition number cond(A) = norminf( |inv(A)||A| )
!>    is computed by computing scaling factors R such that
!>    diag(R)*A*op2(C) is row equilibrated and computing the standard
!>    infinity-norm condition number.
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
!>     The factors L and U from the factorization
!>     A = P*L*U as computed by SGETRF.
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
!>     as computed by SGETRF; row i of the matrix was interchanged
!>     with row IPIV(i).
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
!>     Workspace.2
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
!> \ingroup realGEcomputational
!
!  =====================================================================
      FUNCTION SLA_GERCOND(Trans,N,A,Lda,Af,Ldaf,Ipiv,Cmode,C,Info,Work,&
     &                     Iwork)
      USE S_LSAME
      USE S_SGETRS
      USE S_SLACN2
      USE S_XERBLA
      IMPLICIT NONE
!*--SLA_GERCOND158
      REAL :: SLA_GERCOND
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Trans
      INTEGER :: N
      REAL , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER , INTENT(IN) :: Lda
      REAL , DIMENSION(Ldaf,*) :: Af
      INTEGER :: Ldaf
      INTEGER , DIMENSION(*) :: Ipiv
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
      LOGICAL :: notrans
!
! End of declarations rewritten by SPAG
!
!     ..
!     .. Array Arguments ..
!    ..
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
!     .. Executable Statements ..
!
      SLA_GERCOND = 0.0
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
         CALL XERBLA('SLA_GERCOND',-Info)
         RETURN
      ENDIF
      IF ( N==0 ) THEN
         SLA_GERCOND = 1.0
         RETURN
      ENDIF
!
!     Compute the equilibration matrix R such that
!     inv(R)*A*C has unit 1-norm.
!
      IF ( notrans ) THEN
         DO i = 1 , N
            tmp = 0.0
            IF ( Cmode==1 ) THEN
               DO j = 1 , N
                  tmp = tmp + ABS(A(i,j)*C(j))
               ENDDO
            ELSEIF ( Cmode==0 ) THEN
               DO j = 1 , N
                  tmp = tmp + ABS(A(i,j))
               ENDDO
            ELSE
               DO j = 1 , N
                  tmp = tmp + ABS(A(i,j)/C(j))
               ENDDO
            ENDIF
            Work(2*N+i) = tmp
         ENDDO
      ELSE
         DO i = 1 , N
            tmp = 0.0
            IF ( Cmode==1 ) THEN
               DO j = 1 , N
                  tmp = tmp + ABS(A(j,i)*C(j))
               ENDDO
            ELSEIF ( Cmode==0 ) THEN
               DO j = 1 , N
                  tmp = tmp + ABS(A(j,i))
               ENDDO
            ELSE
               DO j = 1 , N
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
 
               IF ( notrans ) THEN
                  CALL SGETRS('No transpose',N,1,Af,Ldaf,Ipiv,Work,N,   &
     &                        Info)
               ELSE
                  CALL SGETRS('Transpose',N,1,Af,Ldaf,Ipiv,Work,N,Info)
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
 
               IF ( notrans ) THEN
                  CALL SGETRS('Transpose',N,1,Af,Ldaf,Ipiv,Work,N,Info)
               ELSE
                  CALL SGETRS('No transpose',N,1,Af,Ldaf,Ipiv,Work,N,   &
     &                        Info)
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
         IF ( ainvnm/=0.0 ) SLA_GERCOND = (1.0/ainvnm)
         EXIT
      ENDDO
!
!
      END FUNCTION SLA_GERCOND