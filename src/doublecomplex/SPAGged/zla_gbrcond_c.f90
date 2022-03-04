!*==zla_gbrcond_c.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b ZLA_GBRCOND_C computes the infinity norm condition number of op(A)*inv(diag(c)) for general banded matrices.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZLA_GBRCOND_C + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zla_gbrcond_c.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zla_gbrcond_c.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zla_gbrcond_c.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       DOUBLE PRECISION FUNCTION ZLA_GBRCOND_C( TRANS, N, KL, KU, AB,
!                                                LDAB, AFB, LDAFB, IPIV,
!                                                C, CAPPLY, INFO, WORK,
!                                                RWORK )
!
!       .. Scalar Arguments ..
!       CHARACTER          TRANS
!       LOGICAL            CAPPLY
!       INTEGER            N, KL, KU, KD, KE, LDAB, LDAFB, INFO
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * )
!       COMPLEX*16         AB( LDAB, * ), AFB( LDAFB, * ), WORK( * )
!       DOUBLE PRECISION   C( * ), RWORK( * )
!
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    ZLA_GBRCOND_C Computes the infinity norm condition number of
!>    op(A) * inv(diag(C)) where C is a DOUBLE PRECISION vector.
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
!> \param[in] KL
!> \verbatim
!>          KL is INTEGER
!>     The number of subdiagonals within the band of A.  KL >= 0.
!> \endverbatim
!>
!> \param[in] KU
!> \verbatim
!>          KU is INTEGER
!>     The number of superdiagonals within the band of A.  KU >= 0.
!> \endverbatim
!>
!> \param[in] AB
!> \verbatim
!>          AB is COMPLEX*16 array, dimension (LDAB,N)
!>     On entry, the matrix A in band storage, in rows 1 to KL+KU+1.
!>     The j-th column of A is stored in the j-th column of the
!>     array AB as follows:
!>     AB(KU+1+i-j,j) = A(i,j) for max(1,j-KU)<=i<=min(N,j+kl)
!> \endverbatim
!>
!> \param[in] LDAB
!> \verbatim
!>          LDAB is INTEGER
!>     The leading dimension of the array AB.  LDAB >= KL+KU+1.
!> \endverbatim
!>
!> \param[in] AFB
!> \verbatim
!>          AFB is COMPLEX*16 array, dimension (LDAFB,N)
!>     Details of the LU factorization of the band matrix A, as
!>     computed by ZGBTRF.  U is stored as an upper triangular
!>     band matrix with KL+KU superdiagonals in rows 1 to KL+KU+1,
!>     and the multipliers used during the factorization are stored
!>     in rows KL+KU+2 to 2*KL+KU+1.
!> \endverbatim
!>
!> \param[in] LDAFB
!> \verbatim
!>          LDAFB is INTEGER
!>     The leading dimension of the array AFB.  LDAFB >= 2*KL+KU+1.
!> \endverbatim
!>
!> \param[in] IPIV
!> \verbatim
!>          IPIV is INTEGER array, dimension (N)
!>     The pivot indices from the factorization A = P*L*U
!>     as computed by ZGBTRF; row i of the matrix was interchanged
!>     with row IPIV(i).
!> \endverbatim
!>
!> \param[in] C
!> \verbatim
!>          C is DOUBLE PRECISION array, dimension (N)
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
!> \ingroup complex16GBcomputational
!
!  =====================================================================
      FUNCTION ZLA_GBRCOND_C(Trans,N,Kl,Ku,Ab,Ldab,Afb,Ldafb,Ipiv,C,    &
     &                       Capply,Info,Work,Rwork)
      USE F77KINDS                        
      USE S_LSAME
      USE S_XERBLA
      USE S_ZGBTRS
      USE S_ZLACN2
      IMPLICIT NONE
!*--ZLA_GBRCOND_C170
      REAL(R8KIND) :: ZLA_GBRCOND_C
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Trans
      INTEGER :: N
      INTEGER :: Kl
      INTEGER :: Ku
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(Ldab,*) :: Ab
      INTEGER , INTENT(IN) :: Ldab
      COMPLEX(CX16KIND) , DIMENSION(Ldafb,*) :: Afb
      INTEGER :: Ldafb
      INTEGER , DIMENSION(*) :: Ipiv
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: C
      LOGICAL , INTENT(IN) :: Capply
      INTEGER , INTENT(INOUT) :: Info
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Rwork
!
! Local variable declarations rewritten by SPAG
!
      REAL(R8KIND) :: ainvnm , anorm , tmp
      REAL(R8KIND) :: CABS1
      INTEGER :: i , j , kase , kd , ke
      INTEGER , DIMENSION(3) :: isave
      LOGICAL :: notrans
      COMPLEX(CX16KIND) :: zdum
!
! End of declarations rewritten by SPAG
!
!     ..
!     .. Array Arguments ..
!
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
      ZLA_GBRCOND_C = 0.0D+0
!
      Info = 0
      notrans = LSAME(Trans,'N')
      IF ( .NOT.notrans .AND. .NOT.LSAME(Trans,'T') .AND.               &
     &     .NOT.LSAME(Trans,'C') ) THEN
         Info = -1
      ELSEIF ( N<0 ) THEN
         Info = -2
      ELSEIF ( Kl<0 .OR. Kl>N-1 ) THEN
         Info = -3
      ELSEIF ( Ku<0 .OR. Ku>N-1 ) THEN
         Info = -4
      ELSEIF ( Ldab<Kl+Ku+1 ) THEN
         Info = -6
      ELSEIF ( Ldafb<2*Kl+Ku+1 ) THEN
         Info = -8
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('ZLA_GBRCOND_C',-Info)
         RETURN
      ENDIF
!
!     Compute norm of op(A)*op2(C).
!
      anorm = 0.0D+0
      kd = Ku + 1
      ke = Kl + 1
      IF ( notrans ) THEN
         DO i = 1 , N
            tmp = 0.0D+0
            IF ( Capply ) THEN
               DO j = MAX(i-Kl,1) , MIN(i+Ku,N)
                  tmp = tmp + CABS1(Ab(kd+i-j,j))/C(j)
               ENDDO
            ELSE
               DO j = MAX(i-Kl,1) , MIN(i+Ku,N)
                  tmp = tmp + CABS1(Ab(kd+i-j,j))
               ENDDO
            ENDIF
            Rwork(i) = tmp
            anorm = MAX(anorm,tmp)
         ENDDO
      ELSE
         DO i = 1 , N
            tmp = 0.0D+0
            IF ( Capply ) THEN
               DO j = MAX(i-Kl,1) , MIN(i+Ku,N)
                  tmp = tmp + CABS1(Ab(ke-i+j,i))/C(j)
               ENDDO
            ELSE
               DO j = MAX(i-Kl,1) , MIN(i+Ku,N)
                  tmp = tmp + CABS1(Ab(ke-i+j,i))
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
         ZLA_GBRCOND_C = 1.0D+0
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
               IF ( notrans ) THEN
                  CALL ZGBTRS('No transpose',N,Kl,Ku,1,Afb,Ldafb,Ipiv,  &
     &                        Work,N,Info)
               ELSE
                  CALL ZGBTRS('Conjugate transpose',N,Kl,Ku,1,Afb,Ldafb,&
     &                        Ipiv,Work,N,Info)
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
                  CALL ZGBTRS('Conjugate transpose',N,Kl,Ku,1,Afb,Ldafb,&
     &                        Ipiv,Work,N,Info)
               ELSE
                  CALL ZGBTRS('No transpose',N,Kl,Ku,1,Afb,Ldafb,Ipiv,  &
     &                        Work,N,Info)
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
         IF ( ainvnm/=0.0D+0 ) ZLA_GBRCOND_C = 1.0D+0/ainvnm
         EXIT
      ENDDO
!
!
      END FUNCTION ZLA_GBRCOND_C
