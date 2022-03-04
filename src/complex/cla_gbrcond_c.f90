!*==cla_gbrcond_c.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b CLA_GBRCOND_C computes the infinity norm condition number of op(A)*inv(diag(c)) for general banded matrices.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CLA_GBRCOND_C + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cla_gbrcond_c.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cla_gbrcond_c.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cla_gbrcond_c.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       REAL FUNCTION CLA_GBRCOND_C( TRANS, N, KL, KU, AB, LDAB, AFB,
!                                    LDAFB, IPIV, C, CAPPLY, INFO, WORK,
!                                    RWORK )
!
!       .. Scalar Arguments ..
!       CHARACTER          TRANS
!       LOGICAL            CAPPLY
!       INTEGER            N, KL, KU, KD, KE, LDAB, LDAFB, INFO
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * )
!       COMPLEX            AB( LDAB, * ), AFB( LDAFB, * ), WORK( * )
!       REAL               C( * ), RWORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    CLA_GBRCOND_C Computes the infinity norm condition number of
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
!>          AB is COMPLEX array, dimension (LDAB,N)
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
!>          AFB is COMPLEX array, dimension (LDAFB,N)
!>     Details of the LU factorization of the band matrix A, as
!>     computed by CGBTRF.  U is stored as an upper triangular
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
!>     as computed by CGBTRF; row i of the matrix was interchanged
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
!> \ingroup complexGBcomputational
!
!  =====================================================================
      REAL FUNCTION CLA_GBRCOND_C(Trans,N,Kl,Ku,Ab,Ldab,Afb,Ldafb,Ipiv, &
     &                            C,Capply,Info,Work,Rwork)
      IMPLICIT NONE
!*--CLA_GBRCOND_C164
!
!  -- LAPACK computational routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Trans
      LOGICAL Capply
      INTEGER N , Kl , Ku , kd , ke , Ldab , Ldafb , Info
!     ..
!     .. Array Arguments ..
      INTEGER Ipiv(*)
      COMPLEX Ab(Ldab,*) , Afb(Ldafb,*) , Work(*)
      REAL C(*) , Rwork(*)
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
      LOGICAL notrans
      INTEGER kase , i , j
      REAL ainvnm , anorm , tmp
      COMPLEX zdum
!     ..
!     .. Local Arrays ..
      INTEGER isave(3)
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL CLACN2 , CGBTRS , XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , MAX
!     ..
!     .. Statement Functions ..
      REAL CABS1
!     ..
!     .. Statement Function Definitions ..
      CABS1(zdum) = ABS(REAL(zdum)) + ABS(AIMAG(zdum))
!     ..
!     .. Executable Statements ..
      CLA_GBRCOND_C = 0.0E+0
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
         CALL XERBLA('CLA_GBRCOND_C',-Info)
         RETURN
      ENDIF
!
!     Compute norm of op(A)*op2(C).
!
      anorm = 0.0E+0
      kd = Ku + 1
      ke = Kl + 1
      IF ( notrans ) THEN
         DO i = 1 , N
            tmp = 0.0E+0
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
            tmp = 0.0E+0
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
         CLA_GBRCOND_C = 1.0E+0
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
                  CALL CGBTRS('No transpose',N,Kl,Ku,1,Afb,Ldafb,Ipiv,  &
     &                        Work,N,Info)
               ELSE
                  CALL CGBTRS('Conjugate transpose',N,Kl,Ku,1,Afb,Ldafb,&
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
                  CALL CGBTRS('Conjugate transpose',N,Kl,Ku,1,Afb,Ldafb,&
     &                        Ipiv,Work,N,Info)
               ELSE
                  CALL CGBTRS('No transpose',N,Kl,Ku,1,Afb,Ldafb,Ipiv,  &
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
         IF ( ainvnm/=0.0E+0 ) CLA_GBRCOND_C = 1.0E+0/ainvnm
         EXIT
      ENDDO
!
!
      END FUNCTION CLA_GBRCOND_C
