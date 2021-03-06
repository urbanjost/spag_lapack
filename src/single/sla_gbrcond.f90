!*==sla_gbrcond.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b SLA_GBRCOND estimates the Skeel condition number for a general banded matrix.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SLA_GBRCOND + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sla_gbrcond.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sla_gbrcond.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sla_gbrcond.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       REAL FUNCTION SLA_GBRCOND( TRANS, N, KL, KU, AB, LDAB, AFB, LDAFB,
!                                  IPIV, CMODE, C, INFO, WORK, IWORK )
!
!       .. Scalar Arguments ..
!       CHARACTER          TRANS
!       INTEGER            N, LDAB, LDAFB, INFO, KL, KU, CMODE
!       ..
!       .. Array Arguments ..
!       INTEGER            IWORK( * ), IPIV( * )
!       REAL               AB( LDAB, * ), AFB( LDAFB, * ), WORK( * ),
!      $                   C( * )
!      ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    SLA_GBRCOND Estimates the Skeel condition number of  op(A) * op2(C)
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
!>          AB is REAL array, dimension (LDAB,N)
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
!>          AFB is REAL array, dimension (LDAFB,N)
!>     Details of the LU factorization of the band matrix A, as
!>     computed by SGBTRF.  U is stored as an upper triangular
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
!>     as computed by SGBTRF; row i of the matrix was interchanged
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
!>          WORK is REAL array, dimension (5*N).
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
!> \ingroup realGBcomputational
!
!  =====================================================================
      REAL FUNCTION SLA_GBRCOND(Trans,N,Kl,Ku,Ab,Ldab,Afb,Ldafb,Ipiv,   &
     &                          Cmode,C,Info,Work,Iwork)
      IMPLICIT NONE
!*--SLA_GBRCOND172
!
!  -- LAPACK computational routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Trans
      INTEGER N , Ldab , Ldafb , Info , Kl , Ku , Cmode
!     ..
!     .. Array Arguments ..
      INTEGER Iwork(*) , Ipiv(*)
      REAL Ab(Ldab,*) , Afb(Ldafb,*) , Work(*) , C(*)
!    ..
!
!  =====================================================================
!
!     .. Local Scalars ..
      LOGICAL notrans
      INTEGER kase , i , j , kd , ke
      REAL ainvnm , tmp
!     ..
!     .. Local Arrays ..
      INTEGER isave(3)
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL SLACN2 , SGBTRS , XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , MAX
!     ..
!     .. Executable Statements ..
!
      SLA_GBRCOND = 0.0
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
         CALL XERBLA('SLA_GBRCOND',-Info)
         RETURN
      ENDIF
      IF ( N==0 ) THEN
         SLA_GBRCOND = 1.0
         RETURN
      ENDIF
!
!     Compute the equilibration matrix R such that
!     inv(R)*A*C has unit 1-norm.
!
      kd = Ku + 1
      ke = Kl + 1
      IF ( notrans ) THEN
         DO i = 1 , N
            tmp = 0.0
            IF ( Cmode==1 ) THEN
               DO j = MAX(i-Kl,1) , MIN(i+Ku,N)
                  tmp = tmp + ABS(Ab(kd+i-j,j)*C(j))
               ENDDO
            ELSEIF ( Cmode==0 ) THEN
               DO j = MAX(i-Kl,1) , MIN(i+Ku,N)
                  tmp = tmp + ABS(Ab(kd+i-j,j))
               ENDDO
            ELSE
               DO j = MAX(i-Kl,1) , MIN(i+Ku,N)
                  tmp = tmp + ABS(Ab(kd+i-j,j)/C(j))
               ENDDO
            ENDIF
            Work(2*N+i) = tmp
         ENDDO
      ELSE
         DO i = 1 , N
            tmp = 0.0
            IF ( Cmode==1 ) THEN
               DO j = MAX(i-Kl,1) , MIN(i+Ku,N)
                  tmp = tmp + ABS(Ab(ke-i+j,i)*C(j))
               ENDDO
            ELSEIF ( Cmode==0 ) THEN
               DO j = MAX(i-Kl,1) , MIN(i+Ku,N)
                  tmp = tmp + ABS(Ab(ke-i+j,i))
               ENDDO
            ELSE
               DO j = MAX(i-Kl,1) , MIN(i+Ku,N)
                  tmp = tmp + ABS(Ab(ke-i+j,i)/C(j))
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
                  CALL SGBTRS('No transpose',N,Kl,Ku,1,Afb,Ldafb,Ipiv,  &
     &                        Work,N,Info)
               ELSE
                  CALL SGBTRS('Transpose',N,Kl,Ku,1,Afb,Ldafb,Ipiv,Work,&
     &                        N,Info)
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
                  CALL SGBTRS('Transpose',N,Kl,Ku,1,Afb,Ldafb,Ipiv,Work,&
     &                        N,Info)
               ELSE
                  CALL SGBTRS('No transpose',N,Kl,Ku,1,Afb,Ldafb,Ipiv,  &
     &                        Work,N,Info)
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
         IF ( ainvnm/=0.0 ) SLA_GBRCOND = (1.0/ainvnm)
         EXIT
      ENDDO
!
!
      END FUNCTION SLA_GBRCOND
