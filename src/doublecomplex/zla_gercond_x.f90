!*==zla_gercond_x.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b ZLA_GERCOND_X computes the infinity norm condition number of op(A)*diag(x) for general matrices.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZLA_GERCOND_X + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zla_gercond_x.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zla_gercond_x.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zla_gercond_x.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       DOUBLE PRECISION FUNCTION ZLA_GERCOND_X( TRANS, N, A, LDA, AF,
!                                                LDAF, IPIV, X, INFO,
!                                                WORK, RWORK )
!
!       .. Scalar Arguments ..
!       CHARACTER          TRANS
!       INTEGER            N, LDA, LDAF, INFO
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * )
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
!>    ZLA_GERCOND_X computes the infinity norm condition number of
!>    op(A) * diag(X) where X is a COMPLEX*16 vector.
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
!>     The factors L and U from the factorization
!>     A = P*L*U as computed by ZGETRF.
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
!>     as computed by ZGETRF; row i of the matrix was interchanged
!>     with row IPIV(i).
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
!> \ingroup complex16GEcomputational
!
!  =====================================================================
      DOUBLE PRECISION FUNCTION ZLA_GERCOND_X(Trans,N,A,Lda,Af,Ldaf,    &
     &   Ipiv,X,Info,Work,Rwork)
      IMPLICIT NONE
!*--ZLA_GERCOND_X139
!
!  -- LAPACK computational routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Trans
      INTEGER N , Lda , Ldaf , Info
!     ..
!     .. Array Arguments ..
      INTEGER Ipiv(*)
      COMPLEX*16 A(Lda,*) , Af(Ldaf,*) , Work(*) , X(*)
      DOUBLE PRECISION Rwork(*)
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
      LOGICAL notrans
      INTEGER kase
      DOUBLE PRECISION ainvnm , anorm , tmp
      INTEGER i , j
      COMPLEX*16 zdum
!     ..
!     .. Local Arrays ..
      INTEGER isave(3)
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL ZLACN2 , ZGETRS , XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , MAX , REAL , DIMAG
!     ..
!     .. Statement Functions ..
      DOUBLE PRECISION CABS1
!     ..
!     .. Statement Function Definitions ..
      CABS1(zdum) = ABS(DBLE(zdum)) + ABS(DIMAG(zdum))
!     ..
!     .. Executable Statements ..
!
      ZLA_GERCOND_X = 0.0D+0
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
         CALL XERBLA('ZLA_GERCOND_X',-Info)
         RETURN
      ENDIF
!
!     Compute norm of op(A)*op2(C).
!
      anorm = 0.0D+0
      IF ( notrans ) THEN
         DO i = 1 , N
            tmp = 0.0D+0
            DO j = 1 , N
               tmp = tmp + CABS1(A(i,j)*X(j))
            ENDDO
            Rwork(i) = tmp
            anorm = MAX(anorm,tmp)
         ENDDO
      ELSE
         DO i = 1 , N
            tmp = 0.0D+0
            DO j = 1 , N
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
         ZLA_GERCOND_X = 1.0D+0
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
!           Multiply by R.
               DO i = 1 , N
                  Work(i) = Work(i)*Rwork(i)
               ENDDO
!
               IF ( notrans ) THEN
                  CALL ZGETRS('No transpose',N,1,Af,Ldaf,Ipiv,Work,N,   &
     &                        Info)
               ELSE
                  CALL ZGETRS('Conjugate transpose',N,1,Af,Ldaf,Ipiv,   &
     &                        Work,N,Info)
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
               IF ( notrans ) THEN
                  CALL ZGETRS('Conjugate transpose',N,1,Af,Ldaf,Ipiv,   &
     &                        Work,N,Info)
               ELSE
                  CALL ZGETRS('No transpose',N,1,Af,Ldaf,Ipiv,Work,N,   &
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
         IF ( ainvnm/=0.0D+0 ) ZLA_GERCOND_X = 1.0D+0/ainvnm
         EXIT
      ENDDO
!
!
      END FUNCTION ZLA_GERCOND_X
