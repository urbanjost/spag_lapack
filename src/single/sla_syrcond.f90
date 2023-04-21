!*==sla_syrcond.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b SLA_SYRCOND estimates the Skeel condition number for a symmetric indefinite matrix.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SLA_SYRCOND + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sla_syrcond.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sla_syrcond.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sla_syrcond.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       REAL FUNCTION SLA_SYRCOND( UPLO, N, A, LDA, AF, LDAF, IPIV, CMODE,
!                                  C, INFO, WORK, IWORK )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            N, LDA, LDAF, INFO, CMODE
!       ..
!       .. Array Arguments
!       INTEGER            IWORK( * ), IPIV( * )
!       REAL               A( LDA, * ), AF( LDAF, * ), WORK( * ), C( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    SLA_SYRCOND estimates the Skeel condition number of  op(A) * op2(C)
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
!>     The block diagonal matrix D and the multipliers used to
!>     obtain the factor U or L as computed by SSYTRF.
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
!>     Details of the interchanges and the block structure of D
!>     as determined by SSYTRF.
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
!> \ingroup realSYcomputational
!
!  =====================================================================
      REAL FUNCTION SLA_SYRCOND(Uplo,N,A,Lda,Af,Ldaf,Ipiv,Cmode,C,Info, &
     &                          Work,Iwork)
      IMPLICIT NONE
!*--SLA_SYRCOND150
!
!  -- LAPACK computational routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Uplo
      INTEGER N , Lda , Ldaf , Info , Cmode
!     ..
!     .. Array Arguments
      INTEGER Iwork(*) , Ipiv(*)
      REAL A(Lda,*) , Af(Ldaf,*) , Work(*) , C(*)
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
      CHARACTER normin
      INTEGER kase , i , j
      REAL ainvnm , smlnum , tmp
      LOGICAL up
!     ..
!     .. Local Arrays ..
      INTEGER isave(3)
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      REAL SLAMCH
      EXTERNAL LSAME , SLAMCH
!     ..
!     .. External Subroutines ..
      EXTERNAL SLACN2 , XERBLA , SSYTRS
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , MAX
!     ..
!     .. Executable Statements ..
!
      SLA_SYRCOND = 0.0
!
      Info = 0
      IF ( N<0 ) THEN
         Info = -2
      ELSEIF ( Lda<MAX(1,N) ) THEN
         Info = -4
      ELSEIF ( Ldaf<MAX(1,N) ) THEN
         Info = -6
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('SLA_SYRCOND',-Info)
         RETURN
      ENDIF
      IF ( N==0 ) THEN
         SLA_SYRCOND = 1.0
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
      smlnum = SLAMCH('Safe minimum')
      ainvnm = 0.0
      normin = 'N'
 
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
                  CALL SSYTRS('U',N,1,Af,Ldaf,Ipiv,Work,N,Info)
               ELSE
                  CALL SSYTRS('L',N,1,Af,Ldaf,Ipiv,Work,N,Info)
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
                  CALL SSYTRS('U',N,1,Af,Ldaf,Ipiv,Work,N,Info)
               ELSE
                  CALL SSYTRS('L',N,1,Af,Ldaf,Ipiv,Work,N,Info)
               ENDIF
!
!           Multiply by R.
!
               DO i = 1 , N
                  Work(i) = Work(i)*Work(2*N+i)
               ENDDO
            ENDIF
!
            CYCLE
         ENDIF
!
!     Compute the estimate of the reciprocal condition number.
!
         IF ( ainvnm/=0.0 ) SLA_SYRCOND = (1.0/ainvnm)
         EXIT
      ENDDO
!
!
      END FUNCTION SLA_SYRCOND
