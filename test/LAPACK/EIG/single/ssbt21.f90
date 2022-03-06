!*==ssbt21.f90  processed by SPAG 7.51RB at 17:37 on  4 Mar 2022
!> \brief \b SSBT21
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE SSBT21( UPLO, N, KA, KS, A, LDA, D, E, U, LDU, WORK,
!                          RESULT )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            KA, KS, LDA, LDU, N
!       ..
!       .. Array Arguments ..
!       REAL               A( LDA, * ), D( * ), E( * ), RESULT( 2 ),
!      $                   U( LDU, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SSBT21  generally checks a decomposition of the form
!>
!>         A = U S U**T
!>
!> where **T means transpose, A is symmetric banded, U is
!> orthogonal, and S is diagonal (if KS=0) or symmetric
!> tridiagonal (if KS=1).
!>
!> Specifically:
!>
!>         RESULT(1) = | A - U S U**T | / ( |A| n ulp ) and
!>         RESULT(2) = | I - U U**T | / ( n ulp )
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER
!>          If UPLO='U', the upper triangle of A and V will be used and
!>          the (strictly) lower triangle will not be referenced.
!>          If UPLO='L', the lower triangle of A and V will be used and
!>          the (strictly) upper triangle will not be referenced.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The size of the matrix.  If it is zero, SSBT21 does nothing.
!>          It must be at least zero.
!> \endverbatim
!>
!> \param[in] KA
!> \verbatim
!>          KA is INTEGER
!>          The bandwidth of the matrix A.  It must be at least zero.  If
!>          it is larger than N-1, then max( 0, N-1 ) will be used.
!> \endverbatim
!>
!> \param[in] KS
!> \verbatim
!>          KS is INTEGER
!>          The bandwidth of the matrix S.  It may only be zero or one.
!>          If zero, then S is diagonal, and E is not referenced.  If
!>          one, then S is symmetric tri-diagonal.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is REAL array, dimension (LDA, N)
!>          The original (unfactored) matrix.  It is assumed to be
!>          symmetric, and only the upper (UPLO='U') or only the lower
!>          (UPLO='L') will be referenced.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of A.  It must be at least 1
!>          and at least min( KA, N-1 ).
!> \endverbatim
!>
!> \param[in] D
!> \verbatim
!>          D is REAL array, dimension (N)
!>          The diagonal of the (symmetric tri-) diagonal matrix S.
!> \endverbatim
!>
!> \param[in] E
!> \verbatim
!>          E is REAL array, dimension (N-1)
!>          The off-diagonal of the (symmetric tri-) diagonal matrix S.
!>          E(1) is the (1,2) and (2,1) element, E(2) is the (2,3) and
!>          (3,2) element, etc.
!>          Not referenced if KS=0.
!> \endverbatim
!>
!> \param[in] U
!> \verbatim
!>          U is REAL array, dimension (LDU, N)
!>          The orthogonal matrix in the decomposition, expressed as a
!>          dense matrix (i.e., not as a product of Householder
!>          transformations, Givens transformations, etc.)
!> \endverbatim
!>
!> \param[in] LDU
!> \verbatim
!>          LDU is INTEGER
!>          The leading dimension of U.  LDU must be at least N and
!>          at least 1.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is REAL array, dimension (N**2+N)
!> \endverbatim
!>
!> \param[out] RESULT
!> \verbatim
!>          RESULT is REAL array, dimension (2)
!>          The values computed by the two tests described above.  The
!>          values are currently limited to 1/ulp, to avoid overflow.
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
!> \ingroup single_eig
!
!  =====================================================================
      SUBROUTINE SSBT21(Uplo,N,Ka,Ks,A,Lda,D,E,U,Ldu,Work,Result)
      IMPLICIT NONE
!*--SSBT21150
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Uplo
      INTEGER Ka , Ks , Lda , Ldu , N
!     ..
!     .. Array Arguments ..
      REAL A(Lda,*) , D(*) , E(*) , Result(2) , U(Ldu,*) , Work(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ZERO , ONE
      PARAMETER (ZERO=0.0E0,ONE=1.0E0)
!     ..
!     .. Local Scalars ..
      LOGICAL lower
      CHARACTER cuplo
      INTEGER ika , j , jc , jr , lw
      REAL anorm , ulp , unfl , wnorm
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      REAL SLAMCH , SLANGE , SLANSB , SLANSP
      EXTERNAL LSAME , SLAMCH , SLANGE , SLANSB , SLANSP
!     ..
!     .. External Subroutines ..
      EXTERNAL SGEMM , SSPR , SSPR2
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX , MIN , REAL
!     ..
!     .. Executable Statements ..
!
!     Constants
!
      Result(1) = ZERO
      Result(2) = ZERO
      IF ( N<=0 ) RETURN
!
      ika = MAX(0,MIN(N-1,Ka))
      lw = (N*(N+1))/2
!
      IF ( LSAME(Uplo,'U') ) THEN
         lower = .FALSE.
         cuplo = 'U'
      ELSE
         lower = .TRUE.
         cuplo = 'L'
      ENDIF
!
      unfl = SLAMCH('Safe minimum')
      ulp = SLAMCH('Epsilon')*SLAMCH('Base')
!
!     Some Error Checks
!
!     Do Test 1
!
!     Norm of A:
!
      anorm = MAX(SLANSB('1',cuplo,N,ika,A,Lda,Work),unfl)
!
!     Compute error matrix:    Error = A - U S U**T
!
!     Copy A from SB to SP storage format.
!
      j = 0
      DO jc = 1 , N
         IF ( lower ) THEN
            DO jr = 1 , MIN(ika+1,N+1-jc)
               j = j + 1
               Work(j) = A(jr,jc)
            ENDDO
            DO jr = ika + 2 , N + 1 - jc
               j = j + 1
               Work(j) = ZERO
            ENDDO
         ELSE
            DO jr = ika + 2 , jc
               j = j + 1
               Work(j) = ZERO
            ENDDO
            DO jr = MIN(ika,jc-1) , 0 , -1
               j = j + 1
               Work(j) = A(ika+1-jr,jc)
            ENDDO
         ENDIF
      ENDDO
!
      DO j = 1 , N
         CALL SSPR(cuplo,N,-D(j),U(1,j),1,Work)
      ENDDO
!
      IF ( N>1 .AND. Ks==1 ) THEN
         DO j = 1 , N - 1
            CALL SSPR2(cuplo,N,-E(j),U(1,j),1,U(1,j+1),1,Work)
         ENDDO
      ENDIF
      wnorm = SLANSP('1',cuplo,N,Work,Work(lw+1))
!
      IF ( anorm>wnorm ) THEN
         Result(1) = (wnorm/anorm)/(N*ulp)
      ELSEIF ( anorm<ONE ) THEN
         Result(1) = (MIN(wnorm,N*anorm)/anorm)/(N*ulp)
      ELSE
         Result(1) = MIN(wnorm/anorm,REAL(N))/(N*ulp)
      ENDIF
!
!     Do Test 2
!
!     Compute  U U**T - I
!
      CALL SGEMM('N','C',N,N,N,ONE,U,Ldu,U,Ldu,ZERO,Work,N)
!
      DO j = 1 , N
         Work((N+1)*(j-1)+1) = Work((N+1)*(j-1)+1) - ONE
      ENDDO
!
      Result(2) = MIN(SLANGE('1',N,N,Work,N,Work(N**2+1)),REAL(N))      &
     &            /(N*ulp)
!
!
!     End of SSBT21
!
      END SUBROUTINE SSBT21
