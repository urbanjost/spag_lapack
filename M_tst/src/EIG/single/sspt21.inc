!*==sspt21.f90  processed by SPAG 7.51RB at 17:37 on  4 Mar 2022
!> \brief \b sspt21
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE sspt21( ITYPE, UPLO, N, KBAND, AP, D, E, U, LDU, VP,
!                          TAU, WORK, RESULT )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            ITYPE, KBAND, LDU, N
!       ..
!       .. Array Arguments ..
!       REAL               AP( * ), D( * ), E( * ), RESULT( 2 ), TAU( * ),
!      $                   U( LDU, * ), VP( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> sspt21  generally checks a decomposition of the form
!>
!>         A = U S U**T
!>
!> where **T means transpose, A is symmetric (stored in packed format), U
!> is orthogonal, and S is diagonal (if KBAND=0) or symmetric
!> tridiagonal (if KBAND=1).  If ITYPE=1, then U is represented as a
!> dense matrix, otherwise the U is expressed as a product of
!> Householder transformations, whose vectors are stored in the array
!> "V" and whose scaling constants are in "TAU"; we shall use the
!> letter "V" to refer to the product of Householder transformations
!> (which should be equal to U).
!>
!> Specifically, if ITYPE=1, then:
!>
!>         RESULT(1) = | A - U S U**T | / ( |A| n ulp ) and
!>         RESULT(2) = | I - U U**T | / ( n ulp )
!>
!> If ITYPE=2, then:
!>
!>         RESULT(1) = | A - V S V**T | / ( |A| n ulp )
!>
!> If ITYPE=3, then:
!>
!>         RESULT(1) = | I - V U**T | / ( n ulp )
!>
!> Packed storage means that, for example, if UPLO='U', then the columns
!> of the upper triangle of A are stored one after another, so that
!> A(1,j+1) immediately follows A(j,j) in the array AP.  Similarly, if
!> UPLO='L', then the columns of the lower triangle of A are stored one
!> after another in AP, so that A(j+1,j+1) immediately follows A(n,j)
!> in the array AP.  This means that A(i,j) is stored in:
!>
!>    AP( i + j*(j-1)/2 )                 if UPLO='U'
!>
!>    AP( i + (2*n-j)*(j-1)/2 )           if UPLO='L'
!>
!> The array VP bears the same relation to the matrix V that A does to
!> AP.
!>
!> For ITYPE > 1, the transformation U is expressed as a product
!> of Householder transformations:
!>
!>    If UPLO='U', then  V = H(n-1)...H(1),  where
!>
!>        H(j) = I  -  tau(j) v(j) v(j)**T
!>
!>    and the first j-1 elements of v(j) are stored in V(1:j-1,j+1),
!>    (i.e., VP( j*(j+1)/2 + 1 : j*(j+1)/2 + j-1 ) ),
!>    the j-th element is 1, and the last n-j elements are 0.
!>
!>    If UPLO='L', then  V = H(1)...H(n-1),  where
!>
!>        H(j) = I  -  tau(j) v(j) v(j)**T
!>
!>    and the first j elements of v(j) are 0, the (j+1)-st is 1, and the
!>    (j+2)-nd through n-th elements are stored in V(j+2:n,j) (i.e.,
!>    in VP( (2*n-j)*(j-1)/2 + j+2 : (2*n-j)*(j-1)/2 + n ) .)
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] ITYPE
!> \verbatim
!>          ITYPE is INTEGER
!>          Specifies the type of tests to be performed.
!>          1: U expressed as a dense orthogonal matrix:
!>             RESULT(1) = | A - U S U**T | / ( |A| n ulp ) and
!>             RESULT(2) = | I - U U**T | / ( n ulp )
!>
!>          2: U expressed as a product V of Housholder transformations:
!>             RESULT(1) = | A - V S V**T | / ( |A| n ulp )
!>
!>          3: U expressed both as a dense orthogonal matrix and
!>             as a product of Housholder transformations:
!>             RESULT(1) = | I - V U**T | / ( n ulp )
!> \endverbatim
!>
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER
!>          If UPLO='U', AP and VP are considered to contain the upper
!>          triangle of A and V.
!>          If UPLO='L', AP and VP are considered to contain the lower
!>          triangle of A and V.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The size of the matrix.  If it is zero, sspt21 does nothing.
!>          It must be at least zero.
!> \endverbatim
!>
!> \param[in] KBAND
!> \verbatim
!>          KBAND is INTEGER
!>          The bandwidth of the matrix.  It may only be zero or one.
!>          If zero, then S is diagonal, and E is not referenced.  If
!>          one, then S is symmetric tri-diagonal.
!> \endverbatim
!>
!> \param[in] AP
!> \verbatim
!>          AP is REAL array, dimension (N*(N+1)/2)
!>          The original (unfactored) matrix.  It is assumed to be
!>          symmetric, and contains the columns of just the upper
!>          triangle (UPLO='U') or only the lower triangle (UPLO='L'),
!>          packed one after another.
!> \endverbatim
!>
!> \param[in] D
!> \verbatim
!>          D is REAL array, dimension (N)
!>          The diagonal of the (symmetric tri-) diagonal matrix.
!> \endverbatim
!>
!> \param[in] E
!> \verbatim
!>          E is REAL array, dimension (N-1)
!>          The off-diagonal of the (symmetric tri-) diagonal matrix.
!>          E(1) is the (1,2) and (2,1) element, E(2) is the (2,3) and
!>          (3,2) element, etc.
!>          Not referenced if KBAND=0.
!> \endverbatim
!>
!> \param[in] U
!> \verbatim
!>          U is REAL array, dimension (LDU, N)
!>          If ITYPE=1 or 3, this contains the orthogonal matrix in
!>          the decomposition, expressed as a dense matrix.  If ITYPE=2,
!>          then it is not referenced.
!> \endverbatim
!>
!> \param[in] LDU
!> \verbatim
!>          LDU is INTEGER
!>          The leading dimension of U.  LDU must be at least N and
!>          at least 1.
!> \endverbatim
!>
!> \param[in] VP
!> \verbatim
!>          VP is REAL array, dimension (N*(N+1)/2)
!>          If ITYPE=2 or 3, the columns of this array contain the
!>          Householder vectors used to describe the orthogonal matrix
!>          in the decomposition, as described in purpose.
!>          *NOTE* If ITYPE=2 or 3, V is modified and restored.  The
!>          subdiagonal (if UPLO='L') or the superdiagonal (if UPLO='U')
!>          is set to one, and later reset to its original value, during
!>          the course of the calculation.
!>          If ITYPE=1, then it is neither referenced nor modified.
!> \endverbatim
!>
!> \param[in] TAU
!> \verbatim
!>          TAU is REAL array, dimension (N)
!>          If ITYPE >= 2, then TAU(j) is the scalar factor of
!>          v(j) v(j)**T in the Householder transformation H(j) of
!>          the product  U = H(1)...H(n-2)
!>          If ITYPE < 2, then TAU is not referenced.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is REAL array, dimension (N**2+N)
!>          Workspace.
!> \endverbatim
!>
!> \param[out] RESULT
!> \verbatim
!>          RESULT is REAL array, dimension (2)
!>          The values computed by the two tests described above.  The
!>          values are currently limited to 1/ulp, to avoid overflow.
!>          RESULT(1) is always modified.  RESULT(2) is modified only
!>          if ITYPE=1.
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
      SUBROUTINE sspt21(Itype,Uplo,N,Kband,Ap,D,E,U,Ldu,Vp,Tau,Work,    &
     &                  Result)
      IMPLICIT NONE
!*--SSPT21225
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Uplo
      INTEGER Itype , Kband , Ldu , N
!     ..
!     .. Array Arguments ..
      REAL Ap(*) , D(*) , E(*) , Result(2) , Tau(*) , U(Ldu,*) , Vp(*) ,&
     &     Work(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ZERO , ONE , TEN
      PARAMETER (ZERO=0.0E0,ONE=1.0E0,TEN=10.0E0)
      REAL HALF
      PARAMETER (HALF=1.0E+0/2.0E+0)
!     ..
!     .. Local Scalars ..
      LOGICAL lower
      CHARACTER cuplo
      INTEGER iinfo , j , jp , jp1 , jr , lap
      REAL anorm , temp , ulp , unfl , vsave , wnorm
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      REAL SDOT , SLAMCH , SLANGE , SLANSP
      EXTERNAL LSAME , SDOT , SLAMCH , SLANGE , SLANSP
!     ..
!     .. External Subroutines ..
      EXTERNAL SAXPY , SCOPY , SGEMM , SLACPY , SLASET , SOPMTR ,   SSPMV , SSPR , SSPR2
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX , MIN , REAL
!     ..
!     .. Executable Statements ..
!
!     1)      Constants
!
      Result(1) = ZERO
      IF ( Itype==1 ) Result(2) = ZERO
      IF ( N<=0 ) RETURN
!
      lap = (N*(N+1))/2
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
      IF ( Itype<1 .OR. Itype>3 ) THEN
         Result(1) = TEN/ulp
         RETURN
      ENDIF
!
!     Do Test 1
!
!     Norm of A:
!
      IF ( Itype==3 ) THEN
         anorm = ONE
      ELSE
         anorm = MAX(SLANSP('1',cuplo,N,Ap,Work),unfl)
      ENDIF
!
!     Compute error matrix:
!
      IF ( Itype==1 ) THEN
!
!        ITYPE=1: error = A - U S U**T
!
         CALL SLASET('Full',N,N,ZERO,ZERO,Work,N)
         CALL SCOPY(lap,Ap,1,Work,1)
!
         DO j = 1 , N
            CALL SSPR(cuplo,N,-D(j),U(1,j),1,Work)
         ENDDO
!
         IF ( N>1 .AND. Kband==1 ) THEN
            DO j = 1 , N - 1
               CALL SSPR2(cuplo,N,-E(j),U(1,j),1,U(1,j+1),1,Work)
            ENDDO
         ENDIF
         wnorm = SLANSP('1',cuplo,N,Work,Work(N**2+1))
!
      ELSEIF ( Itype==2 ) THEN
!
!        ITYPE=2: error = V S V**T - A
!
         CALL SLASET('Full',N,N,ZERO,ZERO,Work,N)
!
         IF ( lower ) THEN
            Work(lap) = D(N)
            DO j = N - 1 , 1 , -1
               jp = ((2*N-j)*(j-1))/2
               jp1 = jp + N - j
               IF ( Kband==1 ) THEN
                  Work(jp+j+1) = (ONE-Tau(j))*E(j)
                  DO jr = j + 2 , N
                     Work(jp+jr) = -Tau(j)*E(j)*Vp(jp+jr)
                  ENDDO
               ENDIF
!
               IF ( Tau(j)/=ZERO ) THEN
                  vsave = Vp(jp+j+1)
                  Vp(jp+j+1) = ONE
                  CALL SSPMV('L',N-j,ONE,Work(jp1+j+1),Vp(jp+j+1),1,    &
     &                       ZERO,Work(lap+1),1)
                  temp = -HALF*Tau(j)                                   &
     &                   *SDOT(N-j,Work(lap+1),1,Vp(jp+j+1),1)
                  CALL SAXPY(N-j,temp,Vp(jp+j+1),1,Work(lap+1),1)
                  CALL SSPR2('L',N-j,-Tau(j),Vp(jp+j+1),1,Work(lap+1),1,&
     &                       Work(jp1+j+1))
                  Vp(jp+j+1) = vsave
               ENDIF
               Work(jp+j) = D(j)
            ENDDO
         ELSE
            Work(1) = D(1)
            DO j = 1 , N - 1
               jp = (j*(j-1))/2
               jp1 = jp + j
               IF ( Kband==1 ) THEN
                  Work(jp1+j) = (ONE-Tau(j))*E(j)
                  DO jr = 1 , j - 1
                     Work(jp1+jr) = -Tau(j)*E(j)*Vp(jp1+jr)
                  ENDDO
               ENDIF
!
               IF ( Tau(j)/=ZERO ) THEN
                  vsave = Vp(jp1+j)
                  Vp(jp1+j) = ONE
                  CALL SSPMV('U',j,ONE,Work,Vp(jp1+1),1,ZERO,Work(lap+1)&
     &                       ,1)
                  temp = -HALF*Tau(j)*SDOT(j,Work(lap+1),1,Vp(jp1+1),1)
                  CALL SAXPY(j,temp,Vp(jp1+1),1,Work(lap+1),1)
                  CALL SSPR2('U',j,-Tau(j),Vp(jp1+1),1,Work(lap+1),1,   &
     &                       Work)
                  Vp(jp1+j) = vsave
               ENDIF
               Work(jp1+j+1) = D(j+1)
            ENDDO
         ENDIF
!
         DO j = 1 , lap
            Work(j) = Work(j) - Ap(j)
         ENDDO
         wnorm = SLANSP('1',cuplo,N,Work,Work(lap+1))
!
      ELSEIF ( Itype==3 ) THEN
!
!        ITYPE=3: error = U V**T - I
!
         IF ( N<2 ) RETURN
         CALL SLACPY(' ',N,N,U,Ldu,Work,N)
         CALL SOPMTR('R',cuplo,'T',N,N,Vp,Tau,Work,N,Work(N**2+1),iinfo)
         IF ( iinfo/=0 ) THEN
            Result(1) = TEN/ulp
            RETURN
         ENDIF
!
         DO j = 1 , N
            Work((N+1)*(j-1)+1) = Work((N+1)*(j-1)+1) - ONE
         ENDDO
!
         wnorm = SLANGE('1',N,N,Work,N,Work(N**2+1))
      ENDIF
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
      IF ( Itype==1 ) THEN
         CALL SGEMM('N','C',N,N,N,ONE,U,Ldu,U,Ldu,ZERO,Work,N)
!
         DO j = 1 , N
            Work((N+1)*(j-1)+1) = Work((N+1)*(j-1)+1) - ONE
         ENDDO
!
         Result(2) = MIN(SLANGE('1',N,N,Work,N,Work(N**2+1)),REAL(N))   &
     &               /(N*ulp)
      ENDIF
!
!
!     End of sspt21
!
      END SUBROUTINE sspt21
