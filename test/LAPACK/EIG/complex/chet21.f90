!*==chet21.f90  processed by SPAG 7.51RB at 17:37 on  4 Mar 2022
!> \brief \b CHET21
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE CHET21( ITYPE, UPLO, N, KBAND, A, LDA, D, E, U, LDU, V,
!                          LDV, TAU, WORK, RWORK, RESULT )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            ITYPE, KBAND, LDA, LDU, LDV, N
!       ..
!       .. Array Arguments ..
!       REAL               D( * ), E( * ), RESULT( 2 ), RWORK( * )
!       COMPLEX            A( LDA, * ), TAU( * ), U( LDU, * ),
!      $                   V( LDV, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CHET21 generally checks a decomposition of the form
!>
!>    A = U S U**H
!>
!> where **H means conjugate transpose, A is hermitian, U is unitary, and
!> S is diagonal (if KBAND=0) or (real) symmetric tridiagonal (if
!> KBAND=1).
!>
!> If ITYPE=1, then U is represented as a dense matrix; otherwise U is
!> expressed as a product of Householder transformations, whose vectors
!> are stored in the array "V" and whose scaling constants are in "TAU".
!> We shall use the letter "V" to refer to the product of Householder
!> transformations (which should be equal to U).
!>
!> Specifically, if ITYPE=1, then:
!>
!>    RESULT(1) = | A - U S U**H | / ( |A| n ulp ) and
!>    RESULT(2) = | I - U U**H | / ( n ulp )
!>
!> If ITYPE=2, then:
!>
!>    RESULT(1) = | A - V S V**H | / ( |A| n ulp )
!>
!> If ITYPE=3, then:
!>
!>    RESULT(1) = | I - U V**H | / ( n ulp )
!>
!> For ITYPE > 1, the transformation U is expressed as a product
!> V = H(1)...H(n-2),  where H(j) = I  -  tau(j) v(j) v(j)**H and each
!> vector v(j) has its first j elements 0 and the remaining n-j elements
!> stored in V(j+1:n,j).
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] ITYPE
!> \verbatim
!>          ITYPE is INTEGER
!>          Specifies the type of tests to be performed.
!>          1: U expressed as a dense unitary matrix:
!>             RESULT(1) = | A - U S U**H | / ( |A| n ulp ) and
!>             RESULT(2) = | I - U U**H | / ( n ulp )
!>
!>          2: U expressed as a product V of Housholder transformations:
!>             RESULT(1) = | A - V S V**H | / ( |A| n ulp )
!>
!>          3: U expressed both as a dense unitary matrix and
!>             as a product of Housholder transformations:
!>             RESULT(1) = | I - U V**H | / ( n ulp )
!> \endverbatim
!>
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
!>          The size of the matrix.  If it is zero, CHET21 does nothing.
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
!> \param[in] A
!> \verbatim
!>          A is COMPLEX array, dimension (LDA, N)
!>          The original (unfactored) matrix.  It is assumed to be
!>          hermitian, and only the upper (UPLO='U') or only the lower
!>          (UPLO='L') will be referenced.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of A.  It must be at least 1
!>          and at least N.
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
!>          U is COMPLEX array, dimension (LDU, N)
!>          If ITYPE=1 or 3, this contains the unitary matrix in
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
!> \param[in] V
!> \verbatim
!>          V is COMPLEX array, dimension (LDV, N)
!>          If ITYPE=2 or 3, the columns of this array contain the
!>          Householder vectors used to describe the unitary matrix
!>          in the decomposition.  If UPLO='L', then the vectors are in
!>          the lower triangle, if UPLO='U', then in the upper
!>          triangle.
!>          *NOTE* If ITYPE=2 or 3, V is modified and restored.  The
!>          subdiagonal (if UPLO='L') or the superdiagonal (if UPLO='U')
!>          is set to one, and later reset to its original value, during
!>          the course of the calculation.
!>          If ITYPE=1, then it is neither referenced nor modified.
!> \endverbatim
!>
!> \param[in] LDV
!> \verbatim
!>          LDV is INTEGER
!>          The leading dimension of V.  LDV must be at least N and
!>          at least 1.
!> \endverbatim
!>
!> \param[in] TAU
!> \verbatim
!>          TAU is COMPLEX array, dimension (N)
!>          If ITYPE >= 2, then TAU(j) is the scalar factor of
!>          v(j) v(j)**H in the Householder transformation H(j) of
!>          the product  U = H(1)...H(n-2)
!>          If ITYPE < 2, then TAU is not referenced.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX array, dimension (2*N**2)
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is REAL array, dimension (N)
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
!> \ingroup complex_eig
!
!  =====================================================================
      SUBROUTINE CHET21(Itype,Uplo,N,Kband,A,Lda,D,E,U,Ldu,V,Ldv,Tau,   &
     &                  Work,Rwork,Result)
      IMPLICIT NONE
!*--CHET21218
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Uplo
      INTEGER Itype , Kband , Lda , Ldu , Ldv , N
!     ..
!     .. Array Arguments ..
      REAL D(*) , E(*) , Result(2) , Rwork(*)
      COMPLEX A(Lda,*) , Tau(*) , U(Ldu,*) , V(Ldv,*) , Work(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ZERO , ONE , TEN
      PARAMETER (ZERO=0.0E+0,ONE=1.0E+0,TEN=10.0E+0)
      COMPLEX CZERO , CONE
      PARAMETER (CZERO=(0.0E+0,0.0E+0),CONE=(1.0E+0,0.0E+0))
!     ..
!     .. Local Scalars ..
      LOGICAL lower
      CHARACTER cuplo
      INTEGER iinfo , j , jcol , jr , jrow
      REAL anorm , ulp , unfl , wnorm
      COMPLEX vsave
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      REAL CLANGE , CLANHE , SLAMCH
      EXTERNAL LSAME , CLANGE , CLANHE , SLAMCH
!     ..
!     .. External Subroutines ..
      EXTERNAL CGEMM , CHER , CHER2 , CLACPY , CLARFY , CLASET ,        &
     &         CUNM2L , CUNM2R
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC CMPLX , MAX , MIN , REAL
!     ..
!     .. Executable Statements ..
!
      Result(1) = ZERO
      IF ( Itype==1 ) Result(2) = ZERO
      IF ( N<=0 ) RETURN
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
         anorm = MAX(CLANHE('1',cuplo,N,A,Lda,Rwork),unfl)
      ENDIF
!
!     Compute error matrix:
!
      IF ( Itype==1 ) THEN
!
!        ITYPE=1: error = A - U S U**H
!
         CALL CLASET('Full',N,N,CZERO,CZERO,Work,N)
         CALL CLACPY(cuplo,N,N,A,Lda,Work,N)
!
         DO j = 1 , N
            CALL CHER(cuplo,N,-D(j),U(1,j),1,Work,N)
         ENDDO
!
         IF ( N>1 .AND. Kband==1 ) THEN
            DO j = 2 , N - 1
               CALL CHER2(cuplo,N,-CMPLX(E(j)),U(1,j),1,U(1,j-1),1,Work,&
     &                    N)
            ENDDO
         ENDIF
         wnorm = CLANHE('1',cuplo,N,Work,N,Rwork)
!
      ELSEIF ( Itype==2 ) THEN
!
!        ITYPE=2: error = V S V**H - A
!
         CALL CLASET('Full',N,N,CZERO,CZERO,Work,N)
!
         IF ( lower ) THEN
            Work(N**2) = D(N)
            DO j = N - 1 , 1 , -1
               IF ( Kband==1 ) THEN
                  Work((N+1)*(j-1)+2) = (CONE-Tau(j))*E(j)
                  DO jr = j + 2 , N
                     Work((j-1)*N+jr) = -Tau(j)*E(j)*V(jr,j)
                  ENDDO
               ENDIF
!
               vsave = V(j+1,j)
               V(j+1,j) = ONE
               CALL CLARFY('L',N-j,V(j+1,j),1,Tau(j),Work((N+1)*j+1),N, &
     &                     Work(N**2+1))
               V(j+1,j) = vsave
               Work((N+1)*(j-1)+1) = D(j)
            ENDDO
         ELSE
            Work(1) = D(1)
            DO j = 1 , N - 1
               IF ( Kband==1 ) THEN
                  Work((N+1)*j) = (CONE-Tau(j))*E(j)
                  DO jr = 1 , j - 1
                     Work(j*N+jr) = -Tau(j)*E(j)*V(jr,j+1)
                  ENDDO
               ENDIF
!
               vsave = V(j,j+1)
               V(j,j+1) = ONE
               CALL CLARFY('U',j,V(1,j+1),1,Tau(j),Work,N,Work(N**2+1))
               V(j,j+1) = vsave
               Work((N+1)*j+1) = D(j+1)
            ENDDO
         ENDIF
!
         DO jcol = 1 , N
            IF ( lower ) THEN
               DO jrow = jcol , N
                  Work(jrow+N*(jcol-1)) = Work(jrow+N*(jcol-1))         &
     &               - A(jrow,jcol)
               ENDDO
            ELSE
               DO jrow = 1 , jcol
                  Work(jrow+N*(jcol-1)) = Work(jrow+N*(jcol-1))         &
     &               - A(jrow,jcol)
               ENDDO
            ENDIF
         ENDDO
         wnorm = CLANHE('1',cuplo,N,Work,N,Rwork)
!
      ELSEIF ( Itype==3 ) THEN
!
!        ITYPE=3: error = U V**H - I
!
         IF ( N<2 ) RETURN
         CALL CLACPY(' ',N,N,U,Ldu,Work,N)
         IF ( lower ) THEN
            CALL CUNM2R('R','C',N,N-1,N-1,V(2,1),Ldv,Tau,Work(N+1),N,   &
     &                  Work(N**2+1),iinfo)
         ELSE
            CALL CUNM2L('R','C',N,N-1,N-1,V(1,2),Ldv,Tau,Work,N,        &
     &                  Work(N**2+1),iinfo)
         ENDIF
         IF ( iinfo/=0 ) THEN
            Result(1) = TEN/ulp
            RETURN
         ENDIF
!
         DO j = 1 , N
            Work((N+1)*(j-1)+1) = Work((N+1)*(j-1)+1) - CONE
         ENDDO
!
         wnorm = CLANGE('1',N,N,Work,N,Rwork)
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
!     Compute  U U**H - I
!
      IF ( Itype==1 ) THEN
         CALL CGEMM('N','C',N,N,N,CONE,U,Ldu,U,Ldu,CZERO,Work,N)
!
         DO j = 1 , N
            Work((N+1)*(j-1)+1) = Work((N+1)*(j-1)+1) - CONE
         ENDDO
!
         Result(2) = MIN(CLANGE('1',N,N,Work,N,Rwork),REAL(N))/(N*ulp)
      ENDIF
!
!
!     End of CHET21
!
      END SUBROUTINE CHET21
