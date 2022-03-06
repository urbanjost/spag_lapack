!*==dlatm5.f90  processed by SPAG 7.51RB at 17:44 on  4 Mar 2022
!> \brief \b DLATM5
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLATM5( PRTYPE, M, N, A, LDA, B, LDB, C, LDC, D, LDD,
!                          E, LDE, F, LDF, R, LDR, L, LDL, ALPHA, QBLCKA,
!                          QBLCKB )
!
!       .. Scalar Arguments ..
!       INTEGER            LDA, LDB, LDC, LDD, LDE, LDF, LDL, LDR, M, N,
!      $                   PRTYPE, QBLCKA, QBLCKB
!       DOUBLE PRECISION   ALPHA
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), C( LDC, * ),
!      $                   D( LDD, * ), E( LDE, * ), F( LDF, * ),
!      $                   L( LDL, * ), R( LDR, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLATM5 generates matrices involved in the Generalized Sylvester
!> equation:
!>
!>     A * R - L * B = C
!>     D * R - L * E = F
!>
!> They also satisfy (the diagonalization condition)
!>
!>  [ I -L ] ( [ A  -C ], [ D -F ] ) [ I  R ] = ( [ A    ], [ D    ] )
!>  [    I ] ( [     B ]  [    E ] ) [    I ]   ( [    B ]  [    E ] )
!>
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] PRTYPE
!> \verbatim
!>          PRTYPE is INTEGER
!>          "Points" to a certain type of the matrices to generate
!>          (see further details).
!> \endverbatim
!>
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          Specifies the order of A and D and the number of rows in
!>          C, F,  R and L.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          Specifies the order of B and E and the number of columns in
!>          C, F, R and L.
!> \endverbatim
!>
!> \param[out] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (LDA, M).
!>          On exit A M-by-M is initialized according to PRTYPE.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of A.
!> \endverbatim
!>
!> \param[out] B
!> \verbatim
!>          B is DOUBLE PRECISION array, dimension (LDB, N).
!>          On exit B N-by-N is initialized according to PRTYPE.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of B.
!> \endverbatim
!>
!> \param[out] C
!> \verbatim
!>          C is DOUBLE PRECISION array, dimension (LDC, N).
!>          On exit C M-by-N is initialized according to PRTYPE.
!> \endverbatim
!>
!> \param[in] LDC
!> \verbatim
!>          LDC is INTEGER
!>          The leading dimension of C.
!> \endverbatim
!>
!> \param[out] D
!> \verbatim
!>          D is DOUBLE PRECISION array, dimension (LDD, M).
!>          On exit D M-by-M is initialized according to PRTYPE.
!> \endverbatim
!>
!> \param[in] LDD
!> \verbatim
!>          LDD is INTEGER
!>          The leading dimension of D.
!> \endverbatim
!>
!> \param[out] E
!> \verbatim
!>          E is DOUBLE PRECISION array, dimension (LDE, N).
!>          On exit E N-by-N is initialized according to PRTYPE.
!> \endverbatim
!>
!> \param[in] LDE
!> \verbatim
!>          LDE is INTEGER
!>          The leading dimension of E.
!> \endverbatim
!>
!> \param[out] F
!> \verbatim
!>          F is DOUBLE PRECISION array, dimension (LDF, N).
!>          On exit F M-by-N is initialized according to PRTYPE.
!> \endverbatim
!>
!> \param[in] LDF
!> \verbatim
!>          LDF is INTEGER
!>          The leading dimension of F.
!> \endverbatim
!>
!> \param[out] R
!> \verbatim
!>          R is DOUBLE PRECISION array, dimension (LDR, N).
!>          On exit R M-by-N is initialized according to PRTYPE.
!> \endverbatim
!>
!> \param[in] LDR
!> \verbatim
!>          LDR is INTEGER
!>          The leading dimension of R.
!> \endverbatim
!>
!> \param[out] L
!> \verbatim
!>          L is DOUBLE PRECISION array, dimension (LDL, N).
!>          On exit L M-by-N is initialized according to PRTYPE.
!> \endverbatim
!>
!> \param[in] LDL
!> \verbatim
!>          LDL is INTEGER
!>          The leading dimension of L.
!> \endverbatim
!>
!> \param[in] ALPHA
!> \verbatim
!>          ALPHA is DOUBLE PRECISION
!>          Parameter used in generating PRTYPE = 1 and 5 matrices.
!> \endverbatim
!>
!> \param[in] QBLCKA
!> \verbatim
!>          QBLCKA is INTEGER
!>          When PRTYPE = 3, specifies the distance between 2-by-2
!>          blocks on the diagonal in A. Otherwise, QBLCKA is not
!>          referenced. QBLCKA > 1.
!> \endverbatim
!>
!> \param[in] QBLCKB
!> \verbatim
!>          QBLCKB is INTEGER
!>          When PRTYPE = 3, specifies the distance between 2-by-2
!>          blocks on the diagonal in B. Otherwise, QBLCKB is not
!>          referenced. QBLCKB > 1.
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
!> \date June 2016
!
!> \ingroup double_matgen
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  PRTYPE = 1: A and B are Jordan blocks, D and E are identity matrices
!>
!>             A : if (i == j) then A(i, j) = 1.0
!>                 if (j == i + 1) then A(i, j) = -1.0
!>                 else A(i, j) = 0.0,            i, j = 1...M
!>
!>             B : if (i == j) then B(i, j) = 1.0 - ALPHA
!>                 if (j == i + 1) then B(i, j) = 1.0
!>                 else B(i, j) = 0.0,            i, j = 1...N
!>
!>             D : if (i == j) then D(i, j) = 1.0
!>                 else D(i, j) = 0.0,            i, j = 1...M
!>
!>             E : if (i == j) then E(i, j) = 1.0
!>                 else E(i, j) = 0.0,            i, j = 1...N
!>
!>             L =  R are chosen from [-10...10],
!>                  which specifies the right hand sides (C, F).
!>
!>  PRTYPE = 2 or 3: Triangular and/or quasi- triangular.
!>
!>             A : if (i <= j) then A(i, j) = [-1...1]
!>                 else A(i, j) = 0.0,             i, j = 1...M
!>
!>                 if (PRTYPE = 3) then
!>                    A(k + 1, k + 1) = A(k, k)
!>                    A(k + 1, k) = [-1...1]
!>                    sign(A(k, k + 1) = -(sin(A(k + 1, k))
!>                        k = 1, M - 1, QBLCKA
!>
!>             B : if (i <= j) then B(i, j) = [-1...1]
!>                 else B(i, j) = 0.0,            i, j = 1...N
!>
!>                 if (PRTYPE = 3) then
!>                    B(k + 1, k + 1) = B(k, k)
!>                    B(k + 1, k) = [-1...1]
!>                    sign(B(k, k + 1) = -(sign(B(k + 1, k))
!>                        k = 1, N - 1, QBLCKB
!>
!>             D : if (i <= j) then D(i, j) = [-1...1].
!>                 else D(i, j) = 0.0,            i, j = 1...M
!>
!>
!>             E : if (i <= j) then D(i, j) = [-1...1]
!>                 else E(i, j) = 0.0,            i, j = 1...N
!>
!>                 L, R are chosen from [-10...10],
!>                 which specifies the right hand sides (C, F).
!>
!>  PRTYPE = 4 Full
!>             A(i, j) = [-10...10]
!>             D(i, j) = [-1...1]    i,j = 1...M
!>             B(i, j) = [-10...10]
!>             E(i, j) = [-1...1]    i,j = 1...N
!>             R(i, j) = [-10...10]
!>             L(i, j) = [-1...1]    i = 1..M ,j = 1...N
!>
!>             L, R specifies the right hand sides (C, F).
!>
!>  PRTYPE = 5 special case common and/or close eigs.
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE DLATM5(Prtype,M,N,A,Lda,B,Ldb,C,Ldc,D,Ldd,E,Lde,F,Ldf, &
     &                  R,Ldr,L,Ldl,Alpha,Qblcka,Qblckb)
      IMPLICIT NONE
!*--DLATM5271
!
!  -- LAPACK computational routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     June 2016
!
!     .. Scalar Arguments ..
      INTEGER Lda , Ldb , Ldc , Ldd , Lde , Ldf , Ldl , Ldr , M , N ,   &
     &        Prtype , Qblcka , Qblckb
      DOUBLE PRECISION Alpha
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION A(Lda,*) , B(Ldb,*) , C(Ldc,*) , D(Ldd,*) ,      &
     &                 E(Lde,*) , F(Ldf,*) , L(Ldl,*) , R(Ldr,*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ONE , ZERO , TWENTY , HALF , TWO
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0,TWENTY=2.0D+1,HALF=0.5D+0,      &
     &           TWO=2.0D+0)
!     ..
!     .. Local Scalars ..
      INTEGER i , j , k
      DOUBLE PRECISION imeps , reeps
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC DBLE , MOD , SIN
!     ..
!     .. External Subroutines ..
      EXTERNAL DGEMM
!     ..
!     .. Executable Statements ..
!
      IF ( Prtype==1 ) THEN
         DO i = 1 , M
            DO j = 1 , M
               IF ( i==j ) THEN
                  A(i,j) = ONE
                  D(i,j) = ONE
               ELSEIF ( i==j-1 ) THEN
                  A(i,j) = -ONE
                  D(i,j) = ZERO
               ELSE
                  A(i,j) = ZERO
                  D(i,j) = ZERO
               ENDIF
            ENDDO
         ENDDO
!
         DO i = 1 , N
            DO j = 1 , N
               IF ( i==j ) THEN
                  B(i,j) = ONE - Alpha
                  E(i,j) = ONE
               ELSEIF ( i==j-1 ) THEN
                  B(i,j) = ONE
                  E(i,j) = ZERO
               ELSE
                  B(i,j) = ZERO
                  E(i,j) = ZERO
               ENDIF
            ENDDO
         ENDDO
!
         DO i = 1 , M
            DO j = 1 , N
               R(i,j) = (HALF-SIN(DBLE(i/j)))*TWENTY
               L(i,j) = R(i,j)
            ENDDO
         ENDDO
!
      ELSEIF ( Prtype==2 .OR. Prtype==3 ) THEN
         DO i = 1 , M
            DO j = 1 , M
               IF ( i<=j ) THEN
                  A(i,j) = (HALF-SIN(DBLE(i)))*TWO
                  D(i,j) = (HALF-SIN(DBLE(i*j)))*TWO
               ELSE
                  A(i,j) = ZERO
                  D(i,j) = ZERO
               ENDIF
            ENDDO
         ENDDO
!
         DO i = 1 , N
            DO j = 1 , N
               IF ( i<=j ) THEN
                  B(i,j) = (HALF-SIN(DBLE(i+j)))*TWO
                  E(i,j) = (HALF-SIN(DBLE(j)))*TWO
               ELSE
                  B(i,j) = ZERO
                  E(i,j) = ZERO
               ENDIF
            ENDDO
         ENDDO
!
         DO i = 1 , M
            DO j = 1 , N
               R(i,j) = (HALF-SIN(DBLE(i*j)))*TWENTY
               L(i,j) = (HALF-SIN(DBLE(i+j)))*TWENTY
            ENDDO
         ENDDO
!
         IF ( Prtype==3 ) THEN
            IF ( Qblcka<=1 ) Qblcka = 2
            DO k = 1 , M - 1 , Qblcka
               A(k+1,k+1) = A(k,k)
               A(k+1,k) = -SIN(A(k,k+1))
            ENDDO
!
            IF ( Qblckb<=1 ) Qblckb = 2
            DO k = 1 , N - 1 , Qblckb
               B(k+1,k+1) = B(k,k)
               B(k+1,k) = -SIN(B(k,k+1))
            ENDDO
         ENDIF
!
      ELSEIF ( Prtype==4 ) THEN
         DO i = 1 , M
            DO j = 1 , M
               A(i,j) = (HALF-SIN(DBLE(i*j)))*TWENTY
               D(i,j) = (HALF-SIN(DBLE(i+j)))*TWO
            ENDDO
         ENDDO
!
         DO i = 1 , N
            DO j = 1 , N
               B(i,j) = (HALF-SIN(DBLE(i+j)))*TWENTY
               E(i,j) = (HALF-SIN(DBLE(i*j)))*TWO
            ENDDO
         ENDDO
!
         DO i = 1 , M
            DO j = 1 , N
               R(i,j) = (HALF-SIN(DBLE(j/i)))*TWENTY
               L(i,j) = (HALF-SIN(DBLE(i*j)))*TWO
            ENDDO
         ENDDO
!
      ELSEIF ( Prtype>=5 ) THEN
         reeps = HALF*TWO*TWENTY/Alpha
         imeps = (HALF-TWO)/Alpha
         DO i = 1 , M
            DO j = 1 , N
               R(i,j) = (HALF-SIN(DBLE(i*j)))*Alpha/TWENTY
               L(i,j) = (HALF-SIN(DBLE(i+j)))*Alpha/TWENTY
            ENDDO
         ENDDO
!
         DO i = 1 , M
            D(i,i) = ONE
         ENDDO
!
         DO i = 1 , M
            IF ( i<=4 ) THEN
               A(i,i) = ONE
               IF ( i>2 ) A(i,i) = ONE + reeps
               IF ( MOD(i,2)/=0 .AND. i<M ) THEN
                  A(i,i+1) = imeps
               ELSEIF ( i>1 ) THEN
                  A(i,i-1) = -imeps
               ENDIF
            ELSEIF ( i<=8 ) THEN
               IF ( i<=6 ) THEN
                  A(i,i) = reeps
               ELSE
                  A(i,i) = -reeps
               ENDIF
               IF ( MOD(i,2)/=0 .AND. i<M ) THEN
                  A(i,i+1) = ONE
               ELSEIF ( i>1 ) THEN
                  A(i,i-1) = -ONE
               ENDIF
            ELSE
               A(i,i) = ONE
               IF ( MOD(i,2)/=0 .AND. i<M ) THEN
                  A(i,i+1) = imeps*2
               ELSEIF ( i>1 ) THEN
                  A(i,i-1) = -imeps*2
               ENDIF
            ENDIF
         ENDDO
!
         DO i = 1 , N
            E(i,i) = ONE
            IF ( i<=4 ) THEN
               B(i,i) = -ONE
               IF ( i>2 ) B(i,i) = ONE - reeps
               IF ( MOD(i,2)/=0 .AND. i<N ) THEN
                  B(i,i+1) = imeps
               ELSEIF ( i>1 ) THEN
                  B(i,i-1) = -imeps
               ENDIF
            ELSEIF ( i<=8 ) THEN
               IF ( i<=6 ) THEN
                  B(i,i) = reeps
               ELSE
                  B(i,i) = -reeps
               ENDIF
               IF ( MOD(i,2)/=0 .AND. i<N ) THEN
                  B(i,i+1) = ONE + imeps
               ELSEIF ( i>1 ) THEN
                  B(i,i-1) = -ONE - imeps
               ENDIF
            ELSE
               B(i,i) = ONE - reeps
               IF ( MOD(i,2)/=0 .AND. i<N ) THEN
                  B(i,i+1) = imeps*2
               ELSEIF ( i>1 ) THEN
                  B(i,i-1) = -imeps*2
               ENDIF
            ENDIF
         ENDDO
      ENDIF
!
!     Compute rhs (C, F)
!
      CALL DGEMM('N','N',M,N,M,ONE,A,Lda,R,Ldr,ZERO,C,Ldc)
      CALL DGEMM('N','N',M,N,N,-ONE,L,Ldl,B,Ldb,ONE,C,Ldc)
      CALL DGEMM('N','N',M,N,M,ONE,D,Ldd,R,Ldr,ZERO,F,Ldf)
      CALL DGEMM('N','N',M,N,N,-ONE,L,Ldl,E,Lde,ONE,F,Ldf)
!
!     End of DLATM5
!
      END SUBROUTINE DLATM5
