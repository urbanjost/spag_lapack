!*==clatm6.f90  processed by SPAG 7.51RB at 17:44 on  4 Mar 2022
!> \brief \b CLATM6
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE CLATM6( TYPE, N, A, LDA, B, X, LDX, Y, LDY, ALPHA,
!                          BETA, WX, WY, S, DIF )
!
!       .. Scalar Arguments ..
!       INTEGER            LDA, LDX, LDY, N, TYPE
!       COMPLEX            ALPHA, BETA, WX, WY
!       ..
!       .. Array Arguments ..
!       REAL               DIF( * ), S( * )
!       COMPLEX            A( LDA, * ), B( LDA, * ), X( LDX, * ),
!      $                   Y( LDY, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CLATM6 generates test matrices for the generalized eigenvalue
!> problem, their corresponding right and left eigenvector matrices,
!> and also reciprocal condition numbers for all eigenvalues and
!> the reciprocal condition numbers of eigenvectors corresponding to
!> the 1th and 5th eigenvalues.
!>
!> Test Matrices
!> =============
!>
!> Two kinds of test matrix pairs
!>          (A, B) = inverse(YH) * (Da, Db) * inverse(X)
!> are used in the tests:
!>
!> Type 1:
!>    Da = 1+a   0    0    0    0    Db = 1   0   0   0   0
!>          0   2+a   0    0    0         0   1   0   0   0
!>          0    0   3+a   0    0         0   0   1   0   0
!>          0    0    0   4+a   0         0   0   0   1   0
!>          0    0    0    0   5+a ,      0   0   0   0   1
!> and Type 2:
!>    Da = 1+i   0    0       0       0    Db = 1   0   0   0   0
!>          0   1-i   0       0       0         0   1   0   0   0
!>          0    0    1       0       0         0   0   1   0   0
!>          0    0    0 (1+a)+(1+b)i  0         0   0   0   1   0
!>          0    0    0       0 (1+a)-(1+b)i,   0   0   0   0   1 .
!>
!> In both cases the same inverse(YH) and inverse(X) are used to compute
!> (A, B), giving the exact eigenvectors to (A,B) as (YH, X):
!>
!> YH:  =  1    0   -y    y   -y    X =  1   0  -x  -x   x
!>         0    1   -y    y   -y         0   1   x  -x  -x
!>         0    0    1    0    0         0   0   1   0   0
!>         0    0    0    1    0         0   0   0   1   0
!>         0    0    0    0    1,        0   0   0   0   1 , where
!>
!> a, b, x and y will have all values independently of each other.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] TYPE
!> \verbatim
!>          TYPE is INTEGER
!>          Specifies the problem type (see further details).
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          Size of the matrices A and B.
!> \endverbatim
!>
!> \param[out] A
!> \verbatim
!>          A is COMPLEX array, dimension (LDA, N).
!>          On exit A N-by-N is initialized according to TYPE.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of A and of B.
!> \endverbatim
!>
!> \param[out] B
!> \verbatim
!>          B is COMPLEX array, dimension (LDA, N).
!>          On exit B N-by-N is initialized according to TYPE.
!> \endverbatim
!>
!> \param[out] X
!> \verbatim
!>          X is COMPLEX array, dimension (LDX, N).
!>          On exit X is the N-by-N matrix of right eigenvectors.
!> \endverbatim
!>
!> \param[in] LDX
!> \verbatim
!>          LDX is INTEGER
!>          The leading dimension of X.
!> \endverbatim
!>
!> \param[out] Y
!> \verbatim
!>          Y is COMPLEX array, dimension (LDY, N).
!>          On exit Y is the N-by-N matrix of left eigenvectors.
!> \endverbatim
!>
!> \param[in] LDY
!> \verbatim
!>          LDY is INTEGER
!>          The leading dimension of Y.
!> \endverbatim
!>
!> \param[in] ALPHA
!> \verbatim
!>          ALPHA is COMPLEX
!> \endverbatim
!>
!> \param[in] BETA
!> \verbatim
!>          BETA is COMPLEX
!>
!>          Weighting constants for matrix A.
!> \endverbatim
!>
!> \param[in] WX
!> \verbatim
!>          WX is COMPLEX
!>          Constant for right eigenvector matrix.
!> \endverbatim
!>
!> \param[in] WY
!> \verbatim
!>          WY is COMPLEX
!>          Constant for left eigenvector matrix.
!> \endverbatim
!>
!> \param[out] S
!> \verbatim
!>          S is REAL array, dimension (N)
!>          S(i) is the reciprocal condition number for eigenvalue i.
!> \endverbatim
!>
!> \param[out] DIF
!> \verbatim
!>          DIF is REAL array, dimension (N)
!>          DIF(i) is the reciprocal condition number for eigenvector i.
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
!> \ingroup complex_matgen
!
!  =====================================================================
      SUBROUTINE CLATM6(Type,N,A,Lda,B,X,Ldx,Y,Ldy,Alpha,Beta,Wx,Wy,S,  &
     &                  Dif)
      IMPLICIT NONE
!*--CLATM6178
!
!  -- LAPACK computational routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Lda , Ldx , Ldy , N , Type
      COMPLEX Alpha , Beta , Wx , Wy
!     ..
!     .. Array Arguments ..
      REAL Dif(*) , S(*)
      COMPLEX A(Lda,*) , B(Lda,*) , X(Ldx,*) , Y(Ldy,*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL RONE , TWO , THREE
      PARAMETER (RONE=1.0E+0,TWO=2.0E+0,THREE=3.0E+0)
      COMPLEX ZERO , ONE
      PARAMETER (ZERO=(0.0E+0,0.0E+0),ONE=(1.0E+0,0.0E+0))
!     ..
!     .. Local Scalars ..
      INTEGER i , info , j
!     ..
!     .. Local Arrays ..
      REAL rwork(50)
      COMPLEX work(26) , z(8,8)
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC CABS , CMPLX , CONJG , REAL , SQRT
!     ..
!     .. External Subroutines ..
      EXTERNAL CGESVD , CLACPY 
!     ..
!     .. Executable Statements ..
!
!     Generate test problem ...
!     (Da, Db) ...
!
      DO i = 1 , N
         DO j = 1 , N
!
            IF ( i==j ) THEN
               A(i,i) = CMPLX(i) + Alpha
               B(i,i) = ONE
            ELSE
               A(i,j) = ZERO
               B(i,j) = ZERO
            ENDIF
!
         ENDDO
      ENDDO
      IF ( Type==2 ) THEN
         A(1,1) = CMPLX(RONE,RONE)
         A(2,2) = CONJG(A(1,1))
         A(3,3) = ONE
         A(4,4) = CMPLX(REAL(ONE+Alpha),REAL(ONE+Beta))
         A(5,5) = CONJG(A(4,4))
      ENDIF
!
!     Form X and Y
!
      CALL CLACPY('F',N,N,B,Lda,Y,Ldy)
      Y(3,1) = -CONJG(Wy)
      Y(4,1) = CONJG(Wy)
      Y(5,1) = -CONJG(Wy)
      Y(3,2) = -CONJG(Wy)
      Y(4,2) = CONJG(Wy)
      Y(5,2) = -CONJG(Wy)
!
      CALL CLACPY('F',N,N,B,Lda,X,Ldx)
      X(1,3) = -Wx
      X(1,4) = -Wx
      X(1,5) = Wx
      X(2,3) = Wx
      X(2,4) = -Wx
      X(2,5) = -Wx
!
!     Form (A, B)
!
      B(1,3) = Wx + Wy
      B(2,3) = -Wx + Wy
      B(1,4) = Wx - Wy
      B(2,4) = Wx - Wy
      B(1,5) = -Wx + Wy
      B(2,5) = Wx + Wy
      A(1,3) = Wx*A(1,1) + Wy*A(3,3)
      A(2,3) = -Wx*A(2,2) + Wy*A(3,3)
      A(1,4) = Wx*A(1,1) - Wy*A(4,4)
      A(2,4) = Wx*A(2,2) - Wy*A(4,4)
      A(1,5) = -Wx*A(1,1) + Wy*A(5,5)
      A(2,5) = Wx*A(2,2) + Wy*A(5,5)
!
!     Compute condition numbers
!
      S(1) = RONE/SQRT((RONE+THREE*CABS(Wy)*CABS(Wy))                   &
     &       /(RONE+CABS(A(1,1))*CABS(A(1,1))))
      S(2) = RONE/SQRT((RONE+THREE*CABS(Wy)*CABS(Wy))                   &
     &       /(RONE+CABS(A(2,2))*CABS(A(2,2))))
      S(3) = RONE/SQRT((RONE+TWO*CABS(Wx)*CABS(Wx))                     &
     &       /(RONE+CABS(A(3,3))*CABS(A(3,3))))
      S(4) = RONE/SQRT((RONE+TWO*CABS(Wx)*CABS(Wx))                     &
     &       /(RONE+CABS(A(4,4))*CABS(A(4,4))))
      S(5) = RONE/SQRT((RONE+TWO*CABS(Wx)*CABS(Wx))                     &
     &       /(RONE+CABS(A(5,5))*CABS(A(5,5))))
!
      CALL clakf2(1,4,A,Lda,A(2,2),B,B(2,2),z,8)
      CALL CGESVD('N','N',8,8,z,8,rwork,work,1,work(2),1,work(3),24,    &
     &            rwork(9),info)
      Dif(1) = rwork(8)
!
      CALL clakf2(4,1,A,Lda,A(5,5),B,B(5,5),z,8)
      CALL CGESVD('N','N',8,8,z,8,rwork,work,1,work(2),1,work(3),24,    &
     &            rwork(9),info)
      Dif(5) = rwork(8)
!
!
!     End of CLATM6
!
      END SUBROUTINE CLATM6
