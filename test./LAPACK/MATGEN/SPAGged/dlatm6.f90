!*==dlatm6.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b DLATM6
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLATM6( TYPE, N, A, LDA, B, X, LDX, Y, LDY, ALPHA,
!                          BETA, WX, WY, S, DIF )
!
!       .. Scalar Arguments ..
!       INTEGER            LDA, LDX, LDY, N, TYPE
!       DOUBLE PRECISION   ALPHA, BETA, WX, WY
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   A( LDA, * ), B( LDA, * ), DIF( * ), S( * ),
!      $                   X( LDX, * ), Y( LDY, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLATM6 generates test matrices for the generalized eigenvalue
!> problem, their corresponding right and left eigenvector matrices,
!> and also reciprocal condition numbers for all eigenvalues and
!> the reciprocal condition numbers of eigenvectors corresponding to
!> the 1th and 5th eigenvalues.
!>
!> Test Matrices
!> =============
!>
!> Two kinds of test matrix pairs
!>
!>       (A, B) = inverse(YH) * (Da, Db) * inverse(X)
!>
!> are used in the tests:
!>
!> Type 1:
!>    Da = 1+a   0    0    0    0    Db = 1   0   0   0   0
!>          0   2+a   0    0    0         0   1   0   0   0
!>          0    0   3+a   0    0         0   0   1   0   0
!>          0    0    0   4+a   0         0   0   0   1   0
!>          0    0    0    0   5+a ,      0   0   0   0   1 , and
!>
!> Type 2:
!>    Da =  1   -1    0    0    0    Db = 1   0   0   0   0
!>          1    1    0    0    0         0   1   0   0   0
!>          0    0    1    0    0         0   0   1   0   0
!>          0    0    0   1+a  1+b        0   0   0   1   0
!>          0    0    0  -1-b  1+a ,      0   0   0   0   1 .
!>
!> In both cases the same inverse(YH) and inverse(X) are used to compute
!> (A, B), giving the exact eigenvectors to (A,B) as (YH, X):
!>
!> YH:  =  1    0   -y    y   -y    X =  1   0  -x  -x   x
!>         0    1   -y    y   -y         0   1   x  -x  -x
!>         0    0    1    0    0         0   0   1   0   0
!>         0    0    0    1    0         0   0   0   1   0
!>         0    0    0    0    1,        0   0   0   0   1 ,
!>
!> where a, b, x and y will have all values independently of each other.
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
!>          A is DOUBLE PRECISION array, dimension (LDA, N).
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
!>          B is DOUBLE PRECISION array, dimension (LDA, N).
!>          On exit B N-by-N is initialized according to TYPE.
!> \endverbatim
!>
!> \param[out] X
!> \verbatim
!>          X is DOUBLE PRECISION array, dimension (LDX, N).
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
!>          Y is DOUBLE PRECISION array, dimension (LDY, N).
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
!>          ALPHA is DOUBLE PRECISION
!> \endverbatim
!>
!> \param[in] BETA
!> \verbatim
!>          BETA is DOUBLE PRECISION
!>
!>          Weighting constants for matrix A.
!> \endverbatim
!>
!> \param[in] WX
!> \verbatim
!>          WX is DOUBLE PRECISION
!>          Constant for right eigenvector matrix.
!> \endverbatim
!>
!> \param[in] WY
!> \verbatim
!>          WY is DOUBLE PRECISION
!>          Constant for left eigenvector matrix.
!> \endverbatim
!>
!> \param[out] S
!> \verbatim
!>          S is DOUBLE PRECISION array, dimension (N)
!>          S(i) is the reciprocal condition number for eigenvalue i.
!> \endverbatim
!>
!> \param[out] DIF
!> \verbatim
!>          DIF is DOUBLE PRECISION array, dimension (N)
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
!> \ingroup double_matgen
!
!  =====================================================================
      SUBROUTINE DLATM6(Type,N,A,Lda,B,X,Ldx,Y,Ldy,Alpha,Beta,Wx,Wy,S,  &
     &                  Dif)
      IMPLICIT NONE
!*--DLATM6180
!
!  -- LAPACK computational routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Lda , Ldx , Ldy , N , Type
      DOUBLE PRECISION Alpha , Beta , Wx , Wy
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION A(Lda,*) , B(Lda,*) , Dif(*) , S(*) , X(Ldx,*) , &
     &                 Y(Ldy,*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO , ONE , TWO , THREE
      PARAMETER (ZERO=0.0D+0,ONE=1.0D+0,TWO=2.0D+0,THREE=3.0D+0)
!     ..
!     .. Local Scalars ..
      INTEGER i , info , j
!     ..
!     .. Local Arrays ..
      DOUBLE PRECISION work(100) , z(12,12)
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC DBLE , SQRT
!     ..
!     .. External Subroutines ..
      EXTERNAL DGESVD , DLACPY , DLAKF2
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
               A(i,i) = DBLE(i) + Alpha
               B(i,i) = ONE
            ELSE
               A(i,j) = ZERO
               B(i,j) = ZERO
            ENDIF
!
         ENDDO
      ENDDO
!
!     Form X and Y
!
      CALL DLACPY('F',N,N,B,Lda,Y,Ldy)
      Y(3,1) = -Wy
      Y(4,1) = Wy
      Y(5,1) = -Wy
      Y(3,2) = -Wy
      Y(4,2) = Wy
      Y(5,2) = -Wy
!
      CALL DLACPY('F',N,N,B,Lda,X,Ldx)
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
      IF ( Type==1 ) THEN
         A(1,3) = Wx*A(1,1) + Wy*A(3,3)
         A(2,3) = -Wx*A(2,2) + Wy*A(3,3)
         A(1,4) = Wx*A(1,1) - Wy*A(4,4)
         A(2,4) = Wx*A(2,2) - Wy*A(4,4)
         A(1,5) = -Wx*A(1,1) + Wy*A(5,5)
         A(2,5) = Wx*A(2,2) + Wy*A(5,5)
      ELSEIF ( Type==2 ) THEN
         A(1,3) = TWO*Wx + Wy
         A(2,3) = Wy
         A(1,4) = -Wy*(TWO+Alpha+Beta)
         A(2,4) = TWO*Wx - Wy*(TWO+Alpha+Beta)
         A(1,5) = -TWO*Wx + Wy*(Alpha-Beta)
         A(2,5) = Wy*(Alpha-Beta)
         A(1,1) = ONE
         A(1,2) = -ONE
         A(2,1) = ONE
         A(2,2) = A(1,1)
         A(3,3) = ONE
         A(4,4) = ONE + Alpha
         A(4,5) = ONE + Beta
         A(5,4) = -A(4,5)
         A(5,5) = A(4,4)
      ENDIF
!
!     Compute condition numbers
!
      IF ( Type==1 ) THEN
!
         S(1) = ONE/SQRT((ONE+THREE*Wy*Wy)/(ONE+A(1,1)*A(1,1)))
         S(2) = ONE/SQRT((ONE+THREE*Wy*Wy)/(ONE+A(2,2)*A(2,2)))
         S(3) = ONE/SQRT((ONE+TWO*Wx*Wx)/(ONE+A(3,3)*A(3,3)))
         S(4) = ONE/SQRT((ONE+TWO*Wx*Wx)/(ONE+A(4,4)*A(4,4)))
         S(5) = ONE/SQRT((ONE+TWO*Wx*Wx)/(ONE+A(5,5)*A(5,5)))
!
         CALL DLAKF2(1,4,A,Lda,A(2,2),B,B(2,2),z,12)
         CALL DGESVD('N','N',8,8,z,12,work,work(9),1,work(10),1,work(11)&
     &               ,40,info)
         Dif(1) = work(8)
!
         CALL DLAKF2(4,1,A,Lda,A(5,5),B,B(5,5),z,12)
         CALL DGESVD('N','N',8,8,z,12,work,work(9),1,work(10),1,work(11)&
     &               ,40,info)
         Dif(5) = work(8)
!
      ELSEIF ( Type==2 ) THEN
!
         S(1) = ONE/SQRT(ONE/THREE+Wy*Wy)
         S(2) = S(1)
         S(3) = ONE/SQRT(ONE/TWO+Wx*Wx)
         S(4) = ONE/SQRT((ONE+TWO*Wx*Wx)/(ONE+(ONE+Alpha)*(ONE+Alpha)+( &
     &          ONE+Beta)*(ONE+Beta)))
         S(5) = S(4)
!
         CALL DLAKF2(2,3,A,Lda,A(3,3),B,B(3,3),z,12)
         CALL DGESVD('N','N',12,12,z,12,work,work(13),1,work(14),1,     &
     &               work(15),60,info)
         Dif(1) = work(12)
!
         CALL DLAKF2(3,2,A,Lda,A(4,4),B,B(4,4),z,12)
         CALL DGESVD('N','N',12,12,z,12,work,work(13),1,work(14),1,     &
     &               work(15),60,info)
         Dif(5) = work(12)
!
      ENDIF
!
!
!     End of DLATM6
!
      END SUBROUTINE DLATM6
