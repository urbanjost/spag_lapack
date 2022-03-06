!*==slaptm.f90  processed by SPAG 7.51RB at 17:34 on  4 Mar 2022
!> \brief \b slaptm
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE slaptm( N, NRHS, ALPHA, D, E, X, LDX, BETA, B, LDB )
!
!       .. Scalar Arguments ..
!       INTEGER            LDB, LDX, N, NRHS
!       REAL               ALPHA, BETA
!       ..
!       .. Array Arguments ..
!       REAL               B( LDB, * ), D( * ), E( * ), X( LDX, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> slaptm multiplies an N by NRHS matrix X by a symmetric tridiagonal
!> matrix A and stores the result in a matrix B.  The operation has the
!> form
!>
!>    B := alpha * A * X + beta * B
!>
!> where alpha may be either 1. or -1. and beta may be 0., 1., or -1.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in] NRHS
!> \verbatim
!>          NRHS is INTEGER
!>          The number of right hand sides, i.e., the number of columns
!>          of the matrices X and B.
!> \endverbatim
!>
!> \param[in] ALPHA
!> \verbatim
!>          ALPHA is REAL
!>          The scalar alpha.  ALPHA must be 1. or -1.; otherwise,
!>          it is assumed to be 0.
!> \endverbatim
!>
!> \param[in] D
!> \verbatim
!>          D is REAL array, dimension (N)
!>          The n diagonal elements of the tridiagonal matrix A.
!> \endverbatim
!>
!> \param[in] E
!> \verbatim
!>          E is REAL array, dimension (N-1)
!>          The (n-1) subdiagonal or superdiagonal elements of A.
!> \endverbatim
!>
!> \param[in] X
!> \verbatim
!>          X is REAL array, dimension (LDX,NRHS)
!>          The N by NRHS matrix X.
!> \endverbatim
!>
!> \param[in] LDX
!> \verbatim
!>          LDX is INTEGER
!>          The leading dimension of the array X.  LDX >= max(N,1).
!> \endverbatim
!>
!> \param[in] BETA
!> \verbatim
!>          BETA is REAL
!>          The scalar beta.  BETA must be 0., 1., or -1.; otherwise,
!>          it is assumed to be 1.
!> \endverbatim
!>
!> \param[in,out] B
!> \verbatim
!>          B is REAL array, dimension (LDB,NRHS)
!>          On entry, the N by NRHS matrix B.
!>          On exit, B is overwritten by the matrix expression
!>          B := alpha * A * X + beta * B.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of the array B.  LDB >= max(N,1).
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
!> \ingroup single_lin
!
!  =====================================================================
      SUBROUTINE slaptm(N,Nrhs,Alpha,D,E,X,Ldx,Beta,B,Ldb)
      IMPLICIT NONE
!*--SLAPTM120
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Ldb , Ldx , N , Nrhs
      REAL Alpha , Beta
!     ..
!     .. Array Arguments ..
      REAL B(Ldb,*) , D(*) , E(*) , X(Ldx,*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ONE , ZERO
      PARAMETER (ONE=1.0E+0,ZERO=0.0E+0)
!     ..
!     .. Local Scalars ..
      INTEGER i , j
!     ..
!     .. Executable Statements ..
!
      IF ( N==0 ) RETURN
!
!     Multiply B by BETA if BETA.NE.1.
!
      IF ( Beta==ZERO ) THEN
         DO j = 1 , Nrhs
            DO i = 1 , N
               B(i,j) = ZERO
            ENDDO
         ENDDO
      ELSEIF ( Beta==-ONE ) THEN
         DO j = 1 , Nrhs
            DO i = 1 , N
               B(i,j) = -B(i,j)
            ENDDO
         ENDDO
      ENDIF
!
      IF ( Alpha==ONE ) THEN
!
!        Compute B := B + A*X
!
         DO j = 1 , Nrhs
            IF ( N==1 ) THEN
               B(1,j) = B(1,j) + D(1)*X(1,j)
            ELSE
               B(1,j) = B(1,j) + D(1)*X(1,j) + E(1)*X(2,j)
               B(N,j) = B(N,j) + E(N-1)*X(N-1,j) + D(N)*X(N,j)
               DO i = 2 , N - 1
                  B(i,j) = B(i,j) + E(i-1)*X(i-1,j) + D(i)*X(i,j) + E(i)&
     &                     *X(i+1,j)
               ENDDO
            ENDIF
         ENDDO
      ELSEIF ( Alpha==-ONE ) THEN
!
!        Compute B := B - A*X
!
         DO j = 1 , Nrhs
            IF ( N==1 ) THEN
               B(1,j) = B(1,j) - D(1)*X(1,j)
            ELSE
               B(1,j) = B(1,j) - D(1)*X(1,j) - E(1)*X(2,j)
               B(N,j) = B(N,j) - E(N-1)*X(N-1,j) - D(N)*X(N,j)
               DO i = 2 , N - 1
                  B(i,j) = B(i,j) - E(i-1)*X(i-1,j) - D(i)*X(i,j) - E(i)&
     &                     *X(i+1,j)
               ENDDO
            ENDIF
         ENDDO
      ENDIF
!
!     End of slaptm
!
      END SUBROUTINE slaptm
