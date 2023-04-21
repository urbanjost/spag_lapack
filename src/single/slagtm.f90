!*==slagtm.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b SLAGTM performs a matrix-matrix product of the form C = αAB+βC, where A is a tridiagonal matrix, B and C are rectangular matrices, and α and β are scalars, which may be 0, 1, or -1.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SLAGTM + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slagtm.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slagtm.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slagtm.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SLAGTM( TRANS, N, NRHS, ALPHA, DL, D, DU, X, LDX, BETA,
!                          B, LDB )
!
!       .. Scalar Arguments ..
!       CHARACTER          TRANS
!       INTEGER            LDB, LDX, N, NRHS
!       REAL               ALPHA, BETA
!       ..
!       .. Array Arguments ..
!       REAL               B( LDB, * ), D( * ), DL( * ), DU( * ),
!      $                   X( LDX, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SLAGTM performs a matrix-vector product of the form
!>
!>    B := alpha * A * X + beta * B
!>
!> where A is a tridiagonal matrix of order N, B and X are N by NRHS
!> matrices, and alpha and beta are real scalars, each of which may be
!> 0., 1., or -1.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] TRANS
!> \verbatim
!>          TRANS is CHARACTER*1
!>          Specifies the operation applied to A.
!>          = 'N':  No transpose, B := alpha * A * X + beta * B
!>          = 'T':  Transpose,    B := alpha * A'* X + beta * B
!>          = 'C':  Conjugate transpose = Transpose
!> \endverbatim
!>
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
!>          The scalar alpha.  ALPHA must be 0., 1., or -1.; otherwise,
!>          it is assumed to be 0.
!> \endverbatim
!>
!> \param[in] DL
!> \verbatim
!>          DL is REAL array, dimension (N-1)
!>          The (n-1) sub-diagonal elements of T.
!> \endverbatim
!>
!> \param[in] D
!> \verbatim
!>          D is REAL array, dimension (N)
!>          The diagonal elements of T.
!> \endverbatim
!>
!> \param[in] DU
!> \verbatim
!>          DU is REAL array, dimension (N-1)
!>          The (n-1) super-diagonal elements of T.
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
!> \ingroup realOTHERauxiliary
!
!  =====================================================================
      SUBROUTINE SLAGTM(Trans,N,Nrhs,Alpha,Dl,D,Du,X,Ldx,Beta,B,Ldb)
      IMPLICIT NONE
!*--SLAGTM148
!
!  -- LAPACK auxiliary routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Trans
      INTEGER Ldb , Ldx , N , Nrhs
      REAL Alpha , Beta
!     ..
!     .. Array Arguments ..
      REAL B(Ldb,*) , D(*) , Dl(*) , Du(*) , X(Ldx,*)
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
!     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
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
         IF ( LSAME(Trans,'N') ) THEN
!
!           Compute B := B + A*X
!
            DO j = 1 , Nrhs
               IF ( N==1 ) THEN
                  B(1,j) = B(1,j) + D(1)*X(1,j)
               ELSE
                  B(1,j) = B(1,j) + D(1)*X(1,j) + Du(1)*X(2,j)
                  B(N,j) = B(N,j) + Dl(N-1)*X(N-1,j) + D(N)*X(N,j)
                  DO i = 2 , N - 1
                     B(i,j) = B(i,j) + Dl(i-1)*X(i-1,j) + D(i)*X(i,j)   &
     &                        + Du(i)*X(i+1,j)
                  ENDDO
               ENDIF
            ENDDO
         ELSE
!
!           Compute B := B + A**T*X
!
            DO j = 1 , Nrhs
               IF ( N==1 ) THEN
                  B(1,j) = B(1,j) + D(1)*X(1,j)
               ELSE
                  B(1,j) = B(1,j) + D(1)*X(1,j) + Dl(1)*X(2,j)
                  B(N,j) = B(N,j) + Du(N-1)*X(N-1,j) + D(N)*X(N,j)
                  DO i = 2 , N - 1
                     B(i,j) = B(i,j) + Du(i-1)*X(i-1,j) + D(i)*X(i,j)   &
     &                        + Dl(i)*X(i+1,j)
                  ENDDO
               ENDIF
            ENDDO
         ENDIF
      ELSEIF ( Alpha==-ONE ) THEN
         IF ( LSAME(Trans,'N') ) THEN
!
!           Compute B := B - A*X
!
            DO j = 1 , Nrhs
               IF ( N==1 ) THEN
                  B(1,j) = B(1,j) - D(1)*X(1,j)
               ELSE
                  B(1,j) = B(1,j) - D(1)*X(1,j) - Du(1)*X(2,j)
                  B(N,j) = B(N,j) - Dl(N-1)*X(N-1,j) - D(N)*X(N,j)
                  DO i = 2 , N - 1
                     B(i,j) = B(i,j) - Dl(i-1)*X(i-1,j) - D(i)*X(i,j)   &
     &                        - Du(i)*X(i+1,j)
                  ENDDO
               ENDIF
            ENDDO
         ELSE
!
!           Compute B := B - A**T*X
!
            DO j = 1 , Nrhs
               IF ( N==1 ) THEN
                  B(1,j) = B(1,j) - D(1)*X(1,j)
               ELSE
                  B(1,j) = B(1,j) - D(1)*X(1,j) - Dl(1)*X(2,j)
                  B(N,j) = B(N,j) - Du(N-1)*X(N-1,j) - D(N)*X(N,j)
                  DO i = 2 , N - 1
                     B(i,j) = B(i,j) - Du(i-1)*X(i-1,j) - D(i)*X(i,j)   &
     &                        - Dl(i)*X(i+1,j)
                  ENDDO
               ENDIF
            ENDDO
         ENDIF
      ENDIF
!
!     End of SLAGTM
!
      END SUBROUTINE SLAGTM
