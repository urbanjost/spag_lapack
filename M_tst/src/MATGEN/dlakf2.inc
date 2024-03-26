!*==dlakf2.f90  processed by SPAG 7.51RB at 17:44 on  4 Mar 2022
 
!> \brief \b dlakf2
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE dlakf2( M, N, A, LDA, B, D, E, Z, LDZ )
!
!       .. Scalar Arguments ..
!       INTEGER            LDA, LDZ, M, N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   A( LDA, * ), B( LDA, * ), D( LDA, * ),
!      $                   E( LDA, * ), Z( LDZ, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> Form the 2*M*N by 2*M*N matrix
!>
!>        Z = [ kron(In, A)  -kron(B', Im) ]
!>            [ kron(In, D)  -kron(E', Im) ],
!>
!> where In is the identity matrix of size n and X' is the transpose
!> of X. kron(X, Y) is the Kronecker product between the matrices X
!> and Y.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          Size of matrix, must be >= 1.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          Size of matrix, must be >= 1.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is DOUBLE PRECISION, dimension ( LDA, M )
!>          The matrix A in the output matrix Z.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of A, B, D, and E. ( LDA >= M+N )
!> \endverbatim
!>
!> \param[in] B
!> \verbatim
!>          B is DOUBLE PRECISION, dimension ( LDA, N )
!> \endverbatim
!>
!> \param[in] D
!> \verbatim
!>          D is DOUBLE PRECISION, dimension ( LDA, M )
!> \endverbatim
!>
!> \param[in] E
!> \verbatim
!>          E is DOUBLE PRECISION, dimension ( LDA, N )
!>
!>          The matrices used in forming the output matrix Z.
!> \endverbatim
!>
!> \param[out] Z
!> \verbatim
!>          Z is DOUBLE PRECISION, dimension ( LDZ, 2*M*N )
!>          The resultant Kronecker M*N*2 by M*N*2 matrix (see above.)
!> \endverbatim
!>
!> \param[in] LDZ
!> \verbatim
!>          LDZ is INTEGER
!>          The leading dimension of Z. ( LDZ >= 2*M*N )
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
      SUBROUTINE dlakf2(M,N,A,Lda,B,D,E,Z,Ldz)
      IMPLICIT NONE
!*--dlakf2110
!
!  -- LAPACK computational routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Lda , Ldz , M , N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION A(Lda,*) , B(Lda,*) , D(Lda,*) , E(Lda,*) ,      &
     &                 Z(Ldz,*)
!     ..
!
!  ====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D+0)
!     ..
!     .. Local Scalars ..
      INTEGER i , ik , j , jk , l , mn , mn2
!     ..
!     .. External Subroutines ..
      EXTERNAL DLASET
!     ..
!     .. Executable Statements ..
!
!     Initialize Z
!
      mn = M*N
      mn2 = 2*mn
      CALL DLASET('Full',mn2,mn2,ZERO,ZERO,Z,Ldz)
!
      ik = 1
      DO l = 1 , N
!
!        form kron(In, A)
!
         DO i = 1 , M
            DO j = 1 , M
               Z(ik+i-1,ik+j-1) = A(i,j)
            ENDDO
         ENDDO
!
!        form kron(In, D)
!
         DO i = 1 , M
            DO j = 1 , M
               Z(ik+mn+i-1,ik+j-1) = D(i,j)
            ENDDO
         ENDDO
!
         ik = ik + M
      ENDDO
!
      ik = 1
      DO l = 1 , N
         jk = mn + 1
!
         DO j = 1 , N
!
!           form -kron(B', Im)
!
            DO i = 1 , M
               Z(ik+i-1,jk+i-1) = -B(j,l)
            ENDDO
!
!           form -kron(E', Im)
!
            DO i = 1 , M
               Z(ik+mn+i-1,jk+i-1) = -E(j,l)
            ENDDO
!
            jk = jk + M
         ENDDO
!
         ik = ik + M
      ENDDO
!
!
!     End of dlakf2
!
      END SUBROUTINE dlakf2
