!*==sgesc2.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b SGESC2 solves a system of linear equations using the LU factorization with complete pivoting computed by sgetc2.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SGESC2 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgesc2.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgesc2.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgesc2.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SGESC2( N, A, LDA, RHS, IPIV, JPIV, SCALE )
!
!       .. Scalar Arguments ..
!       INTEGER            LDA, N
!       REAL               SCALE
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * ), JPIV( * )
!       REAL               A( LDA, * ), RHS( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SGESC2 solves a system of linear equations
!>
!>           A * X = scale* RHS
!>
!> with a general N-by-N matrix A using the LU factorization with
!> complete pivoting computed by SGETC2.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix A.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is REAL array, dimension (LDA,N)
!>          On entry, the  LU part of the factorization of the n-by-n
!>          matrix A computed by SGETC2:  A = P * L * U * Q
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1, N).
!> \endverbatim
!>
!> \param[in,out] RHS
!> \verbatim
!>          RHS is REAL array, dimension (N).
!>          On entry, the right hand side vector b.
!>          On exit, the solution vector X.
!> \endverbatim
!>
!> \param[in] IPIV
!> \verbatim
!>          IPIV is INTEGER array, dimension (N).
!>          The pivot indices; for 1 <= i <= N, row i of the
!>          matrix has been interchanged with row IPIV(i).
!> \endverbatim
!>
!> \param[in] JPIV
!> \verbatim
!>          JPIV is INTEGER array, dimension (N).
!>          The pivot indices; for 1 <= j <= N, column j of the
!>          matrix has been interchanged with column JPIV(j).
!> \endverbatim
!>
!> \param[out] SCALE
!> \verbatim
!>          SCALE is REAL
!>           On exit, SCALE contains the scale factor. SCALE is chosen
!>           0 <= SCALE <= 1 to prevent overflow in the solution.
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
!> \ingroup realGEauxiliary
!
!> \par Contributors:
!  ==================
!>
!>     Bo Kagstrom and Peter Poromaa, Department of Computing Science,
!>     Umea University, S-901 87 Umea, Sweden.
!
!  =====================================================================
      SUBROUTINE SGESC2(N,A,Lda,Rhs,Ipiv,Jpiv,Scale)
      USE S_ISAMAX
      USE S_SLABAD
      USE S_SLAMCH
      USE S_SLASWP
      USE S_SSCAL
      IMPLICIT NONE
!*--SGESC2123
!
! PARAMETER definitions rewritten by SPAG
!
      REAL , PARAMETER  ::  ONE = 1.0E+0 , TWO = 2.0E+0
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER :: N
      REAL , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , INTENT(INOUT) , DIMENSION(*) :: Rhs
      INTEGER , DIMENSION(*) :: Ipiv
      INTEGER , DIMENSION(*) :: Jpiv
      REAL , INTENT(INOUT) :: Scale
!
! Local variable declarations rewritten by SPAG
!
      REAL :: bignum , eps , smlnum , temp
      INTEGER :: i , j
!
! End of declarations rewritten by SPAG
!
!     ..
!     .. Array Arguments ..
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
!     ..
!     .. Local Scalars ..
!     ..
!     .. External Subroutines ..
!     ..
!     .. External Functions ..
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. Executable Statements ..
!
!      Set constant to control overflow
!
      eps = SLAMCH('P')
      smlnum = SLAMCH('S')/eps
      bignum = ONE/smlnum
      CALL SLABAD(smlnum,bignum)
!
!     Apply permutations IPIV to RHS
!
      CALL SLASWP(1,Rhs,Lda,1,N-1,Ipiv,1)
!
!     Solve for L part
!
      DO i = 1 , N - 1
         DO j = i + 1 , N
            Rhs(j) = Rhs(j) - A(j,i)*Rhs(i)
         ENDDO
      ENDDO
!
!     Solve for U part
!
      Scale = ONE
!
!     Check for scaling
!
      i = ISAMAX(N,Rhs,1)
      IF ( TWO*smlnum*ABS(Rhs(i))>ABS(A(N,N)) ) THEN
         temp = (ONE/TWO)/ABS(Rhs(i))
         CALL SSCAL(N,temp,Rhs(1),1)
         Scale = Scale*temp
      ENDIF
!
      DO i = N , 1 , -1
         temp = ONE/A(i,i)
         Rhs(i) = Rhs(i)*temp
         DO j = i + 1 , N
            Rhs(i) = Rhs(i) - Rhs(j)*(A(i,j)*temp)
         ENDDO
      ENDDO
!
!     Apply permutations JPIV to the solution (RHS)
!
      CALL SLASWP(1,Rhs,Lda,1,N-1,Jpiv,-1)
!
!     End of SGESC2
!
      END SUBROUTINE SGESC2
