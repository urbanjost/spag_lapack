!*==dgesc2.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b DGESC2 solves a system of linear equations using the LU factorization with complete pivoting computed by sgetc2.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DGESC2 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dgesc2.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dgesc2.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dgesc2.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DGESC2( N, A, LDA, RHS, IPIV, JPIV, SCALE )
!
!       .. Scalar Arguments ..
!       INTEGER            LDA, N
!       DOUBLE PRECISION   SCALE
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * ), JPIV( * )
!       DOUBLE PRECISION   A( LDA, * ), RHS( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DGESC2 solves a system of linear equations
!>
!>           A * X = scale* RHS
!>
!> with a general N-by-N matrix A using the LU factorization with
!> complete pivoting computed by DGETC2.
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
!>          A is DOUBLE PRECISION array, dimension (LDA,N)
!>          On entry, the  LU part of the factorization of the n-by-n
!>          matrix A computed by DGETC2:  A = P * L * U * Q
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
!>          RHS is DOUBLE PRECISION array, dimension (N).
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
!>          SCALE is DOUBLE PRECISION
!>          On exit, SCALE contains the scale factor. SCALE is chosen
!>          0 <= SCALE <= 1 to prevent overflow in the solution.
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
!> \date November 2017
!
!> \ingroup doubleGEauxiliary
!
!> \par Contributors:
!  ==================
!>
!>     Bo Kagstrom and Peter Poromaa, Department of Computing Science,
!>     Umea University, S-901 87 Umea, Sweden.
!
!  =====================================================================
      SUBROUTINE DGESC2(N,A,Lda,Rhs,Ipiv,Jpiv,Scale)
      USE F77KINDS                        
      USE S_DLABAD
      USE S_DLAMCH
      USE S_DLASWP
      USE S_DSCAL
      USE S_IDAMAX
      IMPLICIT NONE
!*--DGESC2124
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , TWO = 2.0D+0
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER :: N
      REAL(R8KIND) , INTENT(IN) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Rhs
      INTEGER , DIMENSION(*) :: Ipiv
      INTEGER , DIMENSION(*) :: Jpiv
      REAL(R8KIND) , INTENT(INOUT) :: Scale
!
! Local variable declarations rewritten by SPAG
!
      REAL(R8KIND) :: bignum , eps , smlnum , temp
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
      eps = DLAMCH('P')
      smlnum = DLAMCH('S')/eps
      bignum = ONE/smlnum
      CALL DLABAD(smlnum,bignum)
!
!     Apply permutations IPIV to RHS
!
      CALL DLASWP(1,Rhs,Lda,1,N-1,Ipiv,1)
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
      i = IDAMAX(N,Rhs,1)
      IF ( TWO*smlnum*ABS(Rhs(i))>ABS(A(N,N)) ) THEN
         temp = (ONE/TWO)/ABS(Rhs(i))
         CALL DSCAL(N,temp,Rhs(1),1)
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
      CALL DLASWP(1,Rhs,Lda,1,N-1,Jpiv,-1)
!
!     End of DGESC2
!
      END SUBROUTINE DGESC2
