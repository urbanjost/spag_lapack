!*==cgesc2.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b CGESC2 solves a system of linear equations using the LU factorization with complete pivoting computed by sgetc2.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CGESC2 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgesc2.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgesc2.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgesc2.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CGESC2( N, A, LDA, RHS, IPIV, JPIV, SCALE )
!
!       .. Scalar Arguments ..
!       INTEGER            LDA, N
!       REAL               SCALE
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * ), JPIV( * )
!       COMPLEX            A( LDA, * ), RHS( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CGESC2 solves a system of linear equations
!>
!>           A * X = scale* RHS
!>
!> with a general N-by-N matrix A using the LU factorization with
!> complete pivoting computed by CGETC2.
!>
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix A.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is COMPLEX array, dimension (LDA, N)
!>          On entry, the  LU part of the factorization of the n-by-n
!>          matrix A computed by CGETC2:  A = P * L * U * Q
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
!>          RHS is COMPLEX array, dimension N.
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
!> \ingroup complexGEauxiliary
!
!> \par Contributors:
!  ==================
!>
!>     Bo Kagstrom and Peter Poromaa, Department of Computing Science,
!>     Umea University, S-901 87 Umea, Sweden.
!
!  =====================================================================
      SUBROUTINE CGESC2(N,A,Lda,Rhs,Ipiv,Jpiv,Scale)
      IMPLICIT NONE
!*--CGESC2119
!
!  -- LAPACK auxiliary routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Lda , N
      REAL Scale
!     ..
!     .. Array Arguments ..
      INTEGER Ipiv(*) , Jpiv(*)
      COMPLEX A(Lda,*) , Rhs(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ZERO , ONE , TWO
      PARAMETER (ZERO=0.0E+0,ONE=1.0E+0,TWO=2.0E+0)
!     ..
!     .. Local Scalars ..
      INTEGER i , j
      REAL bignum , eps , smlnum
      COMPLEX temp
!     ..
!     .. External Subroutines ..
      EXTERNAL CLASWP , CSCAL , SLABAD
!     ..
!     .. External Functions ..
      INTEGER ICAMAX
      REAL SLAMCH
      EXTERNAL ICAMAX , SLAMCH
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , CMPLX , REAL
!     ..
!     .. Executable Statements ..
!
!     Set constant to control overflow
!
      eps = SLAMCH('P')
      smlnum = SLAMCH('S')/eps
      bignum = ONE/smlnum
      CALL SLABAD(smlnum,bignum)
!
!     Apply permutations IPIV to RHS
!
      CALL CLASWP(1,Rhs,Lda,1,N-1,Ipiv,1)
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
      i = ICAMAX(N,Rhs,1)
      IF ( TWO*smlnum*ABS(Rhs(i))>ABS(A(N,N)) ) THEN
         temp = CMPLX(ONE/TWO,ZERO)/ABS(Rhs(i))
         CALL CSCAL(N,temp,Rhs(1),1)
         Scale = Scale*REAL(temp)
      ENDIF
      DO i = N , 1 , -1
         temp = CMPLX(ONE,ZERO)/A(i,i)
         Rhs(i) = Rhs(i)*temp
         DO j = i + 1 , N
            Rhs(i) = Rhs(i) - Rhs(j)*(A(i,j)*temp)
         ENDDO
      ENDDO
!
!     Apply permutations JPIV to the solution (RHS)
!
      CALL CLASWP(1,Rhs,Lda,1,N-1,Jpiv,-1)
!
!     End of CGESC2
!
      END SUBROUTINE CGESC2
