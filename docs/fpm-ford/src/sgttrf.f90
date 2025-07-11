!*==sgttrf.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b SGTTRF
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SGTTRF + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sgttrf.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sgttrf.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sgttrf.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SGTTRF( N, DL, D, DU, DU2, IPIV, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, N
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * )
!       REAL               D( * ), DL( * ), DU( * ), DU2( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SGTTRF computes an LU factorization of a real tridiagonal matrix A
!> using elimination with partial pivoting and row interchanges.
!>
!> The factorization has the form
!>    A = L * U
!> where L is a product of permutation and unit lower bidiagonal
!> matrices and U is upper triangular with nonzeros in only the main
!> diagonal and first two superdiagonals.
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
!> \param[in,out] DL
!> \verbatim
!>          DL is REAL array, dimension (N-1)
!>          On entry, DL must contain the (n-1) sub-diagonal elements of
!>          A.
!>
!>          On exit, DL is overwritten by the (n-1) multipliers that
!>          define the matrix L from the LU factorization of A.
!> \endverbatim
!>
!> \param[in,out] D
!> \verbatim
!>          D is REAL array, dimension (N)
!>          On entry, D must contain the diagonal elements of A.
!>
!>          On exit, D is overwritten by the n diagonal elements of the
!>          upper triangular matrix U from the LU factorization of A.
!> \endverbatim
!>
!> \param[in,out] DU
!> \verbatim
!>          DU is REAL array, dimension (N-1)
!>          On entry, DU must contain the (n-1) super-diagonal elements
!>          of A.
!>
!>          On exit, DU is overwritten by the (n-1) elements of the first
!>          super-diagonal of U.
!> \endverbatim
!>
!> \param[out] DU2
!> \verbatim
!>          DU2 is REAL array, dimension (N-2)
!>          On exit, DU2 is overwritten by the (n-2) elements of the
!>          second super-diagonal of U.
!> \endverbatim
!>
!> \param[out] IPIV
!> \verbatim
!>          IPIV is INTEGER array, dimension (N)
!>          The pivot indices; for 1 <= i <= n, row i of the matrix was
!>          interchanged with row IPIV(i).  IPIV(i) will always be either
!>          i or i+1; IPIV(i) = i indicates a row interchange was not
!>          required.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -k, the k-th argument had an illegal value
!>          > 0:  if INFO = k, U(k,k) is exactly zero. The factorization
!>                has been completed, but the factor U is exactly
!>                singular, and division by zero will occur if it is used
!>                to solve a system of equations.
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
!> \ingroup realGTcomputational
!
!  =====================================================================
      SUBROUTINE SGTTRF(N,Dl,D,Du,Du2,Ipiv,Info)
      IMPLICIT NONE
!*--SGTTRF128
!
!  -- LAPACK computational routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Info , N
!     ..
!     .. Array Arguments ..
      INTEGER Ipiv(*)
      REAL D(*) , Dl(*) , Du(*) , Du2(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ZERO
      PARAMETER (ZERO=0.0E+0)
!     ..
!     .. Local Scalars ..
      INTEGER i
      REAL fact , temp
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS
!     ..
!     .. External Subroutines ..
      EXTERNAL XERBLA
!     ..
!     .. Executable Statements ..
!
      Info = 0
      IF ( N<0 ) THEN
         Info = -1
         CALL XERBLA('SGTTRF',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( N==0 ) RETURN
!
!     Initialize IPIV(i) = i and DU2(I) = 0
!
      DO i = 1 , N
         Ipiv(i) = i
      ENDDO
      DO i = 1 , N - 2
         Du2(i) = ZERO
      ENDDO
!
      DO i = 1 , N - 2
         IF ( ABS(D(i))<ABS(Dl(i)) ) THEN
!
!           Interchange rows I and I+1, eliminate DL(I)
!
            fact = D(i)/Dl(i)
            D(i) = Dl(i)
            Dl(i) = fact
            temp = Du(i)
            Du(i) = D(i+1)
            D(i+1) = temp - fact*D(i+1)
            Du2(i) = Du(i+1)
            Du(i+1) = -fact*Du(i+1)
            Ipiv(i) = i + 1
!
!           No row interchange required, eliminate DL(I)
!
         ELSEIF ( D(i)/=ZERO ) THEN
            fact = Dl(i)/D(i)
            Dl(i) = fact
            D(i+1) = D(i+1) - fact*Du(i)
         ENDIF
      ENDDO
      IF ( N>1 ) THEN
         i = N - 1
         IF ( ABS(D(i))<ABS(Dl(i)) ) THEN
            fact = D(i)/Dl(i)
            D(i) = Dl(i)
            Dl(i) = fact
            temp = Du(i)
            Du(i) = D(i+1)
            D(i+1) = temp - fact*D(i+1)
            Ipiv(i) = i + 1
         ELSEIF ( D(i)/=ZERO ) THEN
            fact = Dl(i)/D(i)
            Dl(i) = fact
            D(i+1) = D(i+1) - fact*Du(i)
         ENDIF
      ENDIF
!
!     Check for a zero on the diagonal of U.
!
      DO i = 1 , N
         IF ( D(i)==ZERO ) THEN
            Info = i
            EXIT
         ENDIF
      ENDDO
!
!
!     End of SGTTRF
!
      END SUBROUTINE SGTTRF
