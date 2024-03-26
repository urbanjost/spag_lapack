!*==zlaipd.f90  processed by SPAG 7.51RB at 17:35 on  4 Mar 2022
!> \brief \b zlaipd
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE zlaipd( N, A, INDA, VINDA )
!
!       .. Scalar Arguments ..
!       INTEGER            INDA, N, VINDA
!       ..
!       .. Array Arguments ..
!       COMPLEX*16         A( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> zlaipd sets the imaginary part of the diagonal elements of a complex
!> matrix A to a large value.  This is used to test LAPACK routines for
!> complex Hermitian matrices, which are not supposed to access or use
!> the imaginary parts of the diagonals.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>         The number of diagonal elements of A.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension
!>                        (1+(N-1)*INDA+(N-2)*VINDA)
!>         On entry, the complex (Hermitian) matrix A.
!>         On exit, the imaginary parts of the diagonal elements are set
!>         to BIGNUM = EPS / SAFMIN, where EPS is the machine epsilon and
!>         SAFMIN is the safe minimum.
!> \endverbatim
!>
!> \param[in] INDA
!> \verbatim
!>          INDA is INTEGER
!>         The increment between A(1) and the next diagonal element of A.
!>         Typical values are
!>         = LDA+1:  square matrices with leading dimension LDA
!>         = 2:  packed upper triangular matrix, starting at A(1,1)
!>         = N:  packed lower triangular matrix, starting at A(1,1)
!> \endverbatim
!>
!> \param[in] VINDA
!> \verbatim
!>          VINDA is INTEGER
!>         The change in the diagonal increment between columns of A.
!>         Typical values are
!>         = 0:  no change, the row and column increments in A are fixed
!>         = 1:  packed upper triangular matrix
!>         = -1:  packed lower triangular matrix
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
!> \ingroup complex16_lin
!
!  =====================================================================
      SUBROUTINE zlaipd(N,A,Inda,Vinda)
      IMPLICIT NONE
!*--ZLAIPD87
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Inda , N , Vinda
!     ..
!     .. Array Arguments ..
      COMPLEX*16 A(*)
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
      INTEGER i , ia , ixa
      DOUBLE PRECISION bignum
!     ..
!     .. External Functions ..
      DOUBLE PRECISION DLAMCH
EXTERNAL    DLAMCH
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC DBLE , DCMPLX
!     ..
!     .. Executable Statements ..
!
      bignum = DLAMCH('Epsilon')/DLAMCH('Safe minimum')
      ia = 1
      ixa = Inda
      DO i = 1 , N
         A(ia) = DCMPLX(DBLE(A(ia)),bignum)
         ia = ia + ixa
         ixa = ixa + Vinda
      ENDDO
      END SUBROUTINE zlaipd
