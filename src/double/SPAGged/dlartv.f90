!*==dlartv.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
 
!> \brief \b DLARTV applies a vector of plane rotations with real cosines and real sines to the elements of a pair of vectors.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DLARTV + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlartv.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlartv.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlartv.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLARTV( N, X, INCX, Y, INCY, C, S, INCC )
!
!       .. Scalar Arguments ..
!       INTEGER            INCC, INCX, INCY, N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   C( * ), S( * ), X( * ), Y( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLARTV applies a vector of real plane rotations to elements of the
!> real vectors x and y. For i = 1,2,...,n
!>
!>    ( x(i) ) := (  c(i)  s(i) ) ( x(i) )
!>    ( y(i) )    ( -s(i)  c(i) ) ( y(i) )
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of plane rotations to be applied.
!> \endverbatim
!>
!> \param[in,out] X
!> \verbatim
!>          X is DOUBLE PRECISION array,
!>                         dimension (1+(N-1)*INCX)
!>          The vector x.
!> \endverbatim
!>
!> \param[in] INCX
!> \verbatim
!>          INCX is INTEGER
!>          The increment between elements of X. INCX > 0.
!> \endverbatim
!>
!> \param[in,out] Y
!> \verbatim
!>          Y is DOUBLE PRECISION array,
!>                         dimension (1+(N-1)*INCY)
!>          The vector y.
!> \endverbatim
!>
!> \param[in] INCY
!> \verbatim
!>          INCY is INTEGER
!>          The increment between elements of Y. INCY > 0.
!> \endverbatim
!>
!> \param[in] C
!> \verbatim
!>          C is DOUBLE PRECISION array, dimension (1+(N-1)*INCC)
!>          The cosines of the plane rotations.
!> \endverbatim
!>
!> \param[in] S
!> \verbatim
!>          S is DOUBLE PRECISION array, dimension (1+(N-1)*INCC)
!>          The sines of the plane rotations.
!> \endverbatim
!>
!> \param[in] INCC
!> \verbatim
!>          INCC is INTEGER
!>          The increment between elements of C and S. INCC > 0.
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
!> \ingroup doubleOTHERauxiliary
!
!  =====================================================================
      SUBROUTINE DLARTV(N,X,Incx,Y,Incy,C,S,Incc)
      USE F77KINDS                        
      IMPLICIT NONE
!*--DLARTV114
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: X
      INTEGER , INTENT(IN) :: Incx
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Y
      INTEGER , INTENT(IN) :: Incy
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: C
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: S
      INTEGER , INTENT(IN) :: Incc
!
! Local variable declarations rewritten by SPAG
!
      INTEGER :: i , ic , ix , iy
      REAL(R8KIND) :: xi , yi
!
! End of declarations rewritten by SPAG
!
!     ..
!     .. Array Arguments ..
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
!     ..
!     .. Executable Statements ..
!
      ix = 1
      iy = 1
      ic = 1
      DO i = 1 , N
         xi = X(ix)
         yi = Y(iy)
         X(ix) = C(ic)*xi + S(ic)*yi
         Y(iy) = C(ic)*yi - S(ic)*xi
         ix = ix + Incx
         iy = iy + Incy
         ic = ic + Incc
      ENDDO
!
!     End of DLARTV
!
      END SUBROUTINE DLARTV
