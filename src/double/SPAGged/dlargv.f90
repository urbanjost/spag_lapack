!*==dlargv.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b DLARGV generates a vector of plane rotations with real cosines and real sines.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DLARGV + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlargv.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlargv.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlargv.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLARGV( N, X, INCX, Y, INCY, C, INCC )
!
!       .. Scalar Arguments ..
!       INTEGER            INCC, INCX, INCY, N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   C( * ), X( * ), Y( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLARGV generates a vector of real plane rotations, determined by
!> elements of the real vectors x and y. For i = 1,2,...,n
!>
!>    (  c(i)  s(i) ) ( x(i) ) = ( a(i) )
!>    ( -s(i)  c(i) ) ( y(i) ) = (   0  )
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of plane rotations to be generated.
!> \endverbatim
!>
!> \param[in,out] X
!> \verbatim
!>          X is DOUBLE PRECISION array,
!>                         dimension (1+(N-1)*INCX)
!>          On entry, the vector x.
!>          On exit, x(i) is overwritten by a(i), for i = 1,...,n.
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
!>          On entry, the vector y.
!>          On exit, the sines of the plane rotations.
!> \endverbatim
!>
!> \param[in] INCY
!> \verbatim
!>          INCY is INTEGER
!>          The increment between elements of Y. INCY > 0.
!> \endverbatim
!>
!> \param[out] C
!> \verbatim
!>          C is DOUBLE PRECISION array, dimension (1+(N-1)*INCC)
!>          The cosines of the plane rotations.
!> \endverbatim
!>
!> \param[in] INCC
!> \verbatim
!>          INCC is INTEGER
!>          The increment between elements of C. INCC > 0.
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
      SUBROUTINE DLARGV(N,X,Incx,Y,Incy,C,Incc)
      USE F77KINDS                        
      IMPLICIT NONE
!*--DLARGV109
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0 , ONE = 1.0D+0
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: X
      INTEGER , INTENT(IN) :: Incx
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Y
      INTEGER , INTENT(IN) :: Incy
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: C
      INTEGER , INTENT(IN) :: Incc
!
! Local variable declarations rewritten by SPAG
!
      REAL(R8KIND) :: f , g , t , tt
      INTEGER :: i , ic , ix , iy
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
!     .. Intrinsic Functions ..
!     ..
!     .. Executable Statements ..
!
      ix = 1
      iy = 1
      ic = 1
      DO i = 1 , N
         f = X(ix)
         g = Y(iy)
         IF ( g==ZERO ) THEN
            C(ic) = ONE
         ELSEIF ( f==ZERO ) THEN
            C(ic) = ZERO
            Y(iy) = ONE
            X(ix) = g
         ELSEIF ( ABS(f)>ABS(g) ) THEN
            t = g/f
            tt = SQRT(ONE+t*t)
            C(ic) = ONE/tt
            Y(iy) = t*C(ic)
            X(ix) = f*tt
         ELSE
            t = f/g
            tt = SQRT(ONE+t*t)
            Y(iy) = ONE/tt
            C(ic) = t*Y(iy)
            X(ix) = g*tt
         ENDIF
         ic = ic + Incc
         iy = iy + Incy
         ix = ix + Incx
      ENDDO
!
!     End of DLARGV
!
      END SUBROUTINE DLARGV
