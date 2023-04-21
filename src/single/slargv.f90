!*==slargv.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b SLARGV generates a vector of plane rotations with real cosines and real sines.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SLARGV + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slargv.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slargv.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slargv.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SLARGV( N, X, INCX, Y, INCY, C, INCC )
!
!       .. Scalar Arguments ..
!       INTEGER            INCC, INCX, INCY, N
!       ..
!       .. Array Arguments ..
!       REAL               C( * ), X( * ), Y( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SLARGV generates a vector of real plane rotations, determined by
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
!>          X is REAL array,
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
!>          Y is REAL array,
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
!>          C is REAL array, dimension (1+(N-1)*INCC)
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
!> \ingroup realOTHERauxiliary
!
!  =====================================================================
      SUBROUTINE SLARGV(N,X,Incx,Y,Incy,C,Incc)
      IMPLICIT NONE
!*--SLARGV108
!
!  -- LAPACK auxiliary routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Incc , Incx , Incy , N
!     ..
!     .. Array Arguments ..
      REAL C(*) , X(*) , Y(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ZERO , ONE
      PARAMETER (ZERO=0.0E+0,ONE=1.0E+0)
!     ..
!     .. Local Scalars ..
      INTEGER i , ic , ix , iy
      REAL f , g , t , tt
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , SQRT
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
!     End of SLARGV
!
      END SUBROUTINE SLARGV
