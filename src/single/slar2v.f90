!*==slar2v.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b SLAR2V applies a vector of plane rotations with real cosines and real sines from both sides to a sequence of 2-by-2 symmetric/Hermitian matrices.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SLAR2V + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/slar2v.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/slar2v.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/slar2v.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SLAR2V( N, X, Y, Z, INCX, C, S, INCC )
!
!       .. Scalar Arguments ..
!       INTEGER            INCC, INCX, N
!       ..
!       .. Array Arguments ..
!       REAL               C( * ), S( * ), X( * ), Y( * ), Z( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SLAR2V applies a vector of real plane rotations from both sides to
!> a sequence of 2-by-2 real symmetric matrices, defined by the elements
!> of the vectors x, y and z. For i = 1,2,...,n
!>
!>    ( x(i)  z(i) ) := (  c(i)  s(i) ) ( x(i)  z(i) ) ( c(i) -s(i) )
!>    ( z(i)  y(i) )    ( -s(i)  c(i) ) ( z(i)  y(i) ) ( s(i)  c(i) )
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
!>          X is REAL array,
!>                         dimension (1+(N-1)*INCX)
!>          The vector x.
!> \endverbatim
!>
!> \param[in,out] Y
!> \verbatim
!>          Y is REAL array,
!>                         dimension (1+(N-1)*INCX)
!>          The vector y.
!> \endverbatim
!>
!> \param[in,out] Z
!> \verbatim
!>          Z is REAL array,
!>                         dimension (1+(N-1)*INCX)
!>          The vector z.
!> \endverbatim
!>
!> \param[in] INCX
!> \verbatim
!>          INCX is INTEGER
!>          The increment between elements of X, Y and Z. INCX > 0.
!> \endverbatim
!>
!> \param[in] C
!> \verbatim
!>          C is REAL array, dimension (1+(N-1)*INCC)
!>          The cosines of the plane rotations.
!> \endverbatim
!>
!> \param[in] S
!> \verbatim
!>          S is REAL array, dimension (1+(N-1)*INCC)
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
!> \ingroup realOTHERauxiliary
!
!  =====================================================================
      SUBROUTINE SLAR2V(N,X,Y,Z,Incx,C,S,Incc)
      IMPLICIT NONE
!*--SLAR2V114
!
!  -- LAPACK auxiliary routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Incc , Incx , N
!     ..
!     .. Array Arguments ..
      REAL C(*) , S(*) , X(*) , Y(*) , Z(*)
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
      INTEGER i , ic , ix
      REAL ci , si , t1 , t2 , t3 , t4 , t5 , t6 , xi , yi , zi
!     ..
!     .. Executable Statements ..
!
      ix = 1
      ic = 1
      DO i = 1 , N
         xi = X(ix)
         yi = Y(ix)
         zi = Z(ix)
         ci = C(ic)
         si = S(ic)
         t1 = si*zi
         t2 = ci*zi
         t3 = t2 - si*xi
         t4 = t2 + si*yi
         t5 = ci*xi + t1
         t6 = ci*yi - t1
         X(ix) = ci*t5 + si*t4
         Y(ix) = ci*t6 - si*t3
         Z(ix) = ci*t4 - si*t5
         ix = ix + Incx
         ic = ic + Incc
      ENDDO
!
!     End of SLAR2V
!
      END SUBROUTINE SLAR2V
