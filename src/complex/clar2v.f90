!*==clar2v.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b CLAR2V applies a vector of plane rotations with real cosines and complex sines from both sides to a sequence of 2-by-2 symmetric/Hermitian matrices.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CLAR2V + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/clar2v.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/clar2v.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/clar2v.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CLAR2V( N, X, Y, Z, INCX, C, S, INCC )
!
!       .. Scalar Arguments ..
!       INTEGER            INCC, INCX, N
!       ..
!       .. Array Arguments ..
!       REAL               C( * )
!       COMPLEX            S( * ), X( * ), Y( * ), Z( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CLAR2V applies a vector of complex plane rotations with real cosines
!> from both sides to a sequence of 2-by-2 complex Hermitian matrices,
!> defined by the elements of the vectors x, y and z. For i = 1,2,...,n
!>
!>    (       x(i)  z(i) ) :=
!>    ( conjg(z(i)) y(i) )
!>
!>      (  c(i) conjg(s(i)) ) (       x(i)  z(i) ) ( c(i) -conjg(s(i)) )
!>      ( -s(i)       c(i)  ) ( conjg(z(i)) y(i) ) ( s(i)        c(i)  )
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
!>          X is COMPLEX array, dimension (1+(N-1)*INCX)
!>          The vector x; the elements of x are assumed to be real.
!> \endverbatim
!>
!> \param[in,out] Y
!> \verbatim
!>          Y is COMPLEX array, dimension (1+(N-1)*INCX)
!>          The vector y; the elements of y are assumed to be real.
!> \endverbatim
!>
!> \param[in,out] Z
!> \verbatim
!>          Z is COMPLEX array, dimension (1+(N-1)*INCX)
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
!>          S is COMPLEX array, dimension (1+(N-1)*INCC)
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
!> \ingroup complexOTHERauxiliary
!
!  =====================================================================
      SUBROUTINE CLAR2V(N,X,Y,Z,Incx,C,S,Incc)
      IMPLICIT NONE
!*--CLAR2V115
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
      REAL C(*)
      COMPLEX S(*) , X(*) , Y(*) , Z(*)
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
      INTEGER i , ic , ix
      REAL ci , sii , sir , t1i , t1r , t5 , t6 , xi , yi , zii , zir
      COMPLEX si , t2 , t3 , t4 , zi
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC AIMAG , CMPLX , CONJG , REAL
!     ..
!     .. Executable Statements ..
!
      ix = 1
      ic = 1
      DO i = 1 , N
         xi = REAL(X(ix))
         yi = REAL(Y(ix))
         zi = Z(ix)
         zir = REAL(zi)
         zii = AIMAG(zi)
         ci = C(ic)
         si = S(ic)
         sir = REAL(si)
         sii = AIMAG(si)
         t1r = sir*zir - sii*zii
         t1i = sir*zii + sii*zir
         t2 = ci*zi
         t3 = t2 - CONJG(si)*xi
         t4 = CONJG(t2) + si*yi
         t5 = ci*xi + t1r
         t6 = ci*yi - t1r
         X(ix) = ci*t5 + (sir*REAL(t4)+sii*AIMAG(t4))
         Y(ix) = ci*t6 - (sir*REAL(t3)-sii*AIMAG(t3))
         Z(ix) = ci*t3 + CONJG(si)*CMPLX(t6,t1i)
         ix = ix + Incx
         ic = ic + Incc
      ENDDO
!
!     End of CLAR2V
!
      END SUBROUTINE CLAR2V
