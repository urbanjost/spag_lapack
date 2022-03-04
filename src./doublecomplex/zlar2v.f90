!*==zlar2v.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b ZLAR2V applies a vector of plane rotations with real cosines and complex sines from both sides to a sequence of 2-by-2 symmetric/Hermitian matrices.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZLAR2V + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlar2v.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlar2v.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlar2v.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZLAR2V( N, X, Y, Z, INCX, C, S, INCC )
!
!       .. Scalar Arguments ..
!       INTEGER            INCC, INCX, N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   C( * )
!       COMPLEX*16         S( * ), X( * ), Y( * ), Z( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZLAR2V applies a vector of complex plane rotations with real cosines
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
!>          X is COMPLEX*16 array, dimension (1+(N-1)*INCX)
!>          The vector x; the elements of x are assumed to be real.
!> \endverbatim
!>
!> \param[in,out] Y
!> \verbatim
!>          Y is COMPLEX*16 array, dimension (1+(N-1)*INCX)
!>          The vector y; the elements of y are assumed to be real.
!> \endverbatim
!>
!> \param[in,out] Z
!> \verbatim
!>          Z is COMPLEX*16 array, dimension (1+(N-1)*INCX)
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
!>          C is DOUBLE PRECISION array, dimension (1+(N-1)*INCC)
!>          The cosines of the plane rotations.
!> \endverbatim
!>
!> \param[in] S
!> \verbatim
!>          S is COMPLEX*16 array, dimension (1+(N-1)*INCC)
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
!> \ingroup complex16OTHERauxiliary
!
!  =====================================================================
      SUBROUTINE ZLAR2V(N,X,Y,Z,Incx,C,S,Incc)
      USE F77KINDS                        
      IMPLICIT NONE
!*--ZLAR2V116
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: X
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Y
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Z
      INTEGER , INTENT(IN) :: Incx
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: C
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(*) :: S
      INTEGER , INTENT(IN) :: Incc
!
! Local variable declarations rewritten by SPAG
!
      REAL(R8KIND) :: ci , sii , sir , t1i , t1r , t5 , t6 , xi , yi ,  &
     &                zii , zir
      INTEGER :: i , ic , ix
      COMPLEX(CX16KIND) :: si , t2 , t3 , t4 , zi
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
!     .. Intrinsic Functions ..
!     ..
!     .. Executable Statements ..
!
      ix = 1
      ic = 1
      DO i = 1 , N
         xi = DBLE(X(ix))
         yi = DBLE(Y(ix))
         zi = Z(ix)
         zir = DBLE(zi)
         zii = DIMAG(zi)
         ci = C(ic)
         si = S(ic)
         sir = DBLE(si)
         sii = DIMAG(si)
         t1r = sir*zir - sii*zii
         t1i = sir*zii + sii*zir
         t2 = ci*zi
         t3 = t2 - DCONJG(si)*xi
         t4 = DCONJG(t2) + si*yi
         t5 = ci*xi + t1r
         t6 = ci*yi - t1r
         X(ix) = ci*t5 + (sir*DBLE(t4)+sii*DIMAG(t4))
         Y(ix) = ci*t6 - (sir*DBLE(t3)-sii*DIMAG(t3))
         Z(ix) = ci*t3 + DCONJG(si)*DCMPLX(t6,t1i)
         ix = ix + Incx
         ic = ic + Incc
      ENDDO
!
!     End of ZLAR2V
!
      END SUBROUTINE ZLAR2V
