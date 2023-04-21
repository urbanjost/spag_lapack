!*==zrot.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b ZROT applies a plane rotation with real cosine and complex sine to a pair of complex vectors.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZROT + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zrot.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zrot.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zrot.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZROT( N, CX, INCX, CY, INCY, C, S )
!
!       .. Scalar Arguments ..
!       INTEGER            INCX, INCY, N
!       DOUBLE PRECISION   C
!       COMPLEX*16         S
!       ..
!       .. Array Arguments ..
!       COMPLEX*16         CX( * ), CY( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZROT   applies a plane rotation, where the cos (C) is real and the
!> sin (S) is complex, and the vectors CX and CY are complex.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of elements in the vectors CX and CY.
!> \endverbatim
!>
!> \param[in,out] CX
!> \verbatim
!>          CX is COMPLEX*16 array, dimension (N)
!>          On input, the vector X.
!>          On output, CX is overwritten with C*X + S*Y.
!> \endverbatim
!>
!> \param[in] INCX
!> \verbatim
!>          INCX is INTEGER
!>          The increment between successive values of CY.  INCX <> 0.
!> \endverbatim
!>
!> \param[in,out] CY
!> \verbatim
!>          CY is COMPLEX*16 array, dimension (N)
!>          On input, the vector Y.
!>          On output, CY is overwritten with -CONJG(S)*X + C*Y.
!> \endverbatim
!>
!> \param[in] INCY
!> \verbatim
!>          INCY is INTEGER
!>          The increment between successive values of CY.  INCX <> 0.
!> \endverbatim
!>
!> \param[in] C
!> \verbatim
!>          C is DOUBLE PRECISION
!> \endverbatim
!>
!> \param[in] S
!> \verbatim
!>          S is COMPLEX*16
!>          C and S define a rotation
!>             [  C          S  ]
!>             [ -conjg(S)   C  ]
!>          where C*C + S*CONJG(S) = 1.0.
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
      SUBROUTINE ZROT(N,Cx,Incx,Cy,Incy,C,S)
      IMPLICIT NONE
!*--ZROT107
!
!  -- LAPACK auxiliary routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Incx , Incy , N
      DOUBLE PRECISION C
      COMPLEX*16 S
!     ..
!     .. Array Arguments ..
      COMPLEX*16 Cx(*) , Cy(*)
!     ..
!
! =====================================================================
!
!     .. Local Scalars ..
      INTEGER i , ix , iy
      COMPLEX*16 stemp
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC DCONJG
!     ..
!     .. Executable Statements ..
!
      IF ( N<=0 ) RETURN
      IF ( Incx==1 .AND. Incy==1 ) THEN
!
!     Code for both increments equal to 1
!
         DO i = 1 , N
            stemp = C*Cx(i) + S*Cy(i)
            Cy(i) = C*Cy(i) - DCONJG(S)*Cx(i)
            Cx(i) = stemp
         ENDDO
         GOTO 99999
      ENDIF
!
!     Code for unequal increments or equal increments not equal to 1
!
      ix = 1
      iy = 1
      IF ( Incx<0 ) ix = (-N+1)*Incx + 1
      IF ( Incy<0 ) iy = (-N+1)*Incy + 1
      DO i = 1 , N
         stemp = C*Cx(ix) + S*Cy(iy)
         Cy(iy) = C*Cy(iy) - DCONJG(S)*Cx(ix)
         Cx(ix) = stemp
         ix = ix + Incx
         iy = iy + Incy
      ENDDO
      RETURN
99999 END SUBROUTINE ZROT
