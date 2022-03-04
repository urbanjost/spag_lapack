!*==crot.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b CROT applies a plane rotation with real cosine and complex sine to a pair of complex vectors.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CROT + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/crot.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/crot.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/crot.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CROT( N, CX, INCX, CY, INCY, C, S )
!
!       .. Scalar Arguments ..
!       INTEGER            INCX, INCY, N
!       REAL               C
!       COMPLEX            S
!       ..
!       .. Array Arguments ..
!       COMPLEX            CX( * ), CY( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CROT   applies a plane rotation, where the cos (C) is real and the
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
!>          CX is COMPLEX array, dimension (N)
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
!>          CY is COMPLEX array, dimension (N)
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
!>          C is REAL
!> \endverbatim
!>
!> \param[in] S
!> \verbatim
!>          S is COMPLEX
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
!> \ingroup complexOTHERauxiliary
!
!  =====================================================================
      SUBROUTINE CROT(N,Cx,Incx,Cy,Incy,C,S)
      IMPLICIT NONE
!*--CROT107
!
!  -- LAPACK auxiliary routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Incx , Incy , N
      REAL C
      COMPLEX S
!     ..
!     .. Array Arguments ..
      COMPLEX Cx(*) , Cy(*)
!     ..
!
! =====================================================================
!
!     .. Local Scalars ..
      INTEGER i , ix , iy
      COMPLEX stemp
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC CONJG
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
            Cy(i) = C*Cy(i) - CONJG(S)*Cx(i)
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
         Cy(iy) = C*Cy(iy) - CONJG(S)*Cx(ix)
         Cx(ix) = stemp
         ix = ix + Incx
         iy = iy + Incy
      ENDDO
      RETURN
99999 END SUBROUTINE CROT
