!*==izmax1.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b IZMAX1 finds the index of the first vector element of maximum absolute value.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download IZMAX1 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/izmax1.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/izmax1.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/izmax1.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       INTEGER          FUNCTION IZMAX1( N, ZX, INCX )
!
!       .. Scalar Arguments ..
!       INTEGER            INCX, N
!       ..
!       .. Array Arguments ..
!       COMPLEX*16         ZX( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> IZMAX1 finds the index of the first vector element of maximum absolute value.
!>
!> Based on IZAMAX from Level 1 BLAS.
!> The change is to use the 'genuine' absolute value.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of elements in the vector ZX.
!> \endverbatim
!>
!> \param[in] ZX
!> \verbatim
!>          ZX is COMPLEX*16 array, dimension (N)
!>          The vector ZX. The IZMAX1 function returns the index of its first
!>          element of maximum absolute value.
!> \endverbatim
!>
!> \param[in] INCX
!> \verbatim
!>          INCX is INTEGER
!>          The spacing between successive values of ZX.  INCX >= 1.
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
!> \date February 2014
!
!> \ingroup complexOTHERauxiliary
!
!> \par Contributors:
!  ==================
!>
!> Nick Higham for use with ZLACON.
!
!  =====================================================================
      INTEGER FUNCTION IZMAX1(N,Zx,Incx)
      IMPLICIT NONE
!*--IZMAX185
!
!  -- LAPACK auxiliary routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     February 2014
!
!     .. Scalar Arguments ..
      INTEGER Incx , N
!     ..
!     .. Array Arguments ..
      COMPLEX*16 Zx(*)
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
      DOUBLE PRECISION dmax
      INTEGER i , ix
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS
!     ..
!     .. Executable Statements ..
!
      IZMAX1 = 0
      IF ( N<1 .OR. Incx<=0 ) RETURN
      IZMAX1 = 1
      IF ( N==1 ) RETURN
      IF ( Incx==1 ) THEN
!
!        code for increment equal to 1
!
         dmax = ABS(Zx(1))
         DO i = 2 , N
            IF ( ABS(Zx(i))>dmax ) THEN
               IZMAX1 = i
               dmax = ABS(Zx(i))
            ENDIF
         ENDDO
      ELSE
!
!        code for increment not equal to 1
!
         ix = 1
         dmax = ABS(Zx(1))
         ix = ix + Incx
         DO i = 2 , N
            IF ( ABS(Zx(ix))>dmax ) THEN
               IZMAX1 = i
               dmax = ABS(Zx(ix))
            ENDIF
            ix = ix + Incx
         ENDDO
      ENDIF
!
!     End of IZMAX1
!
      END FUNCTION IZMAX1
