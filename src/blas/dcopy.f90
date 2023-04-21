!*==dcopy.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b DCOPY
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DCOPY(N,DX,INCX,DY,INCY)
!
!       .. Scalar Arguments ..
!       INTEGER INCX,INCY,N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION DX(*),DY(*)
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    DCOPY copies a vector, x, to a vector, y.
!>    uses unrolled loops for increments equal to 1.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>         number of elements in input vector(s)
!> \endverbatim
!>
!> \param[in] DX
!> \verbatim
!>          DX is DOUBLE PRECISION array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
!> \endverbatim
!>
!> \param[in] INCX
!> \verbatim
!>          INCX is INTEGER
!>         storage spacing between elements of DX
!> \endverbatim
!>
!> \param[out] DY
!> \verbatim
!>          DY is DOUBLE PRECISION array, dimension ( 1 + ( N - 1 )*abs( INCY ) )
!> \endverbatim
!>
!> \param[in] INCY
!> \verbatim
!>          INCY is INTEGER
!>         storage spacing between elements of DY
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
!> \date November 2017
!
!> \ingroup double_blas_level1
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>     jack dongarra, linpack, 3/11/78.
!>     modified 12/3/93, array(1) declarations changed to array(*)
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE DCOPY(N,Dx,Incx,Dy,Incy)
      IMPLICIT NONE
!*--DCOPY86
!
!  -- Reference BLAS level1 routine (version 3.8.0) --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2017
!
!     .. Scalar Arguments ..
      INTEGER Incx , Incy , N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION Dx(*) , Dy(*)
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
      INTEGER i , ix , iy , m , mp1
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MOD
!     ..
      IF ( N<=0 ) RETURN
      IF ( Incx==1 .AND. Incy==1 ) THEN
!
!        code for both increments equal to 1
!
!
!        clean-up loop
!
         m = MOD(N,7)
         IF ( m/=0 ) THEN
            DO i = 1 , m
               Dy(i) = Dx(i)
            ENDDO
            IF ( N<7 ) RETURN
         ENDIF
         mp1 = m + 1
         DO i = mp1 , N , 7
            Dy(i) = Dx(i)
            Dy(i+1) = Dx(i+1)
            Dy(i+2) = Dx(i+2)
            Dy(i+3) = Dx(i+3)
            Dy(i+4) = Dx(i+4)
            Dy(i+5) = Dx(i+5)
            Dy(i+6) = Dx(i+6)
         ENDDO
      ELSE
!
!        code for unequal increments or equal increments
!          not equal to 1
!
         ix = 1
         iy = 1
         IF ( Incx<0 ) ix = (-N+1)*Incx + 1
         IF ( Incy<0 ) iy = (-N+1)*Incy + 1
         DO i = 1 , N
            Dy(iy) = Dx(ix)
            ix = ix + Incx
            iy = iy + Incy
         ENDDO
      ENDIF
      END SUBROUTINE DCOPY
