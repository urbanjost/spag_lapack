!*==icopy.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b ICOPY
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ICOPY( N, SX, INCX, SY, INCY )
!
!       .. Scalar Arguments ..
!       INTEGER            INCX, INCY, N
!       ..
!       .. Array Arguments ..
!       INTEGER            SX( * ), SY( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ICOPY copies an integer vector x to an integer vector y.
!> Uses unrolled loops for increments equal to 1.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The length of the vectors SX and SY.
!> \endverbatim
!>
!> \param[in] SX
!> \verbatim
!>          SX is INTEGER array, dimension (1+(N-1)*abs(INCX))
!>          The vector X.
!> \endverbatim
!>
!> \param[in] INCX
!> \verbatim
!>          INCX is INTEGER
!>          The spacing between consecutive elements of SX.
!> \endverbatim
!>
!> \param[out] SY
!> \verbatim
!>          SY is INTEGER array, dimension (1+(N-1)*abs(INCY))
!>          The vector Y.
!> \endverbatim
!>
!> \param[in] INCY
!> \verbatim
!>          INCY is INTEGER
!>          The spacing between consecutive elements of SY.
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
!> \ingroup aux_lin
!
!  =====================================================================
      SUBROUTINE ICOPY(N,Sx,Incx,Sy,Incy)
      IMPLICIT NONE
!*--ICOPY79
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Incx , Incy , N
!     ..
!     .. Array Arguments ..
      INTEGER Sx(*) , Sy(*)
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
!     .. Executable Statements ..
!
      IF ( N<=0 ) RETURN
      IF ( Incx==1 .AND. Incy==1 ) THEN
!
!     Code for both increments equal to 1
!
!     Clean-up loop
!
         m = MOD(N,7)
         IF ( m/=0 ) THEN
            DO i = 1 , m
               Sy(i) = Sx(i)
            ENDDO
            IF ( N<7 ) RETURN
         ENDIF
         mp1 = m + 1
         DO i = mp1 , N , 7
            Sy(i) = Sx(i)
            Sy(i+1) = Sx(i+1)
            Sy(i+2) = Sx(i+2)
            Sy(i+3) = Sx(i+3)
            Sy(i+4) = Sx(i+4)
            Sy(i+5) = Sx(i+5)
            Sy(i+6) = Sx(i+6)
         ENDDO
      ELSE
!
!     Code for unequal increments or equal increments not equal to 1
!
         ix = 1
         iy = 1
         IF ( Incx<0 ) ix = (-N+1)*Incx + 1
         IF ( Incy<0 ) iy = (-N+1)*Incy + 1
         DO i = 1 , N
            Sy(iy) = Sx(ix)
            ix = ix + Incx
            iy = iy + Incy
         ENDDO
         RETURN
      ENDIF
!
!     End of ICOPY
!
      END SUBROUTINE ICOPY
