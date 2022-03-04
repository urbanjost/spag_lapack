!*==dswap.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b DSWAP
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DSWAP(N,DX,INCX,DY,INCY)
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
!>    DSWAP interchanges two vectors.
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
!> \param[in,out] DX
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
!> \param[in,out] DY
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
      SUBROUTINE DSWAP(N,Dx,Incx,Dy,Incy)
      USE F77KINDS                        
      IMPLICIT NONE
!*--DSWAP87
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Dx
      INTEGER , INTENT(IN) :: Incx
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Dy
      INTEGER , INTENT(IN) :: Incy
!
! Local variable declarations rewritten by SPAG
!
      REAL(R8KIND) :: dtemp
      INTEGER :: i , ix , iy , m , mp1
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
      IF ( N<=0 ) RETURN
      IF ( Incx==1 .AND. Incy==1 ) THEN
!
!       code for both increments equal to 1
!
!
!       clean-up loop
!
         m = MOD(N,3)
         IF ( m/=0 ) THEN
            DO i = 1 , m
               dtemp = Dx(i)
               Dx(i) = Dy(i)
               Dy(i) = dtemp
            ENDDO
            IF ( N<3 ) RETURN
         ENDIF
         mp1 = m + 1
         DO i = mp1 , N , 3
            dtemp = Dx(i)
            Dx(i) = Dy(i)
            Dy(i) = dtemp
            dtemp = Dx(i+1)
            Dx(i+1) = Dy(i+1)
            Dy(i+1) = dtemp
            dtemp = Dx(i+2)
            Dx(i+2) = Dy(i+2)
            Dy(i+2) = dtemp
         ENDDO
      ELSE
!
!       code for unequal increments or equal increments not equal
!         to 1
!
         ix = 1
         iy = 1
         IF ( Incx<0 ) ix = (-N+1)*Incx + 1
         IF ( Incy<0 ) iy = (-N+1)*Incy + 1
         DO i = 1 , N
            dtemp = Dx(ix)
            Dx(ix) = Dy(iy)
            Dy(iy) = dtemp
            ix = ix + Incx
            iy = iy + Incy
         ENDDO
      ENDIF
      END SUBROUTINE DSWAP
