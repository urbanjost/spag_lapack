!*==saxpy.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b SAXPY
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE SAXPY(N,SA,SX,INCX,SY,INCY)
!
!       .. Scalar Arguments ..
!       REAL SA
!       INTEGER INCX,INCY,N
!       ..
!       .. Array Arguments ..
!       REAL SX(*),SY(*)
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    SAXPY constant times a vector plus a vector.
!>    uses unrolled loops for increments equal to one.
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
!> \param[in] SA
!> \verbatim
!>          SA is REAL
!>           On entry, SA specifies the scalar alpha.
!> \endverbatim
!>
!> \param[in] SX
!> \verbatim
!>          SX is REAL array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
!> \endverbatim
!>
!> \param[in] INCX
!> \verbatim
!>          INCX is INTEGER
!>         storage spacing between elements of SX
!> \endverbatim
!>
!> \param[in,out] SY
!> \verbatim
!>          SY is REAL array, dimension ( 1 + ( N - 1 )*abs( INCY ) )
!> \endverbatim
!>
!> \param[in] INCY
!> \verbatim
!>          INCY is INTEGER
!>         storage spacing between elements of SY
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
!> \ingroup single_blas_level1
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
      SUBROUTINE SAXPY(N,Sa,Sx,Incx,Sy,Incy)
      IMPLICIT NONE
!*--SAXPY93
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(IN) :: Sa
      REAL , INTENT(IN) , DIMENSION(*) :: Sx
      INTEGER , INTENT(IN) :: Incx
      REAL , INTENT(INOUT) , DIMENSION(*) :: Sy
      INTEGER , INTENT(IN) :: Incy
!
! Local variable declarations rewritten by SPAG
!
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
      IF ( Sa==0.0 ) RETURN
      IF ( Incx==1 .AND. Incy==1 ) THEN
!
!        code for both increments equal to 1
!
!
!        clean-up loop
!
         m = MOD(N,4)
         IF ( m/=0 ) THEN
            DO i = 1 , m
               Sy(i) = Sy(i) + Sa*Sx(i)
            ENDDO
         ENDIF
         IF ( N<4 ) RETURN
         mp1 = m + 1
         DO i = mp1 , N , 4
            Sy(i) = Sy(i) + Sa*Sx(i)
            Sy(i+1) = Sy(i+1) + Sa*Sx(i+1)
            Sy(i+2) = Sy(i+2) + Sa*Sx(i+2)
            Sy(i+3) = Sy(i+3) + Sa*Sx(i+3)
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
            Sy(iy) = Sy(iy) + Sa*Sx(ix)
            ix = ix + Incx
            iy = iy + Incy
         ENDDO
      ENDIF
      END SUBROUTINE SAXPY
