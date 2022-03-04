!*==dscal.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b DSCAL
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DSCAL(N,DA,DX,INCX)
!
!       .. Scalar Arguments ..
!       DOUBLE PRECISION DA
!       INTEGER INCX,N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION DX(*)
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    DSCAL scales a vector by a constant.
!>    uses unrolled loops for increment equal to 1.
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
!> \param[in] DA
!> \verbatim
!>          DA is DOUBLE PRECISION
!>           On entry, DA specifies the scalar alpha.
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
!>     modified 3/93 to return if incx .le. 0.
!>     modified 12/3/93, array(1) declarations changed to array(*)
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE DSCAL(N,Da,Dx,Incx)
      USE F77KINDS                        
      IMPLICIT NONE
!*--DSCAL84
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(IN) :: Da
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Dx
      INTEGER , INTENT(IN) :: Incx
!
! Local variable declarations rewritten by SPAG
!
      INTEGER :: i , m , mp1 , nincx
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
      IF ( N<=0 .OR. Incx<=0 ) RETURN
      IF ( Incx==1 ) THEN
!
!        code for increment equal to 1
!
!
!        clean-up loop
!
         m = MOD(N,5)
         IF ( m/=0 ) THEN
            DO i = 1 , m
               Dx(i) = Da*Dx(i)
            ENDDO
            IF ( N<5 ) RETURN
         ENDIF
         mp1 = m + 1
         DO i = mp1 , N , 5
            Dx(i) = Da*Dx(i)
            Dx(i+1) = Da*Dx(i+1)
            Dx(i+2) = Da*Dx(i+2)
            Dx(i+3) = Da*Dx(i+3)
            Dx(i+4) = Da*Dx(i+4)
         ENDDO
      ELSE
!
!        code for increment not equal to 1
!
         nincx = N*Incx
         DO i = 1 , nincx , Incx
            Dx(i) = Da*Dx(i)
         ENDDO
      ENDIF
      END SUBROUTINE DSCAL
