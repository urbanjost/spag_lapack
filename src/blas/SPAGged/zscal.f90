!*==zscal.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b ZSCAL
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZSCAL(N,ZA,ZX,INCX)
!
!       .. Scalar Arguments ..
!       COMPLEX*16 ZA
!       INTEGER INCX,N
!       ..
!       .. Array Arguments ..
!       COMPLEX*16 ZX(*)
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    ZSCAL scales a vector by a constant.
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
!> \param[in] ZA
!> \verbatim
!>          ZA is COMPLEX*16
!>           On entry, ZA specifies the scalar alpha.
!> \endverbatim
!>
!> \param[in,out] ZX
!> \verbatim
!>          ZX is COMPLEX*16 array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
!> \endverbatim
!>
!> \param[in] INCX
!> \verbatim
!>          INCX is INTEGER
!>         storage spacing between elements of ZX
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
!> \ingroup complex16_blas_level1
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>     jack dongarra, 3/11/78.
!>     modified 3/93 to return if incx .le. 0.
!>     modified 12/3/93, array(1) declarations changed to array(*)
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE ZSCAL(N,Za,Zx,Incx)
      USE F77KINDS                        
      IMPLICIT NONE
!*--ZSCAL83
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(IN) :: Za
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Zx
      INTEGER , INTENT(IN) :: Incx
!
! Local variable declarations rewritten by SPAG
!
      INTEGER :: i , nincx
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
      IF ( N<=0 .OR. Incx<=0 ) RETURN
      IF ( Incx==1 ) THEN
!
!        code for increment equal to 1
!
         DO i = 1 , N
            Zx(i) = Za*Zx(i)
         ENDDO
      ELSE
!
!        code for increment not equal to 1
!
         nincx = N*Incx
         DO i = 1 , nincx , Incx
            Zx(i) = Za*Zx(i)
         ENDDO
      ENDIF
      END SUBROUTINE ZSCAL
