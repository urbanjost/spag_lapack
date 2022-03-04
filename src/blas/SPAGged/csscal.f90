!*==csscal.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b CSSCAL
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE CSSCAL(N,SA,CX,INCX)
!
!       .. Scalar Arguments ..
!       REAL SA
!       INTEGER INCX,N
!       ..
!       .. Array Arguments ..
!       COMPLEX CX(*)
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    CSSCAL scales a complex vector by a real constant.
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
!> \param[in,out] CX
!> \verbatim
!>          CX is COMPLEX array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
!> \endverbatim
!>
!> \param[in] INCX
!> \verbatim
!>          INCX is INTEGER
!>         storage spacing between elements of CX
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
!> \ingroup complex_blas_level1
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
      SUBROUTINE CSSCAL(N,Sa,Cx,Incx)
      IMPLICIT NONE
!*--CSSCAL82
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(IN) :: Sa
      COMPLEX , INTENT(INOUT) , DIMENSION(*) :: Cx
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
!     .. Intrinsic Functions ..
!     ..
      IF ( N<=0 .OR. Incx<=0 ) RETURN
      IF ( Incx==1 ) THEN
!
!        code for increment equal to 1
!
         DO i = 1 , N
            Cx(i) = CMPLX(Sa*REAL(Cx(i)),Sa*AIMAG(Cx(i)))
         ENDDO
      ELSE
!
!        code for increment not equal to 1
!
         nincx = N*Incx
         DO i = 1 , nincx , Incx
            Cx(i) = CMPLX(Sa*REAL(Cx(i)),Sa*AIMAG(Cx(i)))
         ENDDO
      ENDIF
      END SUBROUTINE CSSCAL
