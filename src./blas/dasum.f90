!*==dasum.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
!> \brief \b DASUM
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       DOUBLE PRECISION FUNCTION DASUM(N,DX,INCX)
!
!       .. Scalar Arguments ..
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
!>    DASUM takes the sum of the absolute values.
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
      FUNCTION DASUM(N,Dx,Incx)
      USE F77KINDS
      USE F77KINDS                        
      IMPLICIT NONE
!*--DASUM77
      REAL(r8kind) :: DASUM
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER , INTENT(IN) :: N
      REAL(r8kind) , INTENT(IN) , DIMENSION(*) :: Dx
      INTEGER , INTENT(IN) :: Incx
!
! Local variable declarations rewritten by SPAG
!
      REAL(r8kind) :: dtemp
      INTEGER :: i , m , mp1 , nincx
!
! End of declarations rewritten by SPAG
!
!
! Dummy argument declarations rewritten by SPAG
!
!
! Local variable declarations rewritten by SPAG
!
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
      DASUM = 0.0D0
      dtemp = 0.0D0
      IF ( N<=0 .OR. Incx<=0 ) RETURN
      IF ( Incx==1 ) THEN
!        code for increment equal to 1
!
!
!        clean-up loop
!
         m = MOD(N,6)
         IF ( m/=0 ) THEN
            DO i = 1 , m
               dtemp = dtemp + DABS(Dx(i))
            ENDDO
            IF ( N<6 ) THEN
               DASUM = dtemp
               RETURN
            ENDIF
         ENDIF
         mp1 = m + 1
         DO i = mp1 , N , 6
            dtemp = dtemp + DABS(Dx(i)) + DABS(Dx(i+1)) + DABS(Dx(i+2)) &
     &              + DABS(Dx(i+3)) + DABS(Dx(i+4)) + DABS(Dx(i+5))
         ENDDO
      ELSE
!
!        code for increment not equal to 1
!
         nincx = N*Incx
         DO i = 1 , nincx , Incx
            dtemp = dtemp + DABS(Dx(i))
         ENDDO
      ENDIF
      DASUM = dtemp
      END FUNCTION DASUM
