!*==sscal.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b SSCAL
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE SSCAL(N,SA,SX,INCX)
!
!       .. Scalar Arguments ..
!       REAL SA
!       INTEGER INCX,N
!       ..
!       .. Array Arguments ..
!       REAL SX(*)
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    SSCAL scales a vector by a constant.
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
!> \param[in] SA
!> \verbatim
!>          SA is REAL
!>           On entry, SA specifies the scalar alpha.
!> \endverbatim
!>
!> \param[in,out] SX
!> \verbatim
!>          SX is REAL array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
!> \endverbatim
!>
!> \param[in] INCX
!> \verbatim
!>          INCX is INTEGER
!>         storage spacing between elements of SX
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
!>     modified 3/93 to return if incx .le. 0.
!>     modified 12/3/93, array(1) declarations changed to array(*)
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE SSCAL(N,Sa,Sx,Incx)
      IMPLICIT NONE
!*--SSCAL83
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(IN) :: Sa
      REAL , INTENT(INOUT) , DIMENSION(*) :: Sx
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
               Sx(i) = Sa*Sx(i)
            ENDDO
            IF ( N<5 ) RETURN
         ENDIF
         mp1 = m + 1
         DO i = mp1 , N , 5
            Sx(i) = Sa*Sx(i)
            Sx(i+1) = Sa*Sx(i+1)
            Sx(i+2) = Sa*Sx(i+2)
            Sx(i+3) = Sa*Sx(i+3)
            Sx(i+4) = Sa*Sx(i+4)
         ENDDO
      ELSE
!
!        code for increment not equal to 1
!
         nincx = N*Incx
         DO i = 1 , nincx , Incx
            Sx(i) = Sa*Sx(i)
         ENDDO
      ENDIF
      END SUBROUTINE SSCAL
