!*==idamax.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b IDAMAX
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       INTEGER FUNCTION IDAMAX(N,DX,INCX)
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
!>    IDAMAX finds the index of the first element having maximum absolute value.
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
!> \ingroup aux_blas
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
      FUNCTION IDAMAX(N,Dx,Incx)
      USE F77KINDS                        
      IMPLICIT NONE
!*--IDAMAX76
      INTEGER :: IDAMAX
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: Dx
      INTEGER , INTENT(IN) :: Incx
!
! Local variable declarations rewritten by SPAG
!
      REAL(R8KIND) :: dmax
      INTEGER :: i , ix
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
      IDAMAX = 0
      IF ( N<1 .OR. Incx<=0 ) RETURN
      IDAMAX = 1
      IF ( N==1 ) RETURN
      IF ( Incx==1 ) THEN
!
!        code for increment equal to 1
!
         dmax = DABS(Dx(1))
         DO i = 2 , N
            IF ( DABS(Dx(i))>dmax ) THEN
               IDAMAX = i
               dmax = DABS(Dx(i))
            ENDIF
         ENDDO
      ELSE
!
!        code for increment not equal to 1
!
         ix = 1
         dmax = DABS(Dx(1))
         ix = ix + Incx
         DO i = 2 , N
            IF ( DABS(Dx(ix))>dmax ) THEN
               IDAMAX = i
               dmax = DABS(Dx(ix))
            ENDIF
            ix = ix + Incx
         ENDDO
      ENDIF
      END FUNCTION IDAMAX
