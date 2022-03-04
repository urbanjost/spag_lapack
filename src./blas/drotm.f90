!*==drotm.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
!> \brief \b DROTM
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DROTM(N,DX,INCX,DY,INCY,DPARAM)
!
!       .. Scalar Arguments ..
!       INTEGER INCX,INCY,N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION DPARAM(5),DX(*),DY(*)
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>    APPLY THE MODIFIED GIVENS TRANSFORMATION, H, TO THE 2 BY N MATRIX
!>
!>    (DX**T) , WHERE **T INDICATES TRANSPOSE. THE ELEMENTS OF DX ARE IN
!>    (DY**T)
!>
!>    DX(LX+I*INCX), I = 0 TO N-1, WHERE LX = 1 IF INCX .GE. 0, ELSE
!>    LX = (-INCX)*N, AND SIMILARLY FOR SY USING LY AND INCY.
!>    WITH DPARAM(1)=DFLAG, H HAS ONE OF THE FOLLOWING FORMS..
!>
!>    DFLAG=-1.D0     DFLAG=0.D0        DFLAG=1.D0     DFLAG=-2.D0
!>
!>      (DH11  DH12)    (1.D0  DH12)    (DH11  1.D0)    (1.D0  0.D0)
!>    H=(          )    (          )    (          )    (          )
!>      (DH21  DH22),   (DH21  1.D0),   (-1.D0 DH22),   (0.D0  1.D0).
!>    SEE DROTMG FOR A DESCRIPTION OF DATA STORAGE IN DPARAM.
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
!>
!> \param[in] DPARAM
!> \verbatim
!>          DPARAM is DOUBLE PRECISION array, dimension (5)
!>     DPARAM(1)=DFLAG
!>     DPARAM(2)=DH11
!>     DPARAM(3)=DH21
!>     DPARAM(4)=DH12
!>     DPARAM(5)=DH22
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
!  =====================================================================
      SUBROUTINE DROTM(N,Dx,Incx,Dy,Incy,Dparam)
      USE F77KINDS
      IMPLICIT NONE
!*--DROTM101
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER , INTENT(IN) :: N
      REAL(r8kind) , INTENT(INOUT) , DIMENSION(*) :: Dx
      INTEGER , INTENT(IN) :: Incx
      REAL(r8kind) , INTENT(INOUT) , DIMENSION(*) :: Dy
      INTEGER , INTENT(IN) :: Incy
      REAL(r8kind) , INTENT(IN) , DIMENSION(5) :: Dparam
!
! Local variable declarations rewritten by SPAG
!
      REAL(r8kind) :: dflag , dh11 , dh12 , dh21 , dh22 , w , z
      INTEGER :: i , kx , ky , nsteps
      REAL(r8kind) , SAVE :: two , zero
!
! End of declarations rewritten by SPAG
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
!     .. Data statements ..
      DATA zero , two/0.D0 , 2.D0/
!     ..
!
      dflag = Dparam(1)
      IF ( N<=0 .OR. (dflag+two==zero) ) RETURN
      IF ( Incx==Incy .AND. Incx>0 ) THEN
!
         nsteps = N*Incx
         IF ( dflag<zero ) THEN
            dh11 = Dparam(2)
            dh12 = Dparam(4)
            dh21 = Dparam(3)
            dh22 = Dparam(5)
            DO i = 1 , nsteps , Incx
               w = Dx(i)
               z = Dy(i)
               Dx(i) = w*dh11 + z*dh12
               Dy(i) = w*dh21 + z*dh22
            ENDDO
         ELSEIF ( dflag==zero ) THEN
            dh12 = Dparam(4)
            dh21 = Dparam(3)
            DO i = 1 , nsteps , Incx
               w = Dx(i)
               z = Dy(i)
               Dx(i) = w + z*dh12
               Dy(i) = w*dh21 + z
            ENDDO
         ELSE
            dh11 = Dparam(2)
            dh22 = Dparam(5)
            DO i = 1 , nsteps , Incx
               w = Dx(i)
               z = Dy(i)
               Dx(i) = w*dh11 + z
               Dy(i) = -w + dh22*z
            ENDDO
         ENDIF
      ELSE
         kx = 1
         ky = 1
         IF ( Incx<0 ) kx = 1 + (1-N)*Incx
         IF ( Incy<0 ) ky = 1 + (1-N)*Incy
!
         IF ( dflag<zero ) THEN
            dh11 = Dparam(2)
            dh12 = Dparam(4)
            dh21 = Dparam(3)
            dh22 = Dparam(5)
            DO i = 1 , N
               w = Dx(kx)
               z = Dy(ky)
               Dx(kx) = w*dh11 + z*dh12
               Dy(ky) = w*dh21 + z*dh22
               kx = kx + Incx
               ky = ky + Incy
            ENDDO
         ELSEIF ( dflag==zero ) THEN
            dh12 = Dparam(4)
            dh21 = Dparam(3)
            DO i = 1 , N
               w = Dx(kx)
               z = Dy(ky)
               Dx(kx) = w + z*dh12
               Dy(ky) = w*dh21 + z
               kx = kx + Incx
               ky = ky + Incy
            ENDDO
         ELSE
            dh11 = Dparam(2)
            dh22 = Dparam(5)
            DO i = 1 , N
               w = Dx(kx)
               z = Dy(ky)
               Dx(kx) = w*dh11 + z
               Dy(ky) = -w + dh22*z
               kx = kx + Incx
               ky = ky + Incy
            ENDDO
         ENDIF
      ENDIF
      END SUBROUTINE DROTM
