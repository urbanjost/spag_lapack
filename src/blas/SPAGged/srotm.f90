!*==srotm.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b SROTM
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE SROTM(N,SX,INCX,SY,INCY,SPARAM)
!
!       .. Scalar Arguments ..
!       INTEGER INCX,INCY,N
!       ..
!       .. Array Arguments ..
!       REAL SPARAM(5),SX(*),SY(*)
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
!>    (SX**T) , WHERE **T INDICATES TRANSPOSE. THE ELEMENTS OF SX ARE IN
!>    (SX**T)
!>
!>    SX(LX+I*INCX), I = 0 TO N-1, WHERE LX = 1 IF INCX .GE. 0, ELSE
!>    LX = (-INCX)*N, AND SIMILARLY FOR SY USING USING LY AND INCY.
!>    WITH SPARAM(1)=SFLAG, H HAS ONE OF THE FOLLOWING FORMS..
!>
!>    SFLAG=-1.E0     SFLAG=0.E0        SFLAG=1.E0     SFLAG=-2.E0
!>
!>      (SH11  SH12)    (1.E0  SH12)    (SH11  1.E0)    (1.E0  0.E0)
!>    H=(          )    (          )    (          )    (          )
!>      (SH21  SH22),   (SH21  1.E0),   (-1.E0 SH22),   (0.E0  1.E0).
!>    SEE  SROTMG FOR A DESCRIPTION OF DATA STORAGE IN SPARAM.
!>
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
!>
!> \param[in] SPARAM
!> \verbatim
!>          SPARAM is REAL array, dimension (5)
!>     SPARAM(1)=SFLAG
!>     SPARAM(2)=SH11
!>     SPARAM(3)=SH21
!>     SPARAM(4)=SH12
!>     SPARAM(5)=SH22
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
!  =====================================================================
      SUBROUTINE SROTM(N,Sx,Incx,Sy,Incy,Sparam)
      IMPLICIT NONE
!*--SROTM101
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(INOUT) , DIMENSION(*) :: Sx
      INTEGER , INTENT(IN) :: Incx
      REAL , INTENT(INOUT) , DIMENSION(*) :: Sy
      INTEGER , INTENT(IN) :: Incy
      REAL , INTENT(IN) , DIMENSION(5) :: Sparam
!
! Local variable declarations rewritten by SPAG
!
      INTEGER :: i , kx , ky , nsteps
      REAL :: sflag , sh11 , sh12 , sh21 , sh22 , w , z
      REAL , SAVE :: two , zero
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
      DATA zero , two/0.E0 , 2.E0/
!     ..
!
      sflag = Sparam(1)
      IF ( N<=0 .OR. (sflag+two==zero) ) RETURN
      IF ( Incx==Incy .AND. Incx>0 ) THEN
!
         nsteps = N*Incx
         IF ( sflag<zero ) THEN
            sh11 = Sparam(2)
            sh12 = Sparam(4)
            sh21 = Sparam(3)
            sh22 = Sparam(5)
            DO i = 1 , nsteps , Incx
               w = Sx(i)
               z = Sy(i)
               Sx(i) = w*sh11 + z*sh12
               Sy(i) = w*sh21 + z*sh22
            ENDDO
         ELSEIF ( sflag==zero ) THEN
            sh12 = Sparam(4)
            sh21 = Sparam(3)
            DO i = 1 , nsteps , Incx
               w = Sx(i)
               z = Sy(i)
               Sx(i) = w + z*sh12
               Sy(i) = w*sh21 + z
            ENDDO
         ELSE
            sh11 = Sparam(2)
            sh22 = Sparam(5)
            DO i = 1 , nsteps , Incx
               w = Sx(i)
               z = Sy(i)
               Sx(i) = w*sh11 + z
               Sy(i) = -w + sh22*z
            ENDDO
         ENDIF
      ELSE
         kx = 1
         ky = 1
         IF ( Incx<0 ) kx = 1 + (1-N)*Incx
         IF ( Incy<0 ) ky = 1 + (1-N)*Incy
!
         IF ( sflag<zero ) THEN
            sh11 = Sparam(2)
            sh12 = Sparam(4)
            sh21 = Sparam(3)
            sh22 = Sparam(5)
            DO i = 1 , N
               w = Sx(kx)
               z = Sy(ky)
               Sx(kx) = w*sh11 + z*sh12
               Sy(ky) = w*sh21 + z*sh22
               kx = kx + Incx
               ky = ky + Incy
            ENDDO
         ELSEIF ( sflag==zero ) THEN
            sh12 = Sparam(4)
            sh21 = Sparam(3)
            DO i = 1 , N
               w = Sx(kx)
               z = Sy(ky)
               Sx(kx) = w + z*sh12
               Sy(ky) = w*sh21 + z
               kx = kx + Incx
               ky = ky + Incy
            ENDDO
         ELSE
            sh11 = Sparam(2)
            sh22 = Sparam(5)
            DO i = 1 , N
               w = Sx(kx)
               z = Sy(ky)
               Sx(kx) = w*sh11 + z
               Sy(ky) = -w + sh22*z
               kx = kx + Incx
               ky = ky + Incy
            ENDDO
         ENDIF
      ENDIF
      END SUBROUTINE SROTM
