!*==sdsdot.f90  processed by SPAG 7.51RB at 22:07 on  3 Mar 2022
!> \brief \b SDSDOT
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       REAL FUNCTION SDSDOT(N,SB,SX,INCX,SY,INCY)
!
!       .. Scalar Arguments ..
!       REAL SB
!       INTEGER INCX,INCY,N
!       ..
!       .. Array Arguments ..
!       REAL SX(*),SY(*)
!       ..
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>   Compute the inner product of two vectors with extended
!>   precision accumulation.
!>
!>   Returns S.P. result with dot product accumulated in D.P.
!>   SDSDOT = SB + sum for I = 0 to N-1 of SX(LX+I*INCX)*SY(LY+I*INCY),
!>   where LX = 1 if INCX .GE. 0, else LX = 1+(1-N)*INCX, and LY is
!>   defined in a similar way using INCY.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          number of elements in input vector(s)
!> \endverbatim
!>
!> \param[in] SB
!> \verbatim
!>          SB is REAL
!>          single precision scalar to be added to inner product
!> \endverbatim
!>
!> \param[in] SX
!> \verbatim
!>          SX is REAL array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
!>          single precision vector with N elements
!> \endverbatim
!>
!> \param[in] INCX
!> \verbatim
!>          INCX is INTEGER
!>          storage spacing between elements of SX
!> \endverbatim
!>
!> \param[in] SY
!> \verbatim
!>          SY is REAL array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
!>          single precision vector with N elements
!> \endverbatim
!>
!> \param[in] INCY
!> \verbatim
!>          INCY is INTEGER
!>          storage spacing between elements of SY
!> \endverbatim
!
!  Authors:
!  ========
!
!> \author Lawson, C. L., (JPL), Hanson, R. J., (SNLA),
!> \author Kincaid, D. R., (U. of Texas), Krogh, F. T., (JPL)
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
!>    REFERENCES
!>
!>    C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
!>    Krogh, Basic linear algebra subprograms for Fortran
!>    usage, Algorithm No. 539, Transactions on Mathematical
!>    Software 5, 3 (September 1979), pp. 308-323.
!>
!>    REVISION HISTORY  (YYMMDD)
!>
!>    791001  DATE WRITTEN
!>    890531  Changed all specific intrinsics to generic.  (WRB)
!>    890831  Modified array declarations.  (WRB)
!>    890831  REVISION DATE from Version 3.2
!>    891214  Prologue converted to Version 4.0 format.  (BAB)
!>    920310  Corrected definition of LX in DESCRIPTION.  (WRB)
!>    920501  Reformatted the REFERENCES section.  (WRB)
!>    070118  Reformat to LAPACK coding style
!> \endverbatim
!>
!  =====================================================================
      FUNCTION SDSDOT(N,Sb,Sx,Incx,Sy,Incy)
      USE F77KINDS
      USE F77KINDS                        
      IMPLICIT NONE
!*--SDSDOT119
      REAL :: SDSDOT
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER , INTENT(IN) :: N
      REAL , INTENT(IN) :: Sb
      REAL , INTENT(IN) , DIMENSION(*) :: Sx
      INTEGER , INTENT(IN) :: Incx
      REAL , INTENT(IN) , DIMENSION(*) :: Sy
      INTEGER , INTENT(IN) :: Incy
!
! Local variable declarations rewritten by SPAG
!
      REAL(r8kind) :: dsdot
      INTEGER :: i , kx , ky , ns
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
!     .. Local Scalars ..
!     ..
!     .. Intrinsic Functions ..
!     ..
      dsdot = Sb
      IF ( N<=0 ) THEN
         SDSDOT = dsdot
         RETURN
      ENDIF
      IF ( Incx==Incy .AND. Incx>0 ) THEN
!
!     Code for equal and positive increments.
!
         ns = N*Incx
         DO i = 1 , ns , Incx
            dsdot = dsdot + DBLE(Sx(i))*DBLE(Sy(i))
         ENDDO
      ELSE
!
!     Code for unequal or nonpositive increments.
!
         kx = 1
         ky = 1
         IF ( Incx<0 ) kx = 1 + (1-N)*Incx
         IF ( Incy<0 ) ky = 1 + (1-N)*Incy
         DO i = 1 , N
            dsdot = dsdot + DBLE(Sx(kx))*DBLE(Sy(ky))
            kx = kx + Incx
            ky = ky + Incy
         ENDDO
      ENDIF
      SDSDOT = dsdot
      END FUNCTION SDSDOT
