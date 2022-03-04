!*==icamax.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b ICAMAX
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       INTEGER FUNCTION ICAMAX(N,CX,INCX)
!
!       .. Scalar Arguments ..
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
!>    ICAMAX finds the index of the first element having maximum |Re(.)| + |Im(.)|
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
!> \param[in] CX
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
      INTEGER FUNCTION ICAMAX(N,Cx,Incx)
      IMPLICIT NONE
!*--ICAMAX75
!
!  -- Reference BLAS level1 routine (version 3.8.0) --
!  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2017
!
!     .. Scalar Arguments ..
      INTEGER Incx , N
!     ..
!     .. Array Arguments ..
      COMPLEX Cx(*)
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
      REAL smax
      INTEGER i , ix
!     ..
!     .. External Functions ..
      REAL SCABS1
      EXTERNAL SCABS1
!     ..
      ICAMAX = 0
      IF ( N<1 .OR. Incx<=0 ) RETURN
      ICAMAX = 1
      IF ( N==1 ) RETURN
      IF ( Incx==1 ) THEN
!
!        code for increment equal to 1
!
         smax = SCABS1(Cx(1))
         DO i = 2 , N
            IF ( SCABS1(Cx(i))>smax ) THEN
               ICAMAX = i
               smax = SCABS1(Cx(i))
            ENDIF
         ENDDO
      ELSE
!
!        code for increment not equal to 1
!
         ix = 1
         smax = SCABS1(Cx(1))
         ix = ix + Incx
         DO i = 2 , N
            IF ( SCABS1(Cx(ix))>smax ) THEN
               ICAMAX = i
               smax = SCABS1(Cx(ix))
            ENDIF
            ix = ix + Incx
         ENDDO
      ENDIF
      END FUNCTION ICAMAX
