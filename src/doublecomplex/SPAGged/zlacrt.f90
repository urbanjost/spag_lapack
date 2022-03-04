!*==zlacrt.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b ZLACRT performs a linear transformation of a pair of complex vectors.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZLACRT + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlacrt.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlacrt.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlacrt.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZLACRT( N, CX, INCX, CY, INCY, C, S )
!
!       .. Scalar Arguments ..
!       INTEGER            INCX, INCY, N
!       COMPLEX*16         C, S
!       ..
!       .. Array Arguments ..
!       COMPLEX*16         CX( * ), CY( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZLACRT performs the operation
!>
!>    (  c  s )( x )  ==> ( x )
!>    ( -s  c )( y )      ( y )
!>
!> where c and s are complex and the vectors x and y are complex.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of elements in the vectors CX and CY.
!> \endverbatim
!>
!> \param[in,out] CX
!> \verbatim
!>          CX is COMPLEX*16 array, dimension (N)
!>          On input, the vector x.
!>          On output, CX is overwritten with c*x + s*y.
!> \endverbatim
!>
!> \param[in] INCX
!> \verbatim
!>          INCX is INTEGER
!>          The increment between successive values of CX.  INCX <> 0.
!> \endverbatim
!>
!> \param[in,out] CY
!> \verbatim
!>          CY is COMPLEX*16 array, dimension (N)
!>          On input, the vector y.
!>          On output, CY is overwritten with -s*x + c*y.
!> \endverbatim
!>
!> \param[in] INCY
!> \verbatim
!>          INCY is INTEGER
!>          The increment between successive values of CY.  INCY <> 0.
!> \endverbatim
!>
!> \param[in] C
!> \verbatim
!>          C is COMPLEX*16
!> \endverbatim
!>
!> \param[in] S
!> \verbatim
!>          S is COMPLEX*16
!>          C and S define the matrix
!>             [  C   S  ].
!>             [ -S   C  ]
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
!> \date December 2016
!
!> \ingroup complex16OTHERauxiliary
!
!  =====================================================================
      SUBROUTINE ZLACRT(N,Cx,Incx,Cy,Incy,C,S)
      USE F77KINDS                        
      IMPLICIT NONE
!*--ZLACRT110
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Cx
      INTEGER , INTENT(IN) :: Incx
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Cy
      INTEGER , INTENT(IN) :: Incy
      COMPLEX(CX16KIND) , INTENT(IN) :: C
      COMPLEX(CX16KIND) , INTENT(IN) :: S
!
! Local variable declarations rewritten by SPAG
!
      COMPLEX(CX16KIND) :: ctemp
      INTEGER :: i , ix , iy
!
! End of declarations rewritten by SPAG
!
!     ..
!     .. Array Arguments ..
!     ..
!
! =====================================================================
!
!     .. Local Scalars ..
!     ..
!     .. Executable Statements ..
!
      IF ( N<=0 ) RETURN
      IF ( Incx==1 .AND. Incy==1 ) THEN
!
!     Code for both increments equal to 1
!
         DO i = 1 , N
            ctemp = C*Cx(i) + S*Cy(i)
            Cy(i) = C*Cy(i) - S*Cx(i)
            Cx(i) = ctemp
         ENDDO
         GOTO 99999
      ENDIF
!
!     Code for unequal increments or equal increments not equal to 1
!
      ix = 1
      iy = 1
      IF ( Incx<0 ) ix = (-N+1)*Incx + 1
      IF ( Incy<0 ) iy = (-N+1)*Incy + 1
      DO i = 1 , N
         ctemp = C*Cx(ix) + S*Cy(iy)
         Cy(iy) = C*Cy(iy) - S*Cx(ix)
         Cx(ix) = ctemp
         ix = ix + Incx
         iy = iy + Incy
      ENDDO
99999 END SUBROUTINE ZLACRT
