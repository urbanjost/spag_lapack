!*==icmax1.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b ICMAX1 finds the index of the first vector element of maximum absolute value.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ICMAX1 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/icmax1.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/icmax1.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/icmax1.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       INTEGER          FUNCTION ICMAX1( N, CX, INCX )
!
!       .. Scalar Arguments ..
!       INTEGER            INCX, N
!       ..
!       .. Array Arguments ..
!       COMPLEX            CX( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ICMAX1 finds the index of the first vector element of maximum absolute value.
!>
!> Based on ICAMAX from Level 1 BLAS.
!> The change is to use the 'genuine' absolute value.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of elements in the vector CX.
!> \endverbatim
!>
!> \param[in] CX
!> \verbatim
!>          CX is COMPLEX array, dimension (N)
!>          The vector CX. The ICMAX1 function returns the index of its first
!>          element of maximum absolute value.
!> \endverbatim
!>
!> \param[in] INCX
!> \verbatim
!>          INCX is INTEGER
!>          The spacing between successive values of CX.  INCX >= 1.
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
!> \date February 2014
!
!> \ingroup complexOTHERauxiliary
!
!> \par Contributors:
!  ==================
!>
!> Nick Higham for use with CLACON.
!
!  =====================================================================
      FUNCTION ICMAX1(N,Cx,Incx)
      IMPLICIT NONE
!*--ICMAX185
      INTEGER :: ICMAX1
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(IN) , DIMENSION(*) :: Cx
      INTEGER , INTENT(IN) :: Incx
!
! Local variable declarations rewritten by SPAG
!
      INTEGER :: i , ix
      REAL :: smax
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
!     .. Executable Statements ..
!
      ICMAX1 = 0
      IF ( N<1 .OR. Incx<=0 ) RETURN
      ICMAX1 = 1
      IF ( N==1 ) RETURN
      IF ( Incx==1 ) THEN
!
!        code for increment equal to 1
!
         smax = ABS(Cx(1))
         DO i = 2 , N
            IF ( ABS(Cx(i))>smax ) THEN
               ICMAX1 = i
               smax = ABS(Cx(i))
            ENDIF
         ENDDO
      ELSE
!
!        code for increment not equal to 1
!
         ix = 1
         smax = ABS(Cx(1))
         ix = ix + Incx
         DO i = 2 , N
            IF ( ABS(Cx(ix))>smax ) THEN
               ICMAX1 = i
               smax = ABS(Cx(ix))
            ENDIF
            ix = ix + Incx
         ENDDO
      ENDIF
!
!     End of ICMAX1
!
      END FUNCTION ICMAX1
