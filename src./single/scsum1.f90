!*==scsum1.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b SCSUM1 forms the 1-norm of the complex vector using the true absolute value.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SCSUM1 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/scsum1.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/scsum1.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/scsum1.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       REAL             FUNCTION SCSUM1( N, CX, INCX )
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
!> SCSUM1 takes the sum of the absolute values of a complex
!> vector and returns a single precision result.
!>
!> Based on SCASUM from the Level 1 BLAS.
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
!>          The vector whose elements will be summed.
!> \endverbatim
!>
!> \param[in] INCX
!> \verbatim
!>          INCX is INTEGER
!>          The spacing between successive values of CX.  INCX > 0.
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
!> \ingroup complexOTHERauxiliary
!
!> \par Contributors:
!  ==================
!>
!> Nick Higham for use with CLACON.
!
!  =====================================================================
      FUNCTION SCSUM1(N,Cx,Incx)
      IMPLICIT NONE
!*--SCSUM185
      REAL :: SCSUM1
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER , INTENT(IN) :: N
      COMPLEX , INTENT(IN) , DIMENSION(*) :: Cx
      INTEGER , INTENT(IN) :: Incx
!
! Local variable declarations rewritten by SPAG
!
      INTEGER :: i , nincx
      REAL :: stemp
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
      SCSUM1 = 0.0E0
      stemp = 0.0E0
      IF ( N<=0 ) RETURN
      IF ( Incx==1 ) THEN
!
!     CODE FOR INCREMENT EQUAL TO 1
!
         DO i = 1 , N
!
!        NEXT LINE MODIFIED.
!
            stemp = stemp + ABS(Cx(i))
         ENDDO
         SCSUM1 = stemp
         GOTO 99999
      ENDIF
!
!     CODE FOR INCREMENT NOT EQUAL TO 1
!
      nincx = N*Incx
      DO i = 1 , nincx , Incx
!
!        NEXT LINE MODIFIED.
!
         stemp = stemp + ABS(Cx(i))
      ENDDO
      SCSUM1 = stemp
!
!     End of SCSUM1
!
99999 END FUNCTION SCSUM1
