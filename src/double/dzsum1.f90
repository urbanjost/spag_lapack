!*==dzsum1.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b DZSUM1 forms the 1-norm of the complex vector using the true absolute value.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DZSUM1 + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dzsum1.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dzsum1.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dzsum1.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       DOUBLE PRECISION FUNCTION DZSUM1( N, CX, INCX )
!
!       .. Scalar Arguments ..
!       INTEGER            INCX, N
!       ..
!       .. Array Arguments ..
!       COMPLEX*16         CX( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DZSUM1 takes the sum of the absolute values of a complex
!> vector and returns a double precision result.
!>
!> Based on DZASUM from the Level 1 BLAS.
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
!>          CX is COMPLEX*16 array, dimension (N)
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
!> \ingroup complex16OTHERauxiliary
!
!> \par Contributors:
!  ==================
!>
!> Nick Higham for use with ZLACON.
!
!  =====================================================================
      DOUBLE PRECISION FUNCTION DZSUM1(N,Cx,Incx)
      IMPLICIT NONE
!*--DZSUM185
!
!  -- LAPACK auxiliary routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Incx , N
!     ..
!     .. Array Arguments ..
      COMPLEX*16 Cx(*)
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
      INTEGER i , nincx
      DOUBLE PRECISION stemp
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS
!     ..
!     .. Executable Statements ..
!
      DZSUM1 = 0.0D0
      stemp = 0.0D0
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
         DZSUM1 = stemp
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
      DZSUM1 = stemp
      RETURN
!
!     End of DZSUM1
!
99999 END FUNCTION DZSUM1
