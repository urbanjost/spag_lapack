!*==zlacgv.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b ZLACGV conjugates a complex vector.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZLACGV + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlacgv.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlacgv.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlacgv.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZLACGV( N, X, INCX )
!
!       .. Scalar Arguments ..
!       INTEGER            INCX, N
!       ..
!       .. Array Arguments ..
!       COMPLEX*16         X( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZLACGV conjugates a complex vector of length N.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The length of the vector X.  N >= 0.
!> \endverbatim
!>
!> \param[in,out] X
!> \verbatim
!>          X is COMPLEX*16 array, dimension
!>                         (1+(N-1)*abs(INCX))
!>          On entry, the vector of length N to be conjugated.
!>          On exit, X is overwritten with conjg(X).
!> \endverbatim
!>
!> \param[in] INCX
!> \verbatim
!>          INCX is INTEGER
!>          The spacing between successive elements of X.
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
      SUBROUTINE ZLACGV(N,X,Incx)
      USE F77KINDS                        
      IMPLICIT NONE
!*--ZLACGV79
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER , INTENT(IN) :: N
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: X
      INTEGER , INTENT(IN) :: Incx
!
! Local variable declarations rewritten by SPAG
!
      INTEGER :: i , ioff
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
!     .. Intrinsic Functions ..
!     ..
!     .. Executable Statements ..
!
      IF ( Incx==1 ) THEN
         DO i = 1 , N
            X(i) = DCONJG(X(i))
         ENDDO
      ELSE
         ioff = 1
         IF ( Incx<0 ) ioff = 1 - (N-1)*Incx
         DO i = 1 , N
            X(ioff) = DCONJG(X(ioff))
            ioff = ioff + Incx
         ENDDO
      ENDIF
!
!     End of ZLACGV
!
      END SUBROUTINE ZLACGV
