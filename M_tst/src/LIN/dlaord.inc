!*==dlaord.f90  processed by SPAG 7.51RB at 17:34 on  4 Mar 2022
 
!> \brief \b dlaord
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE dlaord( JOB, N, X, INCX )
!
!       .. Scalar Arguments ..
!       CHARACTER          JOB
!       INTEGER            INCX, N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   X( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> dlaord sorts the elements of a vector x in increasing or decreasing
!> order.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] JOB
!> \verbatim
!>          JOB is CHARACTER
!>          = 'I':  Sort in increasing order
!>          = 'D':  Sort in decreasing order
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The length of the vector X.
!> \endverbatim
!>
!> \param[in,out] X
!> \verbatim
!>          X is DOUBLE PRECISION array, dimension
!>                         (1+(N-1)*INCX)
!>          On entry, the vector of length n to be sorted.
!>          On exit, the vector x is sorted in the prescribed order.
!> \endverbatim
!>
!> \param[in] INCX
!> \verbatim
!>          INCX is INTEGER
!>          The spacing between successive elements of X.  INCX >= 0.
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
!> \ingroup double_lin
!
!  =====================================================================
      SUBROUTINE dlaord(Job,N,X,Incx)
      IMPLICIT NONE
!*--DLAORD78
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Job
      INTEGER Incx , N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION X(*)
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
      INTEGER i , inc , ix , ixnext
      DOUBLE PRECISION temp
!     ..
!     .. External Functions ..
      LOGICAL LSAME
EXTERNAL    LSAME
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS
!     ..
!     .. Executable Statements ..
!
      inc = ABS(Incx)
      IF ( LSAME(Job,'I') ) THEN
!
!        Sort in increasing order
!
         DO i = 2 , N
            ix = 1 + (i-1)*inc
            DO WHILE ( ix/=1 )
               ixnext = ix - inc
               IF ( X(ix)>X(ixnext) ) EXIT
               temp = X(ix)
               X(ix) = X(ixnext)
               X(ixnext) = temp
               ix = ixnext
            ENDDO
         ENDDO
!
      ELSEIF ( LSAME(Job,'D') ) THEN
!
!        Sort in decreasing order
!
         DO i = 2 , N
            ix = 1 + (i-1)*inc
            DO WHILE ( ix/=1 )
               ixnext = ix - inc
               IF ( X(ix)<X(ixnext) ) EXIT
               temp = X(ix)
               X(ix) = X(ixnext)
               X(ixnext) = temp
               ix = ixnext
            ENDDO
         ENDDO
      ENDIF
!
!     End of dlaord
!
      END SUBROUTINE dlaord
