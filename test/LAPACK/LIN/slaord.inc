!*==slaord.f90  processed by SPAG 7.51RB at 17:34 on  4 Mar 2022
 
!> \brief \b slaord
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE slaord( JOB, N, X, INCX )
!
!       .. Scalar Arguments ..
!       CHARACTER          JOB
!       INTEGER            INCX, N
!       ..
!       .. Array Arguments ..
!       REAL               X( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> slaord sorts the elements of a vector x in increasing or decreasing
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
!>          X is REAL array, dimension
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
!> \ingroup single_lin
!
!  =====================================================================
      SUBROUTINE slaord(Job,N,X,Incx)
      IMPLICIT NONE
!*--SLAORD78
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
      REAL X(*)
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
      INTEGER i , inc , ix , ixnext
      REAL temp
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
!     End of slaord
!
      END SUBROUTINE slaord
