!*==ilaprec.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b ILAPREC
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ILAPREC + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ilaprec.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ilaprec.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ilaprec.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       INTEGER FUNCTION ILAPREC( PREC )
!
!       .. Scalar Arguments ..
!       CHARACTER          PREC
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> This subroutine translated from a character string specifying an
!> intermediate precision to the relevant BLAST-specified integer
!> constant.
!>
!> ILAPREC returns an INTEGER.  If ILAPREC < 0, then the input is not a
!> character indicating a supported intermediate precision.  Otherwise
!> ILAPREC returns the constant value corresponding to PREC.
!> \endverbatim
!
!  Arguments:
!  ==========
!
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
!> \ingroup auxOTHERcomputational
!
!  =====================================================================
      FUNCTION ILAPREC(Prec)
      USE S_LSAME
      IMPLICIT NONE
!*--ILAPREC63
!
! PARAMETER definitions rewritten by SPAG
!
      INTEGER , PARAMETER  ::  BLAS_PREC_SINGLE = 211 ,                 &
     &                         BLAS_PREC_DOUBLE = 212 ,                 &
     &                         BLAS_PREC_INDIGENOUS = 213 ,             &
     &                         BLAS_PREC_EXTRA = 214
      INTEGER :: ILAPREC
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Prec
!
! End of declarations rewritten by SPAG
!
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
!     ..
!     .. External Functions ..
!     ..
!     .. Executable Statements ..
      IF ( LSAME(Prec,'S') ) THEN
         ILAPREC = BLAS_PREC_SINGLE
      ELSEIF ( LSAME(Prec,'D') ) THEN
         ILAPREC = BLAS_PREC_DOUBLE
      ELSEIF ( LSAME(Prec,'I') ) THEN
         ILAPREC = BLAS_PREC_INDIGENOUS
      ELSEIF ( LSAME(Prec,'X') .OR. LSAME(Prec,'E') ) THEN
         ILAPREC = BLAS_PREC_EXTRA
      ELSE
         ILAPREC = -1
      ENDIF
!
!     End of ILAPREC
!
      END FUNCTION ILAPREC
