!*==chla_transtype.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b CHLA_TRANSTYPE
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CHLA_TRANSTYPE + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/chla_transtype.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/chla_transtype.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/chla_transtype.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       CHARACTER*1 FUNCTION CHLA_TRANSTYPE( TRANS )
!
!       .. Scalar Arguments ..
!       INTEGER            TRANS
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> This subroutine translates from a BLAST-specified integer constant to
!> the character string specifying a transposition operation.
!>
!> CHLA_TRANSTYPE returns an CHARACTER*1.  If CHLA_TRANSTYPE is 'X',
!> then input is not an integer indicating a transposition operator.
!> Otherwise CHLA_TRANSTYPE returns the constant value corresponding to
!> TRANS.
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
      CHARACTER*1 FUNCTION CHLA_TRANSTYPE(Trans)
      IMPLICIT NONE
!*--CHLA_TRANSTYPE62
!
!  -- LAPACK computational routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Trans
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      INTEGER BLAS_NO_TRANS , BLAS_TRANS , BLAS_CONJ_TRANS
      PARAMETER (BLAS_NO_TRANS=111,BLAS_TRANS=112,BLAS_CONJ_TRANS=113)
!     ..
!     .. Executable Statements ..
      IF ( Trans==BLAS_NO_TRANS ) THEN
         CHLA_TRANSTYPE = 'N'
      ELSEIF ( Trans==BLAS_TRANS ) THEN
         CHLA_TRANSTYPE = 'T'
      ELSEIF ( Trans==BLAS_CONJ_TRANS ) THEN
         CHLA_TRANSTYPE = 'C'
      ELSE
         CHLA_TRANSTYPE = 'X'
      ENDIF
!
!     End of CHLA_TRANSTYPE
!
      END FUNCTION CHLA_TRANSTYPE
