!*==lsamen.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b LSAMEN
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download LSAMEN + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/lsamen.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/lsamen.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/lsamen.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       LOGICAL          FUNCTION LSAMEN( N, CA, CB )
!
!       .. Scalar Arguments ..
!       CHARACTER*( * )    CA, CB
!       INTEGER            N
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> LSAMEN  tests if the first N letters of CA are the same as the
!> first N letters of CB, regardless of case.
!> LSAMEN returns .TRUE. if CA and CB are equivalent except for case
!> and .FALSE. otherwise.  LSAMEN also returns .FALSE. if LEN( CA )
!> or LEN( CB ) is less than N.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of characters in CA and CB to be compared.
!> \endverbatim
!>
!> \param[in] CA
!> \verbatim
!>          CA is CHARACTER*(*)
!> \endverbatim
!>
!> \param[in] CB
!> \verbatim
!>          CB is CHARACTER*(*)
!>          CA and CB specify two character strings of length at least N.
!>          Only the first N characters of each string will be accessed.
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
!> \ingroup OTHERauxiliary
!
!  =====================================================================
      LOGICAL FUNCTION LSAMEN(N,Ca,Cb)
      IMPLICIT NONE
!*--LSAMEN78
!
!  -- LAPACK auxiliary routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER*(*) Ca , Cb
      INTEGER N
!     ..
!
! =====================================================================
!
!     .. Local Scalars ..
      INTEGER i
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC LEN
!     ..
!     .. Executable Statements ..
!
      LSAMEN = .FALSE.
      IF ( LEN(Ca)>=N .AND. LEN(Cb)>=N ) THEN
!
!     Do for each character in the two strings.
!
         DO i = 1 , N
!
!        Test if the characters are equal using LSAME.
!
            IF ( .NOT.LSAME(Ca(i:i),Cb(i:i)) ) GOTO 99999
!
         ENDDO
         LSAMEN = .TRUE.
      ENDIF
!
!
!     End of LSAMEN
!
99999 END FUNCTION LSAMEN
