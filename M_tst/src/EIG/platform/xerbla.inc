!*==xerbla.f90  processed by SPAG 7.51RB at 17:37 on  4 Mar 2022
!> \brief \b xerbla
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE xerbla( SRNAME, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER*(*)      SRNAME
!       INTEGER            INFO
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> This is a special version of xerbla to be used only as part of
!> the test program for testing error exits from the LAPACK routines.
!> Error messages are printed if INFO.NE.INFOT or if SRNAME.NE.SRNAMT,
!> where INFOT and SRNAMT are values stored in COMMON.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] SRNAME
!> \verbatim
!>          SRNAME is CHARACTER*(*)
!>          The name of the subroutine calling xerbla.  This name should
!>          match the COMMON variable SRNAMT.
!> \endverbatim
!>
!> \param[in] INFO
!> \verbatim
!>          INFO is INTEGER
!>          The error return code from the calling subroutine.  INFO
!>          should equal the COMMON variable INFOT.
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
!> \ingroup aux_eig
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  The following variables are passed via the common blocks INFOC and
!>  SRNAMC:
!>
!>  INFOT   INTEGER      Expected integer return code
!>  NOUT    INTEGER      Unit number for printing error messages
!>  OK      LOGICAL      Set to .TRUE. if INFO = INFOT and
!>                       SRNAME = SRNAMT, otherwise set to .FALSE.
!>  LERR    LOGICAL      Set to .TRUE., indicating that xerbla was called
!>  SRNAMT  CHARACTER*(*) Expected name of calling subroutine
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE xerbla(Srname,Info)
      IMPLICIT NONE
!*--XERBLA79
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER*(*) Srname
      INTEGER Info
!     ..
!
!  =====================================================================
!
!     .. Scalars in Common ..
      LOGICAL LERr , OK
      CHARACTER*32 SRNamt
      INTEGER INFot , NOUt
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC LEN_TRIM
!     ..
!     .. Common blocks ..
      COMMON /INFOC / INFot , NOUt , OK , LERr
      COMMON /SRNAMC/ SRNamt
!     ..
!     .. Executable Statements ..
!
      LERr = .TRUE.
      IF ( Info/=INFot ) THEN
         IF ( INFot/=0 ) THEN
            WRITE (NOUt,FMT=99001) SRNamt(1:LEN_TRIM(SRNamt)) , Info ,  &
     &                             INFot
         ELSE
            WRITE (NOUt,FMT=99003) Srname(1:LEN_TRIM(Srname)) , Info
         ENDIF
         OK = .FALSE.
      ENDIF
      IF ( Srname/=SRNamt ) THEN
         WRITE (NOUt,FMT=99002) Srname(1:LEN_TRIM(Srname)) ,            &
     &                          SRNamt(1:LEN_TRIM(SRNamt))
         OK = .FALSE.
      ENDIF
      RETURN
!
99001 FORMAT (' *** xerbla was called from ',A,' with INFO = ',I6,      &
     &        ' instead of ',I2,' ***')
99002 FORMAT (' *** xerbla was called with SRNAME = ',A,' instead of ', &
     &        A6,' ***')
99003 FORMAT (' *** On entry to ',A,' parameter number ',I6,            &
     &        ' had an illegal value ***')
!
!     End of xerbla
!
      END SUBROUTINE xerbla
