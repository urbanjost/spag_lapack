!*==derrrfp.f90  processed by SPAG 7.51RB at 17:34 on  4 Mar 2022
!> \brief \b derrrfp
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE derrrfp( NUNIT )
!
!       .. Scalar Arguments ..
!       INTEGER            NUNIT
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> derrrfp tests the error exits for the DOUBLE PRECISION driver routines
!> for solving linear systems of equations.
!>
!> ddrvrfp tests the DOUBLE PRECISION LAPACK RFP routines:
!>     DTFSM, DTFTRI, DSFRK, DTFTTP, DTFTTR, DPFTRF, DPFTRS, DTPTTF,
!>     DTPTTR, DTRTTF, and DTRTTP
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] NUNIT
!> \verbatim
!>          NUNIT is INTEGER
!>          The unit number for output.
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
      SUBROUTINE derrrfp(Nunit)
      IMPLICIT NONE
!*--DERRRFP56
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Nunit
!     ..
!
!  =====================================================================
!
!     ..
!     .. Local Scalars ..
      INTEGER info
      DOUBLE PRECISION alpha , beta
!     ..
!     .. Local Arrays ..
      DOUBLE PRECISION a(1,1) , b(1,1)
!     ..
!     .. External Subroutines ..
EXTERNAL    DTFSM,   DTFTRI,   DSFRK,   DTFTTP,   DTFTTR
EXTERNAL    DPFTRI,   DPFTRF,   DPFTRS,   DTPTTF,   DTPTTR,   DTRTTF
EXTERNAL    DTRTTP
!     ..
!     .. Scalars in Common ..
      LOGICAL LERr , OK
      CHARACTER*32 SRNamt
      INTEGER INFot , NOUt
!     ..
!     .. Common blocks ..
      COMMON /INFOC / INFot , NOUt , OK , LERr
      COMMON /SRNAMC/ SRNamt
!     ..
!     .. Executable Statements ..
!
      NOUt = Nunit
      OK = .TRUE.
      a(1,1) = 1.0D+0
      b(1,1) = 1.0D+0
      alpha = 1.0D+0
      beta = 1.0D+0
!
      SRNamt = 'DPFTRF'
      INFot = 1
      CALL DPFTRF('/','U',0,a,info)
      CALL CHKXER('DPFTRF',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL DPFTRF('N','/',0,a,info)
      CALL CHKXER('DPFTRF',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL DPFTRF('N','U',-1,a,info)
      CALL CHKXER('DPFTRF',INFot,NOUt,LERr,OK)
!
      SRNamt = 'DPFTRS'
      INFot = 1
      CALL DPFTRS('/','U',0,0,a,b,1,info)
      CALL CHKXER('DPFTRS',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL DPFTRS('N','/',0,0,a,b,1,info)
      CALL CHKXER('DPFTRS',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL DPFTRS('N','U',-1,0,a,b,1,info)
      CALL CHKXER('DPFTRS',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL DPFTRS('N','U',0,-1,a,b,1,info)
      CALL CHKXER('DPFTRS',INFot,NOUt,LERr,OK)
      INFot = 7
      CALL DPFTRS('N','U',0,0,a,b,0,info)
      CALL CHKXER('DPFTRS',INFot,NOUt,LERr,OK)
!
      SRNamt = 'DPFTRI'
      INFot = 1
      CALL DPFTRI('/','U',0,a,info)
      CALL CHKXER('DPFTRI',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL DPFTRI('N','/',0,a,info)
      CALL CHKXER('DPFTRI',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL DPFTRI('N','U',-1,a,info)
      CALL CHKXER('DPFTRI',INFot,NOUt,LERr,OK)
!
      SRNamt = 'DTFSM '
      INFot = 1
      CALL DTFSM('/','L','U','T','U',0,0,alpha,a,b,1)
      CALL CHKXER('DTFSM ',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL DTFSM('N','/','U','T','U',0,0,alpha,a,b,1)
      CALL CHKXER('DTFSM ',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL DTFSM('N','L','/','T','U',0,0,alpha,a,b,1)
      CALL CHKXER('DTFSM ',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL DTFSM('N','L','U','/','U',0,0,alpha,a,b,1)
      CALL CHKXER('DTFSM ',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL DTFSM('N','L','U','T','/',0,0,alpha,a,b,1)
      CALL CHKXER('DTFSM ',INFot,NOUt,LERr,OK)
      INFot = 6
      CALL DTFSM('N','L','U','T','U',-1,0,alpha,a,b,1)
      CALL CHKXER('DTFSM ',INFot,NOUt,LERr,OK)
      INFot = 7
      CALL DTFSM('N','L','U','T','U',0,-1,alpha,a,b,1)
      CALL CHKXER('DTFSM ',INFot,NOUt,LERr,OK)
      INFot = 11
      CALL DTFSM('N','L','U','T','U',0,0,alpha,a,b,0)
      CALL CHKXER('DTFSM ',INFot,NOUt,LERr,OK)
!
      SRNamt = 'DTFTRI'
      INFot = 1
      CALL DTFTRI('/','L','N',0,a,info)
      CALL CHKXER('DTFTRI',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL DTFTRI('N','/','N',0,a,info)
      CALL CHKXER('DTFTRI',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL DTFTRI('N','L','/',0,a,info)
      CALL CHKXER('DTFTRI',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL DTFTRI('N','L','N',-1,a,info)
      CALL CHKXER('DTFTRI',INFot,NOUt,LERr,OK)
!
      SRNamt = 'DTFTTR'
      INFot = 1
      CALL DTFTTR('/','U',0,a,b,1,info)
      CALL CHKXER('DTFTTR',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL DTFTTR('N','/',0,a,b,1,info)
      CALL CHKXER('DTFTTR',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL DTFTTR('N','U',-1,a,b,1,info)
      CALL CHKXER('DTFTTR',INFot,NOUt,LERr,OK)
      INFot = 6
      CALL DTFTTR('N','U',0,a,b,0,info)
      CALL CHKXER('DTFTTR',INFot,NOUt,LERr,OK)
!
      SRNamt = 'DTRTTF'
      INFot = 1
      CALL DTRTTF('/','U',0,a,1,b,info)
      CALL CHKXER('DTRTTF',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL DTRTTF('N','/',0,a,1,b,info)
      CALL CHKXER('DTRTTF',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL DTRTTF('N','U',-1,a,1,b,info)
      CALL CHKXER('DTRTTF',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL DTRTTF('N','U',0,a,0,b,info)
      CALL CHKXER('DTRTTF',INFot,NOUt,LERr,OK)
!
      SRNamt = 'DTFTTP'
      INFot = 1
      CALL DTFTTP('/','U',0,a,b,info)
      CALL CHKXER('DTFTTP',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL DTFTTP('N','/',0,a,b,info)
      CALL CHKXER('DTFTTP',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL DTFTTP('N','U',-1,a,b,info)
      CALL CHKXER('DTFTTP',INFot,NOUt,LERr,OK)
!
      SRNamt = 'DTPTTF'
      INFot = 1
      CALL DTPTTF('/','U',0,a,b,info)
      CALL CHKXER('DTPTTF',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL DTPTTF('N','/',0,a,b,info)
      CALL CHKXER('DTPTTF',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL DTPTTF('N','U',-1,a,b,info)
      CALL CHKXER('DTPTTF',INFot,NOUt,LERr,OK)
!
      SRNamt = 'DTRTTP'
      INFot = 1
      CALL DTRTTP('/',0,a,1,b,info)
      CALL CHKXER('DTRTTP',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL DTRTTP('U',-1,a,1,b,info)
      CALL CHKXER('DTRTTP',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL DTRTTP('U',0,a,0,b,info)
      CALL CHKXER('DTRTTP',INFot,NOUt,LERr,OK)
!
      SRNamt = 'DTPTTR'
      INFot = 1
      CALL DTPTTR('/',0,a,b,1,info)
      CALL CHKXER('DTPTTR',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL DTPTTR('U',-1,a,b,1,info)
      CALL CHKXER('DTPTTR',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL DTPTTR('U',0,a,b,0,info)
      CALL CHKXER('DTPTTR',INFot,NOUt,LERr,OK)
!
      SRNamt = 'DSFRK '
      INFot = 1
      CALL DSFRK('/','U','N',0,0,alpha,a,1,beta,b)
      CALL CHKXER('DSFRK ',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL DSFRK('N','/','N',0,0,alpha,a,1,beta,b)
      CALL CHKXER('DSFRK ',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL DSFRK('N','U','/',0,0,alpha,a,1,beta,b)
      CALL CHKXER('DSFRK ',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL DSFRK('N','U','N',-1,0,alpha,a,1,beta,b)
      CALL CHKXER('DSFRK ',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL DSFRK('N','U','N',0,-1,alpha,a,1,beta,b)
      CALL CHKXER('DSFRK ',INFot,NOUt,LERr,OK)
      INFot = 8
      CALL DSFRK('N','U','N',0,0,alpha,a,0,beta,b)
      CALL CHKXER('DSFRK ',INFot,NOUt,LERr,OK)
!
!     Print a summary line.
!
      IF ( OK ) THEN
         WRITE (NOUt,FMT=99001)
      ELSE
         WRITE (NOUt,FMT=99002)
      ENDIF
!
99001 FORMAT (1X,'DOUBLE PRECISION RFP routines passed the tests of ',  &
     &        'the error exits')
99002 FORMAT (' *** RFP routines failed the tests of the error ',       &
     &        'exits ***')
!
!     End of derrrfp
!
      END SUBROUTINE derrrfp
