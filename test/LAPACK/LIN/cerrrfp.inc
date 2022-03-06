!*==cerrrfp.f90  processed by SPAG 7.51RB at 17:34 on  4 Mar 2022
!> \brief \b cerrrfp
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE cerrrfp( NUNIT )
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
!> cerrrfp tests the error exits for the COMPLEX driver routines
!> for solving linear systems of equations.
!>
!> cdrvrfp tests the COMPLEX LAPACK RFP routines:
!>     CTFSM, CTFTRI, CHFRK, CTFTTP, CTFTTR, CPFTRF, CPFTRS, CTPTTF,
!>     CTPTTR, CTRTTF, and CTRTTP
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
!> \ingroup complex_lin
!
!  =====================================================================
      SUBROUTINE cerrrfp(Nunit)
      IMPLICIT NONE
!*--CERRRFP56
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
      COMPLEX alpha , beta
      real    r_alpha
      real    r_beta
!     ..
!     .. Local Arrays ..
      COMPLEX a(1,1) , b(1,1)
!     ..
!     .. External Subroutines ..
EXTERNAL    CHKXER,   CTFSM,   CTFTRI,   CHFRK,   CTFTTP,   CTFTTR
EXTERNAL    CPFTRI,   CPFTRF,   CPFTRS,   CTPTTF,   CTPTTR,   CTRTTF
EXTERNAL    CTRTTP
!     ..
!     .. Scalars in Common ..
      LOGICAL LERr , OK
      CHARACTER*32 SRNamt
      INTEGER INFot , NOUt
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC CMPLX
!     ..
!     .. Common blocks ..
      COMMON /INFOC / INFot , NOUt , OK , LERr
      COMMON /SRNAMC/ SRNamt
!     ..
!     .. Executable Statements ..
!
      NOUt = Nunit
      OK = .TRUE.
      a(1,1) = CMPLX(1.D0,1.D0)
      b(1,1) = CMPLX(1.D0,1.D0)
      alpha = CMPLX(1.D0,1.D0)
      r_alpha=1.0
      beta = CMPLX(1.D0,1.D0)
      r_beta=1.0
!
      SRNamt = 'CPFTRF'
      INFot = 1
      CALL CPFTRF('/','U',0,a,info)
      CALL CHKXER('CPFTRF',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL CPFTRF('N','/',0,a,info)
      CALL CHKXER('CPFTRF',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL CPFTRF('N','U',-1,a,info)
      CALL CHKXER('CPFTRF',INFot,NOUt,LERr,OK)
!
      SRNamt = 'CPFTRS'
      INFot = 1
      CALL CPFTRS('/','U',0,0,a,b,1,info)
      CALL CHKXER('CPFTRS',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL CPFTRS('N','/',0,0,a,b,1,info)
      CALL CHKXER('CPFTRS',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL CPFTRS('N','U',-1,0,a,b,1,info)
      CALL CHKXER('CPFTRS',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL CPFTRS('N','U',0,-1,a,b,1,info)
      CALL CHKXER('CPFTRS',INFot,NOUt,LERr,OK)
      INFot = 7
      CALL CPFTRS('N','U',0,0,a,b,0,info)
      CALL CHKXER('CPFTRS',INFot,NOUt,LERr,OK)
!
      SRNamt = 'CPFTRI'
      INFot = 1
      CALL CPFTRI('/','U',0,a,info)
      CALL CHKXER('CPFTRI',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL CPFTRI('N','/',0,a,info)
      CALL CHKXER('CPFTRI',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL CPFTRI('N','U',-1,a,info)
      CALL CHKXER('CPFTRI',INFot,NOUt,LERr,OK)
!
      SRNamt = 'CTFSM '
      INFot = 1
      CALL CTFSM('/','L','U','C','U',0,0,alpha,a,b,1)
      CALL CHKXER('CTFSM ',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL CTFSM('N','/','U','C','U',0,0,alpha,a,b,1)
      CALL CHKXER('CTFSM ',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL CTFSM('N','L','/','C','U',0,0,alpha,a,b,1)
      CALL CHKXER('CTFSM ',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL CTFSM('N','L','U','/','U',0,0,alpha,a,b,1)
      CALL CHKXER('CTFSM ',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL CTFSM('N','L','U','C','/',0,0,alpha,a,b,1)
      CALL CHKXER('CTFSM ',INFot,NOUt,LERr,OK)
      INFot = 6
      CALL CTFSM('N','L','U','C','U',-1,0,alpha,a,b,1)
      CALL CHKXER('CTFSM ',INFot,NOUt,LERr,OK)
      INFot = 7
      CALL CTFSM('N','L','U','C','U',0,-1,alpha,a,b,1)
      CALL CHKXER('CTFSM ',INFot,NOUt,LERr,OK)
      INFot = 11
      CALL CTFSM('N','L','U','C','U',0,0,alpha,a,b,0)
      CALL CHKXER('CTFSM ',INFot,NOUt,LERr,OK)
!
      SRNamt = 'CTFTRI'
      INFot = 1
      CALL CTFTRI('/','L','N',0,a,info)
      CALL CHKXER('CTFTRI',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL CTFTRI('N','/','N',0,a,info)
      CALL CHKXER('CTFTRI',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL CTFTRI('N','L','/',0,a,info)
      CALL CHKXER('CTFTRI',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL CTFTRI('N','L','N',-1,a,info)
      CALL CHKXER('CTFTRI',INFot,NOUt,LERr,OK)
!
      SRNamt = 'CTFTTR'
      INFot = 1
      CALL CTFTTR('/','U',0,a,b,1,info)
      CALL CHKXER('CTFTTR',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL CTFTTR('N','/',0,a,b,1,info)
      CALL CHKXER('CTFTTR',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL CTFTTR('N','U',-1,a,b,1,info)
      CALL CHKXER('CTFTTR',INFot,NOUt,LERr,OK)
      INFot = 6
      CALL CTFTTR('N','U',0,a,b,0,info)
      CALL CHKXER('CTFTTR',INFot,NOUt,LERr,OK)
!
      SRNamt = 'CTRTTF'
      INFot = 1
      CALL CTRTTF('/','U',0,a,1,b,info)
      CALL CHKXER('CTRTTF',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL CTRTTF('N','/',0,a,1,b,info)
      CALL CHKXER('CTRTTF',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL CTRTTF('N','U',-1,a,1,b,info)
      CALL CHKXER('CTRTTF',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL CTRTTF('N','U',0,a,0,b,info)
      CALL CHKXER('CTRTTF',INFot,NOUt,LERr,OK)
!
      SRNamt = 'CTFTTP'
      INFot = 1
      CALL CTFTTP('/','U',0,a,b,info)
      CALL CHKXER('CTFTTP',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL CTFTTP('N','/',0,a,b,info)
      CALL CHKXER('CTFTTP',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL CTFTTP('N','U',-1,a,b,info)
      CALL CHKXER('CTFTTP',INFot,NOUt,LERr,OK)
!
      SRNamt = 'CTPTTF'
      INFot = 1
      CALL CTPTTF('/','U',0,a,b,info)
      CALL CHKXER('CTPTTF',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL CTPTTF('N','/',0,a,b,info)
      CALL CHKXER('CTPTTF',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL CTPTTF('N','U',-1,a,b,info)
      CALL CHKXER('CTPTTF',INFot,NOUt,LERr,OK)
!
      SRNamt = 'CTRTTP'
      INFot = 1
      CALL CTRTTP('/',0,a,1,b,info)
      CALL CHKXER('CTRTTP',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL CTRTTP('U',-1,a,1,b,info)
      CALL CHKXER('CTRTTP',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL CTRTTP('U',0,a,0,b,info)
      CALL CHKXER('CTRTTP',INFot,NOUt,LERr,OK)
!
      SRNamt = 'CTPTTR'
      INFot = 1
      CALL CTPTTR('/',0,a,b,1,info)
      CALL CHKXER('CTPTTR',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL CTPTTR('U',-1,a,b,1,info)
      CALL CHKXER('CTPTTR',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL CTPTTR('U',0,a,b,0,info)
      CALL CHKXER('CTPTTR',INFot,NOUt,LERr,OK)
!
      SRNamt = 'CHFRK '
      INFot = 1
      CALL CHFRK('/','U','N',0,0,r_alpha,a,1,r_beta,b)
      CALL CHKXER('CHFRK ',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL CHFRK('N','/','N',0,0,r_alpha,a,1,r_beta,b)
      CALL CHKXER('CHFRK ',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL CHFRK('N','U','/',0,0,r_alpha,a,1,r_beta,b)
      CALL CHKXER('CHFRK ',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL CHFRK('N','U','N',-1,0,r_alpha,a,1,r_beta,b)
      CALL CHKXER('CHFRK ',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL CHFRK('N','U','N',0,-1,r_alpha,a,1,r_beta,b)
      CALL CHKXER('CHFRK ',INFot,NOUt,LERr,OK)
      INFot = 8
      CALL CHFRK('N','U','N',0,0,r_alpha,a,0,r_beta,b)
      CALL CHKXER('CHFRK ',INFot,NOUt,LERr,OK)
!
!     Print a summary line.
!
      IF ( OK ) THEN
         WRITE (NOUt,FMT=99001)
      ELSE
         WRITE (NOUt,FMT=99002)
      ENDIF
!
99001 FORMAT (1X,'COMPLEX RFP routines passed the tests of the ',       &
     &        'error exits')
99002 FORMAT (' *** RFP routines failed the tests of the error ',       &
     &        'exits ***')
!
!     End of cerrrfp
!
      END SUBROUTINE cerrrfp
