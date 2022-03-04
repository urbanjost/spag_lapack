!*==serrrfp.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b SERRRFP
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE SERRRFP( NUNIT )
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
!> SERRRFP tests the error exits for the REAL driver routines
!> for solving linear systems of equations.
!>
!> SDRVRFP tests the REAL LAPACK RFP routines:
!>     STFSM, STFTRI, SSFRK, STFTTP, STFTTR, SPFTRF, SPFTRS, STPTTF,
!>     STPTTR, STRTTF, and STRTTP
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
!> \ingroup single_lin
!
!  =====================================================================
      SUBROUTINE SERRRFP(Nunit)
      IMPLICIT NONE
!*--SERRRFP56
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
      REAL alpha , beta
!     ..
!     .. Local Arrays ..
      REAL a(1,1) , b(1,1)
!     ..
!     .. External Subroutines ..
      EXTERNAL CHKXER , STFSM , STFTRI , SSFRK , STFTTP , STFTTR ,      &
     &         SPFTRI , SPFTRF , SPFTRS , STPTTF , STPTTR , STRTTF ,    &
     &         STRTTP
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
      a(1,1) = 1.0E+0
      b(1,1) = 1.0E+0
      alpha = 1.0E+0
      beta = 1.0E+0
!
      SRNamt = 'SPFTRF'
      INFot = 1
      CALL SPFTRF('/','U',0,a,info)
      CALL CHKXER('SPFTRF',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL SPFTRF('N','/',0,a,info)
      CALL CHKXER('SPFTRF',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL SPFTRF('N','U',-1,a,info)
      CALL CHKXER('SPFTRF',INFot,NOUt,LERr,OK)
!
      SRNamt = 'SPFTRS'
      INFot = 1
      CALL SPFTRS('/','U',0,0,a,b,1,info)
      CALL CHKXER('SPFTRS',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL SPFTRS('N','/',0,0,a,b,1,info)
      CALL CHKXER('SPFTRS',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL SPFTRS('N','U',-1,0,a,b,1,info)
      CALL CHKXER('SPFTRS',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL SPFTRS('N','U',0,-1,a,b,1,info)
      CALL CHKXER('SPFTRS',INFot,NOUt,LERr,OK)
      INFot = 7
      CALL SPFTRS('N','U',0,0,a,b,0,info)
      CALL CHKXER('SPFTRS',INFot,NOUt,LERr,OK)
!
      SRNamt = 'SPFTRI'
      INFot = 1
      CALL SPFTRI('/','U',0,a,info)
      CALL CHKXER('SPFTRI',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL SPFTRI('N','/',0,a,info)
      CALL CHKXER('SPFTRI',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL SPFTRI('N','U',-1,a,info)
      CALL CHKXER('SPFTRI',INFot,NOUt,LERr,OK)
!
      SRNamt = 'STFSM '
      INFot = 1
      CALL STFSM('/','L','U','T','U',0,0,alpha,a,b,1)
      CALL CHKXER('STFSM ',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL STFSM('N','/','U','T','U',0,0,alpha,a,b,1)
      CALL CHKXER('STFSM ',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL STFSM('N','L','/','T','U',0,0,alpha,a,b,1)
      CALL CHKXER('STFSM ',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL STFSM('N','L','U','/','U',0,0,alpha,a,b,1)
      CALL CHKXER('STFSM ',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL STFSM('N','L','U','T','/',0,0,alpha,a,b,1)
      CALL CHKXER('STFSM ',INFot,NOUt,LERr,OK)
      INFot = 6
      CALL STFSM('N','L','U','T','U',-1,0,alpha,a,b,1)
      CALL CHKXER('STFSM ',INFot,NOUt,LERr,OK)
      INFot = 7
      CALL STFSM('N','L','U','T','U',0,-1,alpha,a,b,1)
      CALL CHKXER('STFSM ',INFot,NOUt,LERr,OK)
      INFot = 11
      CALL STFSM('N','L','U','T','U',0,0,alpha,a,b,0)
      CALL CHKXER('STFSM ',INFot,NOUt,LERr,OK)
!
      SRNamt = 'STFTRI'
      INFot = 1
      CALL STFTRI('/','L','N',0,a,info)
      CALL CHKXER('STFTRI',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL STFTRI('N','/','N',0,a,info)
      CALL CHKXER('STFTRI',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL STFTRI('N','L','/',0,a,info)
      CALL CHKXER('STFTRI',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL STFTRI('N','L','N',-1,a,info)
      CALL CHKXER('STFTRI',INFot,NOUt,LERr,OK)
!
      SRNamt = 'STFTTR'
      INFot = 1
      CALL STFTTR('/','U',0,a,b,1,info)
      CALL CHKXER('STFTTR',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL STFTTR('N','/',0,a,b,1,info)
      CALL CHKXER('STFTTR',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL STFTTR('N','U',-1,a,b,1,info)
      CALL CHKXER('STFTTR',INFot,NOUt,LERr,OK)
      INFot = 6
      CALL STFTTR('N','U',0,a,b,0,info)
      CALL CHKXER('STFTTR',INFot,NOUt,LERr,OK)
!
      SRNamt = 'STRTTF'
      INFot = 1
      CALL STRTTF('/','U',0,a,1,b,info)
      CALL CHKXER('STRTTF',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL STRTTF('N','/',0,a,1,b,info)
      CALL CHKXER('STRTTF',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL STRTTF('N','U',-1,a,1,b,info)
      CALL CHKXER('STRTTF',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL STRTTF('N','U',0,a,0,b,info)
      CALL CHKXER('STRTTF',INFot,NOUt,LERr,OK)
!
      SRNamt = 'STFTTP'
      INFot = 1
      CALL STFTTP('/','U',0,a,b,info)
      CALL CHKXER('STFTTP',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL STFTTP('N','/',0,a,b,info)
      CALL CHKXER('STFTTP',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL STFTTP('N','U',-1,a,b,info)
      CALL CHKXER('STFTTP',INFot,NOUt,LERr,OK)
!
      SRNamt = 'STPTTF'
      INFot = 1
      CALL STPTTF('/','U',0,a,b,info)
      CALL CHKXER('STPTTF',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL STPTTF('N','/',0,a,b,info)
      CALL CHKXER('STPTTF',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL STPTTF('N','U',-1,a,b,info)
      CALL CHKXER('STPTTF',INFot,NOUt,LERr,OK)
!
      SRNamt = 'STRTTP'
      INFot = 1
      CALL STRTTP('/',0,a,1,b,info)
      CALL CHKXER('STRTTP',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL STRTTP('U',-1,a,1,b,info)
      CALL CHKXER('STRTTP',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL STRTTP('U',0,a,0,b,info)
      CALL CHKXER('STRTTP',INFot,NOUt,LERr,OK)
!
      SRNamt = 'STPTTR'
      INFot = 1
      CALL STPTTR('/',0,a,b,1,info)
      CALL CHKXER('STPTTR',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL STPTTR('U',-1,a,b,1,info)
      CALL CHKXER('STPTTR',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL STPTTR('U',0,a,b,0,info)
      CALL CHKXER('STPTTR',INFot,NOUt,LERr,OK)
!
      SRNamt = 'SSFRK '
      INFot = 1
      CALL SSFRK('/','U','N',0,0,alpha,a,1,beta,b)
      CALL CHKXER('SSFRK ',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL SSFRK('N','/','N',0,0,alpha,a,1,beta,b)
      CALL CHKXER('SSFRK ',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL SSFRK('N','U','/',0,0,alpha,a,1,beta,b)
      CALL CHKXER('SSFRK ',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL SSFRK('N','U','N',-1,0,alpha,a,1,beta,b)
      CALL CHKXER('SSFRK ',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL SSFRK('N','U','N',0,-1,alpha,a,1,beta,b)
      CALL CHKXER('SSFRK ',INFot,NOUt,LERr,OK)
      INFot = 8
      CALL SSFRK('N','U','N',0,0,alpha,a,0,beta,b)
      CALL CHKXER('SSFRK ',INFot,NOUt,LERr,OK)
!
!     Print a summary line.
!
      IF ( OK ) THEN
         WRITE (NOUt,FMT=99001)
      ELSE
         WRITE (NOUt,FMT=99002)
      ENDIF
!
99001 FORMAT (1X,'REAL RFP routines passed the tests of ',              &
     &        'the error exits')
99002 FORMAT (' *** RFP routines failed the tests of the error ',       &
     &        'exits ***')
!
!     End of SERRRFP
!
      END SUBROUTINE SERRRFP
