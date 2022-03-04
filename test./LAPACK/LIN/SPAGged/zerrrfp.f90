!*==zerrrfp.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b ZERRRFP
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZERRRFP( NUNIT )
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
!> ZERRRFP tests the error exits for the COMPLEX*16 driver routines
!> for solving linear systems of equations.
!>
!> ZDRVRFP tests the COMPLEX*16 LAPACK RFP routines:
!>     ZTFSM, ZTFTRI, ZHFRK, ZTFTTP, ZTFTTR, ZPFTRF, ZPFTRS, ZTPTTF,
!>     ZTPTTR, ZTRTTF, and ZTRTTP
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
!> \ingroup complex16_lin
!
!  =====================================================================
      SUBROUTINE ZERRRFP(Nunit)
      IMPLICIT NONE
!*--ZERRRFP56
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
      COMPLEX*16 calpha
!     ..
!     .. Local Arrays ..
      COMPLEX*16 a(1,1) , b(1,1)
!     ..
!     .. External Subroutines ..
      EXTERNAL CHKXER , ZTFSM , ZTFTRI , ZHFRK , ZTFTTP , ZTFTTR ,      &
     &         ZPFTRI , ZPFTRF , ZPFTRS , ZTPTTF , ZTPTTR , ZTRTTF ,    &
     &         ZTRTTP
!     ..
!     .. Scalars in Common ..
      LOGICAL LERr , OK
      CHARACTER*32 SRNamt
      INTEGER INFot , NOUt
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC DCMPLX
!     ..
!     .. Common blocks ..
      COMMON /INFOC / INFot , NOUt , OK , LERr
      COMMON /SRNAMC/ SRNamt
!     ..
!     .. Executable Statements ..
!
      NOUt = Nunit
      OK = .TRUE.
      a(1,1) = DCMPLX(1.0D0,1.0D0)
      b(1,1) = DCMPLX(1.0D0,1.0D0)
      alpha = 1.0D0
      calpha = DCMPLX(1.0D0,1.0D0)
      beta = 1.0D0
!
      SRNamt = 'ZPFTRF'
      INFot = 1
      CALL ZPFTRF('/','U',0,a,info)
      CALL CHKXER('ZPFTRF',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL ZPFTRF('N','/',0,a,info)
      CALL CHKXER('ZPFTRF',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL ZPFTRF('N','U',-1,a,info)
      CALL CHKXER('ZPFTRF',INFot,NOUt,LERr,OK)
!
      SRNamt = 'ZPFTRS'
      INFot = 1
      CALL ZPFTRS('/','U',0,0,a,b,1,info)
      CALL CHKXER('ZPFTRS',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL ZPFTRS('N','/',0,0,a,b,1,info)
      CALL CHKXER('ZPFTRS',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL ZPFTRS('N','U',-1,0,a,b,1,info)
      CALL CHKXER('ZPFTRS',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL ZPFTRS('N','U',0,-1,a,b,1,info)
      CALL CHKXER('ZPFTRS',INFot,NOUt,LERr,OK)
      INFot = 7
      CALL ZPFTRS('N','U',0,0,a,b,0,info)
      CALL CHKXER('ZPFTRS',INFot,NOUt,LERr,OK)
!
      SRNamt = 'ZPFTRI'
      INFot = 1
      CALL ZPFTRI('/','U',0,a,info)
      CALL CHKXER('ZPFTRI',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL ZPFTRI('N','/',0,a,info)
      CALL CHKXER('ZPFTRI',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL ZPFTRI('N','U',-1,a,info)
      CALL CHKXER('ZPFTRI',INFot,NOUt,LERr,OK)
!
      SRNamt = 'ZTFSM '
      INFot = 1
      CALL ZTFSM('/','L','U','C','U',0,0,calpha,a,b,1)
      CALL CHKXER('ZTFSM ',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL ZTFSM('N','/','U','C','U',0,0,calpha,a,b,1)
      CALL CHKXER('ZTFSM ',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL ZTFSM('N','L','/','C','U',0,0,calpha,a,b,1)
      CALL CHKXER('ZTFSM ',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL ZTFSM('N','L','U','/','U',0,0,calpha,a,b,1)
      CALL CHKXER('ZTFSM ',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL ZTFSM('N','L','U','C','/',0,0,calpha,a,b,1)
      CALL CHKXER('ZTFSM ',INFot,NOUt,LERr,OK)
      INFot = 6
      CALL ZTFSM('N','L','U','C','U',-1,0,calpha,a,b,1)
      CALL CHKXER('ZTFSM ',INFot,NOUt,LERr,OK)
      INFot = 7
      CALL ZTFSM('N','L','U','C','U',0,-1,calpha,a,b,1)
      CALL CHKXER('ZTFSM ',INFot,NOUt,LERr,OK)
      INFot = 11
      CALL ZTFSM('N','L','U','C','U',0,0,calpha,a,b,0)
      CALL CHKXER('ZTFSM ',INFot,NOUt,LERr,OK)
!
      SRNamt = 'ZTFTRI'
      INFot = 1
      CALL ZTFTRI('/','L','N',0,a,info)
      CALL CHKXER('ZTFTRI',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL ZTFTRI('N','/','N',0,a,info)
      CALL CHKXER('ZTFTRI',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL ZTFTRI('N','L','/',0,a,info)
      CALL CHKXER('ZTFTRI',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL ZTFTRI('N','L','N',-1,a,info)
      CALL CHKXER('ZTFTRI',INFot,NOUt,LERr,OK)
!
      SRNamt = 'ZTFTTR'
      INFot = 1
      CALL ZTFTTR('/','U',0,a,b,1,info)
      CALL CHKXER('ZTFTTR',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL ZTFTTR('N','/',0,a,b,1,info)
      CALL CHKXER('ZTFTTR',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL ZTFTTR('N','U',-1,a,b,1,info)
      CALL CHKXER('ZTFTTR',INFot,NOUt,LERr,OK)
      INFot = 6
      CALL ZTFTTR('N','U',0,a,b,0,info)
      CALL CHKXER('ZTFTTR',INFot,NOUt,LERr,OK)
!
      SRNamt = 'ZTRTTF'
      INFot = 1
      CALL ZTRTTF('/','U',0,a,1,b,info)
      CALL CHKXER('ZTRTTF',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL ZTRTTF('N','/',0,a,1,b,info)
      CALL CHKXER('ZTRTTF',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL ZTRTTF('N','U',-1,a,1,b,info)
      CALL CHKXER('ZTRTTF',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL ZTRTTF('N','U',0,a,0,b,info)
      CALL CHKXER('ZTRTTF',INFot,NOUt,LERr,OK)
!
      SRNamt = 'ZTFTTP'
      INFot = 1
      CALL ZTFTTP('/','U',0,a,b,info)
      CALL CHKXER('ZTFTTP',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL ZTFTTP('N','/',0,a,b,info)
      CALL CHKXER('ZTFTTP',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL ZTFTTP('N','U',-1,a,b,info)
      CALL CHKXER('ZTFTTP',INFot,NOUt,LERr,OK)
!
      SRNamt = 'ZTPTTF'
      INFot = 1
      CALL ZTPTTF('/','U',0,a,b,info)
      CALL CHKXER('ZTPTTF',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL ZTPTTF('N','/',0,a,b,info)
      CALL CHKXER('ZTPTTF',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL ZTPTTF('N','U',-1,a,b,info)
      CALL CHKXER('ZTPTTF',INFot,NOUt,LERr,OK)
!
      SRNamt = 'ZTRTTP'
      INFot = 1
      CALL ZTRTTP('/',0,a,1,b,info)
      CALL CHKXER('ZTRTTP',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL ZTRTTP('U',-1,a,1,b,info)
      CALL CHKXER('ZTRTTP',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL ZTRTTP('U',0,a,0,b,info)
      CALL CHKXER('ZTRTTP',INFot,NOUt,LERr,OK)
!
      SRNamt = 'ZTPTTR'
      INFot = 1
      CALL ZTPTTR('/',0,a,b,1,info)
      CALL CHKXER('ZTPTTR',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL ZTPTTR('U',-1,a,b,1,info)
      CALL CHKXER('ZTPTTR',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL ZTPTTR('U',0,a,b,0,info)
      CALL CHKXER('ZTPTTR',INFot,NOUt,LERr,OK)
!
      SRNamt = 'ZHFRK '
      INFot = 1
      CALL ZHFRK('/','U','N',0,0,alpha,a,1,beta,b)
      CALL CHKXER('ZHFRK ',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL ZHFRK('N','/','N',0,0,alpha,a,1,beta,b)
      CALL CHKXER('ZHFRK ',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL ZHFRK('N','U','/',0,0,alpha,a,1,beta,b)
      CALL CHKXER('ZHFRK ',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL ZHFRK('N','U','N',-1,0,alpha,a,1,beta,b)
      CALL CHKXER('ZHFRK ',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL ZHFRK('N','U','N',0,-1,alpha,a,1,beta,b)
      CALL CHKXER('ZHFRK ',INFot,NOUt,LERr,OK)
      INFot = 8
      CALL ZHFRK('N','U','N',0,0,alpha,a,0,beta,b)
      CALL CHKXER('ZHFRK ',INFot,NOUt,LERr,OK)
!
!     Print a summary line.
!
      IF ( OK ) THEN
         WRITE (NOUt,FMT=99001)
      ELSE
         WRITE (NOUt,FMT=99002)
      ENDIF
!
99001 FORMAT (1X,'COMPLEX*16 RFP routines passed the tests of the ',    &
     &        'error exits')
99002 FORMAT (' *** RFP routines failed the tests of the error ',       &
     &        'exits ***')
!
!     End of ZERRRFP
!
      END SUBROUTINE ZERRRFP
