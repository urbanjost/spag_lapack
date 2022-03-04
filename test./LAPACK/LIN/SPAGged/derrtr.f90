!*==derrtr.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b DERRTR
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DERRTR( PATH, NUNIT )
!
!       .. Scalar Arguments ..
!       CHARACTER*3        PATH
!       INTEGER            NUNIT
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DERRTR tests the error exits for the DOUBLE PRECISION triangular
!> routines.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] PATH
!> \verbatim
!>          PATH is CHARACTER*3
!>          The LAPACK path name for the routines to be tested.
!> \endverbatim
!>
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
      SUBROUTINE DERRTR(Path,Nunit)
      IMPLICIT NONE
!*--DERRTR59
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER*3 Path
      INTEGER Nunit
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      INTEGER NMAX
      PARAMETER (NMAX=2)
!     ..
!     .. Local Scalars ..
      CHARACTER*2 c2
      INTEGER info
      DOUBLE PRECISION rcond , scale
!     ..
!     .. Local Arrays ..
      INTEGER iw(NMAX)
      DOUBLE PRECISION a(NMAX,NMAX) , b(NMAX) , r1(NMAX) , r2(NMAX) ,   &
     &                 w(NMAX) , x(NMAX)
!     ..
!     .. External Functions ..
      LOGICAL LSAMEN
      EXTERNAL LSAMEN
!     ..
!     .. External Subroutines ..
      EXTERNAL ALAESM , CHKXER , DLATBS , DLATPS , DLATRS , DTBCON ,    &
     &         DTBRFS , DTBTRS , DTPCON , DTPRFS , DTPTRI , DTPTRS ,    &
     &         DTRCON , DTRRFS , DTRTI2 , DTRTRI , DTRTRS
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
      WRITE (NOUt,FMT=*)
      c2 = Path(2:3)
      a(1,1) = 1.D0
      a(1,2) = 2.D0
      a(2,2) = 3.D0
      a(2,1) = 4.D0
      OK = .TRUE.
!
      IF ( LSAMEN(2,c2,'TR') ) THEN
!
!        Test error exits for the general triangular routines.
!
!        DTRTRI
!
         SRNamt = 'DTRTRI'
         INFot = 1
         CALL DTRTRI('/','N',0,a,1,info)
         CALL CHKXER('DTRTRI',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DTRTRI('U','/',0,a,1,info)
         CALL CHKXER('DTRTRI',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DTRTRI('U','N',-1,a,1,info)
         CALL CHKXER('DTRTRI',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL DTRTRI('U','N',2,a,1,info)
         CALL CHKXER('DTRTRI',INFot,NOUt,LERr,OK)
!
!        DTRTI2
!
         SRNamt = 'DTRTI2'
         INFot = 1
         CALL DTRTI2('/','N',0,a,1,info)
         CALL CHKXER('DTRTI2',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DTRTI2('U','/',0,a,1,info)
         CALL CHKXER('DTRTI2',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DTRTI2('U','N',-1,a,1,info)
         CALL CHKXER('DTRTI2',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL DTRTI2('U','N',2,a,1,info)
         CALL CHKXER('DTRTI2',INFot,NOUt,LERr,OK)
!
!        DTRTRS
!
         SRNamt = 'DTRTRS'
         INFot = 1
         CALL DTRTRS('/','N','N',0,0,a,1,x,1,info)
         CALL CHKXER('DTRTRS',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DTRTRS('U','/','N',0,0,a,1,x,1,info)
         CALL CHKXER('DTRTRS',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DTRTRS('U','N','/',0,0,a,1,x,1,info)
         CALL CHKXER('DTRTRS',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL DTRTRS('U','N','N',-1,0,a,1,x,1,info)
         CALL CHKXER('DTRTRS',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL DTRTRS('U','N','N',0,-1,a,1,x,1,info)
         CALL CHKXER('DTRTRS',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL DTRTRS('U','N','N',2,1,a,1,x,2,info)
         CALL CHKXER('DTRTRS',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL DTRTRS('U','N','N',2,1,a,2,x,1,info)
         CALL CHKXER('DTRTRS',INFot,NOUt,LERr,OK)
!
!        DTRRFS
!
         SRNamt = 'DTRRFS'
         INFot = 1
         CALL DTRRFS('/','N','N',0,0,a,1,b,1,x,1,r1,r2,w,iw,info)
         CALL CHKXER('DTRRFS',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DTRRFS('U','/','N',0,0,a,1,b,1,x,1,r1,r2,w,iw,info)
         CALL CHKXER('DTRRFS',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DTRRFS('U','N','/',0,0,a,1,b,1,x,1,r1,r2,w,iw,info)
         CALL CHKXER('DTRRFS',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL DTRRFS('U','N','N',-1,0,a,1,b,1,x,1,r1,r2,w,iw,info)
         CALL CHKXER('DTRRFS',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL DTRRFS('U','N','N',0,-1,a,1,b,1,x,1,r1,r2,w,iw,info)
         CALL CHKXER('DTRRFS',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL DTRRFS('U','N','N',2,1,a,1,b,2,x,2,r1,r2,w,iw,info)
         CALL CHKXER('DTRRFS',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL DTRRFS('U','N','N',2,1,a,2,b,1,x,2,r1,r2,w,iw,info)
         CALL CHKXER('DTRRFS',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL DTRRFS('U','N','N',2,1,a,2,b,2,x,1,r1,r2,w,iw,info)
         CALL CHKXER('DTRRFS',INFot,NOUt,LERr,OK)
!
!        DTRCON
!
         SRNamt = 'DTRCON'
         INFot = 1
         CALL DTRCON('/','U','N',0,a,1,rcond,w,iw,info)
         CALL CHKXER('DTRCON',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DTRCON('1','/','N',0,a,1,rcond,w,iw,info)
         CALL CHKXER('DTRCON',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DTRCON('1','U','/',0,a,1,rcond,w,iw,info)
         CALL CHKXER('DTRCON',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL DTRCON('1','U','N',-1,a,1,rcond,w,iw,info)
         CALL CHKXER('DTRCON',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL DTRCON('1','U','N',2,a,1,rcond,w,iw,info)
         CALL CHKXER('DTRCON',INFot,NOUt,LERr,OK)
!
!        DLATRS
!
         SRNamt = 'DLATRS'
         INFot = 1
         CALL DLATRS('/','N','N','N',0,a,1,x,scale,w,info)
         CALL CHKXER('DLATRS',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DLATRS('U','/','N','N',0,a,1,x,scale,w,info)
         CALL CHKXER('DLATRS',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DLATRS('U','N','/','N',0,a,1,x,scale,w,info)
         CALL CHKXER('DLATRS',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL DLATRS('U','N','N','/',0,a,1,x,scale,w,info)
         CALL CHKXER('DLATRS',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL DLATRS('U','N','N','N',-1,a,1,x,scale,w,info)
         CALL CHKXER('DLATRS',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL DLATRS('U','N','N','N',2,a,1,x,scale,w,info)
         CALL CHKXER('DLATRS',INFot,NOUt,LERr,OK)
!
      ELSEIF ( LSAMEN(2,c2,'TP') ) THEN
!
!        Test error exits for the packed triangular routines.
!
!        DTPTRI
!
         SRNamt = 'DTPTRI'
         INFot = 1
         CALL DTPTRI('/','N',0,a,info)
         CALL CHKXER('DTPTRI',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DTPTRI('U','/',0,a,info)
         CALL CHKXER('DTPTRI',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DTPTRI('U','N',-1,a,info)
         CALL CHKXER('DTPTRI',INFot,NOUt,LERr,OK)
!
!        DTPTRS
!
         SRNamt = 'DTPTRS'
         INFot = 1
         CALL DTPTRS('/','N','N',0,0,a,x,1,info)
         CALL CHKXER('DTPTRS',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DTPTRS('U','/','N',0,0,a,x,1,info)
         CALL CHKXER('DTPTRS',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DTPTRS('U','N','/',0,0,a,x,1,info)
         CALL CHKXER('DTPTRS',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL DTPTRS('U','N','N',-1,0,a,x,1,info)
         CALL CHKXER('DTPTRS',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL DTPTRS('U','N','N',0,-1,a,x,1,info)
         CALL CHKXER('DTPTRS',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL DTPTRS('U','N','N',2,1,a,x,1,info)
         CALL CHKXER('DTPTRS',INFot,NOUt,LERr,OK)
!
!        DTPRFS
!
         SRNamt = 'DTPRFS'
         INFot = 1
         CALL DTPRFS('/','N','N',0,0,a,b,1,x,1,r1,r2,w,iw,info)
         CALL CHKXER('DTPRFS',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DTPRFS('U','/','N',0,0,a,b,1,x,1,r1,r2,w,iw,info)
         CALL CHKXER('DTPRFS',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DTPRFS('U','N','/',0,0,a,b,1,x,1,r1,r2,w,iw,info)
         CALL CHKXER('DTPRFS',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL DTPRFS('U','N','N',-1,0,a,b,1,x,1,r1,r2,w,iw,info)
         CALL CHKXER('DTPRFS',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL DTPRFS('U','N','N',0,-1,a,b,1,x,1,r1,r2,w,iw,info)
         CALL CHKXER('DTPRFS',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL DTPRFS('U','N','N',2,1,a,b,1,x,2,r1,r2,w,iw,info)
         CALL CHKXER('DTPRFS',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL DTPRFS('U','N','N',2,1,a,b,2,x,1,r1,r2,w,iw,info)
         CALL CHKXER('DTPRFS',INFot,NOUt,LERr,OK)
!
!        DTPCON
!
         SRNamt = 'DTPCON'
         INFot = 1
         CALL DTPCON('/','U','N',0,a,rcond,w,iw,info)
         CALL CHKXER('DTPCON',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DTPCON('1','/','N',0,a,rcond,w,iw,info)
         CALL CHKXER('DTPCON',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DTPCON('1','U','/',0,a,rcond,w,iw,info)
         CALL CHKXER('DTPCON',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL DTPCON('1','U','N',-1,a,rcond,w,iw,info)
         CALL CHKXER('DTPCON',INFot,NOUt,LERr,OK)
!
!        DLATPS
!
         SRNamt = 'DLATPS'
         INFot = 1
         CALL DLATPS('/','N','N','N',0,a,x,scale,w,info)
         CALL CHKXER('DLATPS',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DLATPS('U','/','N','N',0,a,x,scale,w,info)
         CALL CHKXER('DLATPS',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DLATPS('U','N','/','N',0,a,x,scale,w,info)
         CALL CHKXER('DLATPS',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL DLATPS('U','N','N','/',0,a,x,scale,w,info)
         CALL CHKXER('DLATPS',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL DLATPS('U','N','N','N',-1,a,x,scale,w,info)
         CALL CHKXER('DLATPS',INFot,NOUt,LERr,OK)
!
      ELSEIF ( LSAMEN(2,c2,'TB') ) THEN
!
!        Test error exits for the banded triangular routines.
!
!        DTBTRS
!
         SRNamt = 'DTBTRS'
         INFot = 1
         CALL DTBTRS('/','N','N',0,0,0,a,1,x,1,info)
         CALL CHKXER('DTBTRS',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DTBTRS('U','/','N',0,0,0,a,1,x,1,info)
         CALL CHKXER('DTBTRS',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DTBTRS('U','N','/',0,0,0,a,1,x,1,info)
         CALL CHKXER('DTBTRS',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL DTBTRS('U','N','N',-1,0,0,a,1,x,1,info)
         CALL CHKXER('DTBTRS',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL DTBTRS('U','N','N',0,-1,0,a,1,x,1,info)
         CALL CHKXER('DTBTRS',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL DTBTRS('U','N','N',0,0,-1,a,1,x,1,info)
         CALL CHKXER('DTBTRS',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL DTBTRS('U','N','N',2,1,1,a,1,x,2,info)
         CALL CHKXER('DTBTRS',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL DTBTRS('U','N','N',2,0,1,a,1,x,1,info)
         CALL CHKXER('DTBTRS',INFot,NOUt,LERr,OK)
!
!        DTBRFS
!
         SRNamt = 'DTBRFS'
         INFot = 1
         CALL DTBRFS('/','N','N',0,0,0,a,1,b,1,x,1,r1,r2,w,iw,info)
         CALL CHKXER('DTBRFS',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DTBRFS('U','/','N',0,0,0,a,1,b,1,x,1,r1,r2,w,iw,info)
         CALL CHKXER('DTBRFS',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DTBRFS('U','N','/',0,0,0,a,1,b,1,x,1,r1,r2,w,iw,info)
         CALL CHKXER('DTBRFS',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL DTBRFS('U','N','N',-1,0,0,a,1,b,1,x,1,r1,r2,w,iw,info)
         CALL CHKXER('DTBRFS',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL DTBRFS('U','N','N',0,-1,0,a,1,b,1,x,1,r1,r2,w,iw,info)
         CALL CHKXER('DTBRFS',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL DTBRFS('U','N','N',0,0,-1,a,1,b,1,x,1,r1,r2,w,iw,info)
         CALL CHKXER('DTBRFS',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL DTBRFS('U','N','N',2,1,1,a,1,b,2,x,2,r1,r2,w,iw,info)
         CALL CHKXER('DTBRFS',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL DTBRFS('U','N','N',2,1,1,a,2,b,1,x,2,r1,r2,w,iw,info)
         CALL CHKXER('DTBRFS',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL DTBRFS('U','N','N',2,1,1,a,2,b,2,x,1,r1,r2,w,iw,info)
         CALL CHKXER('DTBRFS',INFot,NOUt,LERr,OK)
!
!        DTBCON
!
         SRNamt = 'DTBCON'
         INFot = 1
         CALL DTBCON('/','U','N',0,0,a,1,rcond,w,iw,info)
         CALL CHKXER('DTBCON',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DTBCON('1','/','N',0,0,a,1,rcond,w,iw,info)
         CALL CHKXER('DTBCON',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DTBCON('1','U','/',0,0,a,1,rcond,w,iw,info)
         CALL CHKXER('DTBCON',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL DTBCON('1','U','N',-1,0,a,1,rcond,w,iw,info)
         CALL CHKXER('DTBCON',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL DTBCON('1','U','N',0,-1,a,1,rcond,w,iw,info)
         CALL CHKXER('DTBCON',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL DTBCON('1','U','N',2,1,a,1,rcond,w,iw,info)
         CALL CHKXER('DTBCON',INFot,NOUt,LERr,OK)
!
!        DLATBS
!
         SRNamt = 'DLATBS'
         INFot = 1
         CALL DLATBS('/','N','N','N',0,0,a,1,x,scale,w,info)
         CALL CHKXER('DLATBS',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DLATBS('U','/','N','N',0,0,a,1,x,scale,w,info)
         CALL CHKXER('DLATBS',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DLATBS('U','N','/','N',0,0,a,1,x,scale,w,info)
         CALL CHKXER('DLATBS',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL DLATBS('U','N','N','/',0,0,a,1,x,scale,w,info)
         CALL CHKXER('DLATBS',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL DLATBS('U','N','N','N',-1,0,a,1,x,scale,w,info)
         CALL CHKXER('DLATBS',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL DLATBS('U','N','N','N',1,-1,a,1,x,scale,w,info)
         CALL CHKXER('DLATBS',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL DLATBS('U','N','N','N',2,1,a,1,x,scale,w,info)
         CALL CHKXER('DLATBS',INFot,NOUt,LERr,OK)
      ENDIF
!
!     Print a summary line.
!
      CALL ALAESM(Path,OK,NOUt)
!
!
!     End of DERRTR
!
      END SUBROUTINE DERRTR
