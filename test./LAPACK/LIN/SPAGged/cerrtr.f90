!*==cerrtr.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b CERRTR
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE CERRTR( PATH, NUNIT )
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
!> CERRTR tests the error exits for the COMPLEX triangular routines.
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
!> \ingroup complex_lin
!
!  =====================================================================
      SUBROUTINE CERRTR(Path,Nunit)
      IMPLICIT NONE
!*--CERRTR58
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
      REAL rcond , scale
!     ..
!     .. Local Arrays ..
      REAL r1(NMAX) , r2(NMAX) , rw(NMAX)
      COMPLEX a(NMAX,NMAX) , b(NMAX) , w(NMAX) , x(NMAX)
!     ..
!     .. External Functions ..
      LOGICAL LSAMEN
      EXTERNAL LSAMEN
!     ..
!     .. External Subroutines ..
      EXTERNAL ALAESM , CHKXER , CLATBS , CLATPS , CLATRS , CTBCON ,    &
     &         CTBRFS , CTBTRS , CTPCON , CTPRFS , CTPTRI , CTPTRS ,    &
     &         CTRCON , CTRRFS , CTRTI2 , CTRTRI , CTRTRS
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
      a(1,1) = 1.
      a(1,2) = 2.
      a(2,2) = 3.
      a(2,1) = 4.
      OK = .TRUE.
!
!     Test error exits for the general triangular routines.
!
      IF ( LSAMEN(2,c2,'TR') ) THEN
!
!        CTRTRI
!
         SRNamt = 'CTRTRI'
         INFot = 1
         CALL CTRTRI('/','N',0,a,1,info)
         CALL CHKXER('CTRTRI',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CTRTRI('U','/',0,a,1,info)
         CALL CHKXER('CTRTRI',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CTRTRI('U','N',-1,a,1,info)
         CALL CHKXER('CTRTRI',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL CTRTRI('U','N',2,a,1,info)
         CALL CHKXER('CTRTRI',INFot,NOUt,LERr,OK)
!
!        CTRTI2
!
         SRNamt = 'CTRTI2'
         INFot = 1
         CALL CTRTI2('/','N',0,a,1,info)
         CALL CHKXER('CTRTI2',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CTRTI2('U','/',0,a,1,info)
         CALL CHKXER('CTRTI2',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CTRTI2('U','N',-1,a,1,info)
         CALL CHKXER('CTRTI2',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL CTRTI2('U','N',2,a,1,info)
         CALL CHKXER('CTRTI2',INFot,NOUt,LERr,OK)
!
!
!        CTRTRS
!
         SRNamt = 'CTRTRS'
         INFot = 1
         CALL CTRTRS('/','N','N',0,0,a,1,x,1,info)
         CALL CHKXER('CTRTRS',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CTRTRS('U','/','N',0,0,a,1,x,1,info)
         CALL CHKXER('CTRTRS',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CTRTRS('U','N','/',0,0,a,1,x,1,info)
         CALL CHKXER('CTRTRS',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CTRTRS('U','N','N',-1,0,a,1,x,1,info)
         CALL CHKXER('CTRTRS',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL CTRTRS('U','N','N',0,-1,a,1,x,1,info)
         CALL CHKXER('CTRTRS',INFot,NOUt,LERr,OK)
         INFot = 7
!
!        CTRRFS
!
         SRNamt = 'CTRRFS'
         INFot = 1
         CALL CTRRFS('/','N','N',0,0,a,1,b,1,x,1,r1,r2,w,rw,info)
         CALL CHKXER('CTRRFS',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CTRRFS('U','/','N',0,0,a,1,b,1,x,1,r1,r2,w,rw,info)
         CALL CHKXER('CTRRFS',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CTRRFS('U','N','/',0,0,a,1,b,1,x,1,r1,r2,w,rw,info)
         CALL CHKXER('CTRRFS',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CTRRFS('U','N','N',-1,0,a,1,b,1,x,1,r1,r2,w,rw,info)
         CALL CHKXER('CTRRFS',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL CTRRFS('U','N','N',0,-1,a,1,b,1,x,1,r1,r2,w,rw,info)
         CALL CHKXER('CTRRFS',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL CTRRFS('U','N','N',2,1,a,1,b,2,x,2,r1,r2,w,rw,info)
         CALL CHKXER('CTRRFS',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL CTRRFS('U','N','N',2,1,a,2,b,1,x,2,r1,r2,w,rw,info)
         CALL CHKXER('CTRRFS',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL CTRRFS('U','N','N',2,1,a,2,b,2,x,1,r1,r2,w,rw,info)
         CALL CHKXER('CTRRFS',INFot,NOUt,LERr,OK)
!
!        CTRCON
!
         SRNamt = 'CTRCON'
         INFot = 1
         CALL CTRCON('/','U','N',0,a,1,rcond,w,rw,info)
         CALL CHKXER('CTRCON',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CTRCON('1','/','N',0,a,1,rcond,w,rw,info)
         CALL CHKXER('CTRCON',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CTRCON('1','U','/',0,a,1,rcond,w,rw,info)
         CALL CHKXER('CTRCON',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CTRCON('1','U','N',-1,a,1,rcond,w,rw,info)
         CALL CHKXER('CTRCON',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL CTRCON('1','U','N',2,a,1,rcond,w,rw,info)
         CALL CHKXER('CTRCON',INFot,NOUt,LERr,OK)
!
!        CLATRS
!
         SRNamt = 'CLATRS'
         INFot = 1
         CALL CLATRS('/','N','N','N',0,a,1,x,scale,rw,info)
         CALL CHKXER('CLATRS',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CLATRS('U','/','N','N',0,a,1,x,scale,rw,info)
         CALL CHKXER('CLATRS',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CLATRS('U','N','/','N',0,a,1,x,scale,rw,info)
         CALL CHKXER('CLATRS',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CLATRS('U','N','N','/',0,a,1,x,scale,rw,info)
         CALL CHKXER('CLATRS',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL CLATRS('U','N','N','N',-1,a,1,x,scale,rw,info)
         CALL CHKXER('CLATRS',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL CLATRS('U','N','N','N',2,a,1,x,scale,rw,info)
         CALL CHKXER('CLATRS',INFot,NOUt,LERr,OK)
!
!     Test error exits for the packed triangular routines.
!
      ELSEIF ( LSAMEN(2,c2,'TP') ) THEN
!
!        CTPTRI
!
         SRNamt = 'CTPTRI'
         INFot = 1
         CALL CTPTRI('/','N',0,a,info)
         CALL CHKXER('CTPTRI',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CTPTRI('U','/',0,a,info)
         CALL CHKXER('CTPTRI',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CTPTRI('U','N',-1,a,info)
         CALL CHKXER('CTPTRI',INFot,NOUt,LERr,OK)
!
!        CTPTRS
!
         SRNamt = 'CTPTRS'
         INFot = 1
         CALL CTPTRS('/','N','N',0,0,a,x,1,info)
         CALL CHKXER('CTPTRS',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CTPTRS('U','/','N',0,0,a,x,1,info)
         CALL CHKXER('CTPTRS',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CTPTRS('U','N','/',0,0,a,x,1,info)
         CALL CHKXER('CTPTRS',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CTPTRS('U','N','N',-1,0,a,x,1,info)
         CALL CHKXER('CTPTRS',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL CTPTRS('U','N','N',0,-1,a,x,1,info)
         CALL CHKXER('CTPTRS',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL CTPTRS('U','N','N',2,1,a,x,1,info)
         CALL CHKXER('CTPTRS',INFot,NOUt,LERr,OK)
!
!        CTPRFS
!
         SRNamt = 'CTPRFS'
         INFot = 1
         CALL CTPRFS('/','N','N',0,0,a,b,1,x,1,r1,r2,w,rw,info)
         CALL CHKXER('CTPRFS',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CTPRFS('U','/','N',0,0,a,b,1,x,1,r1,r2,w,rw,info)
         CALL CHKXER('CTPRFS',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CTPRFS('U','N','/',0,0,a,b,1,x,1,r1,r2,w,rw,info)
         CALL CHKXER('CTPRFS',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CTPRFS('U','N','N',-1,0,a,b,1,x,1,r1,r2,w,rw,info)
         CALL CHKXER('CTPRFS',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL CTPRFS('U','N','N',0,-1,a,b,1,x,1,r1,r2,w,rw,info)
         CALL CHKXER('CTPRFS',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL CTPRFS('U','N','N',2,1,a,b,1,x,2,r1,r2,w,rw,info)
         CALL CHKXER('CTPRFS',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL CTPRFS('U','N','N',2,1,a,b,2,x,1,r1,r2,w,rw,info)
         CALL CHKXER('CTPRFS',INFot,NOUt,LERr,OK)
!
!        CTPCON
!
         SRNamt = 'CTPCON'
         INFot = 1
         CALL CTPCON('/','U','N',0,a,rcond,w,rw,info)
         CALL CHKXER('CTPCON',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CTPCON('1','/','N',0,a,rcond,w,rw,info)
         CALL CHKXER('CTPCON',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CTPCON('1','U','/',0,a,rcond,w,rw,info)
         CALL CHKXER('CTPCON',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CTPCON('1','U','N',-1,a,rcond,w,rw,info)
         CALL CHKXER('CTPCON',INFot,NOUt,LERr,OK)
!
!        CLATPS
!
         SRNamt = 'CLATPS'
         INFot = 1
         CALL CLATPS('/','N','N','N',0,a,x,scale,rw,info)
         CALL CHKXER('CLATPS',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CLATPS('U','/','N','N',0,a,x,scale,rw,info)
         CALL CHKXER('CLATPS',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CLATPS('U','N','/','N',0,a,x,scale,rw,info)
         CALL CHKXER('CLATPS',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CLATPS('U','N','N','/',0,a,x,scale,rw,info)
         CALL CHKXER('CLATPS',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL CLATPS('U','N','N','N',-1,a,x,scale,rw,info)
         CALL CHKXER('CLATPS',INFot,NOUt,LERr,OK)
!
!     Test error exits for the banded triangular routines.
!
      ELSEIF ( LSAMEN(2,c2,'TB') ) THEN
!
!        CTBTRS
!
         SRNamt = 'CTBTRS'
         INFot = 1
         CALL CTBTRS('/','N','N',0,0,0,a,1,x,1,info)
         CALL CHKXER('CTBTRS',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CTBTRS('U','/','N',0,0,0,a,1,x,1,info)
         CALL CHKXER('CTBTRS',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CTBTRS('U','N','/',0,0,0,a,1,x,1,info)
         CALL CHKXER('CTBTRS',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CTBTRS('U','N','N',-1,0,0,a,1,x,1,info)
         CALL CHKXER('CTBTRS',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL CTBTRS('U','N','N',0,-1,0,a,1,x,1,info)
         CALL CHKXER('CTBTRS',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL CTBTRS('U','N','N',0,0,-1,a,1,x,1,info)
         CALL CHKXER('CTBTRS',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL CTBTRS('U','N','N',2,1,1,a,1,x,2,info)
         CALL CHKXER('CTBTRS',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL CTBTRS('U','N','N',2,0,1,a,1,x,1,info)
         CALL CHKXER('CTBTRS',INFot,NOUt,LERr,OK)
!
!        CTBRFS
!
         SRNamt = 'CTBRFS'
         INFot = 1
         CALL CTBRFS('/','N','N',0,0,0,a,1,b,1,x,1,r1,r2,w,rw,info)
         CALL CHKXER('CTBRFS',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CTBRFS('U','/','N',0,0,0,a,1,b,1,x,1,r1,r2,w,rw,info)
         CALL CHKXER('CTBRFS',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CTBRFS('U','N','/',0,0,0,a,1,b,1,x,1,r1,r2,w,rw,info)
         CALL CHKXER('CTBRFS',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CTBRFS('U','N','N',-1,0,0,a,1,b,1,x,1,r1,r2,w,rw,info)
         CALL CHKXER('CTBRFS',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL CTBRFS('U','N','N',0,-1,0,a,1,b,1,x,1,r1,r2,w,rw,info)
         CALL CHKXER('CTBRFS',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL CTBRFS('U','N','N',0,0,-1,a,1,b,1,x,1,r1,r2,w,rw,info)
         CALL CHKXER('CTBRFS',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL CTBRFS('U','N','N',2,1,1,a,1,b,2,x,2,r1,r2,w,rw,info)
         CALL CHKXER('CTBRFS',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL CTBRFS('U','N','N',2,1,1,a,2,b,1,x,2,r1,r2,w,rw,info)
         CALL CHKXER('CTBRFS',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL CTBRFS('U','N','N',2,1,1,a,2,b,2,x,1,r1,r2,w,rw,info)
         CALL CHKXER('CTBRFS',INFot,NOUt,LERr,OK)
!
!        CTBCON
!
         SRNamt = 'CTBCON'
         INFot = 1
         CALL CTBCON('/','U','N',0,0,a,1,rcond,w,rw,info)
         CALL CHKXER('CTBCON',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CTBCON('1','/','N',0,0,a,1,rcond,w,rw,info)
         CALL CHKXER('CTBCON',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CTBCON('1','U','/',0,0,a,1,rcond,w,rw,info)
         CALL CHKXER('CTBCON',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CTBCON('1','U','N',-1,0,a,1,rcond,w,rw,info)
         CALL CHKXER('CTBCON',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL CTBCON('1','U','N',0,-1,a,1,rcond,w,rw,info)
         CALL CHKXER('CTBCON',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL CTBCON('1','U','N',2,1,a,1,rcond,w,rw,info)
         CALL CHKXER('CTBCON',INFot,NOUt,LERr,OK)
!
!        CLATBS
!
         SRNamt = 'CLATBS'
         INFot = 1
         CALL CLATBS('/','N','N','N',0,0,a,1,x,scale,rw,info)
         CALL CHKXER('CLATBS',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CLATBS('U','/','N','N',0,0,a,1,x,scale,rw,info)
         CALL CHKXER('CLATBS',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CLATBS('U','N','/','N',0,0,a,1,x,scale,rw,info)
         CALL CHKXER('CLATBS',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CLATBS('U','N','N','/',0,0,a,1,x,scale,rw,info)
         CALL CHKXER('CLATBS',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL CLATBS('U','N','N','N',-1,0,a,1,x,scale,rw,info)
         CALL CHKXER('CLATBS',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL CLATBS('U','N','N','N',1,-1,a,1,x,scale,rw,info)
         CALL CHKXER('CLATBS',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL CLATBS('U','N','N','N',2,1,a,1,x,scale,rw,info)
         CALL CHKXER('CLATBS',INFot,NOUt,LERr,OK)
      ENDIF
!
!     Print a summary line.
!
      CALL ALAESM(Path,OK,NOUt)
!
!
!     End of CERRTR
!
      END SUBROUTINE CERRTR
