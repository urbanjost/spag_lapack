!*==zerrtr.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b ZERRTR
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZERRTR( PATH, NUNIT )
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
!> ZERRTR tests the error exits for the COMPLEX*16 triangular routines.
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
!> \ingroup complex16_lin
!
!  =====================================================================
      SUBROUTINE ZERRTR(Path,Nunit)
      IMPLICIT NONE
!*--ZERRTR58
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
      DOUBLE PRECISION r1(NMAX) , r2(NMAX) , rw(NMAX)
      COMPLEX*16 a(NMAX,NMAX) , b(NMAX) , w(NMAX) , x(NMAX)
!     ..
!     .. External Functions ..
      LOGICAL LSAMEN
      EXTERNAL LSAMEN
!     ..
!     .. External Subroutines ..
      EXTERNAL ALAESM , CHKXER , ZLATBS , ZLATPS , ZLATRS , ZTBCON ,    &
     &         ZTBRFS , ZTBTRS , ZTPCON , ZTPRFS , ZTPTRI , ZTPTRS ,    &
     &         ZTRCON , ZTRRFS , ZTRTI2 , ZTRTRI , ZTRTRS
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
!     Test error exits for the general triangular routines.
!
      IF ( LSAMEN(2,c2,'TR') ) THEN
!
!        ZTRTRI
!
         SRNamt = 'ZTRTRI'
         INFot = 1
         CALL ZTRTRI('/','N',0,a,1,info)
         CALL CHKXER('ZTRTRI',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZTRTRI('U','/',0,a,1,info)
         CALL CHKXER('ZTRTRI',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZTRTRI('U','N',-1,a,1,info)
         CALL CHKXER('ZTRTRI',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL ZTRTRI('U','N',2,a,1,info)
         CALL CHKXER('ZTRTRI',INFot,NOUt,LERr,OK)
!
!        ZTRTI2
!
         SRNamt = 'ZTRTI2'
         INFot = 1
         CALL ZTRTI2('/','N',0,a,1,info)
         CALL CHKXER('ZTRTI2',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZTRTI2('U','/',0,a,1,info)
         CALL CHKXER('ZTRTI2',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZTRTI2('U','N',-1,a,1,info)
         CALL CHKXER('ZTRTI2',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL ZTRTI2('U','N',2,a,1,info)
         CALL CHKXER('ZTRTI2',INFot,NOUt,LERr,OK)
!
!
!        ZTRTRS
!
         SRNamt = 'ZTRTRS'
         INFot = 1
         CALL ZTRTRS('/','N','N',0,0,a,1,x,1,info)
         CALL CHKXER('ZTRTRS',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZTRTRS('U','/','N',0,0,a,1,x,1,info)
         CALL CHKXER('ZTRTRS',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZTRTRS('U','N','/',0,0,a,1,x,1,info)
         CALL CHKXER('ZTRTRS',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZTRTRS('U','N','N',-1,0,a,1,x,1,info)
         CALL CHKXER('ZTRTRS',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL ZTRTRS('U','N','N',0,-1,a,1,x,1,info)
         CALL CHKXER('ZTRTRS',INFot,NOUt,LERr,OK)
         INFot = 7
!
!        ZTRRFS
!
         SRNamt = 'ZTRRFS'
         INFot = 1
         CALL ZTRRFS('/','N','N',0,0,a,1,b,1,x,1,r1,r2,w,rw,info)
         CALL CHKXER('ZTRRFS',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZTRRFS('U','/','N',0,0,a,1,b,1,x,1,r1,r2,w,rw,info)
         CALL CHKXER('ZTRRFS',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZTRRFS('U','N','/',0,0,a,1,b,1,x,1,r1,r2,w,rw,info)
         CALL CHKXER('ZTRRFS',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZTRRFS('U','N','N',-1,0,a,1,b,1,x,1,r1,r2,w,rw,info)
         CALL CHKXER('ZTRRFS',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL ZTRRFS('U','N','N',0,-1,a,1,b,1,x,1,r1,r2,w,rw,info)
         CALL CHKXER('ZTRRFS',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL ZTRRFS('U','N','N',2,1,a,1,b,2,x,2,r1,r2,w,rw,info)
         CALL CHKXER('ZTRRFS',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL ZTRRFS('U','N','N',2,1,a,2,b,1,x,2,r1,r2,w,rw,info)
         CALL CHKXER('ZTRRFS',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL ZTRRFS('U','N','N',2,1,a,2,b,2,x,1,r1,r2,w,rw,info)
         CALL CHKXER('ZTRRFS',INFot,NOUt,LERr,OK)
!
!        ZTRCON
!
         SRNamt = 'ZTRCON'
         INFot = 1
         CALL ZTRCON('/','U','N',0,a,1,rcond,w,rw,info)
         CALL CHKXER('ZTRCON',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZTRCON('1','/','N',0,a,1,rcond,w,rw,info)
         CALL CHKXER('ZTRCON',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZTRCON('1','U','/',0,a,1,rcond,w,rw,info)
         CALL CHKXER('ZTRCON',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZTRCON('1','U','N',-1,a,1,rcond,w,rw,info)
         CALL CHKXER('ZTRCON',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL ZTRCON('1','U','N',2,a,1,rcond,w,rw,info)
         CALL CHKXER('ZTRCON',INFot,NOUt,LERr,OK)
!
!        ZLATRS
!
         SRNamt = 'ZLATRS'
         INFot = 1
         CALL ZLATRS('/','N','N','N',0,a,1,x,scale,rw,info)
         CALL CHKXER('ZLATRS',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZLATRS('U','/','N','N',0,a,1,x,scale,rw,info)
         CALL CHKXER('ZLATRS',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZLATRS('U','N','/','N',0,a,1,x,scale,rw,info)
         CALL CHKXER('ZLATRS',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZLATRS('U','N','N','/',0,a,1,x,scale,rw,info)
         CALL CHKXER('ZLATRS',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL ZLATRS('U','N','N','N',-1,a,1,x,scale,rw,info)
         CALL CHKXER('ZLATRS',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL ZLATRS('U','N','N','N',2,a,1,x,scale,rw,info)
         CALL CHKXER('ZLATRS',INFot,NOUt,LERr,OK)
!
!     Test error exits for the packed triangular routines.
!
      ELSEIF ( LSAMEN(2,c2,'TP') ) THEN
!
!        ZTPTRI
!
         SRNamt = 'ZTPTRI'
         INFot = 1
         CALL ZTPTRI('/','N',0,a,info)
         CALL CHKXER('ZTPTRI',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZTPTRI('U','/',0,a,info)
         CALL CHKXER('ZTPTRI',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZTPTRI('U','N',-1,a,info)
         CALL CHKXER('ZTPTRI',INFot,NOUt,LERr,OK)
!
!        ZTPTRS
!
         SRNamt = 'ZTPTRS'
         INFot = 1
         CALL ZTPTRS('/','N','N',0,0,a,x,1,info)
         CALL CHKXER('ZTPTRS',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZTPTRS('U','/','N',0,0,a,x,1,info)
         CALL CHKXER('ZTPTRS',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZTPTRS('U','N','/',0,0,a,x,1,info)
         CALL CHKXER('ZTPTRS',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZTPTRS('U','N','N',-1,0,a,x,1,info)
         CALL CHKXER('ZTPTRS',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL ZTPTRS('U','N','N',0,-1,a,x,1,info)
         CALL CHKXER('ZTPTRS',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL ZTPTRS('U','N','N',2,1,a,x,1,info)
         CALL CHKXER('ZTPTRS',INFot,NOUt,LERr,OK)
!
!        ZTPRFS
!
         SRNamt = 'ZTPRFS'
         INFot = 1
         CALL ZTPRFS('/','N','N',0,0,a,b,1,x,1,r1,r2,w,rw,info)
         CALL CHKXER('ZTPRFS',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZTPRFS('U','/','N',0,0,a,b,1,x,1,r1,r2,w,rw,info)
         CALL CHKXER('ZTPRFS',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZTPRFS('U','N','/',0,0,a,b,1,x,1,r1,r2,w,rw,info)
         CALL CHKXER('ZTPRFS',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZTPRFS('U','N','N',-1,0,a,b,1,x,1,r1,r2,w,rw,info)
         CALL CHKXER('ZTPRFS',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL ZTPRFS('U','N','N',0,-1,a,b,1,x,1,r1,r2,w,rw,info)
         CALL CHKXER('ZTPRFS',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL ZTPRFS('U','N','N',2,1,a,b,1,x,2,r1,r2,w,rw,info)
         CALL CHKXER('ZTPRFS',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL ZTPRFS('U','N','N',2,1,a,b,2,x,1,r1,r2,w,rw,info)
         CALL CHKXER('ZTPRFS',INFot,NOUt,LERr,OK)
!
!        ZTPCON
!
         SRNamt = 'ZTPCON'
         INFot = 1
         CALL ZTPCON('/','U','N',0,a,rcond,w,rw,info)
         CALL CHKXER('ZTPCON',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZTPCON('1','/','N',0,a,rcond,w,rw,info)
         CALL CHKXER('ZTPCON',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZTPCON('1','U','/',0,a,rcond,w,rw,info)
         CALL CHKXER('ZTPCON',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZTPCON('1','U','N',-1,a,rcond,w,rw,info)
         CALL CHKXER('ZTPCON',INFot,NOUt,LERr,OK)
!
!        ZLATPS
!
         SRNamt = 'ZLATPS'
         INFot = 1
         CALL ZLATPS('/','N','N','N',0,a,x,scale,rw,info)
         CALL CHKXER('ZLATPS',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZLATPS('U','/','N','N',0,a,x,scale,rw,info)
         CALL CHKXER('ZLATPS',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZLATPS('U','N','/','N',0,a,x,scale,rw,info)
         CALL CHKXER('ZLATPS',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZLATPS('U','N','N','/',0,a,x,scale,rw,info)
         CALL CHKXER('ZLATPS',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL ZLATPS('U','N','N','N',-1,a,x,scale,rw,info)
         CALL CHKXER('ZLATPS',INFot,NOUt,LERr,OK)
!
!     Test error exits for the banded triangular routines.
!
      ELSEIF ( LSAMEN(2,c2,'TB') ) THEN
!
!        ZTBTRS
!
         SRNamt = 'ZTBTRS'
         INFot = 1
         CALL ZTBTRS('/','N','N',0,0,0,a,1,x,1,info)
         CALL CHKXER('ZTBTRS',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZTBTRS('U','/','N',0,0,0,a,1,x,1,info)
         CALL CHKXER('ZTBTRS',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZTBTRS('U','N','/',0,0,0,a,1,x,1,info)
         CALL CHKXER('ZTBTRS',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZTBTRS('U','N','N',-1,0,0,a,1,x,1,info)
         CALL CHKXER('ZTBTRS',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL ZTBTRS('U','N','N',0,-1,0,a,1,x,1,info)
         CALL CHKXER('ZTBTRS',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL ZTBTRS('U','N','N',0,0,-1,a,1,x,1,info)
         CALL CHKXER('ZTBTRS',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL ZTBTRS('U','N','N',2,1,1,a,1,x,2,info)
         CALL CHKXER('ZTBTRS',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL ZTBTRS('U','N','N',2,0,1,a,1,x,1,info)
         CALL CHKXER('ZTBTRS',INFot,NOUt,LERr,OK)
!
!        ZTBRFS
!
         SRNamt = 'ZTBRFS'
         INFot = 1
         CALL ZTBRFS('/','N','N',0,0,0,a,1,b,1,x,1,r1,r2,w,rw,info)
         CALL CHKXER('ZTBRFS',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZTBRFS('U','/','N',0,0,0,a,1,b,1,x,1,r1,r2,w,rw,info)
         CALL CHKXER('ZTBRFS',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZTBRFS('U','N','/',0,0,0,a,1,b,1,x,1,r1,r2,w,rw,info)
         CALL CHKXER('ZTBRFS',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZTBRFS('U','N','N',-1,0,0,a,1,b,1,x,1,r1,r2,w,rw,info)
         CALL CHKXER('ZTBRFS',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL ZTBRFS('U','N','N',0,-1,0,a,1,b,1,x,1,r1,r2,w,rw,info)
         CALL CHKXER('ZTBRFS',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL ZTBRFS('U','N','N',0,0,-1,a,1,b,1,x,1,r1,r2,w,rw,info)
         CALL CHKXER('ZTBRFS',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL ZTBRFS('U','N','N',2,1,1,a,1,b,2,x,2,r1,r2,w,rw,info)
         CALL CHKXER('ZTBRFS',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL ZTBRFS('U','N','N',2,1,1,a,2,b,1,x,2,r1,r2,w,rw,info)
         CALL CHKXER('ZTBRFS',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL ZTBRFS('U','N','N',2,1,1,a,2,b,2,x,1,r1,r2,w,rw,info)
         CALL CHKXER('ZTBRFS',INFot,NOUt,LERr,OK)
!
!        ZTBCON
!
         SRNamt = 'ZTBCON'
         INFot = 1
         CALL ZTBCON('/','U','N',0,0,a,1,rcond,w,rw,info)
         CALL CHKXER('ZTBCON',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZTBCON('1','/','N',0,0,a,1,rcond,w,rw,info)
         CALL CHKXER('ZTBCON',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZTBCON('1','U','/',0,0,a,1,rcond,w,rw,info)
         CALL CHKXER('ZTBCON',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZTBCON('1','U','N',-1,0,a,1,rcond,w,rw,info)
         CALL CHKXER('ZTBCON',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL ZTBCON('1','U','N',0,-1,a,1,rcond,w,rw,info)
         CALL CHKXER('ZTBCON',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL ZTBCON('1','U','N',2,1,a,1,rcond,w,rw,info)
         CALL CHKXER('ZTBCON',INFot,NOUt,LERr,OK)
!
!        ZLATBS
!
         SRNamt = 'ZLATBS'
         INFot = 1
         CALL ZLATBS('/','N','N','N',0,0,a,1,x,scale,rw,info)
         CALL CHKXER('ZLATBS',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZLATBS('U','/','N','N',0,0,a,1,x,scale,rw,info)
         CALL CHKXER('ZLATBS',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZLATBS('U','N','/','N',0,0,a,1,x,scale,rw,info)
         CALL CHKXER('ZLATBS',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZLATBS('U','N','N','/',0,0,a,1,x,scale,rw,info)
         CALL CHKXER('ZLATBS',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL ZLATBS('U','N','N','N',-1,0,a,1,x,scale,rw,info)
         CALL CHKXER('ZLATBS',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL ZLATBS('U','N','N','N',1,-1,a,1,x,scale,rw,info)
         CALL CHKXER('ZLATBS',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL ZLATBS('U','N','N','N',2,1,a,1,x,scale,rw,info)
         CALL CHKXER('ZLATBS',INFot,NOUt,LERr,OK)
      ENDIF
!
!     Print a summary line.
!
      CALL ALAESM(Path,OK,NOUt)
!
!
!     End of ZERRTR
!
      END SUBROUTINE ZERRTR
