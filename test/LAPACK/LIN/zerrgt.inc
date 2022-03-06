!*==zerrgt.f90  processed by SPAG 7.51RB at 17:35 on  4 Mar 2022
!> \brief \b zerrgt
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE zerrgt( PATH, NUNIT )
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
!> zerrgt tests the error exits for the COMPLEX*16 tridiagonal
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
!> \ingroup complex16_lin
!
!  =====================================================================
      SUBROUTINE zerrgt(Path,Nunit)
      IMPLICIT NONE
!*--ZERRGT59
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
      INTEGER i , info
      DOUBLE PRECISION anorm , rcond
!     ..
!     .. Local Arrays ..
      INTEGER ip(NMAX)
      DOUBLE PRECISION d(NMAX) , df(NMAX) , r1(NMAX) , r2(NMAX) ,       &
     &                 rw(NMAX)
      COMPLEX*16 b(NMAX) , dl(NMAX) , dlf(NMAX) , du(NMAX) , du2(NMAX) ,&
     &           duf(NMAX) , e(NMAX) , ef(NMAX) , w(NMAX) , x(NMAX)
!     ..
!     .. External Functions ..
      LOGICAL LSAMEN
EXTERNAL    LSAMEN
!     ..
!     .. External Subroutines ..
EXTERNAL    CHKXER,   ZGTCON,   ZGTRFS,   ZGTTRF,   ZGTTRS
EXTERNAL    ZPTCON,   ZPTRFS,   ZPTTRF,   ZPTTRS
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
      DO i = 1 , NMAX
         d(i) = 1.D0
         e(i) = 2.D0
         dl(i) = 3.D0
         du(i) = 4.D0
      ENDDO
      anorm = 1.0D0
      OK = .TRUE.
!
      IF ( LSAMEN(2,c2,'GT') ) THEN
!
!        Test error exits for the general tridiagonal routines.
!
!        ZGTTRF
!
         SRNamt = 'ZGTTRF'
         INFot = 1
         CALL ZGTTRF(-1,dl,e,du,du2,ip,info)
         CALL CHKXER('ZGTTRF',INFot,NOUt,LERr,OK)
!
!        ZGTTRS
!
         SRNamt = 'ZGTTRS'
         INFot = 1
         CALL ZGTTRS('/',0,0,dl,e,du,du2,ip,x,1,info)
         CALL CHKXER('ZGTTRS',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZGTTRS('N',-1,0,dl,e,du,du2,ip,x,1,info)
         CALL CHKXER('ZGTTRS',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZGTTRS('N',0,-1,dl,e,du,du2,ip,x,1,info)
         CALL CHKXER('ZGTTRS',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL ZGTTRS('N',2,1,dl,e,du,du2,ip,x,1,info)
         CALL CHKXER('ZGTTRS',INFot,NOUt,LERr,OK)
!
!        ZGTRFS
!
         SRNamt = 'ZGTRFS'
         INFot = 1
         CALL ZGTRFS('/',0,0,dl,e,du,dlf,ef,duf,du2,ip,b,1,x,1,r1,r2,w, &
     &               rw,info)
         CALL CHKXER('ZGTRFS',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZGTRFS('N',-1,0,dl,e,du,dlf,ef,duf,du2,ip,b,1,x,1,r1,r2,w,&
     &               rw,info)
         CALL CHKXER('ZGTRFS',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZGTRFS('N',0,-1,dl,e,du,dlf,ef,duf,du2,ip,b,1,x,1,r1,r2,w,&
     &               rw,info)
         CALL CHKXER('ZGTRFS',INFot,NOUt,LERr,OK)
         INFot = 13
         CALL ZGTRFS('N',2,1,dl,e,du,dlf,ef,duf,du2,ip,b,1,x,2,r1,r2,w, &
     &               rw,info)
         CALL CHKXER('ZGTRFS',INFot,NOUt,LERr,OK)
         INFot = 15
         CALL ZGTRFS('N',2,1,dl,e,du,dlf,ef,duf,du2,ip,b,2,x,1,r1,r2,w, &
     &               rw,info)
         CALL CHKXER('ZGTRFS',INFot,NOUt,LERr,OK)
!
!        ZGTCON
!
         SRNamt = 'ZGTCON'
         INFot = 1
         CALL ZGTCON('/',0,dl,e,du,du2,ip,anorm,rcond,w,info)
         CALL CHKXER('ZGTCON',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZGTCON('I',-1,dl,e,du,du2,ip,anorm,rcond,w,info)
         CALL CHKXER('ZGTCON',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL ZGTCON('I',0,dl,e,du,du2,ip,-anorm,rcond,w,info)
         CALL CHKXER('ZGTCON',INFot,NOUt,LERr,OK)
!
      ELSEIF ( LSAMEN(2,c2,'PT') ) THEN
!
!        Test error exits for the positive definite tridiagonal
!        routines.
!
!        ZPTTRF
!
         SRNamt = 'ZPTTRF'
         INFot = 1
         CALL ZPTTRF(-1,d,e,info)
         CALL CHKXER('ZPTTRF',INFot,NOUt,LERr,OK)
!
!        ZPTTRS
!
         SRNamt = 'ZPTTRS'
         INFot = 1
         CALL ZPTTRS('/',1,0,d,e,x,1,info)
         CALL CHKXER('ZPTTRS',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZPTTRS('U',-1,0,d,e,x,1,info)
         CALL CHKXER('ZPTTRS',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZPTTRS('U',0,-1,d,e,x,1,info)
         CALL CHKXER('ZPTTRS',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL ZPTTRS('U',2,1,d,e,x,1,info)
         CALL CHKXER('ZPTTRS',INFot,NOUt,LERr,OK)
!
!        ZPTRFS
!
         SRNamt = 'ZPTRFS'
         INFot = 1
         CALL ZPTRFS('/',1,0,d,e,df,ef,b,1,x,1,r1,r2,w,rw,info)
         CALL CHKXER('ZPTRFS',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZPTRFS('U',-1,0,d,e,df,ef,b,1,x,1,r1,r2,w,rw,info)
         CALL CHKXER('ZPTRFS',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZPTRFS('U',0,-1,d,e,df,ef,b,1,x,1,r1,r2,w,rw,info)
         CALL CHKXER('ZPTRFS',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL ZPTRFS('U',2,1,d,e,df,ef,b,1,x,2,r1,r2,w,rw,info)
         CALL CHKXER('ZPTRFS',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL ZPTRFS('U',2,1,d,e,df,ef,b,2,x,1,r1,r2,w,rw,info)
         CALL CHKXER('ZPTRFS',INFot,NOUt,LERr,OK)
!
!        ZPTCON
!
         SRNamt = 'ZPTCON'
         INFot = 1
         CALL ZPTCON(-1,d,e,anorm,rcond,rw,info)
         CALL CHKXER('ZPTCON',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZPTCON(0,d,e,-anorm,rcond,rw,info)
         CALL CHKXER('ZPTCON',INFot,NOUt,LERr,OK)
      ENDIF
!
!     Print a summary line.
!
      CALL alaesm(Path,OK,NOUt)
!
!
!     End of zerrgt
!
      END SUBROUTINE zerrgt
