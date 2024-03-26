!*==derrgt.f90  processed by SPAG 7.51RB at 17:34 on  4 Mar 2022
!> \brief \b derrgt
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE derrgt( PATH, NUNIT )
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
!> derrgt tests the error exits for the DOUBLE PRECISION tridiagonal
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
      SUBROUTINE derrgt(Path,Nunit)
      IMPLICIT NONE
!*--DERRGT59
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
      DOUBLE PRECISION anorm , rcond
!     ..
!     .. Local Arrays ..
      INTEGER ip(NMAX) , iw(NMAX)
      DOUBLE PRECISION b(NMAX) , c(NMAX) , cf(NMAX) , d(NMAX) , df(NMAX)&
     &                 , e(NMAX) , ef(NMAX) , f(NMAX) , r1(NMAX) ,      &
     &                 r2(NMAX) , w(NMAX) , x(NMAX)
!     ..
!     .. External Functions ..
      LOGICAL LSAMEN
EXTERNAL    LSAMEN
!     ..
!     .. External Subroutines ..
EXTERNAL    DGTCON,   DGTRFS,   DGTTRF,   DGTTRS
EXTERNAL    DPTCON,   DPTRFS,   DPTTRF,   DPTTRS
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
      d(1) = 1.D0
      d(2) = 2.D0
      df(1) = 1.D0
      df(2) = 2.D0
      e(1) = 3.D0
      e(2) = 4.D0
      ef(1) = 3.D0
      ef(2) = 4.D0
      anorm = 1.0D0
      OK = .TRUE.
!
      IF ( LSAMEN(2,c2,'GT') ) THEN
!
!        Test error exits for the general tridiagonal routines.
!
!        DGTTRF
!
         SRNamt = 'DGTTRF'
         INFot = 1
         CALL DGTTRF(-1,c,d,e,f,ip,info)
         CALL CHKXER('DGTTRF',INFot,NOUt,LERr,OK)
!
!        DGTTRS
!
         SRNamt = 'DGTTRS'
         INFot = 1
         CALL DGTTRS('/',0,0,c,d,e,f,ip,x,1,info)
         CALL CHKXER('DGTTRS',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DGTTRS('N',-1,0,c,d,e,f,ip,x,1,info)
         CALL CHKXER('DGTTRS',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DGTTRS('N',0,-1,c,d,e,f,ip,x,1,info)
         CALL CHKXER('DGTTRS',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL DGTTRS('N',2,1,c,d,e,f,ip,x,1,info)
         CALL CHKXER('DGTTRS',INFot,NOUt,LERr,OK)
!
!        DGTRFS
!
         SRNamt = 'DGTRFS'
         INFot = 1
         CALL DGTRFS('/',0,0,c,d,e,cf,df,ef,f,ip,b,1,x,1,r1,r2,w,iw,    &
     &               info)
         CALL CHKXER('DGTRFS',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DGTRFS('N',-1,0,c,d,e,cf,df,ef,f,ip,b,1,x,1,r1,r2,w,iw,   &
     &               info)
         CALL CHKXER('DGTRFS',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DGTRFS('N',0,-1,c,d,e,cf,df,ef,f,ip,b,1,x,1,r1,r2,w,iw,   &
     &               info)
         CALL CHKXER('DGTRFS',INFot,NOUt,LERr,OK)
         INFot = 13
         CALL DGTRFS('N',2,1,c,d,e,cf,df,ef,f,ip,b,1,x,2,r1,r2,w,iw,    &
     &               info)
         CALL CHKXER('DGTRFS',INFot,NOUt,LERr,OK)
         INFot = 15
         CALL DGTRFS('N',2,1,c,d,e,cf,df,ef,f,ip,b,2,x,1,r1,r2,w,iw,    &
     &               info)
         CALL CHKXER('DGTRFS',INFot,NOUt,LERr,OK)
!
!        DGTCON
!
         SRNamt = 'DGTCON'
         INFot = 1
         CALL DGTCON('/',0,c,d,e,f,ip,anorm,rcond,w,iw,info)
         CALL CHKXER('DGTCON',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DGTCON('I',-1,c,d,e,f,ip,anorm,rcond,w,iw,info)
         CALL CHKXER('DGTCON',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL DGTCON('I',0,c,d,e,f,ip,-anorm,rcond,w,iw,info)
         CALL CHKXER('DGTCON',INFot,NOUt,LERr,OK)
!
      ELSEIF ( LSAMEN(2,c2,'PT') ) THEN
!
!        Test error exits for the positive definite tridiagonal
!        routines.
!
!        DPTTRF
!
         SRNamt = 'DPTTRF'
         INFot = 1
         CALL DPTTRF(-1,d,e,info)
         CALL CHKXER('DPTTRF',INFot,NOUt,LERr,OK)
!
!        DPTTRS
!
         SRNamt = 'DPTTRS'
         INFot = 1
         CALL DPTTRS(-1,0,d,e,x,1,info)
         CALL CHKXER('DPTTRS',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DPTTRS(0,-1,d,e,x,1,info)
         CALL CHKXER('DPTTRS',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL DPTTRS(2,1,d,e,x,1,info)
         CALL CHKXER('DPTTRS',INFot,NOUt,LERr,OK)
!
!        DPTRFS
!
         SRNamt = 'DPTRFS'
         INFot = 1
         CALL DPTRFS(-1,0,d,e,df,ef,b,1,x,1,r1,r2,w,info)
         CALL CHKXER('DPTRFS',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DPTRFS(0,-1,d,e,df,ef,b,1,x,1,r1,r2,w,info)
         CALL CHKXER('DPTRFS',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL DPTRFS(2,1,d,e,df,ef,b,1,x,2,r1,r2,w,info)
         CALL CHKXER('DPTRFS',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL DPTRFS(2,1,d,e,df,ef,b,2,x,1,r1,r2,w,info)
         CALL CHKXER('DPTRFS',INFot,NOUt,LERr,OK)
!
!        DPTCON
!
         SRNamt = 'DPTCON'
         INFot = 1
         CALL DPTCON(-1,d,e,anorm,rcond,w,info)
         CALL CHKXER('DPTCON',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL DPTCON(0,d,e,-anorm,rcond,w,info)
         CALL CHKXER('DPTCON',INFot,NOUt,LERr,OK)
      ENDIF
!
!     Print a summary line.
!
      CALL alaesm(Path,OK,NOUt)
!
!
!     End of derrgt
!
      END SUBROUTINE derrgt
