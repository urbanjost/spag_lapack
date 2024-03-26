!*==cerrgt.f90  processed by SPAG 7.51RB at 17:34 on  4 Mar 2022
!> \brief \b cerrgt
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE cerrgt( PATH, NUNIT )
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
!> cerrgt tests the error exits for the COMPLEX tridiagonal
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
!> \ingroup complex_lin
!
!  =====================================================================
      SUBROUTINE cerrgt(Path,Nunit)
      IMPLICIT NONE
!*--CERRGT59
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
      REAL anorm , rcond
!     ..
!     .. Local Arrays ..
      INTEGER ip(NMAX)
      REAL d(NMAX) , df(NMAX) , r1(NMAX) , r2(NMAX) , rw(NMAX)
      COMPLEX b(NMAX) , dl(NMAX) , dlf(NMAX) , du(NMAX) , du2(NMAX) ,   &
     &        duf(NMAX) , e(NMAX) , ef(NMAX) , w(NMAX) , x(NMAX)
!     ..
!     .. External Functions ..
      LOGICAL LSAMEN
EXTERNAL    LSAMEN
!     ..
!     .. External Subroutines ..
EXTERNAL    CGTCON,   CGTRFS,   CGTTRF,   CGTTRS   
EXTERNAL    CPTCON,   CPTRFS,   CPTTRF,   CPTTRS
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
         d(i) = 1.
         e(i) = 2.
         dl(i) = 3.
         du(i) = 4.
      ENDDO
      anorm = 1.0
      OK = .TRUE.
!
      IF ( LSAMEN(2,c2,'GT') ) THEN
!
!        Test error exits for the general tridiagonal routines.
!
!        CGTTRF
!
         SRNamt = 'CGTTRF'
         INFot = 1
         CALL CGTTRF(-1,dl,e,du,du2,ip,info)
         CALL CHKXER('CGTTRF',INFot,NOUt,LERr,OK)
!
!        CGTTRS
!
         SRNamt = 'CGTTRS'
         INFot = 1
         CALL CGTTRS('/',0,0,dl,e,du,du2,ip,x,1,info)
         CALL CHKXER('CGTTRS',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CGTTRS('N',-1,0,dl,e,du,du2,ip,x,1,info)
         CALL CHKXER('CGTTRS',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CGTTRS('N',0,-1,dl,e,du,du2,ip,x,1,info)
         CALL CHKXER('CGTTRS',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL CGTTRS('N',2,1,dl,e,du,du2,ip,x,1,info)
         CALL CHKXER('CGTTRS',INFot,NOUt,LERr,OK)
!
!        CGTRFS
!
         SRNamt = 'CGTRFS'
         INFot = 1
         CALL CGTRFS('/',0,0,dl,e,du,dlf,ef,duf,du2,ip,b,1,x,1,r1,r2,w, &
     &               rw,info)
         CALL CHKXER('CGTRFS',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CGTRFS('N',-1,0,dl,e,du,dlf,ef,duf,du2,ip,b,1,x,1,r1,r2,w,&
     &               rw,info)
         CALL CHKXER('CGTRFS',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CGTRFS('N',0,-1,dl,e,du,dlf,ef,duf,du2,ip,b,1,x,1,r1,r2,w,&
     &               rw,info)
         CALL CHKXER('CGTRFS',INFot,NOUt,LERr,OK)
         INFot = 13
         CALL CGTRFS('N',2,1,dl,e,du,dlf,ef,duf,du2,ip,b,1,x,2,r1,r2,w, &
     &               rw,info)
         CALL CHKXER('CGTRFS',INFot,NOUt,LERr,OK)
         INFot = 15
         CALL CGTRFS('N',2,1,dl,e,du,dlf,ef,duf,du2,ip,b,2,x,1,r1,r2,w, &
     &               rw,info)
         CALL CHKXER('CGTRFS',INFot,NOUt,LERr,OK)
!
!        CGTCON
!
         SRNamt = 'CGTCON'
         INFot = 1
         CALL CGTCON('/',0,dl,e,du,du2,ip,anorm,rcond,w,info)
         CALL CHKXER('CGTCON',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CGTCON('I',-1,dl,e,du,du2,ip,anorm,rcond,w,info)
         CALL CHKXER('CGTCON',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL CGTCON('I',0,dl,e,du,du2,ip,-anorm,rcond,w,info)
         CALL CHKXER('CGTCON',INFot,NOUt,LERr,OK)
!
      ELSEIF ( LSAMEN(2,c2,'PT') ) THEN
!
!        Test error exits for the positive definite tridiagonal
!        routines.
!
!        CPTTRF
!
         SRNamt = 'CPTTRF'
         INFot = 1
         CALL CPTTRF(-1,d,e,info)
         CALL CHKXER('CPTTRF',INFot,NOUt,LERr,OK)
!
!        CPTTRS
!
         SRNamt = 'CPTTRS'
         INFot = 1
         CALL CPTTRS('/',1,0,d,e,x,1,info)
         CALL CHKXER('CPTTRS',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CPTTRS('U',-1,0,d,e,x,1,info)
         CALL CHKXER('CPTTRS',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CPTTRS('U',0,-1,d,e,x,1,info)
         CALL CHKXER('CPTTRS',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL CPTTRS('U',2,1,d,e,x,1,info)
         CALL CHKXER('CPTTRS',INFot,NOUt,LERr,OK)
!
!        CPTRFS
!
         SRNamt = 'CPTRFS'
         INFot = 1
         CALL CPTRFS('/',1,0,d,e,df,ef,b,1,x,1,r1,r2,w,rw,info)
         CALL CHKXER('CPTRFS',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CPTRFS('U',-1,0,d,e,df,ef,b,1,x,1,r1,r2,w,rw,info)
         CALL CHKXER('CPTRFS',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CPTRFS('U',0,-1,d,e,df,ef,b,1,x,1,r1,r2,w,rw,info)
         CALL CHKXER('CPTRFS',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL CPTRFS('U',2,1,d,e,df,ef,b,1,x,2,r1,r2,w,rw,info)
         CALL CHKXER('CPTRFS',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL CPTRFS('U',2,1,d,e,df,ef,b,2,x,1,r1,r2,w,rw,info)
         CALL CHKXER('CPTRFS',INFot,NOUt,LERr,OK)
!
!        CPTCON
!
         SRNamt = 'CPTCON'
         INFot = 1
         CALL CPTCON(-1,d,e,anorm,rcond,rw,info)
         CALL CHKXER('CPTCON',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CPTCON(0,d,e,-anorm,rcond,rw,info)
         CALL CHKXER('CPTCON',INFot,NOUt,LERr,OK)
      ENDIF
!
!     Print a summary line.
!
      CALL alaesm(Path,OK,NOUt)
!
!
!     End of cerrgt
!
      END SUBROUTINE cerrgt
