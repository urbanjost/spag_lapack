!*==serrgt.f90  processed by SPAG 7.51RB at 17:34 on  4 Mar 2022
!> \brief \b serrgt
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE serrgt( PATH, NUNIT )
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
!> serrgt tests the error exits for the REAL tridiagonal
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
!> \ingroup single_lin
!
!  =====================================================================
      SUBROUTINE serrgt(Path,Nunit)
      IMPLICIT NONE
!*--SERRGT59
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
      REAL anorm , rcond
!     ..
!     .. Local Arrays ..
      INTEGER ip(NMAX) , iw(NMAX)
      REAL b(NMAX) , c(NMAX) , cf(NMAX) , d(NMAX) , df(NMAX) , e(NMAX) ,&
     &     ef(NMAX) , f(NMAX) , r1(NMAX) , r2(NMAX) , w(NMAX) , x(NMAX)
!     ..
!     .. External Functions ..
      LOGICAL LSAMEN
EXTERNAL    LSAMEN
!     ..
!     .. External Subroutines ..
EXTERNAL    SGTCON,   SGTRFS,   SGTTRF,   SGTTRS
EXTERNAL    SPTCON,   SPTRFS,   SPTTRF,   SPTTRS
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
      d(1) = 1.
      d(2) = 2.
      df(1) = 1.
      df(2) = 2.
      e(1) = 3.
      e(2) = 4.
      ef(1) = 3.
      ef(2) = 4.
      anorm = 1.0
      OK = .TRUE.
!
      IF ( LSAMEN(2,c2,'GT') ) THEN
!
!        Test error exits for the general tridiagonal routines.
!
!        SGTTRF
!
         SRNamt = 'SGTTRF'
         INFot = 1
         CALL SGTTRF(-1,c,d,e,f,ip,info)
         CALL CHKXER('SGTTRF',INFot,NOUt,LERr,OK)
!
!        SGTTRS
!
         SRNamt = 'SGTTRS'
         INFot = 1
         CALL SGTTRS('/',0,0,c,d,e,f,ip,x,1,info)
         CALL CHKXER('SGTTRS',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SGTTRS('N',-1,0,c,d,e,f,ip,x,1,info)
         CALL CHKXER('SGTTRS',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SGTTRS('N',0,-1,c,d,e,f,ip,x,1,info)
         CALL CHKXER('SGTTRS',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL SGTTRS('N',2,1,c,d,e,f,ip,x,1,info)
         CALL CHKXER('SGTTRS',INFot,NOUt,LERr,OK)
!
!        SGTRFS
!
         SRNamt = 'SGTRFS'
         INFot = 1
         CALL SGTRFS('/',0,0,c,d,e,cf,df,ef,f,ip,b,1,x,1,r1,r2,w,iw,    &
     &               info)
         CALL CHKXER('SGTRFS',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SGTRFS('N',-1,0,c,d,e,cf,df,ef,f,ip,b,1,x,1,r1,r2,w,iw,   &
     &               info)
         CALL CHKXER('SGTRFS',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SGTRFS('N',0,-1,c,d,e,cf,df,ef,f,ip,b,1,x,1,r1,r2,w,iw,   &
     &               info)
         CALL CHKXER('SGTRFS',INFot,NOUt,LERr,OK)
         INFot = 13
         CALL SGTRFS('N',2,1,c,d,e,cf,df,ef,f,ip,b,1,x,2,r1,r2,w,iw,    &
     &               info)
         CALL CHKXER('SGTRFS',INFot,NOUt,LERr,OK)
         INFot = 15
         CALL SGTRFS('N',2,1,c,d,e,cf,df,ef,f,ip,b,2,x,1,r1,r2,w,iw,    &
     &               info)
         CALL CHKXER('SGTRFS',INFot,NOUt,LERr,OK)
!
!        SGTCON
!
         SRNamt = 'SGTCON'
         INFot = 1
         CALL SGTCON('/',0,c,d,e,f,ip,anorm,rcond,w,iw,info)
         CALL CHKXER('SGTCON',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SGTCON('I',-1,c,d,e,f,ip,anorm,rcond,w,iw,info)
         CALL CHKXER('SGTCON',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL SGTCON('I',0,c,d,e,f,ip,-anorm,rcond,w,iw,info)
         CALL CHKXER('SGTCON',INFot,NOUt,LERr,OK)
!
      ELSEIF ( LSAMEN(2,c2,'PT') ) THEN
!
!        Test error exits for the positive definite tridiagonal
!        routines.
!
!        SPTTRF
!
         SRNamt = 'SPTTRF'
         INFot = 1
         CALL SPTTRF(-1,d,e,info)
         CALL CHKXER('SPTTRF',INFot,NOUt,LERr,OK)
!
!        SPTTRS
!
         SRNamt = 'SPTTRS'
         INFot = 1
         CALL SPTTRS(-1,0,d,e,x,1,info)
         CALL CHKXER('SPTTRS',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SPTTRS(0,-1,d,e,x,1,info)
         CALL CHKXER('SPTTRS',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL SPTTRS(2,1,d,e,x,1,info)
         CALL CHKXER('SPTTRS',INFot,NOUt,LERr,OK)
!
!        SPTRFS
!
         SRNamt = 'SPTRFS'
         INFot = 1
         CALL SPTRFS(-1,0,d,e,df,ef,b,1,x,1,r1,r2,w,info)
         CALL CHKXER('SPTRFS',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SPTRFS(0,-1,d,e,df,ef,b,1,x,1,r1,r2,w,info)
         CALL CHKXER('SPTRFS',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL SPTRFS(2,1,d,e,df,ef,b,1,x,2,r1,r2,w,info)
         CALL CHKXER('SPTRFS',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL SPTRFS(2,1,d,e,df,ef,b,2,x,1,r1,r2,w,info)
         CALL CHKXER('SPTRFS',INFot,NOUt,LERr,OK)
!
!        SPTCON
!
         SRNamt = 'SPTCON'
         INFot = 1
         CALL SPTCON(-1,d,e,anorm,rcond,w,info)
         CALL CHKXER('SPTCON',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL SPTCON(0,d,e,-anorm,rcond,w,info)
         CALL CHKXER('SPTCON',INFot,NOUt,LERr,OK)
      ENDIF
!
!     Print a summary line.
!
      CALL alaesm(Path,OK,NOUt)
!
!
!     End of serrgt
!
      END SUBROUTINE serrgt
