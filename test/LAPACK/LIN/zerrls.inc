!*==zerrls.f90  processed by SPAG 7.51RB at 17:35 on  4 Mar 2022
!> \brief \b zerrls
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE zerrls( PATH, NUNIT )
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
!> zerrls tests the error exits for the COMPLEX*16 least squares
!> driver routines (ZGELS, CGELSS, CGELSY, CGELSD).
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
      SUBROUTINE zerrls(Path,Nunit)
      IMPLICIT NONE
!*--ZERRLS59
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
      INTEGER info , irnk
      DOUBLE PRECISION rcond
!     ..
!     .. Local Arrays ..
      INTEGER ip(NMAX)
      DOUBLE PRECISION rw(NMAX) , s(NMAX)
      COMPLEX*16 a(NMAX,NMAX) , b(NMAX,NMAX) , w(NMAX)
!     ..
!     .. External Functions ..
      LOGICAL LSAMEN
EXTERNAL    LSAMEN
!     ..
!     .. External Subroutines ..
EXTERNAL    CHKXER,   ZGELS,   ZGELSD,   ZGELSS,   ZGELSY
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
      c2 = Path(2:3)
      a(1,1) = (1.0D+0,0.0D+0)
      a(1,2) = (2.0D+0,0.0D+0)
      a(2,2) = (3.0D+0,0.0D+0)
      a(2,1) = (4.0D+0,0.0D+0)
      OK = .TRUE.
      WRITE (NOUt,FMT=*)
!
!     Test error exits for the least squares driver routines.
!
      IF ( LSAMEN(2,c2,'LS') ) THEN
!
!        ZGELS
!
         SRNamt = 'ZGELS '
         INFot = 1
         CALL ZGELS('/',0,0,0,a,1,b,1,w,1,info)
         CALL CHKXER('ZGELS ',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZGELS('N',-1,0,0,a,1,b,1,w,1,info)
         CALL CHKXER('ZGELS ',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZGELS('N',0,-1,0,a,1,b,1,w,1,info)
         CALL CHKXER('ZGELS ',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZGELS('N',0,0,-1,a,1,b,1,w,1,info)
         CALL CHKXER('ZGELS ',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL ZGELS('N',2,0,0,a,1,b,2,w,2,info)
         CALL CHKXER('ZGELS ',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL ZGELS('N',2,0,0,a,2,b,1,w,2,info)
         CALL CHKXER('ZGELS ',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL ZGELS('N',1,1,0,a,1,b,1,w,1,info)
         CALL CHKXER('ZGELS ',INFot,NOUt,LERr,OK)
!
!        ZGELSS
!
         SRNamt = 'ZGELSS'
         INFot = 1
         CALL ZGELSS(-1,0,0,a,1,b,1,s,rcond,irnk,w,1,rw,info)
         CALL CHKXER('ZGELSS',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZGELSS(0,-1,0,a,1,b,1,s,rcond,irnk,w,1,rw,info)
         CALL CHKXER('ZGELSS',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZGELSS(0,0,-1,a,1,b,1,s,rcond,irnk,w,1,rw,info)
         CALL CHKXER('ZGELSS',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL ZGELSS(2,0,0,a,1,b,2,s,rcond,irnk,w,2,rw,info)
         CALL CHKXER('ZGELSS',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL ZGELSS(2,0,0,a,2,b,1,s,rcond,irnk,w,2,rw,info)
         CALL CHKXER('ZGELSS',INFot,NOUt,LERr,OK)
!
!        ZGELSY
!
         SRNamt = 'ZGELSY'
         INFot = 1
         CALL ZGELSY(-1,0,0,a,1,b,1,ip,rcond,irnk,w,10,rw,info)
         CALL CHKXER('ZGELSY',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZGELSY(0,-1,0,a,1,b,1,ip,rcond,irnk,w,10,rw,info)
         CALL CHKXER('ZGELSY',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZGELSY(0,0,-1,a,1,b,1,ip,rcond,irnk,w,10,rw,info)
         CALL CHKXER('ZGELSY',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL ZGELSY(2,0,0,a,1,b,2,ip,rcond,irnk,w,10,rw,info)
         CALL CHKXER('ZGELSY',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL ZGELSY(2,0,0,a,2,b,1,ip,rcond,irnk,w,10,rw,info)
         CALL CHKXER('ZGELSY',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL ZGELSY(0,3,0,a,1,b,3,ip,rcond,irnk,w,1,rw,info)
         CALL CHKXER('ZGELSY',INFot,NOUt,LERr,OK)
!
!        ZGELSD
!
         SRNamt = 'ZGELSD'
         INFot = 1
         CALL ZGELSD(-1,0,0,a,1,b,1,s,rcond,irnk,w,10,rw,ip,info)
         CALL CHKXER('ZGELSD',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZGELSD(0,-1,0,a,1,b,1,s,rcond,irnk,w,10,rw,ip,info)
         CALL CHKXER('ZGELSD',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZGELSD(0,0,-1,a,1,b,1,s,rcond,irnk,w,10,rw,ip,info)
         CALL CHKXER('ZGELSD',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL ZGELSD(2,0,0,a,1,b,2,s,rcond,irnk,w,10,rw,ip,info)
         CALL CHKXER('ZGELSD',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL ZGELSD(2,0,0,a,2,b,1,s,rcond,irnk,w,10,rw,ip,info)
         CALL CHKXER('ZGELSD',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL ZGELSD(2,2,1,a,2,b,2,s,rcond,irnk,w,1,rw,ip,info)
         CALL CHKXER('ZGELSD',INFot,NOUt,LERr,OK)
      ENDIF
!
!     Print a summary line.
!
      CALL alaesm(Path,OK,NOUt)
!
!
!     End of zerrls
!
      END SUBROUTINE zerrls
