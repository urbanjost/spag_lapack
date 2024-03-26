!*==cerrls.f90  processed by SPAG 7.51RB at 17:34 on  4 Mar 2022
!> \brief \b cerrls
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE cerrls( PATH, NUNIT )
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
!> cerrls tests the error exits for the COMPLEX least squares
!> driver routines (CGELS, CGELSS, CGELSY, CGELSD).
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
      SUBROUTINE cerrls(Path,Nunit)
      IMPLICIT NONE
!*--CERRLS59
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
      REAL rcond
!     ..
!     .. Local Arrays ..
      INTEGER ip(NMAX)
      REAL rw(NMAX) , s(NMAX)
      COMPLEX a(NMAX,NMAX) , b(NMAX,NMAX) , w(NMAX)
!     ..
!     .. External Functions ..
      LOGICAL LSAMEN
EXTERNAL    LSAMEN
!     ..
!     .. External Subroutines ..
EXTERNAL    CGELS,   CGELSD,   CGELSS,   CGELSY
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
      a(1,1) = (1.0E+0,0.0E+0)
      a(1,2) = (2.0E+0,0.0E+0)
      a(2,2) = (3.0E+0,0.0E+0)
      a(2,1) = (4.0E+0,0.0E+0)
      OK = .TRUE.
      WRITE (NOUt,FMT=*)
!
!     Test error exits for the least squares driver routines.
!
      IF ( LSAMEN(2,c2,'LS') ) THEN
!
!        CGELS
!
         SRNamt = 'CGELS '
         INFot = 1
         CALL CGELS('/',0,0,0,a,1,b,1,w,1,info)
         CALL CHKXER('CGELS ',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CGELS('N',-1,0,0,a,1,b,1,w,1,info)
         CALL CHKXER('CGELS ',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CGELS('N',0,-1,0,a,1,b,1,w,1,info)
         CALL CHKXER('CGELS ',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CGELS('N',0,0,-1,a,1,b,1,w,1,info)
         CALL CHKXER('CGELS ',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL CGELS('N',2,0,0,a,1,b,2,w,2,info)
         CALL CHKXER('CGELS ',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL CGELS('N',2,0,0,a,2,b,1,w,2,info)
         CALL CHKXER('CGELS ',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL CGELS('N',1,1,0,a,1,b,1,w,1,info)
         CALL CHKXER('CGELS ',INFot,NOUt,LERr,OK)
!
!        CGELSS
!
         SRNamt = 'CGELSS'
         INFot = 1
         CALL CGELSS(-1,0,0,a,1,b,1,s,rcond,irnk,w,1,rw,info)
         CALL CHKXER('CGELSS',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CGELSS(0,-1,0,a,1,b,1,s,rcond,irnk,w,1,rw,info)
         CALL CHKXER('CGELSS',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CGELSS(0,0,-1,a,1,b,1,s,rcond,irnk,w,1,rw,info)
         CALL CHKXER('CGELSS',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL CGELSS(2,0,0,a,1,b,2,s,rcond,irnk,w,2,rw,info)
         CALL CHKXER('CGELSS',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL CGELSS(2,0,0,a,2,b,1,s,rcond,irnk,w,2,rw,info)
         CALL CHKXER('CGELSS',INFot,NOUt,LERr,OK)
!
!        CGELSY
!
         SRNamt = 'CGELSY'
         INFot = 1
         CALL CGELSY(-1,0,0,a,1,b,1,ip,rcond,irnk,w,10,rw,info)
         CALL CHKXER('CGELSY',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CGELSY(0,-1,0,a,1,b,1,ip,rcond,irnk,w,10,rw,info)
         CALL CHKXER('CGELSY',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CGELSY(0,0,-1,a,1,b,1,ip,rcond,irnk,w,10,rw,info)
         CALL CHKXER('CGELSY',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL CGELSY(2,0,0,a,1,b,2,ip,rcond,irnk,w,10,rw,info)
         CALL CHKXER('CGELSY',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL CGELSY(2,0,0,a,2,b,1,ip,rcond,irnk,w,10,rw,info)
         CALL CHKXER('CGELSY',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL CGELSY(0,3,0,a,1,b,3,ip,rcond,irnk,w,1,rw,info)
         CALL CHKXER('CGELSY',INFot,NOUt,LERr,OK)
!
!        CGELSD
!
         SRNamt = 'CGELSD'
         INFot = 1
         CALL CGELSD(-1,0,0,a,1,b,1,s,rcond,irnk,w,10,rw,ip,info)
         CALL CHKXER('CGELSD',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CGELSD(0,-1,0,a,1,b,1,s,rcond,irnk,w,10,rw,ip,info)
         CALL CHKXER('CGELSD',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CGELSD(0,0,-1,a,1,b,1,s,rcond,irnk,w,10,rw,ip,info)
         CALL CHKXER('CGELSD',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL CGELSD(2,0,0,a,1,b,2,s,rcond,irnk,w,10,rw,ip,info)
         CALL CHKXER('CGELSD',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL CGELSD(2,0,0,a,2,b,1,s,rcond,irnk,w,10,rw,ip,info)
         CALL CHKXER('CGELSD',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL CGELSD(2,2,1,a,2,b,2,s,rcond,irnk,w,1,rw,ip,info)
         CALL CHKXER('CGELSD',INFot,NOUt,LERr,OK)
      ENDIF
!
!     Print a summary line.
!
      CALL alaesm(Path,OK,NOUt)
!
!
!     End of cerrls
!
      END SUBROUTINE cerrls
