!*==serrls.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b SERRLS
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE SERRLS( PATH, NUNIT )
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
!> SERRLS tests the error exits for the REAL least squares
!> driver routines (SGELS, SGELSS, SGELSY, SGELSD).
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
      SUBROUTINE SERRLS(Path,Nunit)
      IMPLICIT NONE
!*--SERRLS59
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
      REAL a(NMAX,NMAX) , b(NMAX,NMAX) , s(NMAX) , w(NMAX)
!     ..
!     .. External Functions ..
      LOGICAL LSAMEN
      EXTERNAL LSAMEN
!     ..
!     .. External Subroutines ..
      EXTERNAL ALAESM , CHKXER , SGELS , SGELSD , SGELSS , SGELSY
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
      a(1,1) = 1.0E+0
      a(1,2) = 2.0E+0
      a(2,2) = 3.0E+0
      a(2,1) = 4.0E+0
      OK = .TRUE.
!
      IF ( LSAMEN(2,c2,'LS') ) THEN
!
!        Test error exits for the least squares driver routines.
!
!        SGELS
!
         SRNamt = 'SGELS '
         INFot = 1
         CALL SGELS('/',0,0,0,a,1,b,1,w,1,info)
         CALL CHKXER('SGELS ',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SGELS('N',-1,0,0,a,1,b,1,w,1,info)
         CALL CHKXER('SGELS ',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SGELS('N',0,-1,0,a,1,b,1,w,1,info)
         CALL CHKXER('SGELS ',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL SGELS('N',0,0,-1,a,1,b,1,w,1,info)
         CALL CHKXER('SGELS ',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL SGELS('N',2,0,0,a,1,b,2,w,2,info)
         CALL CHKXER('SGELS ',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL SGELS('N',2,0,0,a,2,b,1,w,2,info)
         CALL CHKXER('SGELS ',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL SGELS('N',1,1,0,a,1,b,1,w,1,info)
         CALL CHKXER('SGELS ',INFot,NOUt,LERr,OK)
!
!        SGELSS
!
         SRNamt = 'SGELSS'
         INFot = 1
         CALL SGELSS(-1,0,0,a,1,b,1,s,rcond,irnk,w,1,info)
         CALL CHKXER('SGELSS',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SGELSS(0,-1,0,a,1,b,1,s,rcond,irnk,w,1,info)
         CALL CHKXER('SGELSS',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SGELSS(0,0,-1,a,1,b,1,s,rcond,irnk,w,1,info)
         CALL CHKXER('SGELSS',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL SGELSS(2,0,0,a,1,b,2,s,rcond,irnk,w,2,info)
         CALL CHKXER('SGELSS',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL SGELSS(2,0,0,a,2,b,1,s,rcond,irnk,w,2,info)
         CALL CHKXER('SGELSS',INFot,NOUt,LERr,OK)
!
!        SGELSY
!
         SRNamt = 'SGELSY'
         INFot = 1
         CALL SGELSY(-1,0,0,a,1,b,1,ip,rcond,irnk,w,10,info)
         CALL CHKXER('SGELSY',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SGELSY(0,-1,0,a,1,b,1,ip,rcond,irnk,w,10,info)
         CALL CHKXER('SGELSY',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SGELSY(0,0,-1,a,1,b,1,ip,rcond,irnk,w,10,info)
         CALL CHKXER('SGELSY',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL SGELSY(2,0,0,a,1,b,2,ip,rcond,irnk,w,10,info)
         CALL CHKXER('SGELSY',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL SGELSY(2,0,0,a,2,b,1,ip,rcond,irnk,w,10,info)
         CALL CHKXER('SGELSY',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL SGELSY(2,2,1,a,2,b,2,ip,rcond,irnk,w,1,info)
         CALL CHKXER('SGELSY',INFot,NOUt,LERr,OK)
!
!        SGELSD
!
         SRNamt = 'SGELSD'
         INFot = 1
         CALL SGELSD(-1,0,0,a,1,b,1,s,rcond,irnk,w,10,ip,info)
         CALL CHKXER('SGELSD',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SGELSD(0,-1,0,a,1,b,1,s,rcond,irnk,w,10,ip,info)
         CALL CHKXER('SGELSD',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SGELSD(0,0,-1,a,1,b,1,s,rcond,irnk,w,10,ip,info)
         CALL CHKXER('SGELSD',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL SGELSD(2,0,0,a,1,b,2,s,rcond,irnk,w,10,ip,info)
         CALL CHKXER('SGELSD',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL SGELSD(2,0,0,a,2,b,1,s,rcond,irnk,w,10,ip,info)
         CALL CHKXER('SGELSD',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL SGELSD(2,2,1,a,2,b,2,s,rcond,irnk,w,1,ip,info)
         CALL CHKXER('SGELSD',INFot,NOUt,LERr,OK)
      ENDIF
!
!     Print a summary line.
!
      CALL ALAESM(Path,OK,NOUt)
!
!
!     End of SERRLS
!
      END SUBROUTINE SERRLS
