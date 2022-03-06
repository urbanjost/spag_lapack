!*==derrls.f90  processed by SPAG 7.51RB at 17:34 on  4 Mar 2022
!> \brief \b derrls
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE derrls( PATH, NUNIT )
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
!> derrls tests the error exits for the DOUBLE PRECISION least squares
!> driver routines (DGELS, SGELSS, SGELSY, SGELSD).
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
      SUBROUTINE derrls(Path,Nunit)
      IMPLICIT NONE
!*--DERRLS59
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
      DOUBLE PRECISION a(NMAX,NMAX) , b(NMAX,NMAX) , s(NMAX) , w(NMAX)
!     ..
!     .. External Functions ..
      LOGICAL LSAMEN
EXTERNAL    LSAMEN
!     ..
!     .. External Subroutines ..
EXTERNAL    CHKXER,   DGELS,   DGELSD,   DGELSS,   DGELSY
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
      a(1,1) = 1.0D+0
      a(1,2) = 2.0D+0
      a(2,2) = 3.0D+0
      a(2,1) = 4.0D+0
      OK = .TRUE.
!
      IF ( LSAMEN(2,c2,'LS') ) THEN
!
!        Test error exits for the least squares driver routines.
!
!        DGELS
!
         SRNamt = 'DGELS '
         INFot = 1
         CALL DGELS('/',0,0,0,a,1,b,1,w,1,info)
         CALL CHKXER('DGELS ',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DGELS('N',-1,0,0,a,1,b,1,w,1,info)
         CALL CHKXER('DGELS ',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DGELS('N',0,-1,0,a,1,b,1,w,1,info)
         CALL CHKXER('DGELS ',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL DGELS('N',0,0,-1,a,1,b,1,w,1,info)
         CALL CHKXER('DGELS ',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL DGELS('N',2,0,0,a,1,b,2,w,2,info)
         CALL CHKXER('DGELS ',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL DGELS('N',2,0,0,a,2,b,1,w,2,info)
         CALL CHKXER('DGELS ',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL DGELS('N',1,1,0,a,1,b,1,w,1,info)
         CALL CHKXER('DGELS ',INFot,NOUt,LERr,OK)
!
!        DGELSS
!
         SRNamt = 'DGELSS'
         INFot = 1
         CALL DGELSS(-1,0,0,a,1,b,1,s,rcond,irnk,w,1,info)
         CALL CHKXER('DGELSS',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DGELSS(0,-1,0,a,1,b,1,s,rcond,irnk,w,1,info)
         CALL CHKXER('DGELSS',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DGELSS(0,0,-1,a,1,b,1,s,rcond,irnk,w,1,info)
         CALL CHKXER('DGELSS',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL DGELSS(2,0,0,a,1,b,2,s,rcond,irnk,w,2,info)
         CALL CHKXER('DGELSS',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL DGELSS(2,0,0,a,2,b,1,s,rcond,irnk,w,2,info)
         CALL CHKXER('DGELSS',INFot,NOUt,LERr,OK)
!
!        DGELSY
!
         SRNamt = 'DGELSY'
         INFot = 1
         CALL DGELSY(-1,0,0,a,1,b,1,ip,rcond,irnk,w,10,info)
         CALL CHKXER('DGELSY',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DGELSY(0,-1,0,a,1,b,1,ip,rcond,irnk,w,10,info)
         CALL CHKXER('DGELSY',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DGELSY(0,0,-1,a,1,b,1,ip,rcond,irnk,w,10,info)
         CALL CHKXER('DGELSY',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL DGELSY(2,0,0,a,1,b,2,ip,rcond,irnk,w,10,info)
         CALL CHKXER('DGELSY',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL DGELSY(2,0,0,a,2,b,1,ip,rcond,irnk,w,10,info)
         CALL CHKXER('DGELSY',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL DGELSY(2,2,1,a,2,b,2,ip,rcond,irnk,w,1,info)
         CALL CHKXER('DGELSY',INFot,NOUt,LERr,OK)
!
!        DGELSD
!
         SRNamt = 'DGELSD'
         INFot = 1
         CALL DGELSD(-1,0,0,a,1,b,1,s,rcond,irnk,w,10,ip,info)
         CALL CHKXER('DGELSD',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DGELSD(0,-1,0,a,1,b,1,s,rcond,irnk,w,10,ip,info)
         CALL CHKXER('DGELSD',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DGELSD(0,0,-1,a,1,b,1,s,rcond,irnk,w,10,ip,info)
         CALL CHKXER('DGELSD',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL DGELSD(2,0,0,a,1,b,2,s,rcond,irnk,w,10,ip,info)
         CALL CHKXER('DGELSD',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL DGELSD(2,0,0,a,2,b,1,s,rcond,irnk,w,10,ip,info)
         CALL CHKXER('DGELSD',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL DGELSD(2,2,1,a,2,b,2,s,rcond,irnk,w,1,ip,info)
         CALL CHKXER('DGELSD',INFot,NOUt,LERr,OK)
      ENDIF
!
!     Print a summary line.
!
      CALL alaesm(Path,OK,NOUt)
!
!
!     End of derrls
!
      END SUBROUTINE derrls
