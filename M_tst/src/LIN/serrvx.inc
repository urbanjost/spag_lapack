!*==serrvx.f90  processed by SPAG 7.51RB at 17:34 on  4 Mar 2022
!> \brief \b serrvx
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE serrvx( PATH, NUNIT )
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
!> serrvx tests the error exits for the REAL driver routines
!> for solving linear systems of equations.
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
!> \date November 2017
!
!> \ingroup single_lin
!
!  =====================================================================
      SUBROUTINE serrvx(Path,Nunit)
      IMPLICIT NONE
!*--SERRVX59
!
!  -- LAPACK test routine (version 3.8.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2017
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
      PARAMETER (NMAX=4)
!     ..
!     .. Local Scalars ..
      CHARACTER eq
      CHARACTER*2 c2
      INTEGER i , info , j
      REAL rcond
!     ..
!     .. Local Arrays ..
      INTEGER ip(NMAX) , iw(NMAX)
      REAL a(NMAX,NMAX) , af(NMAX,NMAX) , b(NMAX) , c(NMAX) , e(NMAX) , &
     &     r(NMAX) , r1(NMAX) , r2(NMAX) , w(2*NMAX) , x(NMAX)
!     ..
!     .. External Functions ..
      LOGICAL LSAMEN
EXTERNAL    LSAMEN
!     ..
!     .. External Subroutines ..
EXTERNAL    SGBSV,   SGBSVX,   SGESV,   SGESVX,   SGTSV
EXTERNAL    SGTSVX,   SPBSV,   SPBSVX,   SPOSV,   SPOSVX,   SPPSV
EXTERNAL    SPPSVX,   SPTSV,   SPTSVX,   SSPSV,   SSPSVX,   SSYSV
EXTERNAL    SSYSV_AA,   SSYSV_RK,   SSYSV_ROOK,   SSYSVX
EXTERNAL    SSYSV_AA_2STAGE
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
!     .. Intrinsic Functions ..
      INTRINSIC REAL
!     ..
!     .. Executable Statements ..
!
      NOUt = Nunit
      WRITE (NOUt,FMT=*)
      c2 = Path(2:3)
!
!     Set the variables to innocuous values.
!
      DO j = 1 , NMAX
         DO i = 1 , NMAX
            a(i,j) = 1./REAL(i+j)
            af(i,j) = 1./REAL(i+j)
         ENDDO
         b(j) = 0.E+0
         e(j) = 0.E+0
         r1(j) = 0.E+0
         r2(j) = 0.E+0
         w(j) = 0.E+0
         x(j) = 0.E+0
         c(j) = 0.E+0
         r(j) = 0.E+0
         ip(j) = j
      ENDDO
      eq = ' '
      OK = .TRUE.
!
      IF ( LSAMEN(2,c2,'GE') ) THEN
!
!        SGESV
!
         SRNamt = 'SGESV '
         INFot = 1
         CALL SGESV(-1,0,a,1,ip,b,1,info)
         CALL CHKXER('SGESV ',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SGESV(0,-1,a,1,ip,b,1,info)
         CALL CHKXER('SGESV ',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL SGESV(2,1,a,1,ip,b,2,info)
         CALL CHKXER('SGESV ',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL SGESV(2,1,a,2,ip,b,1,info)
         CALL CHKXER('SGESV ',INFot,NOUt,LERr,OK)
!
!        SGESVX
!
         SRNamt = 'SGESVX'
         INFot = 1
         CALL SGESVX('/','N',0,0,a,1,af,1,ip,eq,r,c,b,1,x,1,rcond,r1,r2,&
     &               w,iw,info)
         CALL CHKXER('SGESVX',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SGESVX('N','/',0,0,a,1,af,1,ip,eq,r,c,b,1,x,1,rcond,r1,r2,&
     &               w,iw,info)
         CALL CHKXER('SGESVX',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SGESVX('N','N',-1,0,a,1,af,1,ip,eq,r,c,b,1,x,1,rcond,r1,  &
     &               r2,w,iw,info)
         CALL CHKXER('SGESVX',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL SGESVX('N','N',0,-1,a,1,af,1,ip,eq,r,c,b,1,x,1,rcond,r1,  &
     &               r2,w,iw,info)
         CALL CHKXER('SGESVX',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL SGESVX('N','N',2,1,a,1,af,2,ip,eq,r,c,b,2,x,2,rcond,r1,r2,&
     &               w,iw,info)
         CALL CHKXER('SGESVX',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL SGESVX('N','N',2,1,a,2,af,1,ip,eq,r,c,b,2,x,2,rcond,r1,r2,&
     &               w,iw,info)
         CALL CHKXER('SGESVX',INFot,NOUt,LERr,OK)
         INFot = 10
         eq = '/'
         CALL SGESVX('F','N',0,0,a,1,af,1,ip,eq,r,c,b,1,x,1,rcond,r1,r2,&
     &               w,iw,info)
         CALL CHKXER('SGESVX',INFot,NOUt,LERr,OK)
         INFot = 11
         eq = 'R'
         CALL SGESVX('F','N',1,0,a,1,af,1,ip,eq,r,c,b,1,x,1,rcond,r1,r2,&
     &               w,iw,info)
         CALL CHKXER('SGESVX',INFot,NOUt,LERr,OK)
         INFot = 12
         eq = 'C'
         CALL SGESVX('F','N',1,0,a,1,af,1,ip,eq,r,c,b,1,x,1,rcond,r1,r2,&
     &               w,iw,info)
         CALL CHKXER('SGESVX',INFot,NOUt,LERr,OK)
         INFot = 14
         CALL SGESVX('N','N',2,1,a,2,af,2,ip,eq,r,c,b,1,x,2,rcond,r1,r2,&
     &               w,iw,info)
         CALL CHKXER('SGESVX',INFot,NOUt,LERr,OK)
         INFot = 16
         CALL SGESVX('N','N',2,1,a,2,af,2,ip,eq,r,c,b,2,x,1,rcond,r1,r2,&
     &               w,iw,info)
         CALL CHKXER('SGESVX',INFot,NOUt,LERr,OK)
!
      ELSEIF ( LSAMEN(2,c2,'GB') ) THEN
!
!        SGBSV
!
         SRNamt = 'SGBSV '
         INFot = 1
         CALL SGBSV(-1,0,0,0,a,1,ip,b,1,info)
         CALL CHKXER('SGBSV ',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SGBSV(1,-1,0,0,a,1,ip,b,1,info)
         CALL CHKXER('SGBSV ',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SGBSV(1,0,-1,0,a,1,ip,b,1,info)
         CALL CHKXER('SGBSV ',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL SGBSV(0,0,0,-1,a,1,ip,b,1,info)
         CALL CHKXER('SGBSV ',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL SGBSV(1,1,1,0,a,3,ip,b,1,info)
         CALL CHKXER('SGBSV ',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL SGBSV(2,0,0,0,a,1,ip,b,1,info)
         CALL CHKXER('SGBSV ',INFot,NOUt,LERr,OK)
!
!        SGBSVX
!
         SRNamt = 'SGBSVX'
         INFot = 1
         CALL SGBSVX('/','N',0,0,0,0,a,1,af,1,ip,eq,r,c,b,1,x,1,rcond,  &
     &               r1,r2,w,iw,info)
         CALL CHKXER('SGBSVX',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SGBSVX('N','/',0,0,0,0,a,1,af,1,ip,eq,r,c,b,1,x,1,rcond,  &
     &               r1,r2,w,iw,info)
         CALL CHKXER('SGBSVX',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SGBSVX('N','N',-1,0,0,0,a,1,af,1,ip,eq,r,c,b,1,x,1,rcond, &
     &               r1,r2,w,iw,info)
         CALL CHKXER('SGBSVX',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL SGBSVX('N','N',1,-1,0,0,a,1,af,1,ip,eq,r,c,b,1,x,1,rcond, &
     &               r1,r2,w,iw,info)
         CALL CHKXER('SGBSVX',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL SGBSVX('N','N',1,0,-1,0,a,1,af,1,ip,eq,r,c,b,1,x,1,rcond, &
     &               r1,r2,w,iw,info)
         CALL CHKXER('SGBSVX',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL SGBSVX('N','N',0,0,0,-1,a,1,af,1,ip,eq,r,c,b,1,x,1,rcond, &
     &               r1,r2,w,iw,info)
         CALL CHKXER('SGBSVX',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL SGBSVX('N','N',1,1,1,0,a,2,af,4,ip,eq,r,c,b,1,x,1,rcond,  &
     &               r1,r2,w,iw,info)
         CALL CHKXER('SGBSVX',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL SGBSVX('N','N',1,1,1,0,a,3,af,3,ip,eq,r,c,b,1,x,1,rcond,  &
     &               r1,r2,w,iw,info)
         CALL CHKXER('SGBSVX',INFot,NOUt,LERr,OK)
         INFot = 12
         eq = '/'
         CALL SGBSVX('F','N',0,0,0,0,a,1,af,1,ip,eq,r,c,b,1,x,1,rcond,  &
     &               r1,r2,w,iw,info)
         CALL CHKXER('SGBSVX',INFot,NOUt,LERr,OK)
         INFot = 13
         eq = 'R'
         CALL SGBSVX('F','N',1,0,0,0,a,1,af,1,ip,eq,r,c,b,1,x,1,rcond,  &
     &               r1,r2,w,iw,info)
         CALL CHKXER('SGBSVX',INFot,NOUt,LERr,OK)
         INFot = 14
         eq = 'C'
         CALL SGBSVX('F','N',1,0,0,0,a,1,af,1,ip,eq,r,c,b,1,x,1,rcond,  &
     &               r1,r2,w,iw,info)
         CALL CHKXER('SGBSVX',INFot,NOUt,LERr,OK)
         INFot = 16
         CALL SGBSVX('N','N',2,0,0,0,a,1,af,1,ip,eq,r,c,b,1,x,2,rcond,  &
     &               r1,r2,w,iw,info)
         CALL CHKXER('SGBSVX',INFot,NOUt,LERr,OK)
         INFot = 18
         CALL SGBSVX('N','N',2,0,0,0,a,1,af,1,ip,eq,r,c,b,2,x,1,rcond,  &
     &               r1,r2,w,iw,info)
         CALL CHKXER('SGBSVX',INFot,NOUt,LERr,OK)
!
      ELSEIF ( LSAMEN(2,c2,'GT') ) THEN
!
!        SGTSV
!
         SRNamt = 'SGTSV '
         INFot = 1
         CALL SGTSV(-1,0,a(1,1),a(1,2),a(1,3),b,1,info)
         CALL CHKXER('SGTSV ',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SGTSV(0,-1,a(1,1),a(1,2),a(1,3),b,1,info)
         CALL CHKXER('SGTSV ',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL SGTSV(2,0,a(1,1),a(1,2),a(1,3),b,1,info)
         CALL CHKXER('SGTSV ',INFot,NOUt,LERr,OK)
!
!        SGTSVX
!
         SRNamt = 'SGTSVX'
         INFot = 1
         CALL SGTSVX('/','N',0,0,a(1,1),a(1,2),a(1,3),af(1,1),af(1,2),  &
     &               af(1,3),af(1,4),ip,b,1,x,1,rcond,r1,r2,w,iw,info)
         CALL CHKXER('SGTSVX',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SGTSVX('N','/',0,0,a(1,1),a(1,2),a(1,3),af(1,1),af(1,2),  &
     &               af(1,3),af(1,4),ip,b,1,x,1,rcond,r1,r2,w,iw,info)
         CALL CHKXER('SGTSVX',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SGTSVX('N','N',-1,0,a(1,1),a(1,2),a(1,3),af(1,1),af(1,2), &
     &               af(1,3),af(1,4),ip,b,1,x,1,rcond,r1,r2,w,iw,info)
         CALL CHKXER('SGTSVX',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL SGTSVX('N','N',0,-1,a(1,1),a(1,2),a(1,3),af(1,1),af(1,2), &
     &               af(1,3),af(1,4),ip,b,1,x,1,rcond,r1,r2,w,iw,info)
         CALL CHKXER('SGTSVX',INFot,NOUt,LERr,OK)
         INFot = 14
         CALL SGTSVX('N','N',2,0,a(1,1),a(1,2),a(1,3),af(1,1),af(1,2),  &
     &               af(1,3),af(1,4),ip,b,1,x,2,rcond,r1,r2,w,iw,info)
         CALL CHKXER('SGTSVX',INFot,NOUt,LERr,OK)
         INFot = 16
         CALL SGTSVX('N','N',2,0,a(1,1),a(1,2),a(1,3),af(1,1),af(1,2),  &
     &               af(1,3),af(1,4),ip,b,2,x,1,rcond,r1,r2,w,iw,info)
         CALL CHKXER('SGTSVX',INFot,NOUt,LERr,OK)
!
      ELSEIF ( LSAMEN(2,c2,'PO') ) THEN
!
!        SPOSV
!
         SRNamt = 'SPOSV '
         INFot = 1
         CALL SPOSV('/',0,0,a,1,b,1,info)
         CALL CHKXER('SPOSV ',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SPOSV('U',-1,0,a,1,b,1,info)
         CALL CHKXER('SPOSV ',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SPOSV('U',0,-1,a,1,b,1,info)
         CALL CHKXER('SPOSV ',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL SPOSV('U',2,0,a,1,b,2,info)
         CALL CHKXER('SPOSV ',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL SPOSV('U',2,0,a,2,b,1,info)
         CALL CHKXER('SPOSV ',INFot,NOUt,LERr,OK)
!
!        SPOSVX
!
         SRNamt = 'SPOSVX'
         INFot = 1
         CALL SPOSVX('/','U',0,0,a,1,af,1,eq,c,b,1,x,1,rcond,r1,r2,w,iw,&
     &               info)
         CALL CHKXER('SPOSVX',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SPOSVX('N','/',0,0,a,1,af,1,eq,c,b,1,x,1,rcond,r1,r2,w,iw,&
     &               info)
         CALL CHKXER('SPOSVX',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SPOSVX('N','U',-1,0,a,1,af,1,eq,c,b,1,x,1,rcond,r1,r2,w,  &
     &               iw,info)
         CALL CHKXER('SPOSVX',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL SPOSVX('N','U',0,-1,a,1,af,1,eq,c,b,1,x,1,rcond,r1,r2,w,  &
     &               iw,info)
         CALL CHKXER('SPOSVX',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL SPOSVX('N','U',2,0,a,1,af,2,eq,c,b,2,x,2,rcond,r1,r2,w,iw,&
     &               info)
         CALL CHKXER('SPOSVX',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL SPOSVX('N','U',2,0,a,2,af,1,eq,c,b,2,x,2,rcond,r1,r2,w,iw,&
     &               info)
         CALL CHKXER('SPOSVX',INFot,NOUt,LERr,OK)
         INFot = 9
         eq = '/'
         CALL SPOSVX('F','U',0,0,a,1,af,1,eq,c,b,1,x,1,rcond,r1,r2,w,iw,&
     &               info)
         CALL CHKXER('SPOSVX',INFot,NOUt,LERr,OK)
         INFot = 10
         eq = 'Y'
         CALL SPOSVX('F','U',1,0,a,1,af,1,eq,c,b,1,x,1,rcond,r1,r2,w,iw,&
     &               info)
         CALL CHKXER('SPOSVX',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL SPOSVX('N','U',2,0,a,2,af,2,eq,c,b,1,x,2,rcond,r1,r2,w,iw,&
     &               info)
         CALL CHKXER('SPOSVX',INFot,NOUt,LERr,OK)
         INFot = 14
         CALL SPOSVX('N','U',2,0,a,2,af,2,eq,c,b,2,x,1,rcond,r1,r2,w,iw,&
     &               info)
         CALL CHKXER('SPOSVX',INFot,NOUt,LERr,OK)
!
      ELSEIF ( LSAMEN(2,c2,'PP') ) THEN
!
!        SPPSV
!
         SRNamt = 'SPPSV '
         INFot = 1
         CALL SPPSV('/',0,0,a,b,1,info)
         CALL CHKXER('SPPSV ',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SPPSV('U',-1,0,a,b,1,info)
         CALL CHKXER('SPPSV ',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SPPSV('U',0,-1,a,b,1,info)
         CALL CHKXER('SPPSV ',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL SPPSV('U',2,0,a,b,1,info)
         CALL CHKXER('SPPSV ',INFot,NOUt,LERr,OK)
!
!        SPPSVX
!
         SRNamt = 'SPPSVX'
         INFot = 1
         CALL SPPSVX('/','U',0,0,a,af,eq,c,b,1,x,1,rcond,r1,r2,w,iw,    &
     &               info)
         CALL CHKXER('SPPSVX',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SPPSVX('N','/',0,0,a,af,eq,c,b,1,x,1,rcond,r1,r2,w,iw,    &
     &               info)
         CALL CHKXER('SPPSVX',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SPPSVX('N','U',-1,0,a,af,eq,c,b,1,x,1,rcond,r1,r2,w,iw,   &
     &               info)
         CALL CHKXER('SPPSVX',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL SPPSVX('N','U',0,-1,a,af,eq,c,b,1,x,1,rcond,r1,r2,w,iw,   &
     &               info)
         CALL CHKXER('SPPSVX',INFot,NOUt,LERr,OK)
         INFot = 7
         eq = '/'
         CALL SPPSVX('F','U',0,0,a,af,eq,c,b,1,x,1,rcond,r1,r2,w,iw,    &
     &               info)
         CALL CHKXER('SPPSVX',INFot,NOUt,LERr,OK)
         INFot = 8
         eq = 'Y'
         CALL SPPSVX('F','U',1,0,a,af,eq,c,b,1,x,1,rcond,r1,r2,w,iw,    &
     &               info)
         CALL CHKXER('SPPSVX',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL SPPSVX('N','U',2,0,a,af,eq,c,b,1,x,2,rcond,r1,r2,w,iw,    &
     &               info)
         CALL CHKXER('SPPSVX',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL SPPSVX('N','U',2,0,a,af,eq,c,b,2,x,1,rcond,r1,r2,w,iw,    &
     &               info)
         CALL CHKXER('SPPSVX',INFot,NOUt,LERr,OK)
!
      ELSEIF ( LSAMEN(2,c2,'PB') ) THEN
!
!        SPBSV
!
         SRNamt = 'SPBSV '
         INFot = 1
         CALL SPBSV('/',0,0,0,a,1,b,1,info)
         CALL CHKXER('SPBSV ',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SPBSV('U',-1,0,0,a,1,b,1,info)
         CALL CHKXER('SPBSV ',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SPBSV('U',1,-1,0,a,1,b,1,info)
         CALL CHKXER('SPBSV ',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL SPBSV('U',0,0,-1,a,1,b,1,info)
         CALL CHKXER('SPBSV ',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL SPBSV('U',1,1,0,a,1,b,2,info)
         CALL CHKXER('SPBSV ',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL SPBSV('U',2,0,0,a,1,b,1,info)
         CALL CHKXER('SPBSV ',INFot,NOUt,LERr,OK)
!
!        SPBSVX
!
         SRNamt = 'SPBSVX'
         INFot = 1
         CALL SPBSVX('/','U',0,0,0,a,1,af,1,eq,c,b,1,x,1,rcond,r1,r2,w, &
     &               iw,info)
         CALL CHKXER('SPBSVX',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SPBSVX('N','/',0,0,0,a,1,af,1,eq,c,b,1,x,1,rcond,r1,r2,w, &
     &               iw,info)
         CALL CHKXER('SPBSVX',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SPBSVX('N','U',-1,0,0,a,1,af,1,eq,c,b,1,x,1,rcond,r1,r2,w,&
     &               iw,info)
         CALL CHKXER('SPBSVX',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL SPBSVX('N','U',1,-1,0,a,1,af,1,eq,c,b,1,x,1,rcond,r1,r2,w,&
     &               iw,info)
         CALL CHKXER('SPBSVX',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL SPBSVX('N','U',0,0,-1,a,1,af,1,eq,c,b,1,x,1,rcond,r1,r2,w,&
     &               iw,info)
         CALL CHKXER('SPBSVX',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL SPBSVX('N','U',1,1,0,a,1,af,2,eq,c,b,2,x,2,rcond,r1,r2,w, &
     &               iw,info)
         CALL CHKXER('SPBSVX',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL SPBSVX('N','U',1,1,0,a,2,af,1,eq,c,b,2,x,2,rcond,r1,r2,w, &
     &               iw,info)
         CALL CHKXER('SPBSVX',INFot,NOUt,LERr,OK)
         INFot = 10
         eq = '/'
         CALL SPBSVX('F','U',0,0,0,a,1,af,1,eq,c,b,1,x,1,rcond,r1,r2,w, &
     &               iw,info)
         CALL CHKXER('SPBSVX',INFot,NOUt,LERr,OK)
         INFot = 11
         eq = 'Y'
         CALL SPBSVX('F','U',1,0,0,a,1,af,1,eq,c,b,1,x,1,rcond,r1,r2,w, &
     &               iw,info)
         CALL CHKXER('SPBSVX',INFot,NOUt,LERr,OK)
         INFot = 13
         CALL SPBSVX('N','U',2,0,0,a,1,af,1,eq,c,b,1,x,2,rcond,r1,r2,w, &
     &               iw,info)
         CALL CHKXER('SPBSVX',INFot,NOUt,LERr,OK)
         INFot = 15
         CALL SPBSVX('N','U',2,0,0,a,1,af,1,eq,c,b,2,x,1,rcond,r1,r2,w, &
     &               iw,info)
         CALL CHKXER('SPBSVX',INFot,NOUt,LERr,OK)
!
      ELSEIF ( LSAMEN(2,c2,'PT') ) THEN
!
!        SPTSV
!
         SRNamt = 'SPTSV '
         INFot = 1
         CALL SPTSV(-1,0,a(1,1),a(1,2),b,1,info)
         CALL CHKXER('SPTSV ',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SPTSV(0,-1,a(1,1),a(1,2),b,1,info)
         CALL CHKXER('SPTSV ',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL SPTSV(2,0,a(1,1),a(1,2),b,1,info)
         CALL CHKXER('SPTSV ',INFot,NOUt,LERr,OK)
!
!        SPTSVX
!
         SRNamt = 'SPTSVX'
         INFot = 1
         CALL SPTSVX('/',0,0,a(1,1),a(1,2),af(1,1),af(1,2),b,1,x,1,     &
     &               rcond,r1,r2,w,info)
         CALL CHKXER('SPTSVX',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SPTSVX('N',-1,0,a(1,1),a(1,2),af(1,1),af(1,2),b,1,x,1,    &
     &               rcond,r1,r2,w,info)
         CALL CHKXER('SPTSVX',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SPTSVX('N',0,-1,a(1,1),a(1,2),af(1,1),af(1,2),b,1,x,1,    &
     &               rcond,r1,r2,w,info)
         CALL CHKXER('SPTSVX',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL SPTSVX('N',2,0,a(1,1),a(1,2),af(1,1),af(1,2),b,1,x,2,     &
     &               rcond,r1,r2,w,info)
         CALL CHKXER('SPTSVX',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL SPTSVX('N',2,0,a(1,1),a(1,2),af(1,1),af(1,2),b,2,x,1,     &
     &               rcond,r1,r2,w,info)
         CALL CHKXER('SPTSVX',INFot,NOUt,LERr,OK)
!
      ELSEIF ( LSAMEN(2,c2,'SY') ) THEN
!
!        SSYSV
!
         SRNamt = 'SSYSV '
         INFot = 1
         CALL SSYSV('/',0,0,a,1,ip,b,1,w,1,info)
         CALL CHKXER('SSYSV ',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SSYSV('U',-1,0,a,1,ip,b,1,w,1,info)
         CALL CHKXER('SSYSV ',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SSYSV('U',0,-1,a,1,ip,b,1,w,1,info)
         CALL CHKXER('SSYSV ',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL SSYSV('U',2,0,a,2,ip,b,1,w,1,info)
         CALL CHKXER('SSYSV ',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL SSYSV('U',0,0,a,1,ip,b,1,w,0,info)
         CALL CHKXER('SSYSV ',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL SSYSV('U',0,0,a,1,ip,b,1,w,-2,info)
         CALL CHKXER('SSYSV ',INFot,NOUt,LERr,OK)
!
!        SSYSVX
!
         SRNamt = 'SSYSVX'
         INFot = 1
         CALL SSYSVX('/','U',0,0,a,1,af,1,ip,b,1,x,1,rcond,r1,r2,w,1,iw,&
     &               info)
         CALL CHKXER('SSYSVX',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SSYSVX('N','/',0,0,a,1,af,1,ip,b,1,x,1,rcond,r1,r2,w,1,iw,&
     &               info)
         CALL CHKXER('SSYSVX',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SSYSVX('N','U',-1,0,a,1,af,1,ip,b,1,x,1,rcond,r1,r2,w,1,  &
     &               iw,info)
         CALL CHKXER('SSYSVX',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL SSYSVX('N','U',0,-1,a,1,af,1,ip,b,1,x,1,rcond,r1,r2,w,1,  &
     &               iw,info)
         CALL CHKXER('SSYSVX',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL SSYSVX('N','U',2,0,a,1,af,2,ip,b,2,x,2,rcond,r1,r2,w,4,iw,&
     &               info)
         CALL CHKXER('SSYSVX',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL SSYSVX('N','U',2,0,a,2,af,1,ip,b,2,x,2,rcond,r1,r2,w,4,iw,&
     &               info)
         CALL CHKXER('SSYSVX',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL SSYSVX('N','U',2,0,a,2,af,2,ip,b,1,x,2,rcond,r1,r2,w,4,iw,&
     &               info)
         CALL CHKXER('SSYSVX',INFot,NOUt,LERr,OK)
         INFot = 13
         CALL SSYSVX('N','U',2,0,a,2,af,2,ip,b,2,x,1,rcond,r1,r2,w,4,iw,&
     &               info)
         CALL CHKXER('SSYSVX',INFot,NOUt,LERr,OK)
         INFot = 18
         CALL SSYSVX('N','U',2,0,a,2,af,2,ip,b,2,x,2,rcond,r1,r2,w,3,iw,&
     &               info)
         CALL CHKXER('SSYSVX',INFot,NOUt,LERr,OK)
!
!
      ELSEIF ( LSAMEN(2,c2,'SR') ) THEN
!
!        SSYSV_ROOK
!
         SRNamt = 'SSYSV_ROOK'
         INFot = 1
         CALL SSYSV_ROOK('/',0,0,a,1,ip,b,1,w,1,info)
         CALL CHKXER('SSYSV_ROOK',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SSYSV_ROOK('U',-1,0,a,1,ip,b,1,w,1,info)
         CALL CHKXER('SSYSV_ROOK',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SSYSV_ROOK('U',0,-1,a,1,ip,b,1,w,1,info)
         CALL CHKXER('SSYSV_ROOK',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL SSYSV_ROOK('U',2,0,a,2,ip,b,1,w,1,info)
         CALL CHKXER('SSYSV_ROOK',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL SSYSV_ROOK('U',0,0,a,1,ip,b,1,w,0,info)
         CALL CHKXER('SSYSV_ROOK',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL SSYSV_ROOK('U',0,0,a,1,ip,b,1,w,-2,info)
         CALL CHKXER('SSYSV_ROOK',INFot,NOUt,LERr,OK)
!
      ELSEIF ( LSAMEN(2,c2,'SK') ) THEN
!
!        SSYSV_RK
!
!        Test error exits of the driver that uses factorization
!        of a symmetric indefinite matrix with rook
!        (bounded Bunch-Kaufman) pivoting with the new storage
!        format for factors L ( or U) and D.
!
!        L (or U) is stored in A, diagonal of D is stored on the
!        diagonal of A, subdiagonal of D is stored in a separate array E.
!
         SRNamt = 'SSYSV_RK'
         INFot = 1
         CALL SSYSV_RK('/',0,0,a,1,e,ip,b,1,w,1,info)
         CALL CHKXER('SSYSV_RK',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SSYSV_RK('U',-1,0,a,1,e,ip,b,1,w,1,info)
         CALL CHKXER('SSYSV_RK',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SSYSV_RK('U',0,-1,a,1,e,ip,b,1,w,1,info)
         CALL CHKXER('SSYSV_RK',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL SSYSV_RK('U',2,0,a,1,e,ip,b,2,w,1,info)
         CALL CHKXER('SSYSV_RK',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL SSYSV_RK('U',2,0,a,2,e,ip,b,1,w,1,info)
         CALL CHKXER('SSYSV_RK',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL SSYSV_RK('U',0,0,a,1,e,ip,b,1,w,0,info)
         CALL CHKXER('SSYSV_RK',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL SSYSV_RK('U',0,0,a,1,e,ip,b,1,w,-2,info)
         CALL CHKXER('SSYSV_RK',INFot,NOUt,LERr,OK)
!
      ELSEIF ( LSAMEN(2,c2,'SA') ) THEN
!
!        SSYSV_AA
!
         SRNamt = 'SSYSV_AA'
         INFot = 1
         CALL SSYSV_AA('/',0,0,a,1,ip,b,1,w,1,info)
         CALL CHKXER('SSYSV_AA',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SSYSV_AA('U',-1,0,a,1,ip,b,1,w,1,info)
         CALL CHKXER('SSYSV_AA',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SSYSV_AA('U',0,-1,a,1,ip,b,1,w,1,info)
         CALL CHKXER('SSYSV_AA',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL SSYSV_AA('U',2,0,a,2,ip,b,1,w,1,info)
         CALL CHKXER('SSYSV_AA',INFot,NOUt,LERr,OK)
!
      ELSEIF ( LSAMEN(2,c2,'S2') ) THEN
!
!        DSYSV_AASEN_2STAGE
!
         SRNamt = 'SSYSV_AA_2STAGE'
         INFot = 1
         CALL SSYSV_AA_2STAGE('/',0,0,a,1,a,1,ip,ip,b,1,w,1,info)
         CALL CHKXER('SSYSV_AA_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SSYSV_AA_2STAGE('U',-1,0,a,1,a,1,ip,ip,b,1,w,1,info)
         CALL CHKXER('SSYSV_AA_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SSYSV_AA_2STAGE('U',0,-1,a,1,a,1,ip,ip,b,1,w,1,info)
         CALL CHKXER('SSYSV_AA_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL SSYSV_AA_2STAGE('U',2,1,a,1,a,1,ip,ip,b,1,w,1,info)
         CALL CHKXER('SSYSV_AA_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL SSYSV_AA_2STAGE('U',2,1,a,2,a,8,ip,ip,b,1,w,1,info)
         CALL CHKXER('SSYSV_AA_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL SSYSV_AA_2STAGE('U',2,1,a,2,a,1,ip,ip,b,2,w,1,info)
         CALL CHKXER('SSYSV_AA_2STAGE',INFot,NOUt,LERr,OK)
!
      ELSEIF ( LSAMEN(2,c2,'SP') ) THEN
!
!        SSPSV
!
         SRNamt = 'SSPSV '
         INFot = 1
         CALL SSPSV('/',0,0,a,ip,b,1,info)
         CALL CHKXER('SSPSV ',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SSPSV('U',-1,0,a,ip,b,1,info)
         CALL CHKXER('SSPSV ',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SSPSV('U',0,-1,a,ip,b,1,info)
         CALL CHKXER('SSPSV ',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL SSPSV('U',2,0,a,ip,b,1,info)
         CALL CHKXER('SSPSV ',INFot,NOUt,LERr,OK)
!
!        SSPSVX
!
         SRNamt = 'SSPSVX'
         INFot = 1
         CALL SSPSVX('/','U',0,0,a,af,ip,b,1,x,1,rcond,r1,r2,w,iw,info)
         CALL CHKXER('SSPSVX',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SSPSVX('N','/',0,0,a,af,ip,b,1,x,1,rcond,r1,r2,w,iw,info)
         CALL CHKXER('SSPSVX',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SSPSVX('N','U',-1,0,a,af,ip,b,1,x,1,rcond,r1,r2,w,iw,info)
         CALL CHKXER('SSPSVX',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL SSPSVX('N','U',0,-1,a,af,ip,b,1,x,1,rcond,r1,r2,w,iw,info)
         CALL CHKXER('SSPSVX',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL SSPSVX('N','U',2,0,a,af,ip,b,1,x,2,rcond,r1,r2,w,iw,info)
         CALL CHKXER('SSPSVX',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL SSPSVX('N','U',2,0,a,af,ip,b,2,x,1,rcond,r1,r2,w,iw,info)
         CALL CHKXER('SSPSVX',INFot,NOUt,LERr,OK)
      ENDIF
!
!     Print a summary line.
!
      IF ( OK ) THEN
         WRITE (NOUt,FMT=99001) Path
      ELSE
         WRITE (NOUt,FMT=99002) Path
      ENDIF
!
99001 FORMAT (1X,A3,' drivers passed the tests of the error exits')
99002 FORMAT (' *** ',A3,' drivers failed the tests of the error ',     &
     &        'exits ***')
!
!
!     End of serrvx
!
      END SUBROUTINE serrvx
