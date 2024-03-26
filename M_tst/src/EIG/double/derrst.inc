!*==derrst.f90  processed by SPAG 7.51RB at 17:37 on  4 Mar 2022
!> \brief \b derrst
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE derrst( PATH, NUNIT )
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
!> derrst tests the error exits for DSYTRD, DORGTR, DORMTR, DSPTRD,
!> DOPGTR, DOPMTR, DSTEQR, SSTERF, SSTEBZ, SSTEIN, DPTEQR, DSBTRD,
!> DSYEV, SSYEVX, SSYEVD, DSBEV, SSBEVX, SSBEVD,
!> DSPEV, SSPEVX, SSPEVD, DSTEV, SSTEVX, SSTEVD, and SSTEDC.
!> DSYEVD_2STAGE, DSYEVR_2STAGE, DSYEVX_2STAGE,
!> DSYEV_2STAGE, DSBEV_2STAGE, DSBEVD_2STAGE,
!> DSBEVX_2STAGE, DSYTRD_2STAGE, DSYTRD_SY2SB,
!> DSYTRD_SB2ST
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
!> \ingroup double_eig
!
!  =====================================================================
      SUBROUTINE derrst(Path,Nunit)
      IMPLICIT NONE
!*--DERRST65
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
!     NMAX has to be at least 3 or LIW may be too small
!     .. Parameters ..
      INTEGER NMAX , LIW , LW
      PARAMETER (NMAX=3,LIW=12*NMAX,LW=20*NMAX)
!     ..
!     .. Local Scalars ..
      CHARACTER*2 c2
      INTEGER i , info , j , m , n , nsplit , nt
!     ..
!     .. Local Arrays ..
      INTEGER i1(NMAX) , i2(NMAX) , i3(NMAX) , iw(LIW)
      DOUBLE PRECISION a(NMAX,NMAX) , c(NMAX,NMAX) , d(NMAX) , e(NMAX) ,&
     &                 q(NMAX,NMAX) , r(NMAX) , tau(NMAX) , w(LW) ,     &
     &                 x(NMAX) , z(NMAX,NMAX)
!     ..
!     .. External Functions ..
      LOGICAL LSAMEN
      EXTERNAL LSAMEN
!     ..
!     .. External Subroutines ..
      EXTERNAL DOPGTR , DOPMTR , DORGTR , DORMTR , DPTEQR 
     EXTERNAL  DSBEV , DSBEVD , DSBEVX , DSBTRD , DSPEV , DSPEVD 
     EXTERNAL  DSPEVX , DSPTRD , DSTEBZ , DSTEDC , DSTEIN , DSTEQR 
     EXTERNAL  DSTERF , DSTEV , DSTEVD , DSTEVR , DSTEVX , DSYEV 
     EXTERNAL  DSYEVD , DSYEVR , DSYEVX , DSYTRD , DSYEVD_2STAGE 
     EXTERNAL  DSYEVR_2STAGE , DSYEVX_2STAGE , DSYEV_2STAGE 
     EXTERNAL  DSBEV_2STAGE , DSBEVD_2STAGE , DSBEVX_2STAGE 
     EXTERNAL  DSYTRD_2STAGE , DSYTRD_SY2SB , DSYTRD_SB2ST
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
      INTRINSIC DBLE
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
            a(i,j) = 1.D0/DBLE(i+j)
         ENDDO
      ENDDO
      DO j = 1 , NMAX
         d(j) = DBLE(j)
         e(j) = 0.0D0
         i1(j) = j
         i2(j) = j
         tau(j) = 1.D0
      ENDDO
      OK = .TRUE.
      nt = 0
!
!     Test error exits for the ST path.
!
      IF ( LSAMEN(2,c2,'ST') ) THEN
!
!        DSYTRD
!
         SRNamt = 'DSYTRD'
         INFot = 1
         CALL DSYTRD('/',0,a,1,d,e,tau,w,1,info)
         CALL chkxer('DSYTRD',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DSYTRD('U',-1,a,1,d,e,tau,w,1,info)
         CALL chkxer('DSYTRD',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL DSYTRD('U',2,a,1,d,e,tau,w,1,info)
         CALL chkxer('DSYTRD',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL DSYTRD('U',0,a,1,d,e,tau,w,0,info)
         CALL chkxer('DSYTRD',INFot,NOUt,LERr,OK)
         nt = nt + 4
!
!        DSYTRD_2STAGE
!
         SRNamt = 'DSYTRD_2STAGE'
         INFot = 1
         CALL DSYTRD_2STAGE('/','U',0,a,1,d,e,tau,c,1,w,1,info)
         CALL chkxer('DSYTRD_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 1
         CALL DSYTRD_2STAGE('H','U',0,a,1,d,e,tau,c,1,w,1,info)
         CALL chkxer('DSYTRD_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DSYTRD_2STAGE('N','/',0,a,1,d,e,tau,c,1,w,1,info)
         CALL chkxer('DSYTRD_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DSYTRD_2STAGE('N','U',-1,a,1,d,e,tau,c,1,w,1,info)
         CALL chkxer('DSYTRD_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL DSYTRD_2STAGE('N','U',2,a,1,d,e,tau,c,1,w,1,info)
         CALL chkxer('DSYTRD_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL DSYTRD_2STAGE('N','U',0,a,1,d,e,tau,c,0,w,1,info)
         CALL chkxer('DSYTRD_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL DSYTRD_2STAGE('N','U',0,a,1,d,e,tau,c,1,w,0,info)
         CALL chkxer('DSYTRD_2STAGE',INFot,NOUt,LERr,OK)
         nt = nt + 7
!
!        DSYTRD_SY2SB
!
         SRNamt = 'DSYTRD_SY2SB'
         INFot = 1
         CALL DSYTRD_SY2SB('/',0,0,a,1,c,1,tau,w,1,info)
         CALL chkxer('DSYTRD_SY2SB',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DSYTRD_SY2SB('U',-1,0,a,1,c,1,tau,w,1,info)
         CALL chkxer('DSYTRD_SY2SB',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DSYTRD_SY2SB('U',0,-1,a,1,c,1,tau,w,1,info)
         CALL chkxer('DSYTRD_SY2SB',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL DSYTRD_SY2SB('U',2,0,a,1,c,1,tau,w,1,info)
         CALL chkxer('DSYTRD_SY2SB',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL DSYTRD_SY2SB('U',0,2,a,1,c,1,tau,w,1,info)
         CALL chkxer('DSYTRD_SY2SB',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL DSYTRD_SY2SB('U',0,0,a,1,c,1,tau,w,0,info)
         CALL chkxer('DSYTRD_SY2SB',INFot,NOUt,LERr,OK)
         nt = nt + 6
!
!        DSYTRD_SB2ST
!
         SRNamt = 'DSYTRD_SB2ST'
         INFot = 1
         CALL DSYTRD_SB2ST('/','N','U',0,0,a,1,d,e,c,1,w,1,info)
         CALL chkxer('DSYTRD_SB2ST',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DSYTRD_SB2ST('Y','/','U',0,0,a,1,d,e,c,1,w,1,info)
         CALL chkxer('DSYTRD_SB2ST',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DSYTRD_SB2ST('Y','H','U',0,0,a,1,d,e,c,1,w,1,info)
         CALL chkxer('DSYTRD_SB2ST',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DSYTRD_SB2ST('Y','N','/',0,0,a,1,d,e,c,1,w,1,info)
         CALL chkxer('DSYTRD_SB2ST',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL DSYTRD_SB2ST('Y','N','U',-1,0,a,1,d,e,c,1,w,1,info)
         CALL chkxer('DSYTRD_SB2ST',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL DSYTRD_SB2ST('Y','N','U',0,-1,a,1,d,e,c,1,w,1,info)
         CALL chkxer('DSYTRD_SB2ST',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL DSYTRD_SB2ST('Y','N','U',0,1,a,1,d,e,c,1,w,1,info)
         CALL chkxer('DSYTRD_SB2ST',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL DSYTRD_SB2ST('Y','N','U',0,0,a,1,d,e,c,0,w,1,info)
         CALL chkxer('DSYTRD_SB2ST',INFot,NOUt,LERr,OK)
         INFot = 13
         CALL DSYTRD_SB2ST('Y','N','U',0,0,a,1,d,e,c,1,w,0,info)
         CALL chkxer('DSYTRD_SB2ST',INFot,NOUt,LERr,OK)
         nt = nt + 9
!
!        DORGTR
!
         SRNamt = 'DORGTR'
         INFot = 1
         CALL DORGTR('/',0,a,1,tau,w,1,info)
         CALL chkxer('DORGTR',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DORGTR('U',-1,a,1,tau,w,1,info)
         CALL chkxer('DORGTR',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL DORGTR('U',2,a,1,tau,w,1,info)
         CALL chkxer('DORGTR',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL DORGTR('U',3,a,3,tau,w,1,info)
         CALL chkxer('DORGTR',INFot,NOUt,LERr,OK)
         nt = nt + 4
!
!        DORMTR
!
         SRNamt = 'DORMTR'
         INFot = 1
         CALL DORMTR('/','U','N',0,0,a,1,tau,c,1,w,1,info)
         CALL chkxer('DORMTR',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DORMTR('L','/','N',0,0,a,1,tau,c,1,w,1,info)
         CALL chkxer('DORMTR',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DORMTR('L','U','/',0,0,a,1,tau,c,1,w,1,info)
         CALL chkxer('DORMTR',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL DORMTR('L','U','N',-1,0,a,1,tau,c,1,w,1,info)
         CALL chkxer('DORMTR',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL DORMTR('L','U','N',0,-1,a,1,tau,c,1,w,1,info)
         CALL chkxer('DORMTR',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL DORMTR('L','U','N',2,0,a,1,tau,c,2,w,1,info)
         CALL chkxer('DORMTR',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL DORMTR('R','U','N',0,2,a,1,tau,c,1,w,1,info)
         CALL chkxer('DORMTR',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL DORMTR('L','U','N',2,0,a,2,tau,c,1,w,1,info)
         CALL chkxer('DORMTR',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL DORMTR('L','U','N',0,2,a,1,tau,c,1,w,1,info)
         CALL chkxer('DORMTR',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL DORMTR('R','U','N',2,0,a,1,tau,c,2,w,1,info)
         CALL chkxer('DORMTR',INFot,NOUt,LERr,OK)
         nt = nt + 10
!
!        DSPTRD
!
         SRNamt = 'DSPTRD'
         INFot = 1
         CALL DSPTRD('/',0,a,d,e,tau,info)
         CALL chkxer('DSPTRD',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DSPTRD('U',-1,a,d,e,tau,info)
         CALL chkxer('DSPTRD',INFot,NOUt,LERr,OK)
         nt = nt + 2
!
!        DOPGTR
!
         SRNamt = 'DOPGTR'
         INFot = 1
         CALL DOPGTR('/',0,a,tau,z,1,w,info)
         CALL chkxer('DOPGTR',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DOPGTR('U',-1,a,tau,z,1,w,info)
         CALL chkxer('DOPGTR',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL DOPGTR('U',2,a,tau,z,1,w,info)
         CALL chkxer('DOPGTR',INFot,NOUt,LERr,OK)
         nt = nt + 3
!
!        DOPMTR
!
         SRNamt = 'DOPMTR'
         INFot = 1
         CALL DOPMTR('/','U','N',0,0,a,tau,c,1,w,info)
         CALL chkxer('DOPMTR',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DOPMTR('L','/','N',0,0,a,tau,c,1,w,info)
         CALL chkxer('DOPMTR',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DOPMTR('L','U','/',0,0,a,tau,c,1,w,info)
         CALL chkxer('DOPMTR',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL DOPMTR('L','U','N',-1,0,a,tau,c,1,w,info)
         CALL chkxer('DOPMTR',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL DOPMTR('L','U','N',0,-1,a,tau,c,1,w,info)
         CALL chkxer('DOPMTR',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL DOPMTR('L','U','N',2,0,a,tau,c,1,w,info)
         CALL chkxer('DOPMTR',INFot,NOUt,LERr,OK)
         nt = nt + 6
!
!        DPTEQR
!
         SRNamt = 'DPTEQR'
         INFot = 1
         CALL DPTEQR('/',0,d,e,z,1,w,info)
         CALL chkxer('DPTEQR',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DPTEQR('N',-1,d,e,z,1,w,info)
         CALL chkxer('DPTEQR',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL DPTEQR('V',2,d,e,z,1,w,info)
         CALL chkxer('DPTEQR',INFot,NOUt,LERr,OK)
         nt = nt + 3
!
!        DSTEBZ
!
         SRNamt = 'DSTEBZ'
         INFot = 1
         CALL DSTEBZ('/','E',0,0.0D0,1.0D0,1,0,0.0D0,d,e,m,nsplit,x,i1, &
     &               i2,w,iw,info)
         CALL chkxer('DSTEBZ',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DSTEBZ('A','/',0,0.0D0,0.0D0,0,0,0.0D0,d,e,m,nsplit,x,i1, &
     &               i2,w,iw,info)
         CALL chkxer('DSTEBZ',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DSTEBZ('A','E',-1,0.0D0,0.0D0,0,0,0.0D0,d,e,m,nsplit,x,i1,&
     &               i2,w,iw,info)
         CALL chkxer('DSTEBZ',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL DSTEBZ('V','E',0,0.0D0,0.0D0,0,0,0.0D0,d,e,m,nsplit,x,i1, &
     &               i2,w,iw,info)
         CALL chkxer('DSTEBZ',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL DSTEBZ('I','E',0,0.0D0,0.0D0,0,0,0.0D0,d,e,m,nsplit,x,i1, &
     &               i2,w,iw,info)
         CALL chkxer('DSTEBZ',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL DSTEBZ('I','E',1,0.0D0,0.0D0,2,1,0.0D0,d,e,m,nsplit,x,i1, &
     &               i2,w,iw,info)
         CALL chkxer('DSTEBZ',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL DSTEBZ('I','E',1,0.0D0,0.0D0,1,0,0.0D0,d,e,m,nsplit,x,i1, &
     &               i2,w,iw,info)
         CALL chkxer('DSTEBZ',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL DSTEBZ('I','E',1,0.0D0,0.0D0,1,2,0.0D0,d,e,m,nsplit,x,i1, &
     &               i2,w,iw,info)
         CALL chkxer('DSTEBZ',INFot,NOUt,LERr,OK)
         nt = nt + 8
!
!        DSTEIN
!
         SRNamt = 'DSTEIN'
         INFot = 1
         CALL DSTEIN(-1,d,e,0,x,i1,i2,z,1,w,iw,i3,info)
         CALL chkxer('DSTEIN',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL DSTEIN(0,d,e,-1,x,i1,i2,z,1,w,iw,i3,info)
         CALL chkxer('DSTEIN',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL DSTEIN(0,d,e,1,x,i1,i2,z,1,w,iw,i3,info)
         CALL chkxer('DSTEIN',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL DSTEIN(2,d,e,0,x,i1,i2,z,1,w,iw,i3,info)
         CALL chkxer('DSTEIN',INFot,NOUt,LERr,OK)
         nt = nt + 4
!
!        DSTEQR
!
         SRNamt = 'DSTEQR'
         INFot = 1
         CALL DSTEQR('/',0,d,e,z,1,w,info)
         CALL chkxer('DSTEQR',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DSTEQR('N',-1,d,e,z,1,w,info)
         CALL chkxer('DSTEQR',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL DSTEQR('V',2,d,e,z,1,w,info)
         CALL chkxer('DSTEQR',INFot,NOUt,LERr,OK)
         nt = nt + 3
!
!        DSTERF
!
         SRNamt = 'DSTERF'
         INFot = 1
         CALL DSTERF(-1,d,e,info)
         CALL chkxer('DSTERF',INFot,NOUt,LERr,OK)
         nt = nt + 1
!
!        DSTEDC
!
         SRNamt = 'DSTEDC'
         INFot = 1
         CALL DSTEDC('/',0,d,e,z,1,w,1,iw,1,info)
         CALL chkxer('DSTEDC',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DSTEDC('N',-1,d,e,z,1,w,1,iw,1,info)
         CALL chkxer('DSTEDC',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL DSTEDC('V',2,d,e,z,1,w,23,iw,28,info)
         CALL chkxer('DSTEDC',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL DSTEDC('N',1,d,e,z,1,w,0,iw,1,info)
         CALL chkxer('DSTEDC',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL DSTEDC('I',2,d,e,z,2,w,0,iw,12,info)
         CALL chkxer('DSTEDC',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL DSTEDC('V',2,d,e,z,2,w,0,iw,28,info)
         CALL chkxer('DSTEDC',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL DSTEDC('N',1,d,e,z,1,w,1,iw,0,info)
         CALL chkxer('DSTEDC',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL DSTEDC('I',2,d,e,z,2,w,19,iw,0,info)
         CALL chkxer('DSTEDC',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL DSTEDC('V',2,d,e,z,2,w,23,iw,0,info)
         CALL chkxer('DSTEDC',INFot,NOUt,LERr,OK)
         nt = nt + 9
!
!        DSTEVD
!
         SRNamt = 'DSTEVD'
         INFot = 1
         CALL DSTEVD('/',0,d,e,z,1,w,1,iw,1,info)
         CALL chkxer('DSTEVD',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DSTEVD('N',-1,d,e,z,1,w,1,iw,1,info)
         CALL chkxer('DSTEVD',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL DSTEVD('V',2,d,e,z,1,w,19,iw,12,info)
         CALL chkxer('DSTEVD',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL DSTEVD('N',1,d,e,z,1,w,0,iw,1,info)
         CALL chkxer('DSTEVD',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL DSTEVD('V',2,d,e,z,2,w,12,iw,12,info)
         CALL chkxer('DSTEVD',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL DSTEVD('N',0,d,e,z,1,w,1,iw,0,info)
         CALL chkxer('DSTEVD',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL DSTEVD('V',2,d,e,z,2,w,19,iw,11,info)
         CALL chkxer('DSTEVD',INFot,NOUt,LERr,OK)
         nt = nt + 7
!
!        DSTEV
!
         SRNamt = 'DSTEV '
         INFot = 1
         CALL DSTEV('/',0,d,e,z,1,w,info)
         CALL chkxer('DSTEV ',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DSTEV('N',-1,d,e,z,1,w,info)
         CALL chkxer('DSTEV ',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL DSTEV('V',2,d,e,z,1,w,info)
         CALL chkxer('DSTEV ',INFot,NOUt,LERr,OK)
         nt = nt + 3
!
!        DSTEVX
!
         SRNamt = 'DSTEVX'
         INFot = 1
         CALL DSTEVX('/','A',0,d,e,0.0D0,0.0D0,0,0,0.0D0,m,x,z,1,w,iw,  &
     &               i3,info)
         CALL chkxer('DSTEVX',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DSTEVX('N','/',0,d,e,0.0D0,1.0D0,1,0,0.0D0,m,x,z,1,w,iw,  &
     &               i3,info)
         CALL chkxer('DSTEVX',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DSTEVX('N','A',-1,d,e,0.0D0,0.0D0,0,0,0.0D0,m,x,z,1,w,iw, &
     &               i3,info)
         CALL chkxer('DSTEVX',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL DSTEVX('N','V',1,d,e,0.0D0,0.0D0,0,0,0.0D0,m,x,z,1,w,iw,  &
     &               i3,info)
         CALL chkxer('DSTEVX',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL DSTEVX('N','I',1,d,e,0.0D0,0.0D0,0,0,0.0D0,m,x,z,1,w,iw,  &
     &               i3,info)
         CALL chkxer('DSTEVX',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL DSTEVX('N','I',1,d,e,0.0D0,0.0D0,2,1,0.0D0,m,x,z,1,w,iw,  &
     &               i3,info)
         CALL chkxer('DSTEVX',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL DSTEVX('N','I',2,d,e,0.0D0,0.0D0,2,1,0.0D0,m,x,z,1,w,iw,  &
     &               i3,info)
         CALL chkxer('DSTEVX',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL DSTEVX('N','I',1,d,e,0.0D0,0.0D0,1,2,0.0D0,m,x,z,1,w,iw,  &
     &               i3,info)
         CALL chkxer('DSTEVX',INFot,NOUt,LERr,OK)
         INFot = 14
         CALL DSTEVX('V','A',2,d,e,0.0D0,0.0D0,0,0,0.0D0,m,x,z,1,w,iw,  &
     &               i3,info)
         CALL chkxer('DSTEVX',INFot,NOUt,LERr,OK)
         nt = nt + 9
!
!        DSTEVR
!
         n = 1
         SRNamt = 'DSTEVR'
         INFot = 1
         CALL DSTEVR('/','A',0,d,e,0.0D0,0.0D0,1,1,0.0D0,m,r,z,1,iw,x,  &
     &               20*n,iw(2*n+1),10*n,info)
         CALL chkxer('DSTEVR',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DSTEVR('V','/',0,d,e,0.0D0,0.0D0,1,1,0.0D0,m,r,z,1,iw,x,  &
     &               20*n,iw(2*n+1),10*n,info)
         CALL chkxer('DSTEVR',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DSTEVR('V','A',-1,d,e,0.0D0,0.0D0,1,1,0.0D0,m,r,z,1,iw,x, &
     &               20*n,iw(2*n+1),10*n,info)
         CALL chkxer('DSTEVR',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL DSTEVR('V','V',1,d,e,0.0D0,0.0D0,1,1,0.0D0,m,r,z,1,iw,x,  &
     &               20*n,iw(2*n+1),10*n,info)
         CALL chkxer('DSTEVR',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL DSTEVR('V','I',1,d,e,0.0D0,0.0D0,0,1,0.0D0,m,w,z,1,iw,x,  &
     &               20*n,iw(2*n+1),10*n,info)
         CALL chkxer('DSTEVR',INFot,NOUt,LERr,OK)
         INFot = 9
         n = 2
         CALL DSTEVR('V','I',2,d,e,0.0D0,0.0D0,2,1,0.0D0,m,w,z,1,iw,x,  &
     &               20*n,iw(2*n+1),10*n,info)
         CALL chkxer('DSTEVR',INFot,NOUt,LERr,OK)
         INFot = 14
         n = 1
         CALL DSTEVR('V','I',1,d,e,0.0D0,0.0D0,1,1,0.0D0,m,w,z,0,iw,x,  &
     &               20*n,iw(2*n+1),10*n,info)
         CALL chkxer('DSTEVR',INFot,NOUt,LERr,OK)
         INFot = 17
         CALL DSTEVR('V','I',1,d,e,0.0D0,0.0D0,1,1,0.0D0,m,w,z,1,iw,x,  &
     &               20*n-1,iw(2*n+1),10*n,info)
         CALL chkxer('DSTEVR',INFot,NOUt,LERr,OK)
         INFot = 19
         CALL DSTEVR('V','I',1,d,e,0.0D0,0.0D0,1,1,0.0D0,m,w,z,1,iw,x,  &
     &               20*n,iw(2*n+1),10*n-1,info)
         CALL chkxer('DSTEVR',INFot,NOUt,LERr,OK)
         nt = nt + 9
!
!        DSYEVD
!
         SRNamt = 'DSYEVD'
         INFot = 1
         CALL DSYEVD('/','U',0,a,1,x,w,1,iw,1,info)
         CALL chkxer('DSYEVD',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DSYEVD('N','/',0,a,1,x,w,1,iw,1,info)
         CALL chkxer('DSYEVD',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DSYEVD('N','U',-1,a,1,x,w,1,iw,1,info)
         CALL chkxer('DSYEVD',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL DSYEVD('N','U',2,a,1,x,w,3,iw,1,info)
         CALL chkxer('DSYEVD',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL DSYEVD('N','U',1,a,1,x,w,0,iw,1,info)
         CALL chkxer('DSYEVD',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL DSYEVD('N','U',2,a,2,x,w,4,iw,1,info)
         CALL chkxer('DSYEVD',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL DSYEVD('V','U',2,a,2,x,w,20,iw,12,info)
         CALL chkxer('DSYEVD',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL DSYEVD('N','U',1,a,1,x,w,1,iw,0,info)
         CALL chkxer('DSYEVD',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL DSYEVD('N','U',2,a,2,x,w,5,iw,0,info)
         CALL chkxer('DSYEVD',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL DSYEVD('V','U',2,a,2,x,w,27,iw,11,info)
         CALL chkxer('DSYEVD',INFot,NOUt,LERr,OK)
         nt = nt + 10
!
!        DSYEVD_2STAGE
!
         SRNamt = 'DSYEVD_2STAGE'
         INFot = 1
         CALL DSYEVD_2STAGE('/','U',0,a,1,x,w,1,iw,1,info)
         CALL chkxer('DSYEVD_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 1
         CALL DSYEVD_2STAGE('V','U',0,a,1,x,w,1,iw,1,info)
         CALL chkxer('DSYEVD_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DSYEVD_2STAGE('N','/',0,a,1,x,w,1,iw,1,info)
         CALL chkxer('DSYEVD_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DSYEVD_2STAGE('N','U',-1,a,1,x,w,1,iw,1,info)
         CALL chkxer('DSYEVD_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL DSYEVD_2STAGE('N','U',2,a,1,x,w,3,iw,1,info)
         CALL chkxer('DSYEVD_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL DSYEVD_2STAGE('N','U',1,a,1,x,w,0,iw,1,info)
         CALL chkxer('DSYEVD_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL DSYEVD_2STAGE('N','U',2,a,2,x,w,4,iw,1,info)
         CALL chkxer('DSYEVD_2STAGE',INFot,NOUt,LERr,OK)
!         INFOT = 8
!         CALL DSYEVD_2STAGE( 'V', 'U', 2, A, 2, X, W, 20, IW, 12, INFO )
!         CALL chkxer( 'DSYEVD_2STAGE', INFOT, NOUT, LERR, OK )
         INFot = 10
         CALL DSYEVD_2STAGE('N','U',1,a,1,x,w,1,iw,0,info)
         CALL chkxer('DSYEVD_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL DSYEVD_2STAGE('N','U',2,a,2,x,w,25,iw,0,info)
         CALL chkxer('DSYEVD_2STAGE',INFot,NOUt,LERr,OK)
!         INFOT = 10
!         CALL DSYEVD_2STAGE( 'V', 'U', 2, A, 2, X, W, 27, IW, 11, INFO )
!         CALL chkxer( 'DSYEVD_2STAGE', INFOT, NOUT, LERR, OK )
         nt = nt + 9
!
!        DSYEVR
!
         SRNamt = 'DSYEVR'
         n = 1
         INFot = 1
         CALL DSYEVR('/','A','U',0,a,1,0.0D0,0.0D0,1,1,0.0D0,m,r,z,1,iw,&
     &               q,26*n,iw(2*n+1),10*n,info)
         CALL chkxer('DSYEVR',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DSYEVR('V','/','U',0,a,1,0.0D0,0.0D0,1,1,0.0D0,m,r,z,1,iw,&
     &               q,26*n,iw(2*n+1),10*n,info)
         CALL chkxer('DSYEVR',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DSYEVR('V','A','/',-1,a,1,0.0D0,0.0D0,1,1,0.0D0,m,r,z,1,  &
     &               iw,q,26*n,iw(2*n+1),10*n,info)
         CALL chkxer('DSYEVR',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL DSYEVR('V','A','U',-1,a,1,0.0D0,0.0D0,1,1,0.0D0,m,r,z,1,  &
     &               iw,q,26*n,iw(2*n+1),10*n,info)
         CALL chkxer('DSYEVR',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL DSYEVR('V','A','U',2,a,1,0.0D0,0.0D0,1,1,0.0D0,m,r,z,1,iw,&
     &               q,26*n,iw(2*n+1),10*n,info)
         CALL chkxer('DSYEVR',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL DSYEVR('V','V','U',1,a,1,0.0D0,0.0D0,1,1,0.0D0,m,r,z,1,iw,&
     &               q,26*n,iw(2*n+1),10*n,info)
         CALL chkxer('DSYEVR',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL DSYEVR('V','I','U',1,a,1,0.0D0,0.0D0,0,1,0.0D0,m,r,z,1,iw,&
     &               q,26*n,iw(2*n+1),10*n,info)
         CALL chkxer('DSYEVR',INFot,NOUt,LERr,OK)
         INFot = 10
!
         CALL DSYEVR('V','I','U',2,a,2,0.0D0,0.0D0,2,1,0.0D0,m,r,z,1,iw,&
     &               q,26*n,iw(2*n+1),10*n,info)
         CALL chkxer('DSYEVR',INFot,NOUt,LERr,OK)
         INFot = 15
         CALL DSYEVR('V','I','U',1,a,1,0.0D0,0.0D0,1,1,0.0D0,m,r,z,0,iw,&
     &               q,26*n,iw(2*n+1),10*n,info)
         CALL chkxer('DSYEVR',INFot,NOUt,LERr,OK)
         INFot = 18
         CALL DSYEVR('V','I','U',1,a,1,0.0D0,0.0D0,1,1,0.0D0,m,r,z,1,iw,&
     &               q,26*n-1,iw(2*n+1),10*n,info)
         CALL chkxer('DSYEVR',INFot,NOUt,LERr,OK)
         INFot = 20
         CALL DSYEVR('V','I','U',1,a,1,0.0D0,0.0D0,1,1,0.0D0,m,r,z,1,iw,&
     &               q,26*n,iw(2*n+1),10*n-1,info)
         CALL chkxer('DSYEVR',INFot,NOUt,LERr,OK)
         nt = nt + 11
!
!        DSYEVR_2STAGE
!
         SRNamt = 'DSYEVR_2STAGE'
         n = 1
         INFot = 1
         CALL DSYEVR_2STAGE('/','A','U',0,a,1,0.0D0,0.0D0,1,1,0.0D0,m,r,&
     &                      z,1,iw,q,26*n,iw(2*n+1),10*n,info)
         CALL chkxer('DSYEVR_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 1
         CALL DSYEVR_2STAGE('V','A','U',0,a,1,0.0D0,0.0D0,1,1,0.0D0,m,r,&
     &                      z,1,iw,q,26*n,iw(2*n+1),10*n,info)
         CALL chkxer('DSYEVR_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DSYEVR_2STAGE('N','/','U',0,a,1,0.0D0,0.0D0,1,1,0.0D0,m,r,&
     &                      z,1,iw,q,26*n,iw(2*n+1),10*n,info)
         CALL chkxer('DSYEVR_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DSYEVR_2STAGE('N','A','/',-1,a,1,0.0D0,0.0D0,1,1,0.0D0,m, &
     &                      r,z,1,iw,q,26*n,iw(2*n+1),10*n,info)
         CALL chkxer('DSYEVR_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL DSYEVR_2STAGE('N','A','U',-1,a,1,0.0D0,0.0D0,1,1,0.0D0,m, &
     &                      r,z,1,iw,q,26*n,iw(2*n+1),10*n,info)
         CALL chkxer('DSYEVR_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL DSYEVR_2STAGE('N','A','U',2,a,1,0.0D0,0.0D0,1,1,0.0D0,m,r,&
     &                      z,1,iw,q,26*n,iw(2*n+1),10*n,info)
         CALL chkxer('DSYEVR_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL DSYEVR_2STAGE('N','V','U',1,a,1,0.0D0,0.0D0,1,1,0.0D0,m,r,&
     &                      z,1,iw,q,26*n,iw(2*n+1),10*n,info)
         CALL chkxer('DSYEVR_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL DSYEVR_2STAGE('N','I','U',1,a,1,0.0D0,0.0D0,0,1,0.0D0,m,r,&
     &                      z,1,iw,q,26*n,iw(2*n+1),10*n,info)
         CALL chkxer('DSYEVR_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL DSYEVR_2STAGE('N','I','U',2,a,2,0.0D0,0.0D0,2,1,0.0D0,m,r,&
     &                      z,1,iw,q,26*n,iw(2*n+1),10*n,info)
         CALL chkxer('DSYEVR_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 15
         CALL DSYEVR_2STAGE('N','I','U',1,a,1,0.0D0,0.0D0,1,1,0.0D0,m,r,&
     &                      z,0,iw,q,26*n,iw(2*n+1),10*n,info)
         CALL chkxer('DSYEVR_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 18
         CALL DSYEVR_2STAGE('N','I','U',1,a,1,0.0D0,0.0D0,1,1,0.0D0,m,r,&
     &                      z,1,iw,q,0,iw(2*n+1),10*n,info)
         CALL chkxer('DSYEVR_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 20
         CALL DSYEVR_2STAGE('N','I','U',1,a,1,0.0D0,0.0D0,1,1,0.0D0,m,r,&
     &                      z,1,iw,q,26*n,iw(2*n+1),0,info)
         CALL chkxer('DSYEVR_2STAGE',INFot,NOUt,LERr,OK)
         nt = nt + 12
!
!        DSYEV
!
         SRNamt = 'DSYEV '
         INFot = 1
         CALL DSYEV('/','U',0,a,1,x,w,1,info)
         CALL chkxer('DSYEV ',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DSYEV('N','/',0,a,1,x,w,1,info)
         CALL chkxer('DSYEV ',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DSYEV('N','U',-1,a,1,x,w,1,info)
         CALL chkxer('DSYEV ',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL DSYEV('N','U',2,a,1,x,w,3,info)
         CALL chkxer('DSYEV ',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL DSYEV('N','U',1,a,1,x,w,1,info)
         CALL chkxer('DSYEV ',INFot,NOUt,LERr,OK)
         nt = nt + 5
!
!        DSYEV_2STAGE
!
         SRNamt = 'DSYEV_2STAGE '
         INFot = 1
         CALL DSYEV_2STAGE('/','U',0,a,1,x,w,1,info)
         CALL chkxer('DSYEV_2STAGE ',INFot,NOUt,LERr,OK)
         INFot = 1
         CALL DSYEV_2STAGE('V','U',0,a,1,x,w,1,info)
         CALL chkxer('DSYEV_2STAGE ',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DSYEV_2STAGE('N','/',0,a,1,x,w,1,info)
         CALL chkxer('DSYEV_2STAGE ',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DSYEV_2STAGE('N','U',-1,a,1,x,w,1,info)
         CALL chkxer('DSYEV_2STAGE ',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL DSYEV_2STAGE('N','U',2,a,1,x,w,3,info)
         CALL chkxer('DSYEV_2STAGE ',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL DSYEV_2STAGE('N','U',1,a,1,x,w,1,info)
         CALL chkxer('DSYEV_2STAGE ',INFot,NOUt,LERr,OK)
         nt = nt + 6
!
!        DSYEVX
!
         SRNamt = 'DSYEVX'
         INFot = 1
         CALL DSYEVX('/','A','U',0,a,1,0.0D0,0.0D0,0,0,0.0D0,m,x,z,1,w, &
     &               1,iw,i3,info)
         CALL chkxer('DSYEVX',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DSYEVX('N','/','U',0,a,1,0.0D0,1.0D0,1,0,0.0D0,m,x,z,1,w, &
     &               1,iw,i3,info)
         CALL chkxer('DSYEVX',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DSYEVX('N','A','/',0,a,1,0.0D0,0.0D0,0,0,0.0D0,m,x,z,1,w, &
     &               1,iw,i3,info)
         INFot = 4
         CALL DSYEVX('N','A','U',-1,a,1,0.0D0,0.0D0,0,0,0.0D0,m,x,z,1,w,&
     &               1,iw,i3,info)
         CALL chkxer('DSYEVX',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL DSYEVX('N','A','U',2,a,1,0.0D0,0.0D0,0,0,0.0D0,m,x,z,1,w, &
     &               16,iw,i3,info)
         CALL chkxer('DSYEVX',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL DSYEVX('N','V','U',1,a,1,0.0D0,0.0D0,0,0,0.0D0,m,x,z,1,w, &
     &               8,iw,i3,info)
         CALL chkxer('DSYEVX',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL DSYEVX('N','I','U',1,a,1,0.0D0,0.0D0,0,0,0.0D0,m,x,z,1,w, &
     &               8,iw,i3,info)
         CALL chkxer('DSYEVX',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL DSYEVX('N','I','U',1,a,1,0.0D0,0.0D0,2,1,0.0D0,m,x,z,1,w, &
     &               8,iw,i3,info)
         CALL chkxer('DSYEVX',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL DSYEVX('N','I','U',2,a,2,0.0D0,0.0D0,2,1,0.0D0,m,x,z,1,w, &
     &               16,iw,i3,info)
         CALL chkxer('DSYEVX',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL DSYEVX('N','I','U',1,a,1,0.0D0,0.0D0,1,2,0.0D0,m,x,z,1,w, &
     &               8,iw,i3,info)
         CALL chkxer('DSYEVX',INFot,NOUt,LERr,OK)
         INFot = 15
         CALL DSYEVX('V','A','U',2,a,2,0.0D0,0.0D0,0,0,0.0D0,m,x,z,1,w, &
     &               16,iw,i3,info)
         CALL chkxer('DSYEVX',INFot,NOUt,LERr,OK)
         INFot = 17
         CALL DSYEVX('V','A','U',1,a,1,0.0D0,0.0D0,0,0,0.0D0,m,x,z,1,w, &
     &               0,iw,i3,info)
         CALL chkxer('DSYEVX',INFot,NOUt,LERr,OK)
         nt = nt + 12
!
!        DSYEVX_2STAGE
!
         SRNamt = 'DSYEVX_2STAGE'
         INFot = 1
         CALL DSYEVX_2STAGE('/','A','U',0,a,1,0.0D0,0.0D0,0,0,0.0D0,m,x,&
     &                      z,1,w,1,iw,i3,info)
         CALL chkxer('DSYEVX_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 1
         CALL DSYEVX_2STAGE('V','A','U',0,a,1,0.0D0,0.0D0,0,0,0.0D0,m,x,&
     &                      z,1,w,1,iw,i3,info)
         CALL chkxer('DSYEVX_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DSYEVX_2STAGE('N','/','U',0,a,1,0.0D0,1.0D0,1,0,0.0D0,m,x,&
     &                      z,1,w,1,iw,i3,info)
         CALL chkxer('DSYEVX_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DSYEVX_2STAGE('N','A','/',0,a,1,0.0D0,0.0D0,0,0,0.0D0,m,x,&
     &                      z,1,w,1,iw,i3,info)
         INFot = 4
         CALL DSYEVX_2STAGE('N','A','U',-1,a,1,0.0D0,0.0D0,0,0,0.0D0,m, &
     &                      x,z,1,w,1,iw,i3,info)
         CALL chkxer('DSYEVX_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL DSYEVX_2STAGE('N','A','U',2,a,1,0.0D0,0.0D0,0,0,0.0D0,m,x,&
     &                      z,1,w,16,iw,i3,info)
         CALL chkxer('DSYEVX_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL DSYEVX_2STAGE('N','V','U',1,a,1,0.0D0,0.0D0,0,0,0.0D0,m,x,&
     &                      z,1,w,8,iw,i3,info)
         CALL chkxer('DSYEVX_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL DSYEVX_2STAGE('N','I','U',1,a,1,0.0D0,0.0D0,0,0,0.0D0,m,x,&
     &                      z,1,w,8,iw,i3,info)
         CALL chkxer('DSYEVX_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL DSYEVX_2STAGE('N','I','U',1,a,1,0.0D0,0.0D0,2,1,0.0D0,m,x,&
     &                      z,1,w,8,iw,i3,info)
         CALL chkxer('DSYEVX_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL DSYEVX_2STAGE('N','I','U',2,a,2,0.0D0,0.0D0,2,1,0.0D0,m,x,&
     &                      z,1,w,16,iw,i3,info)
         CALL chkxer('DSYEVX_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL DSYEVX_2STAGE('N','I','U',1,a,1,0.0D0,0.0D0,1,2,0.0D0,m,x,&
     &                      z,1,w,8,iw,i3,info)
         CALL chkxer('DSYEVX_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 15
         CALL DSYEVX_2STAGE('N','A','U',2,a,2,0.0D0,0.0D0,0,0,0.0D0,m,x,&
     &                      z,0,w,16,iw,i3,info)
         CALL chkxer('DSYEVX_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 17
         CALL DSYEVX_2STAGE('N','A','U',1,a,1,0.0D0,0.0D0,0,0,0.0D0,m,x,&
     &                      z,1,w,0,iw,i3,info)
         CALL chkxer('DSYEVX_2STAGE',INFot,NOUt,LERr,OK)
         nt = nt + 13
!
!        DSPEVD
!
         SRNamt = 'DSPEVD'
         INFot = 1
         CALL DSPEVD('/','U',0,a,x,z,1,w,1,iw,1,info)
         CALL chkxer('DSPEVD',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DSPEVD('N','/',0,a,x,z,1,w,1,iw,1,info)
         CALL chkxer('DSPEVD',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DSPEVD('N','U',-1,a,x,z,1,w,1,iw,1,info)
         CALL chkxer('DSPEVD',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL DSPEVD('V','U',2,a,x,z,1,w,23,iw,12,info)
         CALL chkxer('DSPEVD',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL DSPEVD('N','U',1,a,x,z,1,w,0,iw,1,info)
         CALL chkxer('DSPEVD',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL DSPEVD('N','U',2,a,x,z,1,w,3,iw,1,info)
         CALL chkxer('DSPEVD',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL DSPEVD('V','U',2,a,x,z,2,w,16,iw,12,info)
         CALL chkxer('DSPEVD',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL DSPEVD('N','U',1,a,x,z,1,w,1,iw,0,info)
         CALL chkxer('DSPEVD',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL DSPEVD('N','U',2,a,x,z,1,w,4,iw,0,info)
         CALL chkxer('DSPEVD',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL DSPEVD('V','U',2,a,x,z,2,w,23,iw,11,info)
         CALL chkxer('DSPEVD',INFot,NOUt,LERr,OK)
         nt = nt + 10
!
!        DSPEV
!
         SRNamt = 'DSPEV '
         INFot = 1
         CALL DSPEV('/','U',0,a,w,z,1,x,info)
         CALL chkxer('DSPEV ',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DSPEV('N','/',0,a,w,z,1,x,info)
         CALL chkxer('DSPEV ',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DSPEV('N','U',-1,a,w,z,1,x,info)
         CALL chkxer('DSPEV ',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL DSPEV('V','U',2,a,w,z,1,x,info)
         CALL chkxer('DSPEV ',INFot,NOUt,LERr,OK)
         nt = nt + 4
!
!        DSPEVX
!
         SRNamt = 'DSPEVX'
         INFot = 1
         CALL DSPEVX('/','A','U',0,a,0.0D0,0.0D0,0,0,0.0D0,m,x,z,1,w,iw,&
     &               i3,info)
         CALL chkxer('DSPEVX',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DSPEVX('N','/','U',0,a,0.0D0,0.0D0,0,0,0.0D0,m,x,z,1,w,iw,&
     &               i3,info)
         CALL chkxer('DSPEVX',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DSPEVX('N','A','/',0,a,0.0D0,0.0D0,0,0,0.0D0,m,x,z,1,w,iw,&
     &               i3,info)
         INFot = 4
         CALL DSPEVX('N','A','U',-1,a,0.0D0,0.0D0,0,0,0.0D0,m,x,z,1,w,  &
     &               iw,i3,info)
         CALL chkxer('DSPEVX',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL DSPEVX('N','V','U',1,a,0.0D0,0.0D0,0,0,0.0D0,m,x,z,1,w,iw,&
     &               i3,info)
         CALL chkxer('DSPEVX',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL DSPEVX('N','I','U',1,a,0.0D0,0.0D0,0,0,0.0D0,m,x,z,1,w,iw,&
     &               i3,info)
         CALL chkxer('DSPEVX',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL DSPEVX('N','I','U',1,a,0.0D0,0.0D0,2,1,0.0D0,m,x,z,1,w,iw,&
     &               i3,info)
         CALL chkxer('DSPEVX',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL DSPEVX('N','I','U',2,a,0.0D0,0.0D0,2,1,0.0D0,m,x,z,1,w,iw,&
     &               i3,info)
         CALL chkxer('DSPEVX',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL DSPEVX('N','I','U',1,a,0.0D0,0.0D0,1,2,0.0D0,m,x,z,1,w,iw,&
     &               i3,info)
         CALL chkxer('DSPEVX',INFot,NOUt,LERr,OK)
         INFot = 14
         CALL DSPEVX('V','A','U',2,a,0.0D0,0.0D0,0,0,0.0D0,m,x,z,1,w,iw,&
     &               i3,info)
         CALL chkxer('DSPEVX',INFot,NOUt,LERr,OK)
         nt = nt + 10
!
!     Test error exits for the SB path.
!
      ELSEIF ( LSAMEN(2,c2,'SB') ) THEN
!
!        DSBTRD
!
         SRNamt = 'DSBTRD'
         INFot = 1
         CALL DSBTRD('/','U',0,0,a,1,d,e,z,1,w,info)
         CALL chkxer('DSBTRD',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DSBTRD('N','/',0,0,a,1,d,e,z,1,w,info)
         CALL chkxer('DSBTRD',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DSBTRD('N','U',-1,0,a,1,d,e,z,1,w,info)
         CALL chkxer('DSBTRD',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL DSBTRD('N','U',0,-1,a,1,d,e,z,1,w,info)
         CALL chkxer('DSBTRD',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL DSBTRD('N','U',1,1,a,1,d,e,z,1,w,info)
         CALL chkxer('DSBTRD',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL DSBTRD('V','U',2,0,a,1,d,e,z,1,w,info)
         CALL chkxer('DSBTRD',INFot,NOUt,LERr,OK)
         nt = nt + 6
!
!        DSYTRD_SB2ST
!
         SRNamt = 'DSYTRD_SB2ST'
         INFot = 1
         CALL DSYTRD_SB2ST('/','N','U',0,0,a,1,d,e,c,1,w,1,info)
         CALL chkxer('DSYTRD_SB2ST',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DSYTRD_SB2ST('N','/','U',0,0,a,1,d,e,c,1,w,1,info)
         CALL chkxer('DSYTRD_SB2ST',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DSYTRD_SB2ST('N','H','U',0,0,a,1,d,e,c,1,w,1,info)
         CALL chkxer('DSYTRD_SB2ST',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DSYTRD_SB2ST('N','N','/',0,0,a,1,d,e,c,1,w,1,info)
         CALL chkxer('DSYTRD_SB2ST',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL DSYTRD_SB2ST('N','N','U',-1,0,a,1,d,e,c,1,w,1,info)
         CALL chkxer('DSYTRD_SB2ST',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL DSYTRD_SB2ST('N','N','U',0,-1,a,1,d,e,c,1,w,1,info)
         CALL chkxer('DSYTRD_SB2ST',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL DSYTRD_SB2ST('N','N','U',0,1,a,1,d,e,c,1,w,1,info)
         CALL chkxer('DSYTRD_SB2ST',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL DSYTRD_SB2ST('N','N','U',0,0,a,1,d,e,c,0,w,1,info)
         CALL chkxer('DSYTRD_SB2ST',INFot,NOUt,LERr,OK)
         INFot = 13
         CALL DSYTRD_SB2ST('N','N','U',0,0,a,1,d,e,c,1,w,0,info)
         CALL chkxer('DSYTRD_SB2ST',INFot,NOUt,LERr,OK)
         nt = nt + 9
!
!        DSBEVD
!
         SRNamt = 'DSBEVD'
         INFot = 1
         CALL DSBEVD('/','U',0,0,a,1,x,z,1,w,1,iw,1,info)
         CALL chkxer('DSBEVD',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DSBEVD('N','/',0,0,a,1,x,z,1,w,1,iw,1,info)
         CALL chkxer('DSBEVD',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DSBEVD('N','U',-1,0,a,1,x,z,1,w,1,iw,1,info)
         CALL chkxer('DSBEVD',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL DSBEVD('N','U',0,-1,a,1,x,z,1,w,1,iw,1,info)
         CALL chkxer('DSBEVD',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL DSBEVD('N','U',2,1,a,1,x,z,1,w,4,iw,1,info)
         CALL chkxer('DSBEVD',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL DSBEVD('V','U',2,1,a,2,x,z,1,w,25,iw,12,info)
         CALL chkxer('DSBEVD',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL DSBEVD('N','U',1,0,a,1,x,z,1,w,0,iw,1,info)
         CALL chkxer('DSBEVD',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL DSBEVD('N','U',2,0,a,1,x,z,1,w,3,iw,1,info)
         CALL chkxer('DSBEVD',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL DSBEVD('V','U',2,0,a,1,x,z,2,w,18,iw,12,info)
         CALL chkxer('DSBEVD',INFot,NOUt,LERr,OK)
         INFot = 13
         CALL DSBEVD('N','U',1,0,a,1,x,z,1,w,1,iw,0,info)
         CALL chkxer('DSBEVD',INFot,NOUt,LERr,OK)
         INFot = 13
         CALL DSBEVD('V','U',2,0,a,1,x,z,2,w,25,iw,11,info)
         CALL chkxer('DSBEVD',INFot,NOUt,LERr,OK)
         nt = nt + 11
!
!        DSBEVD_2STAGE
!
         SRNamt = 'DSBEVD_2STAGE'
         INFot = 1
         CALL DSBEVD_2STAGE('/','U',0,0,a,1,x,z,1,w,1,iw,1,info)
         CALL chkxer('DSBEVD_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 1
         CALL DSBEVD_2STAGE('V','U',0,0,a,1,x,z,1,w,1,iw,1,info)
         CALL chkxer('DSBEVD_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DSBEVD_2STAGE('N','/',0,0,a,1,x,z,1,w,1,iw,1,info)
         CALL chkxer('DSBEVD_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DSBEVD_2STAGE('N','U',-1,0,a,1,x,z,1,w,1,iw,1,info)
         CALL chkxer('DSBEVD_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL DSBEVD_2STAGE('N','U',0,-1,a,1,x,z,1,w,1,iw,1,info)
         CALL chkxer('DSBEVD_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL DSBEVD_2STAGE('N','U',2,1,a,1,x,z,1,w,4,iw,1,info)
         CALL chkxer('DSBEVD_2STAGE',INFot,NOUt,LERr,OK)
!         INFOT = 9
!         CALL DSBEVD_2STAGE( 'V', 'U', 2, 1, A, 2, X, Z, 1, W,
!     $                                      25, IW, 12, INFO )
!         CALL chkxer( 'DSBEVD_2STAGE', INFOT, NOUT, LERR, OK )
         INFot = 11
         CALL DSBEVD_2STAGE('N','U',1,0,a,1,x,z,1,w,0,iw,1,info)
         CALL chkxer('DSBEVD_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL DSBEVD_2STAGE('N','U',2,0,a,1,x,z,1,w,3,iw,1,info)
         CALL chkxer('DSBEVD_2STAGE',INFot,NOUt,LERr,OK)
!         INFOT = 11
!         CALL DSBEVD_2STAGE( 'V', 'U', 2, 0, A, 1, X, Z, 2, W,
!     $                                      18, IW, 12, INFO )
!         CALL chkxer( 'DSBEVD_2STAGE', INFOT, NOUT, LERR, OK )
         INFot = 13
         CALL DSBEVD_2STAGE('N','U',1,0,a,1,x,z,1,w,1,iw,0,info)
         CALL chkxer('DSBEVD_2STAGE',INFot,NOUt,LERr,OK)
!         INFOT = 13
!         CALL DSBEVD_2STAGE( 'V', 'U', 2, 0, A, 1, X, Z, 2, W,
!     $                                      25, IW, 11, INFO )
!         CALL chkxer( 'DSBEVD_2STAGE', INFOT, NOUT, LERR, OK )
!         NT = NT + 12
         nt = nt + 9
!
!        DSBEV
!
         SRNamt = 'DSBEV '
         INFot = 1
         CALL DSBEV('/','U',0,0,a,1,x,z,1,w,info)
         CALL chkxer('DSBEV ',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DSBEV('N','/',0,0,a,1,x,z,1,w,info)
         CALL chkxer('DSBEV ',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DSBEV('N','U',-1,0,a,1,x,z,1,w,info)
         CALL chkxer('DSBEV ',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL DSBEV('N','U',0,-1,a,1,x,z,1,w,info)
         CALL chkxer('DSBEV ',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL DSBEV('N','U',2,1,a,1,x,z,1,w,info)
         CALL chkxer('DSBEV ',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL DSBEV('V','U',2,0,a,1,x,z,1,w,info)
         CALL chkxer('DSBEV ',INFot,NOUt,LERr,OK)
         nt = nt + 6
!
!        DSBEV_2STAGE
!
         SRNamt = 'DSBEV_2STAGE '
         INFot = 1
         CALL DSBEV_2STAGE('/','U',0,0,a,1,x,z,1,w,0,info)
         CALL chkxer('DSBEV_2STAGE ',INFot,NOUt,LERr,OK)
         INFot = 1
         CALL DSBEV_2STAGE('V','U',0,0,a,1,x,z,1,w,0,info)
         CALL chkxer('DSBEV_2STAGE ',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DSBEV_2STAGE('N','/',0,0,a,1,x,z,1,w,0,info)
         CALL chkxer('DSBEV_2STAGE ',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DSBEV_2STAGE('N','U',-1,0,a,1,x,z,1,w,0,info)
         CALL chkxer('DSBEV_2STAGE ',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL DSBEV_2STAGE('N','U',0,-1,a,1,x,z,1,w,0,info)
         CALL chkxer('DSBEV_2STAGE ',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL DSBEV_2STAGE('N','U',2,1,a,1,x,z,1,w,0,info)
         CALL chkxer('DSBEV_2STAGE ',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL DSBEV_2STAGE('N','U',2,0,a,1,x,z,0,w,0,info)
         CALL chkxer('DSBEV_2STAGE ',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL DSBEV_2STAGE('N','U',0,0,a,1,x,z,1,w,0,info)
         CALL chkxer('DSBEV_2STAGE ',INFot,NOUt,LERr,OK)
         nt = nt + 8
!
!        DSBEVX
!
         SRNamt = 'DSBEVX'
         INFot = 1
         CALL DSBEVX('/','A','U',0,0,a,1,q,1,0.0D0,0.0D0,0,0,0.0D0,m,x, &
     &               z,1,w,iw,i3,info)
         CALL chkxer('DSBEVX',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DSBEVX('N','/','U',0,0,a,1,q,1,0.0D0,0.0D0,0,0,0.0D0,m,x, &
     &               z,1,w,iw,i3,info)
         CALL chkxer('DSBEVX',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DSBEVX('N','A','/',0,0,a,1,q,1,0.0D0,0.0D0,0,0,0.0D0,m,x, &
     &               z,1,w,iw,i3,info)
         CALL chkxer('DSBEVX',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL DSBEVX('N','A','U',-1,0,a,1,q,1,0.0D0,0.0D0,0,0,0.0D0,m,x,&
     &               z,1,w,iw,i3,info)
         CALL chkxer('DSBEVX',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL DSBEVX('N','A','U',0,-1,a,1,q,1,0.0D0,0.0D0,0,0,0.0D0,m,x,&
     &               z,1,w,iw,i3,info)
         CALL chkxer('DSBEVX',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL DSBEVX('N','A','U',2,1,a,1,q,1,0.0D0,0.0D0,0,0,0.0D0,m,x, &
     &               z,1,w,iw,i3,info)
         CALL chkxer('DSBEVX',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL DSBEVX('V','A','U',2,0,a,1,q,1,0.0D0,0.0D0,0,0,0.0D0,m,x, &
     &               z,2,w,iw,i3,info)
         CALL chkxer('DSBEVX',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL DSBEVX('N','V','U',1,0,a,1,q,1,0.0D0,0.0D0,0,0,0.0D0,m,x, &
     &               z,1,w,iw,i3,info)
         CALL chkxer('DSBEVX',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL DSBEVX('N','I','U',1,0,a,1,q,1,0.0D0,0.0D0,0,0,0.0D0,m,x, &
     &               z,1,w,iw,i3,info)
         CALL chkxer('DSBEVX',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL DSBEVX('N','I','U',1,0,a,1,q,1,0.0D0,0.0D0,2,1,0.0D0,m,x, &
     &               z,1,w,iw,i3,info)
         CALL chkxer('DSBEVX',INFot,NOUt,LERr,OK)
         INFot = 13
         CALL DSBEVX('N','I','U',2,0,a,1,q,1,0.0D0,0.0D0,2,1,0.0D0,m,x, &
     &               z,1,w,iw,i3,info)
         CALL chkxer('DSBEVX',INFot,NOUt,LERr,OK)
         INFot = 13
         CALL DSBEVX('N','I','U',1,0,a,1,q,1,0.0D0,0.0D0,1,2,0.0D0,m,x, &
     &               z,1,w,iw,i3,info)
         CALL chkxer('DSBEVX',INFot,NOUt,LERr,OK)
         INFot = 18
         CALL DSBEVX('V','A','U',2,0,a,1,q,2,0.0D0,0.0D0,0,0,0.0D0,m,x, &
     &               z,1,w,iw,i3,info)
         CALL chkxer('DSBEVX',INFot,NOUt,LERr,OK)
         nt = nt + 13
!
!        DSBEVX_2STAGE
!
         SRNamt = 'DSBEVX_2STAGE'
         INFot = 1
         CALL DSBEVX_2STAGE('/','A','U',0,0,a,1,q,1,0.0D0,0.0D0,0,0,    &
     &                      0.0D0,m,x,z,1,w,0,iw,i3,info)
         CALL chkxer('DSBEVX_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 1
         CALL DSBEVX_2STAGE('V','A','U',0,0,a,1,q,1,0.0D0,0.0D0,0,0,    &
     &                      0.0D0,m,x,z,1,w,0,iw,i3,info)
         CALL chkxer('DSBEVX_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DSBEVX_2STAGE('N','/','U',0,0,a,1,q,1,0.0D0,0.0D0,0,0,    &
     &                      0.0D0,m,x,z,1,w,0,iw,i3,info)
         CALL chkxer('DSBEVX_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DSBEVX_2STAGE('N','A','/',0,0,a,1,q,1,0.0D0,0.0D0,0,0,    &
     &                      0.0D0,m,x,z,1,w,0,iw,i3,info)
         CALL chkxer('DSBEVX_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL DSBEVX_2STAGE('N','A','U',-1,0,a,1,q,1,0.0D0,0.0D0,0,0,   &
     &                      0.0D0,m,x,z,1,w,0,iw,i3,info)
         CALL chkxer('DSBEVX_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL DSBEVX_2STAGE('N','A','U',0,-1,a,1,q,1,0.0D0,0.0D0,0,0,   &
     &                      0.0D0,m,x,z,1,w,0,iw,i3,info)
         CALL chkxer('DSBEVX_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL DSBEVX_2STAGE('N','A','U',2,1,a,1,q,1,0.0D0,0.0D0,0,0,    &
     &                      0.0D0,m,x,z,1,w,0,iw,i3,info)
         CALL chkxer('DSBEVX_2STAGE',INFot,NOUt,LERr,OK)
!         INFOT = 9
!         CALL DSBEVX_2STAGE( 'V', 'A', 'U', 2, 0, A, 1, Q, 1, 0.0D0,
!     $          0.0D0, 0, 0, 0.0D0, M, X, Z, 2, W, 0, IW, I3, INFO )
!         CALL chkxer( 'DSBEVX_2STAGE', INFOT, NOUT, LERR, OK )
         INFot = 11
         CALL DSBEVX_2STAGE('N','V','U',1,0,a,1,q,1,0.0D0,0.0D0,0,0,    &
     &                      0.0D0,m,x,z,1,w,0,iw,i3,info)
         CALL chkxer('DSBEVX_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL DSBEVX_2STAGE('N','I','U',1,0,a,1,q,1,0.0D0,0.0D0,0,0,    &
     &                      0.0D0,m,x,z,1,w,0,iw,i3,info)
         CALL chkxer('DSBEVX_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL DSBEVX_2STAGE('N','I','U',1,0,a,1,q,1,0.0D0,0.0D0,2,1,    &
     &                      0.0D0,m,x,z,1,w,0,iw,i3,info)
         CALL chkxer('DSBEVX_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 13
         CALL DSBEVX_2STAGE('N','I','U',2,0,a,1,q,1,0.0D0,0.0D0,2,1,    &
     &                      0.0D0,m,x,z,1,w,0,iw,i3,info)
         CALL chkxer('DSBEVX_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 13
         CALL DSBEVX_2STAGE('N','I','U',1,0,a,1,q,1,0.0D0,0.0D0,1,2,    &
     &                      0.0D0,m,x,z,1,w,0,iw,i3,info)
         CALL chkxer('DSBEVX_2STAGE',INFot,NOUt,LERr,OK)
!         INFOT = 18
!         CALL DSBEVX_2STAGE( 'V', 'A', 'U', 2, 0, A, 1, Q, 2, 0.0D0,
!     $          0.0D0, 0, 0, 0.0D0, M, X, Z, 1, W, 0, IW, I3, INFO )
!         CALL chkxer( 'DSBEVX_2STAGE', INFOT, NOUT, LERR, OK )
         INFot = 20
         CALL DSBEVX_2STAGE('N','A','U',0,0,a,1,q,1,0.0D0,0.0D0,0,0,    &
     &                      0.0D0,m,x,z,1,w,0,iw,i3,info)
         CALL chkxer('DSBEVX_2STAGE',INFot,NOUt,LERr,OK)
!         NT = NT + 15
         nt = nt + 13
      ENDIF
!
!     Print a summary line.
!
      IF ( OK ) THEN
         WRITE (NOUt,FMT=99001) Path , nt
      ELSE
         WRITE (NOUt,FMT=99002) Path
      ENDIF
!
99001 FORMAT (1X,A3,' routines passed the tests of the error exits',    &
     &        ' (',I3,' tests done)')
99002 FORMAT (' *** ',A3,' routines failed the tests of the error ',    &
     &        'exits ***')
!
!
!     End of derrst
!
      END SUBROUTINE derrst
