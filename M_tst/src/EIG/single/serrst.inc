!*==serrst.f90  processed by SPAG 7.51RB at 17:37 on  4 Mar 2022
!> \brief \b serrst
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE serrst( PATH, NUNIT )
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
!> serrst tests the error exits for SSYTRD, SORGTR, SORMTR, SSPTRD,
!> SOPGTR, SOPMTR, SSTEQR, SSTERF, SSTEBZ, SSTEIN, SPTEQR, SSBTRD,
!> SSYEV, SSYEVX, SSYEVD, SSBEV, SSBEVX, SSBEVD,
!> SSPEV, SSPEVX, SSPEVD, SSTEV, SSTEVX, SSTEVD, and SSTEDC.
!> SSYEVD_2STAGE, SSYEVR_2STAGE, SSYEVX_2STAGE,
!> SSYEV_2STAGE, SSBEV_2STAGE, SSBEVD_2STAGE,
!> SSBEVX_2STAGE, SSYTRD_2STAGE, SSYTRD_SY2SB,
!> SSYTRD_SB2ST
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
!> \ingroup single_eig
!
!  =====================================================================
      SUBROUTINE serrst(Path,Nunit)
      IMPLICIT NONE
!*--SERRST65
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
      REAL a(NMAX,NMAX) , c(NMAX,NMAX) , d(NMAX) , e(NMAX) ,            &
     &     q(NMAX,NMAX) , r(NMAX) , tau(NMAX) , w(LW) , x(NMAX) ,       &
     &     z(NMAX,NMAX)
!     ..
!     .. External Functions ..
      LOGICAL LSAMEN
      EXTERNAL LSAMEN
!     ..
!     .. External Subroutines ..
      EXTERNAL SOPGTR , SOPMTR , SORGTR , SORMTR , SPTEQR 
     EXTERNAL         SSBEV , SSBEVD , SSBEVX , SSBTRD , SSPEV , SSPEVD 
     EXTERNAL         SSPEVX , SSPTRD , SSTEBZ , SSTEDC , SSTEIN , SSTEQR 
     EXTERNAL         SSTERF , SSTEV , SSTEVD , SSTEVR , SSTEVX , SSYEV 
     EXTERNAL         SSYEVD , SSYEVR , SSYEVX , SSYTRD , SSYEVD_2STAGE 
     EXTERNAL         SSYEVR_2STAGE , SSYEVX_2STAGE , SSYEV_2STAGE 
     EXTERNAL         SSBEV_2STAGE , SSBEVD_2STAGE , SSBEVX_2STAGE 
     EXTERNAL         SSYTRD_2STAGE , SSYTRD_SY2SB , SSYTRD_SB2ST
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
         ENDDO
      ENDDO
      DO j = 1 , NMAX
         d(j) = REAL(j)
         e(j) = 0.0
         i1(j) = j
         i2(j) = j
         tau(j) = 1.
      ENDDO
      OK = .TRUE.
      nt = 0
!
!     Test error exits for the ST path.
!
      IF ( LSAMEN(2,c2,'ST') ) THEN
!
!        SSYTRD
!
         SRNamt = 'SSYTRD'
         INFot = 1
         CALL SSYTRD('/',0,a,1,d,e,tau,w,1,info)
         CALL chkxer('SSYTRD',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SSYTRD('U',-1,a,1,d,e,tau,w,1,info)
         CALL chkxer('SSYTRD',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL SSYTRD('U',2,a,1,d,e,tau,w,1,info)
         CALL chkxer('SSYTRD',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL SSYTRD('U',0,a,1,d,e,tau,w,0,info)
         CALL chkxer('SSYTRD',INFot,NOUt,LERr,OK)
         nt = nt + 4
!
!        SSYTRD_2STAGE
!
         SRNamt = 'SSYTRD_2STAGE'
         INFot = 1
         CALL SSYTRD_2STAGE('/','U',0,a,1,d,e,tau,c,1,w,1,info)
         CALL chkxer('SSYTRD_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 1
         CALL SSYTRD_2STAGE('H','U',0,a,1,d,e,tau,c,1,w,1,info)
         CALL chkxer('SSYTRD_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SSYTRD_2STAGE('N','/',0,a,1,d,e,tau,c,1,w,1,info)
         CALL chkxer('SSYTRD_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SSYTRD_2STAGE('N','U',-1,a,1,d,e,tau,c,1,w,1,info)
         CALL chkxer('SSYTRD_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL SSYTRD_2STAGE('N','U',2,a,1,d,e,tau,c,1,w,1,info)
         CALL chkxer('SSYTRD_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL SSYTRD_2STAGE('N','U',0,a,1,d,e,tau,c,0,w,1,info)
         CALL chkxer('SSYTRD_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL SSYTRD_2STAGE('N','U',0,a,1,d,e,tau,c,1,w,0,info)
         CALL chkxer('SSYTRD_2STAGE',INFot,NOUt,LERr,OK)
         nt = nt + 7
!
!        SSYTRD_SY2SB
!
         SRNamt = 'SSYTRD_SY2SB'
         INFot = 1
         CALL SSYTRD_SY2SB('/',0,0,a,1,c,1,tau,w,1,info)
         CALL chkxer('SSYTRD_SY2SB',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SSYTRD_SY2SB('U',-1,0,a,1,c,1,tau,w,1,info)
         CALL chkxer('SSYTRD_SY2SB',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SSYTRD_SY2SB('U',0,-1,a,1,c,1,tau,w,1,info)
         CALL chkxer('SSYTRD_SY2SB',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL SSYTRD_SY2SB('U',2,0,a,1,c,1,tau,w,1,info)
         CALL chkxer('SSYTRD_SY2SB',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL SSYTRD_SY2SB('U',0,2,a,1,c,1,tau,w,1,info)
         CALL chkxer('SSYTRD_SY2SB',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL SSYTRD_SY2SB('U',0,0,a,1,c,1,tau,w,0,info)
         CALL chkxer('SSYTRD_SY2SB',INFot,NOUt,LERr,OK)
         nt = nt + 6
!
!        SSYTRD_SB2ST
!
         SRNamt = 'SSYTRD_SB2ST'
         INFot = 1
         CALL SSYTRD_SB2ST('/','N','U',0,0,a,1,d,e,c,1,w,1,info)
         CALL chkxer('SSYTRD_SB2ST',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SSYTRD_SB2ST('Y','/','U',0,0,a,1,d,e,c,1,w,1,info)
         CALL chkxer('SSYTRD_SB2ST',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SSYTRD_SB2ST('Y','H','U',0,0,a,1,d,e,c,1,w,1,info)
         CALL chkxer('SSYTRD_SB2ST',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SSYTRD_SB2ST('Y','N','/',0,0,a,1,d,e,c,1,w,1,info)
         CALL chkxer('SSYTRD_SB2ST',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL SSYTRD_SB2ST('Y','N','U',-1,0,a,1,d,e,c,1,w,1,info)
         CALL chkxer('SSYTRD_SB2ST',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL SSYTRD_SB2ST('Y','N','U',0,-1,a,1,d,e,c,1,w,1,info)
         CALL chkxer('SSYTRD_SB2ST',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL SSYTRD_SB2ST('Y','N','U',0,1,a,1,d,e,c,1,w,1,info)
         CALL chkxer('SSYTRD_SB2ST',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL SSYTRD_SB2ST('Y','N','U',0,0,a,1,d,e,c,0,w,1,info)
         CALL chkxer('SSYTRD_SB2ST',INFot,NOUt,LERr,OK)
         INFot = 13
         CALL SSYTRD_SB2ST('Y','N','U',0,0,a,1,d,e,c,1,w,0,info)
         CALL chkxer('SSYTRD_SB2ST',INFot,NOUt,LERr,OK)
         nt = nt + 9
!
!        SORGTR
!
         SRNamt = 'SORGTR'
         INFot = 1
         CALL SORGTR('/',0,a,1,tau,w,1,info)
         CALL chkxer('SORGTR',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SORGTR('U',-1,a,1,tau,w,1,info)
         CALL chkxer('SORGTR',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL SORGTR('U',2,a,1,tau,w,1,info)
         CALL chkxer('SORGTR',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL SORGTR('U',3,a,3,tau,w,1,info)
         CALL chkxer('SORGTR',INFot,NOUt,LERr,OK)
         nt = nt + 4
!
!        SORMTR
!
         SRNamt = 'SORMTR'
         INFot = 1
         CALL SORMTR('/','U','N',0,0,a,1,tau,c,1,w,1,info)
         CALL chkxer('SORMTR',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SORMTR('L','/','N',0,0,a,1,tau,c,1,w,1,info)
         CALL chkxer('SORMTR',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SORMTR('L','U','/',0,0,a,1,tau,c,1,w,1,info)
         CALL chkxer('SORMTR',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL SORMTR('L','U','N',-1,0,a,1,tau,c,1,w,1,info)
         CALL chkxer('SORMTR',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL SORMTR('L','U','N',0,-1,a,1,tau,c,1,w,1,info)
         CALL chkxer('SORMTR',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL SORMTR('L','U','N',2,0,a,1,tau,c,2,w,1,info)
         CALL chkxer('SORMTR',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL SORMTR('R','U','N',0,2,a,1,tau,c,1,w,1,info)
         CALL chkxer('SORMTR',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL SORMTR('L','U','N',2,0,a,2,tau,c,1,w,1,info)
         CALL chkxer('SORMTR',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL SORMTR('L','U','N',0,2,a,1,tau,c,1,w,1,info)
         CALL chkxer('SORMTR',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL SORMTR('R','U','N',2,0,a,1,tau,c,2,w,1,info)
         CALL chkxer('SORMTR',INFot,NOUt,LERr,OK)
         nt = nt + 10
!
!        SSPTRD
!
         SRNamt = 'SSPTRD'
         INFot = 1
         CALL SSPTRD('/',0,a,d,e,tau,info)
         CALL chkxer('SSPTRD',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SSPTRD('U',-1,a,d,e,tau,info)
         CALL chkxer('SSPTRD',INFot,NOUt,LERr,OK)
         nt = nt + 2
!
!        SOPGTR
!
         SRNamt = 'SOPGTR'
         INFot = 1
         CALL SOPGTR('/',0,a,tau,z,1,w,info)
         CALL chkxer('SOPGTR',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SOPGTR('U',-1,a,tau,z,1,w,info)
         CALL chkxer('SOPGTR',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL SOPGTR('U',2,a,tau,z,1,w,info)
         CALL chkxer('SOPGTR',INFot,NOUt,LERr,OK)
         nt = nt + 3
!
!        SOPMTR
!
         SRNamt = 'SOPMTR'
         INFot = 1
         CALL SOPMTR('/','U','N',0,0,a,tau,c,1,w,info)
         CALL chkxer('SOPMTR',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SOPMTR('L','/','N',0,0,a,tau,c,1,w,info)
         CALL chkxer('SOPMTR',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SOPMTR('L','U','/',0,0,a,tau,c,1,w,info)
         CALL chkxer('SOPMTR',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL SOPMTR('L','U','N',-1,0,a,tau,c,1,w,info)
         CALL chkxer('SOPMTR',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL SOPMTR('L','U','N',0,-1,a,tau,c,1,w,info)
         CALL chkxer('SOPMTR',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL SOPMTR('L','U','N',2,0,a,tau,c,1,w,info)
         CALL chkxer('SOPMTR',INFot,NOUt,LERr,OK)
         nt = nt + 6
!
!        SPTEQR
!
         SRNamt = 'SPTEQR'
         INFot = 1
         CALL SPTEQR('/',0,d,e,z,1,w,info)
         CALL chkxer('SPTEQR',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SPTEQR('N',-1,d,e,z,1,w,info)
         CALL chkxer('SPTEQR',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL SPTEQR('V',2,d,e,z,1,w,info)
         CALL chkxer('SPTEQR',INFot,NOUt,LERr,OK)
         nt = nt + 3
!
!        SSTEBZ
!
         SRNamt = 'SSTEBZ'
         INFot = 1
         CALL SSTEBZ('/','E',0,0.0,1.0,1,0,0.0,d,e,m,nsplit,x,i1,i2,w,  &
     &               iw,info)
         CALL chkxer('SSTEBZ',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SSTEBZ('A','/',0,0.0,0.0,0,0,0.0,d,e,m,nsplit,x,i1,i2,w,  &
     &               iw,info)
         CALL chkxer('SSTEBZ',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SSTEBZ('A','E',-1,0.0,0.0,0,0,0.0,d,e,m,nsplit,x,i1,i2,w, &
     &               iw,info)
         CALL chkxer('SSTEBZ',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL SSTEBZ('V','E',0,0.0,0.0,0,0,0.0,d,e,m,nsplit,x,i1,i2,w,  &
     &               iw,info)
         CALL chkxer('SSTEBZ',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL SSTEBZ('I','E',0,0.0,0.0,0,0,0.0,d,e,m,nsplit,x,i1,i2,w,  &
     &               iw,info)
         CALL chkxer('SSTEBZ',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL SSTEBZ('I','E',1,0.0,0.0,2,1,0.0,d,e,m,nsplit,x,i1,i2,w,  &
     &               iw,info)
         CALL chkxer('SSTEBZ',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL SSTEBZ('I','E',1,0.0,0.0,1,0,0.0,d,e,m,nsplit,x,i1,i2,w,  &
     &               iw,info)
         CALL chkxer('SSTEBZ',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL SSTEBZ('I','E',1,0.0,0.0,1,2,0.0,d,e,m,nsplit,x,i1,i2,w,  &
     &               iw,info)
         CALL chkxer('SSTEBZ',INFot,NOUt,LERr,OK)
         nt = nt + 8
!
!        SSTEIN
!
         SRNamt = 'SSTEIN'
         INFot = 1
         CALL SSTEIN(-1,d,e,0,x,i1,i2,z,1,w,iw,i3,info)
         CALL chkxer('SSTEIN',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL SSTEIN(0,d,e,-1,x,i1,i2,z,1,w,iw,i3,info)
         CALL chkxer('SSTEIN',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL SSTEIN(0,d,e,1,x,i1,i2,z,1,w,iw,i3,info)
         CALL chkxer('SSTEIN',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL SSTEIN(2,d,e,0,x,i1,i2,z,1,w,iw,i3,info)
         CALL chkxer('SSTEIN',INFot,NOUt,LERr,OK)
         nt = nt + 4
!
!        SSTEQR
!
         SRNamt = 'SSTEQR'
         INFot = 1
         CALL SSTEQR('/',0,d,e,z,1,w,info)
         CALL chkxer('SSTEQR',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SSTEQR('N',-1,d,e,z,1,w,info)
         CALL chkxer('SSTEQR',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL SSTEQR('V',2,d,e,z,1,w,info)
         CALL chkxer('SSTEQR',INFot,NOUt,LERr,OK)
         nt = nt + 3
!
!        SSTERF
!
         SRNamt = 'SSTERF'
         INFot = 1
         CALL SSTERF(-1,d,e,info)
         CALL chkxer('SSTERF',INFot,NOUt,LERr,OK)
         nt = nt + 1
!
!        SSTEDC
!
         SRNamt = 'SSTEDC'
         INFot = 1
         CALL SSTEDC('/',0,d,e,z,1,w,1,iw,1,info)
         CALL chkxer('SSTEDC',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SSTEDC('N',-1,d,e,z,1,w,1,iw,1,info)
         CALL chkxer('SSTEDC',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL SSTEDC('V',2,d,e,z,1,w,23,iw,28,info)
         CALL chkxer('SSTEDC',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL SSTEDC('N',1,d,e,z,1,w,0,iw,1,info)
         CALL chkxer('SSTEDC',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL SSTEDC('I',2,d,e,z,2,w,0,iw,12,info)
         CALL chkxer('SSTEDC',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL SSTEDC('V',2,d,e,z,2,w,0,iw,28,info)
         CALL chkxer('SSTEDC',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL SSTEDC('N',1,d,e,z,1,w,1,iw,0,info)
         CALL chkxer('SSTEDC',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL SSTEDC('I',2,d,e,z,2,w,19,iw,0,info)
         CALL chkxer('SSTEDC',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL SSTEDC('V',2,d,e,z,2,w,23,iw,0,info)
         CALL chkxer('SSTEDC',INFot,NOUt,LERr,OK)
         nt = nt + 9
!
!        SSTEVD
!
         SRNamt = 'SSTEVD'
         INFot = 1
         CALL SSTEVD('/',0,d,e,z,1,w,1,iw,1,info)
         CALL chkxer('SSTEVD',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SSTEVD('N',-1,d,e,z,1,w,1,iw,1,info)
         CALL chkxer('SSTEVD',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL SSTEVD('V',2,d,e,z,1,w,19,iw,12,info)
         CALL chkxer('SSTEVD',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL SSTEVD('N',1,d,e,z,1,w,0,iw,1,info)
         CALL chkxer('SSTEVD',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL SSTEVD('V',2,d,e,z,2,w,12,iw,12,info)
         CALL chkxer('SSTEVD',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL SSTEVD('N',0,d,e,z,1,w,1,iw,0,info)
         CALL chkxer('SSTEVD',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL SSTEVD('V',2,d,e,z,2,w,19,iw,11,info)
         CALL chkxer('SSTEVD',INFot,NOUt,LERr,OK)
         nt = nt + 7
!
!        SSTEV
!
         SRNamt = 'SSTEV '
         INFot = 1
         CALL SSTEV('/',0,d,e,z,1,w,info)
         CALL chkxer('SSTEV ',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SSTEV('N',-1,d,e,z,1,w,info)
         CALL chkxer('SSTEV ',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL SSTEV('V',2,d,e,z,1,w,info)
         CALL chkxer('SSTEV ',INFot,NOUt,LERr,OK)
         nt = nt + 3
!
!        SSTEVX
!
         SRNamt = 'SSTEVX'
         INFot = 1
         CALL SSTEVX('/','A',0,d,e,0.0,0.0,0,0,0.0,m,x,z,1,w,iw,i3,info)
         CALL chkxer('SSTEVX',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SSTEVX('N','/',0,d,e,0.0,1.0,1,0,0.0,m,x,z,1,w,iw,i3,info)
         CALL chkxer('SSTEVX',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SSTEVX('N','A',-1,d,e,0.0,0.0,0,0,0.0,m,x,z,1,w,iw,i3,    &
     &               info)
         CALL chkxer('SSTEVX',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL SSTEVX('N','V',1,d,e,0.0,0.0,0,0,0.0,m,x,z,1,w,iw,i3,info)
         CALL chkxer('SSTEVX',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL SSTEVX('N','I',1,d,e,0.0,0.0,0,0,0.0,m,x,z,1,w,iw,i3,info)
         CALL chkxer('SSTEVX',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL SSTEVX('N','I',1,d,e,0.0,0.0,2,1,0.0,m,x,z,1,w,iw,i3,info)
         CALL chkxer('SSTEVX',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL SSTEVX('N','I',2,d,e,0.0,0.0,2,1,0.0,m,x,z,1,w,iw,i3,info)
         CALL chkxer('SSTEVX',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL SSTEVX('N','I',1,d,e,0.0,0.0,1,2,0.0,m,x,z,1,w,iw,i3,info)
         CALL chkxer('SSTEVX',INFot,NOUt,LERr,OK)
         INFot = 14
         CALL SSTEVX('V','A',2,d,e,0.0,0.0,0,0,0.0,m,x,z,1,w,iw,i3,info)
         CALL chkxer('SSTEVX',INFot,NOUt,LERr,OK)
         nt = nt + 9
!
!        SSTEVR
!
         n = 1
         SRNamt = 'SSTEVR'
         INFot = 1
         CALL SSTEVR('/','A',0,d,e,0.0,0.0,1,1,0.0,m,r,z,1,iw,x,20*n,   &
     &               iw(2*n+1),10*n,info)
         CALL chkxer('SSTEVR',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SSTEVR('V','/',0,d,e,0.0,0.0,1,1,0.0,m,r,z,1,iw,x,20*n,   &
     &               iw(2*n+1),10*n,info)
         CALL chkxer('SSTEVR',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SSTEVR('V','A',-1,d,e,0.0,0.0,1,1,0.0,m,r,z,1,iw,x,20*n,  &
     &               iw(2*n+1),10*n,info)
         CALL chkxer('SSTEVR',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL SSTEVR('V','V',1,d,e,0.0,0.0,1,1,0.0,m,r,z,1,iw,x,20*n,   &
     &               iw(2*n+1),10*n,info)
         CALL chkxer('SSTEVR',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL SSTEVR('V','I',1,d,e,0.0,0.0,0,1,0.0,m,w,z,1,iw,x,20*n,   &
     &               iw(2*n+1),10*n,info)
         CALL chkxer('SSTEVR',INFot,NOUt,LERr,OK)
         INFot = 9
         n = 2
         CALL SSTEVR('V','I',2,d,e,0.0,0.0,2,1,0.0,m,w,z,1,iw,x,20*n,   &
     &               iw(2*n+1),10*n,info)
         CALL chkxer('SSTEVR',INFot,NOUt,LERr,OK)
         INFot = 14
         n = 1
         CALL SSTEVR('V','I',1,d,e,0.0,0.0,1,1,0.0,m,w,z,0,iw,x,20*n,   &
     &               iw(2*n+1),10*n,info)
         CALL chkxer('SSTEVR',INFot,NOUt,LERr,OK)
         INFot = 17
         CALL SSTEVR('V','I',1,d,e,0.0,0.0,1,1,0.0,m,w,z,1,iw,x,20*n-1, &
     &               iw(2*n+1),10*n,info)
         CALL chkxer('SSTEVR',INFot,NOUt,LERr,OK)
         INFot = 19
         CALL SSTEVR('V','I',1,d,e,0.0,0.0,1,1,0.0,m,w,z,1,iw,x,20*n,   &
     &               iw(2*n+1),10*n-1,info)
         CALL chkxer('SSTEVR',INFot,NOUt,LERr,OK)
         nt = nt + 9
!
!        SSYEVD
!
         SRNamt = 'SSYEVD'
         INFot = 1
         CALL SSYEVD('/','U',0,a,1,x,w,1,iw,1,info)
         CALL chkxer('SSYEVD',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SSYEVD('N','/',0,a,1,x,w,1,iw,1,info)
         CALL chkxer('SSYEVD',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SSYEVD('N','U',-1,a,1,x,w,1,iw,1,info)
         CALL chkxer('SSYEVD',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL SSYEVD('N','U',2,a,1,x,w,3,iw,1,info)
         CALL chkxer('SSYEVD',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL SSYEVD('N','U',1,a,1,x,w,0,iw,1,info)
         CALL chkxer('SSYEVD',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL SSYEVD('N','U',2,a,2,x,w,4,iw,1,info)
         CALL chkxer('SSYEVD',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL SSYEVD('V','U',2,a,2,x,w,20,iw,12,info)
         CALL chkxer('SSYEVD',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL SSYEVD('N','U',1,a,1,x,w,1,iw,0,info)
         CALL chkxer('SSYEVD',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL SSYEVD('N','U',2,a,2,x,w,5,iw,0,info)
         CALL chkxer('SSYEVD',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL SSYEVD('V','U',2,a,2,x,w,27,iw,11,info)
         CALL chkxer('SSYEVD',INFot,NOUt,LERr,OK)
         nt = nt + 10
!
!        SSYEVD_2STAGE
!
         SRNamt = 'SSYEVD_2STAGE'
         INFot = 1
         CALL SSYEVD_2STAGE('/','U',0,a,1,x,w,1,iw,1,info)
         CALL chkxer('SSYEVD_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 1
         CALL SSYEVD_2STAGE('V','U',0,a,1,x,w,1,iw,1,info)
         CALL chkxer('SSYEVD_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SSYEVD_2STAGE('N','/',0,a,1,x,w,1,iw,1,info)
         CALL chkxer('SSYEVD_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SSYEVD_2STAGE('N','U',-1,a,1,x,w,1,iw,1,info)
         CALL chkxer('SSYEVD_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL SSYEVD_2STAGE('N','U',2,a,1,x,w,3,iw,1,info)
         CALL chkxer('SSYEVD_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL SSYEVD_2STAGE('N','U',1,a,1,x,w,0,iw,1,info)
         CALL chkxer('SSYEVD_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL SSYEVD_2STAGE('N','U',2,a,2,x,w,4,iw,1,info)
         CALL chkxer('SSYEVD_2STAGE',INFot,NOUt,LERr,OK)
!         INFOT = 8
!         CALL SSYEVD_2STAGE( 'V', 'U', 2, A, 2, X, W, 20, IW, 12, INFO )
!         CALL chkxer( 'SSYEVD_2STAGE', INFOT, NOUT, LERR, OK )
         INFot = 10
         CALL SSYEVD_2STAGE('N','U',1,a,1,x,w,1,iw,0,info)
         CALL chkxer('SSYEVD_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL SSYEVD_2STAGE('N','U',2,a,2,x,w,25,iw,0,info)
         CALL chkxer('SSYEVD_2STAGE',INFot,NOUt,LERr,OK)
!         INFOT = 10
!         CALL SSYEVD_2STAGE( 'V', 'U', 2, A, 2, X, W, 27, IW, 11, INFO )
!         CALL chkxer( 'SSYEVD_2STAGE', INFOT, NOUT, LERR, OK )
         nt = nt + 9
!
!        SSYEVR
!
         SRNamt = 'SSYEVR'
         n = 1
         INFot = 1
         CALL SSYEVR('/','A','U',0,a,1,0.0,0.0,1,1,0.0,m,r,z,1,iw,q,    &
     &               26*n,iw(2*n+1),10*n,info)
         CALL chkxer('SSYEVR',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SSYEVR('V','/','U',0,a,1,0.0,0.0,1,1,0.0,m,r,z,1,iw,q,    &
     &               26*n,iw(2*n+1),10*n,info)
         CALL chkxer('SSYEVR',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SSYEVR('V','A','/',-1,a,1,0.0,0.0,1,1,0.0,m,r,z,1,iw,q,   &
     &               26*n,iw(2*n+1),10*n,info)
         CALL chkxer('SSYEVR',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL SSYEVR('V','A','U',-1,a,1,0.0,0.0,1,1,0.0,m,r,z,1,iw,q,   &
     &               26*n,iw(2*n+1),10*n,info)
         CALL chkxer('SSYEVR',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL SSYEVR('V','A','U',2,a,1,0.0,0.0,1,1,0.0,m,r,z,1,iw,q,    &
     &               26*n,iw(2*n+1),10*n,info)
         CALL chkxer('SSYEVR',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL SSYEVR('V','V','U',1,a,1,0.0E0,0.0E0,1,1,0.0,m,r,z,1,iw,q,&
     &               26*n,iw(2*n+1),10*n,info)
         CALL chkxer('SSYEVR',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL SSYEVR('V','I','U',1,a,1,0.0E0,0.0E0,0,1,0.0,m,r,z,1,iw,q,&
     &               26*n,iw(2*n+1),10*n,info)
         CALL chkxer('SSYEVR',INFot,NOUt,LERr,OK)
         INFot = 10
!
         CALL SSYEVR('V','I','U',2,a,2,0.0E0,0.0E0,2,1,0.0,m,r,z,1,iw,q,&
     &               26*n,iw(2*n+1),10*n,info)
         CALL chkxer('SSYEVR',INFot,NOUt,LERr,OK)
         INFot = 15
         CALL SSYEVR('V','I','U',1,a,1,0.0E0,0.0E0,1,1,0.0,m,r,z,0,iw,q,&
     &               26*n,iw(2*n+1),10*n,info)
         CALL chkxer('SSYEVR',INFot,NOUt,LERr,OK)
         INFot = 18
         CALL SSYEVR('V','I','U',1,a,1,0.0E0,0.0E0,1,1,0.0,m,r,z,1,iw,q,&
     &               26*n-1,iw(2*n+1),10*n,info)
         CALL chkxer('SSYEVR',INFot,NOUt,LERr,OK)
         INFot = 20
         CALL SSYEVR('V','I','U',1,a,1,0.0E0,0.0E0,1,1,0.0,m,r,z,1,iw,q,&
     &               26*n,iw(2*n+1),10*n-1,info)
         CALL chkxer('SSYEVR',INFot,NOUt,LERr,OK)
         nt = nt + 11
!
!        SSYEVR_2STAGE
!
         SRNamt = 'SSYEVR_2STAGE'
         n = 1
         INFot = 1
         CALL SSYEVR_2STAGE('/','A','U',0,a,1,0.0E0,0.0E0,1,1,0.0E0,m,r,&
     &                      z,1,iw,q,26*n,iw(2*n+1),10*n,info)
         CALL chkxer('SSYEVR_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 1
         CALL SSYEVR_2STAGE('V','A','U',0,a,1,0.0E0,0.0E0,1,1,0.0E0,m,r,&
     &                      z,1,iw,q,26*n,iw(2*n+1),10*n,info)
         CALL chkxer('SSYEVR_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SSYEVR_2STAGE('N','/','U',0,a,1,0.0E0,0.0E0,1,1,0.0E0,m,r,&
     &                      z,1,iw,q,26*n,iw(2*n+1),10*n,info)
         CALL chkxer('SSYEVR_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SSYEVR_2STAGE('N','A','/',-1,a,1,0.0E0,0.0E0,1,1,0.0E0,m, &
     &                      r,z,1,iw,q,26*n,iw(2*n+1),10*n,info)
         CALL chkxer('SSYEVR_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL SSYEVR_2STAGE('N','A','U',-1,a,1,0.0E0,0.0E0,1,1,0.0E0,m, &
     &                      r,z,1,iw,q,26*n,iw(2*n+1),10*n,info)
         CALL chkxer('SSYEVR_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL SSYEVR_2STAGE('N','A','U',2,a,1,0.0E0,0.0E0,1,1,0.0E0,m,r,&
     &                      z,1,iw,q,26*n,iw(2*n+1),10*n,info)
         CALL chkxer('SSYEVR_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL SSYEVR_2STAGE('N','V','U',1,a,1,0.0E0,0.0E0,1,1,0.0E0,m,r,&
     &                      z,1,iw,q,26*n,iw(2*n+1),10*n,info)
         CALL chkxer('SSYEVR_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL SSYEVR_2STAGE('N','I','U',1,a,1,0.0E0,0.0E0,0,1,0.0E0,m,r,&
     &                      z,1,iw,q,26*n,iw(2*n+1),10*n,info)
         CALL chkxer('SSYEVR_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL SSYEVR_2STAGE('N','I','U',2,a,2,0.0E0,0.0E0,2,1,0.0E0,m,r,&
     &                      z,1,iw,q,26*n,iw(2*n+1),10*n,info)
         CALL chkxer('SSYEVR_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 15
         CALL SSYEVR_2STAGE('N','I','U',1,a,1,0.0E0,0.0E0,1,1,0.0E0,m,r,&
     &                      z,0,iw,q,26*n,iw(2*n+1),10*n,info)
         CALL chkxer('SSYEVR_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 18
         CALL SSYEVR_2STAGE('N','I','U',1,a,1,0.0E0,0.0E0,1,1,0.0E0,m,r,&
     &                      z,1,iw,q,0,iw(2*n+1),10*n,info)
         CALL chkxer('SSYEVR_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 20
         CALL SSYEVR_2STAGE('N','I','U',1,a,1,0.0E0,0.0E0,1,1,0.0E0,m,r,&
     &                      z,1,iw,q,26*n,iw(2*n+1),0,info)
         CALL chkxer('SSYEVR_2STAGE',INFot,NOUt,LERr,OK)
         nt = nt + 12
!
!        SSYEV
!
         SRNamt = 'SSYEV '
         INFot = 1
         CALL SSYEV('/','U',0,a,1,x,w,1,info)
         CALL chkxer('SSYEV ',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SSYEV('N','/',0,a,1,x,w,1,info)
         CALL chkxer('SSYEV ',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SSYEV('N','U',-1,a,1,x,w,1,info)
         CALL chkxer('SSYEV ',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL SSYEV('N','U',2,a,1,x,w,3,info)
         CALL chkxer('SSYEV ',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL SSYEV('N','U',1,a,1,x,w,1,info)
         CALL chkxer('SSYEV ',INFot,NOUt,LERr,OK)
         nt = nt + 5
!
!        SSYEV_2STAGE
!
         SRNamt = 'SSYEV_2STAGE '
         INFot = 1
         CALL SSYEV_2STAGE('/','U',0,a,1,x,w,1,info)
         CALL chkxer('SSYEV_2STAGE ',INFot,NOUt,LERr,OK)
         INFot = 1
         CALL SSYEV_2STAGE('V','U',0,a,1,x,w,1,info)
         CALL chkxer('SSYEV_2STAGE ',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SSYEV_2STAGE('N','/',0,a,1,x,w,1,info)
         CALL chkxer('SSYEV_2STAGE ',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SSYEV_2STAGE('N','U',-1,a,1,x,w,1,info)
         CALL chkxer('SSYEV_2STAGE ',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL SSYEV_2STAGE('N','U',2,a,1,x,w,3,info)
         CALL chkxer('SSYEV_2STAGE ',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL SSYEV_2STAGE('N','U',1,a,1,x,w,1,info)
         CALL chkxer('SSYEV_2STAGE ',INFot,NOUt,LERr,OK)
         nt = nt + 6
!
!        SSYEVX
!
         SRNamt = 'SSYEVX'
         INFot = 1
         CALL SSYEVX('/','A','U',0,a,1,0.0,0.0,0,0,0.0,m,x,z,1,w,1,iw,  &
     &               i3,info)
         CALL chkxer('SSYEVX',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SSYEVX('N','/','U',0,a,1,0.0,1.0,1,0,0.0,m,x,z,1,w,1,iw,  &
     &               i3,info)
         CALL chkxer('SSYEVX',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SSYEVX('N','A','/',0,a,1,0.0,0.0,0,0,0.0,m,x,z,1,w,1,iw,  &
     &               i3,info)
         INFot = 4
         CALL SSYEVX('N','A','U',-1,a,1,0.0,0.0,0,0,0.0,m,x,z,1,w,1,iw, &
     &               i3,info)
         CALL chkxer('SSYEVX',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL SSYEVX('N','A','U',2,a,1,0.0,0.0,0,0,0.0,m,x,z,1,w,16,iw, &
     &               i3,info)
         CALL chkxer('SSYEVX',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL SSYEVX('N','V','U',1,a,1,0.0,0.0,0,0,0.0,m,x,z,1,w,8,iw,  &
     &               i3,info)
         CALL chkxer('SSYEVX',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL SSYEVX('N','I','U',1,a,1,0.0,0.0,0,0,0.0,m,x,z,1,w,8,iw,  &
     &               i3,info)
         CALL chkxer('SSYEVX',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL SSYEVX('N','I','U',1,a,1,0.0,0.0,2,1,0.0,m,x,z,1,w,8,iw,  &
     &               i3,info)
         CALL chkxer('SSYEVX',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL SSYEVX('N','I','U',2,a,2,0.0,0.0,2,1,0.0,m,x,z,1,w,16,iw, &
     &               i3,info)
         CALL chkxer('SSYEVX',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL SSYEVX('N','I','U',1,a,1,0.0,0.0,1,2,0.0,m,x,z,1,w,8,iw,  &
     &               i3,info)
         CALL chkxer('SSYEVX',INFot,NOUt,LERr,OK)
         INFot = 15
         CALL SSYEVX('V','A','U',2,a,2,0.0,0.0,0,0,0.0,m,x,z,1,w,16,iw, &
     &               i3,info)
         CALL chkxer('SSYEVX',INFot,NOUt,LERr,OK)
         INFot = 17
         CALL SSYEVX('V','A','U',1,a,1,0.0,0.0,0,0,0.0,m,x,z,1,w,0,iw,  &
     &               i3,info)
         CALL chkxer('SSYEVX',INFot,NOUt,LERr,OK)
         nt = nt + 12
!
!        SSYEVX_2STAGE
!
         SRNamt = 'SSYEVX_2STAGE'
         INFot = 1
         CALL SSYEVX_2STAGE('/','A','U',0,a,1,0.0E0,0.0E0,0,0,0.0E0,m,x,&
     &                      z,1,w,1,iw,i3,info)
         CALL chkxer('SSYEVX_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 1
         CALL SSYEVX_2STAGE('V','A','U',0,a,1,0.0E0,0.0E0,0,0,0.0E0,m,x,&
     &                      z,1,w,1,iw,i3,info)
         CALL chkxer('SSYEVX_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SSYEVX_2STAGE('N','/','U',0,a,1,0.0E0,1.0E0,1,0,0.0E0,m,x,&
     &                      z,1,w,1,iw,i3,info)
         CALL chkxer('SSYEVX_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SSYEVX_2STAGE('N','A','/',0,a,1,0.0E0,0.0E0,0,0,0.0E0,m,x,&
     &                      z,1,w,1,iw,i3,info)
         INFot = 4
         CALL SSYEVX_2STAGE('N','A','U',-1,a,1,0.0E0,0.0E0,0,0,0.0E0,m, &
     &                      x,z,1,w,1,iw,i3,info)
         CALL chkxer('SSYEVX_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL SSYEVX_2STAGE('N','A','U',2,a,1,0.0E0,0.0E0,0,0,0.0E0,m,x,&
     &                      z,1,w,16,iw,i3,info)
         CALL chkxer('SSYEVX_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL SSYEVX_2STAGE('N','V','U',1,a,1,0.0E0,0.0E0,0,0,0.0E0,m,x,&
     &                      z,1,w,8,iw,i3,info)
         CALL chkxer('SSYEVX_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL SSYEVX_2STAGE('N','I','U',1,a,1,0.0E0,0.0E0,0,0,0.0E0,m,x,&
     &                      z,1,w,8,iw,i3,info)
         CALL chkxer('SSYEVX_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL SSYEVX_2STAGE('N','I','U',1,a,1,0.0E0,0.0E0,2,1,0.0E0,m,x,&
     &                      z,1,w,8,iw,i3,info)
         CALL chkxer('SSYEVX_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL SSYEVX_2STAGE('N','I','U',2,a,2,0.0E0,0.0E0,2,1,0.0E0,m,x,&
     &                      z,1,w,16,iw,i3,info)
         CALL chkxer('SSYEVX_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL SSYEVX_2STAGE('N','I','U',1,a,1,0.0E0,0.0E0,1,2,0.0E0,m,x,&
     &                      z,1,w,8,iw,i3,info)
         CALL chkxer('SSYEVX_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 15
         CALL SSYEVX_2STAGE('N','A','U',2,a,2,0.0E0,0.0E0,0,0,0.0E0,m,x,&
     &                      z,0,w,16,iw,i3,info)
         CALL chkxer('SSYEVX_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 17
         CALL SSYEVX_2STAGE('N','A','U',1,a,1,0.0E0,0.0E0,0,0,0.0E0,m,x,&
     &                      z,1,w,0,iw,i3,info)
         CALL chkxer('SSYEVX_2STAGE',INFot,NOUt,LERr,OK)
         nt = nt + 13
!
!        SSPEVD
!
         SRNamt = 'SSPEVD'
         INFot = 1
         CALL SSPEVD('/','U',0,a,x,z,1,w,1,iw,1,info)
         CALL chkxer('SSPEVD',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SSPEVD('N','/',0,a,x,z,1,w,1,iw,1,info)
         CALL chkxer('SSPEVD',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SSPEVD('N','U',-1,a,x,z,1,w,1,iw,1,info)
         CALL chkxer('SSPEVD',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL SSPEVD('V','U',2,a,x,z,1,w,23,iw,12,info)
         CALL chkxer('SSPEVD',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL SSPEVD('N','U',1,a,x,z,1,w,0,iw,1,info)
         CALL chkxer('SSPEVD',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL SSPEVD('N','U',2,a,x,z,1,w,3,iw,1,info)
         CALL chkxer('SSPEVD',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL SSPEVD('V','U',2,a,x,z,2,w,16,iw,12,info)
         CALL chkxer('SSPEVD',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL SSPEVD('N','U',1,a,x,z,1,w,1,iw,0,info)
         CALL chkxer('SSPEVD',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL SSPEVD('N','U',2,a,x,z,1,w,4,iw,0,info)
         CALL chkxer('SSPEVD',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL SSPEVD('V','U',2,a,x,z,2,w,23,iw,11,info)
         CALL chkxer('SSPEVD',INFot,NOUt,LERr,OK)
         nt = nt + 10
!
!        SSPEV
!
         SRNamt = 'SSPEV '
         INFot = 1
         CALL SSPEV('/','U',0,a,w,z,1,x,info)
         CALL chkxer('SSPEV ',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SSPEV('N','/',0,a,w,z,1,x,info)
         CALL chkxer('SSPEV ',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SSPEV('N','U',-1,a,w,z,1,x,info)
         CALL chkxer('SSPEV ',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL SSPEV('V','U',2,a,w,z,1,x,info)
         CALL chkxer('SSPEV ',INFot,NOUt,LERr,OK)
         nt = nt + 4
!
!        SSPEVX
!
         SRNamt = 'SSPEVX'
         INFot = 1
         CALL SSPEVX('/','A','U',0,a,0.0,0.0,0,0,0.0,m,x,z,1,w,iw,i3,   &
     &               info)
         CALL chkxer('SSPEVX',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SSPEVX('N','/','U',0,a,0.0,0.0,0,0,0.0,m,x,z,1,w,iw,i3,   &
     &               info)
         CALL chkxer('SSPEVX',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SSPEVX('N','A','/',0,a,0.0,0.0,0,0,0.0,m,x,z,1,w,iw,i3,   &
     &               info)
         INFot = 4
         CALL SSPEVX('N','A','U',-1,a,0.0,0.0,0,0,0.0,m,x,z,1,w,iw,i3,  &
     &               info)
         CALL chkxer('SSPEVX',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL SSPEVX('N','V','U',1,a,0.0,0.0,0,0,0.0,m,x,z,1,w,iw,i3,   &
     &               info)
         CALL chkxer('SSPEVX',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL SSPEVX('N','I','U',1,a,0.0,0.0,0,0,0.0,m,x,z,1,w,iw,i3,   &
     &               info)
         CALL chkxer('SSPEVX',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL SSPEVX('N','I','U',1,a,0.0,0.0,2,1,0.0,m,x,z,1,w,iw,i3,   &
     &               info)
         CALL chkxer('SSPEVX',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL SSPEVX('N','I','U',2,a,0.0,0.0,2,1,0.0,m,x,z,1,w,iw,i3,   &
     &               info)
         CALL chkxer('SSPEVX',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL SSPEVX('N','I','U',1,a,0.0,0.0,1,2,0.0,m,x,z,1,w,iw,i3,   &
     &               info)
         CALL chkxer('SSPEVX',INFot,NOUt,LERr,OK)
         INFot = 14
         CALL SSPEVX('V','A','U',2,a,0.0,0.0,0,0,0.0,m,x,z,1,w,iw,i3,   &
     &               info)
         CALL chkxer('SSPEVX',INFot,NOUt,LERr,OK)
         nt = nt + 10
!
!     Test error exits for the SB path.
!
      ELSEIF ( LSAMEN(2,c2,'SB') ) THEN
!
!        SSBTRD
!
         SRNamt = 'SSBTRD'
         INFot = 1
         CALL SSBTRD('/','U',0,0,a,1,d,e,z,1,w,info)
         CALL chkxer('SSBTRD',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SSBTRD('N','/',0,0,a,1,d,e,z,1,w,info)
         CALL chkxer('SSBTRD',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SSBTRD('N','U',-1,0,a,1,d,e,z,1,w,info)
         CALL chkxer('SSBTRD',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL SSBTRD('N','U',0,-1,a,1,d,e,z,1,w,info)
         CALL chkxer('SSBTRD',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL SSBTRD('N','U',1,1,a,1,d,e,z,1,w,info)
         CALL chkxer('SSBTRD',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL SSBTRD('V','U',2,0,a,1,d,e,z,1,w,info)
         CALL chkxer('SSBTRD',INFot,NOUt,LERr,OK)
         nt = nt + 6
!
!        SSYTRD_SB2ST
!
         SRNamt = 'SSYTRD_SB2ST'
         INFot = 1
         CALL SSYTRD_SB2ST('/','N','U',0,0,a,1,d,e,c,1,w,1,info)
         CALL chkxer('SSYTRD_SB2ST',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SSYTRD_SB2ST('N','/','U',0,0,a,1,d,e,c,1,w,1,info)
         CALL chkxer('SSYTRD_SB2ST',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SSYTRD_SB2ST('N','H','U',0,0,a,1,d,e,c,1,w,1,info)
         CALL chkxer('SSYTRD_SB2ST',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SSYTRD_SB2ST('N','N','/',0,0,a,1,d,e,c,1,w,1,info)
         CALL chkxer('SSYTRD_SB2ST',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL SSYTRD_SB2ST('N','N','U',-1,0,a,1,d,e,c,1,w,1,info)
         CALL chkxer('SSYTRD_SB2ST',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL SSYTRD_SB2ST('N','N','U',0,-1,a,1,d,e,c,1,w,1,info)
         CALL chkxer('SSYTRD_SB2ST',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL SSYTRD_SB2ST('N','N','U',0,1,a,1,d,e,c,1,w,1,info)
         CALL chkxer('SSYTRD_SB2ST',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL SSYTRD_SB2ST('N','N','U',0,0,a,1,d,e,c,0,w,1,info)
         CALL chkxer('SSYTRD_SB2ST',INFot,NOUt,LERr,OK)
         INFot = 13
         CALL SSYTRD_SB2ST('N','N','U',0,0,a,1,d,e,c,1,w,0,info)
         CALL chkxer('SSYTRD_SB2ST',INFot,NOUt,LERr,OK)
         nt = nt + 9
!
!        SSBEVD
!
         SRNamt = 'SSBEVD'
         INFot = 1
         CALL SSBEVD('/','U',0,0,a,1,x,z,1,w,1,iw,1,info)
         CALL chkxer('SSBEVD',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SSBEVD('N','/',0,0,a,1,x,z,1,w,1,iw,1,info)
         CALL chkxer('SSBEVD',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SSBEVD('N','U',-1,0,a,1,x,z,1,w,1,iw,1,info)
         CALL chkxer('SSBEVD',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL SSBEVD('N','U',0,-1,a,1,x,z,1,w,1,iw,1,info)
         CALL chkxer('SSBEVD',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL SSBEVD('N','U',2,1,a,1,x,z,1,w,4,iw,1,info)
         CALL chkxer('SSBEVD',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL SSBEVD('V','U',2,1,a,2,x,z,1,w,25,iw,12,info)
         CALL chkxer('SSBEVD',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL SSBEVD('N','U',1,0,a,1,x,z,1,w,0,iw,1,info)
         CALL chkxer('SSBEVD',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL SSBEVD('N','U',2,0,a,1,x,z,1,w,3,iw,1,info)
         CALL chkxer('SSBEVD',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL SSBEVD('V','U',2,0,a,1,x,z,2,w,18,iw,12,info)
         CALL chkxer('SSBEVD',INFot,NOUt,LERr,OK)
         INFot = 13
         CALL SSBEVD('N','U',1,0,a,1,x,z,1,w,1,iw,0,info)
         CALL chkxer('SSBEVD',INFot,NOUt,LERr,OK)
         INFot = 13
         CALL SSBEVD('V','U',2,0,a,1,x,z,2,w,25,iw,11,info)
         CALL chkxer('SSBEVD',INFot,NOUt,LERr,OK)
         nt = nt + 11
!
!        SSBEVD_2STAGE
!
         SRNamt = 'SSBEVD_2STAGE'
         INFot = 1
         CALL SSBEVD_2STAGE('/','U',0,0,a,1,x,z,1,w,1,iw,1,info)
         CALL chkxer('SSBEVD_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 1
         CALL SSBEVD_2STAGE('V','U',0,0,a,1,x,z,1,w,1,iw,1,info)
         CALL chkxer('SSBEVD_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SSBEVD_2STAGE('N','/',0,0,a,1,x,z,1,w,1,iw,1,info)
         CALL chkxer('SSBEVD_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SSBEVD_2STAGE('N','U',-1,0,a,1,x,z,1,w,1,iw,1,info)
         CALL chkxer('SSBEVD_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL SSBEVD_2STAGE('N','U',0,-1,a,1,x,z,1,w,1,iw,1,info)
         CALL chkxer('SSBEVD_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL SSBEVD_2STAGE('N','U',2,1,a,1,x,z,1,w,4,iw,1,info)
         CALL chkxer('SSBEVD_2STAGE',INFot,NOUt,LERr,OK)
!         INFOT = 9
!         CALL SSBEVD_2STAGE( 'V', 'U', 2, 1, A, 2, X, Z, 1, W,
!     $                                      25, IW, 12, INFO )
!         CALL chkxer( 'SSBEVD_2STAGE', INFOT, NOUT, LERR, OK )
         INFot = 11
         CALL SSBEVD_2STAGE('N','U',1,0,a,1,x,z,1,w,0,iw,1,info)
         CALL chkxer('SSBEVD_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL SSBEVD_2STAGE('N','U',2,0,a,1,x,z,1,w,3,iw,1,info)
         CALL chkxer('SSBEVD_2STAGE',INFot,NOUt,LERr,OK)
!         INFOT = 11
!         CALL SSBEVD_2STAGE( 'V', 'U', 2, 0, A, 1, X, Z, 2, W,
!     $                                      18, IW, 12, INFO )
!         CALL chkxer( 'SSBEVD_2STAGE', INFOT, NOUT, LERR, OK )
         INFot = 13
         CALL SSBEVD_2STAGE('N','U',1,0,a,1,x,z,1,w,1,iw,0,info)
         CALL chkxer('SSBEVD_2STAGE',INFot,NOUt,LERr,OK)
!         INFOT = 13
!         CALL SSBEVD_2STAGE( 'V', 'U', 2, 0, A, 1, X, Z, 2, W,
!     $                                      25, IW, 11, INFO )
!         CALL chkxer( 'SSBEVD_2STAGE', INFOT, NOUT, LERR, OK )
!         NT = NT + 12
         nt = nt + 9
!
!        SSBEV
!
         SRNamt = 'SSBEV '
         INFot = 1
         CALL SSBEV('/','U',0,0,a,1,x,z,1,w,info)
         CALL chkxer('SSBEV ',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SSBEV('N','/',0,0,a,1,x,z,1,w,info)
         CALL chkxer('SSBEV ',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SSBEV('N','U',-1,0,a,1,x,z,1,w,info)
         CALL chkxer('SSBEV ',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL SSBEV('N','U',0,-1,a,1,x,z,1,w,info)
         CALL chkxer('SSBEV ',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL SSBEV('N','U',2,1,a,1,x,z,1,w,info)
         CALL chkxer('SSBEV ',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL SSBEV('V','U',2,0,a,1,x,z,1,w,info)
         CALL chkxer('SSBEV ',INFot,NOUt,LERr,OK)
         nt = nt + 6
!
!        SSBEV_2STAGE
!
         SRNamt = 'SSBEV_2STAGE '
         INFot = 1
         CALL SSBEV_2STAGE('/','U',0,0,a,1,x,z,1,w,0,info)
         CALL chkxer('SSBEV_2STAGE ',INFot,NOUt,LERr,OK)
         INFot = 1
         CALL SSBEV_2STAGE('V','U',0,0,a,1,x,z,1,w,0,info)
         CALL chkxer('SSBEV_2STAGE ',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SSBEV_2STAGE('N','/',0,0,a,1,x,z,1,w,0,info)
         CALL chkxer('SSBEV_2STAGE ',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SSBEV_2STAGE('N','U',-1,0,a,1,x,z,1,w,0,info)
         CALL chkxer('SSBEV_2STAGE ',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL SSBEV_2STAGE('N','U',0,-1,a,1,x,z,1,w,0,info)
         CALL chkxer('SSBEV_2STAGE ',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL SSBEV_2STAGE('N','U',2,1,a,1,x,z,1,w,0,info)
         CALL chkxer('SSBEV_2STAGE ',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL SSBEV_2STAGE('N','U',2,0,a,1,x,z,0,w,0,info)
         CALL chkxer('SSBEV_2STAGE ',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL SSBEV_2STAGE('N','U',0,0,a,1,x,z,1,w,0,info)
         CALL chkxer('SSBEV_2STAGE ',INFot,NOUt,LERr,OK)
         nt = nt + 8
!
!        SSBEVX
!
         SRNamt = 'SSBEVX'
         INFot = 1
         CALL SSBEVX('/','A','U',0,0,a,1,q,1,0.0,0.0,0,0,0.0,m,x,z,1,w, &
     &               iw,i3,info)
         CALL chkxer('SSBEVX',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SSBEVX('N','/','U',0,0,a,1,q,1,0.0,0.0,0,0,0.0,m,x,z,1,w, &
     &               iw,i3,info)
         CALL chkxer('SSBEVX',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SSBEVX('N','A','/',0,0,a,1,q,1,0.0,0.0,0,0,0.0,m,x,z,1,w, &
     &               iw,i3,info)
         CALL chkxer('SSBEVX',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL SSBEVX('N','A','U',-1,0,a,1,q,1,0.0,0.0,0,0,0.0,m,x,z,1,w,&
     &               iw,i3,info)
         CALL chkxer('SSBEVX',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL SSBEVX('N','A','U',0,-1,a,1,q,1,0.0,0.0,0,0,0.0,m,x,z,1,w,&
     &               iw,i3,info)
         CALL chkxer('SSBEVX',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL SSBEVX('N','A','U',2,1,a,1,q,1,0.0,0.0,0,0,0.0,m,x,z,1,w, &
     &               iw,i3,info)
         CALL chkxer('SSBEVX',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL SSBEVX('V','A','U',2,0,a,1,q,1,0.0,0.0,0,0,0.0,m,x,z,2,w, &
     &               iw,i3,info)
         CALL chkxer('SSBEVX',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL SSBEVX('N','V','U',1,0,a,1,q,1,0.0,0.0,0,0,0.0,m,x,z,1,w, &
     &               iw,i3,info)
         CALL chkxer('SSBEVX',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL SSBEVX('N','I','U',1,0,a,1,q,1,0.0,0.0,0,0,0.0,m,x,z,1,w, &
     &               iw,i3,info)
         CALL chkxer('SSBEVX',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL SSBEVX('N','I','U',1,0,a,1,q,1,0.0,0.0,2,1,0.0,m,x,z,1,w, &
     &               iw,i3,info)
         CALL chkxer('SSBEVX',INFot,NOUt,LERr,OK)
         INFot = 13
         CALL SSBEVX('N','I','U',2,0,a,1,q,1,0.0,0.0,2,1,0.0,m,x,z,1,w, &
     &               iw,i3,info)
         CALL chkxer('SSBEVX',INFot,NOUt,LERr,OK)
         INFot = 13
         CALL SSBEVX('N','I','U',1,0,a,1,q,1,0.0,0.0,1,2,0.0,m,x,z,1,w, &
     &               iw,i3,info)
         CALL chkxer('SSBEVX',INFot,NOUt,LERr,OK)
         INFot = 18
         CALL SSBEVX('V','A','U',2,0,a,1,q,2,0.0,0.0,0,0,0.0,m,x,z,1,w, &
     &               iw,i3,info)
         CALL chkxer('SSBEVX',INFot,NOUt,LERr,OK)
         nt = nt + 13
!
!        SSBEVX_2STAGE
!
         SRNamt = 'SSBEVX_2STAGE'
         INFot = 1
         CALL SSBEVX_2STAGE('/','A','U',0,0,a,1,q,1,0.0E0,0.0E0,0,0,    &
     &                      0.0E0,m,x,z,1,w,0,iw,i3,info)
         CALL chkxer('SSBEVX_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 1
         CALL SSBEVX_2STAGE('V','A','U',0,0,a,1,q,1,0.0E0,0.0E0,0,0,    &
     &                      0.0E0,m,x,z,1,w,0,iw,i3,info)
         CALL chkxer('SSBEVX_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SSBEVX_2STAGE('N','/','U',0,0,a,1,q,1,0.0E0,0.0E0,0,0,    &
     &                      0.0E0,m,x,z,1,w,0,iw,i3,info)
         CALL chkxer('SSBEVX_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SSBEVX_2STAGE('N','A','/',0,0,a,1,q,1,0.0E0,0.0E0,0,0,    &
     &                      0.0E0,m,x,z,1,w,0,iw,i3,info)
         CALL chkxer('SSBEVX_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL SSBEVX_2STAGE('N','A','U',-1,0,a,1,q,1,0.0E0,0.0E0,0,0,   &
     &                      0.0E0,m,x,z,1,w,0,iw,i3,info)
         CALL chkxer('SSBEVX_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL SSBEVX_2STAGE('N','A','U',0,-1,a,1,q,1,0.0E0,0.0E0,0,0,   &
     &                      0.0E0,m,x,z,1,w,0,iw,i3,info)
         CALL chkxer('SSBEVX_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL SSBEVX_2STAGE('N','A','U',2,1,a,1,q,1,0.0E0,0.0E0,0,0,    &
     &                      0.0E0,m,x,z,1,w,0,iw,i3,info)
         CALL chkxer('SSBEVX_2STAGE',INFot,NOUt,LERr,OK)
!         INFOT = 9
!         CALL SSBEVX_2STAGE( 'V', 'A', 'U', 2, 0, A, 1, Q, 1, 0.0E0,
!     $          0.0E0, 0, 0, 0.0E0, M, X, Z, 2, W, 0, IW, I3, INFO )
!         CALL chkxer( 'SSBEVX_2STAGE', INFOT, NOUT, LERR, OK )
         INFot = 11
         CALL SSBEVX_2STAGE('N','V','U',1,0,a,1,q,1,0.0E0,0.0E0,0,0,    &
     &                      0.0E0,m,x,z,1,w,0,iw,i3,info)
         CALL chkxer('SSBEVX_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL SSBEVX_2STAGE('N','I','U',1,0,a,1,q,1,0.0E0,0.0E0,0,0,    &
     &                      0.0E0,m,x,z,1,w,0,iw,i3,info)
         CALL chkxer('SSBEVX_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL SSBEVX_2STAGE('N','I','U',1,0,a,1,q,1,0.0E0,0.0E0,2,1,    &
     &                      0.0E0,m,x,z,1,w,0,iw,i3,info)
         CALL chkxer('SSBEVX_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 13
         CALL SSBEVX_2STAGE('N','I','U',2,0,a,1,q,1,0.0E0,0.0E0,2,1,    &
     &                      0.0E0,m,x,z,1,w,0,iw,i3,info)
         CALL chkxer('SSBEVX_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 13
         CALL SSBEVX_2STAGE('N','I','U',1,0,a,1,q,1,0.0E0,0.0E0,1,2,    &
     &                      0.0E0,m,x,z,1,w,0,iw,i3,info)
         CALL chkxer('SSBEVX_2STAGE',INFot,NOUt,LERr,OK)
!         INFOT = 18
!         CALL SSBEVX_2STAGE( 'V', 'A', 'U', 2, 0, A, 1, Q, 2, 0.0E0,
!     $          0.0E0, 0, 0, 0.0E0, M, X, Z, 1, W, 0, IW, I3, INFO )
!         CALL chkxer( 'SSBEVX_2STAGE', INFOT, NOUT, LERR, OK )
         INFot = 20
         CALL SSBEVX_2STAGE('N','A','U',0,0,a,1,q,1,0.0E0,0.0E0,0,0,    &
     &                      0.0E0,m,x,z,1,w,0,iw,i3,info)
         CALL chkxer('SSBEVX_2STAGE',INFot,NOUt,LERr,OK)
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
!     End of serrst
!
      END SUBROUTINE serrst
