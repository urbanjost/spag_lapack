!*==zerrst.f90  processed by SPAG 7.51RB at 17:37 on  4 Mar 2022
!> \brief \b ZERRST
!
!  @precisions fortran z -> c
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZERRST( PATH, NUNIT )
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
!> ZERRST tests the error exits for ZHETRD, ZUNGTR, CUNMTR, ZHPTRD,
!> ZUNGTR, ZUPMTR, ZSTEQR, CSTEIN, ZPTEQR, ZHBTRD,
!> ZHEEV, CHEEVX, CHEEVD, ZHBEV, CHBEVX, CHBEVD,
!> ZHPEV, CHPEVX, CHPEVD, and ZSTEDC.
!> ZHEEVD_2STAGE, ZHEEVR_2STAGE, ZHEEVX_2STAGE,
!> ZHEEV_2STAGE, ZHBEV_2STAGE, ZHBEVD_2STAGE,
!> ZHBEVX_2STAGE, ZHETRD_2STAGE
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
!> \date June 2017
!
!> \ingroup complex16_eig
!
!  =====================================================================
      SUBROUTINE ZERRST(Path,Nunit)
      IMPLICIT NONE
!*--ZERRST66
!
!  -- LAPACK test routine (version 3.7.1) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     June 2017
!
!     .. Scalar Arguments ..
      CHARACTER*3 Path
      INTEGER Nunit
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      INTEGER NMAX , LIW , LW
      PARAMETER (NMAX=3,LIW=12*NMAX,LW=20*NMAX)
!     ..
!     .. Local Scalars ..
      CHARACTER*2 c2
      INTEGER i , info , j , m , n , nt
!     ..
!     .. Local Arrays ..
      INTEGER i1(NMAX) , i2(NMAX) , i3(NMAX) , iw(LIW)
      DOUBLE PRECISION d(NMAX) , e(NMAX) , r(LW) , rw(LW) , x(NMAX)
      COMPLEX*16 a(NMAX,NMAX) , c(NMAX,NMAX) , q(NMAX,NMAX) , tau(NMAX) &
     &           , w(LW) , z(NMAX,NMAX)
!     ..
!     .. External Functions ..
      LOGICAL LSAMEN
      EXTERNAL LSAMEN
!     ..
!     .. External Subroutines ..
      EXTERNAL CHKXER , ZHBEV , ZHBEVD , ZHBEVX , ZHBTRD , ZHEEV ,      &
     &         ZHEEVD , ZHEEVR , ZHEEVX , ZHETRD , ZHPEV , ZHPEVD ,     &
     &         ZHPEVX , ZHPTRD , ZPTEQR , ZSTEDC , ZSTEIN , ZSTEQR ,    &
     &         ZUNGTR , ZUNMTR , ZUPGTR , ZUPMTR , ZHEEVD_2STAGE ,      &
     &         ZHEEVR_2STAGE , ZHEEVX_2STAGE , ZHEEV_2STAGE ,           &
     &         ZHBEV_2STAGE , ZHBEVD_2STAGE , ZHBEVX_2STAGE ,           &
     &         ZHETRD_2STAGE
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
!        ZHETRD
!
         SRNamt = 'ZHETRD'
         INFot = 1
         CALL ZHETRD('/',0,a,1,d,e,tau,w,1,info)
         CALL CHKXER('ZHETRD',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZHETRD('U',-1,a,1,d,e,tau,w,1,info)
         CALL CHKXER('ZHETRD',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZHETRD('U',2,a,1,d,e,tau,w,1,info)
         CALL CHKXER('ZHETRD',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL ZHETRD('U',0,a,1,d,e,tau,w,0,info)
         CALL CHKXER('ZHETRD',INFot,NOUt,LERr,OK)
         nt = nt + 4
!
!        ZHETRD_2STAGE
!
         SRNamt = 'ZHETRD_2STAGE'
         INFot = 1
         CALL ZHETRD_2STAGE('/','U',0,a,1,d,e,tau,c,1,w,1,info)
         CALL CHKXER('ZHETRD_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 1
         CALL ZHETRD_2STAGE('H','U',0,a,1,d,e,tau,c,1,w,1,info)
         CALL CHKXER('ZHETRD_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZHETRD_2STAGE('N','/',0,a,1,d,e,tau,c,1,w,1,info)
         CALL CHKXER('ZHETRD_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZHETRD_2STAGE('N','U',-1,a,1,d,e,tau,c,1,w,1,info)
         CALL CHKXER('ZHETRD_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL ZHETRD_2STAGE('N','U',2,a,1,d,e,tau,c,1,w,1,info)
         CALL CHKXER('ZHETRD_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL ZHETRD_2STAGE('N','U',0,a,1,d,e,tau,c,0,w,1,info)
         CALL CHKXER('ZHETRD_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL ZHETRD_2STAGE('N','U',0,a,1,d,e,tau,c,1,w,0,info)
         CALL CHKXER('ZHETRD_2STAGE',INFot,NOUt,LERr,OK)
         nt = nt + 7
!
!        ZHETRD_HE2HB
!
         SRNamt = 'ZHETRD_HE2HB'
         INFot = 1
         CALL ZHETRD_HE2HB('/',0,0,a,1,c,1,tau,w,1,info)
         CALL CHKXER('ZHETRD_HE2HB',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZHETRD_HE2HB('U',-1,0,a,1,c,1,tau,w,1,info)
         CALL CHKXER('ZHETRD_HE2HB',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZHETRD_HE2HB('U',0,-1,a,1,c,1,tau,w,1,info)
         CALL CHKXER('ZHETRD_HE2HB',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL ZHETRD_HE2HB('U',2,0,a,1,c,1,tau,w,1,info)
         CALL CHKXER('ZHETRD_HE2HB',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL ZHETRD_HE2HB('U',0,2,a,1,c,1,tau,w,1,info)
         CALL CHKXER('ZHETRD_HE2HB',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL ZHETRD_HE2HB('U',0,0,a,1,c,1,tau,w,0,info)
         CALL CHKXER('ZHETRD_HE2HB',INFot,NOUt,LERr,OK)
         nt = nt + 6
!
!        ZHETRD_HB2ST
!
         SRNamt = 'ZHETRD_HB2ST'
         INFot = 1
         CALL ZHETRD_HB2ST('/','N','U',0,0,a,1,d,e,c,1,w,1,info)
         CALL CHKXER('ZHETRD_HB2ST',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZHETRD_HB2ST('Y','/','U',0,0,a,1,d,e,c,1,w,1,info)
         CALL CHKXER('ZHETRD_HB2ST',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZHETRD_HB2ST('Y','H','U',0,0,a,1,d,e,c,1,w,1,info)
         CALL CHKXER('ZHETRD_HB2ST',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZHETRD_HB2ST('Y','N','/',0,0,a,1,d,e,c,1,w,1,info)
         CALL CHKXER('ZHETRD_HB2ST',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZHETRD_HB2ST('Y','N','U',-1,0,a,1,d,e,c,1,w,1,info)
         CALL CHKXER('ZHETRD_HB2ST',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL ZHETRD_HB2ST('Y','N','U',0,-1,a,1,d,e,c,1,w,1,info)
         CALL CHKXER('ZHETRD_HB2ST',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL ZHETRD_HB2ST('Y','N','U',0,1,a,1,d,e,c,1,w,1,info)
         CALL CHKXER('ZHETRD_HB2ST',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL ZHETRD_HB2ST('Y','N','U',0,0,a,1,d,e,c,0,w,1,info)
         CALL CHKXER('ZHETRD_HB2ST',INFot,NOUt,LERr,OK)
         INFot = 13
         CALL ZHETRD_HB2ST('Y','N','U',0,0,a,1,d,e,c,1,w,0,info)
         CALL CHKXER('ZHETRD_HB2ST',INFot,NOUt,LERr,OK)
         nt = nt + 9
!
!        ZUNGTR
!
         SRNamt = 'ZUNGTR'
         INFot = 1
         CALL ZUNGTR('/',0,a,1,tau,w,1,info)
         CALL CHKXER('ZUNGTR',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZUNGTR('U',-1,a,1,tau,w,1,info)
         CALL CHKXER('ZUNGTR',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZUNGTR('U',2,a,1,tau,w,1,info)
         CALL CHKXER('ZUNGTR',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL ZUNGTR('U',3,a,3,tau,w,1,info)
         CALL CHKXER('ZUNGTR',INFot,NOUt,LERr,OK)
         nt = nt + 4
!
!        ZUNMTR
!
         SRNamt = 'ZUNMTR'
         INFot = 1
         CALL ZUNMTR('/','U','N',0,0,a,1,tau,c,1,w,1,info)
         CALL CHKXER('ZUNMTR',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZUNMTR('L','/','N',0,0,a,1,tau,c,1,w,1,info)
         CALL CHKXER('ZUNMTR',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZUNMTR('L','U','/',0,0,a,1,tau,c,1,w,1,info)
         CALL CHKXER('ZUNMTR',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZUNMTR('L','U','N',-1,0,a,1,tau,c,1,w,1,info)
         CALL CHKXER('ZUNMTR',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL ZUNMTR('L','U','N',0,-1,a,1,tau,c,1,w,1,info)
         CALL CHKXER('ZUNMTR',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL ZUNMTR('L','U','N',2,0,a,1,tau,c,2,w,1,info)
         CALL CHKXER('ZUNMTR',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL ZUNMTR('R','U','N',0,2,a,1,tau,c,1,w,1,info)
         CALL CHKXER('ZUNMTR',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL ZUNMTR('L','U','N',2,0,a,2,tau,c,1,w,1,info)
         CALL CHKXER('ZUNMTR',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL ZUNMTR('L','U','N',0,2,a,1,tau,c,1,w,1,info)
         CALL CHKXER('ZUNMTR',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL ZUNMTR('R','U','N',2,0,a,1,tau,c,2,w,1,info)
         CALL CHKXER('ZUNMTR',INFot,NOUt,LERr,OK)
         nt = nt + 10
!
!        ZHPTRD
!
         SRNamt = 'ZHPTRD'
         INFot = 1
         CALL ZHPTRD('/',0,a,d,e,tau,info)
         CALL CHKXER('ZHPTRD',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZHPTRD('U',-1,a,d,e,tau,info)
         CALL CHKXER('ZHPTRD',INFot,NOUt,LERr,OK)
         nt = nt + 2
!
!        ZUPGTR
!
         SRNamt = 'ZUPGTR'
         INFot = 1
         CALL ZUPGTR('/',0,a,tau,z,1,w,info)
         CALL CHKXER('ZUPGTR',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZUPGTR('U',-1,a,tau,z,1,w,info)
         CALL CHKXER('ZUPGTR',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL ZUPGTR('U',2,a,tau,z,1,w,info)
         CALL CHKXER('ZUPGTR',INFot,NOUt,LERr,OK)
         nt = nt + 3
!
!        ZUPMTR
!
         SRNamt = 'ZUPMTR'
         INFot = 1
         CALL ZUPMTR('/','U','N',0,0,a,tau,c,1,w,info)
         CALL CHKXER('ZUPMTR',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZUPMTR('L','/','N',0,0,a,tau,c,1,w,info)
         CALL CHKXER('ZUPMTR',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZUPMTR('L','U','/',0,0,a,tau,c,1,w,info)
         CALL CHKXER('ZUPMTR',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZUPMTR('L','U','N',-1,0,a,tau,c,1,w,info)
         CALL CHKXER('ZUPMTR',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL ZUPMTR('L','U','N',0,-1,a,tau,c,1,w,info)
         CALL CHKXER('ZUPMTR',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL ZUPMTR('L','U','N',2,0,a,tau,c,1,w,info)
         CALL CHKXER('ZUPMTR',INFot,NOUt,LERr,OK)
         nt = nt + 6
!
!        ZPTEQR
!
         SRNamt = 'ZPTEQR'
         INFot = 1
         CALL ZPTEQR('/',0,d,e,z,1,rw,info)
         CALL CHKXER('ZPTEQR',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZPTEQR('N',-1,d,e,z,1,rw,info)
         CALL CHKXER('ZPTEQR',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL ZPTEQR('V',2,d,e,z,1,rw,info)
         CALL CHKXER('ZPTEQR',INFot,NOUt,LERr,OK)
         nt = nt + 3
!
!        ZSTEIN
!
         SRNamt = 'ZSTEIN'
         INFot = 1
         CALL ZSTEIN(-1,d,e,0,x,i1,i2,z,1,rw,iw,i3,info)
         CALL CHKXER('ZSTEIN',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZSTEIN(0,d,e,-1,x,i1,i2,z,1,rw,iw,i3,info)
         CALL CHKXER('ZSTEIN',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZSTEIN(0,d,e,1,x,i1,i2,z,1,rw,iw,i3,info)
         CALL CHKXER('ZSTEIN',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL ZSTEIN(2,d,e,0,x,i1,i2,z,1,rw,iw,i3,info)
         CALL CHKXER('ZSTEIN',INFot,NOUt,LERr,OK)
         nt = nt + 4
!
!        ZSTEQR
!
         SRNamt = 'ZSTEQR'
         INFot = 1
         CALL ZSTEQR('/',0,d,e,z,1,rw,info)
         CALL CHKXER('ZSTEQR',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZSTEQR('N',-1,d,e,z,1,rw,info)
         CALL CHKXER('ZSTEQR',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL ZSTEQR('V',2,d,e,z,1,rw,info)
         CALL CHKXER('ZSTEQR',INFot,NOUt,LERr,OK)
         nt = nt + 3
!
!        ZSTEDC
!
         SRNamt = 'ZSTEDC'
         INFot = 1
         CALL ZSTEDC('/',0,d,e,z,1,w,1,rw,1,iw,1,info)
         CALL CHKXER('ZSTEDC',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZSTEDC('N',-1,d,e,z,1,w,1,rw,1,iw,1,info)
         CALL CHKXER('ZSTEDC',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL ZSTEDC('V',2,d,e,z,1,w,4,rw,23,iw,28,info)
         CALL CHKXER('ZSTEDC',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL ZSTEDC('N',2,d,e,z,1,w,0,rw,1,iw,1,info)
         CALL CHKXER('ZSTEDC',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL ZSTEDC('V',2,d,e,z,2,w,0,rw,23,iw,28,info)
         CALL CHKXER('ZSTEDC',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL ZSTEDC('N',2,d,e,z,1,w,1,rw,0,iw,1,info)
         CALL CHKXER('ZSTEDC',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL ZSTEDC('I',2,d,e,z,2,w,1,rw,1,iw,12,info)
         CALL CHKXER('ZSTEDC',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL ZSTEDC('V',2,d,e,z,2,w,4,rw,1,iw,28,info)
         CALL CHKXER('ZSTEDC',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL ZSTEDC('N',2,d,e,z,1,w,1,rw,1,iw,0,info)
         CALL CHKXER('ZSTEDC',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL ZSTEDC('I',2,d,e,z,2,w,1,rw,23,iw,0,info)
         CALL CHKXER('ZSTEDC',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL ZSTEDC('V',2,d,e,z,2,w,4,rw,23,iw,0,info)
         CALL CHKXER('ZSTEDC',INFot,NOUt,LERr,OK)
         nt = nt + 11
!
!        ZHEEVD
!
         SRNamt = 'ZHEEVD'
         INFot = 1
         CALL ZHEEVD('/','U',0,a,1,x,w,1,rw,1,iw,1,info)
         CALL CHKXER('ZHEEVD',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZHEEVD('N','/',0,a,1,x,w,1,rw,1,iw,1,info)
         CALL CHKXER('ZHEEVD',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZHEEVD('N','U',-1,a,1,x,w,1,rw,1,iw,1,info)
         CALL CHKXER('ZHEEVD',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL ZHEEVD('N','U',2,a,1,x,w,3,rw,2,iw,1,info)
         CALL CHKXER('ZHEEVD',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL ZHEEVD('N','U',1,a,1,x,w,0,rw,1,iw,1,info)
         CALL CHKXER('ZHEEVD',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL ZHEEVD('N','U',2,a,2,x,w,2,rw,2,iw,1,info)
         CALL CHKXER('ZHEEVD',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL ZHEEVD('V','U',2,a,2,x,w,3,rw,25,iw,12,info)
         CALL CHKXER('ZHEEVD',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL ZHEEVD('N','U',1,a,1,x,w,1,rw,0,iw,1,info)
         CALL CHKXER('ZHEEVD',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL ZHEEVD('N','U',2,a,2,x,w,3,rw,1,iw,1,info)
         CALL CHKXER('ZHEEVD',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL ZHEEVD('V','U',2,a,2,x,w,8,rw,18,iw,12,info)
         CALL CHKXER('ZHEEVD',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL ZHEEVD('N','U',1,a,1,x,w,1,rw,1,iw,0,info)
         CALL CHKXER('ZHEEVD',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL ZHEEVD('V','U',2,a,2,x,w,8,rw,25,iw,11,info)
         CALL CHKXER('ZHEEVD',INFot,NOUt,LERr,OK)
         nt = nt + 12
!
!        ZHEEVD_2STAGE
!
         SRNamt = 'ZHEEVD_2STAGE'
         INFot = 1
         CALL ZHEEVD_2STAGE('/','U',0,a,1,x,w,1,rw,1,iw,1,info)
         CALL CHKXER('ZHEEVD_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 1
         CALL ZHEEVD_2STAGE('V','U',0,a,1,x,w,1,rw,1,iw,1,info)
         CALL CHKXER('ZHEEVD_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZHEEVD_2STAGE('N','/',0,a,1,x,w,1,rw,1,iw,1,info)
         CALL CHKXER('ZHEEVD_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZHEEVD_2STAGE('N','U',-1,a,1,x,w,1,rw,1,iw,1,info)
         CALL CHKXER('ZHEEVD_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL ZHEEVD_2STAGE('N','U',2,a,1,x,w,3,rw,2,iw,1,info)
         CALL CHKXER('ZHEEVD_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL ZHEEVD_2STAGE('N','U',1,a,1,x,w,0,rw,1,iw,1,info)
         CALL CHKXER('ZHEEVD_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL ZHEEVD_2STAGE('N','U',2,a,2,x,w,2,rw,2,iw,1,info)
         CALL CHKXER('ZHEEVD_2STAGE',INFot,NOUt,LERr,OK)
!         INFOT = 8
!         CALL ZHEEVD_2STAGE( 'V', 'U', 2, A, 2, X, W, 3,
!     $                            RW, 25, IW, 12, INFO )
!         CALL CHKXER( 'ZHEEVD_2STAGE', INFOT, NOUT, LERR, OK )
         INFot = 10
         CALL ZHEEVD_2STAGE('N','U',1,a,1,x,w,1,rw,0,iw,1,info)
         CALL CHKXER('ZHEEVD_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL ZHEEVD_2STAGE('N','U',2,a,2,x,w,25,rw,1,iw,1,info)
         CALL CHKXER('ZHEEVD_2STAGE',INFot,NOUt,LERr,OK)
!         INFOT = 10
!         CALL ZHEEVD_2STAGE( 'V', 'U', 2, A, 2, X, W, 8,
!     $                            RW, 18, IW, 12, INFO )
!         CALL CHKXER( 'ZHEEVD_2STAGE', INFOT, NOUT, LERR, OK )
         INFot = 12
         CALL ZHEEVD_2STAGE('N','U',1,a,1,x,w,1,rw,1,iw,0,info)
         CALL CHKXER('ZHEEVD_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 12
!         CALL ZHEEVD_2STAGE( 'V', 'U', 2, A, 2, X, W, 8,
!     $                            RW, 25, IW, 11, INFO )
!         CALL CHKXER( 'ZHEEVD_2STAGE', INFOT, NOUT, LERR, OK )
         nt = nt + 10
!
!        ZHEEV
!
         SRNamt = 'ZHEEV '
         INFot = 1
         CALL ZHEEV('/','U',0,a,1,x,w,1,rw,info)
         CALL CHKXER('ZHEEV ',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZHEEV('N','/',0,a,1,x,w,1,rw,info)
         CALL CHKXER('ZHEEV ',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZHEEV('N','U',-1,a,1,x,w,1,rw,info)
         CALL CHKXER('ZHEEV ',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL ZHEEV('N','U',2,a,1,x,w,3,rw,info)
         CALL CHKXER('ZHEEV ',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL ZHEEV('N','U',2,a,2,x,w,2,rw,info)
         CALL CHKXER('ZHEEV ',INFot,NOUt,LERr,OK)
         nt = nt + 5
!
!        ZHEEV_2STAGE
!
         SRNamt = 'ZHEEV_2STAGE '
         INFot = 1
         CALL ZHEEV_2STAGE('/','U',0,a,1,x,w,1,rw,info)
         CALL CHKXER('ZHEEV_2STAGE ',INFot,NOUt,LERr,OK)
         INFot = 1
         CALL ZHEEV_2STAGE('V','U',0,a,1,x,w,1,rw,info)
         CALL CHKXER('ZHEEV_2STAGE ',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZHEEV_2STAGE('N','/',0,a,1,x,w,1,rw,info)
         CALL CHKXER('ZHEEV_2STAGE ',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZHEEV_2STAGE('N','U',-1,a,1,x,w,1,rw,info)
         CALL CHKXER('ZHEEV_2STAGE ',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL ZHEEV_2STAGE('N','U',2,a,1,x,w,3,rw,info)
         CALL CHKXER('ZHEEV_2STAGE ',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL ZHEEV_2STAGE('N','U',2,a,2,x,w,2,rw,info)
         CALL CHKXER('ZHEEV_2STAGE ',INFot,NOUt,LERr,OK)
         nt = nt + 6
!
!        ZHEEVX
!
         SRNamt = 'ZHEEVX'
         INFot = 1
         CALL ZHEEVX('/','A','U',0,a,1,0.0D0,0.0D0,0,0,0.0D0,m,x,z,1,w, &
     &               1,rw,iw,i3,info)
         CALL CHKXER('ZHEEVX',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZHEEVX('V','/','U',0,a,1,0.0D0,1.0D0,1,0,0.0D0,m,x,z,1,w, &
     &               1,rw,iw,i3,info)
         CALL CHKXER('ZHEEVX',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZHEEVX('V','A','/',0,a,1,0.0D0,0.0D0,0,0,0.0D0,m,x,z,1,w, &
     &               1,rw,iw,i3,info)
         INFot = 4
         CALL ZHEEVX('V','A','U',-1,a,1,0.0D0,0.0D0,0,0,0.0D0,m,x,z,1,w,&
     &               1,rw,iw,i3,info)
         CALL CHKXER('ZHEEVX',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL ZHEEVX('V','A','U',2,a,1,0.0D0,0.0D0,0,0,0.0D0,m,x,z,2,w, &
     &               3,rw,iw,i3,info)
         CALL CHKXER('ZHEEVX',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL ZHEEVX('V','V','U',1,a,1,0.0D0,0.0D0,0,0,0.0D0,m,x,z,1,w, &
     &               1,rw,iw,i3,info)
         CALL CHKXER('ZHEEVX',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL ZHEEVX('V','I','U',1,a,1,0.0D0,0.0D0,0,0,0.0D0,m,x,z,1,w, &
     &               1,rw,iw,i3,info)
         CALL CHKXER('ZHEEVX',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL ZHEEVX('V','I','U',2,a,2,0.0D0,0.0D0,2,1,0.0D0,m,x,z,2,w, &
     &               3,rw,iw,i3,info)
         CALL CHKXER('ZHEEVX',INFot,NOUt,LERr,OK)
         INFot = 15
         CALL ZHEEVX('V','A','U',2,a,2,0.0D0,0.0D0,0,0,0.0D0,m,x,z,1,w, &
     &               3,rw,iw,i3,info)
         CALL CHKXER('ZHEEVX',INFot,NOUt,LERr,OK)
         INFot = 17
         CALL ZHEEVX('V','A','U',2,a,2,0.0D0,0.0D0,0,0,0.0D0,m,x,z,2,w, &
     &               2,rw,iw,i1,info)
         CALL CHKXER('ZHEEVX',INFot,NOUt,LERr,OK)
         nt = nt + 10
!
!        ZHEEVX_2STAGE
!
         SRNamt = 'ZHEEVX_2STAGE'
         INFot = 1
         CALL ZHEEVX_2STAGE('/','A','U',0,a,1,0.0D0,0.0D0,0,0,0.0D0,m,x,&
     &                      z,1,w,1,rw,iw,i3,info)
         CALL CHKXER('ZHEEVX_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 1
         CALL ZHEEVX_2STAGE('V','A','U',0,a,1,0.0D0,0.0D0,0,0,0.0D0,m,x,&
     &                      z,1,w,1,rw,iw,i3,info)
         CALL CHKXER('ZHEEVX_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZHEEVX_2STAGE('N','/','U',0,a,1,0.0D0,1.0D0,1,0,0.0D0,m,x,&
     &                      z,1,w,1,rw,iw,i3,info)
         CALL CHKXER('ZHEEVX_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZHEEVX_2STAGE('N','A','/',0,a,1,0.0D0,0.0D0,0,0,0.0D0,m,x,&
     &                      z,1,w,1,rw,iw,i3,info)
         INFot = 4
         CALL ZHEEVX_2STAGE('N','A','U',-1,a,1,0.0D0,0.0D0,0,0,0.0D0,m, &
     &                      x,z,1,w,1,rw,iw,i3,info)
         CALL CHKXER('ZHEEVX_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL ZHEEVX_2STAGE('N','A','U',2,a,1,0.0D0,0.0D0,0,0,0.0D0,m,x,&
     &                      z,2,w,3,rw,iw,i3,info)
         CALL CHKXER('ZHEEVX_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL ZHEEVX_2STAGE('N','V','U',1,a,1,0.0D0,0.0D0,0,0,0.0D0,m,x,&
     &                      z,1,w,1,rw,iw,i3,info)
         CALL CHKXER('ZHEEVX_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL ZHEEVX_2STAGE('N','I','U',1,a,1,0.0D0,0.0D0,0,0,0.0D0,m,x,&
     &                      z,1,w,1,rw,iw,i3,info)
         CALL CHKXER('ZHEEVX_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL ZHEEVX_2STAGE('N','I','U',2,a,2,0.0D0,0.0D0,2,1,0.0D0,m,x,&
     &                      z,2,w,3,rw,iw,i3,info)
         CALL CHKXER('ZHEEVX_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 15
         CALL ZHEEVX_2STAGE('N','A','U',2,a,2,0.0D0,0.0D0,0,0,0.0D0,m,x,&
     &                      z,0,w,3,rw,iw,i3,info)
         CALL CHKXER('ZHEEVX_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 17
         CALL ZHEEVX_2STAGE('N','A','U',2,a,2,0.0D0,0.0D0,0,0,0.0D0,m,x,&
     &                      z,2,w,0,rw,iw,i1,info)
         CALL CHKXER('ZHEEVX_2STAGE',INFot,NOUt,LERr,OK)
         nt = nt + 11
!
!        ZHEEVR
!
         SRNamt = 'ZHEEVR'
         n = 1
         INFot = 1
         CALL ZHEEVR('/','A','U',0,a,1,0.0D0,0.0D0,1,1,0.0D0,m,r,z,1,iw,&
     &               q,2*n,rw,24*n,iw(2*n+1),10*n,info)
         CALL CHKXER('ZHEEVR',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZHEEVR('V','/','U',0,a,1,0.0D0,0.0D0,1,1,0.0D0,m,r,z,1,iw,&
     &               q,2*n,rw,24*n,iw(2*n+1),10*n,info)
         CALL CHKXER('ZHEEVR',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZHEEVR('V','A','/',-1,a,1,0.0D0,0.0D0,1,1,0.0D0,m,r,z,1,  &
     &               iw,q,2*n,rw,24*n,iw(2*n+1),10*n,info)
         CALL CHKXER('ZHEEVR',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZHEEVR('V','A','U',-1,a,1,0.0D0,0.0D0,1,1,0.0D0,m,r,z,1,  &
     &               iw,q,2*n,rw,24*n,iw(2*n+1),10*n,info)
         CALL CHKXER('ZHEEVR',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL ZHEEVR('V','A','U',2,a,1,0.0D0,0.0D0,1,1,0.0D0,m,r,z,1,iw,&
     &               q,2*n,rw,24*n,iw(2*n+1),10*n,info)
         CALL CHKXER('ZHEEVR',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL ZHEEVR('V','V','U',1,a,1,0.0D0,0.0D0,1,1,0.0D0,m,r,z,1,iw,&
     &               q,2*n,rw,24*n,iw(2*n+1),10*n,info)
         CALL CHKXER('ZHEEVR',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL ZHEEVR('V','I','U',1,a,1,0.0D0,0.0D0,0,1,0.0D0,m,r,z,1,iw,&
     &               q,2*n,rw,24*n,iw(2*n+1),10*n,info)
         CALL CHKXER('ZHEEVR',INFot,NOUt,LERr,OK)
         INFot = 10
!
         CALL ZHEEVR('V','I','U',2,a,2,0.0D0,0.0D0,2,1,0.0D0,m,r,z,1,iw,&
     &               q,2*n,rw,24*n,iw(2*n+1),10*n,info)
         CALL CHKXER('ZHEEVR',INFot,NOUt,LERr,OK)
         INFot = 15
         CALL ZHEEVR('V','I','U',1,a,1,0.0D0,0.0D0,1,1,0.0D0,m,r,z,0,iw,&
     &               q,2*n,rw,24*n,iw(2*n+1),10*n,info)
         CALL CHKXER('ZHEEVR',INFot,NOUt,LERr,OK)
         INFot = 18
         CALL ZHEEVR('V','I','U',1,a,1,0.0D0,0.0D0,1,1,0.0D0,m,r,z,1,iw,&
     &               q,2*n-1,rw,24*n,iw(2*n+1),10*n,info)
         CALL CHKXER('ZHEEVR',INFot,NOUt,LERr,OK)
         INFot = 20
         CALL ZHEEVR('V','I','U',1,a,1,0.0D0,0.0D0,1,1,0.0D0,m,r,z,1,iw,&
     &               q,2*n,rw,24*n-1,iw(2*n-1),10*n,info)
         CALL CHKXER('ZHEEVR',INFot,NOUt,LERr,OK)
         INFot = 22
         CALL ZHEEVR('V','I','U',1,a,1,0.0D0,0.0D0,1,1,0.0D0,m,r,z,1,iw,&
     &               q,2*n,rw,24*n,iw,10*n-1,info)
         CALL CHKXER('ZHEEVR',INFot,NOUt,LERr,OK)
         nt = nt + 12
!
!        ZHEEVR_2STAGE
!
         SRNamt = 'ZHEEVR_2STAGE'
         n = 1
         INFot = 1
         CALL ZHEEVR_2STAGE('/','A','U',0,a,1,0.0D0,0.0D0,1,1,0.0D0,m,r,&
     &                      z,1,iw,q,2*n,rw,24*n,iw(2*n+1),10*n,info)
         CALL CHKXER('ZHEEVR_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 1
         CALL ZHEEVR_2STAGE('V','A','U',0,a,1,0.0D0,0.0D0,1,1,0.0D0,m,r,&
     &                      z,1,iw,q,2*n,rw,24*n,iw(2*n+1),10*n,info)
         CALL CHKXER('ZHEEVR_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZHEEVR_2STAGE('N','/','U',0,a,1,0.0D0,0.0D0,1,1,0.0D0,m,r,&
     &                      z,1,iw,q,2*n,rw,24*n,iw(2*n+1),10*n,info)
         CALL CHKXER('ZHEEVR_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZHEEVR_2STAGE('N','A','/',-1,a,1,0.0D0,0.0D0,1,1,0.0D0,m, &
     &                      r,z,1,iw,q,2*n,rw,24*n,iw(2*n+1),10*n,info)
         CALL CHKXER('ZHEEVR_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZHEEVR_2STAGE('N','A','U',-1,a,1,0.0D0,0.0D0,1,1,0.0D0,m, &
     &                      r,z,1,iw,q,2*n,rw,24*n,iw(2*n+1),10*n,info)
         CALL CHKXER('ZHEEVR_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL ZHEEVR_2STAGE('N','A','U',2,a,1,0.0D0,0.0D0,1,1,0.0D0,m,r,&
     &                      z,1,iw,q,2*n,rw,24*n,iw(2*n+1),10*n,info)
         CALL CHKXER('ZHEEVR_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL ZHEEVR_2STAGE('N','V','U',1,a,1,0.0D0,0.0D0,1,1,0.0D0,m,r,&
     &                      z,1,iw,q,2*n,rw,24*n,iw(2*n+1),10*n,info)
         CALL CHKXER('ZHEEVR_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL ZHEEVR_2STAGE('N','I','U',1,a,1,0.0D0,0.0D0,0,1,0.0D0,m,r,&
     &                      z,1,iw,q,2*n,rw,24*n,iw(2*n+1),10*n,info)
         CALL CHKXER('ZHEEVR_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL ZHEEVR_2STAGE('N','I','U',2,a,2,0.0D0,0.0D0,2,1,0.0D0,m,r,&
     &                      z,1,iw,q,2*n,rw,24*n,iw(2*n+1),10*n,info)
         CALL CHKXER('ZHEEVR_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 15
         CALL ZHEEVR_2STAGE('N','I','U',1,a,1,0.0D0,0.0D0,1,1,0.0D0,m,r,&
     &                      z,0,iw,q,2*n,rw,24*n,iw(2*n+1),10*n,info)
         CALL CHKXER('ZHEEVR_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 18
         CALL ZHEEVR_2STAGE('N','I','U',1,a,1,0.0D0,0.0D0,1,1,0.0D0,m,r,&
     &                      z,1,iw,q,2*n-1,rw,24*n,iw(2*n+1),10*n,info)
         CALL CHKXER('ZHEEVR_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 20
         CALL ZHEEVR_2STAGE('N','I','U',1,a,1,0.0D0,0.0D0,1,1,0.0D0,m,r,&
     &                      z,1,iw,q,26*n,rw,24*n-1,iw(2*n-1),10*n,info)
         CALL CHKXER('ZHEEVR_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 22
         CALL ZHEEVR_2STAGE('N','I','U',1,a,1,0.0D0,0.0D0,1,1,0.0D0,m,r,&
     &                      z,1,iw,q,26*n,rw,24*n,iw,10*n-1,info)
         CALL CHKXER('ZHEEVR_2STAGE',INFot,NOUt,LERr,OK)
         nt = nt + 13
!
!        ZHPEVD
!
         SRNamt = 'ZHPEVD'
         INFot = 1
         CALL ZHPEVD('/','U',0,a,x,z,1,w,1,rw,1,iw,1,info)
         CALL CHKXER('ZHPEVD',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZHPEVD('N','/',0,a,x,z,1,w,1,rw,1,iw,1,info)
         CALL CHKXER('ZHPEVD',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZHPEVD('N','U',-1,a,x,z,1,w,1,rw,1,iw,1,info)
         CALL CHKXER('ZHPEVD',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL ZHPEVD('V','U',2,a,x,z,1,w,4,rw,25,iw,12,info)
         CALL CHKXER('ZHPEVD',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL ZHPEVD('N','U',1,a,x,z,1,w,0,rw,1,iw,1,info)
         CALL CHKXER('ZHPEVD',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL ZHPEVD('N','U',2,a,x,z,2,w,1,rw,2,iw,1,info)
         CALL CHKXER('ZHPEVD',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL ZHPEVD('V','U',2,a,x,z,2,w,2,rw,25,iw,12,info)
         CALL CHKXER('ZHPEVD',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL ZHPEVD('N','U',1,a,x,z,1,w,1,rw,0,iw,1,info)
         CALL CHKXER('ZHPEVD',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL ZHPEVD('N','U',2,a,x,z,2,w,2,rw,1,iw,1,info)
         CALL CHKXER('ZHPEVD',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL ZHPEVD('V','U',2,a,x,z,2,w,4,rw,18,iw,12,info)
         CALL CHKXER('ZHPEVD',INFot,NOUt,LERr,OK)
         INFot = 13
         CALL ZHPEVD('N','U',1,a,x,z,1,w,1,rw,1,iw,0,info)
         CALL CHKXER('ZHPEVD',INFot,NOUt,LERr,OK)
         INFot = 13
         CALL ZHPEVD('N','U',2,a,x,z,2,w,2,rw,2,iw,0,info)
         CALL CHKXER('ZHPEVD',INFot,NOUt,LERr,OK)
         INFot = 13
         CALL ZHPEVD('V','U',2,a,x,z,2,w,4,rw,25,iw,2,info)
         CALL CHKXER('ZHPEVD',INFot,NOUt,LERr,OK)
         nt = nt + 13
!
!        ZHPEV
!
         SRNamt = 'ZHPEV '
         INFot = 1
         CALL ZHPEV('/','U',0,a,x,z,1,w,rw,info)
         CALL CHKXER('ZHPEV ',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZHPEV('N','/',0,a,x,z,1,w,rw,info)
         CALL CHKXER('ZHPEV ',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZHPEV('N','U',-1,a,x,z,1,w,rw,info)
         CALL CHKXER('ZHPEV ',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL ZHPEV('V','U',2,a,x,z,1,w,rw,info)
         CALL CHKXER('ZHPEV ',INFot,NOUt,LERr,OK)
         nt = nt + 4
!
!        ZHPEVX
!
         SRNamt = 'ZHPEVX'
         INFot = 1
         CALL ZHPEVX('/','A','U',0,a,0.0D0,0.0D0,0,0,0.0D0,m,x,z,1,w,rw,&
     &               iw,i3,info)
         CALL CHKXER('ZHPEVX',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZHPEVX('V','/','U',0,a,0.0D0,1.0D0,1,0,0.0D0,m,x,z,1,w,rw,&
     &               iw,i3,info)
         CALL CHKXER('ZHPEVX',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZHPEVX('V','A','/',0,a,0.0D0,0.0D0,0,0,0.0D0,m,x,z,1,w,rw,&
     &               iw,i3,info)
         CALL CHKXER('ZHPEVX',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZHPEVX('V','A','U',-1,a,0.0D0,0.0D0,0,0,0.0D0,m,x,z,1,w,  &
     &               rw,iw,i3,info)
         CALL CHKXER('ZHPEVX',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL ZHPEVX('V','V','U',1,a,0.0D0,0.0D0,0,0,0.0D0,m,x,z,1,w,rw,&
     &               iw,i3,info)
         CALL CHKXER('ZHPEVX',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL ZHPEVX('V','I','U',1,a,0.0D0,0.0D0,0,0,0.0D0,m,x,z,1,w,rw,&
     &               iw,i3,info)
         CALL CHKXER('ZHPEVX',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL ZHPEVX('V','I','U',2,a,0.0D0,0.0D0,2,1,0.0D0,m,x,z,2,w,rw,&
     &               iw,i3,info)
         CALL CHKXER('ZHPEVX',INFot,NOUt,LERr,OK)
         INFot = 14
         CALL ZHPEVX('V','A','U',2,a,0.0D0,0.0D0,0,0,0.0D0,m,x,z,1,w,rw,&
     &               iw,i3,info)
         CALL CHKXER('ZHPEVX',INFot,NOUt,LERr,OK)
         nt = nt + 8
!
!     Test error exits for the HB path.
!
      ELSEIF ( LSAMEN(2,c2,'HB') ) THEN
!
!        ZHBTRD
!
         SRNamt = 'ZHBTRD'
         INFot = 1
         CALL ZHBTRD('/','U',0,0,a,1,d,e,z,1,w,info)
         CALL CHKXER('ZHBTRD',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZHBTRD('N','/',0,0,a,1,d,e,z,1,w,info)
         CALL CHKXER('ZHBTRD',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZHBTRD('N','U',-1,0,a,1,d,e,z,1,w,info)
         CALL CHKXER('ZHBTRD',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZHBTRD('N','U',0,-1,a,1,d,e,z,1,w,info)
         CALL CHKXER('ZHBTRD',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL ZHBTRD('N','U',1,1,a,1,d,e,z,1,w,info)
         CALL CHKXER('ZHBTRD',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL ZHBTRD('V','U',2,0,a,1,d,e,z,1,w,info)
         CALL CHKXER('ZHBTRD',INFot,NOUt,LERr,OK)
         nt = nt + 6
!
!        ZHETRD_HB2ST
!
         SRNamt = 'ZHETRD_HB2ST'
         INFot = 1
         CALL ZHETRD_HB2ST('/','N','U',0,0,a,1,d,e,c,1,w,1,info)
         CALL CHKXER('ZHETRD_HB2ST',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZHETRD_HB2ST('N','/','U',0,0,a,1,d,e,c,1,w,1,info)
         CALL CHKXER('ZHETRD_HB2ST',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZHETRD_HB2ST('N','H','U',0,0,a,1,d,e,c,1,w,1,info)
         CALL CHKXER('ZHETRD_HB2ST',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZHETRD_HB2ST('N','N','/',0,0,a,1,d,e,c,1,w,1,info)
         CALL CHKXER('ZHETRD_HB2ST',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZHETRD_HB2ST('N','N','U',-1,0,a,1,d,e,c,1,w,1,info)
         CALL CHKXER('ZHETRD_HB2ST',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL ZHETRD_HB2ST('N','N','U',0,-1,a,1,d,e,c,1,w,1,info)
         CALL CHKXER('ZHETRD_HB2ST',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL ZHETRD_HB2ST('N','N','U',0,1,a,1,d,e,c,1,w,1,info)
         CALL CHKXER('ZHETRD_HB2ST',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL ZHETRD_HB2ST('N','N','U',0,0,a,1,d,e,c,0,w,1,info)
         CALL CHKXER('ZHETRD_HB2ST',INFot,NOUt,LERr,OK)
         INFot = 13
         CALL ZHETRD_HB2ST('N','N','U',0,0,a,1,d,e,c,1,w,0,info)
         CALL CHKXER('ZHETRD_HB2ST',INFot,NOUt,LERr,OK)
         nt = nt + 9
!
!        ZHBEVD
!
         SRNamt = 'ZHBEVD'
         INFot = 1
         CALL ZHBEVD('/','U',0,0,a,1,x,z,1,w,1,rw,1,iw,1,info)
         CALL CHKXER('ZHBEVD',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZHBEVD('N','/',0,0,a,1,x,z,1,w,1,rw,1,iw,1,info)
         CALL CHKXER('ZHBEVD',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZHBEVD('N','U',-1,0,a,1,x,z,1,w,1,rw,1,iw,1,info)
         CALL CHKXER('ZHBEVD',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZHBEVD('N','U',0,-1,a,1,x,z,1,w,1,rw,1,iw,1,info)
         CALL CHKXER('ZHBEVD',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL ZHBEVD('N','U',2,1,a,1,x,z,1,w,2,rw,2,iw,1,info)
         CALL CHKXER('ZHBEVD',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL ZHBEVD('V','U',2,1,a,2,x,z,1,w,8,rw,25,iw,12,info)
         CALL CHKXER('ZHBEVD',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL ZHBEVD('N','U',1,0,a,1,x,z,1,w,0,rw,1,iw,1,info)
         CALL CHKXER('ZHBEVD',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL ZHBEVD('N','U',2,1,a,2,x,z,2,w,1,rw,2,iw,1,info)
         CALL CHKXER('ZHBEVD',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL ZHBEVD('V','U',2,1,a,2,x,z,2,w,2,rw,25,iw,12,info)
         CALL CHKXER('ZHBEVD',INFot,NOUt,LERr,OK)
         INFot = 13
         CALL ZHBEVD('N','U',1,0,a,1,x,z,1,w,1,rw,0,iw,1,info)
         CALL CHKXER('ZHBEVD',INFot,NOUt,LERr,OK)
         INFot = 13
         CALL ZHBEVD('N','U',2,1,a,2,x,z,2,w,2,rw,1,iw,1,info)
         CALL CHKXER('ZHBEVD',INFot,NOUt,LERr,OK)
         INFot = 13
         CALL ZHBEVD('V','U',2,1,a,2,x,z,2,w,8,rw,2,iw,12,info)
         CALL CHKXER('ZHBEVD',INFot,NOUt,LERr,OK)
         INFot = 15
         CALL ZHBEVD('N','U',1,0,a,1,x,z,1,w,1,rw,1,iw,0,info)
         CALL CHKXER('ZHBEVD',INFot,NOUt,LERr,OK)
         INFot = 15
         CALL ZHBEVD('N','U',2,1,a,2,x,z,2,w,2,rw,2,iw,0,info)
         CALL CHKXER('ZHBEVD',INFot,NOUt,LERr,OK)
         INFot = 15
         CALL ZHBEVD('V','U',2,1,a,2,x,z,2,w,8,rw,25,iw,2,info)
         CALL CHKXER('ZHBEVD',INFot,NOUt,LERr,OK)
         nt = nt + 15
!
!        ZHBEVD_2STAGE
!
         SRNamt = 'ZHBEVD_2STAGE'
         INFot = 1
         CALL ZHBEVD_2STAGE('/','U',0,0,a,1,x,z,1,w,1,rw,1,iw,1,info)
         CALL CHKXER('ZHBEVD_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 1
         CALL ZHBEVD_2STAGE('V','U',0,0,a,1,x,z,1,w,1,rw,1,iw,1,info)
         CALL CHKXER('ZHBEVD_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZHBEVD_2STAGE('N','/',0,0,a,1,x,z,1,w,1,rw,1,iw,1,info)
         CALL CHKXER('ZHBEVD_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZHBEVD_2STAGE('N','U',-1,0,a,1,x,z,1,w,1,rw,1,iw,1,info)
         CALL CHKXER('ZHBEVD_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZHBEVD_2STAGE('N','U',0,-1,a,1,x,z,1,w,1,rw,1,iw,1,info)
         CALL CHKXER('ZHBEVD_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL ZHBEVD_2STAGE('N','U',2,1,a,1,x,z,1,w,2,rw,2,iw,1,info)
         CALL CHKXER('ZHBEVD_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL ZHBEVD_2STAGE('N','U',2,1,a,2,x,z,0,w,8,rw,25,iw,12,info)
         CALL CHKXER('ZHBEVD_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL ZHBEVD_2STAGE('N','U',1,0,a,1,x,z,1,w,0,rw,1,iw,1,info)
         CALL CHKXER('ZHBEVD_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL ZHBEVD_2STAGE('N','U',2,1,a,2,x,z,2,w,1,rw,2,iw,1,info)
         CALL CHKXER('ZHBEVD_2STAGE',INFot,NOUt,LERr,OK)
!         INFOT = 11
!         CALL ZHBEVD_2STAGE( 'V', 'U', 2, 1, A, 2, X, Z, 2,
!     $                         W, 2, RW, 25, IW, 12, INFO )
!         CALL CHKXER( 'ZHBEVD_2STAGE', INFOT, NOUT, LERR, OK )
         INFot = 13
         CALL ZHBEVD_2STAGE('N','U',1,0,a,1,x,z,1,w,1,rw,0,iw,1,info)
         CALL CHKXER('ZHBEVD_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 13
         CALL ZHBEVD_2STAGE('N','U',2,1,a,2,x,z,2,w,25,rw,1,iw,1,info)
         CALL CHKXER('ZHBEVD_2STAGE',INFot,NOUt,LERr,OK)
!         INFOT = 13
!         CALL ZHBEVD_2STAGE( 'V', 'U', 2, 1, A, 2, X, Z, 2,
!     $                          W, 25, RW, 2, IW, 12, INFO )
!         CALL CHKXER( 'ZHBEVD_2STAGE', INFOT, NOUT, LERR, OK )
         INFot = 15
         CALL ZHBEVD_2STAGE('N','U',1,0,a,1,x,z,1,w,1,rw,1,iw,0,info)
         CALL CHKXER('ZHBEVD_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 15
         CALL ZHBEVD_2STAGE('N','U',2,1,a,2,x,z,2,w,25,rw,2,iw,0,info)
         CALL CHKXER('ZHBEVD_2STAGE',INFot,NOUt,LERr,OK)
!         INFOT = 15
!         CALL ZHBEVD_2STAGE( 'V', 'U', 2, 1, A, 2, X, Z, 2,
!     $                          W, 25, RW, 25, IW, 2, INFO )
!         CALL CHKXER( 'ZHBEVD_2STAGE', INFOT, NOUT, LERR, OK )
         nt = nt + 13
!
!        ZHBEV
!
         SRNamt = 'ZHBEV '
         INFot = 1
         CALL ZHBEV('/','U',0,0,a,1,x,z,1,w,rw,info)
         CALL CHKXER('ZHBEV ',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZHBEV('N','/',0,0,a,1,x,z,1,w,rw,info)
         CALL CHKXER('ZHBEV ',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZHBEV('N','U',-1,0,a,1,x,z,1,w,rw,info)
         CALL CHKXER('ZHBEV ',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZHBEV('N','U',0,-1,a,1,x,z,1,w,rw,info)
         CALL CHKXER('ZHBEV ',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL ZHBEV('N','U',2,1,a,1,x,z,1,w,rw,info)
         CALL CHKXER('ZHBEV ',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL ZHBEV('V','U',2,0,a,1,x,z,1,w,rw,info)
         CALL CHKXER('ZHBEV ',INFot,NOUt,LERr,OK)
         nt = nt + 6
!
!        ZHBEV_2STAGE
!
         SRNamt = 'ZHBEV_2STAGE '
         INFot = 1
         CALL ZHBEV_2STAGE('/','U',0,0,a,1,x,z,1,w,0,rw,info)
         CALL CHKXER('ZHBEV_2STAGE ',INFot,NOUt,LERr,OK)
         INFot = 1
         CALL ZHBEV_2STAGE('V','U',0,0,a,1,x,z,1,w,0,rw,info)
         CALL CHKXER('ZHBEV_2STAGE ',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZHBEV_2STAGE('N','/',0,0,a,1,x,z,1,w,0,rw,info)
         CALL CHKXER('ZHBEV_2STAGE ',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZHBEV_2STAGE('N','U',-1,0,a,1,x,z,1,w,0,rw,info)
         CALL CHKXER('ZHBEV_2STAGE ',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZHBEV_2STAGE('N','U',0,-1,a,1,x,z,1,w,0,rw,info)
         CALL CHKXER('ZHBEV_2STAGE ',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL ZHBEV_2STAGE('N','U',2,1,a,1,x,z,1,w,0,rw,info)
         CALL CHKXER('ZHBEV_2STAGE ',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL ZHBEV_2STAGE('N','U',2,0,a,1,x,z,0,w,0,rw,info)
         CALL CHKXER('ZHBEV_2STAGE ',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL ZHBEV_2STAGE('N','U',2,0,a,1,x,z,1,w,0,rw,info)
         CALL CHKXER('ZHBEV_2STAGE ',INFot,NOUt,LERr,OK)
         nt = nt + 8
!
!        ZHBEVX
!
         SRNamt = 'ZHBEVX'
         INFot = 1
         CALL ZHBEVX('/','A','U',0,0,a,1,q,1,0.0D0,0.0D0,0,0,0.0D0,m,x, &
     &               z,1,w,rw,iw,i3,info)
         CALL CHKXER('ZHBEVX',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZHBEVX('V','/','U',0,0,a,1,q,1,0.0D0,1.0D0,1,0,0.0D0,m,x, &
     &               z,1,w,rw,iw,i3,info)
         CALL CHKXER('ZHBEVX',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZHBEVX('V','A','/',0,0,a,1,q,1,0.0D0,0.0D0,0,0,0.0D0,m,x, &
     &               z,1,w,rw,iw,i3,info)
         INFot = 4
         CALL ZHBEVX('V','A','U',-1,0,a,1,q,1,0.0D0,0.0D0,0,0,0.0D0,m,x,&
     &               z,1,w,rw,iw,i3,info)
         CALL CHKXER('ZHBEVX',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL ZHBEVX('V','A','U',0,-1,a,1,q,1,0.0D0,0.0D0,0,0,0.0D0,m,x,&
     &               z,1,w,rw,iw,i3,info)
         CALL CHKXER('ZHBEVX',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL ZHBEVX('V','A','U',2,1,a,1,q,2,0.0D0,0.0D0,0,0,0.0D0,m,x, &
     &               z,2,w,rw,iw,i3,info)
         CALL CHKXER('ZHBEVX',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL ZHBEVX('V','A','U',2,0,a,1,q,1,0.0D0,0.0D0,0,0,0.0D0,m,x, &
     &               z,2,w,rw,iw,i3,info)
         CALL CHKXER('ZHBEVX',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL ZHBEVX('V','V','U',1,0,a,1,q,1,0.0D0,0.0D0,0,0,0.0D0,m,x, &
     &               z,1,w,rw,iw,i3,info)
         CALL CHKXER('ZHBEVX',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL ZHBEVX('V','I','U',1,0,a,1,q,1,0.0D0,0.0D0,0,0,0.0D0,m,x, &
     &               z,1,w,rw,iw,i3,info)
         CALL CHKXER('ZHBEVX',INFot,NOUt,LERr,OK)
         INFot = 13
         CALL ZHBEVX('V','I','U',1,0,a,1,q,1,0.0D0,0.0D0,1,2,0.0D0,m,x, &
     &               z,1,w,rw,iw,i3,info)
         CALL CHKXER('ZHBEVX',INFot,NOUt,LERr,OK)
         INFot = 18
         CALL ZHBEVX('V','A','U',2,0,a,1,q,2,0.0D0,0.0D0,0,0,0.0D0,m,x, &
     &               z,1,w,rw,iw,i3,info)
         CALL CHKXER('ZHBEVX',INFot,NOUt,LERr,OK)
         nt = nt + 11
!
!        ZHBEVX_2STAGE
!
         SRNamt = 'ZHBEVX_2STAGE'
         INFot = 1
         CALL ZHBEVX_2STAGE('/','A','U',0,0,a,1,q,1,0.0D0,0.0D0,0,0,    &
     &                      0.0D0,m,x,z,1,w,0,rw,iw,i3,info)
         INFot = 1
         CALL ZHBEVX_2STAGE('V','A','U',0,0,a,1,q,1,0.0D0,0.0D0,0,0,    &
     &                      0.0D0,m,x,z,1,w,0,rw,iw,i3,info)
         CALL CHKXER('ZHBEVX_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZHBEVX_2STAGE('N','/','U',0,0,a,1,q,1,0.0D0,1.0D0,1,0,    &
     &                      0.0D0,m,x,z,1,w,0,rw,iw,i3,info)
         CALL CHKXER('ZHBEVX_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZHBEVX_2STAGE('N','A','/',0,0,a,1,q,1,0.0D0,0.0D0,0,0,    &
     &                      0.0D0,m,x,z,1,w,0,rw,iw,i3,info)
         INFot = 4
         CALL ZHBEVX_2STAGE('N','A','U',-1,0,a,1,q,1,0.0D0,0.0D0,0,0,   &
     &                      0.0D0,m,x,z,1,w,0,rw,iw,i3,info)
         CALL CHKXER('ZHBEVX_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL ZHBEVX_2STAGE('N','A','U',0,-1,a,1,q,1,0.0D0,0.0D0,0,0,   &
     &                      0.0D0,m,x,z,1,w,0,rw,iw,i3,info)
         CALL CHKXER('ZHBEVX_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL ZHBEVX_2STAGE('N','A','U',2,1,a,1,q,2,0.0D0,0.0D0,0,0,    &
     &                      0.0D0,m,x,z,2,w,0,rw,iw,i3,info)
         CALL CHKXER('ZHBEVX_2STAGE',INFot,NOUt,LERr,OK)
!         INFOT = 9
!         CALL ZHBEVX_2STAGE( 'V', 'A', 'U', 2, 0, A, 1, Q, 1,
!     $                       0.0D0, 0.0D0, 0, 0, 0.0D0,
!     $                       M, X, Z, 2, W, 0, RW, IW, I3, INFO )
!         CALL CHKXER( 'ZHBEVX_2STAGE', INFOT, NOUT, LERR, OK )
         INFot = 11
         CALL ZHBEVX_2STAGE('N','V','U',1,0,a,1,q,1,0.0D0,0.0D0,0,0,    &
     &                      0.0D0,m,x,z,1,w,0,rw,iw,i3,info)
         CALL CHKXER('ZHBEVX_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL ZHBEVX_2STAGE('N','I','U',1,0,a,1,q,1,0.0D0,0.0D0,0,0,    &
     &                      0.0D0,m,x,z,1,w,0,rw,iw,i3,info)
         CALL CHKXER('ZHBEVX_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 13
         CALL ZHBEVX_2STAGE('N','I','U',1,0,a,1,q,1,0.0D0,0.0D0,1,2,    &
     &                      0.0D0,m,x,z,1,w,0,rw,iw,i3,info)
         CALL CHKXER('ZHBEVX_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 18
         CALL ZHBEVX_2STAGE('N','A','U',2,0,a,1,q,2,0.0D0,0.0D0,0,0,    &
     &                      0.0D0,m,x,z,0,w,0,rw,iw,i3,info)
         CALL CHKXER('ZHBEVX_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 20
         CALL ZHBEVX_2STAGE('N','A','U',2,0,a,1,q,2,0.0D0,0.0D0,0,0,    &
     &                      0.0D0,m,x,z,1,w,0,rw,iw,i3,info)
         CALL CHKXER('ZHBEVX_2STAGE',INFot,NOUt,LERr,OK)
         nt = nt + 12
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
!     End of ZERRST
!
      END SUBROUTINE ZERRST
