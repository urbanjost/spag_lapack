!*==cerrst.f90  processed by SPAG 7.51RB at 17:37 on  4 Mar 2022
!> \brief \b CERRST
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE CERRST( PATH, NUNIT )
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
!> CERRST tests the error exits for CHETRD, CUNGTR, CUNMTR, CHPTRD,
!> CUNGTR, CUPMTR, CSTEQR, CSTEIN, CPTEQR, CHBTRD,
!> CHEEV, CHEEVX, CHEEVD, CHBEV, CHBEVX, CHBEVD,
!> CHPEV, CHPEVX, CHPEVD, and CSTEDC.
!> CHEEVD_2STAGE, CHEEVR_2STAGE, CHEEVX_2STAGE,
!> CHEEV_2STAGE, CHBEV_2STAGE, CHBEVD_2STAGE,
!> CHBEVX_2STAGE, CHETRD_2STAGE, CHETRD_HE2HB,
!> CHETRD_HB2ST
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
!> \ingroup complex_eig
!
!  =====================================================================
      SUBROUTINE CERRST(Path,Nunit)
      IMPLICIT NONE
!*--CERRST65
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
      REAL d(NMAX) , e(NMAX) , r(LW) , rw(LW) , x(NMAX)
      COMPLEX a(NMAX,NMAX) , c(NMAX,NMAX) , q(NMAX,NMAX) , tau(NMAX) ,  &
     &        w(LW) , z(NMAX,NMAX)
!     ..
!     .. External Functions ..
      LOGICAL LSAMEN
      EXTERNAL LSAMEN
!     ..
!     .. External Subroutines ..
      EXTERNAL CHBEV , CHBEVD , CHBEVX , CHBTRD , CHEEV , CHEEVD ,      &
     &         CHEEVR , CHEEVX , CHETRD , CHKXER , CHPEV , CHPEVD ,     &
     &         CHPEVX , CHPTRD , CPTEQR , CSTEDC , CSTEIN , CSTEQR ,    &
     &         CUNGTR , CUNMTR , CUPGTR , CUPMTR , CHEEVD_2STAGE ,      &
     &         CHEEVR_2STAGE , CHEEVX_2STAGE , CHEEV_2STAGE ,           &
     &         CHBEV_2STAGE , CHBEVD_2STAGE , CHBEVX_2STAGE ,           &
     &         CHETRD_2STAGE , CHETRD_HE2HB , CHETRD_HB2ST
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
!        CHETRD
!
         SRNamt = 'CHETRD'
         INFot = 1
         CALL CHETRD('/',0,a,1,d,e,tau,w,1,info)
         CALL CHKXER('CHETRD',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CHETRD('U',-1,a,1,d,e,tau,w,1,info)
         CALL CHKXER('CHETRD',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CHETRD('U',2,a,1,d,e,tau,w,1,info)
         CALL CHKXER('CHETRD',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL CHETRD('U',0,a,1,d,e,tau,w,0,info)
         CALL CHKXER('CHETRD',INFot,NOUt,LERr,OK)
         nt = nt + 4
!
!        CHETRD_2STAGE
!
         SRNamt = 'CHETRD_2STAGE'
         INFot = 1
         CALL CHETRD_2STAGE('/','U',0,a,1,d,e,tau,c,1,w,1,info)
         CALL CHKXER('CHETRD_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 1
         CALL CHETRD_2STAGE('H','U',0,a,1,d,e,tau,c,1,w,1,info)
         CALL CHKXER('CHETRD_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CHETRD_2STAGE('N','/',0,a,1,d,e,tau,c,1,w,1,info)
         CALL CHKXER('CHETRD_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CHETRD_2STAGE('N','U',-1,a,1,d,e,tau,c,1,w,1,info)
         CALL CHKXER('CHETRD_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL CHETRD_2STAGE('N','U',2,a,1,d,e,tau,c,1,w,1,info)
         CALL CHKXER('CHETRD_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL CHETRD_2STAGE('N','U',0,a,1,d,e,tau,c,0,w,1,info)
         CALL CHKXER('CHETRD_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL CHETRD_2STAGE('N','U',0,a,1,d,e,tau,c,1,w,0,info)
         CALL CHKXER('CHETRD_2STAGE',INFot,NOUt,LERr,OK)
         nt = nt + 7
!
!        CHETRD_HE2HB
!
         SRNamt = 'CHETRD_HE2HB'
         INFot = 1
         CALL CHETRD_HE2HB('/',0,0,a,1,c,1,tau,w,1,info)
         CALL CHKXER('CHETRD_HE2HB',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CHETRD_HE2HB('U',-1,0,a,1,c,1,tau,w,1,info)
         CALL CHKXER('CHETRD_HE2HB',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CHETRD_HE2HB('U',0,-1,a,1,c,1,tau,w,1,info)
         CALL CHKXER('CHETRD_HE2HB',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL CHETRD_HE2HB('U',2,0,a,1,c,1,tau,w,1,info)
         CALL CHKXER('CHETRD_HE2HB',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL CHETRD_HE2HB('U',0,2,a,1,c,1,tau,w,1,info)
         CALL CHKXER('CHETRD_HE2HB',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL CHETRD_HE2HB('U',0,0,a,1,c,1,tau,w,0,info)
         CALL CHKXER('CHETRD_HE2HB',INFot,NOUt,LERr,OK)
         nt = nt + 6
!
!        CHETRD_HB2ST
!
         SRNamt = 'CHETRD_HB2ST'
         INFot = 1
         CALL CHETRD_HB2ST('/','N','U',0,0,a,1,d,e,c,1,w,1,info)
         CALL CHKXER('CHETRD_HB2ST',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CHETRD_HB2ST('Y','/','U',0,0,a,1,d,e,c,1,w,1,info)
         CALL CHKXER('CHETRD_HB2ST',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CHETRD_HB2ST('Y','H','U',0,0,a,1,d,e,c,1,w,1,info)
         CALL CHKXER('CHETRD_HB2ST',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CHETRD_HB2ST('Y','N','/',0,0,a,1,d,e,c,1,w,1,info)
         CALL CHKXER('CHETRD_HB2ST',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CHETRD_HB2ST('Y','N','U',-1,0,a,1,d,e,c,1,w,1,info)
         CALL CHKXER('CHETRD_HB2ST',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL CHETRD_HB2ST('Y','N','U',0,-1,a,1,d,e,c,1,w,1,info)
         CALL CHKXER('CHETRD_HB2ST',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL CHETRD_HB2ST('Y','N','U',0,1,a,1,d,e,c,1,w,1,info)
         CALL CHKXER('CHETRD_HB2ST',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL CHETRD_HB2ST('Y','N','U',0,0,a,1,d,e,c,0,w,1,info)
         CALL CHKXER('CHETRD_HB2ST',INFot,NOUt,LERr,OK)
         INFot = 13
         CALL CHETRD_HB2ST('Y','N','U',0,0,a,1,d,e,c,1,w,0,info)
         CALL CHKXER('CHETRD_HB2ST',INFot,NOUt,LERr,OK)
         nt = nt + 9
!
!        CUNGTR
!
         SRNamt = 'CUNGTR'
         INFot = 1
         CALL CUNGTR('/',0,a,1,tau,w,1,info)
         CALL CHKXER('CUNGTR',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CUNGTR('U',-1,a,1,tau,w,1,info)
         CALL CHKXER('CUNGTR',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CUNGTR('U',2,a,1,tau,w,1,info)
         CALL CHKXER('CUNGTR',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL CUNGTR('U',3,a,3,tau,w,1,info)
         CALL CHKXER('CUNGTR',INFot,NOUt,LERr,OK)
         nt = nt + 4
!
!        CUNMTR
!
         SRNamt = 'CUNMTR'
         INFot = 1
         CALL CUNMTR('/','U','N',0,0,a,1,tau,c,1,w,1,info)
         CALL CHKXER('CUNMTR',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CUNMTR('L','/','N',0,0,a,1,tau,c,1,w,1,info)
         CALL CHKXER('CUNMTR',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CUNMTR('L','U','/',0,0,a,1,tau,c,1,w,1,info)
         CALL CHKXER('CUNMTR',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CUNMTR('L','U','N',-1,0,a,1,tau,c,1,w,1,info)
         CALL CHKXER('CUNMTR',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL CUNMTR('L','U','N',0,-1,a,1,tau,c,1,w,1,info)
         CALL CHKXER('CUNMTR',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL CUNMTR('L','U','N',2,0,a,1,tau,c,2,w,1,info)
         CALL CHKXER('CUNMTR',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL CUNMTR('R','U','N',0,2,a,1,tau,c,1,w,1,info)
         CALL CHKXER('CUNMTR',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL CUNMTR('L','U','N',2,0,a,2,tau,c,1,w,1,info)
         CALL CHKXER('CUNMTR',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL CUNMTR('L','U','N',0,2,a,1,tau,c,1,w,1,info)
         CALL CHKXER('CUNMTR',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL CUNMTR('R','U','N',2,0,a,1,tau,c,2,w,1,info)
         CALL CHKXER('CUNMTR',INFot,NOUt,LERr,OK)
         nt = nt + 10
!
!        CHPTRD
!
         SRNamt = 'CHPTRD'
         INFot = 1
         CALL CHPTRD('/',0,a,d,e,tau,info)
         CALL CHKXER('CHPTRD',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CHPTRD('U',-1,a,d,e,tau,info)
         CALL CHKXER('CHPTRD',INFot,NOUt,LERr,OK)
         nt = nt + 2
!
!        CUPGTR
!
         SRNamt = 'CUPGTR'
         INFot = 1
         CALL CUPGTR('/',0,a,tau,z,1,w,info)
         CALL CHKXER('CUPGTR',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CUPGTR('U',-1,a,tau,z,1,w,info)
         CALL CHKXER('CUPGTR',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL CUPGTR('U',2,a,tau,z,1,w,info)
         CALL CHKXER('CUPGTR',INFot,NOUt,LERr,OK)
         nt = nt + 3
!
!        CUPMTR
!
         SRNamt = 'CUPMTR'
         INFot = 1
         CALL CUPMTR('/','U','N',0,0,a,tau,c,1,w,info)
         CALL CHKXER('CUPMTR',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CUPMTR('L','/','N',0,0,a,tau,c,1,w,info)
         CALL CHKXER('CUPMTR',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CUPMTR('L','U','/',0,0,a,tau,c,1,w,info)
         CALL CHKXER('CUPMTR',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CUPMTR('L','U','N',-1,0,a,tau,c,1,w,info)
         CALL CHKXER('CUPMTR',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL CUPMTR('L','U','N',0,-1,a,tau,c,1,w,info)
         CALL CHKXER('CUPMTR',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL CUPMTR('L','U','N',2,0,a,tau,c,1,w,info)
         CALL CHKXER('CUPMTR',INFot,NOUt,LERr,OK)
         nt = nt + 6
!
!        CPTEQR
!
         SRNamt = 'CPTEQR'
         INFot = 1
         CALL CPTEQR('/',0,d,e,z,1,rw,info)
         CALL CHKXER('CPTEQR',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CPTEQR('N',-1,d,e,z,1,rw,info)
         CALL CHKXER('CPTEQR',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL CPTEQR('V',2,d,e,z,1,rw,info)
         CALL CHKXER('CPTEQR',INFot,NOUt,LERr,OK)
         nt = nt + 3
!
!        CSTEIN
!
         SRNamt = 'CSTEIN'
         INFot = 1
         CALL CSTEIN(-1,d,e,0,x,i1,i2,z,1,rw,iw,i3,info)
         CALL CHKXER('CSTEIN',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CSTEIN(0,d,e,-1,x,i1,i2,z,1,rw,iw,i3,info)
         CALL CHKXER('CSTEIN',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CSTEIN(0,d,e,1,x,i1,i2,z,1,rw,iw,i3,info)
         CALL CHKXER('CSTEIN',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL CSTEIN(2,d,e,0,x,i1,i2,z,1,rw,iw,i3,info)
         CALL CHKXER('CSTEIN',INFot,NOUt,LERr,OK)
         nt = nt + 4
!
!        CSTEQR
!
         SRNamt = 'CSTEQR'
         INFot = 1
         CALL CSTEQR('/',0,d,e,z,1,rw,info)
         CALL CHKXER('CSTEQR',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CSTEQR('N',-1,d,e,z,1,rw,info)
         CALL CHKXER('CSTEQR',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL CSTEQR('V',2,d,e,z,1,rw,info)
         CALL CHKXER('CSTEQR',INFot,NOUt,LERr,OK)
         nt = nt + 3
!
!        CSTEDC
!
         SRNamt = 'CSTEDC'
         INFot = 1
         CALL CSTEDC('/',0,d,e,z,1,w,1,rw,1,iw,1,info)
         CALL CHKXER('CSTEDC',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CSTEDC('N',-1,d,e,z,1,w,1,rw,1,iw,1,info)
         CALL CHKXER('CSTEDC',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL CSTEDC('V',2,d,e,z,1,w,4,rw,23,iw,28,info)
         CALL CHKXER('CSTEDC',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL CSTEDC('N',2,d,e,z,1,w,0,rw,1,iw,1,info)
         CALL CHKXER('CSTEDC',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL CSTEDC('V',2,d,e,z,2,w,0,rw,23,iw,28,info)
         CALL CHKXER('CSTEDC',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL CSTEDC('N',2,d,e,z,1,w,1,rw,0,iw,1,info)
         CALL CHKXER('CSTEDC',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL CSTEDC('I',2,d,e,z,2,w,1,rw,1,iw,12,info)
         CALL CHKXER('CSTEDC',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL CSTEDC('V',2,d,e,z,2,w,4,rw,1,iw,28,info)
         CALL CHKXER('CSTEDC',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL CSTEDC('N',2,d,e,z,1,w,1,rw,1,iw,0,info)
         CALL CHKXER('CSTEDC',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL CSTEDC('I',2,d,e,z,2,w,1,rw,23,iw,0,info)
         CALL CHKXER('CSTEDC',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL CSTEDC('V',2,d,e,z,2,w,4,rw,23,iw,0,info)
         CALL CHKXER('CSTEDC',INFot,NOUt,LERr,OK)
         nt = nt + 11
!
!        CHEEVD
!
         SRNamt = 'CHEEVD'
         INFot = 1
         CALL CHEEVD('/','U',0,a,1,x,w,1,rw,1,iw,1,info)
         CALL CHKXER('CHEEVD',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CHEEVD('N','/',0,a,1,x,w,1,rw,1,iw,1,info)
         CALL CHKXER('CHEEVD',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CHEEVD('N','U',-1,a,1,x,w,1,rw,1,iw,1,info)
         CALL CHKXER('CHEEVD',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL CHEEVD('N','U',2,a,1,x,w,3,rw,2,iw,1,info)
         CALL CHKXER('CHEEVD',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL CHEEVD('N','U',1,a,1,x,w,0,rw,1,iw,1,info)
         CALL CHKXER('CHEEVD',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL CHEEVD('N','U',2,a,2,x,w,2,rw,2,iw,1,info)
         CALL CHKXER('CHEEVD',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL CHEEVD('V','U',2,a,2,x,w,3,rw,25,iw,12,info)
         CALL CHKXER('CHEEVD',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL CHEEVD('N','U',1,a,1,x,w,1,rw,0,iw,1,info)
         CALL CHKXER('CHEEVD',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL CHEEVD('N','U',2,a,2,x,w,3,rw,1,iw,1,info)
         CALL CHKXER('CHEEVD',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL CHEEVD('V','U',2,a,2,x,w,8,rw,18,iw,12,info)
         CALL CHKXER('CHEEVD',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL CHEEVD('N','U',1,a,1,x,w,1,rw,1,iw,0,info)
         CALL CHKXER('CHEEVD',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL CHEEVD('V','U',2,a,2,x,w,8,rw,25,iw,11,info)
         CALL CHKXER('CHEEVD',INFot,NOUt,LERr,OK)
         nt = nt + 12
!
!        CHEEVD_2STAGE
!
         SRNamt = 'CHEEVD_2STAGE'
         INFot = 1
         CALL CHEEVD_2STAGE('/','U',0,a,1,x,w,1,rw,1,iw,1,info)
         CALL CHKXER('CHEEVD_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 1
         CALL CHEEVD_2STAGE('V','U',0,a,1,x,w,1,rw,1,iw,1,info)
         CALL CHKXER('CHEEVD_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CHEEVD_2STAGE('N','/',0,a,1,x,w,1,rw,1,iw,1,info)
         CALL CHKXER('CHEEVD_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CHEEVD_2STAGE('N','U',-1,a,1,x,w,1,rw,1,iw,1,info)
         CALL CHKXER('CHEEVD_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL CHEEVD_2STAGE('N','U',2,a,1,x,w,3,rw,2,iw,1,info)
         CALL CHKXER('CHEEVD_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL CHEEVD_2STAGE('N','U',1,a,1,x,w,0,rw,1,iw,1,info)
         CALL CHKXER('CHEEVD_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL CHEEVD_2STAGE('N','U',2,a,2,x,w,2,rw,2,iw,1,info)
         CALL CHKXER('CHEEVD_2STAGE',INFot,NOUt,LERr,OK)
!         INFOT = 8
!         CALL CHEEVD_2STAGE( 'V', 'U', 2, A, 2, X, W, 3,
!     $                            RW, 25, IW, 12, INFO )
!         CALL CHKXER( 'CHEEVD_2STAGE', INFOT, NOUT, LERR, OK )
         INFot = 10
         CALL CHEEVD_2STAGE('N','U',1,a,1,x,w,1,rw,0,iw,1,info)
         CALL CHKXER('CHEEVD_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL CHEEVD_2STAGE('N','U',2,a,2,x,w,25,rw,1,iw,1,info)
         CALL CHKXER('CHEEVD_2STAGE',INFot,NOUt,LERr,OK)
!         INFOT = 10
!         CALL CHEEVD_2STAGE( 'V', 'U', 2, A, 2, X, W, 8,
!     $                            RW, 18, IW, 12, INFO )
!         CALL CHKXER( 'CHEEVD_2STAGE', INFOT, NOUT, LERR, OK )
         INFot = 12
         CALL CHEEVD_2STAGE('N','U',1,a,1,x,w,1,rw,1,iw,0,info)
         CALL CHKXER('CHEEVD_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 12
!         CALL CHEEVD_2STAGE( 'V', 'U', 2, A, 2, X, W, 8,
!     $                            RW, 25, IW, 11, INFO )
!         CALL CHKXER( 'CHEEVD_2STAGE', INFOT, NOUT, LERR, OK )
         nt = nt + 10
!
!        CHEEV
!
         SRNamt = 'CHEEV '
         INFot = 1
         CALL CHEEV('/','U',0,a,1,x,w,1,rw,info)
         CALL CHKXER('CHEEV ',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CHEEV('N','/',0,a,1,x,w,1,rw,info)
         CALL CHKXER('CHEEV ',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CHEEV('N','U',-1,a,1,x,w,1,rw,info)
         CALL CHKXER('CHEEV ',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL CHEEV('N','U',2,a,1,x,w,3,rw,info)
         CALL CHKXER('CHEEV ',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL CHEEV('N','U',2,a,2,x,w,2,rw,info)
         CALL CHKXER('CHEEV ',INFot,NOUt,LERr,OK)
         nt = nt + 5
!
!        CHEEV_2STAGE
!
         SRNamt = 'CHEEV_2STAGE '
         INFot = 1
         CALL CHEEV_2STAGE('/','U',0,a,1,x,w,1,rw,info)
         CALL CHKXER('CHEEV_2STAGE ',INFot,NOUt,LERr,OK)
         INFot = 1
         CALL CHEEV_2STAGE('V','U',0,a,1,x,w,1,rw,info)
         CALL CHKXER('CHEEV_2STAGE ',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CHEEV_2STAGE('N','/',0,a,1,x,w,1,rw,info)
         CALL CHKXER('CHEEV_2STAGE ',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CHEEV_2STAGE('N','U',-1,a,1,x,w,1,rw,info)
         CALL CHKXER('CHEEV_2STAGE ',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL CHEEV_2STAGE('N','U',2,a,1,x,w,3,rw,info)
         CALL CHKXER('CHEEV_2STAGE ',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL CHEEV_2STAGE('N','U',2,a,2,x,w,2,rw,info)
         CALL CHKXER('CHEEV_2STAGE ',INFot,NOUt,LERr,OK)
         nt = nt + 6
!
!        CHEEVX
!
         SRNamt = 'CHEEVX'
         INFot = 1
         CALL CHEEVX('/','A','U',0,a,1,0.0,0.0,0,0,0.0,m,x,z,1,w,1,rw,  &
     &               iw,i3,info)
         CALL CHKXER('CHEEVX',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CHEEVX('V','/','U',0,a,1,0.0,1.0,1,0,0.0,m,x,z,1,w,1,rw,  &
     &               iw,i3,info)
         CALL CHKXER('CHEEVX',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CHEEVX('V','A','/',0,a,1,0.0,0.0,0,0,0.0,m,x,z,1,w,1,rw,  &
     &               iw,i3,info)
         INFot = 4
         CALL CHEEVX('V','A','U',-1,a,1,0.0,0.0,0,0,0.0,m,x,z,1,w,1,rw, &
     &               iw,i3,info)
         CALL CHKXER('CHEEVX',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL CHEEVX('V','A','U',2,a,1,0.0,0.0,0,0,0.0,m,x,z,2,w,3,rw,  &
     &               iw,i3,info)
         CALL CHKXER('CHEEVX',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL CHEEVX('V','V','U',1,a,1,0.0,0.0,0,0,0.0,m,x,z,1,w,1,rw,  &
     &               iw,i3,info)
         CALL CHKXER('CHEEVX',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL CHEEVX('V','I','U',1,a,1,0.0,0.0,0,0,0.0,m,x,z,1,w,1,rw,  &
     &               iw,i3,info)
         CALL CHKXER('CHEEVX',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL CHEEVX('V','I','U',2,a,2,0.0,0.0,2,1,0.0,m,x,z,2,w,3,rw,  &
     &               iw,i3,info)
         CALL CHKXER('CHEEVX',INFot,NOUt,LERr,OK)
         INFot = 15
         CALL CHEEVX('V','A','U',2,a,2,0.0,0.0,0,0,0.0,m,x,z,1,w,3,rw,  &
     &               iw,i3,info)
         CALL CHKXER('CHEEVX',INFot,NOUt,LERr,OK)
         INFot = 17
         CALL CHEEVX('V','A','U',2,a,2,0.0,0.0,0,0,0.0,m,x,z,2,w,2,rw,  &
     &               iw,i1,info)
         CALL CHKXER('CHEEVX',INFot,NOUt,LERr,OK)
         nt = nt + 10
!
!        CHEEVX_2STAGE
!
         SRNamt = 'CHEEVX_2STAGE'
         INFot = 1
         CALL CHEEVX_2STAGE('/','A','U',0,a,1,0.0E0,0.0E0,0,0,0.0E0,m,x,&
     &                      z,1,w,1,rw,iw,i3,info)
         CALL CHKXER('CHEEVX_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 1
         CALL CHEEVX_2STAGE('V','A','U',0,a,1,0.0E0,0.0E0,0,0,0.0E0,m,x,&
     &                      z,1,w,1,rw,iw,i3,info)
         CALL CHKXER('CHEEVX_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CHEEVX_2STAGE('N','/','U',0,a,1,0.0E0,1.0E0,1,0,0.0E0,m,x,&
     &                      z,1,w,1,rw,iw,i3,info)
         CALL CHKXER('CHEEVX_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CHEEVX_2STAGE('N','A','/',0,a,1,0.0E0,0.0E0,0,0,0.0E0,m,x,&
     &                      z,1,w,1,rw,iw,i3,info)
         INFot = 4
         CALL CHEEVX_2STAGE('N','A','U',-1,a,1,0.0E0,0.0E0,0,0,0.0E0,m, &
     &                      x,z,1,w,1,rw,iw,i3,info)
         CALL CHKXER('CHEEVX_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL CHEEVX_2STAGE('N','A','U',2,a,1,0.0E0,0.0E0,0,0,0.0E0,m,x,&
     &                      z,2,w,3,rw,iw,i3,info)
         CALL CHKXER('CHEEVX_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL CHEEVX_2STAGE('N','V','U',1,a,1,0.0E0,0.0E0,0,0,0.0E0,m,x,&
     &                      z,1,w,1,rw,iw,i3,info)
         CALL CHKXER('CHEEVX_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL CHEEVX_2STAGE('N','I','U',1,a,1,0.0E0,0.0E0,0,0,0.0E0,m,x,&
     &                      z,1,w,1,rw,iw,i3,info)
         CALL CHKXER('CHEEVX_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL CHEEVX_2STAGE('N','I','U',2,a,2,0.0E0,0.0E0,2,1,0.0E0,m,x,&
     &                      z,2,w,3,rw,iw,i3,info)
         CALL CHKXER('CHEEVX_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 15
         CALL CHEEVX_2STAGE('N','A','U',2,a,2,0.0E0,0.0E0,0,0,0.0E0,m,x,&
     &                      z,0,w,3,rw,iw,i3,info)
         CALL CHKXER('CHEEVX_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 17
         CALL CHEEVX_2STAGE('N','A','U',2,a,2,0.0E0,0.0E0,0,0,0.0E0,m,x,&
     &                      z,2,w,0,rw,iw,i1,info)
         CALL CHKXER('CHEEVX_2STAGE',INFot,NOUt,LERr,OK)
         nt = nt + 11
!
!        CHEEVR
!
         SRNamt = 'CHEEVR'
         n = 1
         INFot = 1
         CALL CHEEVR('/','A','U',0,a,1,0.0,0.0,1,1,0.0,m,r,z,1,iw,q,2*n,&
     &               rw,24*n,iw(2*n+1),10*n,info)
         CALL CHKXER('CHEEVR',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CHEEVR('V','/','U',0,a,1,0.0,0.0,1,1,0.0,m,r,z,1,iw,q,2*n,&
     &               rw,24*n,iw(2*n+1),10*n,info)
         CALL CHKXER('CHEEVR',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CHEEVR('V','A','/',-1,a,1,0.0,0.0,1,1,0.0,m,r,z,1,iw,q,   &
     &               2*n,rw,24*n,iw(2*n+1),10*n,info)
         CALL CHKXER('CHEEVR',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CHEEVR('V','A','U',-1,a,1,0.0,0.0,1,1,0.0,m,r,z,1,iw,q,   &
     &               2*n,rw,24*n,iw(2*n+1),10*n,info)
         CALL CHKXER('CHEEVR',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL CHEEVR('V','A','U',2,a,1,0.0,0.0,1,1,0.0,m,r,z,1,iw,q,2*n,&
     &               rw,24*n,iw(2*n+1),10*n,info)
         CALL CHKXER('CHEEVR',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL CHEEVR('V','V','U',1,a,1,0.0E0,0.0E0,1,1,0.0,m,r,z,1,iw,q,&
     &               2*n,rw,24*n,iw(2*n+1),10*n,info)
         CALL CHKXER('CHEEVR',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL CHEEVR('V','I','U',1,a,1,0.0E0,0.0E0,0,1,0.0,m,r,z,1,iw,q,&
     &               2*n,rw,24*n,iw(2*n+1),10*n,info)
         CALL CHKXER('CHEEVR',INFot,NOUt,LERr,OK)
         INFot = 10
!
         CALL CHEEVR('V','I','U',2,a,2,0.0E0,0.0E0,2,1,0.0,m,r,z,1,iw,q,&
     &               2*n,rw,24*n,iw(2*n+1),10*n,info)
         CALL CHKXER('CHEEVR',INFot,NOUt,LERr,OK)
         INFot = 15
         CALL CHEEVR('V','I','U',1,a,1,0.0E0,0.0E0,1,1,0.0,m,r,z,0,iw,q,&
     &               2*n,rw,24*n,iw(2*n+1),10*n,info)
         CALL CHKXER('CHEEVR',INFot,NOUt,LERr,OK)
         INFot = 18
         CALL CHEEVR('V','I','U',1,a,1,0.0E0,0.0E0,1,1,0.0,m,r,z,1,iw,q,&
     &               2*n-1,rw,24*n,iw(2*n+1),10*n,info)
         CALL CHKXER('CHEEVR',INFot,NOUt,LERr,OK)
         INFot = 20
         CALL CHEEVR('V','I','U',1,a,1,0.0E0,0.0E0,1,1,0.0,m,r,z,1,iw,q,&
     &               2*n,rw,24*n-1,iw(2*n-1),10*n,info)
         CALL CHKXER('CHEEVR',INFot,NOUt,LERr,OK)
         INFot = 22
         CALL CHEEVR('V','I','U',1,a,1,0.0E0,0.0E0,1,1,0.0,m,r,z,1,iw,q,&
     &               2*n,rw,24*n,iw,10*n-1,info)
         CALL CHKXER('CHEEVR',INFot,NOUt,LERr,OK)
         nt = nt + 12
!
!        CHEEVR_2STAGE
!
         SRNamt = 'CHEEVR_2STAGE'
         n = 1
         INFot = 1
         CALL CHEEVR_2STAGE('/','A','U',0,a,1,0.0D0,0.0D0,1,1,0.0D0,m,r,&
     &                      z,1,iw,q,2*n,rw,24*n,iw(2*n+1),10*n,info)
         CALL CHKXER('CHEEVR_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 1
         CALL CHEEVR_2STAGE('V','A','U',0,a,1,0.0D0,0.0D0,1,1,0.0D0,m,r,&
     &                      z,1,iw,q,2*n,rw,24*n,iw(2*n+1),10*n,info)
         CALL CHKXER('CHEEVR_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CHEEVR_2STAGE('N','/','U',0,a,1,0.0D0,0.0D0,1,1,0.0D0,m,r,&
     &                      z,1,iw,q,2*n,rw,24*n,iw(2*n+1),10*n,info)
         CALL CHKXER('CHEEVR_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CHEEVR_2STAGE('N','A','/',-1,a,1,0.0D0,0.0D0,1,1,0.0D0,m, &
     &                      r,z,1,iw,q,2*n,rw,24*n,iw(2*n+1),10*n,info)
         CALL CHKXER('CHEEVR_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CHEEVR_2STAGE('N','A','U',-1,a,1,0.0D0,0.0D0,1,1,0.0D0,m, &
     &                      r,z,1,iw,q,2*n,rw,24*n,iw(2*n+1),10*n,info)
         CALL CHKXER('CHEEVR_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL CHEEVR_2STAGE('N','A','U',2,a,1,0.0D0,0.0D0,1,1,0.0D0,m,r,&
     &                      z,1,iw,q,2*n,rw,24*n,iw(2*n+1),10*n,info)
         CALL CHKXER('CHEEVR_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL CHEEVR_2STAGE('N','V','U',1,a,1,0.0D0,0.0D0,1,1,0.0D0,m,r,&
     &                      z,1,iw,q,2*n,rw,24*n,iw(2*n+1),10*n,info)
         CALL CHKXER('CHEEVR_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL CHEEVR_2STAGE('N','I','U',1,a,1,0.0D0,0.0D0,0,1,0.0D0,m,r,&
     &                      z,1,iw,q,2*n,rw,24*n,iw(2*n+1),10*n,info)
         CALL CHKXER('CHEEVR_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL CHEEVR_2STAGE('N','I','U',2,a,2,0.0D0,0.0D0,2,1,0.0D0,m,r,&
     &                      z,1,iw,q,2*n,rw,24*n,iw(2*n+1),10*n,info)
         CALL CHKXER('CHEEVR_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 15
         CALL CHEEVR_2STAGE('N','I','U',1,a,1,0.0D0,0.0D0,1,1,0.0D0,m,r,&
     &                      z,0,iw,q,2*n,rw,24*n,iw(2*n+1),10*n,info)
         CALL CHKXER('CHEEVR_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 18
         CALL CHEEVR_2STAGE('N','I','U',1,a,1,0.0D0,0.0D0,1,1,0.0D0,m,r,&
     &                      z,1,iw,q,2*n-1,rw,24*n,iw(2*n+1),10*n,info)
         CALL CHKXER('CHEEVR_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 20
         CALL CHEEVR_2STAGE('N','I','U',1,a,1,0.0D0,0.0D0,1,1,0.0D0,m,r,&
     &                      z,1,iw,q,26*n,rw,24*n-1,iw(2*n-1),10*n,info)
         CALL CHKXER('CHEEVR_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 22
         CALL CHEEVR_2STAGE('N','I','U',1,a,1,0.0D0,0.0D0,1,1,0.0D0,m,r,&
     &                      z,1,iw,q,26*n,rw,24*n,iw,10*n-1,info)
         CALL CHKXER('CHEEVR_2STAGE',INFot,NOUt,LERr,OK)
         nt = nt + 13
!
!        CHPEVD
!
         SRNamt = 'CHPEVD'
         INFot = 1
         CALL CHPEVD('/','U',0,a,x,z,1,w,1,rw,1,iw,1,info)
         CALL CHKXER('CHPEVD',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CHPEVD('N','/',0,a,x,z,1,w,1,rw,1,iw,1,info)
         CALL CHKXER('CHPEVD',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CHPEVD('N','U',-1,a,x,z,1,w,1,rw,1,iw,1,info)
         CALL CHKXER('CHPEVD',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL CHPEVD('V','U',2,a,x,z,1,w,4,rw,25,iw,12,info)
         CALL CHKXER('CHPEVD',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL CHPEVD('N','U',1,a,x,z,1,w,0,rw,1,iw,1,info)
         CALL CHKXER('CHPEVD',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL CHPEVD('N','U',2,a,x,z,2,w,1,rw,2,iw,1,info)
         CALL CHKXER('CHPEVD',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL CHPEVD('V','U',2,a,x,z,2,w,2,rw,25,iw,12,info)
         CALL CHKXER('CHPEVD',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL CHPEVD('N','U',1,a,x,z,1,w,1,rw,0,iw,1,info)
         CALL CHKXER('CHPEVD',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL CHPEVD('N','U',2,a,x,z,2,w,2,rw,1,iw,1,info)
         CALL CHKXER('CHPEVD',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL CHPEVD('V','U',2,a,x,z,2,w,4,rw,18,iw,12,info)
         CALL CHKXER('CHPEVD',INFot,NOUt,LERr,OK)
         INFot = 13
         CALL CHPEVD('N','U',1,a,x,z,1,w,1,rw,1,iw,0,info)
         CALL CHKXER('CHPEVD',INFot,NOUt,LERr,OK)
         INFot = 13
         CALL CHPEVD('N','U',2,a,x,z,2,w,2,rw,2,iw,0,info)
         CALL CHKXER('CHPEVD',INFot,NOUt,LERr,OK)
         INFot = 13
         CALL CHPEVD('V','U',2,a,x,z,2,w,4,rw,25,iw,2,info)
         CALL CHKXER('CHPEVD',INFot,NOUt,LERr,OK)
         nt = nt + 13
!
!        CHPEV
!
         SRNamt = 'CHPEV '
         INFot = 1
         CALL CHPEV('/','U',0,a,x,z,1,w,rw,info)
         CALL CHKXER('CHPEV ',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CHPEV('N','/',0,a,x,z,1,w,rw,info)
         CALL CHKXER('CHPEV ',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CHPEV('N','U',-1,a,x,z,1,w,rw,info)
         CALL CHKXER('CHPEV ',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL CHPEV('V','U',2,a,x,z,1,w,rw,info)
         CALL CHKXER('CHPEV ',INFot,NOUt,LERr,OK)
         nt = nt + 4
!
!        CHPEVX
!
         SRNamt = 'CHPEVX'
         INFot = 1
         CALL CHPEVX('/','A','U',0,a,0.0,0.0,0,0,0.0,m,x,z,1,w,rw,iw,i3,&
     &               info)
         CALL CHKXER('CHPEVX',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CHPEVX('V','/','U',0,a,0.0,1.0,1,0,0.0,m,x,z,1,w,rw,iw,i3,&
     &               info)
         CALL CHKXER('CHPEVX',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CHPEVX('V','A','/',0,a,0.0,0.0,0,0,0.0,m,x,z,1,w,rw,iw,i3,&
     &               info)
         CALL CHKXER('CHPEVX',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CHPEVX('V','A','U',-1,a,0.0,0.0,0,0,0.0,m,x,z,1,w,rw,iw,  &
     &               i3,info)
         CALL CHKXER('CHPEVX',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL CHPEVX('V','V','U',1,a,0.0,0.0,0,0,0.0,m,x,z,1,w,rw,iw,i3,&
     &               info)
         CALL CHKXER('CHPEVX',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL CHPEVX('V','I','U',1,a,0.0,0.0,0,0,0.0,m,x,z,1,w,rw,iw,i3,&
     &               info)
         CALL CHKXER('CHPEVX',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL CHPEVX('V','I','U',2,a,0.0,0.0,2,1,0.0,m,x,z,2,w,rw,iw,i3,&
     &               info)
         CALL CHKXER('CHPEVX',INFot,NOUt,LERr,OK)
         INFot = 14
         CALL CHPEVX('V','A','U',2,a,0.0,0.0,0,0,0.0,m,x,z,1,w,rw,iw,i3,&
     &               info)
         CALL CHKXER('CHPEVX',INFot,NOUt,LERr,OK)
         nt = nt + 8
!
!     Test error exits for the HB path.
!
      ELSEIF ( LSAMEN(2,c2,'HB') ) THEN
!
!        CHBTRD
!
         SRNamt = 'CHBTRD'
         INFot = 1
         CALL CHBTRD('/','U',0,0,a,1,d,e,z,1,w,info)
         CALL CHKXER('CHBTRD',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CHBTRD('N','/',0,0,a,1,d,e,z,1,w,info)
         CALL CHKXER('CHBTRD',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CHBTRD('N','U',-1,0,a,1,d,e,z,1,w,info)
         CALL CHKXER('CHBTRD',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CHBTRD('N','U',0,-1,a,1,d,e,z,1,w,info)
         CALL CHKXER('CHBTRD',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL CHBTRD('N','U',1,1,a,1,d,e,z,1,w,info)
         CALL CHKXER('CHBTRD',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL CHBTRD('V','U',2,0,a,1,d,e,z,1,w,info)
         CALL CHKXER('CHBTRD',INFot,NOUt,LERr,OK)
         nt = nt + 6
!
!        CHETRD_HB2ST
!
         SRNamt = 'CHETRD_HB2ST'
         INFot = 1
         CALL CHETRD_HB2ST('/','N','U',0,0,a,1,d,e,c,1,w,1,info)
         CALL CHKXER('CHETRD_HB2ST',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CHETRD_HB2ST('N','/','U',0,0,a,1,d,e,c,1,w,1,info)
         CALL CHKXER('CHETRD_HB2ST',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CHETRD_HB2ST('N','H','U',0,0,a,1,d,e,c,1,w,1,info)
         CALL CHKXER('CHETRD_HB2ST',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CHETRD_HB2ST('N','N','/',0,0,a,1,d,e,c,1,w,1,info)
         CALL CHKXER('CHETRD_HB2ST',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CHETRD_HB2ST('N','N','U',-1,0,a,1,d,e,c,1,w,1,info)
         CALL CHKXER('CHETRD_HB2ST',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL CHETRD_HB2ST('N','N','U',0,-1,a,1,d,e,c,1,w,1,info)
         CALL CHKXER('CHETRD_HB2ST',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL CHETRD_HB2ST('N','N','U',0,1,a,1,d,e,c,1,w,1,info)
         CALL CHKXER('CHETRD_HB2ST',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL CHETRD_HB2ST('N','N','U',0,0,a,1,d,e,c,0,w,1,info)
         CALL CHKXER('CHETRD_HB2ST',INFot,NOUt,LERr,OK)
         INFot = 13
         CALL CHETRD_HB2ST('N','N','U',0,0,a,1,d,e,c,1,w,0,info)
         CALL CHKXER('CHETRD_HB2ST',INFot,NOUt,LERr,OK)
         nt = nt + 9
!
!        CHBEVD
!
         SRNamt = 'CHBEVD'
         INFot = 1
         CALL CHBEVD('/','U',0,0,a,1,x,z,1,w,1,rw,1,iw,1,info)
         CALL CHKXER('CHBEVD',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CHBEVD('N','/',0,0,a,1,x,z,1,w,1,rw,1,iw,1,info)
         CALL CHKXER('CHBEVD',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CHBEVD('N','U',-1,0,a,1,x,z,1,w,1,rw,1,iw,1,info)
         CALL CHKXER('CHBEVD',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CHBEVD('N','U',0,-1,a,1,x,z,1,w,1,rw,1,iw,1,info)
         CALL CHKXER('CHBEVD',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL CHBEVD('N','U',2,1,a,1,x,z,1,w,2,rw,2,iw,1,info)
         CALL CHKXER('CHBEVD',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL CHBEVD('V','U',2,1,a,2,x,z,1,w,8,rw,25,iw,12,info)
         CALL CHKXER('CHBEVD',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL CHBEVD('N','U',1,0,a,1,x,z,1,w,0,rw,1,iw,1,info)
         CALL CHKXER('CHBEVD',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL CHBEVD('N','U',2,1,a,2,x,z,2,w,1,rw,2,iw,1,info)
         CALL CHKXER('CHBEVD',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL CHBEVD('V','U',2,1,a,2,x,z,2,w,2,rw,25,iw,12,info)
         CALL CHKXER('CHBEVD',INFot,NOUt,LERr,OK)
         INFot = 13
         CALL CHBEVD('N','U',1,0,a,1,x,z,1,w,1,rw,0,iw,1,info)
         CALL CHKXER('CHBEVD',INFot,NOUt,LERr,OK)
         INFot = 13
         CALL CHBEVD('N','U',2,1,a,2,x,z,2,w,2,rw,1,iw,1,info)
         CALL CHKXER('CHBEVD',INFot,NOUt,LERr,OK)
         INFot = 13
         CALL CHBEVD('V','U',2,1,a,2,x,z,2,w,8,rw,2,iw,12,info)
         CALL CHKXER('CHBEVD',INFot,NOUt,LERr,OK)
         INFot = 15
         CALL CHBEVD('N','U',1,0,a,1,x,z,1,w,1,rw,1,iw,0,info)
         CALL CHKXER('CHBEVD',INFot,NOUt,LERr,OK)
         INFot = 15
         CALL CHBEVD('N','U',2,1,a,2,x,z,2,w,2,rw,2,iw,0,info)
         CALL CHKXER('CHBEVD',INFot,NOUt,LERr,OK)
         INFot = 15
         CALL CHBEVD('V','U',2,1,a,2,x,z,2,w,8,rw,25,iw,2,info)
         CALL CHKXER('CHBEVD',INFot,NOUt,LERr,OK)
         nt = nt + 15
!
!        CHBEVD_2STAGE
!
         SRNamt = 'CHBEVD_2STAGE'
         INFot = 1
         CALL CHBEVD_2STAGE('/','U',0,0,a,1,x,z,1,w,1,rw,1,iw,1,info)
         CALL CHKXER('CHBEVD_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 1
         CALL CHBEVD_2STAGE('V','U',0,0,a,1,x,z,1,w,1,rw,1,iw,1,info)
         CALL CHKXER('CHBEVD_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CHBEVD_2STAGE('N','/',0,0,a,1,x,z,1,w,1,rw,1,iw,1,info)
         CALL CHKXER('CHBEVD_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CHBEVD_2STAGE('N','U',-1,0,a,1,x,z,1,w,1,rw,1,iw,1,info)
         CALL CHKXER('CHBEVD_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CHBEVD_2STAGE('N','U',0,-1,a,1,x,z,1,w,1,rw,1,iw,1,info)
         CALL CHKXER('CHBEVD_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL CHBEVD_2STAGE('N','U',2,1,a,1,x,z,1,w,2,rw,2,iw,1,info)
         CALL CHKXER('CHBEVD_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL CHBEVD_2STAGE('N','U',2,1,a,2,x,z,0,w,8,rw,25,iw,12,info)
         CALL CHKXER('CHBEVD_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL CHBEVD_2STAGE('N','U',1,0,a,1,x,z,1,w,0,rw,1,iw,1,info)
         CALL CHKXER('CHBEVD_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL CHBEVD_2STAGE('N','U',2,1,a,2,x,z,2,w,1,rw,2,iw,1,info)
         CALL CHKXER('CHBEVD_2STAGE',INFot,NOUt,LERr,OK)
!         INFOT = 11
!         CALL CHBEVD_2STAGE( 'V', 'U', 2, 1, A, 2, X, Z, 2,
!     $                         W, 2, RW, 25, IW, 12, INFO )
!         CALL CHKXER( 'CHBEVD_2STAGE', INFOT, NOUT, LERR, OK )
         INFot = 13
         CALL CHBEVD_2STAGE('N','U',1,0,a,1,x,z,1,w,1,rw,0,iw,1,info)
         CALL CHKXER('CHBEVD_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 13
         CALL CHBEVD_2STAGE('N','U',2,1,a,2,x,z,2,w,25,rw,1,iw,1,info)
         CALL CHKXER('CHBEVD_2STAGE',INFot,NOUt,LERr,OK)
!         INFOT = 13
!         CALL CHBEVD_2STAGE( 'V', 'U', 2, 1, A, 2, X, Z, 2,
!     $                          W, 25, RW, 2, IW, 12, INFO )
!         CALL CHKXER( 'CHBEVD_2STAGE', INFOT, NOUT, LERR, OK )
         INFot = 15
         CALL CHBEVD_2STAGE('N','U',1,0,a,1,x,z,1,w,1,rw,1,iw,0,info)
         CALL CHKXER('CHBEVD_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 15
         CALL CHBEVD_2STAGE('N','U',2,1,a,2,x,z,2,w,25,rw,2,iw,0,info)
         CALL CHKXER('CHBEVD_2STAGE',INFot,NOUt,LERr,OK)
!         INFOT = 15
!         CALL CHBEVD_2STAGE( 'V', 'U', 2, 1, A, 2, X, Z, 2,
!     $                          W, 25, RW, 25, IW, 2, INFO )
!         CALL CHKXER( 'CHBEVD_2STAGE', INFOT, NOUT, LERR, OK )
         nt = nt + 13
!
!        CHBEV
!
         SRNamt = 'CHBEV '
         INFot = 1
         CALL CHBEV('/','U',0,0,a,1,x,z,1,w,rw,info)
         CALL CHKXER('CHBEV ',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CHBEV('N','/',0,0,a,1,x,z,1,w,rw,info)
         CALL CHKXER('CHBEV ',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CHBEV('N','U',-1,0,a,1,x,z,1,w,rw,info)
         CALL CHKXER('CHBEV ',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CHBEV('N','U',0,-1,a,1,x,z,1,w,rw,info)
         CALL CHKXER('CHBEV ',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL CHBEV('N','U',2,1,a,1,x,z,1,w,rw,info)
         CALL CHKXER('CHBEV ',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL CHBEV('V','U',2,0,a,1,x,z,1,w,rw,info)
         CALL CHKXER('CHBEV ',INFot,NOUt,LERr,OK)
         nt = nt + 6
!
!        CHBEV_2STAGE
!
         SRNamt = 'CHBEV_2STAGE '
         INFot = 1
         CALL CHBEV_2STAGE('/','U',0,0,a,1,x,z,1,w,0,rw,info)
         CALL CHKXER('CHBEV_2STAGE ',INFot,NOUt,LERr,OK)
         INFot = 1
         CALL CHBEV_2STAGE('V','U',0,0,a,1,x,z,1,w,0,rw,info)
         CALL CHKXER('CHBEV_2STAGE ',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CHBEV_2STAGE('N','/',0,0,a,1,x,z,1,w,0,rw,info)
         CALL CHKXER('CHBEV_2STAGE ',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CHBEV_2STAGE('N','U',-1,0,a,1,x,z,1,w,0,rw,info)
         CALL CHKXER('CHBEV_2STAGE ',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CHBEV_2STAGE('N','U',0,-1,a,1,x,z,1,w,0,rw,info)
         CALL CHKXER('CHBEV_2STAGE ',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL CHBEV_2STAGE('N','U',2,1,a,1,x,z,1,w,0,rw,info)
         CALL CHKXER('CHBEV_2STAGE ',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL CHBEV_2STAGE('N','U',2,0,a,1,x,z,0,w,0,rw,info)
         CALL CHKXER('CHBEV_2STAGE ',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL CHBEV_2STAGE('N','U',2,0,a,1,x,z,1,w,0,rw,info)
         CALL CHKXER('CHBEV_2STAGE ',INFot,NOUt,LERr,OK)
         nt = nt + 8
!
!        CHBEVX
!
         SRNamt = 'CHBEVX'
         INFot = 1
         CALL CHBEVX('/','A','U',0,0,a,1,q,1,0.0,0.0,0,0,0.0,m,x,z,1,w, &
     &               rw,iw,i3,info)
         CALL CHKXER('CHBEVX',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CHBEVX('V','/','U',0,0,a,1,q,1,0.0,1.0,1,0,0.0,m,x,z,1,w, &
     &               rw,iw,i3,info)
         CALL CHKXER('CHBEVX',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CHBEVX('V','A','/',0,0,a,1,q,1,0.0,0.0,0,0,0.0,m,x,z,1,w, &
     &               rw,iw,i3,info)
         INFot = 4
         CALL CHBEVX('V','A','U',-1,0,a,1,q,1,0.0,0.0,0,0,0.0,m,x,z,1,w,&
     &               rw,iw,i3,info)
         CALL CHKXER('CHBEVX',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL CHBEVX('V','A','U',0,-1,a,1,q,1,0.0,0.0,0,0,0.0,m,x,z,1,w,&
     &               rw,iw,i3,info)
         CALL CHKXER('CHBEVX',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL CHBEVX('V','A','U',2,1,a,1,q,2,0.0,0.0,0,0,0.0,m,x,z,2,w, &
     &               rw,iw,i3,info)
         CALL CHKXER('CHBEVX',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL CHBEVX('V','A','U',2,0,a,1,q,1,0.0,0.0,0,0,0.0,m,x,z,2,w, &
     &               rw,iw,i3,info)
         CALL CHKXER('CHBEVX',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL CHBEVX('V','V','U',1,0,a,1,q,1,0.0,0.0,0,0,0.0,m,x,z,1,w, &
     &               rw,iw,i3,info)
         CALL CHKXER('CHBEVX',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL CHBEVX('V','I','U',1,0,a,1,q,1,0.0,0.0,0,0,0.0,m,x,z,1,w, &
     &               rw,iw,i3,info)
         CALL CHKXER('CHBEVX',INFot,NOUt,LERr,OK)
         INFot = 13
         CALL CHBEVX('V','I','U',1,0,a,1,q,1,0.0,0.0,1,2,0.0,m,x,z,1,w, &
     &               rw,iw,i3,info)
         CALL CHKXER('CHBEVX',INFot,NOUt,LERr,OK)
         INFot = 18
         CALL CHBEVX('V','A','U',2,0,a,1,q,2,0.0,0.0,0,0,0.0,m,x,z,1,w, &
     &               rw,iw,i3,info)
         CALL CHKXER('CHBEVX',INFot,NOUt,LERr,OK)
         nt = nt + 11
!
!        CHBEVX_2STAGE
!
         SRNamt = 'CHBEVX_2STAGE'
         INFot = 1
         CALL CHBEVX_2STAGE('/','A','U',0,0,a,1,q,1,0.0D0,0.0D0,0,0,    &
     &                      0.0D0,m,x,z,1,w,0,rw,iw,i3,info)
         INFot = 1
         CALL CHBEVX_2STAGE('V','A','U',0,0,a,1,q,1,0.0D0,0.0D0,0,0,    &
     &                      0.0D0,m,x,z,1,w,0,rw,iw,i3,info)
         CALL CHKXER('CHBEVX_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CHBEVX_2STAGE('N','/','U',0,0,a,1,q,1,0.0D0,1.0D0,1,0,    &
     &                      0.0D0,m,x,z,1,w,0,rw,iw,i3,info)
         CALL CHKXER('CHBEVX_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CHBEVX_2STAGE('N','A','/',0,0,a,1,q,1,0.0D0,0.0D0,0,0,    &
     &                      0.0D0,m,x,z,1,w,0,rw,iw,i3,info)
         INFot = 4
         CALL CHBEVX_2STAGE('N','A','U',-1,0,a,1,q,1,0.0D0,0.0D0,0,0,   &
     &                      0.0D0,m,x,z,1,w,0,rw,iw,i3,info)
         CALL CHKXER('CHBEVX_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL CHBEVX_2STAGE('N','A','U',0,-1,a,1,q,1,0.0D0,0.0D0,0,0,   &
     &                      0.0D0,m,x,z,1,w,0,rw,iw,i3,info)
         CALL CHKXER('CHBEVX_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL CHBEVX_2STAGE('N','A','U',2,1,a,1,q,2,0.0D0,0.0D0,0,0,    &
     &                      0.0D0,m,x,z,2,w,0,rw,iw,i3,info)
         CALL CHKXER('CHBEVX_2STAGE',INFot,NOUt,LERr,OK)
!         INFOT = 9
!         CALL CHBEVX_2STAGE( 'V', 'A', 'U', 2, 0, A, 1, Q, 1,
!     $                       0.0D0, 0.0D0, 0, 0, 0.0D0,
!     $                       M, X, Z, 2, W, 0, RW, IW, I3, INFO )
!         CALL CHKXER( 'CHBEVX_2STAGE', INFOT, NOUT, LERR, OK )
         INFot = 11
         CALL CHBEVX_2STAGE('N','V','U',1,0,a,1,q,1,0.0D0,0.0D0,0,0,    &
     &                      0.0D0,m,x,z,1,w,0,rw,iw,i3,info)
         CALL CHKXER('CHBEVX_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL CHBEVX_2STAGE('N','I','U',1,0,a,1,q,1,0.0D0,0.0D0,0,0,    &
     &                      0.0D0,m,x,z,1,w,0,rw,iw,i3,info)
         CALL CHKXER('CHBEVX_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 13
         CALL CHBEVX_2STAGE('N','I','U',1,0,a,1,q,1,0.0D0,0.0D0,1,2,    &
     &                      0.0D0,m,x,z,1,w,0,rw,iw,i3,info)
         CALL CHKXER('CHBEVX_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 18
         CALL CHBEVX_2STAGE('N','A','U',2,0,a,1,q,2,0.0D0,0.0D0,0,0,    &
     &                      0.0D0,m,x,z,0,w,0,rw,iw,i3,info)
         CALL CHKXER('CHBEVX_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 20
         CALL CHBEVX_2STAGE('N','A','U',2,0,a,1,q,2,0.0D0,0.0D0,0,0,    &
     &                      0.0D0,m,x,z,1,w,0,rw,iw,i3,info)
         CALL CHKXER('CHBEVX_2STAGE',INFot,NOUt,LERr,OK)
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
!     End of CERRST
!
      END SUBROUTINE CERRST
