!*==derrvx.f90  processed by SPAG 7.51RB at 17:34 on  4 Mar 2022
!> \brief \b DERRVXX
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE derrvx( PATH, NUNIT )
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
!> derrvx tests the error exits for the DOUBLE PRECISION driver routines
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
!> \date December 2016
!
!> \ingroup double_lin
!
!  =====================================================================
      SUBROUTINE derrvx(Path,Nunit)
      IMPLICIT NONE
!*--DERRVX59
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
      PARAMETER (NMAX=4)
      REAL ONE
      PARAMETER (ONE=1.0D+0)
!     ..
!     .. Local Scalars ..
      CHARACTER eq
      CHARACTER*2 c2
      INTEGER i , info , j , n_err_bnds , nparams
      DOUBLE PRECISION rcond , rpvgrw , berr
!     ..
!     .. Local Arrays ..
      INTEGER ip(NMAX) , iw(NMAX)
      DOUBLE PRECISION a(NMAX,NMAX) , af(NMAX,NMAX) , b(NMAX) , c(NMAX) &
     &                 , e(NMAX) , r(NMAX) , r1(NMAX) , r2(NMAX) ,      &
     &                 w(2*NMAX) , x(NMAX) , err_bnds_n(NMAX,3) ,       &
     &                 err_bnds_c(NMAX,3) , params(1)
!     ..
!     .. External Functions ..
      LOGICAL LSAMEN
EXTERNAL    LSAMEN
!     ..
!     .. External Subroutines ..
EXTERNAL    DGBSV,   DGBSVX,   DGESV,   DGESVX,   DGTSV
EXTERNAL    DGTSVX,   DPBSV,   DPBSVX,   DPOSV,   DPOSVX,   DPPSV
EXTERNAL    DPPSVX,   DPTSV,   DPTSVX,   DSPSV,   DSPSVX,   DSYSV
EXTERNAL    DSYSV_RK,   DSYSV_ROOK,   DSYSVX,   DGESVXX,   DSYSVXX
EXTERNAL    DPOSVXX,   DGBSVXX
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
            af(i,j) = 1.D0/DBLE(i+j)
         ENDDO
         b(j) = 0.D+0
         e(j) = 0.D+0
         r1(j) = 0.D+0
         r2(j) = 0.D+0
         w(j) = 0.D+0
         x(j) = 0.D+0
         c(j) = 0.D+0
         r(j) = 0.D+0
         ip(j) = j
      ENDDO
      eq = ' '
      OK = .TRUE.
!
      IF ( LSAMEN(2,c2,'GE') ) THEN
!
!        DGESV
!
         SRNamt = 'DGESV '
         INFot = 1
         CALL DGESV(-1,0,a,1,ip,b,1,info)
         CALL CHKXER('DGESV ',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DGESV(0,-1,a,1,ip,b,1,info)
         CALL CHKXER('DGESV ',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL DGESV(2,1,a,1,ip,b,2,info)
         CALL CHKXER('DGESV ',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL DGESV(2,1,a,2,ip,b,1,info)
         CALL CHKXER('DGESV ',INFot,NOUt,LERr,OK)
!
!        DGESVX
!
         SRNamt = 'DGESVX'
         INFot = 1
         CALL DGESVX('/','N',0,0,a,1,af,1,ip,eq,r,c,b,1,x,1,rcond,r1,r2,&
     &               w,iw,info)
         CALL CHKXER('DGESVX',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DGESVX('N','/',0,0,a,1,af,1,ip,eq,r,c,b,1,x,1,rcond,r1,r2,&
     &               w,iw,info)
         CALL CHKXER('DGESVX',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DGESVX('N','N',-1,0,a,1,af,1,ip,eq,r,c,b,1,x,1,rcond,r1,  &
     &               r2,w,iw,info)
         CALL CHKXER('DGESVX',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL DGESVX('N','N',0,-1,a,1,af,1,ip,eq,r,c,b,1,x,1,rcond,r1,  &
     &               r2,w,iw,info)
         CALL CHKXER('DGESVX',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL DGESVX('N','N',2,1,a,1,af,2,ip,eq,r,c,b,2,x,2,rcond,r1,r2,&
     &               w,iw,info)
         CALL CHKXER('DGESVX',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL DGESVX('N','N',2,1,a,2,af,1,ip,eq,r,c,b,2,x,2,rcond,r1,r2,&
     &               w,iw,info)
         CALL CHKXER('DGESVX',INFot,NOUt,LERr,OK)
         INFot = 10
         eq = '/'
         CALL DGESVX('F','N',0,0,a,1,af,1,ip,eq,r,c,b,1,x,1,rcond,r1,r2,&
     &               w,iw,info)
         CALL CHKXER('DGESVX',INFot,NOUt,LERr,OK)
         INFot = 11
         eq = 'R'
         CALL DGESVX('F','N',1,0,a,1,af,1,ip,eq,r,c,b,1,x,1,rcond,r1,r2,&
     &               w,iw,info)
         CALL CHKXER('DGESVX',INFot,NOUt,LERr,OK)
         INFot = 12
         eq = 'C'
         CALL DGESVX('F','N',1,0,a,1,af,1,ip,eq,r,c,b,1,x,1,rcond,r1,r2,&
     &               w,iw,info)
         CALL CHKXER('DGESVX',INFot,NOUt,LERr,OK)
         INFot = 14
         CALL DGESVX('N','N',2,1,a,2,af,2,ip,eq,r,c,b,1,x,2,rcond,r1,r2,&
     &               w,iw,info)
         CALL CHKXER('DGESVX',INFot,NOUt,LERr,OK)
         INFot = 16
         CALL DGESVX('N','N',2,1,a,2,af,2,ip,eq,r,c,b,2,x,1,rcond,r1,r2,&
     &               w,iw,info)
         CALL CHKXER('DGESVX',INFot,NOUt,LERr,OK)
!
!        DGESVXX
!
         n_err_bnds = 3
         nparams = 1
         SRNamt = 'DGESVXX'
         INFot = 1
         CALL DGESVXX('/','N',0,0,a,1,af,1,ip,eq,r,c,b,1,x,1,rcond,     &
     &                rpvgrw,berr,n_err_bnds,err_bnds_n,err_bnds_c,     &
     &                nparams,params,w,iw,info)
         CALL CHKXER('DGESVXX',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DGESVXX('N','/',0,0,a,1,af,1,ip,eq,r,c,b,1,x,1,rcond,     &
     &                rpvgrw,berr,n_err_bnds,err_bnds_n,err_bnds_c,     &
     &                nparams,params,w,iw,info)
         CALL CHKXER('DGESVXX',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DGESVXX('N','N',-1,0,a,1,af,1,ip,eq,r,c,b,1,x,1,rcond,    &
     &                rpvgrw,berr,n_err_bnds,err_bnds_n,err_bnds_c,     &
     &                nparams,params,w,iw,info)
         CALL CHKXER('DGESVXX',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL DGESVXX('N','N',0,-1,a,1,af,1,ip,eq,r,c,b,1,x,1,rcond,    &
     &                rpvgrw,berr,n_err_bnds,err_bnds_n,err_bnds_c,     &
     &                nparams,params,w,iw,info)
         CALL CHKXER('DGESVXX',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL DGESVXX('N','N',2,1,a,1,af,2,ip,eq,r,c,b,2,x,2,rcond,     &
     &                rpvgrw,berr,n_err_bnds,err_bnds_n,err_bnds_c,     &
     &                nparams,params,w,iw,info)
         CALL CHKXER('DGESVXX',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL DGESVXX('N','N',2,1,a,2,af,1,ip,eq,r,c,b,2,x,2,rcond,     &
     &                rpvgrw,berr,n_err_bnds,err_bnds_n,err_bnds_c,     &
     &                nparams,params,w,iw,info)
         CALL CHKXER('DGESVXX',INFot,NOUt,LERr,OK)
         INFot = 10
         eq = '/'
         CALL DGESVXX('F','N',0,0,a,1,af,1,ip,eq,r,c,b,1,x,1,rcond,     &
     &                rpvgrw,berr,n_err_bnds,err_bnds_n,err_bnds_c,     &
     &                nparams,params,w,iw,info)
         CALL CHKXER('DGESVXX',INFot,NOUt,LERr,OK)
         INFot = 11
         eq = 'R'
         CALL DGESVXX('F','N',1,0,a,1,af,1,ip,eq,r,c,b,1,x,1,rcond,     &
     &                rpvgrw,berr,n_err_bnds,err_bnds_n,err_bnds_c,     &
     &                nparams,params,w,iw,info)
         CALL CHKXER('DGESVXX',INFot,NOUt,LERr,OK)
         INFot = 12
         eq = 'C'
         CALL DGESVXX('F','N',1,0,a,1,af,1,ip,eq,r,c,b,1,x,1,rcond,     &
     &                rpvgrw,berr,n_err_bnds,err_bnds_n,err_bnds_c,     &
     &                nparams,params,w,iw,info)
         CALL CHKXER('DGESVXX',INFot,NOUt,LERr,OK)
         INFot = 14
         CALL DGESVXX('N','N',2,1,a,2,af,2,ip,eq,r,c,b,1,x,2,rcond,     &
     &                rpvgrw,berr,n_err_bnds,err_bnds_n,err_bnds_c,     &
     &                nparams,params,w,iw,info)
         CALL CHKXER('DGESVXX',INFot,NOUt,LERr,OK)
         INFot = 16
         CALL DGESVXX('N','N',2,1,a,2,af,2,ip,eq,r,c,b,2,x,1,rcond,     &
     &                rpvgrw,berr,n_err_bnds,err_bnds_n,err_bnds_c,     &
     &                nparams,params,w,iw,info)
         CALL CHKXER('DGESVXX',INFot,NOUt,LERr,OK)
!
      ELSEIF ( LSAMEN(2,c2,'GB') ) THEN
!
!        DGBSV
!
         SRNamt = 'DGBSV '
         INFot = 1
         CALL DGBSV(-1,0,0,0,a,1,ip,b,1,info)
         CALL CHKXER('DGBSV ',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DGBSV(1,-1,0,0,a,1,ip,b,1,info)
         CALL CHKXER('DGBSV ',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DGBSV(1,0,-1,0,a,1,ip,b,1,info)
         CALL CHKXER('DGBSV ',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL DGBSV(0,0,0,-1,a,1,ip,b,1,info)
         CALL CHKXER('DGBSV ',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL DGBSV(1,1,1,0,a,3,ip,b,1,info)
         CALL CHKXER('DGBSV ',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL DGBSV(2,0,0,0,a,1,ip,b,1,info)
         CALL CHKXER('DGBSV ',INFot,NOUt,LERr,OK)
!
!        DGBSVX
!
         SRNamt = 'DGBSVX'
         INFot = 1
         CALL DGBSVX('/','N',0,0,0,0,a,1,af,1,ip,eq,r,c,b,1,x,1,rcond,  &
     &               r1,r2,w,iw,info)
         CALL CHKXER('DGBSVX',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DGBSVX('N','/',0,0,0,0,a,1,af,1,ip,eq,r,c,b,1,x,1,rcond,  &
     &               r1,r2,w,iw,info)
         CALL CHKXER('DGBSVX',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DGBSVX('N','N',-1,0,0,0,a,1,af,1,ip,eq,r,c,b,1,x,1,rcond, &
     &               r1,r2,w,iw,info)
         CALL CHKXER('DGBSVX',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL DGBSVX('N','N',1,-1,0,0,a,1,af,1,ip,eq,r,c,b,1,x,1,rcond, &
     &               r1,r2,w,iw,info)
         CALL CHKXER('DGBSVX',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL DGBSVX('N','N',1,0,-1,0,a,1,af,1,ip,eq,r,c,b,1,x,1,rcond, &
     &               r1,r2,w,iw,info)
         CALL CHKXER('DGBSVX',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL DGBSVX('N','N',0,0,0,-1,a,1,af,1,ip,eq,r,c,b,1,x,1,rcond, &
     &               r1,r2,w,iw,info)
         CALL CHKXER('DGBSVX',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL DGBSVX('N','N',1,1,1,0,a,2,af,4,ip,eq,r,c,b,1,x,1,rcond,  &
     &               r1,r2,w,iw,info)
         CALL CHKXER('DGBSVX',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL DGBSVX('N','N',1,1,1,0,a,3,af,3,ip,eq,r,c,b,1,x,1,rcond,  &
     &               r1,r2,w,iw,info)
         CALL CHKXER('DGBSVX',INFot,NOUt,LERr,OK)
         INFot = 12
         eq = '/'
         CALL DGBSVX('F','N',0,0,0,0,a,1,af,1,ip,eq,r,c,b,1,x,1,rcond,  &
     &               r1,r2,w,iw,info)
         CALL CHKXER('DGBSVX',INFot,NOUt,LERr,OK)
         INFot = 13
         eq = 'R'
         CALL DGBSVX('F','N',1,0,0,0,a,1,af,1,ip,eq,r,c,b,1,x,1,rcond,  &
     &               r1,r2,w,iw,info)
         CALL CHKXER('DGBSVX',INFot,NOUt,LERr,OK)
         INFot = 14
         eq = 'C'
         CALL DGBSVX('F','N',1,0,0,0,a,1,af,1,ip,eq,r,c,b,1,x,1,rcond,  &
     &               r1,r2,w,iw,info)
         CALL CHKXER('DGBSVX',INFot,NOUt,LERr,OK)
         INFot = 16
         CALL DGBSVX('N','N',2,0,0,0,a,1,af,1,ip,eq,r,c,b,1,x,2,rcond,  &
     &               r1,r2,w,iw,info)
         CALL CHKXER('DGBSVX',INFot,NOUt,LERr,OK)
         INFot = 18
         CALL DGBSVX('N','N',2,0,0,0,a,1,af,1,ip,eq,r,c,b,2,x,1,rcond,  &
     &               r1,r2,w,iw,info)
         CALL CHKXER('DGBSVX',INFot,NOUt,LERr,OK)
!
!        DGBSVXX
!
         n_err_bnds = 3
         nparams = 1
         SRNamt = 'DGBSVXX'
         INFot = 1
         CALL DGBSVXX('/','N',0,0,0,0,a,1,af,1,ip,eq,r,c,b,1,x,1,rcond, &
     &                rpvgrw,berr,n_err_bnds,err_bnds_n,err_bnds_c,     &
     &                nparams,params,w,iw,info)
         CALL CHKXER('DGBSVXX',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DGBSVXX('N','/',0,1,1,0,a,1,af,1,ip,eq,r,c,b,1,x,1,rcond, &
     &                rpvgrw,berr,n_err_bnds,err_bnds_n,err_bnds_c,     &
     &                nparams,params,w,iw,info)
         CALL CHKXER('DGBSVXX',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DGBSVXX('N','N',-1,1,1,0,a,1,af,1,ip,eq,r,c,b,1,x,1,rcond,&
     &                rpvgrw,berr,n_err_bnds,err_bnds_n,err_bnds_c,     &
     &                nparams,params,w,iw,info)
         CALL CHKXER('DGBSVXX',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL DGBSVXX('N','N',2,-1,1,0,a,1,af,1,ip,eq,r,c,b,1,x,1,rcond,&
     &                rpvgrw,berr,n_err_bnds,err_bnds_n,err_bnds_c,     &
     &                nparams,params,w,iw,info)
         CALL CHKXER('DGBSVXX',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL DGBSVXX('N','N',2,1,-1,0,a,1,af,1,ip,eq,r,c,b,1,x,1,rcond,&
     &                rpvgrw,berr,n_err_bnds,err_bnds_n,err_bnds_c,     &
     &                nparams,params,w,iw,info)
         CALL CHKXER('DGBSVXX',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL DGBSVXX('N','N',0,1,1,-1,a,1,af,1,ip,eq,r,c,b,1,x,1,rcond,&
     &                rpvgrw,berr,n_err_bnds,err_bnds_n,err_bnds_c,     &
     &                nparams,params,w,iw,info)
         CALL CHKXER('DGBSVXX',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL DGBSVXX('N','N',2,1,1,1,a,2,af,2,ip,eq,r,c,b,2,x,2,rcond, &
     &                rpvgrw,berr,n_err_bnds,err_bnds_n,err_bnds_c,     &
     &                nparams,params,w,iw,info)
         CALL CHKXER('DGBSVXX',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL DGBSVXX('N','N',2,1,1,1,a,3,af,3,ip,eq,r,c,b,2,x,2,rcond, &
     &                rpvgrw,berr,n_err_bnds,err_bnds_n,err_bnds_c,     &
     &                nparams,params,w,iw,info)
         CALL CHKXER('DGBSVXX',INFot,NOUt,LERr,OK)
         INFot = 12
         eq = '/'
         CALL DGBSVXX('F','N',0,1,1,0,a,3,af,4,ip,eq,r,c,b,1,x,1,rcond, &
     &                rpvgrw,berr,n_err_bnds,err_bnds_n,err_bnds_c,     &
     &                nparams,params,w,iw,info)
         CALL CHKXER('DGBSVXX',INFot,NOUt,LERr,OK)
         INFot = 13
         eq = 'R'
         CALL DGBSVXX('F','N',1,1,1,0,a,3,af,4,ip,eq,r,c,b,1,x,1,rcond, &
     &                rpvgrw,berr,n_err_bnds,err_bnds_n,err_bnds_c,     &
     &                nparams,params,w,iw,info)
         CALL CHKXER('DGBSVXX',INFot,NOUt,LERr,OK)
         INFot = 14
         eq = 'C'
         CALL DGBSVXX('F','N',1,1,1,0,a,3,af,4,ip,eq,r,c,b,1,x,1,rcond, &
     &                rpvgrw,berr,n_err_bnds,err_bnds_n,err_bnds_c,     &
     &                nparams,params,w,iw,info)
         CALL CHKXER('DGBSVXX',INFot,NOUt,LERr,OK)
         INFot = 15
         CALL DGBSVXX('N','N',2,1,1,1,a,3,af,4,ip,eq,r,c,b,1,x,2,rcond, &
     &                rpvgrw,berr,n_err_bnds,err_bnds_n,err_bnds_c,     &
     &                nparams,params,w,iw,info)
         CALL CHKXER('DGBSVXX',INFot,NOUt,LERr,OK)
         INFot = 16
         CALL DGBSVXX('N','N',2,1,1,1,a,3,af,4,ip,eq,r,c,b,2,x,1,rcond, &
     &                rpvgrw,berr,n_err_bnds,err_bnds_n,err_bnds_c,     &
     &                nparams,params,w,iw,info)
         CALL CHKXER('DGBSVXX',INFot,NOUt,LERr,OK)
!
      ELSEIF ( LSAMEN(2,c2,'GT') ) THEN
!
!        DGTSV
!
         SRNamt = 'DGTSV '
         INFot = 1
         CALL DGTSV(-1,0,a(1,1),a(1,2),a(1,3),b,1,info)
         CALL CHKXER('DGTSV ',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DGTSV(0,-1,a(1,1),a(1,2),a(1,3),b,1,info)
         CALL CHKXER('DGTSV ',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL DGTSV(2,0,a(1,1),a(1,2),a(1,3),b,1,info)
         CALL CHKXER('DGTSV ',INFot,NOUt,LERr,OK)
!
!        DGTSVX
!
         SRNamt = 'DGTSVX'
         INFot = 1
         CALL DGTSVX('/','N',0,0,a(1,1),a(1,2),a(1,3),af(1,1),af(1,2),  &
     &               af(1,3),af(1,4),ip,b,1,x,1,rcond,r1,r2,w,iw,info)
         CALL CHKXER('DGTSVX',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DGTSVX('N','/',0,0,a(1,1),a(1,2),a(1,3),af(1,1),af(1,2),  &
     &               af(1,3),af(1,4),ip,b,1,x,1,rcond,r1,r2,w,iw,info)
         CALL CHKXER('DGTSVX',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DGTSVX('N','N',-1,0,a(1,1),a(1,2),a(1,3),af(1,1),af(1,2), &
     &               af(1,3),af(1,4),ip,b,1,x,1,rcond,r1,r2,w,iw,info)
         CALL CHKXER('DGTSVX',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL DGTSVX('N','N',0,-1,a(1,1),a(1,2),a(1,3),af(1,1),af(1,2), &
     &               af(1,3),af(1,4),ip,b,1,x,1,rcond,r1,r2,w,iw,info)
         CALL CHKXER('DGTSVX',INFot,NOUt,LERr,OK)
         INFot = 14
         CALL DGTSVX('N','N',2,0,a(1,1),a(1,2),a(1,3),af(1,1),af(1,2),  &
     &               af(1,3),af(1,4),ip,b,1,x,2,rcond,r1,r2,w,iw,info)
         CALL CHKXER('DGTSVX',INFot,NOUt,LERr,OK)
         INFot = 16
         CALL DGTSVX('N','N',2,0,a(1,1),a(1,2),a(1,3),af(1,1),af(1,2),  &
     &               af(1,3),af(1,4),ip,b,2,x,1,rcond,r1,r2,w,iw,info)
         CALL CHKXER('DGTSVX',INFot,NOUt,LERr,OK)
!
      ELSEIF ( LSAMEN(2,c2,'PO') ) THEN
!
!        DPOSV
!
         SRNamt = 'DPOSV '
         INFot = 1
         CALL DPOSV('/',0,0,a,1,b,1,info)
         CALL CHKXER('DPOSV ',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DPOSV('U',-1,0,a,1,b,1,info)
         CALL CHKXER('DPOSV ',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DPOSV('U',0,-1,a,1,b,1,info)
         CALL CHKXER('DPOSV ',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL DPOSV('U',2,0,a,1,b,2,info)
         CALL CHKXER('DPOSV ',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL DPOSV('U',2,0,a,2,b,1,info)
         CALL CHKXER('DPOSV ',INFot,NOUt,LERr,OK)
!
!        DPOSVX
!
         SRNamt = 'DPOSVX'
         INFot = 1
         CALL DPOSVX('/','U',0,0,a,1,af,1,eq,c,b,1,x,1,rcond,r1,r2,w,iw,&
     &               info)
         CALL CHKXER('DPOSVX',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DPOSVX('N','/',0,0,a,1,af,1,eq,c,b,1,x,1,rcond,r1,r2,w,iw,&
     &               info)
         CALL CHKXER('DPOSVX',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DPOSVX('N','U',-1,0,a,1,af,1,eq,c,b,1,x,1,rcond,r1,r2,w,  &
     &               iw,info)
         CALL CHKXER('DPOSVX',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL DPOSVX('N','U',0,-1,a,1,af,1,eq,c,b,1,x,1,rcond,r1,r2,w,  &
     &               iw,info)
         CALL CHKXER('DPOSVX',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL DPOSVX('N','U',2,0,a,1,af,2,eq,c,b,2,x,2,rcond,r1,r2,w,iw,&
     &               info)
         CALL CHKXER('DPOSVX',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL DPOSVX('N','U',2,0,a,2,af,1,eq,c,b,2,x,2,rcond,r1,r2,w,iw,&
     &               info)
         CALL CHKXER('DPOSVX',INFot,NOUt,LERr,OK)
         INFot = 9
         eq = '/'
         CALL DPOSVX('F','U',0,0,a,1,af,1,eq,c,b,1,x,1,rcond,r1,r2,w,iw,&
     &               info)
         CALL CHKXER('DPOSVX',INFot,NOUt,LERr,OK)
         INFot = 10
         eq = 'Y'
         CALL DPOSVX('F','U',1,0,a,1,af,1,eq,c,b,1,x,1,rcond,r1,r2,w,iw,&
     &               info)
         CALL CHKXER('DPOSVX',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL DPOSVX('N','U',2,0,a,2,af,2,eq,c,b,1,x,2,rcond,r1,r2,w,iw,&
     &               info)
         CALL CHKXER('DPOSVX',INFot,NOUt,LERr,OK)
         INFot = 14
         CALL DPOSVX('N','U',2,0,a,2,af,2,eq,c,b,2,x,1,rcond,r1,r2,w,iw,&
     &               info)
         CALL CHKXER('DPOSVX',INFot,NOUt,LERr,OK)
!
!        DPOSVXX
!
         n_err_bnds = 3
         nparams = 1
         SRNamt = 'DPOSVXX'
         INFot = 1
         CALL DPOSVXX('/','U',0,0,a,1,af,1,eq,c,b,1,x,1,rcond,rpvgrw,   &
     &                berr,n_err_bnds,err_bnds_n,err_bnds_c,nparams,    &
     &                params,w,iw,info)
         CALL CHKXER('DPOSVXX',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DPOSVXX('N','/',0,0,a,1,af,1,eq,c,b,1,x,1,rcond,rpvgrw,   &
     &                berr,n_err_bnds,err_bnds_n,err_bnds_c,nparams,    &
     &                params,w,iw,info)
         CALL CHKXER('DPOSVXX',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DPOSVXX('N','U',-1,0,a,1,af,1,eq,c,b,1,x,1,rcond,rpvgrw,  &
     &                berr,n_err_bnds,err_bnds_n,err_bnds_c,nparams,    &
     &                params,w,iw,info)
         CALL CHKXER('DPOSVXX',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL DPOSVXX('N','U',0,-1,a,1,af,1,eq,c,b,1,x,1,rcond,rpvgrw,  &
     &                berr,n_err_bnds,err_bnds_n,err_bnds_c,nparams,    &
     &                params,w,iw,info)
         CALL CHKXER('DPOSVXX',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL DPOSVXX('N','U',2,0,a,1,af,2,eq,c,b,2,x,2,rcond,rpvgrw,   &
     &                berr,n_err_bnds,err_bnds_n,err_bnds_c,nparams,    &
     &                params,w,iw,info)
         CALL CHKXER('DPOSVXX',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL DPOSVXX('N','U',2,0,a,2,af,1,eq,c,b,2,x,2,rcond,rpvgrw,   &
     &                berr,n_err_bnds,err_bnds_n,err_bnds_c,nparams,    &
     &                params,w,iw,info)
         CALL CHKXER('DPOSVXX',INFot,NOUt,LERr,OK)
         INFot = 9
         eq = '/'
         CALL DPOSVXX('F','U',0,0,a,1,af,1,eq,c,b,1,x,1,rcond,rpvgrw,   &
     &                berr,n_err_bnds,err_bnds_n,err_bnds_c,nparams,    &
     &                params,w,iw,info)
         CALL CHKXER('DPOSVXX',INFot,NOUt,LERr,OK)
         INFot = 10
         eq = 'Y'
         CALL DPOSVXX('F','U',1,0,a,1,af,1,eq,c,b,1,x,1,rcond,rpvgrw,   &
     &                berr,n_err_bnds,err_bnds_n,err_bnds_c,nparams,    &
     &                params,w,iw,info)
         CALL CHKXER('DPOSVXX',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL DPOSVXX('N','U',2,0,a,2,af,2,eq,c,b,1,x,2,rcond,rpvgrw,   &
     &                berr,n_err_bnds,err_bnds_n,err_bnds_c,nparams,    &
     &                params,w,iw,info)
         CALL CHKXER('DPOSVXX',INFot,NOUt,LERr,OK)
         INFot = 14
         CALL DPOSVXX('N','U',2,0,a,2,af,2,eq,c,b,2,x,1,rcond,rpvgrw,   &
     &                berr,n_err_bnds,err_bnds_n,err_bnds_c,nparams,    &
     &                params,w,iw,info)
         CALL CHKXER('DPOSVXX',INFot,NOUt,LERr,OK)
!
      ELSEIF ( LSAMEN(2,c2,'PP') ) THEN
!
!        DPPSV
!
         SRNamt = 'DPPSV '
         INFot = 1
         CALL DPPSV('/',0,0,a,b,1,info)
         CALL CHKXER('DPPSV ',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DPPSV('U',-1,0,a,b,1,info)
         CALL CHKXER('DPPSV ',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DPPSV('U',0,-1,a,b,1,info)
         CALL CHKXER('DPPSV ',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL DPPSV('U',2,0,a,b,1,info)
         CALL CHKXER('DPPSV ',INFot,NOUt,LERr,OK)
!
!        DPPSVX
!
         SRNamt = 'DPPSVX'
         INFot = 1
         CALL DPPSVX('/','U',0,0,a,af,eq,c,b,1,x,1,rcond,r1,r2,w,iw,    &
     &               info)
         CALL CHKXER('DPPSVX',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DPPSVX('N','/',0,0,a,af,eq,c,b,1,x,1,rcond,r1,r2,w,iw,    &
     &               info)
         CALL CHKXER('DPPSVX',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DPPSVX('N','U',-1,0,a,af,eq,c,b,1,x,1,rcond,r1,r2,w,iw,   &
     &               info)
         CALL CHKXER('DPPSVX',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL DPPSVX('N','U',0,-1,a,af,eq,c,b,1,x,1,rcond,r1,r2,w,iw,   &
     &               info)
         CALL CHKXER('DPPSVX',INFot,NOUt,LERr,OK)
         INFot = 7
         eq = '/'
         CALL DPPSVX('F','U',0,0,a,af,eq,c,b,1,x,1,rcond,r1,r2,w,iw,    &
     &               info)
         CALL CHKXER('DPPSVX',INFot,NOUt,LERr,OK)
         INFot = 8
         eq = 'Y'
         CALL DPPSVX('F','U',1,0,a,af,eq,c,b,1,x,1,rcond,r1,r2,w,iw,    &
     &               info)
         CALL CHKXER('DPPSVX',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL DPPSVX('N','U',2,0,a,af,eq,c,b,1,x,2,rcond,r1,r2,w,iw,    &
     &               info)
         CALL CHKXER('DPPSVX',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL DPPSVX('N','U',2,0,a,af,eq,c,b,2,x,1,rcond,r1,r2,w,iw,    &
     &               info)
         CALL CHKXER('DPPSVX',INFot,NOUt,LERr,OK)
!
      ELSEIF ( LSAMEN(2,c2,'PB') ) THEN
!
!        DPBSV
!
         SRNamt = 'DPBSV '
         INFot = 1
         CALL DPBSV('/',0,0,0,a,1,b,1,info)
         CALL CHKXER('DPBSV ',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DPBSV('U',-1,0,0,a,1,b,1,info)
         CALL CHKXER('DPBSV ',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DPBSV('U',1,-1,0,a,1,b,1,info)
         CALL CHKXER('DPBSV ',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL DPBSV('U',0,0,-1,a,1,b,1,info)
         CALL CHKXER('DPBSV ',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL DPBSV('U',1,1,0,a,1,b,2,info)
         CALL CHKXER('DPBSV ',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL DPBSV('U',2,0,0,a,1,b,1,info)
         CALL CHKXER('DPBSV ',INFot,NOUt,LERr,OK)
!
!        DPBSVX
!
         SRNamt = 'DPBSVX'
         INFot = 1
         CALL DPBSVX('/','U',0,0,0,a,1,af,1,eq,c,b,1,x,1,rcond,r1,r2,w, &
     &               iw,info)
         CALL CHKXER('DPBSVX',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DPBSVX('N','/',0,0,0,a,1,af,1,eq,c,b,1,x,1,rcond,r1,r2,w, &
     &               iw,info)
         CALL CHKXER('DPBSVX',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DPBSVX('N','U',-1,0,0,a,1,af,1,eq,c,b,1,x,1,rcond,r1,r2,w,&
     &               iw,info)
         CALL CHKXER('DPBSVX',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL DPBSVX('N','U',1,-1,0,a,1,af,1,eq,c,b,1,x,1,rcond,r1,r2,w,&
     &               iw,info)
         CALL CHKXER('DPBSVX',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL DPBSVX('N','U',0,0,-1,a,1,af,1,eq,c,b,1,x,1,rcond,r1,r2,w,&
     &               iw,info)
         CALL CHKXER('DPBSVX',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL DPBSVX('N','U',1,1,0,a,1,af,2,eq,c,b,2,x,2,rcond,r1,r2,w, &
     &               iw,info)
         CALL CHKXER('DPBSVX',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL DPBSVX('N','U',1,1,0,a,2,af,1,eq,c,b,2,x,2,rcond,r1,r2,w, &
     &               iw,info)
         CALL CHKXER('DPBSVX',INFot,NOUt,LERr,OK)
         INFot = 10
         eq = '/'
         CALL DPBSVX('F','U',0,0,0,a,1,af,1,eq,c,b,1,x,1,rcond,r1,r2,w, &
     &               iw,info)
         CALL CHKXER('DPBSVX',INFot,NOUt,LERr,OK)
         INFot = 11
         eq = 'Y'
         CALL DPBSVX('F','U',1,0,0,a,1,af,1,eq,c,b,1,x,1,rcond,r1,r2,w, &
     &               iw,info)
         CALL CHKXER('DPBSVX',INFot,NOUt,LERr,OK)
         INFot = 13
         CALL DPBSVX('N','U',2,0,0,a,1,af,1,eq,c,b,1,x,2,rcond,r1,r2,w, &
     &               iw,info)
         CALL CHKXER('DPBSVX',INFot,NOUt,LERr,OK)
         INFot = 15
         CALL DPBSVX('N','U',2,0,0,a,1,af,1,eq,c,b,2,x,1,rcond,r1,r2,w, &
     &               iw,info)
         CALL CHKXER('DPBSVX',INFot,NOUt,LERr,OK)
!
      ELSEIF ( LSAMEN(2,c2,'PT') ) THEN
!
!        DPTSV
!
         SRNamt = 'DPTSV '
         INFot = 1
         CALL DPTSV(-1,0,a(1,1),a(1,2),b,1,info)
         CALL CHKXER('DPTSV ',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DPTSV(0,-1,a(1,1),a(1,2),b,1,info)
         CALL CHKXER('DPTSV ',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL DPTSV(2,0,a(1,1),a(1,2),b,1,info)
         CALL CHKXER('DPTSV ',INFot,NOUt,LERr,OK)
!
!        DPTSVX
!
         SRNamt = 'DPTSVX'
         INFot = 1
         CALL DPTSVX('/',0,0,a(1,1),a(1,2),af(1,1),af(1,2),b,1,x,1,     &
     &               rcond,r1,r2,w,info)
         CALL CHKXER('DPTSVX',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DPTSVX('N',-1,0,a(1,1),a(1,2),af(1,1),af(1,2),b,1,x,1,    &
     &               rcond,r1,r2,w,info)
         CALL CHKXER('DPTSVX',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DPTSVX('N',0,-1,a(1,1),a(1,2),af(1,1),af(1,2),b,1,x,1,    &
     &               rcond,r1,r2,w,info)
         CALL CHKXER('DPTSVX',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL DPTSVX('N',2,0,a(1,1),a(1,2),af(1,1),af(1,2),b,1,x,2,     &
     &               rcond,r1,r2,w,info)
         CALL CHKXER('DPTSVX',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL DPTSVX('N',2,0,a(1,1),a(1,2),af(1,1),af(1,2),b,2,x,1,     &
     &               rcond,r1,r2,w,info)
         CALL CHKXER('DPTSVX',INFot,NOUt,LERr,OK)
!
      ELSEIF ( LSAMEN(2,c2,'SY') ) THEN
!
!        DSYSV
!
         SRNamt = 'DSYSV '
         INFot = 1
         CALL DSYSV('/',0,0,a,1,ip,b,1,w,1,info)
         CALL CHKXER('DSYSV ',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DSYSV('U',-1,0,a,1,ip,b,1,w,1,info)
         CALL CHKXER('DSYSV ',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DSYSV('U',0,-1,a,1,ip,b,1,w,1,info)
         CALL CHKXER('DSYSV ',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL DSYSV('U',2,0,a,1,ip,b,2,w,1,info)
         CALL CHKXER('DSYSV_ROOK',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL DSYSV('U',2,0,a,2,ip,b,1,w,1,info)
         CALL CHKXER('DSYSV ',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL DSYSV('U',0,0,a,1,ip,b,1,w,0,info)
         CALL CHKXER('DSYSV ',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL DSYSV('U',0,0,a,1,ip,b,1,w,-2,info)
         CALL CHKXER('DSYSV ',INFot,NOUt,LERr,OK)
!
!        DSYSVX
!
         SRNamt = 'DSYSVX'
         INFot = 1
         CALL DSYSVX('/','U',0,0,a,1,af,1,ip,b,1,x,1,rcond,r1,r2,w,1,iw,&
     &               info)
         CALL CHKXER('DSYSVX',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DSYSVX('N','/',0,0,a,1,af,1,ip,b,1,x,1,rcond,r1,r2,w,1,iw,&
     &               info)
         CALL CHKXER('DSYSVX',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DSYSVX('N','U',-1,0,a,1,af,1,ip,b,1,x,1,rcond,r1,r2,w,1,  &
     &               iw,info)
         CALL CHKXER('DSYSVX',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL DSYSVX('N','U',0,-1,a,1,af,1,ip,b,1,x,1,rcond,r1,r2,w,1,  &
     &               iw,info)
         CALL CHKXER('DSYSVX',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL DSYSVX('N','U',2,0,a,1,af,2,ip,b,2,x,2,rcond,r1,r2,w,4,iw,&
     &               info)
         CALL CHKXER('DSYSVX',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL DSYSVX('N','U',2,0,a,2,af,1,ip,b,2,x,2,rcond,r1,r2,w,4,iw,&
     &               info)
         CALL CHKXER('DSYSVX',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL DSYSVX('N','U',2,0,a,2,af,2,ip,b,1,x,2,rcond,r1,r2,w,4,iw,&
     &               info)
         CALL CHKXER('DSYSVX',INFot,NOUt,LERr,OK)
         INFot = 13
         CALL DSYSVX('N','U',2,0,a,2,af,2,ip,b,2,x,1,rcond,r1,r2,w,4,iw,&
     &               info)
         CALL CHKXER('DSYSVX',INFot,NOUt,LERr,OK)
         INFot = 18
         CALL DSYSVX('N','U',2,0,a,2,af,2,ip,b,2,x,2,rcond,r1,r2,w,3,iw,&
     &               info)
         CALL CHKXER('DSYSVX',INFot,NOUt,LERr,OK)
!
!        DSYSVXX
!
         n_err_bnds = 3
         nparams = 1
         SRNamt = 'DSYSVXX'
         INFot = 1
         eq = 'N'
         CALL DSYSVXX('/','U',0,0,a,1,af,1,ip,eq,r,b,1,x,1,rcond,rpvgrw,&
     &                berr,n_err_bnds,err_bnds_n,err_bnds_c,nparams,    &
     &                params,w,iw,info)
         CALL CHKXER('DSYSVXX',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DSYSVXX('N','/',0,0,a,1,af,1,ip,eq,r,b,1,x,1,rcond,rpvgrw,&
     &                berr,n_err_bnds,err_bnds_n,err_bnds_c,nparams,    &
     &                params,w,iw,info)
         CALL CHKXER('DSYSVXX',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DSYSVXX('N','U',-1,0,a,1,af,1,ip,eq,r,b,1,x,1,rcond,      &
     &                rpvgrw,berr,n_err_bnds,err_bnds_n,err_bnds_c,     &
     &                nparams,params,w,iw,info)
         CALL CHKXER('DSYSVXX',INFot,NOUt,LERr,OK)
         INFot = 4
         eq = '/'
         CALL DSYSVXX('N','U',0,-1,a,1,af,1,ip,eq,r,b,1,x,1,rcond,      &
     &                rpvgrw,berr,n_err_bnds,err_bnds_n,err_bnds_c,     &
     &                nparams,params,w,iw,info)
         CALL CHKXER('DSYSVXX',INFot,NOUt,LERr,OK)
         eq = 'Y'
         INFot = 6
         CALL DSYSVXX('N','U',2,0,a,1,af,2,ip,eq,r,b,2,x,2,rcond,rpvgrw,&
     &                berr,n_err_bnds,err_bnds_n,err_bnds_c,nparams,    &
     &                params,w,iw,info)
         CALL CHKXER('DSYSVXX',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL DSYSVXX('N','U',2,0,a,2,af,1,ip,eq,r,b,2,x,2,rcond,rpvgrw,&
     &                berr,n_err_bnds,err_bnds_n,err_bnds_c,nparams,    &
     &                params,w,iw,info)
         CALL CHKXER('DSYSVXX',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL DSYSVXX('F','U',2,0,a,2,af,2,ip,'A',r,b,2,x,2,rcond,      &
     &                rpvgrw,berr,n_err_bnds,err_bnds_n,err_bnds_c,     &
     &                nparams,params,w,iw,info)
         CALL CHKXER('DSYSVXX',INFot,NOUt,LERr,OK)
         INFot = 11
         eq = 'Y'
         CALL DSYSVXX('F','U',2,0,a,2,af,2,ip,eq,r,b,2,x,2,rcond,rpvgrw,&
     &                berr,n_err_bnds,err_bnds_n,err_bnds_c,nparams,    &
     &                params,w,iw,info)
         CALL CHKXER('DSYSVXX',INFot,NOUt,LERr,OK)
         INFot = 11
         eq = 'Y'
         r(1) = -ONE
         CALL DSYSVXX('F','U',2,0,a,2,af,2,ip,eq,r,b,2,x,2,rcond,rpvgrw,&
     &                berr,n_err_bnds,err_bnds_n,err_bnds_c,nparams,    &
     &                params,w,iw,info)
         CALL CHKXER('DSYSVXX',INFot,NOUt,LERr,OK)
         INFot = 13
         eq = 'N'
         CALL DSYSVXX('N','U',2,0,a,2,af,2,ip,eq,r,b,1,x,2,rcond,rpvgrw,&
     &                berr,n_err_bnds,err_bnds_n,err_bnds_c,nparams,    &
     &                params,w,iw,info)
         CALL CHKXER('DSYSVXX',INFot,NOUt,LERr,OK)
         INFot = 15
         CALL DSYSVXX('N','U',2,0,a,2,af,2,ip,eq,r,b,2,x,1,rcond,rpvgrw,&
     &                berr,n_err_bnds,err_bnds_n,err_bnds_c,nparams,    &
     &                params,w,iw,info)
         CALL CHKXER('DSYSVXX',INFot,NOUt,LERr,OK)
!
      ELSEIF ( LSAMEN(2,c2,'SR') ) THEN
!
!        DSYSV_ROOK
!
         SRNamt = 'DSYSV_ROOK'
         INFot = 1
         CALL DSYSV_ROOK('/',0,0,a,1,ip,b,1,w,1,info)
         CALL CHKXER('DSYSV_ROOK',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DSYSV_ROOK('U',-1,0,a,1,ip,b,1,w,1,info)
         CALL CHKXER('DSYSV_ROOK',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DSYSV_ROOK('U',0,-1,a,1,ip,b,1,w,1,info)
         CALL CHKXER('DSYSV_ROOK',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL DSYSV_ROOK('U',2,0,a,1,ip,b,2,w,1,info)
         CALL CHKXER('DSYSV_ROOK',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL DSYSV_ROOK('U',2,0,a,2,ip,b,1,w,1,info)
         CALL CHKXER('DSYSV_ROOK',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL DSYSV_ROOK('U',0,0,a,1,ip,b,1,w,0,info)
         CALL CHKXER('DSYSV_ROOK',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL DSYSV_ROOK('U',0,0,a,1,ip,b,1,w,-2,info)
         CALL CHKXER('DSYSV_ROOK',INFot,NOUt,LERr,OK)
!
      ELSEIF ( LSAMEN(2,c2,'SK') ) THEN
!
!        DSYSV_RK
!
!        Test error exits of the driver that uses factorization
!        of a symmetric indefinite matrix with rook
!        (bounded Bunch-Kaufman) pivoting with the new storage
!        format for factors L ( or U) and D.
!
!        L (or U) is stored in A, diagonal of D is stored on the
!        diagonal of A, subdiagonal of D is stored in a separate array E.
!
         SRNamt = 'DSYSV_RK'
         INFot = 1
         CALL DSYSV_RK('/',0,0,a,1,e,ip,b,1,w,1,info)
         CALL CHKXER('DSYSV_RK',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DSYSV_RK('U',-1,0,a,1,e,ip,b,1,w,1,info)
         CALL CHKXER('DSYSV_RK',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DSYSV_RK('U',0,-1,a,1,e,ip,b,1,w,1,info)
         CALL CHKXER('DSYSV_RK',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL DSYSV_RK('U',2,0,a,1,e,ip,b,2,w,1,info)
         CALL CHKXER('DSYSV_RK',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL DSYSV_RK('U',2,0,a,2,e,ip,b,1,w,1,info)
         CALL CHKXER('DSYSV_RK',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL DSYSV_RK('U',0,0,a,1,e,ip,b,1,w,0,info)
         CALL CHKXER('DSYSV_RK',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL DSYSV_RK('U',0,0,a,1,e,ip,b,1,w,-2,info)
         CALL CHKXER('DSYSV_RK',INFot,NOUt,LERr,OK)
!
      ELSEIF ( LSAMEN(2,c2,'SP') ) THEN
!
!        DSPSV
!
         SRNamt = 'DSPSV '
         INFot = 1
         CALL DSPSV('/',0,0,a,ip,b,1,info)
         CALL CHKXER('DSPSV ',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DSPSV('U',-1,0,a,ip,b,1,info)
         CALL CHKXER('DSPSV ',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DSPSV('U',0,-1,a,ip,b,1,info)
         CALL CHKXER('DSPSV ',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL DSPSV('U',2,0,a,ip,b,1,info)
         CALL CHKXER('DSPSV ',INFot,NOUt,LERr,OK)
!
!        DSPSVX
!
         SRNamt = 'DSPSVX'
         INFot = 1
         CALL DSPSVX('/','U',0,0,a,af,ip,b,1,x,1,rcond,r1,r2,w,iw,info)
         CALL CHKXER('DSPSVX',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DSPSVX('N','/',0,0,a,af,ip,b,1,x,1,rcond,r1,r2,w,iw,info)
         CALL CHKXER('DSPSVX',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DSPSVX('N','U',-1,0,a,af,ip,b,1,x,1,rcond,r1,r2,w,iw,info)
         CALL CHKXER('DSPSVX',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL DSPSVX('N','U',0,-1,a,af,ip,b,1,x,1,rcond,r1,r2,w,iw,info)
         CALL CHKXER('DSPSVX',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL DSPSVX('N','U',2,0,a,af,ip,b,1,x,2,rcond,r1,r2,w,iw,info)
         CALL CHKXER('DSPSVX',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL DSPSVX('N','U',2,0,a,af,ip,b,2,x,1,rcond,r1,r2,w,iw,info)
         CALL CHKXER('DSPSVX',INFot,NOUt,LERr,OK)
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
!     End of derrvx
!
      END SUBROUTINE derrvx
