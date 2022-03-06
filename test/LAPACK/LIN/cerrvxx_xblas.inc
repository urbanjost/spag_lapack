!*==cerrvx.f90  processed by SPAG 7.51RB at 17:34 on  4 Mar 2022
!> \brief \b CERRVXX
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE cerrvx( PATH, NUNIT )
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
!> cerrvx tests the error exits for the COMPLEX driver routines
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
!> \ingroup complex_lin
!
!  =====================================================================
      SUBROUTINE cerrvx(Path,Nunit)
      IMPLICIT NONE
!*--CERRVX59
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
      PARAMETER (ONE=1.0E+0)
!     ..
!     .. Local Scalars ..
      CHARACTER eq
      CHARACTER*2 c2
      INTEGER i , info , j , n_err_bnds , nparams
      REAL rcond , rpvgrw , berr
!     ..
!     .. Local Arrays ..
      INTEGER ip(NMAX)
      REAL c(NMAX) , r(NMAX) , r1(NMAX) , r2(NMAX) , rf(NMAX) , rw(NMAX)&
     &     , err_bnds_n(NMAX,3) , err_bnds_c(NMAX,3) , params(1)
      COMPLEX a(NMAX,NMAX) , af(NMAX,NMAX) , b(NMAX) , e(NMAX) ,        &
     &        w(2*NMAX) , x(NMAX)
!     ..
!     .. External Functions ..
      LOGICAL LSAMEN
EXTERNAL    LSAMEN
!     ..
!     .. External Subroutines ..
EXTERNAL    CGBSV,   CGBSVX,   CGESV,   CGESVX,   CGTSV,   CGTSVX
EXTERNAL    CHESV,   CHESV_RK,   CHESV_ROOK,   CHESVX,   CHKXER,   CHPSV
EXTERNAL    CHPSVX,   CPBSV,   CPBSVX,   CPOSV,   CPOSVX,   CPPSV
EXTERNAL    CPPSVX,   CPTSV,   CPTSVX,   CSPSV,   CSPSVX,   CSYSV
EXTERNAL    CSYSV_RK,   CSYSV_ROOK,   CSYSVX,   CGESVXX,   CPOSVXX
EXTERNAL    CSYSVXX,   CHESVXX,   CGBSVXX
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
      INTRINSIC CMPLX , REAL
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
            a(i,j) = CMPLX(1./REAL(i+j),-1./REAL(i+j))
            af(i,j) = CMPLX(1./REAL(i+j),-1./REAL(i+j))
         ENDDO
         b(j) = 0.E+0
         e(j) = 0E+0
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
!        CGESV
!
         SRNamt = 'CGESV '
         INFot = 1
         CALL CGESV(-1,0,a,1,ip,b,1,info)
         CALL CHKXER('CGESV ',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CGESV(0,-1,a,1,ip,b,1,info)
         CALL CHKXER('CGESV ',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CGESV(2,1,a,1,ip,b,2,info)
         CALL CHKXER('CGESV ',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL CGESV(2,1,a,2,ip,b,1,info)
         CALL CHKXER('CGESV ',INFot,NOUt,LERr,OK)
!
!        CGESVX
!
         SRNamt = 'CGESVX'
         INFot = 1
         CALL CGESVX('/','N',0,0,a,1,af,1,ip,eq,r,c,b,1,x,1,rcond,r1,r2,&
     &               w,rw,info)
         CALL CHKXER('CGESVX',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CGESVX('N','/',0,0,a,1,af,1,ip,eq,r,c,b,1,x,1,rcond,r1,r2,&
     &               w,rw,info)
         CALL CHKXER('CGESVX',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CGESVX('N','N',-1,0,a,1,af,1,ip,eq,r,c,b,1,x,1,rcond,r1,  &
     &               r2,w,rw,info)
         CALL CHKXER('CGESVX',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CGESVX('N','N',0,-1,a,1,af,1,ip,eq,r,c,b,1,x,1,rcond,r1,  &
     &               r2,w,rw,info)
         CALL CHKXER('CGESVX',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL CGESVX('N','N',2,1,a,1,af,2,ip,eq,r,c,b,2,x,2,rcond,r1,r2,&
     &               w,rw,info)
         CALL CHKXER('CGESVX',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL CGESVX('N','N',2,1,a,2,af,1,ip,eq,r,c,b,2,x,2,rcond,r1,r2,&
     &               w,rw,info)
         CALL CHKXER('CGESVX',INFot,NOUt,LERr,OK)
         INFot = 10
         eq = '/'
         CALL CGESVX('F','N',0,0,a,1,af,1,ip,eq,r,c,b,1,x,1,rcond,r1,r2,&
     &               w,rw,info)
         CALL CHKXER('CGESVX',INFot,NOUt,LERr,OK)
         INFot = 11
         eq = 'R'
         CALL CGESVX('F','N',1,0,a,1,af,1,ip,eq,r,c,b,1,x,1,rcond,r1,r2,&
     &               w,rw,info)
         CALL CHKXER('CGESVX',INFot,NOUt,LERr,OK)
         INFot = 12
         eq = 'C'
         CALL CGESVX('F','N',1,0,a,1,af,1,ip,eq,r,c,b,1,x,1,rcond,r1,r2,&
     &               w,rw,info)
         CALL CHKXER('CGESVX',INFot,NOUt,LERr,OK)
         INFot = 14
         CALL CGESVX('N','N',2,1,a,2,af,2,ip,eq,r,c,b,1,x,2,rcond,r1,r2,&
     &               w,rw,info)
         CALL CHKXER('CGESVX',INFot,NOUt,LERr,OK)
         INFot = 16
         CALL CGESVX('N','N',2,1,a,2,af,2,ip,eq,r,c,b,2,x,1,rcond,r1,r2,&
     &               w,rw,info)
         CALL CHKXER('CGESVX',INFot,NOUt,LERr,OK)
!
!        CGESVXX
!
         n_err_bnds = 3
         nparams = 1
         SRNamt = 'CGESVXX'
         INFot = 1
         CALL CGESVXX('/','N',0,0,a,1,af,1,ip,eq,r,c,b,1,x,1,rcond,     &
     &                rpvgrw,berr,n_err_bnds,err_bnds_n,err_bnds_c,     &
     &                nparams,params,w,rw,info)
         CALL CHKXER('CGESVXX',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CGESVXX('N','/',0,0,a,1,af,1,ip,eq,r,c,b,1,x,1,rcond,     &
     &                rpvgrw,berr,n_err_bnds,err_bnds_n,err_bnds_c,     &
     &                nparams,params,w,rw,info)
         CALL CHKXER('CGESVXX',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CGESVXX('N','N',-1,0,a,1,af,1,ip,eq,r,c,b,1,x,1,rcond,    &
     &                rpvgrw,berr,n_err_bnds,err_bnds_n,err_bnds_c,     &
     &                nparams,params,w,rw,info)
         CALL CHKXER('CGESVXX',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CGESVXX('N','N',0,-1,a,1,af,1,ip,eq,r,c,b,1,x,1,rcond,    &
     &                rpvgrw,berr,n_err_bnds,err_bnds_n,err_bnds_c,     &
     &                nparams,params,w,rw,info)
         CALL CHKXER('CGESVXX',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL CGESVXX('N','N',2,1,a,1,af,2,ip,eq,r,c,b,2,x,2,rcond,     &
     &                rpvgrw,berr,n_err_bnds,err_bnds_n,err_bnds_c,     &
     &                nparams,params,w,rw,info)
         CALL CHKXER('CGESVXX',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL CGESVXX('N','N',2,1,a,2,af,1,ip,eq,r,c,b,2,x,2,rcond,     &
     &                rpvgrw,berr,n_err_bnds,err_bnds_n,err_bnds_c,     &
     &                nparams,params,w,rw,info)
         CALL CHKXER('CGESVXX',INFot,NOUt,LERr,OK)
         INFot = 10
         eq = '/'
         CALL CGESVXX('F','N',0,0,a,1,af,1,ip,eq,r,c,b,1,x,1,rcond,     &
     &                rpvgrw,berr,n_err_bnds,err_bnds_n,err_bnds_c,     &
     &                nparams,params,w,rw,info)
         CALL CHKXER('CGESVXX',INFot,NOUt,LERr,OK)
         INFot = 11
         eq = 'R'
         CALL CGESVXX('F','N',1,0,a,1,af,1,ip,eq,r,c,b,1,x,1,rcond,     &
     &                rpvgrw,berr,n_err_bnds,err_bnds_n,err_bnds_c,     &
     &                nparams,params,w,rw,info)
         CALL CHKXER('CGESVXX',INFot,NOUt,LERr,OK)
         INFot = 12
         eq = 'C'
         CALL CGESVXX('F','N',1,0,a,1,af,1,ip,eq,r,c,b,1,x,1,rcond,     &
     &                rpvgrw,berr,n_err_bnds,err_bnds_n,err_bnds_c,     &
     &                nparams,params,w,rw,info)
         CALL CHKXER('CGESVXX',INFot,NOUt,LERr,OK)
         INFot = 14
         CALL CGESVXX('N','N',2,1,a,2,af,2,ip,eq,r,c,b,1,x,2,rcond,     &
     &                rpvgrw,berr,n_err_bnds,err_bnds_n,err_bnds_c,     &
     &                nparams,params,w,rw,info)
         CALL CHKXER('CGESVXX',INFot,NOUt,LERr,OK)
         INFot = 16
         CALL CGESVXX('N','N',2,1,a,2,af,2,ip,eq,r,c,b,2,x,1,rcond,     &
     &                rpvgrw,berr,n_err_bnds,err_bnds_n,err_bnds_c,     &
     &                nparams,params,w,rw,info)
         CALL CHKXER('CGESVXX',INFot,NOUt,LERr,OK)
!
      ELSEIF ( LSAMEN(2,c2,'GB') ) THEN
!
!        CGBSV
!
         SRNamt = 'CGBSV '
         INFot = 1
         CALL CGBSV(-1,0,0,0,a,1,ip,b,1,info)
         CALL CHKXER('CGBSV ',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CGBSV(1,-1,0,0,a,1,ip,b,1,info)
         CALL CHKXER('CGBSV ',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CGBSV(1,0,-1,0,a,1,ip,b,1,info)
         CALL CHKXER('CGBSV ',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CGBSV(0,0,0,-1,a,1,ip,b,1,info)
         CALL CHKXER('CGBSV ',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL CGBSV(1,1,1,0,a,3,ip,b,1,info)
         CALL CHKXER('CGBSV ',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL CGBSV(2,0,0,0,a,1,ip,b,1,info)
         CALL CHKXER('CGBSV ',INFot,NOUt,LERr,OK)
!
!        CGBSVX
!
         SRNamt = 'CGBSVX'
         INFot = 1
         CALL CGBSVX('/','N',0,0,0,0,a,1,af,1,ip,eq,r,c,b,1,x,1,rcond,  &
     &               r1,r2,w,rw,info)
         CALL CHKXER('CGBSVX',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CGBSVX('N','/',0,0,0,0,a,1,af,1,ip,eq,r,c,b,1,x,1,rcond,  &
     &               r1,r2,w,rw,info)
         CALL CHKXER('CGBSVX',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CGBSVX('N','N',-1,0,0,0,a,1,af,1,ip,eq,r,c,b,1,x,1,rcond, &
     &               r1,r2,w,rw,info)
         CALL CHKXER('CGBSVX',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CGBSVX('N','N',1,-1,0,0,a,1,af,1,ip,eq,r,c,b,1,x,1,rcond, &
     &               r1,r2,w,rw,info)
         CALL CHKXER('CGBSVX',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL CGBSVX('N','N',1,0,-1,0,a,1,af,1,ip,eq,r,c,b,1,x,1,rcond, &
     &               r1,r2,w,rw,info)
         CALL CHKXER('CGBSVX',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL CGBSVX('N','N',0,0,0,-1,a,1,af,1,ip,eq,r,c,b,1,x,1,rcond, &
     &               r1,r2,w,rw,info)
         CALL CHKXER('CGBSVX',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL CGBSVX('N','N',1,1,1,0,a,2,af,4,ip,eq,r,c,b,1,x,1,rcond,  &
     &               r1,r2,w,rw,info)
         CALL CHKXER('CGBSVX',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL CGBSVX('N','N',1,1,1,0,a,3,af,3,ip,eq,r,c,b,1,x,1,rcond,  &
     &               r1,r2,w,rw,info)
         CALL CHKXER('CGBSVX',INFot,NOUt,LERr,OK)
         INFot = 12
         eq = '/'
         CALL CGBSVX('F','N',0,0,0,0,a,1,af,1,ip,eq,r,c,b,1,x,1,rcond,  &
     &               r1,r2,w,rw,info)
         CALL CHKXER('CGBSVX',INFot,NOUt,LERr,OK)
         INFot = 13
         eq = 'R'
         CALL CGBSVX('F','N',1,0,0,0,a,1,af,1,ip,eq,r,c,b,1,x,1,rcond,  &
     &               r1,r2,w,rw,info)
         CALL CHKXER('CGBSVX',INFot,NOUt,LERr,OK)
         INFot = 14
         eq = 'C'
         CALL CGBSVX('F','N',1,0,0,0,a,1,af,1,ip,eq,r,c,b,1,x,1,rcond,  &
     &               r1,r2,w,rw,info)
         CALL CHKXER('CGBSVX',INFot,NOUt,LERr,OK)
         INFot = 16
         CALL CGBSVX('N','N',2,0,0,0,a,1,af,1,ip,eq,r,c,b,1,x,2,rcond,  &
     &               r1,r2,w,rw,info)
         CALL CHKXER('CGBSVX',INFot,NOUt,LERr,OK)
         INFot = 18
         CALL CGBSVX('N','N',2,0,0,0,a,1,af,1,ip,eq,r,c,b,2,x,1,rcond,  &
     &               r1,r2,w,rw,info)
         CALL CHKXER('CGBSVX',INFot,NOUt,LERr,OK)
!
!        CGBSVXX
!
         n_err_bnds = 3
         nparams = 1
         SRNamt = 'CGBSVXX'
         INFot = 1
         CALL CGBSVXX('/','N',0,0,0,0,a,1,af,1,ip,eq,r,c,b,1,x,1,rcond, &
     &                rpvgrw,berr,n_err_bnds,err_bnds_n,err_bnds_c,     &
     &                nparams,params,w,rw,info)
         CALL CHKXER('CGBSVXX',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CGBSVXX('N','/',0,1,1,0,a,1,af,1,ip,eq,r,c,b,1,x,1,rcond, &
     &                rpvgrw,berr,n_err_bnds,err_bnds_n,err_bnds_c,     &
     &                nparams,params,w,rw,info)
         CALL CHKXER('CGBSVXX',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CGBSVXX('N','N',-1,1,1,0,a,1,af,1,ip,eq,r,c,b,1,x,1,rcond,&
     &                rpvgrw,berr,n_err_bnds,err_bnds_n,err_bnds_c,     &
     &                nparams,params,w,rw,info)
         CALL CHKXER('CGBSVXX',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CGBSVXX('N','N',2,-1,1,0,a,1,af,1,ip,eq,r,c,b,1,x,1,rcond,&
     &                rpvgrw,berr,n_err_bnds,err_bnds_n,err_bnds_c,     &
     &                nparams,params,w,rw,info)
         CALL CHKXER('CGBSVXX',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL CGBSVXX('N','N',2,1,-1,0,a,1,af,1,ip,eq,r,c,b,1,x,1,rcond,&
     &                rpvgrw,berr,n_err_bnds,err_bnds_n,err_bnds_c,     &
     &                nparams,params,w,rw,info)
         CALL CHKXER('CGBSVXX',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL CGBSVXX('N','N',0,1,1,-1,a,1,af,1,ip,eq,r,c,b,1,x,1,rcond,&
     &                rpvgrw,berr,n_err_bnds,err_bnds_n,err_bnds_c,     &
     &                nparams,params,w,rw,info)
         CALL CHKXER('CGBSVXX',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL CGBSVXX('N','N',2,1,1,1,a,2,af,2,ip,eq,r,c,b,2,x,2,rcond, &
     &                rpvgrw,berr,n_err_bnds,err_bnds_n,err_bnds_c,     &
     &                nparams,params,w,rw,info)
         CALL CHKXER('CGBSVXX',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL CGBSVXX('N','N',2,1,1,1,a,3,af,3,ip,eq,r,c,b,2,x,2,rcond, &
     &                rpvgrw,berr,n_err_bnds,err_bnds_n,err_bnds_c,     &
     &                nparams,params,w,rw,info)
         CALL CHKXER('CGBSVXX',INFot,NOUt,LERr,OK)
         INFot = 12
         eq = '/'
         CALL CGBSVXX('F','N',0,1,1,0,a,3,af,4,ip,eq,r,c,b,1,x,1,rcond, &
     &                rpvgrw,berr,n_err_bnds,err_bnds_n,err_bnds_c,     &
     &                nparams,params,w,rw,info)
         CALL CHKXER('CGBSVXX',INFot,NOUt,LERr,OK)
         INFot = 13
         eq = 'R'
         CALL CGBSVXX('F','N',1,1,1,0,a,3,af,4,ip,eq,r,c,b,1,x,1,rcond, &
     &                rpvgrw,berr,n_err_bnds,err_bnds_n,err_bnds_c,     &
     &                nparams,params,w,rw,info)
         CALL CHKXER('CGBSVXX',INFot,NOUt,LERr,OK)
         INFot = 14
         eq = 'C'
         CALL CGBSVXX('F','N',1,1,1,0,a,3,af,4,ip,eq,r,c,b,1,x,1,rcond, &
     &                rpvgrw,berr,n_err_bnds,err_bnds_n,err_bnds_c,     &
     &                nparams,params,w,rw,info)
         CALL CHKXER('CGBSVXX',INFot,NOUt,LERr,OK)
         INFot = 15
         CALL CGBSVXX('N','N',2,1,1,1,a,3,af,4,ip,eq,r,c,b,1,x,2,rcond, &
     &                rpvgrw,berr,n_err_bnds,err_bnds_n,err_bnds_c,     &
     &                nparams,params,w,rw,info)
         CALL CHKXER('CGBSVXX',INFot,NOUt,LERr,OK)
         INFot = 16
         CALL CGBSVXX('N','N',2,1,1,1,a,3,af,4,ip,eq,r,c,b,2,x,1,rcond, &
     &                rpvgrw,berr,n_err_bnds,err_bnds_n,err_bnds_c,     &
     &                nparams,params,w,rw,info)
         CALL CHKXER('CGBSVXX',INFot,NOUt,LERr,OK)
!
      ELSEIF ( LSAMEN(2,c2,'GT') ) THEN
!
!        CGTSV
!
         SRNamt = 'CGTSV '
         INFot = 1
         CALL CGTSV(-1,0,a(1,1),a(1,2),a(1,3),b,1,info)
         CALL CHKXER('CGTSV ',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CGTSV(0,-1,a(1,1),a(1,2),a(1,3),b,1,info)
         CALL CHKXER('CGTSV ',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL CGTSV(2,0,a(1,1),a(1,2),a(1,3),b,1,info)
         CALL CHKXER('CGTSV ',INFot,NOUt,LERr,OK)
!
!        CGTSVX
!
         SRNamt = 'CGTSVX'
         INFot = 1
         CALL CGTSVX('/','N',0,0,a(1,1),a(1,2),a(1,3),af(1,1),af(1,2),  &
     &               af(1,3),af(1,4),ip,b,1,x,1,rcond,r1,r2,w,rw,info)
         CALL CHKXER('CGTSVX',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CGTSVX('N','/',0,0,a(1,1),a(1,2),a(1,3),af(1,1),af(1,2),  &
     &               af(1,3),af(1,4),ip,b,1,x,1,rcond,r1,r2,w,rw,info)
         CALL CHKXER('CGTSVX',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CGTSVX('N','N',-1,0,a(1,1),a(1,2),a(1,3),af(1,1),af(1,2), &
     &               af(1,3),af(1,4),ip,b,1,x,1,rcond,r1,r2,w,rw,info)
         CALL CHKXER('CGTSVX',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CGTSVX('N','N',0,-1,a(1,1),a(1,2),a(1,3),af(1,1),af(1,2), &
     &               af(1,3),af(1,4),ip,b,1,x,1,rcond,r1,r2,w,rw,info)
         CALL CHKXER('CGTSVX',INFot,NOUt,LERr,OK)
         INFot = 14
         CALL CGTSVX('N','N',2,0,a(1,1),a(1,2),a(1,3),af(1,1),af(1,2),  &
     &               af(1,3),af(1,4),ip,b,1,x,2,rcond,r1,r2,w,rw,info)
         CALL CHKXER('CGTSVX',INFot,NOUt,LERr,OK)
         INFot = 16
         CALL CGTSVX('N','N',2,0,a(1,1),a(1,2),a(1,3),af(1,1),af(1,2),  &
     &               af(1,3),af(1,4),ip,b,2,x,1,rcond,r1,r2,w,rw,info)
         CALL CHKXER('CGTSVX',INFot,NOUt,LERr,OK)
!
      ELSEIF ( LSAMEN(2,c2,'PO') ) THEN
!
!        CPOSV
!
         SRNamt = 'CPOSV '
         INFot = 1
         CALL CPOSV('/',0,0,a,1,b,1,info)
         CALL CHKXER('CPOSV ',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CPOSV('U',-1,0,a,1,b,1,info)
         CALL CHKXER('CPOSV ',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CPOSV('U',0,-1,a,1,b,1,info)
         CALL CHKXER('CPOSV ',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL CPOSV('U',2,0,a,1,b,2,info)
         CALL CHKXER('CPOSV ',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL CPOSV('U',2,0,a,2,b,1,info)
         CALL CHKXER('CPOSV ',INFot,NOUt,LERr,OK)
!
!        CPOSVX
!
         SRNamt = 'CPOSVX'
         INFot = 1
         CALL CPOSVX('/','U',0,0,a,1,af,1,eq,c,b,1,x,1,rcond,r1,r2,w,rw,&
     &               info)
         CALL CHKXER('CPOSVX',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CPOSVX('N','/',0,0,a,1,af,1,eq,c,b,1,x,1,rcond,r1,r2,w,rw,&
     &               info)
         CALL CHKXER('CPOSVX',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CPOSVX('N','U',-1,0,a,1,af,1,eq,c,b,1,x,1,rcond,r1,r2,w,  &
     &               rw,info)
         CALL CHKXER('CPOSVX',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CPOSVX('N','U',0,-1,a,1,af,1,eq,c,b,1,x,1,rcond,r1,r2,w,  &
     &               rw,info)
         CALL CHKXER('CPOSVX',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL CPOSVX('N','U',2,0,a,1,af,2,eq,c,b,2,x,2,rcond,r1,r2,w,rw,&
     &               info)
         CALL CHKXER('CPOSVX',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL CPOSVX('N','U',2,0,a,2,af,1,eq,c,b,2,x,2,rcond,r1,r2,w,rw,&
     &               info)
         CALL CHKXER('CPOSVX',INFot,NOUt,LERr,OK)
         INFot = 9
         eq = '/'
         CALL CPOSVX('F','U',0,0,a,1,af,1,eq,c,b,1,x,1,rcond,r1,r2,w,rw,&
     &               info)
         CALL CHKXER('CPOSVX',INFot,NOUt,LERr,OK)
         INFot = 10
         eq = 'Y'
         CALL CPOSVX('F','U',1,0,a,1,af,1,eq,c,b,1,x,1,rcond,r1,r2,w,rw,&
     &               info)
         CALL CHKXER('CPOSVX',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL CPOSVX('N','U',2,0,a,2,af,2,eq,c,b,1,x,2,rcond,r1,r2,w,rw,&
     &               info)
         CALL CHKXER('CPOSVX',INFot,NOUt,LERr,OK)
         INFot = 14
         CALL CPOSVX('N','U',2,0,a,2,af,2,eq,c,b,2,x,1,rcond,r1,r2,w,rw,&
     &               info)
         CALL CHKXER('CPOSVX',INFot,NOUt,LERr,OK)
!
!        CPOSVXX
!
         n_err_bnds = 3
         nparams = 1
         SRNamt = 'CPOSVXX'
         INFot = 1
         CALL CPOSVXX('/','U',0,0,a,1,af,1,eq,c,b,1,x,1,rcond,rpvgrw,   &
     &                berr,n_err_bnds,err_bnds_n,err_bnds_c,nparams,    &
     &                params,w,rw,info)
         CALL CHKXER('CPOSVXX',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CPOSVXX('N','/',0,0,a,1,af,1,eq,c,b,1,x,1,rcond,rpvgrw,   &
     &                berr,n_err_bnds,err_bnds_n,err_bnds_c,nparams,    &
     &                params,w,rw,info)
         CALL CHKXER('CPOSVXX',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CPOSVXX('N','U',-1,0,a,1,af,1,eq,c,b,1,x,1,rcond,rpvgrw,  &
     &                berr,n_err_bnds,err_bnds_n,err_bnds_c,nparams,    &
     &                params,w,rw,info)
         CALL CHKXER('CPOSVXX',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CPOSVXX('N','U',0,-1,a,1,af,1,eq,c,b,1,x,1,rcond,rpvgrw,  &
     &                berr,n_err_bnds,err_bnds_n,err_bnds_c,nparams,    &
     &                params,w,rw,info)
         CALL CHKXER('CPOSVXX',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL CPOSVXX('N','U',2,0,a,1,af,2,eq,c,b,2,x,2,rcond,rpvgrw,   &
     &                berr,n_err_bnds,err_bnds_n,err_bnds_c,nparams,    &
     &                params,w,rw,info)
         CALL CHKXER('CPOSVXX',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL CPOSVXX('N','U',2,0,a,2,af,1,eq,c,b,2,x,2,rcond,rpvgrw,   &
     &                berr,n_err_bnds,err_bnds_n,err_bnds_c,nparams,    &
     &                params,w,rw,info)
         CALL CHKXER('CPOSVXX',INFot,NOUt,LERr,OK)
         INFot = 9
         eq = '/'
         CALL CPOSVXX('F','U',0,0,a,1,af,1,eq,c,b,1,x,1,rcond,rpvgrw,   &
     &                berr,n_err_bnds,err_bnds_n,err_bnds_c,nparams,    &
     &                params,w,rw,info)
         CALL CHKXER('CPOSVXX',INFot,NOUt,LERr,OK)
         INFot = 10
         eq = 'Y'
         CALL CPOSVXX('F','U',1,0,a,1,af,1,eq,c,b,1,x,1,rcond,rpvgrw,   &
     &                berr,n_err_bnds,err_bnds_n,err_bnds_c,nparams,    &
     &                params,w,rw,info)
         CALL CHKXER('CPOSVXX',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL CPOSVXX('N','U',2,0,a,2,af,2,eq,c,b,1,x,2,rcond,rpvgrw,   &
     &                berr,n_err_bnds,err_bnds_n,err_bnds_c,nparams,    &
     &                params,w,rw,info)
         CALL CHKXER('CPOSVXX',INFot,NOUt,LERr,OK)
         INFot = 14
         CALL CPOSVXX('N','U',2,0,a,2,af,2,eq,c,b,2,x,1,rcond,rpvgrw,   &
     &                berr,n_err_bnds,err_bnds_n,err_bnds_c,nparams,    &
     &                params,w,rw,info)
         CALL CHKXER('CPOSVXX',INFot,NOUt,LERr,OK)
!
      ELSEIF ( LSAMEN(2,c2,'PP') ) THEN
!
!        CPPSV
!
         SRNamt = 'CPPSV '
         INFot = 1
         CALL CPPSV('/',0,0,a,b,1,info)
         CALL CHKXER('CPPSV ',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CPPSV('U',-1,0,a,b,1,info)
         CALL CHKXER('CPPSV ',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CPPSV('U',0,-1,a,b,1,info)
         CALL CHKXER('CPPSV ',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL CPPSV('U',2,0,a,b,1,info)
         CALL CHKXER('CPPSV ',INFot,NOUt,LERr,OK)
!
!        CPPSVX
!
         SRNamt = 'CPPSVX'
         INFot = 1
         CALL CPPSVX('/','U',0,0,a,af,eq,c,b,1,x,1,rcond,r1,r2,w,rw,    &
     &               info)
         CALL CHKXER('CPPSVX',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CPPSVX('N','/',0,0,a,af,eq,c,b,1,x,1,rcond,r1,r2,w,rw,    &
     &               info)
         CALL CHKXER('CPPSVX',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CPPSVX('N','U',-1,0,a,af,eq,c,b,1,x,1,rcond,r1,r2,w,rw,   &
     &               info)
         CALL CHKXER('CPPSVX',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CPPSVX('N','U',0,-1,a,af,eq,c,b,1,x,1,rcond,r1,r2,w,rw,   &
     &               info)
         CALL CHKXER('CPPSVX',INFot,NOUt,LERr,OK)
         INFot = 7
         eq = '/'
         CALL CPPSVX('F','U',0,0,a,af,eq,c,b,1,x,1,rcond,r1,r2,w,rw,    &
     &               info)
         CALL CHKXER('CPPSVX',INFot,NOUt,LERr,OK)
         INFot = 8
         eq = 'Y'
         CALL CPPSVX('F','U',1,0,a,af,eq,c,b,1,x,1,rcond,r1,r2,w,rw,    &
     &               info)
         CALL CHKXER('CPPSVX',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL CPPSVX('N','U',2,0,a,af,eq,c,b,1,x,2,rcond,r1,r2,w,rw,    &
     &               info)
         CALL CHKXER('CPPSVX',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL CPPSVX('N','U',2,0,a,af,eq,c,b,2,x,1,rcond,r1,r2,w,rw,    &
     &               info)
         CALL CHKXER('CPPSVX',INFot,NOUt,LERr,OK)
!
      ELSEIF ( LSAMEN(2,c2,'PB') ) THEN
!
!        CPBSV
!
         SRNamt = 'CPBSV '
         INFot = 1
         CALL CPBSV('/',0,0,0,a,1,b,1,info)
         CALL CHKXER('CPBSV ',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CPBSV('U',-1,0,0,a,1,b,1,info)
         CALL CHKXER('CPBSV ',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CPBSV('U',1,-1,0,a,1,b,1,info)
         CALL CHKXER('CPBSV ',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CPBSV('U',0,0,-1,a,1,b,1,info)
         CALL CHKXER('CPBSV ',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL CPBSV('U',1,1,0,a,1,b,2,info)
         CALL CHKXER('CPBSV ',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL CPBSV('U',2,0,0,a,1,b,1,info)
         CALL CHKXER('CPBSV ',INFot,NOUt,LERr,OK)
!
!        CPBSVX
!
         SRNamt = 'CPBSVX'
         INFot = 1
         CALL CPBSVX('/','U',0,0,0,a,1,af,1,eq,c,b,1,x,1,rcond,r1,r2,w, &
     &               rw,info)
         CALL CHKXER('CPBSVX',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CPBSVX('N','/',0,0,0,a,1,af,1,eq,c,b,1,x,1,rcond,r1,r2,w, &
     &               rw,info)
         CALL CHKXER('CPBSVX',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CPBSVX('N','U',-1,0,0,a,1,af,1,eq,c,b,1,x,1,rcond,r1,r2,w,&
     &               rw,info)
         CALL CHKXER('CPBSVX',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CPBSVX('N','U',1,-1,0,a,1,af,1,eq,c,b,1,x,1,rcond,r1,r2,w,&
     &               rw,info)
         CALL CHKXER('CPBSVX',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL CPBSVX('N','U',0,0,-1,a,1,af,1,eq,c,b,1,x,1,rcond,r1,r2,w,&
     &               rw,info)
         CALL CHKXER('CPBSVX',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL CPBSVX('N','U',1,1,0,a,1,af,2,eq,c,b,2,x,2,rcond,r1,r2,w, &
     &               rw,info)
         CALL CHKXER('CPBSVX',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL CPBSVX('N','U',1,1,0,a,2,af,1,eq,c,b,2,x,2,rcond,r1,r2,w, &
     &               rw,info)
         CALL CHKXER('CPBSVX',INFot,NOUt,LERr,OK)
         INFot = 10
         eq = '/'
         CALL CPBSVX('F','U',0,0,0,a,1,af,1,eq,c,b,1,x,1,rcond,r1,r2,w, &
     &               rw,info)
         CALL CHKXER('CPBSVX',INFot,NOUt,LERr,OK)
         INFot = 11
         eq = 'Y'
         CALL CPBSVX('F','U',1,0,0,a,1,af,1,eq,c,b,1,x,1,rcond,r1,r2,w, &
     &               rw,info)
         CALL CHKXER('CPBSVX',INFot,NOUt,LERr,OK)
         INFot = 13
         CALL CPBSVX('N','U',2,0,0,a,1,af,1,eq,c,b,1,x,2,rcond,r1,r2,w, &
     &               rw,info)
         CALL CHKXER('CPBSVX',INFot,NOUt,LERr,OK)
         INFot = 15
         CALL CPBSVX('N','U',2,0,0,a,1,af,1,eq,c,b,2,x,1,rcond,r1,r2,w, &
     &               rw,info)
         CALL CHKXER('CPBSVX',INFot,NOUt,LERr,OK)
!
      ELSEIF ( LSAMEN(2,c2,'PT') ) THEN
!
!        CPTSV
!
         SRNamt = 'CPTSV '
         INFot = 1
         CALL CPTSV(-1,0,r,a(1,1),b,1,info)
         CALL CHKXER('CPTSV ',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CPTSV(0,-1,r,a(1,1),b,1,info)
         CALL CHKXER('CPTSV ',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL CPTSV(2,0,r,a(1,1),b,1,info)
         CALL CHKXER('CPTSV ',INFot,NOUt,LERr,OK)
!
!        CPTSVX
!
         SRNamt = 'CPTSVX'
         INFot = 1
         CALL CPTSVX('/',0,0,r,a(1,1),rf,af(1,1),b,1,x,1,rcond,r1,r2,w, &
     &               rw,info)
         CALL CHKXER('CPTSVX',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CPTSVX('N',-1,0,r,a(1,1),rf,af(1,1),b,1,x,1,rcond,r1,r2,w,&
     &               rw,info)
         CALL CHKXER('CPTSVX',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CPTSVX('N',0,-1,r,a(1,1),rf,af(1,1),b,1,x,1,rcond,r1,r2,w,&
     &               rw,info)
         CALL CHKXER('CPTSVX',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL CPTSVX('N',2,0,r,a(1,1),rf,af(1,1),b,1,x,2,rcond,r1,r2,w, &
     &               rw,info)
         CALL CHKXER('CPTSVX',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL CPTSVX('N',2,0,r,a(1,1),rf,af(1,1),b,2,x,1,rcond,r1,r2,w, &
     &               rw,info)
         CALL CHKXER('CPTSVX',INFot,NOUt,LERr,OK)
!
      ELSEIF ( LSAMEN(2,c2,'HE') ) THEN
!
!        CHESV
!
         SRNamt = 'CHESV '
         INFot = 1
         CALL CHESV('/',0,0,a,1,ip,b,1,w,1,info)
         CALL CHKXER('CHESV ',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CHESV('U',-1,0,a,1,ip,b,1,w,1,info)
         CALL CHKXER('CHESV ',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CHESV('U',0,-1,a,1,ip,b,1,w,1,info)
         CALL CHKXER('CHESV ',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL CHESV('U',2,0,a,1,ip,b,2,w,1,info)
         CALL CHKXER('CHESV ',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL CHESV('U',2,0,a,2,ip,b,1,w,1,info)
         CALL CHKXER('CHESV ',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL CHESV('U',0,0,a,1,ip,b,1,w,0,info)
         CALL CHKXER('CHESV ',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL CHESV('U',0,0,a,1,ip,b,1,w,-2,info)
         CALL CHKXER('CHESV ',INFot,NOUt,LERr,OK)
!
!        CHESVX
!
         SRNamt = 'CHESVX'
         INFot = 1
         CALL CHESVX('/','U',0,0,a,1,af,1,ip,b,1,x,1,rcond,r1,r2,w,1,rw,&
     &               info)
         CALL CHKXER('CHESVX',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CHESVX('N','/',0,0,a,1,af,1,ip,b,1,x,1,rcond,r1,r2,w,1,rw,&
     &               info)
         CALL CHKXER('CHESVX',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CHESVX('N','U',-1,0,a,1,af,1,ip,b,1,x,1,rcond,r1,r2,w,1,  &
     &               rw,info)
         CALL CHKXER('CHESVX',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CHESVX('N','U',0,-1,a,1,af,1,ip,b,1,x,1,rcond,r1,r2,w,1,  &
     &               rw,info)
         CALL CHKXER('CHESVX',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL CHESVX('N','U',2,0,a,1,af,2,ip,b,2,x,2,rcond,r1,r2,w,4,rw,&
     &               info)
         CALL CHKXER('CHESVX',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL CHESVX('N','U',2,0,a,2,af,1,ip,b,2,x,2,rcond,r1,r2,w,4,rw,&
     &               info)
         CALL CHKXER('CHESVX',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL CHESVX('N','U',2,0,a,2,af,2,ip,b,1,x,2,rcond,r1,r2,w,4,rw,&
     &               info)
         CALL CHKXER('CHESVX',INFot,NOUt,LERr,OK)
         INFot = 13
         CALL CHESVX('N','U',2,0,a,2,af,2,ip,b,2,x,1,rcond,r1,r2,w,4,rw,&
     &               info)
         CALL CHKXER('CHESVX',INFot,NOUt,LERr,OK)
         INFot = 18
         CALL CHESVX('N','U',2,0,a,2,af,2,ip,b,2,x,2,rcond,r1,r2,w,3,rw,&
     &               info)
         CALL CHKXER('CHESVX',INFot,NOUt,LERr,OK)
!
!        CHESVXX
!
         n_err_bnds = 3
         nparams = 1
         SRNamt = 'CHESVXX'
         INFot = 1
         CALL CHESVXX('/','U',0,0,a,1,af,1,ip,eq,c,b,1,x,1,rcond,rpvgrw,&
     &                berr,n_err_bnds,err_bnds_n,err_bnds_c,nparams,    &
     &                params,w,rw,info)
         CALL CHKXER('CHESVXX',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CHESVXX('N','/',0,0,a,1,af,1,ip,eq,c,b,1,x,1,rcond,rpvgrw,&
     &                berr,n_err_bnds,err_bnds_n,err_bnds_c,nparams,    &
     &                params,w,rw,info)
         CALL CHKXER('CHESVXX',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CHESVXX('N','U',-1,0,a,1,af,1,ip,eq,c,b,1,x,1,rcond,      &
     &                rpvgrw,berr,n_err_bnds,err_bnds_n,err_bnds_c,     &
     &                nparams,params,w,rw,info)
         CALL CHKXER('CHESVXX',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CHESVXX('N','U',0,-1,a,1,af,1,ip,eq,c,b,1,x,1,rcond,      &
     &                rpvgrw,berr,n_err_bnds,err_bnds_n,err_bnds_c,     &
     &                nparams,params,w,rw,info)
         CALL CHKXER('CHESVXX',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL CHESVXX('N','U',2,0,a,1,af,2,ip,eq,c,b,2,x,2,rcond,rpvgrw,&
     &                berr,n_err_bnds,err_bnds_n,err_bnds_c,nparams,    &
     &                params,w,rw,info)
         CALL CHKXER('CHESVXX',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL CHESVXX('N','U',2,0,a,2,af,1,ip,eq,c,b,2,x,2,rcond,rpvgrw,&
     &                berr,n_err_bnds,err_bnds_n,err_bnds_c,nparams,    &
     &                params,w,rw,info)
         CALL CHKXER('CHESVXX',INFot,NOUt,LERr,OK)
         INFot = 9
         eq = '/'
         CALL CHESVXX('F','U',0,0,a,1,af,1,ip,eq,c,b,1,x,1,rcond,rpvgrw,&
     &                berr,n_err_bnds,err_bnds_n,err_bnds_c,nparams,    &
     &                params,w,rw,info)
         CALL CHKXER('CHESVXX',INFot,NOUt,LERr,OK)
         INFot = 10
         eq = 'Y'
         CALL CHESVXX('F','U',1,0,a,1,af,1,ip,eq,c,b,1,x,1,rcond,rpvgrw,&
     &                berr,n_err_bnds,err_bnds_n,err_bnds_c,nparams,    &
     &                params,w,rw,info)
         CALL CHKXER('CHESVXX',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL CHESVXX('N','U',2,0,a,2,af,2,ip,eq,c,b,1,x,2,rcond,rpvgrw,&
     &                berr,n_err_bnds,err_bnds_n,err_bnds_c,nparams,    &
     &                params,w,rw,info)
         CALL CHKXER('CHESVXX',INFot,NOUt,LERr,OK)
         INFot = 14
         CALL CHESVXX('N','U',2,0,a,2,af,2,ip,eq,c,b,2,x,1,rcond,rpvgrw,&
     &                berr,n_err_bnds,err_bnds_n,err_bnds_c,nparams,    &
     &                params,w,rw,info)
         CALL CHKXER('CHESVXX',INFot,NOUt,LERr,OK)
!
      ELSEIF ( LSAMEN(2,c2,'HR') ) THEN
!
!        CHESV_ROOK
!
         SRNamt = 'CHESV_ROOK'
         INFot = 1
         CALL CHESV_ROOK('/',0,0,a,1,ip,b,1,w,1,info)
         CALL CHKXER('CHESV_ROOK',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CHESV_ROOK('U',-1,0,a,1,ip,b,1,w,1,info)
         CALL CHKXER('CHESV_ROOK',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CHESV_ROOK('U',0,-1,a,1,ip,b,1,w,1,info)
         CALL CHKXER('CHESV_ROOK',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL CHESV_ROOK('U',2,0,a,2,ip,b,1,w,1,info)
         CALL CHKXER('CHESV_ROOK',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL CHESV_ROOK('U',0,0,a,1,ip,b,1,w,0,info)
         CALL CHKXER('CHESV_ROOK',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL CHESV_ROOK('U',0,0,a,1,ip,b,1,w,-2,info)
         CALL CHKXER('CHESV_ROOK',INFot,NOUt,LERr,OK)
!
      ELSEIF ( LSAMEN(2,c2,'HK') ) THEN
!
!        CHESV_RK
!
!        Test error exits of the driver that uses factorization
!        of a symmetric indefinite matrix with rook
!        (bounded Bunch-Kaufman) pivoting with the new storage
!        format for factors L ( or U) and D.
!
!        L (or U) is stored in A, diagonal of D is stored on the
!        diagonal of A, subdiagonal of D is stored in a separate array E.
!
         SRNamt = 'CHESV_RK'
         INFot = 1
         CALL CHESV_RK('/',0,0,a,1,e,ip,b,1,w,1,info)
         CALL CHKXER('CHESV_RK',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CHESV_RK('U',-1,0,a,1,e,ip,b,1,w,1,info)
         CALL CHKXER('CHESV_RK',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CHESV_RK('U',0,-1,a,1,e,ip,b,1,w,1,info)
         CALL CHKXER('CHESV_RK',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL CHESV_RK('U',2,0,a,1,e,ip,b,2,w,1,info)
         CALL CHKXER('CHESV_RK',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL CHESV_RK('U',2,0,a,2,e,ip,b,1,w,1,info)
         CALL CHKXER('CHESV_RK',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL CHESV_RK('U',0,0,a,1,e,ip,b,1,w,0,info)
         CALL CHKXER('CHESV_RK',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL CHESV_RK('U',0,0,a,1,e,ip,b,1,w,-2,info)
         CALL CHKXER('CHESV_RK',INFot,NOUt,LERr,OK)
!
      ELSEIF ( LSAMEN(2,c2,'HP') ) THEN
!
!        CHPSV
!
         SRNamt = 'CHPSV '
         INFot = 1
         CALL CHPSV('/',0,0,a,ip,b,1,info)
         CALL CHKXER('CHPSV ',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CHPSV('U',-1,0,a,ip,b,1,info)
         CALL CHKXER('CHPSV ',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CHPSV('U',0,-1,a,ip,b,1,info)
         CALL CHKXER('CHPSV ',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL CHPSV('U',2,0,a,ip,b,1,info)
         CALL CHKXER('CHPSV ',INFot,NOUt,LERr,OK)
!
!        CHPSVX
!
         SRNamt = 'CHPSVX'
         INFot = 1
         CALL CHPSVX('/','U',0,0,a,af,ip,b,1,x,1,rcond,r1,r2,w,rw,info)
         CALL CHKXER('CHPSVX',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CHPSVX('N','/',0,0,a,af,ip,b,1,x,1,rcond,r1,r2,w,rw,info)
         CALL CHKXER('CHPSVX',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CHPSVX('N','U',-1,0,a,af,ip,b,1,x,1,rcond,r1,r2,w,rw,info)
         CALL CHKXER('CHPSVX',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CHPSVX('N','U',0,-1,a,af,ip,b,1,x,1,rcond,r1,r2,w,rw,info)
         CALL CHKXER('CHPSVX',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL CHPSVX('N','U',2,0,a,af,ip,b,1,x,2,rcond,r1,r2,w,rw,info)
         CALL CHKXER('CHPSVX',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL CHPSVX('N','U',2,0,a,af,ip,b,2,x,1,rcond,r1,r2,w,rw,info)
         CALL CHKXER('CHPSVX',INFot,NOUt,LERr,OK)
!
      ELSEIF ( LSAMEN(2,c2,'SY') ) THEN
!
!        CSYSV
!
         SRNamt = 'CSYSV '
         INFot = 1
         CALL CSYSV('/',0,0,a,1,ip,b,1,w,1,info)
         CALL CHKXER('CSYSV ',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CSYSV('U',-1,0,a,1,ip,b,1,w,1,info)
         CALL CHKXER('CSYSV ',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CSYSV('U',0,-1,a,1,ip,b,1,w,1,info)
         CALL CHKXER('CSYSV ',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL CSYSV('U',2,0,a,2,ip,b,1,w,1,info)
         CALL CHKXER('CSYSV ',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL CSYSV('U',0,0,a,1,ip,b,1,w,0,info)
         CALL CHKXER('CSYSV ',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL CSYSV('U',0,0,a,1,ip,b,1,w,-2,info)
         CALL CHKXER('CSYSV ',INFot,NOUt,LERr,OK)
!
!        CSYSVX
!
         SRNamt = 'CSYSVX'
         INFot = 1
         CALL CSYSVX('/','U',0,0,a,1,af,1,ip,b,1,x,1,rcond,r1,r2,w,1,rw,&
     &               info)
         CALL CHKXER('CSYSVX',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CSYSVX('N','/',0,0,a,1,af,1,ip,b,1,x,1,rcond,r1,r2,w,1,rw,&
     &               info)
         CALL CHKXER('CSYSVX',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CSYSVX('N','U',-1,0,a,1,af,1,ip,b,1,x,1,rcond,r1,r2,w,1,  &
     &               rw,info)
         CALL CHKXER('CSYSVX',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CSYSVX('N','U',0,-1,a,1,af,1,ip,b,1,x,1,rcond,r1,r2,w,1,  &
     &               rw,info)
         CALL CHKXER('CSYSVX',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL CSYSVX('N','U',2,0,a,1,af,2,ip,b,2,x,2,rcond,r1,r2,w,4,rw,&
     &               info)
         CALL CHKXER('CSYSVX',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL CSYSVX('N','U',2,0,a,2,af,1,ip,b,2,x,2,rcond,r1,r2,w,4,rw,&
     &               info)
         CALL CHKXER('CSYSVX',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL CSYSVX('N','U',2,0,a,2,af,2,ip,b,1,x,2,rcond,r1,r2,w,4,rw,&
     &               info)
         CALL CHKXER('CSYSVX',INFot,NOUt,LERr,OK)
         INFot = 13
         CALL CSYSVX('N','U',2,0,a,2,af,2,ip,b,2,x,1,rcond,r1,r2,w,4,rw,&
     &               info)
         CALL CHKXER('CSYSVX',INFot,NOUt,LERr,OK)
         INFot = 18
         CALL CSYSVX('N','U',2,0,a,2,af,2,ip,b,2,x,2,rcond,r1,r2,w,3,rw,&
     &               info)
         CALL CHKXER('CSYSVX',INFot,NOUt,LERr,OK)
!
!        CSYSVXX
!
         n_err_bnds = 3
         nparams = 1
         SRNamt = 'CSYSVXX'
         INFot = 1
         eq = 'N'
         CALL CSYSVXX('/','U',0,0,a,1,af,1,ip,eq,r,b,1,x,1,rcond,rpvgrw,&
     &                berr,n_err_bnds,err_bnds_n,err_bnds_c,nparams,    &
     &                params,w,rw,info)
         CALL CHKXER('CSYSVXX',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CSYSVXX('N','/',0,0,a,1,af,1,ip,eq,r,b,1,x,1,rcond,rpvgrw,&
     &                berr,n_err_bnds,err_bnds_n,err_bnds_c,nparams,    &
     &                params,w,rw,info)
         CALL CHKXER('CSYSVXX',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CSYSVXX('N','U',-1,0,a,1,af,1,ip,eq,r,b,1,x,1,rcond,      &
     &                rpvgrw,berr,n_err_bnds,err_bnds_n,err_bnds_c,     &
     &                nparams,params,w,rw,info)
         CALL CHKXER('CSYSVXX',INFot,NOUt,LERr,OK)
         INFot = 4
         eq = '/'
         CALL CSYSVXX('N','U',0,-1,a,1,af,1,ip,eq,r,b,1,x,1,rcond,      &
     &                rpvgrw,berr,n_err_bnds,err_bnds_n,err_bnds_c,     &
     &                nparams,params,w,rw,info)
         CALL CHKXER('CSYSVXX',INFot,NOUt,LERr,OK)
         eq = 'Y'
         INFot = 6
         CALL CSYSVXX('N','U',2,0,a,1,af,2,ip,eq,r,b,2,x,2,rcond,rpvgrw,&
     &                berr,n_err_bnds,err_bnds_n,err_bnds_c,nparams,    &
     &                params,w,rw,info)
         CALL CHKXER('CSYSVXX',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL CSYSVXX('N','U',2,0,a,2,af,1,ip,eq,r,b,2,x,2,rcond,rpvgrw,&
     &                berr,n_err_bnds,err_bnds_n,err_bnds_c,nparams,    &
     &                params,w,rw,info)
         CALL CHKXER('CSYSVXX',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL CSYSVXX('F','U',2,0,a,2,af,2,ip,'A',r,b,2,x,2,rcond,      &
     &                rpvgrw,berr,n_err_bnds,err_bnds_n,err_bnds_c,     &
     &                nparams,params,w,rw,info)
         CALL CHKXER('CSYSVXX',INFot,NOUt,LERr,OK)
         INFot = 11
         eq = 'Y'
         CALL CSYSVXX('F','U',2,0,a,2,af,2,ip,eq,r,b,2,x,2,rcond,rpvgrw,&
     &                berr,n_err_bnds,err_bnds_n,err_bnds_c,nparams,    &
     &                params,w,rw,info)
         CALL CHKXER('CSYSVXX',INFot,NOUt,LERr,OK)
         INFot = 11
         eq = 'Y'
         r(1) = -ONE
         CALL CSYSVXX('F','U',2,0,a,2,af,2,ip,eq,r,b,2,x,2,rcond,rpvgrw,&
     &                berr,n_err_bnds,err_bnds_n,err_bnds_c,nparams,    &
     &                params,w,rw,info)
         CALL CHKXER('CSYSVXX',INFot,NOUt,LERr,OK)
         INFot = 13
         eq = 'N'
         CALL CSYSVXX('N','U',2,0,a,2,af,2,ip,eq,r,b,1,x,2,rcond,rpvgrw,&
     &                berr,n_err_bnds,err_bnds_n,err_bnds_c,nparams,    &
     &                params,w,rw,info)
         CALL CHKXER('CSYSVXX',INFot,NOUt,LERr,OK)
         INFot = 15
         CALL CSYSVXX('N','U',2,0,a,2,af,2,ip,eq,r,b,2,x,1,rcond,rpvgrw,&
     &                berr,n_err_bnds,err_bnds_n,err_bnds_c,nparams,    &
     &                params,w,rw,info)
         CALL CHKXER('CSYSVXX',INFot,NOUt,LERr,OK)
!
      ELSEIF ( LSAMEN(2,c2,'SR') ) THEN
!
!        CSYSV_ROOK
!
         SRNamt = 'CSYSV_ROOK'
         INFot = 1
         CALL CSYSV_ROOK('/',0,0,a,1,ip,b,1,w,1,info)
         CALL CHKXER('CSYSV_ROOK',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CSYSV_ROOK('U',-1,0,a,1,ip,b,1,w,1,info)
         CALL CHKXER('CSYSV_ROOK',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CSYSV_ROOK('U',0,-1,a,1,ip,b,1,w,1,info)
         CALL CHKXER('CSYSV_ROOK',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL CSYSV_ROOK('U',2,0,a,2,ip,b,1,w,1,info)
         CALL CHKXER('CSYSV_ROOK',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL CSYSV_ROOK('U',0,0,a,1,ip,b,1,w,0,info)
         CALL CHKXER('CSYSV_ROOK',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL CSYSV_ROOK('U',0,0,a,1,ip,b,1,w,-2,info)
         CALL CHKXER('CSYSV_ROOK',INFot,NOUt,LERr,OK)
!
      ELSEIF ( LSAMEN(2,c2,'SK') ) THEN
!
!        CSYSV_RK
!
!        Test error exits of the driver that uses factorization
!        of a symmetric indefinite matrix with rook
!        (bounded Bunch-Kaufman) pivoting with the new storage
!        format for factors L ( or U) and D.
!
!        L (or U) is stored in A, diagonal of D is stored on the
!        diagonal of A, subdiagonal of D is stored in a separate array E.
!
         SRNamt = 'CSYSV_RK'
         INFot = 1
         CALL CSYSV_RK('/',0,0,a,1,e,ip,b,1,w,1,info)
         CALL CHKXER('CSYSV_RK',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CSYSV_RK('U',-1,0,a,1,e,ip,b,1,w,1,info)
         CALL CHKXER('CSYSV_RK',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CSYSV_RK('U',0,-1,a,1,e,ip,b,1,w,1,info)
         CALL CHKXER('CSYSV_RK',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL CSYSV_RK('U',2,0,a,1,e,ip,b,2,w,1,info)
         CALL CHKXER('CSYSV_RK',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL CSYSV_RK('U',2,0,a,2,e,ip,b,1,w,1,info)
         CALL CHKXER('CSYSV_RK',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL CSYSV_RK('U',0,0,a,1,e,ip,b,1,w,0,info)
         CALL CHKXER('CSYSV_RK',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL CSYSV_RK('U',0,0,a,1,e,ip,b,1,w,-2,info)
         CALL CHKXER('CSYSV_RK',INFot,NOUt,LERr,OK)
!
      ELSEIF ( LSAMEN(2,c2,'SP') ) THEN
!
!        CSPSV
!
         SRNamt = 'CSPSV '
         INFot = 1
         CALL CSPSV('/',0,0,a,ip,b,1,info)
         CALL CHKXER('CSPSV ',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CSPSV('U',-1,0,a,ip,b,1,info)
         CALL CHKXER('CSPSV ',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CSPSV('U',0,-1,a,ip,b,1,info)
         CALL CHKXER('CSPSV ',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL CSPSV('U',2,0,a,ip,b,1,info)
         CALL CHKXER('CSPSV ',INFot,NOUt,LERr,OK)
!
!        CSPSVX
!
         SRNamt = 'CSPSVX'
         INFot = 1
         CALL CSPSVX('/','U',0,0,a,af,ip,b,1,x,1,rcond,r1,r2,w,rw,info)
         CALL CHKXER('CSPSVX',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CSPSVX('N','/',0,0,a,af,ip,b,1,x,1,rcond,r1,r2,w,rw,info)
         CALL CHKXER('CSPSVX',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CSPSVX('N','U',-1,0,a,af,ip,b,1,x,1,rcond,r1,r2,w,rw,info)
         CALL CHKXER('CSPSVX',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CSPSVX('N','U',0,-1,a,af,ip,b,1,x,1,rcond,r1,r2,w,rw,info)
         CALL CHKXER('CSPSVX',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL CSPSVX('N','U',2,0,a,af,ip,b,1,x,2,rcond,r1,r2,w,rw,info)
         CALL CHKXER('CSPSVX',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL CSPSVX('N','U',2,0,a,af,ip,b,2,x,1,rcond,r1,r2,w,rw,info)
         CALL CHKXER('CSPSVX',INFot,NOUt,LERr,OK)
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
!     End of cerrvx
!
      END SUBROUTINE cerrvx
