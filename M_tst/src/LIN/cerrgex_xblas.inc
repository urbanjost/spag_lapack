!*==cerrge.f90  processed by SPAG 7.51RB at 17:34 on  4 Mar 2022
!> \brief \b CERRGEX
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE cerrge( PATH, NUNIT )
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
!> cerrge tests the error exits for the COMPLEX routines
!> for general matrices.
!>
!> Note that this file is used only when the XBLAS are available,
!> otherwise cerrge.f defines this subroutine.
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
      SUBROUTINE cerrge(Path,Nunit)
      IMPLICIT NONE
!*--CERRGE62
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
!     ..
!     .. Local Scalars ..
      CHARACTER eq
      CHARACTER*2 c2
      INTEGER i , info , j , n_err_bnds , nparams
      REAL anrm , ccond , rcond , berr
!     ..
!     .. Local Arrays ..
      INTEGER ip(NMAX)
      REAL r(NMAX) , r1(NMAX) , r2(NMAX) , cs(NMAX) , rs(NMAX)
      COMPLEX a(NMAX,NMAX) , af(NMAX,NMAX) , b(NMAX) , w(2*NMAX) ,      &
     &        x(NMAX) , err_bnds_n(NMAX,3) , err_bnds_c(NMAX,3) ,       &
     &        params(1)
!     ..
!     .. External Functions ..
      LOGICAL LSAMEN
EXTERNAL    LSAMEN
!     ..
!     .. External Subroutines ..
EXTERNAL    CGBCON,   CGBEQU,   CGBRFS,   CGBTF2,   CGBTRF
EXTERNAL    CGBTRS,   CGECON,   CGEEQU,   CGERFS,   CGETF2,   CGETRF
EXTERNAL    CGETRI,   CGETRS,   CGEEQUB,   CGERFSX,   CGBEQUB
EXTERNAL    CGBRFSX
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
         b(j) = 0.
         r1(j) = 0.
         r2(j) = 0.
         w(j) = 0.
         x(j) = 0.
         cs(j) = 0.
         rs(j) = 0.
         ip(j) = j
      ENDDO
      OK = .TRUE.
!
!     Test error exits of the routines that use the LU decomposition
!     of a general matrix.
!
      IF ( LSAMEN(2,c2,'GE') ) THEN
!
!        CGETRF
!
         SRNamt = 'CGETRF'
         INFot = 1
         CALL CGETRF(-1,0,a,1,ip,info)
         CALL CHKXER('CGETRF',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CGETRF(0,-1,a,1,ip,info)
         CALL CHKXER('CGETRF',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CGETRF(2,1,a,1,ip,info)
         CALL CHKXER('CGETRF',INFot,NOUt,LERr,OK)
!
!        CGETF2
!
         SRNamt = 'CGETF2'
         INFot = 1
         CALL CGETF2(-1,0,a,1,ip,info)
         CALL CHKXER('CGETF2',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CGETF2(0,-1,a,1,ip,info)
         CALL CHKXER('CGETF2',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CGETF2(2,1,a,1,ip,info)
         CALL CHKXER('CGETF2',INFot,NOUt,LERr,OK)
!
!        CGETRI
!
         SRNamt = 'CGETRI'
         INFot = 1
         CALL CGETRI(-1,a,1,ip,w,1,info)
         CALL CHKXER('CGETRI',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CGETRI(2,a,1,ip,w,2,info)
         CALL CHKXER('CGETRI',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL CGETRI(2,a,2,ip,w,1,info)
         CALL CHKXER('CGETRI',INFot,NOUt,LERr,OK)
!
!        CGETRS
!
         SRNamt = 'CGETRS'
         INFot = 1
         CALL CGETRS('/',0,0,a,1,ip,b,1,info)
         CALL CHKXER('CGETRS',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CGETRS('N',-1,0,a,1,ip,b,1,info)
         CALL CHKXER('CGETRS',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CGETRS('N',0,-1,a,1,ip,b,1,info)
         CALL CHKXER('CGETRS',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL CGETRS('N',2,1,a,1,ip,b,2,info)
         CALL CHKXER('CGETRS',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL CGETRS('N',2,1,a,2,ip,b,1,info)
         CALL CHKXER('CGETRS',INFot,NOUt,LERr,OK)
!
!        CGERFS
!
         SRNamt = 'CGERFS'
         INFot = 1
         CALL CGERFS('/',0,0,a,1,af,1,ip,b,1,x,1,r1,r2,w,r,info)
         CALL CHKXER('CGERFS',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CGERFS('N',-1,0,a,1,af,1,ip,b,1,x,1,r1,r2,w,r,info)
         CALL CHKXER('CGERFS',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CGERFS('N',0,-1,a,1,af,1,ip,b,1,x,1,r1,r2,w,r,info)
         CALL CHKXER('CGERFS',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL CGERFS('N',2,1,a,1,af,2,ip,b,2,x,2,r1,r2,w,r,info)
         CALL CHKXER('CGERFS',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL CGERFS('N',2,1,a,2,af,1,ip,b,2,x,2,r1,r2,w,r,info)
         CALL CHKXER('CGERFS',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL CGERFS('N',2,1,a,2,af,2,ip,b,1,x,2,r1,r2,w,r,info)
         CALL CHKXER('CGERFS',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL CGERFS('N',2,1,a,2,af,2,ip,b,2,x,1,r1,r2,w,r,info)
         CALL CHKXER('CGERFS',INFot,NOUt,LERr,OK)
!
!        CGERFSX
!
         n_err_bnds = 3
         nparams = 0
         SRNamt = 'CGERFSX'
         INFot = 1
         CALL CGERFSX('/',eq,0,0,a,1,af,1,ip,rs,cs,b,1,x,1,rcond,berr,  &
     &                n_err_bnds,err_bnds_n,err_bnds_c,nparams,params,w,&
     &                r,info)
         CALL CHKXER('CGERFSX',INFot,NOUt,LERr,OK)
         INFot = 2
         eq = '/'
         CALL CGERFSX('N',eq,2,1,a,1,af,2,ip,rs,cs,b,2,x,2,rcond,berr,  &
     &                n_err_bnds,err_bnds_n,err_bnds_c,nparams,params,w,&
     &                r,info)
         CALL CHKXER('CGERFSX',INFot,NOUt,LERr,OK)
         INFot = 3
         eq = 'R'
         CALL CGERFSX('N',eq,-1,0,a,1,af,1,ip,rs,cs,b,1,x,1,rcond,berr, &
     &                n_err_bnds,err_bnds_n,err_bnds_c,nparams,params,w,&
     &                r,info)
         CALL CHKXER('CGERFSX',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CGERFSX('N',eq,0,-1,a,1,af,1,ip,rs,cs,b,1,x,1,rcond,berr, &
     &                n_err_bnds,err_bnds_n,err_bnds_c,nparams,params,w,&
     &                r,info)
         CALL CHKXER('CGERFSX',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL CGERFSX('N',eq,2,1,a,1,af,2,ip,rs,cs,b,2,x,2,rcond,berr,  &
     &                n_err_bnds,err_bnds_n,err_bnds_c,nparams,params,w,&
     &                r,info)
         CALL CHKXER('CGERFSX',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL CGERFSX('N',eq,2,1,a,2,af,1,ip,rs,cs,b,2,x,2,rcond,berr,  &
     &                n_err_bnds,err_bnds_n,err_bnds_c,nparams,params,w,&
     &                r,info)
         CALL CHKXER('CGERFSX',INFot,NOUt,LERr,OK)
         INFot = 13
         eq = 'C'
         CALL CGERFSX('N',eq,2,1,a,2,af,2,ip,rs,cs,b,1,x,2,rcond,berr,  &
     &                n_err_bnds,err_bnds_n,err_bnds_c,nparams,params,w,&
     &                r,info)
         CALL CHKXER('CGERFSX',INFot,NOUt,LERr,OK)
         INFot = 15
         CALL CGERFSX('N',eq,2,1,a,2,af,2,ip,rs,cs,b,2,x,1,rcond,berr,  &
     &                n_err_bnds,err_bnds_n,err_bnds_c,nparams,params,w,&
     &                r,info)
         CALL CHKXER('CGERFSX',INFot,NOUt,LERr,OK)
!
!        CGECON
!
         SRNamt = 'CGECON'
         INFot = 1
         CALL CGECON('/',0,a,1,anrm,rcond,w,r,info)
         CALL CHKXER('CGECON',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CGECON('1',-1,a,1,anrm,rcond,w,r,info)
         CALL CHKXER('CGECON',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CGECON('1',2,a,1,anrm,rcond,w,r,info)
         CALL CHKXER('CGECON',INFot,NOUt,LERr,OK)
!
!        CGEEQU
!
         SRNamt = 'CGEEQU'
         INFot = 1
         CALL CGEEQU(-1,0,a,1,r1,r2,rcond,ccond,anrm,info)
         CALL CHKXER('CGEEQU',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CGEEQU(0,-1,a,1,r1,r2,rcond,ccond,anrm,info)
         CALL CHKXER('CGEEQU',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CGEEQU(2,2,a,1,r1,r2,rcond,ccond,anrm,info)
         CALL CHKXER('CGEEQU',INFot,NOUt,LERr,OK)
!
!        CGEEQUB
!
         SRNamt = 'CGEEQUB'
         INFot = 1
         CALL CGEEQUB(-1,0,a,1,r1,r2,rcond,ccond,anrm,info)
         CALL CHKXER('CGEEQUB',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CGEEQUB(0,-1,a,1,r1,r2,rcond,ccond,anrm,info)
         CALL CHKXER('CGEEQUB',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CGEEQUB(2,2,a,1,r1,r2,rcond,ccond,anrm,info)
         CALL CHKXER('CGEEQUB',INFot,NOUt,LERr,OK)
!
!     Test error exits of the routines that use the LU decomposition
!     of a general band matrix.
!
      ELSEIF ( LSAMEN(2,c2,'GB') ) THEN
!
!        CGBTRF
!
         SRNamt = 'CGBTRF'
         INFot = 1
         CALL CGBTRF(-1,0,0,0,a,1,ip,info)
         CALL CHKXER('CGBTRF',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CGBTRF(0,-1,0,0,a,1,ip,info)
         CALL CHKXER('CGBTRF',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CGBTRF(1,1,-1,0,a,1,ip,info)
         CALL CHKXER('CGBTRF',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CGBTRF(1,1,0,-1,a,1,ip,info)
         CALL CHKXER('CGBTRF',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL CGBTRF(2,2,1,1,a,3,ip,info)
         CALL CHKXER('CGBTRF',INFot,NOUt,LERr,OK)
!
!        CGBTF2
!
         SRNamt = 'CGBTF2'
         INFot = 1
         CALL CGBTF2(-1,0,0,0,a,1,ip,info)
         CALL CHKXER('CGBTF2',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CGBTF2(0,-1,0,0,a,1,ip,info)
         CALL CHKXER('CGBTF2',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CGBTF2(1,1,-1,0,a,1,ip,info)
         CALL CHKXER('CGBTF2',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CGBTF2(1,1,0,-1,a,1,ip,info)
         CALL CHKXER('CGBTF2',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL CGBTF2(2,2,1,1,a,3,ip,info)
         CALL CHKXER('CGBTF2',INFot,NOUt,LERr,OK)
!
!        CGBTRS
!
         SRNamt = 'CGBTRS'
         INFot = 1
         CALL CGBTRS('/',0,0,0,1,a,1,ip,b,1,info)
         CALL CHKXER('CGBTRS',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CGBTRS('N',-1,0,0,1,a,1,ip,b,1,info)
         CALL CHKXER('CGBTRS',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CGBTRS('N',1,-1,0,1,a,1,ip,b,1,info)
         CALL CHKXER('CGBTRS',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CGBTRS('N',1,0,-1,1,a,1,ip,b,1,info)
         CALL CHKXER('CGBTRS',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL CGBTRS('N',1,0,0,-1,a,1,ip,b,1,info)
         CALL CHKXER('CGBTRS',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL CGBTRS('N',2,1,1,1,a,3,ip,b,2,info)
         CALL CHKXER('CGBTRS',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL CGBTRS('N',2,0,0,1,a,1,ip,b,1,info)
         CALL CHKXER('CGBTRS',INFot,NOUt,LERr,OK)
!
!        CGBRFS
!
         SRNamt = 'CGBRFS'
         INFot = 1
         CALL CGBRFS('/',0,0,0,0,a,1,af,1,ip,b,1,x,1,r1,r2,w,r,info)
         CALL CHKXER('CGBRFS',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CGBRFS('N',-1,0,0,0,a,1,af,1,ip,b,1,x,1,r1,r2,w,r,info)
         CALL CHKXER('CGBRFS',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CGBRFS('N',1,-1,0,0,a,1,af,1,ip,b,1,x,1,r1,r2,w,r,info)
         CALL CHKXER('CGBRFS',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CGBRFS('N',1,0,-1,0,a,1,af,1,ip,b,1,x,1,r1,r2,w,r,info)
         CALL CHKXER('CGBRFS',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL CGBRFS('N',1,0,0,-1,a,1,af,1,ip,b,1,x,1,r1,r2,w,r,info)
         CALL CHKXER('CGBRFS',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL CGBRFS('N',2,1,1,1,a,2,af,4,ip,b,2,x,2,r1,r2,w,r,info)
         CALL CHKXER('CGBRFS',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL CGBRFS('N',2,1,1,1,a,3,af,3,ip,b,2,x,2,r1,r2,w,r,info)
         CALL CHKXER('CGBRFS',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL CGBRFS('N',2,0,0,1,a,1,af,1,ip,b,1,x,2,r1,r2,w,r,info)
         CALL CHKXER('CGBRFS',INFot,NOUt,LERr,OK)
         INFot = 14
         CALL CGBRFS('N',2,0,0,1,a,1,af,1,ip,b,2,x,1,r1,r2,w,r,info)
         CALL CHKXER('CGBRFS',INFot,NOUt,LERr,OK)
!
!        CGBRFSX
!
         n_err_bnds = 3
         nparams = 0
         SRNamt = 'CGBRFSX'
         INFot = 1
         CALL CGBRFSX('/',eq,0,0,0,0,a,1,af,1,ip,rs,cs,b,1,x,1,rcond,   &
     &                berr,n_err_bnds,err_bnds_n,err_bnds_c,nparams,    &
     &                params,w,r,info)
         CALL CHKXER('CGBRFSX',INFot,NOUt,LERr,OK)
         INFot = 2
         eq = '/'
         CALL CGBRFSX('N',eq,2,1,1,1,a,1,af,2,ip,rs,cs,b,2,x,2,rcond,   &
     &                berr,n_err_bnds,err_bnds_n,err_bnds_c,nparams,    &
     &                params,w,r,info)
         CALL CHKXER('CGBRFSX',INFot,NOUt,LERr,OK)
         INFot = 3
         eq = 'R'
         CALL CGBRFSX('N',eq,-1,1,1,0,a,1,af,1,ip,rs,cs,b,1,x,1,rcond,  &
     &                berr,n_err_bnds,err_bnds_n,err_bnds_c,nparams,    &
     &                params,w,r,info)
         CALL CHKXER('CGBRFSX',INFot,NOUt,LERr,OK)
         INFot = 4
         eq = 'R'
         CALL CGBRFSX('N',eq,2,-1,1,1,a,3,af,4,ip,rs,cs,b,1,x,1,rcond,  &
     &                berr,n_err_bnds,err_bnds_n,err_bnds_c,nparams,    &
     &                params,w,r,info)
         CALL CHKXER('CGBRFSX',INFot,NOUt,LERr,OK)
         INFot = 5
         eq = 'R'
         CALL CGBRFSX('N',eq,2,1,-1,1,a,3,af,4,ip,rs,cs,b,1,x,1,rcond,  &
     &                berr,n_err_bnds,err_bnds_n,err_bnds_c,nparams,    &
     &                params,w,r,info)
         CALL CHKXER('CGBRFSX',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL CGBRFSX('N',eq,0,0,0,-1,a,1,af,1,ip,rs,cs,b,1,x,1,rcond,  &
     &                berr,n_err_bnds,err_bnds_n,err_bnds_c,nparams,    &
     &                params,w,r,info)
         CALL CHKXER('CGBRFSX',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL CGBRFSX('N',eq,2,1,1,1,a,1,af,2,ip,rs,cs,b,2,x,2,rcond,   &
     &                berr,n_err_bnds,err_bnds_n,err_bnds_c,nparams,    &
     &                params,w,r,info)
         CALL CHKXER('CGBRFSX',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL CGBRFSX('N',eq,2,1,1,1,a,3,af,3,ip,rs,cs,b,2,x,2,rcond,   &
     &                berr,n_err_bnds,err_bnds_n,err_bnds_c,nparams,    &
     &                params,w,r,info)
         CALL CHKXER('CGBRFSX',INFot,NOUt,LERr,OK)
         INFot = 13
         eq = 'C'
         CALL CGBRFSX('N',eq,2,1,1,1,a,3,af,5,ip,rs,cs,b,1,x,2,rcond,   &
     &                berr,n_err_bnds,err_bnds_n,err_bnds_c,nparams,    &
     &                params,w,r,info)
         CALL CHKXER('CGBRFSX',INFot,NOUt,LERr,OK)
         INFot = 15
         CALL CGBRFSX('N',eq,2,1,1,1,a,3,af,5,ip,rs,cs,b,2,x,1,rcond,   &
     &                berr,n_err_bnds,err_bnds_n,err_bnds_c,nparams,    &
     &                params,w,r,info)
         CALL CHKXER('CGBRFSX',INFot,NOUt,LERr,OK)
!
!        CGBCON
!
         SRNamt = 'CGBCON'
         INFot = 1
         CALL CGBCON('/',0,0,0,a,1,ip,anrm,rcond,w,r,info)
         CALL CHKXER('CGBCON',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CGBCON('1',-1,0,0,a,1,ip,anrm,rcond,w,r,info)
         CALL CHKXER('CGBCON',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CGBCON('1',1,-1,0,a,1,ip,anrm,rcond,w,r,info)
         CALL CHKXER('CGBCON',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CGBCON('1',1,0,-1,a,1,ip,anrm,rcond,w,r,info)
         CALL CHKXER('CGBCON',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL CGBCON('1',2,1,1,a,3,ip,anrm,rcond,w,r,info)
         CALL CHKXER('CGBCON',INFot,NOUt,LERr,OK)
!
!        CGBEQU
!
         SRNamt = 'CGBEQU'
         INFot = 1
         CALL CGBEQU(-1,0,0,0,a,1,r1,r2,rcond,ccond,anrm,info)
         CALL CHKXER('CGBEQU',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CGBEQU(0,-1,0,0,a,1,r1,r2,rcond,ccond,anrm,info)
         CALL CHKXER('CGBEQU',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CGBEQU(1,1,-1,0,a,1,r1,r2,rcond,ccond,anrm,info)
         CALL CHKXER('CGBEQU',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CGBEQU(1,1,0,-1,a,1,r1,r2,rcond,ccond,anrm,info)
         CALL CHKXER('CGBEQU',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL CGBEQU(2,2,1,1,a,2,r1,r2,rcond,ccond,anrm,info)
         CALL CHKXER('CGBEQU',INFot,NOUt,LERr,OK)
!
!        CGBEQUB
!
         SRNamt = 'CGBEQUB'
         INFot = 1
         CALL CGBEQUB(-1,0,0,0,a,1,r1,r2,rcond,ccond,anrm,info)
         CALL CHKXER('CGBEQUB',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CGBEQUB(0,-1,0,0,a,1,r1,r2,rcond,ccond,anrm,info)
         CALL CHKXER('CGBEQUB',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CGBEQUB(1,1,-1,0,a,1,r1,r2,rcond,ccond,anrm,info)
         CALL CHKXER('CGBEQUB',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CGBEQUB(1,1,0,-1,a,1,r1,r2,rcond,ccond,anrm,info)
         CALL CHKXER('CGBEQUB',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL CGBEQUB(2,2,1,1,a,2,r1,r2,rcond,ccond,anrm,info)
         CALL CHKXER('CGBEQUB',INFot,NOUt,LERr,OK)
      ENDIF
!
!     Print a summary line.
!
      CALL alaesm(Path,OK,NOUt)
!
!
!     End of cerrge
!
      END SUBROUTINE cerrge
