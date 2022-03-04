!*==serrge.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b SERRGEX
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE SERRGE( PATH, NUNIT )
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
!> SERRGE tests the error exits for the REAL routines
!> for general matrices.
!>
!> Note that this file is used only when the XBLAS are available,
!> otherwise serrge.f defines this subroutine.
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
      SUBROUTINE SERRGE(Path,Nunit)
      IMPLICIT NONE
!*--SERRGE62
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
      INTEGER NMAX , LW
      PARAMETER (NMAX=4,LW=3*NMAX)
!     ..
!     .. Local Scalars ..
      CHARACTER eq
      CHARACTER*2 c2
      INTEGER i , info , j , n_err_bnds , nparams
      REAL anrm , ccond , rcond , berr
!     ..
!     .. Local Arrays ..
      INTEGER ip(NMAX) , iw(NMAX)
      REAL a(NMAX,NMAX) , af(NMAX,NMAX) , b(NMAX) , c(NMAX) , r(NMAX) , &
     &     r1(NMAX) , r2(NMAX) , w(LW) , x(NMAX) , err_bnds_n(NMAX,3) , &
     &     err_bnds_c(NMAX,3) , params(1)
!     ..
!     .. External Functions ..
      LOGICAL LSAMEN
      EXTERNAL LSAMEN
!     ..
!     .. External Subroutines ..
      EXTERNAL ALAESM , CHKXER , SGBCON , SGBEQU , SGBRFS , SGBTF2 ,    &
     &         SGBTRF , SGBTRS , SGECON , SGEEQU , SGERFS , SGETF2 ,    &
     &         SGETRF , SGETRI , SGETRS , SGEEQUB , SGERFSX , SGBEQUB , &
     &         SGBRFSX
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
         b(j) = 0.
         r1(j) = 0.
         r2(j) = 0.
         w(j) = 0.
         x(j) = 0.
         c(j) = 0.
         r(j) = 0.
         ip(j) = j
         iw(j) = j
      ENDDO
      OK = .TRUE.
!
      IF ( LSAMEN(2,c2,'GE') ) THEN
!
!        Test error exits of the routines that use the LU decomposition
!        of a general matrix.
!
!        SGETRF
!
         SRNamt = 'SGETRF'
         INFot = 1
         CALL SGETRF(-1,0,a,1,ip,info)
         CALL CHKXER('SGETRF',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SGETRF(0,-1,a,1,ip,info)
         CALL CHKXER('SGETRF',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL SGETRF(2,1,a,1,ip,info)
         CALL CHKXER('SGETRF',INFot,NOUt,LERr,OK)
!
!        SGETF2
!
         SRNamt = 'SGETF2'
         INFot = 1
         CALL SGETF2(-1,0,a,1,ip,info)
         CALL CHKXER('SGETF2',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SGETF2(0,-1,a,1,ip,info)
         CALL CHKXER('SGETF2',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL SGETF2(2,1,a,1,ip,info)
         CALL CHKXER('SGETF2',INFot,NOUt,LERr,OK)
!
!        SGETRI
!
         SRNamt = 'SGETRI'
         INFot = 1
         CALL SGETRI(-1,a,1,ip,w,LW,info)
         CALL CHKXER('SGETRI',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SGETRI(2,a,1,ip,w,LW,info)
         CALL CHKXER('SGETRI',INFot,NOUt,LERr,OK)
!
!        SGETRS
!
         SRNamt = 'SGETRS'
         INFot = 1
         CALL SGETRS('/',0,0,a,1,ip,b,1,info)
         CALL CHKXER('SGETRS',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SGETRS('N',-1,0,a,1,ip,b,1,info)
         CALL CHKXER('SGETRS',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SGETRS('N',0,-1,a,1,ip,b,1,info)
         CALL CHKXER('SGETRS',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL SGETRS('N',2,1,a,1,ip,b,2,info)
         CALL CHKXER('SGETRS',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL SGETRS('N',2,1,a,2,ip,b,1,info)
         CALL CHKXER('SGETRS',INFot,NOUt,LERr,OK)
!
!        SGERFS
!
         SRNamt = 'SGERFS'
         INFot = 1
         CALL SGERFS('/',0,0,a,1,af,1,ip,b,1,x,1,r1,r2,w,iw,info)
         CALL CHKXER('SGERFS',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SGERFS('N',-1,0,a,1,af,1,ip,b,1,x,1,r1,r2,w,iw,info)
         CALL CHKXER('SGERFS',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SGERFS('N',0,-1,a,1,af,1,ip,b,1,x,1,r1,r2,w,iw,info)
         CALL CHKXER('SGERFS',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL SGERFS('N',2,1,a,1,af,2,ip,b,2,x,2,r1,r2,w,iw,info)
         CALL CHKXER('SGERFS',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL SGERFS('N',2,1,a,2,af,1,ip,b,2,x,2,r1,r2,w,iw,info)
         CALL CHKXER('SGERFS',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL SGERFS('N',2,1,a,2,af,2,ip,b,1,x,2,r1,r2,w,iw,info)
         CALL CHKXER('SGERFS',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL SGERFS('N',2,1,a,2,af,2,ip,b,2,x,1,r1,r2,w,iw,info)
         CALL CHKXER('SGERFS',INFot,NOUt,LERr,OK)
!
!        SGERFSX
!
         n_err_bnds = 3
         nparams = 0
         SRNamt = 'SGERFSX'
         INFot = 1
         CALL SGERFSX('/',eq,0,0,a,1,af,1,ip,r,c,b,1,x,1,rcond,berr,    &
     &                n_err_bnds,err_bnds_n,err_bnds_c,nparams,params,w,&
     &                iw,info)
         CALL CHKXER('SGERFSX',INFot,NOUt,LERr,OK)
         INFot = 2
         eq = '/'
         CALL SGERFSX('N',eq,2,1,a,1,af,2,ip,r,c,b,2,x,2,rcond,berr,    &
     &                n_err_bnds,err_bnds_n,err_bnds_c,nparams,params,w,&
     &                iw,info)
         CALL CHKXER('SGERFSX',INFot,NOUt,LERr,OK)
         INFot = 3
         eq = 'R'
         CALL SGERFSX('N',eq,-1,0,a,1,af,1,ip,r,c,b,1,x,1,rcond,berr,   &
     &                n_err_bnds,err_bnds_n,err_bnds_c,nparams,params,w,&
     &                iw,info)
         CALL CHKXER('SGERFSX',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL SGERFSX('N',eq,0,-1,a,1,af,1,ip,r,c,b,1,x,1,rcond,berr,   &
     &                n_err_bnds,err_bnds_n,err_bnds_c,nparams,params,w,&
     &                iw,info)
         CALL CHKXER('SGERFSX',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL SGERFSX('N',eq,2,1,a,1,af,2,ip,r,c,b,2,x,2,rcond,berr,    &
     &                n_err_bnds,err_bnds_n,err_bnds_c,nparams,params,w,&
     &                iw,info)
         CALL CHKXER('SGERFSX',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL SGERFSX('N',eq,2,1,a,2,af,1,ip,r,c,b,2,x,2,rcond,berr,    &
     &                n_err_bnds,err_bnds_n,err_bnds_c,nparams,params,w,&
     &                iw,info)
         CALL CHKXER('SGERFSX',INFot,NOUt,LERr,OK)
         INFot = 13
         eq = 'C'
         CALL SGERFSX('N',eq,2,1,a,2,af,2,ip,r,c,b,1,x,2,rcond,berr,    &
     &                n_err_bnds,err_bnds_n,err_bnds_c,nparams,params,w,&
     &                iw,info)
         CALL CHKXER('SGERFSX',INFot,NOUt,LERr,OK)
         INFot = 15
         CALL SGERFSX('N',eq,2,1,a,2,af,2,ip,r,c,b,2,x,1,rcond,berr,    &
     &                n_err_bnds,err_bnds_n,err_bnds_c,nparams,params,w,&
     &                iw,info)
         CALL CHKXER('SGERFSX',INFot,NOUt,LERr,OK)
!
!        SGECON
!
         SRNamt = 'SGECON'
         INFot = 1
         CALL SGECON('/',0,a,1,anrm,rcond,w,iw,info)
         CALL CHKXER('SGECON',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SGECON('1',-1,a,1,anrm,rcond,w,iw,info)
         CALL CHKXER('SGECON',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL SGECON('1',2,a,1,anrm,rcond,w,iw,info)
         CALL CHKXER('SGECON',INFot,NOUt,LERr,OK)
!
!        SGEEQU
!
         SRNamt = 'SGEEQU'
         INFot = 1
         CALL SGEEQU(-1,0,a,1,r1,r2,rcond,ccond,anrm,info)
         CALL CHKXER('SGEEQU',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SGEEQU(0,-1,a,1,r1,r2,rcond,ccond,anrm,info)
         CALL CHKXER('SGEEQU',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL SGEEQU(2,2,a,1,r1,r2,rcond,ccond,anrm,info)
         CALL CHKXER('SGEEQU',INFot,NOUt,LERr,OK)
!
!        SGEEQUB
!
         SRNamt = 'SGEEQUB'
         INFot = 1
         CALL SGEEQUB(-1,0,a,1,r1,r2,rcond,ccond,anrm,info)
         CALL CHKXER('SGEEQUB',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SGEEQUB(0,-1,a,1,r1,r2,rcond,ccond,anrm,info)
         CALL CHKXER('SGEEQUB',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL SGEEQUB(2,2,a,1,r1,r2,rcond,ccond,anrm,info)
         CALL CHKXER('SGEEQUB',INFot,NOUt,LERr,OK)
!
      ELSEIF ( LSAMEN(2,c2,'GB') ) THEN
!
!        Test error exits of the routines that use the LU decomposition
!        of a general band matrix.
!
!        SGBTRF
!
         SRNamt = 'SGBTRF'
         INFot = 1
         CALL SGBTRF(-1,0,0,0,a,1,ip,info)
         CALL CHKXER('SGBTRF',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SGBTRF(0,-1,0,0,a,1,ip,info)
         CALL CHKXER('SGBTRF',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SGBTRF(1,1,-1,0,a,1,ip,info)
         CALL CHKXER('SGBTRF',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL SGBTRF(1,1,0,-1,a,1,ip,info)
         CALL CHKXER('SGBTRF',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL SGBTRF(2,2,1,1,a,3,ip,info)
         CALL CHKXER('SGBTRF',INFot,NOUt,LERr,OK)
!
!        SGBTF2
!
         SRNamt = 'SGBTF2'
         INFot = 1
         CALL SGBTF2(-1,0,0,0,a,1,ip,info)
         CALL CHKXER('SGBTF2',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SGBTF2(0,-1,0,0,a,1,ip,info)
         CALL CHKXER('SGBTF2',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SGBTF2(1,1,-1,0,a,1,ip,info)
         CALL CHKXER('SGBTF2',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL SGBTF2(1,1,0,-1,a,1,ip,info)
         CALL CHKXER('SGBTF2',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL SGBTF2(2,2,1,1,a,3,ip,info)
         CALL CHKXER('SGBTF2',INFot,NOUt,LERr,OK)
!
!        SGBTRS
!
         SRNamt = 'SGBTRS'
         INFot = 1
         CALL SGBTRS('/',0,0,0,1,a,1,ip,b,1,info)
         CALL CHKXER('SGBTRS',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SGBTRS('N',-1,0,0,1,a,1,ip,b,1,info)
         CALL CHKXER('SGBTRS',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SGBTRS('N',1,-1,0,1,a,1,ip,b,1,info)
         CALL CHKXER('SGBTRS',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL SGBTRS('N',1,0,-1,1,a,1,ip,b,1,info)
         CALL CHKXER('SGBTRS',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL SGBTRS('N',1,0,0,-1,a,1,ip,b,1,info)
         CALL CHKXER('SGBTRS',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL SGBTRS('N',2,1,1,1,a,3,ip,b,2,info)
         CALL CHKXER('SGBTRS',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL SGBTRS('N',2,0,0,1,a,1,ip,b,1,info)
         CALL CHKXER('SGBTRS',INFot,NOUt,LERr,OK)
!
!        SGBRFS
!
         SRNamt = 'SGBRFS'
         INFot = 1
         CALL SGBRFS('/',0,0,0,0,a,1,af,1,ip,b,1,x,1,r1,r2,w,iw,info)
         CALL CHKXER('SGBRFS',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SGBRFS('N',-1,0,0,0,a,1,af,1,ip,b,1,x,1,r1,r2,w,iw,info)
         CALL CHKXER('SGBRFS',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SGBRFS('N',1,-1,0,0,a,1,af,1,ip,b,1,x,1,r1,r2,w,iw,info)
         CALL CHKXER('SGBRFS',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL SGBRFS('N',1,0,-1,0,a,1,af,1,ip,b,1,x,1,r1,r2,w,iw,info)
         CALL CHKXER('SGBRFS',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL SGBRFS('N',1,0,0,-1,a,1,af,1,ip,b,1,x,1,r1,r2,w,iw,info)
         CALL CHKXER('SGBRFS',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL SGBRFS('N',2,1,1,1,a,2,af,4,ip,b,2,x,2,r1,r2,w,iw,info)
         CALL CHKXER('SGBRFS',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL SGBRFS('N',2,1,1,1,a,3,af,3,ip,b,2,x,2,r1,r2,w,iw,info)
         CALL CHKXER('SGBRFS',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL SGBRFS('N',2,0,0,1,a,1,af,1,ip,b,1,x,2,r1,r2,w,iw,info)
         CALL CHKXER('SGBRFS',INFot,NOUt,LERr,OK)
         INFot = 14
         CALL SGBRFS('N',2,0,0,1,a,1,af,1,ip,b,2,x,1,r1,r2,w,iw,info)
         CALL CHKXER('SGBRFS',INFot,NOUt,LERr,OK)
!
!        SGBRFSX
!
         n_err_bnds = 3
         nparams = 0
         SRNamt = 'SGBRFSX'
         INFot = 1
         CALL SGBRFSX('/',eq,0,0,0,0,a,1,af,1,ip,r,c,b,1,x,1,rcond,berr,&
     &                n_err_bnds,err_bnds_n,err_bnds_c,nparams,params,w,&
     &                iw,info)
         CALL CHKXER('SGBRFSX',INFot,NOUt,LERr,OK)
         INFot = 2
         eq = '/'
         CALL SGBRFSX('N',eq,2,1,1,1,a,1,af,2,ip,r,c,b,2,x,2,rcond,berr,&
     &                n_err_bnds,err_bnds_n,err_bnds_c,nparams,params,w,&
     &                iw,info)
         CALL CHKXER('SGBRFSX',INFot,NOUt,LERr,OK)
         INFot = 3
         eq = 'R'
         CALL SGBRFSX('N',eq,-1,1,1,0,a,1,af,1,ip,r,c,b,1,x,1,rcond,    &
     &                berr,n_err_bnds,err_bnds_n,err_bnds_c,nparams,    &
     &                params,w,iw,info)
         CALL CHKXER('SGBRFSX',INFot,NOUt,LERr,OK)
         INFot = 4
         eq = 'R'
         CALL SGBRFSX('N',eq,2,-1,1,1,a,3,af,4,ip,r,c,b,1,x,1,rcond,    &
     &                berr,n_err_bnds,err_bnds_n,err_bnds_c,nparams,    &
     &                params,w,iw,info)
         CALL CHKXER('SGBRFSX',INFot,NOUt,LERr,OK)
         INFot = 5
         eq = 'R'
         CALL SGBRFSX('N',eq,2,1,-1,1,a,3,af,4,ip,r,c,b,1,x,1,rcond,    &
     &                berr,n_err_bnds,err_bnds_n,err_bnds_c,nparams,    &
     &                params,w,iw,info)
         CALL CHKXER('SGBRFSX',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL SGBRFSX('N',eq,0,0,0,-1,a,1,af,1,ip,r,c,b,1,x,1,rcond,    &
     &                berr,n_err_bnds,err_bnds_n,err_bnds_c,nparams,    &
     &                params,w,iw,info)
         CALL CHKXER('SGBRFSX',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL SGBRFSX('N',eq,2,1,1,1,a,1,af,2,ip,r,c,b,2,x,2,rcond,berr,&
     &                n_err_bnds,err_bnds_n,err_bnds_c,nparams,params,w,&
     &                iw,info)
         CALL CHKXER('SGBRFSX',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL SGBRFSX('N',eq,2,1,1,1,a,3,af,3,ip,r,c,b,2,x,2,rcond,berr,&
     &                n_err_bnds,err_bnds_n,err_bnds_c,nparams,params,w,&
     &                iw,info)
         CALL CHKXER('SGBRFSX',INFot,NOUt,LERr,OK)
         INFot = 13
         eq = 'C'
         CALL SGBRFSX('N',eq,2,1,1,1,a,3,af,5,ip,r,c,b,1,x,2,rcond,berr,&
     &                n_err_bnds,err_bnds_n,err_bnds_c,nparams,params,w,&
     &                iw,info)
         CALL CHKXER('SGBRFSX',INFot,NOUt,LERr,OK)
         INFot = 15
         CALL SGBRFSX('N',eq,2,1,1,1,a,3,af,5,ip,r,c,b,2,x,1,rcond,berr,&
     &                n_err_bnds,err_bnds_n,err_bnds_c,nparams,params,w,&
     &                iw,info)
         CALL CHKXER('SGBRFSX',INFot,NOUt,LERr,OK)
!
!        SGBCON
!
         SRNamt = 'SGBCON'
         INFot = 1
         CALL SGBCON('/',0,0,0,a,1,ip,anrm,rcond,w,iw,info)
         CALL CHKXER('SGBCON',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SGBCON('1',-1,0,0,a,1,ip,anrm,rcond,w,iw,info)
         CALL CHKXER('SGBCON',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SGBCON('1',1,-1,0,a,1,ip,anrm,rcond,w,iw,info)
         CALL CHKXER('SGBCON',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL SGBCON('1',1,0,-1,a,1,ip,anrm,rcond,w,iw,info)
         CALL CHKXER('SGBCON',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL SGBCON('1',2,1,1,a,3,ip,anrm,rcond,w,iw,info)
         CALL CHKXER('SGBCON',INFot,NOUt,LERr,OK)
!
!        SGBEQU
!
         SRNamt = 'SGBEQU'
         INFot = 1
         CALL SGBEQU(-1,0,0,0,a,1,r1,r2,rcond,ccond,anrm,info)
         CALL CHKXER('SGBEQU',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SGBEQU(0,-1,0,0,a,1,r1,r2,rcond,ccond,anrm,info)
         CALL CHKXER('SGBEQU',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SGBEQU(1,1,-1,0,a,1,r1,r2,rcond,ccond,anrm,info)
         CALL CHKXER('SGBEQU',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL SGBEQU(1,1,0,-1,a,1,r1,r2,rcond,ccond,anrm,info)
         CALL CHKXER('SGBEQU',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL SGBEQU(2,2,1,1,a,2,r1,r2,rcond,ccond,anrm,info)
         CALL CHKXER('SGBEQU',INFot,NOUt,LERr,OK)
!
!        SGBEQUB
!
         SRNamt = 'SGBEQUB'
         INFot = 1
         CALL SGBEQUB(-1,0,0,0,a,1,r1,r2,rcond,ccond,anrm,info)
         CALL CHKXER('SGBEQUB',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SGBEQUB(0,-1,0,0,a,1,r1,r2,rcond,ccond,anrm,info)
         CALL CHKXER('SGBEQUB',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SGBEQUB(1,1,-1,0,a,1,r1,r2,rcond,ccond,anrm,info)
         CALL CHKXER('SGBEQUB',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL SGBEQUB(1,1,0,-1,a,1,r1,r2,rcond,ccond,anrm,info)
         CALL CHKXER('SGBEQUB',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL SGBEQUB(2,2,1,1,a,2,r1,r2,rcond,ccond,anrm,info)
         CALL CHKXER('SGBEQUB',INFot,NOUt,LERr,OK)
      ENDIF
!
!     Print a summary line.
!
      CALL ALAESM(Path,OK,NOUt)
!
!
!     End of SERRGE
!
      END SUBROUTINE SERRGE
