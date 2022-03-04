!*==derrge.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b DERRGEX
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DERRGE( PATH, NUNIT )
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
!> DERRGE tests the error exits for the DOUBLE PRECISION routines
!> for general matrices.
!>
!> Note that this file is used only when the XBLAS are available,
!> otherwise derrge.f defines this subroutine.
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
      SUBROUTINE DERRGE(Path,Nunit)
      IMPLICIT NONE
!*--DERRGE62
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
      DOUBLE PRECISION anrm , ccond , rcond , berr
!     ..
!     .. Local Arrays ..
      INTEGER ip(NMAX) , iw(NMAX)
      DOUBLE PRECISION a(NMAX,NMAX) , af(NMAX,NMAX) , b(NMAX) , c(NMAX) &
     &                 , r(NMAX) , r1(NMAX) , r2(NMAX) , w(LW) , x(NMAX)&
     &                 , err_bnds_n(NMAX,3) , err_bnds_c(NMAX,3) ,      &
     &                 params(1)
!     ..
!     .. External Functions ..
      LOGICAL LSAMEN
      EXTERNAL LSAMEN
!     ..
!     .. External Subroutines ..
      EXTERNAL ALAESM , CHKXER , DGBCON , DGBEQU , DGBRFS , DGBTF2 ,    &
     &         DGBTRF , DGBTRS , DGECON , DGEEQU , DGERFS , DGETF2 ,    &
     &         DGETRF , DGETRI , DGETRS , DGEEQUB , DGERFSX , DGBEQUB , &
     &         DGBRFSX
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
         b(j) = 0.D0
         r1(j) = 0.D0
         r2(j) = 0.D0
         w(j) = 0.D0
         x(j) = 0.D0
         c(j) = 0.D0
         r(j) = 0.D0
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
!        DGETRF
!
         SRNamt = 'DGETRF'
         INFot = 1
         CALL DGETRF(-1,0,a,1,ip,info)
         CALL CHKXER('DGETRF',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DGETRF(0,-1,a,1,ip,info)
         CALL CHKXER('DGETRF',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL DGETRF(2,1,a,1,ip,info)
         CALL CHKXER('DGETRF',INFot,NOUt,LERr,OK)
!
!        DGETF2
!
         SRNamt = 'DGETF2'
         INFot = 1
         CALL DGETF2(-1,0,a,1,ip,info)
         CALL CHKXER('DGETF2',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DGETF2(0,-1,a,1,ip,info)
         CALL CHKXER('DGETF2',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL DGETF2(2,1,a,1,ip,info)
         CALL CHKXER('DGETF2',INFot,NOUt,LERr,OK)
!
!        DGETRI
!
         SRNamt = 'DGETRI'
         INFot = 1
         CALL DGETRI(-1,a,1,ip,w,LW,info)
         CALL CHKXER('DGETRI',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DGETRI(2,a,1,ip,w,LW,info)
         CALL CHKXER('DGETRI',INFot,NOUt,LERr,OK)
!
!        DGETRS
!
         SRNamt = 'DGETRS'
         INFot = 1
         CALL DGETRS('/',0,0,a,1,ip,b,1,info)
         CALL CHKXER('DGETRS',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DGETRS('N',-1,0,a,1,ip,b,1,info)
         CALL CHKXER('DGETRS',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DGETRS('N',0,-1,a,1,ip,b,1,info)
         CALL CHKXER('DGETRS',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL DGETRS('N',2,1,a,1,ip,b,2,info)
         CALL CHKXER('DGETRS',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL DGETRS('N',2,1,a,2,ip,b,1,info)
         CALL CHKXER('DGETRS',INFot,NOUt,LERr,OK)
!
!        DGERFS
!
         SRNamt = 'DGERFS'
         INFot = 1
         CALL DGERFS('/',0,0,a,1,af,1,ip,b,1,x,1,r1,r2,w,iw,info)
         CALL CHKXER('DGERFS',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DGERFS('N',-1,0,a,1,af,1,ip,b,1,x,1,r1,r2,w,iw,info)
         CALL CHKXER('DGERFS',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DGERFS('N',0,-1,a,1,af,1,ip,b,1,x,1,r1,r2,w,iw,info)
         CALL CHKXER('DGERFS',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL DGERFS('N',2,1,a,1,af,2,ip,b,2,x,2,r1,r2,w,iw,info)
         CALL CHKXER('DGERFS',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL DGERFS('N',2,1,a,2,af,1,ip,b,2,x,2,r1,r2,w,iw,info)
         CALL CHKXER('DGERFS',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL DGERFS('N',2,1,a,2,af,2,ip,b,1,x,2,r1,r2,w,iw,info)
         CALL CHKXER('DGERFS',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL DGERFS('N',2,1,a,2,af,2,ip,b,2,x,1,r1,r2,w,iw,info)
         CALL CHKXER('DGERFS',INFot,NOUt,LERr,OK)
!
!        DGERFSX
!
         n_err_bnds = 3
         nparams = 0
         SRNamt = 'DGERFSX'
         INFot = 1
         CALL DGERFSX('/',eq,0,0,a,1,af,1,ip,r,c,b,1,x,1,rcond,berr,    &
     &                n_err_bnds,err_bnds_n,err_bnds_c,nparams,params,w,&
     &                iw,info)
         CALL CHKXER('DGERFSX',INFot,NOUt,LERr,OK)
         INFot = 2
         eq = '/'
         CALL DGERFSX('N',eq,2,1,a,1,af,2,ip,r,c,b,2,x,2,rcond,berr,    &
     &                n_err_bnds,err_bnds_n,err_bnds_c,nparams,params,w,&
     &                iw,info)
         CALL CHKXER('DGERFSX',INFot,NOUt,LERr,OK)
         INFot = 3
         eq = 'R'
         CALL DGERFSX('N',eq,-1,0,a,1,af,1,ip,r,c,b,1,x,1,rcond,berr,   &
     &                n_err_bnds,err_bnds_n,err_bnds_c,nparams,params,w,&
     &                iw,info)
         CALL CHKXER('DGERFSX',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL DGERFSX('N',eq,0,-1,a,1,af,1,ip,r,c,b,1,x,1,rcond,berr,   &
     &                n_err_bnds,err_bnds_n,err_bnds_c,nparams,params,w,&
     &                iw,info)
         CALL CHKXER('DGERFSX',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL DGERFSX('N',eq,2,1,a,1,af,2,ip,r,c,b,2,x,2,rcond,berr,    &
     &                n_err_bnds,err_bnds_n,err_bnds_c,nparams,params,w,&
     &                iw,info)
         CALL CHKXER('DGERFSX',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL DGERFSX('N',eq,2,1,a,2,af,1,ip,r,c,b,2,x,2,rcond,berr,    &
     &                n_err_bnds,err_bnds_n,err_bnds_c,nparams,params,w,&
     &                iw,info)
         CALL CHKXER('DGERFSX',INFot,NOUt,LERr,OK)
         INFot = 13
         eq = 'C'
         CALL DGERFSX('N',eq,2,1,a,2,af,2,ip,r,c,b,1,x,2,rcond,berr,    &
     &                n_err_bnds,err_bnds_n,err_bnds_c,nparams,params,w,&
     &                iw,info)
         CALL CHKXER('DGERFSX',INFot,NOUt,LERr,OK)
         INFot = 15
         CALL DGERFSX('N',eq,2,1,a,2,af,2,ip,r,c,b,2,x,1,rcond,berr,    &
     &                n_err_bnds,err_bnds_n,err_bnds_c,nparams,params,w,&
     &                iw,info)
         CALL CHKXER('DGERFSX',INFot,NOUt,LERr,OK)
!
!        DGECON
!
         SRNamt = 'DGECON'
         INFot = 1
         CALL DGECON('/',0,a,1,anrm,rcond,w,iw,info)
         CALL CHKXER('DGECON',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DGECON('1',-1,a,1,anrm,rcond,w,iw,info)
         CALL CHKXER('DGECON',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL DGECON('1',2,a,1,anrm,rcond,w,iw,info)
         CALL CHKXER('DGECON',INFot,NOUt,LERr,OK)
!
!        DGEEQU
!
         SRNamt = 'DGEEQU'
         INFot = 1
         CALL DGEEQU(-1,0,a,1,r1,r2,rcond,ccond,anrm,info)
         CALL CHKXER('DGEEQU',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DGEEQU(0,-1,a,1,r1,r2,rcond,ccond,anrm,info)
         CALL CHKXER('DGEEQU',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL DGEEQU(2,2,a,1,r1,r2,rcond,ccond,anrm,info)
         CALL CHKXER('DGEEQU',INFot,NOUt,LERr,OK)
!
!        DGEEQUB
!
         SRNamt = 'DGEEQUB'
         INFot = 1
         CALL DGEEQUB(-1,0,a,1,r1,r2,rcond,ccond,anrm,info)
         CALL CHKXER('DGEEQUB',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DGEEQUB(0,-1,a,1,r1,r2,rcond,ccond,anrm,info)
         CALL CHKXER('DGEEQUB',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL DGEEQUB(2,2,a,1,r1,r2,rcond,ccond,anrm,info)
         CALL CHKXER('DGEEQUB',INFot,NOUt,LERr,OK)
!
      ELSEIF ( LSAMEN(2,c2,'GB') ) THEN
!
!        Test error exits of the routines that use the LU decomposition
!        of a general band matrix.
!
!        DGBTRF
!
         SRNamt = 'DGBTRF'
         INFot = 1
         CALL DGBTRF(-1,0,0,0,a,1,ip,info)
         CALL CHKXER('DGBTRF',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DGBTRF(0,-1,0,0,a,1,ip,info)
         CALL CHKXER('DGBTRF',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DGBTRF(1,1,-1,0,a,1,ip,info)
         CALL CHKXER('DGBTRF',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL DGBTRF(1,1,0,-1,a,1,ip,info)
         CALL CHKXER('DGBTRF',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL DGBTRF(2,2,1,1,a,3,ip,info)
         CALL CHKXER('DGBTRF',INFot,NOUt,LERr,OK)
!
!        DGBTF2
!
         SRNamt = 'DGBTF2'
         INFot = 1
         CALL DGBTF2(-1,0,0,0,a,1,ip,info)
         CALL CHKXER('DGBTF2',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DGBTF2(0,-1,0,0,a,1,ip,info)
         CALL CHKXER('DGBTF2',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DGBTF2(1,1,-1,0,a,1,ip,info)
         CALL CHKXER('DGBTF2',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL DGBTF2(1,1,0,-1,a,1,ip,info)
         CALL CHKXER('DGBTF2',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL DGBTF2(2,2,1,1,a,3,ip,info)
         CALL CHKXER('DGBTF2',INFot,NOUt,LERr,OK)
!
!        DGBTRS
!
         SRNamt = 'DGBTRS'
         INFot = 1
         CALL DGBTRS('/',0,0,0,1,a,1,ip,b,1,info)
         CALL CHKXER('DGBTRS',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DGBTRS('N',-1,0,0,1,a,1,ip,b,1,info)
         CALL CHKXER('DGBTRS',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DGBTRS('N',1,-1,0,1,a,1,ip,b,1,info)
         CALL CHKXER('DGBTRS',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL DGBTRS('N',1,0,-1,1,a,1,ip,b,1,info)
         CALL CHKXER('DGBTRS',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL DGBTRS('N',1,0,0,-1,a,1,ip,b,1,info)
         CALL CHKXER('DGBTRS',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL DGBTRS('N',2,1,1,1,a,3,ip,b,2,info)
         CALL CHKXER('DGBTRS',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL DGBTRS('N',2,0,0,1,a,1,ip,b,1,info)
         CALL CHKXER('DGBTRS',INFot,NOUt,LERr,OK)
!
!        DGBRFS
!
         SRNamt = 'DGBRFS'
         INFot = 1
         CALL DGBRFS('/',0,0,0,0,a,1,af,1,ip,b,1,x,1,r1,r2,w,iw,info)
         CALL CHKXER('DGBRFS',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DGBRFS('N',-1,0,0,0,a,1,af,1,ip,b,1,x,1,r1,r2,w,iw,info)
         CALL CHKXER('DGBRFS',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DGBRFS('N',1,-1,0,0,a,1,af,1,ip,b,1,x,1,r1,r2,w,iw,info)
         CALL CHKXER('DGBRFS',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL DGBRFS('N',1,0,-1,0,a,1,af,1,ip,b,1,x,1,r1,r2,w,iw,info)
         CALL CHKXER('DGBRFS',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL DGBRFS('N',1,0,0,-1,a,1,af,1,ip,b,1,x,1,r1,r2,w,iw,info)
         CALL CHKXER('DGBRFS',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL DGBRFS('N',2,1,1,1,a,2,af,4,ip,b,2,x,2,r1,r2,w,iw,info)
         CALL CHKXER('DGBRFS',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL DGBRFS('N',2,1,1,1,a,3,af,3,ip,b,2,x,2,r1,r2,w,iw,info)
         CALL CHKXER('DGBRFS',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL DGBRFS('N',2,0,0,1,a,1,af,1,ip,b,1,x,2,r1,r2,w,iw,info)
         CALL CHKXER('DGBRFS',INFot,NOUt,LERr,OK)
         INFot = 14
         CALL DGBRFS('N',2,0,0,1,a,1,af,1,ip,b,2,x,1,r1,r2,w,iw,info)
         CALL CHKXER('DGBRFS',INFot,NOUt,LERr,OK)
!
!        DGBRFSX
!
         n_err_bnds = 3
         nparams = 0
         SRNamt = 'DGBRFSX'
         INFot = 1
         CALL DGBRFSX('/',eq,0,0,0,0,a,1,af,1,ip,r,c,b,1,x,1,rcond,berr,&
     &                n_err_bnds,err_bnds_n,err_bnds_c,nparams,params,w,&
     &                iw,info)
         CALL CHKXER('DGBRFSX',INFot,NOUt,LERr,OK)
         INFot = 2
         eq = '/'
         CALL DGBRFSX('N',eq,2,1,1,1,a,1,af,2,ip,r,c,b,2,x,2,rcond,berr,&
     &                n_err_bnds,err_bnds_n,err_bnds_c,nparams,params,w,&
     &                iw,info)
         CALL CHKXER('DGBRFSX',INFot,NOUt,LERr,OK)
         INFot = 3
         eq = 'R'
         CALL DGBRFSX('N',eq,-1,1,1,0,a,1,af,1,ip,r,c,b,1,x,1,rcond,    &
     &                berr,n_err_bnds,err_bnds_n,err_bnds_c,nparams,    &
     &                params,w,iw,info)
         CALL CHKXER('DGBRFSX',INFot,NOUt,LERr,OK)
         INFot = 4
         eq = 'R'
         CALL DGBRFSX('N',eq,2,-1,1,1,a,3,af,4,ip,r,c,b,1,x,1,rcond,    &
     &                berr,n_err_bnds,err_bnds_n,err_bnds_c,nparams,    &
     &                params,w,iw,info)
         CALL CHKXER('DGBRFSX',INFot,NOUt,LERr,OK)
         INFot = 5
         eq = 'R'
         CALL DGBRFSX('N',eq,2,1,-1,1,a,3,af,4,ip,r,c,b,1,x,1,rcond,    &
     &                berr,n_err_bnds,err_bnds_n,err_bnds_c,nparams,    &
     &                params,w,iw,info)
         CALL CHKXER('DGBRFSX',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL DGBRFSX('N',eq,0,0,0,-1,a,1,af,1,ip,r,c,b,1,x,1,rcond,    &
     &                berr,n_err_bnds,err_bnds_n,err_bnds_c,nparams,    &
     &                params,w,iw,info)
         CALL CHKXER('DGBRFSX',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL DGBRFSX('N',eq,2,1,1,1,a,1,af,2,ip,r,c,b,2,x,2,rcond,berr,&
     &                n_err_bnds,err_bnds_n,err_bnds_c,nparams,params,w,&
     &                iw,info)
         CALL CHKXER('DGBRFSX',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL DGBRFSX('N',eq,2,1,1,1,a,3,af,3,ip,r,c,b,2,x,2,rcond,berr,&
     &                n_err_bnds,err_bnds_n,err_bnds_c,nparams,params,w,&
     &                iw,info)
         CALL CHKXER('DGBRFSX',INFot,NOUt,LERr,OK)
         INFot = 13
         eq = 'C'
         CALL DGBRFSX('N',eq,2,1,1,1,a,3,af,5,ip,r,c,b,1,x,2,rcond,berr,&
     &                n_err_bnds,err_bnds_n,err_bnds_c,nparams,params,w,&
     &                iw,info)
         CALL CHKXER('DGBRFSX',INFot,NOUt,LERr,OK)
         INFot = 15
         CALL DGBRFSX('N',eq,2,1,1,1,a,3,af,5,ip,r,c,b,2,x,1,rcond,berr,&
     &                n_err_bnds,err_bnds_n,err_bnds_c,nparams,params,w,&
     &                iw,info)
         CALL CHKXER('DGBRFSX',INFot,NOUt,LERr,OK)
!
!        DGBCON
!
         SRNamt = 'DGBCON'
         INFot = 1
         CALL DGBCON('/',0,0,0,a,1,ip,anrm,rcond,w,iw,info)
         CALL CHKXER('DGBCON',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DGBCON('1',-1,0,0,a,1,ip,anrm,rcond,w,iw,info)
         CALL CHKXER('DGBCON',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DGBCON('1',1,-1,0,a,1,ip,anrm,rcond,w,iw,info)
         CALL CHKXER('DGBCON',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL DGBCON('1',1,0,-1,a,1,ip,anrm,rcond,w,iw,info)
         CALL CHKXER('DGBCON',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL DGBCON('1',2,1,1,a,3,ip,anrm,rcond,w,iw,info)
         CALL CHKXER('DGBCON',INFot,NOUt,LERr,OK)
!
!        DGBEQU
!
         SRNamt = 'DGBEQU'
         INFot = 1
         CALL DGBEQU(-1,0,0,0,a,1,r1,r2,rcond,ccond,anrm,info)
         CALL CHKXER('DGBEQU',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DGBEQU(0,-1,0,0,a,1,r1,r2,rcond,ccond,anrm,info)
         CALL CHKXER('DGBEQU',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DGBEQU(1,1,-1,0,a,1,r1,r2,rcond,ccond,anrm,info)
         CALL CHKXER('DGBEQU',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL DGBEQU(1,1,0,-1,a,1,r1,r2,rcond,ccond,anrm,info)
         CALL CHKXER('DGBEQU',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL DGBEQU(2,2,1,1,a,2,r1,r2,rcond,ccond,anrm,info)
         CALL CHKXER('DGBEQU',INFot,NOUt,LERr,OK)
!
!        DGBEQUB
!
         SRNamt = 'DGBEQUB'
         INFot = 1
         CALL DGBEQUB(-1,0,0,0,a,1,r1,r2,rcond,ccond,anrm,info)
         CALL CHKXER('DGBEQUB',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DGBEQUB(0,-1,0,0,a,1,r1,r2,rcond,ccond,anrm,info)
         CALL CHKXER('DGBEQUB',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DGBEQUB(1,1,-1,0,a,1,r1,r2,rcond,ccond,anrm,info)
         CALL CHKXER('DGBEQUB',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL DGBEQUB(1,1,0,-1,a,1,r1,r2,rcond,ccond,anrm,info)
         CALL CHKXER('DGBEQUB',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL DGBEQUB(2,2,1,1,a,2,r1,r2,rcond,ccond,anrm,info)
         CALL CHKXER('DGBEQUB',INFot,NOUt,LERr,OK)
      ENDIF
!
!     Print a summary line.
!
      CALL ALAESM(Path,OK,NOUt)
!
!
!     End of DERRGE
!
      END SUBROUTINE DERRGE
