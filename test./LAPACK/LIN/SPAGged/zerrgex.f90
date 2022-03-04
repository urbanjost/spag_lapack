!*==zerrge.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b ZERRGEX
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZERRGE( PATH, NUNIT )
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
!> ZERRGE tests the error exits for the COMPLEX*16 routines
!> for general matrices.
!>
!> Note that this file is used only when the XBLAS are available,
!> otherwise zerrge.f defines this subroutine.
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
!> \ingroup complex16_lin
!
!  =====================================================================
      SUBROUTINE ZERRGE(Path,Nunit)
      IMPLICIT NONE
!*--ZERRGE62
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
      DOUBLE PRECISION anrm , ccond , rcond , berr
!     ..
!     .. Local Arrays ..
      INTEGER ip(NMAX)
      DOUBLE PRECISION r(NMAX) , r1(NMAX) , r2(NMAX) , cs(NMAX) ,       &
     &                 rs(NMAX)
      COMPLEX*16 a(NMAX,NMAX) , af(NMAX,NMAX) , b(NMAX) , w(2*NMAX) ,   &
     &           x(NMAX) , err_bnds_n(NMAX,3) , err_bnds_c(NMAX,3) ,    &
     &           params
!     ..
!     .. External Functions ..
      LOGICAL LSAMEN
      EXTERNAL LSAMEN
!     ..
!     .. External Subroutines ..
      EXTERNAL ALAESM , CHKXER , ZGBCON , ZGBEQU , ZGBRFS , ZGBTF2 ,    &
     &         ZGBTRF , ZGBTRS , ZGECON , ZGEEQU , ZGERFS , ZGETF2 ,    &
     &         ZGETRF , ZGETRI , ZGETRS , ZGEEQUB , ZGERFSX , ZGBEQUB , &
     &         ZGBRFSX
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
      INTRINSIC DBLE , DCMPLX
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
            a(i,j) = DCMPLX(1.D0/DBLE(i+j),-1.D0/DBLE(i+j))
            af(i,j) = DCMPLX(1.D0/DBLE(i+j),-1.D0/DBLE(i+j))
         ENDDO
         b(j) = 0.D0
         r1(j) = 0.D0
         r2(j) = 0.D0
         w(j) = 0.D0
         x(j) = 0.D0
         cs(j) = 0.D0
         rs(j) = 0.D0
         ip(j) = j
      ENDDO
      OK = .TRUE.
!
!     Test error exits of the routines that use the LU decomposition
!     of a general matrix.
!
      IF ( LSAMEN(2,c2,'GE') ) THEN
!
!        ZGETRF
!
         SRNamt = 'ZGETRF'
         INFot = 1
         CALL ZGETRF(-1,0,a,1,ip,info)
         CALL CHKXER('ZGETRF',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZGETRF(0,-1,a,1,ip,info)
         CALL CHKXER('ZGETRF',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZGETRF(2,1,a,1,ip,info)
         CALL CHKXER('ZGETRF',INFot,NOUt,LERr,OK)
!
!        ZGETF2
!
         SRNamt = 'ZGETF2'
         INFot = 1
         CALL ZGETF2(-1,0,a,1,ip,info)
         CALL CHKXER('ZGETF2',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZGETF2(0,-1,a,1,ip,info)
         CALL CHKXER('ZGETF2',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZGETF2(2,1,a,1,ip,info)
         CALL CHKXER('ZGETF2',INFot,NOUt,LERr,OK)
!
!        ZGETRI
!
         SRNamt = 'ZGETRI'
         INFot = 1
         CALL ZGETRI(-1,a,1,ip,w,1,info)
         CALL CHKXER('ZGETRI',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZGETRI(2,a,1,ip,w,2,info)
         CALL CHKXER('ZGETRI',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL ZGETRI(2,a,2,ip,w,1,info)
         CALL CHKXER('ZGETRI',INFot,NOUt,LERr,OK)
!
!        ZGETRS
!
         SRNamt = 'ZGETRS'
         INFot = 1
         CALL ZGETRS('/',0,0,a,1,ip,b,1,info)
         CALL CHKXER('ZGETRS',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZGETRS('N',-1,0,a,1,ip,b,1,info)
         CALL CHKXER('ZGETRS',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZGETRS('N',0,-1,a,1,ip,b,1,info)
         CALL CHKXER('ZGETRS',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL ZGETRS('N',2,1,a,1,ip,b,2,info)
         CALL CHKXER('ZGETRS',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL ZGETRS('N',2,1,a,2,ip,b,1,info)
         CALL CHKXER('ZGETRS',INFot,NOUt,LERr,OK)
!
!        ZGERFS
!
         SRNamt = 'ZGERFS'
         INFot = 1
         CALL ZGERFS('/',0,0,a,1,af,1,ip,b,1,x,1,r1,r2,w,r,info)
         CALL CHKXER('ZGERFS',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZGERFS('N',-1,0,a,1,af,1,ip,b,1,x,1,r1,r2,w,r,info)
         CALL CHKXER('ZGERFS',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZGERFS('N',0,-1,a,1,af,1,ip,b,1,x,1,r1,r2,w,r,info)
         CALL CHKXER('ZGERFS',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL ZGERFS('N',2,1,a,1,af,2,ip,b,2,x,2,r1,r2,w,r,info)
         CALL CHKXER('ZGERFS',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL ZGERFS('N',2,1,a,2,af,1,ip,b,2,x,2,r1,r2,w,r,info)
         CALL CHKXER('ZGERFS',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL ZGERFS('N',2,1,a,2,af,2,ip,b,1,x,2,r1,r2,w,r,info)
         CALL CHKXER('ZGERFS',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL ZGERFS('N',2,1,a,2,af,2,ip,b,2,x,1,r1,r2,w,r,info)
         CALL CHKXER('ZGERFS',INFot,NOUt,LERr,OK)
!
!        ZGERFSX
!
         n_err_bnds = 3
         nparams = 0
         SRNamt = 'ZGERFSX'
         INFot = 1
         CALL ZGERFSX('/',eq,0,0,a,1,af,1,ip,rs,cs,b,1,x,1,rcond,berr,  &
     &                n_err_bnds,err_bnds_n,err_bnds_c,nparams,params,w,&
     &                r,info)
         CALL CHKXER('ZGERFSX',INFot,NOUt,LERr,OK)
         INFot = 2
         eq = '/'
         CALL ZGERFSX('N',eq,2,1,a,1,af,2,ip,rs,cs,b,2,x,2,rcond,berr,  &
     &                n_err_bnds,err_bnds_n,err_bnds_c,nparams,params,w,&
     &                r,info)
         CALL CHKXER('ZGERFSX',INFot,NOUt,LERr,OK)
         INFot = 3
         eq = 'R'
         CALL ZGERFSX('N',eq,-1,0,a,1,af,1,ip,rs,cs,b,1,x,1,rcond,berr, &
     &                n_err_bnds,err_bnds_n,err_bnds_c,nparams,params,w,&
     &                r,info)
         CALL CHKXER('ZGERFSX',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZGERFSX('N',eq,0,-1,a,1,af,1,ip,rs,cs,b,1,x,1,rcond,berr, &
     &                n_err_bnds,err_bnds_n,err_bnds_c,nparams,params,w,&
     &                r,info)
         CALL CHKXER('ZGERFSX',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL ZGERFSX('N',eq,2,1,a,1,af,2,ip,rs,cs,b,2,x,2,rcond,berr,  &
     &                n_err_bnds,err_bnds_n,err_bnds_c,nparams,params,w,&
     &                r,info)
         CALL CHKXER('ZGERFSX',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL ZGERFSX('N',eq,2,1,a,2,af,1,ip,rs,cs,b,2,x,2,rcond,berr,  &
     &                n_err_bnds,err_bnds_n,err_bnds_c,nparams,params,w,&
     &                r,info)
         CALL CHKXER('ZGERFSX',INFot,NOUt,LERr,OK)
         INFot = 13
         eq = 'C'
         CALL ZGERFSX('N',eq,2,1,a,2,af,2,ip,rs,cs,b,1,x,2,rcond,berr,  &
     &                n_err_bnds,err_bnds_n,err_bnds_c,nparams,params,w,&
     &                r,info)
         CALL CHKXER('ZGERFSX',INFot,NOUt,LERr,OK)
         INFot = 15
         CALL ZGERFSX('N',eq,2,1,a,2,af,2,ip,rs,cs,b,2,x,1,rcond,berr,  &
     &                n_err_bnds,err_bnds_n,err_bnds_c,nparams,params,w,&
     &                r,info)
         CALL CHKXER('ZGERFSX',INFot,NOUt,LERr,OK)
!
!        ZGECON
!
         SRNamt = 'ZGECON'
         INFot = 1
         CALL ZGECON('/',0,a,1,anrm,rcond,w,r,info)
         CALL CHKXER('ZGECON',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZGECON('1',-1,a,1,anrm,rcond,w,r,info)
         CALL CHKXER('ZGECON',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZGECON('1',2,a,1,anrm,rcond,w,r,info)
         CALL CHKXER('ZGECON',INFot,NOUt,LERr,OK)
!
!        ZGEEQU
!
         SRNamt = 'ZGEEQU'
         INFot = 1
         CALL ZGEEQU(-1,0,a,1,r1,r2,rcond,ccond,anrm,info)
         CALL CHKXER('ZGEEQU',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZGEEQU(0,-1,a,1,r1,r2,rcond,ccond,anrm,info)
         CALL CHKXER('ZGEEQU',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZGEEQU(2,2,a,1,r1,r2,rcond,ccond,anrm,info)
         CALL CHKXER('ZGEEQU',INFot,NOUt,LERr,OK)
!
!        ZGEEQUB
!
         SRNamt = 'ZGEEQUB'
         INFot = 1
         CALL ZGEEQUB(-1,0,a,1,r1,r2,rcond,ccond,anrm,info)
         CALL CHKXER('ZGEEQUB',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZGEEQUB(0,-1,a,1,r1,r2,rcond,ccond,anrm,info)
         CALL CHKXER('ZGEEQUB',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZGEEQUB(2,2,a,1,r1,r2,rcond,ccond,anrm,info)
         CALL CHKXER('ZGEEQUB',INFot,NOUt,LERr,OK)
!
!     Test error exits of the routines that use the LU decomposition
!     of a general band matrix.
!
      ELSEIF ( LSAMEN(2,c2,'GB') ) THEN
!
!        ZGBTRF
!
         SRNamt = 'ZGBTRF'
         INFot = 1
         CALL ZGBTRF(-1,0,0,0,a,1,ip,info)
         CALL CHKXER('ZGBTRF',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZGBTRF(0,-1,0,0,a,1,ip,info)
         CALL CHKXER('ZGBTRF',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZGBTRF(1,1,-1,0,a,1,ip,info)
         CALL CHKXER('ZGBTRF',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZGBTRF(1,1,0,-1,a,1,ip,info)
         CALL CHKXER('ZGBTRF',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL ZGBTRF(2,2,1,1,a,3,ip,info)
         CALL CHKXER('ZGBTRF',INFot,NOUt,LERr,OK)
!
!        ZGBTF2
!
         SRNamt = 'ZGBTF2'
         INFot = 1
         CALL ZGBTF2(-1,0,0,0,a,1,ip,info)
         CALL CHKXER('ZGBTF2',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZGBTF2(0,-1,0,0,a,1,ip,info)
         CALL CHKXER('ZGBTF2',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZGBTF2(1,1,-1,0,a,1,ip,info)
         CALL CHKXER('ZGBTF2',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZGBTF2(1,1,0,-1,a,1,ip,info)
         CALL CHKXER('ZGBTF2',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL ZGBTF2(2,2,1,1,a,3,ip,info)
         CALL CHKXER('ZGBTF2',INFot,NOUt,LERr,OK)
!
!        ZGBTRS
!
         SRNamt = 'ZGBTRS'
         INFot = 1
         CALL ZGBTRS('/',0,0,0,1,a,1,ip,b,1,info)
         CALL CHKXER('ZGBTRS',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZGBTRS('N',-1,0,0,1,a,1,ip,b,1,info)
         CALL CHKXER('ZGBTRS',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZGBTRS('N',1,-1,0,1,a,1,ip,b,1,info)
         CALL CHKXER('ZGBTRS',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZGBTRS('N',1,0,-1,1,a,1,ip,b,1,info)
         CALL CHKXER('ZGBTRS',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL ZGBTRS('N',1,0,0,-1,a,1,ip,b,1,info)
         CALL CHKXER('ZGBTRS',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL ZGBTRS('N',2,1,1,1,a,3,ip,b,2,info)
         CALL CHKXER('ZGBTRS',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL ZGBTRS('N',2,0,0,1,a,1,ip,b,1,info)
         CALL CHKXER('ZGBTRS',INFot,NOUt,LERr,OK)
!
!        ZGBRFS
!
         SRNamt = 'ZGBRFS'
         INFot = 1
         CALL ZGBRFS('/',0,0,0,0,a,1,af,1,ip,b,1,x,1,r1,r2,w,r,info)
         CALL CHKXER('ZGBRFS',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZGBRFS('N',-1,0,0,0,a,1,af,1,ip,b,1,x,1,r1,r2,w,r,info)
         CALL CHKXER('ZGBRFS',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZGBRFS('N',1,-1,0,0,a,1,af,1,ip,b,1,x,1,r1,r2,w,r,info)
         CALL CHKXER('ZGBRFS',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZGBRFS('N',1,0,-1,0,a,1,af,1,ip,b,1,x,1,r1,r2,w,r,info)
         CALL CHKXER('ZGBRFS',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL ZGBRFS('N',1,0,0,-1,a,1,af,1,ip,b,1,x,1,r1,r2,w,r,info)
         CALL CHKXER('ZGBRFS',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL ZGBRFS('N',2,1,1,1,a,2,af,4,ip,b,2,x,2,r1,r2,w,r,info)
         CALL CHKXER('ZGBRFS',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL ZGBRFS('N',2,1,1,1,a,3,af,3,ip,b,2,x,2,r1,r2,w,r,info)
         CALL CHKXER('ZGBRFS',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL ZGBRFS('N',2,0,0,1,a,1,af,1,ip,b,1,x,2,r1,r2,w,r,info)
         CALL CHKXER('ZGBRFS',INFot,NOUt,LERr,OK)
         INFot = 14
         CALL ZGBRFS('N',2,0,0,1,a,1,af,1,ip,b,2,x,1,r1,r2,w,r,info)
         CALL CHKXER('ZGBRFS',INFot,NOUt,LERr,OK)
!
!        ZGBRFSX
!
         n_err_bnds = 3
         nparams = 0
         SRNamt = 'ZGBRFSX'
         INFot = 1
         CALL ZGBRFSX('/',eq,0,0,0,0,a,1,af,1,ip,rs,cs,b,1,x,1,rcond,   &
     &                berr,n_err_bnds,err_bnds_n,err_bnds_c,nparams,    &
     &                params,w,r,info)
         CALL CHKXER('ZGBRFSX',INFot,NOUt,LERr,OK)
         INFot = 2
         eq = '/'
         CALL ZGBRFSX('N',eq,2,1,1,1,a,1,af,2,ip,rs,cs,b,2,x,2,rcond,   &
     &                berr,n_err_bnds,err_bnds_n,err_bnds_c,nparams,    &
     &                params,w,r,info)
         CALL CHKXER('ZGBRFSX',INFot,NOUt,LERr,OK)
         INFot = 3
         eq = 'R'
         CALL ZGBRFSX('N',eq,-1,1,1,0,a,1,af,1,ip,rs,cs,b,1,x,1,rcond,  &
     &                berr,n_err_bnds,err_bnds_n,err_bnds_c,nparams,    &
     &                params,w,r,info)
         CALL CHKXER('ZGBRFSX',INFot,NOUt,LERr,OK)
         INFot = 4
         eq = 'R'
         CALL ZGBRFSX('N',eq,2,-1,1,1,a,3,af,4,ip,rs,cs,b,1,x,1,rcond,  &
     &                berr,n_err_bnds,err_bnds_n,err_bnds_c,nparams,    &
     &                params,w,r,info)
         CALL CHKXER('ZGBRFSX',INFot,NOUt,LERr,OK)
         INFot = 5
         eq = 'R'
         CALL ZGBRFSX('N',eq,2,1,-1,1,a,3,af,4,ip,rs,cs,b,1,x,1,rcond,  &
     &                berr,n_err_bnds,err_bnds_n,err_bnds_c,nparams,    &
     &                params,w,r,info)
         CALL CHKXER('ZGBRFSX',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL ZGBRFSX('N',eq,0,0,0,-1,a,1,af,1,ip,rs,cs,b,1,x,1,rcond,  &
     &                berr,n_err_bnds,err_bnds_n,err_bnds_c,nparams,    &
     &                params,w,r,info)
         CALL CHKXER('ZGBRFSX',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL ZGBRFSX('N',eq,2,1,1,1,a,1,af,2,ip,rs,cs,b,2,x,2,rcond,   &
     &                berr,n_err_bnds,err_bnds_n,err_bnds_c,nparams,    &
     &                params,w,r,info)
         CALL CHKXER('ZGBRFSX',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL ZGBRFSX('N',eq,2,1,1,1,a,3,af,3,ip,rs,cs,b,2,x,2,rcond,   &
     &                berr,n_err_bnds,err_bnds_n,err_bnds_c,nparams,    &
     &                params,w,r,info)
         CALL CHKXER('ZGBRFSX',INFot,NOUt,LERr,OK)
         INFot = 13
         eq = 'C'
         CALL ZGBRFSX('N',eq,2,1,1,1,a,3,af,5,ip,rs,cs,b,1,x,2,rcond,   &
     &                berr,n_err_bnds,err_bnds_n,err_bnds_c,nparams,    &
     &                params,w,r,info)
         CALL CHKXER('ZGBRFSX',INFot,NOUt,LERr,OK)
         INFot = 15
         CALL ZGBRFSX('N',eq,2,1,1,1,a,3,af,5,ip,rs,cs,b,2,x,1,rcond,   &
     &                berr,n_err_bnds,err_bnds_n,err_bnds_c,nparams,    &
     &                params,w,r,info)
         CALL CHKXER('ZGBRFSX',INFot,NOUt,LERr,OK)
!
!        ZGBCON
!
         SRNamt = 'ZGBCON'
         INFot = 1
         CALL ZGBCON('/',0,0,0,a,1,ip,anrm,rcond,w,r,info)
         CALL CHKXER('ZGBCON',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZGBCON('1',-1,0,0,a,1,ip,anrm,rcond,w,r,info)
         CALL CHKXER('ZGBCON',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZGBCON('1',1,-1,0,a,1,ip,anrm,rcond,w,r,info)
         CALL CHKXER('ZGBCON',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZGBCON('1',1,0,-1,a,1,ip,anrm,rcond,w,r,info)
         CALL CHKXER('ZGBCON',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL ZGBCON('1',2,1,1,a,3,ip,anrm,rcond,w,r,info)
         CALL CHKXER('ZGBCON',INFot,NOUt,LERr,OK)
!
!        ZGBEQU
!
         SRNamt = 'ZGBEQU'
         INFot = 1
         CALL ZGBEQU(-1,0,0,0,a,1,r1,r2,rcond,ccond,anrm,info)
         CALL CHKXER('ZGBEQU',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZGBEQU(0,-1,0,0,a,1,r1,r2,rcond,ccond,anrm,info)
         CALL CHKXER('ZGBEQU',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZGBEQU(1,1,-1,0,a,1,r1,r2,rcond,ccond,anrm,info)
         CALL CHKXER('ZGBEQU',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZGBEQU(1,1,0,-1,a,1,r1,r2,rcond,ccond,anrm,info)
         CALL CHKXER('ZGBEQU',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL ZGBEQU(2,2,1,1,a,2,r1,r2,rcond,ccond,anrm,info)
         CALL CHKXER('ZGBEQU',INFot,NOUt,LERr,OK)
!
!        ZGBEQUB
!
         SRNamt = 'ZGBEQUB'
         INFot = 1
         CALL ZGBEQUB(-1,0,0,0,a,1,r1,r2,rcond,ccond,anrm,info)
         CALL CHKXER('ZGBEQUB',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZGBEQUB(0,-1,0,0,a,1,r1,r2,rcond,ccond,anrm,info)
         CALL CHKXER('ZGBEQUB',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZGBEQUB(1,1,-1,0,a,1,r1,r2,rcond,ccond,anrm,info)
         CALL CHKXER('ZGBEQUB',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZGBEQUB(1,1,0,-1,a,1,r1,r2,rcond,ccond,anrm,info)
         CALL CHKXER('ZGBEQUB',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL ZGBEQUB(2,2,1,1,a,2,r1,r2,rcond,ccond,anrm,info)
         CALL CHKXER('ZGBEQUB',INFot,NOUt,LERr,OK)
      ENDIF
!
!     Print a summary line.
!
      CALL ALAESM(Path,OK,NOUt)
!
!
!     End of ZERRGE
!
      END SUBROUTINE ZERRGE
