!*==zerrpo.f90  processed by SPAG 7.51RB at 17:35 on  4 Mar 2022
!> \brief \b ZERRPOX
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE zerrpo( PATH, NUNIT )
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
!> zerrpo tests the error exits for the COMPLEX*16 routines
!> for Hermitian positive definite matrices.
!>
!> Note that this file is used only when the XBLAS are available,
!> otherwise zerrpo.f defines this subroutine.
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
      SUBROUTINE zerrpo(Path,Nunit)
      IMPLICIT NONE
!*--ZERRPO62
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
      DOUBLE PRECISION anrm , rcond , berr
!     ..
!     .. Local Arrays ..
      DOUBLE PRECISION s(NMAX) , r(NMAX) , r1(NMAX) , r2(NMAX) ,        &
     &                 err_bnds_n(NMAX,3) , err_bnds_c(NMAX,3) ,        &
     &                 params(1)
      COMPLEX*16 a(NMAX,NMAX) , af(NMAX,NMAX) , b(NMAX) , w(2*NMAX) ,   &
     &           x(NMAX)
!     ..
!     .. External Functions ..
      LOGICAL LSAMEN
EXTERNAL    LSAMEN
!     ..
!     .. External Subroutines ..
EXTERNAL    CHKXER,   ZPBCON,   ZPBEQU,   ZPBRFS,   ZPBTF2
EXTERNAL    ZPBTRF,   ZPBTRS,   ZPOCON,   ZPOEQU,   ZPORFS,   ZPOTF2
EXTERNAL    ZPOTRF,   ZPOTRI,   ZPOTRS,   ZPPCON,   ZPPEQU,   ZPPRFS
EXTERNAL    ZPPTRF,   ZPPTRI,   ZPPTRS,   ZPOEQUB,   ZPORFSX
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
         s(j) = 0.D0
      ENDDO
      anrm = 1.D0
      OK = .TRUE.
!
!     Test error exits of the routines that use the Cholesky
!     decomposition of a Hermitian positive definite matrix.
!
      IF ( LSAMEN(2,c2,'PO') ) THEN
!
!        ZPOTRF
!
         SRNamt = 'ZPOTRF'
         INFot = 1
         CALL ZPOTRF('/',0,a,1,info)
         CALL CHKXER('ZPOTRF',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZPOTRF('U',-1,a,1,info)
         CALL CHKXER('ZPOTRF',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZPOTRF('U',2,a,1,info)
         CALL CHKXER('ZPOTRF',INFot,NOUt,LERr,OK)
!
!        ZPOTF2
!
         SRNamt = 'ZPOTF2'
         INFot = 1
         CALL ZPOTF2('/',0,a,1,info)
         CALL CHKXER('ZPOTF2',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZPOTF2('U',-1,a,1,info)
         CALL CHKXER('ZPOTF2',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZPOTF2('U',2,a,1,info)
         CALL CHKXER('ZPOTF2',INFot,NOUt,LERr,OK)
!
!        ZPOTRI
!
         SRNamt = 'ZPOTRI'
         INFot = 1
         CALL ZPOTRI('/',0,a,1,info)
         CALL CHKXER('ZPOTRI',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZPOTRI('U',-1,a,1,info)
         CALL CHKXER('ZPOTRI',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZPOTRI('U',2,a,1,info)
         CALL CHKXER('ZPOTRI',INFot,NOUt,LERr,OK)
!
!        ZPOTRS
!
         SRNamt = 'ZPOTRS'
         INFot = 1
         CALL ZPOTRS('/',0,0,a,1,b,1,info)
         CALL CHKXER('ZPOTRS',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZPOTRS('U',-1,0,a,1,b,1,info)
         CALL CHKXER('ZPOTRS',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZPOTRS('U',0,-1,a,1,b,1,info)
         CALL CHKXER('ZPOTRS',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL ZPOTRS('U',2,1,a,1,b,2,info)
         CALL CHKXER('ZPOTRS',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL ZPOTRS('U',2,1,a,2,b,1,info)
         CALL CHKXER('ZPOTRS',INFot,NOUt,LERr,OK)
!
!        ZPORFS
!
         SRNamt = 'ZPORFS'
         INFot = 1
         CALL ZPORFS('/',0,0,a,1,af,1,b,1,x,1,r1,r2,w,r,info)
         CALL CHKXER('ZPORFS',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZPORFS('U',-1,0,a,1,af,1,b,1,x,1,r1,r2,w,r,info)
         CALL CHKXER('ZPORFS',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZPORFS('U',0,-1,a,1,af,1,b,1,x,1,r1,r2,w,r,info)
         CALL CHKXER('ZPORFS',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL ZPORFS('U',2,1,a,1,af,2,b,2,x,2,r1,r2,w,r,info)
         CALL CHKXER('ZPORFS',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL ZPORFS('U',2,1,a,2,af,1,b,2,x,2,r1,r2,w,r,info)
         CALL CHKXER('ZPORFS',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL ZPORFS('U',2,1,a,2,af,2,b,1,x,2,r1,r2,w,r,info)
         CALL CHKXER('ZPORFS',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL ZPORFS('U',2,1,a,2,af,2,b,2,x,1,r1,r2,w,r,info)
         CALL CHKXER('ZPORFS',INFot,NOUt,LERr,OK)
!
!        ZPORFSX
!
         n_err_bnds = 3
         nparams = 0
         SRNamt = 'ZPORFSX'
         INFot = 1
         CALL ZPORFSX('/',eq,0,0,a,1,af,1,s,b,1,x,1,rcond,berr,         &
     &                n_err_bnds,err_bnds_n,err_bnds_c,nparams,params,w,&
     &                r,info)
         CALL CHKXER('ZPORFSX',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZPORFSX('U',"/",-1,0,a,1,af,1,s,b,1,x,1,rcond,berr,       &
     &                n_err_bnds,err_bnds_n,err_bnds_c,nparams,params,w,&
     &                r,info)
         CALL CHKXER('ZPORFSX',INFot,NOUt,LERr,OK)
         eq = 'N'
         INFot = 3
         CALL ZPORFSX('U',eq,-1,0,a,1,af,1,s,b,1,x,1,rcond,berr,        &
     &                n_err_bnds,err_bnds_n,err_bnds_c,nparams,params,w,&
     &                r,info)
         CALL CHKXER('ZPORFSX',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZPORFSX('U',eq,0,-1,a,1,af,1,s,b,1,x,1,rcond,berr,        &
     &                n_err_bnds,err_bnds_n,err_bnds_c,nparams,params,w,&
     &                r,info)
         CALL CHKXER('ZPORFSX',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL ZPORFSX('U',eq,2,1,a,1,af,2,s,b,2,x,2,rcond,berr,         &
     &                n_err_bnds,err_bnds_n,err_bnds_c,nparams,params,w,&
     &                r,info)
         CALL CHKXER('ZPORFSX',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL ZPORFSX('U',eq,2,1,a,2,af,1,s,b,2,x,2,rcond,berr,         &
     &                n_err_bnds,err_bnds_n,err_bnds_c,nparams,params,w,&
     &                r,info)
         CALL CHKXER('ZPORFSX',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL ZPORFSX('U',eq,2,1,a,2,af,2,s,b,1,x,2,rcond,berr,         &
     &                n_err_bnds,err_bnds_n,err_bnds_c,nparams,params,w,&
     &                r,info)
         CALL CHKXER('ZPORFSX',INFot,NOUt,LERr,OK)
         INFot = 13
         CALL ZPORFSX('U',eq,2,1,a,2,af,2,s,b,2,x,1,rcond,berr,         &
     &                n_err_bnds,err_bnds_n,err_bnds_c,nparams,params,w,&
     &                r,info)
         CALL CHKXER('ZPORFSX',INFot,NOUt,LERr,OK)
!
!        ZPOCON
!
         SRNamt = 'ZPOCON'
         INFot = 1
         CALL ZPOCON('/',0,a,1,anrm,rcond,w,r,info)
         CALL CHKXER('ZPOCON',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZPOCON('U',-1,a,1,anrm,rcond,w,r,info)
         CALL CHKXER('ZPOCON',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZPOCON('U',2,a,1,anrm,rcond,w,r,info)
         CALL CHKXER('ZPOCON',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL ZPOCON('U',1,a,1,-anrm,rcond,w,r,info)
         CALL CHKXER('ZPOCON',INFot,NOUt,LERr,OK)
!
!        ZPOEQU
!
         SRNamt = 'ZPOEQU'
         INFot = 1
         CALL ZPOEQU(-1,a,1,r1,rcond,anrm,info)
         CALL CHKXER('ZPOEQU',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZPOEQU(2,a,1,r1,rcond,anrm,info)
         CALL CHKXER('ZPOEQU',INFot,NOUt,LERr,OK)
!
!        ZPOEQUB
!
         SRNamt = 'ZPOEQUB'
         INFot = 1
         CALL ZPOEQUB(-1,a,1,r1,rcond,anrm,info)
         CALL CHKXER('ZPOEQUB',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZPOEQUB(2,a,1,r1,rcond,anrm,info)
         CALL CHKXER('ZPOEQUB',INFot,NOUt,LERr,OK)
!
!     Test error exits of the routines that use the Cholesky
!     decomposition of a Hermitian positive definite packed matrix.
!
      ELSEIF ( LSAMEN(2,c2,'PP') ) THEN
!
!        ZPPTRF
!
         SRNamt = 'ZPPTRF'
         INFot = 1
         CALL ZPPTRF('/',0,a,info)
         CALL CHKXER('ZPPTRF',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZPPTRF('U',-1,a,info)
         CALL CHKXER('ZPPTRF',INFot,NOUt,LERr,OK)
!
!        ZPPTRI
!
         SRNamt = 'ZPPTRI'
         INFot = 1
         CALL ZPPTRI('/',0,a,info)
         CALL CHKXER('ZPPTRI',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZPPTRI('U',-1,a,info)
         CALL CHKXER('ZPPTRI',INFot,NOUt,LERr,OK)
!
!        ZPPTRS
!
         SRNamt = 'ZPPTRS'
         INFot = 1
         CALL ZPPTRS('/',0,0,a,b,1,info)
         CALL CHKXER('ZPPTRS',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZPPTRS('U',-1,0,a,b,1,info)
         CALL CHKXER('ZPPTRS',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZPPTRS('U',0,-1,a,b,1,info)
         CALL CHKXER('ZPPTRS',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL ZPPTRS('U',2,1,a,b,1,info)
         CALL CHKXER('ZPPTRS',INFot,NOUt,LERr,OK)
!
!        ZPPRFS
!
         SRNamt = 'ZPPRFS'
         INFot = 1
         CALL ZPPRFS('/',0,0,a,af,b,1,x,1,r1,r2,w,r,info)
         CALL CHKXER('ZPPRFS',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZPPRFS('U',-1,0,a,af,b,1,x,1,r1,r2,w,r,info)
         CALL CHKXER('ZPPRFS',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZPPRFS('U',0,-1,a,af,b,1,x,1,r1,r2,w,r,info)
         CALL CHKXER('ZPPRFS',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL ZPPRFS('U',2,1,a,af,b,1,x,2,r1,r2,w,r,info)
         CALL CHKXER('ZPPRFS',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL ZPPRFS('U',2,1,a,af,b,2,x,1,r1,r2,w,r,info)
         CALL CHKXER('ZPPRFS',INFot,NOUt,LERr,OK)
!
!        ZPPCON
!
         SRNamt = 'ZPPCON'
         INFot = 1
         CALL ZPPCON('/',0,a,anrm,rcond,w,r,info)
         CALL CHKXER('ZPPCON',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZPPCON('U',-1,a,anrm,rcond,w,r,info)
         CALL CHKXER('ZPPCON',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZPPCON('U',1,a,-anrm,rcond,w,r,info)
         CALL CHKXER('ZPPCON',INFot,NOUt,LERr,OK)
!
!        ZPPEQU
!
         SRNamt = 'ZPPEQU'
         INFot = 1
         CALL ZPPEQU('/',0,a,r1,rcond,anrm,info)
         CALL CHKXER('ZPPEQU',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZPPEQU('U',-1,a,r1,rcond,anrm,info)
         CALL CHKXER('ZPPEQU',INFot,NOUt,LERr,OK)
!
!     Test error exits of the routines that use the Cholesky
!     decomposition of a Hermitian positive definite band matrix.
!
      ELSEIF ( LSAMEN(2,c2,'PB') ) THEN
!
!        ZPBTRF
!
         SRNamt = 'ZPBTRF'
         INFot = 1
         CALL ZPBTRF('/',0,0,a,1,info)
         CALL CHKXER('ZPBTRF',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZPBTRF('U',-1,0,a,1,info)
         CALL CHKXER('ZPBTRF',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZPBTRF('U',1,-1,a,1,info)
         CALL CHKXER('ZPBTRF',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL ZPBTRF('U',2,1,a,1,info)
         CALL CHKXER('ZPBTRF',INFot,NOUt,LERr,OK)
!
!        ZPBTF2
!
         SRNamt = 'ZPBTF2'
         INFot = 1
         CALL ZPBTF2('/',0,0,a,1,info)
         CALL CHKXER('ZPBTF2',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZPBTF2('U',-1,0,a,1,info)
         CALL CHKXER('ZPBTF2',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZPBTF2('U',1,-1,a,1,info)
         CALL CHKXER('ZPBTF2',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL ZPBTF2('U',2,1,a,1,info)
         CALL CHKXER('ZPBTF2',INFot,NOUt,LERr,OK)
!
!        ZPBTRS
!
         SRNamt = 'ZPBTRS'
         INFot = 1
         CALL ZPBTRS('/',0,0,0,a,1,b,1,info)
         CALL CHKXER('ZPBTRS',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZPBTRS('U',-1,0,0,a,1,b,1,info)
         CALL CHKXER('ZPBTRS',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZPBTRS('U',1,-1,0,a,1,b,1,info)
         CALL CHKXER('ZPBTRS',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZPBTRS('U',0,0,-1,a,1,b,1,info)
         CALL CHKXER('ZPBTRS',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL ZPBTRS('U',2,1,1,a,1,b,1,info)
         CALL CHKXER('ZPBTRS',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL ZPBTRS('U',2,0,1,a,1,b,1,info)
         CALL CHKXER('ZPBTRS',INFot,NOUt,LERr,OK)
!
!        ZPBRFS
!
         SRNamt = 'ZPBRFS'
         INFot = 1
         CALL ZPBRFS('/',0,0,0,a,1,af,1,b,1,x,1,r1,r2,w,r,info)
         CALL CHKXER('ZPBRFS',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZPBRFS('U',-1,0,0,a,1,af,1,b,1,x,1,r1,r2,w,r,info)
         CALL CHKXER('ZPBRFS',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZPBRFS('U',1,-1,0,a,1,af,1,b,1,x,1,r1,r2,w,r,info)
         CALL CHKXER('ZPBRFS',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZPBRFS('U',0,0,-1,a,1,af,1,b,1,x,1,r1,r2,w,r,info)
         CALL CHKXER('ZPBRFS',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL ZPBRFS('U',2,1,1,a,1,af,2,b,2,x,2,r1,r2,w,r,info)
         CALL CHKXER('ZPBRFS',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL ZPBRFS('U',2,1,1,a,2,af,1,b,2,x,2,r1,r2,w,r,info)
         CALL CHKXER('ZPBRFS',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL ZPBRFS('U',2,0,1,a,1,af,1,b,1,x,2,r1,r2,w,r,info)
         CALL CHKXER('ZPBRFS',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL ZPBRFS('U',2,0,1,a,1,af,1,b,2,x,1,r1,r2,w,r,info)
         CALL CHKXER('ZPBRFS',INFot,NOUt,LERr,OK)
!
!        ZPBCON
!
         SRNamt = 'ZPBCON'
         INFot = 1
         CALL ZPBCON('/',0,0,a,1,anrm,rcond,w,r,info)
         CALL CHKXER('ZPBCON',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZPBCON('U',-1,0,a,1,anrm,rcond,w,r,info)
         CALL CHKXER('ZPBCON',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZPBCON('U',1,-1,a,1,anrm,rcond,w,r,info)
         CALL CHKXER('ZPBCON',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL ZPBCON('U',2,1,a,1,anrm,rcond,w,r,info)
         CALL CHKXER('ZPBCON',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL ZPBCON('U',1,0,a,1,-anrm,rcond,w,r,info)
         CALL CHKXER('ZPBCON',INFot,NOUt,LERr,OK)
!
!        ZPBEQU
!
         SRNamt = 'ZPBEQU'
         INFot = 1
         CALL ZPBEQU('/',0,0,a,1,r1,rcond,anrm,info)
         CALL CHKXER('ZPBEQU',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZPBEQU('U',-1,0,a,1,r1,rcond,anrm,info)
         CALL CHKXER('ZPBEQU',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZPBEQU('U',1,-1,a,1,r1,rcond,anrm,info)
         CALL CHKXER('ZPBEQU',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL ZPBEQU('U',2,1,a,1,r1,rcond,anrm,info)
         CALL CHKXER('ZPBEQU',INFot,NOUt,LERr,OK)
      ENDIF
!
!     Print a summary line.
!
      CALL alaesm(Path,OK,NOUt)
!
!
!     End of zerrpo
!
      END SUBROUTINE zerrpo
