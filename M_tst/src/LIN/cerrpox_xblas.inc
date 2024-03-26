!*==cerrpo.f90  processed by SPAG 7.51RB at 17:34 on  4 Mar 2022
!> \brief \b CERRPOX
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE cerrpo( PATH, NUNIT )
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
!> cerrpo tests the error exits for the COMPLEX routines
!> for Hermitian positive definite matrices.
!>
!> Note that this file is used only when the XBLAS are available,
!> otherwise cerrpo.f defines this subroutine.
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
      SUBROUTINE cerrpo(Path,Nunit)
      IMPLICIT NONE
!*--CERRPO62
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
      REAL anrm , rcond , berr
!     ..
!     .. Local Arrays ..
      REAL s(NMAX) , r(NMAX) , r1(NMAX) , r2(NMAX) , err_bnds_n(NMAX,3) &
     &     , err_bnds_c(NMAX,3) , params(1)
      COMPLEX a(NMAX,NMAX) , af(NMAX,NMAX) , b(NMAX) , w(2*NMAX) ,      &
     &        x(NMAX)
!     ..
!     .. External Functions ..
      LOGICAL LSAMEN
EXTERNAL    LSAMEN
!     ..
!     .. External Subroutines ..
EXTERNAL    CPBCON,   CPBEQU,   CPBRFS,   CPBTF2
EXTERNAL    CPBTRF,   CPBTRS,   CPOCON,   CPOEQU,   CPORFS,   CPOTF2
EXTERNAL    CPOTRF,   CPOTRI,   CPOTRS,   CPPCON,   CPPEQU,   CPPRFS
EXTERNAL    CPPTRF,   CPPTRI,   CPPTRS,   CPOEQUB,   CPORFSX
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
         s(j) = 0.
      ENDDO
      anrm = 1.
      OK = .TRUE.
!
!     Test error exits of the routines that use the Cholesky
!     decomposition of a Hermitian positive definite matrix.
!
      IF ( LSAMEN(2,c2,'PO') ) THEN
!
!        CPOTRF
!
         SRNamt = 'CPOTRF'
         INFot = 1
         CALL CPOTRF('/',0,a,1,info)
         CALL CHKXER('CPOTRF',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CPOTRF('U',-1,a,1,info)
         CALL CHKXER('CPOTRF',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CPOTRF('U',2,a,1,info)
         CALL CHKXER('CPOTRF',INFot,NOUt,LERr,OK)
!
!        CPOTF2
!
         SRNamt = 'CPOTF2'
         INFot = 1
         CALL CPOTF2('/',0,a,1,info)
         CALL CHKXER('CPOTF2',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CPOTF2('U',-1,a,1,info)
         CALL CHKXER('CPOTF2',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CPOTF2('U',2,a,1,info)
         CALL CHKXER('CPOTF2',INFot,NOUt,LERr,OK)
!
!        CPOTRI
!
         SRNamt = 'CPOTRI'
         INFot = 1
         CALL CPOTRI('/',0,a,1,info)
         CALL CHKXER('CPOTRI',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CPOTRI('U',-1,a,1,info)
         CALL CHKXER('CPOTRI',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CPOTRI('U',2,a,1,info)
         CALL CHKXER('CPOTRI',INFot,NOUt,LERr,OK)
!
!        CPOTRS
!
         SRNamt = 'CPOTRS'
         INFot = 1
         CALL CPOTRS('/',0,0,a,1,b,1,info)
         CALL CHKXER('CPOTRS',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CPOTRS('U',-1,0,a,1,b,1,info)
         CALL CHKXER('CPOTRS',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CPOTRS('U',0,-1,a,1,b,1,info)
         CALL CHKXER('CPOTRS',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL CPOTRS('U',2,1,a,1,b,2,info)
         CALL CHKXER('CPOTRS',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL CPOTRS('U',2,1,a,2,b,1,info)
         CALL CHKXER('CPOTRS',INFot,NOUt,LERr,OK)
!
!        CPORFS
!
         SRNamt = 'CPORFS'
         INFot = 1
         CALL CPORFS('/',0,0,a,1,af,1,b,1,x,1,r1,r2,w,r,info)
         CALL CHKXER('CPORFS',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CPORFS('U',-1,0,a,1,af,1,b,1,x,1,r1,r2,w,r,info)
         CALL CHKXER('CPORFS',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CPORFS('U',0,-1,a,1,af,1,b,1,x,1,r1,r2,w,r,info)
         CALL CHKXER('CPORFS',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL CPORFS('U',2,1,a,1,af,2,b,2,x,2,r1,r2,w,r,info)
         CALL CHKXER('CPORFS',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL CPORFS('U',2,1,a,2,af,1,b,2,x,2,r1,r2,w,r,info)
         CALL CHKXER('CPORFS',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL CPORFS('U',2,1,a,2,af,2,b,1,x,2,r1,r2,w,r,info)
         CALL CHKXER('CPORFS',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL CPORFS('U',2,1,a,2,af,2,b,2,x,1,r1,r2,w,r,info)
         CALL CHKXER('CPORFS',INFot,NOUt,LERr,OK)
!
!        CPORFSX
!
         n_err_bnds = 3
         nparams = 0
         SRNamt = 'CPORFSX'
         INFot = 1
         CALL CPORFSX('/',eq,0,0,a,1,af,1,s,b,1,x,1,rcond,berr,         &
     &                n_err_bnds,err_bnds_n,err_bnds_c,nparams,params,w,&
     &                r,info)
         CALL CHKXER('CPORFSX',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CPORFSX('U','/',-1,0,a,1,af,1,s,b,1,x,1,rcond,berr,       &
     &                n_err_bnds,err_bnds_n,err_bnds_c,nparams,params,w,&
     &                r,info)
         CALL CHKXER('CPORFSX',INFot,NOUt,LERr,OK)
         eq = 'N'
         INFot = 3
         CALL CPORFSX('U',eq,-1,0,a,1,af,1,s,b,1,x,1,rcond,berr,        &
     &                n_err_bnds,err_bnds_n,err_bnds_c,nparams,params,w,&
     &                r,info)
         CALL CHKXER('CPORFSX',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CPORFSX('U',eq,0,-1,a,1,af,1,s,b,1,x,1,rcond,berr,        &
     &                n_err_bnds,err_bnds_n,err_bnds_c,nparams,params,w,&
     &                r,info)
         CALL CHKXER('CPORFSX',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL CPORFSX('U',eq,2,1,a,1,af,2,s,b,2,x,2,rcond,berr,         &
     &                n_err_bnds,err_bnds_n,err_bnds_c,nparams,params,w,&
     &                r,info)
         CALL CHKXER('CPORFSX',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL CPORFSX('U',eq,2,1,a,2,af,1,s,b,2,x,2,rcond,berr,         &
     &                n_err_bnds,err_bnds_n,err_bnds_c,nparams,params,w,&
     &                r,info)
         CALL CHKXER('CPORFSX',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL CPORFSX('U',eq,2,1,a,2,af,2,s,b,1,x,2,rcond,berr,         &
     &                n_err_bnds,err_bnds_n,err_bnds_c,nparams,params,w,&
     &                r,info)
         CALL CHKXER('CPORFSX',INFot,NOUt,LERr,OK)
         INFot = 13
         CALL CPORFSX('U',eq,2,1,a,2,af,2,s,b,2,x,1,rcond,berr,         &
     &                n_err_bnds,err_bnds_n,err_bnds_c,nparams,params,w,&
     &                r,info)
         CALL CHKXER('CPORFSX',INFot,NOUt,LERr,OK)
!
!        CPOCON
!
         SRNamt = 'CPOCON'
         INFot = 1
         CALL CPOCON('/',0,a,1,anrm,rcond,w,r,info)
         CALL CHKXER('CPOCON',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CPOCON('U',-1,a,1,anrm,rcond,w,r,info)
         CALL CHKXER('CPOCON',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CPOCON('U',2,a,1,anrm,rcond,w,r,info)
         CALL CHKXER('CPOCON',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL CPOCON('U',1,a,1,-anrm,rcond,w,r,info)
         CALL CHKXER('CPOCON',INFot,NOUt,LERr,OK)
!
!        CPOEQU
!
         SRNamt = 'CPOEQU'
         INFot = 1
         CALL CPOEQU(-1,a,1,r1,rcond,anrm,info)
         CALL CHKXER('CPOEQU',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CPOEQU(2,a,1,r1,rcond,anrm,info)
         CALL CHKXER('CPOEQU',INFot,NOUt,LERr,OK)
!
!        CPOEQUB
!
         SRNamt = 'CPOEQUB'
         INFot = 1
         CALL CPOEQUB(-1,a,1,r1,rcond,anrm,info)
         CALL CHKXER('CPOEQUB',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CPOEQUB(2,a,1,r1,rcond,anrm,info)
         CALL CHKXER('CPOEQUB',INFot,NOUt,LERr,OK)
!
!     Test error exits of the routines that use the Cholesky
!     decomposition of a Hermitian positive definite packed matrix.
!
      ELSEIF ( LSAMEN(2,c2,'PP') ) THEN
!
!        CPPTRF
!
         SRNamt = 'CPPTRF'
         INFot = 1
         CALL CPPTRF('/',0,a,info)
         CALL CHKXER('CPPTRF',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CPPTRF('U',-1,a,info)
         CALL CHKXER('CPPTRF',INFot,NOUt,LERr,OK)
!
!        CPPTRI
!
         SRNamt = 'CPPTRI'
         INFot = 1
         CALL CPPTRI('/',0,a,info)
         CALL CHKXER('CPPTRI',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CPPTRI('U',-1,a,info)
         CALL CHKXER('CPPTRI',INFot,NOUt,LERr,OK)
!
!        CPPTRS
!
         SRNamt = 'CPPTRS'
         INFot = 1
         CALL CPPTRS('/',0,0,a,b,1,info)
         CALL CHKXER('CPPTRS',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CPPTRS('U',-1,0,a,b,1,info)
         CALL CHKXER('CPPTRS',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CPPTRS('U',0,-1,a,b,1,info)
         CALL CHKXER('CPPTRS',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL CPPTRS('U',2,1,a,b,1,info)
         CALL CHKXER('CPPTRS',INFot,NOUt,LERr,OK)
!
!        CPPRFS
!
         SRNamt = 'CPPRFS'
         INFot = 1
         CALL CPPRFS('/',0,0,a,af,b,1,x,1,r1,r2,w,r,info)
         CALL CHKXER('CPPRFS',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CPPRFS('U',-1,0,a,af,b,1,x,1,r1,r2,w,r,info)
         CALL CHKXER('CPPRFS',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CPPRFS('U',0,-1,a,af,b,1,x,1,r1,r2,w,r,info)
         CALL CHKXER('CPPRFS',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL CPPRFS('U',2,1,a,af,b,1,x,2,r1,r2,w,r,info)
         CALL CHKXER('CPPRFS',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL CPPRFS('U',2,1,a,af,b,2,x,1,r1,r2,w,r,info)
         CALL CHKXER('CPPRFS',INFot,NOUt,LERr,OK)
!
!        CPPCON
!
         SRNamt = 'CPPCON'
         INFot = 1
         CALL CPPCON('/',0,a,anrm,rcond,w,r,info)
         CALL CHKXER('CPPCON',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CPPCON('U',-1,a,anrm,rcond,w,r,info)
         CALL CHKXER('CPPCON',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CPPCON('U',1,a,-anrm,rcond,w,r,info)
         CALL CHKXER('CPPCON',INFot,NOUt,LERr,OK)
!
!        CPPEQU
!
         SRNamt = 'CPPEQU'
         INFot = 1
         CALL CPPEQU('/',0,a,r1,rcond,anrm,info)
         CALL CHKXER('CPPEQU',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CPPEQU('U',-1,a,r1,rcond,anrm,info)
         CALL CHKXER('CPPEQU',INFot,NOUt,LERr,OK)
!
!     Test error exits of the routines that use the Cholesky
!     decomposition of a Hermitian positive definite band matrix.
!
      ELSEIF ( LSAMEN(2,c2,'PB') ) THEN
!
!        CPBTRF
!
         SRNamt = 'CPBTRF'
         INFot = 1
         CALL CPBTRF('/',0,0,a,1,info)
         CALL CHKXER('CPBTRF',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CPBTRF('U',-1,0,a,1,info)
         CALL CHKXER('CPBTRF',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CPBTRF('U',1,-1,a,1,info)
         CALL CHKXER('CPBTRF',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL CPBTRF('U',2,1,a,1,info)
         CALL CHKXER('CPBTRF',INFot,NOUt,LERr,OK)
!
!        CPBTF2
!
         SRNamt = 'CPBTF2'
         INFot = 1
         CALL CPBTF2('/',0,0,a,1,info)
         CALL CHKXER('CPBTF2',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CPBTF2('U',-1,0,a,1,info)
         CALL CHKXER('CPBTF2',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CPBTF2('U',1,-1,a,1,info)
         CALL CHKXER('CPBTF2',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL CPBTF2('U',2,1,a,1,info)
         CALL CHKXER('CPBTF2',INFot,NOUt,LERr,OK)
!
!        CPBTRS
!
         SRNamt = 'CPBTRS'
         INFot = 1
         CALL CPBTRS('/',0,0,0,a,1,b,1,info)
         CALL CHKXER('CPBTRS',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CPBTRS('U',-1,0,0,a,1,b,1,info)
         CALL CHKXER('CPBTRS',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CPBTRS('U',1,-1,0,a,1,b,1,info)
         CALL CHKXER('CPBTRS',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CPBTRS('U',0,0,-1,a,1,b,1,info)
         CALL CHKXER('CPBTRS',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL CPBTRS('U',2,1,1,a,1,b,1,info)
         CALL CHKXER('CPBTRS',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL CPBTRS('U',2,0,1,a,1,b,1,info)
         CALL CHKXER('CPBTRS',INFot,NOUt,LERr,OK)
!
!        CPBRFS
!
         SRNamt = 'CPBRFS'
         INFot = 1
         CALL CPBRFS('/',0,0,0,a,1,af,1,b,1,x,1,r1,r2,w,r,info)
         CALL CHKXER('CPBRFS',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CPBRFS('U',-1,0,0,a,1,af,1,b,1,x,1,r1,r2,w,r,info)
         CALL CHKXER('CPBRFS',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CPBRFS('U',1,-1,0,a,1,af,1,b,1,x,1,r1,r2,w,r,info)
         CALL CHKXER('CPBRFS',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CPBRFS('U',0,0,-1,a,1,af,1,b,1,x,1,r1,r2,w,r,info)
         CALL CHKXER('CPBRFS',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL CPBRFS('U',2,1,1,a,1,af,2,b,2,x,2,r1,r2,w,r,info)
         CALL CHKXER('CPBRFS',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL CPBRFS('U',2,1,1,a,2,af,1,b,2,x,2,r1,r2,w,r,info)
         CALL CHKXER('CPBRFS',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL CPBRFS('U',2,0,1,a,1,af,1,b,1,x,2,r1,r2,w,r,info)
         CALL CHKXER('CPBRFS',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL CPBRFS('U',2,0,1,a,1,af,1,b,2,x,1,r1,r2,w,r,info)
         CALL CHKXER('CPBRFS',INFot,NOUt,LERr,OK)
!
!        CPBCON
!
         SRNamt = 'CPBCON'
         INFot = 1
         CALL CPBCON('/',0,0,a,1,anrm,rcond,w,r,info)
         CALL CHKXER('CPBCON',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CPBCON('U',-1,0,a,1,anrm,rcond,w,r,info)
         CALL CHKXER('CPBCON',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CPBCON('U',1,-1,a,1,anrm,rcond,w,r,info)
         CALL CHKXER('CPBCON',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL CPBCON('U',2,1,a,1,anrm,rcond,w,r,info)
         CALL CHKXER('CPBCON',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL CPBCON('U',1,0,a,1,-anrm,rcond,w,r,info)
         CALL CHKXER('CPBCON',INFot,NOUt,LERr,OK)
!
!        CPBEQU
!
         SRNamt = 'CPBEQU'
         INFot = 1
         CALL CPBEQU('/',0,0,a,1,r1,rcond,anrm,info)
         CALL CHKXER('CPBEQU',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CPBEQU('U',-1,0,a,1,r1,rcond,anrm,info)
         CALL CHKXER('CPBEQU',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CPBEQU('U',1,-1,a,1,r1,rcond,anrm,info)
         CALL CHKXER('CPBEQU',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL CPBEQU('U',2,1,a,1,r1,rcond,anrm,info)
         CALL CHKXER('CPBEQU',INFot,NOUt,LERr,OK)
      ENDIF
!
!     Print a summary line.
!
      CALL alaesm(Path,OK,NOUt)
!
!
!     End of cerrpo
!
      END SUBROUTINE cerrpo
