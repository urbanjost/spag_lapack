!*==derrpo.f90  processed by SPAG 7.51RB at 17:34 on  4 Mar 2022
!> \brief \b DERRPOX
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE derrpo( PATH, NUNIT )
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
!> derrpo tests the error exits for the DOUBLE PRECISION routines
!> for symmetric positive definite matrices.
!>
!> Note that this file is used only when the XBLAS are available,
!> otherwise derrpo.f defines this subroutine.
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
      SUBROUTINE derrpo(Path,Nunit)
      IMPLICIT NONE
!*--DERRPO62
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
      INTEGER iw(NMAX)
      DOUBLE PRECISION a(NMAX,NMAX) , af(NMAX,NMAX) , b(NMAX) , r1(NMAX)&
     &                 , r2(NMAX) , w(3*NMAX) , x(NMAX) , s(NMAX) ,     &
     &                 err_bnds_n(NMAX,3) , err_bnds_c(NMAX,3) ,        &
     &                 params(1)
!     ..
!     .. External Functions ..
      LOGICAL LSAMEN
EXTERNAL    LSAMEN
!     ..
!     .. External Subroutines ..
EXTERNAL    DPBCON,   DPBEQU,   DPBRFS,   DPBTF2
EXTERNAL    DPBTRF,   DPBTRS,   DPOCON,   DPOEQU,   DPORFS,   DPOTF2
EXTERNAL    DPOTRF,   DPOTRI,   DPOTRS,   DPPCON,   DPPEQU,   DPPRFS
EXTERNAL    DPPTRF,   DPPTRI,   DPPTRS,   DPOEQUB,   DPORFSX
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
         s(j) = 0.D0
         iw(j) = j
      ENDDO
      OK = .TRUE.
!
      IF ( LSAMEN(2,c2,'PO') ) THEN
!
!        Test error exits of the routines that use the Cholesky
!        decomposition of a symmetric positive definite matrix.
!
!        DPOTRF
!
         SRNamt = 'DPOTRF'
         INFot = 1
         CALL DPOTRF('/',0,a,1,info)
         CALL CHKXER('DPOTRF',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DPOTRF('U',-1,a,1,info)
         CALL CHKXER('DPOTRF',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL DPOTRF('U',2,a,1,info)
         CALL CHKXER('DPOTRF',INFot,NOUt,LERr,OK)
!
!        DPOTF2
!
         SRNamt = 'DPOTF2'
         INFot = 1
         CALL DPOTF2('/',0,a,1,info)
         CALL CHKXER('DPOTF2',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DPOTF2('U',-1,a,1,info)
         CALL CHKXER('DPOTF2',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL DPOTF2('U',2,a,1,info)
         CALL CHKXER('DPOTF2',INFot,NOUt,LERr,OK)
!
!        DPOTRI
!
         SRNamt = 'DPOTRI'
         INFot = 1
         CALL DPOTRI('/',0,a,1,info)
         CALL CHKXER('DPOTRI',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DPOTRI('U',-1,a,1,info)
         CALL CHKXER('DPOTRI',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL DPOTRI('U',2,a,1,info)
         CALL CHKXER('DPOTRI',INFot,NOUt,LERr,OK)
!
!        DPOTRS
!
         SRNamt = 'DPOTRS'
         INFot = 1
         CALL DPOTRS('/',0,0,a,1,b,1,info)
         CALL CHKXER('DPOTRS',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DPOTRS('U',-1,0,a,1,b,1,info)
         CALL CHKXER('DPOTRS',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DPOTRS('U',0,-1,a,1,b,1,info)
         CALL CHKXER('DPOTRS',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL DPOTRS('U',2,1,a,1,b,2,info)
         CALL CHKXER('DPOTRS',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL DPOTRS('U',2,1,a,2,b,1,info)
         CALL CHKXER('DPOTRS',INFot,NOUt,LERr,OK)
!
!        DPORFS
!
         SRNamt = 'DPORFS'
         INFot = 1
         CALL DPORFS('/',0,0,a,1,af,1,b,1,x,1,r1,r2,w,iw,info)
         CALL CHKXER('DPORFS',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DPORFS('U',-1,0,a,1,af,1,b,1,x,1,r1,r2,w,iw,info)
         CALL CHKXER('DPORFS',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DPORFS('U',0,-1,a,1,af,1,b,1,x,1,r1,r2,w,iw,info)
         CALL CHKXER('DPORFS',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL DPORFS('U',2,1,a,1,af,2,b,2,x,2,r1,r2,w,iw,info)
         CALL CHKXER('DPORFS',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL DPORFS('U',2,1,a,2,af,1,b,2,x,2,r1,r2,w,iw,info)
         CALL CHKXER('DPORFS',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL DPORFS('U',2,1,a,2,af,2,b,1,x,2,r1,r2,w,iw,info)
         CALL CHKXER('DPORFS',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL DPORFS('U',2,1,a,2,af,2,b,2,x,1,r1,r2,w,iw,info)
         CALL CHKXER('DPORFS',INFot,NOUt,LERr,OK)
!
!        DPORFSX
!
         n_err_bnds = 3
         nparams = 0
         SRNamt = 'DPORFSX'
         INFot = 1
         CALL DPORFSX('/',eq,0,0,a,1,af,1,s,b,1,x,1,rcond,berr,         &
     &                n_err_bnds,err_bnds_n,err_bnds_c,nparams,params,w,&
     &                iw,info)
         CALL CHKXER('DPORFSX',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DPORFSX('U',"/",-1,0,a,1,af,1,s,b,1,x,1,rcond,berr,       &
     &                n_err_bnds,err_bnds_n,err_bnds_c,nparams,params,w,&
     &                iw,info)
         CALL CHKXER('DPORFSX',INFot,NOUt,LERr,OK)
         eq = 'N'
         INFot = 3
         CALL DPORFSX('U',eq,-1,0,a,1,af,1,s,b,1,x,1,rcond,berr,        &
     &                n_err_bnds,err_bnds_n,err_bnds_c,nparams,params,w,&
     &                iw,info)
         CALL CHKXER('DPORFSX',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL DPORFSX('U',eq,0,-1,a,1,af,1,s,b,1,x,1,rcond,berr,        &
     &                n_err_bnds,err_bnds_n,err_bnds_c,nparams,params,w,&
     &                iw,info)
         CALL CHKXER('DPORFSX',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL DPORFSX('U',eq,2,1,a,1,af,2,s,b,2,x,2,rcond,berr,         &
     &                n_err_bnds,err_bnds_n,err_bnds_c,nparams,params,w,&
     &                iw,info)
         CALL CHKXER('DPORFSX',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL DPORFSX('U',eq,2,1,a,2,af,1,s,b,2,x,2,rcond,berr,         &
     &                n_err_bnds,err_bnds_n,err_bnds_c,nparams,params,w,&
     &                iw,info)
         CALL CHKXER('DPORFSX',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL DPORFSX('U',eq,2,1,a,2,af,2,s,b,1,x,2,rcond,berr,         &
     &                n_err_bnds,err_bnds_n,err_bnds_c,nparams,params,w,&
     &                iw,info)
         CALL CHKXER('DPORFSX',INFot,NOUt,LERr,OK)
         INFot = 13
         CALL DPORFSX('U',eq,2,1,a,2,af,2,s,b,2,x,1,rcond,berr,         &
     &                n_err_bnds,err_bnds_n,err_bnds_c,nparams,params,w,&
     &                iw,info)
         CALL CHKXER('DPORFSX',INFot,NOUt,LERr,OK)
!
!        DPOCON
!
         SRNamt = 'DPOCON'
         INFot = 1
         CALL DPOCON('/',0,a,1,anrm,rcond,w,iw,info)
         CALL CHKXER('DPOCON',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DPOCON('U',-1,a,1,anrm,rcond,w,iw,info)
         CALL CHKXER('DPOCON',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL DPOCON('U',2,a,1,anrm,rcond,w,iw,info)
         CALL CHKXER('DPOCON',INFot,NOUt,LERr,OK)
!
!        DPOEQU
!
         SRNamt = 'DPOEQU'
         INFot = 1
         CALL DPOEQU(-1,a,1,r1,rcond,anrm,info)
         CALL CHKXER('DPOEQU',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DPOEQU(2,a,1,r1,rcond,anrm,info)
         CALL CHKXER('DPOEQU',INFot,NOUt,LERr,OK)
!
!        DPOEQUB
!
         SRNamt = 'DPOEQUB'
         INFot = 1
         CALL DPOEQUB(-1,a,1,r1,rcond,anrm,info)
         CALL CHKXER('DPOEQUB',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DPOEQUB(2,a,1,r1,rcond,anrm,info)
         CALL CHKXER('DPOEQUB',INFot,NOUt,LERr,OK)
!
      ELSEIF ( LSAMEN(2,c2,'PP') ) THEN
!
!        Test error exits of the routines that use the Cholesky
!        decomposition of a symmetric positive definite packed matrix.
!
!        DPPTRF
!
         SRNamt = 'DPPTRF'
         INFot = 1
         CALL DPPTRF('/',0,a,info)
         CALL CHKXER('DPPTRF',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DPPTRF('U',-1,a,info)
         CALL CHKXER('DPPTRF',INFot,NOUt,LERr,OK)
!
!        DPPTRI
!
         SRNamt = 'DPPTRI'
         INFot = 1
         CALL DPPTRI('/',0,a,info)
         CALL CHKXER('DPPTRI',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DPPTRI('U',-1,a,info)
         CALL CHKXER('DPPTRI',INFot,NOUt,LERr,OK)
!
!        DPPTRS
!
         SRNamt = 'DPPTRS'
         INFot = 1
         CALL DPPTRS('/',0,0,a,b,1,info)
         CALL CHKXER('DPPTRS',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DPPTRS('U',-1,0,a,b,1,info)
         CALL CHKXER('DPPTRS',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DPPTRS('U',0,-1,a,b,1,info)
         CALL CHKXER('DPPTRS',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL DPPTRS('U',2,1,a,b,1,info)
         CALL CHKXER('DPPTRS',INFot,NOUt,LERr,OK)
!
!        DPPRFS
!
         SRNamt = 'DPPRFS'
         INFot = 1
         CALL DPPRFS('/',0,0,a,af,b,1,x,1,r1,r2,w,iw,info)
         CALL CHKXER('DPPRFS',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DPPRFS('U',-1,0,a,af,b,1,x,1,r1,r2,w,iw,info)
         CALL CHKXER('DPPRFS',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DPPRFS('U',0,-1,a,af,b,1,x,1,r1,r2,w,iw,info)
         CALL CHKXER('DPPRFS',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL DPPRFS('U',2,1,a,af,b,1,x,2,r1,r2,w,iw,info)
         CALL CHKXER('DPPRFS',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL DPPRFS('U',2,1,a,af,b,2,x,1,r1,r2,w,iw,info)
         CALL CHKXER('DPPRFS',INFot,NOUt,LERr,OK)
!
!        DPPCON
!
         SRNamt = 'DPPCON'
         INFot = 1
         CALL DPPCON('/',0,a,anrm,rcond,w,iw,info)
         CALL CHKXER('DPPCON',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DPPCON('U',-1,a,anrm,rcond,w,iw,info)
         CALL CHKXER('DPPCON',INFot,NOUt,LERr,OK)
!
!        DPPEQU
!
         SRNamt = 'DPPEQU'
         INFot = 1
         CALL DPPEQU('/',0,a,r1,rcond,anrm,info)
         CALL CHKXER('DPPEQU',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DPPEQU('U',-1,a,r1,rcond,anrm,info)
         CALL CHKXER('DPPEQU',INFot,NOUt,LERr,OK)
!
      ELSEIF ( LSAMEN(2,c2,'PB') ) THEN
!
!        Test error exits of the routines that use the Cholesky
!        decomposition of a symmetric positive definite band matrix.
!
!        DPBTRF
!
         SRNamt = 'DPBTRF'
         INFot = 1
         CALL DPBTRF('/',0,0,a,1,info)
         CALL CHKXER('DPBTRF',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DPBTRF('U',-1,0,a,1,info)
         CALL CHKXER('DPBTRF',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DPBTRF('U',1,-1,a,1,info)
         CALL CHKXER('DPBTRF',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL DPBTRF('U',2,1,a,1,info)
         CALL CHKXER('DPBTRF',INFot,NOUt,LERr,OK)
!
!        DPBTF2
!
         SRNamt = 'DPBTF2'
         INFot = 1
         CALL DPBTF2('/',0,0,a,1,info)
         CALL CHKXER('DPBTF2',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DPBTF2('U',-1,0,a,1,info)
         CALL CHKXER('DPBTF2',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DPBTF2('U',1,-1,a,1,info)
         CALL CHKXER('DPBTF2',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL DPBTF2('U',2,1,a,1,info)
         CALL CHKXER('DPBTF2',INFot,NOUt,LERr,OK)
!
!        DPBTRS
!
         SRNamt = 'DPBTRS'
         INFot = 1
         CALL DPBTRS('/',0,0,0,a,1,b,1,info)
         CALL CHKXER('DPBTRS',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DPBTRS('U',-1,0,0,a,1,b,1,info)
         CALL CHKXER('DPBTRS',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DPBTRS('U',1,-1,0,a,1,b,1,info)
         CALL CHKXER('DPBTRS',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL DPBTRS('U',0,0,-1,a,1,b,1,info)
         CALL CHKXER('DPBTRS',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL DPBTRS('U',2,1,1,a,1,b,1,info)
         CALL CHKXER('DPBTRS',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL DPBTRS('U',2,0,1,a,1,b,1,info)
         CALL CHKXER('DPBTRS',INFot,NOUt,LERr,OK)
!
!        DPBRFS
!
         SRNamt = 'DPBRFS'
         INFot = 1
         CALL DPBRFS('/',0,0,0,a,1,af,1,b,1,x,1,r1,r2,w,iw,info)
         CALL CHKXER('DPBRFS',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DPBRFS('U',-1,0,0,a,1,af,1,b,1,x,1,r1,r2,w,iw,info)
         CALL CHKXER('DPBRFS',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DPBRFS('U',1,-1,0,a,1,af,1,b,1,x,1,r1,r2,w,iw,info)
         CALL CHKXER('DPBRFS',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL DPBRFS('U',0,0,-1,a,1,af,1,b,1,x,1,r1,r2,w,iw,info)
         CALL CHKXER('DPBRFS',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL DPBRFS('U',2,1,1,a,1,af,2,b,2,x,2,r1,r2,w,iw,info)
         CALL CHKXER('DPBRFS',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL DPBRFS('U',2,1,1,a,2,af,1,b,2,x,2,r1,r2,w,iw,info)
         CALL CHKXER('DPBRFS',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL DPBRFS('U',2,0,1,a,1,af,1,b,1,x,2,r1,r2,w,iw,info)
         CALL CHKXER('DPBRFS',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL DPBRFS('U',2,0,1,a,1,af,1,b,2,x,1,r1,r2,w,iw,info)
         CALL CHKXER('DPBRFS',INFot,NOUt,LERr,OK)
!
!        DPBCON
!
         SRNamt = 'DPBCON'
         INFot = 1
         CALL DPBCON('/',0,0,a,1,anrm,rcond,w,iw,info)
         CALL CHKXER('DPBCON',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DPBCON('U',-1,0,a,1,anrm,rcond,w,iw,info)
         CALL CHKXER('DPBCON',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DPBCON('U',1,-1,a,1,anrm,rcond,w,iw,info)
         CALL CHKXER('DPBCON',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL DPBCON('U',2,1,a,1,anrm,rcond,w,iw,info)
         CALL CHKXER('DPBCON',INFot,NOUt,LERr,OK)
!
!        DPBEQU
!
         SRNamt = 'DPBEQU'
         INFot = 1
         CALL DPBEQU('/',0,0,a,1,r1,rcond,anrm,info)
         CALL CHKXER('DPBEQU',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DPBEQU('U',-1,0,a,1,r1,rcond,anrm,info)
         CALL CHKXER('DPBEQU',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DPBEQU('U',1,-1,a,1,r1,rcond,anrm,info)
         CALL CHKXER('DPBEQU',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL DPBEQU('U',2,1,a,1,r1,rcond,anrm,info)
         CALL CHKXER('DPBEQU',INFot,NOUt,LERr,OK)
      ENDIF
!
!     Print a summary line.
!
      CALL alaesm(Path,OK,NOUt)
!
!
!     End of derrpo
!
      END SUBROUTINE derrpo
