!*==serrpo.f90  processed by SPAG 7.51RB at 17:34 on  4 Mar 2022
!> \brief \b SERRPOX
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE serrpo( PATH, NUNIT )
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
!> serrpo tests the error exits for the REAL routines
!> for symmetric positive definite matrices.
!>
!> Note that this file is used only when the XBLAS are available,
!> otherwise serrpo.f defines this subroutine.
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
      SUBROUTINE serrpo(Path,Nunit)
      IMPLICIT NONE
!*--SERRPO62
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
      INTEGER iw(NMAX)
      REAL a(NMAX,NMAX) , af(NMAX,NMAX) , b(NMAX) , r1(NMAX) , r2(NMAX) &
     &     , w(3*NMAX) , x(NMAX) , s(NMAX) , err_bnds_n(NMAX,3) ,       &
     &     err_bnds_c(NMAX,3) , params(1)
!     ..
!     .. External Functions ..
      LOGICAL LSAMEN
EXTERNAL    LSAMEN
!     ..
!     .. External Subroutines ..
EXTERNAL    CHKXER,   SPBCON,   SPBEQU,   SPBRFS,   SPBTF2
EXTERNAL    SPBTRF,   SPBTRS,   SPOCON,   SPOEQU,   SPORFS,   SPOTF2
EXTERNAL    SPOTRF,   SPOTRI,   SPOTRS,   SPPCON,   SPPEQU,   SPPRFS
EXTERNAL    SPPTRF,   SPPTRI,   SPPTRS,   SPOEQUB,   SPORFSX
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
         s(j) = 0.
         iw(j) = j
      ENDDO
      OK = .TRUE.
!
      IF ( LSAMEN(2,c2,'PO') ) THEN
!
!        Test error exits of the routines that use the Cholesky
!        decomposition of a symmetric positive definite matrix.
!
!        SPOTRF
!
         SRNamt = 'SPOTRF'
         INFot = 1
         CALL SPOTRF('/',0,a,1,info)
         CALL CHKXER('SPOTRF',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SPOTRF('U',-1,a,1,info)
         CALL CHKXER('SPOTRF',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL SPOTRF('U',2,a,1,info)
         CALL CHKXER('SPOTRF',INFot,NOUt,LERr,OK)
!
!        SPOTF2
!
         SRNamt = 'SPOTF2'
         INFot = 1
         CALL SPOTF2('/',0,a,1,info)
         CALL CHKXER('SPOTF2',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SPOTF2('U',-1,a,1,info)
         CALL CHKXER('SPOTF2',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL SPOTF2('U',2,a,1,info)
         CALL CHKXER('SPOTF2',INFot,NOUt,LERr,OK)
!
!        SPOTRI
!
         SRNamt = 'SPOTRI'
         INFot = 1
         CALL SPOTRI('/',0,a,1,info)
         CALL CHKXER('SPOTRI',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SPOTRI('U',-1,a,1,info)
         CALL CHKXER('SPOTRI',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL SPOTRI('U',2,a,1,info)
         CALL CHKXER('SPOTRI',INFot,NOUt,LERr,OK)
!
!        SPOTRS
!
         SRNamt = 'SPOTRS'
         INFot = 1
         CALL SPOTRS('/',0,0,a,1,b,1,info)
         CALL CHKXER('SPOTRS',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SPOTRS('U',-1,0,a,1,b,1,info)
         CALL CHKXER('SPOTRS',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SPOTRS('U',0,-1,a,1,b,1,info)
         CALL CHKXER('SPOTRS',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL SPOTRS('U',2,1,a,1,b,2,info)
         CALL CHKXER('SPOTRS',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL SPOTRS('U',2,1,a,2,b,1,info)
         CALL CHKXER('SPOTRS',INFot,NOUt,LERr,OK)
!
!        SPORFS
!
         SRNamt = 'SPORFS'
         INFot = 1
         CALL SPORFS('/',0,0,a,1,af,1,b,1,x,1,r1,r2,w,iw,info)
         CALL CHKXER('SPORFS',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SPORFS('U',-1,0,a,1,af,1,b,1,x,1,r1,r2,w,iw,info)
         CALL CHKXER('SPORFS',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SPORFS('U',0,-1,a,1,af,1,b,1,x,1,r1,r2,w,iw,info)
         CALL CHKXER('SPORFS',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL SPORFS('U',2,1,a,1,af,2,b,2,x,2,r1,r2,w,iw,info)
         CALL CHKXER('SPORFS',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL SPORFS('U',2,1,a,2,af,1,b,2,x,2,r1,r2,w,iw,info)
         CALL CHKXER('SPORFS',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL SPORFS('U',2,1,a,2,af,2,b,1,x,2,r1,r2,w,iw,info)
         CALL CHKXER('SPORFS',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL SPORFS('U',2,1,a,2,af,2,b,2,x,1,r1,r2,w,iw,info)
         CALL CHKXER('SPORFS',INFot,NOUt,LERr,OK)
!
!        SPORFSX
!
         n_err_bnds = 3
         nparams = 0
         SRNamt = 'SPORFSX'
         INFot = 1
         CALL SPORFSX('/',eq,0,0,a,1,af,1,s,b,1,x,1,rcond,berr,         &
     &                n_err_bnds,err_bnds_n,err_bnds_c,nparams,params,w,&
     &                iw,info)
         CALL CHKXER('SPORFSX',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SPORFSX('U',"/",-1,0,a,1,af,1,s,b,1,x,1,rcond,berr,       &
     &                n_err_bnds,err_bnds_n,err_bnds_c,nparams,params,w,&
     &                iw,info)
         CALL CHKXER('SPORFSX',INFot,NOUt,LERr,OK)
         eq = 'N'
         INFot = 3
         CALL SPORFSX('U',eq,-1,0,a,1,af,1,s,b,1,x,1,rcond,berr,        &
     &                n_err_bnds,err_bnds_n,err_bnds_c,nparams,params,w,&
     &                iw,info)
         CALL CHKXER('SPORFSX',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL SPORFSX('U',eq,0,-1,a,1,af,1,s,b,1,x,1,rcond,berr,        &
     &                n_err_bnds,err_bnds_n,err_bnds_c,nparams,params,w,&
     &                iw,info)
         CALL CHKXER('SPORFSX',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL SPORFSX('U',eq,2,1,a,1,af,2,s,b,2,x,2,rcond,berr,         &
     &                n_err_bnds,err_bnds_n,err_bnds_c,nparams,params,w,&
     &                iw,info)
         CALL CHKXER('SPORFSX',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL SPORFSX('U',eq,2,1,a,2,af,1,s,b,2,x,2,rcond,berr,         &
     &                n_err_bnds,err_bnds_n,err_bnds_c,nparams,params,w,&
     &                iw,info)
         CALL CHKXER('SPORFSX',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL SPORFSX('U',eq,2,1,a,2,af,2,s,b,1,x,2,rcond,berr,         &
     &                n_err_bnds,err_bnds_n,err_bnds_c,nparams,params,w,&
     &                iw,info)
         CALL CHKXER('SPORFSX',INFot,NOUt,LERr,OK)
         INFot = 13
         CALL SPORFSX('U',eq,2,1,a,2,af,2,s,b,2,x,1,rcond,berr,         &
     &                n_err_bnds,err_bnds_n,err_bnds_c,nparams,params,w,&
     &                iw,info)
         CALL CHKXER('SPORFSX',INFot,NOUt,LERr,OK)
!
!        SPOCON
!
         SRNamt = 'SPOCON'
         INFot = 1
         CALL SPOCON('/',0,a,1,anrm,rcond,w,iw,info)
         CALL CHKXER('SPOCON',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SPOCON('U',-1,a,1,anrm,rcond,w,iw,info)
         CALL CHKXER('SPOCON',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL SPOCON('U',2,a,1,anrm,rcond,w,iw,info)
         CALL CHKXER('SPOCON',INFot,NOUt,LERr,OK)
!
!        SPOEQU
!
         SRNamt = 'SPOEQU'
         INFot = 1
         CALL SPOEQU(-1,a,1,r1,rcond,anrm,info)
         CALL CHKXER('SPOEQU',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SPOEQU(2,a,1,r1,rcond,anrm,info)
         CALL CHKXER('SPOEQU',INFot,NOUt,LERr,OK)
!
!        SPOEQUB
!
         SRNamt = 'SPOEQUB'
         INFot = 1
         CALL SPOEQUB(-1,a,1,r1,rcond,anrm,info)
         CALL CHKXER('SPOEQUB',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SPOEQUB(2,a,1,r1,rcond,anrm,info)
         CALL CHKXER('SPOEQUB',INFot,NOUt,LERr,OK)
!
      ELSEIF ( LSAMEN(2,c2,'PP') ) THEN
!
!        Test error exits of the routines that use the Cholesky
!        decomposition of a symmetric positive definite packed matrix.
!
!        SPPTRF
!
         SRNamt = 'SPPTRF'
         INFot = 1
         CALL SPPTRF('/',0,a,info)
         CALL CHKXER('SPPTRF',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SPPTRF('U',-1,a,info)
         CALL CHKXER('SPPTRF',INFot,NOUt,LERr,OK)
!
!        SPPTRI
!
         SRNamt = 'SPPTRI'
         INFot = 1
         CALL SPPTRI('/',0,a,info)
         CALL CHKXER('SPPTRI',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SPPTRI('U',-1,a,info)
         CALL CHKXER('SPPTRI',INFot,NOUt,LERr,OK)
!
!        SPPTRS
!
         SRNamt = 'SPPTRS'
         INFot = 1
         CALL SPPTRS('/',0,0,a,b,1,info)
         CALL CHKXER('SPPTRS',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SPPTRS('U',-1,0,a,b,1,info)
         CALL CHKXER('SPPTRS',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SPPTRS('U',0,-1,a,b,1,info)
         CALL CHKXER('SPPTRS',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL SPPTRS('U',2,1,a,b,1,info)
         CALL CHKXER('SPPTRS',INFot,NOUt,LERr,OK)
!
!        SPPRFS
!
         SRNamt = 'SPPRFS'
         INFot = 1
         CALL SPPRFS('/',0,0,a,af,b,1,x,1,r1,r2,w,iw,info)
         CALL CHKXER('SPPRFS',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SPPRFS('U',-1,0,a,af,b,1,x,1,r1,r2,w,iw,info)
         CALL CHKXER('SPPRFS',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SPPRFS('U',0,-1,a,af,b,1,x,1,r1,r2,w,iw,info)
         CALL CHKXER('SPPRFS',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL SPPRFS('U',2,1,a,af,b,1,x,2,r1,r2,w,iw,info)
         CALL CHKXER('SPPRFS',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL SPPRFS('U',2,1,a,af,b,2,x,1,r1,r2,w,iw,info)
         CALL CHKXER('SPPRFS',INFot,NOUt,LERr,OK)
!
!        SPPCON
!
         SRNamt = 'SPPCON'
         INFot = 1
         CALL SPPCON('/',0,a,anrm,rcond,w,iw,info)
         CALL CHKXER('SPPCON',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SPPCON('U',-1,a,anrm,rcond,w,iw,info)
         CALL CHKXER('SPPCON',INFot,NOUt,LERr,OK)
!
!        SPPEQU
!
         SRNamt = 'SPPEQU'
         INFot = 1
         CALL SPPEQU('/',0,a,r1,rcond,anrm,info)
         CALL CHKXER('SPPEQU',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SPPEQU('U',-1,a,r1,rcond,anrm,info)
         CALL CHKXER('SPPEQU',INFot,NOUt,LERr,OK)
!
      ELSEIF ( LSAMEN(2,c2,'PB') ) THEN
!
!        Test error exits of the routines that use the Cholesky
!        decomposition of a symmetric positive definite band matrix.
!
!        SPBTRF
!
         SRNamt = 'SPBTRF'
         INFot = 1
         CALL SPBTRF('/',0,0,a,1,info)
         CALL CHKXER('SPBTRF',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SPBTRF('U',-1,0,a,1,info)
         CALL CHKXER('SPBTRF',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SPBTRF('U',1,-1,a,1,info)
         CALL CHKXER('SPBTRF',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL SPBTRF('U',2,1,a,1,info)
         CALL CHKXER('SPBTRF',INFot,NOUt,LERr,OK)
!
!        SPBTF2
!
         SRNamt = 'SPBTF2'
         INFot = 1
         CALL SPBTF2('/',0,0,a,1,info)
         CALL CHKXER('SPBTF2',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SPBTF2('U',-1,0,a,1,info)
         CALL CHKXER('SPBTF2',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SPBTF2('U',1,-1,a,1,info)
         CALL CHKXER('SPBTF2',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL SPBTF2('U',2,1,a,1,info)
         CALL CHKXER('SPBTF2',INFot,NOUt,LERr,OK)
!
!        SPBTRS
!
         SRNamt = 'SPBTRS'
         INFot = 1
         CALL SPBTRS('/',0,0,0,a,1,b,1,info)
         CALL CHKXER('SPBTRS',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SPBTRS('U',-1,0,0,a,1,b,1,info)
         CALL CHKXER('SPBTRS',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SPBTRS('U',1,-1,0,a,1,b,1,info)
         CALL CHKXER('SPBTRS',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL SPBTRS('U',0,0,-1,a,1,b,1,info)
         CALL CHKXER('SPBTRS',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL SPBTRS('U',2,1,1,a,1,b,1,info)
         CALL CHKXER('SPBTRS',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL SPBTRS('U',2,0,1,a,1,b,1,info)
         CALL CHKXER('SPBTRS',INFot,NOUt,LERr,OK)
!
!        SPBRFS
!
         SRNamt = 'SPBRFS'
         INFot = 1
         CALL SPBRFS('/',0,0,0,a,1,af,1,b,1,x,1,r1,r2,w,iw,info)
         CALL CHKXER('SPBRFS',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SPBRFS('U',-1,0,0,a,1,af,1,b,1,x,1,r1,r2,w,iw,info)
         CALL CHKXER('SPBRFS',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SPBRFS('U',1,-1,0,a,1,af,1,b,1,x,1,r1,r2,w,iw,info)
         CALL CHKXER('SPBRFS',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL SPBRFS('U',0,0,-1,a,1,af,1,b,1,x,1,r1,r2,w,iw,info)
         CALL CHKXER('SPBRFS',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL SPBRFS('U',2,1,1,a,1,af,2,b,2,x,2,r1,r2,w,iw,info)
         CALL CHKXER('SPBRFS',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL SPBRFS('U',2,1,1,a,2,af,1,b,2,x,2,r1,r2,w,iw,info)
         CALL CHKXER('SPBRFS',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL SPBRFS('U',2,0,1,a,1,af,1,b,1,x,2,r1,r2,w,iw,info)
         CALL CHKXER('SPBRFS',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL SPBRFS('U',2,0,1,a,1,af,1,b,2,x,1,r1,r2,w,iw,info)
         CALL CHKXER('SPBRFS',INFot,NOUt,LERr,OK)
!
!        SPBCON
!
         SRNamt = 'SPBCON'
         INFot = 1
         CALL SPBCON('/',0,0,a,1,anrm,rcond,w,iw,info)
         CALL CHKXER('SPBCON',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SPBCON('U',-1,0,a,1,anrm,rcond,w,iw,info)
         CALL CHKXER('SPBCON',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SPBCON('U',1,-1,a,1,anrm,rcond,w,iw,info)
         CALL CHKXER('SPBCON',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL SPBCON('U',2,1,a,1,anrm,rcond,w,iw,info)
         CALL CHKXER('SPBCON',INFot,NOUt,LERr,OK)
!
!        SPBEQU
!
         SRNamt = 'SPBEQU'
         INFot = 1
         CALL SPBEQU('/',0,0,a,1,r1,rcond,anrm,info)
         CALL CHKXER('SPBEQU',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SPBEQU('U',-1,0,a,1,r1,rcond,anrm,info)
         CALL CHKXER('SPBEQU',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SPBEQU('U',1,-1,a,1,r1,rcond,anrm,info)
         CALL CHKXER('SPBEQU',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL SPBEQU('U',2,1,a,1,r1,rcond,anrm,info)
         CALL CHKXER('SPBEQU',INFot,NOUt,LERr,OK)
      ENDIF
!
!     Print a summary line.
!
      CALL alaesm(Path,OK,NOUt)
!
!
!     End of serrpo
!
      END SUBROUTINE serrpo
