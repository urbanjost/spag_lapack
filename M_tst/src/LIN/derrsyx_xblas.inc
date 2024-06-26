!*==derrsy.f90  processed by SPAG 7.51RB at 17:34 on  4 Mar 2022
!> \brief \b DERRSYX
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE derrsy( PATH, NUNIT )
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
!> derrsy tests the error exits for the DOUBLE PRECISION routines
!> for symmetric indefinite matrices.
!>
!> Note that this file is used only when the XBLAS are available,
!> otherwise derrsy.f defines this subroutine.
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
      SUBROUTINE derrsy(Path,Nunit)
      IMPLICIT NONE
!*--DERRSY62
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
      INTEGER ip(NMAX) , iw(NMAX)
      DOUBLE PRECISION a(NMAX,NMAX) , af(NMAX,NMAX) , b(NMAX) , e(NMAX) &
     &                 , r1(NMAX) , r2(NMAX) , w(3*NMAX) , x(NMAX) ,    &
     &                 s(NMAX) , err_bnds_n(NMAX,3) , err_bnds_c(NMAX,3)&
     &                 , params(1)
!     ..
!     .. External Functions ..
      LOGICAL LSAMEN
EXTERNAL    LSAMEN
!     ..
!     .. External Subroutines ..
EXTERNAL    DSPCON,   DSPRFS,   DSPTRF,   DSPTRI
EXTERNAL    DSPTRS,   DSYCON,   DSYCON_3,   DSYCON_ROOK,   DSYRFS
EXTERNAL    DSYTF2,   DSYTF2_RK,   DSYTF2_ROOK,   DSYTRF,   DSYTRF_RK
EXTERNAL    DSYTRF_ROOK,   DSYTRI,   DSYTRI_3,   DSYTRI_3X
EXTERNAL    DSYTRI_ROOK,   DSYTRI2,   DSYTRI2X,   DSYTRS,   DSYTRS_3
EXTERNAL    DSYTRS_ROOK,   DSYRFSX
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
         e(j) = 0.D0
         r1(j) = 0.D0
         r2(j) = 0.D0
         w(j) = 0.D0
         x(j) = 0.D0
         s(j) = 0.D0
         ip(j) = j
         iw(j) = j
      ENDDO
      anrm = 1.0D0
      rcond = 1.0D0
      OK = .TRUE.
!
      IF ( LSAMEN(2,c2,'SY') ) THEN
!
!        Test error exits of the routines that use factorization
!        of a symmetric indefinite matrix with patrial
!        (Bunch-Kaufman) pivoting.
!
!        DSYTRF
!
         SRNamt = 'DSYTRF'
         INFot = 1
         CALL DSYTRF('/',0,a,1,ip,w,1,info)
         CALL CHKXER('DSYTRF',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DSYTRF('U',-1,a,1,ip,w,1,info)
         CALL CHKXER('DSYTRF',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL DSYTRF('U',2,a,1,ip,w,4,info)
         CALL CHKXER('DSYTRF',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL DSYTRF('U',0,a,1,ip,w,0,info)
         CALL CHKXER('DSYTRF',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL DSYTRF('U',0,a,1,ip,w,-2,info)
         CALL CHKXER('DSYTRF',INFot,NOUt,LERr,OK)
!
!        DSYTF2
!
         SRNamt = 'DSYTF2'
         INFot = 1
         CALL DSYTF2('/',0,a,1,ip,info)
         CALL CHKXER('DSYTF2',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DSYTF2('U',-1,a,1,ip,info)
         CALL CHKXER('DSYTF2',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL DSYTF2('U',2,a,1,ip,info)
         CALL CHKXER('DSYTF2',INFot,NOUt,LERr,OK)
!
!        DSYTRI
!
         SRNamt = 'DSYTRI'
         INFot = 1
         CALL DSYTRI('/',0,a,1,ip,w,info)
         CALL CHKXER('DSYTRI',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DSYTRI('U',-1,a,1,ip,w,info)
         CALL CHKXER('DSYTRI',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL DSYTRI('U',2,a,1,ip,w,info)
         CALL CHKXER('DSYTRI',INFot,NOUt,LERr,OK)
!
!        DSYTRI2
!
         SRNamt = 'DSYTRI2'
         INFot = 1
         CALL DSYTRI2('/',0,a,1,ip,w,iw,info)
         CALL CHKXER('DSYTRI2',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DSYTRI2('U',-1,a,1,ip,w,iw,info)
         CALL CHKXER('DSYTRI2',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL DSYTRI2('U',2,a,1,ip,w,iw,info)
         CALL CHKXER('DSYTRI2',INFot,NOUt,LERr,OK)
!
!        DSYTRI2X
!
         SRNamt = 'DSYTRI2X'
         INFot = 1
         CALL DSYTRI2X('/',0,a,1,ip,w,1,info)
         CALL CHKXER('DSYTRI2X',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DSYTRI2X('U',-1,a,1,ip,w,1,info)
         CALL CHKXER('DSYTRI2X',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL DSYTRI2X('U',2,a,1,ip,w,1,info)
         CALL CHKXER('DSYTRI2X',INFot,NOUt,LERr,OK)
!
!        DSYTRS
!
         SRNamt = 'DSYTRS'
         INFot = 1
         CALL DSYTRS('/',0,0,a,1,ip,b,1,info)
         CALL CHKXER('DSYTRS',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DSYTRS('U',-1,0,a,1,ip,b,1,info)
         CALL CHKXER('DSYTRS',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DSYTRS('U',0,-1,a,1,ip,b,1,info)
         CALL CHKXER('DSYTRS',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL DSYTRS('U',2,1,a,1,ip,b,2,info)
         CALL CHKXER('DSYTRS',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL DSYTRS('U',2,1,a,2,ip,b,1,info)
         CALL CHKXER('DSYTRS',INFot,NOUt,LERr,OK)
!
!        DSYRFS
!
         SRNamt = 'DSYRFS'
         INFot = 1
         CALL DSYRFS('/',0,0,a,1,af,1,ip,b,1,x,1,r1,r2,w,iw,info)
         CALL CHKXER('DSYRFS',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DSYRFS('U',-1,0,a,1,af,1,ip,b,1,x,1,r1,r2,w,iw,info)
         CALL CHKXER('DSYRFS',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DSYRFS('U',0,-1,a,1,af,1,ip,b,1,x,1,r1,r2,w,iw,info)
         CALL CHKXER('DSYRFS',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL DSYRFS('U',2,1,a,1,af,2,ip,b,2,x,2,r1,r2,w,iw,info)
         CALL CHKXER('DSYRFS',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL DSYRFS('U',2,1,a,2,af,1,ip,b,2,x,2,r1,r2,w,iw,info)
         CALL CHKXER('DSYRFS',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL DSYRFS('U',2,1,a,2,af,2,ip,b,1,x,2,r1,r2,w,iw,info)
         CALL CHKXER('DSYRFS',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL DSYRFS('U',2,1,a,2,af,2,ip,b,2,x,1,r1,r2,w,iw,info)
         CALL CHKXER('DSYRFS',INFot,NOUt,LERr,OK)
!
!        DSYRFSX
!
         n_err_bnds = 3
         nparams = 0
         SRNamt = 'DSYRFSX'
         INFot = 1
         CALL DSYRFSX('/',eq,0,0,a,1,af,1,ip,s,b,1,x,1,rcond,berr,      &
     &                n_err_bnds,err_bnds_n,err_bnds_c,nparams,params,w,&
     &                iw,info)
         CALL CHKXER('DSYRFSX',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DSYRFSX('U',eq,-1,0,a,1,af,1,ip,s,b,1,x,1,rcond,berr,     &
     &                n_err_bnds,err_bnds_n,err_bnds_c,nparams,params,w,&
     &                iw,info)
         CALL CHKXER('DSYRFSX',INFot,NOUt,LERr,OK)
         eq = 'N'
         INFot = 3
         CALL DSYRFSX('U',eq,-1,0,a,1,af,1,ip,s,b,1,x,1,rcond,berr,     &
     &                n_err_bnds,err_bnds_n,err_bnds_c,nparams,params,w,&
     &                iw,info)
         CALL CHKXER('DSYRFSX',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL DSYRFSX('U',eq,0,-1,a,1,af,1,ip,s,b,1,x,1,rcond,berr,     &
     &                n_err_bnds,err_bnds_n,err_bnds_c,nparams,params,w,&
     &                iw,info)
         CALL CHKXER('DSYRFSX',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL DSYRFSX('U',eq,2,1,a,1,af,2,ip,s,b,2,x,2,rcond,berr,      &
     &                n_err_bnds,err_bnds_n,err_bnds_c,nparams,params,w,&
     &                iw,info)
         CALL CHKXER('DSYRFSX',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL DSYRFSX('U',eq,2,1,a,2,af,1,ip,s,b,2,x,2,rcond,berr,      &
     &                n_err_bnds,err_bnds_n,err_bnds_c,nparams,params,w,&
     &                iw,info)
         CALL CHKXER('DSYRFSX',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL DSYRFSX('U',eq,2,1,a,2,af,2,ip,s,b,1,x,2,rcond,berr,      &
     &                n_err_bnds,err_bnds_n,err_bnds_c,nparams,params,w,&
     &                iw,info)
         CALL CHKXER('DSYRFSX',INFot,NOUt,LERr,OK)
         INFot = 14
         CALL DSYRFSX('U',eq,2,1,a,2,af,2,ip,s,b,2,x,1,rcond,berr,      &
     &                n_err_bnds,err_bnds_n,err_bnds_c,nparams,params,w,&
     &                iw,info)
         CALL CHKXER('DSYRFSX',INFot,NOUt,LERr,OK)
!
!        DSYCON
!
         SRNamt = 'DSYCON'
         INFot = 1
         CALL DSYCON('/',0,a,1,ip,anrm,rcond,w,iw,info)
         CALL CHKXER('DSYCON',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DSYCON('U',-1,a,1,ip,anrm,rcond,w,iw,info)
         CALL CHKXER('DSYCON',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL DSYCON('U',2,a,1,ip,anrm,rcond,w,iw,info)
         CALL CHKXER('DSYCON',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL DSYCON('U',1,a,1,ip,-1.0D0,rcond,w,iw,info)
         CALL CHKXER('DSYCON',INFot,NOUt,LERr,OK)
!
      ELSEIF ( LSAMEN(2,c2,'SR') ) THEN
!
!        Test error exits of the routines that use factorization
!        of a symmetric indefinite matrix with rook
!        (bounded Bunch-Kaufman) pivoting.
!
!        DSYTRF_ROOK
!
         SRNamt = 'DSYTRF_ROOK'
         INFot = 1
         CALL DSYTRF_ROOK('/',0,a,1,ip,w,1,info)
         CALL CHKXER('DSYTRF_ROOK',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DSYTRF_ROOK('U',-1,a,1,ip,w,1,info)
         CALL CHKXER('DSYTRF_ROOK',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL DSYTRF_ROOK('U',2,a,1,ip,w,4,info)
         CALL CHKXER('DSYTRF_ROOK',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL DSYTRF_ROOK('U',0,a,1,ip,w,0,info)
         CALL CHKXER('DSYTRF_ROOK',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL DSYTRF_ROOK('U',0,a,1,ip,w,-2,info)
         CALL CHKXER('DSYTRF_ROOK',INFot,NOUt,LERr,OK)
!
!        DSYTF2_ROOK
!
         SRNamt = 'DSYTF2_ROOK'
         INFot = 1
         CALL DSYTF2_ROOK('/',0,a,1,ip,info)
         CALL CHKXER('DSYTF2_ROOK',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DSYTF2_ROOK('U',-1,a,1,ip,info)
         CALL CHKXER('DSYTF2_ROOK',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL DSYTF2_ROOK('U',2,a,1,ip,info)
         CALL CHKXER('DSYTF2_ROOK',INFot,NOUt,LERr,OK)
!
!        DSYTRI_ROOK
!
         SRNamt = 'DSYTRI_ROOK'
         INFot = 1
         CALL DSYTRI_ROOK('/',0,a,1,ip,w,info)
         CALL CHKXER('DSYTRI_ROOK',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DSYTRI_ROOK('U',-1,a,1,ip,w,info)
         CALL CHKXER('DSYTRI_ROOK',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL DSYTRI_ROOK('U',2,a,1,ip,w,info)
         CALL CHKXER('DSYTRI_ROOK',INFot,NOUt,LERr,OK)
!
!        DSYTRS_ROOK
!
         SRNamt = 'DSYTRS_ROOK'
         INFot = 1
         CALL DSYTRS_ROOK('/',0,0,a,1,ip,b,1,info)
         CALL CHKXER('DSYTRS_ROOK',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DSYTRS_ROOK('U',-1,0,a,1,ip,b,1,info)
         CALL CHKXER('DSYTRS_ROOK',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DSYTRS_ROOK('U',0,-1,a,1,ip,b,1,info)
         CALL CHKXER('DSYTRS_ROOK',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL DSYTRS_ROOK('U',2,1,a,1,ip,b,2,info)
         CALL CHKXER('DSYTRS_ROOK',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL DSYTRS_ROOK('U',2,1,a,2,ip,b,1,info)
         CALL CHKXER('DSYTRS_ROOK',INFot,NOUt,LERr,OK)
!
!        DSYCON_ROOK
!
         SRNamt = 'DSYCON_ROOK'
         INFot = 1
         CALL DSYCON_ROOK('/',0,a,1,ip,anrm,rcond,w,iw,info)
         CALL CHKXER('DSYCON_ROOK',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DSYCON_ROOK('U',-1,a,1,ip,anrm,rcond,w,iw,info)
         CALL CHKXER('DSYCON_ROOK',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL DSYCON_ROOK('U',2,a,1,ip,anrm,rcond,w,iw,info)
         CALL CHKXER('DSYCON_ROOK',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL DSYCON_ROOK('U',1,a,1,ip,-1.0D0,rcond,w,iw,info)
         CALL CHKXER('DSYCON_ROOK',INFot,NOUt,LERr,OK)
!
      ELSEIF ( LSAMEN(2,c2,'SK') ) THEN
!
!        Test error exits of the routines that use factorization
!        of a symmetric indefinite matrix with rook
!        (bounded Bunch-Kaufman) pivoting with the new storage
!        format for factors L ( or U) and D.
!
!        L (or U) is stored in A, diagonal of D is stored on the
!        diagonal of A, subdiagonal of D is stored in a separate array E.
!
!        DSYTRF_RK
!
         SRNamt = 'DSYTRF_RK'
         INFot = 1
         CALL DSYTRF_RK('/',0,a,1,e,ip,w,1,info)
         CALL CHKXER('DSYTRF_RK',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DSYTRF_RK('U',-1,a,1,e,ip,w,1,info)
         CALL CHKXER('DSYTRF_RK',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL DSYTRF_RK('U',2,a,1,e,ip,w,1,info)
         CALL CHKXER('DSYTRF_RK',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL DSYTRF_RK('U',0,a,1,e,ip,w,0,info)
         CALL CHKXER('DSYTRF_RK',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL DSYTRF_RK('U',0,a,1,e,ip,w,-2,info)
         CALL CHKXER('DSYTRF_RK',INFot,NOUt,LERr,OK)
!
!        DSYTF2_RK
!
         SRNamt = 'DSYTF2_RK'
         INFot = 1
         CALL DSYTF2_RK('/',0,a,1,e,ip,info)
         CALL CHKXER('DSYTF2_RK',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DSYTF2_RK('U',-1,a,1,e,ip,info)
         CALL CHKXER('DSYTF2_RK',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL DSYTF2_RK('U',2,a,1,e,ip,info)
         CALL CHKXER('DSYTF2_RK',INFot,NOUt,LERr,OK)
!
!        DSYTRI_3
!
         SRNamt = 'DSYTRI_3'
         INFot = 1
         CALL DSYTRI_3('/',0,a,1,e,ip,w,1,info)
         CALL CHKXER('DSYTRI_3',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DSYTRI_3('U',-1,a,1,e,ip,w,1,info)
         CALL CHKXER('DSYTRI_3',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL DSYTRI_3('U',2,a,1,e,ip,w,1,info)
         CALL CHKXER('DSYTRI_3',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL DSYTRI_3('U',0,a,1,e,ip,w,0,info)
         CALL CHKXER('DSYTRI_3',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL DSYTRI_3('U',0,a,1,e,ip,w,-2,info)
         CALL CHKXER('DSYTRI_3',INFot,NOUt,LERr,OK)
!
!        DSYTRI_3X
!
         SRNamt = 'DSYTRI_3X'
         INFot = 1
         CALL DSYTRI_3X('/',0,a,1,e,ip,w,1,info)
         CALL CHKXER('DSYTRI_3X',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DSYTRI_3X('U',-1,a,1,e,ip,w,1,info)
         CALL CHKXER('DSYTRI_3X',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL DSYTRI_3X('U',2,a,1,e,ip,w,1,info)
         CALL CHKXER('DSYTRI_3X',INFot,NOUt,LERr,OK)
!
!        DSYTRS_3
!
         SRNamt = 'DSYTRS_3'
         INFot = 1
         CALL DSYTRS_3('/',0,0,a,1,e,ip,b,1,info)
         CALL CHKXER('DSYTRS_3',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DSYTRS_3('U',-1,0,a,1,e,ip,b,1,info)
         CALL CHKXER('DSYTRS_3',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DSYTRS_3('U',0,-1,a,1,e,ip,b,1,info)
         CALL CHKXER('DSYTRS_3',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL DSYTRS_3('U',2,1,a,1,e,ip,b,2,info)
         CALL CHKXER('DSYTRS_3',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL DSYTRS_3('U',2,1,a,2,e,ip,b,1,info)
         CALL CHKXER('DSYTRS_3',INFot,NOUt,LERr,OK)
!
!        DSYCON_3
!
         SRNamt = 'DSYCON_3'
         INFot = 1
         CALL DSYCON_3('/',0,a,1,e,ip,anrm,rcond,w,iw,info)
         CALL CHKXER('DSYCON_3',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DSYCON_3('U',-1,a,1,e,ip,anrm,rcond,w,iw,info)
         CALL CHKXER('DSYCON_3',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL DSYCON_3('U',2,a,1,e,ip,anrm,rcond,w,iw,info)
         CALL CHKXER('DSYCON_3',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL DSYCON_3('U',1,a,1,e,ip,-1.0D0,rcond,w,iw,info)
         CALL CHKXER('DSYCON_3',INFot,NOUt,LERr,OK)
!
      ELSEIF ( LSAMEN(2,c2,'SP') ) THEN
!
!        Test error exits of the routines that use factorization
!        of a symmetric indefinite packed matrix with patrial
!        (Bunch-Kaufman) pivoting.
!
!        DSPTRF
!
         SRNamt = 'DSPTRF'
         INFot = 1
         CALL DSPTRF('/',0,a,ip,info)
         CALL CHKXER('DSPTRF',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DSPTRF('U',-1,a,ip,info)
         CALL CHKXER('DSPTRF',INFot,NOUt,LERr,OK)
!
!        DSPTRI
!
         SRNamt = 'DSPTRI'
         INFot = 1
         CALL DSPTRI('/',0,a,ip,w,info)
         CALL CHKXER('DSPTRI',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DSPTRI('U',-1,a,ip,w,info)
         CALL CHKXER('DSPTRI',INFot,NOUt,LERr,OK)
!
!        DSPTRS
!
         SRNamt = 'DSPTRS'
         INFot = 1
         CALL DSPTRS('/',0,0,a,ip,b,1,info)
         CALL CHKXER('DSPTRS',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DSPTRS('U',-1,0,a,ip,b,1,info)
         CALL CHKXER('DSPTRS',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DSPTRS('U',0,-1,a,ip,b,1,info)
         CALL CHKXER('DSPTRS',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL DSPTRS('U',2,1,a,ip,b,1,info)
         CALL CHKXER('DSPTRS',INFot,NOUt,LERr,OK)
!
!        DSPRFS
!
         SRNamt = 'DSPRFS'
         INFot = 1
         CALL DSPRFS('/',0,0,a,af,ip,b,1,x,1,r1,r2,w,iw,info)
         CALL CHKXER('DSPRFS',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DSPRFS('U',-1,0,a,af,ip,b,1,x,1,r1,r2,w,iw,info)
         CALL CHKXER('DSPRFS',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DSPRFS('U',0,-1,a,af,ip,b,1,x,1,r1,r2,w,iw,info)
         CALL CHKXER('DSPRFS',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL DSPRFS('U',2,1,a,af,ip,b,1,x,2,r1,r2,w,iw,info)
         CALL CHKXER('DSPRFS',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL DSPRFS('U',2,1,a,af,ip,b,2,x,1,r1,r2,w,iw,info)
         CALL CHKXER('DSPRFS',INFot,NOUt,LERr,OK)
!
!        DSPCON
!
         SRNamt = 'DSPCON'
         INFot = 1
         CALL DSPCON('/',0,a,ip,anrm,rcond,w,iw,info)
         CALL CHKXER('DSPCON',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DSPCON('U',-1,a,ip,anrm,rcond,w,iw,info)
         CALL CHKXER('DSPCON',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL DSPCON('U',1,a,ip,-1.0D0,rcond,w,iw,info)
         CALL CHKXER('DSPCON',INFot,NOUt,LERr,OK)
      ENDIF
!
!     Print a summary line.
!
      CALL alaesm(Path,OK,NOUt)
!
!
!     End of derrsy
!
      END SUBROUTINE derrsy
