!*==serrsy.f90  processed by SPAG 7.51RB at 17:34 on  4 Mar 2022
!> \brief \b SERRSYX
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE serrsy( PATH, NUNIT )
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
!> serrsy tests the error exits for the REAL routines
!> for symmetric indefinite matrices.
!>
!> Note that this file is used only when the XBLAS are available,
!> otherwise serrsy.f defines this subroutine.
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
      SUBROUTINE serrsy(Path,Nunit)
      IMPLICIT NONE
!*--SERRSY62
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
      INTEGER ip(NMAX) , iw(NMAX)
      REAL a(NMAX,NMAX) , af(NMAX,NMAX) , b(NMAX) , e(NMAX) , r1(NMAX) ,&
     &     r2(NMAX) , w(3*NMAX) , x(NMAX) , s(NMAX) , err_bnds_n(NMAX,3)&
     &     , err_bnds_c(NMAX,3) , params(1)
!     ..
!     .. External Functions ..
      LOGICAL LSAMEN
EXTERNAL    LSAMEN
!     ..
!     .. External Subroutines ..
EXTERNAL    SSPCON,   SSPRFS,   SSPTRF,   SSPTRI
EXTERNAL    SSPTRS,   SSYCON,   SSYCON_3,   SSYCON_ROOK,   SSYRFS
EXTERNAL    SSYTF2,   SSYTF2_RK,   SSYTF2_ROOK,   SSYTRF,   SSYTRF_RK
EXTERNAL    SSYTRF_ROOK,   SSYTRI,   SSYTRI_3,   SSYTRI_3X
EXTERNAL    SSYTRI_ROOK,   SSYTRI2,   SSYTRI2X,   SSYTRS,   SSYTRS_3
EXTERNAL    SSYTRS_ROOK,   SSYRFSX
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
         b(j) = 0.E+0
         e(j) = 0.E+0
         r1(j) = 0.E+0
         r2(j) = 0.E+0
         w(j) = 0.E+0
         x(j) = 0.E+0
         ip(j) = j
         iw(j) = j
      ENDDO
      anrm = 1.0
      rcond = 1.0
      OK = .TRUE.
!
      IF ( LSAMEN(2,c2,'SY') ) THEN
!
!        Test error exits of the routines that use factorization
!        of a symmetric indefinite matrix with patrial
!        (Bunch-Kaufman) pivoting.
!
!        SSYTRF
!
         SRNamt = 'SSYTRF'
         INFot = 1
         CALL SSYTRF('/',0,a,1,ip,w,1,info)
         CALL CHKXER('SSYTRF',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SSYTRF('U',-1,a,1,ip,w,1,info)
         CALL CHKXER('SSYTRF',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL SSYTRF('U',2,a,1,ip,w,4,info)
         CALL CHKXER('SSYTRF',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL SSYTRF('U',0,a,1,ip,w,0,info)
         CALL CHKXER('SSYTRF',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL SSYTRF('U',0,a,1,ip,w,-2,info)
         CALL CHKXER('SSYTRF',INFot,NOUt,LERr,OK)
!
!        SSYTF2
!
         SRNamt = 'SSYTF2'
         INFot = 1
         CALL SSYTF2('/',0,a,1,ip,info)
         CALL CHKXER('SSYTF2',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SSYTF2('U',-1,a,1,ip,info)
         CALL CHKXER('SSYTF2',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL SSYTF2('U',2,a,1,ip,info)
         CALL CHKXER('SSYTF2',INFot,NOUt,LERr,OK)
!
!        SSYTRI
!
         SRNamt = 'SSYTRI'
         INFot = 1
         CALL SSYTRI('/',0,a,1,ip,w,info)
         CALL CHKXER('SSYTRI',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SSYTRI('U',-1,a,1,ip,w,info)
         CALL CHKXER('SSYTRI',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL SSYTRI('U',2,a,1,ip,w,info)
         CALL CHKXER('SSYTRI',INFot,NOUt,LERr,OK)
!
!        SSYTRI2
!
         SRNamt = 'SSYTRI2'
         INFot = 1
         CALL SSYTRI2('/',0,a,1,ip,w,iw,info)
         CALL CHKXER('SSYTRI',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SSYTRI2('U',-1,a,1,ip,w,iw,info)
         CALL CHKXER('SSYTRI',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL SSYTRI2('U',2,a,1,ip,w,iw,info)
         CALL CHKXER('SSYTRI',INFot,NOUt,LERr,OK)
!
!        SSYTRI2X
!
         SRNamt = 'SSYTRI2X'
         INFot = 1
         CALL SSYTRI2X('/',0,a,1,ip,w,1,info)
         CALL CHKXER('SSYTRI2X',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SSYTRI2X('U',-1,a,1,ip,w,1,info)
         CALL CHKXER('SSYTRI2X',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL SSYTRI2X('U',2,a,1,ip,w,1,info)
         CALL CHKXER('SSYTRI2X',INFot,NOUt,LERr,OK)
!
!        SSYTRS
!
         SRNamt = 'SSYTRS'
         INFot = 1
         CALL SSYTRS('/',0,0,a,1,ip,b,1,info)
         CALL CHKXER('SSYTRS',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SSYTRS('U',-1,0,a,1,ip,b,1,info)
         CALL CHKXER('SSYTRS',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SSYTRS('U',0,-1,a,1,ip,b,1,info)
         CALL CHKXER('SSYTRS',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL SSYTRS('U',2,1,a,1,ip,b,2,info)
         CALL CHKXER('SSYTRS',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL SSYTRS('U',2,1,a,2,ip,b,1,info)
         CALL CHKXER('SSYTRS',INFot,NOUt,LERr,OK)
!
!        SSYRFS
!
         SRNamt = 'SSYRFS'
         INFot = 1
         CALL SSYRFS('/',0,0,a,1,af,1,ip,b,1,x,1,r1,r2,w,iw,info)
         CALL CHKXER('SSYRFS',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SSYRFS('U',-1,0,a,1,af,1,ip,b,1,x,1,r1,r2,w,iw,info)
         CALL CHKXER('SSYRFS',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SSYRFS('U',0,-1,a,1,af,1,ip,b,1,x,1,r1,r2,w,iw,info)
         CALL CHKXER('SSYRFS',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL SSYRFS('U',2,1,a,1,af,2,ip,b,2,x,2,r1,r2,w,iw,info)
         CALL CHKXER('SSYRFS',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL SSYRFS('U',2,1,a,2,af,1,ip,b,2,x,2,r1,r2,w,iw,info)
         CALL CHKXER('SSYRFS',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL SSYRFS('U',2,1,a,2,af,2,ip,b,1,x,2,r1,r2,w,iw,info)
         CALL CHKXER('SSYRFS',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL SSYRFS('U',2,1,a,2,af,2,ip,b,2,x,1,r1,r2,w,iw,info)
         CALL CHKXER('SSYRFS',INFot,NOUt,LERr,OK)
!
!        SSYRFSX
!
         n_err_bnds = 3
         nparams = 0
         SRNamt = 'SSYRFSX'
         INFot = 1
         CALL SSYRFSX('/',eq,0,0,a,1,af,1,ip,s,b,1,x,1,rcond,berr,      &
     &                n_err_bnds,err_bnds_n,err_bnds_c,nparams,params,w,&
     &                iw,info)
         CALL CHKXER('SSYRFSX',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SSYRFSX('U',eq,-1,0,a,1,af,1,ip,s,b,1,x,1,rcond,berr,     &
     &                n_err_bnds,err_bnds_n,err_bnds_c,nparams,params,w,&
     &                iw,info)
         CALL CHKXER('SSYRFSX',INFot,NOUt,LERr,OK)
         eq = 'N'
         INFot = 3
         CALL SSYRFSX('U',eq,-1,0,a,1,af,1,ip,s,b,1,x,1,rcond,berr,     &
     &                n_err_bnds,err_bnds_n,err_bnds_c,nparams,params,w,&
     &                iw,info)
         CALL CHKXER('SSYRFSX',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL SSYRFSX('U',eq,0,-1,a,1,af,1,ip,s,b,1,x,1,rcond,berr,     &
     &                n_err_bnds,err_bnds_n,err_bnds_c,nparams,params,w,&
     &                iw,info)
         CALL CHKXER('SSYRFSX',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL SSYRFSX('U',eq,2,1,a,1,af,2,ip,s,b,2,x,2,rcond,berr,      &
     &                n_err_bnds,err_bnds_n,err_bnds_c,nparams,params,w,&
     &                iw,info)
         CALL CHKXER('SSYRFSX',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL SSYRFSX('U',eq,2,1,a,2,af,1,ip,s,b,2,x,2,rcond,berr,      &
     &                n_err_bnds,err_bnds_n,err_bnds_c,nparams,params,w,&
     &                iw,info)
         CALL CHKXER('SSYRFSX',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL SSYRFSX('U',eq,2,1,a,2,af,2,ip,s,b,1,x,2,rcond,berr,      &
     &                n_err_bnds,err_bnds_n,err_bnds_c,nparams,params,w,&
     &                iw,info)
         CALL CHKXER('SSYRFSX',INFot,NOUt,LERr,OK)
         INFot = 14
         CALL SSYRFSX('U',eq,2,1,a,2,af,2,ip,s,b,2,x,1,rcond,berr,      &
     &                n_err_bnds,err_bnds_n,err_bnds_c,nparams,params,w,&
     &                iw,info)
         CALL CHKXER('SSYRFSX',INFot,NOUt,LERr,OK)
!
!        SSYCON
!
         SRNamt = 'SSYCON'
         INFot = 1
         CALL SSYCON('/',0,a,1,ip,anrm,rcond,w,iw,info)
         CALL CHKXER('SSYCON',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SSYCON('U',-1,a,1,ip,anrm,rcond,w,iw,info)
         CALL CHKXER('SSYCON',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL SSYCON('U',2,a,1,ip,anrm,rcond,w,iw,info)
         CALL CHKXER('SSYCON',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL SSYCON('U',1,a,1,ip,-1.0,rcond,w,iw,info)
         CALL CHKXER('SSYCON',INFot,NOUt,LERr,OK)
!
      ELSEIF ( LSAMEN(2,c2,'SR') ) THEN
!
!        Test error exits of the routines that use factorization
!        of a symmetric indefinite matrix with rook
!        (bounded Bunch-Kaufman) pivoting.
!
!        SSYTRF_ROOK
!
         SRNamt = 'SSYTRF_ROOK'
         INFot = 1
         CALL SSYTRF_ROOK('/',0,a,1,ip,w,1,info)
         CALL CHKXER('SSYTRF_ROOK',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SSYTRF_ROOK('U',-1,a,1,ip,w,1,info)
         CALL CHKXER('SSYTRF_ROOK',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL SSYTRF_ROOK('U',2,a,1,ip,w,4,info)
         CALL CHKXER('SSYTRF_ROOK',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL SSYTRF_ROOK('U',0,a,1,ip,w,0,info)
         CALL CHKXER('SSYTRF_ROOK',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL SSYTRF_ROOK('U',0,a,1,ip,w,-2,info)
         CALL CHKXER('SSYTRF_ROOK',INFot,NOUt,LERr,OK)
!
!        SSYTF2_ROOK
!
         SRNamt = 'SSYTF2_ROOK'
         INFot = 1
         CALL SSYTF2_ROOK('/',0,a,1,ip,info)
         CALL CHKXER('SSYTF2_ROOK',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SSYTF2_ROOK('U',-1,a,1,ip,info)
         CALL CHKXER('SSYTF2_ROOK',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL SSYTF2_ROOK('U',2,a,1,ip,info)
         CALL CHKXER('SSYTF2_ROOK',INFot,NOUt,LERr,OK)
!
!        SSYTRI_ROOK
!
         SRNamt = 'SSYTRI_ROOK'
         INFot = 1
         CALL SSYTRI_ROOK('/',0,a,1,ip,w,info)
         CALL CHKXER('SSYTRI_ROOK',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SSYTRI_ROOK('U',-1,a,1,ip,w,info)
         CALL CHKXER('SSYTRI_ROOK',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL SSYTRI_ROOK('U',2,a,1,ip,w,info)
         CALL CHKXER('SSYTRI_ROOK',INFot,NOUt,LERr,OK)
!
!        SSYTRS_ROOK
!
         SRNamt = 'SSYTRS_ROOK'
         INFot = 1
         CALL SSYTRS_ROOK('/',0,0,a,1,ip,b,1,info)
         CALL CHKXER('SSYTRS_ROOK',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SSYTRS_ROOK('U',-1,0,a,1,ip,b,1,info)
         CALL CHKXER('SSYTRS_ROOK',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SSYTRS_ROOK('U',0,-1,a,1,ip,b,1,info)
         CALL CHKXER('SSYTRS_ROOK',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL SSYTRS_ROOK('U',2,1,a,1,ip,b,2,info)
         CALL CHKXER('SSYTRS_ROOK',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL SSYTRS_ROOK('U',2,1,a,2,ip,b,1,info)
         CALL CHKXER('SSYTRS_ROOK',INFot,NOUt,LERr,OK)
!
!        SSYCON_ROOK
!
         SRNamt = 'SSYCON_ROOK'
         INFot = 1
         CALL SSYCON_ROOK('/',0,a,1,ip,anrm,rcond,w,iw,info)
         CALL CHKXER('SSYCON_ROOK',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SSYCON_ROOK('U',-1,a,1,ip,anrm,rcond,w,iw,info)
         CALL CHKXER('SSYCON_ROOK',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL SSYCON_ROOK('U',2,a,1,ip,anrm,rcond,w,iw,info)
         CALL CHKXER('SSYCON_ROOK',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL SSYCON_ROOK('U',1,a,1,ip,-1.0,rcond,w,iw,info)
         CALL CHKXER('SSYCON_ROOK',INFot,NOUt,LERr,OK)
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
!        SSYTRF_RK
!
         SRNamt = 'SSYTRF_RK'
         INFot = 1
         CALL SSYTRF_RK('/',0,a,1,e,ip,w,1,info)
         CALL CHKXER('SSYTRF_RK',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SSYTRF_RK('U',-1,a,1,e,ip,w,1,info)
         CALL CHKXER('SSYTRF_RK',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL SSYTRF_RK('U',2,a,1,e,ip,w,4,info)
         CALL CHKXER('SSYTRF_RK',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL SSYTRF_RK('U',0,a,1,e,ip,w,0,info)
         CALL CHKXER('SSYTRF_RK',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL SSYTRF_RK('U',0,a,1,e,ip,w,-2,info)
         CALL CHKXER('SSYTRF_RK',INFot,NOUt,LERr,OK)
!
!        SSYTF2_RK
!
         SRNamt = 'SSYTF2_RK'
         INFot = 1
         CALL SSYTF2_RK('/',0,a,1,e,ip,info)
         CALL CHKXER('SSYTF2_RK',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SSYTF2_RK('U',-1,a,1,e,ip,info)
         CALL CHKXER('SSYTF2_RK',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL SSYTF2_RK('U',2,a,1,e,ip,info)
         CALL CHKXER('SSYTF2_RK',INFot,NOUt,LERr,OK)
!
!        SSYTRI_3
!
         SRNamt = 'SSYTRI_3'
         INFot = 1
         CALL SSYTRI_3('/',0,a,1,e,ip,w,1,info)
         CALL CHKXER('SSYTRI_3',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SSYTRI_3('U',-1,a,1,e,ip,w,1,info)
         CALL CHKXER('SSYTRI_3',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL SSYTRI_3('U',2,a,1,e,ip,w,1,info)
         CALL CHKXER('SSYTRI_3',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL SSYTRI_3('U',0,a,1,e,ip,w,0,info)
         CALL CHKXER('SSYTRI_3',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL SSYTRI_3('U',0,a,1,e,ip,w,-2,info)
         CALL CHKXER('SSYTRI_3',INFot,NOUt,LERr,OK)
!
!        SSYTRI_3X
!
         SRNamt = 'SSYTRI_3X'
         INFot = 1
         CALL SSYTRI_3X('/',0,a,1,e,ip,w,1,info)
         CALL CHKXER('SSYTRI_3X',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SSYTRI_3X('U',-1,a,1,e,ip,w,1,info)
         CALL CHKXER('SSYTRI_3X',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL SSYTRI_3X('U',2,a,1,e,ip,w,1,info)
         CALL CHKXER('SSYTRI_3X',INFot,NOUt,LERr,OK)
!
!        SSYTRS_3
!
         SRNamt = 'SSYTRS_3'
         INFot = 1
         CALL SSYTRS_3('/',0,0,a,1,e,ip,b,1,info)
         CALL CHKXER('SSYTRS_3',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SSYTRS_3('U',-1,0,a,1,e,ip,b,1,info)
         CALL CHKXER('SSYTRS_3',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SSYTRS_3('U',0,-1,a,1,e,ip,b,1,info)
         CALL CHKXER('SSYTRS_3',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL SSYTRS_3('U',2,1,a,1,e,ip,b,2,info)
         CALL CHKXER('SSYTRS_3',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL SSYTRS_3('U',2,1,a,2,e,ip,b,1,info)
         CALL CHKXER('SSYTRS_3',INFot,NOUt,LERr,OK)
!
!        SSYCON_3
!
         SRNamt = 'SSYCON_3'
         INFot = 1
         CALL SSYCON_3('/',0,a,1,e,ip,anrm,rcond,w,iw,info)
         CALL CHKXER('SSYCON_3',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SSYCON_3('U',-1,a,1,e,ip,anrm,rcond,w,iw,info)
         CALL CHKXER('SSYCON_3',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL SSYCON_3('U',2,a,1,e,ip,anrm,rcond,w,iw,info)
         CALL CHKXER('SSYCON_3',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL SSYCON_3('U',1,a,1,e,ip,-1.0E0,rcond,w,iw,info)
         CALL CHKXER('SSYCON_3',INFot,NOUt,LERr,OK)
!
      ELSEIF ( LSAMEN(2,c2,'SP') ) THEN
!
!        Test error exits of the routines that use factorization
!        of a symmetric indefinite packed matrix with patrial
!        (Bunch-Kaufman) pivoting.
!
!        SSPTRF
!
         SRNamt = 'SSPTRF'
         INFot = 1
         CALL SSPTRF('/',0,a,ip,info)
         CALL CHKXER('SSPTRF',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SSPTRF('U',-1,a,ip,info)
         CALL CHKXER('SSPTRF',INFot,NOUt,LERr,OK)
!
!        SSPTRI
!
         SRNamt = 'SSPTRI'
         INFot = 1
         CALL SSPTRI('/',0,a,ip,w,info)
         CALL CHKXER('SSPTRI',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SSPTRI('U',-1,a,ip,w,info)
         CALL CHKXER('SSPTRI',INFot,NOUt,LERr,OK)
!
!        SSPTRS
!
         SRNamt = 'SSPTRS'
         INFot = 1
         CALL SSPTRS('/',0,0,a,ip,b,1,info)
         CALL CHKXER('SSPTRS',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SSPTRS('U',-1,0,a,ip,b,1,info)
         CALL CHKXER('SSPTRS',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SSPTRS('U',0,-1,a,ip,b,1,info)
         CALL CHKXER('SSPTRS',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL SSPTRS('U',2,1,a,ip,b,1,info)
         CALL CHKXER('SSPTRS',INFot,NOUt,LERr,OK)
!
!        SSPRFS
!
         SRNamt = 'SSPRFS'
         INFot = 1
         CALL SSPRFS('/',0,0,a,af,ip,b,1,x,1,r1,r2,w,iw,info)
         CALL CHKXER('SSPRFS',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SSPRFS('U',-1,0,a,af,ip,b,1,x,1,r1,r2,w,iw,info)
         CALL CHKXER('SSPRFS',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SSPRFS('U',0,-1,a,af,ip,b,1,x,1,r1,r2,w,iw,info)
         CALL CHKXER('SSPRFS',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL SSPRFS('U',2,1,a,af,ip,b,1,x,2,r1,r2,w,iw,info)
         CALL CHKXER('SSPRFS',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL SSPRFS('U',2,1,a,af,ip,b,2,x,1,r1,r2,w,iw,info)
         CALL CHKXER('SSPRFS',INFot,NOUt,LERr,OK)
!
!        SSPCON
!
         SRNamt = 'SSPCON'
         INFot = 1
         CALL SSPCON('/',0,a,ip,anrm,rcond,w,iw,info)
         CALL CHKXER('SSPCON',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SSPCON('U',-1,a,ip,anrm,rcond,w,iw,info)
         CALL CHKXER('SSPCON',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL SSPCON('U',1,a,ip,-1.0,rcond,w,iw,info)
         CALL CHKXER('SSPCON',INFot,NOUt,LERr,OK)
      ENDIF
!
!     Print a summary line.
!
      CALL alaesm(Path,OK,NOUt)
!
!
!     End of serrsy
!
      END SUBROUTINE serrsy
