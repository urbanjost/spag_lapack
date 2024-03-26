!*==zerrsy.f90  processed by SPAG 7.51RB at 17:35 on  4 Mar 2022
!> \brief \b ZERRSYX
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE zerrsy( PATH, NUNIT )
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
!> zerrsy tests the error exits for the COMPLEX*16 routines
!> for symmetric indefinite matrices.
!>
!> Note that this file is used only when the XBLAS are available,
!> otherwise zerrsy.f defines this subroutine.
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
      SUBROUTINE zerrsy(Path,Nunit)
      IMPLICIT NONE
!*--ZERRSY62
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
      INTEGER ip(NMAX)
      DOUBLE PRECISION r(NMAX) , r1(NMAX) , r2(NMAX) , s(NMAX) ,        &
     &                 err_bnds_n(NMAX,3) , err_bnds_c(NMAX,3) ,        &
     &                 params(1)
      COMPLEX*16 a(NMAX,NMAX) , af(NMAX,NMAX) , b(NMAX) , e(NMAX) ,     &
     &           w(2*NMAX) , x(NMAX)
!     ..
!     .. External Functions ..
      LOGICAL LSAMEN
EXTERNAL    LSAMEN
!     ..
!     .. External Subroutines ..
EXTERNAL    ZSPCON,   ZSPRFS,   ZSPTRF,   ZSPTRI
EXTERNAL    ZSPTRS,   ZSYCON,   ZSYCON_3,   ZSYCON_ROOK,   ZSYRFS
EXTERNAL    ZSYTF2,   ZSYTF2_RK,   ZSYTF2_ROOK,   ZSYTRF,   ZSYTRF_RK
EXTERNAL    ZSYTRF_ROOK,   ZSYTRI,   ZSYTRI_3,   ZSYTRI_3X
EXTERNAL    ZSYTRI_ROOK,   ZSYTRI2,   ZSYTRI2X,   ZSYTRS,   ZSYTRS_3
EXTERNAL    ZSYTRS_ROOK,   ZSYRFSX
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
         e(j) = 0.D0
         r1(j) = 0.D0
         r2(j) = 0.D0
         w(j) = 0.D0
         x(j) = 0.D0
         s(j) = 0.D0
         ip(j) = j
      ENDDO
      anrm = 1.0D0
      OK = .TRUE.
!
      IF ( LSAMEN(2,c2,'SY') ) THEN
!
!        Test error exits of the routines that use factorization
!        of a symmetric indefinite matrix with patrial
!        (Bunch-Kaufman) diagonal pivoting method.
!
!        ZSYTRF
!
         SRNamt = 'ZSYTRF'
         INFot = 1
         CALL ZSYTRF('/',0,a,1,ip,w,1,info)
         CALL CHKXER('ZSYTRF',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZSYTRF('U',-1,a,1,ip,w,1,info)
         CALL CHKXER('ZSYTRF',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZSYTRF('U',2,a,1,ip,w,4,info)
         CALL CHKXER('ZSYTRF',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL ZSYTRF('U',0,a,1,ip,w,0,info)
         CALL CHKXER('ZSYTRF',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL ZSYTRF('U',0,a,1,ip,w,-2,info)
         CALL CHKXER('ZSYTRF',INFot,NOUt,LERr,OK)
!
!        ZSYTF2
!
         SRNamt = 'ZSYTF2'
         INFot = 1
         CALL ZSYTF2('/',0,a,1,ip,info)
         CALL CHKXER('ZSYTF2',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZSYTF2('U',-1,a,1,ip,info)
         CALL CHKXER('ZSYTF2',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZSYTF2('U',2,a,1,ip,info)
         CALL CHKXER('ZSYTF2',INFot,NOUt,LERr,OK)
!
!        ZSYTRI
!
         SRNamt = 'ZSYTRI'
         INFot = 1
         CALL ZSYTRI('/',0,a,1,ip,w,info)
         CALL CHKXER('ZSYTRI',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZSYTRI('U',-1,a,1,ip,w,info)
         CALL CHKXER('ZSYTRI',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZSYTRI('U',2,a,1,ip,w,info)
         CALL CHKXER('ZSYTRI',INFot,NOUt,LERr,OK)
!
!        ZSYTRI2
!
         SRNamt = 'ZSYTRI2'
         INFot = 1
         CALL ZSYTRI2('/',0,a,1,ip,w,1,info)
         CALL CHKXER('ZSYTRI2',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZSYTRI2('U',-1,a,1,ip,w,1,info)
         CALL CHKXER('ZSYTRI2',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZSYTRI2('U',2,a,1,ip,w,1,info)
         CALL CHKXER('ZSYTRI2',INFot,NOUt,LERr,OK)
!
!        ZSYTRI2X
!
         SRNamt = 'ZSYTRI2X'
         INFot = 1
         CALL ZSYTRI2X('/',0,a,1,ip,w,1,info)
         CALL CHKXER('ZSYTRI2X',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZSYTRI2X('U',-1,a,1,ip,w,1,info)
         CALL CHKXER('ZSYTRI2X',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZSYTRI2X('U',2,a,1,ip,w,1,info)
         CALL CHKXER('ZSYTRI2X',INFot,NOUt,LERr,OK)
!
!        ZSYTRS
!
         SRNamt = 'ZSYTRS'
         INFot = 1
         CALL ZSYTRS('/',0,0,a,1,ip,b,1,info)
         CALL CHKXER('ZSYTRS',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZSYTRS('U',-1,0,a,1,ip,b,1,info)
         CALL CHKXER('ZSYTRS',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZSYTRS('U',0,-1,a,1,ip,b,1,info)
         CALL CHKXER('ZSYTRS',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL ZSYTRS('U',2,1,a,1,ip,b,2,info)
         CALL CHKXER('ZSYTRS',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL ZSYTRS('U',2,1,a,2,ip,b,1,info)
         CALL CHKXER('ZSYTRS',INFot,NOUt,LERr,OK)
!
!        ZSYRFS
!
         SRNamt = 'ZSYRFS'
         INFot = 1
         CALL ZSYRFS('/',0,0,a,1,af,1,ip,b,1,x,1,r1,r2,w,r,info)
         CALL CHKXER('ZSYRFS',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZSYRFS('U',-1,0,a,1,af,1,ip,b,1,x,1,r1,r2,w,r,info)
         CALL CHKXER('ZSYRFS',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZSYRFS('U',0,-1,a,1,af,1,ip,b,1,x,1,r1,r2,w,r,info)
         CALL CHKXER('ZSYRFS',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL ZSYRFS('U',2,1,a,1,af,2,ip,b,2,x,2,r1,r2,w,r,info)
         CALL CHKXER('ZSYRFS',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL ZSYRFS('U',2,1,a,2,af,1,ip,b,2,x,2,r1,r2,w,r,info)
         CALL CHKXER('ZSYRFS',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL ZSYRFS('U',2,1,a,2,af,2,ip,b,1,x,2,r1,r2,w,r,info)
         CALL CHKXER('ZSYRFS',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL ZSYRFS('U',2,1,a,2,af,2,ip,b,2,x,1,r1,r2,w,r,info)
         CALL CHKXER('ZSYRFS',INFot,NOUt,LERr,OK)
!
!        ZSYRFSX
!
         n_err_bnds = 3
         nparams = 0
         SRNamt = 'ZSYRFSX'
         INFot = 1
         CALL ZSYRFSX('/',eq,0,0,a,1,af,1,ip,s,b,1,x,1,rcond,berr,      &
     &                n_err_bnds,err_bnds_n,err_bnds_c,nparams,params,w,&
     &                r,info)
         CALL CHKXER('ZSYRFSX',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZSYRFSX('U',eq,-1,0,a,1,af,1,ip,s,b,1,x,1,rcond,berr,     &
     &                n_err_bnds,err_bnds_n,err_bnds_c,nparams,params,w,&
     &                r,info)
         CALL CHKXER('ZSYRFSX',INFot,NOUt,LERr,OK)
         eq = 'N'
         INFot = 3
         CALL ZSYRFSX('U',eq,-1,0,a,1,af,1,ip,s,b,1,x,1,rcond,berr,     &
     &                n_err_bnds,err_bnds_n,err_bnds_c,nparams,params,w,&
     &                r,info)
         CALL CHKXER('ZSYRFSX',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZSYRFSX('U',eq,0,-1,a,1,af,1,ip,s,b,1,x,1,rcond,berr,     &
     &                n_err_bnds,err_bnds_n,err_bnds_c,nparams,params,w,&
     &                r,info)
         CALL CHKXER('ZSYRFSX',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL ZSYRFSX('U',eq,2,1,a,1,af,2,ip,s,b,2,x,2,rcond,berr,      &
     &                n_err_bnds,err_bnds_n,err_bnds_c,nparams,params,w,&
     &                r,info)
         CALL CHKXER('ZSYRFSX',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL ZSYRFSX('U',eq,2,1,a,2,af,1,ip,s,b,2,x,2,rcond,berr,      &
     &                n_err_bnds,err_bnds_n,err_bnds_c,nparams,params,w,&
     &                r,info)
         CALL CHKXER('ZSYRFSX',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL ZSYRFSX('U',eq,2,1,a,2,af,2,ip,s,b,1,x,2,rcond,berr,      &
     &                n_err_bnds,err_bnds_n,err_bnds_c,nparams,params,w,&
     &                r,info)
         CALL CHKXER('ZSYRFSX',INFot,NOUt,LERr,OK)
         INFot = 14
         CALL ZSYRFSX('U',eq,2,1,a,2,af,2,ip,s,b,2,x,1,rcond,berr,      &
     &                n_err_bnds,err_bnds_n,err_bnds_c,nparams,params,w,&
     &                r,info)
         CALL CHKXER('ZSYRFSX',INFot,NOUt,LERr,OK)
!
!        ZSYCON
!
         SRNamt = 'ZSYCON'
         INFot = 1
         CALL ZSYCON('/',0,a,1,ip,anrm,rcond,w,info)
         CALL CHKXER('ZSYCON',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZSYCON('U',-1,a,1,ip,anrm,rcond,w,info)
         CALL CHKXER('ZSYCON',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZSYCON('U',2,a,1,ip,anrm,rcond,w,info)
         CALL CHKXER('ZSYCON',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL ZSYCON('U',1,a,1,ip,-anrm,rcond,w,info)
         CALL CHKXER('ZSYCON',INFot,NOUt,LERr,OK)
!
      ELSEIF ( LSAMEN(2,c2,'SR') ) THEN
!
!        Test error exits of the routines that use factorization
!        of a symmetric indefinite matrix with rook
!        (bounded Bunch-Kaufman) diagonal pivoting method.
!
!        ZSYTRF_ROOK
!
         SRNamt = 'ZSYTRF_ROOK'
         INFot = 1
         CALL ZSYTRF_ROOK('/',0,a,1,ip,w,1,info)
         CALL CHKXER('ZSYTRF_ROOK',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZSYTRF_ROOK('U',-1,a,1,ip,w,1,info)
         CALL CHKXER('ZSYTRF_ROOK',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZSYTRF_ROOK('U',2,a,1,ip,w,4,info)
         CALL CHKXER('ZSYTRF_ROOK',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL ZSYTRF_ROOK('U',0,a,1,ip,w,0,info)
         CALL CHKXER('ZSYTRF_ROOK',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL ZSYTRF_ROOK('U',0,a,1,ip,w,-2,info)
         CALL CHKXER('ZSYTRF_ROOK',INFot,NOUt,LERr,OK)
!
!        ZSYTF2_ROOK
!
         SRNamt = 'ZSYTF2_ROOK'
         INFot = 1
         CALL ZSYTF2_ROOK('/',0,a,1,ip,info)
         CALL CHKXER('ZSYTF2_ROOK',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZSYTF2_ROOK('U',-1,a,1,ip,info)
         CALL CHKXER('ZSYTF2_ROOK',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZSYTF2_ROOK('U',2,a,1,ip,info)
         CALL CHKXER('ZSYTF2_ROOK',INFot,NOUt,LERr,OK)
!
!        ZSYTRI_ROOK
!
         SRNamt = 'ZSYTRI_ROOK'
         INFot = 1
         CALL ZSYTRI_ROOK('/',0,a,1,ip,w,info)
         CALL CHKXER('ZSYTRI_ROOK',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZSYTRI_ROOK('U',-1,a,1,ip,w,info)
         CALL CHKXER('ZSYTRI_ROOK',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZSYTRI_ROOK('U',2,a,1,ip,w,info)
         CALL CHKXER('ZSYTRI_ROOK',INFot,NOUt,LERr,OK)
!
!        ZSYTRS_ROOK
!
         SRNamt = 'ZSYTRS_ROOK'
         INFot = 1
         CALL ZSYTRS_ROOK('/',0,0,a,1,ip,b,1,info)
         CALL CHKXER('ZSYTRS_ROOK',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZSYTRS_ROOK('U',-1,0,a,1,ip,b,1,info)
         CALL CHKXER('ZSYTRS_ROOK',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZSYTRS_ROOK('U',0,-1,a,1,ip,b,1,info)
         CALL CHKXER('ZSYTRS_ROOK',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL ZSYTRS_ROOK('U',2,1,a,1,ip,b,2,info)
         CALL CHKXER('ZSYTRS_ROOK',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL ZSYTRS_ROOK('U',2,1,a,2,ip,b,1,info)
         CALL CHKXER('ZSYTRS_ROOK',INFot,NOUt,LERr,OK)
!
!        ZSYCON_ROOK
!
         SRNamt = 'ZSYCON_ROOK'
         INFot = 1
         CALL ZSYCON_ROOK('/',0,a,1,ip,anrm,rcond,w,info)
         CALL CHKXER('ZSYCON_ROOK',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZSYCON_ROOK('U',-1,a,1,ip,anrm,rcond,w,info)
         CALL CHKXER('ZSYCON_ROOK',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZSYCON_ROOK('U',2,a,1,ip,anrm,rcond,w,info)
         CALL CHKXER('ZSYCON_ROOK',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL ZSYCON_ROOK('U',1,a,1,ip,-anrm,rcond,w,info)
         CALL CHKXER('ZSYCON_ROOK',INFot,NOUt,LERr,OK)
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
!        ZSYTRF_RK
!
         SRNamt = 'ZSYTRF_RK'
         INFot = 1
         CALL ZSYTRF_RK('/',0,a,1,e,ip,w,1,info)
         CALL CHKXER('ZSYTRF_RK',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZSYTRF_RK('U',-1,a,1,e,ip,w,1,info)
         CALL CHKXER('ZSYTRF_RK',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZSYTRF_RK('U',2,a,1,e,ip,w,4,info)
         CALL CHKXER('ZSYTRF_RK',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL ZSYTRF_RK('U',0,a,1,e,ip,w,0,info)
         CALL CHKXER('ZSYTRF_RK',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL ZSYTRF_RK('U',0,a,1,e,ip,w,-2,info)
         CALL CHKXER('ZSYTRF_RK',INFot,NOUt,LERr,OK)
!
!        ZSYTF2_RK
!
         SRNamt = 'ZSYTF2_RK'
         INFot = 1
         CALL ZSYTF2_RK('/',0,a,1,e,ip,info)
         CALL CHKXER('ZSYTF2_RK',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZSYTF2_RK('U',-1,a,1,e,ip,info)
         CALL CHKXER('ZSYTF2_RK',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZSYTF2_RK('U',2,a,1,e,ip,info)
         CALL CHKXER('ZSYTF2_RK',INFot,NOUt,LERr,OK)
!
!        ZSYTRI_3
!
         SRNamt = 'ZSYTRI_3'
         INFot = 1
         CALL ZSYTRI_3('/',0,a,1,e,ip,w,1,info)
         CALL CHKXER('ZSYTRI_3',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZSYTRI_3('U',-1,a,1,e,ip,w,1,info)
         CALL CHKXER('ZSYTRI_3',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZSYTRI_3('U',2,a,1,e,ip,w,1,info)
         CALL CHKXER('ZSYTRI_3',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL ZSYTRI_3('U',0,a,1,e,ip,w,0,info)
         CALL CHKXER('ZSYTRI_3',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL ZSYTRI_3('U',0,a,1,e,ip,w,-2,info)
         CALL CHKXER('ZSYTRI_3',INFot,NOUt,LERr,OK)
!
!        ZSYTRI_3X
!
         SRNamt = 'ZSYTRI_3X'
         INFot = 1
         CALL ZSYTRI_3X('/',0,a,1,e,ip,w,1,info)
         CALL CHKXER('ZSYTRI_3X',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZSYTRI_3X('U',-1,a,1,e,ip,w,1,info)
         CALL CHKXER('ZSYTRI_3X',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZSYTRI_3X('U',2,a,1,e,ip,w,1,info)
         CALL CHKXER('ZSYTRI_3X',INFot,NOUt,LERr,OK)
!
!        ZSYTRS_3
!
         SRNamt = 'ZSYTRS_3'
         INFot = 1
         CALL ZSYTRS_3('/',0,0,a,1,e,ip,b,1,info)
         CALL CHKXER('ZSYTRS_3',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZSYTRS_3('U',-1,0,a,1,e,ip,b,1,info)
         CALL CHKXER('ZSYTRS_3',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZSYTRS_3('U',0,-1,a,1,e,ip,b,1,info)
         CALL CHKXER('ZSYTRS_3',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL ZSYTRS_3('U',2,1,a,1,e,ip,b,2,info)
         CALL CHKXER('ZSYTRS_3',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL ZSYTRS_3('U',2,1,a,2,e,ip,b,1,info)
         CALL CHKXER('ZSYTRS_3',INFot,NOUt,LERr,OK)
!
!        ZSYCON_3
!
         SRNamt = 'ZSYCON_3'
         INFot = 1
         CALL ZSYCON_3('/',0,a,1,e,ip,anrm,rcond,w,info)
         CALL CHKXER('ZSYCON_3',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZSYCON_3('U',-1,a,1,e,ip,anrm,rcond,w,info)
         CALL CHKXER('ZSYCON_3',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZSYCON_3('U',2,a,1,e,ip,anrm,rcond,w,info)
         CALL CHKXER('ZSYCON_3',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL ZSYCON_3('U',1,a,1,e,ip,-1.0D0,rcond,w,info)
         CALL CHKXER('ZSYCON_3',INFot,NOUt,LERr,OK)
!
      ELSEIF ( LSAMEN(2,c2,'SP') ) THEN
!
!        Test error exits of the routines that use factorization
!        of a symmetric indefinite packed matrix with patrial
!        (Bunch-Kaufman) pivoting.
!
!        ZSPTRF
!
         SRNamt = 'ZSPTRF'
         INFot = 1
         CALL ZSPTRF('/',0,a,ip,info)
         CALL CHKXER('ZSPTRF',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZSPTRF('U',-1,a,ip,info)
         CALL CHKXER('ZSPTRF',INFot,NOUt,LERr,OK)
!
!        ZSPTRI
!
         SRNamt = 'ZSPTRI'
         INFot = 1
         CALL ZSPTRI('/',0,a,ip,w,info)
         CALL CHKXER('ZSPTRI',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZSPTRI('U',-1,a,ip,w,info)
         CALL CHKXER('ZSPTRI',INFot,NOUt,LERr,OK)
!
!        ZSPTRS
!
         SRNamt = 'ZSPTRS'
         INFot = 1
         CALL ZSPTRS('/',0,0,a,ip,b,1,info)
         CALL CHKXER('ZSPTRS',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZSPTRS('U',-1,0,a,ip,b,1,info)
         CALL CHKXER('ZSPTRS',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZSPTRS('U',0,-1,a,ip,b,1,info)
         CALL CHKXER('ZSPTRS',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL ZSPTRS('U',2,1,a,ip,b,1,info)
         CALL CHKXER('ZSPTRS',INFot,NOUt,LERr,OK)
!
!        ZSPRFS
!
         SRNamt = 'ZSPRFS'
         INFot = 1
         CALL ZSPRFS('/',0,0,a,af,ip,b,1,x,1,r1,r2,w,r,info)
         CALL CHKXER('ZSPRFS',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZSPRFS('U',-1,0,a,af,ip,b,1,x,1,r1,r2,w,r,info)
         CALL CHKXER('ZSPRFS',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZSPRFS('U',0,-1,a,af,ip,b,1,x,1,r1,r2,w,r,info)
         CALL CHKXER('ZSPRFS',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL ZSPRFS('U',2,1,a,af,ip,b,1,x,2,r1,r2,w,r,info)
         CALL CHKXER('ZSPRFS',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL ZSPRFS('U',2,1,a,af,ip,b,2,x,1,r1,r2,w,r,info)
         CALL CHKXER('ZSPRFS',INFot,NOUt,LERr,OK)
!
!        ZSPCON
!
         SRNamt = 'ZSPCON'
         INFot = 1
         CALL ZSPCON('/',0,a,ip,anrm,rcond,w,info)
         CALL CHKXER('ZSPCON',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZSPCON('U',-1,a,ip,anrm,rcond,w,info)
         CALL CHKXER('ZSPCON',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL ZSPCON('U',1,a,ip,-anrm,rcond,w,info)
         CALL CHKXER('ZSPCON',INFot,NOUt,LERr,OK)
      ENDIF
!
!     Print a summary line.
!
      CALL alaesm(Path,OK,NOUt)
!
!
!     End of zerrsy
!
      END SUBROUTINE zerrsy
