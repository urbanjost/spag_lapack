!*==cerrhe.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b CERRHEX
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE CERRHE( PATH, NUNIT )
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
!> CERRHE tests the error exits for the COMPLEX routines
!> for Hermitian indefinite matrices.
!>
!> Note that this file is used only when the XBLAS are available,
!> otherwise cerrhe.f defines this subroutine.
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
      SUBROUTINE CERRHE(Path,Nunit)
      IMPLICIT NONE
!*--CERRHE62
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
      INTEGER ip(NMAX)
      REAL r(NMAX) , r1(NMAX) , r2(NMAX) , s(NMAX) , err_bnds_n(NMAX,3) &
     &     , err_bnds_c(NMAX,3) , params(1)
      COMPLEX a(NMAX,NMAX) , af(NMAX,NMAX) , b(NMAX) , e(NMAX) ,        &
     &        w(2*NMAX) , x(NMAX)
!     ..
!     .. External Functions ..
      LOGICAL LSAMEN
      EXTERNAL LSAMEN
!     ..
!     .. External Subroutines ..
      EXTERNAL ALAESM , CHECON , CHECON_3 , CHECON_ROOK , CHERFS ,      &
     &         CHETF2 , CHETF2_RK , CHETF2_ROOK , CHETRF , CHETRF_RK ,  &
     &         CHETRF_ROOK , CHETRI , CHETRI_3 , CHETRI_3X ,            &
     &         CHETRI_ROOK , CHETRI2 , CHETRI2X , CHETRS , CHETRS_3 ,   &
     &         CHETRS_ROOK , CHKXER , CHPCON , CHPRFS , CHPTRF ,        &
     &         CHPTRI , CHPTRS , CHERFSX
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
         b(j) = 0.E+0
         e(j) = 0.E+0
         r1(j) = 0.E+0
         r2(j) = 0.E+0
         w(j) = 0.E+0
         x(j) = 0.E+0
         ip(j) = j
      ENDDO
      anrm = 1.0
      OK = .TRUE.
!
      IF ( LSAMEN(2,c2,'HE') ) THEN
!
!        Test error exits of the routines that use factorization
!        of a Hermitian indefinite matrix with patrial
!        (Bunch-Kaufman) diagonal pivoting method.
!
!        CHETRF
!
         SRNamt = 'CHETRF'
         INFot = 1
         CALL CHETRF('/',0,a,1,ip,w,1,info)
         CALL CHKXER('CHETRF',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CHETRF('U',-1,a,1,ip,w,1,info)
         CALL CHKXER('CHETRF',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CHETRF('U',2,a,1,ip,w,4,info)
         CALL CHKXER('CHETRF',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL CHETRF('U',0,a,1,ip,w,0,info)
         CALL CHKXER('CHETRF',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL CHETRF('U',0,a,1,ip,w,-2,info)
         CALL CHKXER('CHETRF',INFot,NOUt,LERr,OK)
!
!        CHETF2
!
         SRNamt = 'CHETF2'
         INFot = 1
         CALL CHETF2('/',0,a,1,ip,info)
         CALL CHKXER('CHETF2',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CHETF2('U',-1,a,1,ip,info)
         CALL CHKXER('CHETF2',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CHETF2('U',2,a,1,ip,info)
         CALL CHKXER('CHETF2',INFot,NOUt,LERr,OK)
!
!        CHETRI
!
         SRNamt = 'CHETRI'
         INFot = 1
         CALL CHETRI('/',0,a,1,ip,w,info)
         CALL CHKXER('CHETRI',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CHETRI('U',-1,a,1,ip,w,info)
         CALL CHKXER('CHETRI',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CHETRI('U',2,a,1,ip,w,info)
         CALL CHKXER('CHETRI',INFot,NOUt,LERr,OK)
!
!        CHETRI2
!
         SRNamt = 'CHETRI2'
         INFot = 1
         CALL CHETRI2('/',0,a,1,ip,w,1,info)
         CALL CHKXER('CHETRI2',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CHETRI2('U',-1,a,1,ip,w,1,info)
         CALL CHKXER('CHETRI2',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CHETRI2('U',2,a,1,ip,w,1,info)
         CALL CHKXER('CHETRI2',INFot,NOUt,LERr,OK)
!
!        CHETRI2X
!
         SRNamt = 'CHETRI2X'
         INFot = 1
         CALL CHETRI2X('/',0,a,1,ip,w,1,info)
         CALL CHKXER('CHETRI2X',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CHETRI2X('U',-1,a,1,ip,w,1,info)
         CALL CHKXER('CHETRI2X',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CHETRI2X('U',2,a,1,ip,w,1,info)
         CALL CHKXER('CHETRI2X',INFot,NOUt,LERr,OK)
!
!        CHETRS
!
         SRNamt = 'CHETRS'
         INFot = 1
         CALL CHETRS('/',0,0,a,1,ip,b,1,info)
         CALL CHKXER('CHETRS',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CHETRS('U',-1,0,a,1,ip,b,1,info)
         CALL CHKXER('CHETRS',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CHETRS('U',0,-1,a,1,ip,b,1,info)
         CALL CHKXER('CHETRS',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL CHETRS('U',2,1,a,1,ip,b,2,info)
         CALL CHKXER('CHETRS',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL CHETRS('U',2,1,a,2,ip,b,1,info)
         CALL CHKXER('CHETRS',INFot,NOUt,LERr,OK)
!
!        CHERFS
!
         SRNamt = 'CHERFS'
         INFot = 1
         CALL CHERFS('/',0,0,a,1,af,1,ip,b,1,x,1,r1,r2,w,r,info)
         CALL CHKXER('CHERFS',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CHERFS('U',-1,0,a,1,af,1,ip,b,1,x,1,r1,r2,w,r,info)
         CALL CHKXER('CHERFS',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CHERFS('U',0,-1,a,1,af,1,ip,b,1,x,1,r1,r2,w,r,info)
         CALL CHKXER('CHERFS',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL CHERFS('U',2,1,a,1,af,2,ip,b,2,x,2,r1,r2,w,r,info)
         CALL CHKXER('CHERFS',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL CHERFS('U',2,1,a,2,af,1,ip,b,2,x,2,r1,r2,w,r,info)
         CALL CHKXER('CHERFS',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL CHERFS('U',2,1,a,2,af,2,ip,b,1,x,2,r1,r2,w,r,info)
         CALL CHKXER('CHERFS',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL CHERFS('U',2,1,a,2,af,2,ip,b,2,x,1,r1,r2,w,r,info)
         CALL CHKXER('CHERFS',INFot,NOUt,LERr,OK)
!
!        CHECON
!
         SRNamt = 'CHECON'
         INFot = 1
         CALL CHECON('/',0,a,1,ip,anrm,rcond,w,info)
         CALL CHKXER('CHECON',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CHECON('U',-1,a,1,ip,anrm,rcond,w,info)
         CALL CHKXER('CHECON',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CHECON('U',2,a,1,ip,anrm,rcond,w,info)
         CALL CHKXER('CHECON',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL CHECON('U',1,a,1,ip,-anrm,rcond,w,info)
         CALL CHKXER('CHECON',INFot,NOUt,LERr,OK)
!
!        CHERFSX
!
         n_err_bnds = 3
         nparams = 0
         SRNamt = 'CHERFSX'
         INFot = 1
         CALL CHERFSX('/',eq,0,0,a,1,af,1,ip,s,b,1,x,1,rcond,berr,      &
     &                n_err_bnds,err_bnds_n,err_bnds_c,nparams,params,w,&
     &                r,info)
         CALL CHKXER('CHERFSX',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CHERFSX('U',eq,-1,0,a,1,af,1,ip,s,b,1,x,1,rcond,berr,     &
     &                n_err_bnds,err_bnds_n,err_bnds_c,nparams,params,w,&
     &                r,info)
         CALL CHKXER('CHERFSX',INFot,NOUt,LERr,OK)
         eq = 'N'
         INFot = 3
         CALL CHERFSX('U',eq,-1,0,a,1,af,1,ip,s,b,1,x,1,rcond,berr,     &
     &                n_err_bnds,err_bnds_n,err_bnds_c,nparams,params,w,&
     &                r,info)
         CALL CHKXER('CHERFSX',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CHERFSX('U',eq,0,-1,a,1,af,1,ip,s,b,1,x,1,rcond,berr,     &
     &                n_err_bnds,err_bnds_n,err_bnds_c,nparams,params,w,&
     &                r,info)
         CALL CHKXER('CHERFSX',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL CHERFSX('U',eq,2,1,a,1,af,2,ip,s,b,2,x,2,rcond,berr,      &
     &                n_err_bnds,err_bnds_n,err_bnds_c,nparams,params,w,&
     &                r,info)
         CALL CHKXER('CHERFSX',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL CHERFSX('U',eq,2,1,a,2,af,1,ip,s,b,2,x,2,rcond,berr,      &
     &                n_err_bnds,err_bnds_n,err_bnds_c,nparams,params,w,&
     &                r,info)
         CALL CHKXER('CHERFSX',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL CHERFSX('U',eq,2,1,a,2,af,2,ip,s,b,1,x,2,rcond,berr,      &
     &                n_err_bnds,err_bnds_n,err_bnds_c,nparams,params,w,&
     &                r,info)
         CALL CHKXER('CHERFSX',INFot,NOUt,LERr,OK)
         INFot = 14
         CALL CHERFSX('U',eq,2,1,a,2,af,2,ip,s,b,2,x,1,rcond,berr,      &
     &                n_err_bnds,err_bnds_n,err_bnds_c,nparams,params,w,&
     &                r,info)
         CALL CHKXER('CHERFSX',INFot,NOUt,LERr,OK)
!
      ELSEIF ( LSAMEN(2,c2,'HR') ) THEN
!
!        Test error exits of the routines that use factorization
!        of a Hermitian indefinite matrix with rook
!        (bounded Bunch-Kaufman) diagonal pivoting method.
!
!        CHETRF_ROOK
!
         SRNamt = 'CHETRF_ROOK'
         INFot = 1
         CALL CHETRF_ROOK('/',0,a,1,ip,w,1,info)
         CALL CHKXER('CHETRF_ROOK',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CHETRF_ROOK('U',-1,a,1,ip,w,1,info)
         CALL CHKXER('CHETRF_ROOK',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CHETRF_ROOK('U',2,a,1,ip,w,4,info)
         CALL CHKXER('CHETRF_ROOK',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL CHETRF_ROOK('U',0,a,1,ip,w,0,info)
         CALL CHKXER('CHETRF_ROOK',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL CHETRF_ROOK('U',0,a,1,ip,w,-2,info)
         CALL CHKXER('CHETRF_ROOK',INFot,NOUt,LERr,OK)
!
!        CHETF2_ROOK
!
         SRNamt = 'CHETF2_ROOK'
         INFot = 1
         CALL CHETF2_ROOK('/',0,a,1,ip,info)
         CALL CHKXER('CHETF2_ROOK',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CHETF2_ROOK('U',-1,a,1,ip,info)
         CALL CHKXER('CHETF2_ROOK',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CHETF2_ROOK('U',2,a,1,ip,info)
         CALL CHKXER('CHETF2_ROOK',INFot,NOUt,LERr,OK)
!
!        CHETRI_ROOK
!
         SRNamt = 'CHETRI_ROOK'
         INFot = 1
         CALL CHETRI_ROOK('/',0,a,1,ip,w,info)
         CALL CHKXER('CHETRI_ROOK',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CHETRI_ROOK('U',-1,a,1,ip,w,info)
         CALL CHKXER('CHETRI_ROOK',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CHETRI_ROOK('U',2,a,1,ip,w,info)
         CALL CHKXER('CHETRI_ROOK',INFot,NOUt,LERr,OK)
!
!        CHETRS_ROOK
!
         SRNamt = 'CHETRS_ROOK'
         INFot = 1
         CALL CHETRS_ROOK('/',0,0,a,1,ip,b,1,info)
         CALL CHKXER('CHETRS_ROOK',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CHETRS_ROOK('U',-1,0,a,1,ip,b,1,info)
         CALL CHKXER('CHETRS_ROOK',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CHETRS_ROOK('U',0,-1,a,1,ip,b,1,info)
         CALL CHKXER('CHETRS_ROOK',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL CHETRS_ROOK('U',2,1,a,1,ip,b,2,info)
         CALL CHKXER('CHETRS_ROOK',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL CHETRS_ROOK('U',2,1,a,2,ip,b,1,info)
         CALL CHKXER('CHETRS_ROOK',INFot,NOUt,LERr,OK)
!
!        CHECON_ROOK
!
         SRNamt = 'CHECON_ROOK'
         INFot = 1
         CALL CHECON_ROOK('/',0,a,1,ip,anrm,rcond,w,info)
         CALL CHKXER('CHECON_ROOK',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CHECON_ROOK('U',-1,a,1,ip,anrm,rcond,w,info)
         CALL CHKXER('CHECON_ROOK',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CHECON_ROOK('U',2,a,1,ip,anrm,rcond,w,info)
         CALL CHKXER('CHECON_ROOK',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL CHECON_ROOK('U',1,a,1,ip,-anrm,rcond,w,info)
         CALL CHKXER('CHECON_ROOK',INFot,NOUt,LERr,OK)
!
      ELSEIF ( LSAMEN(2,c2,'HK') ) THEN
!
!        Test error exits of the routines that use factorization
!        of a Hermitian indefinite matrix with rook
!        (bounded Bunch-Kaufman) pivoting with the new storage
!        format for factors L ( or U) and D.
!
!        L (or U) is stored in A, diagonal of D is stored on the
!        diagonal of A, subdiagonal of D is stored in a separate array E.
!
!        CHETRF_RK
!
         SRNamt = 'CHETRF_RK'
         INFot = 1
         CALL CHETRF_RK('/',0,a,1,e,ip,w,1,info)
         CALL CHKXER('CHETRF_RK',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CHETRF_RK('U',-1,a,1,e,ip,w,1,info)
         CALL CHKXER('CHETRF_RK',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CHETRF_RK('U',2,a,1,e,ip,w,4,info)
         CALL CHKXER('CHETRF_RK',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL CHETRF_RK('U',0,a,1,e,ip,w,0,info)
         CALL CHKXER('CHETRF_RK',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL CHETRF_RK('U',0,a,1,e,ip,w,-2,info)
         CALL CHKXER('CHETRF_RK',INFot,NOUt,LERr,OK)
!
!        CHETF2_RK
!
         SRNamt = 'CHETF2_RK'
         INFot = 1
         CALL CHETF2_RK('/',0,a,1,e,ip,info)
         CALL CHKXER('CHETF2_RK',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CHETF2_RK('U',-1,a,1,e,ip,info)
         CALL CHKXER('CHETF2_RK',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CHETF2_RK('U',2,a,1,e,ip,info)
         CALL CHKXER('CHETF2_RK',INFot,NOUt,LERr,OK)
!
!        CHETRI_3
!
         SRNamt = 'CHETRI_3'
         INFot = 1
         CALL CHETRI_3('/',0,a,1,e,ip,w,1,info)
         CALL CHKXER('CHETRI_3',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CHETRI_3('U',-1,a,1,e,ip,w,1,info)
         CALL CHKXER('CHETRI_3',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CHETRI_3('U',2,a,1,e,ip,w,1,info)
         CALL CHKXER('CHETRI_3',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL CHETRI_3('U',0,a,1,e,ip,w,0,info)
         CALL CHKXER('CHETRI_3',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL CHETRI_3('U',0,a,1,e,ip,w,-2,info)
         CALL CHKXER('CHETRI_3',INFot,NOUt,LERr,OK)
!
!        CHETRI_3X
!
         SRNamt = 'CHETRI_3X'
         INFot = 1
         CALL CHETRI_3X('/',0,a,1,e,ip,w,1,info)
         CALL CHKXER('CHETRI_3X',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CHETRI_3X('U',-1,a,1,e,ip,w,1,info)
         CALL CHKXER('CHETRI_3X',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CHETRI_3X('U',2,a,1,e,ip,w,1,info)
         CALL CHKXER('CHETRI_3X',INFot,NOUt,LERr,OK)
!
!        CHETRS_3
!
         SRNamt = 'CHETRS_3'
         INFot = 1
         CALL CHETRS_3('/',0,0,a,1,e,ip,b,1,info)
         CALL CHKXER('CHETRS_3',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CHETRS_3('U',-1,0,a,1,e,ip,b,1,info)
         CALL CHKXER('CHETRS_3',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CHETRS_3('U',0,-1,a,1,e,ip,b,1,info)
         CALL CHKXER('CHETRS_3',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL CHETRS_3('U',2,1,a,1,e,ip,b,2,info)
         CALL CHKXER('CHETRS_3',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL CHETRS_3('U',2,1,a,2,e,ip,b,1,info)
         CALL CHKXER('CHETRS_3',INFot,NOUt,LERr,OK)
!
!        CHECON_3
!
         SRNamt = 'CHECON_3'
         INFot = 1
         CALL CHECON_3('/',0,a,1,e,ip,anrm,rcond,w,info)
         CALL CHKXER('CHECON_3',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CHECON_3('U',-1,a,1,e,ip,anrm,rcond,w,info)
         CALL CHKXER('CHECON_3',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CHECON_3('U',2,a,1,e,ip,anrm,rcond,w,info)
         CALL CHKXER('CHECON_3',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL CHECON_3('U',1,a,1,e,ip,-1.0E0,rcond,w,info)
         CALL CHKXER('CHECON_3',INFot,NOUt,LERr,OK)
!
      ELSEIF ( LSAMEN(2,c2,'HP') ) THEN
!
!     Test error exits of the routines that use factorization
!     of a Hermitian indefinite packed matrix with patrial
!     (Bunch-Kaufman) diagonal pivoting method.
!
!        CHPTRF
!
         SRNamt = 'CHPTRF'
         INFot = 1
         CALL CHPTRF('/',0,a,ip,info)
         CALL CHKXER('CHPTRF',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CHPTRF('U',-1,a,ip,info)
         CALL CHKXER('CHPTRF',INFot,NOUt,LERr,OK)
!
!        CHPTRI
!
         SRNamt = 'CHPTRI'
         INFot = 1
         CALL CHPTRI('/',0,a,ip,w,info)
         CALL CHKXER('CHPTRI',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CHPTRI('U',-1,a,ip,w,info)
         CALL CHKXER('CHPTRI',INFot,NOUt,LERr,OK)
!
!        CHPTRS
!
         SRNamt = 'CHPTRS'
         INFot = 1
         CALL CHPTRS('/',0,0,a,ip,b,1,info)
         CALL CHKXER('CHPTRS',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CHPTRS('U',-1,0,a,ip,b,1,info)
         CALL CHKXER('CHPTRS',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CHPTRS('U',0,-1,a,ip,b,1,info)
         CALL CHKXER('CHPTRS',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL CHPTRS('U',2,1,a,ip,b,1,info)
         CALL CHKXER('CHPTRS',INFot,NOUt,LERr,OK)
!
!        CHPRFS
!
         SRNamt = 'CHPRFS'
         INFot = 1
         CALL CHPRFS('/',0,0,a,af,ip,b,1,x,1,r1,r2,w,r,info)
         CALL CHKXER('CHPRFS',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CHPRFS('U',-1,0,a,af,ip,b,1,x,1,r1,r2,w,r,info)
         CALL CHKXER('CHPRFS',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CHPRFS('U',0,-1,a,af,ip,b,1,x,1,r1,r2,w,r,info)
         CALL CHKXER('CHPRFS',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL CHPRFS('U',2,1,a,af,ip,b,1,x,2,r1,r2,w,r,info)
         CALL CHKXER('CHPRFS',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL CHPRFS('U',2,1,a,af,ip,b,2,x,1,r1,r2,w,r,info)
         CALL CHKXER('CHPRFS',INFot,NOUt,LERr,OK)
!
!        CHPCON
!
         SRNamt = 'CHPCON'
         INFot = 1
         CALL CHPCON('/',0,a,ip,anrm,rcond,w,info)
         CALL CHKXER('CHPCON',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CHPCON('U',-1,a,ip,anrm,rcond,w,info)
         CALL CHKXER('CHPCON',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL CHPCON('U',1,a,ip,-anrm,rcond,w,info)
         CALL CHKXER('CHPCON',INFot,NOUt,LERr,OK)
      ENDIF
!
!     Print a summary line.
!
      CALL ALAESM(Path,OK,NOUt)
!
!
!     End of CERRHE
!
      END SUBROUTINE CERRHE
