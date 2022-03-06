!*==zerrhe.f90  processed by SPAG 7.51RB at 17:35 on  4 Mar 2022
!> \brief \b zerrhe
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE zerrhe( PATH, NUNIT )
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
!> zerrhe tests the error exits for the COMPLEX*16 routines
!> for Hermitian indefinite matrices.
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
!> \date November 2017
!
!> \ingroup complex16_lin
!
!  =====================================================================
      SUBROUTINE zerrhe(Path,Nunit)
      IMPLICIT NONE
!*--ZERRHE59
!
!  -- LAPACK test routine (version 3.8.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2017
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
      CHARACTER*2 c2
      INTEGER i , info , j
      DOUBLE PRECISION anrm , rcond
!     ..
!     .. Local Arrays ..
      INTEGER ip(NMAX)
      DOUBLE PRECISION r(NMAX) , r1(NMAX) , r2(NMAX)
      COMPLEX*16 a(NMAX,NMAX) , af(NMAX,NMAX) , b(NMAX) , e(NMAX) ,     &
     &           w(2*NMAX) , x(NMAX)
!     ..
!     .. External Functions ..
      LOGICAL LSAMEN
EXTERNAL    LSAMEN
!     ..
!     .. External Subroutines ..
EXTERNAL    CHKXER,   ZHECON,   ZHECON_3,   ZHECON_ROOK
EXTERNAL    ZHERFS,   ZHETF2,   ZHETF2_RK,   ZHETF2_ROOK,   ZHETRF
EXTERNAL    ZHETRF_RK,   ZHETRF_ROOK,   ZHETRF_AA,   ZHETRF_AA_2STAGE
EXTERNAL    ZHETRI,   ZHETRI_3,   ZHETRI_3X,   ZHETRI_ROOK,   ZHETRI2
EXTERNAL    ZHETRI2X,   ZHETRS,   ZHETRS_3,   ZHETRS_ROOK,   ZHETRS_AA
EXTERNAL    ZHETRS_AA_2STAGE,   ZHPCON,   ZHPRFS,   ZHPTRF,   ZHPTRI
EXTERNAL    ZHPTRS
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
         ip(j) = j
      ENDDO
      anrm = 1.0D0
      OK = .TRUE.
!
      IF ( LSAMEN(2,c2,'HE') ) THEN
!
!        Test error exits of the routines that use factorization
!        of a Hermitian indefinite matrix with patrial
!        (Bunch-Kaufman) diagonal pivoting method.
!
!        ZHETRF
!
         SRNamt = 'ZHETRF'
         INFot = 1
         CALL ZHETRF('/',0,a,1,ip,w,1,info)
         CALL CHKXER('ZHETRF',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZHETRF('U',-1,a,1,ip,w,1,info)
         CALL CHKXER('ZHETRF',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZHETRF('U',2,a,1,ip,w,4,info)
         CALL CHKXER('ZHETRF',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL ZHETRF('U',0,a,1,ip,w,0,info)
         CALL CHKXER('ZHETRF',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL ZHETRF('U',0,a,1,ip,w,-2,info)
         CALL CHKXER('ZHETRF',INFot,NOUt,LERr,OK)
!
!        ZHETF2
!
         SRNamt = 'ZHETF2'
         INFot = 1
         CALL ZHETF2('/',0,a,1,ip,info)
         CALL CHKXER('ZHETF2',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZHETF2('U',-1,a,1,ip,info)
         CALL CHKXER('ZHETF2',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZHETF2('U',2,a,1,ip,info)
         CALL CHKXER('ZHETF2',INFot,NOUt,LERr,OK)
!
!        ZHETRI
!
         SRNamt = 'ZHETRI'
         INFot = 1
         CALL ZHETRI('/',0,a,1,ip,w,info)
         CALL CHKXER('ZHETRI',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZHETRI('U',-1,a,1,ip,w,info)
         CALL CHKXER('ZHETRI',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZHETRI('U',2,a,1,ip,w,info)
         CALL CHKXER('ZHETRI',INFot,NOUt,LERr,OK)
!
!        ZHETRI2
!
         SRNamt = 'ZHETRI2'
         INFot = 1
         CALL ZHETRI2('/',0,a,1,ip,w,1,info)
         CALL CHKXER('ZHETRI2',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZHETRI2('U',-1,a,1,ip,w,1,info)
         CALL CHKXER('ZHETRI2',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZHETRI2('U',2,a,1,ip,w,1,info)
         CALL CHKXER('ZHETRI2',INFot,NOUt,LERr,OK)
!
!        ZHETRI2X
!
         SRNamt = 'ZHETRI2X'
         INFot = 1
         CALL ZHETRI2X('/',0,a,1,ip,w,1,info)
         CALL CHKXER('ZHETRI2X',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZHETRI2X('U',-1,a,1,ip,w,1,info)
         CALL CHKXER('ZHETRI2X',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZHETRI2X('U',2,a,1,ip,w,1,info)
         CALL CHKXER('ZHETRI2X',INFot,NOUt,LERr,OK)
!
!        ZHETRS
!
         SRNamt = 'ZHETRS'
         INFot = 1
         CALL ZHETRS('/',0,0,a,1,ip,b,1,info)
         CALL CHKXER('ZHETRS',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZHETRS('U',-1,0,a,1,ip,b,1,info)
         CALL CHKXER('ZHETRS',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZHETRS('U',0,-1,a,1,ip,b,1,info)
         CALL CHKXER('ZHETRS',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL ZHETRS('U',2,1,a,1,ip,b,2,info)
         CALL CHKXER('ZHETRS',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL ZHETRS('U',2,1,a,2,ip,b,1,info)
         CALL CHKXER('ZHETRS',INFot,NOUt,LERr,OK)
!
!        ZHERFS
!
         SRNamt = 'ZHERFS'
         INFot = 1
         CALL ZHERFS('/',0,0,a,1,af,1,ip,b,1,x,1,r1,r2,w,r,info)
         CALL CHKXER('ZHERFS',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZHERFS('U',-1,0,a,1,af,1,ip,b,1,x,1,r1,r2,w,r,info)
         CALL CHKXER('ZHERFS',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZHERFS('U',0,-1,a,1,af,1,ip,b,1,x,1,r1,r2,w,r,info)
         CALL CHKXER('ZHERFS',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL ZHERFS('U',2,1,a,1,af,2,ip,b,2,x,2,r1,r2,w,r,info)
         CALL CHKXER('ZHERFS',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL ZHERFS('U',2,1,a,2,af,1,ip,b,2,x,2,r1,r2,w,r,info)
         CALL CHKXER('ZHERFS',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL ZHERFS('U',2,1,a,2,af,2,ip,b,1,x,2,r1,r2,w,r,info)
         CALL CHKXER('ZHERFS',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL ZHERFS('U',2,1,a,2,af,2,ip,b,2,x,1,r1,r2,w,r,info)
         CALL CHKXER('ZHERFS',INFot,NOUt,LERr,OK)
!
!        ZHECON
!
         SRNamt = 'ZHECON'
         INFot = 1
         CALL ZHECON('/',0,a,1,ip,anrm,rcond,w,info)
         CALL CHKXER('ZHECON',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZHECON('U',-1,a,1,ip,anrm,rcond,w,info)
         CALL CHKXER('ZHECON',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZHECON('U',2,a,1,ip,anrm,rcond,w,info)
         CALL CHKXER('ZHECON',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL ZHECON('U',1,a,1,ip,-anrm,rcond,w,info)
         CALL CHKXER('ZHECON',INFot,NOUt,LERr,OK)
!
      ELSEIF ( LSAMEN(2,c2,'HR') ) THEN
!
!        Test error exits of the routines that use factorization
!        of a Hermitian indefinite matrix with rook
!        (bounded Bunch-Kaufman) diagonal pivoting method.
!
!        ZHETRF_ROOK
!
         SRNamt = 'ZHETRF_ROOK'
         INFot = 1
         CALL ZHETRF_ROOK('/',0,a,1,ip,w,1,info)
         CALL CHKXER('ZHETRF_ROOK',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZHETRF_ROOK('U',-1,a,1,ip,w,1,info)
         CALL CHKXER('ZHETRF_ROOK',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZHETRF_ROOK('U',2,a,1,ip,w,4,info)
         CALL CHKXER('ZHETRF_ROOK',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL ZHETRF_ROOK('U',0,a,1,ip,w,0,info)
         CALL CHKXER('ZHETRF_ROOK',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL ZHETRF_ROOK('U',0,a,1,ip,w,-2,info)
         CALL CHKXER('ZHETRF_ROOK',INFot,NOUt,LERr,OK)
!
!        ZHETF2_ROOK
!
         SRNamt = 'ZHETF2_ROOK'
         INFot = 1
         CALL ZHETF2_ROOK('/',0,a,1,ip,info)
         CALL CHKXER('ZHETF2_ROOK',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZHETF2_ROOK('U',-1,a,1,ip,info)
         CALL CHKXER('ZHETF2_ROOK',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZHETF2_ROOK('U',2,a,1,ip,info)
         CALL CHKXER('ZHETF2_ROOK',INFot,NOUt,LERr,OK)
!
!        ZHETRI_ROOK
!
         SRNamt = 'ZHETRI_ROOK'
         INFot = 1
         CALL ZHETRI_ROOK('/',0,a,1,ip,w,info)
         CALL CHKXER('ZHETRI_ROOK',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZHETRI_ROOK('U',-1,a,1,ip,w,info)
         CALL CHKXER('ZHETRI_ROOK',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZHETRI_ROOK('U',2,a,1,ip,w,info)
         CALL CHKXER('ZHETRI_ROOK',INFot,NOUt,LERr,OK)
!
!        ZHETRS_ROOK
!
         SRNamt = 'ZHETRS_ROOK'
         INFot = 1
         CALL ZHETRS_ROOK('/',0,0,a,1,ip,b,1,info)
         CALL CHKXER('ZHETRS_ROOK',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZHETRS_ROOK('U',-1,0,a,1,ip,b,1,info)
         CALL CHKXER('ZHETRS_ROOK',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZHETRS_ROOK('U',0,-1,a,1,ip,b,1,info)
         CALL CHKXER('ZHETRS_ROOK',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL ZHETRS_ROOK('U',2,1,a,1,ip,b,2,info)
         CALL CHKXER('ZHETRS_ROOK',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL ZHETRS_ROOK('U',2,1,a,2,ip,b,1,info)
         CALL CHKXER('ZHETRS_ROOK',INFot,NOUt,LERr,OK)
!
!        ZHECON_ROOK
!
         SRNamt = 'ZHECON_ROOK'
         INFot = 1
         CALL ZHECON_ROOK('/',0,a,1,ip,anrm,rcond,w,info)
         CALL CHKXER('ZHECON_ROOK',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZHECON_ROOK('U',-1,a,1,ip,anrm,rcond,w,info)
         CALL CHKXER('ZHECON_ROOK',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZHECON_ROOK('U',2,a,1,ip,anrm,rcond,w,info)
         CALL CHKXER('ZHECON_ROOK',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL ZHECON_ROOK('U',1,a,1,ip,-anrm,rcond,w,info)
         CALL CHKXER('ZHECON_ROOK',INFot,NOUt,LERr,OK)
!
      ELSEIF ( LSAMEN(2,c2,'HK') ) THEN
!
!        Test error exits of the routines that use factorization
!        of a symmetric indefinite matrix with rook
!        (bounded Bunch-Kaufman) pivoting with the new storage
!        format for factors L ( or U) and D.
!
!        L (or U) is stored in A, diagonal of D is stored on the
!        diagonal of A, subdiagonal of D is stored in a separate array E.
!
!        ZHETRF_RK
!
         SRNamt = 'ZHETRF_RK'
         INFot = 1
         CALL ZHETRF_RK('/',0,a,1,e,ip,w,1,info)
         CALL CHKXER('ZHETRF_RK',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZHETRF_RK('U',-1,a,1,e,ip,w,1,info)
         CALL CHKXER('ZHETRF_RK',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZHETRF_RK('U',2,a,1,e,ip,w,4,info)
         CALL CHKXER('ZHETRF_RK',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL ZHETRF_RK('U',0,a,1,e,ip,w,0,info)
         CALL CHKXER('ZHETRF_RK',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL ZHETRF_RK('U',0,a,1,e,ip,w,-2,info)
         CALL CHKXER('ZHETRF_RK',INFot,NOUt,LERr,OK)
!
!        ZHETF2_RK
!
         SRNamt = 'ZHETF2_RK'
         INFot = 1
         CALL ZHETF2_RK('/',0,a,1,e,ip,info)
         CALL CHKXER('ZHETF2_RK',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZHETF2_RK('U',-1,a,1,e,ip,info)
         CALL CHKXER('ZHETF2_RK',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZHETF2_RK('U',2,a,1,e,ip,info)
         CALL CHKXER('ZHETF2_RK',INFot,NOUt,LERr,OK)
!
!        ZHETRI_3
!
         SRNamt = 'ZHETRI_3'
         INFot = 1
         CALL ZHETRI_3('/',0,a,1,e,ip,w,1,info)
         CALL CHKXER('ZHETRI_3',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZHETRI_3('U',-1,a,1,e,ip,w,1,info)
         CALL CHKXER('ZHETRI_3',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZHETRI_3('U',2,a,1,e,ip,w,1,info)
         CALL CHKXER('ZHETRI_3',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL ZHETRI_3('U',0,a,1,e,ip,w,0,info)
         CALL CHKXER('ZHETRI_3',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL ZHETRI_3('U',0,a,1,e,ip,w,-2,info)
         CALL CHKXER('ZHETRI_3',INFot,NOUt,LERr,OK)
!
!        ZHETRI_3X
!
         SRNamt = 'ZHETRI_3X'
         INFot = 1
         CALL ZHETRI_3X('/',0,a,1,e,ip,w,1,info)
         CALL CHKXER('ZHETRI_3X',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZHETRI_3X('U',-1,a,1,e,ip,w,1,info)
         CALL CHKXER('ZHETRI_3X',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZHETRI_3X('U',2,a,1,e,ip,w,1,info)
         CALL CHKXER('ZHETRI_3X',INFot,NOUt,LERr,OK)
!
!        ZHETRS_3
!
         SRNamt = 'ZHETRS_3'
         INFot = 1
         CALL ZHETRS_3('/',0,0,a,1,e,ip,b,1,info)
         CALL CHKXER('ZHETRS_3',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZHETRS_3('U',-1,0,a,1,e,ip,b,1,info)
         CALL CHKXER('ZHETRS_3',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZHETRS_3('U',0,-1,a,1,e,ip,b,1,info)
         CALL CHKXER('ZHETRS_3',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL ZHETRS_3('U',2,1,a,1,e,ip,b,2,info)
         CALL CHKXER('ZHETRS_3',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL ZHETRS_3('U',2,1,a,2,e,ip,b,1,info)
         CALL CHKXER('ZHETRS_3',INFot,NOUt,LERr,OK)
!
!        ZHECON_3
!
         SRNamt = 'ZHECON_3'
         INFot = 1
         CALL ZHECON_3('/',0,a,1,e,ip,anrm,rcond,w,info)
         CALL CHKXER('ZHECON_3',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZHECON_3('U',-1,a,1,e,ip,anrm,rcond,w,info)
         CALL CHKXER('ZHECON_3',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZHECON_3('U',2,a,1,e,ip,anrm,rcond,w,info)
         CALL CHKXER('ZHECON_3',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL ZHECON_3('U',1,a,1,e,ip,-1.0D0,rcond,w,info)
         CALL CHKXER('ZHECON_3',INFot,NOUt,LERr,OK)
!
!        Test error exits of the routines that use factorization
!        of a Hermitian indefinite matrix with Aasen's algorithm.
!
      ELSEIF ( LSAMEN(2,c2,'HA') ) THEN
!
!        ZHETRF_AA
!
         SRNamt = 'ZHETRF_AA'
         INFot = 1
         CALL ZHETRF_AA('/',0,a,1,ip,w,1,info)
         CALL CHKXER('ZHETRF_AA',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZHETRF_AA('U',-1,a,1,ip,w,1,info)
         CALL CHKXER('ZHETRF_AA',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZHETRF_AA('U',2,a,1,ip,w,4,info)
         CALL CHKXER('ZHETRF_AA',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL ZHETRF_AA('U',0,a,1,ip,w,0,info)
         CALL CHKXER('ZHETRF_AA',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL ZHETRF_AA('U',0,a,1,ip,w,-2,info)
         CALL CHKXER('ZHETRF_AA',INFot,NOUt,LERr,OK)
!
!        ZHETRS_AA
!
         SRNamt = 'ZHETRS_AA'
         INFot = 1
         CALL ZHETRS_AA('/',0,0,a,1,ip,b,1,w,1,info)
         CALL CHKXER('ZHETRS_AA',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZHETRS_AA('U',-1,0,a,1,ip,b,1,w,1,info)
         CALL CHKXER('ZHETRS_AA',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZHETRS_AA('U',0,-1,a,1,ip,b,1,w,1,info)
         CALL CHKXER('ZHETRS_AA',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL ZHETRS_AA('U',2,1,a,1,ip,b,2,w,1,info)
         CALL CHKXER('ZHETRS_AA',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL ZHETRS_AA('U',2,1,a,2,ip,b,1,w,1,info)
         CALL CHKXER('ZHETRS_AA',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL ZHETRS_AA('U',0,1,a,1,ip,b,1,w,0,info)
         CALL CHKXER('ZHETRS_AA',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL ZHETRS_AA('U',0,1,a,1,ip,b,1,w,-2,info)
         CALL CHKXER('ZHETRS_AA',INFot,NOUt,LERr,OK)
!
      ELSEIF ( LSAMEN(2,c2,'S2') ) THEN
!
!        Test error exits of the routines that use factorization
!        of a symmetric indefinite matrix with Aasen's algorithm.
!
!        ZHETRF_AA_2STAGE
!
         SRNamt = 'ZHETRF_AA_2STAGE'
         INFot = 1
         CALL ZHETRF_AA_2STAGE('/',0,a,1,a,1,ip,ip,w,1,info)
         CALL CHKXER('ZHETRF_AA_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZHETRF_AA_2STAGE('U',-1,a,1,a,1,ip,ip,w,1,info)
         CALL CHKXER('ZHETRF_AA_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZHETRF_AA_2STAGE('U',2,a,1,a,2,ip,ip,w,1,info)
         CALL CHKXER('ZHETRF_AA_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL ZHETRF_AA_2STAGE('U',2,a,2,a,1,ip,ip,w,1,info)
         CALL CHKXER('ZHETRF_AA_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL ZHETRF_AA_2STAGE('U',2,a,2,a,8,ip,ip,w,0,info)
         CALL CHKXER('ZHETRF_AA_2STAGE',INFot,NOUt,LERr,OK)
!
!        ZHETRS_AA_2STAGE
!
         SRNamt = 'ZHETRS_AA_2STAGE'
         INFot = 1
         CALL ZHETRS_AA_2STAGE('/',0,0,a,1,a,1,ip,ip,b,1,info)
         CALL CHKXER('ZHETRS_AA_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZHETRS_AA_2STAGE('U',-1,0,a,1,a,1,ip,ip,b,1,info)
         CALL CHKXER('ZHETRS_AA_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZHETRS_AA_2STAGE('U',0,-1,a,1,a,1,ip,ip,b,1,info)
         CALL CHKXER('ZHETRS_AA_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL ZHETRS_AA_2STAGE('U',2,1,a,1,a,1,ip,ip,b,1,info)
         CALL CHKXER('ZHETRS_AA_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL ZHETRS_AA_2STAGE('U',2,1,a,2,a,1,ip,ip,b,1,info)
         CALL CHKXER('ZHETRS_AA_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL ZHETRS_AA_2STAGE('U',2,1,a,2,a,8,ip,ip,b,1,info)
         CALL CHKXER('ZHETRS_AA_STAGE',INFot,NOUt,LERr,OK)
!
      ELSEIF ( LSAMEN(2,c2,'HP') ) THEN
!
!        Test error exits of the routines that use factorization
!        of a Hermitian indefinite packed matrix with patrial
!        (Bunch-Kaufman) diagonal pivoting method.
!
!        ZHPTRF
!
         SRNamt = 'ZHPTRF'
         INFot = 1
         CALL ZHPTRF('/',0,a,ip,info)
         CALL CHKXER('ZHPTRF',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZHPTRF('U',-1,a,ip,info)
         CALL CHKXER('ZHPTRF',INFot,NOUt,LERr,OK)
!
!        ZHPTRI
!
         SRNamt = 'ZHPTRI'
         INFot = 1
         CALL ZHPTRI('/',0,a,ip,w,info)
         CALL CHKXER('ZHPTRI',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZHPTRI('U',-1,a,ip,w,info)
         CALL CHKXER('ZHPTRI',INFot,NOUt,LERr,OK)
!
!        ZHPTRS
!
         SRNamt = 'ZHPTRS'
         INFot = 1
         CALL ZHPTRS('/',0,0,a,ip,b,1,info)
         CALL CHKXER('ZHPTRS',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZHPTRS('U',-1,0,a,ip,b,1,info)
         CALL CHKXER('ZHPTRS',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZHPTRS('U',0,-1,a,ip,b,1,info)
         CALL CHKXER('ZHPTRS',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL ZHPTRS('U',2,1,a,ip,b,1,info)
         CALL CHKXER('ZHPTRS',INFot,NOUt,LERr,OK)
!
!        ZHPRFS
!
         SRNamt = 'ZHPRFS'
         INFot = 1
         CALL ZHPRFS('/',0,0,a,af,ip,b,1,x,1,r1,r2,w,r,info)
         CALL CHKXER('ZHPRFS',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZHPRFS('U',-1,0,a,af,ip,b,1,x,1,r1,r2,w,r,info)
         CALL CHKXER('ZHPRFS',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZHPRFS('U',0,-1,a,af,ip,b,1,x,1,r1,r2,w,r,info)
         CALL CHKXER('ZHPRFS',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL ZHPRFS('U',2,1,a,af,ip,b,1,x,2,r1,r2,w,r,info)
         CALL CHKXER('ZHPRFS',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL ZHPRFS('U',2,1,a,af,ip,b,2,x,1,r1,r2,w,r,info)
         CALL CHKXER('ZHPRFS',INFot,NOUt,LERr,OK)
!
!        ZHPCON
!
         SRNamt = 'ZHPCON'
         INFot = 1
         CALL ZHPCON('/',0,a,ip,anrm,rcond,w,info)
         CALL CHKXER('ZHPCON',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZHPCON('U',-1,a,ip,anrm,rcond,w,info)
         CALL CHKXER('ZHPCON',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL ZHPCON('U',1,a,ip,-anrm,rcond,w,info)
         CALL CHKXER('ZHPCON',INFot,NOUt,LERr,OK)
      ENDIF
!
!     Print a summary line.
!
      CALL alaesm(Path,OK,NOUt)
!
!
!     End of zerrhe
!
      END SUBROUTINE zerrhe
