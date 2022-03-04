!*==cerrsy.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b CERRSY
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE CERRSY( PATH, NUNIT )
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
!> CERRSY tests the error exits for the COMPLEX routines
!> for symmetric indefinite matrices.
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
!> \ingroup complex_lin
!
!  =====================================================================
      SUBROUTINE CERRSY(Path,Nunit)
      IMPLICIT NONE
!*--CERRSY59
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
!     .. Parameters ..
      INTEGER NMAX
      PARAMETER (NMAX=4)
!     ..
!     .. Local Scalars ..
      CHARACTER*2 c2
      INTEGER i , info , j
      REAL anrm , rcond
!     ..
!     .. Local Arrays ..
      INTEGER ip(NMAX)
      REAL r(NMAX) , r1(NMAX) , r2(NMAX)
      COMPLEX a(NMAX,NMAX) , af(NMAX,NMAX) , b(NMAX) , e(NMAX) ,        &
     &        w(2*NMAX) , x(NMAX)
!     ..
!     .. External Functions ..
      LOGICAL LSAMEN
      EXTERNAL LSAMEN
!     ..
!     .. External Subroutines ..
      EXTERNAL ALAESM , CHKXER , CSPCON , CSPRFS , CSPTRF , CSPTRI ,    &
     &         CSPTRS , CSYCON , CSYCON_3 , CSYCON_ROOK , CSYRFS ,      &
     &         CSYTF2 , CSYTF2_RK , CSYTF2_ROOK , CSYTRF , CSYTRF_RK ,  &
     &         CSYTRF_ROOK , CSYTRI , CSYTRI_3 , CSYTRI_3X ,            &
     &         CSYTRI_ROOK , CSYTRI2 , CSYTRI2X , CSYTRS , CSYTRS_3 ,   &
     &         CSYTRS_ROOK
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
         b(j) = 0.E0
         e(j) = 0.E0
         r1(j) = 0.E0
         r2(j) = 0.E0
         w(j) = 0.E0
         x(j) = 0.E0
         ip(j) = j
      ENDDO
      anrm = 1.0
      OK = .TRUE.
!
      IF ( LSAMEN(2,c2,'SY') ) THEN
!
!        Test error exits of the routines that use factorization
!        of a symmetric indefinite matrix with patrial
!        (Bunch-Kaufman) diagonal pivoting method.
!
!        CSYTRF
!
         SRNamt = 'CSYTRF'
         INFot = 1
         CALL CSYTRF('/',0,a,1,ip,w,1,info)
         CALL CHKXER('CSYTRF',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CSYTRF('U',-1,a,1,ip,w,1,info)
         CALL CHKXER('CSYTRF',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CSYTRF('U',2,a,1,ip,w,4,info)
         CALL CHKXER('CSYTRF',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL CSYTRF('U',0,a,1,ip,w,0,info)
         CALL CHKXER('CSYTRF',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL CSYTRF('U',0,a,1,ip,w,-2,info)
         CALL CHKXER('CSYTRF',INFot,NOUt,LERr,OK)
!
!        CSYTF2
!
         SRNamt = 'CSYTF2'
         INFot = 1
         CALL CSYTF2('/',0,a,1,ip,info)
         CALL CHKXER('CSYTF2',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CSYTF2('U',-1,a,1,ip,info)
         CALL CHKXER('CSYTF2',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CSYTF2('U',2,a,1,ip,info)
         CALL CHKXER('CSYTF2',INFot,NOUt,LERr,OK)
!
!        CSYTRI
!
         SRNamt = 'CSYTRI'
         INFot = 1
         CALL CSYTRI('/',0,a,1,ip,w,info)
         CALL CHKXER('CSYTRI',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CSYTRI('U',-1,a,1,ip,w,info)
         CALL CHKXER('CSYTRI',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CSYTRI('U',2,a,1,ip,w,info)
         CALL CHKXER('CSYTRI',INFot,NOUt,LERr,OK)
!
!        CSYTRI2
!
         SRNamt = 'CSYTRI2'
         INFot = 1
         CALL CSYTRI2('/',0,a,1,ip,w,1,info)
         CALL CHKXER('CSYTRI2',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CSYTRI2('U',-1,a,1,ip,w,1,info)
         CALL CHKXER('CSYTRI2',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CSYTRI2('U',2,a,1,ip,w,1,info)
         CALL CHKXER('CSYTRI2',INFot,NOUt,LERr,OK)
!
!        CSYTRI2X
!
         SRNamt = 'CSYTRI2X'
         INFot = 1
         CALL CSYTRI2X('/',0,a,1,ip,w,1,info)
         CALL CHKXER('CSYTRI2X',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CSYTRI2X('U',-1,a,1,ip,w,1,info)
         CALL CHKXER('CSYTRI2X',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CSYTRI2X('U',2,a,1,ip,w,1,info)
         CALL CHKXER('CSYTRI2X',INFot,NOUt,LERr,OK)
!
!        CSYTRS
!
         SRNamt = 'CSYTRS'
         INFot = 1
         CALL CSYTRS('/',0,0,a,1,ip,b,1,info)
         CALL CHKXER('CSYTRS',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CSYTRS('U',-1,0,a,1,ip,b,1,info)
         CALL CHKXER('CSYTRS',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CSYTRS('U',0,-1,a,1,ip,b,1,info)
         CALL CHKXER('CSYTRS',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL CSYTRS('U',2,1,a,1,ip,b,2,info)
         CALL CHKXER('CSYTRS',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL CSYTRS('U',2,1,a,2,ip,b,1,info)
         CALL CHKXER('CSYTRS',INFot,NOUt,LERr,OK)
!
!        CSYRFS
!
         SRNamt = 'CSYRFS'
         INFot = 1
         CALL CSYRFS('/',0,0,a,1,af,1,ip,b,1,x,1,r1,r2,w,r,info)
         CALL CHKXER('CSYRFS',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CSYRFS('U',-1,0,a,1,af,1,ip,b,1,x,1,r1,r2,w,r,info)
         CALL CHKXER('CSYRFS',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CSYRFS('U',0,-1,a,1,af,1,ip,b,1,x,1,r1,r2,w,r,info)
         CALL CHKXER('CSYRFS',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL CSYRFS('U',2,1,a,1,af,2,ip,b,2,x,2,r1,r2,w,r,info)
         CALL CHKXER('CSYRFS',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL CSYRFS('U',2,1,a,2,af,1,ip,b,2,x,2,r1,r2,w,r,info)
         CALL CHKXER('CSYRFS',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL CSYRFS('U',2,1,a,2,af,2,ip,b,1,x,2,r1,r2,w,r,info)
         CALL CHKXER('CSYRFS',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL CSYRFS('U',2,1,a,2,af,2,ip,b,2,x,1,r1,r2,w,r,info)
         CALL CHKXER('CSYRFS',INFot,NOUt,LERr,OK)
!
!        CSYCON
!
         SRNamt = 'CSYCON'
         INFot = 1
         CALL CSYCON('/',0,a,1,ip,anrm,rcond,w,info)
         CALL CHKXER('CSYCON',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CSYCON('U',-1,a,1,ip,anrm,rcond,w,info)
         CALL CHKXER('CSYCON',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CSYCON('U',2,a,1,ip,anrm,rcond,w,info)
         CALL CHKXER('CSYCON',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL CSYCON('U',1,a,1,ip,-anrm,rcond,w,info)
         CALL CHKXER('CSYCON',INFot,NOUt,LERr,OK)
!
      ELSEIF ( LSAMEN(2,c2,'SR') ) THEN
!
!        Test error exits of the routines that use factorization
!        of a symmetric indefinite matrix with rook
!        (bounded Bunch-Kaufman) diagonal pivoting method.
!
!        CSYTRF_ROOK
!
         SRNamt = 'CSYTRF_ROOK'
         INFot = 1
         CALL CSYTRF_ROOK('/',0,a,1,ip,w,1,info)
         CALL CHKXER('CSYTRF_ROOK',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CSYTRF_ROOK('U',-1,a,1,ip,w,1,info)
         CALL CHKXER('CSYTRF_ROOK',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CSYTRF_ROOK('U',2,a,1,ip,w,4,info)
         CALL CHKXER('CSYTRF_ROOK',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL CSYTRF_ROOK('U',0,a,1,ip,w,0,info)
         CALL CHKXER('CSYTRF_ROOK',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL CSYTRF_ROOK('U',0,a,1,ip,w,-2,info)
         CALL CHKXER('CSYTRF_ROOK',INFot,NOUt,LERr,OK)
!
!        CSYTF2_ROOK
!
         SRNamt = 'CSYTF2_ROOK'
         INFot = 1
         CALL CSYTF2_ROOK('/',0,a,1,ip,info)
         CALL CHKXER('CSYTF2_ROOK',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CSYTF2_ROOK('U',-1,a,1,ip,info)
         CALL CHKXER('CSYTF2_ROOK',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CSYTF2_ROOK('U',2,a,1,ip,info)
         CALL CHKXER('CSYTF2_ROOK',INFot,NOUt,LERr,OK)
!
!        CSYTRI_ROOK
!
         SRNamt = 'CSYTRI_ROOK'
         INFot = 1
         CALL CSYTRI_ROOK('/',0,a,1,ip,w,info)
         CALL CHKXER('CSYTRI_ROOK',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CSYTRI_ROOK('U',-1,a,1,ip,w,info)
         CALL CHKXER('CSYTRI_ROOK',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CSYTRI_ROOK('U',2,a,1,ip,w,info)
         CALL CHKXER('CSYTRI_ROOK',INFot,NOUt,LERr,OK)
!
!        CSYTRS_ROOK
!
         SRNamt = 'CSYTRS_ROOK'
         INFot = 1
         CALL CSYTRS_ROOK('/',0,0,a,1,ip,b,1,info)
         CALL CHKXER('CSYTRS_ROOK',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CSYTRS_ROOK('U',-1,0,a,1,ip,b,1,info)
         CALL CHKXER('CSYTRS_ROOK',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CSYTRS_ROOK('U',0,-1,a,1,ip,b,1,info)
         CALL CHKXER('CSYTRS_ROOK',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL CSYTRS_ROOK('U',2,1,a,1,ip,b,2,info)
         CALL CHKXER('CSYTRS_ROOK',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL CSYTRS_ROOK('U',2,1,a,2,ip,b,1,info)
         CALL CHKXER('CSYTRS_ROOK',INFot,NOUt,LERr,OK)
!
!        CSYCON_ROOK
!
         SRNamt = 'CSYCON_ROOK'
         INFot = 1
         CALL CSYCON_ROOK('/',0,a,1,ip,anrm,rcond,w,info)
         CALL CHKXER('CSYCON_ROOK',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CSYCON_ROOK('U',-1,a,1,ip,anrm,rcond,w,info)
         CALL CHKXER('CSYCON_ROOK',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CSYCON_ROOK('U',2,a,1,ip,anrm,rcond,w,info)
         CALL CHKXER('CSYCON_ROOK',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL CSYCON_ROOK('U',1,a,1,ip,-anrm,rcond,w,info)
         CALL CHKXER('CSYCON_ROOK',INFot,NOUt,LERr,OK)
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
!        CSYTRF_RK
!
         SRNamt = 'CSYTRF_RK'
         INFot = 1
         CALL CSYTRF_RK('/',0,a,1,e,ip,w,1,info)
         CALL CHKXER('CSYTRF_RK',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CSYTRF_RK('U',-1,a,1,e,ip,w,1,info)
         CALL CHKXER('CSYTRF_RK',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CSYTRF_RK('U',2,a,1,e,ip,w,4,info)
         CALL CHKXER('CSYTRF_RK',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL CSYTRF_RK('U',0,a,1,e,ip,w,0,info)
         CALL CHKXER('CSYTRF_RK',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL CSYTRF_RK('U',0,a,1,e,ip,w,-2,info)
         CALL CHKXER('CSYTRF_RK',INFot,NOUt,LERr,OK)
!
!        CSYTF2_RK
!
         SRNamt = 'CSYTF2_RK'
         INFot = 1
         CALL CSYTF2_RK('/',0,a,1,e,ip,info)
         CALL CHKXER('CSYTF2_RK',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CSYTF2_RK('U',-1,a,1,e,ip,info)
         CALL CHKXER('CSYTF2_RK',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CSYTF2_RK('U',2,a,1,e,ip,info)
         CALL CHKXER('CSYTF2_RK',INFot,NOUt,LERr,OK)
!
!        CSYTRI_3
!
         SRNamt = 'CSYTRI_3'
         INFot = 1
         CALL CSYTRI_3('/',0,a,1,e,ip,w,1,info)
         CALL CHKXER('CSYTRI_3',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CSYTRI_3('U',-1,a,1,e,ip,w,1,info)
         CALL CHKXER('CSYTRI_3',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CSYTRI_3('U',2,a,1,e,ip,w,1,info)
         CALL CHKXER('CSYTRI_3',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL CSYTRI_3('U',0,a,1,e,ip,w,0,info)
         CALL CHKXER('CSYTRI_3',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL CSYTRI_3('U',0,a,1,e,ip,w,-2,info)
         CALL CHKXER('CSYTRI_3',INFot,NOUt,LERr,OK)
!
!        CSYTRI_3X
!
         SRNamt = 'CSYTRI_3X'
         INFot = 1
         CALL CSYTRI_3X('/',0,a,1,e,ip,w,1,info)
         CALL CHKXER('CSYTRI_3X',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CSYTRI_3X('U',-1,a,1,e,ip,w,1,info)
         CALL CHKXER('CSYTRI_3X',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CSYTRI_3X('U',2,a,1,e,ip,w,1,info)
         CALL CHKXER('CSYTRI_3X',INFot,NOUt,LERr,OK)
!
!        CSYTRS_3
!
         SRNamt = 'CSYTRS_3'
         INFot = 1
         CALL CSYTRS_3('/',0,0,a,1,e,ip,b,1,info)
         CALL CHKXER('CSYTRS_3',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CSYTRS_3('U',-1,0,a,1,e,ip,b,1,info)
         CALL CHKXER('CSYTRS_3',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CSYTRS_3('U',0,-1,a,1,e,ip,b,1,info)
         CALL CHKXER('CSYTRS_3',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL CSYTRS_3('U',2,1,a,1,e,ip,b,2,info)
         CALL CHKXER('CSYTRS_3',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL CSYTRS_3('U',2,1,a,2,e,ip,b,1,info)
         CALL CHKXER('CSYTRS_3',INFot,NOUt,LERr,OK)
!
!        CSYCON_3
!
         SRNamt = 'CSYCON_3'
         INFot = 1
         CALL CSYCON_3('/',0,a,1,e,ip,anrm,rcond,w,info)
         CALL CHKXER('CSYCON_3',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CSYCON_3('U',-1,a,1,e,ip,anrm,rcond,w,info)
         CALL CHKXER('CSYCON_3',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CSYCON_3('U',2,a,1,e,ip,anrm,rcond,w,info)
         CALL CHKXER('CSYCON_3',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL CSYCON_3('U',1,a,1,e,ip,-1.0E0,rcond,w,info)
         CALL CHKXER('CSYCON_3',INFot,NOUt,LERr,OK)
!
      ELSEIF ( LSAMEN(2,c2,'SP') ) THEN
!
!        Test error exits of the routines that use factorization
!        of a symmetric indefinite packed matrix with patrial
!        (Bunch-Kaufman) diagonal pivoting method.
!
!        CSPTRF
!
         SRNamt = 'CSPTRF'
         INFot = 1
         CALL CSPTRF('/',0,a,ip,info)
         CALL CHKXER('CSPTRF',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CSPTRF('U',-1,a,ip,info)
         CALL CHKXER('CSPTRF',INFot,NOUt,LERr,OK)
!
!        CSPTRI
!
         SRNamt = 'CSPTRI'
         INFot = 1
         CALL CSPTRI('/',0,a,ip,w,info)
         CALL CHKXER('CSPTRI',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CSPTRI('U',-1,a,ip,w,info)
         CALL CHKXER('CSPTRI',INFot,NOUt,LERr,OK)
!
!        CSPTRS
!
         SRNamt = 'CSPTRS'
         INFot = 1
         CALL CSPTRS('/',0,0,a,ip,b,1,info)
         CALL CHKXER('CSPTRS',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CSPTRS('U',-1,0,a,ip,b,1,info)
         CALL CHKXER('CSPTRS',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CSPTRS('U',0,-1,a,ip,b,1,info)
         CALL CHKXER('CSPTRS',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL CSPTRS('U',2,1,a,ip,b,1,info)
         CALL CHKXER('CSPTRS',INFot,NOUt,LERr,OK)
!
!        CSPRFS
!
         SRNamt = 'CSPRFS'
         INFot = 1
         CALL CSPRFS('/',0,0,a,af,ip,b,1,x,1,r1,r2,w,r,info)
         CALL CHKXER('CSPRFS',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CSPRFS('U',-1,0,a,af,ip,b,1,x,1,r1,r2,w,r,info)
         CALL CHKXER('CSPRFS',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CSPRFS('U',0,-1,a,af,ip,b,1,x,1,r1,r2,w,r,info)
         CALL CHKXER('CSPRFS',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL CSPRFS('U',2,1,a,af,ip,b,1,x,2,r1,r2,w,r,info)
         CALL CHKXER('CSPRFS',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL CSPRFS('U',2,1,a,af,ip,b,2,x,1,r1,r2,w,r,info)
         CALL CHKXER('CSPRFS',INFot,NOUt,LERr,OK)
!
!        CSPCON
!
         SRNamt = 'CSPCON'
         INFot = 1
         CALL CSPCON('/',0,a,ip,anrm,rcond,w,info)
         CALL CHKXER('CSPCON',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CSPCON('U',-1,a,ip,anrm,rcond,w,info)
         CALL CHKXER('CSPCON',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL CSPCON('U',1,a,ip,-anrm,rcond,w,info)
         CALL CHKXER('CSPCON',INFot,NOUt,LERr,OK)
!
      ELSEIF ( LSAMEN(2,c2,'SA') ) THEN
!
!        Test error exits of the routines that use factorization
!        of a symmetric indefinite matrix with Aasen's algorithm
!
!        CSYTRF_AA
!
         SRNamt = 'CSYTRF_AA'
         INFot = 1
         CALL CSYTRF_AA('/',0,a,1,ip,w,1,info)
         CALL CHKXER('CSYTRF_AA',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CSYTRF_AA('U',-1,a,1,ip,w,1,info)
         CALL CHKXER('CSYTRF_AA',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CSYTRF_AA('U',2,a,1,ip,w,4,info)
         CALL CHKXER('CSYTRF_AA',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL CSYTRF_AA('U',0,a,1,ip,w,0,info)
         CALL CHKXER('CSYTRF_AA',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL CSYTRF_AA('U',0,a,1,ip,w,-2,info)
         CALL CHKXER('CSYTRF_AA',INFot,NOUt,LERr,OK)
!
!        CSYTRS_AA
!
         SRNamt = 'CSYTRS_AA'
         INFot = 1
         CALL CSYTRS_AA('/',0,0,a,1,ip,b,1,w,1,info)
         CALL CHKXER('CSYTRS_AA',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CSYTRS_AA('U',-1,0,a,1,ip,b,1,w,1,info)
         CALL CHKXER('CSYTRS_AA',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CSYTRS_AA('U',0,-1,a,1,ip,b,1,w,1,info)
         CALL CHKXER('CSYTRS_AA',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL CSYTRS_AA('U',2,1,a,1,ip,b,2,w,1,info)
         CALL CHKXER('CSYTRS_AA',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL CSYTRS_AA('U',2,1,a,2,ip,b,1,w,1,info)
         CALL CHKXER('CSYTRS_AA',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL CSYTRS_AA('U',0,1,a,1,ip,b,1,w,0,info)
         CALL CHKXER('CSYTRS_AA',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL CSYTRS_AA('U',0,1,a,1,ip,b,1,w,-2,info)
         CALL CHKXER('CSYTRS_AA',INFot,NOUt,LERr,OK)
!
      ELSEIF ( LSAMEN(2,c2,'S2') ) THEN
!
!        Test error exits of the routines that use factorization
!        of a symmetric indefinite matrix with Aasen's algorithm.
!
!        CSYTRF_AA_2STAGE
!
         SRNamt = 'CSYTRF_AA_2STAGE'
         INFot = 1
         CALL CSYTRF_AA_2STAGE('/',0,a,1,a,1,ip,ip,w,1,info)
         CALL CHKXER('CSYTRF_AA_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CSYTRF_AA_2STAGE('U',-1,a,1,a,1,ip,ip,w,1,info)
         CALL CHKXER('CSYTRF_AA_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CSYTRF_AA_2STAGE('U',2,a,1,a,2,ip,ip,w,1,info)
         CALL CHKXER('CSYTRF_AA_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL CSYTRF_AA_2STAGE('U',2,a,2,a,1,ip,ip,w,1,info)
         CALL CHKXER('CSYTRF_AA_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL CSYTRF_AA_2STAGE('U',2,a,2,a,8,ip,ip,w,0,info)
         CALL CHKXER('CSYTRF_AA_2STAGE',INFot,NOUt,LERr,OK)
!
!        CHETRS_AA_2STAGE
!
         SRNamt = 'CSYTRS_AA_2STAGE'
         INFot = 1
         CALL CSYTRS_AA_2STAGE('/',0,0,a,1,a,1,ip,ip,b,1,info)
         CALL CHKXER('CSYTRS_AA_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CSYTRS_AA_2STAGE('U',-1,0,a,1,a,1,ip,ip,b,1,info)
         CALL CHKXER('CSYTRS_AA_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CSYTRS_AA_2STAGE('U',0,-1,a,1,a,1,ip,ip,b,1,info)
         CALL CHKXER('CSYTRS_AA_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL CSYTRS_AA_2STAGE('U',2,1,a,1,a,1,ip,ip,b,1,info)
         CALL CHKXER('CSYTRS_AA_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL CSYTRS_AA_2STAGE('U',2,1,a,2,a,1,ip,ip,b,1,info)
         CALL CHKXER('CSYTRS_AA_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL CSYTRS_AA_2STAGE('U',2,1,a,2,a,8,ip,ip,b,1,info)
         CALL CHKXER('CSYTRS_AA_STAGE',INFot,NOUt,LERr,OK)
!
      ENDIF
!
!     Print a summary line.
!
      CALL ALAESM(Path,OK,NOUt)
!
!
!     End of CERRSY
!
      END SUBROUTINE CERRSY
