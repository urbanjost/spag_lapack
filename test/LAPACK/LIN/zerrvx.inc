!*==zerrvx.f90  processed by SPAG 7.51RB at 17:35 on  4 Mar 2022
!> \brief \b zerrvx
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE zerrvx( PATH, NUNIT )
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
!> zerrvx tests the error exits for the COMPLEX*16 driver routines
!> for solving linear systems of equations.
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
      SUBROUTINE zerrvx(Path,Nunit)
      IMPLICIT NONE
!*--ZERRVX59
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
      CHARACTER eq
      CHARACTER*2 c2
      INTEGER i , info , j
      DOUBLE PRECISION rcond
!     ..
!     .. Local Arrays ..
      INTEGER ip(NMAX)
      DOUBLE PRECISION c(NMAX) , r(NMAX) , r1(NMAX) , r2(NMAX) ,        &
     &                 rf(NMAX) , rw(NMAX)
      COMPLEX*16 a(NMAX,NMAX) , af(NMAX,NMAX) , b(NMAX) , e(NMAX) ,     &
     &           w(2*NMAX) , x(NMAX)
!     ..
!     .. External Functions ..
      LOGICAL LSAMEN
EXTERNAL    LSAMEN
!     ..
!     .. External Subroutines ..
EXTERNAL    CHKXER,   ZGBSV,   ZGBSVX,   ZGESV,   ZGESVX,   ZGTSV
EXTERNAL    ZGTSVX,   ZHESV,   ZHESV_RK,   ZHESV_ROOK,   ZHESVX,   ZHPSV
EXTERNAL    ZHPSVX,   ZPBSV,   ZPBSVX,   ZPOSV,   ZPOSVX,   ZPPSV
EXTERNAL    ZPPSVX,   ZPTSV,   ZPTSVX,   ZSPSV,   ZSPSVX,   ZSYSV
EXTERNAL    ZSYSV_AA,   ZSYSV_RK,   ZSYSV_ROOK,   ZSYSVX
EXTERNAL    ZHESV_AA_2STAGE
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
         c(j) = 0.D0
         r(j) = 0.D0
         ip(j) = j
      ENDDO
      eq = ' '
      OK = .TRUE.
!
      IF ( LSAMEN(2,c2,'GE') ) THEN
!
!        ZGESV
!
         SRNamt = 'ZGESV '
         INFot = 1
         CALL ZGESV(-1,0,a,1,ip,b,1,info)
         CALL CHKXER('ZGESV ',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZGESV(0,-1,a,1,ip,b,1,info)
         CALL CHKXER('ZGESV ',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZGESV(2,1,a,1,ip,b,2,info)
         CALL CHKXER('ZGESV ',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL ZGESV(2,1,a,2,ip,b,1,info)
         CALL CHKXER('ZGESV ',INFot,NOUt,LERr,OK)
!
!        ZGESVX
!
         SRNamt = 'ZGESVX'
         INFot = 1
         CALL ZGESVX('/','N',0,0,a,1,af,1,ip,eq,r,c,b,1,x,1,rcond,r1,r2,&
     &               w,rw,info)
         CALL CHKXER('ZGESVX',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZGESVX('N','/',0,0,a,1,af,1,ip,eq,r,c,b,1,x,1,rcond,r1,r2,&
     &               w,rw,info)
         CALL CHKXER('ZGESVX',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZGESVX('N','N',-1,0,a,1,af,1,ip,eq,r,c,b,1,x,1,rcond,r1,  &
     &               r2,w,rw,info)
         CALL CHKXER('ZGESVX',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZGESVX('N','N',0,-1,a,1,af,1,ip,eq,r,c,b,1,x,1,rcond,r1,  &
     &               r2,w,rw,info)
         CALL CHKXER('ZGESVX',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL ZGESVX('N','N',2,1,a,1,af,2,ip,eq,r,c,b,2,x,2,rcond,r1,r2,&
     &               w,rw,info)
         CALL CHKXER('ZGESVX',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL ZGESVX('N','N',2,1,a,2,af,1,ip,eq,r,c,b,2,x,2,rcond,r1,r2,&
     &               w,rw,info)
         CALL CHKXER('ZGESVX',INFot,NOUt,LERr,OK)
         INFot = 10
         eq = '/'
         CALL ZGESVX('F','N',0,0,a,1,af,1,ip,eq,r,c,b,1,x,1,rcond,r1,r2,&
     &               w,rw,info)
         CALL CHKXER('ZGESVX',INFot,NOUt,LERr,OK)
         INFot = 11
         eq = 'R'
         CALL ZGESVX('F','N',1,0,a,1,af,1,ip,eq,r,c,b,1,x,1,rcond,r1,r2,&
     &               w,rw,info)
         CALL CHKXER('ZGESVX',INFot,NOUt,LERr,OK)
         INFot = 12
         eq = 'C'
         CALL ZGESVX('F','N',1,0,a,1,af,1,ip,eq,r,c,b,1,x,1,rcond,r1,r2,&
     &               w,rw,info)
         CALL CHKXER('ZGESVX',INFot,NOUt,LERr,OK)
         INFot = 14
         CALL ZGESVX('N','N',2,1,a,2,af,2,ip,eq,r,c,b,1,x,2,rcond,r1,r2,&
     &               w,rw,info)
         CALL CHKXER('ZGESVX',INFot,NOUt,LERr,OK)
         INFot = 16
         CALL ZGESVX('N','N',2,1,a,2,af,2,ip,eq,r,c,b,2,x,1,rcond,r1,r2,&
     &               w,rw,info)
         CALL CHKXER('ZGESVX',INFot,NOUt,LERr,OK)
!
      ELSEIF ( LSAMEN(2,c2,'GB') ) THEN
!
!        ZGBSV
!
         SRNamt = 'ZGBSV '
         INFot = 1
         CALL ZGBSV(-1,0,0,0,a,1,ip,b,1,info)
         CALL CHKXER('ZGBSV ',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZGBSV(1,-1,0,0,a,1,ip,b,1,info)
         CALL CHKXER('ZGBSV ',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZGBSV(1,0,-1,0,a,1,ip,b,1,info)
         CALL CHKXER('ZGBSV ',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZGBSV(0,0,0,-1,a,1,ip,b,1,info)
         CALL CHKXER('ZGBSV ',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL ZGBSV(1,1,1,0,a,3,ip,b,1,info)
         CALL CHKXER('ZGBSV ',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL ZGBSV(2,0,0,0,a,1,ip,b,1,info)
         CALL CHKXER('ZGBSV ',INFot,NOUt,LERr,OK)
!
!        ZGBSVX
!
         SRNamt = 'ZGBSVX'
         INFot = 1
         CALL ZGBSVX('/','N',0,0,0,0,a,1,af,1,ip,eq,r,c,b,1,x,1,rcond,  &
     &               r1,r2,w,rw,info)
         CALL CHKXER('ZGBSVX',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZGBSVX('N','/',0,0,0,0,a,1,af,1,ip,eq,r,c,b,1,x,1,rcond,  &
     &               r1,r2,w,rw,info)
         CALL CHKXER('ZGBSVX',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZGBSVX('N','N',-1,0,0,0,a,1,af,1,ip,eq,r,c,b,1,x,1,rcond, &
     &               r1,r2,w,rw,info)
         CALL CHKXER('ZGBSVX',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZGBSVX('N','N',1,-1,0,0,a,1,af,1,ip,eq,r,c,b,1,x,1,rcond, &
     &               r1,r2,w,rw,info)
         CALL CHKXER('ZGBSVX',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL ZGBSVX('N','N',1,0,-1,0,a,1,af,1,ip,eq,r,c,b,1,x,1,rcond, &
     &               r1,r2,w,rw,info)
         CALL CHKXER('ZGBSVX',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL ZGBSVX('N','N',0,0,0,-1,a,1,af,1,ip,eq,r,c,b,1,x,1,rcond, &
     &               r1,r2,w,rw,info)
         CALL CHKXER('ZGBSVX',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL ZGBSVX('N','N',1,1,1,0,a,2,af,4,ip,eq,r,c,b,1,x,1,rcond,  &
     &               r1,r2,w,rw,info)
         CALL CHKXER('ZGBSVX',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL ZGBSVX('N','N',1,1,1,0,a,3,af,3,ip,eq,r,c,b,1,x,1,rcond,  &
     &               r1,r2,w,rw,info)
         CALL CHKXER('ZGBSVX',INFot,NOUt,LERr,OK)
         INFot = 12
         eq = '/'
         CALL ZGBSVX('F','N',0,0,0,0,a,1,af,1,ip,eq,r,c,b,1,x,1,rcond,  &
     &               r1,r2,w,rw,info)
         CALL CHKXER('ZGBSVX',INFot,NOUt,LERr,OK)
         INFot = 13
         eq = 'R'
         CALL ZGBSVX('F','N',1,0,0,0,a,1,af,1,ip,eq,r,c,b,1,x,1,rcond,  &
     &               r1,r2,w,rw,info)
         CALL CHKXER('ZGBSVX',INFot,NOUt,LERr,OK)
         INFot = 14
         eq = 'C'
         CALL ZGBSVX('F','N',1,0,0,0,a,1,af,1,ip,eq,r,c,b,1,x,1,rcond,  &
     &               r1,r2,w,rw,info)
         CALL CHKXER('ZGBSVX',INFot,NOUt,LERr,OK)
         INFot = 16
         CALL ZGBSVX('N','N',2,0,0,0,a,1,af,1,ip,eq,r,c,b,1,x,2,rcond,  &
     &               r1,r2,w,rw,info)
         CALL CHKXER('ZGBSVX',INFot,NOUt,LERr,OK)
         INFot = 18
         CALL ZGBSVX('N','N',2,0,0,0,a,1,af,1,ip,eq,r,c,b,2,x,1,rcond,  &
     &               r1,r2,w,rw,info)
         CALL CHKXER('ZGBSVX',INFot,NOUt,LERr,OK)
!
      ELSEIF ( LSAMEN(2,c2,'GT') ) THEN
!
!        ZGTSV
!
         SRNamt = 'ZGTSV '
         INFot = 1
         CALL ZGTSV(-1,0,a(1,1),a(1,2),a(1,3),b,1,info)
         CALL CHKXER('ZGTSV ',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZGTSV(0,-1,a(1,1),a(1,2),a(1,3),b,1,info)
         CALL CHKXER('ZGTSV ',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL ZGTSV(2,0,a(1,1),a(1,2),a(1,3),b,1,info)
         CALL CHKXER('ZGTSV ',INFot,NOUt,LERr,OK)
!
!        ZGTSVX
!
         SRNamt = 'ZGTSVX'
         INFot = 1
         CALL ZGTSVX('/','N',0,0,a(1,1),a(1,2),a(1,3),af(1,1),af(1,2),  &
     &               af(1,3),af(1,4),ip,b,1,x,1,rcond,r1,r2,w,rw,info)
         CALL CHKXER('ZGTSVX',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZGTSVX('N','/',0,0,a(1,1),a(1,2),a(1,3),af(1,1),af(1,2),  &
     &               af(1,3),af(1,4),ip,b,1,x,1,rcond,r1,r2,w,rw,info)
         CALL CHKXER('ZGTSVX',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZGTSVX('N','N',-1,0,a(1,1),a(1,2),a(1,3),af(1,1),af(1,2), &
     &               af(1,3),af(1,4),ip,b,1,x,1,rcond,r1,r2,w,rw,info)
         CALL CHKXER('ZGTSVX',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZGTSVX('N','N',0,-1,a(1,1),a(1,2),a(1,3),af(1,1),af(1,2), &
     &               af(1,3),af(1,4),ip,b,1,x,1,rcond,r1,r2,w,rw,info)
         CALL CHKXER('ZGTSVX',INFot,NOUt,LERr,OK)
         INFot = 14
         CALL ZGTSVX('N','N',2,0,a(1,1),a(1,2),a(1,3),af(1,1),af(1,2),  &
     &               af(1,3),af(1,4),ip,b,1,x,2,rcond,r1,r2,w,rw,info)
         CALL CHKXER('ZGTSVX',INFot,NOUt,LERr,OK)
         INFot = 16
         CALL ZGTSVX('N','N',2,0,a(1,1),a(1,2),a(1,3),af(1,1),af(1,2),  &
     &               af(1,3),af(1,4),ip,b,2,x,1,rcond,r1,r2,w,rw,info)
         CALL CHKXER('ZGTSVX',INFot,NOUt,LERr,OK)
!
      ELSEIF ( LSAMEN(2,c2,'PO') ) THEN
!
!        ZPOSV
!
         SRNamt = 'ZPOSV '
         INFot = 1
         CALL ZPOSV('/',0,0,a,1,b,1,info)
         CALL CHKXER('ZPOSV ',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZPOSV('U',-1,0,a,1,b,1,info)
         CALL CHKXER('ZPOSV ',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZPOSV('U',0,-1,a,1,b,1,info)
         CALL CHKXER('ZPOSV ',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL ZPOSV('U',2,0,a,1,b,2,info)
         CALL CHKXER('ZPOSV ',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL ZPOSV('U',2,0,a,2,b,1,info)
         CALL CHKXER('ZPOSV ',INFot,NOUt,LERr,OK)
!
!        ZPOSVX
!
         SRNamt = 'ZPOSVX'
         INFot = 1
         CALL ZPOSVX('/','U',0,0,a,1,af,1,eq,c,b,1,x,1,rcond,r1,r2,w,rw,&
     &               info)
         CALL CHKXER('ZPOSVX',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZPOSVX('N','/',0,0,a,1,af,1,eq,c,b,1,x,1,rcond,r1,r2,w,rw,&
     &               info)
         CALL CHKXER('ZPOSVX',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZPOSVX('N','U',-1,0,a,1,af,1,eq,c,b,1,x,1,rcond,r1,r2,w,  &
     &               rw,info)
         CALL CHKXER('ZPOSVX',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZPOSVX('N','U',0,-1,a,1,af,1,eq,c,b,1,x,1,rcond,r1,r2,w,  &
     &               rw,info)
         CALL CHKXER('ZPOSVX',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL ZPOSVX('N','U',2,0,a,1,af,2,eq,c,b,2,x,2,rcond,r1,r2,w,rw,&
     &               info)
         CALL CHKXER('ZPOSVX',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL ZPOSVX('N','U',2,0,a,2,af,1,eq,c,b,2,x,2,rcond,r1,r2,w,rw,&
     &               info)
         CALL CHKXER('ZPOSVX',INFot,NOUt,LERr,OK)
         INFot = 9
         eq = '/'
         CALL ZPOSVX('F','U',0,0,a,1,af,1,eq,c,b,1,x,1,rcond,r1,r2,w,rw,&
     &               info)
         CALL CHKXER('ZPOSVX',INFot,NOUt,LERr,OK)
         INFot = 10
         eq = 'Y'
         CALL ZPOSVX('F','U',1,0,a,1,af,1,eq,c,b,1,x,1,rcond,r1,r2,w,rw,&
     &               info)
         CALL CHKXER('ZPOSVX',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL ZPOSVX('N','U',2,0,a,2,af,2,eq,c,b,1,x,2,rcond,r1,r2,w,rw,&
     &               info)
         CALL CHKXER('ZPOSVX',INFot,NOUt,LERr,OK)
         INFot = 14
         CALL ZPOSVX('N','U',2,0,a,2,af,2,eq,c,b,2,x,1,rcond,r1,r2,w,rw,&
     &               info)
         CALL CHKXER('ZPOSVX',INFot,NOUt,LERr,OK)
!
      ELSEIF ( LSAMEN(2,c2,'PP') ) THEN
!
!        ZPPSV
!
         SRNamt = 'ZPPSV '
         INFot = 1
         CALL ZPPSV('/',0,0,a,b,1,info)
         CALL CHKXER('ZPPSV ',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZPPSV('U',-1,0,a,b,1,info)
         CALL CHKXER('ZPPSV ',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZPPSV('U',0,-1,a,b,1,info)
         CALL CHKXER('ZPPSV ',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL ZPPSV('U',2,0,a,b,1,info)
         CALL CHKXER('ZPPSV ',INFot,NOUt,LERr,OK)
!
!        ZPPSVX
!
         SRNamt = 'ZPPSVX'
         INFot = 1
         CALL ZPPSVX('/','U',0,0,a,af,eq,c,b,1,x,1,rcond,r1,r2,w,rw,    &
     &               info)
         CALL CHKXER('ZPPSVX',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZPPSVX('N','/',0,0,a,af,eq,c,b,1,x,1,rcond,r1,r2,w,rw,    &
     &               info)
         CALL CHKXER('ZPPSVX',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZPPSVX('N','U',-1,0,a,af,eq,c,b,1,x,1,rcond,r1,r2,w,rw,   &
     &               info)
         CALL CHKXER('ZPPSVX',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZPPSVX('N','U',0,-1,a,af,eq,c,b,1,x,1,rcond,r1,r2,w,rw,   &
     &               info)
         CALL CHKXER('ZPPSVX',INFot,NOUt,LERr,OK)
         INFot = 7
         eq = '/'
         CALL ZPPSVX('F','U',0,0,a,af,eq,c,b,1,x,1,rcond,r1,r2,w,rw,    &
     &               info)
         CALL CHKXER('ZPPSVX',INFot,NOUt,LERr,OK)
         INFot = 8
         eq = 'Y'
         CALL ZPPSVX('F','U',1,0,a,af,eq,c,b,1,x,1,rcond,r1,r2,w,rw,    &
     &               info)
         CALL CHKXER('ZPPSVX',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL ZPPSVX('N','U',2,0,a,af,eq,c,b,1,x,2,rcond,r1,r2,w,rw,    &
     &               info)
         CALL CHKXER('ZPPSVX',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL ZPPSVX('N','U',2,0,a,af,eq,c,b,2,x,1,rcond,r1,r2,w,rw,    &
     &               info)
         CALL CHKXER('ZPPSVX',INFot,NOUt,LERr,OK)
!
      ELSEIF ( LSAMEN(2,c2,'PB') ) THEN
!
!        ZPBSV
!
         SRNamt = 'ZPBSV '
         INFot = 1
         CALL ZPBSV('/',0,0,0,a,1,b,1,info)
         CALL CHKXER('ZPBSV ',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZPBSV('U',-1,0,0,a,1,b,1,info)
         CALL CHKXER('ZPBSV ',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZPBSV('U',1,-1,0,a,1,b,1,info)
         CALL CHKXER('ZPBSV ',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZPBSV('U',0,0,-1,a,1,b,1,info)
         CALL CHKXER('ZPBSV ',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL ZPBSV('U',1,1,0,a,1,b,2,info)
         CALL CHKXER('ZPBSV ',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL ZPBSV('U',2,0,0,a,1,b,1,info)
         CALL CHKXER('ZPBSV ',INFot,NOUt,LERr,OK)
!
!        ZPBSVX
!
         SRNamt = 'ZPBSVX'
         INFot = 1
         CALL ZPBSVX('/','U',0,0,0,a,1,af,1,eq,c,b,1,x,1,rcond,r1,r2,w, &
     &               rw,info)
         CALL CHKXER('ZPBSVX',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZPBSVX('N','/',0,0,0,a,1,af,1,eq,c,b,1,x,1,rcond,r1,r2,w, &
     &               rw,info)
         CALL CHKXER('ZPBSVX',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZPBSVX('N','U',-1,0,0,a,1,af,1,eq,c,b,1,x,1,rcond,r1,r2,w,&
     &               rw,info)
         CALL CHKXER('ZPBSVX',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZPBSVX('N','U',1,-1,0,a,1,af,1,eq,c,b,1,x,1,rcond,r1,r2,w,&
     &               rw,info)
         CALL CHKXER('ZPBSVX',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL ZPBSVX('N','U',0,0,-1,a,1,af,1,eq,c,b,1,x,1,rcond,r1,r2,w,&
     &               rw,info)
         CALL CHKXER('ZPBSVX',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL ZPBSVX('N','U',1,1,0,a,1,af,2,eq,c,b,2,x,2,rcond,r1,r2,w, &
     &               rw,info)
         CALL CHKXER('ZPBSVX',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL ZPBSVX('N','U',1,1,0,a,2,af,1,eq,c,b,2,x,2,rcond,r1,r2,w, &
     &               rw,info)
         CALL CHKXER('ZPBSVX',INFot,NOUt,LERr,OK)
         INFot = 10
         eq = '/'
         CALL ZPBSVX('F','U',0,0,0,a,1,af,1,eq,c,b,1,x,1,rcond,r1,r2,w, &
     &               rw,info)
         CALL CHKXER('ZPBSVX',INFot,NOUt,LERr,OK)
         INFot = 11
         eq = 'Y'
         CALL ZPBSVX('F','U',1,0,0,a,1,af,1,eq,c,b,1,x,1,rcond,r1,r2,w, &
     &               rw,info)
         CALL CHKXER('ZPBSVX',INFot,NOUt,LERr,OK)
         INFot = 13
         CALL ZPBSVX('N','U',2,0,0,a,1,af,1,eq,c,b,1,x,2,rcond,r1,r2,w, &
     &               rw,info)
         CALL CHKXER('ZPBSVX',INFot,NOUt,LERr,OK)
         INFot = 15
         CALL ZPBSVX('N','U',2,0,0,a,1,af,1,eq,c,b,2,x,1,rcond,r1,r2,w, &
     &               rw,info)
         CALL CHKXER('ZPBSVX',INFot,NOUt,LERr,OK)
!
      ELSEIF ( LSAMEN(2,c2,'PT') ) THEN
!
!        ZPTSV
!
         SRNamt = 'ZPTSV '
         INFot = 1
         CALL ZPTSV(-1,0,r,a(1,1),b,1,info)
         CALL CHKXER('ZPTSV ',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZPTSV(0,-1,r,a(1,1),b,1,info)
         CALL CHKXER('ZPTSV ',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL ZPTSV(2,0,r,a(1,1),b,1,info)
         CALL CHKXER('ZPTSV ',INFot,NOUt,LERr,OK)
!
!        ZPTSVX
!
         SRNamt = 'ZPTSVX'
         INFot = 1
         CALL ZPTSVX('/',0,0,r,a(1,1),rf,af(1,1),b,1,x,1,rcond,r1,r2,w, &
     &               rw,info)
         CALL CHKXER('ZPTSVX',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZPTSVX('N',-1,0,r,a(1,1),rf,af(1,1),b,1,x,1,rcond,r1,r2,w,&
     &               rw,info)
         CALL CHKXER('ZPTSVX',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZPTSVX('N',0,-1,r,a(1,1),rf,af(1,1),b,1,x,1,rcond,r1,r2,w,&
     &               rw,info)
         CALL CHKXER('ZPTSVX',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL ZPTSVX('N',2,0,r,a(1,1),rf,af(1,1),b,1,x,2,rcond,r1,r2,w, &
     &               rw,info)
         CALL CHKXER('ZPTSVX',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL ZPTSVX('N',2,0,r,a(1,1),rf,af(1,1),b,2,x,1,rcond,r1,r2,w, &
     &               rw,info)
         CALL CHKXER('ZPTSVX',INFot,NOUt,LERr,OK)
!
      ELSEIF ( LSAMEN(2,c2,'HE') ) THEN
!
!        ZHESV
!
         SRNamt = 'ZHESV '
         INFot = 1
         CALL ZHESV('/',0,0,a,1,ip,b,1,w,1,info)
         CALL CHKXER('ZHESV ',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZHESV('U',-1,0,a,1,ip,b,1,w,1,info)
         CALL CHKXER('ZHESV ',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZHESV('U',0,-1,a,1,ip,b,1,w,1,info)
         CALL CHKXER('ZHESV ',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL ZHESV('U',2,0,a,1,ip,b,2,w,1,info)
         CALL CHKXER('ZHESV ',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL ZHESV('U',2,0,a,2,ip,b,1,w,1,info)
         CALL CHKXER('ZHESV ',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL ZHESV('U',0,0,a,1,ip,b,1,w,0,info)
         CALL CHKXER('ZHESV ',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL ZHESV('U',0,0,a,1,ip,b,1,w,-2,info)
         CALL CHKXER('ZHESV ',INFot,NOUt,LERr,OK)
!
!        ZHESVX
!
         SRNamt = 'ZHESVX'
         INFot = 1
         CALL ZHESVX('/','U',0,0,a,1,af,1,ip,b,1,x,1,rcond,r1,r2,w,1,rw,&
     &               info)
         CALL CHKXER('ZHESVX',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZHESVX('N','/',0,0,a,1,af,1,ip,b,1,x,1,rcond,r1,r2,w,1,rw,&
     &               info)
         CALL CHKXER('ZHESVX',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZHESVX('N','U',-1,0,a,1,af,1,ip,b,1,x,1,rcond,r1,r2,w,1,  &
     &               rw,info)
         CALL CHKXER('ZHESVX',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZHESVX('N','U',0,-1,a,1,af,1,ip,b,1,x,1,rcond,r1,r2,w,1,  &
     &               rw,info)
         CALL CHKXER('ZHESVX',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL ZHESVX('N','U',2,0,a,1,af,2,ip,b,2,x,2,rcond,r1,r2,w,4,rw,&
     &               info)
         CALL CHKXER('ZHESVX',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL ZHESVX('N','U',2,0,a,2,af,1,ip,b,2,x,2,rcond,r1,r2,w,4,rw,&
     &               info)
         CALL CHKXER('ZHESVX',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL ZHESVX('N','U',2,0,a,2,af,2,ip,b,1,x,2,rcond,r1,r2,w,4,rw,&
     &               info)
         CALL CHKXER('ZHESVX',INFot,NOUt,LERr,OK)
         INFot = 13
         CALL ZHESVX('N','U',2,0,a,2,af,2,ip,b,2,x,1,rcond,r1,r2,w,4,rw,&
     &               info)
         CALL CHKXER('ZHESVX',INFot,NOUt,LERr,OK)
         INFot = 18
         CALL ZHESVX('N','U',2,0,a,2,af,2,ip,b,2,x,2,rcond,r1,r2,w,3,rw,&
     &               info)
         CALL CHKXER('ZHESVX',INFot,NOUt,LERr,OK)
!
      ELSEIF ( LSAMEN(2,c2,'HR') ) THEN
!
!        ZHESV_ROOK
!
         SRNamt = 'ZHESV_ROOK'
         INFot = 1
         CALL ZHESV_ROOK('/',0,0,a,1,ip,b,1,w,1,info)
         CALL CHKXER('ZHESV_ROOK',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZHESV_ROOK('U',-1,0,a,1,ip,b,1,w,1,info)
         CALL CHKXER('ZHESV_ROOK',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZHESV_ROOK('U',0,-1,a,1,ip,b,1,w,1,info)
         CALL CHKXER('ZHESV_ROOK',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL ZHESV_ROOK('U',2,0,a,2,ip,b,1,w,1,info)
         CALL CHKXER('ZHESV_ROOK',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL ZHESV_ROOK('U',0,0,a,1,ip,b,1,w,0,info)
         CALL CHKXER('ZHESV_ROOK',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL ZHESV_ROOK('U',0,0,a,1,ip,b,1,w,-2,info)
         CALL CHKXER('ZHESV_ROOK',INFot,NOUt,LERr,OK)
!
      ELSEIF ( LSAMEN(2,c2,'HK') ) THEN
!
!        ZSYSV_RK
!
!        Test error exits of the driver that uses factorization
!        of a Hermitian indefinite matrix with rook
!        (bounded Bunch-Kaufman) pivoting with the new storage
!        format for factors L ( or U) and D.
!
!        L (or U) is stored in A, diagonal of D is stored on the
!        diagonal of A, subdiagonal of D is stored in a separate array E.
!
         SRNamt = 'ZHESV_RK'
         INFot = 1
         CALL ZHESV_RK('/',0,0,a,1,e,ip,b,1,w,1,info)
         CALL CHKXER('ZHESV_RK',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZHESV_RK('U',-1,0,a,1,e,ip,b,1,w,1,info)
         CALL CHKXER('ZHESV_RK',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZHESV_RK('U',0,-1,a,1,e,ip,b,1,w,1,info)
         CALL CHKXER('ZHESV_RK',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL ZHESV_RK('U',2,0,a,1,e,ip,b,2,w,1,info)
         CALL CHKXER('ZHESV_RK',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL ZHESV_RK('U',2,0,a,2,e,ip,b,1,w,1,info)
         CALL CHKXER('ZHESV_RK',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL ZHESV_RK('U',0,0,a,1,e,ip,b,1,w,0,info)
         CALL CHKXER('ZHESV_RK',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL ZHESV_RK('U',0,0,a,1,e,ip,b,1,w,-2,info)
         CALL CHKXER('ZHESV_RK',INFot,NOUt,LERr,OK)
!
      ELSEIF ( LSAMEN(2,c2,'HA') ) THEN
!
!        ZHESV_AA
!
         SRNamt = 'ZHESV_AA'
         INFot = 1
         CALL ZHESV_AA('/',0,0,a,1,ip,b,1,w,1,info)
         CALL CHKXER('ZHESV_AA',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZHESV_AA('U',-1,0,a,1,ip,b,1,w,1,info)
         CALL CHKXER('ZHESV_AA',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZHESV_AA('U',0,-1,a,1,ip,b,1,w,1,info)
         CALL CHKXER('ZHESV_AA',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL ZHESV_AA('U',2,0,a,2,ip,b,1,w,1,info)
         CALL CHKXER('ZHESV_AA',INFot,NOUt,LERr,OK)
!
      ELSEIF ( LSAMEN(2,c2,'H2') ) THEN
!
!        ZHESV_AASEN_2STAGE
!
         SRNamt = 'ZHESV_AA_2STAGE'
         INFot = 1
         CALL ZHESV_AA_2STAGE('/',0,0,a,1,a,1,ip,ip,b,1,w,1,info)
         CALL CHKXER('ZHESV_AA_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZHESV_AA_2STAGE('U',-1,0,a,1,a,1,ip,ip,b,1,w,1,info)
         CALL CHKXER('ZHESV_AA_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZHESV_AA_2STAGE('U',0,-1,a,1,a,1,ip,ip,b,1,w,1,info)
         CALL CHKXER('ZHESV_AA_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL ZHESV_AA_2STAGE('U',2,1,a,1,a,1,ip,ip,b,1,w,1,info)
         CALL CHKXER('ZHESV_AA_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL ZHESV_AA_2STAGE('U',2,1,a,2,a,8,ip,ip,b,1,w,1,info)
         CALL CHKXER('ZHESV_AA_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL ZHESV_AA_2STAGE('U',2,1,a,2,a,1,ip,ip,b,2,w,1,info)
         CALL CHKXER('ZHESV_AA_2STAGE',INFot,NOUt,LERr,OK)
!
      ELSEIF ( LSAMEN(2,c2,'S2') ) THEN
!
!        ZSYSV_AASEN_2STAGE
!
         SRNamt = 'ZSYSV_AA_2STAGE'
         INFot = 1
         CALL ZSYSV_AA_2STAGE('/',0,0,a,1,a,1,ip,ip,b,1,w,1,info)
         CALL CHKXER('ZSYSV_AA_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZSYSV_AA_2STAGE('U',-1,0,a,1,a,1,ip,ip,b,1,w,1,info)
         CALL CHKXER('ZSYSV_AA_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZSYSV_AA_2STAGE('U',0,-1,a,1,a,1,ip,ip,b,1,w,1,info)
         CALL CHKXER('ZSYSV_AA_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL ZSYSV_AA_2STAGE('U',2,1,a,1,a,1,ip,ip,b,1,w,1,info)
         CALL CHKXER('ZSYSV_AA_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL ZSYSV_AA_2STAGE('U',2,1,a,2,a,8,ip,ip,b,1,w,1,info)
         CALL CHKXER('ZSYSV_AA_2STAGE',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL ZSYSV_AA_2STAGE('U',2,1,a,2,a,1,ip,ip,b,2,w,1,info)
         CALL CHKXER('ZSYSV_AA_2STAGE',INFot,NOUt,LERr,OK)
!*
      ELSEIF ( LSAMEN(2,c2,'HP') ) THEN
!
!        ZHPSV
!
         SRNamt = 'ZHPSV '
         INFot = 1
         CALL ZHPSV('/',0,0,a,ip,b,1,info)
         CALL CHKXER('ZHPSV ',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZHPSV('U',-1,0,a,ip,b,1,info)
         CALL CHKXER('ZHPSV ',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZHPSV('U',0,-1,a,ip,b,1,info)
         CALL CHKXER('ZHPSV ',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL ZHPSV('U',2,0,a,ip,b,1,info)
         CALL CHKXER('ZHPSV ',INFot,NOUt,LERr,OK)
!
!        ZHPSVX
!
         SRNamt = 'ZHPSVX'
         INFot = 1
         CALL ZHPSVX('/','U',0,0,a,af,ip,b,1,x,1,rcond,r1,r2,w,rw,info)
         CALL CHKXER('ZHPSVX',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZHPSVX('N','/',0,0,a,af,ip,b,1,x,1,rcond,r1,r2,w,rw,info)
         CALL CHKXER('ZHPSVX',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZHPSVX('N','U',-1,0,a,af,ip,b,1,x,1,rcond,r1,r2,w,rw,info)
         CALL CHKXER('ZHPSVX',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZHPSVX('N','U',0,-1,a,af,ip,b,1,x,1,rcond,r1,r2,w,rw,info)
         CALL CHKXER('ZHPSVX',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL ZHPSVX('N','U',2,0,a,af,ip,b,1,x,2,rcond,r1,r2,w,rw,info)
         CALL CHKXER('ZHPSVX',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL ZHPSVX('N','U',2,0,a,af,ip,b,2,x,1,rcond,r1,r2,w,rw,info)
         CALL CHKXER('ZHPSVX',INFot,NOUt,LERr,OK)
!
      ELSEIF ( LSAMEN(2,c2,'SY') ) THEN
!
!        ZSYSV
!
         SRNamt = 'ZSYSV '
         INFot = 1
         CALL ZSYSV('/',0,0,a,1,ip,b,1,w,1,info)
         CALL CHKXER('ZSYSV ',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZSYSV('U',-1,0,a,1,ip,b,1,w,1,info)
         CALL CHKXER('ZSYSV ',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZSYSV('U',0,-1,a,1,ip,b,1,w,1,info)
         CALL CHKXER('ZSYSV ',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL ZSYSV('U',2,0,a,2,ip,b,1,w,1,info)
         CALL CHKXER('ZSYSV ',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL ZSYSV('U',0,0,a,1,ip,b,1,w,0,info)
         CALL CHKXER('ZSYSV ',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL ZSYSV('U',0,0,a,1,ip,b,1,w,-2,info)
         CALL CHKXER('ZSYSV ',INFot,NOUt,LERr,OK)
!
!        ZSYSVX
!
         SRNamt = 'ZSYSVX'
         INFot = 1
         CALL ZSYSVX('/','U',0,0,a,1,af,1,ip,b,1,x,1,rcond,r1,r2,w,1,rw,&
     &               info)
         CALL CHKXER('ZSYSVX',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZSYSVX('N','/',0,0,a,1,af,1,ip,b,1,x,1,rcond,r1,r2,w,1,rw,&
     &               info)
         CALL CHKXER('ZSYSVX',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZSYSVX('N','U',-1,0,a,1,af,1,ip,b,1,x,1,rcond,r1,r2,w,1,  &
     &               rw,info)
         CALL CHKXER('ZSYSVX',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZSYSVX('N','U',0,-1,a,1,af,1,ip,b,1,x,1,rcond,r1,r2,w,1,  &
     &               rw,info)
         CALL CHKXER('ZSYSVX',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL ZSYSVX('N','U',2,0,a,1,af,2,ip,b,2,x,2,rcond,r1,r2,w,4,rw,&
     &               info)
         CALL CHKXER('ZSYSVX',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL ZSYSVX('N','U',2,0,a,2,af,1,ip,b,2,x,2,rcond,r1,r2,w,4,rw,&
     &               info)
         CALL CHKXER('ZSYSVX',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL ZSYSVX('N','U',2,0,a,2,af,2,ip,b,1,x,2,rcond,r1,r2,w,4,rw,&
     &               info)
         CALL CHKXER('ZSYSVX',INFot,NOUt,LERr,OK)
         INFot = 13
         CALL ZSYSVX('N','U',2,0,a,2,af,2,ip,b,2,x,1,rcond,r1,r2,w,4,rw,&
     &               info)
         CALL CHKXER('ZSYSVX',INFot,NOUt,LERr,OK)
         INFot = 18
         CALL ZSYSVX('N','U',2,0,a,2,af,2,ip,b,2,x,2,rcond,r1,r2,w,3,rw,&
     &               info)
         CALL CHKXER('ZSYSVX',INFot,NOUt,LERr,OK)
!
      ELSEIF ( LSAMEN(2,c2,'SR') ) THEN
!
!        ZSYSV_ROOK
!
         SRNamt = 'ZSYSV_ROOK'
         INFot = 1
         CALL ZSYSV_ROOK('/',0,0,a,1,ip,b,1,w,1,info)
         CALL CHKXER('ZSYSV_ROOK',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZSYSV_ROOK('U',-1,0,a,1,ip,b,1,w,1,info)
         CALL CHKXER('ZSYSV_ROOK',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZSYSV_ROOK('U',0,-1,a,1,ip,b,1,w,1,info)
         CALL CHKXER('ZSYSV_ROOK',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL ZSYSV_ROOK('U',2,0,a,2,ip,b,1,w,1,info)
         CALL CHKXER('ZSYSV_ROOK',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL ZSYSV_ROOK('U',0,0,a,1,ip,b,1,w,0,info)
         CALL CHKXER('ZSYSV_ROOK',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL ZSYSV_ROOK('U',0,0,a,1,ip,b,1,w,-2,info)
!
      ELSEIF ( LSAMEN(2,c2,'SK') ) THEN
!
!        ZSYSV_RK
!
!        Test error exits of the driver that uses factorization
!        of a symmetric indefinite matrix with rook
!        (bounded Bunch-Kaufman) pivoting with the new storage
!        format for factors L ( or U) and D.
!
!        L (or U) is stored in A, diagonal of D is stored on the
!        diagonal of A, subdiagonal of D is stored in a separate array E.
!
         SRNamt = 'ZSYSV_RK'
         INFot = 1
         CALL ZSYSV_RK('/',0,0,a,1,e,ip,b,1,w,1,info)
         CALL CHKXER('ZSYSV_RK',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZSYSV_RK('U',-1,0,a,1,e,ip,b,1,w,1,info)
         CALL CHKXER('ZSYSV_RK',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZSYSV_RK('U',0,-1,a,1,e,ip,b,1,w,1,info)
         CALL CHKXER('ZSYSV_RK',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL ZSYSV_RK('U',2,0,a,1,e,ip,b,2,w,1,info)
         CALL CHKXER('ZSYSV_RK',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL ZSYSV_RK('U',2,0,a,2,e,ip,b,1,w,1,info)
         CALL CHKXER('ZSYSV_RK',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL ZSYSV_RK('U',0,0,a,1,e,ip,b,1,w,0,info)
         CALL CHKXER('ZSYSV_RK',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL ZSYSV_RK('U',0,0,a,1,e,ip,b,1,w,-2,info)
         CALL CHKXER('ZSYSV_RK',INFot,NOUt,LERr,OK)
!
      ELSEIF ( LSAMEN(2,c2,'SP') ) THEN
!
!        ZSPSV
!
         SRNamt = 'ZSPSV '
         INFot = 1
         CALL ZSPSV('/',0,0,a,ip,b,1,info)
         CALL CHKXER('ZSPSV ',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZSPSV('U',-1,0,a,ip,b,1,info)
         CALL CHKXER('ZSPSV ',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZSPSV('U',0,-1,a,ip,b,1,info)
         CALL CHKXER('ZSPSV ',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL ZSPSV('U',2,0,a,ip,b,1,info)
         CALL CHKXER('ZSPSV ',INFot,NOUt,LERr,OK)
!
!        ZSPSVX
!
         SRNamt = 'ZSPSVX'
         INFot = 1
         CALL ZSPSVX('/','U',0,0,a,af,ip,b,1,x,1,rcond,r1,r2,w,rw,info)
         CALL CHKXER('ZSPSVX',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZSPSVX('N','/',0,0,a,af,ip,b,1,x,1,rcond,r1,r2,w,rw,info)
         CALL CHKXER('ZSPSVX',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZSPSVX('N','U',-1,0,a,af,ip,b,1,x,1,rcond,r1,r2,w,rw,info)
         CALL CHKXER('ZSPSVX',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZSPSVX('N','U',0,-1,a,af,ip,b,1,x,1,rcond,r1,r2,w,rw,info)
         CALL CHKXER('ZSPSVX',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL ZSPSVX('N','U',2,0,a,af,ip,b,1,x,2,rcond,r1,r2,w,rw,info)
         CALL CHKXER('ZSPSVX',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL ZSPSVX('N','U',2,0,a,af,ip,b,2,x,1,rcond,r1,r2,w,rw,info)
         CALL CHKXER('ZSPSVX',INFot,NOUt,LERr,OK)
      ENDIF
!
!     Print a summary line.
!
      IF ( OK ) THEN
         WRITE (NOUt,FMT=99001) Path
      ELSE
         WRITE (NOUt,FMT=99002) Path
      ENDIF
!
99001 FORMAT (1X,A3,' drivers passed the tests of the error exits')
99002 FORMAT (' *** ',A3,' drivers failed the tests of the error ',     &
     &        'exits ***')
!
!
!     End of zerrvx
!
      END SUBROUTINE zerrvx
