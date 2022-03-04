!*==zerred.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b ZERRED
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZERRED( PATH, NUNIT )
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
!> ZERRED tests the error exits for the eigenvalue driver routines for
!> DOUBLE COMPLEX PRECISION matrices:
!>
!> PATH  driver   description
!> ----  ------   -----------
!> ZEV   ZGEEV    find eigenvalues/eigenvectors for nonsymmetric A
!> ZES   ZGEES    find eigenvalues/Schur form for nonsymmetric A
!> ZVX   ZGEEVX   ZGEEV + balancing and condition estimation
!> ZSX   ZGEESX   ZGEES + balancing and condition estimation
!> ZBD   ZGESVD   compute SVD of an M-by-N matrix A
!>       ZGESDD   compute SVD of an M-by-N matrix A(by divide and
!>                conquer)
!>       ZGEJSV   compute SVD of an M-by-N matrix A where M >= N
!>       ZGESVDX  compute SVD of an M-by-N matrix A(by bisection
!>                and inverse iteration)
!>       ZGESVDQ  compute SVD of an M-by-N matrix A(with a
!>                QR-Preconditioned )
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
!> \date June 2016
!
!> \ingroup complex16_eig
!
!  =====================================================================
      SUBROUTINE ZERRED(Path,Nunit)
      IMPLICIT NONE
!*--ZERRED74
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     June 2016
!
!     .. Scalar Arguments ..
      CHARACTER*3 Path
      INTEGER Nunit
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      INTEGER NMAX , LW
      PARAMETER (NMAX=4,LW=5*NMAX)
      DOUBLE PRECISION ONE , ZERO
      PARAMETER (ONE=1.0D0,ZERO=0.0D0)
!     ..
!     .. Local Scalars ..
      CHARACTER*2 c2
      INTEGER i , ihi , ilo , info , j , ns , nt , sdim
      DOUBLE PRECISION abnrm
!     ..
!     .. Local Arrays ..
      LOGICAL b(NMAX)
      INTEGER iw(4*NMAX)
      DOUBLE PRECISION r1(NMAX) , r2(NMAX) , rw(LW) , s(NMAX)
      COMPLEX*16 a(NMAX,NMAX) , u(NMAX,NMAX) , vl(NMAX,NMAX) ,          &
     &           vr(NMAX,NMAX) , vt(NMAX,NMAX) , w(10*NMAX) , x(NMAX)
!     ..
!     .. External Subroutines ..
      EXTERNAL CHKXER , ZGEES , ZGEESX , ZGEEV , ZGEEVX , ZGESVJ ,      &
     &         ZGESDD , ZGESVD , ZGESVDX , ZGESVQ
!     ..
!     .. External Functions ..
      LOGICAL LSAMEN , ZSLECT
      EXTERNAL LSAMEN , ZSLECT
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC LEN_TRIM
!     ..
!     .. Arrays in Common ..
      LOGICAL SELval(20)
      DOUBLE PRECISION SELwi(20) , SELwr(20)
!     ..
!     .. Scalars in Common ..
      LOGICAL LERr , OK
      CHARACTER*32 SRNamt
      INTEGER INFot , NOUt , SELdim , SELopt
!     ..
!     .. Common blocks ..
      COMMON /INFOC / INFot , NOUt , OK , LERr
      COMMON /SRNAMC/ SRNamt
      COMMON /SSLCT / SELopt , SELdim , SELval , SELwr , SELwi
!     ..
!     .. Executable Statements ..
!
      NOUt = Nunit
      WRITE (NOUt,FMT=*)
      c2 = Path(2:3)
!
!     Initialize A
!
      DO j = 1 , NMAX
         DO i = 1 , NMAX
            a(i,j) = ZERO
         ENDDO
      ENDDO
      DO i = 1 , NMAX
         a(i,i) = ONE
      ENDDO
      OK = .TRUE.
      nt = 0
!
      IF ( LSAMEN(2,c2,'EV') ) THEN
!
!        Test ZGEEV
!
         SRNamt = 'ZGEEV '
         INFot = 1
         CALL ZGEEV('X','N',0,a,1,x,vl,1,vr,1,w,1,rw,info)
         CALL CHKXER('ZGEEV ',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZGEEV('N','X',0,a,1,x,vl,1,vr,1,w,1,rw,info)
         CALL CHKXER('ZGEEV ',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZGEEV('N','N',-1,a,1,x,vl,1,vr,1,w,1,rw,info)
         CALL CHKXER('ZGEEV ',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL ZGEEV('N','N',2,a,1,x,vl,1,vr,1,w,4,rw,info)
         CALL CHKXER('ZGEEV ',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL ZGEEV('V','N',2,a,2,x,vl,1,vr,1,w,4,rw,info)
         CALL CHKXER('ZGEEV ',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL ZGEEV('N','V',2,a,2,x,vl,1,vr,1,w,4,rw,info)
         CALL CHKXER('ZGEEV ',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL ZGEEV('V','V',1,a,1,x,vl,1,vr,1,w,1,rw,info)
         CALL CHKXER('ZGEEV ',INFot,NOUt,LERr,OK)
         nt = nt + 7
!
      ELSEIF ( LSAMEN(2,c2,'ES') ) THEN
!
!        Test ZGEES
!
         SRNamt = 'ZGEES '
         INFot = 1
         CALL ZGEES('X','N',ZSLECT,0,a,1,sdim,x,vl,1,w,1,rw,b,info)
         CALL CHKXER('ZGEES ',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZGEES('N','X',ZSLECT,0,a,1,sdim,x,vl,1,w,1,rw,b,info)
         CALL CHKXER('ZGEES ',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZGEES('N','S',ZSLECT,-1,a,1,sdim,x,vl,1,w,1,rw,b,info)
         CALL CHKXER('ZGEES ',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL ZGEES('N','S',ZSLECT,2,a,1,sdim,x,vl,1,w,4,rw,b,info)
         CALL CHKXER('ZGEES ',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL ZGEES('V','S',ZSLECT,2,a,2,sdim,x,vl,1,w,4,rw,b,info)
         CALL CHKXER('ZGEES ',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL ZGEES('N','S',ZSLECT,1,a,1,sdim,x,vl,1,w,1,rw,b,info)
         CALL CHKXER('ZGEES ',INFot,NOUt,LERr,OK)
         nt = nt + 6
!
      ELSEIF ( LSAMEN(2,c2,'VX') ) THEN
!
!        Test ZGEEVX
!
         SRNamt = 'ZGEEVX'
         INFot = 1
         CALL ZGEEVX('X','N','N','N',0,a,1,x,vl,1,vr,1,ilo,ihi,s,abnrm, &
     &               r1,r2,w,1,rw,info)
         CALL CHKXER('ZGEEVX',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZGEEVX('N','X','N','N',0,a,1,x,vl,1,vr,1,ilo,ihi,s,abnrm, &
     &               r1,r2,w,1,rw,info)
         CALL CHKXER('ZGEEVX',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZGEEVX('N','N','X','N',0,a,1,x,vl,1,vr,1,ilo,ihi,s,abnrm, &
     &               r1,r2,w,1,rw,info)
         CALL CHKXER('ZGEEVX',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZGEEVX('N','N','N','X',0,a,1,x,vl,1,vr,1,ilo,ihi,s,abnrm, &
     &               r1,r2,w,1,rw,info)
         CALL CHKXER('ZGEEVX',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL ZGEEVX('N','N','N','N',-1,a,1,x,vl,1,vr,1,ilo,ihi,s,abnrm,&
     &               r1,r2,w,1,rw,info)
         CALL CHKXER('ZGEEVX',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL ZGEEVX('N','N','N','N',2,a,1,x,vl,1,vr,1,ilo,ihi,s,abnrm, &
     &               r1,r2,w,4,rw,info)
         CALL CHKXER('ZGEEVX',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL ZGEEVX('N','V','N','N',2,a,2,x,vl,1,vr,1,ilo,ihi,s,abnrm, &
     &               r1,r2,w,4,rw,info)
         CALL CHKXER('ZGEEVX',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL ZGEEVX('N','N','V','N',2,a,2,x,vl,1,vr,1,ilo,ihi,s,abnrm, &
     &               r1,r2,w,4,rw,info)
         CALL CHKXER('ZGEEVX',INFot,NOUt,LERr,OK)
         INFot = 20
         CALL ZGEEVX('N','N','N','N',1,a,1,x,vl,1,vr,1,ilo,ihi,s,abnrm, &
     &               r1,r2,w,1,rw,info)
         CALL CHKXER('ZGEEVX',INFot,NOUt,LERr,OK)
         INFot = 20
         CALL ZGEEVX('N','N','V','V',1,a,1,x,vl,1,vr,1,ilo,ihi,s,abnrm, &
     &               r1,r2,w,2,rw,info)
         CALL CHKXER('ZGEEVX',INFot,NOUt,LERr,OK)
         nt = nt + 10
!
      ELSEIF ( LSAMEN(2,c2,'SX') ) THEN
!
!        Test ZGEESX
!
         SRNamt = 'ZGEESX'
         INFot = 1
         CALL ZGEESX('X','N',ZSLECT,'N',0,a,1,sdim,x,vl,1,r1(1),r2(1),w,&
     &               1,rw,b,info)
         CALL CHKXER('ZGEESX',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZGEESX('N','X',ZSLECT,'N',0,a,1,sdim,x,vl,1,r1(1),r2(1),w,&
     &               1,rw,b,info)
         CALL CHKXER('ZGEESX',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZGEESX('N','N',ZSLECT,'X',0,a,1,sdim,x,vl,1,r1(1),r2(1),w,&
     &               1,rw,b,info)
         CALL CHKXER('ZGEESX',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL ZGEESX('N','N',ZSLECT,'N',-1,a,1,sdim,x,vl,1,r1(1),r2(1), &
     &               w,1,rw,b,info)
         CALL CHKXER('ZGEESX',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL ZGEESX('N','N',ZSLECT,'N',2,a,1,sdim,x,vl,1,r1(1),r2(1),w,&
     &               4,rw,b,info)
         CALL CHKXER('ZGEESX',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL ZGEESX('V','N',ZSLECT,'N',2,a,2,sdim,x,vl,1,r1(1),r2(1),w,&
     &               4,rw,b,info)
         CALL CHKXER('ZGEESX',INFot,NOUt,LERr,OK)
         INFot = 15
         CALL ZGEESX('N','N',ZSLECT,'N',1,a,1,sdim,x,vl,1,r1(1),r2(1),w,&
     &               1,rw,b,info)
         CALL CHKXER('ZGEESX',INFot,NOUt,LERr,OK)
         nt = nt + 7
!
      ELSEIF ( LSAMEN(2,c2,'BD') ) THEN
!
!        Test ZGESVD
!
         SRNamt = 'ZGESVD'
         INFot = 1
         CALL ZGESVD('X','N',0,0,a,1,s,u,1,vt,1,w,1,rw,info)
         CALL CHKXER('ZGESVD',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZGESVD('N','X',0,0,a,1,s,u,1,vt,1,w,1,rw,info)
         CALL CHKXER('ZGESVD',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZGESVD('O','O',0,0,a,1,s,u,1,vt,1,w,1,rw,info)
         CALL CHKXER('ZGESVD',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZGESVD('N','N',-1,0,a,1,s,u,1,vt,1,w,1,rw,info)
         CALL CHKXER('ZGESVD',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZGESVD('N','N',0,-1,a,1,s,u,1,vt,1,w,1,rw,info)
         CALL CHKXER('ZGESVD',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL ZGESVD('N','N',2,1,a,1,s,u,1,vt,1,w,5,rw,info)
         CALL CHKXER('ZGESVD',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL ZGESVD('A','N',2,1,a,2,s,u,1,vt,1,w,5,rw,info)
         CALL CHKXER('ZGESVD',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL ZGESVD('N','A',1,2,a,1,s,u,1,vt,1,w,5,rw,info)
         CALL CHKXER('ZGESVD',INFot,NOUt,LERr,OK)
         nt = nt + 8
         IF ( OK ) THEN
            WRITE (NOUt,FMT=99001) SRNamt(1:LEN_TRIM(SRNamt)) , nt
         ELSE
            WRITE (NOUt,FMT=99002)
         ENDIF
!
!        Test ZGESDD
!
         SRNamt = 'ZGESDD'
         INFot = 1
         CALL ZGESDD('X',0,0,a,1,s,u,1,vt,1,w,1,rw,iw,info)
         CALL CHKXER('ZGESDD',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZGESDD('N',-1,0,a,1,s,u,1,vt,1,w,1,rw,iw,info)
         CALL CHKXER('ZGESDD',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZGESDD('N',0,-1,a,1,s,u,1,vt,1,w,1,rw,iw,info)
         CALL CHKXER('ZGESDD',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL ZGESDD('N',2,1,a,1,s,u,1,vt,1,w,5,rw,iw,info)
         CALL CHKXER('ZGESDD',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL ZGESDD('A',2,1,a,2,s,u,1,vt,1,w,5,rw,iw,info)
         CALL CHKXER('ZGESDD',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL ZGESDD('A',1,2,a,1,s,u,1,vt,1,w,5,rw,iw,info)
         CALL CHKXER('ZGESDD',INFot,NOUt,LERr,OK)
         nt = nt - 2
         IF ( OK ) THEN
            WRITE (NOUt,FMT=99001) SRNamt(1:LEN_TRIM(SRNamt)) , nt
         ELSE
            WRITE (NOUt,FMT=99002)
         ENDIF
!
!        Test ZGEJSV
!
         SRNamt = 'ZGEJSV'
         INFot = 1
         CALL ZGEJSV('X','U','V','R','N','N',0,0,a,1,s,u,1,vt,1,w,1,rw, &
     &               1,iw,info)
         CALL CHKXER('ZGEJSV',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZGEJSV('G','X','V','R','N','N',0,0,a,1,s,u,1,vt,1,w,1,rw, &
     &               1,iw,info)
         CALL CHKXER('ZGEJSV',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZGEJSV('G','U','X','R','N','N',0,0,a,1,s,u,1,vt,1,w,1,rw, &
     &               1,iw,info)
         CALL CHKXER('ZGEJSV',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZGEJSV('G','U','V','X','N','N',0,0,a,1,s,u,1,vt,1,w,1,rw, &
     &               1,iw,info)
         CALL CHKXER('ZGEJSV',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL ZGEJSV('G','U','V','R','X','N',0,0,a,1,s,u,1,vt,1,w,1,rw, &
     &               1,iw,info)
         CALL CHKXER('ZGEJSV',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL ZGEJSV('G','U','V','R','N','X',0,0,a,1,s,u,1,vt,1,w,1,rw, &
     &               1,iw,info)
         CALL CHKXER('ZGEJSV',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL ZGEJSV('G','U','V','R','N','N',-1,0,a,1,s,u,1,vt,1,w,1,rw,&
     &               1,iw,info)
         CALL CHKXER('ZGEJSV',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL ZGEJSV('G','U','V','R','N','N',0,-1,a,1,s,u,1,vt,1,w,1,rw,&
     &               1,iw,info)
         CALL CHKXER('ZGEJSV',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL ZGEJSV('G','U','V','R','N','N',2,1,a,1,s,u,1,vt,1,w,1,rw, &
     &               1,iw,info)
         CALL CHKXER('ZGEJSV',INFot,NOUt,LERr,OK)
         INFot = 13
         CALL ZGEJSV('G','U','V','R','N','N',2,2,a,2,s,u,1,vt,2,w,1,rw, &
     &               1,iw,info)
         CALL CHKXER('ZGEJSV',INFot,NOUt,LERr,OK)
         INFot = 15
         CALL ZGEJSV('G','U','V','R','N','N',2,2,a,2,s,u,2,vt,1,w,1,rw, &
     &               1,iw,info)
         CALL CHKXER('ZGEJSV',INFot,NOUt,LERr,OK)
         nt = 11
         IF ( OK ) THEN
            WRITE (NOUt,FMT=99001) SRNamt(1:LEN_TRIM(SRNamt)) , nt
         ELSE
            WRITE (NOUt,FMT=99002)
         ENDIF
!
!        Test ZGESVDX
!
         SRNamt = 'ZGESVDX'
         INFot = 1
         CALL ZGESVDX('X','N','A',0,0,a,1,ZERO,ZERO,0,0,ns,s,u,1,vt,1,w,&
     &                1,rw,iw,info)
         CALL CHKXER('ZGESVDX',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZGESVDX('N','X','A',0,0,a,1,ZERO,ZERO,0,0,ns,s,u,1,vt,1,w,&
     &                1,rw,iw,info)
         CALL CHKXER('ZGESVDX',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZGESVDX('N','N','X',0,0,a,1,ZERO,ZERO,0,0,ns,s,u,1,vt,1,w,&
     &                1,rw,iw,info)
         CALL CHKXER('ZGESVDX',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZGESVDX('N','N','A',-1,0,a,1,ZERO,ZERO,0,0,ns,s,u,1,vt,1, &
     &                w,1,rw,iw,info)
         CALL CHKXER('ZGESVDX',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL ZGESVDX('N','N','A',0,-1,a,1,ZERO,ZERO,0,0,ns,s,u,1,vt,1, &
     &                w,1,rw,iw,info)
         CALL CHKXER('ZGESVDX',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL ZGESVDX('N','N','A',2,1,a,1,ZERO,ZERO,0,0,ns,s,u,1,vt,1,w,&
     &                1,rw,iw,info)
         CALL CHKXER('ZGESVDX',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL ZGESVDX('N','N','V',2,1,a,2,-ONE,ZERO,0,0,ns,s,u,1,vt,1,w,&
     &                1,rw,iw,info)
         CALL CHKXER('ZGESVDX',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL ZGESVDX('N','N','V',2,1,a,2,ONE,ZERO,0,0,ns,s,u,1,vt,1,w, &
     &                1,rw,iw,info)
         CALL CHKXER('ZGESVDX',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL ZGESVDX('N','N','I',2,2,a,2,ZERO,ZERO,0,1,ns,s,u,1,vt,1,w,&
     &                1,rw,iw,info)
         CALL CHKXER('ZGESVDX',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL ZGESVDX('V','N','I',2,2,a,2,ZERO,ZERO,1,0,ns,s,u,1,vt,1,w,&
     &                1,rw,iw,info)
         CALL CHKXER('ZGESVDX',INFot,NOUt,LERr,OK)
         INFot = 15
         CALL ZGESVDX('V','N','A',2,2,a,2,ZERO,ZERO,0,0,ns,s,u,1,vt,1,w,&
     &                1,rw,iw,info)
         CALL CHKXER('ZGESVDX',INFot,NOUt,LERr,OK)
         INFot = 17
         CALL ZGESVDX('N','V','A',2,2,a,2,ZERO,ZERO,0,0,ns,s,u,1,vt,1,w,&
     &                1,rw,iw,info)
         CALL CHKXER('ZGESVDX',INFot,NOUt,LERr,OK)
         nt = 12
         IF ( OK ) THEN
            WRITE (NOUt,FMT=99001) SRNamt(1:LEN_TRIM(SRNamt)) , nt
         ELSE
            WRITE (NOUt,FMT=99002)
         ENDIF
!
!        Test ZGESVDQ
!
         SRNamt = 'ZGESVDQ'
         INFot = 1
         CALL ZGESVDQ('X','P','T','A','A',0,0,a,1,s,u,0,vt,0,ns,iw,1,w, &
     &                1,rw,1,info)
         CALL CHKXER('ZGESVDQ',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZGESVDQ('A','X','T','A','A',0,0,a,1,s,u,0,vt,0,ns,iw,1,w, &
     &                1,rw,1,info)
         CALL CHKXER('ZGESVDQ',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZGESVDQ('A','P','X','A','A',0,0,a,1,s,u,0,vt,0,ns,iw,1,w, &
     &                1,rw,1,info)
         CALL CHKXER('ZGESVDQ',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZGESVDQ('A','P','T','X','A',0,0,a,1,s,u,0,vt,0,ns,iw,1,w, &
     &                1,rw,1,info)
         CALL CHKXER('ZGESVDQ',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL ZGESVDQ('A','P','T','A','X',0,0,a,1,s,u,0,vt,0,ns,iw,1,w, &
     &                1,rw,1,info)
         CALL CHKXER('ZGESVDQ',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL ZGESVDQ('A','P','T','A','A',-1,0,a,1,s,u,0,vt,0,ns,iw,1,w,&
     &                1,rw,1,info)
         CALL CHKXER('ZGESVDQ',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL ZGESVDQ('A','P','T','A','A',0,1,a,1,s,u,0,vt,0,ns,iw,1,w, &
     &                1,rw,1,info)
         CALL CHKXER('ZGESVDQ',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL ZGESVDQ('A','P','T','A','A',1,1,a,0,s,u,0,vt,0,ns,iw,1,w, &
     &                1,rw,1,info)
         CALL CHKXER('ZGESVDQ',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL ZGESVDQ('A','P','T','A','A',1,1,a,1,s,u,-1,vt,0,ns,iw,1,w,&
     &                1,rw,1,info)
         CALL CHKXER('ZGESVDQ',INFot,NOUt,LERr,OK)
         INFot = 14
         CALL ZGESVDQ('A','P','T','A','A',1,1,a,1,s,u,1,vt,-1,ns,iw,1,w,&
     &                1,rw,1,info)
         CALL CHKXER('ZGESVDQ',INFot,NOUt,LERr,OK)
         INFot = 17
         CALL ZGESVDQ('A','P','T','A','A',1,1,a,1,s,u,1,vt,1,ns,iw,-5,w,&
     &                1,rw,1,info)
         CALL CHKXER('ZGESVDQ',INFot,NOUt,LERr,OK)
         nt = 11
         IF ( OK ) THEN
            WRITE (NOUt,FMT=99001) SRNamt(1:LEN_TRIM(SRNamt)) , nt
         ELSE
            WRITE (NOUt,FMT=99002)
         ENDIF
      ENDIF
!
!     Print a summary line.
!
      IF ( .NOT.LSAMEN(2,c2,'BD') ) THEN
         IF ( OK ) THEN
            WRITE (NOUt,FMT=99001) SRNamt(1:LEN_TRIM(SRNamt)) , nt
         ELSE
            WRITE (NOUt,FMT=99002)
         ENDIF
      ENDIF
!
99001 FORMAT (1X,A,' passed the tests of the error exits (',I3,         &
     &        ' tests done)')
99002 FORMAT (' *** ',A,' failed the tests of the error exits ***')
!
!     End of ZERRED
!
      END SUBROUTINE ZERRED
