!*==derred.f90  processed by SPAG 7.51RB at 17:37 on  4 Mar 2022
!> \brief \b DERRED
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DERRED( PATH, NUNIT )
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
!> DERRED tests the error exits for the eigenvalue driver routines for
!> DOUBLE PRECISION matrices:
!>
!> PATH  driver   description
!> ----  ------   -----------
!> SEV   DGEEV    find eigenvalues/eigenvectors for nonsymmetric A
!> SES   DGEES    find eigenvalues/Schur form for nonsymmetric A
!> SVX   DGEEVX   SGEEV + balancing and condition estimation
!> SSX   DGEESX   SGEES + balancing and condition estimation
!> DBD   DGESVD   compute SVD of an M-by-N matrix A
!>       DGESDD   compute SVD of an M-by-N matrix A (by divide and
!>                conquer)
!>       DGEJSV   compute SVD of an M-by-N matrix A where M >= N
!>       DGESVDX  compute SVD of an M-by-N matrix A(by bisection
!>                and inverse iteration)
!>       DGESVDQ  compute SVD of an M-by-N matrix A(with a
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
!> \ingroup double_eig
!
!  =====================================================================
      SUBROUTINE DERRED(Path,Nunit)
      IMPLICIT NONE
!*--DERRED74
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
      INTEGER NMAX
      DOUBLE PRECISION ONE , ZERO
      PARAMETER (NMAX=4,ONE=1.0D0,ZERO=0.0D0)
!     ..
!     .. Local Scalars ..
      CHARACTER*2 c2
      INTEGER i , ihi , ilo , info , j , ns , nt , sdim
      DOUBLE PRECISION abnrm
!     ..
!     .. Local Arrays ..
      LOGICAL b(NMAX)
      INTEGER iw(2*NMAX)
      DOUBLE PRECISION a(NMAX,NMAX) , r1(NMAX) , r2(NMAX) , s(NMAX) ,   &
     &                 u(NMAX,NMAX) , vl(NMAX,NMAX) , vr(NMAX,NMAX) ,   &
     &                 vt(NMAX,NMAX) , w(10*NMAX) , wi(NMAX) , wr(NMAX)
!     ..
!     .. External Subroutines ..
      EXTERNAL CHKXER , DGEES , DGEESX , DGEEV , DGEEVX , DGEJSV ,      &
     &         DGESDD , DGESVD , DGESVDX , DGESVQ
!     ..
!     .. External Functions ..
      LOGICAL DSLECT , LSAMEN
      EXTERNAL DSLECT , LSAMEN
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
!        Test DGEEV
!
         SRNamt = 'DGEEV '
         INFot = 1
         CALL DGEEV('X','N',0,a,1,wr,wi,vl,1,vr,1,w,1,info)
         CALL CHKXER('DGEEV ',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DGEEV('N','X',0,a,1,wr,wi,vl,1,vr,1,w,1,info)
         CALL CHKXER('DGEEV ',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DGEEV('N','N',-1,a,1,wr,wi,vl,1,vr,1,w,1,info)
         CALL CHKXER('DGEEV ',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL DGEEV('N','N',2,a,1,wr,wi,vl,1,vr,1,w,6,info)
         CALL CHKXER('DGEEV ',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL DGEEV('V','N',2,a,2,wr,wi,vl,1,vr,1,w,8,info)
         CALL CHKXER('DGEEV ',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL DGEEV('N','V',2,a,2,wr,wi,vl,1,vr,1,w,8,info)
         CALL CHKXER('DGEEV ',INFot,NOUt,LERr,OK)
         INFot = 13
         CALL DGEEV('V','V',1,a,1,wr,wi,vl,1,vr,1,w,3,info)
         CALL CHKXER('DGEEV ',INFot,NOUt,LERr,OK)
         nt = nt + 7
!
      ELSEIF ( LSAMEN(2,c2,'ES') ) THEN
!
!        Test DGEES
!
         SRNamt = 'DGEES '
         INFot = 1
         CALL DGEES('X','N',DSLECT,0,a,1,sdim,wr,wi,vl,1,w,1,b,info)
         CALL CHKXER('DGEES ',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DGEES('N','X',DSLECT,0,a,1,sdim,wr,wi,vl,1,w,1,b,info)
         CALL CHKXER('DGEES ',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL DGEES('N','S',DSLECT,-1,a,1,sdim,wr,wi,vl,1,w,1,b,info)
         CALL CHKXER('DGEES ',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL DGEES('N','S',DSLECT,2,a,1,sdim,wr,wi,vl,1,w,6,b,info)
         CALL CHKXER('DGEES ',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL DGEES('V','S',DSLECT,2,a,2,sdim,wr,wi,vl,1,w,6,b,info)
         CALL CHKXER('DGEES ',INFot,NOUt,LERr,OK)
         INFot = 13
         CALL DGEES('N','S',DSLECT,1,a,1,sdim,wr,wi,vl,1,w,2,b,info)
         CALL CHKXER('DGEES ',INFot,NOUt,LERr,OK)
         nt = nt + 6
!
      ELSEIF ( LSAMEN(2,c2,'VX') ) THEN
!
!        Test DGEEVX
!
         SRNamt = 'DGEEVX'
         INFot = 1
         CALL DGEEVX('X','N','N','N',0,a,1,wr,wi,vl,1,vr,1,ilo,ihi,s,   &
     &               abnrm,r1,r2,w,1,iw,info)
         CALL CHKXER('DGEEVX',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DGEEVX('N','X','N','N',0,a,1,wr,wi,vl,1,vr,1,ilo,ihi,s,   &
     &               abnrm,r1,r2,w,1,iw,info)
         CALL CHKXER('DGEEVX',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DGEEVX('N','N','X','N',0,a,1,wr,wi,vl,1,vr,1,ilo,ihi,s,   &
     &               abnrm,r1,r2,w,1,iw,info)
         CALL CHKXER('DGEEVX',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL DGEEVX('N','N','N','X',0,a,1,wr,wi,vl,1,vr,1,ilo,ihi,s,   &
     &               abnrm,r1,r2,w,1,iw,info)
         CALL CHKXER('DGEEVX',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL DGEEVX('N','N','N','N',-1,a,1,wr,wi,vl,1,vr,1,ilo,ihi,s,  &
     &               abnrm,r1,r2,w,1,iw,info)
         CALL CHKXER('DGEEVX',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL DGEEVX('N','N','N','N',2,a,1,wr,wi,vl,1,vr,1,ilo,ihi,s,   &
     &               abnrm,r1,r2,w,1,iw,info)
         CALL CHKXER('DGEEVX',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL DGEEVX('N','V','N','N',2,a,2,wr,wi,vl,1,vr,1,ilo,ihi,s,   &
     &               abnrm,r1,r2,w,6,iw,info)
         CALL CHKXER('DGEEVX',INFot,NOUt,LERr,OK)
         INFot = 13
         CALL DGEEVX('N','N','V','N',2,a,2,wr,wi,vl,1,vr,1,ilo,ihi,s,   &
     &               abnrm,r1,r2,w,6,iw,info)
         CALL CHKXER('DGEEVX',INFot,NOUt,LERr,OK)
         INFot = 21
         CALL DGEEVX('N','N','N','N',1,a,1,wr,wi,vl,1,vr,1,ilo,ihi,s,   &
     &               abnrm,r1,r2,w,1,iw,info)
         CALL CHKXER('DGEEVX',INFot,NOUt,LERr,OK)
         INFot = 21
         CALL DGEEVX('N','V','N','N',1,a,1,wr,wi,vl,1,vr,1,ilo,ihi,s,   &
     &               abnrm,r1,r2,w,2,iw,info)
         CALL CHKXER('DGEEVX',INFot,NOUt,LERr,OK)
         INFot = 21
         CALL DGEEVX('N','N','V','V',1,a,1,wr,wi,vl,1,vr,1,ilo,ihi,s,   &
     &               abnrm,r1,r2,w,3,iw,info)
         CALL CHKXER('DGEEVX',INFot,NOUt,LERr,OK)
         nt = nt + 11
!
      ELSEIF ( LSAMEN(2,c2,'SX') ) THEN
!
!        Test DGEESX
!
         SRNamt = 'DGEESX'
         INFot = 1
         CALL DGEESX('X','N',DSLECT,'N',0,a,1,sdim,wr,wi,vl,1,r1(1),    &
     &               r2(1),w,1,iw,1,b,info)
         CALL CHKXER('DGEESX',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DGEESX('N','X',DSLECT,'N',0,a,1,sdim,wr,wi,vl,1,r1(1),    &
     &               r2(1),w,1,iw,1,b,info)
         CALL CHKXER('DGEESX',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL DGEESX('N','N',DSLECT,'X',0,a,1,sdim,wr,wi,vl,1,r1(1),    &
     &               r2(1),w,1,iw,1,b,info)
         CALL CHKXER('DGEESX',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL DGEESX('N','N',DSLECT,'N',-1,a,1,sdim,wr,wi,vl,1,r1(1),   &
     &               r2(1),w,1,iw,1,b,info)
         CALL CHKXER('DGEESX',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL DGEESX('N','N',DSLECT,'N',2,a,1,sdim,wr,wi,vl,1,r1(1),    &
     &               r2(1),w,6,iw,1,b,info)
         CALL CHKXER('DGEESX',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL DGEESX('V','N',DSLECT,'N',2,a,2,sdim,wr,wi,vl,1,r1(1),    &
     &               r2(1),w,6,iw,1,b,info)
         CALL CHKXER('DGEESX',INFot,NOUt,LERr,OK)
         INFot = 16
         CALL DGEESX('N','N',DSLECT,'N',1,a,1,sdim,wr,wi,vl,1,r1(1),    &
     &               r2(1),w,2,iw,1,b,info)
         CALL CHKXER('DGEESX',INFot,NOUt,LERr,OK)
         nt = nt + 7
!
      ELSEIF ( LSAMEN(2,c2,'BD') ) THEN
!
!        Test DGESVD
!
         SRNamt = 'DGESVD'
         INFot = 1
         CALL DGESVD('X','N',0,0,a,1,s,u,1,vt,1,w,1,info)
         CALL CHKXER('DGESVD',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DGESVD('N','X',0,0,a,1,s,u,1,vt,1,w,1,info)
         CALL CHKXER('DGESVD',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DGESVD('O','O',0,0,a,1,s,u,1,vt,1,w,1,info)
         CALL CHKXER('DGESVD',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DGESVD('N','N',-1,0,a,1,s,u,1,vt,1,w,1,info)
         CALL CHKXER('DGESVD',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL DGESVD('N','N',0,-1,a,1,s,u,1,vt,1,w,1,info)
         CALL CHKXER('DGESVD',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL DGESVD('N','N',2,1,a,1,s,u,1,vt,1,w,5,info)
         CALL CHKXER('DGESVD',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL DGESVD('A','N',2,1,a,2,s,u,1,vt,1,w,5,info)
         CALL CHKXER('DGESVD',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL DGESVD('N','A',1,2,a,1,s,u,1,vt,1,w,5,info)
         CALL CHKXER('DGESVD',INFot,NOUt,LERr,OK)
         nt = 8
         IF ( OK ) THEN
            WRITE (NOUt,FMT=99001) SRNamt(1:LEN_TRIM(SRNamt)) , nt
         ELSE
            WRITE (NOUt,FMT=99002)
         ENDIF
!
!        Test DGESDD
!
         SRNamt = 'DGESDD'
         INFot = 1
         CALL DGESDD('X',0,0,a,1,s,u,1,vt,1,w,1,iw,info)
         CALL CHKXER('DGESDD',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DGESDD('N',-1,0,a,1,s,u,1,vt,1,w,1,iw,info)
         CALL CHKXER('DGESDD',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DGESDD('N',0,-1,a,1,s,u,1,vt,1,w,1,iw,info)
         CALL CHKXER('DGESDD',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL DGESDD('N',2,1,a,1,s,u,1,vt,1,w,5,iw,info)
         CALL CHKXER('DGESDD',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL DGESDD('A',2,1,a,2,s,u,1,vt,1,w,5,iw,info)
         CALL CHKXER('DGESDD',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL DGESDD('A',1,2,a,1,s,u,1,vt,1,w,5,iw,info)
         CALL CHKXER('DGESDD',INFot,NOUt,LERr,OK)
         nt = 6
         IF ( OK ) THEN
            WRITE (NOUt,FMT=99001) SRNamt(1:LEN_TRIM(SRNamt)) , nt
         ELSE
            WRITE (NOUt,FMT=99002)
         ENDIF
!
!        Test DGEJSV
!
         SRNamt = 'DGEJSV'
         INFot = 1
         CALL DGEJSV('X','U','V','R','N','N',0,0,a,1,s,u,1,vt,1,w,1,iw, &
     &               info)
         CALL CHKXER('DGEJSV',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DGEJSV('G','X','V','R','N','N',0,0,a,1,s,u,1,vt,1,w,1,iw, &
     &               info)
         CALL CHKXER('DGEJSV',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DGEJSV('G','U','X','R','N','N',0,0,a,1,s,u,1,vt,1,w,1,iw, &
     &               info)
         CALL CHKXER('DGEJSV',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL DGEJSV('G','U','V','X','N','N',0,0,a,1,s,u,1,vt,1,w,1,iw, &
     &               info)
         CALL CHKXER('DGEJSV',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL DGEJSV('G','U','V','R','X','N',0,0,a,1,s,u,1,vt,1,w,1,iw, &
     &               info)
         CALL CHKXER('DGEJSV',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL DGEJSV('G','U','V','R','N','X',0,0,a,1,s,u,1,vt,1,w,1,iw, &
     &               info)
         CALL CHKXER('DGEJSV',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL DGEJSV('G','U','V','R','N','N',-1,0,a,1,s,u,1,vt,1,w,1,iw,&
     &               info)
         CALL CHKXER('DGEJSV',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL DGEJSV('G','U','V','R','N','N',0,-1,a,1,s,u,1,vt,1,w,1,iw,&
     &               info)
         CALL CHKXER('DGEJSV',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL DGEJSV('G','U','V','R','N','N',2,1,a,1,s,u,1,vt,1,w,1,iw, &
     &               info)
         CALL CHKXER('DGEJSV',INFot,NOUt,LERr,OK)
         INFot = 13
         CALL DGEJSV('G','U','V','R','N','N',2,2,a,2,s,u,1,vt,2,w,1,iw, &
     &               info)
         CALL CHKXER('DGEJSV',INFot,NOUt,LERr,OK)
         INFot = 15
         CALL DGEJSV('G','U','V','R','N','N',2,2,a,2,s,u,2,vt,1,w,1,iw, &
     &               info)
         CALL CHKXER('DGEJSV',INFot,NOUt,LERr,OK)
         nt = 11
         IF ( OK ) THEN
            WRITE (NOUt,FMT=99001) SRNamt(1:LEN_TRIM(SRNamt)) , nt
         ELSE
            WRITE (NOUt,FMT=99002)
         ENDIF
!
!        Test DGESVDX
!
         SRNamt = 'DGESVDX'
         INFot = 1
         CALL DGESVDX('X','N','A',0,0,a,1,ZERO,ZERO,0,0,ns,s,u,1,vt,1,w,&
     &                1,iw,info)
         CALL CHKXER('DGESVDX',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DGESVDX('N','X','A',0,0,a,1,ZERO,ZERO,0,0,ns,s,u,1,vt,1,w,&
     &                1,iw,info)
         CALL CHKXER('DGESVDX',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DGESVDX('N','N','X',0,0,a,1,ZERO,ZERO,0,0,ns,s,u,1,vt,1,w,&
     &                1,iw,info)
         CALL CHKXER('DGESVDX',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL DGESVDX('N','N','A',-1,0,a,1,ZERO,ZERO,0,0,ns,s,u,1,vt,1, &
     &                w,1,iw,info)
         CALL CHKXER('DGESVDX',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL DGESVDX('N','N','A',0,-1,a,1,ZERO,ZERO,0,0,ns,s,u,1,vt,1, &
     &                w,1,iw,info)
         CALL CHKXER('DGESVDX',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL DGESVDX('N','N','A',2,1,a,1,ZERO,ZERO,0,0,ns,s,u,1,vt,1,w,&
     &                1,iw,info)
         CALL CHKXER('DGESVDX',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL DGESVDX('N','N','V',2,1,a,2,-ONE,ZERO,0,0,ns,s,u,1,vt,1,w,&
     &                1,iw,info)
         CALL CHKXER('DGESVDX',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL DGESVDX('N','N','V',2,1,a,2,ONE,ZERO,0,0,ns,s,u,1,vt,1,w, &
     &                1,iw,info)
         CALL CHKXER('DGESVDX',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL DGESVDX('N','N','I',2,2,a,2,ZERO,ZERO,0,1,ns,s,u,1,vt,1,w,&
     &                1,iw,info)
         CALL CHKXER('DGESVDX',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL DGESVDX('V','N','I',2,2,a,2,ZERO,ZERO,1,0,ns,s,u,1,vt,1,w,&
     &                1,iw,info)
         CALL CHKXER('DGESVDX',INFot,NOUt,LERr,OK)
         INFot = 15
         CALL DGESVDX('V','N','A',2,2,a,2,ZERO,ZERO,0,0,ns,s,u,1,vt,1,w,&
     &                1,iw,info)
         CALL CHKXER('DGESVDX',INFot,NOUt,LERr,OK)
         INFot = 17
         CALL DGESVDX('N','V','A',2,2,a,2,ZERO,ZERO,0,0,ns,s,u,1,vt,1,w,&
     &                1,iw,info)
         CALL CHKXER('DGESVDX',INFot,NOUt,LERr,OK)
         nt = 12
         IF ( OK ) THEN
            WRITE (NOUt,FMT=99001) SRNamt(1:LEN_TRIM(SRNamt)) , nt
         ELSE
            WRITE (NOUt,FMT=99002)
         ENDIF
!
!        Test DGESVDQ
!
         SRNamt = 'DGESVDQ'
         INFot = 1
         CALL DGESVDQ('X','P','T','A','A',0,0,a,1,s,u,0,vt,0,ns,iw,1,w, &
     &                1,w,1,info)
         CALL CHKXER('DGESVDQ',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DGESVDQ('A','X','T','A','A',0,0,a,1,s,u,0,vt,0,ns,iw,1,w, &
     &                1,w,1,info)
         CALL CHKXER('DGESVDQ',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DGESVDQ('A','P','X','A','A',0,0,a,1,s,u,0,vt,0,ns,iw,1,w, &
     &                1,w,1,info)
         CALL CHKXER('DGESVDQ',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL DGESVDQ('A','P','T','X','A',0,0,a,1,s,u,0,vt,0,ns,iw,1,w, &
     &                1,w,1,info)
         CALL CHKXER('DGESVDQ',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL DGESVDQ('A','P','T','A','X',0,0,a,1,s,u,0,vt,0,ns,iw,1,w, &
     &                1,w,1,info)
         CALL CHKXER('DGESVDQ',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL DGESVDQ('A','P','T','A','A',-1,0,a,1,s,u,0,vt,0,ns,iw,1,w,&
     &                1,w,1,info)
         CALL CHKXER('DGESVDQ',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL DGESVDQ('A','P','T','A','A',0,1,a,1,s,u,0,vt,0,ns,iw,1,w, &
     &                1,w,1,info)
         CALL CHKXER('DGESVDQ',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL DGESVDQ('A','P','T','A','A',1,1,a,0,s,u,0,vt,0,ns,iw,1,w, &
     &                1,w,1,info)
         CALL CHKXER('DGESVDQ',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL DGESVDQ('A','P','T','A','A',1,1,a,1,s,u,-1,vt,0,ns,iw,1,w,&
     &                1,w,1,info)
         CALL CHKXER('DGESVDQ',INFot,NOUt,LERr,OK)
         INFot = 14
         CALL DGESVDQ('A','P','T','A','A',1,1,a,1,s,u,1,vt,-1,ns,iw,1,w,&
     &                1,w,1,info)
         CALL CHKXER('DGESVDQ',INFot,NOUt,LERr,OK)
         INFot = 17
         CALL DGESVDQ('A','P','T','A','A',1,1,a,1,s,u,1,vt,1,ns,iw,-5,w,&
     &                1,w,1,info)
         CALL CHKXER('DGESVDQ',INFot,NOUt,LERr,OK)
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
!     End of DERRED
      END SUBROUTINE DERRED
