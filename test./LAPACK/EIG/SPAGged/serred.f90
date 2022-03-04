!*==serred.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b SERRED
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE SERRED( PATH, NUNIT )
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
!> SERRED tests the error exits for the eigenvalue driver routines for
!> REAL matrices:
!>
!> PATH  driver   description
!> ----  ------   -----------
!> SEV   SGEEV    find eigenvalues/eigenvectors for nonsymmetric A
!> SES   SGEES    find eigenvalues/Schur form for nonsymmetric A
!> SVX   SGEEVX   SGEEV + balancing and condition estimation
!> SSX   SGEESX   SGEES + balancing and condition estimation
!> SBD   SGESVD   compute SVD of an M-by-N matrix A
!>       SGESDD   compute SVD of an M-by-N matrix A (by divide and
!>                conquer)
!>       SGEJSV   compute SVD of an M-by-N matrix A where M >= N
!>       SGESVDX  compute SVD of an M-by-N matrix A(by bisection
!>                and inverse iteration)
!>       SGESVDQ  compute SVD of an M-by-N matrix A(with a
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
!> \ingroup single_eig
!
!  =====================================================================
      SUBROUTINE SERRED(Path,Nunit)
      IMPLICIT NONE
!*--SERRED74
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
      REAL ONE , ZERO
      PARAMETER (NMAX=4,ONE=1.0E0,ZERO=0.0E0)
!     ..
!     .. Local Scalars ..
      CHARACTER*2 c2
      INTEGER i , ihi , ilo , info , j , ns , nt , sdim
      REAL abnrm
!     ..
!     .. Local Arrays ..
      LOGICAL b(NMAX)
      INTEGER iw(2*NMAX)
      REAL a(NMAX,NMAX) , r1(NMAX) , r2(NMAX) , s(NMAX) , u(NMAX,NMAX) ,&
     &     vl(NMAX,NMAX) , vr(NMAX,NMAX) , vt(NMAX,NMAX) , w(10*NMAX) , &
     &     wi(NMAX) , wr(NMAX)
!     ..
!     .. External Subroutines ..
      EXTERNAL CHKXER , SGEES , SGEESX , SGEEV , SGEEVX , SGEJSV ,      &
     &         SGESDD , SGESVD , SGESVDX , SGESVDQ
!     ..
!     .. External Functions ..
      LOGICAL SSLECT , LSAMEN
      EXTERNAL SSLECT , LSAMEN
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC LEN_TRIM
!     ..
!     .. Arrays in Common ..
      LOGICAL SELval(20)
      REAL SELwi(20) , SELwr(20)
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
!        Test SGEEV
!
         SRNamt = 'SGEEV '
         INFot = 1
         CALL SGEEV('X','N',0,a,1,wr,wi,vl,1,vr,1,w,1,info)
         CALL CHKXER('SGEEV ',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SGEEV('N','X',0,a,1,wr,wi,vl,1,vr,1,w,1,info)
         CALL CHKXER('SGEEV ',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SGEEV('N','N',-1,a,1,wr,wi,vl,1,vr,1,w,1,info)
         CALL CHKXER('SGEEV ',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL SGEEV('N','N',2,a,1,wr,wi,vl,1,vr,1,w,6,info)
         CALL CHKXER('SGEEV ',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL SGEEV('V','N',2,a,2,wr,wi,vl,1,vr,1,w,8,info)
         CALL CHKXER('SGEEV ',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL SGEEV('N','V',2,a,2,wr,wi,vl,1,vr,1,w,8,info)
         CALL CHKXER('SGEEV ',INFot,NOUt,LERr,OK)
         INFot = 13
         CALL SGEEV('V','V',1,a,1,wr,wi,vl,1,vr,1,w,3,info)
         CALL CHKXER('SGEEV ',INFot,NOUt,LERr,OK)
         nt = nt + 7
!
      ELSEIF ( LSAMEN(2,c2,'ES') ) THEN
!
!        Test SGEES
!
         SRNamt = 'SGEES '
         INFot = 1
         CALL SGEES('X','N',SSLECT,0,a,1,sdim,wr,wi,vl,1,w,1,b,info)
         CALL CHKXER('SGEES ',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SGEES('N','X',SSLECT,0,a,1,sdim,wr,wi,vl,1,w,1,b,info)
         CALL CHKXER('SGEES ',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL SGEES('N','S',SSLECT,-1,a,1,sdim,wr,wi,vl,1,w,1,b,info)
         CALL CHKXER('SGEES ',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL SGEES('N','S',SSLECT,2,a,1,sdim,wr,wi,vl,1,w,6,b,info)
         CALL CHKXER('SGEES ',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL SGEES('V','S',SSLECT,2,a,2,sdim,wr,wi,vl,1,w,6,b,info)
         CALL CHKXER('SGEES ',INFot,NOUt,LERr,OK)
         INFot = 13
         CALL SGEES('N','S',SSLECT,1,a,1,sdim,wr,wi,vl,1,w,2,b,info)
         CALL CHKXER('SGEES ',INFot,NOUt,LERr,OK)
         nt = nt + 6
!
      ELSEIF ( LSAMEN(2,c2,'VX') ) THEN
!
!        Test SGEEVX
!
         SRNamt = 'SGEEVX'
         INFot = 1
         CALL SGEEVX('X','N','N','N',0,a,1,wr,wi,vl,1,vr,1,ilo,ihi,s,   &
     &               abnrm,r1,r2,w,1,iw,info)
         CALL CHKXER('SGEEVX',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SGEEVX('N','X','N','N',0,a,1,wr,wi,vl,1,vr,1,ilo,ihi,s,   &
     &               abnrm,r1,r2,w,1,iw,info)
         CALL CHKXER('SGEEVX',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SGEEVX('N','N','X','N',0,a,1,wr,wi,vl,1,vr,1,ilo,ihi,s,   &
     &               abnrm,r1,r2,w,1,iw,info)
         CALL CHKXER('SGEEVX',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL SGEEVX('N','N','N','X',0,a,1,wr,wi,vl,1,vr,1,ilo,ihi,s,   &
     &               abnrm,r1,r2,w,1,iw,info)
         CALL CHKXER('SGEEVX',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL SGEEVX('N','N','N','N',-1,a,1,wr,wi,vl,1,vr,1,ilo,ihi,s,  &
     &               abnrm,r1,r2,w,1,iw,info)
         CALL CHKXER('SGEEVX',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL SGEEVX('N','N','N','N',2,a,1,wr,wi,vl,1,vr,1,ilo,ihi,s,   &
     &               abnrm,r1,r2,w,1,iw,info)
         CALL CHKXER('SGEEVX',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL SGEEVX('N','V','N','N',2,a,2,wr,wi,vl,1,vr,1,ilo,ihi,s,   &
     &               abnrm,r1,r2,w,6,iw,info)
         CALL CHKXER('SGEEVX',INFot,NOUt,LERr,OK)
         INFot = 13
         CALL SGEEVX('N','N','V','N',2,a,2,wr,wi,vl,1,vr,1,ilo,ihi,s,   &
     &               abnrm,r1,r2,w,6,iw,info)
         CALL CHKXER('SGEEVX',INFot,NOUt,LERr,OK)
         INFot = 21
         CALL SGEEVX('N','N','N','N',1,a,1,wr,wi,vl,1,vr,1,ilo,ihi,s,   &
     &               abnrm,r1,r2,w,1,iw,info)
         CALL CHKXER('SGEEVX',INFot,NOUt,LERr,OK)
         INFot = 21
         CALL SGEEVX('N','V','N','N',1,a,1,wr,wi,vl,1,vr,1,ilo,ihi,s,   &
     &               abnrm,r1,r2,w,2,iw,info)
         CALL CHKXER('SGEEVX',INFot,NOUt,LERr,OK)
         INFot = 21
         CALL SGEEVX('N','N','V','V',1,a,1,wr,wi,vl,1,vr,1,ilo,ihi,s,   &
     &               abnrm,r1,r2,w,3,iw,info)
         CALL CHKXER('SGEEVX',INFot,NOUt,LERr,OK)
         nt = nt + 11
!
      ELSEIF ( LSAMEN(2,c2,'SX') ) THEN
!
!        Test SGEESX
!
         SRNamt = 'SGEESX'
         INFot = 1
         CALL SGEESX('X','N',SSLECT,'N',0,a,1,sdim,wr,wi,vl,1,r1(1),    &
     &               r2(1),w,1,iw,1,b,info)
         CALL CHKXER('SGEESX',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SGEESX('N','X',SSLECT,'N',0,a,1,sdim,wr,wi,vl,1,r1(1),    &
     &               r2(1),w,1,iw,1,b,info)
         CALL CHKXER('SGEESX',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL SGEESX('N','N',SSLECT,'X',0,a,1,sdim,wr,wi,vl,1,r1(1),    &
     &               r2(1),w,1,iw,1,b,info)
         CALL CHKXER('SGEESX',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL SGEESX('N','N',SSLECT,'N',-1,a,1,sdim,wr,wi,vl,1,r1(1),   &
     &               r2(1),w,1,iw,1,b,info)
         CALL CHKXER('SGEESX',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL SGEESX('N','N',SSLECT,'N',2,a,1,sdim,wr,wi,vl,1,r1(1),    &
     &               r2(1),w,6,iw,1,b,info)
         CALL CHKXER('SGEESX',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL SGEESX('V','N',SSLECT,'N',2,a,2,sdim,wr,wi,vl,1,r1(1),    &
     &               r2(1),w,6,iw,1,b,info)
         CALL CHKXER('SGEESX',INFot,NOUt,LERr,OK)
         INFot = 16
         CALL SGEESX('N','N',SSLECT,'N',1,a,1,sdim,wr,wi,vl,1,r1(1),    &
     &               r2(1),w,2,iw,1,b,info)
         CALL CHKXER('SGEESX',INFot,NOUt,LERr,OK)
         nt = nt + 7
!
      ELSEIF ( LSAMEN(2,c2,'BD') ) THEN
!
!        Test SGESVD
!
         SRNamt = 'SGESVD'
         INFot = 1
         CALL SGESVD('X','N',0,0,a,1,s,u,1,vt,1,w,1,info)
         CALL CHKXER('SGESVD',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SGESVD('N','X',0,0,a,1,s,u,1,vt,1,w,1,info)
         CALL CHKXER('SGESVD',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SGESVD('O','O',0,0,a,1,s,u,1,vt,1,w,1,info)
         CALL CHKXER('SGESVD',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SGESVD('N','N',-1,0,a,1,s,u,1,vt,1,w,1,info)
         CALL CHKXER('SGESVD',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL SGESVD('N','N',0,-1,a,1,s,u,1,vt,1,w,1,info)
         CALL CHKXER('SGESVD',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL SGESVD('N','N',2,1,a,1,s,u,1,vt,1,w,5,info)
         CALL CHKXER('SGESVD',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL SGESVD('A','N',2,1,a,2,s,u,1,vt,1,w,5,info)
         CALL CHKXER('SGESVD',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL SGESVD('N','A',1,2,a,1,s,u,1,vt,1,w,5,info)
         CALL CHKXER('SGESVD',INFot,NOUt,LERr,OK)
         nt = 8
         IF ( OK ) THEN
            WRITE (NOUt,FMT=99001) SRNamt(1:LEN_TRIM(SRNamt)) , nt
         ELSE
            WRITE (NOUt,FMT=99002)
         ENDIF
!
!        Test SGESDD
!
         SRNamt = 'SGESDD'
         INFot = 1
         CALL SGESDD('X',0,0,a,1,s,u,1,vt,1,w,1,iw,info)
         CALL CHKXER('SGESDD',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SGESDD('N',-1,0,a,1,s,u,1,vt,1,w,1,iw,info)
         CALL CHKXER('SGESDD',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SGESDD('N',0,-1,a,1,s,u,1,vt,1,w,1,iw,info)
         CALL CHKXER('SGESDD',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL SGESDD('N',2,1,a,1,s,u,1,vt,1,w,5,iw,info)
         CALL CHKXER('SGESDD',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL SGESDD('A',2,1,a,2,s,u,1,vt,1,w,5,iw,info)
         CALL CHKXER('SGESDD',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL SGESDD('A',1,2,a,1,s,u,1,vt,1,w,5,iw,info)
         CALL CHKXER('SGESDD',INFot,NOUt,LERr,OK)
         nt = 6
         IF ( OK ) THEN
            WRITE (NOUt,FMT=99001) SRNamt(1:LEN_TRIM(SRNamt)) , nt
         ELSE
            WRITE (NOUt,FMT=99002)
         ENDIF
!
!        Test SGEJSV
!
         SRNamt = 'SGEJSV'
         INFot = 1
         CALL SGEJSV('X','U','V','R','N','N',0,0,a,1,s,u,1,vt,1,w,1,iw, &
     &               info)
         CALL CHKXER('SGEJSV',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SGEJSV('G','X','V','R','N','N',0,0,a,1,s,u,1,vt,1,w,1,iw, &
     &               info)
         CALL CHKXER('SGEJSV',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SGEJSV('G','U','X','R','N','N',0,0,a,1,s,u,1,vt,1,w,1,iw, &
     &               info)
         CALL CHKXER('SGEJSV',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL SGEJSV('G','U','V','X','N','N',0,0,a,1,s,u,1,vt,1,w,1,iw, &
     &               info)
         CALL CHKXER('SGEJSV',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL SGEJSV('G','U','V','R','X','N',0,0,a,1,s,u,1,vt,1,w,1,iw, &
     &               info)
         CALL CHKXER('SGEJSV',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL SGEJSV('G','U','V','R','N','X',0,0,a,1,s,u,1,vt,1,w,1,iw, &
     &               info)
         CALL CHKXER('SGEJSV',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL SGEJSV('G','U','V','R','N','N',-1,0,a,1,s,u,1,vt,1,w,1,iw,&
     &               info)
         CALL CHKXER('SGEJSV',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL SGEJSV('G','U','V','R','N','N',0,-1,a,1,s,u,1,vt,1,w,1,iw,&
     &               info)
         CALL CHKXER('SGEJSV',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL SGEJSV('G','U','V','R','N','N',2,1,a,1,s,u,1,vt,1,w,1,iw, &
     &               info)
         CALL CHKXER('SGEJSV',INFot,NOUt,LERr,OK)
         INFot = 13
         CALL SGEJSV('G','U','V','R','N','N',2,2,a,2,s,u,1,vt,2,w,1,iw, &
     &               info)
         CALL CHKXER('SGEJSV',INFot,NOUt,LERr,OK)
         INFot = 15
         CALL SGEJSV('G','U','V','R','N','N',2,2,a,2,s,u,2,vt,1,w,1,iw, &
     &               info)
         CALL CHKXER('SGEJSV',INFot,NOUt,LERr,OK)
         nt = 11
         IF ( OK ) THEN
            WRITE (NOUt,FMT=99001) SRNamt(1:LEN_TRIM(SRNamt)) , nt
         ELSE
            WRITE (NOUt,FMT=99002)
         ENDIF
!
!        Test SGESVDX
!
         SRNamt = 'SGESVDX'
         INFot = 1
         CALL SGESVDX('X','N','A',0,0,a,1,ZERO,ZERO,0,0,ns,s,u,1,vt,1,w,&
     &                1,iw,info)
         CALL CHKXER('SGESVDX',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SGESVDX('N','X','A',0,0,a,1,ZERO,ZERO,0,0,ns,s,u,1,vt,1,w,&
     &                1,iw,info)
         CALL CHKXER('SGESVDX',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SGESVDX('N','N','X',0,0,a,1,ZERO,ZERO,0,0,ns,s,u,1,vt,1,w,&
     &                1,iw,info)
         CALL CHKXER('SGESVDX',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL SGESVDX('N','N','A',-1,0,a,1,ZERO,ZERO,0,0,ns,s,u,1,vt,1, &
     &                w,1,iw,info)
         CALL CHKXER('SGESVDX',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL SGESVDX('N','N','A',0,-1,a,1,ZERO,ZERO,0,0,ns,s,u,1,vt,1, &
     &                w,1,iw,info)
         CALL CHKXER('SGESVDX',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL SGESVDX('N','N','A',2,1,a,1,ZERO,ZERO,0,0,ns,s,u,1,vt,1,w,&
     &                1,iw,info)
         CALL CHKXER('SGESVDX',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL SGESVDX('N','N','V',2,1,a,2,-ONE,ZERO,0,0,ns,s,u,1,vt,1,w,&
     &                1,iw,info)
         CALL CHKXER('SGESVDX',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL SGESVDX('N','N','V',2,1,a,2,ONE,ZERO,0,0,ns,s,u,1,vt,1,w, &
     &                1,iw,info)
         CALL CHKXER('SGESVDX',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL SGESVDX('N','N','I',2,2,a,2,ZERO,ZERO,0,1,ns,s,u,1,vt,1,w,&
     &                1,iw,info)
         CALL CHKXER('SGESVDX',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL SGESVDX('V','N','I',2,2,a,2,ZERO,ZERO,1,0,ns,s,u,1,vt,1,w,&
     &                1,iw,info)
         CALL CHKXER('SGESVDX',INFot,NOUt,LERr,OK)
         INFot = 15
         CALL SGESVDX('V','N','A',2,2,a,2,ZERO,ZERO,0,0,ns,s,u,1,vt,1,w,&
     &                1,iw,info)
         CALL CHKXER('SGESVDX',INFot,NOUt,LERr,OK)
         INFot = 17
         CALL SGESVDX('N','V','A',2,2,a,2,ZERO,ZERO,0,0,ns,s,u,1,vt,1,w,&
     &                1,iw,info)
         CALL CHKXER('SGESVDX',INFot,NOUt,LERr,OK)
         nt = 12
         IF ( OK ) THEN
            WRITE (NOUt,FMT=99001) SRNamt(1:LEN_TRIM(SRNamt)) , nt
         ELSE
            WRITE (NOUt,FMT=99002)
         ENDIF
!
!        Test SGESVDQ
!
         SRNamt = 'SGESVDQ'
         INFot = 1
         CALL SGESVDQ('X','P','T','A','A',0,0,a,1,s,u,0,vt,0,ns,iw,1,w, &
     &                1,w,1,info)
         CALL CHKXER('SGESVDQ',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SGESVDQ('A','X','T','A','A',0,0,a,1,s,u,0,vt,0,ns,iw,1,w, &
     &                1,w,1,info)
         CALL CHKXER('SGESVDQ',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SGESVDQ('A','P','X','A','A',0,0,a,1,s,u,0,vt,0,ns,iw,1,w, &
     &                1,w,1,info)
         CALL CHKXER('SGESVDQ',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL SGESVDQ('A','P','T','X','A',0,0,a,1,s,u,0,vt,0,ns,iw,1,w, &
     &                1,w,1,info)
         CALL CHKXER('SGESVDQ',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL SGESVDQ('A','P','T','A','X',0,0,a,1,s,u,0,vt,0,ns,iw,1,w, &
     &                1,w,1,info)
         CALL CHKXER('SGESVDQ',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL SGESVDQ('A','P','T','A','A',-1,0,a,1,s,u,0,vt,0,ns,iw,1,w,&
     &                1,w,1,info)
         CALL CHKXER('SGESVDQ',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL SGESVDQ('A','P','T','A','A',0,1,a,1,s,u,0,vt,0,ns,iw,1,w, &
     &                1,w,1,info)
         CALL CHKXER('SGESVDQ',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL SGESVDQ('A','P','T','A','A',1,1,a,0,s,u,0,vt,0,ns,iw,1,w, &
     &                1,w,1,info)
         CALL CHKXER('SGESVDQ',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL SGESVDQ('A','P','T','A','A',1,1,a,1,s,u,-1,vt,0,ns,iw,1,w,&
     &                1,w,1,info)
         CALL CHKXER('SGESVDQ',INFot,NOUt,LERr,OK)
         INFot = 14
         CALL SGESVDQ('A','P','T','A','A',1,1,a,1,s,u,1,vt,-1,ns,iw,1,w,&
     &                1,w,1,info)
         CALL CHKXER('SGESVDQ',INFot,NOUt,LERr,OK)
         INFot = 17
         CALL SGESVDQ('A','P','T','A','A',1,1,a,1,s,u,1,vt,1,ns,iw,-5,w,&
     &                1,w,1,info)
         CALL CHKXER('SGESVDQ',INFot,NOUt,LERr,OK)
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
!     End of SERRED
!
      END SUBROUTINE SERRED
