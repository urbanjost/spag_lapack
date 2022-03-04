!*==cerred.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b CERRED
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE CERRED( PATH, NUNIT )
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
!> CERRED tests the error exits for the eigenvalue driver routines for
!> REAL matrices:
!>
!> PATH  driver   description
!> ----  ------   -----------
!> CEV   CGEEV    find eigenvalues/eigenvectors for nonsymmetric A
!> CES   CGEES    find eigenvalues/Schur form for nonsymmetric A
!> CVX   CGEEVX   CGEEV + balancing and condition estimation
!> CSX   CGEESX   CGEES + balancing and condition estimation
!> CBD   CGESVD   compute SVD of an M-by-N matrix A
!>       CGESDD   compute SVD of an M-by-N matrix A(by divide and
!>                conquer)
!>       CGEJSV   compute SVD of an M-by-N matrix A where M >= N
!>       CGESVDX  compute SVD of an M-by-N matrix A(by bisection
!>                and inverse iteration)
!>       CGESVDQ  compute SVD of an M-by-N matrix A(with a
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
!> \ingroup complex_eig
!
!  =====================================================================
      SUBROUTINE CERRED(Path,Nunit)
      IMPLICIT NONE
!*--CERRED74
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
      REAL ONE , ZERO
      PARAMETER (ONE=1.0E0,ZERO=0.0E0)
!     ..
!     .. Local Scalars ..
      CHARACTER*2 c2
      INTEGER i , ihi , ilo , info , j , ns , nt , sdim
      REAL abnrm
!     ..
!     .. Local Arrays ..
      LOGICAL b(NMAX)
      INTEGER iw(4*NMAX)
      REAL r1(NMAX) , r2(NMAX) , rw(LW) , s(NMAX)
      COMPLEX a(NMAX,NMAX) , u(NMAX,NMAX) , vl(NMAX,NMAX) ,             &
     &        vr(NMAX,NMAX) , vt(NMAX,NMAX) , w(10*NMAX) , x(NMAX)
!     ..
!     .. External Subroutines ..
      EXTERNAL CHKXER , CGEES , CGEESX , CGEEV , CGEEVX , CGEJSV ,      &
     &         CGESDD , CGESVD , CGESVDX , CGESVDQ
!     ..
!     .. External Functions ..
      LOGICAL LSAMEN , CSLECT
      EXTERNAL LSAMEN , CSLECT
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
!        Test CGEEV
!
         SRNamt = 'CGEEV '
         INFot = 1
         CALL CGEEV('X','N',0,a,1,x,vl,1,vr,1,w,1,rw,info)
         CALL CHKXER('CGEEV ',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CGEEV('N','X',0,a,1,x,vl,1,vr,1,w,1,rw,info)
         CALL CHKXER('CGEEV ',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CGEEV('N','N',-1,a,1,x,vl,1,vr,1,w,1,rw,info)
         CALL CHKXER('CGEEV ',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL CGEEV('N','N',2,a,1,x,vl,1,vr,1,w,4,rw,info)
         CALL CHKXER('CGEEV ',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL CGEEV('V','N',2,a,2,x,vl,1,vr,1,w,4,rw,info)
         CALL CHKXER('CGEEV ',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL CGEEV('N','V',2,a,2,x,vl,1,vr,1,w,4,rw,info)
         CALL CHKXER('CGEEV ',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL CGEEV('V','V',1,a,1,x,vl,1,vr,1,w,1,rw,info)
         CALL CHKXER('CGEEV ',INFot,NOUt,LERr,OK)
         nt = nt + 7
!
      ELSEIF ( LSAMEN(2,c2,'ES') ) THEN
!
!        Test CGEES
!
         SRNamt = 'CGEES '
         INFot = 1
         CALL CGEES('X','N',CSLECT,0,a,1,sdim,x,vl,1,w,1,rw,b,info)
         CALL CHKXER('CGEES ',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CGEES('N','X',CSLECT,0,a,1,sdim,x,vl,1,w,1,rw,b,info)
         CALL CHKXER('CGEES ',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CGEES('N','S',CSLECT,-1,a,1,sdim,x,vl,1,w,1,rw,b,info)
         CALL CHKXER('CGEES ',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL CGEES('N','S',CSLECT,2,a,1,sdim,x,vl,1,w,4,rw,b,info)
         CALL CHKXER('CGEES ',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL CGEES('V','S',CSLECT,2,a,2,sdim,x,vl,1,w,4,rw,b,info)
         CALL CHKXER('CGEES ',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL CGEES('N','S',CSLECT,1,a,1,sdim,x,vl,1,w,1,rw,b,info)
         CALL CHKXER('CGEES ',INFot,NOUt,LERr,OK)
         nt = nt + 6
!
      ELSEIF ( LSAMEN(2,c2,'VX') ) THEN
!
!        Test CGEEVX
!
         SRNamt = 'CGEEVX'
         INFot = 1
         CALL CGEEVX('X','N','N','N',0,a,1,x,vl,1,vr,1,ilo,ihi,s,abnrm, &
     &               r1,r2,w,1,rw,info)
         CALL CHKXER('CGEEVX',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CGEEVX('N','X','N','N',0,a,1,x,vl,1,vr,1,ilo,ihi,s,abnrm, &
     &               r1,r2,w,1,rw,info)
         CALL CHKXER('CGEEVX',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CGEEVX('N','N','X','N',0,a,1,x,vl,1,vr,1,ilo,ihi,s,abnrm, &
     &               r1,r2,w,1,rw,info)
         CALL CHKXER('CGEEVX',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CGEEVX('N','N','N','X',0,a,1,x,vl,1,vr,1,ilo,ihi,s,abnrm, &
     &               r1,r2,w,1,rw,info)
         CALL CHKXER('CGEEVX',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL CGEEVX('N','N','N','N',-1,a,1,x,vl,1,vr,1,ilo,ihi,s,abnrm,&
     &               r1,r2,w,1,rw,info)
         CALL CHKXER('CGEEVX',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL CGEEVX('N','N','N','N',2,a,1,x,vl,1,vr,1,ilo,ihi,s,abnrm, &
     &               r1,r2,w,4,rw,info)
         CALL CHKXER('CGEEVX',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL CGEEVX('N','V','N','N',2,a,2,x,vl,1,vr,1,ilo,ihi,s,abnrm, &
     &               r1,r2,w,4,rw,info)
         CALL CHKXER('CGEEVX',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL CGEEVX('N','N','V','N',2,a,2,x,vl,1,vr,1,ilo,ihi,s,abnrm, &
     &               r1,r2,w,4,rw,info)
         CALL CHKXER('CGEEVX',INFot,NOUt,LERr,OK)
         INFot = 20
         CALL CGEEVX('N','N','N','N',1,a,1,x,vl,1,vr,1,ilo,ihi,s,abnrm, &
     &               r1,r2,w,1,rw,info)
         CALL CHKXER('CGEEVX',INFot,NOUt,LERr,OK)
         INFot = 20
         CALL CGEEVX('N','N','V','V',1,a,1,x,vl,1,vr,1,ilo,ihi,s,abnrm, &
     &               r1,r2,w,2,rw,info)
         CALL CHKXER('CGEEVX',INFot,NOUt,LERr,OK)
         nt = nt + 10
!
      ELSEIF ( LSAMEN(2,c2,'SX') ) THEN
!
!        Test CGEESX
!
         SRNamt = 'CGEESX'
         INFot = 1
         CALL CGEESX('X','N',CSLECT,'N',0,a,1,sdim,x,vl,1,r1(1),r2(1),w,&
     &               1,rw,b,info)
         CALL CHKXER('CGEESX',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CGEESX('N','X',CSLECT,'N',0,a,1,sdim,x,vl,1,r1(1),r2(1),w,&
     &               1,rw,b,info)
         CALL CHKXER('CGEESX',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CGEESX('N','N',CSLECT,'X',0,a,1,sdim,x,vl,1,r1(1),r2(1),w,&
     &               1,rw,b,info)
         CALL CHKXER('CGEESX',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL CGEESX('N','N',CSLECT,'N',-1,a,1,sdim,x,vl,1,r1(1),r2(1), &
     &               w,1,rw,b,info)
         CALL CHKXER('CGEESX',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL CGEESX('N','N',CSLECT,'N',2,a,1,sdim,x,vl,1,r1(1),r2(1),w,&
     &               4,rw,b,info)
         CALL CHKXER('CGEESX',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL CGEESX('V','N',CSLECT,'N',2,a,2,sdim,x,vl,1,r1(1),r2(1),w,&
     &               4,rw,b,info)
         CALL CHKXER('CGEESX',INFot,NOUt,LERr,OK)
         INFot = 15
         CALL CGEESX('N','N',CSLECT,'N',1,a,1,sdim,x,vl,1,r1(1),r2(1),w,&
     &               1,rw,b,info)
         CALL CHKXER('CGEESX',INFot,NOUt,LERr,OK)
         nt = nt + 7
!
      ELSEIF ( LSAMEN(2,c2,'BD') ) THEN
!
!        Test CGESVD
!
         SRNamt = 'CGESVD'
         INFot = 1
         CALL CGESVD('X','N',0,0,a,1,s,u,1,vt,1,w,1,rw,info)
         CALL CHKXER('CGESVD',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CGESVD('N','X',0,0,a,1,s,u,1,vt,1,w,1,rw,info)
         CALL CHKXER('CGESVD',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CGESVD('O','O',0,0,a,1,s,u,1,vt,1,w,1,rw,info)
         CALL CHKXER('CGESVD',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CGESVD('N','N',-1,0,a,1,s,u,1,vt,1,w,1,rw,info)
         CALL CHKXER('CGESVD',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CGESVD('N','N',0,-1,a,1,s,u,1,vt,1,w,1,rw,info)
         CALL CHKXER('CGESVD',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL CGESVD('N','N',2,1,a,1,s,u,1,vt,1,w,5,rw,info)
         CALL CHKXER('CGESVD',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL CGESVD('A','N',2,1,a,2,s,u,1,vt,1,w,5,rw,info)
         CALL CHKXER('CGESVD',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL CGESVD('N','A',1,2,a,1,s,u,1,vt,1,w,5,rw,info)
         CALL CHKXER('CGESVD',INFot,NOUt,LERr,OK)
         nt = nt + 8
         IF ( OK ) THEN
            WRITE (NOUt,FMT=99001) SRNamt(1:LEN_TRIM(SRNamt)) , nt
         ELSE
            WRITE (NOUt,FMT=99002)
         ENDIF
!
!        Test CGESDD
!
         SRNamt = 'CGESDD'
         INFot = 1
         CALL CGESDD('X',0,0,a,1,s,u,1,vt,1,w,1,rw,iw,info)
         CALL CHKXER('CGESDD',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CGESDD('N',-1,0,a,1,s,u,1,vt,1,w,1,rw,iw,info)
         CALL CHKXER('CGESDD',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CGESDD('N',0,-1,a,1,s,u,1,vt,1,w,1,rw,iw,info)
         CALL CHKXER('CGESDD',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL CGESDD('N',2,1,a,1,s,u,1,vt,1,w,5,rw,iw,info)
         CALL CHKXER('CGESDD',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL CGESDD('A',2,1,a,2,s,u,1,vt,1,w,5,rw,iw,info)
         CALL CHKXER('CGESDD',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL CGESDD('A',1,2,a,1,s,u,1,vt,1,w,5,rw,iw,info)
         CALL CHKXER('CGESDD',INFot,NOUt,LERr,OK)
         nt = nt - 2
         IF ( OK ) THEN
            WRITE (NOUt,FMT=99001) SRNamt(1:LEN_TRIM(SRNamt)) , nt
         ELSE
            WRITE (NOUt,FMT=99002)
         ENDIF
!
!        Test CGEJSV
!
         SRNamt = 'CGEJSV'
         INFot = 1
         CALL CGEJSV('X','U','V','R','N','N',0,0,a,1,s,u,1,vt,1,w,1,rw, &
     &               1,iw,info)
         CALL CHKXER('CGEJSV',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CGEJSV('G','X','V','R','N','N',0,0,a,1,s,u,1,vt,1,w,1,rw, &
     &               1,iw,info)
         CALL CHKXER('CGEJSV',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CGEJSV('G','U','X','R','N','N',0,0,a,1,s,u,1,vt,1,w,1,rw, &
     &               1,iw,info)
         CALL CHKXER('CGEJSV',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CGEJSV('G','U','V','X','N','N',0,0,a,1,s,u,1,vt,1,w,1,rw, &
     &               1,iw,info)
         CALL CHKXER('CGEJSV',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL CGEJSV('G','U','V','R','X','N',0,0,a,1,s,u,1,vt,1,w,1,rw, &
     &               1,iw,info)
         CALL CHKXER('CGEJSV',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL CGEJSV('G','U','V','R','N','X',0,0,a,1,s,u,1,vt,1,w,1,rw, &
     &               1,iw,info)
         CALL CHKXER('CGEJSV',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL CGEJSV('G','U','V','R','N','N',-1,0,a,1,s,u,1,vt,1,w,1,rw,&
     &               1,iw,info)
         CALL CHKXER('CGEJSV',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL CGEJSV('G','U','V','R','N','N',0,-1,a,1,s,u,1,vt,1,w,1,rw,&
     &               1,iw,info)
         CALL CHKXER('CGEJSV',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL CGEJSV('G','U','V','R','N','N',2,1,a,1,s,u,1,vt,1,w,1,rw, &
     &               1,iw,info)
         CALL CHKXER('CGEJSV',INFot,NOUt,LERr,OK)
         INFot = 13
         CALL CGEJSV('G','U','V','R','N','N',2,2,a,2,s,u,1,vt,2,w,1,rw, &
     &               1,iw,info)
         CALL CHKXER('CGEJSV',INFot,NOUt,LERr,OK)
         INFot = 15
         CALL CGEJSV('G','U','V','R','N','N',2,2,a,2,s,u,2,vt,1,w,1,rw, &
     &               1,iw,info)
         CALL CHKXER('CGEJSV',INFot,NOUt,LERr,OK)
         nt = 11
         IF ( OK ) THEN
            WRITE (NOUt,FMT=99001) SRNamt(1:LEN_TRIM(SRNamt)) , nt
         ELSE
            WRITE (NOUt,FMT=99002)
         ENDIF
!
!        Test CGESVDX
!
         SRNamt = 'CGESVDX'
         INFot = 1
         CALL CGESVDX('X','N','A',0,0,a,1,ZERO,ZERO,0,0,ns,s,u,1,vt,1,w,&
     &                1,rw,iw,info)
         CALL CHKXER('CGESVDX',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CGESVDX('N','X','A',0,0,a,1,ZERO,ZERO,0,0,ns,s,u,1,vt,1,w,&
     &                1,rw,iw,info)
         CALL CHKXER('CGESVDX',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CGESVDX('N','N','X',0,0,a,1,ZERO,ZERO,0,0,ns,s,u,1,vt,1,w,&
     &                1,rw,iw,info)
         CALL CHKXER('CGESVDX',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CGESVDX('N','N','A',-1,0,a,1,ZERO,ZERO,0,0,ns,s,u,1,vt,1, &
     &                w,1,rw,iw,info)
         CALL CHKXER('CGESVDX',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL CGESVDX('N','N','A',0,-1,a,1,ZERO,ZERO,0,0,ns,s,u,1,vt,1, &
     &                w,1,rw,iw,info)
         CALL CHKXER('CGESVDX',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL CGESVDX('N','N','A',2,1,a,1,ZERO,ZERO,0,0,ns,s,u,1,vt,1,w,&
     &                1,rw,iw,info)
         CALL CHKXER('CGESVDX',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL CGESVDX('N','N','V',2,1,a,2,-ONE,ZERO,0,0,ns,s,u,1,vt,1,w,&
     &                1,rw,iw,info)
         CALL CHKXER('CGESVDX',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL CGESVDX('N','N','V',2,1,a,2,ONE,ZERO,0,0,ns,s,u,1,vt,1,w, &
     &                1,rw,iw,info)
         CALL CHKXER('CGESVDX',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL CGESVDX('N','N','I',2,2,a,2,ZERO,ZERO,0,1,ns,s,u,1,vt,1,w,&
     &                1,rw,iw,info)
         CALL CHKXER('CGESVDX',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL CGESVDX('V','N','I',2,2,a,2,ZERO,ZERO,1,0,ns,s,u,1,vt,1,w,&
     &                1,rw,iw,info)
         CALL CHKXER('CGESVDX',INFot,NOUt,LERr,OK)
         INFot = 15
         CALL CGESVDX('V','N','A',2,2,a,2,ZERO,ZERO,0,0,ns,s,u,1,vt,1,w,&
     &                1,rw,iw,info)
         CALL CHKXER('CGESVDX',INFot,NOUt,LERr,OK)
         INFot = 17
         CALL CGESVDX('N','V','A',2,2,a,2,ZERO,ZERO,0,0,ns,s,u,1,vt,1,w,&
     &                1,rw,iw,info)
         CALL CHKXER('CGESVDX',INFot,NOUt,LERr,OK)
         nt = 12
         IF ( OK ) THEN
            WRITE (NOUt,FMT=99001) SRNamt(1:LEN_TRIM(SRNamt)) , nt
         ELSE
            WRITE (NOUt,FMT=99002)
         ENDIF
!
!        Test CGESVDQ
!
         SRNamt = 'CGESVDQ'
         INFot = 1
         CALL CGESVDQ('X','P','T','A','A',0,0,a,1,s,u,0,vt,0,ns,iw,1,w, &
     &                1,rw,1,info)
         CALL CHKXER('CGESVDQ',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CGESVDQ('A','X','T','A','A',0,0,a,1,s,u,0,vt,0,ns,iw,1,w, &
     &                1,rw,1,info)
         CALL CHKXER('CGESVDQ',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CGESVDQ('A','P','X','A','A',0,0,a,1,s,u,0,vt,0,ns,iw,1,w, &
     &                1,rw,1,info)
         CALL CHKXER('CGESVDQ',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CGESVDQ('A','P','T','X','A',0,0,a,1,s,u,0,vt,0,ns,iw,1,w, &
     &                1,rw,1,info)
         CALL CHKXER('CGESVDQ',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL CGESVDQ('A','P','T','A','X',0,0,a,1,s,u,0,vt,0,ns,iw,1,w, &
     &                1,rw,1,info)
         CALL CHKXER('CGESVDQ',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL CGESVDQ('A','P','T','A','A',-1,0,a,1,s,u,0,vt,0,ns,iw,1,w,&
     &                1,rw,1,info)
         CALL CHKXER('CGESVDQ',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL CGESVDQ('A','P','T','A','A',0,1,a,1,s,u,0,vt,0,ns,iw,1,w, &
     &                1,rw,1,info)
         CALL CHKXER('CGESVDQ',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL CGESVDQ('A','P','T','A','A',1,1,a,0,s,u,0,vt,0,ns,iw,1,w, &
     &                1,rw,1,info)
         CALL CHKXER('CGESVDQ',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL CGESVDQ('A','P','T','A','A',1,1,a,1,s,u,-1,vt,0,ns,iw,1,w,&
     &                1,rw,1,info)
         CALL CHKXER('CGESVDQ',INFot,NOUt,LERr,OK)
         INFot = 14
         CALL CGESVDQ('A','P','T','A','A',1,1,a,1,s,u,1,vt,-1,ns,iw,1,w,&
     &                1,rw,1,info)
         CALL CHKXER('CGESVDQ',INFot,NOUt,LERr,OK)
         INFot = 17
         CALL CGESVDQ('A','P','T','A','A',1,1,a,1,s,u,1,vt,1,ns,iw,-5,w,&
     &                1,rw,1,info)
         CALL CHKXER('CGESVDQ',INFot,NOUt,LERr,OK)
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
!     End of CERRED
!
      END SUBROUTINE CERRED
