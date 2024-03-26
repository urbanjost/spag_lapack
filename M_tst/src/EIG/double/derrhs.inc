!*==derrhs.f90  processed by SPAG 7.51RB at 17:37 on  4 Mar 2022
!> \brief \b derrhs
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE derrhs( PATH, NUNIT )
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
!> derrhs tests the error exits for DGEBAK, SGEBAL, SGEHRD, DORGHR,
!> DORMHR, DHSEQR, SHSEIN, and DTREVC.
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
!> \ingroup double_eig
!
!  =====================================================================
      SUBROUTINE derrhs(Path,Nunit)
      IMPLICIT NONE
!*--DERRHS59
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
      INTEGER NMAX , LW
      PARAMETER (NMAX=3,LW=(NMAX+2)*(NMAX+2)+NMAX)
!     ..
!     .. Local Scalars ..
      CHARACTER*2 c2
      INTEGER i , ihi , ilo , info , j , m , nt
!     ..
!     .. Local Arrays ..
      LOGICAL sel(NMAX)
      INTEGER ifaill(NMAX) , ifailr(NMAX)
      DOUBLE PRECISION a(NMAX,NMAX) , c(NMAX,NMAX) , s(NMAX) , tau(NMAX)&
     &                 , vl(NMAX,NMAX) , vr(NMAX,NMAX) , w(LW) ,        &
     &                 wi(NMAX) , wr(NMAX)
!     ..
!     .. External Functions ..
      LOGICAL LSAMEN
      EXTERNAL LSAMEN
!     ..
!     .. External Subroutines ..
      EXTERNAL DGEBAK , DGEBAL , DGEHRD , DHSEIN , DHSEQR ,    DORGHR , DORMHR , DTREVC
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC DBLE
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
            a(i,j) = 1.D0/DBLE(i+j)
         ENDDO
         wi(j) = DBLE(j)
         sel(j) = .TRUE.
      ENDDO
      OK = .TRUE.
      nt = 0
!
!     Test error exits of the nonsymmetric eigenvalue routines.
!
      IF ( LSAMEN(2,c2,'HS') ) THEN
!
!        DGEBAL
!
         SRNamt = 'DGEBAL'
         INFot = 1
         CALL DGEBAL('/',0,a,1,ilo,ihi,s,info)
         CALL chkxer('DGEBAL',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DGEBAL('N',-1,a,1,ilo,ihi,s,info)
         CALL chkxer('DGEBAL',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL DGEBAL('N',2,a,1,ilo,ihi,s,info)
         CALL chkxer('DGEBAL',INFot,NOUt,LERr,OK)
         nt = nt + 3
!
!        DGEBAK
!
         SRNamt = 'DGEBAK'
         INFot = 1
         CALL DGEBAK('/','R',0,1,0,s,0,a,1,info)
         CALL chkxer('DGEBAK',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DGEBAK('N','/',0,1,0,s,0,a,1,info)
         CALL chkxer('DGEBAK',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DGEBAK('N','R',-1,1,0,s,0,a,1,info)
         CALL chkxer('DGEBAK',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL DGEBAK('N','R',0,0,0,s,0,a,1,info)
         CALL chkxer('DGEBAK',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL DGEBAK('N','R',0,2,0,s,0,a,1,info)
         CALL chkxer('DGEBAK',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL DGEBAK('N','R',2,2,1,s,0,a,2,info)
         CALL chkxer('DGEBAK',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL DGEBAK('N','R',0,1,1,s,0,a,1,info)
         CALL chkxer('DGEBAK',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL DGEBAK('N','R',0,1,0,s,-1,a,1,info)
         CALL chkxer('DGEBAK',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL DGEBAK('N','R',2,1,2,s,0,a,1,info)
         CALL chkxer('DGEBAK',INFot,NOUt,LERr,OK)
         nt = nt + 9
!
!        DGEHRD
!
         SRNamt = 'DGEHRD'
         INFot = 1
         CALL DGEHRD(-1,1,1,a,1,tau,w,1,info)
         CALL chkxer('DGEHRD',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DGEHRD(0,0,0,a,1,tau,w,1,info)
         CALL chkxer('DGEHRD',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DGEHRD(0,2,0,a,1,tau,w,1,info)
         CALL chkxer('DGEHRD',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DGEHRD(1,1,0,a,1,tau,w,1,info)
         CALL chkxer('DGEHRD',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DGEHRD(0,1,1,a,1,tau,w,1,info)
         CALL chkxer('DGEHRD',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL DGEHRD(2,1,1,a,1,tau,w,2,info)
         CALL chkxer('DGEHRD',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL DGEHRD(2,1,2,a,2,tau,w,1,info)
         CALL chkxer('DGEHRD',INFot,NOUt,LERr,OK)
         nt = nt + 7
!
!        DORGHR
!
         SRNamt = 'DORGHR'
         INFot = 1
         CALL DORGHR(-1,1,1,a,1,tau,w,1,info)
         CALL chkxer('DORGHR',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DORGHR(0,0,0,a,1,tau,w,1,info)
         CALL chkxer('DORGHR',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DORGHR(0,2,0,a,1,tau,w,1,info)
         CALL chkxer('DORGHR',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DORGHR(1,1,0,a,1,tau,w,1,info)
         CALL chkxer('DORGHR',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DORGHR(0,1,1,a,1,tau,w,1,info)
         CALL chkxer('DORGHR',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL DORGHR(2,1,1,a,1,tau,w,1,info)
         CALL chkxer('DORGHR',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL DORGHR(3,1,3,a,3,tau,w,1,info)
         CALL chkxer('DORGHR',INFot,NOUt,LERr,OK)
         nt = nt + 7
!
!        DORMHR
!
         SRNamt = 'DORMHR'
         INFot = 1
         CALL DORMHR('/','N',0,0,1,0,a,1,tau,c,1,w,1,info)
         CALL chkxer('DORMHR',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DORMHR('L','/',0,0,1,0,a,1,tau,c,1,w,1,info)
         CALL chkxer('DORMHR',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DORMHR('L','N',-1,0,1,0,a,1,tau,c,1,w,1,info)
         CALL chkxer('DORMHR',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL DORMHR('L','N',0,-1,1,0,a,1,tau,c,1,w,1,info)
         CALL chkxer('DORMHR',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL DORMHR('L','N',0,0,0,0,a,1,tau,c,1,w,1,info)
         CALL chkxer('DORMHR',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL DORMHR('L','N',0,0,2,0,a,1,tau,c,1,w,1,info)
         CALL chkxer('DORMHR',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL DORMHR('L','N',1,2,2,1,a,1,tau,c,1,w,2,info)
         CALL chkxer('DORMHR',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL DORMHR('R','N',2,1,2,1,a,1,tau,c,2,w,2,info)
         CALL chkxer('DORMHR',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL DORMHR('L','N',1,1,1,0,a,1,tau,c,1,w,1,info)
         CALL chkxer('DORMHR',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL DORMHR('L','N',0,1,1,1,a,1,tau,c,1,w,1,info)
         CALL chkxer('DORMHR',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL DORMHR('R','N',1,0,1,1,a,1,tau,c,1,w,1,info)
         CALL chkxer('DORMHR',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL DORMHR('L','N',2,1,1,1,a,1,tau,c,2,w,1,info)
         CALL chkxer('DORMHR',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL DORMHR('R','N',1,2,1,1,a,1,tau,c,1,w,1,info)
         CALL chkxer('DORMHR',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL DORMHR('L','N',2,1,1,1,a,2,tau,c,1,w,1,info)
         CALL chkxer('DORMHR',INFot,NOUt,LERr,OK)
         INFot = 13
         CALL DORMHR('L','N',1,2,1,1,a,1,tau,c,1,w,1,info)
         CALL chkxer('DORMHR',INFot,NOUt,LERr,OK)
         INFot = 13
         CALL DORMHR('R','N',2,1,1,1,a,1,tau,c,2,w,1,info)
         CALL chkxer('DORMHR',INFot,NOUt,LERr,OK)
         nt = nt + 16
!
!        DHSEQR
!
         SRNamt = 'DHSEQR'
         INFot = 1
         CALL DHSEQR('/','N',0,1,0,a,1,wr,wi,c,1,w,1,info)
         CALL chkxer('DHSEQR',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DHSEQR('E','/',0,1,0,a,1,wr,wi,c,1,w,1,info)
         CALL chkxer('DHSEQR',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DHSEQR('E','N',-1,1,0,a,1,wr,wi,c,1,w,1,info)
         CALL chkxer('DHSEQR',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL DHSEQR('E','N',0,0,0,a,1,wr,wi,c,1,w,1,info)
         CALL chkxer('DHSEQR',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL DHSEQR('E','N',0,2,0,a,1,wr,wi,c,1,w,1,info)
         CALL chkxer('DHSEQR',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL DHSEQR('E','N',1,1,0,a,1,wr,wi,c,1,w,1,info)
         CALL chkxer('DHSEQR',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL DHSEQR('E','N',1,1,2,a,1,wr,wi,c,1,w,1,info)
         CALL chkxer('DHSEQR',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL DHSEQR('E','N',2,1,2,a,1,wr,wi,c,2,w,1,info)
         CALL chkxer('DHSEQR',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL DHSEQR('E','V',2,1,2,a,2,wr,wi,c,1,w,1,info)
         CALL chkxer('DHSEQR',INFot,NOUt,LERr,OK)
         nt = nt + 9
!
!        DHSEIN
!
         SRNamt = 'DHSEIN'
         INFot = 1
         CALL DHSEIN('/','N','N',sel,0,a,1,wr,wi,vl,1,vr,1,0,m,w,ifaill,&
     &               ifailr,info)
         CALL chkxer('DHSEIN',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DHSEIN('R','/','N',sel,0,a,1,wr,wi,vl,1,vr,1,0,m,w,ifaill,&
     &               ifailr,info)
         CALL chkxer('DHSEIN',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DHSEIN('R','N','/',sel,0,a,1,wr,wi,vl,1,vr,1,0,m,w,ifaill,&
     &               ifailr,info)
         CALL chkxer('DHSEIN',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL DHSEIN('R','N','N',sel,-1,a,1,wr,wi,vl,1,vr,1,0,m,w,      &
     &               ifaill,ifailr,info)
         CALL chkxer('DHSEIN',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL DHSEIN('R','N','N',sel,2,a,1,wr,wi,vl,1,vr,2,4,m,w,ifaill,&
     &               ifailr,info)
         CALL chkxer('DHSEIN',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL DHSEIN('L','N','N',sel,2,a,2,wr,wi,vl,1,vr,1,4,m,w,ifaill,&
     &               ifailr,info)
         CALL chkxer('DHSEIN',INFot,NOUt,LERr,OK)
         INFot = 13
         CALL DHSEIN('R','N','N',sel,2,a,2,wr,wi,vl,1,vr,1,4,m,w,ifaill,&
     &               ifailr,info)
         CALL chkxer('DHSEIN',INFot,NOUt,LERr,OK)
         INFot = 14
         CALL DHSEIN('R','N','N',sel,2,a,2,wr,wi,vl,1,vr,2,1,m,w,ifaill,&
     &               ifailr,info)
         CALL chkxer('DHSEIN',INFot,NOUt,LERr,OK)
         nt = nt + 8
!
!        DTREVC
!
         SRNamt = 'DTREVC'
         INFot = 1
         CALL DTREVC('/','A',sel,0,a,1,vl,1,vr,1,0,m,w,info)
         CALL chkxer('DTREVC',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DTREVC('L','/',sel,0,a,1,vl,1,vr,1,0,m,w,info)
         CALL chkxer('DTREVC',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL DTREVC('L','A',sel,-1,a,1,vl,1,vr,1,0,m,w,info)
         CALL chkxer('DTREVC',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL DTREVC('L','A',sel,2,a,1,vl,2,vr,1,4,m,w,info)
         CALL chkxer('DTREVC',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL DTREVC('L','A',sel,2,a,2,vl,1,vr,1,4,m,w,info)
         CALL chkxer('DTREVC',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL DTREVC('R','A',sel,2,a,2,vl,1,vr,1,4,m,w,info)
         CALL chkxer('DTREVC',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL DTREVC('L','A',sel,2,a,2,vl,2,vr,1,1,m,w,info)
         CALL chkxer('DTREVC',INFot,NOUt,LERr,OK)
         nt = nt + 7
      ENDIF
!
!     Print a summary line.
!
      IF ( OK ) THEN
         WRITE (NOUt,FMT=99001) Path , nt
      ELSE
         WRITE (NOUt,FMT=99002) Path
      ENDIF
!
99001 FORMAT (1X,A3,' routines passed the tests of the error exits',    &
     &        ' (',I3,' tests done)')
99002 FORMAT (' *** ',A3,' routines failed the tests of the error ',    &
     &        'exits ***')
!
!
!     End of derrhs
!
      END SUBROUTINE derrhs
