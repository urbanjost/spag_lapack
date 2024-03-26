!*==serrhs.f90  processed by SPAG 7.51RB at 17:37 on  4 Mar 2022
!> \brief \b serrhs
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE serrhs( PATH, NUNIT )
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
!> serrhs tests the error exits for SGEBAK, SGEBAL, SGEHRD, SORGHR,
!> SORMHR, SHSEQR, SHSEIN, and STREVC.
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
!> \ingroup single_eig
!
!  =====================================================================
      SUBROUTINE serrhs(Path,Nunit)
      IMPLICIT NONE
!*--SERRHS59
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
      INTEGER i , ilo , ihi , info , j , m , nt
!     ..
!     .. Local Arrays ..
      LOGICAL sel(NMAX)
      INTEGER ifaill(NMAX) , ifailr(NMAX)
      REAL a(NMAX,NMAX) , c(NMAX,NMAX) , tau(NMAX) , vl(NMAX,NMAX) ,    &
     &     vr(NMAX,NMAX) , w(LW) , wi(NMAX) , wr(NMAX) , s(NMAX)
!     ..
!     .. External Functions ..
      LOGICAL LSAMEN
      EXTERNAL LSAMEN
!     ..
!     .. External Subroutines ..
      EXTERNAL SGEBAK , SGEBAL , SGEHRD , SHSEIN , SHSEQR ,  SORGHR , SORMHR , STREVC
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC REAL
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
            a(i,j) = 1./REAL(i+j)
         ENDDO
         wi(j) = REAL(j)
         sel(j) = .TRUE.
      ENDDO
      OK = .TRUE.
      nt = 0
!
!     Test error exits of the nonsymmetric eigenvalue routines.
!
      IF ( LSAMEN(2,c2,'HS') ) THEN
!
!        SGEBAL
!
         SRNamt = 'SGEBAL'
         INFot = 1
         CALL SGEBAL('/',0,a,1,ilo,ihi,s,info)
         CALL chkxer('SGEBAL',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SGEBAL('N',-1,a,1,ilo,ihi,s,info)
         CALL chkxer('SGEBAL',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL SGEBAL('N',2,a,1,ilo,ihi,s,info)
         CALL chkxer('SGEBAL',INFot,NOUt,LERr,OK)
         nt = nt + 3
!
!        SGEBAK
!
         SRNamt = 'SGEBAK'
         INFot = 1
         CALL SGEBAK('/','R',0,1,0,s,0,a,1,info)
         CALL chkxer('SGEBAK',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SGEBAK('N','/',0,1,0,s,0,a,1,info)
         CALL chkxer('SGEBAK',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SGEBAK('N','R',-1,1,0,s,0,a,1,info)
         CALL chkxer('SGEBAK',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL SGEBAK('N','R',0,0,0,s,0,a,1,info)
         CALL chkxer('SGEBAK',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL SGEBAK('N','R',0,2,0,s,0,a,1,info)
         CALL chkxer('SGEBAK',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL SGEBAK('N','R',2,2,1,s,0,a,2,info)
         CALL chkxer('SGEBAK',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL SGEBAK('N','R',0,1,1,s,0,a,1,info)
         CALL chkxer('SGEBAK',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL SGEBAK('N','R',0,1,0,s,-1,a,1,info)
         CALL chkxer('SGEBAK',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL SGEBAK('N','R',2,1,2,s,0,a,1,info)
         CALL chkxer('SGEBAK',INFot,NOUt,LERr,OK)
         nt = nt + 9
!
!        SGEHRD
!
         SRNamt = 'SGEHRD'
         INFot = 1
         CALL SGEHRD(-1,1,1,a,1,tau,w,1,info)
         CALL chkxer('SGEHRD',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SGEHRD(0,0,0,a,1,tau,w,1,info)
         CALL chkxer('SGEHRD',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SGEHRD(0,2,0,a,1,tau,w,1,info)
         CALL chkxer('SGEHRD',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SGEHRD(1,1,0,a,1,tau,w,1,info)
         CALL chkxer('SGEHRD',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SGEHRD(0,1,1,a,1,tau,w,1,info)
         CALL chkxer('SGEHRD',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL SGEHRD(2,1,1,a,1,tau,w,2,info)
         CALL chkxer('SGEHRD',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL SGEHRD(2,1,2,a,2,tau,w,1,info)
         CALL chkxer('SGEHRD',INFot,NOUt,LERr,OK)
         nt = nt + 7
!
!        SORGHR
!
         SRNamt = 'SORGHR'
         INFot = 1
         CALL SORGHR(-1,1,1,a,1,tau,w,1,info)
         CALL chkxer('SORGHR',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SORGHR(0,0,0,a,1,tau,w,1,info)
         CALL chkxer('SORGHR',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SORGHR(0,2,0,a,1,tau,w,1,info)
         CALL chkxer('SORGHR',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SORGHR(1,1,0,a,1,tau,w,1,info)
         CALL chkxer('SORGHR',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SORGHR(0,1,1,a,1,tau,w,1,info)
         CALL chkxer('SORGHR',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL SORGHR(2,1,1,a,1,tau,w,1,info)
         CALL chkxer('SORGHR',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL SORGHR(3,1,3,a,3,tau,w,1,info)
         CALL chkxer('SORGHR',INFot,NOUt,LERr,OK)
         nt = nt + 7
!
!        SORMHR
!
         SRNamt = 'SORMHR'
         INFot = 1
         CALL SORMHR('/','N',0,0,1,0,a,1,tau,c,1,w,1,info)
         CALL chkxer('SORMHR',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SORMHR('L','/',0,0,1,0,a,1,tau,c,1,w,1,info)
         CALL chkxer('SORMHR',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SORMHR('L','N',-1,0,1,0,a,1,tau,c,1,w,1,info)
         CALL chkxer('SORMHR',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL SORMHR('L','N',0,-1,1,0,a,1,tau,c,1,w,1,info)
         CALL chkxer('SORMHR',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL SORMHR('L','N',0,0,0,0,a,1,tau,c,1,w,1,info)
         CALL chkxer('SORMHR',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL SORMHR('L','N',0,0,2,0,a,1,tau,c,1,w,1,info)
         CALL chkxer('SORMHR',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL SORMHR('L','N',1,2,2,1,a,1,tau,c,1,w,2,info)
         CALL chkxer('SORMHR',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL SORMHR('R','N',2,1,2,1,a,1,tau,c,2,w,2,info)
         CALL chkxer('SORMHR',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL SORMHR('L','N',1,1,1,0,a,1,tau,c,1,w,1,info)
         CALL chkxer('SORMHR',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL SORMHR('L','N',0,1,1,1,a,1,tau,c,1,w,1,info)
         CALL chkxer('SORMHR',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL SORMHR('R','N',1,0,1,1,a,1,tau,c,1,w,1,info)
         CALL chkxer('SORMHR',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL SORMHR('L','N',2,1,1,1,a,1,tau,c,2,w,1,info)
         CALL chkxer('SORMHR',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL SORMHR('R','N',1,2,1,1,a,1,tau,c,1,w,1,info)
         CALL chkxer('SORMHR',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL SORMHR('L','N',2,1,1,1,a,2,tau,c,1,w,1,info)
         CALL chkxer('SORMHR',INFot,NOUt,LERr,OK)
         INFot = 13
         CALL SORMHR('L','N',1,2,1,1,a,1,tau,c,1,w,1,info)
         CALL chkxer('SORMHR',INFot,NOUt,LERr,OK)
         INFot = 13
         CALL SORMHR('R','N',2,1,1,1,a,1,tau,c,2,w,1,info)
         CALL chkxer('SORMHR',INFot,NOUt,LERr,OK)
         nt = nt + 16
!
!        SHSEQR
!
         SRNamt = 'SHSEQR'
         INFot = 1
         CALL SHSEQR('/','N',0,1,0,a,1,wr,wi,c,1,w,1,info)
         CALL chkxer('SHSEQR',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SHSEQR('E','/',0,1,0,a,1,wr,wi,c,1,w,1,info)
         CALL chkxer('SHSEQR',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SHSEQR('E','N',-1,1,0,a,1,wr,wi,c,1,w,1,info)
         CALL chkxer('SHSEQR',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL SHSEQR('E','N',0,0,0,a,1,wr,wi,c,1,w,1,info)
         CALL chkxer('SHSEQR',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL SHSEQR('E','N',0,2,0,a,1,wr,wi,c,1,w,1,info)
         CALL chkxer('SHSEQR',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL SHSEQR('E','N',1,1,0,a,1,wr,wi,c,1,w,1,info)
         CALL chkxer('SHSEQR',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL SHSEQR('E','N',1,1,2,a,1,wr,wi,c,1,w,1,info)
         CALL chkxer('SHSEQR',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL SHSEQR('E','N',2,1,2,a,1,wr,wi,c,2,w,1,info)
         CALL chkxer('SHSEQR',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL SHSEQR('E','V',2,1,2,a,2,wr,wi,c,1,w,1,info)
         CALL chkxer('SHSEQR',INFot,NOUt,LERr,OK)
         nt = nt + 9
!
!        SHSEIN
!
         SRNamt = 'SHSEIN'
         INFot = 1
         CALL SHSEIN('/','N','N',sel,0,a,1,wr,wi,vl,1,vr,1,0,m,w,ifaill,&
     &               ifailr,info)
         CALL chkxer('SHSEIN',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SHSEIN('R','/','N',sel,0,a,1,wr,wi,vl,1,vr,1,0,m,w,ifaill,&
     &               ifailr,info)
         CALL chkxer('SHSEIN',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SHSEIN('R','N','/',sel,0,a,1,wr,wi,vl,1,vr,1,0,m,w,ifaill,&
     &               ifailr,info)
         CALL chkxer('SHSEIN',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL SHSEIN('R','N','N',sel,-1,a,1,wr,wi,vl,1,vr,1,0,m,w,      &
     &               ifaill,ifailr,info)
         CALL chkxer('SHSEIN',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL SHSEIN('R','N','N',sel,2,a,1,wr,wi,vl,1,vr,2,4,m,w,ifaill,&
     &               ifailr,info)
         CALL chkxer('SHSEIN',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL SHSEIN('L','N','N',sel,2,a,2,wr,wi,vl,1,vr,1,4,m,w,ifaill,&
     &               ifailr,info)
         CALL chkxer('SHSEIN',INFot,NOUt,LERr,OK)
         INFot = 13
         CALL SHSEIN('R','N','N',sel,2,a,2,wr,wi,vl,1,vr,1,4,m,w,ifaill,&
     &               ifailr,info)
         CALL chkxer('SHSEIN',INFot,NOUt,LERr,OK)
         INFot = 14
         CALL SHSEIN('R','N','N',sel,2,a,2,wr,wi,vl,1,vr,2,1,m,w,ifaill,&
     &               ifailr,info)
         CALL chkxer('SHSEIN',INFot,NOUt,LERr,OK)
         nt = nt + 8
!
!        STREVC
!
         SRNamt = 'STREVC'
         INFot = 1
         CALL STREVC('/','A',sel,0,a,1,vl,1,vr,1,0,m,w,info)
         CALL chkxer('STREVC',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL STREVC('L','/',sel,0,a,1,vl,1,vr,1,0,m,w,info)
         CALL chkxer('STREVC',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL STREVC('L','A',sel,-1,a,1,vl,1,vr,1,0,m,w,info)
         CALL chkxer('STREVC',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL STREVC('L','A',sel,2,a,1,vl,2,vr,1,4,m,w,info)
         CALL chkxer('STREVC',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL STREVC('L','A',sel,2,a,2,vl,1,vr,1,4,m,w,info)
         CALL chkxer('STREVC',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL STREVC('R','A',sel,2,a,2,vl,1,vr,1,4,m,w,info)
         CALL chkxer('STREVC',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL STREVC('L','A',sel,2,a,2,vl,2,vr,1,1,m,w,info)
         CALL chkxer('STREVC',INFot,NOUt,LERr,OK)
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
!     End of serrhs
!
      END SUBROUTINE serrhs
