!*==zerrhs.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b ZERRHS
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZERRHS( PATH, NUNIT )
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
!> ZERRHS tests the error exits for ZGEBAK, CGEBAL, CGEHRD, ZUNGHR,
!> ZUNMHR, ZHSEQR, CHSEIN, and ZTREVC.
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
!> \ingroup complex16_eig
!
!  =====================================================================
      SUBROUTINE ZERRHS(Path,Nunit)
      IMPLICIT NONE
!*--ZERRHS59
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
      PARAMETER (NMAX=3,LW=NMAX*NMAX)
!     ..
!     .. Local Scalars ..
      CHARACTER*2 c2
      INTEGER i , ihi , ilo , info , j , m , nt
!     ..
!     .. Local Arrays ..
      LOGICAL sel(NMAX)
      INTEGER ifaill(NMAX) , ifailr(NMAX)
      DOUBLE PRECISION rw(NMAX) , s(NMAX)
      COMPLEX*16 a(NMAX,NMAX) , c(NMAX,NMAX) , tau(NMAX) , vl(NMAX,NMAX)&
     &           , vr(NMAX,NMAX) , w(LW) , x(NMAX)
!     ..
!     .. External Functions ..
      LOGICAL LSAMEN
      EXTERNAL LSAMEN
!     ..
!     .. External Subroutines ..
      EXTERNAL CHKXER , ZGEBAK , ZGEBAL , ZGEHRD , ZHSEIN , ZHSEQR ,    &
     &         ZTREVC , ZUNGHR , ZUNMHR
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
         sel(j) = .TRUE.
      ENDDO
      OK = .TRUE.
      nt = 0
!
!     Test error exits of the nonsymmetric eigenvalue routines.
!
      IF ( LSAMEN(2,c2,'HS') ) THEN
!
!        ZGEBAL
!
         SRNamt = 'ZGEBAL'
         INFot = 1
         CALL ZGEBAL('/',0,a,1,ilo,ihi,s,info)
         CALL CHKXER('ZGEBAL',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZGEBAL('N',-1,a,1,ilo,ihi,s,info)
         CALL CHKXER('ZGEBAL',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZGEBAL('N',2,a,1,ilo,ihi,s,info)
         CALL CHKXER('ZGEBAL',INFot,NOUt,LERr,OK)
         nt = nt + 3
!
!        ZGEBAK
!
         SRNamt = 'ZGEBAK'
         INFot = 1
         CALL ZGEBAK('/','R',0,1,0,s,0,a,1,info)
         CALL CHKXER('ZGEBAK',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZGEBAK('N','/',0,1,0,s,0,a,1,info)
         CALL CHKXER('ZGEBAK',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZGEBAK('N','R',-1,1,0,s,0,a,1,info)
         CALL CHKXER('ZGEBAK',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZGEBAK('N','R',0,0,0,s,0,a,1,info)
         CALL CHKXER('ZGEBAK',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZGEBAK('N','R',0,2,0,s,0,a,1,info)
         CALL CHKXER('ZGEBAK',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL ZGEBAK('N','R',2,2,1,s,0,a,2,info)
         CALL CHKXER('ZGEBAK',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL ZGEBAK('N','R',0,1,1,s,0,a,1,info)
         CALL CHKXER('ZGEBAK',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL ZGEBAK('N','R',0,1,0,s,-1,a,1,info)
         CALL CHKXER('ZGEBAK',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL ZGEBAK('N','R',2,1,2,s,0,a,1,info)
         CALL CHKXER('ZGEBAK',INFot,NOUt,LERr,OK)
         nt = nt + 9
!
!        ZGEHRD
!
         SRNamt = 'ZGEHRD'
         INFot = 1
         CALL ZGEHRD(-1,1,1,a,1,tau,w,1,info)
         CALL CHKXER('ZGEHRD',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZGEHRD(0,0,0,a,1,tau,w,1,info)
         CALL CHKXER('ZGEHRD',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZGEHRD(0,2,0,a,1,tau,w,1,info)
         CALL CHKXER('ZGEHRD',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZGEHRD(1,1,0,a,1,tau,w,1,info)
         CALL CHKXER('ZGEHRD',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZGEHRD(0,1,1,a,1,tau,w,1,info)
         CALL CHKXER('ZGEHRD',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL ZGEHRD(2,1,1,a,1,tau,w,2,info)
         CALL CHKXER('ZGEHRD',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL ZGEHRD(2,1,2,a,2,tau,w,1,info)
         CALL CHKXER('ZGEHRD',INFot,NOUt,LERr,OK)
         nt = nt + 7
!
!        ZUNGHR
!
         SRNamt = 'ZUNGHR'
         INFot = 1
         CALL ZUNGHR(-1,1,1,a,1,tau,w,1,info)
         CALL CHKXER('ZUNGHR',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZUNGHR(0,0,0,a,1,tau,w,1,info)
         CALL CHKXER('ZUNGHR',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZUNGHR(0,2,0,a,1,tau,w,1,info)
         CALL CHKXER('ZUNGHR',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZUNGHR(1,1,0,a,1,tau,w,1,info)
         CALL CHKXER('ZUNGHR',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZUNGHR(0,1,1,a,1,tau,w,1,info)
         CALL CHKXER('ZUNGHR',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL ZUNGHR(2,1,1,a,1,tau,w,1,info)
         CALL CHKXER('ZUNGHR',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL ZUNGHR(3,1,3,a,3,tau,w,1,info)
         CALL CHKXER('ZUNGHR',INFot,NOUt,LERr,OK)
         nt = nt + 7
!
!        ZUNMHR
!
         SRNamt = 'ZUNMHR'
         INFot = 1
         CALL ZUNMHR('/','N',0,0,1,0,a,1,tau,c,1,w,1,info)
         CALL CHKXER('ZUNMHR',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZUNMHR('L','/',0,0,1,0,a,1,tau,c,1,w,1,info)
         CALL CHKXER('ZUNMHR',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZUNMHR('L','N',-1,0,1,0,a,1,tau,c,1,w,1,info)
         CALL CHKXER('ZUNMHR',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZUNMHR('L','N',0,-1,1,0,a,1,tau,c,1,w,1,info)
         CALL CHKXER('ZUNMHR',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL ZUNMHR('L','N',0,0,0,0,a,1,tau,c,1,w,1,info)
         CALL CHKXER('ZUNMHR',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL ZUNMHR('L','N',0,0,2,0,a,1,tau,c,1,w,1,info)
         CALL CHKXER('ZUNMHR',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL ZUNMHR('L','N',1,2,2,1,a,1,tau,c,1,w,2,info)
         CALL CHKXER('ZUNMHR',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL ZUNMHR('R','N',2,1,2,1,a,1,tau,c,2,w,2,info)
         CALL CHKXER('ZUNMHR',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL ZUNMHR('L','N',1,1,1,0,a,1,tau,c,1,w,1,info)
         CALL CHKXER('ZUNMHR',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL ZUNMHR('L','N',0,1,1,1,a,1,tau,c,1,w,1,info)
         CALL CHKXER('ZUNMHR',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL ZUNMHR('R','N',1,0,1,1,a,1,tau,c,1,w,1,info)
         CALL CHKXER('ZUNMHR',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL ZUNMHR('L','N',2,1,1,1,a,1,tau,c,2,w,1,info)
         CALL CHKXER('ZUNMHR',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL ZUNMHR('R','N',1,2,1,1,a,1,tau,c,1,w,1,info)
         CALL CHKXER('ZUNMHR',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL ZUNMHR('L','N',2,1,1,1,a,2,tau,c,1,w,1,info)
         CALL CHKXER('ZUNMHR',INFot,NOUt,LERr,OK)
         INFot = 13
         CALL ZUNMHR('L','N',1,2,1,1,a,1,tau,c,1,w,1,info)
         CALL CHKXER('ZUNMHR',INFot,NOUt,LERr,OK)
         INFot = 13
         CALL ZUNMHR('R','N',2,1,1,1,a,1,tau,c,2,w,1,info)
         CALL CHKXER('ZUNMHR',INFot,NOUt,LERr,OK)
         nt = nt + 16
!
!        ZHSEQR
!
         SRNamt = 'ZHSEQR'
         INFot = 1
         CALL ZHSEQR('/','N',0,1,0,a,1,x,c,1,w,1,info)
         CALL CHKXER('ZHSEQR',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZHSEQR('E','/',0,1,0,a,1,x,c,1,w,1,info)
         CALL CHKXER('ZHSEQR',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZHSEQR('E','N',-1,1,0,a,1,x,c,1,w,1,info)
         CALL CHKXER('ZHSEQR',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZHSEQR('E','N',0,0,0,a,1,x,c,1,w,1,info)
         CALL CHKXER('ZHSEQR',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZHSEQR('E','N',0,2,0,a,1,x,c,1,w,1,info)
         CALL CHKXER('ZHSEQR',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL ZHSEQR('E','N',1,1,0,a,1,x,c,1,w,1,info)
         CALL CHKXER('ZHSEQR',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL ZHSEQR('E','N',1,1,2,a,1,x,c,1,w,1,info)
         CALL CHKXER('ZHSEQR',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL ZHSEQR('E','N',2,1,2,a,1,x,c,2,w,1,info)
         CALL CHKXER('ZHSEQR',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL ZHSEQR('E','V',2,1,2,a,2,x,c,1,w,1,info)
         CALL CHKXER('ZHSEQR',INFot,NOUt,LERr,OK)
         nt = nt + 9
!
!        ZHSEIN
!
         SRNamt = 'ZHSEIN'
         INFot = 1
         CALL ZHSEIN('/','N','N',sel,0,a,1,x,vl,1,vr,1,0,m,w,rw,ifaill, &
     &               ifailr,info)
         CALL CHKXER('ZHSEIN',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZHSEIN('R','/','N',sel,0,a,1,x,vl,1,vr,1,0,m,w,rw,ifaill, &
     &               ifailr,info)
         CALL CHKXER('ZHSEIN',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZHSEIN('R','N','/',sel,0,a,1,x,vl,1,vr,1,0,m,w,rw,ifaill, &
     &               ifailr,info)
         CALL CHKXER('ZHSEIN',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL ZHSEIN('R','N','N',sel,-1,a,1,x,vl,1,vr,1,0,m,w,rw,ifaill,&
     &               ifailr,info)
         CALL CHKXER('ZHSEIN',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL ZHSEIN('R','N','N',sel,2,a,1,x,vl,1,vr,2,4,m,w,rw,ifaill, &
     &               ifailr,info)
         CALL CHKXER('ZHSEIN',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL ZHSEIN('L','N','N',sel,2,a,2,x,vl,1,vr,1,4,m,w,rw,ifaill, &
     &               ifailr,info)
         CALL CHKXER('ZHSEIN',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL ZHSEIN('R','N','N',sel,2,a,2,x,vl,1,vr,1,4,m,w,rw,ifaill, &
     &               ifailr,info)
         CALL CHKXER('ZHSEIN',INFot,NOUt,LERr,OK)
         INFot = 13
         CALL ZHSEIN('R','N','N',sel,2,a,2,x,vl,1,vr,2,1,m,w,rw,ifaill, &
     &               ifailr,info)
         CALL CHKXER('ZHSEIN',INFot,NOUt,LERr,OK)
         nt = nt + 8
!
!        ZTREVC
!
         SRNamt = 'ZTREVC'
         INFot = 1
         CALL ZTREVC('/','A',sel,0,a,1,vl,1,vr,1,0,m,w,rw,info)
         CALL CHKXER('ZTREVC',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZTREVC('L','/',sel,0,a,1,vl,1,vr,1,0,m,w,rw,info)
         CALL CHKXER('ZTREVC',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZTREVC('L','A',sel,-1,a,1,vl,1,vr,1,0,m,w,rw,info)
         CALL CHKXER('ZTREVC',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL ZTREVC('L','A',sel,2,a,1,vl,2,vr,1,4,m,w,rw,info)
         CALL CHKXER('ZTREVC',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL ZTREVC('L','A',sel,2,a,2,vl,1,vr,1,4,m,w,rw,info)
         CALL CHKXER('ZTREVC',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL ZTREVC('R','A',sel,2,a,2,vl,1,vr,1,4,m,w,rw,info)
         CALL CHKXER('ZTREVC',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL ZTREVC('L','A',sel,2,a,2,vl,2,vr,1,1,m,w,rw,info)
         CALL CHKXER('ZTREVC',INFot,NOUt,LERr,OK)
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
!     End of ZERRHS
!
      END SUBROUTINE ZERRHS
