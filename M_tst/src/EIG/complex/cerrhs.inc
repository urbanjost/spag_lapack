!*==cerrhs.f90  processed by SPAG 7.51RB at 17:37 on  4 Mar 2022
!> \brief \b cerrhs
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE cerrhs( PATH, NUNIT )
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
!> cerrhs tests the error exits for CGEBAK, CGEBAL, CGEHRD, CUNGHR,
!> CUNMHR, CHSEQR, CHSEIN, and CTREVC.
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
!> \ingroup complex_eig
!
!  =====================================================================
      SUBROUTINE cerrhs(Path,Nunit)
      IMPLICIT NONE
!*--CERRHS59
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
      REAL rw(NMAX) , s(NMAX)
      COMPLEX a(NMAX,NMAX) , c(NMAX,NMAX) , tau(NMAX) , vl(NMAX,NMAX) , &
     &        vr(NMAX,NMAX) , w(LW) , x(NMAX)
!     ..
!     .. External Functions ..
      LOGICAL LSAMEN
      EXTERNAL LSAMEN
!     ..
!     .. External Subroutines ..
      EXTERNAL CGEBAK , CGEBAL , CGEHRD , CHSEIN , CHSEQR ,    CUNGHR , CUNMHR , CTREVC
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
         sel(j) = .TRUE.
      ENDDO
      OK = .TRUE.
      nt = 0
!
!     Test error exits of the nonsymmetric eigenvalue routines.
!
      IF ( LSAMEN(2,c2,'HS') ) THEN
!
!        CGEBAL
!
         SRNamt = 'CGEBAL'
         INFot = 1
         CALL CGEBAL('/',0,a,1,ilo,ihi,s,info)
         CALL chkxer('CGEBAL',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CGEBAL('N',-1,a,1,ilo,ihi,s,info)
         CALL chkxer('CGEBAL',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CGEBAL('N',2,a,1,ilo,ihi,s,info)
         CALL chkxer('CGEBAL',INFot,NOUt,LERr,OK)
         nt = nt + 3
!
!        CGEBAK
!
         SRNamt = 'CGEBAK'
         INFot = 1
         CALL CGEBAK('/','R',0,1,0,s,0,a,1,info)
         CALL chkxer('CGEBAK',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CGEBAK('N','/',0,1,0,s,0,a,1,info)
         CALL chkxer('CGEBAK',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CGEBAK('N','R',-1,1,0,s,0,a,1,info)
         CALL chkxer('CGEBAK',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CGEBAK('N','R',0,0,0,s,0,a,1,info)
         CALL chkxer('CGEBAK',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CGEBAK('N','R',0,2,0,s,0,a,1,info)
         CALL chkxer('CGEBAK',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL CGEBAK('N','R',2,2,1,s,0,a,2,info)
         CALL chkxer('CGEBAK',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL CGEBAK('N','R',0,1,1,s,0,a,1,info)
         CALL chkxer('CGEBAK',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL CGEBAK('N','R',0,1,0,s,-1,a,1,info)
         CALL chkxer('CGEBAK',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL CGEBAK('N','R',2,1,2,s,0,a,1,info)
         CALL chkxer('CGEBAK',INFot,NOUt,LERr,OK)
         nt = nt + 9
!
!        CGEHRD
!
         SRNamt = 'CGEHRD'
         INFot = 1
         CALL CGEHRD(-1,1,1,a,1,tau,w,1,info)
         CALL chkxer('CGEHRD',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CGEHRD(0,0,0,a,1,tau,w,1,info)
         CALL chkxer('CGEHRD',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CGEHRD(0,2,0,a,1,tau,w,1,info)
         CALL chkxer('CGEHRD',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CGEHRD(1,1,0,a,1,tau,w,1,info)
         CALL chkxer('CGEHRD',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CGEHRD(0,1,1,a,1,tau,w,1,info)
         CALL chkxer('CGEHRD',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL CGEHRD(2,1,1,a,1,tau,w,2,info)
         CALL chkxer('CGEHRD',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL CGEHRD(2,1,2,a,2,tau,w,1,info)
         CALL chkxer('CGEHRD',INFot,NOUt,LERr,OK)
         nt = nt + 7
!
!        CUNGHR
!
         SRNamt = 'CUNGHR'
         INFot = 1
         CALL CUNGHR(-1,1,1,a,1,tau,w,1,info)
         CALL chkxer('CUNGHR',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CUNGHR(0,0,0,a,1,tau,w,1,info)
         CALL chkxer('CUNGHR',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CUNGHR(0,2,0,a,1,tau,w,1,info)
         CALL chkxer('CUNGHR',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CUNGHR(1,1,0,a,1,tau,w,1,info)
         CALL chkxer('CUNGHR',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CUNGHR(0,1,1,a,1,tau,w,1,info)
         CALL chkxer('CUNGHR',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL CUNGHR(2,1,1,a,1,tau,w,1,info)
         CALL chkxer('CUNGHR',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL CUNGHR(3,1,3,a,3,tau,w,1,info)
         CALL chkxer('CUNGHR',INFot,NOUt,LERr,OK)
         nt = nt + 7
!
!        CUNMHR
!
         SRNamt = 'CUNMHR'
         INFot = 1
         CALL CUNMHR('/','N',0,0,1,0,a,1,tau,c,1,w,1,info)
         CALL chkxer('CUNMHR',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CUNMHR('L','/',0,0,1,0,a,1,tau,c,1,w,1,info)
         CALL chkxer('CUNMHR',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CUNMHR('L','N',-1,0,1,0,a,1,tau,c,1,w,1,info)
         CALL chkxer('CUNMHR',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CUNMHR('L','N',0,-1,1,0,a,1,tau,c,1,w,1,info)
         CALL chkxer('CUNMHR',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL CUNMHR('L','N',0,0,0,0,a,1,tau,c,1,w,1,info)
         CALL chkxer('CUNMHR',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL CUNMHR('L','N',0,0,2,0,a,1,tau,c,1,w,1,info)
         CALL chkxer('CUNMHR',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL CUNMHR('L','N',1,2,2,1,a,1,tau,c,1,w,2,info)
         CALL chkxer('CUNMHR',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL CUNMHR('R','N',2,1,2,1,a,1,tau,c,2,w,2,info)
         CALL chkxer('CUNMHR',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL CUNMHR('L','N',1,1,1,0,a,1,tau,c,1,w,1,info)
         CALL chkxer('CUNMHR',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL CUNMHR('L','N',0,1,1,1,a,1,tau,c,1,w,1,info)
         CALL chkxer('CUNMHR',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL CUNMHR('R','N',1,0,1,1,a,1,tau,c,1,w,1,info)
         CALL chkxer('CUNMHR',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL CUNMHR('L','N',2,1,1,1,a,1,tau,c,2,w,1,info)
         CALL chkxer('CUNMHR',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL CUNMHR('R','N',1,2,1,1,a,1,tau,c,1,w,1,info)
         CALL chkxer('CUNMHR',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL CUNMHR('L','N',2,1,1,1,a,2,tau,c,1,w,1,info)
         CALL chkxer('CUNMHR',INFot,NOUt,LERr,OK)
         INFot = 13
         CALL CUNMHR('L','N',1,2,1,1,a,1,tau,c,1,w,1,info)
         CALL chkxer('CUNMHR',INFot,NOUt,LERr,OK)
         INFot = 13
         CALL CUNMHR('R','N',2,1,1,1,a,1,tau,c,2,w,1,info)
         CALL chkxer('CUNMHR',INFot,NOUt,LERr,OK)
         nt = nt + 16
!
!        CHSEQR
!
         SRNamt = 'CHSEQR'
         INFot = 1
         CALL CHSEQR('/','N',0,1,0,a,1,x,c,1,w,1,info)
         CALL chkxer('CHSEQR',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CHSEQR('E','/',0,1,0,a,1,x,c,1,w,1,info)
         CALL chkxer('CHSEQR',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CHSEQR('E','N',-1,1,0,a,1,x,c,1,w,1,info)
         CALL chkxer('CHSEQR',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CHSEQR('E','N',0,0,0,a,1,x,c,1,w,1,info)
         CALL chkxer('CHSEQR',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CHSEQR('E','N',0,2,0,a,1,x,c,1,w,1,info)
         CALL chkxer('CHSEQR',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL CHSEQR('E','N',1,1,0,a,1,x,c,1,w,1,info)
         CALL chkxer('CHSEQR',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL CHSEQR('E','N',1,1,2,a,1,x,c,1,w,1,info)
         CALL chkxer('CHSEQR',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL CHSEQR('E','N',2,1,2,a,1,x,c,2,w,1,info)
         CALL chkxer('CHSEQR',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL CHSEQR('E','V',2,1,2,a,2,x,c,1,w,1,info)
         CALL chkxer('CHSEQR',INFot,NOUt,LERr,OK)
         nt = nt + 9
!
!        CHSEIN
!
         SRNamt = 'CHSEIN'
         INFot = 1
         CALL CHSEIN('/','N','N',sel,0,a,1,x,vl,1,vr,1,0,m,w,rw,ifaill, &
     &               ifailr,info)
         CALL chkxer('CHSEIN',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CHSEIN('R','/','N',sel,0,a,1,x,vl,1,vr,1,0,m,w,rw,ifaill, &
     &               ifailr,info)
         CALL chkxer('CHSEIN',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CHSEIN('R','N','/',sel,0,a,1,x,vl,1,vr,1,0,m,w,rw,ifaill, &
     &               ifailr,info)
         CALL chkxer('CHSEIN',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL CHSEIN('R','N','N',sel,-1,a,1,x,vl,1,vr,1,0,m,w,rw,ifaill,&
     &               ifailr,info)
         CALL chkxer('CHSEIN',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL CHSEIN('R','N','N',sel,2,a,1,x,vl,1,vr,2,4,m,w,rw,ifaill, &
     &               ifailr,info)
         CALL chkxer('CHSEIN',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL CHSEIN('L','N','N',sel,2,a,2,x,vl,1,vr,1,4,m,w,rw,ifaill, &
     &               ifailr,info)
         CALL chkxer('CHSEIN',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL CHSEIN('R','N','N',sel,2,a,2,x,vl,1,vr,1,4,m,w,rw,ifaill, &
     &               ifailr,info)
         CALL chkxer('CHSEIN',INFot,NOUt,LERr,OK)
         INFot = 13
         CALL CHSEIN('R','N','N',sel,2,a,2,x,vl,1,vr,2,1,m,w,rw,ifaill, &
     &               ifailr,info)
         CALL chkxer('CHSEIN',INFot,NOUt,LERr,OK)
         nt = nt + 8
!
!        CTREVC
!
         SRNamt = 'CTREVC'
         INFot = 1
         CALL CTREVC('/','A',sel,0,a,1,vl,1,vr,1,0,m,w,rw,info)
         CALL chkxer('CTREVC',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CTREVC('L','/',sel,0,a,1,vl,1,vr,1,0,m,w,rw,info)
         CALL chkxer('CTREVC',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CTREVC('L','A',sel,-1,a,1,vl,1,vr,1,0,m,w,rw,info)
         CALL chkxer('CTREVC',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL CTREVC('L','A',sel,2,a,1,vl,2,vr,1,4,m,w,rw,info)
         CALL chkxer('CTREVC',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL CTREVC('L','A',sel,2,a,2,vl,1,vr,1,4,m,w,rw,info)
         CALL chkxer('CTREVC',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL CTREVC('R','A',sel,2,a,2,vl,1,vr,1,4,m,w,rw,info)
         CALL chkxer('CTREVC',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL CTREVC('L','A',sel,2,a,2,vl,2,vr,1,1,m,w,rw,info)
         CALL chkxer('CTREVC',INFot,NOUt,LERr,OK)
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
!     End of cerrhs
!
      END SUBROUTINE cerrhs
