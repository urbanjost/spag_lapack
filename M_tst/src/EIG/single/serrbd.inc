!*==serrbd.f90  processed by SPAG 7.51RB at 17:37 on  4 Mar 2022
!> \brief \b serrbd
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE serrbd( PATH, NUNIT )
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
!> serrbd tests the error exits for SGEBD2, SGEBRD, SORGBR, SORMBR,
!> SBDSQR, SBDSDC and SBDSVDX.
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
      SUBROUTINE serrbd(Path,Nunit)
      IMPLICIT NONE
!*--SERRBD59
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
      PARAMETER (NMAX=4,LW=NMAX)
      REAL ZERO , ONE
      PARAMETER (ZERO=0.0E0,ONE=1.0E0)
!     ..
!     .. Local Scalars ..
      CHARACTER*2 c2
      INTEGER i , info , j , ns , nt
!     ..
!     .. Local Arrays ..
      INTEGER iq(NMAX,NMAX) , iw(NMAX)
      REAL a(NMAX,NMAX) , d(NMAX) , e(NMAX) , q(NMAX,NMAX) , s(NMAX) ,  &
     &     tp(NMAX) , tq(NMAX) , u(NMAX,NMAX) , v(NMAX,NMAX) , w(LW)
!     ..
!     .. External Functions ..
      LOGICAL LSAMEN
      EXTERNAL LSAMEN
!     ..
!     .. External Subroutines ..
      EXTERNAL SBDSDC , SBDSQR , SBDSVDX , SGEBD2 , SGEBRD , SORGBR , SORMBR
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
      INTRINSIC REAL
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
            a(i,j) = 1.D0/REAL(i+j)
         ENDDO
      ENDDO
      OK = .TRUE.
      nt = 0
!
!     Test error exits of the SVD routines.
!
      IF ( LSAMEN(2,c2,'BD') ) THEN
!
!        SGEBRD
!
         SRNamt = 'SGEBRD'
         INFot = 1
         CALL SGEBRD(-1,0,a,1,d,e,tq,tp,w,1,info)
         CALL chkxer('SGEBRD',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SGEBRD(0,-1,a,1,d,e,tq,tp,w,1,info)
         CALL chkxer('SGEBRD',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL SGEBRD(2,1,a,1,d,e,tq,tp,w,2,info)
         CALL chkxer('SGEBRD',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL SGEBRD(2,1,a,2,d,e,tq,tp,w,1,info)
         CALL chkxer('SGEBRD',INFot,NOUt,LERr,OK)
         nt = nt + 4
!
!        SGEBD2
!
         SRNamt = 'SGEBD2'
         INFot = 1
         CALL SGEBD2(-1,0,a,1,d,e,tq,tp,w,info)
         CALL chkxer('SGEBD2',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SGEBD2(0,-1,a,1,d,e,tq,tp,w,info)
         CALL chkxer('SGEBD2',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL SGEBD2(2,1,a,1,d,e,tq,tp,w,info)
         CALL chkxer('SGEBD2',INFot,NOUt,LERr,OK)
         nt = nt + 3
!
!        SORGBR
!
         SRNamt = 'SORGBR'
         INFot = 1
         CALL SORGBR('/',0,0,0,a,1,tq,w,1,info)
         CALL chkxer('SORGBR',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SORGBR('Q',-1,0,0,a,1,tq,w,1,info)
         CALL chkxer('SORGBR',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SORGBR('Q',0,-1,0,a,1,tq,w,1,info)
         CALL chkxer('SORGBR',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SORGBR('Q',0,1,0,a,1,tq,w,1,info)
         CALL chkxer('SORGBR',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SORGBR('Q',1,0,1,a,1,tq,w,1,info)
         CALL chkxer('SORGBR',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SORGBR('P',1,0,0,a,1,tq,w,1,info)
         CALL chkxer('SORGBR',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SORGBR('P',0,1,1,a,1,tq,w,1,info)
         CALL chkxer('SORGBR',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL SORGBR('Q',0,0,-1,a,1,tq,w,1,info)
         CALL chkxer('SORGBR',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL SORGBR('Q',2,1,1,a,1,tq,w,1,info)
         CALL chkxer('SORGBR',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL SORGBR('Q',2,2,1,a,2,tq,w,1,info)
         CALL chkxer('SORGBR',INFot,NOUt,LERr,OK)
         nt = nt + 10
!
!        SORMBR
!
         SRNamt = 'SORMBR'
         INFot = 1
         CALL SORMBR('/','L','T',0,0,0,a,1,tq,u,1,w,1,info)
         CALL chkxer('SORMBR',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SORMBR('Q','/','T',0,0,0,a,1,tq,u,1,w,1,info)
         CALL chkxer('SORMBR',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SORMBR('Q','L','/',0,0,0,a,1,tq,u,1,w,1,info)
         CALL chkxer('SORMBR',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL SORMBR('Q','L','T',-1,0,0,a,1,tq,u,1,w,1,info)
         CALL chkxer('SORMBR',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL SORMBR('Q','L','T',0,-1,0,a,1,tq,u,1,w,1,info)
         CALL chkxer('SORMBR',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL SORMBR('Q','L','T',0,0,-1,a,1,tq,u,1,w,1,info)
         CALL chkxer('SORMBR',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL SORMBR('Q','L','T',2,0,0,a,1,tq,u,2,w,1,info)
         CALL chkxer('SORMBR',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL SORMBR('Q','R','T',0,2,0,a,1,tq,u,1,w,1,info)
         CALL chkxer('SORMBR',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL SORMBR('P','L','T',2,0,2,a,1,tq,u,2,w,1,info)
         CALL chkxer('SORMBR',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL SORMBR('P','R','T',0,2,2,a,1,tq,u,1,w,1,info)
         CALL chkxer('SORMBR',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL SORMBR('Q','R','T',2,0,0,a,1,tq,u,1,w,1,info)
         CALL chkxer('SORMBR',INFot,NOUt,LERr,OK)
         INFot = 13
         CALL SORMBR('Q','L','T',0,2,0,a,1,tq,u,1,w,1,info)
         CALL chkxer('SORMBR',INFot,NOUt,LERr,OK)
         INFot = 13
         CALL SORMBR('Q','R','T',2,0,0,a,1,tq,u,2,w,1,info)
         CALL chkxer('SORMBR',INFot,NOUt,LERr,OK)
         nt = nt + 13
!
!        SBDSQR
!
         SRNamt = 'SBDSQR'
         INFot = 1
         CALL SBDSQR('/',0,0,0,0,d,e,v,1,u,1,a,1,w,info)
         CALL chkxer('SBDSQR',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SBDSQR('U',-1,0,0,0,d,e,v,1,u,1,a,1,w,info)
         CALL chkxer('SBDSQR',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SBDSQR('U',0,-1,0,0,d,e,v,1,u,1,a,1,w,info)
         CALL chkxer('SBDSQR',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL SBDSQR('U',0,0,-1,0,d,e,v,1,u,1,a,1,w,info)
         CALL chkxer('SBDSQR',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL SBDSQR('U',0,0,0,-1,d,e,v,1,u,1,a,1,w,info)
         CALL chkxer('SBDSQR',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL SBDSQR('U',2,1,0,0,d,e,v,1,u,1,a,1,w,info)
         CALL chkxer('SBDSQR',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL SBDSQR('U',0,0,2,0,d,e,v,1,u,1,a,1,w,info)
         CALL chkxer('SBDSQR',INFot,NOUt,LERr,OK)
         INFot = 13
         CALL SBDSQR('U',2,0,0,1,d,e,v,1,u,1,a,1,w,info)
         CALL chkxer('SBDSQR',INFot,NOUt,LERr,OK)
         nt = nt + 8
!
!        SBDSDC
!
         SRNamt = 'SBDSDC'
         INFot = 1
         CALL SBDSDC('/','N',0,d,e,u,1,v,1,q,iq,w,iw,info)
         CALL chkxer('SBDSDC',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SBDSDC('U','/',0,d,e,u,1,v,1,q,iq,w,iw,info)
         CALL chkxer('SBDSDC',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SBDSDC('U','N',-1,d,e,u,1,v,1,q,iq,w,iw,info)
         CALL chkxer('SBDSDC',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL SBDSDC('U','I',2,d,e,u,1,v,1,q,iq,w,iw,info)
         CALL chkxer('SBDSDC',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL SBDSDC('U','I',2,d,e,u,2,v,1,q,iq,w,iw,info)
         CALL chkxer('SBDSDC',INFot,NOUt,LERr,OK)
         nt = nt + 5
!
!        SBDSVDX
!
         SRNamt = 'SBDSVDX'
         INFot = 1
         CALL SBDSVDX('X','N','A',1,d,e,ZERO,ONE,0,0,ns,s,q,1,w,iw,info)
         CALL chkxer('SBDSVDX',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SBDSVDX('U','X','A',1,d,e,ZERO,ONE,0,0,ns,s,q,1,w,iw,info)
         CALL chkxer('SBDSVDX',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SBDSVDX('U','V','X',1,d,e,ZERO,ONE,0,0,ns,s,q,1,w,iw,info)
         CALL chkxer('SBDSVDX',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL SBDSVDX('U','V','A',-1,d,e,ZERO,ONE,0,0,ns,s,q,1,w,iw,    &
     &                info)
         CALL chkxer('SBDSVDX',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL SBDSVDX('U','V','V',2,d,e,-ONE,ZERO,0,0,ns,s,q,1,w,iw,    &
     &                info)
         CALL chkxer('SBDSVDX',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL SBDSVDX('U','V','V',2,d,e,ONE,ZERO,0,0,ns,s,q,1,w,iw,info)
         CALL chkxer('SBDSVDX',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL SBDSVDX('L','V','I',2,d,e,ZERO,ZERO,0,2,ns,s,q,1,w,iw,    &
     &                info)
         CALL chkxer('SBDSVDX',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL SBDSVDX('L','V','I',4,d,e,ZERO,ZERO,5,2,ns,s,q,1,w,iw,    &
     &                info)
         CALL chkxer('SBDSVDX',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL SBDSVDX('L','V','I',4,d,e,ZERO,ZERO,3,2,ns,s,q,1,w,iw,    &
     &                info)
         CALL chkxer('SBDSVDX',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL SBDSVDX('L','V','I',4,d,e,ZERO,ZERO,3,5,ns,s,q,1,w,iw,    &
     &                info)
         CALL chkxer('SBDSVDX',INFot,NOUt,LERr,OK)
         INFot = 14
         CALL SBDSVDX('L','V','A',4,d,e,ZERO,ZERO,0,0,ns,s,q,0,w,iw,    &
     &                info)
         CALL chkxer('SBDSVDX',INFot,NOUt,LERr,OK)
         INFot = 14
         CALL SBDSVDX('L','V','A',4,d,e,ZERO,ZERO,0,0,ns,s,q,2,w,iw,    &
     &                info)
         CALL chkxer('SBDSVDX',INFot,NOUt,LERr,OK)
         nt = nt + 12
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
!     End of serrbd
!
      END SUBROUTINE serrbd
