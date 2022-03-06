!*==derrbd.f90  processed by SPAG 7.51RB at 17:37 on  4 Mar 2022
!> \brief \b DERRBD
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DERRBD( PATH, NUNIT )
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
!> DERRBD tests the error exits for DGEBD2, DGEBRD, DORGBR, DORMBR,
!> DBDSQR, DBDSDC and DBDSVDX.
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
      SUBROUTINE DERRBD(Path,Nunit)
      IMPLICIT NONE
!*--DERRBD59
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
      DOUBLE PRECISION ZERO , ONE
      PARAMETER (ZERO=0.0D0,ONE=1.0D0)
!     ..
!     .. Local Scalars ..
      CHARACTER*2 c2
      INTEGER i , info , j , ns , nt
!     ..
!     .. Local Arrays ..
      INTEGER iq(NMAX,NMAX) , iw(NMAX)
      DOUBLE PRECISION a(NMAX,NMAX) , d(NMAX) , e(NMAX) , q(NMAX,NMAX) ,&
     &                 s(NMAX) , tp(NMAX) , tq(NMAX) , u(NMAX,NMAX) ,   &
     &                 v(NMAX,NMAX) , w(LW)
!     ..
!     .. External Functions ..
      LOGICAL LSAMEN
      EXTERNAL LSAMEN
!     ..
!     .. External Subroutines ..
      EXTERNAL CHKXER , DBDSDC , DBDSQR , DBDSVDX , DGEBD2 , DGEBRD ,   &
     &         DORGBR , DORMBR
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
      INTRINSIC DBLE
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
      ENDDO
      OK = .TRUE.
      nt = 0
!
!     Test error exits of the SVD routines.
!
      IF ( LSAMEN(2,c2,'BD') ) THEN
!
!        DGEBRD
!
         SRNamt = 'DGEBRD'
         INFot = 1
         CALL DGEBRD(-1,0,a,1,d,e,tq,tp,w,1,info)
         CALL CHKXER('DGEBRD',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DGEBRD(0,-1,a,1,d,e,tq,tp,w,1,info)
         CALL CHKXER('DGEBRD',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL DGEBRD(2,1,a,1,d,e,tq,tp,w,2,info)
         CALL CHKXER('DGEBRD',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL DGEBRD(2,1,a,2,d,e,tq,tp,w,1,info)
         CALL CHKXER('DGEBRD',INFot,NOUt,LERr,OK)
         nt = nt + 4
!
!        DGEBD2
!
         SRNamt = 'DGEBD2'
         INFot = 1
         CALL DGEBD2(-1,0,a,1,d,e,tq,tp,w,info)
         CALL CHKXER('DGEBD2',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DGEBD2(0,-1,a,1,d,e,tq,tp,w,info)
         CALL CHKXER('DGEBD2',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL DGEBD2(2,1,a,1,d,e,tq,tp,w,info)
         CALL CHKXER('DGEBD2',INFot,NOUt,LERr,OK)
         nt = nt + 3
!
!        DORGBR
!
         SRNamt = 'DORGBR'
         INFot = 1
         CALL DORGBR('/',0,0,0,a,1,tq,w,1,info)
         CALL CHKXER('DORGBR',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DORGBR('Q',-1,0,0,a,1,tq,w,1,info)
         CALL CHKXER('DORGBR',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DORGBR('Q',0,-1,0,a,1,tq,w,1,info)
         CALL CHKXER('DORGBR',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DORGBR('Q',0,1,0,a,1,tq,w,1,info)
         CALL CHKXER('DORGBR',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DORGBR('Q',1,0,1,a,1,tq,w,1,info)
         CALL CHKXER('DORGBR',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DORGBR('P',1,0,0,a,1,tq,w,1,info)
         CALL CHKXER('DORGBR',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DORGBR('P',0,1,1,a,1,tq,w,1,info)
         CALL CHKXER('DORGBR',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL DORGBR('Q',0,0,-1,a,1,tq,w,1,info)
         CALL CHKXER('DORGBR',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL DORGBR('Q',2,1,1,a,1,tq,w,1,info)
         CALL CHKXER('DORGBR',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL DORGBR('Q',2,2,1,a,2,tq,w,1,info)
         CALL CHKXER('DORGBR',INFot,NOUt,LERr,OK)
         nt = nt + 10
!
!        DORMBR
!
         SRNamt = 'DORMBR'
         INFot = 1
         CALL DORMBR('/','L','T',0,0,0,a,1,tq,u,1,w,1,info)
         CALL CHKXER('DORMBR',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DORMBR('Q','/','T',0,0,0,a,1,tq,u,1,w,1,info)
         CALL CHKXER('DORMBR',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DORMBR('Q','L','/',0,0,0,a,1,tq,u,1,w,1,info)
         CALL CHKXER('DORMBR',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL DORMBR('Q','L','T',-1,0,0,a,1,tq,u,1,w,1,info)
         CALL CHKXER('DORMBR',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL DORMBR('Q','L','T',0,-1,0,a,1,tq,u,1,w,1,info)
         CALL CHKXER('DORMBR',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL DORMBR('Q','L','T',0,0,-1,a,1,tq,u,1,w,1,info)
         CALL CHKXER('DORMBR',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL DORMBR('Q','L','T',2,0,0,a,1,tq,u,2,w,1,info)
         CALL CHKXER('DORMBR',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL DORMBR('Q','R','T',0,2,0,a,1,tq,u,1,w,1,info)
         CALL CHKXER('DORMBR',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL DORMBR('P','L','T',2,0,2,a,1,tq,u,2,w,1,info)
         CALL CHKXER('DORMBR',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL DORMBR('P','R','T',0,2,2,a,1,tq,u,1,w,1,info)
         CALL CHKXER('DORMBR',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL DORMBR('Q','R','T',2,0,0,a,1,tq,u,1,w,1,info)
         CALL CHKXER('DORMBR',INFot,NOUt,LERr,OK)
         INFot = 13
         CALL DORMBR('Q','L','T',0,2,0,a,1,tq,u,1,w,1,info)
         CALL CHKXER('DORMBR',INFot,NOUt,LERr,OK)
         INFot = 13
         CALL DORMBR('Q','R','T',2,0,0,a,1,tq,u,2,w,1,info)
         CALL CHKXER('DORMBR',INFot,NOUt,LERr,OK)
         nt = nt + 13
!
!        DBDSQR
!
         SRNamt = 'DBDSQR'
         INFot = 1
         CALL DBDSQR('/',0,0,0,0,d,e,v,1,u,1,a,1,w,info)
         CALL CHKXER('DBDSQR',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DBDSQR('U',-1,0,0,0,d,e,v,1,u,1,a,1,w,info)
         CALL CHKXER('DBDSQR',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DBDSQR('U',0,-1,0,0,d,e,v,1,u,1,a,1,w,info)
         CALL CHKXER('DBDSQR',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL DBDSQR('U',0,0,-1,0,d,e,v,1,u,1,a,1,w,info)
         CALL CHKXER('DBDSQR',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL DBDSQR('U',0,0,0,-1,d,e,v,1,u,1,a,1,w,info)
         CALL CHKXER('DBDSQR',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL DBDSQR('U',2,1,0,0,d,e,v,1,u,1,a,1,w,info)
         CALL CHKXER('DBDSQR',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL DBDSQR('U',0,0,2,0,d,e,v,1,u,1,a,1,w,info)
         CALL CHKXER('DBDSQR',INFot,NOUt,LERr,OK)
         INFot = 13
         CALL DBDSQR('U',2,0,0,1,d,e,v,1,u,1,a,1,w,info)
         CALL CHKXER('DBDSQR',INFot,NOUt,LERr,OK)
         nt = nt + 8
!
!        DBDSDC
!
         SRNamt = 'DBDSDC'
         INFot = 1
         CALL DBDSDC('/','N',0,d,e,u,1,v,1,q,iq,w,iw,info)
         CALL CHKXER('DBDSDC',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DBDSDC('U','/',0,d,e,u,1,v,1,q,iq,w,iw,info)
         CALL CHKXER('DBDSDC',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DBDSDC('U','N',-1,d,e,u,1,v,1,q,iq,w,iw,info)
         CALL CHKXER('DBDSDC',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL DBDSDC('U','I',2,d,e,u,1,v,1,q,iq,w,iw,info)
         CALL CHKXER('DBDSDC',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL DBDSDC('U','I',2,d,e,u,2,v,1,q,iq,w,iw,info)
         CALL CHKXER('DBDSDC',INFot,NOUt,LERr,OK)
         nt = nt + 5
!
!        DBDSVDX
!
         SRNamt = 'DBDSVDX'
         INFot = 1
         CALL DBDSVDX('X','N','A',1,d,e,ZERO,ONE,0,0,ns,s,q,1,w,iw,info)
         CALL CHKXER('DBDSVDX',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DBDSVDX('U','X','A',1,d,e,ZERO,ONE,0,0,ns,s,q,1,w,iw,info)
         CALL CHKXER('DBDSVDX',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DBDSVDX('U','V','X',1,d,e,ZERO,ONE,0,0,ns,s,q,1,w,iw,info)
         CALL CHKXER('DBDSVDX',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL DBDSVDX('U','V','A',-1,d,e,ZERO,ONE,0,0,ns,s,q,1,w,iw,    &
     &                info)
         CALL CHKXER('DBDSVDX',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL DBDSVDX('U','V','V',2,d,e,-ONE,ZERO,0,0,ns,s,q,1,w,iw,    &
     &                info)
         CALL CHKXER('DBDSVDX',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL DBDSVDX('U','V','V',2,d,e,ONE,ZERO,0,0,ns,s,q,1,w,iw,info)
         CALL CHKXER('DBDSVDX',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL DBDSVDX('L','V','I',2,d,e,ZERO,ZERO,0,2,ns,s,q,1,w,iw,    &
     &                info)
         CALL CHKXER('DBDSVDX',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL DBDSVDX('L','V','I',4,d,e,ZERO,ZERO,5,2,ns,s,q,1,w,iw,    &
     &                info)
         CALL CHKXER('DBDSVDX',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL DBDSVDX('L','V','I',4,d,e,ZERO,ZERO,3,2,ns,s,q,1,w,iw,    &
     &                info)
         CALL CHKXER('DBDSVDX',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL DBDSVDX('L','V','I',4,d,e,ZERO,ZERO,3,5,ns,s,q,1,w,iw,    &
     &                info)
         CALL CHKXER('DBDSVDX',INFot,NOUt,LERr,OK)
         INFot = 14
         CALL DBDSVDX('L','V','A',4,d,e,ZERO,ZERO,0,0,ns,s,q,0,w,iw,    &
     &                info)
         CALL CHKXER('DBDSVDX',INFot,NOUt,LERr,OK)
         INFot = 14
         CALL DBDSVDX('L','V','A',4,d,e,ZERO,ZERO,0,0,ns,s,q,2,w,iw,    &
     &                info)
         CALL CHKXER('DBDSVDX',INFot,NOUt,LERr,OK)
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
!     End of DERRBD
!
      END SUBROUTINE DERRBD
