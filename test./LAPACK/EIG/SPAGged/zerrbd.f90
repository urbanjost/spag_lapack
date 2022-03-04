!*==zerrbd.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b ZERRBD
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZERRBD( PATH, NUNIT )
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
!> ZERRBD tests the error exits for ZGEBRD, ZUNGBR, ZUNMBR, and ZBDSQR.
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
      SUBROUTINE ZERRBD(Path,Nunit)
      IMPLICIT NONE
!*--ZERRBD58
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
      PARAMETER (NMAX=4,LW=NMAX)
!     ..
!     .. Local Scalars ..
      CHARACTER*2 c2
      INTEGER i , info , j , nt
!     ..
!     .. Local Arrays ..
      DOUBLE PRECISION d(NMAX) , e(NMAX) , rw(4*NMAX)
      COMPLEX*16 a(NMAX,NMAX) , tp(NMAX) , tq(NMAX) , u(NMAX,NMAX) ,    &
     &           v(NMAX,NMAX) , w(LW)
!     ..
!     .. External Functions ..
      LOGICAL LSAMEN
      EXTERNAL LSAMEN
!     ..
!     .. External Subroutines ..
      EXTERNAL CHKXER , ZBDSQR , ZGEBRD , ZUNGBR , ZUNMBR
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
!        ZGEBRD
!
         SRNamt = 'ZGEBRD'
         INFot = 1
         CALL ZGEBRD(-1,0,a,1,d,e,tq,tp,w,1,info)
         CALL CHKXER('ZGEBRD',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZGEBRD(0,-1,a,1,d,e,tq,tp,w,1,info)
         CALL CHKXER('ZGEBRD',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZGEBRD(2,1,a,1,d,e,tq,tp,w,2,info)
         CALL CHKXER('ZGEBRD',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL ZGEBRD(2,1,a,2,d,e,tq,tp,w,1,info)
         CALL CHKXER('ZGEBRD',INFot,NOUt,LERr,OK)
         nt = nt + 4
!
!        ZUNGBR
!
         SRNamt = 'ZUNGBR'
         INFot = 1
         CALL ZUNGBR('/',0,0,0,a,1,tq,w,1,info)
         CALL CHKXER('ZUNGBR',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZUNGBR('Q',-1,0,0,a,1,tq,w,1,info)
         CALL CHKXER('ZUNGBR',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZUNGBR('Q',0,-1,0,a,1,tq,w,1,info)
         CALL CHKXER('ZUNGBR',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZUNGBR('Q',0,1,0,a,1,tq,w,1,info)
         CALL CHKXER('ZUNGBR',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZUNGBR('Q',1,0,1,a,1,tq,w,1,info)
         CALL CHKXER('ZUNGBR',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZUNGBR('P',1,0,0,a,1,tq,w,1,info)
         CALL CHKXER('ZUNGBR',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZUNGBR('P',0,1,1,a,1,tq,w,1,info)
         CALL CHKXER('ZUNGBR',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZUNGBR('Q',0,0,-1,a,1,tq,w,1,info)
         CALL CHKXER('ZUNGBR',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL ZUNGBR('Q',2,1,1,a,1,tq,w,1,info)
         CALL CHKXER('ZUNGBR',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL ZUNGBR('Q',2,2,1,a,2,tq,w,1,info)
         CALL CHKXER('ZUNGBR',INFot,NOUt,LERr,OK)
         nt = nt + 10
!
!        ZUNMBR
!
         SRNamt = 'ZUNMBR'
         INFot = 1
         CALL ZUNMBR('/','L','T',0,0,0,a,1,tq,u,1,w,1,info)
         CALL CHKXER('ZUNMBR',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZUNMBR('Q','/','T',0,0,0,a,1,tq,u,1,w,1,info)
         CALL CHKXER('ZUNMBR',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZUNMBR('Q','L','/',0,0,0,a,1,tq,u,1,w,1,info)
         CALL CHKXER('ZUNMBR',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZUNMBR('Q','L','C',-1,0,0,a,1,tq,u,1,w,1,info)
         CALL CHKXER('ZUNMBR',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL ZUNMBR('Q','L','C',0,-1,0,a,1,tq,u,1,w,1,info)
         CALL CHKXER('ZUNMBR',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL ZUNMBR('Q','L','C',0,0,-1,a,1,tq,u,1,w,1,info)
         CALL CHKXER('ZUNMBR',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL ZUNMBR('Q','L','C',2,0,0,a,1,tq,u,2,w,1,info)
         CALL CHKXER('ZUNMBR',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL ZUNMBR('Q','R','C',0,2,0,a,1,tq,u,1,w,1,info)
         CALL CHKXER('ZUNMBR',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL ZUNMBR('P','L','C',2,0,2,a,1,tq,u,2,w,1,info)
         CALL CHKXER('ZUNMBR',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL ZUNMBR('P','R','C',0,2,2,a,1,tq,u,1,w,1,info)
         CALL CHKXER('ZUNMBR',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL ZUNMBR('Q','R','C',2,0,0,a,1,tq,u,1,w,1,info)
         CALL CHKXER('ZUNMBR',INFot,NOUt,LERr,OK)
         INFot = 13
         CALL ZUNMBR('Q','L','C',0,2,0,a,1,tq,u,1,w,0,info)
         CALL CHKXER('ZUNMBR',INFot,NOUt,LERr,OK)
         INFot = 13
         CALL ZUNMBR('Q','R','C',2,0,0,a,1,tq,u,2,w,0,info)
         CALL CHKXER('ZUNMBR',INFot,NOUt,LERr,OK)
         nt = nt + 13
!
!        ZBDSQR
!
         SRNamt = 'ZBDSQR'
         INFot = 1
         CALL ZBDSQR('/',0,0,0,0,d,e,v,1,u,1,a,1,rw,info)
         CALL CHKXER('ZBDSQR',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZBDSQR('U',-1,0,0,0,d,e,v,1,u,1,a,1,rw,info)
         CALL CHKXER('ZBDSQR',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZBDSQR('U',0,-1,0,0,d,e,v,1,u,1,a,1,rw,info)
         CALL CHKXER('ZBDSQR',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZBDSQR('U',0,0,-1,0,d,e,v,1,u,1,a,1,rw,info)
         CALL CHKXER('ZBDSQR',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL ZBDSQR('U',0,0,0,-1,d,e,v,1,u,1,a,1,rw,info)
         CALL CHKXER('ZBDSQR',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL ZBDSQR('U',2,1,0,0,d,e,v,1,u,1,a,1,rw,info)
         CALL CHKXER('ZBDSQR',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL ZBDSQR('U',0,0,2,0,d,e,v,1,u,1,a,1,rw,info)
         CALL CHKXER('ZBDSQR',INFot,NOUt,LERr,OK)
         INFot = 13
         CALL ZBDSQR('U',2,0,0,1,d,e,v,1,u,1,a,1,rw,info)
         CALL CHKXER('ZBDSQR',INFot,NOUt,LERr,OK)
         nt = nt + 8
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
99001 FORMAT (1X,A3,' routines passed the tests of the error exits (',  &
     &        I3,' tests done)')
99002 FORMAT (' *** ',A3,' routines failed the tests of the error ',    &
     &        'exits ***')
!
!
!     End of ZERRBD
!
      END SUBROUTINE ZERRBD
