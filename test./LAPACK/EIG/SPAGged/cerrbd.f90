!*==cerrbd.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b CERRBD
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE CERRBD( PATH, NUNIT )
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
!> CERRBD tests the error exits for CGEBRD, CUNGBR, CUNMBR, and CBDSQR.
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
      SUBROUTINE CERRBD(Path,Nunit)
      IMPLICIT NONE
!*--CERRBD58
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
      REAL d(NMAX) , e(NMAX) , rw(4*NMAX)
      COMPLEX a(NMAX,NMAX) , tp(NMAX) , tq(NMAX) , u(NMAX,NMAX) ,       &
     &        v(NMAX,NMAX) , w(LW)
!     ..
!     .. External Functions ..
      LOGICAL LSAMEN
      EXTERNAL LSAMEN
!     ..
!     .. External Subroutines ..
      EXTERNAL CBDSQR , CGEBRD , CHKXER , CUNGBR , CUNMBR
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
            a(i,j) = 1./REAL(i+j)
         ENDDO
      ENDDO
      OK = .TRUE.
      nt = 0
!
!     Test error exits of the SVD routines.
!
      IF ( LSAMEN(2,c2,'BD') ) THEN
!
!        CGEBRD
!
         SRNamt = 'CGEBRD'
         INFot = 1
         CALL CGEBRD(-1,0,a,1,d,e,tq,tp,w,1,info)
         CALL CHKXER('CGEBRD',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CGEBRD(0,-1,a,1,d,e,tq,tp,w,1,info)
         CALL CHKXER('CGEBRD',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CGEBRD(2,1,a,1,d,e,tq,tp,w,2,info)
         CALL CHKXER('CGEBRD',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL CGEBRD(2,1,a,2,d,e,tq,tp,w,1,info)
         CALL CHKXER('CGEBRD',INFot,NOUt,LERr,OK)
         nt = nt + 4
!
!        CUNGBR
!
         SRNamt = 'CUNGBR'
         INFot = 1
         CALL CUNGBR('/',0,0,0,a,1,tq,w,1,info)
         CALL CHKXER('CUNGBR',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CUNGBR('Q',-1,0,0,a,1,tq,w,1,info)
         CALL CHKXER('CUNGBR',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CUNGBR('Q',0,-1,0,a,1,tq,w,1,info)
         CALL CHKXER('CUNGBR',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CUNGBR('Q',0,1,0,a,1,tq,w,1,info)
         CALL CHKXER('CUNGBR',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CUNGBR('Q',1,0,1,a,1,tq,w,1,info)
         CALL CHKXER('CUNGBR',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CUNGBR('P',1,0,0,a,1,tq,w,1,info)
         CALL CHKXER('CUNGBR',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CUNGBR('P',0,1,1,a,1,tq,w,1,info)
         CALL CHKXER('CUNGBR',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CUNGBR('Q',0,0,-1,a,1,tq,w,1,info)
         CALL CHKXER('CUNGBR',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL CUNGBR('Q',2,1,1,a,1,tq,w,1,info)
         CALL CHKXER('CUNGBR',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL CUNGBR('Q',2,2,1,a,2,tq,w,1,info)
         CALL CHKXER('CUNGBR',INFot,NOUt,LERr,OK)
         nt = nt + 10
!
!        CUNMBR
!
         SRNamt = 'CUNMBR'
         INFot = 1
         CALL CUNMBR('/','L','T',0,0,0,a,1,tq,u,1,w,1,info)
         CALL CHKXER('CUNMBR',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CUNMBR('Q','/','T',0,0,0,a,1,tq,u,1,w,1,info)
         CALL CHKXER('CUNMBR',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CUNMBR('Q','L','/',0,0,0,a,1,tq,u,1,w,1,info)
         CALL CHKXER('CUNMBR',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CUNMBR('Q','L','C',-1,0,0,a,1,tq,u,1,w,1,info)
         CALL CHKXER('CUNMBR',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL CUNMBR('Q','L','C',0,-1,0,a,1,tq,u,1,w,1,info)
         CALL CHKXER('CUNMBR',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL CUNMBR('Q','L','C',0,0,-1,a,1,tq,u,1,w,1,info)
         CALL CHKXER('CUNMBR',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL CUNMBR('Q','L','C',2,0,0,a,1,tq,u,2,w,1,info)
         CALL CHKXER('CUNMBR',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL CUNMBR('Q','R','C',0,2,0,a,1,tq,u,1,w,1,info)
         CALL CHKXER('CUNMBR',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL CUNMBR('P','L','C',2,0,2,a,1,tq,u,2,w,1,info)
         CALL CHKXER('CUNMBR',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL CUNMBR('P','R','C',0,2,2,a,1,tq,u,1,w,1,info)
         CALL CHKXER('CUNMBR',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL CUNMBR('Q','R','C',2,0,0,a,1,tq,u,1,w,1,info)
         CALL CHKXER('CUNMBR',INFot,NOUt,LERr,OK)
         INFot = 13
         CALL CUNMBR('Q','L','C',0,2,0,a,1,tq,u,1,w,0,info)
         CALL CHKXER('CUNMBR',INFot,NOUt,LERr,OK)
         INFot = 13
         CALL CUNMBR('Q','R','C',2,0,0,a,1,tq,u,2,w,0,info)
         CALL CHKXER('CUNMBR',INFot,NOUt,LERr,OK)
         nt = nt + 13
!
!        CBDSQR
!
         SRNamt = 'CBDSQR'
         INFot = 1
         CALL CBDSQR('/',0,0,0,0,d,e,v,1,u,1,a,1,rw,info)
         CALL CHKXER('CBDSQR',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CBDSQR('U',-1,0,0,0,d,e,v,1,u,1,a,1,rw,info)
         CALL CHKXER('CBDSQR',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CBDSQR('U',0,-1,0,0,d,e,v,1,u,1,a,1,rw,info)
         CALL CHKXER('CBDSQR',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CBDSQR('U',0,0,-1,0,d,e,v,1,u,1,a,1,rw,info)
         CALL CHKXER('CBDSQR',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL CBDSQR('U',0,0,0,-1,d,e,v,1,u,1,a,1,rw,info)
         CALL CHKXER('CBDSQR',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL CBDSQR('U',2,1,0,0,d,e,v,1,u,1,a,1,rw,info)
         CALL CHKXER('CBDSQR',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL CBDSQR('U',0,0,2,0,d,e,v,1,u,1,a,1,rw,info)
         CALL CHKXER('CBDSQR',INFot,NOUt,LERr,OK)
         INFot = 13
         CALL CBDSQR('U',2,0,0,1,d,e,v,1,u,1,a,1,rw,info)
         CALL CHKXER('CBDSQR',INFot,NOUt,LERr,OK)
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
!     End of CERRBD
!
      END SUBROUTINE CERRBD
