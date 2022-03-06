!*==zerrec.f90  processed by SPAG 7.51RB at 17:37 on  4 Mar 2022
!> \brief \b ZERREC
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZERREC( PATH, NUNIT )
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
!> ZERREC tests the error exits for the routines for eigen- condition
!> estimation for DOUBLE PRECISION matrices:
!>    ZTRSYL, ZTREXC, ZTRSNA and ZTRSEN.
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
      SUBROUTINE ZERREC(Path,Nunit)
      IMPLICIT NONE
!*--ZERREC60
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
      PARAMETER (NMAX=4,LW=NMAX*(NMAX+2))
      DOUBLE PRECISION ONE , ZERO
      PARAMETER (ONE=1.0D0,ZERO=0.0D0)
!     ..
!     .. Local Scalars ..
      INTEGER i , ifst , ilst , info , j , m , nt
      DOUBLE PRECISION scale
!     ..
!     .. Local Arrays ..
      LOGICAL sel(NMAX)
      DOUBLE PRECISION rw(LW) , s(NMAX) , sep(NMAX)
      COMPLEX*16 a(NMAX,NMAX) , b(NMAX,NMAX) , c(NMAX,NMAX) , work(LW) ,&
     &           x(NMAX)
!     ..
!     .. External Subroutines ..
      EXTERNAL CHKXER , ZTREXC , ZTRSEN , ZTRSNA , ZTRSYL
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
      OK = .TRUE.
      nt = 0
!
!     Initialize A, B and SEL
!
      DO j = 1 , NMAX
         DO i = 1 , NMAX
            a(i,j) = ZERO
            b(i,j) = ZERO
         ENDDO
      ENDDO
      DO i = 1 , NMAX
         a(i,i) = ONE
         sel(i) = .TRUE.
      ENDDO
!
!     Test ZTRSYL
!
      SRNamt = 'ZTRSYL'
      INFot = 1
      CALL ZTRSYL('X','N',1,0,0,a,1,b,1,c,1,scale,info)
      CALL CHKXER('ZTRSYL',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL ZTRSYL('N','X',1,0,0,a,1,b,1,c,1,scale,info)
      CALL CHKXER('ZTRSYL',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL ZTRSYL('N','N',0,0,0,a,1,b,1,c,1,scale,info)
      CALL CHKXER('ZTRSYL',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL ZTRSYL('N','N',1,-1,0,a,1,b,1,c,1,scale,info)
      CALL CHKXER('ZTRSYL',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL ZTRSYL('N','N',1,0,-1,a,1,b,1,c,1,scale,info)
      CALL CHKXER('ZTRSYL',INFot,NOUt,LERr,OK)
      INFot = 7
      CALL ZTRSYL('N','N',1,2,0,a,1,b,1,c,2,scale,info)
      CALL CHKXER('ZTRSYL',INFot,NOUt,LERr,OK)
      INFot = 9
      CALL ZTRSYL('N','N',1,0,2,a,1,b,1,c,1,scale,info)
      CALL CHKXER('ZTRSYL',INFot,NOUt,LERr,OK)
      INFot = 11
      CALL ZTRSYL('N','N',1,2,0,a,2,b,1,c,1,scale,info)
      CALL CHKXER('ZTRSYL',INFot,NOUt,LERr,OK)
      nt = nt + 8
!
!     Test ZTREXC
!
      SRNamt = 'ZTREXC'
      ifst = 1
      ilst = 1
      INFot = 1
      CALL ZTREXC('X',1,a,1,b,1,ifst,ilst,info)
      CALL CHKXER('ZTREXC',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL ZTREXC('N',-1,a,1,b,1,ifst,ilst,info)
      CALL CHKXER('ZTREXC',INFot,NOUt,LERr,OK)
      INFot = 4
      ilst = 2
      CALL ZTREXC('N',2,a,1,b,1,ifst,ilst,info)
      CALL CHKXER('ZTREXC',INFot,NOUt,LERr,OK)
      INFot = 6
      CALL ZTREXC('V',2,a,2,b,1,ifst,ilst,info)
      CALL CHKXER('ZTREXC',INFot,NOUt,LERr,OK)
      INFot = 7
      ifst = 0
      ilst = 1
      CALL ZTREXC('V',1,a,1,b,1,ifst,ilst,info)
      CALL CHKXER('ZTREXC',INFot,NOUt,LERr,OK)
      INFot = 7
      ifst = 2
      CALL ZTREXC('V',1,a,1,b,1,ifst,ilst,info)
      CALL CHKXER('ZTREXC',INFot,NOUt,LERr,OK)
      INFot = 8
      ifst = 1
      ilst = 0
      CALL ZTREXC('V',1,a,1,b,1,ifst,ilst,info)
      CALL CHKXER('ZTREXC',INFot,NOUt,LERr,OK)
      INFot = 8
      ilst = 2
      CALL ZTREXC('V',1,a,1,b,1,ifst,ilst,info)
      CALL CHKXER('ZTREXC',INFot,NOUt,LERr,OK)
      nt = nt + 8
!
!     Test ZTRSNA
!
      SRNamt = 'ZTRSNA'
      INFot = 1
      CALL ZTRSNA('X','A',sel,0,a,1,b,1,c,1,s,sep,1,m,work,1,rw,info)
      CALL CHKXER('ZTRSNA',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL ZTRSNA('B','X',sel,0,a,1,b,1,c,1,s,sep,1,m,work,1,rw,info)
      CALL CHKXER('ZTRSNA',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL ZTRSNA('B','A',sel,-1,a,1,b,1,c,1,s,sep,1,m,work,1,rw,info)
      CALL CHKXER('ZTRSNA',INFot,NOUt,LERr,OK)
      INFot = 6
      CALL ZTRSNA('V','A',sel,2,a,1,b,1,c,1,s,sep,2,m,work,2,rw,info)
      CALL CHKXER('ZTRSNA',INFot,NOUt,LERr,OK)
      INFot = 8
      CALL ZTRSNA('B','A',sel,2,a,2,b,1,c,2,s,sep,2,m,work,2,rw,info)
      CALL CHKXER('ZTRSNA',INFot,NOUt,LERr,OK)
      INFot = 10
      CALL ZTRSNA('B','A',sel,2,a,2,b,2,c,1,s,sep,2,m,work,2,rw,info)
      CALL CHKXER('ZTRSNA',INFot,NOUt,LERr,OK)
      INFot = 13
      CALL ZTRSNA('B','A',sel,1,a,1,b,1,c,1,s,sep,0,m,work,1,rw,info)
      CALL CHKXER('ZTRSNA',INFot,NOUt,LERr,OK)
      INFot = 13
      CALL ZTRSNA('B','S',sel,2,a,2,b,2,c,2,s,sep,1,m,work,1,rw,info)
      CALL CHKXER('ZTRSNA',INFot,NOUt,LERr,OK)
      INFot = 16
      CALL ZTRSNA('B','A',sel,2,a,2,b,2,c,2,s,sep,2,m,work,1,rw,info)
      CALL CHKXER('ZTRSNA',INFot,NOUt,LERr,OK)
      nt = nt + 9
!
!     Test ZTRSEN
!
      sel(1) = .FALSE.
      SRNamt = 'ZTRSEN'
      INFot = 1
      CALL ZTRSEN('X','N',sel,0,a,1,b,1,x,m,s(1),sep(1),work,1,info)
      CALL CHKXER('ZTRSEN',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL ZTRSEN('N','X',sel,0,a,1,b,1,x,m,s(1),sep(1),work,1,info)
      CALL CHKXER('ZTRSEN',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL ZTRSEN('N','N',sel,-1,a,1,b,1,x,m,s(1),sep(1),work,1,info)
      CALL CHKXER('ZTRSEN',INFot,NOUt,LERr,OK)
      INFot = 6
      CALL ZTRSEN('N','N',sel,2,a,1,b,1,x,m,s(1),sep(1),work,2,info)
      CALL CHKXER('ZTRSEN',INFot,NOUt,LERr,OK)
      INFot = 8
      CALL ZTRSEN('N','V',sel,2,a,2,b,1,x,m,s(1),sep(1),work,1,info)
      CALL CHKXER('ZTRSEN',INFot,NOUt,LERr,OK)
      INFot = 14
      CALL ZTRSEN('N','V',sel,2,a,2,b,2,x,m,s(1),sep(1),work,0,info)
      CALL CHKXER('ZTRSEN',INFot,NOUt,LERr,OK)
      INFot = 14
      CALL ZTRSEN('E','V',sel,3,a,3,b,3,x,m,s(1),sep(1),work,1,info)
      CALL CHKXER('ZTRSEN',INFot,NOUt,LERr,OK)
      INFot = 14
      CALL ZTRSEN('V','V',sel,3,a,3,b,3,x,m,s(1),sep(1),work,3,info)
      CALL CHKXER('ZTRSEN',INFot,NOUt,LERr,OK)
      nt = nt + 8
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
!     End of ZERREC
!
      END SUBROUTINE ZERREC
