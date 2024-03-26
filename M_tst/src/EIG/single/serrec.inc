!*==serrec.f90  processed by SPAG 7.51RB at 17:37 on  4 Mar 2022
!> \brief \b serrec
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE serrec( PATH, NUNIT )
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
!> serrec tests the error exits for the routines for eigen- condition
!> estimation for REAL matrices:
!>    STRSYL, STREXC, STRSNA and STRSEN.
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
      SUBROUTINE serrec(Path,Nunit)
      IMPLICIT NONE
!*--SERREC60
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
      INTEGER NMAX
      REAL ONE , ZERO
      PARAMETER (NMAX=4,ONE=1.0E0,ZERO=0.0E0)
!     ..
!     .. Local Scalars ..
      INTEGER i , ifst , ilst , info , j , m , nt
      REAL scale
!     ..
!     .. Local Arrays ..
      LOGICAL sel(NMAX)
      INTEGER iwork(NMAX)
      REAL a(NMAX,NMAX) , b(NMAX,NMAX) , c(NMAX,NMAX) , s(NMAX) ,       &
     &     sep(NMAX) , wi(NMAX) , work(NMAX) , wr(NMAX)
!     ..
!     .. External Subroutines ..
      EXTERNAL STREXC , STRSEN , STRSNA , STRSYL
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
!     Test STRSYL
!
      SRNamt = 'STRSYL'
      INFot = 1
      CALL STRSYL('X','N',1,0,0,a,1,b,1,c,1,scale,info)
      CALL chkxer('STRSYL',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL STRSYL('N','X',1,0,0,a,1,b,1,c,1,scale,info)
      CALL chkxer('STRSYL',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL STRSYL('N','N',0,0,0,a,1,b,1,c,1,scale,info)
      CALL chkxer('STRSYL',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL STRSYL('N','N',1,-1,0,a,1,b,1,c,1,scale,info)
      CALL chkxer('STRSYL',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL STRSYL('N','N',1,0,-1,a,1,b,1,c,1,scale,info)
      CALL chkxer('STRSYL',INFot,NOUt,LERr,OK)
      INFot = 7
      CALL STRSYL('N','N',1,2,0,a,1,b,1,c,2,scale,info)
      CALL chkxer('STRSYL',INFot,NOUt,LERr,OK)
      INFot = 9
      CALL STRSYL('N','N',1,0,2,a,1,b,1,c,1,scale,info)
      CALL chkxer('STRSYL',INFot,NOUt,LERr,OK)
      INFot = 11
      CALL STRSYL('N','N',1,2,0,a,2,b,1,c,1,scale,info)
      CALL chkxer('STRSYL',INFot,NOUt,LERr,OK)
      nt = nt + 8
!
!     Test STREXC
!
      SRNamt = 'STREXC'
      ifst = 1
      ilst = 1
      INFot = 1
      CALL STREXC('X',1,a,1,b,1,ifst,ilst,work,info)
      CALL chkxer('STREXC',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL STREXC('N',-1,a,1,b,1,ifst,ilst,work,info)
      CALL chkxer('STREXC',INFot,NOUt,LERr,OK)
      INFot = 4
      ilst = 2
      CALL STREXC('N',2,a,1,b,1,ifst,ilst,work,info)
      CALL chkxer('STREXC',INFot,NOUt,LERr,OK)
      INFot = 6
      CALL STREXC('V',2,a,2,b,1,ifst,ilst,work,info)
      CALL chkxer('STREXC',INFot,NOUt,LERr,OK)
      INFot = 7
      ifst = 0
      ilst = 1
      CALL STREXC('V',1,a,1,b,1,ifst,ilst,work,info)
      CALL chkxer('STREXC',INFot,NOUt,LERr,OK)
      INFot = 7
      ifst = 2
      CALL STREXC('V',1,a,1,b,1,ifst,ilst,work,info)
      CALL chkxer('STREXC',INFot,NOUt,LERr,OK)
      INFot = 8
      ifst = 1
      ilst = 0
      CALL STREXC('V',1,a,1,b,1,ifst,ilst,work,info)
      CALL chkxer('STREXC',INFot,NOUt,LERr,OK)
      INFot = 8
      ilst = 2
      CALL STREXC('V',1,a,1,b,1,ifst,ilst,work,info)
      CALL chkxer('STREXC',INFot,NOUt,LERr,OK)
      nt = nt + 8
!
!     Test STRSNA
!
      SRNamt = 'STRSNA'
      INFot = 1
      CALL STRSNA('X','A',sel,0,a,1,b,1,c,1,s,sep,1,m,work,1,iwork,info)
      CALL chkxer('STRSNA',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL STRSNA('B','X',sel,0,a,1,b,1,c,1,s,sep,1,m,work,1,iwork,info)
      CALL chkxer('STRSNA',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL STRSNA('B','A',sel,-1,a,1,b,1,c,1,s,sep,1,m,work,1,iwork,    &
     &            info)
      CALL chkxer('STRSNA',INFot,NOUt,LERr,OK)
      INFot = 6
      CALL STRSNA('V','A',sel,2,a,1,b,1,c,1,s,sep,2,m,work,2,iwork,info)
      CALL chkxer('STRSNA',INFot,NOUt,LERr,OK)
      INFot = 8
      CALL STRSNA('B','A',sel,2,a,2,b,1,c,2,s,sep,2,m,work,2,iwork,info)
      CALL chkxer('STRSNA',INFot,NOUt,LERr,OK)
      INFot = 10
      CALL STRSNA('B','A',sel,2,a,2,b,2,c,1,s,sep,2,m,work,2,iwork,info)
      CALL chkxer('STRSNA',INFot,NOUt,LERr,OK)
      INFot = 13
      CALL STRSNA('B','A',sel,1,a,1,b,1,c,1,s,sep,0,m,work,1,iwork,info)
      CALL chkxer('STRSNA',INFot,NOUt,LERr,OK)
      INFot = 13
      CALL STRSNA('B','S',sel,2,a,2,b,2,c,2,s,sep,1,m,work,2,iwork,info)
      CALL chkxer('STRSNA',INFot,NOUt,LERr,OK)
      INFot = 16
      CALL STRSNA('B','A',sel,2,a,2,b,2,c,2,s,sep,2,m,work,1,iwork,info)
      CALL chkxer('STRSNA',INFot,NOUt,LERr,OK)
      nt = nt + 9
!
!     Test STRSEN
!
      sel(1) = .FALSE.
      SRNamt = 'STRSEN'
      INFot = 1
      CALL STRSEN('X','N',sel,0,a,1,b,1,wr,wi,m,s(1),sep(1),work,1,     &
     &            iwork,1,info)
      CALL chkxer('STRSEN',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL STRSEN('N','X',sel,0,a,1,b,1,wr,wi,m,s(1),sep(1),work,1,     &
     &            iwork,1,info)
      CALL chkxer('STRSEN',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL STRSEN('N','N',sel,-1,a,1,b,1,wr,wi,m,s(1),sep(1),work,1,    &
     &            iwork,1,info)
      CALL chkxer('STRSEN',INFot,NOUt,LERr,OK)
      INFot = 6
      CALL STRSEN('N','N',sel,2,a,1,b,1,wr,wi,m,s(1),sep(1),work,2,     &
     &            iwork,1,info)
      CALL chkxer('STRSEN',INFot,NOUt,LERr,OK)
      INFot = 8
      CALL STRSEN('N','V',sel,2,a,2,b,1,wr,wi,m,s(1),sep(1),work,1,     &
     &            iwork,1,info)
      CALL chkxer('STRSEN',INFot,NOUt,LERr,OK)
      INFot = 15
      CALL STRSEN('N','V',sel,2,a,2,b,2,wr,wi,m,s(1),sep(1),work,0,     &
     &            iwork,1,info)
      CALL chkxer('STRSEN',INFot,NOUt,LERr,OK)
      INFot = 15
      CALL STRSEN('E','V',sel,3,a,3,b,3,wr,wi,m,s(1),sep(1),work,1,     &
     &            iwork,1,info)
      CALL chkxer('STRSEN',INFot,NOUt,LERr,OK)
      INFot = 15
      CALL STRSEN('V','V',sel,3,a,3,b,3,wr,wi,m,s(1),sep(1),work,3,     &
     &            iwork,2,info)
      CALL chkxer('STRSEN',INFot,NOUt,LERr,OK)
      INFot = 17
      CALL STRSEN('E','V',sel,2,a,2,b,2,wr,wi,m,s(1),sep(1),work,1,     &
     &            iwork,0,info)
      CALL chkxer('STRSEN',INFot,NOUt,LERr,OK)
      INFot = 17
      CALL STRSEN('V','V',sel,3,a,3,b,3,wr,wi,m,s(1),sep(1),work,4,     &
     &            iwork,1,info)
      CALL chkxer('STRSEN',INFot,NOUt,LERr,OK)
      nt = nt + 10
!
!     Print a summary line.
!
      IF ( OK ) THEN
         WRITE (NOUt,FMT=99001) Path , nt
      ELSE
         WRITE (NOUt,FMT=99002) Path
      ENDIF
!
      RETURN
99001 FORMAT (1X,A3,' routines passed the tests of the error exits (',  &
     &        I3,' tests done)')
99002 FORMAT (' *** ',A3,' routines failed the tests of the error ex',  &
     &        'its ***')
!
!     End of serrec
!
      END SUBROUTINE serrec
