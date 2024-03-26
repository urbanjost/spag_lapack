!*==derrec.f90  processed by SPAG 7.51RB at 17:37 on  4 Mar 2022
!> \brief \b derrec
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE derrec( PATH, NUNIT )
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
!> derrec tests the error exits for the routines for eigen- condition
!> estimation for DOUBLE PRECISION matrices:
!>    DTRSYL, DTREXC, DTRSNA and DTRSEN.
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
      SUBROUTINE derrec(Path,Nunit)
      IMPLICIT NONE
!*--DERREC60
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
      DOUBLE PRECISION ONE , ZERO
      PARAMETER (NMAX=4,ONE=1.0D0,ZERO=0.0D0)
!     ..
!     .. Local Scalars ..
      INTEGER i , ifst , ilst , info , j , m , nt
      DOUBLE PRECISION scale
!     ..
!     .. Local Arrays ..
      LOGICAL sel(NMAX)
      INTEGER iwork(NMAX)
      DOUBLE PRECISION a(NMAX,NMAX) , b(NMAX,NMAX) , c(NMAX,NMAX) ,     &
     &                 s(NMAX) , sep(NMAX) , wi(NMAX) , work(NMAX) ,    &
     &                 wr(NMAX)
!     ..
!     .. External Subroutines ..
      EXTERNAL DTREXC , DTRSEN , DTRSNA , DTRSYL
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
!     Test DTRSYL
!
      SRNamt = 'DTRSYL'
      INFot = 1
      CALL DTRSYL('X','N',1,0,0,a,1,b,1,c,1,scale,info)
      CALL chkxer('DTRSYL',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL DTRSYL('N','X',1,0,0,a,1,b,1,c,1,scale,info)
      CALL chkxer('DTRSYL',INFot,NOUt,LERr,OK)
      INFot = 3
      CALL DTRSYL('N','N',0,0,0,a,1,b,1,c,1,scale,info)
      CALL chkxer('DTRSYL',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL DTRSYL('N','N',1,-1,0,a,1,b,1,c,1,scale,info)
      CALL chkxer('DTRSYL',INFot,NOUt,LERr,OK)
      INFot = 5
      CALL DTRSYL('N','N',1,0,-1,a,1,b,1,c,1,scale,info)
      CALL chkxer('DTRSYL',INFot,NOUt,LERr,OK)
      INFot = 7
      CALL DTRSYL('N','N',1,2,0,a,1,b,1,c,2,scale,info)
      CALL chkxer('DTRSYL',INFot,NOUt,LERr,OK)
      INFot = 9
      CALL DTRSYL('N','N',1,0,2,a,1,b,1,c,1,scale,info)
      CALL chkxer('DTRSYL',INFot,NOUt,LERr,OK)
      INFot = 11
      CALL DTRSYL('N','N',1,2,0,a,2,b,1,c,1,scale,info)
      CALL chkxer('DTRSYL',INFot,NOUt,LERr,OK)
      nt = nt + 8
!
!     Test DTREXC
!
      SRNamt = 'DTREXC'
      ifst = 1
      ilst = 1
      INFot = 1
      CALL DTREXC('X',1,a,1,b,1,ifst,ilst,work,info)
      CALL chkxer('DTREXC',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL DTREXC('N',-1,a,1,b,1,ifst,ilst,work,info)
      CALL chkxer('DTREXC',INFot,NOUt,LERr,OK)
      INFot = 4
      ilst = 2
      CALL DTREXC('N',2,a,1,b,1,ifst,ilst,work,info)
      CALL chkxer('DTREXC',INFot,NOUt,LERr,OK)
      INFot = 6
      CALL DTREXC('V',2,a,2,b,1,ifst,ilst,work,info)
      CALL chkxer('DTREXC',INFot,NOUt,LERr,OK)
      INFot = 7
      ifst = 0
      ilst = 1
      CALL DTREXC('V',1,a,1,b,1,ifst,ilst,work,info)
      CALL chkxer('DTREXC',INFot,NOUt,LERr,OK)
      INFot = 7
      ifst = 2
      CALL DTREXC('V',1,a,1,b,1,ifst,ilst,work,info)
      CALL chkxer('DTREXC',INFot,NOUt,LERr,OK)
      INFot = 8
      ifst = 1
      ilst = 0
      CALL DTREXC('V',1,a,1,b,1,ifst,ilst,work,info)
      CALL chkxer('DTREXC',INFot,NOUt,LERr,OK)
      INFot = 8
      ilst = 2
      CALL DTREXC('V',1,a,1,b,1,ifst,ilst,work,info)
      CALL chkxer('DTREXC',INFot,NOUt,LERr,OK)
      nt = nt + 8
!
!     Test DTRSNA
!
      SRNamt = 'DTRSNA'
      INFot = 1
      CALL DTRSNA('X','A',sel,0,a,1,b,1,c,1,s,sep,1,m,work,1,iwork,info)
      CALL chkxer('DTRSNA',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL DTRSNA('B','X',sel,0,a,1,b,1,c,1,s,sep,1,m,work,1,iwork,info)
      CALL chkxer('DTRSNA',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL DTRSNA('B','A',sel,-1,a,1,b,1,c,1,s,sep,1,m,work,1,iwork,    &
     &            info)
      CALL chkxer('DTRSNA',INFot,NOUt,LERr,OK)
      INFot = 6
      CALL DTRSNA('V','A',sel,2,a,1,b,1,c,1,s,sep,2,m,work,2,iwork,info)
      CALL chkxer('DTRSNA',INFot,NOUt,LERr,OK)
      INFot = 8
      CALL DTRSNA('B','A',sel,2,a,2,b,1,c,2,s,sep,2,m,work,2,iwork,info)
      CALL chkxer('DTRSNA',INFot,NOUt,LERr,OK)
      INFot = 10
      CALL DTRSNA('B','A',sel,2,a,2,b,2,c,1,s,sep,2,m,work,2,iwork,info)
      CALL chkxer('DTRSNA',INFot,NOUt,LERr,OK)
      INFot = 13
      CALL DTRSNA('B','A',sel,1,a,1,b,1,c,1,s,sep,0,m,work,1,iwork,info)
      CALL chkxer('DTRSNA',INFot,NOUt,LERr,OK)
      INFot = 13
      CALL DTRSNA('B','S',sel,2,a,2,b,2,c,2,s,sep,1,m,work,2,iwork,info)
      CALL chkxer('DTRSNA',INFot,NOUt,LERr,OK)
      INFot = 16
      CALL DTRSNA('B','A',sel,2,a,2,b,2,c,2,s,sep,2,m,work,1,iwork,info)
      CALL chkxer('DTRSNA',INFot,NOUt,LERr,OK)
      nt = nt + 9
!
!     Test DTRSEN
!
      sel(1) = .FALSE.
      SRNamt = 'DTRSEN'
      INFot = 1
      CALL DTRSEN('X','N',sel,0,a,1,b,1,wr,wi,m,s(1),sep(1),work,1,     &
     &            iwork,1,info)
      CALL chkxer('DTRSEN',INFot,NOUt,LERr,OK)
      INFot = 2
      CALL DTRSEN('N','X',sel,0,a,1,b,1,wr,wi,m,s(1),sep(1),work,1,     &
     &            iwork,1,info)
      CALL chkxer('DTRSEN',INFot,NOUt,LERr,OK)
      INFot = 4
      CALL DTRSEN('N','N',sel,-1,a,1,b,1,wr,wi,m,s(1),sep(1),work,1,    &
     &            iwork,1,info)
      CALL chkxer('DTRSEN',INFot,NOUt,LERr,OK)
      INFot = 6
      CALL DTRSEN('N','N',sel,2,a,1,b,1,wr,wi,m,s(1),sep(1),work,2,     &
     &            iwork,1,info)
      CALL chkxer('DTRSEN',INFot,NOUt,LERr,OK)
      INFot = 8
      CALL DTRSEN('N','V',sel,2,a,2,b,1,wr,wi,m,s(1),sep(1),work,1,     &
     &            iwork,1,info)
      CALL chkxer('DTRSEN',INFot,NOUt,LERr,OK)
      INFot = 15
      CALL DTRSEN('N','V',sel,2,a,2,b,2,wr,wi,m,s(1),sep(1),work,0,     &
     &            iwork,1,info)
      CALL chkxer('DTRSEN',INFot,NOUt,LERr,OK)
      INFot = 15
      CALL DTRSEN('E','V',sel,3,a,3,b,3,wr,wi,m,s(1),sep(1),work,1,     &
     &            iwork,1,info)
      CALL chkxer('DTRSEN',INFot,NOUt,LERr,OK)
      INFot = 15
      CALL DTRSEN('V','V',sel,3,a,3,b,3,wr,wi,m,s(1),sep(1),work,3,     &
     &            iwork,2,info)
      CALL chkxer('DTRSEN',INFot,NOUt,LERr,OK)
      INFot = 17
      CALL DTRSEN('E','V',sel,2,a,2,b,2,wr,wi,m,s(1),sep(1),work,1,     &
     &            iwork,0,info)
      CALL chkxer('DTRSEN',INFot,NOUt,LERr,OK)
      INFot = 17
      CALL DTRSEN('V','V',sel,3,a,3,b,3,wr,wi,m,s(1),sep(1),work,4,     &
     &            iwork,1,info)
      CALL chkxer('DTRSEN',INFot,NOUt,LERr,OK)
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
!     End of derrec
!
      END SUBROUTINE derrec
