!*==derrgg.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b DERRGG
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DERRGG( PATH, NUNIT )
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
!> DERRGG tests the error exits for DGGES, DGGESX, DGGEV,  DGGEVX,
!> DGGGLM, DGGHRD, DGGLSE, DGGQRF, DGGRQF, DGGSVD3,
!> DGGSVP3, DHGEQZ, DORCSD, DTGEVC, DTGEXC, DTGSEN, DTGSJA, DTGSNA,
!> DGGES3, DGGEV3, and DTGSYL.
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
      SUBROUTINE DERRGG(Path,Nunit)
      IMPLICIT NONE
!*--DERRGG61
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
      PARAMETER (NMAX=3,LW=6*NMAX)
      DOUBLE PRECISION ONE , ZERO
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
!     ..
!     .. Local Scalars ..
      CHARACTER*2 c2
      INTEGER dummyk , dummyl , i , ifst , ilo , ihi , ilst , info , j ,&
     &        m , ncycle , nt , sdim , lwork
      DOUBLE PRECISION anrm , bnrm , dif , scale , tola , tolb
!     ..
!     .. Local Arrays ..
      LOGICAL bw(NMAX) , sel(NMAX)
      INTEGER iw(NMAX) , idum(NMAX)
      DOUBLE PRECISION a(NMAX,NMAX) , b(NMAX,NMAX) , ls(NMAX) ,         &
     &                 q(NMAX,NMAX) , r1(NMAX) , r2(NMAX) , r3(NMAX) ,  &
     &                 rce(2) , rcv(2) , rs(NMAX) , tau(NMAX) ,         &
     &                 u(NMAX,NMAX) , v(NMAX,NMAX) , w(LW) ,            &
     &                 z(NMAX,NMAX)
!     ..
!     .. External Functions ..
      LOGICAL DLCTES , DLCTSX , LSAMEN
      EXTERNAL DLCTES , DLCTSX , LSAMEN
!     ..
!     .. External Subroutines ..
      EXTERNAL CHKXER , DGGES , DGGESX , DGGEV , DGGEVX , DGGGLM ,      &
     &         DGGHRD , DGGLSE , DGGQRF , DGGRQF , DHGEQZ , DORCSD ,    &
     &         DTGEVC , DTGEXC , DTGSEN , DTGSJA , DTGSNA , DTGSYL ,    &
     &         DGGHD3 , DGGES3 , DGGEV3 , DGGSVD3 , DGGSVP3
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
         sel(j) = .TRUE.
         DO i = 1 , NMAX
            a(i,j) = ZERO
            b(i,j) = ZERO
         ENDDO
      ENDDO
      DO i = 1 , NMAX
         a(i,i) = ONE
         b(i,i) = ONE
      ENDDO
      OK = .TRUE.
      tola = 1.0D0
      tolb = 1.0D0
      ifst = 1
      ilst = 1
      nt = 0
      lwork = 1
!
!     Test error exits for the GG path.
!
      IF ( LSAMEN(2,c2,'GG') ) THEN
!
!        DGGHRD
!
         SRNamt = 'DGGHRD'
         INFot = 1
         CALL DGGHRD('/','N',0,1,0,a,1,b,1,q,1,z,1,info)
         CALL CHKXER('DGGHRD',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DGGHRD('N','/',0,1,0,a,1,b,1,q,1,z,1,info)
         CALL CHKXER('DGGHRD',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DGGHRD('N','N',-1,0,0,a,1,b,1,q,1,z,1,info)
         CALL CHKXER('DGGHRD',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL DGGHRD('N','N',0,0,0,a,1,b,1,q,1,z,1,info)
         CALL CHKXER('DGGHRD',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL DGGHRD('N','N',0,1,1,a,1,b,1,q,1,z,1,info)
         CALL CHKXER('DGGHRD',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL DGGHRD('N','N',2,1,1,a,1,b,2,q,1,z,1,info)
         CALL CHKXER('DGGHRD',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL DGGHRD('N','N',2,1,1,a,2,b,1,q,1,z,1,info)
         CALL CHKXER('DGGHRD',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL DGGHRD('V','N',2,1,1,a,2,b,2,q,1,z,1,info)
         CALL CHKXER('DGGHRD',INFot,NOUt,LERr,OK)
         INFot = 13
         CALL DGGHRD('N','V',2,1,1,a,2,b,2,q,1,z,1,info)
         CALL CHKXER('DGGHRD',INFot,NOUt,LERr,OK)
         nt = nt + 9
!
!        DGGHD3
!
         SRNamt = 'DGGHD3'
         INFot = 1
         CALL DGGHD3('/','N',0,1,0,a,1,b,1,q,1,z,1,w,LW,info)
         CALL CHKXER('DGGHD3',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DGGHD3('N','/',0,1,0,a,1,b,1,q,1,z,1,w,LW,info)
         CALL CHKXER('DGGHD3',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DGGHD3('N','N',-1,0,0,a,1,b,1,q,1,z,1,w,LW,info)
         CALL CHKXER('DGGHD3',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL DGGHD3('N','N',0,0,0,a,1,b,1,q,1,z,1,w,LW,info)
         CALL CHKXER('DGGHD3',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL DGGHD3('N','N',0,1,1,a,1,b,1,q,1,z,1,w,LW,info)
         CALL CHKXER('DGGHD3',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL DGGHD3('N','N',2,1,1,a,1,b,2,q,1,z,1,w,LW,info)
         CALL CHKXER('DGGHD3',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL DGGHD3('N','N',2,1,1,a,2,b,1,q,1,z,1,w,LW,info)
         CALL CHKXER('DGGHD3',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL DGGHD3('V','N',2,1,1,a,2,b,2,q,1,z,1,w,LW,info)
         CALL CHKXER('DGGHD3',INFot,NOUt,LERr,OK)
         INFot = 13
         CALL DGGHD3('N','V',2,1,1,a,2,b,2,q,1,z,1,w,LW,info)
         CALL CHKXER('DGGHD3',INFot,NOUt,LERr,OK)
         nt = nt + 9
!
!        DHGEQZ
!
         SRNamt = 'DHGEQZ'
         INFot = 1
         CALL DHGEQZ('/','N','N',0,1,0,a,1,b,1,r1,r2,r3,q,1,z,1,w,LW,   &
     &               info)
         CALL CHKXER('DHGEQZ',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DHGEQZ('E','/','N',0,1,0,a,1,b,1,r1,r2,r3,q,1,z,1,w,LW,   &
     &               info)
         CALL CHKXER('DHGEQZ',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DHGEQZ('E','N','/',0,1,0,a,1,b,1,r1,r2,r3,q,1,z,1,w,LW,   &
     &               info)
         CALL CHKXER('DHGEQZ',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL DHGEQZ('E','N','N',-1,0,0,a,1,b,1,r1,r2,r3,q,1,z,1,w,LW,  &
     &               info)
         CALL CHKXER('DHGEQZ',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL DHGEQZ('E','N','N',0,0,0,a,1,b,1,r1,r2,r3,q,1,z,1,w,LW,   &
     &               info)
         CALL CHKXER('DHGEQZ',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL DHGEQZ('E','N','N',0,1,1,a,1,b,1,r1,r2,r3,q,1,z,1,w,LW,   &
     &               info)
         CALL CHKXER('DHGEQZ',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL DHGEQZ('E','N','N',2,1,1,a,1,b,2,r1,r2,r3,q,1,z,1,w,LW,   &
     &               info)
         CALL CHKXER('DHGEQZ',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL DHGEQZ('E','N','N',2,1,1,a,2,b,1,r1,r2,r3,q,1,z,1,w,LW,   &
     &               info)
         CALL CHKXER('DHGEQZ',INFot,NOUt,LERr,OK)
         INFot = 15
         CALL DHGEQZ('E','V','N',2,1,1,a,2,b,2,r1,r2,r3,q,1,z,1,w,LW,   &
     &               info)
         CALL CHKXER('DHGEQZ',INFot,NOUt,LERr,OK)
         INFot = 17
         CALL DHGEQZ('E','N','V',2,1,1,a,2,b,2,r1,r2,r3,q,1,z,1,w,LW,   &
     &               info)
         CALL CHKXER('DHGEQZ',INFot,NOUt,LERr,OK)
         nt = nt + 10
!
!        DTGEVC
!
         SRNamt = 'DTGEVC'
         INFot = 1
         CALL DTGEVC('/','A',sel,0,a,1,b,1,q,1,z,1,0,m,w,info)
         CALL CHKXER('DTGEVC',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DTGEVC('R','/',sel,0,a,1,b,1,q,1,z,1,0,m,w,info)
         CALL CHKXER('DTGEVC',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL DTGEVC('R','A',sel,-1,a,1,b,1,q,1,z,1,0,m,w,info)
         CALL CHKXER('DTGEVC',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL DTGEVC('R','A',sel,2,a,1,b,2,q,1,z,2,0,m,w,info)
         CALL CHKXER('DTGEVC',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL DTGEVC('R','A',sel,2,a,2,b,1,q,1,z,2,0,m,w,info)
         CALL CHKXER('DTGEVC',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL DTGEVC('L','A',sel,2,a,2,b,2,q,1,z,1,0,m,w,info)
         CALL CHKXER('DTGEVC',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL DTGEVC('R','A',sel,2,a,2,b,2,q,1,z,1,0,m,w,info)
         CALL CHKXER('DTGEVC',INFot,NOUt,LERr,OK)
         INFot = 13
         CALL DTGEVC('R','A',sel,2,a,2,b,2,q,1,z,2,1,m,w,info)
         CALL CHKXER('DTGEVC',INFot,NOUt,LERr,OK)
         nt = nt + 8
!
!     Test error exits for the GSV path.
!
      ELSEIF ( LSAMEN(3,Path,'GSV') ) THEN
!
!        DGGSVD3
!
         SRNamt = 'DGGSVD3'
         INFot = 1
         CALL DGGSVD3('/','N','N',0,0,0,dummyk,dummyl,a,1,b,1,r1,r2,u,1,&
     &                v,1,q,1,w,lwork,idum,info)
         CALL CHKXER('DGGSVD3',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DGGSVD3('N','/','N',0,0,0,dummyk,dummyl,a,1,b,1,r1,r2,u,1,&
     &                v,1,q,1,w,lwork,idum,info)
         CALL CHKXER('DGGSVD3',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DGGSVD3('N','N','/',0,0,0,dummyk,dummyl,a,1,b,1,r1,r2,u,1,&
     &                v,1,q,1,w,lwork,idum,info)
         CALL CHKXER('DGGSVD3',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL DGGSVD3('N','N','N',-1,0,0,dummyk,dummyl,a,1,b,1,r1,r2,u, &
     &                1,v,1,q,1,w,lwork,idum,info)
         CALL CHKXER('DGGSVD3',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL DGGSVD3('N','N','N',0,-1,0,dummyk,dummyl,a,1,b,1,r1,r2,u, &
     &                1,v,1,q,1,w,lwork,idum,info)
         CALL CHKXER('DGGSVD3',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL DGGSVD3('N','N','N',0,0,-1,dummyk,dummyl,a,1,b,1,r1,r2,u, &
     &                1,v,1,q,1,w,lwork,idum,info)
         CALL CHKXER('DGGSVD3',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL DGGSVD3('N','N','N',2,1,1,dummyk,dummyl,a,1,b,1,r1,r2,u,1,&
     &                v,1,q,1,w,lwork,idum,info)
         CALL CHKXER('DGGSVD3',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL DGGSVD3('N','N','N',1,1,2,dummyk,dummyl,a,1,b,1,r1,r2,u,1,&
     &                v,1,q,1,w,lwork,idum,info)
         CALL CHKXER('DGGSVD3',INFot,NOUt,LERr,OK)
         INFot = 16
         CALL DGGSVD3('U','N','N',2,2,2,dummyk,dummyl,a,2,b,2,r1,r2,u,1,&
     &                v,1,q,1,w,lwork,idum,info)
         CALL CHKXER('DGGSVD3',INFot,NOUt,LERr,OK)
         INFot = 18
         CALL DGGSVD3('N','V','N',1,1,2,dummyk,dummyl,a,1,b,2,r1,r2,u,1,&
     &                v,1,q,1,w,lwork,idum,info)
         CALL CHKXER('DGGSVD3',INFot,NOUt,LERr,OK)
         INFot = 20
         CALL DGGSVD3('N','N','Q',1,2,1,dummyk,dummyl,a,1,b,1,r1,r2,u,1,&
     &                v,1,q,1,w,lwork,idum,info)
         CALL CHKXER('DGGSVD3',INFot,NOUt,LERr,OK)
         nt = nt + 11
!
!        DGGSVP3
!
         SRNamt = 'DGGSVP3'
         INFot = 1
         CALL DGGSVP3('/','N','N',0,0,0,a,1,b,1,tola,tolb,dummyk,dummyl,&
     &                u,1,v,1,q,1,iw,tau,w,lwork,info)
         CALL CHKXER('DGGSVP3',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DGGSVP3('N','/','N',0,0,0,a,1,b,1,tola,tolb,dummyk,dummyl,&
     &                u,1,v,1,q,1,iw,tau,w,lwork,info)
         CALL CHKXER('DGGSVP3',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DGGSVP3('N','N','/',0,0,0,a,1,b,1,tola,tolb,dummyk,dummyl,&
     &                u,1,v,1,q,1,iw,tau,w,lwork,info)
         CALL CHKXER('DGGSVP3',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL DGGSVP3('N','N','N',-1,0,0,a,1,b,1,tola,tolb,dummyk,      &
     &                dummyl,u,1,v,1,q,1,iw,tau,w,lwork,info)
         CALL CHKXER('DGGSVP3',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL DGGSVP3('N','N','N',0,-1,0,a,1,b,1,tola,tolb,dummyk,      &
     &                dummyl,u,1,v,1,q,1,iw,tau,w,lwork,info)
         CALL CHKXER('DGGSVP3',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL DGGSVP3('N','N','N',0,0,-1,a,1,b,1,tola,tolb,dummyk,      &
     &                dummyl,u,1,v,1,q,1,iw,tau,w,lwork,info)
         CALL CHKXER('DGGSVP3',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL DGGSVP3('N','N','N',2,1,1,a,1,b,1,tola,tolb,dummyk,dummyl,&
     &                u,1,v,1,q,1,iw,tau,w,lwork,info)
         CALL CHKXER('DGGSVP3',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL DGGSVP3('N','N','N',1,2,1,a,1,b,1,tola,tolb,dummyk,dummyl,&
     &                u,1,v,1,q,1,iw,tau,w,lwork,info)
         CALL CHKXER('DGGSVP3',INFot,NOUt,LERr,OK)
         INFot = 16
         CALL DGGSVP3('U','N','N',2,2,2,a,2,b,2,tola,tolb,dummyk,dummyl,&
     &                u,1,v,1,q,1,iw,tau,w,lwork,info)
         CALL CHKXER('DGGSVP3',INFot,NOUt,LERr,OK)
         INFot = 18
         CALL DGGSVP3('N','V','N',1,2,1,a,1,b,2,tola,tolb,dummyk,dummyl,&
     &                u,1,v,1,q,1,iw,tau,w,lwork,info)
         CALL CHKXER('DGGSVP3',INFot,NOUt,LERr,OK)
         INFot = 20
         CALL DGGSVP3('N','N','Q',1,1,2,a,1,b,1,tola,tolb,dummyk,dummyl,&
     &                u,1,v,1,q,1,iw,tau,w,lwork,info)
         CALL CHKXER('DGGSVP3',INFot,NOUt,LERr,OK)
         nt = nt + 11
!
!        DTGSJA
!
         SRNamt = 'DTGSJA'
         INFot = 1
         CALL DTGSJA('/','N','N',0,0,0,dummyk,dummyl,a,1,b,1,tola,tolb, &
     &               r1,r2,u,1,v,1,q,1,w,ncycle,info)
         CALL CHKXER('DTGSJA',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DTGSJA('N','/','N',0,0,0,dummyk,dummyl,a,1,b,1,tola,tolb, &
     &               r1,r2,u,1,v,1,q,1,w,ncycle,info)
         CALL CHKXER('DTGSJA',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DTGSJA('N','N','/',0,0,0,dummyk,dummyl,a,1,b,1,tola,tolb, &
     &               r1,r2,u,1,v,1,q,1,w,ncycle,info)
         CALL CHKXER('DTGSJA',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL DTGSJA('N','N','N',-1,0,0,dummyk,dummyl,a,1,b,1,tola,tolb,&
     &               r1,r2,u,1,v,1,q,1,w,ncycle,info)
         CALL CHKXER('DTGSJA',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL DTGSJA('N','N','N',0,-1,0,dummyk,dummyl,a,1,b,1,tola,tolb,&
     &               r1,r2,u,1,v,1,q,1,w,ncycle,info)
         CALL CHKXER('DTGSJA',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL DTGSJA('N','N','N',0,0,-1,dummyk,dummyl,a,1,b,1,tola,tolb,&
     &               r1,r2,u,1,v,1,q,1,w,ncycle,info)
         CALL CHKXER('DTGSJA',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL DTGSJA('N','N','N',0,0,0,dummyk,dummyl,a,0,b,1,tola,tolb, &
     &               r1,r2,u,1,v,1,q,1,w,ncycle,info)
         CALL CHKXER('DTGSJA',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL DTGSJA('N','N','N',0,0,0,dummyk,dummyl,a,1,b,0,tola,tolb, &
     &               r1,r2,u,1,v,1,q,1,w,ncycle,info)
         CALL CHKXER('DTGSJA',INFot,NOUt,LERr,OK)
         INFot = 18
         CALL DTGSJA('U','N','N',0,0,0,dummyk,dummyl,a,1,b,1,tola,tolb, &
     &               r1,r2,u,0,v,1,q,1,w,ncycle,info)
         CALL CHKXER('DTGSJA',INFot,NOUt,LERr,OK)
         INFot = 20
         CALL DTGSJA('N','V','N',0,0,0,dummyk,dummyl,a,1,b,1,tola,tolb, &
     &               r1,r2,u,1,v,0,q,1,w,ncycle,info)
         CALL CHKXER('DTGSJA',INFot,NOUt,LERr,OK)
         INFot = 22
         CALL DTGSJA('N','N','Q',0,0,0,dummyk,dummyl,a,1,b,1,tola,tolb, &
     &               r1,r2,u,1,v,1,q,0,w,ncycle,info)
         CALL CHKXER('DTGSJA',INFot,NOUt,LERr,OK)
         nt = nt + 11
!
!     Test error exits for the GLM path.
!
      ELSEIF ( LSAMEN(3,Path,'GLM') ) THEN
!
!        DGGGLM
!
         SRNamt = 'DGGGLM'
         INFot = 1
         CALL DGGGLM(-1,0,0,a,1,b,1,r1,r2,r3,w,LW,info)
         CALL CHKXER('DGGGLM',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DGGGLM(0,-1,0,a,1,b,1,r1,r2,r3,w,LW,info)
         CALL CHKXER('DGGGLM',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DGGGLM(0,1,0,a,1,b,1,r1,r2,r3,w,LW,info)
         CALL CHKXER('DGGGLM',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DGGGLM(0,0,-1,a,1,b,1,r1,r2,r3,w,LW,info)
         CALL CHKXER('DGGGLM',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DGGGLM(1,0,0,a,1,b,1,r1,r2,r3,w,LW,info)
         CALL CHKXER('DGGGLM',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL DGGGLM(0,0,0,a,0,b,1,r1,r2,r3,w,LW,info)
         CALL CHKXER('DGGGLM',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL DGGGLM(0,0,0,a,1,b,0,r1,r2,r3,w,LW,info)
         CALL CHKXER('DGGGLM',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL DGGGLM(1,1,1,a,1,b,1,r1,r2,r3,w,1,info)
         CALL CHKXER('DGGGLM',INFot,NOUt,LERr,OK)
         nt = nt + 8
!
!     Test error exits for the LSE path.
!
      ELSEIF ( LSAMEN(3,Path,'LSE') ) THEN
!
!        DGGLSE
!
         SRNamt = 'DGGLSE'
         INFot = 1
         CALL DGGLSE(-1,0,0,a,1,b,1,r1,r2,r3,w,LW,info)
         CALL CHKXER('DGGLSE',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DGGLSE(0,-1,0,a,1,b,1,r1,r2,r3,w,LW,info)
         CALL CHKXER('DGGLSE',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DGGLSE(0,0,-1,a,1,b,1,r1,r2,r3,w,LW,info)
         CALL CHKXER('DGGLSE',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DGGLSE(0,0,1,a,1,b,1,r1,r2,r3,w,LW,info)
         CALL CHKXER('DGGLSE',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DGGLSE(0,1,0,a,1,b,1,r1,r2,r3,w,LW,info)
         CALL CHKXER('DGGLSE',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL DGGLSE(0,0,0,a,0,b,1,r1,r2,r3,w,LW,info)
         CALL CHKXER('DGGLSE',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL DGGLSE(0,0,0,a,1,b,0,r1,r2,r3,w,LW,info)
         CALL CHKXER('DGGLSE',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL DGGLSE(1,1,1,a,1,b,1,r1,r2,r3,w,1,info)
         CALL CHKXER('DGGLSE',INFot,NOUt,LERr,OK)
         nt = nt + 8
!
!     Test error exits for the CSD path.
!
      ELSEIF ( LSAMEN(3,Path,'CSD') ) THEN
!
!        DORCSD
!
         SRNamt = 'DORCSD'
         INFot = 7
         CALL DORCSD('Y','Y','Y','Y','N','N',-1,0,0,a,1,a,1,a,1,a,1,a,a,&
     &               1,a,1,a,1,a,1,w,LW,iw,info)
         CALL CHKXER('DORCSD',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL DORCSD('Y','Y','Y','Y','N','N',1,-1,0,a,1,a,1,a,1,a,1,a,a,&
     &               1,a,1,a,1,a,1,w,LW,iw,info)
         CALL CHKXER('DORCSD',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL DORCSD('Y','Y','Y','Y','N','N',1,1,-1,a,1,a,1,a,1,a,1,a,a,&
     &               1,a,1,a,1,a,1,w,LW,iw,info)
         CALL CHKXER('DORCSD',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL DORCSD('Y','Y','Y','Y','N','N',1,1,1,a,-1,a,1,a,1,a,1,a,a,&
     &               1,a,1,a,1,a,1,w,LW,iw,info)
         CALL CHKXER('DORCSD',INFot,NOUt,LERr,OK)
         INFot = 20
         CALL DORCSD('Y','Y','Y','Y','N','N',1,1,1,a,1,a,1,a,1,a,1,a,a, &
     &               -1,a,1,a,1,a,1,w,LW,iw,info)
         CALL CHKXER('DORCSD',INFot,NOUt,LERr,OK)
         INFot = 22
         CALL DORCSD('Y','Y','Y','Y','N','N',1,1,1,a,1,a,1,a,1,a,1,a,a, &
     &               1,a,-1,a,1,a,1,w,LW,iw,info)
         CALL CHKXER('DORCSD',INFot,NOUt,LERr,OK)
         INFot = 24
         CALL DORCSD('Y','Y','Y','Y','N','N',1,1,1,a,1,a,1,a,1,a,1,a,a, &
     &               1,a,1,a,-1,a,1,w,LW,iw,info)
         CALL CHKXER('DORCSD',INFot,NOUt,LERr,OK)
         INFot = 26
         CALL DORCSD('Y','Y','Y','Y','N','N',1,1,1,a,1,a,1,a,1,a,1,a,a, &
     &               1,a,1,a,1,a,-1,w,LW,iw,info)
         CALL CHKXER('DORCSD',INFot,NOUt,LERr,OK)
         nt = nt + 8
!
!     Test error exits for the GQR path.
!
      ELSEIF ( LSAMEN(3,Path,'GQR') ) THEN
!
!        DGGQRF
!
         SRNamt = 'DGGQRF'
         INFot = 1
         CALL DGGQRF(-1,0,0,a,1,r1,b,1,r2,w,LW,info)
         CALL CHKXER('DGGQRF',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DGGQRF(0,-1,0,a,1,r1,b,1,r2,w,LW,info)
         CALL CHKXER('DGGQRF',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DGGQRF(0,0,-1,a,1,r1,b,1,r2,w,LW,info)
         CALL CHKXER('DGGQRF',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL DGGQRF(0,0,0,a,0,r1,b,1,r2,w,LW,info)
         CALL CHKXER('DGGQRF',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL DGGQRF(0,0,0,a,1,r1,b,0,r2,w,LW,info)
         CALL CHKXER('DGGQRF',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL DGGQRF(1,1,2,a,1,r1,b,1,r2,w,1,info)
         CALL CHKXER('DGGQRF',INFot,NOUt,LERr,OK)
         nt = nt + 6
!
!        DGGRQF
!
         SRNamt = 'DGGRQF'
         INFot = 1
         CALL DGGRQF(-1,0,0,a,1,r1,b,1,r2,w,LW,info)
         CALL CHKXER('DGGRQF',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DGGRQF(0,-1,0,a,1,r1,b,1,r2,w,LW,info)
         CALL CHKXER('DGGRQF',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DGGRQF(0,0,-1,a,1,r1,b,1,r2,w,LW,info)
         CALL CHKXER('DGGRQF',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL DGGRQF(0,0,0,a,0,r1,b,1,r2,w,LW,info)
         CALL CHKXER('DGGRQF',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL DGGRQF(0,0,0,a,1,r1,b,0,r2,w,LW,info)
         CALL CHKXER('DGGRQF',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL DGGRQF(1,1,2,a,1,r1,b,1,r2,w,1,info)
         CALL CHKXER('DGGRQF',INFot,NOUt,LERr,OK)
         nt = nt + 6
!
!     Test error exits for the DGS, DGV, DGX, and DXV paths.
!
      ELSEIF ( LSAMEN(3,Path,'DGS') .OR. LSAMEN(3,Path,'DGV') .OR.      &
     &         LSAMEN(3,Path,'DGX') .OR. LSAMEN(3,Path,'DXV') ) THEN
!
!        DGGES
!
         SRNamt = 'DGGES '
         INFot = 1
         CALL DGGES('/','N','S',DLCTES,1,a,1,b,1,sdim,r1,r2,r3,q,1,u,1, &
     &              w,1,bw,info)
         CALL CHKXER('DGGES ',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DGGES('N','/','S',DLCTES,1,a,1,b,1,sdim,r1,r2,r3,q,1,u,1, &
     &              w,1,bw,info)
         CALL CHKXER('DGGES ',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DGGES('N','V','/',DLCTES,1,a,1,b,1,sdim,r1,r2,r3,q,1,u,1, &
     &              w,1,bw,info)
         CALL CHKXER('DGGES ',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL DGGES('N','V','S',DLCTES,-1,a,1,b,1,sdim,r1,r2,r3,q,1,u,1,&
     &              w,1,bw,info)
         CALL CHKXER('DGGES ',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL DGGES('N','V','S',DLCTES,1,a,0,b,1,sdim,r1,r2,r3,q,1,u,1, &
     &              w,1,bw,info)
         CALL CHKXER('DGGES ',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL DGGES('N','V','S',DLCTES,1,a,1,b,0,sdim,r1,r2,r3,q,1,u,1, &
     &              w,1,bw,info)
         CALL CHKXER('DGGES ',INFot,NOUt,LERr,OK)
         INFot = 15
         CALL DGGES('N','V','S',DLCTES,1,a,1,b,1,sdim,r1,r2,r3,q,0,u,1, &
     &              w,1,bw,info)
         CALL CHKXER('DGGES ',INFot,NOUt,LERr,OK)
         INFot = 15
         CALL DGGES('V','V','S',DLCTES,2,a,2,b,2,sdim,r1,r2,r3,q,1,u,2, &
     &              w,1,bw,info)
         CALL CHKXER('DGGES ',INFot,NOUt,LERr,OK)
         INFot = 17
         CALL DGGES('N','V','S',DLCTES,1,a,1,b,1,sdim,r1,r2,r3,q,1,u,0, &
     &              w,1,bw,info)
         CALL CHKXER('DGGES ',INFot,NOUt,LERr,OK)
         INFot = 17
         CALL DGGES('V','V','S',DLCTES,2,a,2,b,2,sdim,r1,r2,r3,q,2,u,1, &
     &              w,1,bw,info)
         CALL CHKXER('DGGES ',INFot,NOUt,LERr,OK)
         INFot = 19
         CALL DGGES('V','V','S',DLCTES,2,a,2,b,2,sdim,r1,r2,r3,q,2,u,2, &
     &              w,1,bw,info)
         CALL CHKXER('DGGES ',INFot,NOUt,LERr,OK)
         nt = nt + 11
!
!        DGGES3
!
         SRNamt = 'DGGES3 '
         INFot = 1
         CALL DGGES3('/','N','S',DLCTES,1,a,1,b,1,sdim,r1,r2,r3,q,1,u,1,&
     &               w,1,bw,info)
         CALL CHKXER('DGGES3 ',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DGGES3('N','/','S',DLCTES,1,a,1,b,1,sdim,r1,r2,r3,q,1,u,1,&
     &               w,1,bw,info)
         CALL CHKXER('DGGES3 ',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DGGES3('N','V','/',DLCTES,1,a,1,b,1,sdim,r1,r2,r3,q,1,u,1,&
     &               w,1,bw,info)
         CALL CHKXER('DGGES3 ',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL DGGES3('N','V','S',DLCTES,-1,a,1,b,1,sdim,r1,r2,r3,q,1,u, &
     &               1,w,1,bw,info)
         CALL CHKXER('DGGES3 ',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL DGGES3('N','V','S',DLCTES,1,a,0,b,1,sdim,r1,r2,r3,q,1,u,1,&
     &               w,1,bw,info)
         CALL CHKXER('DGGES3 ',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL DGGES3('N','V','S',DLCTES,1,a,1,b,0,sdim,r1,r2,r3,q,1,u,1,&
     &               w,1,bw,info)
         CALL CHKXER('DGGES3 ',INFot,NOUt,LERr,OK)
         INFot = 15
         CALL DGGES3('N','V','S',DLCTES,1,a,1,b,1,sdim,r1,r2,r3,q,0,u,1,&
     &               w,1,bw,info)
         CALL CHKXER('DGGES3 ',INFot,NOUt,LERr,OK)
         INFot = 15
         CALL DGGES3('V','V','S',DLCTES,2,a,2,b,2,sdim,r1,r2,r3,q,1,u,2,&
     &               w,1,bw,info)
         CALL CHKXER('DGGES3 ',INFot,NOUt,LERr,OK)
         INFot = 17
         CALL DGGES3('N','V','S',DLCTES,1,a,1,b,1,sdim,r1,r2,r3,q,1,u,0,&
     &               w,1,bw,info)
         CALL CHKXER('DGGES3 ',INFot,NOUt,LERr,OK)
         INFot = 17
         CALL DGGES3('V','V','S',DLCTES,2,a,2,b,2,sdim,r1,r2,r3,q,2,u,1,&
     &               w,1,bw,info)
         CALL CHKXER('DGGES3 ',INFot,NOUt,LERr,OK)
         INFot = 19
         CALL DGGES3('V','V','S',DLCTES,2,a,2,b,2,sdim,r1,r2,r3,q,2,u,2,&
     &               w,1,bw,info)
         CALL CHKXER('DGGES3 ',INFot,NOUt,LERr,OK)
         nt = nt + 11
!
!        DGGESX
!
         SRNamt = 'DGGESX'
         INFot = 1
         CALL DGGESX('/','N','S',DLCTSX,'N',1,a,1,b,1,sdim,r1,r2,r3,q,1,&
     &               u,1,rce,rcv,w,1,iw,1,bw,info)
         CALL CHKXER('DGGESX',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DGGESX('N','/','S',DLCTSX,'N',1,a,1,b,1,sdim,r1,r2,r3,q,1,&
     &               u,1,rce,rcv,w,1,iw,1,bw,info)
         CALL CHKXER('DGGESX',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DGGESX('V','V','/',DLCTSX,'N',1,a,1,b,1,sdim,r1,r2,r3,q,1,&
     &               u,1,rce,rcv,w,1,iw,1,bw,info)
         CALL CHKXER('DGGESX',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL DGGESX('V','V','S',DLCTSX,'/',1,a,1,b,1,sdim,r1,r2,r3,q,1,&
     &               u,1,rce,rcv,w,1,iw,1,bw,info)
         CALL CHKXER('DGGESX',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL DGGESX('V','V','S',DLCTSX,'B',-1,a,1,b,1,sdim,r1,r2,r3,q, &
     &               1,u,1,rce,rcv,w,1,iw,1,bw,info)
         CALL CHKXER('DGGESX',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL DGGESX('V','V','S',DLCTSX,'B',1,a,0,b,1,sdim,r1,r2,r3,q,1,&
     &               u,1,rce,rcv,w,1,iw,1,bw,info)
         CALL CHKXER('DGGESX',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL DGGESX('V','V','S',DLCTSX,'B',1,a,1,b,0,sdim,r1,r2,r3,q,1,&
     &               u,1,rce,rcv,w,1,iw,1,bw,info)
         CALL CHKXER('DGGESX',INFot,NOUt,LERr,OK)
         INFot = 16
         CALL DGGESX('V','V','S',DLCTSX,'B',1,a,1,b,1,sdim,r1,r2,r3,q,0,&
     &               u,1,rce,rcv,w,1,iw,1,bw,info)
         CALL CHKXER('DGGESX',INFot,NOUt,LERr,OK)
         INFot = 16
         CALL DGGESX('V','V','S',DLCTSX,'B',2,a,2,b,2,sdim,r1,r2,r3,q,1,&
     &               u,1,rce,rcv,w,1,iw,1,bw,info)
         CALL CHKXER('DGGESX',INFot,NOUt,LERr,OK)
         INFot = 18
         CALL DGGESX('V','V','S',DLCTSX,'B',1,a,1,b,1,sdim,r1,r2,r3,q,1,&
     &               u,0,rce,rcv,w,1,iw,1,bw,info)
         CALL CHKXER('DGGESX',INFot,NOUt,LERr,OK)
         INFot = 18
         CALL DGGESX('V','V','S',DLCTSX,'B',2,a,2,b,2,sdim,r1,r2,r3,q,2,&
     &               u,1,rce,rcv,w,1,iw,1,bw,info)
         CALL CHKXER('DGGESX',INFot,NOUt,LERr,OK)
         INFot = 22
         CALL DGGESX('V','V','S',DLCTSX,'B',2,a,2,b,2,sdim,r1,r2,r3,q,2,&
     &               u,2,rce,rcv,w,1,iw,1,bw,info)
         CALL CHKXER('DGGESX',INFot,NOUt,LERr,OK)
         INFot = 24
         CALL DGGESX('V','V','S',DLCTSX,'V',1,a,1,b,1,sdim,r1,r2,r3,q,1,&
     &               u,1,rce,rcv,w,32,iw,0,bw,info)
         CALL CHKXER('DGGESX',INFot,NOUt,LERr,OK)
         nt = nt + 13
!
!        DGGEV
!
         SRNamt = 'DGGEV '
         INFot = 1
         CALL DGGEV('/','N',1,a,1,b,1,r1,r2,r3,q,1,u,1,w,1,info)
         CALL CHKXER('DGGEV ',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DGGEV('N','/',1,a,1,b,1,r1,r2,r3,q,1,u,1,w,1,info)
         CALL CHKXER('DGGEV ',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DGGEV('V','V',-1,a,1,b,1,r1,r2,r3,q,1,u,1,w,1,info)
         CALL CHKXER('DGGEV ',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL DGGEV('V','V',1,a,0,b,1,r1,r2,r3,q,1,u,1,w,1,info)
         CALL CHKXER('DGGEV ',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL DGGEV('V','V',1,a,1,b,0,r1,r2,r3,q,1,u,1,w,1,info)
         CALL CHKXER('DGGEV ',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL DGGEV('N','V',1,a,1,b,1,r1,r2,r3,q,0,u,1,w,1,info)
         CALL CHKXER('DGGEV ',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL DGGEV('V','V',2,a,2,b,2,r1,r2,r3,q,1,u,2,w,1,info)
         CALL CHKXER('DGGEV ',INFot,NOUt,LERr,OK)
         INFot = 14
         CALL DGGEV('V','N',2,a,2,b,2,r1,r2,r3,q,2,u,0,w,1,info)
         CALL CHKXER('DGGEV ',INFot,NOUt,LERr,OK)
         INFot = 14
         CALL DGGEV('V','V',2,a,2,b,2,r1,r2,r3,q,2,u,1,w,1,info)
         CALL CHKXER('DGGEV ',INFot,NOUt,LERr,OK)
         INFot = 16
         CALL DGGEV('V','V',1,a,1,b,1,r1,r2,r3,q,1,u,1,w,1,info)
         CALL CHKXER('DGGEV ',INFot,NOUt,LERr,OK)
         nt = nt + 10
!
!        DGGEV3
!
         SRNamt = 'DGGEV3 '
         INFot = 1
         CALL DGGEV3('/','N',1,a,1,b,1,r1,r2,r3,q,1,u,1,w,1,info)
         CALL CHKXER('DGGEV3 ',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DGGEV3('N','/',1,a,1,b,1,r1,r2,r3,q,1,u,1,w,1,info)
         CALL CHKXER('DGGEV3 ',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DGGEV3('V','V',-1,a,1,b,1,r1,r2,r3,q,1,u,1,w,1,info)
         CALL CHKXER('DGGEV3 ',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL DGGEV3('V','V',1,a,0,b,1,r1,r2,r3,q,1,u,1,w,1,info)
         CALL CHKXER('DGGEV3 ',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL DGGEV3('V','V',1,a,1,b,0,r1,r2,r3,q,1,u,1,w,1,info)
         CALL CHKXER('DGGEV3 ',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL DGGEV3('N','V',1,a,1,b,1,r1,r2,r3,q,0,u,1,w,1,info)
         CALL CHKXER('DGGEV3 ',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL DGGEV3('V','V',2,a,2,b,2,r1,r2,r3,q,1,u,2,w,1,info)
         CALL CHKXER('DGGEV3 ',INFot,NOUt,LERr,OK)
         INFot = 14
         CALL DGGEV3('V','N',2,a,2,b,2,r1,r2,r3,q,2,u,0,w,1,info)
         CALL CHKXER('DGGEV3 ',INFot,NOUt,LERr,OK)
         INFot = 14
         CALL DGGEV3('V','V',2,a,2,b,2,r1,r2,r3,q,2,u,1,w,1,info)
         CALL CHKXER('DGGEV3 ',INFot,NOUt,LERr,OK)
         INFot = 16
         CALL DGGEV3('V','V',1,a,1,b,1,r1,r2,r3,q,1,u,1,w,1,info)
         CALL CHKXER('DGGEV3 ',INFot,NOUt,LERr,OK)
         nt = nt + 10
!
!        DGGEVX
!
         SRNamt = 'DGGEVX'
         INFot = 1
         CALL DGGEVX('/','N','N','N',1,a,1,b,1,r1,r2,r3,q,1,u,1,ilo,ihi,&
     &               ls,rs,anrm,bnrm,rce,rcv,w,1,iw,bw,info)
         CALL CHKXER('DGGEVX',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DGGEVX('N','/','N','N',1,a,1,b,1,r1,r2,r3,q,1,u,1,ilo,ihi,&
     &               ls,rs,anrm,bnrm,rce,rcv,w,1,iw,bw,info)
         CALL CHKXER('DGGEVX',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DGGEVX('N','N','/','N',1,a,1,b,1,r1,r2,r3,q,1,u,1,ilo,ihi,&
     &               ls,rs,anrm,bnrm,rce,rcv,w,1,iw,bw,info)
         CALL CHKXER('DGGEVX',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL DGGEVX('N','N','N','/',1,a,1,b,1,r1,r2,r3,q,1,u,1,ilo,ihi,&
     &               ls,rs,anrm,bnrm,rce,rcv,w,1,iw,bw,info)
         CALL CHKXER('DGGEVX',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL DGGEVX('N','N','N','N',-1,a,1,b,1,r1,r2,r3,q,1,u,1,ilo,   &
     &               ihi,ls,rs,anrm,bnrm,rce,rcv,w,1,iw,bw,info)
         CALL CHKXER('DGGEVX',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL DGGEVX('N','N','N','N',1,a,0,b,1,r1,r2,r3,q,1,u,1,ilo,ihi,&
     &               ls,rs,anrm,bnrm,rce,rcv,w,1,iw,bw,info)
         CALL CHKXER('DGGEVX',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL DGGEVX('N','N','N','N',1,a,1,b,0,r1,r2,r3,q,1,u,1,ilo,ihi,&
     &               ls,rs,anrm,bnrm,rce,rcv,w,1,iw,bw,info)
         CALL CHKXER('DGGEVX',INFot,NOUt,LERr,OK)
         INFot = 14
         CALL DGGEVX('N','N','N','N',1,a,1,b,1,r1,r2,r3,q,0,u,1,ilo,ihi,&
     &               ls,rs,anrm,bnrm,rce,rcv,w,1,iw,bw,info)
         CALL CHKXER('DGGEVX',INFot,NOUt,LERr,OK)
         INFot = 14
         CALL DGGEVX('N','V','N','N',2,a,2,b,2,r1,r2,r3,q,1,u,2,ilo,ihi,&
     &               ls,rs,anrm,bnrm,rce,rcv,w,1,iw,bw,info)
         CALL CHKXER('DGGEVX',INFot,NOUt,LERr,OK)
         INFot = 16
         CALL DGGEVX('N','N','N','N',1,a,1,b,1,r1,r2,r3,q,1,u,0,ilo,ihi,&
     &               ls,rs,anrm,bnrm,rce,rcv,w,1,iw,bw,info)
         CALL CHKXER('DGGEVX',INFot,NOUt,LERr,OK)
         INFot = 16
         CALL DGGEVX('N','N','V','N',2,a,2,b,2,r1,r2,r3,q,2,u,1,ilo,ihi,&
     &               ls,rs,anrm,bnrm,rce,rcv,w,1,iw,bw,info)
         CALL CHKXER('DGGEVX',INFot,NOUt,LERr,OK)
         INFot = 26
         CALL DGGEVX('N','N','V','N',2,a,2,b,2,r1,r2,r3,q,2,u,2,ilo,ihi,&
     &               ls,rs,anrm,bnrm,rce,rcv,w,1,iw,bw,info)
         CALL CHKXER('DGGEVX',INFot,NOUt,LERr,OK)
         nt = nt + 12
!
!        DTGEXC
!
         SRNamt = 'DTGEXC'
         INFot = 3
         CALL DTGEXC(.TRUE.,.TRUE.,-1,a,1,b,1,q,1,z,1,ifst,ilst,w,1,    &
     &               info)
         CALL CHKXER('DTGEXC',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL DTGEXC(.TRUE.,.TRUE.,1,a,0,b,1,q,1,z,1,ifst,ilst,w,1,info)
         CALL CHKXER('DTGEXC',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL DTGEXC(.TRUE.,.TRUE.,1,a,1,b,0,q,1,z,1,ifst,ilst,w,1,info)
         CALL CHKXER('DTGEXC',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL DTGEXC(.FALSE.,.TRUE.,1,a,1,b,1,q,0,z,1,ifst,ilst,w,1,    &
     &               info)
         CALL CHKXER('DTGEXC',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL DTGEXC(.TRUE.,.TRUE.,1,a,1,b,1,q,0,z,1,ifst,ilst,w,1,info)
         CALL CHKXER('DTGEXC',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL DTGEXC(.TRUE.,.FALSE.,1,a,1,b,1,q,1,z,0,ifst,ilst,w,1,    &
     &               info)
         CALL CHKXER('DTGEXC',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL DTGEXC(.TRUE.,.TRUE.,1,a,1,b,1,q,1,z,0,ifst,ilst,w,1,info)
         CALL CHKXER('DTGEXC',INFot,NOUt,LERr,OK)
         INFot = 15
         CALL DTGEXC(.TRUE.,.TRUE.,1,a,1,b,1,q,1,z,1,ifst,ilst,w,0,info)
         CALL CHKXER('DTGEXC',INFot,NOUt,LERr,OK)
         nt = nt + 8
!
!        DTGSEN
!
         SRNamt = 'DTGSEN'
         INFot = 1
         CALL DTGSEN(-1,.TRUE.,.TRUE.,sel,1,a,1,b,1,r1,r2,r3,q,1,z,1,m, &
     &               tola,tolb,rcv,w,1,iw,1,info)
         CALL CHKXER('DTGSEN',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL DTGSEN(1,.TRUE.,.TRUE.,sel,-1,a,1,b,1,r1,r2,r3,q,1,z,1,m, &
     &               tola,tolb,rcv,w,1,iw,1,info)
         CALL CHKXER('DTGSEN',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL DTGSEN(1,.TRUE.,.TRUE.,sel,1,a,0,b,1,r1,r2,r3,q,1,z,1,m,  &
     &               tola,tolb,rcv,w,1,iw,1,info)
         CALL CHKXER('DTGSEN',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL DTGSEN(1,.TRUE.,.TRUE.,sel,1,a,1,b,0,r1,r2,r3,q,1,z,1,m,  &
     &               tola,tolb,rcv,w,1,iw,1,info)
         CALL CHKXER('DTGSEN',INFot,NOUt,LERr,OK)
         INFot = 14
         CALL DTGSEN(1,.TRUE.,.TRUE.,sel,1,a,1,b,1,r1,r2,r3,q,0,z,1,m,  &
     &               tola,tolb,rcv,w,1,iw,1,info)
         CALL CHKXER('DTGSEN',INFot,NOUt,LERr,OK)
         INFot = 16
         CALL DTGSEN(1,.TRUE.,.TRUE.,sel,1,a,1,b,1,r1,r2,r3,q,1,z,0,m,  &
     &               tola,tolb,rcv,w,1,iw,1,info)
         CALL CHKXER('DTGSEN',INFot,NOUt,LERr,OK)
         INFot = 22
         CALL DTGSEN(0,.TRUE.,.TRUE.,sel,1,a,1,b,1,r1,r2,r3,q,1,z,1,m,  &
     &               tola,tolb,rcv,w,1,iw,1,info)
         CALL CHKXER('DTGSEN',INFot,NOUt,LERr,OK)
         INFot = 22
         CALL DTGSEN(1,.TRUE.,.TRUE.,sel,1,a,1,b,1,r1,r2,r3,q,1,z,1,m,  &
     &               tola,tolb,rcv,w,1,iw,1,info)
         CALL CHKXER('DTGSEN',INFot,NOUt,LERr,OK)
         INFot = 22
         CALL DTGSEN(2,.TRUE.,.TRUE.,sel,1,a,1,b,1,r1,r2,r3,q,1,z,1,m,  &
     &               tola,tolb,rcv,w,1,iw,1,info)
         CALL CHKXER('DTGSEN',INFot,NOUt,LERr,OK)
         INFot = 24
         CALL DTGSEN(0,.TRUE.,.TRUE.,sel,1,a,1,b,1,r1,r2,r3,q,1,z,1,m,  &
     &               tola,tolb,rcv,w,20,iw,0,info)
         CALL CHKXER('DTGSEN',INFot,NOUt,LERr,OK)
         INFot = 24
         CALL DTGSEN(1,.TRUE.,.TRUE.,sel,1,a,1,b,1,r1,r2,r3,q,1,z,1,m,  &
     &               tola,tolb,rcv,w,20,iw,0,info)
         CALL CHKXER('DTGSEN',INFot,NOUt,LERr,OK)
         INFot = 24
         CALL DTGSEN(2,.TRUE.,.TRUE.,sel,1,a,1,b,1,r1,r2,r3,q,1,z,1,m,  &
     &               tola,tolb,rcv,w,20,iw,1,info)
         CALL CHKXER('DTGSEN',INFot,NOUt,LERr,OK)
         nt = nt + 12
!
!        DTGSNA
!
         SRNamt = 'DTGSNA'
         INFot = 1
         CALL DTGSNA('/','A',sel,1,a,1,b,1,q,1,u,1,r1,r2,1,m,w,1,iw,    &
     &               info)
         CALL CHKXER('DTGSNA',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DTGSNA('B','/',sel,1,a,1,b,1,q,1,u,1,r1,r2,1,m,w,1,iw,    &
     &               info)
         CALL CHKXER('DTGSNA',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL DTGSNA('B','A',sel,-1,a,1,b,1,q,1,u,1,r1,r2,1,m,w,1,iw,   &
     &               info)
         CALL CHKXER('DTGSNA',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL DTGSNA('B','A',sel,1,a,0,b,1,q,1,u,1,r1,r2,1,m,w,1,iw,    &
     &               info)
         CALL CHKXER('DTGSNA',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL DTGSNA('B','A',sel,1,a,1,b,0,q,1,u,1,r1,r2,1,m,w,1,iw,    &
     &               info)
         CALL CHKXER('DTGSNA',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL DTGSNA('E','A',sel,1,a,1,b,1,q,0,u,1,r1,r2,1,m,w,1,iw,    &
     &               info)
         CALL CHKXER('DTGSNA',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL DTGSNA('E','A',sel,1,a,1,b,1,q,1,u,0,r1,r2,1,m,w,1,iw,    &
     &               info)
         CALL CHKXER('DTGSNA',INFot,NOUt,LERr,OK)
         INFot = 15
         CALL DTGSNA('E','A',sel,1,a,1,b,1,q,1,u,1,r1,r2,0,m,w,1,iw,    &
     &               info)
         CALL CHKXER('DTGSNA',INFot,NOUt,LERr,OK)
         INFot = 18
         CALL DTGSNA('E','A',sel,1,a,1,b,1,q,1,u,1,r1,r2,1,m,w,0,iw,    &
     &               info)
         CALL CHKXER('DTGSNA',INFot,NOUt,LERr,OK)
         nt = nt + 9
!
!        DTGSYL
!
         SRNamt = 'DTGSYL'
         INFot = 1
         CALL DTGSYL('/',0,1,1,a,1,b,1,q,1,u,1,v,1,z,1,scale,dif,w,1,iw,&
     &               info)
         CALL CHKXER('DTGSYL',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL DTGSYL('N',-1,1,1,a,1,b,1,q,1,u,1,v,1,z,1,scale,dif,w,1,  &
     &               iw,info)
         CALL CHKXER('DTGSYL',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL DTGSYL('N',0,0,1,a,1,b,1,q,1,u,1,v,1,z,1,scale,dif,w,1,iw,&
     &               info)
         CALL CHKXER('DTGSYL',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL DTGSYL('N',0,1,0,a,1,b,1,q,1,u,1,v,1,z,1,scale,dif,w,1,iw,&
     &               info)
         CALL CHKXER('DTGSYL',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL DTGSYL('N',0,1,1,a,0,b,1,q,1,u,1,v,1,z,1,scale,dif,w,1,iw,&
     &               info)
         CALL CHKXER('DTGSYL',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL DTGSYL('N',0,1,1,a,1,b,0,q,1,u,1,v,1,z,1,scale,dif,w,1,iw,&
     &               info)
         CALL CHKXER('DTGSYL',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL DTGSYL('N',0,1,1,a,1,b,1,q,0,u,1,v,1,z,1,scale,dif,w,1,iw,&
     &               info)
         CALL CHKXER('DTGSYL',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL DTGSYL('N',0,1,1,a,1,b,1,q,1,u,0,v,1,z,1,scale,dif,w,1,iw,&
     &               info)
         CALL CHKXER('DTGSYL',INFot,NOUt,LERr,OK)
         INFot = 14
         CALL DTGSYL('N',0,1,1,a,1,b,1,q,1,u,1,v,0,z,1,scale,dif,w,1,iw,&
     &               info)
         CALL CHKXER('DTGSYL',INFot,NOUt,LERr,OK)
         INFot = 16
         CALL DTGSYL('N',0,1,1,a,1,b,1,q,1,u,1,v,1,z,0,scale,dif,w,1,iw,&
     &               info)
         CALL CHKXER('DTGSYL',INFot,NOUt,LERr,OK)
         INFot = 20
         CALL DTGSYL('N',1,1,1,a,1,b,1,q,1,u,1,v,1,z,1,scale,dif,w,1,iw,&
     &               info)
         CALL CHKXER('DTGSYL',INFot,NOUt,LERr,OK)
         INFot = 20
         CALL DTGSYL('N',2,1,1,a,1,b,1,q,1,u,1,v,1,z,1,scale,dif,w,1,iw,&
     &               info)
         CALL CHKXER('DTGSYL',INFot,NOUt,LERr,OK)
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
99001 FORMAT (1X,A3,' routines passed the tests of the error exits (',  &
     &        I3,' tests done)')
99002 FORMAT (' *** ',A3,' routines failed the tests of the error ',    &
     &        'exits ***')
!
!
!     End of DERRGG
!
      END SUBROUTINE DERRGG
