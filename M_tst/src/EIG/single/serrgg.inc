!*==serrgg.f90  processed by SPAG 7.51RB at 17:37 on  4 Mar 2022
!> \brief \b serrgg
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE serrgg( PATH, NUNIT )
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
!> serrgg tests the error exits for SGGES, SGGESX, SGGEV, SGGEVX,
!> SGGES3, SGGEV3, SGGGLM, SGGHRD, SGGLSE, SGGQRF, SGGRQF,
!> SGGSVD3, SGGSVP3, SHGEQZ, SORCSD, STGEVC, STGEXC, STGSEN,
!> STGSJA, STGSNA, and STGSYL.
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
      SUBROUTINE serrgg(Path,Nunit)
      IMPLICIT NONE
!*--SERRGG61
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
      REAL ONE , ZERO
      PARAMETER (ONE=1.0E+0,ZERO=0.0E+0)
!     ..
!     .. Local Scalars ..
      CHARACTER*2 c2
      INTEGER dummyk , dummyl , i , ifst , ilo , ihi , ilst , info , j ,&
     &        m , ncycle , nt , sdim , lwork
      REAL anrm , bnrm , dif , scale , tola , tolb
!     ..
!     .. Local Arrays ..
      LOGICAL bw(NMAX) , sel(NMAX)
      INTEGER iw(NMAX) , idum(NMAX)
      REAL a(NMAX,NMAX) , b(NMAX,NMAX) , ls(NMAX) , q(NMAX,NMAX) ,      &
     &     r1(NMAX) , r2(NMAX) , r3(NMAX) , rce(2) , rcv(2) , rs(NMAX) ,&
     &     tau(NMAX) , u(NMAX,NMAX) , v(NMAX,NMAX) , w(LW) ,            &
     &     z(NMAX,NMAX)
!     ..
!     .. External Functions ..
      LOGICAL LSAMEN 
      EXTERNAL LSAMEN 
!     ..
!     .. External Subroutines ..
      EXTERNAL SGGES , SGGESX , SGGEV , SGGEVX , SGGGLM 
     EXTERNAL         SGGHRD , SGGLSE , SGGQRF , SGGRQF , SHGEQZ , SORCSD 
     EXTERNAL         STGEVC , STGEXC , STGSEN , STGSJA , STGSNA , STGSYL 
     EXTERNAL         SGGES3 , SGGEV3 , SGGHD3 , SGGSVD3 , SGGSVP3
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
      tola = 1.0E0
      tolb = 1.0E0
      ifst = 1
      ilst = 1
      nt = 0
      lwork = 1
!
!     Test error exits for the GG path.
!
      IF ( LSAMEN(2,c2,'GG') ) THEN
!
!        SGGHRD
!
         SRNamt = 'SGGHRD'
         INFot = 1
         CALL SGGHRD('/','N',0,1,0,a,1,b,1,q,1,z,1,info)
         CALL chkxer('SGGHRD',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SGGHRD('N','/',0,1,0,a,1,b,1,q,1,z,1,info)
         CALL chkxer('SGGHRD',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SGGHRD('N','N',-1,0,0,a,1,b,1,q,1,z,1,info)
         CALL chkxer('SGGHRD',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL SGGHRD('N','N',0,0,0,a,1,b,1,q,1,z,1,info)
         CALL chkxer('SGGHRD',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL SGGHRD('N','N',0,1,1,a,1,b,1,q,1,z,1,info)
         CALL chkxer('SGGHRD',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL SGGHRD('N','N',2,1,1,a,1,b,2,q,1,z,1,info)
         CALL chkxer('SGGHRD',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL SGGHRD('N','N',2,1,1,a,2,b,1,q,1,z,1,info)
         CALL chkxer('SGGHRD',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL SGGHRD('V','N',2,1,1,a,2,b,2,q,1,z,1,info)
         CALL chkxer('SGGHRD',INFot,NOUt,LERr,OK)
         INFot = 13
         CALL SGGHRD('N','V',2,1,1,a,2,b,2,q,1,z,1,info)
         CALL chkxer('SGGHRD',INFot,NOUt,LERr,OK)
         nt = nt + 9
!
!        SGGHD3
!
         SRNamt = 'SGGHD3'
         INFot = 1
         CALL SGGHD3('/','N',0,1,0,a,1,b,1,q,1,z,1,w,LW,info)
         CALL chkxer('SGGHD3',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SGGHD3('N','/',0,1,0,a,1,b,1,q,1,z,1,w,LW,info)
         CALL chkxer('SGGHD3',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SGGHD3('N','N',-1,0,0,a,1,b,1,q,1,z,1,w,LW,info)
         CALL chkxer('SGGHD3',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL SGGHD3('N','N',0,0,0,a,1,b,1,q,1,z,1,w,LW,info)
         CALL chkxer('SGGHD3',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL SGGHD3('N','N',0,1,1,a,1,b,1,q,1,z,1,w,LW,info)
         CALL chkxer('SGGHD3',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL SGGHD3('N','N',2,1,1,a,1,b,2,q,1,z,1,w,LW,info)
         CALL chkxer('SGGHD3',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL SGGHD3('N','N',2,1,1,a,2,b,1,q,1,z,1,w,LW,info)
         CALL chkxer('SGGHD3',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL SGGHD3('V','N',2,1,1,a,2,b,2,q,1,z,1,w,LW,info)
         CALL chkxer('SGGHD3',INFot,NOUt,LERr,OK)
         INFot = 13
         CALL SGGHD3('N','V',2,1,1,a,2,b,2,q,1,z,1,w,LW,info)
         CALL chkxer('SGGHD3',INFot,NOUt,LERr,OK)
         nt = nt + 9
!
!        SHGEQZ
!
         SRNamt = 'SHGEQZ'
         INFot = 1
         CALL SHGEQZ('/','N','N',0,1,0,a,1,b,1,r1,r2,r3,q,1,z,1,w,LW,   &
     &               info)
         CALL chkxer('SHGEQZ',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SHGEQZ('E','/','N',0,1,0,a,1,b,1,r1,r2,r3,q,1,z,1,w,LW,   &
     &               info)
         CALL chkxer('SHGEQZ',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SHGEQZ('E','N','/',0,1,0,a,1,b,1,r1,r2,r3,q,1,z,1,w,LW,   &
     &               info)
         CALL chkxer('SHGEQZ',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL SHGEQZ('E','N','N',-1,0,0,a,1,b,1,r1,r2,r3,q,1,z,1,w,LW,  &
     &               info)
         CALL chkxer('SHGEQZ',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL SHGEQZ('E','N','N',0,0,0,a,1,b,1,r1,r2,r3,q,1,z,1,w,LW,   &
     &               info)
         CALL chkxer('SHGEQZ',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL SHGEQZ('E','N','N',0,1,1,a,1,b,1,r1,r2,r3,q,1,z,1,w,LW,   &
     &               info)
         CALL chkxer('SHGEQZ',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL SHGEQZ('E','N','N',2,1,1,a,1,b,2,r1,r2,r3,q,1,z,1,w,LW,   &
     &               info)
         CALL chkxer('SHGEQZ',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL SHGEQZ('E','N','N',2,1,1,a,2,b,1,r1,r2,r3,q,1,z,1,w,LW,   &
     &               info)
         CALL chkxer('SHGEQZ',INFot,NOUt,LERr,OK)
         INFot = 15
         CALL SHGEQZ('E','V','N',2,1,1,a,2,b,2,r1,r2,r3,q,1,z,1,w,LW,   &
     &               info)
         CALL chkxer('SHGEQZ',INFot,NOUt,LERr,OK)
         INFot = 17
         CALL SHGEQZ('E','N','V',2,1,1,a,2,b,2,r1,r2,r3,q,1,z,1,w,LW,   &
     &               info)
         CALL chkxer('SHGEQZ',INFot,NOUt,LERr,OK)
         nt = nt + 10
!
!        STGEVC
!
         SRNamt = 'STGEVC'
         INFot = 1
         CALL STGEVC('/','A',sel,0,a,1,b,1,q,1,z,1,0,m,w,info)
         CALL chkxer('STGEVC',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL STGEVC('R','/',sel,0,a,1,b,1,q,1,z,1,0,m,w,info)
         CALL chkxer('STGEVC',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL STGEVC('R','A',sel,-1,a,1,b,1,q,1,z,1,0,m,w,info)
         CALL chkxer('STGEVC',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL STGEVC('R','A',sel,2,a,1,b,2,q,1,z,2,0,m,w,info)
         CALL chkxer('STGEVC',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL STGEVC('R','A',sel,2,a,2,b,1,q,1,z,2,0,m,w,info)
         CALL chkxer('STGEVC',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL STGEVC('L','A',sel,2,a,2,b,2,q,1,z,1,0,m,w,info)
         CALL chkxer('STGEVC',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL STGEVC('R','A',sel,2,a,2,b,2,q,1,z,1,0,m,w,info)
         CALL chkxer('STGEVC',INFot,NOUt,LERr,OK)
         INFot = 13
         CALL STGEVC('R','A',sel,2,a,2,b,2,q,1,z,2,1,m,w,info)
         CALL chkxer('STGEVC',INFot,NOUt,LERr,OK)
         nt = nt + 8
!
!     Test error exits for the GSV path.
!
      ELSEIF ( LSAMEN(3,Path,'GSV') ) THEN
!
!        SGGSVD3
!
         SRNamt = 'SGGSVD3'
         INFot = 1
         CALL SGGSVD3('/','N','N',0,0,0,dummyk,dummyl,a,1,b,1,r1,r2,u,1,&
     &                v,1,q,1,w,lwork,idum,info)
         CALL chkxer('SGGSVD3',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SGGSVD3('N','/','N',0,0,0,dummyk,dummyl,a,1,b,1,r1,r2,u,1,&
     &                v,1,q,1,w,lwork,idum,info)
         CALL chkxer('SGGSVD3',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SGGSVD3('N','N','/',0,0,0,dummyk,dummyl,a,1,b,1,r1,r2,u,1,&
     &                v,1,q,1,w,lwork,idum,info)
         CALL chkxer('SGGSVD3',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL SGGSVD3('N','N','N',-1,0,0,dummyk,dummyl,a,1,b,1,r1,r2,u, &
     &                1,v,1,q,1,w,lwork,idum,info)
         CALL chkxer('SGGSVD3',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL SGGSVD3('N','N','N',0,-1,0,dummyk,dummyl,a,1,b,1,r1,r2,u, &
     &                1,v,1,q,1,w,lwork,idum,info)
         CALL chkxer('SGGSVD3',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL SGGSVD3('N','N','N',0,0,-1,dummyk,dummyl,a,1,b,1,r1,r2,u, &
     &                1,v,1,q,1,w,lwork,idum,info)
         CALL chkxer('SGGSVD3',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL SGGSVD3('N','N','N',2,1,1,dummyk,dummyl,a,1,b,1,r1,r2,u,1,&
     &                v,1,q,1,w,lwork,idum,info)
         CALL chkxer('SGGSVD3',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL SGGSVD3('N','N','N',1,1,2,dummyk,dummyl,a,1,b,1,r1,r2,u,1,&
     &                v,1,q,1,w,lwork,idum,info)
         CALL chkxer('SGGSVD3',INFot,NOUt,LERr,OK)
         INFot = 16
         CALL SGGSVD3('U','N','N',2,2,2,dummyk,dummyl,a,2,b,2,r1,r2,u,1,&
     &                v,1,q,1,w,lwork,idum,info)
         CALL chkxer('SGGSVD3',INFot,NOUt,LERr,OK)
         INFot = 18
         CALL SGGSVD3('N','V','N',1,1,2,dummyk,dummyl,a,1,b,2,r1,r2,u,1,&
     &                v,1,q,1,w,lwork,idum,info)
         CALL chkxer('SGGSVD3',INFot,NOUt,LERr,OK)
         INFot = 20
         CALL SGGSVD3('N','N','Q',1,2,1,dummyk,dummyl,a,1,b,1,r1,r2,u,1,&
     &                v,1,q,1,w,lwork,idum,info)
         CALL chkxer('SGGSVD3',INFot,NOUt,LERr,OK)
         nt = nt + 11
!
!        SGGSVP3
!
         SRNamt = 'SGGSVP3'
         INFot = 1
         CALL SGGSVP3('/','N','N',0,0,0,a,1,b,1,tola,tolb,dummyk,dummyl,&
     &                u,1,v,1,q,1,iw,tau,w,lwork,info)
         CALL chkxer('SGGSVP3',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SGGSVP3('N','/','N',0,0,0,a,1,b,1,tola,tolb,dummyk,dummyl,&
     &                u,1,v,1,q,1,iw,tau,w,lwork,info)
         CALL chkxer('SGGSVP3',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SGGSVP3('N','N','/',0,0,0,a,1,b,1,tola,tolb,dummyk,dummyl,&
     &                u,1,v,1,q,1,iw,tau,w,lwork,info)
         CALL chkxer('SGGSVP3',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL SGGSVP3('N','N','N',-1,0,0,a,1,b,1,tola,tolb,dummyk,      &
     &                dummyl,u,1,v,1,q,1,iw,tau,w,lwork,info)
         CALL chkxer('SGGSVP3',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL SGGSVP3('N','N','N',0,-1,0,a,1,b,1,tola,tolb,dummyk,      &
     &                dummyl,u,1,v,1,q,1,iw,tau,w,lwork,info)
         CALL chkxer('SGGSVP3',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL SGGSVP3('N','N','N',0,0,-1,a,1,b,1,tola,tolb,dummyk,      &
     &                dummyl,u,1,v,1,q,1,iw,tau,w,lwork,info)
         CALL chkxer('SGGSVP3',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL SGGSVP3('N','N','N',2,1,1,a,1,b,1,tola,tolb,dummyk,dummyl,&
     &                u,1,v,1,q,1,iw,tau,w,lwork,info)
         CALL chkxer('SGGSVP3',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL SGGSVP3('N','N','N',1,2,1,a,1,b,1,tola,tolb,dummyk,dummyl,&
     &                u,1,v,1,q,1,iw,tau,w,lwork,info)
         CALL chkxer('SGGSVP3',INFot,NOUt,LERr,OK)
         INFot = 16
         CALL SGGSVP3('U','N','N',2,2,2,a,2,b,2,tola,tolb,dummyk,dummyl,&
     &                u,1,v,1,q,1,iw,tau,w,lwork,info)
         CALL chkxer('SGGSVP3',INFot,NOUt,LERr,OK)
         INFot = 18
         CALL SGGSVP3('N','V','N',1,2,1,a,1,b,2,tola,tolb,dummyk,dummyl,&
     &                u,1,v,1,q,1,iw,tau,w,lwork,info)
         CALL chkxer('SGGSVP3',INFot,NOUt,LERr,OK)
         INFot = 20
         CALL SGGSVP3('N','N','Q',1,1,2,a,1,b,1,tola,tolb,dummyk,dummyl,&
     &                u,1,v,1,q,1,iw,tau,w,lwork,info)
         CALL chkxer('SGGSVP3',INFot,NOUt,LERr,OK)
         nt = nt + 11
!
!        STGSJA
!
         SRNamt = 'STGSJA'
         INFot = 1
         CALL STGSJA('/','N','N',0,0,0,dummyk,dummyl,a,1,b,1,tola,tolb, &
     &               r1,r2,u,1,v,1,q,1,w,ncycle,info)
         CALL chkxer('STGSJA',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL STGSJA('N','/','N',0,0,0,dummyk,dummyl,a,1,b,1,tola,tolb, &
     &               r1,r2,u,1,v,1,q,1,w,ncycle,info)
         CALL chkxer('STGSJA',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL STGSJA('N','N','/',0,0,0,dummyk,dummyl,a,1,b,1,tola,tolb, &
     &               r1,r2,u,1,v,1,q,1,w,ncycle,info)
         CALL chkxer('STGSJA',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL STGSJA('N','N','N',-1,0,0,dummyk,dummyl,a,1,b,1,tola,tolb,&
     &               r1,r2,u,1,v,1,q,1,w,ncycle,info)
         CALL chkxer('STGSJA',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL STGSJA('N','N','N',0,-1,0,dummyk,dummyl,a,1,b,1,tola,tolb,&
     &               r1,r2,u,1,v,1,q,1,w,ncycle,info)
         CALL chkxer('STGSJA',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL STGSJA('N','N','N',0,0,-1,dummyk,dummyl,a,1,b,1,tola,tolb,&
     &               r1,r2,u,1,v,1,q,1,w,ncycle,info)
         CALL chkxer('STGSJA',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL STGSJA('N','N','N',0,0,0,dummyk,dummyl,a,0,b,1,tola,tolb, &
     &               r1,r2,u,1,v,1,q,1,w,ncycle,info)
         CALL chkxer('STGSJA',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL STGSJA('N','N','N',0,0,0,dummyk,dummyl,a,1,b,0,tola,tolb, &
     &               r1,r2,u,1,v,1,q,1,w,ncycle,info)
         CALL chkxer('STGSJA',INFot,NOUt,LERr,OK)
         INFot = 18
         CALL STGSJA('U','N','N',0,0,0,dummyk,dummyl,a,1,b,1,tola,tolb, &
     &               r1,r2,u,0,v,1,q,1,w,ncycle,info)
         CALL chkxer('STGSJA',INFot,NOUt,LERr,OK)
         INFot = 20
         CALL STGSJA('N','V','N',0,0,0,dummyk,dummyl,a,1,b,1,tola,tolb, &
     &               r1,r2,u,1,v,0,q,1,w,ncycle,info)
         CALL chkxer('STGSJA',INFot,NOUt,LERr,OK)
         INFot = 22
         CALL STGSJA('N','N','Q',0,0,0,dummyk,dummyl,a,1,b,1,tola,tolb, &
     &               r1,r2,u,1,v,1,q,0,w,ncycle,info)
         CALL chkxer('STGSJA',INFot,NOUt,LERr,OK)
         nt = nt + 11
!
!     Test error exits for the GLM path.
!
      ELSEIF ( LSAMEN(3,Path,'GLM') ) THEN
!
!        SGGGLM
!
         SRNamt = 'SGGGLM'
         INFot = 1
         CALL SGGGLM(-1,0,0,a,1,b,1,r1,r2,r3,w,LW,info)
         CALL chkxer('SGGGLM',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SGGGLM(0,-1,0,a,1,b,1,r1,r2,r3,w,LW,info)
         CALL chkxer('SGGGLM',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SGGGLM(0,1,0,a,1,b,1,r1,r2,r3,w,LW,info)
         CALL chkxer('SGGGLM',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SGGGLM(0,0,-1,a,1,b,1,r1,r2,r3,w,LW,info)
         CALL chkxer('SGGGLM',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SGGGLM(1,0,0,a,1,b,1,r1,r2,r3,w,LW,info)
         CALL chkxer('SGGGLM',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL SGGGLM(0,0,0,a,0,b,1,r1,r2,r3,w,LW,info)
         CALL chkxer('SGGGLM',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL SGGGLM(0,0,0,a,1,b,0,r1,r2,r3,w,LW,info)
         CALL chkxer('SGGGLM',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL SGGGLM(1,1,1,a,1,b,1,r1,r2,r3,w,1,info)
         CALL chkxer('SGGGLM',INFot,NOUt,LERr,OK)
         nt = nt + 8
!
!     Test error exits for the LSE path.
!
      ELSEIF ( LSAMEN(3,Path,'LSE') ) THEN
!
!        SGGLSE
!
         SRNamt = 'SGGLSE'
         INFot = 1
         CALL SGGLSE(-1,0,0,a,1,b,1,r1,r2,r3,w,LW,info)
         CALL chkxer('SGGLSE',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SGGLSE(0,-1,0,a,1,b,1,r1,r2,r3,w,LW,info)
         CALL chkxer('SGGLSE',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SGGLSE(0,0,-1,a,1,b,1,r1,r2,r3,w,LW,info)
         CALL chkxer('SGGLSE',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SGGLSE(0,0,1,a,1,b,1,r1,r2,r3,w,LW,info)
         CALL chkxer('SGGLSE',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SGGLSE(0,1,0,a,1,b,1,r1,r2,r3,w,LW,info)
         CALL chkxer('SGGLSE',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL SGGLSE(0,0,0,a,0,b,1,r1,r2,r3,w,LW,info)
         CALL chkxer('SGGLSE',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL SGGLSE(0,0,0,a,1,b,0,r1,r2,r3,w,LW,info)
         CALL chkxer('SGGLSE',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL SGGLSE(1,1,1,a,1,b,1,r1,r2,r3,w,1,info)
         CALL chkxer('SGGLSE',INFot,NOUt,LERr,OK)
         nt = nt + 8
!
!     Test error exits for the CSD path.
!
      ELSEIF ( LSAMEN(3,Path,'CSD') ) THEN
!
!        SORCSD
!
         SRNamt = 'SORCSD'
         INFot = 7
         CALL SORCSD('Y','Y','Y','Y','N','N',-1,0,0,a,1,a,1,a,1,a,1,a,a,&
     &               1,a,1,a,1,a,1,w,LW,iw,info)
         CALL chkxer('SORCSD',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL SORCSD('Y','Y','Y','Y','N','N',1,-1,0,a,1,a,1,a,1,a,1,a,a,&
     &               1,a,1,a,1,a,1,w,LW,iw,info)
         CALL chkxer('SORCSD',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL SORCSD('Y','Y','Y','Y','N','N',1,1,-1,a,1,a,1,a,1,a,1,a,a,&
     &               1,a,1,a,1,a,1,w,LW,iw,info)
         CALL chkxer('SORCSD',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL SORCSD('Y','Y','Y','Y','N','N',1,1,1,a,-1,a,1,a,1,a,1,a,a,&
     &               1,a,1,a,1,a,1,w,LW,iw,info)
         CALL chkxer('SORCSD',INFot,NOUt,LERr,OK)
         INFot = 20
         CALL SORCSD('Y','Y','Y','Y','N','N',1,1,1,a,1,a,1,a,1,a,1,a,a, &
     &               -1,a,1,a,1,a,1,w,LW,iw,info)
         CALL chkxer('SORCSD',INFot,NOUt,LERr,OK)
         INFot = 22
         CALL SORCSD('Y','Y','Y','Y','N','N',1,1,1,a,1,a,1,a,1,a,1,a,a, &
     &               1,a,-1,a,1,a,1,w,LW,iw,info)
         CALL chkxer('SORCSD',INFot,NOUt,LERr,OK)
         INFot = 24
         CALL SORCSD('Y','Y','Y','Y','N','N',1,1,1,a,1,a,1,a,1,a,1,a,a, &
     &               1,a,1,a,-1,a,1,w,LW,iw,info)
         CALL chkxer('SORCSD',INFot,NOUt,LERr,OK)
         INFot = 26
         CALL SORCSD('Y','Y','Y','Y','N','N',1,1,1,a,1,a,1,a,1,a,1,a,a, &
     &               1,a,1,a,1,a,-1,w,LW,iw,info)
         CALL chkxer('SORCSD',INFot,NOUt,LERr,OK)
         nt = nt + 8
!
!     Test error exits for the GQR path.
!
      ELSEIF ( LSAMEN(3,Path,'GQR') ) THEN
!
!        SGGQRF
!
         SRNamt = 'SGGQRF'
         INFot = 1
         CALL SGGQRF(-1,0,0,a,1,r1,b,1,r2,w,LW,info)
         CALL chkxer('SGGQRF',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SGGQRF(0,-1,0,a,1,r1,b,1,r2,w,LW,info)
         CALL chkxer('SGGQRF',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SGGQRF(0,0,-1,a,1,r1,b,1,r2,w,LW,info)
         CALL chkxer('SGGQRF',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL SGGQRF(0,0,0,a,0,r1,b,1,r2,w,LW,info)
         CALL chkxer('SGGQRF',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL SGGQRF(0,0,0,a,1,r1,b,0,r2,w,LW,info)
         CALL chkxer('SGGQRF',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL SGGQRF(1,1,2,a,1,r1,b,1,r2,w,1,info)
         CALL chkxer('SGGQRF',INFot,NOUt,LERr,OK)
         nt = nt + 6
!
!        SGGRQF
!
         SRNamt = 'SGGRQF'
         INFot = 1
         CALL SGGRQF(-1,0,0,a,1,r1,b,1,r2,w,LW,info)
         CALL chkxer('SGGRQF',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SGGRQF(0,-1,0,a,1,r1,b,1,r2,w,LW,info)
         CALL chkxer('SGGRQF',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SGGRQF(0,0,-1,a,1,r1,b,1,r2,w,LW,info)
         CALL chkxer('SGGRQF',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL SGGRQF(0,0,0,a,0,r1,b,1,r2,w,LW,info)
         CALL chkxer('SGGRQF',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL SGGRQF(0,0,0,a,1,r1,b,0,r2,w,LW,info)
         CALL chkxer('SGGRQF',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL SGGRQF(1,1,2,a,1,r1,b,1,r2,w,1,info)
         CALL chkxer('SGGRQF',INFot,NOUt,LERr,OK)
         nt = nt + 6
!
!     Test error exits for the SGS, SGV, SGX, and SXV paths.
!
      ELSEIF ( LSAMEN(3,Path,'SGS') .OR. LSAMEN(3,Path,'SGV') .OR.      &
     &         LSAMEN(3,Path,'SGX') .OR. LSAMEN(3,Path,'SXV') ) THEN
!
!        SGGES
!
         SRNamt = 'SGGES '
         INFot = 1
         CALL SGGES('/','N','S',slctes,1,a,1,b,1,sdim,r1,r2,r3,q,1,u,1, &
     &              w,1,bw,info)
         CALL chkxer('SGGES ',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SGGES('N','/','S',slctes,1,a,1,b,1,sdim,r1,r2,r3,q,1,u,1, &
     &              w,1,bw,info)
         CALL chkxer('SGGES ',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SGGES('N','V','/',slctes,1,a,1,b,1,sdim,r1,r2,r3,q,1,u,1, &
     &              w,1,bw,info)
         CALL chkxer('SGGES ',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL SGGES('N','V','S',slctes,-1,a,1,b,1,sdim,r1,r2,r3,q,1,u,1,&
     &              w,1,bw,info)
         CALL chkxer('SGGES ',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL SGGES('N','V','S',slctes,1,a,0,b,1,sdim,r1,r2,r3,q,1,u,1, &
     &              w,1,bw,info)
         CALL chkxer('SGGES ',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL SGGES('N','V','S',slctes,1,a,1,b,0,sdim,r1,r2,r3,q,1,u,1, &
     &              w,1,bw,info)
         CALL chkxer('SGGES ',INFot,NOUt,LERr,OK)
         INFot = 15
         CALL SGGES('N','V','S',slctes,1,a,1,b,1,sdim,r1,r2,r3,q,0,u,1, &
     &              w,1,bw,info)
         CALL chkxer('SGGES ',INFot,NOUt,LERr,OK)
         INFot = 15
         CALL SGGES('V','V','S',slctes,2,a,2,b,2,sdim,r1,r2,r3,q,1,u,2, &
     &              w,1,bw,info)
         CALL chkxer('SGGES ',INFot,NOUt,LERr,OK)
         INFot = 17
         CALL SGGES('N','V','S',slctes,1,a,1,b,1,sdim,r1,r2,r3,q,1,u,0, &
     &              w,1,bw,info)
         CALL chkxer('SGGES ',INFot,NOUt,LERr,OK)
         INFot = 17
         CALL SGGES('V','V','S',slctes,2,a,2,b,2,sdim,r1,r2,r3,q,2,u,1, &
     &              w,1,bw,info)
         CALL chkxer('SGGES ',INFot,NOUt,LERr,OK)
         INFot = 19
         CALL SGGES('V','V','S',slctes,2,a,2,b,2,sdim,r1,r2,r3,q,2,u,2, &
     &              w,1,bw,info)
         CALL chkxer('SGGES ',INFot,NOUt,LERr,OK)
         nt = nt + 11
!
!        SGGES3
!
         SRNamt = 'SGGES3'
         INFot = 1
         CALL SGGES3('/','N','S',slctes,1,a,1,b,1,sdim,r1,r2,r3,q,1,u,1,&
     &               w,1,bw,info)
         CALL chkxer('SGGES3 ',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SGGES3('N','/','S',slctes,1,a,1,b,1,sdim,r1,r2,r3,q,1,u,1,&
     &               w,1,bw,info)
         CALL chkxer('SGGES3 ',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SGGES3('N','V','/',slctes,1,a,1,b,1,sdim,r1,r2,r3,q,1,u,1,&
     &               w,1,bw,info)
         CALL chkxer('SGGES3 ',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL SGGES3('N','V','S',slctes,-1,a,1,b,1,sdim,r1,r2,r3,q,1,u, &
     &               1,w,1,bw,info)
         CALL chkxer('SGGES3 ',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL SGGES3('N','V','S',slctes,1,a,0,b,1,sdim,r1,r2,r3,q,1,u,1,&
     &               w,1,bw,info)
         CALL chkxer('SGGES3 ',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL SGGES3('N','V','S',slctes,1,a,1,b,0,sdim,r1,r2,r3,q,1,u,1,&
     &               w,1,bw,info)
         CALL chkxer('SGGES3 ',INFot,NOUt,LERr,OK)
         INFot = 15
         CALL SGGES3('N','V','S',slctes,1,a,1,b,1,sdim,r1,r2,r3,q,0,u,1,&
     &               w,1,bw,info)
         CALL chkxer('SGGES3 ',INFot,NOUt,LERr,OK)
         INFot = 15
         CALL SGGES3('V','V','S',slctes,2,a,2,b,2,sdim,r1,r2,r3,q,1,u,2,&
     &               w,1,bw,info)
         CALL chkxer('SGGES3 ',INFot,NOUt,LERr,OK)
         INFot = 17
         CALL SGGES3('N','V','S',slctes,1,a,1,b,1,sdim,r1,r2,r3,q,1,u,0,&
     &               w,1,bw,info)
         CALL chkxer('SGGES3 ',INFot,NOUt,LERr,OK)
         INFot = 17
         CALL SGGES3('V','V','S',slctes,2,a,2,b,2,sdim,r1,r2,r3,q,2,u,1,&
     &               w,1,bw,info)
         CALL chkxer('SGGES3 ',INFot,NOUt,LERr,OK)
         INFot = 19
         CALL SGGES3('V','V','S',slctes,2,a,2,b,2,sdim,r1,r2,r3,q,2,u,2,&
     &               w,1,bw,info)
         CALL chkxer('SGGES3 ',INFot,NOUt,LERr,OK)
         nt = nt + 11
!
!        SGGESX
!
         SRNamt = 'SGGESX'
         INFot = 1
         CALL SGGESX('/','N','S',slctsx,'N',1,a,1,b,1,sdim,r1,r2,r3,q,1,&
     &               u,1,rce,rcv,w,1,iw,1,bw,info)
         CALL chkxer('SGGESX',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SGGESX('N','/','S',slctsx,'N',1,a,1,b,1,sdim,r1,r2,r3,q,1,&
     &               u,1,rce,rcv,w,1,iw,1,bw,info)
         CALL chkxer('SGGESX',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SGGESX('V','V','/',slctsx,'N',1,a,1,b,1,sdim,r1,r2,r3,q,1,&
     &               u,1,rce,rcv,w,1,iw,1,bw,info)
         CALL chkxer('SGGESX',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL SGGESX('V','V','S',slctsx,'/',1,a,1,b,1,sdim,r1,r2,r3,q,1,&
     &               u,1,rce,rcv,w,1,iw,1,bw,info)
         CALL chkxer('SGGESX',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL SGGESX('V','V','S',slctsx,'B',-1,a,1,b,1,sdim,r1,r2,r3,q, &
     &               1,u,1,rce,rcv,w,1,iw,1,bw,info)
         CALL chkxer('SGGESX',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL SGGESX('V','V','S',slctsx,'B',1,a,0,b,1,sdim,r1,r2,r3,q,1,&
     &               u,1,rce,rcv,w,1,iw,1,bw,info)
         CALL chkxer('SGGESX',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL SGGESX('V','V','S',slctsx,'B',1,a,1,b,0,sdim,r1,r2,r3,q,1,&
     &               u,1,rce,rcv,w,1,iw,1,bw,info)
         CALL chkxer('SGGESX',INFot,NOUt,LERr,OK)
         INFot = 16
         CALL SGGESX('V','V','S',slctsx,'B',1,a,1,b,1,sdim,r1,r2,r3,q,0,&
     &               u,1,rce,rcv,w,1,iw,1,bw,info)
         CALL chkxer('SGGESX',INFot,NOUt,LERr,OK)
         INFot = 16
         CALL SGGESX('V','V','S',slctsx,'B',2,a,2,b,2,sdim,r1,r2,r3,q,1,&
     &               u,1,rce,rcv,w,1,iw,1,bw,info)
         CALL chkxer('SGGESX',INFot,NOUt,LERr,OK)
         INFot = 18
         CALL SGGESX('V','V','S',slctsx,'B',1,a,1,b,1,sdim,r1,r2,r3,q,1,&
     &               u,0,rce,rcv,w,1,iw,1,bw,info)
         CALL chkxer('SGGESX',INFot,NOUt,LERr,OK)
         INFot = 18
         CALL SGGESX('V','V','S',slctsx,'B',2,a,2,b,2,sdim,r1,r2,r3,q,2,&
     &               u,1,rce,rcv,w,1,iw,1,bw,info)
         CALL chkxer('SGGESX',INFot,NOUt,LERr,OK)
         INFot = 22
         CALL SGGESX('V','V','S',slctsx,'B',2,a,2,b,2,sdim,r1,r2,r3,q,2,&
     &               u,2,rce,rcv,w,1,iw,1,bw,info)
         CALL chkxer('SGGESX',INFot,NOUt,LERr,OK)
         INFot = 24
         CALL SGGESX('V','V','S',slctsx,'V',1,a,1,b,1,sdim,r1,r2,r3,q,1,&
     &               u,1,rce,rcv,w,32,iw,0,bw,info)
         CALL chkxer('SGGESX',INFot,NOUt,LERr,OK)
         nt = nt + 13
!
!        SGGEV
!
         SRNamt = 'SGGEV '
         INFot = 1
         CALL SGGEV('/','N',1,a,1,b,1,r1,r2,r3,q,1,u,1,w,1,info)
         CALL chkxer('SGGEV ',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SGGEV('N','/',1,a,1,b,1,r1,r2,r3,q,1,u,1,w,1,info)
         CALL chkxer('SGGEV ',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SGGEV('V','V',-1,a,1,b,1,r1,r2,r3,q,1,u,1,w,1,info)
         CALL chkxer('SGGEV ',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL SGGEV('V','V',1,a,0,b,1,r1,r2,r3,q,1,u,1,w,1,info)
         CALL chkxer('SGGEV ',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL SGGEV('V','V',1,a,1,b,0,r1,r2,r3,q,1,u,1,w,1,info)
         CALL chkxer('SGGEV ',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL SGGEV('N','V',1,a,1,b,1,r1,r2,r3,q,0,u,1,w,1,info)
         CALL chkxer('SGGEV ',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL SGGEV('V','V',2,a,2,b,2,r1,r2,r3,q,1,u,2,w,1,info)
         CALL chkxer('SGGEV ',INFot,NOUt,LERr,OK)
         INFot = 14
         CALL SGGEV('V','N',2,a,2,b,2,r1,r2,r3,q,2,u,0,w,1,info)
         CALL chkxer('SGGEV ',INFot,NOUt,LERr,OK)
         INFot = 14
         CALL SGGEV('V','V',2,a,2,b,2,r1,r2,r3,q,2,u,1,w,1,info)
         CALL chkxer('SGGEV ',INFot,NOUt,LERr,OK)
         INFot = 16
         CALL SGGEV('V','V',1,a,1,b,1,r1,r2,r3,q,1,u,1,w,1,info)
         CALL chkxer('SGGEV ',INFot,NOUt,LERr,OK)
         nt = nt + 10
!
!        SGGEV3
!
         SRNamt = 'SGGEV3 '
         INFot = 1
         CALL SGGEV3('/','N',1,a,1,b,1,r1,r2,r3,q,1,u,1,w,1,info)
         CALL chkxer('SGGEV3 ',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SGGEV3('N','/',1,a,1,b,1,r1,r2,r3,q,1,u,1,w,1,info)
         CALL chkxer('SGGEV3 ',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SGGEV3('V','V',-1,a,1,b,1,r1,r2,r3,q,1,u,1,w,1,info)
         CALL chkxer('SGGEV3 ',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL SGGEV3('V','V',1,a,0,b,1,r1,r2,r3,q,1,u,1,w,1,info)
         CALL chkxer('SGGEV3 ',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL SGGEV3('V','V',1,a,1,b,0,r1,r2,r3,q,1,u,1,w,1,info)
         CALL chkxer('SGGEV3 ',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL SGGEV3('N','V',1,a,1,b,1,r1,r2,r3,q,0,u,1,w,1,info)
         CALL chkxer('SGGEV3 ',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL SGGEV3('V','V',2,a,2,b,2,r1,r2,r3,q,1,u,2,w,1,info)
         CALL chkxer('SGGEV3 ',INFot,NOUt,LERr,OK)
         INFot = 14
         CALL SGGEV3('V','N',2,a,2,b,2,r1,r2,r3,q,2,u,0,w,1,info)
         CALL chkxer('SGGEV3 ',INFot,NOUt,LERr,OK)
         INFot = 14
         CALL SGGEV3('V','V',2,a,2,b,2,r1,r2,r3,q,2,u,1,w,1,info)
         CALL chkxer('SGGEV3 ',INFot,NOUt,LERr,OK)
         INFot = 16
         CALL SGGEV3('V','V',1,a,1,b,1,r1,r2,r3,q,1,u,1,w,1,info)
         CALL chkxer('SGGEV3 ',INFot,NOUt,LERr,OK)
         nt = nt + 10
!
!        SGGEVX
!
         SRNamt = 'SGGEVX'
         INFot = 1
         CALL SGGEVX('/','N','N','N',1,a,1,b,1,r1,r2,r3,q,1,u,1,ilo,ihi,&
     &               ls,rs,anrm,bnrm,rce,rcv,w,1,iw,bw,info)
         CALL chkxer('SGGEVX',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL SGGEVX('N','/','N','N',1,a,1,b,1,r1,r2,r3,q,1,u,1,ilo,ihi,&
     &               ls,rs,anrm,bnrm,rce,rcv,w,1,iw,bw,info)
         CALL chkxer('SGGEVX',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL SGGEVX('N','N','/','N',1,a,1,b,1,r1,r2,r3,q,1,u,1,ilo,ihi,&
     &               ls,rs,anrm,bnrm,rce,rcv,w,1,iw,bw,info)
         CALL chkxer('SGGEVX',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL SGGEVX('N','N','N','/',1,a,1,b,1,r1,r2,r3,q,1,u,1,ilo,ihi,&
     &               ls,rs,anrm,bnrm,rce,rcv,w,1,iw,bw,info)
         CALL chkxer('SGGEVX',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL SGGEVX('N','N','N','N',-1,a,1,b,1,r1,r2,r3,q,1,u,1,ilo,   &
     &               ihi,ls,rs,anrm,bnrm,rce,rcv,w,1,iw,bw,info)
         CALL chkxer('SGGEVX',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL SGGEVX('N','N','N','N',1,a,0,b,1,r1,r2,r3,q,1,u,1,ilo,ihi,&
     &               ls,rs,anrm,bnrm,rce,rcv,w,1,iw,bw,info)
         CALL chkxer('SGGEVX',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL SGGEVX('N','N','N','N',1,a,1,b,0,r1,r2,r3,q,1,u,1,ilo,ihi,&
     &               ls,rs,anrm,bnrm,rce,rcv,w,1,iw,bw,info)
         CALL chkxer('SGGEVX',INFot,NOUt,LERr,OK)
         INFot = 14
         CALL SGGEVX('N','N','N','N',1,a,1,b,1,r1,r2,r3,q,0,u,1,ilo,ihi,&
     &               ls,rs,anrm,bnrm,rce,rcv,w,1,iw,bw,info)
         CALL chkxer('SGGEVX',INFot,NOUt,LERr,OK)
         INFot = 14
         CALL SGGEVX('N','V','N','N',2,a,2,b,2,r1,r2,r3,q,1,u,2,ilo,ihi,&
     &               ls,rs,anrm,bnrm,rce,rcv,w,1,iw,bw,info)
         CALL chkxer('SGGEVX',INFot,NOUt,LERr,OK)
         INFot = 16
         CALL SGGEVX('N','N','N','N',1,a,1,b,1,r1,r2,r3,q,1,u,0,ilo,ihi,&
     &               ls,rs,anrm,bnrm,rce,rcv,w,1,iw,bw,info)
         CALL chkxer('SGGEVX',INFot,NOUt,LERr,OK)
         INFot = 16
         CALL SGGEVX('N','N','V','N',2,a,2,b,2,r1,r2,r3,q,2,u,1,ilo,ihi,&
     &               ls,rs,anrm,bnrm,rce,rcv,w,1,iw,bw,info)
         CALL chkxer('SGGEVX',INFot,NOUt,LERr,OK)
         INFot = 26
         CALL SGGEVX('N','N','V','N',2,a,2,b,2,r1,r2,r3,q,2,u,2,ilo,ihi,&
     &               ls,rs,anrm,bnrm,rce,rcv,w,1,iw,bw,info)
         CALL chkxer('SGGEVX',INFot,NOUt,LERr,OK)
         nt = nt + 12
!
!        STGEXC
!
         SRNamt = 'STGEXC'
         INFot = 3
         CALL STGEXC(.TRUE.,.TRUE.,-1,a,1,b,1,q,1,z,1,ifst,ilst,w,1,    &
     &               info)
         CALL chkxer('STGEXC',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL STGEXC(.TRUE.,.TRUE.,1,a,0,b,1,q,1,z,1,ifst,ilst,w,1,info)
         CALL chkxer('STGEXC',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL STGEXC(.TRUE.,.TRUE.,1,a,1,b,0,q,1,z,1,ifst,ilst,w,1,info)
         CALL chkxer('STGEXC',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL STGEXC(.FALSE.,.TRUE.,1,a,1,b,1,q,0,z,1,ifst,ilst,w,1,    &
     &               info)
         CALL chkxer('STGEXC',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL STGEXC(.TRUE.,.TRUE.,1,a,1,b,1,q,0,z,1,ifst,ilst,w,1,info)
         CALL chkxer('STGEXC',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL STGEXC(.TRUE.,.FALSE.,1,a,1,b,1,q,1,z,0,ifst,ilst,w,1,    &
     &               info)
         CALL chkxer('STGEXC',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL STGEXC(.TRUE.,.TRUE.,1,a,1,b,1,q,1,z,0,ifst,ilst,w,1,info)
         CALL chkxer('STGEXC',INFot,NOUt,LERr,OK)
         INFot = 15
         CALL STGEXC(.TRUE.,.TRUE.,1,a,1,b,1,q,1,z,1,ifst,ilst,w,0,info)
         CALL chkxer('STGEXC',INFot,NOUt,LERr,OK)
         nt = nt + 8
!
!        STGSEN
!
         SRNamt = 'STGSEN'
         INFot = 1
         CALL STGSEN(-1,.TRUE.,.TRUE.,sel,1,a,1,b,1,r1,r2,r3,q,1,z,1,m, &
     &               tola,tolb,rcv,w,1,iw,1,info)
         CALL chkxer('STGSEN',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL STGSEN(1,.TRUE.,.TRUE.,sel,-1,a,1,b,1,r1,r2,r3,q,1,z,1,m, &
     &               tola,tolb,rcv,w,1,iw,1,info)
         CALL chkxer('STGSEN',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL STGSEN(1,.TRUE.,.TRUE.,sel,1,a,0,b,1,r1,r2,r3,q,1,z,1,m,  &
     &               tola,tolb,rcv,w,1,iw,1,info)
         CALL chkxer('STGSEN',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL STGSEN(1,.TRUE.,.TRUE.,sel,1,a,1,b,0,r1,r2,r3,q,1,z,1,m,  &
     &               tola,tolb,rcv,w,1,iw,1,info)
         CALL chkxer('STGSEN',INFot,NOUt,LERr,OK)
         INFot = 14
         CALL STGSEN(1,.TRUE.,.TRUE.,sel,1,a,1,b,1,r1,r2,r3,q,0,z,1,m,  &
     &               tola,tolb,rcv,w,1,iw,1,info)
         CALL chkxer('STGSEN',INFot,NOUt,LERr,OK)
         INFot = 16
         CALL STGSEN(1,.TRUE.,.TRUE.,sel,1,a,1,b,1,r1,r2,r3,q,1,z,0,m,  &
     &               tola,tolb,rcv,w,1,iw,1,info)
         CALL chkxer('STGSEN',INFot,NOUt,LERr,OK)
         INFot = 22
         CALL STGSEN(0,.TRUE.,.TRUE.,sel,1,a,1,b,1,r1,r2,r3,q,1,z,1,m,  &
     &               tola,tolb,rcv,w,1,iw,1,info)
         CALL chkxer('STGSEN',INFot,NOUt,LERr,OK)
         INFot = 22
         CALL STGSEN(1,.TRUE.,.TRUE.,sel,1,a,1,b,1,r1,r2,r3,q,1,z,1,m,  &
     &               tola,tolb,rcv,w,1,iw,1,info)
         CALL chkxer('STGSEN',INFot,NOUt,LERr,OK)
         INFot = 22
         CALL STGSEN(2,.TRUE.,.TRUE.,sel,1,a,1,b,1,r1,r2,r3,q,1,z,1,m,  &
     &               tola,tolb,rcv,w,1,iw,1,info)
         CALL chkxer('STGSEN',INFot,NOUt,LERr,OK)
         INFot = 24
         CALL STGSEN(0,.TRUE.,.TRUE.,sel,1,a,1,b,1,r1,r2,r3,q,1,z,1,m,  &
     &               tola,tolb,rcv,w,20,iw,0,info)
         CALL chkxer('STGSEN',INFot,NOUt,LERr,OK)
         INFot = 24
         CALL STGSEN(1,.TRUE.,.TRUE.,sel,1,a,1,b,1,r1,r2,r3,q,1,z,1,m,  &
     &               tola,tolb,rcv,w,20,iw,0,info)
         CALL chkxer('STGSEN',INFot,NOUt,LERr,OK)
         INFot = 24
         CALL STGSEN(2,.TRUE.,.TRUE.,sel,1,a,1,b,1,r1,r2,r3,q,1,z,1,m,  &
     &               tola,tolb,rcv,w,20,iw,1,info)
         CALL chkxer('STGSEN',INFot,NOUt,LERr,OK)
         nt = nt + 12
!
!        STGSNA
!
         SRNamt = 'STGSNA'
         INFot = 1
         CALL STGSNA('/','A',sel,1,a,1,b,1,q,1,u,1,r1,r2,1,m,w,1,iw,    &
     &               info)
         CALL chkxer('STGSNA',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL STGSNA('B','/',sel,1,a,1,b,1,q,1,u,1,r1,r2,1,m,w,1,iw,    &
     &               info)
         CALL chkxer('STGSNA',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL STGSNA('B','A',sel,-1,a,1,b,1,q,1,u,1,r1,r2,1,m,w,1,iw,   &
     &               info)
         CALL chkxer('STGSNA',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL STGSNA('B','A',sel,1,a,0,b,1,q,1,u,1,r1,r2,1,m,w,1,iw,    &
     &               info)
         CALL chkxer('STGSNA',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL STGSNA('B','A',sel,1,a,1,b,0,q,1,u,1,r1,r2,1,m,w,1,iw,    &
     &               info)
         CALL chkxer('STGSNA',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL STGSNA('E','A',sel,1,a,1,b,1,q,0,u,1,r1,r2,1,m,w,1,iw,    &
     &               info)
         CALL chkxer('STGSNA',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL STGSNA('E','A',sel,1,a,1,b,1,q,1,u,0,r1,r2,1,m,w,1,iw,    &
     &               info)
         CALL chkxer('STGSNA',INFot,NOUt,LERr,OK)
         INFot = 15
         CALL STGSNA('E','A',sel,1,a,1,b,1,q,1,u,1,r1,r2,0,m,w,1,iw,    &
     &               info)
         CALL chkxer('STGSNA',INFot,NOUt,LERr,OK)
         INFot = 18
         CALL STGSNA('E','A',sel,1,a,1,b,1,q,1,u,1,r1,r2,1,m,w,0,iw,    &
     &               info)
         CALL chkxer('STGSNA',INFot,NOUt,LERr,OK)
         nt = nt + 9
!
!        STGSYL
!
         SRNamt = 'STGSYL'
         INFot = 1
         CALL STGSYL('/',0,1,1,a,1,b,1,q,1,u,1,v,1,z,1,scale,dif,w,1,iw,&
     &               info)
         CALL chkxer('STGSYL',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL STGSYL('N',-1,1,1,a,1,b,1,q,1,u,1,v,1,z,1,scale,dif,w,1,  &
     &               iw,info)
         CALL chkxer('STGSYL',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL STGSYL('N',0,0,1,a,1,b,1,q,1,u,1,v,1,z,1,scale,dif,w,1,iw,&
     &               info)
         CALL chkxer('STGSYL',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL STGSYL('N',0,1,0,a,1,b,1,q,1,u,1,v,1,z,1,scale,dif,w,1,iw,&
     &               info)
         CALL chkxer('STGSYL',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL STGSYL('N',0,1,1,a,0,b,1,q,1,u,1,v,1,z,1,scale,dif,w,1,iw,&
     &               info)
         CALL chkxer('STGSYL',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL STGSYL('N',0,1,1,a,1,b,0,q,1,u,1,v,1,z,1,scale,dif,w,1,iw,&
     &               info)
         CALL chkxer('STGSYL',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL STGSYL('N',0,1,1,a,1,b,1,q,0,u,1,v,1,z,1,scale,dif,w,1,iw,&
     &               info)
         CALL chkxer('STGSYL',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL STGSYL('N',0,1,1,a,1,b,1,q,1,u,0,v,1,z,1,scale,dif,w,1,iw,&
     &               info)
         CALL chkxer('STGSYL',INFot,NOUt,LERr,OK)
         INFot = 14
         CALL STGSYL('N',0,1,1,a,1,b,1,q,1,u,1,v,0,z,1,scale,dif,w,1,iw,&
     &               info)
         CALL chkxer('STGSYL',INFot,NOUt,LERr,OK)
         INFot = 16
         CALL STGSYL('N',0,1,1,a,1,b,1,q,1,u,1,v,1,z,0,scale,dif,w,1,iw,&
     &               info)
         CALL chkxer('STGSYL',INFot,NOUt,LERr,OK)
         INFot = 20
         CALL STGSYL('N',1,1,1,a,1,b,1,q,1,u,1,v,1,z,1,scale,dif,w,1,iw,&
     &               info)
         CALL chkxer('STGSYL',INFot,NOUt,LERr,OK)
         INFot = 20
         CALL STGSYL('N',2,1,1,a,1,b,1,q,1,u,1,v,1,z,1,scale,dif,w,1,iw,&
     &               info)
         CALL chkxer('STGSYL',INFot,NOUt,LERr,OK)
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
!     End of serrgg
!
      END SUBROUTINE serrgg
