!*==cerrgg.f90  processed by SPAG 7.51RB at 17:37 on  4 Mar 2022
!> \brief \b cerrgg
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE cerrgg( PATH, NUNIT )
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
!> cerrgg tests the error exits for CGGES, CGGESX, CGGEV, CGGEVX,
!> CGGES3, CGGEV3, CGGGLM, CGGHRD, CGGLSE, CGGQRF, CGGRQF,
!> CGGSVD3, CGGSVP3, CHGEQZ, CTGEVC, CTGEXC, CTGSEN, CTGSJA,
!> CTGSNA, CTGSYL, and CUNCSD.
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
!> \ingroup complex_eig
!
!  =====================================================================
      SUBROUTINE cerrgg(Path,Nunit)
      IMPLICIT NONE
!*--CERRGG61
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
      INTEGER dummyk , dummyl , i , ifst , ihi , ilo , ilst , info , j ,&
     &        m , ncycle , nt , sdim , lwork
      REAL anrm , bnrm , dif , scale , tola , tolb
!     ..
!     .. Local Arrays ..
      LOGICAL bw(NMAX) , sel(NMAX)
      INTEGER iw(LW) , idum(NMAX)
      REAL ls(NMAX) , r1(NMAX) , r2(NMAX) , rce(NMAX) , rcv(NMAX) ,     &
     &     rs(NMAX) , rw(LW)
      COMPLEX a(NMAX,NMAX) , alpha(NMAX) , b(NMAX,NMAX) , beta(NMAX) ,  &
     &        q(NMAX,NMAX) , tau(NMAX) , u(NMAX,NMAX) , v(NMAX,NMAX) ,  &
     &        w(LW) , z(NMAX,NMAX)
!     ..
!     .. External Functions ..
      LOGICAL LSAMEN
      EXTERNAL LSAMEN
!     ..
!     .. External Subroutines ..
      EXTERNAL CGGES , CGGESX , CGGEV , CGGEVX , CGGGLM , CGGHRD 
     EXTERNAL         CGGLSE , CGGQRF , CGGRQF , CHGEQZ , CTGEVC 
     EXTERNAL         CTGEXC , CTGSEN , CTGSJA , CTGSNA , CTGSYL , CUNCSD 
     EXTERNAL         CGGES3 , CGGEV3 , CGGHD3 , CGGSVD3 , CGGSVP3
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
!        CGGHRD
!
         SRNamt = 'CGGHRD'
         INFot = 1
         CALL CGGHRD('/','N',0,1,0,a,1,b,1,q,1,z,1,info)
         CALL chkxer('CGGHRD',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CGGHRD('N','/',0,1,0,a,1,b,1,q,1,z,1,info)
         CALL chkxer('CGGHRD',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CGGHRD('N','N',-1,0,0,a,1,b,1,q,1,z,1,info)
         CALL chkxer('CGGHRD',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CGGHRD('N','N',0,0,0,a,1,b,1,q,1,z,1,info)
         CALL chkxer('CGGHRD',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL CGGHRD('N','N',0,1,1,a,1,b,1,q,1,z,1,info)
         CALL chkxer('CGGHRD',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL CGGHRD('N','N',2,1,1,a,1,b,2,q,1,z,1,info)
         CALL chkxer('CGGHRD',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL CGGHRD('N','N',2,1,1,a,2,b,1,q,1,z,1,info)
         CALL chkxer('CGGHRD',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL CGGHRD('V','N',2,1,1,a,2,b,2,q,1,z,1,info)
         CALL chkxer('CGGHRD',INFot,NOUt,LERr,OK)
         INFot = 13
         CALL CGGHRD('N','V',2,1,1,a,2,b,2,q,1,z,1,info)
         CALL chkxer('CGGHRD',INFot,NOUt,LERr,OK)
         nt = nt + 9
!
!        CGGHD3
!
         SRNamt = 'CGGHD3'
         INFot = 1
         CALL CGGHD3('/','N',0,1,0,a,1,b,1,q,1,z,1,w,LW,info)
         CALL chkxer('CGGHD3',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CGGHD3('N','/',0,1,0,a,1,b,1,q,1,z,1,w,LW,info)
         CALL chkxer('CGGHD3',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CGGHD3('N','N',-1,0,0,a,1,b,1,q,1,z,1,w,LW,info)
         CALL chkxer('CGGHD3',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CGGHD3('N','N',0,0,0,a,1,b,1,q,1,z,1,w,LW,info)
         CALL chkxer('CGGHD3',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL CGGHD3('N','N',0,1,1,a,1,b,1,q,1,z,1,w,LW,info)
         CALL chkxer('CGGHD3',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL CGGHD3('N','N',2,1,1,a,1,b,2,q,1,z,1,w,LW,info)
         CALL chkxer('CGGHD3',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL CGGHD3('N','N',2,1,1,a,2,b,1,q,1,z,1,w,LW,info)
         CALL chkxer('CGGHD3',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL CGGHD3('V','N',2,1,1,a,2,b,2,q,1,z,1,w,LW,info)
         CALL chkxer('CGGHD3',INFot,NOUt,LERr,OK)
         INFot = 13
         CALL CGGHD3('N','V',2,1,1,a,2,b,2,q,1,z,1,w,LW,info)
         CALL chkxer('CGGHD3',INFot,NOUt,LERr,OK)
         nt = nt + 9
!
!        CHGEQZ
!
         SRNamt = 'CHGEQZ'
         INFot = 1
         CALL CHGEQZ('/','N','N',0,1,0,a,1,b,1,alpha,beta,q,1,z,1,w,1,  &
     &               rw,info)
         CALL chkxer('CHGEQZ',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CHGEQZ('E','/','N',0,1,0,a,1,b,1,alpha,beta,q,1,z,1,w,1,  &
     &               rw,info)
         CALL chkxer('CHGEQZ',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CHGEQZ('E','N','/',0,1,0,a,1,b,1,alpha,beta,q,1,z,1,w,1,  &
     &               rw,info)
         CALL chkxer('CHGEQZ',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CHGEQZ('E','N','N',-1,0,0,a,1,b,1,alpha,beta,q,1,z,1,w,1, &
     &               rw,info)
         CALL chkxer('CHGEQZ',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL CHGEQZ('E','N','N',0,0,0,a,1,b,1,alpha,beta,q,1,z,1,w,1,  &
     &               rw,info)
         CALL chkxer('CHGEQZ',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL CHGEQZ('E','N','N',0,1,1,a,1,b,1,alpha,beta,q,1,z,1,w,1,  &
     &               rw,info)
         CALL chkxer('CHGEQZ',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL CHGEQZ('E','N','N',2,1,1,a,1,b,2,alpha,beta,q,1,z,1,w,1,  &
     &               rw,info)
         CALL chkxer('CHGEQZ',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL CHGEQZ('E','N','N',2,1,1,a,2,b,1,alpha,beta,q,1,z,1,w,1,  &
     &               rw,info)
         CALL chkxer('CHGEQZ',INFot,NOUt,LERr,OK)
         INFot = 14
         CALL CHGEQZ('E','V','N',2,1,1,a,2,b,2,alpha,beta,q,1,z,1,w,1,  &
     &               rw,info)
         CALL chkxer('CHGEQZ',INFot,NOUt,LERr,OK)
         INFot = 16
         CALL CHGEQZ('E','N','V',2,1,1,a,2,b,2,alpha,beta,q,1,z,1,w,1,  &
     &               rw,info)
         CALL chkxer('CHGEQZ',INFot,NOUt,LERr,OK)
         nt = nt + 10
!
!        CTGEVC
!
         SRNamt = 'CTGEVC'
         INFot = 1
         CALL CTGEVC('/','A',sel,0,a,1,b,1,q,1,z,1,0,m,w,rw,info)
         CALL chkxer('CTGEVC',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CTGEVC('R','/',sel,0,a,1,b,1,q,1,z,1,0,m,w,rw,info)
         CALL chkxer('CTGEVC',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CTGEVC('R','A',sel,-1,a,1,b,1,q,1,z,1,0,m,w,rw,info)
         CALL chkxer('CTGEVC',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL CTGEVC('R','A',sel,2,a,1,b,2,q,1,z,2,0,m,w,rw,info)
         CALL chkxer('CTGEVC',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL CTGEVC('R','A',sel,2,a,2,b,1,q,1,z,2,0,m,w,rw,info)
         CALL chkxer('CTGEVC',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL CTGEVC('L','A',sel,2,a,2,b,2,q,1,z,1,0,m,w,rw,info)
         CALL chkxer('CTGEVC',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL CTGEVC('R','A',sel,2,a,2,b,2,q,1,z,1,0,m,w,rw,info)
         CALL chkxer('CTGEVC',INFot,NOUt,LERr,OK)
         INFot = 13
         CALL CTGEVC('R','A',sel,2,a,2,b,2,q,1,z,2,1,m,w,rw,info)
         CALL chkxer('CTGEVC',INFot,NOUt,LERr,OK)
         nt = nt + 8
!
!     Test error exits for the GSV path.
!
      ELSEIF ( LSAMEN(3,Path,'GSV') ) THEN
!
!        CGGSVD3
!
         SRNamt = 'CGGSVD3'
         INFot = 1
         CALL CGGSVD3('/','N','N',0,0,0,dummyk,dummyl,a,1,b,1,r1,r2,u,1,&
     &                v,1,q,1,w,lwork,rw,idum,info)
         CALL chkxer('CGGSVD3',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CGGSVD3('N','/','N',0,0,0,dummyk,dummyl,a,1,b,1,r1,r2,u,1,&
     &                v,1,q,1,w,lwork,rw,idum,info)
         CALL chkxer('CGGSVD3',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CGGSVD3('N','N','/',0,0,0,dummyk,dummyl,a,1,b,1,r1,r2,u,1,&
     &                v,1,q,1,w,lwork,rw,idum,info)
         CALL chkxer('CGGSVD3',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CGGSVD3('N','N','N',-1,0,0,dummyk,dummyl,a,1,b,1,r1,r2,u, &
     &                1,v,1,q,1,w,lwork,rw,idum,info)
         CALL chkxer('CGGSVD3',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL CGGSVD3('N','N','N',0,-1,0,dummyk,dummyl,a,1,b,1,r1,r2,u, &
     &                1,v,1,q,1,w,lwork,rw,idum,info)
         CALL chkxer('CGGSVD3',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL CGGSVD3('N','N','N',0,0,-1,dummyk,dummyl,a,1,b,1,r1,r2,u, &
     &                1,v,1,q,1,w,lwork,rw,idum,info)
         CALL chkxer('CGGSVD3',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL CGGSVD3('N','N','N',2,1,1,dummyk,dummyl,a,1,b,1,r1,r2,u,1,&
     &                v,1,q,1,w,lwork,rw,idum,info)
         CALL chkxer('CGGSVD3',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL CGGSVD3('N','N','N',1,1,2,dummyk,dummyl,a,1,b,1,r1,r2,u,1,&
     &                v,1,q,1,w,lwork,rw,idum,info)
         CALL chkxer('CGGSVD3',INFot,NOUt,LERr,OK)
         INFot = 16
         CALL CGGSVD3('U','N','N',2,2,2,dummyk,dummyl,a,2,b,2,r1,r2,u,1,&
     &                v,1,q,1,w,lwork,rw,idum,info)
         CALL chkxer('CGGSVD3',INFot,NOUt,LERr,OK)
         INFot = 18
         CALL CGGSVD3('N','V','N',2,2,2,dummyk,dummyl,a,2,b,2,r1,r2,u,2,&
     &                v,1,q,1,w,lwork,rw,idum,info)
         CALL chkxer('CGGSVD3',INFot,NOUt,LERr,OK)
         INFot = 20
         CALL CGGSVD3('N','N','Q',2,2,2,dummyk,dummyl,a,2,b,2,r1,r2,u,2,&
     &                v,2,q,1,w,lwork,rw,idum,info)
         CALL chkxer('CGGSVD3',INFot,NOUt,LERr,OK)
         nt = nt + 11
!
!        CGGSVP3
!
         SRNamt = 'CGGSVP3'
         INFot = 1
         CALL CGGSVP3('/','N','N',0,0,0,a,1,b,1,tola,tolb,dummyk,dummyl,&
     &                u,1,v,1,q,1,iw,rw,tau,w,lwork,info)
         CALL chkxer('CGGSVP3',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CGGSVP3('N','/','N',0,0,0,a,1,b,1,tola,tolb,dummyk,dummyl,&
     &                u,1,v,1,q,1,iw,rw,tau,w,lwork,info)
         CALL chkxer('CGGSVP3',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CGGSVP3('N','N','/',0,0,0,a,1,b,1,tola,tolb,dummyk,dummyl,&
     &                u,1,v,1,q,1,iw,rw,tau,w,lwork,info)
         CALL chkxer('CGGSVP3',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CGGSVP3('N','N','N',-1,0,0,a,1,b,1,tola,tolb,dummyk,      &
     &                dummyl,u,1,v,1,q,1,iw,rw,tau,w,lwork,info)
         CALL chkxer('CGGSVP3',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL CGGSVP3('N','N','N',0,-1,0,a,1,b,1,tola,tolb,dummyk,      &
     &                dummyl,u,1,v,1,q,1,iw,rw,tau,w,lwork,info)
         CALL chkxer('CGGSVP3',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL CGGSVP3('N','N','N',0,0,-1,a,1,b,1,tola,tolb,dummyk,      &
     &                dummyl,u,1,v,1,q,1,iw,rw,tau,w,lwork,info)
         CALL chkxer('CGGSVP3',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL CGGSVP3('N','N','N',2,1,1,a,1,b,1,tola,tolb,dummyk,dummyl,&
     &                u,1,v,1,q,1,iw,rw,tau,w,lwork,info)
         CALL chkxer('CGGSVP3',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL CGGSVP3('N','N','N',1,2,1,a,1,b,1,tola,tolb,dummyk,dummyl,&
     &                u,1,v,1,q,1,iw,rw,tau,w,lwork,info)
         CALL chkxer('CGGSVP3',INFot,NOUt,LERr,OK)
         INFot = 16
         CALL CGGSVP3('U','N','N',2,2,2,a,2,b,2,tola,tolb,dummyk,dummyl,&
     &                u,1,v,1,q,1,iw,rw,tau,w,lwork,info)
         CALL chkxer('CGGSVP3',INFot,NOUt,LERr,OK)
         INFot = 18
         CALL CGGSVP3('N','V','N',2,2,2,a,2,b,2,tola,tolb,dummyk,dummyl,&
     &                u,2,v,1,q,1,iw,rw,tau,w,lwork,info)
         CALL chkxer('CGGSVP3',INFot,NOUt,LERr,OK)
         INFot = 20
         CALL CGGSVP3('N','N','Q',2,2,2,a,2,b,2,tola,tolb,dummyk,dummyl,&
     &                u,2,v,2,q,1,iw,rw,tau,w,lwork,info)
         CALL chkxer('CGGSVP3',INFot,NOUt,LERr,OK)
         nt = nt + 11
!
!        CTGSJA
!
         SRNamt = 'CTGSJA'
         INFot = 1
         CALL CTGSJA('/','N','N',0,0,0,dummyk,dummyl,a,1,b,1,tola,tolb, &
     &               r1,r2,u,1,v,1,q,1,w,ncycle,info)
         CALL chkxer('CTGSJA',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CTGSJA('N','/','N',0,0,0,dummyk,dummyl,a,1,b,1,tola,tolb, &
     &               r1,r2,u,1,v,1,q,1,w,ncycle,info)
         CALL chkxer('CTGSJA',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CTGSJA('N','N','/',0,0,0,dummyk,dummyl,a,1,b,1,tola,tolb, &
     &               r1,r2,u,1,v,1,q,1,w,ncycle,info)
         CALL chkxer('CTGSJA',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CTGSJA('N','N','N',-1,0,0,dummyk,dummyl,a,1,b,1,tola,tolb,&
     &               r1,r2,u,1,v,1,q,1,w,ncycle,info)
         CALL chkxer('CTGSJA',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL CTGSJA('N','N','N',0,-1,0,dummyk,dummyl,a,1,b,1,tola,tolb,&
     &               r1,r2,u,1,v,1,q,1,w,ncycle,info)
         CALL chkxer('CTGSJA',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL CTGSJA('N','N','N',0,0,-1,dummyk,dummyl,a,1,b,1,tola,tolb,&
     &               r1,r2,u,1,v,1,q,1,w,ncycle,info)
         CALL chkxer('CTGSJA',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL CTGSJA('N','N','N',0,0,0,dummyk,dummyl,a,0,b,1,tola,tolb, &
     &               r1,r2,u,1,v,1,q,1,w,ncycle,info)
         CALL chkxer('CTGSJA',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL CTGSJA('N','N','N',0,0,0,dummyk,dummyl,a,1,b,0,tola,tolb, &
     &               r1,r2,u,1,v,1,q,1,w,ncycle,info)
         CALL chkxer('CTGSJA',INFot,NOUt,LERr,OK)
         INFot = 18
         CALL CTGSJA('U','N','N',0,0,0,dummyk,dummyl,a,1,b,1,tola,tolb, &
     &               r1,r2,u,0,v,1,q,1,w,ncycle,info)
         CALL chkxer('CTGSJA',INFot,NOUt,LERr,OK)
         INFot = 20
         CALL CTGSJA('N','V','N',0,0,0,dummyk,dummyl,a,1,b,1,tola,tolb, &
     &               r1,r2,u,1,v,0,q,1,w,ncycle,info)
         CALL chkxer('CTGSJA',INFot,NOUt,LERr,OK)
         INFot = 22
         CALL CTGSJA('N','N','Q',0,0,0,dummyk,dummyl,a,1,b,1,tola,tolb, &
     &               r1,r2,u,1,v,1,q,0,w,ncycle,info)
         CALL chkxer('CTGSJA',INFot,NOUt,LERr,OK)
         nt = nt + 11
!
!     Test error exits for the GLM path.
!
      ELSEIF ( LSAMEN(3,Path,'GLM') ) THEN
!
!        CGGGLM
!
         SRNamt = 'CGGGLM'
         INFot = 1
         CALL CGGGLM(-1,0,0,a,1,b,1,tau,alpha,beta,w,LW,info)
         CALL chkxer('CGGGLM',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CGGGLM(0,-1,0,a,1,b,1,tau,alpha,beta,w,LW,info)
         CALL chkxer('CGGGLM',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CGGGLM(0,1,0,a,1,b,1,tau,alpha,beta,w,LW,info)
         CALL chkxer('CGGGLM',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CGGGLM(0,0,-1,a,1,b,1,tau,alpha,beta,w,LW,info)
         CALL chkxer('CGGGLM',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CGGGLM(1,0,0,a,1,b,1,tau,alpha,beta,w,LW,info)
         CALL chkxer('CGGGLM',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL CGGGLM(0,0,0,a,0,b,1,tau,alpha,beta,w,LW,info)
         CALL chkxer('CGGGLM',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL CGGGLM(0,0,0,a,1,b,0,tau,alpha,beta,w,LW,info)
         CALL chkxer('CGGGLM',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL CGGGLM(1,1,1,a,1,b,1,tau,alpha,beta,w,1,info)
         CALL chkxer('CGGGLM',INFot,NOUt,LERr,OK)
         nt = nt + 8
!
!     Test error exits for the LSE path.
!
      ELSEIF ( LSAMEN(3,Path,'LSE') ) THEN
!
!        CGGLSE
!
         SRNamt = 'CGGLSE'
         INFot = 1
         CALL CGGLSE(-1,0,0,a,1,b,1,tau,alpha,beta,w,LW,info)
         CALL chkxer('CGGLSE',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CGGLSE(0,-1,0,a,1,b,1,tau,alpha,beta,w,LW,info)
         CALL chkxer('CGGLSE',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CGGLSE(0,0,-1,a,1,b,1,tau,alpha,beta,w,LW,info)
         CALL chkxer('CGGLSE',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CGGLSE(0,0,1,a,1,b,1,tau,alpha,beta,w,LW,info)
         CALL chkxer('CGGLSE',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CGGLSE(0,1,0,a,1,b,1,tau,alpha,beta,w,LW,info)
         CALL chkxer('CGGLSE',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL CGGLSE(0,0,0,a,0,b,1,tau,alpha,beta,w,LW,info)
         CALL chkxer('CGGLSE',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL CGGLSE(0,0,0,a,1,b,0,tau,alpha,beta,w,LW,info)
         CALL chkxer('CGGLSE',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL CGGLSE(1,1,1,a,1,b,1,tau,alpha,beta,w,1,info)
         CALL chkxer('CGGLSE',INFot,NOUt,LERr,OK)
         nt = nt + 8
!
!     Test error exits for the CSD path.
!
      ELSEIF ( LSAMEN(3,Path,'CSD') ) THEN
!
!        CUNCSD
!
         SRNamt = 'CUNCSD'
         INFot = 7
         CALL CUNCSD('Y','Y','Y','Y','N','N',-1,0,0,a,1,a,1,a,1,a,1,rs, &
     &               a,1,a,1,a,1,a,1,w,LW,rw,LW,iw,info)
         CALL chkxer('CUNCSD',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL CUNCSD('Y','Y','Y','Y','N','N',1,-1,0,a,1,a,1,a,1,a,1,rs, &
     &               a,1,a,1,a,1,a,1,w,LW,rw,LW,iw,info)
         CALL chkxer('CUNCSD',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL CUNCSD('Y','Y','Y','Y','N','N',1,1,-1,a,1,a,1,a,1,a,1,rs, &
     &               a,1,a,1,a,1,a,1,w,LW,rw,LW,iw,info)
         CALL chkxer('CUNCSD',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL CUNCSD('Y','Y','Y','Y','N','N',1,1,1,a,-1,a,1,a,1,a,1,rs, &
     &               a,1,a,1,a,1,a,1,w,LW,rw,LW,iw,info)
         CALL chkxer('CUNCSD',INFot,NOUt,LERr,OK)
         INFot = 20
         CALL CUNCSD('Y','Y','Y','Y','N','N',1,1,1,a,1,a,1,a,1,a,1,rs,a,&
     &               -1,a,1,a,1,a,1,w,LW,rw,LW,iw,info)
         CALL chkxer('CUNCSD',INFot,NOUt,LERr,OK)
         INFot = 22
         CALL CUNCSD('Y','Y','Y','Y','N','N',1,1,1,a,1,a,1,a,1,a,1,rs,a,&
     &               1,a,-1,a,1,a,1,w,LW,rw,LW,iw,info)
         CALL chkxer('CUNCSD',INFot,NOUt,LERr,OK)
         INFot = 24
         CALL CUNCSD('Y','Y','Y','Y','N','N',1,1,1,a,1,a,1,a,1,a,1,rs,a,&
     &               1,a,1,a,-1,a,1,w,LW,rw,LW,iw,info)
         CALL chkxer('CUNCSD',INFot,NOUt,LERr,OK)
         INFot = 26
         CALL CUNCSD('Y','Y','Y','Y','N','N',1,1,1,a,1,a,1,a,1,a,1,rs,a,&
     &               1,a,1,a,1,a,-1,w,LW,rw,LW,iw,info)
         CALL chkxer('CUNCSD',INFot,NOUt,LERr,OK)
         nt = nt + 8
!
!     Test error exits for the GQR path.
!
      ELSEIF ( LSAMEN(3,Path,'GQR') ) THEN
!
!        CGGQRF
!
         SRNamt = 'CGGQRF'
         INFot = 1
         CALL CGGQRF(-1,0,0,a,1,alpha,b,1,beta,w,LW,info)
         CALL chkxer('CGGQRF',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CGGQRF(0,-1,0,a,1,alpha,b,1,beta,w,LW,info)
         CALL chkxer('CGGQRF',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CGGQRF(0,0,-1,a,1,alpha,b,1,beta,w,LW,info)
         CALL chkxer('CGGQRF',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL CGGQRF(0,0,0,a,0,alpha,b,1,beta,w,LW,info)
         CALL chkxer('CGGQRF',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL CGGQRF(0,0,0,a,1,alpha,b,0,beta,w,LW,info)
         CALL chkxer('CGGQRF',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL CGGQRF(1,1,2,a,1,alpha,b,1,beta,w,1,info)
         CALL chkxer('CGGQRF',INFot,NOUt,LERr,OK)
         nt = nt + 6
!
!        CGGRQF
!
         SRNamt = 'CGGRQF'
         INFot = 1
         CALL CGGRQF(-1,0,0,a,1,alpha,b,1,beta,w,LW,info)
         CALL chkxer('CGGRQF',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CGGRQF(0,-1,0,a,1,alpha,b,1,beta,w,LW,info)
         CALL chkxer('CGGRQF',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CGGRQF(0,0,-1,a,1,alpha,b,1,beta,w,LW,info)
         CALL chkxer('CGGRQF',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL CGGRQF(0,0,0,a,0,alpha,b,1,beta,w,LW,info)
         CALL chkxer('CGGRQF',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL CGGRQF(0,0,0,a,1,alpha,b,0,beta,w,LW,info)
         CALL chkxer('CGGRQF',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL CGGRQF(1,1,2,a,1,alpha,b,1,beta,w,1,info)
         CALL chkxer('CGGRQF',INFot,NOUt,LERr,OK)
         nt = nt + 6
!
!     Test error exits for the CGS, CGV, CGX, and CXV paths.
!
      ELSEIF ( LSAMEN(3,Path,'CGS') .OR. LSAMEN(3,Path,'CGV') .OR.      &
     &         LSAMEN(3,Path,'CGX') .OR. LSAMEN(3,Path,'CXV') ) THEN
!
!        CGGES
!
         SRNamt = 'CGGES '
         INFot = 1
         CALL CGGES('/','N','S',clctes,1,a,1,b,1,sdim,alpha,beta,q,1,u, &
     &              1,w,1,rw,bw,info)
         CALL chkxer('CGGES ',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CGGES('N','/','S',clctes,1,a,1,b,1,sdim,alpha,beta,q,1,u, &
     &              1,w,1,rw,bw,info)
         CALL chkxer('CGGES ',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CGGES('N','V','/',clctes,1,a,1,b,1,sdim,alpha,beta,q,1,u, &
     &              1,w,1,rw,bw,info)
         CALL chkxer('CGGES ',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL CGGES('N','V','S',clctes,-1,a,1,b,1,sdim,alpha,beta,q,1,u,&
     &              1,w,1,rw,bw,info)
         CALL chkxer('CGGES ',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL CGGES('N','V','S',clctes,1,a,0,b,1,sdim,alpha,beta,q,1,u, &
     &              1,w,1,rw,bw,info)
         CALL chkxer('CGGES ',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL CGGES('N','V','S',clctes,1,a,1,b,0,sdim,alpha,beta,q,1,u, &
     &              1,w,1,rw,bw,info)
         CALL chkxer('CGGES ',INFot,NOUt,LERr,OK)
         INFot = 14
         CALL CGGES('N','V','S',clctes,1,a,1,b,1,sdim,alpha,beta,q,0,u, &
     &              1,w,1,rw,bw,info)
         CALL chkxer('CGGES ',INFot,NOUt,LERr,OK)
         INFot = 14
         CALL CGGES('V','V','S',clctes,2,a,2,b,2,sdim,alpha,beta,q,1,u, &
     &              2,w,1,rw,bw,info)
         CALL chkxer('CGGES ',INFot,NOUt,LERr,OK)
         INFot = 16
         CALL CGGES('N','V','S',clctes,1,a,1,b,1,sdim,alpha,beta,q,1,u, &
     &              0,w,1,rw,bw,info)
         CALL chkxer('CGGES ',INFot,NOUt,LERr,OK)
         INFot = 16
         CALL CGGES('V','V','S',clctes,2,a,2,b,2,sdim,alpha,beta,q,2,u, &
     &              1,w,1,rw,bw,info)
         CALL chkxer('CGGES ',INFot,NOUt,LERr,OK)
         INFot = 18
         CALL CGGES('V','V','S',clctes,2,a,2,b,2,sdim,alpha,beta,q,2,u, &
     &              2,w,1,rw,bw,info)
         CALL chkxer('CGGES ',INFot,NOUt,LERr,OK)
         nt = nt + 11
!
!        CGGES3
!
         SRNamt = 'CGGES3'
         INFot = 1
         CALL CGGES3('/','N','S',clctes,1,a,1,b,1,sdim,alpha,beta,q,1,u,&
     &               1,w,1,rw,bw,info)
         CALL chkxer('CGGES3',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CGGES3('N','/','S',clctes,1,a,1,b,1,sdim,alpha,beta,q,1,u,&
     &               1,w,1,rw,bw,info)
         CALL chkxer('CGGES3',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CGGES3('N','V','/',clctes,1,a,1,b,1,sdim,alpha,beta,q,1,u,&
     &               1,w,1,rw,bw,info)
         CALL chkxer('CGGES3',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL CGGES3('N','V','S',clctes,-1,a,1,b,1,sdim,alpha,beta,q,1, &
     &               u,1,w,1,rw,bw,info)
         CALL chkxer('CGGES3',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL CGGES3('N','V','S',clctes,1,a,0,b,1,sdim,alpha,beta,q,1,u,&
     &               1,w,1,rw,bw,info)
         CALL chkxer('CGGES3',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL CGGES3('N','V','S',clctes,1,a,1,b,0,sdim,alpha,beta,q,1,u,&
     &               1,w,1,rw,bw,info)
         CALL chkxer('CGGES3',INFot,NOUt,LERr,OK)
         INFot = 14
         CALL CGGES3('N','V','S',clctes,1,a,1,b,1,sdim,alpha,beta,q,0,u,&
     &               1,w,1,rw,bw,info)
         CALL chkxer('CGGES3',INFot,NOUt,LERr,OK)
         INFot = 14
         CALL CGGES3('V','V','S',clctes,2,a,2,b,2,sdim,alpha,beta,q,1,u,&
     &               2,w,1,rw,bw,info)
         CALL chkxer('CGGES3',INFot,NOUt,LERr,OK)
         INFot = 16
         CALL CGGES3('N','V','S',clctes,1,a,1,b,1,sdim,alpha,beta,q,1,u,&
     &               0,w,1,rw,bw,info)
         CALL chkxer('CGGES3',INFot,NOUt,LERr,OK)
         INFot = 16
         CALL CGGES3('V','V','S',clctes,2,a,2,b,2,sdim,alpha,beta,q,2,u,&
     &               1,w,1,rw,bw,info)
         CALL chkxer('CGGES3',INFot,NOUt,LERr,OK)
         INFot = 18
         CALL CGGES3('V','V','S',clctes,2,a,2,b,2,sdim,alpha,beta,q,2,u,&
     &               2,w,1,rw,bw,info)
         CALL chkxer('CGGES3',INFot,NOUt,LERr,OK)
         nt = nt + 11
!
!        CGGESX
!
         SRNamt = 'CGGESX'
         INFot = 1
         CALL CGGESX('/','N','S',clctsx,'N',1,a,1,b,1,sdim,alpha,beta,q,&
     &               1,u,1,rce,rcv,w,1,rw,iw,1,bw,info)
         CALL chkxer('CGGESX',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CGGESX('N','/','S',clctsx,'N',1,a,1,b,1,sdim,alpha,beta,q,&
     &               1,u,1,rce,rcv,w,1,rw,iw,1,bw,info)
         CALL chkxer('CGGESX',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CGGESX('V','V','/',clctsx,'N',1,a,1,b,1,sdim,alpha,beta,q,&
     &               1,u,1,rce,rcv,w,1,rw,iw,1,bw,info)
         CALL chkxer('CGGESX',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL CGGESX('V','V','S',clctsx,'/',1,a,1,b,1,sdim,alpha,beta,q,&
     &               1,u,1,rce,rcv,w,1,rw,iw,1,bw,info)
         CALL chkxer('CGGESX',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL CGGESX('V','V','S',clctsx,'B',-1,a,1,b,1,sdim,alpha,beta, &
     &               q,1,u,1,rce,rcv,w,1,rw,iw,1,bw,info)
         CALL chkxer('CGGESX',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL CGGESX('V','V','S',clctsx,'B',1,a,0,b,1,sdim,alpha,beta,q,&
     &               1,u,1,rce,rcv,w,1,rw,iw,1,bw,info)
         CALL chkxer('CGGESX',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL CGGESX('V','V','S',clctsx,'B',1,a,1,b,0,sdim,alpha,beta,q,&
     &               1,u,1,rce,rcv,w,1,rw,iw,1,bw,info)
         CALL chkxer('CGGESX',INFot,NOUt,LERr,OK)
         INFot = 15
         CALL CGGESX('V','V','S',clctsx,'B',1,a,1,b,1,sdim,alpha,beta,q,&
     &               0,u,1,rce,rcv,w,1,rw,iw,1,bw,info)
         CALL chkxer('CGGESX',INFot,NOUt,LERr,OK)
         INFot = 15
         CALL CGGESX('V','V','S',clctsx,'B',2,a,2,b,2,sdim,alpha,beta,q,&
     &               1,u,1,rce,rcv,w,1,rw,iw,1,bw,info)
         CALL chkxer('CGGESX',INFot,NOUt,LERr,OK)
         INFot = 17
         CALL CGGESX('V','V','S',clctsx,'B',1,a,1,b,1,sdim,alpha,beta,q,&
     &               1,u,0,rce,rcv,w,1,rw,iw,1,bw,info)
         CALL chkxer('CGGESX',INFot,NOUt,LERr,OK)
         INFot = 17
         CALL CGGESX('V','V','S',clctsx,'B',2,a,2,b,2,sdim,alpha,beta,q,&
     &               2,u,1,rce,rcv,w,1,rw,iw,1,bw,info)
         CALL chkxer('CGGESX',INFot,NOUt,LERr,OK)
         INFot = 21
         CALL CGGESX('V','V','S',clctsx,'B',2,a,2,b,2,sdim,alpha,beta,q,&
     &               2,u,2,rce,rcv,w,1,rw,iw,1,bw,info)
         CALL chkxer('CGGESX',INFot,NOUt,LERr,OK)
         INFot = 24
         CALL CGGESX('V','V','S',clctsx,'V',1,a,1,b,1,sdim,alpha,beta,q,&
     &               1,u,1,rce,rcv,w,32,rw,iw,0,bw,info)
         CALL chkxer('CGGESX',INFot,NOUt,LERr,OK)
         nt = nt + 13
!
!        CGGEV
!
         SRNamt = 'CGGEV '
         INFot = 1
         CALL CGGEV('/','N',1,a,1,b,1,alpha,beta,q,1,u,1,w,1,rw,info)
         CALL chkxer('CGGEV ',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CGGEV('N','/',1,a,1,b,1,alpha,beta,q,1,u,1,w,1,rw,info)
         CALL chkxer('CGGEV ',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CGGEV('V','V',-1,a,1,b,1,alpha,beta,q,1,u,1,w,1,rw,info)
         CALL chkxer('CGGEV ',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL CGGEV('V','V',1,a,0,b,1,alpha,beta,q,1,u,1,w,1,rw,info)
         CALL chkxer('CGGEV ',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL CGGEV('V','V',1,a,1,b,0,alpha,beta,q,1,u,1,w,1,rw,info)
         CALL chkxer('CGGEV ',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL CGGEV('N','V',1,a,1,b,1,alpha,beta,q,0,u,1,w,1,rw,info)
         CALL chkxer('CGGEV ',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL CGGEV('V','V',2,a,2,b,2,alpha,beta,q,1,u,2,w,1,rw,info)
         CALL chkxer('CGGEV ',INFot,NOUt,LERr,OK)
         INFot = 13
         CALL CGGEV('V','N',2,a,2,b,2,alpha,beta,q,2,u,0,w,1,rw,info)
         CALL chkxer('CGGEV ',INFot,NOUt,LERr,OK)
         INFot = 13
         CALL CGGEV('V','V',2,a,2,b,2,alpha,beta,q,2,u,1,w,1,rw,info)
         CALL chkxer('CGGEV ',INFot,NOUt,LERr,OK)
         INFot = 15
         CALL CGGEV('V','V',1,a,1,b,1,alpha,beta,q,1,u,1,w,1,rw,info)
         CALL chkxer('CGGEV ',INFot,NOUt,LERr,OK)
         nt = nt + 10
!
!        CGGEV3
!
         SRNamt = 'CGGEV3'
         INFot = 1
         CALL CGGEV3('/','N',1,a,1,b,1,alpha,beta,q,1,u,1,w,1,rw,info)
         CALL chkxer('CGGEV3',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CGGEV3('N','/',1,a,1,b,1,alpha,beta,q,1,u,1,w,1,rw,info)
         CALL chkxer('CGGEV3',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CGGEV3('V','V',-1,a,1,b,1,alpha,beta,q,1,u,1,w,1,rw,info)
         CALL chkxer('CGGEV3',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL CGGEV3('V','V',1,a,0,b,1,alpha,beta,q,1,u,1,w,1,rw,info)
         CALL chkxer('CGGEV3',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL CGGEV3('V','V',1,a,1,b,0,alpha,beta,q,1,u,1,w,1,rw,info)
         CALL chkxer('CGGEV3',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL CGGEV3('N','V',1,a,1,b,1,alpha,beta,q,0,u,1,w,1,rw,info)
         CALL chkxer('CGGEV3',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL CGGEV3('V','V',2,a,2,b,2,alpha,beta,q,1,u,2,w,1,rw,info)
         CALL chkxer('CGGEV3',INFot,NOUt,LERr,OK)
         INFot = 13
         CALL CGGEV3('V','N',2,a,2,b,2,alpha,beta,q,2,u,0,w,1,rw,info)
         CALL chkxer('CGGEV3',INFot,NOUt,LERr,OK)
         INFot = 13
         CALL CGGEV3('V','V',2,a,2,b,2,alpha,beta,q,2,u,1,w,1,rw,info)
         CALL chkxer('CGGEV3',INFot,NOUt,LERr,OK)
         INFot = 15
         CALL CGGEV3('V','V',1,a,1,b,1,alpha,beta,q,1,u,1,w,1,rw,info)
         CALL chkxer('CGGEV3',INFot,NOUt,LERr,OK)
         nt = nt + 10
!
!        CGGEVX
!
         SRNamt = 'CGGEVX'
         INFot = 1
         CALL CGGEVX('/','N','N','N',1,a,1,b,1,alpha,beta,q,1,u,1,ilo,  &
     &               ihi,ls,rs,anrm,bnrm,rce,rcv,w,1,rw,iw,bw,info)
         CALL chkxer('CGGEVX',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CGGEVX('N','/','N','N',1,a,1,b,1,alpha,beta,q,1,u,1,ilo,  &
     &               ihi,ls,rs,anrm,bnrm,rce,rcv,w,1,rw,iw,bw,info)
         CALL chkxer('CGGEVX',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CGGEVX('N','N','/','N',1,a,1,b,1,alpha,beta,q,1,u,1,ilo,  &
     &               ihi,ls,rs,anrm,bnrm,rce,rcv,w,1,rw,iw,bw,info)
         CALL chkxer('CGGEVX',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CGGEVX('N','N','N','/',1,a,1,b,1,alpha,beta,q,1,u,1,ilo,  &
     &               ihi,ls,rs,anrm,bnrm,rce,rcv,w,1,rw,iw,bw,info)
         CALL chkxer('CGGEVX',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL CGGEVX('N','N','N','N',-1,a,1,b,1,alpha,beta,q,1,u,1,ilo, &
     &               ihi,ls,rs,anrm,bnrm,rce,rcv,w,1,rw,iw,bw,info)
         CALL chkxer('CGGEVX',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL CGGEVX('N','N','N','N',1,a,0,b,1,alpha,beta,q,1,u,1,ilo,  &
     &               ihi,ls,rs,anrm,bnrm,rce,rcv,w,1,rw,iw,bw,info)
         CALL chkxer('CGGEVX',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL CGGEVX('N','N','N','N',1,a,1,b,0,alpha,beta,q,1,u,1,ilo,  &
     &               ihi,ls,rs,anrm,bnrm,rce,rcv,w,1,rw,iw,bw,info)
         CALL chkxer('CGGEVX',INFot,NOUt,LERr,OK)
         INFot = 13
         CALL CGGEVX('N','N','N','N',1,a,1,b,1,alpha,beta,q,0,u,1,ilo,  &
     &               ihi,ls,rs,anrm,bnrm,rce,rcv,w,1,rw,iw,bw,info)
         CALL chkxer('CGGEVX',INFot,NOUt,LERr,OK)
         INFot = 13
         CALL CGGEVX('N','V','N','N',2,a,2,b,2,alpha,beta,q,1,u,2,ilo,  &
     &               ihi,ls,rs,anrm,bnrm,rce,rcv,w,1,rw,iw,bw,info)
         CALL chkxer('CGGEVX',INFot,NOUt,LERr,OK)
         INFot = 15
         CALL CGGEVX('N','N','N','N',1,a,1,b,1,alpha,beta,q,1,u,0,ilo,  &
     &               ihi,ls,rs,anrm,bnrm,rce,rcv,w,1,rw,iw,bw,info)
         CALL chkxer('CGGEVX',INFot,NOUt,LERr,OK)
         INFot = 15
         CALL CGGEVX('N','N','V','N',2,a,2,b,2,alpha,beta,q,2,u,1,ilo,  &
     &               ihi,ls,rs,anrm,bnrm,rce,rcv,w,1,rw,iw,bw,info)
         CALL chkxer('CGGEVX',INFot,NOUt,LERr,OK)
         INFot = 25
         CALL CGGEVX('N','N','V','N',2,a,2,b,2,alpha,beta,q,2,u,2,ilo,  &
     &               ihi,ls,rs,anrm,bnrm,rce,rcv,w,0,rw,iw,bw,info)
         CALL chkxer('CGGEVX',INFot,NOUt,LERr,OK)
         nt = nt + 12
!
!        CTGEXC
!
         SRNamt = 'CTGEXC'
         INFot = 3
         CALL CTGEXC(.TRUE.,.TRUE.,-1,a,1,b,1,q,1,z,1,ifst,ilst,info)
         CALL chkxer('CTGEXC',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL CTGEXC(.TRUE.,.TRUE.,1,a,0,b,1,q,1,z,1,ifst,ilst,info)
         CALL chkxer('CTGEXC',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL CTGEXC(.TRUE.,.TRUE.,1,a,1,b,0,q,1,z,1,ifst,ilst,info)
         CALL chkxer('CTGEXC',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL CTGEXC(.FALSE.,.TRUE.,1,a,1,b,1,q,0,z,1,ifst,ilst,info)
         CALL chkxer('CTGEXC',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL CTGEXC(.TRUE.,.TRUE.,1,a,1,b,1,q,0,z,1,ifst,ilst,info)
         CALL chkxer('CTGEXC',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL CTGEXC(.TRUE.,.FALSE.,1,a,1,b,1,q,1,z,0,ifst,ilst,info)
         CALL chkxer('CTGEXC',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL CTGEXC(.TRUE.,.TRUE.,1,a,1,b,1,q,1,z,0,ifst,ilst,info)
         CALL chkxer('CTGEXC',INFot,NOUt,LERr,OK)
         nt = nt + 7
!
!        CTGSEN
!
         SRNamt = 'CTGSEN'
         INFot = 1
         CALL CTGSEN(-1,.TRUE.,.TRUE.,sel,1,a,1,b,1,alpha,beta,q,1,z,1, &
     &               m,tola,tolb,rcv,w,1,iw,1,info)
         CALL chkxer('CTGSEN',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL CTGSEN(1,.TRUE.,.TRUE.,sel,-1,a,1,b,1,alpha,beta,q,1,z,1, &
     &               m,tola,tolb,rcv,w,1,iw,1,info)
         CALL chkxer('CTGSEN',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL CTGSEN(1,.TRUE.,.TRUE.,sel,1,a,0,b,1,alpha,beta,q,1,z,1,m,&
     &               tola,tolb,rcv,w,1,iw,1,info)
         CALL chkxer('CTGSEN',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL CTGSEN(1,.TRUE.,.TRUE.,sel,1,a,1,b,0,alpha,beta,q,1,z,1,m,&
     &               tola,tolb,rcv,w,1,iw,1,info)
         CALL chkxer('CTGSEN',INFot,NOUt,LERr,OK)
         INFot = 13
         CALL CTGSEN(1,.TRUE.,.TRUE.,sel,1,a,1,b,1,alpha,beta,q,0,z,1,m,&
     &               tola,tolb,rcv,w,1,iw,1,info)
         CALL chkxer('CTGSEN',INFot,NOUt,LERr,OK)
         INFot = 15
         CALL CTGSEN(1,.TRUE.,.TRUE.,sel,1,a,1,b,1,alpha,beta,q,1,z,0,m,&
     &               tola,tolb,rcv,w,1,iw,1,info)
         CALL chkxer('CTGSEN',INFot,NOUt,LERr,OK)
         INFot = 21
         CALL CTGSEN(3,.TRUE.,.TRUE.,sel,1,a,1,b,1,alpha,beta,q,1,z,1,m,&
     &               tola,tolb,rcv,w,-5,iw,1,info)
         CALL chkxer('CTGSEN',INFot,NOUt,LERr,OK)
         INFot = 23
         CALL CTGSEN(0,.TRUE.,.TRUE.,sel,1,a,1,b,1,alpha,beta,q,1,z,1,m,&
     &               tola,tolb,rcv,w,20,iw,0,info)
         CALL chkxer('CTGSEN',INFot,NOUt,LERr,OK)
         INFot = 23
         CALL CTGSEN(1,.TRUE.,.TRUE.,sel,1,a,1,b,1,alpha,beta,q,1,z,1,m,&
     &               tola,tolb,rcv,w,20,iw,0,info)
         CALL chkxer('CTGSEN',INFot,NOUt,LERr,OK)
         INFot = 23
         CALL CTGSEN(5,.TRUE.,.TRUE.,sel,1,a,1,b,1,alpha,beta,q,1,z,1,m,&
     &               tola,tolb,rcv,w,20,iw,1,info)
         CALL chkxer('CTGSEN',INFot,NOUt,LERr,OK)
         nt = nt + 11
!
!        CTGSNA
!
         SRNamt = 'CTGSNA'
         INFot = 1
         CALL CTGSNA('/','A',sel,1,a,1,b,1,q,1,u,1,r1,r2,1,m,w,1,iw,    &
     &               info)
         CALL chkxer('CTGSNA',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CTGSNA('B','/',sel,1,a,1,b,1,q,1,u,1,r1,r2,1,m,w,1,iw,    &
     &               info)
         CALL chkxer('CTGSNA',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CTGSNA('B','A',sel,-1,a,1,b,1,q,1,u,1,r1,r2,1,m,w,1,iw,   &
     &               info)
         CALL chkxer('CTGSNA',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL CTGSNA('B','A',sel,1,a,0,b,1,q,1,u,1,r1,r2,1,m,w,1,iw,    &
     &               info)
         CALL chkxer('CTGSNA',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL CTGSNA('B','A',sel,1,a,1,b,0,q,1,u,1,r1,r2,1,m,w,1,iw,    &
     &               info)
         CALL chkxer('CTGSNA',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL CTGSNA('E','A',sel,1,a,1,b,1,q,0,u,1,r1,r2,1,m,w,1,iw,    &
     &               info)
         CALL chkxer('CTGSNA',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL CTGSNA('E','A',sel,1,a,1,b,1,q,1,u,0,r1,r2,1,m,w,1,iw,    &
     &               info)
         CALL chkxer('CTGSNA',INFot,NOUt,LERr,OK)
         INFot = 15
         CALL CTGSNA('E','A',sel,1,a,1,b,1,q,1,u,1,r1,r2,0,m,w,1,iw,    &
     &               info)
         CALL chkxer('CTGSNA',INFot,NOUt,LERr,OK)
         INFot = 18
         CALL CTGSNA('E','A',sel,1,a,1,b,1,q,1,u,1,r1,r2,1,m,w,0,iw,    &
     &               info)
         CALL chkxer('CTGSNA',INFot,NOUt,LERr,OK)
         nt = nt + 9
!
!        CTGSYL
!
         SRNamt = 'CTGSYL'
         INFot = 1
         CALL CTGSYL('/',0,1,1,a,1,b,1,q,1,u,1,v,1,z,1,scale,dif,w,1,iw,&
     &               info)
         CALL chkxer('CTGSYL',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL CTGSYL('N',-1,1,1,a,1,b,1,q,1,u,1,v,1,z,1,scale,dif,w,1,  &
     &               iw,info)
         CALL chkxer('CTGSYL',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL CTGSYL('N',0,0,1,a,1,b,1,q,1,u,1,v,1,z,1,scale,dif,w,1,iw,&
     &               info)
         CALL chkxer('CTGSYL',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL CTGSYL('N',0,1,0,a,1,b,1,q,1,u,1,v,1,z,1,scale,dif,w,1,iw,&
     &               info)
         CALL chkxer('CTGSYL',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL CTGSYL('N',0,1,1,a,0,b,1,q,1,u,1,v,1,z,1,scale,dif,w,1,iw,&
     &               info)
         CALL chkxer('CTGSYL',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL CTGSYL('N',0,1,1,a,1,b,0,q,1,u,1,v,1,z,1,scale,dif,w,1,iw,&
     &               info)
         CALL chkxer('CTGSYL',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL CTGSYL('N',0,1,1,a,1,b,1,q,0,u,1,v,1,z,1,scale,dif,w,1,iw,&
     &               info)
         CALL chkxer('CTGSYL',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL CTGSYL('N',0,1,1,a,1,b,1,q,1,u,0,v,1,z,1,scale,dif,w,1,iw,&
     &               info)
         CALL chkxer('CTGSYL',INFot,NOUt,LERr,OK)
         INFot = 14
         CALL CTGSYL('N',0,1,1,a,1,b,1,q,1,u,1,v,0,z,1,scale,dif,w,1,iw,&
     &               info)
         CALL chkxer('CTGSYL',INFot,NOUt,LERr,OK)
         INFot = 16
         CALL CTGSYL('N',0,1,1,a,1,b,1,q,1,u,1,v,1,z,0,scale,dif,w,1,iw,&
     &               info)
         CALL chkxer('CTGSYL',INFot,NOUt,LERr,OK)
         INFot = 20
         CALL CTGSYL('N',1,1,1,a,1,b,1,q,1,u,1,v,1,z,1,scale,dif,w,1,iw,&
     &               info)
         CALL chkxer('CTGSYL',INFot,NOUt,LERr,OK)
         INFot = 20
         CALL CTGSYL('N',2,1,1,a,1,b,1,q,1,u,1,v,1,z,1,scale,dif,w,1,iw,&
     &               info)
         CALL chkxer('CTGSYL',INFot,NOUt,LERr,OK)
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
!     End of cerrgg
!
      END SUBROUTINE cerrgg
