!*==zerrgg.f90  processed by SPAG 7.51RB at 17:37 on  4 Mar 2022
!> \brief \b ZERRGG
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZERRGG( PATH, NUNIT )
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
!> ZERRGG tests the error exits for ZGGES, ZGGESX, ZGGEV, ZGGEVX,
!> ZGGES3, ZGGEV3, ZGGGLM, ZGGHRD, ZGGLSE, ZGGQRF, ZGGRQF,
!> ZGGSVD3, ZGGSVP3, ZHGEQZ, ZTGEVC, ZTGEXC, ZTGSEN, ZTGSJA,
!> ZTGSNA, ZTGSYL, and ZUNCSD.
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
!> \ingroup complex16_eig
!
!  =====================================================================
      SUBROUTINE ZERRGG(Path,Nunit)
      IMPLICIT NONE
!*--ZERRGG61
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
      INTEGER dummyk , dummyl , i , ifst , ihi , ilo , ilst , info , j ,&
     &        m , ncycle , nt , sdim , lwork
      DOUBLE PRECISION anrm , bnrm , dif , scale , tola , tolb
!     ..
!     .. Local Arrays ..
      LOGICAL bw(NMAX) , sel(NMAX)
      INTEGER iw(LW) , idum(NMAX)
      DOUBLE PRECISION ls(NMAX) , r1(NMAX) , r2(NMAX) , rce(NMAX) ,     &
     &                 rcv(NMAX) , rs(NMAX) , rw(LW)
      COMPLEX*16 a(NMAX,NMAX) , alpha(NMAX) , b(NMAX,NMAX) , beta(NMAX) &
     &           , q(NMAX,NMAX) , tau(NMAX) , u(NMAX,NMAX) ,            &
     &           v(NMAX,NMAX) , w(LW) , z(NMAX,NMAX)
!     ..
!     .. External Functions ..
      LOGICAL LSAMEN , ZLCTES , ZLCTSX
      EXTERNAL LSAMEN , ZLCTES , ZLCTSX
!     ..
!     .. External Subroutines ..
      EXTERNAL CHKXER , ZGGES , ZGGESX , ZGGEV , ZGGEVX , ZGGGLM ,      &
     &         ZGGHRD , ZGGLSE , ZGGQRF , ZGGRQF , ZHGEQZ , ZTGEVC ,    &
     &         ZTGEXC , ZTGSEN , ZTGSJA , ZTGSNA , ZTGSYL , ZUNCSD ,    &
     &         ZGGES3 , ZGGEV3 , ZGGHD3 , ZGGSVD3 , ZGGSVP3
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
!        ZGGHRD
!
         SRNamt = 'ZGGHRD'
         INFot = 1
         CALL ZGGHRD('/','N',0,1,0,a,1,b,1,q,1,z,1,info)
         CALL CHKXER('ZGGHRD',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZGGHRD('N','/',0,1,0,a,1,b,1,q,1,z,1,info)
         CALL CHKXER('ZGGHRD',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZGGHRD('N','N',-1,0,0,a,1,b,1,q,1,z,1,info)
         CALL CHKXER('ZGGHRD',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZGGHRD('N','N',0,0,0,a,1,b,1,q,1,z,1,info)
         CALL CHKXER('ZGGHRD',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL ZGGHRD('N','N',0,1,1,a,1,b,1,q,1,z,1,info)
         CALL CHKXER('ZGGHRD',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL ZGGHRD('N','N',2,1,1,a,1,b,2,q,1,z,1,info)
         CALL CHKXER('ZGGHRD',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL ZGGHRD('N','N',2,1,1,a,2,b,1,q,1,z,1,info)
         CALL CHKXER('ZGGHRD',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL ZGGHRD('V','N',2,1,1,a,2,b,2,q,1,z,1,info)
         CALL CHKXER('ZGGHRD',INFot,NOUt,LERr,OK)
         INFot = 13
         CALL ZGGHRD('N','V',2,1,1,a,2,b,2,q,1,z,1,info)
         CALL CHKXER('ZGGHRD',INFot,NOUt,LERr,OK)
         nt = nt + 9
!
!        ZGGHD3
!
         SRNamt = 'ZGGHD3'
         INFot = 1
         CALL ZGGHD3('/','N',0,1,0,a,1,b,1,q,1,z,1,w,LW,info)
         CALL CHKXER('ZGGHD3',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZGGHD3('N','/',0,1,0,a,1,b,1,q,1,z,1,w,LW,info)
         CALL CHKXER('ZGGHD3',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZGGHD3('N','N',-1,0,0,a,1,b,1,q,1,z,1,w,LW,info)
         CALL CHKXER('ZGGHD3',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZGGHD3('N','N',0,0,0,a,1,b,1,q,1,z,1,w,LW,info)
         CALL CHKXER('ZGGHD3',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL ZGGHD3('N','N',0,1,1,a,1,b,1,q,1,z,1,w,LW,info)
         CALL CHKXER('ZGGHD3',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL ZGGHD3('N','N',2,1,1,a,1,b,2,q,1,z,1,w,LW,info)
         CALL CHKXER('ZGGHD3',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL ZGGHD3('N','N',2,1,1,a,2,b,1,q,1,z,1,w,LW,info)
         CALL CHKXER('ZGGHD3',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL ZGGHD3('V','N',2,1,1,a,2,b,2,q,1,z,1,w,LW,info)
         CALL CHKXER('ZGGHD3',INFot,NOUt,LERr,OK)
         INFot = 13
         CALL ZGGHD3('N','V',2,1,1,a,2,b,2,q,1,z,1,w,LW,info)
         CALL CHKXER('ZGGHD3',INFot,NOUt,LERr,OK)
         nt = nt + 9
!
!        ZHGEQZ
!
         SRNamt = 'ZHGEQZ'
         INFot = 1
         CALL ZHGEQZ('/','N','N',0,1,0,a,1,b,1,alpha,beta,q,1,z,1,w,1,  &
     &               rw,info)
         CALL CHKXER('ZHGEQZ',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZHGEQZ('E','/','N',0,1,0,a,1,b,1,alpha,beta,q,1,z,1,w,1,  &
     &               rw,info)
         CALL CHKXER('ZHGEQZ',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZHGEQZ('E','N','/',0,1,0,a,1,b,1,alpha,beta,q,1,z,1,w,1,  &
     &               rw,info)
         CALL CHKXER('ZHGEQZ',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZHGEQZ('E','N','N',-1,0,0,a,1,b,1,alpha,beta,q,1,z,1,w,1, &
     &               rw,info)
         CALL CHKXER('ZHGEQZ',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL ZHGEQZ('E','N','N',0,0,0,a,1,b,1,alpha,beta,q,1,z,1,w,1,  &
     &               rw,info)
         CALL CHKXER('ZHGEQZ',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL ZHGEQZ('E','N','N',0,1,1,a,1,b,1,alpha,beta,q,1,z,1,w,1,  &
     &               rw,info)
         CALL CHKXER('ZHGEQZ',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL ZHGEQZ('E','N','N',2,1,1,a,1,b,2,alpha,beta,q,1,z,1,w,1,  &
     &               rw,info)
         CALL CHKXER('ZHGEQZ',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL ZHGEQZ('E','N','N',2,1,1,a,2,b,1,alpha,beta,q,1,z,1,w,1,  &
     &               rw,info)
         CALL CHKXER('ZHGEQZ',INFot,NOUt,LERr,OK)
         INFot = 14
         CALL ZHGEQZ('E','V','N',2,1,1,a,2,b,2,alpha,beta,q,1,z,1,w,1,  &
     &               rw,info)
         CALL CHKXER('ZHGEQZ',INFot,NOUt,LERr,OK)
         INFot = 16
         CALL ZHGEQZ('E','N','V',2,1,1,a,2,b,2,alpha,beta,q,1,z,1,w,1,  &
     &               rw,info)
         CALL CHKXER('ZHGEQZ',INFot,NOUt,LERr,OK)
         nt = nt + 10
!
!        ZTGEVC
!
         SRNamt = 'ZTGEVC'
         INFot = 1
         CALL ZTGEVC('/','A',sel,0,a,1,b,1,q,1,z,1,0,m,w,rw,info)
         CALL CHKXER('ZTGEVC',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZTGEVC('R','/',sel,0,a,1,b,1,q,1,z,1,0,m,w,rw,info)
         CALL CHKXER('ZTGEVC',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZTGEVC('R','A',sel,-1,a,1,b,1,q,1,z,1,0,m,w,rw,info)
         CALL CHKXER('ZTGEVC',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL ZTGEVC('R','A',sel,2,a,1,b,2,q,1,z,2,0,m,w,rw,info)
         CALL CHKXER('ZTGEVC',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL ZTGEVC('R','A',sel,2,a,2,b,1,q,1,z,2,0,m,w,rw,info)
         CALL CHKXER('ZTGEVC',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL ZTGEVC('L','A',sel,2,a,2,b,2,q,1,z,1,0,m,w,rw,info)
         CALL CHKXER('ZTGEVC',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL ZTGEVC('R','A',sel,2,a,2,b,2,q,1,z,1,0,m,w,rw,info)
         CALL CHKXER('ZTGEVC',INFot,NOUt,LERr,OK)
         INFot = 13
         CALL ZTGEVC('R','A',sel,2,a,2,b,2,q,1,z,2,1,m,w,rw,info)
         CALL CHKXER('ZTGEVC',INFot,NOUt,LERr,OK)
         nt = nt + 8
!
!     Test error exits for the GSV path.
!
      ELSEIF ( LSAMEN(3,Path,'GSV') ) THEN
!
!        ZGGSVD3
!
         SRNamt = 'ZGGSVD3'
         INFot = 1
         CALL ZGGSVD3('/','N','N',0,0,0,dummyk,dummyl,a,1,b,1,r1,r2,u,1,&
     &                v,1,q,1,w,lwork,rw,idum,info)
         CALL CHKXER('ZGGSVD3',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZGGSVD3('N','/','N',0,0,0,dummyk,dummyl,a,1,b,1,r1,r2,u,1,&
     &                v,1,q,1,w,lwork,rw,idum,info)
         CALL CHKXER('ZGGSVD3',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZGGSVD3('N','N','/',0,0,0,dummyk,dummyl,a,1,b,1,r1,r2,u,1,&
     &                v,1,q,1,w,lwork,rw,idum,info)
         CALL CHKXER('ZGGSVD3',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZGGSVD3('N','N','N',-1,0,0,dummyk,dummyl,a,1,b,1,r1,r2,u, &
     &                1,v,1,q,1,w,lwork,rw,idum,info)
         CALL CHKXER('ZGGSVD3',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL ZGGSVD3('N','N','N',0,-1,0,dummyk,dummyl,a,1,b,1,r1,r2,u, &
     &                1,v,1,q,1,w,lwork,rw,idum,info)
         CALL CHKXER('ZGGSVD3',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL ZGGSVD3('N','N','N',0,0,-1,dummyk,dummyl,a,1,b,1,r1,r2,u, &
     &                1,v,1,q,1,w,lwork,rw,idum,info)
         CALL CHKXER('ZGGSVD3',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL ZGGSVD3('N','N','N',2,1,1,dummyk,dummyl,a,1,b,1,r1,r2,u,1,&
     &                v,1,q,1,w,lwork,rw,idum,info)
         CALL CHKXER('ZGGSVD3',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL ZGGSVD3('N','N','N',1,1,2,dummyk,dummyl,a,1,b,1,r1,r2,u,1,&
     &                v,1,q,1,w,lwork,rw,idum,info)
         CALL CHKXER('ZGGSVD3',INFot,NOUt,LERr,OK)
         INFot = 16
         CALL ZGGSVD3('U','N','N',2,2,2,dummyk,dummyl,a,2,b,2,r1,r2,u,1,&
     &                v,1,q,1,w,lwork,rw,idum,info)
         CALL CHKXER('ZGGSVD3',INFot,NOUt,LERr,OK)
         INFot = 18
         CALL ZGGSVD3('N','V','N',2,2,2,dummyk,dummyl,a,2,b,2,r1,r2,u,2,&
     &                v,1,q,1,w,lwork,rw,idum,info)
         CALL CHKXER('ZGGSVD3',INFot,NOUt,LERr,OK)
         INFot = 20
         CALL ZGGSVD3('N','N','Q',2,2,2,dummyk,dummyl,a,2,b,2,r1,r2,u,2,&
     &                v,2,q,1,w,lwork,rw,idum,info)
         CALL CHKXER('ZGGSVD3',INFot,NOUt,LERr,OK)
         nt = nt + 11
!
!        ZGGSVP3
!
         SRNamt = 'ZGGSVP3'
         INFot = 1
         CALL ZGGSVP3('/','N','N',0,0,0,a,1,b,1,tola,tolb,dummyk,dummyl,&
     &                u,1,v,1,q,1,iw,rw,tau,w,lwork,info)
         CALL CHKXER('ZGGSVP3',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZGGSVP3('N','/','N',0,0,0,a,1,b,1,tola,tolb,dummyk,dummyl,&
     &                u,1,v,1,q,1,iw,rw,tau,w,lwork,info)
         CALL CHKXER('ZGGSVP3',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZGGSVP3('N','N','/',0,0,0,a,1,b,1,tola,tolb,dummyk,dummyl,&
     &                u,1,v,1,q,1,iw,rw,tau,w,lwork,info)
         CALL CHKXER('ZGGSVP3',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZGGSVP3('N','N','N',-1,0,0,a,1,b,1,tola,tolb,dummyk,      &
     &                dummyl,u,1,v,1,q,1,iw,rw,tau,w,lwork,info)
         CALL CHKXER('ZGGSVP3',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL ZGGSVP3('N','N','N',0,-1,0,a,1,b,1,tola,tolb,dummyk,      &
     &                dummyl,u,1,v,1,q,1,iw,rw,tau,w,lwork,info)
         CALL CHKXER('ZGGSVP3',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL ZGGSVP3('N','N','N',0,0,-1,a,1,b,1,tola,tolb,dummyk,      &
     &                dummyl,u,1,v,1,q,1,iw,rw,tau,w,lwork,info)
         CALL CHKXER('ZGGSVP3',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL ZGGSVP3('N','N','N',2,1,1,a,1,b,1,tola,tolb,dummyk,dummyl,&
     &                u,1,v,1,q,1,iw,rw,tau,w,lwork,info)
         CALL CHKXER('ZGGSVP3',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL ZGGSVP3('N','N','N',1,2,1,a,1,b,1,tola,tolb,dummyk,dummyl,&
     &                u,1,v,1,q,1,iw,rw,tau,w,lwork,info)
         CALL CHKXER('ZGGSVP3',INFot,NOUt,LERr,OK)
         INFot = 16
         CALL ZGGSVP3('U','N','N',2,2,2,a,2,b,2,tola,tolb,dummyk,dummyl,&
     &                u,1,v,1,q,1,iw,rw,tau,w,lwork,info)
         CALL CHKXER('ZGGSVP3',INFot,NOUt,LERr,OK)
         INFot = 18
         CALL ZGGSVP3('N','V','N',2,2,2,a,2,b,2,tola,tolb,dummyk,dummyl,&
     &                u,2,v,1,q,1,iw,rw,tau,w,lwork,info)
         CALL CHKXER('ZGGSVP3',INFot,NOUt,LERr,OK)
         INFot = 20
         CALL ZGGSVP3('N','N','Q',2,2,2,a,2,b,2,tola,tolb,dummyk,dummyl,&
     &                u,2,v,2,q,1,iw,rw,tau,w,lwork,info)
         CALL CHKXER('ZGGSVP3',INFot,NOUt,LERr,OK)
         nt = nt + 11
!
!        ZTGSJA
!
         SRNamt = 'ZTGSJA'
         INFot = 1
         CALL ZTGSJA('/','N','N',0,0,0,dummyk,dummyl,a,1,b,1,tola,tolb, &
     &               r1,r2,u,1,v,1,q,1,w,ncycle,info)
         CALL CHKXER('ZTGSJA',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZTGSJA('N','/','N',0,0,0,dummyk,dummyl,a,1,b,1,tola,tolb, &
     &               r1,r2,u,1,v,1,q,1,w,ncycle,info)
         CALL CHKXER('ZTGSJA',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZTGSJA('N','N','/',0,0,0,dummyk,dummyl,a,1,b,1,tola,tolb, &
     &               r1,r2,u,1,v,1,q,1,w,ncycle,info)
         CALL CHKXER('ZTGSJA',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZTGSJA('N','N','N',-1,0,0,dummyk,dummyl,a,1,b,1,tola,tolb,&
     &               r1,r2,u,1,v,1,q,1,w,ncycle,info)
         CALL CHKXER('ZTGSJA',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL ZTGSJA('N','N','N',0,-1,0,dummyk,dummyl,a,1,b,1,tola,tolb,&
     &               r1,r2,u,1,v,1,q,1,w,ncycle,info)
         CALL CHKXER('ZTGSJA',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL ZTGSJA('N','N','N',0,0,-1,dummyk,dummyl,a,1,b,1,tola,tolb,&
     &               r1,r2,u,1,v,1,q,1,w,ncycle,info)
         CALL CHKXER('ZTGSJA',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL ZTGSJA('N','N','N',0,0,0,dummyk,dummyl,a,0,b,1,tola,tolb, &
     &               r1,r2,u,1,v,1,q,1,w,ncycle,info)
         CALL CHKXER('ZTGSJA',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL ZTGSJA('N','N','N',0,0,0,dummyk,dummyl,a,1,b,0,tola,tolb, &
     &               r1,r2,u,1,v,1,q,1,w,ncycle,info)
         CALL CHKXER('ZTGSJA',INFot,NOUt,LERr,OK)
         INFot = 18
         CALL ZTGSJA('U','N','N',0,0,0,dummyk,dummyl,a,1,b,1,tola,tolb, &
     &               r1,r2,u,0,v,1,q,1,w,ncycle,info)
         CALL CHKXER('ZTGSJA',INFot,NOUt,LERr,OK)
         INFot = 20
         CALL ZTGSJA('N','V','N',0,0,0,dummyk,dummyl,a,1,b,1,tola,tolb, &
     &               r1,r2,u,1,v,0,q,1,w,ncycle,info)
         CALL CHKXER('ZTGSJA',INFot,NOUt,LERr,OK)
         INFot = 22
         CALL ZTGSJA('N','N','Q',0,0,0,dummyk,dummyl,a,1,b,1,tola,tolb, &
     &               r1,r2,u,1,v,1,q,0,w,ncycle,info)
         CALL CHKXER('ZTGSJA',INFot,NOUt,LERr,OK)
         nt = nt + 11
!
!     Test error exits for the GLM path.
!
      ELSEIF ( LSAMEN(3,Path,'GLM') ) THEN
!
!        ZGGGLM
!
         SRNamt = 'ZGGGLM'
         INFot = 1
         CALL ZGGGLM(-1,0,0,a,1,b,1,tau,alpha,beta,w,LW,info)
         CALL CHKXER('ZGGGLM',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZGGGLM(0,-1,0,a,1,b,1,tau,alpha,beta,w,LW,info)
         CALL CHKXER('ZGGGLM',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZGGGLM(0,1,0,a,1,b,1,tau,alpha,beta,w,LW,info)
         CALL CHKXER('ZGGGLM',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZGGGLM(0,0,-1,a,1,b,1,tau,alpha,beta,w,LW,info)
         CALL CHKXER('ZGGGLM',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZGGGLM(1,0,0,a,1,b,1,tau,alpha,beta,w,LW,info)
         CALL CHKXER('ZGGGLM',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL ZGGGLM(0,0,0,a,0,b,1,tau,alpha,beta,w,LW,info)
         CALL CHKXER('ZGGGLM',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL ZGGGLM(0,0,0,a,1,b,0,tau,alpha,beta,w,LW,info)
         CALL CHKXER('ZGGGLM',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL ZGGGLM(1,1,1,a,1,b,1,tau,alpha,beta,w,1,info)
         CALL CHKXER('ZGGGLM',INFot,NOUt,LERr,OK)
         nt = nt + 8
!
!     Test error exits for the LSE path.
!
      ELSEIF ( LSAMEN(3,Path,'LSE') ) THEN
!
!        ZGGLSE
!
         SRNamt = 'ZGGLSE'
         INFot = 1
         CALL ZGGLSE(-1,0,0,a,1,b,1,tau,alpha,beta,w,LW,info)
         CALL CHKXER('ZGGLSE',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZGGLSE(0,-1,0,a,1,b,1,tau,alpha,beta,w,LW,info)
         CALL CHKXER('ZGGLSE',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZGGLSE(0,0,-1,a,1,b,1,tau,alpha,beta,w,LW,info)
         CALL CHKXER('ZGGLSE',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZGGLSE(0,0,1,a,1,b,1,tau,alpha,beta,w,LW,info)
         CALL CHKXER('ZGGLSE',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZGGLSE(0,1,0,a,1,b,1,tau,alpha,beta,w,LW,info)
         CALL CHKXER('ZGGLSE',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL ZGGLSE(0,0,0,a,0,b,1,tau,alpha,beta,w,LW,info)
         CALL CHKXER('ZGGLSE',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL ZGGLSE(0,0,0,a,1,b,0,tau,alpha,beta,w,LW,info)
         CALL CHKXER('ZGGLSE',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL ZGGLSE(1,1,1,a,1,b,1,tau,alpha,beta,w,1,info)
         CALL CHKXER('ZGGLSE',INFot,NOUt,LERr,OK)
         nt = nt + 8
!
!     Test error exits for the CSD path.
!
      ELSEIF ( LSAMEN(3,Path,'CSD') ) THEN
!
!        ZUNCSD
!
         SRNamt = 'ZUNCSD'
         INFot = 7
         CALL ZUNCSD('Y','Y','Y','Y','N','N',-1,0,0,a,1,a,1,a,1,a,1,rs, &
     &               a,1,a,1,a,1,a,1,w,LW,rw,LW,iw,info)
         CALL CHKXER('ZUNCSD',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL ZUNCSD('Y','Y','Y','Y','N','N',1,-1,0,a,1,a,1,a,1,a,1,rs, &
     &               a,1,a,1,a,1,a,1,w,LW,rw,LW,iw,info)
         CALL CHKXER('ZUNCSD',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL ZUNCSD('Y','Y','Y','Y','N','N',1,1,-1,a,1,a,1,a,1,a,1,rs, &
     &               a,1,a,1,a,1,a,1,w,LW,rw,LW,iw,info)
         CALL CHKXER('ZUNCSD',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL ZUNCSD('Y','Y','Y','Y','N','N',1,1,1,a,-1,a,1,a,1,a,1,rs, &
     &               a,1,a,1,a,1,a,1,w,LW,rw,LW,iw,info)
         CALL CHKXER('ZUNCSD',INFot,NOUt,LERr,OK)
         INFot = 20
         CALL ZUNCSD('Y','Y','Y','Y','N','N',1,1,1,a,1,a,1,a,1,a,1,rs,a,&
     &               -1,a,1,a,1,a,1,w,LW,rw,LW,iw,info)
         CALL CHKXER('ZUNCSD',INFot,NOUt,LERr,OK)
         INFot = 22
         CALL ZUNCSD('Y','Y','Y','Y','N','N',1,1,1,a,1,a,1,a,1,a,1,rs,a,&
     &               1,a,-1,a,1,a,1,w,LW,rw,LW,iw,info)
         CALL CHKXER('ZUNCSD',INFot,NOUt,LERr,OK)
         INFot = 24
         CALL ZUNCSD('Y','Y','Y','Y','N','N',1,1,1,a,1,a,1,a,1,a,1,rs,a,&
     &               1,a,1,a,-1,a,1,w,LW,rw,LW,iw,info)
         CALL CHKXER('ZUNCSD',INFot,NOUt,LERr,OK)
         INFot = 26
         CALL ZUNCSD('Y','Y','Y','Y','N','N',1,1,1,a,1,a,1,a,1,a,1,rs,a,&
     &               1,a,1,a,1,a,-1,w,LW,rw,LW,iw,info)
         CALL CHKXER('ZUNCSD',INFot,NOUt,LERr,OK)
         nt = nt + 8
!
!     Test error exits for the GQR path.
!
      ELSEIF ( LSAMEN(3,Path,'GQR') ) THEN
!
!        ZGGQRF
!
         SRNamt = 'ZGGQRF'
         INFot = 1
         CALL ZGGQRF(-1,0,0,a,1,alpha,b,1,beta,w,LW,info)
         CALL CHKXER('ZGGQRF',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZGGQRF(0,-1,0,a,1,alpha,b,1,beta,w,LW,info)
         CALL CHKXER('ZGGQRF',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZGGQRF(0,0,-1,a,1,alpha,b,1,beta,w,LW,info)
         CALL CHKXER('ZGGQRF',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL ZGGQRF(0,0,0,a,0,alpha,b,1,beta,w,LW,info)
         CALL CHKXER('ZGGQRF',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL ZGGQRF(0,0,0,a,1,alpha,b,0,beta,w,LW,info)
         CALL CHKXER('ZGGQRF',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL ZGGQRF(1,1,2,a,1,alpha,b,1,beta,w,1,info)
         CALL CHKXER('ZGGQRF',INFot,NOUt,LERr,OK)
         nt = nt + 6
!
!        ZGGRQF
!
         SRNamt = 'ZGGRQF'
         INFot = 1
         CALL ZGGRQF(-1,0,0,a,1,alpha,b,1,beta,w,LW,info)
         CALL CHKXER('ZGGRQF',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZGGRQF(0,-1,0,a,1,alpha,b,1,beta,w,LW,info)
         CALL CHKXER('ZGGRQF',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZGGRQF(0,0,-1,a,1,alpha,b,1,beta,w,LW,info)
         CALL CHKXER('ZGGRQF',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL ZGGRQF(0,0,0,a,0,alpha,b,1,beta,w,LW,info)
         CALL CHKXER('ZGGRQF',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL ZGGRQF(0,0,0,a,1,alpha,b,0,beta,w,LW,info)
         CALL CHKXER('ZGGRQF',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL ZGGRQF(1,1,2,a,1,alpha,b,1,beta,w,1,info)
         CALL CHKXER('ZGGRQF',INFot,NOUt,LERr,OK)
         nt = nt + 6
!
!     Test error exits for the ZGS, ZGV, ZGX, and ZXV paths.
!
      ELSEIF ( LSAMEN(3,Path,'ZGS') .OR. LSAMEN(3,Path,'ZGV') .OR.      &
     &         LSAMEN(3,Path,'ZGX') .OR. LSAMEN(3,Path,'ZXV') ) THEN
!
!        ZGGES
!
         SRNamt = 'ZGGES '
         INFot = 1
         CALL ZGGES('/','N','S',ZLCTES,1,a,1,b,1,sdim,alpha,beta,q,1,u, &
     &              1,w,1,rw,bw,info)
         CALL CHKXER('ZGGES ',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZGGES('N','/','S',ZLCTES,1,a,1,b,1,sdim,alpha,beta,q,1,u, &
     &              1,w,1,rw,bw,info)
         CALL CHKXER('ZGGES ',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZGGES('N','V','/',ZLCTES,1,a,1,b,1,sdim,alpha,beta,q,1,u, &
     &              1,w,1,rw,bw,info)
         CALL CHKXER('ZGGES ',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL ZGGES('N','V','S',ZLCTES,-1,a,1,b,1,sdim,alpha,beta,q,1,u,&
     &              1,w,1,rw,bw,info)
         CALL CHKXER('ZGGES ',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL ZGGES('N','V','S',ZLCTES,1,a,0,b,1,sdim,alpha,beta,q,1,u, &
     &              1,w,1,rw,bw,info)
         CALL CHKXER('ZGGES ',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL ZGGES('N','V','S',ZLCTES,1,a,1,b,0,sdim,alpha,beta,q,1,u, &
     &              1,w,1,rw,bw,info)
         CALL CHKXER('ZGGES ',INFot,NOUt,LERr,OK)
         INFot = 14
         CALL ZGGES('N','V','S',ZLCTES,1,a,1,b,1,sdim,alpha,beta,q,0,u, &
     &              1,w,1,rw,bw,info)
         CALL CHKXER('ZGGES ',INFot,NOUt,LERr,OK)
         INFot = 14
         CALL ZGGES('V','V','S',ZLCTES,2,a,2,b,2,sdim,alpha,beta,q,1,u, &
     &              2,w,1,rw,bw,info)
         CALL CHKXER('ZGGES ',INFot,NOUt,LERr,OK)
         INFot = 16
         CALL ZGGES('N','V','S',ZLCTES,1,a,1,b,1,sdim,alpha,beta,q,1,u, &
     &              0,w,1,rw,bw,info)
         CALL CHKXER('ZGGES ',INFot,NOUt,LERr,OK)
         INFot = 16
         CALL ZGGES('V','V','S',ZLCTES,2,a,2,b,2,sdim,alpha,beta,q,2,u, &
     &              1,w,1,rw,bw,info)
         CALL CHKXER('ZGGES ',INFot,NOUt,LERr,OK)
         INFot = 18
         CALL ZGGES('V','V','S',ZLCTES,2,a,2,b,2,sdim,alpha,beta,q,2,u, &
     &              2,w,1,rw,bw,info)
         CALL CHKXER('ZGGES ',INFot,NOUt,LERr,OK)
         nt = nt + 11
!
!        ZGGES3
!
         SRNamt = 'ZGGES3'
         INFot = 1
         CALL ZGGES3('/','N','S',ZLCTES,1,a,1,b,1,sdim,alpha,beta,q,1,u,&
     &               1,w,1,rw,bw,info)
         CALL CHKXER('ZGGES3',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZGGES3('N','/','S',ZLCTES,1,a,1,b,1,sdim,alpha,beta,q,1,u,&
     &               1,w,1,rw,bw,info)
         CALL CHKXER('ZGGES3',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZGGES3('N','V','/',ZLCTES,1,a,1,b,1,sdim,alpha,beta,q,1,u,&
     &               1,w,1,rw,bw,info)
         CALL CHKXER('ZGGES3',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL ZGGES3('N','V','S',ZLCTES,-1,a,1,b,1,sdim,alpha,beta,q,1, &
     &               u,1,w,1,rw,bw,info)
         CALL CHKXER('ZGGES3',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL ZGGES3('N','V','S',ZLCTES,1,a,0,b,1,sdim,alpha,beta,q,1,u,&
     &               1,w,1,rw,bw,info)
         CALL CHKXER('ZGGES3',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL ZGGES3('N','V','S',ZLCTES,1,a,1,b,0,sdim,alpha,beta,q,1,u,&
     &               1,w,1,rw,bw,info)
         CALL CHKXER('ZGGES3',INFot,NOUt,LERr,OK)
         INFot = 14
         CALL ZGGES3('N','V','S',ZLCTES,1,a,1,b,1,sdim,alpha,beta,q,0,u,&
     &               1,w,1,rw,bw,info)
         CALL CHKXER('ZGGES3',INFot,NOUt,LERr,OK)
         INFot = 14
         CALL ZGGES3('V','V','S',ZLCTES,2,a,2,b,2,sdim,alpha,beta,q,1,u,&
     &               2,w,1,rw,bw,info)
         CALL CHKXER('ZGGES3',INFot,NOUt,LERr,OK)
         INFot = 16
         CALL ZGGES3('N','V','S',ZLCTES,1,a,1,b,1,sdim,alpha,beta,q,1,u,&
     &               0,w,1,rw,bw,info)
         CALL CHKXER('ZGGES3',INFot,NOUt,LERr,OK)
         INFot = 16
         CALL ZGGES3('V','V','S',ZLCTES,2,a,2,b,2,sdim,alpha,beta,q,2,u,&
     &               1,w,1,rw,bw,info)
         CALL CHKXER('ZGGES3',INFot,NOUt,LERr,OK)
         INFot = 18
         CALL ZGGES3('V','V','S',ZLCTES,2,a,2,b,2,sdim,alpha,beta,q,2,u,&
     &               2,w,1,rw,bw,info)
         CALL CHKXER('ZGGES3',INFot,NOUt,LERr,OK)
         nt = nt + 11
!
!        ZGGESX
!
         SRNamt = 'ZGGESX'
         INFot = 1
         CALL ZGGESX('/','N','S',ZLCTSX,'N',1,a,1,b,1,sdim,alpha,beta,q,&
     &               1,u,1,rce,rcv,w,1,rw,iw,1,bw,info)
         CALL CHKXER('ZGGESX',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZGGESX('N','/','S',ZLCTSX,'N',1,a,1,b,1,sdim,alpha,beta,q,&
     &               1,u,1,rce,rcv,w,1,rw,iw,1,bw,info)
         CALL CHKXER('ZGGESX',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZGGESX('V','V','/',ZLCTSX,'N',1,a,1,b,1,sdim,alpha,beta,q,&
     &               1,u,1,rce,rcv,w,1,rw,iw,1,bw,info)
         CALL CHKXER('ZGGESX',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL ZGGESX('V','V','S',ZLCTSX,'/',1,a,1,b,1,sdim,alpha,beta,q,&
     &               1,u,1,rce,rcv,w,1,rw,iw,1,bw,info)
         CALL CHKXER('ZGGESX',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL ZGGESX('V','V','S',ZLCTSX,'B',-1,a,1,b,1,sdim,alpha,beta, &
     &               q,1,u,1,rce,rcv,w,1,rw,iw,1,bw,info)
         CALL CHKXER('ZGGESX',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL ZGGESX('V','V','S',ZLCTSX,'B',1,a,0,b,1,sdim,alpha,beta,q,&
     &               1,u,1,rce,rcv,w,1,rw,iw,1,bw,info)
         CALL CHKXER('ZGGESX',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL ZGGESX('V','V','S',ZLCTSX,'B',1,a,1,b,0,sdim,alpha,beta,q,&
     &               1,u,1,rce,rcv,w,1,rw,iw,1,bw,info)
         CALL CHKXER('ZGGESX',INFot,NOUt,LERr,OK)
         INFot = 15
         CALL ZGGESX('V','V','S',ZLCTSX,'B',1,a,1,b,1,sdim,alpha,beta,q,&
     &               0,u,1,rce,rcv,w,1,rw,iw,1,bw,info)
         CALL CHKXER('ZGGESX',INFot,NOUt,LERr,OK)
         INFot = 15
         CALL ZGGESX('V','V','S',ZLCTSX,'B',2,a,2,b,2,sdim,alpha,beta,q,&
     &               1,u,1,rce,rcv,w,1,rw,iw,1,bw,info)
         CALL CHKXER('ZGGESX',INFot,NOUt,LERr,OK)
         INFot = 17
         CALL ZGGESX('V','V','S',ZLCTSX,'B',1,a,1,b,1,sdim,alpha,beta,q,&
     &               1,u,0,rce,rcv,w,1,rw,iw,1,bw,info)
         CALL CHKXER('ZGGESX',INFot,NOUt,LERr,OK)
         INFot = 17
         CALL ZGGESX('V','V','S',ZLCTSX,'B',2,a,2,b,2,sdim,alpha,beta,q,&
     &               2,u,1,rce,rcv,w,1,rw,iw,1,bw,info)
         CALL CHKXER('ZGGESX',INFot,NOUt,LERr,OK)
         INFot = 21
         CALL ZGGESX('V','V','S',ZLCTSX,'B',2,a,2,b,2,sdim,alpha,beta,q,&
     &               2,u,2,rce,rcv,w,1,rw,iw,1,bw,info)
         CALL CHKXER('ZGGESX',INFot,NOUt,LERr,OK)
         INFot = 24
         CALL ZGGESX('V','V','S',ZLCTSX,'V',1,a,1,b,1,sdim,alpha,beta,q,&
     &               1,u,1,rce,rcv,w,32,rw,iw,0,bw,info)
         CALL CHKXER('ZGGESX',INFot,NOUt,LERr,OK)
         nt = nt + 13
!
!        ZGGEV
!
         SRNamt = 'ZGGEV '
         INFot = 1
         CALL ZGGEV('/','N',1,a,1,b,1,alpha,beta,q,1,u,1,w,1,rw,info)
         CALL CHKXER('ZGGEV ',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZGGEV('N','/',1,a,1,b,1,alpha,beta,q,1,u,1,w,1,rw,info)
         CALL CHKXER('ZGGEV ',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZGGEV('V','V',-1,a,1,b,1,alpha,beta,q,1,u,1,w,1,rw,info)
         CALL CHKXER('ZGGEV ',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL ZGGEV('V','V',1,a,0,b,1,alpha,beta,q,1,u,1,w,1,rw,info)
         CALL CHKXER('ZGGEV ',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL ZGGEV('V','V',1,a,1,b,0,alpha,beta,q,1,u,1,w,1,rw,info)
         CALL CHKXER('ZGGEV ',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL ZGGEV('N','V',1,a,1,b,1,alpha,beta,q,0,u,1,w,1,rw,info)
         CALL CHKXER('ZGGEV ',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL ZGGEV('V','V',2,a,2,b,2,alpha,beta,q,1,u,2,w,1,rw,info)
         CALL CHKXER('ZGGEV ',INFot,NOUt,LERr,OK)
         INFot = 13
         CALL ZGGEV('V','N',2,a,2,b,2,alpha,beta,q,2,u,0,w,1,rw,info)
         CALL CHKXER('ZGGEV ',INFot,NOUt,LERr,OK)
         INFot = 13
         CALL ZGGEV('V','V',2,a,2,b,2,alpha,beta,q,2,u,1,w,1,rw,info)
         CALL CHKXER('ZGGEV ',INFot,NOUt,LERr,OK)
         INFot = 15
         CALL ZGGEV('V','V',1,a,1,b,1,alpha,beta,q,1,u,1,w,1,rw,info)
         CALL CHKXER('ZGGEV ',INFot,NOUt,LERr,OK)
         nt = nt + 10
!
!        ZGGEV3
!
         SRNamt = 'ZGGEV3'
         INFot = 1
         CALL ZGGEV3('/','N',1,a,1,b,1,alpha,beta,q,1,u,1,w,1,rw,info)
         CALL CHKXER('ZGGEV3',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZGGEV3('N','/',1,a,1,b,1,alpha,beta,q,1,u,1,w,1,rw,info)
         CALL CHKXER('ZGGEV3',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZGGEV3('V','V',-1,a,1,b,1,alpha,beta,q,1,u,1,w,1,rw,info)
         CALL CHKXER('ZGGEV3',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL ZGGEV3('V','V',1,a,0,b,1,alpha,beta,q,1,u,1,w,1,rw,info)
         CALL CHKXER('ZGGEV3',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL ZGGEV3('V','V',1,a,1,b,0,alpha,beta,q,1,u,1,w,1,rw,info)
         CALL CHKXER('ZGGEV3',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL ZGGEV3('N','V',1,a,1,b,1,alpha,beta,q,0,u,1,w,1,rw,info)
         CALL CHKXER('ZGGEV3',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL ZGGEV3('V','V',2,a,2,b,2,alpha,beta,q,1,u,2,w,1,rw,info)
         CALL CHKXER('ZGGEV3',INFot,NOUt,LERr,OK)
         INFot = 13
         CALL ZGGEV3('V','N',2,a,2,b,2,alpha,beta,q,2,u,0,w,1,rw,info)
         CALL CHKXER('ZGGEV3',INFot,NOUt,LERr,OK)
         INFot = 13
         CALL ZGGEV3('V','V',2,a,2,b,2,alpha,beta,q,2,u,1,w,1,rw,info)
         CALL CHKXER('ZGGEV3',INFot,NOUt,LERr,OK)
         INFot = 15
         CALL ZGGEV3('V','V',1,a,1,b,1,alpha,beta,q,1,u,1,w,1,rw,info)
         CALL CHKXER('ZGGEV3',INFot,NOUt,LERr,OK)
         nt = nt + 10
!
!        ZGGEVX
!
         SRNamt = 'ZGGEVX'
         INFot = 1
         CALL ZGGEVX('/','N','N','N',1,a,1,b,1,alpha,beta,q,1,u,1,ilo,  &
     &               ihi,ls,rs,anrm,bnrm,rce,rcv,w,1,rw,iw,bw,info)
         CALL CHKXER('ZGGEVX',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZGGEVX('N','/','N','N',1,a,1,b,1,alpha,beta,q,1,u,1,ilo,  &
     &               ihi,ls,rs,anrm,bnrm,rce,rcv,w,1,rw,iw,bw,info)
         CALL CHKXER('ZGGEVX',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZGGEVX('N','N','/','N',1,a,1,b,1,alpha,beta,q,1,u,1,ilo,  &
     &               ihi,ls,rs,anrm,bnrm,rce,rcv,w,1,rw,iw,bw,info)
         CALL CHKXER('ZGGEVX',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZGGEVX('N','N','N','/',1,a,1,b,1,alpha,beta,q,1,u,1,ilo,  &
     &               ihi,ls,rs,anrm,bnrm,rce,rcv,w,1,rw,iw,bw,info)
         CALL CHKXER('ZGGEVX',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL ZGGEVX('N','N','N','N',-1,a,1,b,1,alpha,beta,q,1,u,1,ilo, &
     &               ihi,ls,rs,anrm,bnrm,rce,rcv,w,1,rw,iw,bw,info)
         CALL CHKXER('ZGGEVX',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL ZGGEVX('N','N','N','N',1,a,0,b,1,alpha,beta,q,1,u,1,ilo,  &
     &               ihi,ls,rs,anrm,bnrm,rce,rcv,w,1,rw,iw,bw,info)
         CALL CHKXER('ZGGEVX',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL ZGGEVX('N','N','N','N',1,a,1,b,0,alpha,beta,q,1,u,1,ilo,  &
     &               ihi,ls,rs,anrm,bnrm,rce,rcv,w,1,rw,iw,bw,info)
         CALL CHKXER('ZGGEVX',INFot,NOUt,LERr,OK)
         INFot = 13
         CALL ZGGEVX('N','N','N','N',1,a,1,b,1,alpha,beta,q,0,u,1,ilo,  &
     &               ihi,ls,rs,anrm,bnrm,rce,rcv,w,1,rw,iw,bw,info)
         CALL CHKXER('ZGGEVX',INFot,NOUt,LERr,OK)
         INFot = 13
         CALL ZGGEVX('N','V','N','N',2,a,2,b,2,alpha,beta,q,1,u,2,ilo,  &
     &               ihi,ls,rs,anrm,bnrm,rce,rcv,w,1,rw,iw,bw,info)
         CALL CHKXER('ZGGEVX',INFot,NOUt,LERr,OK)
         INFot = 15
         CALL ZGGEVX('N','N','N','N',1,a,1,b,1,alpha,beta,q,1,u,0,ilo,  &
     &               ihi,ls,rs,anrm,bnrm,rce,rcv,w,1,rw,iw,bw,info)
         CALL CHKXER('ZGGEVX',INFot,NOUt,LERr,OK)
         INFot = 15
         CALL ZGGEVX('N','N','V','N',2,a,2,b,2,alpha,beta,q,2,u,1,ilo,  &
     &               ihi,ls,rs,anrm,bnrm,rce,rcv,w,1,rw,iw,bw,info)
         CALL CHKXER('ZGGEVX',INFot,NOUt,LERr,OK)
         INFot = 25
         CALL ZGGEVX('N','N','V','N',2,a,2,b,2,alpha,beta,q,2,u,2,ilo,  &
     &               ihi,ls,rs,anrm,bnrm,rce,rcv,w,0,rw,iw,bw,info)
         CALL CHKXER('ZGGEVX',INFot,NOUt,LERr,OK)
         nt = nt + 12
!
!        ZTGEXC
!
         SRNamt = 'ZTGEXC'
         INFot = 3
         CALL ZTGEXC(.TRUE.,.TRUE.,-1,a,1,b,1,q,1,z,1,ifst,ilst,info)
         CALL CHKXER('ZTGEXC',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL ZTGEXC(.TRUE.,.TRUE.,1,a,0,b,1,q,1,z,1,ifst,ilst,info)
         CALL CHKXER('ZTGEXC',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL ZTGEXC(.TRUE.,.TRUE.,1,a,1,b,0,q,1,z,1,ifst,ilst,info)
         CALL CHKXER('ZTGEXC',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL ZTGEXC(.FALSE.,.TRUE.,1,a,1,b,1,q,0,z,1,ifst,ilst,info)
         CALL CHKXER('ZTGEXC',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL ZTGEXC(.TRUE.,.TRUE.,1,a,1,b,1,q,0,z,1,ifst,ilst,info)
         CALL CHKXER('ZTGEXC',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL ZTGEXC(.TRUE.,.FALSE.,1,a,1,b,1,q,1,z,0,ifst,ilst,info)
         CALL CHKXER('ZTGEXC',INFot,NOUt,LERr,OK)
         INFot = 11
         CALL ZTGEXC(.TRUE.,.TRUE.,1,a,1,b,1,q,1,z,0,ifst,ilst,info)
         CALL CHKXER('ZTGEXC',INFot,NOUt,LERr,OK)
         nt = nt + 7
!
!        ZTGSEN
!
         SRNamt = 'ZTGSEN'
         INFot = 1
         CALL ZTGSEN(-1,.TRUE.,.TRUE.,sel,1,a,1,b,1,alpha,beta,q,1,z,1, &
     &               m,tola,tolb,rcv,w,1,iw,1,info)
         CALL CHKXER('ZTGSEN',INFot,NOUt,LERr,OK)
         INFot = 5
         CALL ZTGSEN(1,.TRUE.,.TRUE.,sel,-1,a,1,b,1,alpha,beta,q,1,z,1, &
     &               m,tola,tolb,rcv,w,1,iw,1,info)
         CALL CHKXER('ZTGSEN',INFot,NOUt,LERr,OK)
         INFot = 7
         CALL ZTGSEN(1,.TRUE.,.TRUE.,sel,1,a,0,b,1,alpha,beta,q,1,z,1,m,&
     &               tola,tolb,rcv,w,1,iw,1,info)
         CALL CHKXER('ZTGSEN',INFot,NOUt,LERr,OK)
         INFot = 9
         CALL ZTGSEN(1,.TRUE.,.TRUE.,sel,1,a,1,b,0,alpha,beta,q,1,z,1,m,&
     &               tola,tolb,rcv,w,1,iw,1,info)
         CALL CHKXER('ZTGSEN',INFot,NOUt,LERr,OK)
         INFot = 13
         CALL ZTGSEN(1,.TRUE.,.TRUE.,sel,1,a,1,b,1,alpha,beta,q,0,z,1,m,&
     &               tola,tolb,rcv,w,1,iw,1,info)
         CALL CHKXER('ZTGSEN',INFot,NOUt,LERr,OK)
         INFot = 15
         CALL ZTGSEN(1,.TRUE.,.TRUE.,sel,1,a,1,b,1,alpha,beta,q,1,z,0,m,&
     &               tola,tolb,rcv,w,1,iw,1,info)
         CALL CHKXER('ZTGSEN',INFot,NOUt,LERr,OK)
         INFot = 21
         CALL ZTGSEN(3,.TRUE.,.TRUE.,sel,1,a,1,b,1,alpha,beta,q,1,z,1,m,&
     &               tola,tolb,rcv,w,-5,iw,1,info)
         CALL CHKXER('ZTGSEN',INFot,NOUt,LERr,OK)
         INFot = 23
         CALL ZTGSEN(0,.TRUE.,.TRUE.,sel,1,a,1,b,1,alpha,beta,q,1,z,1,m,&
     &               tola,tolb,rcv,w,20,iw,0,info)
         CALL CHKXER('ZTGSEN',INFot,NOUt,LERr,OK)
         INFot = 23
         CALL ZTGSEN(1,.TRUE.,.TRUE.,sel,1,a,1,b,1,alpha,beta,q,1,z,1,m,&
     &               tola,tolb,rcv,w,20,iw,0,info)
         CALL CHKXER('ZTGSEN',INFot,NOUt,LERr,OK)
         INFot = 23
         CALL ZTGSEN(5,.TRUE.,.TRUE.,sel,1,a,1,b,1,alpha,beta,q,1,z,1,m,&
     &               tola,tolb,rcv,w,20,iw,1,info)
         CALL CHKXER('ZTGSEN',INFot,NOUt,LERr,OK)
         nt = nt + 11
!
!        ZTGSNA
!
         SRNamt = 'ZTGSNA'
         INFot = 1
         CALL ZTGSNA('/','A',sel,1,a,1,b,1,q,1,u,1,r1,r2,1,m,w,1,iw,    &
     &               info)
         CALL CHKXER('ZTGSNA',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZTGSNA('B','/',sel,1,a,1,b,1,q,1,u,1,r1,r2,1,m,w,1,iw,    &
     &               info)
         CALL CHKXER('ZTGSNA',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZTGSNA('B','A',sel,-1,a,1,b,1,q,1,u,1,r1,r2,1,m,w,1,iw,   &
     &               info)
         CALL CHKXER('ZTGSNA',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL ZTGSNA('B','A',sel,1,a,0,b,1,q,1,u,1,r1,r2,1,m,w,1,iw,    &
     &               info)
         CALL CHKXER('ZTGSNA',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL ZTGSNA('B','A',sel,1,a,1,b,0,q,1,u,1,r1,r2,1,m,w,1,iw,    &
     &               info)
         CALL CHKXER('ZTGSNA',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL ZTGSNA('E','A',sel,1,a,1,b,1,q,0,u,1,r1,r2,1,m,w,1,iw,    &
     &               info)
         CALL CHKXER('ZTGSNA',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL ZTGSNA('E','A',sel,1,a,1,b,1,q,1,u,0,r1,r2,1,m,w,1,iw,    &
     &               info)
         CALL CHKXER('ZTGSNA',INFot,NOUt,LERr,OK)
         INFot = 15
         CALL ZTGSNA('E','A',sel,1,a,1,b,1,q,1,u,1,r1,r2,0,m,w,1,iw,    &
     &               info)
         CALL CHKXER('ZTGSNA',INFot,NOUt,LERr,OK)
         INFot = 18
         CALL ZTGSNA('E','A',sel,1,a,1,b,1,q,1,u,1,r1,r2,1,m,w,0,iw,    &
     &               info)
         CALL CHKXER('ZTGSNA',INFot,NOUt,LERr,OK)
         nt = nt + 9
!
!        ZTGSYL
!
         SRNamt = 'ZTGSYL'
         INFot = 1
         CALL ZTGSYL('/',0,1,1,a,1,b,1,q,1,u,1,v,1,z,1,scale,dif,w,1,iw,&
     &               info)
         CALL CHKXER('ZTGSYL',INFot,NOUt,LERr,OK)
         INFot = 2
         CALL ZTGSYL('N',-1,1,1,a,1,b,1,q,1,u,1,v,1,z,1,scale,dif,w,1,  &
     &               iw,info)
         CALL CHKXER('ZTGSYL',INFot,NOUt,LERr,OK)
         INFot = 3
         CALL ZTGSYL('N',0,0,1,a,1,b,1,q,1,u,1,v,1,z,1,scale,dif,w,1,iw,&
     &               info)
         CALL CHKXER('ZTGSYL',INFot,NOUt,LERr,OK)
         INFot = 4
         CALL ZTGSYL('N',0,1,0,a,1,b,1,q,1,u,1,v,1,z,1,scale,dif,w,1,iw,&
     &               info)
         CALL CHKXER('ZTGSYL',INFot,NOUt,LERr,OK)
         INFot = 6
         CALL ZTGSYL('N',0,1,1,a,0,b,1,q,1,u,1,v,1,z,1,scale,dif,w,1,iw,&
     &               info)
         CALL CHKXER('ZTGSYL',INFot,NOUt,LERr,OK)
         INFot = 8
         CALL ZTGSYL('N',0,1,1,a,1,b,0,q,1,u,1,v,1,z,1,scale,dif,w,1,iw,&
     &               info)
         CALL CHKXER('ZTGSYL',INFot,NOUt,LERr,OK)
         INFot = 10
         CALL ZTGSYL('N',0,1,1,a,1,b,1,q,0,u,1,v,1,z,1,scale,dif,w,1,iw,&
     &               info)
         CALL CHKXER('ZTGSYL',INFot,NOUt,LERr,OK)
         INFot = 12
         CALL ZTGSYL('N',0,1,1,a,1,b,1,q,1,u,0,v,1,z,1,scale,dif,w,1,iw,&
     &               info)
         CALL CHKXER('ZTGSYL',INFot,NOUt,LERr,OK)
         INFot = 14
         CALL ZTGSYL('N',0,1,1,a,1,b,1,q,1,u,1,v,0,z,1,scale,dif,w,1,iw,&
     &               info)
         CALL CHKXER('ZTGSYL',INFot,NOUt,LERr,OK)
         INFot = 16
         CALL ZTGSYL('N',0,1,1,a,1,b,1,q,1,u,1,v,1,z,0,scale,dif,w,1,iw,&
     &               info)
         CALL CHKXER('ZTGSYL',INFot,NOUt,LERr,OK)
         INFot = 20
         CALL ZTGSYL('N',1,1,1,a,1,b,1,q,1,u,1,v,1,z,1,scale,dif,w,1,iw,&
     &               info)
         CALL CHKXER('ZTGSYL',INFot,NOUt,LERr,OK)
         INFot = 20
         CALL ZTGSYL('N',2,1,1,a,1,b,1,q,1,u,1,v,1,z,1,scale,dif,w,1,iw,&
     &               info)
         CALL CHKXER('ZTGSYL',INFot,NOUt,LERr,OK)
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
!     End of ZERRGG
!
      END SUBROUTINE ZERRGG
