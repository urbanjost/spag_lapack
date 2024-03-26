!*==ztsqr01.f90  processed by SPAG 7.51RB at 17:35 on  4 Mar 2022
!> \brief \b ztsqr01
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ztsqr01(TSSW, M,N, MB, NB, RESULT)
!
!       .. Scalar Arguments ..
!       INTEGER M, N, MB
!       .. Return values ..
!       DOUBLE PRECISION RESULT(6)
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ztsqr01 tests ZGEQR , ZGELQ, ZGEMLQ and ZGEMQR.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] TSSW
!> \verbatim
!>          TSSW is CHARACTER
!>          'TS' for testing tall skinny QR
!>               and anything else for testing short wide LQ
!> \endverbatim
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          Number of rows in test matrix.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          Number of columns in test matrix.
!> \endverbatim
!> \param[in] MB
!> \verbatim
!>          MB is INTEGER
!>          Number of row in row block in test matrix.
!> \endverbatim
!>
!> \param[in] NB
!> \verbatim
!>          NB is INTEGER
!>          Number of columns in column block test matrix.
!> \endverbatim
!>
!> \param[out] RESULT
!> \verbatim
!>          RESULT is DOUBLE PRECISION array, dimension (6)
!>          Results of each of the six tests below.
!>
!>          RESULT(1) = | A - Q R | or | A - L Q |
!>          RESULT(2) = | I - Q^H Q | or | I - Q Q^H |
!>          RESULT(3) = | Q C - Q C |
!>          RESULT(4) = | Q^H C - Q^H C |
!>          RESULT(5) = | C Q - C Q |
!>          RESULT(6) = | C Q^H - C Q^H |
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
!> \date April 2012
!
!  =====================================================================
      SUBROUTINE ztsqr01(Tssw,M,N,Mb,Nb,Result)
      IMPLICIT NONE
!*--ZTSQR0186
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     April 2012
!
!     .. Scalar Arguments ..
      CHARACTER Tssw
      INTEGER M , N , Mb , Nb
!     .. Return values ..
      DOUBLE PRECISION Result(6)
!
!  =====================================================================
!
!     ..
!     .. Local allocatable arrays
      COMPLEX*16 , ALLOCATABLE  ::  af(:,:) , q(:,:) , r(:,:) &
     &                              , work(:) , t(:) , cf(:,:) , df(:,:)&
     &                              , a(:,:) , c(:,:) , d(:,:) , lq(:,:)
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO
      double precision, allocatable :: dpwork(:)
      COMPLEX*16 ONE , CZERO
      PARAMETER (ZERO=0.0,ONE=(1.0,0.0),CZERO=(0.0,0.0))
!     ..
!     .. Local Scalars ..
      LOGICAL testzeros , ts
      INTEGER info , j , k , l , lwork , tsize , mnb
      DOUBLE PRECISION anorm , eps , resid , cnorm , dnorm
!     ..
!     .. Local Arrays ..
      INTEGER iseed(4)
      COMPLEX*16 tquery(5) , workquery(1)
!     ..
!     .. External Functions ..
      DOUBLE PRECISION DLAMCH , ZLANGE , ZLANSY
      LOGICAL LSAME
      INTEGER ilaenv
EXTERNAL    DLAMCH,   ZLANGE,   ZLANSY,   LSAME
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX , MIN
!     .. Scalars in Common ..
      CHARACTER*32 SRNamt
!     ..
!     .. Common blocks ..
      COMMON /SRNAMC/ SRNamt
!     ..
!     .. Data statements ..
      DATA iseed/1988 , 1989 , 1990 , 1991/
!
!     TEST TALL SKINNY OR SHORT WIDE
!
      ts = LSAME(Tssw,'TS')
!
!     TEST MATRICES WITH HALF OF MATRIX BEING ZEROS
!
      testzeros = .FALSE.
!
      eps = DLAMCH('Epsilon')
      k = MIN(M,N)
      l = MAX(M,N,1)
      mnb = MAX(Mb,Nb)
      lwork = MAX(3,l)*mnb
!
!     Dynamically allocate local arrays
!
      ALLOCATE (a(M,N),af(M,N),q(l,l),r(M,l),dpwork(l),c(M,N),cf(M,N),   &
     &          d(N,M),df(N,M),lq(l,N))
!
!     Put random numbers into A and copy to AF
!
      DO j = 1 , N
         CALL ZLARNV(2,iseed,M,a(1,j))
      ENDDO
      IF ( testzeros ) THEN
         IF ( M>=4 ) THEN
            DO j = 1 , N
               CALL ZLARNV(2,iseed,M/2,a(M/4,j))
            ENDDO
         ENDIF
      ENDIF
      CALL ZLACPY('Full',M,N,a,M,af,M)
!
      IF ( ts ) THEN
!
!     Factor the matrix A in the array AF.
!
         CALL ZGEQR(M,N,af,M,tquery,-1,workquery,-1,info)
         tsize = INT(tquery(1))
         lwork = INT(workquery(1))
         CALL ZGEMQR('L','N',M,M,k,af,M,tquery,tsize,cf,M,workquery,-1, &
     &               info)
         lwork = MAX(lwork,INT(workquery(1)))
         CALL ZGEMQR('L','N',M,N,k,af,M,tquery,tsize,cf,M,workquery,-1, &
     &               info)
         lwork = MAX(lwork,INT(workquery(1)))
         CALL ZGEMQR('L','C',M,N,k,af,M,tquery,tsize,cf,M,workquery,-1, &
     &               info)
         lwork = MAX(lwork,INT(workquery(1)))
         CALL ZGEMQR('R','N',N,M,k,af,M,tquery,tsize,df,N,workquery,-1, &
     &               info)
         lwork = MAX(lwork,INT(workquery(1)))
         CALL ZGEMQR('R','C',N,M,k,af,M,tquery,tsize,df,N,workquery,-1, &
     &               info)
         lwork = MAX(lwork,INT(workquery(1)))
         ALLOCATE (t(tsize))
         ALLOCATE (work(lwork))
         SRNamt = 'ZGEQR'
         CALL ZGEQR(M,N,af,M,t,tsize,work,lwork,info)
!
!     Generate the m-by-m matrix Q
!
         CALL ZLASET('Full',M,M,CZERO,ONE,q,M)
         SRNamt = 'ZGEMQR'
         CALL ZGEMQR('L','N',M,M,k,af,M,t,tsize,q,M,work,lwork,info)
!
!     Copy R
!
         CALL ZLASET('Full',M,N,CZERO,CZERO,r,M)
         CALL ZLACPY('Upper',M,N,af,M,r,M)
!
!     Compute |R - Q'*A| / |A| and store in RESULT(1)
!
         CALL ZGEMM('C','N',M,N,M,-ONE,q,M,a,M,ONE,r,M)
!================================================================
!      CHARACTER :: Norm
!      INTEGER :: M
!      INTEGER , INTENT(IN) :: N
!      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
!      INTEGER , INTENT(IN) :: Lda
!      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
!================================================================
         anorm = ZLANGE('1',M,N,a,M,G_DPNULL)
         resid = ZLANGE('1',M,N,r,M,G_DPNULL)
         IF ( anorm>ZERO ) THEN
            Result(1) = resid/(eps*MAX(1,M)*anorm)
         ELSE
            Result(1) = ZERO
         ENDIF
!
!     Compute |I - Q'*Q| and store in RESULT(2)
!
         CALL ZLASET('Full',M,M,CZERO,ONE,r,M)
         CALL ZHERK('U','C',M,M,DREAL(-ONE),q,M,DREAL(ONE),r,M)
         resid = ZLANSY('1','Upper',M,r,M,dpwork)
         Result(2) = resid/(eps*MAX(1,M))
!
!     Generate random m-by-n matrix C and a copy CF
!
         DO j = 1 , N
            CALL ZLARNV(2,iseed,M,c(1,j))
         ENDDO
         cnorm = ZLANGE('1',M,N,c,M,G_DPNULL)
         CALL ZLACPY('Full',M,N,c,M,cf,M)
!
!     Apply Q to C as Q*C
!
         SRNamt = 'ZGEMQR'
         CALL ZGEMQR('L','N',M,N,k,af,M,t,tsize,cf,M,work,lwork,info)
!
!     Compute |Q*C - Q*C| / |C|
!
         CALL ZGEMM('N','N',M,N,M,-ONE,q,M,c,M,ONE,cf,M)
         resid = ZLANGE('1',M,N,cf,M,G_DPNULL)
         IF ( cnorm>ZERO ) THEN
            Result(3) = resid/(eps*MAX(1,M)*cnorm)
         ELSE
            Result(3) = ZERO
         ENDIF
!
!     Copy C into CF again
!
         CALL ZLACPY('Full',M,N,c,M,cf,M)
!
!     Apply Q to C as QT*C
!
         SRNamt = 'ZGEMQR'
         CALL ZGEMQR('L','C',M,N,k,af,M,t,tsize,cf,M,work,lwork,info)
!
!     Compute |QT*C - QT*C| / |C|
!
         CALL ZGEMM('C','N',M,N,M,-ONE,q,M,c,M,ONE,cf,M)
         resid = ZLANGE('1',M,N,cf,M,G_DPNULL)
         IF ( cnorm>ZERO ) THEN
            Result(4) = resid/(eps*MAX(1,M)*cnorm)
         ELSE
            Result(4) = ZERO
         ENDIF
!
!     Generate random n-by-m matrix D and a copy DF
!
         DO j = 1 , M
            CALL ZLARNV(2,iseed,N,d(1,j))
         ENDDO
         dnorm = ZLANGE('1',N,M,d,N,G_DPNULL)
         CALL ZLACPY('Full',N,M,d,N,df,N)
!
!     Apply Q to D as D*Q
!
         SRNamt = 'ZGEMQR'
         CALL ZGEMQR('R','N',N,M,k,af,M,t,tsize,df,N,work,lwork,info)
!
!     Compute |D*Q - D*Q| / |D|
!
         CALL ZGEMM('N','N',N,M,M,-ONE,d,N,q,M,ONE,df,N)
         resid = ZLANGE('1',N,M,df,N,G_DPNULL)
         IF ( dnorm>ZERO ) THEN
            Result(5) = resid/(eps*MAX(1,M)*dnorm)
         ELSE
            Result(5) = ZERO
         ENDIF
!
!     Copy D into DF again
!
         CALL ZLACPY('Full',N,M,d,N,df,N)
!
!     Apply Q to D as D*QT
!
         CALL ZGEMQR('R','C',N,M,k,af,M,t,tsize,df,N,work,lwork,info)
!
!     Compute |D*QT - D*QT| / |D|
!
         CALL ZGEMM('N','C',N,M,M,-ONE,d,N,q,M,ONE,df,N)
         resid = ZLANGE('1',N,M,df,N,G_DPNULL)
         IF ( cnorm>ZERO ) THEN
            Result(6) = resid/(eps*MAX(1,M)*dnorm)
         ELSE
            Result(6) = ZERO
         ENDIF
!
!     Short and wide
!
      ELSE
         CALL ZGELQ(M,N,af,M,tquery,-1,workquery,-1,info)
         tsize = INT(tquery(1))
         lwork = INT(workquery(1))
         CALL ZGEMLQ('R','N',N,N,k,af,M,tquery,tsize,q,N,workquery,-1,  &
     &               info)
         lwork = MAX(lwork,INT(workquery(1)))
         CALL ZGEMLQ('L','N',N,M,k,af,M,tquery,tsize,df,N,workquery,-1, &
     &               info)
         lwork = MAX(lwork,INT(workquery(1)))
         CALL ZGEMLQ('L','C',N,M,k,af,M,tquery,tsize,df,N,workquery,-1, &
     &               info)
         lwork = MAX(lwork,INT(workquery(1)))
         CALL ZGEMLQ('R','N',M,N,k,af,M,tquery,tsize,cf,M,workquery,-1, &
     &               info)
         lwork = MAX(lwork,INT(workquery(1)))
         CALL ZGEMLQ('R','C',M,N,k,af,M,tquery,tsize,cf,M,workquery,-1, &
     &               info)
         lwork = MAX(lwork,INT(workquery(1)))
         ALLOCATE (t(tsize))
         ALLOCATE (work(lwork))
         SRNamt = 'ZGELQ'
         CALL ZGELQ(M,N,af,M,t,tsize,work,lwork,info)
!
!
!     Generate the n-by-n matrix Q
!
         CALL ZLASET('Full',N,N,CZERO,ONE,q,N)
         SRNamt = 'ZGEMLQ'
         CALL ZGEMLQ('R','N',N,N,k,af,M,t,tsize,q,N,work,lwork,info)
!
!     Copy R
!
         CALL ZLASET('Full',M,N,CZERO,CZERO,lq,l)
         CALL ZLACPY('Lower',M,N,af,M,lq,l)
!
!     Compute |L - A*Q'| / |A| and store in RESULT(1)
!
         CALL ZGEMM('N','C',M,N,N,-ONE,a,M,q,N,ONE,lq,l)
         anorm = ZLANGE('1',M,N,a,M,G_DPNULL)
         resid = ZLANGE('1',M,N,lq,l,G_DPNULL)
         IF ( anorm>ZERO ) THEN
            Result(1) = resid/(eps*MAX(1,N)*anorm)
         ELSE
            Result(1) = ZERO
         ENDIF
!
!     Compute |I - Q'*Q| and store in RESULT(2)
!
         CALL ZLASET('Full',N,N,CZERO,ONE,lq,l)
         CALL ZHERK('U','C',N,N,DREAL(-ONE),q,N,DREAL(ONE),lq,l)
         resid = ZLANSY('1','Upper',N,lq,l,dpwork)
         Result(2) = resid/(eps*MAX(1,N))
!
!     Generate random m-by-n matrix C and a copy CF
!
         DO j = 1 , M
            CALL ZLARNV(2,iseed,N,d(1,j))
         ENDDO
         dnorm = ZLANGE('1',N,M,d,N,G_DPNULL)
         CALL ZLACPY('Full',N,M,d,N,df,N)
!
!     Apply Q to C as Q*C
!
         CALL ZGEMLQ('L','N',N,M,k,af,M,t,tsize,df,N,work,lwork,info)
!
!     Compute |Q*D - Q*D| / |D|
!
         CALL ZGEMM('N','N',N,M,N,-ONE,q,N,d,N,ONE,df,N)
         resid = ZLANGE('1',N,M,df,N,G_DPNULL)
         IF ( dnorm>ZERO ) THEN
            Result(3) = resid/(eps*MAX(1,N)*dnorm)
         ELSE
            Result(3) = ZERO
         ENDIF
!
!     Copy D into DF again
!
         CALL ZLACPY('Full',N,M,d,N,df,N)
!
!     Apply Q to D as QT*D
!
         CALL ZGEMLQ('L','C',N,M,k,af,M,t,tsize,df,N,work,lwork,info)
!
!     Compute |QT*D - QT*D| / |D|
!
         CALL ZGEMM('C','N',N,M,N,-ONE,q,N,d,N,ONE,df,N)
         resid = ZLANGE('1',N,M,df,N,G_DPNULL)
         IF ( dnorm>ZERO ) THEN
            Result(4) = resid/(eps*MAX(1,N)*dnorm)
         ELSE
            Result(4) = ZERO
         ENDIF
!
!     Generate random n-by-m matrix D and a copy DF
!
         DO j = 1 , N
            CALL ZLARNV(2,iseed,M,c(1,j))
         ENDDO
         cnorm = ZLANGE('1',M,N,c,M,G_DPNULL)
         CALL ZLACPY('Full',M,N,c,M,cf,M)
!
!     Apply Q to C as C*Q
!
         CALL ZGEMLQ('R','N',M,N,k,af,M,t,tsize,cf,M,work,lwork,info)
!
!     Compute |C*Q - C*Q| / |C|
!
         CALL ZGEMM('N','N',M,N,N,-ONE,c,M,q,N,ONE,cf,M)
         resid = ZLANGE('1',N,M,df,N,G_DPNULL)
         IF ( cnorm>ZERO ) THEN
            Result(5) = resid/(eps*MAX(1,N)*cnorm)
         ELSE
            Result(5) = ZERO
         ENDIF
!
!     Copy C into CF again
!
         CALL ZLACPY('Full',M,N,c,M,cf,M)
!
!     Apply Q to D as D*QT
!
         CALL ZGEMLQ('R','C',M,N,k,af,M,t,tsize,cf,M,work,lwork,info)
!
!     Compute |C*QT - C*QT| / |C|
!
         CALL ZGEMM('N','C',M,N,N,-ONE,c,M,q,N,ONE,cf,M)
         resid = ZLANGE('1',M,N,cf,M,G_DPNULL)
         IF ( cnorm>ZERO ) THEN
            Result(6) = resid/(eps*MAX(1,N)*cnorm)
         ELSE
            Result(6) = ZERO
         ENDIF
!
      ENDIF
!
!     Deallocate all arrays
!
      DEALLOCATE (a,af,q,r,dpwork,work,t,c,d,cf,df)
!
      END SUBROUTINE ztsqr01
