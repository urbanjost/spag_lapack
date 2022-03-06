!*==zlqt05.f90  processed by SPAG 7.51RB at 17:35 on  4 Mar 2022
 
!> \brief \b zlqt05
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE zlqt05(M,N,L,NB,RESULT)
!
!       .. Scalar Arguments ..
!       INTEGER LWORK, M, N, L, NB, LDT
!       .. Return values ..
!       DOUBLE PRECISION RESULT(6)
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> zqrt05 tests ZTPLQT and ZTPMLQT.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          Number of rows in lower part of the test matrix.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          Number of columns in test matrix.
!> \endverbatim
!>
!> \param[in] L
!> \verbatim
!>          L is INTEGER
!>          The number of rows of the upper trapezoidal part the
!>          lower test matrix.  0 <= L <= M.
!> \endverbatim
!>
!> \param[in] NB
!> \verbatim
!>          NB is INTEGER
!>          Block size of test matrix.  NB <= N.
!> \endverbatim
!>
!> \param[out] RESULT
!> \verbatim
!>          RESULT is DOUBLE PRECISION array, dimension (6)
!>          Results of each of the six tests below.
!>
!>          RESULT(1) = | A - Q R |
!>          RESULT(2) = | I - Q^H Q |
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
!> \ingroup double_lin
!
!  =====================================================================
      SUBROUTINE zlqt05(M,N,L,Nb,Result)
      IMPLICIT NONE
!*--ZLQT0585
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     April 2012
!
!     .. Scalar Arguments ..
      INTEGER lwork , M , N , L , Nb , ldt
!     .. Return values ..
      DOUBLE PRECISION Result(6)
!
!  =====================================================================
!
!     ..
!     .. Local allocatable arrays
      COMPLEX*16 , ALLOCATABLE  ::  af(:,:) , q(:,:) , r(:,:) , rwork(:)&
     &                              , work(:) , t(:,:) , cf(:,:) ,      &
     &                              df(:,:) , a(:,:) , c(:,:) , d(:,:)
     double precision, allocatable :: dbwork(:)
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO
      COMPLEX*16 ONE , CZERO
      PARAMETER (ZERO=0.0,ONE=(1.0,0.0),CZERO=(0.0,0.0))
!     ..
!     .. Local Scalars ..
      INTEGER info , j , k , n2 , np1
      DOUBLE PRECISION anorm , eps , resid , cnorm , dnorm
!     ..
!     .. Local Arrays ..
      INTEGER iseed(4)
!     ..
!     .. External Functions ..
      DOUBLE PRECISION DLAMCH
      DOUBLE PRECISION ZLANGE , ZLANSY
      LOGICAL LSAME
EXTERNAL    DLAMCH,   ZLANGE,   ZLANSY,   LSAME
!     ..
!     .. Data statements ..
      DATA iseed/1988 , 1989 , 1990 , 1991/
!
      eps = DLAMCH('Epsilon')
      k = M
      n2 = M + N
      IF ( N>0 ) THEN
         np1 = M + 1
      ELSE
         np1 = 1
      ENDIF
      lwork = n2*n2*Nb
!
!     Dynamically allocate all arrays
!
      ALLOCATE (a(M,n2),af(M,n2),q(n2,n2),r(n2,n2),rwork(n2),dbwork(n2),work(lwork)&
     &          ,t(Nb,M),c(n2,M),cf(n2,M),d(M,n2),df(M,n2))
!
!     Put random stuff into A
!
      ldt = Nb
      CALL ZLASET('Full',M,n2,CZERO,CZERO,a,M)
      CALL ZLASET('Full',Nb,M,CZERO,CZERO,t,Nb)
      DO j = 1 , M
         CALL ZLARNV(2,iseed,M-j+1,a(j,j))
      ENDDO
      IF ( N>0 ) THEN
         DO j = 1 , N - L
            CALL ZLARNV(2,iseed,M,a(1,MIN(N+M,M+1)+j-1))
         ENDDO
      ENDIF
      IF ( L>0 ) THEN
         DO j = 1 , L
            CALL ZLARNV(2,iseed,M-j+1,a(j,MIN(N+M,N+M-L+1)+j-1))
         ENDDO
      ENDIF
!
!     Copy the matrix A to the array AF.
!
      CALL ZLACPY('Full',M,n2,a,M,af,M)
!
!     Factor the matrix A in the array AF.
!
      CALL ZTPLQT(M,N,L,Nb,af,M,af(1,np1),M,t,ldt,work,info)
!
!     Generate the (M+N)-by-(M+N) matrix Q by applying H to I
!
      CALL ZLASET('Full',n2,n2,CZERO,ONE,q,n2)
      CALL ZGEMLQT('L','N',n2,n2,k,Nb,af,M,t,ldt,q,n2,work,info)
!
!     Copy L
!
      CALL ZLASET('Full',n2,n2,CZERO,CZERO,r,n2)
      CALL ZLACPY('Lower',M,n2,af,M,r,n2)
!
!     Compute |L - A*Q*C| / |A| and store in RESULT(1)
!
      CALL ZGEMM('N','C',M,n2,n2,-ONE,a,M,q,n2,ONE,r,n2)
      anorm = ZLANGE('1',M,n2,a,M,dbwork)
      resid = ZLANGE('1',M,n2,r,n2,dbwork)
      IF ( anorm>ZERO ) THEN
         Result(1) = resid/(eps*anorm*MAX(1,n2))
      ELSE
         Result(1) = ZERO
      ENDIF
!
!     Compute |I - Q*Q'| and store in RESULT(2)
!
      CALL ZLASET('Full',n2,n2,CZERO,ONE,r,n2)
      CALL ZHERK('U','N',n2,n2,DREAL(-ONE),q,n2,DREAL(ONE),r,n2)
      resid = ZLANSY('1','Upper',n2,r,n2,dbwork)
      Result(2) = resid/(eps*MAX(1,n2))
!
!     Generate random m-by-n matrix C and a copy CF
!
      CALL ZLASET('Full',n2,M,CZERO,ONE,c,n2)
      DO j = 1 , M
         CALL ZLARNV(2,iseed,n2,c(1,j))
      ENDDO
      cnorm = ZLANGE('1',n2,M,c,n2,dbwork)
      CALL ZLACPY('Full',n2,M,c,n2,cf,n2)
!
!     Apply Q to C as Q*C
!
      CALL ZTPMLQT('L','N',N,M,k,L,Nb,af(1,np1),M,t,ldt,cf,n2,cf(np1,1),&
     &             n2,work,info)
!
!     Compute |Q*C - Q*C| / |C|
!
      CALL ZGEMM('N','N',n2,M,n2,-ONE,q,n2,c,n2,ONE,cf,n2)
      resid = ZLANGE('1',n2,M,cf,n2,dbwork)
      IF ( cnorm>ZERO ) THEN
         Result(3) = resid/(eps*MAX(1,n2)*cnorm)
      ELSE
         Result(3) = ZERO
      ENDIF
 
!
!     Copy C into CF again
!
      CALL ZLACPY('Full',n2,M,c,n2,cf,n2)
!
!     Apply Q to C as QT*C
!
      CALL ZTPMLQT('L','C',N,M,k,L,Nb,af(1,np1),M,t,ldt,cf,n2,cf(np1,1),&
     &             n2,work,info)
!
!     Compute |QT*C - QT*C| / |C|
!
      CALL ZGEMM('C','N',n2,M,n2,-ONE,q,n2,c,n2,ONE,cf,n2)
      resid = ZLANGE('1',n2,M,cf,n2,dbwork)
 
      IF ( cnorm>ZERO ) THEN
         Result(4) = resid/(eps*MAX(1,n2)*cnorm)
      ELSE
         Result(4) = ZERO
      ENDIF
!
!     Generate random m-by-n matrix D and a copy DF
!
      DO j = 1 , n2
         CALL ZLARNV(2,iseed,M,d(1,j))
      ENDDO
      dnorm = ZLANGE('1',M,n2,d,M,dbwork)
      CALL ZLACPY('Full',M,n2,d,M,df,M)
!
!     Apply Q to D as D*Q
!
      CALL ZTPMLQT('R','N',M,N,k,L,Nb,af(1,np1),M,t,ldt,df,M,df(1,np1), &
     &             M,work,info)
!
!     Compute |D*Q - D*Q| / |D|
!
      CALL ZGEMM('N','N',M,n2,n2,-ONE,d,M,q,n2,ONE,df,M)
      resid = ZLANGE('1',M,n2,df,M,dbwork)
      IF ( cnorm>ZERO ) THEN
         Result(5) = resid/(eps*MAX(1,n2)*dnorm)
      ELSE
         Result(5) = ZERO
      ENDIF
!
!     Copy D into DF again
!
      CALL ZLACPY('Full',M,n2,d,M,df,M)
!
!     Apply Q to D as D*QT
!
      CALL ZTPMLQT('R','C',M,N,k,L,Nb,af(1,np1),M,t,ldt,df,M,df(1,np1), &
     &             M,work,info)
 
!
!     Compute |D*QT - D*QT| / |D|
!
      CALL ZGEMM('N','C',M,n2,n2,-ONE,d,M,q,n2,ONE,df,M)
      resid = ZLANGE('1',M,n2,df,M,dbwork)
      IF ( cnorm>ZERO ) THEN
         Result(6) = resid/(eps*MAX(1,n2)*dnorm)
      ELSE
         Result(6) = ZERO
      ENDIF
!
!     Deallocate all arrays
!
      DEALLOCATE (a,af,q,r,dbwork,rwork,work,t,c,d,cf,df)
      END SUBROUTINE zlqt05
