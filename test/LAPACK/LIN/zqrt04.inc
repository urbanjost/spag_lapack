!*==zqrt04.f90  processed by SPAG 7.51RB at 17:35 on  4 Mar 2022
!> \brief \b zqrt04
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE zqrt04(M,N,NB,RESULT)
!
!       .. Scalar Arguments ..
!       INTEGER M, N, NB, LDT
!       .. Return values ..
!       DOUBLE PRECISION RESULT(6)
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> zqrt04 tests ZGEQRT and ZGEMQRT.
!> \endverbatim
!
!  Arguments:
!  ==========
!
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
!>
!> \param[in] NB
!> \verbatim
!>          NB is INTEGER
!>          Block size of test matrix.  NB <= Min(M,N).
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
!> \ingroup complex16_lin
!
!  =====================================================================
      SUBROUTINE zqrt04(M,N,Nb,Result)
      IMPLICIT NONE
!*--ZQRT0477
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     April 2012
!
!     .. Scalar Arguments ..
      INTEGER M , N , Nb , ldt
!     .. Return values ..
      DOUBLE PRECISION Result(6)
!
!  =====================================================================
!
!     ..
!     .. Local allocatable arrays
      COMPLEX*16 , ALLOCATABLE  ::  af(:,:) , q(:,:) , r(:,:) , work(:) &
     &                              , t(:,:) , cf(:,:) , df(:,:) ,      &
     &                              a(:,:) , c(:,:) , d(:,:)
      DOUBLE PRECISION , ALLOCATABLE  ::  rwork(:)
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO
      COMPLEX*16 ONE , CZERO
      PARAMETER (ZERO=0.0,ONE=(1.0,0.0),CZERO=(0.0,0.0))
!     ..
!     .. Local Scalars ..
      INTEGER info , j , k , l , lwork
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
!     .. Intrinsic Functions ..
      INTRINSIC MAX , MIN
!     ..
!     .. Data statements ..
      DATA iseed/1988 , 1989 , 1990 , 1991/
!
      eps = DLAMCH('Epsilon')
      k = MIN(M,N)
      l = MAX(M,N)
      lwork = MAX(2,l)*MAX(2,l)*Nb
!
!     Dynamically allocate local arrays
!
      ALLOCATE (a(M,N),af(M,N),q(M,M),r(M,l),rwork(l),work(lwork),      &
     &          t(Nb,N),c(M,N),cf(M,N),d(N,M),df(N,M))
!
!     Put random numbers into A and copy to AF
!
      ldt = Nb
      DO j = 1 , N
         CALL ZLARNV(2,iseed,M,a(1,j))
      ENDDO
      CALL ZLACPY('Full',M,N,a,M,af,M)
!
!     Factor the matrix A in the array AF.
!
      CALL ZGEQRT(M,N,Nb,af,M,t,ldt,work,info)
!
!     Generate the m-by-m matrix Q
!
      CALL ZLASET('Full',M,M,CZERO,ONE,q,M)
      CALL ZGEMQRT('R','N',M,M,k,Nb,af,M,t,ldt,q,M,work,info)
!
!     Copy R
!
      CALL ZLASET('Full',M,N,CZERO,CZERO,r,M)
      CALL ZLACPY('Upper',M,N,af,M,r,M)
!
!     Compute |R - Q'*A| / |A| and store in RESULT(1)
!
      CALL ZGEMM('C','N',M,N,M,-ONE,q,M,a,M,ONE,r,M)
      anorm = ZLANGE('1',M,N,a,M,rwork)
      resid = ZLANGE('1',M,N,r,M,rwork)
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
      resid = ZLANSY('1','Upper',M,r,M,rwork)
      Result(2) = resid/(eps*MAX(1,M))
!
!     Generate random m-by-n matrix C and a copy CF
!
      DO j = 1 , N
         CALL ZLARNV(2,iseed,M,c(1,j))
      ENDDO
      cnorm = ZLANGE('1',M,N,c,M,rwork)
      CALL ZLACPY('Full',M,N,c,M,cf,M)
!
!     Apply Q to C as Q*C
!
      CALL ZGEMQRT('L','N',M,N,k,Nb,af,M,t,Nb,cf,M,work,info)
!
!     Compute |Q*C - Q*C| / |C|
!
      CALL ZGEMM('N','N',M,N,M,-ONE,q,M,c,M,ONE,cf,M)
      resid = ZLANGE('1',M,N,cf,M,rwork)
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
      CALL ZGEMQRT('L','C',M,N,k,Nb,af,M,t,Nb,cf,M,work,info)
!
!     Compute |QT*C - QT*C| / |C|
!
      CALL ZGEMM('C','N',M,N,M,-ONE,q,M,c,M,ONE,cf,M)
      resid = ZLANGE('1',M,N,cf,M,rwork)
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
      dnorm = ZLANGE('1',N,M,d,N,rwork)
      CALL ZLACPY('Full',N,M,d,N,df,N)
!
!     Apply Q to D as D*Q
!
      CALL ZGEMQRT('R','N',N,M,k,Nb,af,M,t,Nb,df,N,work,info)
!
!     Compute |D*Q - D*Q| / |D|
!
      CALL ZGEMM('N','N',N,M,M,-ONE,d,N,q,M,ONE,df,N)
      resid = ZLANGE('1',N,M,df,N,rwork)
      IF ( cnorm>ZERO ) THEN
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
      CALL ZGEMQRT('R','C',N,M,k,Nb,af,M,t,Nb,df,N,work,info)
!
!     Compute |D*QT - D*QT| / |D|
!
      CALL ZGEMM('N','C',N,M,M,-ONE,d,N,q,M,ONE,df,N)
      resid = ZLANGE('1',N,M,df,N,rwork)
      IF ( cnorm>ZERO ) THEN
         Result(6) = resid/(eps*MAX(1,M)*dnorm)
      ELSE
         Result(6) = ZERO
      ENDIF
!
!     Deallocate all arrays
!
      DEALLOCATE (a,af,q,r,rwork,work,t,c,d,cf,df)
!
      END SUBROUTINE zqrt04
