!*==dqrt04.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b DQRT04
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DQRT04(M,N,NB,RESULT)
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
!> DQRT04 tests DGEQRT and DGEMQRT.
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
!> \ingroup double_lin
!
!  =====================================================================
      SUBROUTINE DQRT04(M,N,Nb,Result)
      IMPLICIT NONE
!*--DQRT0477
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
      DOUBLE PRECISION , ALLOCATABLE  ::  af(:,:) , q(:,:) , r(:,:) ,   &
     &   rwork(:) , work(:) , t(:,:) , cf(:,:) , df(:,:) , a(:,:) ,     &
     &   c(:,:) , d(:,:)
!
!     .. Parameters ..
      DOUBLE PRECISION ONE , ZERO
      PARAMETER (ZERO=0.0,ONE=1.0)
!     ..
!     .. Local Scalars ..
      INTEGER info , j , k , l , lwork
      DOUBLE PRECISION anorm , eps , resid , cnorm , dnorm
!     ..
!     .. Local Arrays ..
      INTEGER iseed(4)
!     ..
!     .. External Functions ..
      DOUBLE PRECISION DLAMCH , DLANGE , DLANSY
      LOGICAL LSAME
      EXTERNAL DLAMCH , DLANGE , DLANSY , LSAME
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
         CALL DLARNV(2,iseed,M,a(1,j))
      ENDDO
      CALL DLACPY('Full',M,N,a,M,af,M)
!
!     Factor the matrix A in the array AF.
!
      CALL DGEQRT(M,N,Nb,af,M,t,ldt,work,info)
!
!     Generate the m-by-m matrix Q
!
      CALL DLASET('Full',M,M,ZERO,ONE,q,M)
      CALL DGEMQRT('R','N',M,M,k,Nb,af,M,t,ldt,q,M,work,info)
!
!     Copy R
!
      CALL DLASET('Full',M,N,ZERO,ZERO,r,M)
      CALL DLACPY('Upper',M,N,af,M,r,M)
!
!     Compute |R - Q'*A| / |A| and store in RESULT(1)
!
      CALL DGEMM('T','N',M,N,M,-ONE,q,M,a,M,ONE,r,M)
      anorm = DLANGE('1',M,N,a,M,rwork)
      resid = DLANGE('1',M,N,r,M,rwork)
      IF ( anorm>ZERO ) THEN
         Result(1) = resid/(eps*MAX(1,M)*anorm)
      ELSE
         Result(1) = ZERO
      ENDIF
!
!     Compute |I - Q'*Q| and store in RESULT(2)
!
      CALL DLASET('Full',M,M,ZERO,ONE,r,M)
      CALL DSYRK('U','C',M,M,-ONE,q,M,ONE,r,M)
      resid = DLANSY('1','Upper',M,r,M,rwork)
      Result(2) = resid/(eps*MAX(1,M))
!
!     Generate random m-by-n matrix C and a copy CF
!
      DO j = 1 , N
         CALL DLARNV(2,iseed,M,c(1,j))
      ENDDO
      cnorm = DLANGE('1',M,N,c,M,rwork)
      CALL DLACPY('Full',M,N,c,M,cf,M)
!
!     Apply Q to C as Q*C
!
      CALL DGEMQRT('L','N',M,N,k,Nb,af,M,t,Nb,cf,M,work,info)
!
!     Compute |Q*C - Q*C| / |C|
!
      CALL DGEMM('N','N',M,N,M,-ONE,q,M,c,M,ONE,cf,M)
      resid = DLANGE('1',M,N,cf,M,rwork)
      IF ( cnorm>ZERO ) THEN
         Result(3) = resid/(eps*MAX(1,M)*cnorm)
      ELSE
         Result(3) = ZERO
      ENDIF
!
!     Copy C into CF again
!
      CALL DLACPY('Full',M,N,c,M,cf,M)
!
!     Apply Q to C as QT*C
!
      CALL DGEMQRT('L','T',M,N,k,Nb,af,M,t,Nb,cf,M,work,info)
!
!     Compute |QT*C - QT*C| / |C|
!
      CALL DGEMM('T','N',M,N,M,-ONE,q,M,c,M,ONE,cf,M)
      resid = DLANGE('1',M,N,cf,M,rwork)
      IF ( cnorm>ZERO ) THEN
         Result(4) = resid/(eps*MAX(1,M)*cnorm)
      ELSE
         Result(4) = ZERO
      ENDIF
!
!     Generate random n-by-m matrix D and a copy DF
!
      DO j = 1 , M
         CALL DLARNV(2,iseed,N,d(1,j))
      ENDDO
      dnorm = DLANGE('1',N,M,d,N,rwork)
      CALL DLACPY('Full',N,M,d,N,df,N)
!
!     Apply Q to D as D*Q
!
      CALL DGEMQRT('R','N',N,M,k,Nb,af,M,t,Nb,df,N,work,info)
!
!     Compute |D*Q - D*Q| / |D|
!
      CALL DGEMM('N','N',N,M,M,-ONE,d,N,q,M,ONE,df,N)
      resid = DLANGE('1',N,M,df,N,rwork)
      IF ( cnorm>ZERO ) THEN
         Result(5) = resid/(eps*MAX(1,M)*dnorm)
      ELSE
         Result(5) = ZERO
      ENDIF
!
!     Copy D into DF again
!
      CALL DLACPY('Full',N,M,d,N,df,N)
!
!     Apply Q to D as D*QT
!
      CALL DGEMQRT('R','T',N,M,k,Nb,af,M,t,Nb,df,N,work,info)
!
!     Compute |D*QT - D*QT| / |D|
!
      CALL DGEMM('N','T',N,M,M,-ONE,d,N,q,M,ONE,df,N)
      resid = DLANGE('1',N,M,df,N,rwork)
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
      END SUBROUTINE DQRT04
