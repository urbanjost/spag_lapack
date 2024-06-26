!*==dlqt04.f90  processed by SPAG 7.51RB at 17:34 on  4 Mar 2022
!> \brief \b dlqt04
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE dlqt04(M,N,NB,RESULT)
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
!> dlqt04 tests DGELQT and DGEMLQT.
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
!>          RESULT(1) = | A - L Q |
!>          RESULT(2) = | I - Q Q^H |
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
      SUBROUTINE dlqt04(M,N,Nb,Result)
      IMPLICIT NONE
!*--DLQT0477
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
      DOUBLE PRECISION , ALLOCATABLE  ::  af(:,:) , q(:,:) , l(:,:) ,   &
     &   rwork(:) , work(:) , t(:,:) , cf(:,:) , df(:,:) , a(:,:) ,     &
     &   c(:,:) , d(:,:)
!
!     .. Parameters ..
      DOUBLE PRECISION ONE , ZERO
      PARAMETER (ZERO=0.0,ONE=1.0)
!     ..
!     .. Local Scalars ..
      INTEGER info , j , k , ll , lwork
      DOUBLE PRECISION anorm , eps , resid , cnorm , dnorm
!     ..
!     .. Local Arrays ..
      INTEGER iseed(4)
!     ..
!     .. External Functions ..
      DOUBLE PRECISION DLAMCH , DLANGE , DLANSY
      LOGICAL LSAME
EXTERNAL    DLAMCH,   DLANGE,   DLANSY,   LSAME
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX , MIN
!     ..
!     .. Data statements ..
      DATA iseed/1988 , 1989 , 1990 , 1991/
!
      eps = DLAMCH('Epsilon')
      k = MIN(M,N)
      ll = MAX(M,N)
      lwork = MAX(2,ll)*MAX(2,ll)*Nb
!
!     Dynamically allocate local arrays
!
      ALLOCATE (a(M,N),af(M,N),q(N,N),l(ll,N),rwork(ll),work(lwork),    &
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
      CALL DGELQT(M,N,Nb,af,M,t,ldt,work,info)
!
!     Generate the n-by-n matrix Q
!
      CALL DLASET('Full',N,N,ZERO,ONE,q,N)
      CALL DGEMLQT('R','N',N,N,k,Nb,af,M,t,ldt,q,N,work,info)
!
!     Copy R
!
      CALL DLASET('Full',M,N,ZERO,ZERO,l,ll)
      CALL DLACPY('Lower',M,N,af,M,l,ll)
!
!     Compute |L - A*Q'| / |A| and store in RESULT(1)
!
      CALL DGEMM('N','T',M,N,N,-ONE,a,M,q,N,ONE,l,ll)
      anorm = DLANGE('1',M,N,a,M,rwork)
      resid = DLANGE('1',M,N,l,ll,rwork)
      IF ( anorm>ZERO ) THEN
         Result(1) = resid/(eps*MAX(1,M)*anorm)
      ELSE
         Result(1) = ZERO
      ENDIF
!
!     Compute |I - Q'*Q| and store in RESULT(2)
!
      CALL DLASET('Full',N,N,ZERO,ONE,l,ll)
      CALL DSYRK('U','C',N,N,-ONE,q,N,ONE,l,ll)
      resid = DLANSY('1','Upper',N,l,ll,rwork)
      Result(2) = resid/(eps*MAX(1,N))
!
!     Generate random m-by-n matrix C and a copy CF
!
      DO j = 1 , M
         CALL DLARNV(2,iseed,N,d(1,j))
      ENDDO
      dnorm = DLANGE('1',N,M,d,N,rwork)
      CALL DLACPY('Full',N,M,d,N,df,N)
!
!     Apply Q to C as Q*C
!
      CALL DGEMLQT('L','N',N,M,k,Nb,af,M,t,Nb,df,N,work,info)
!
!     Compute |Q*D - Q*D| / |D|
!
      CALL DGEMM('N','N',N,M,N,-ONE,q,N,d,N,ONE,df,N)
      resid = DLANGE('1',N,M,df,N,rwork)
      IF ( dnorm>ZERO ) THEN
         Result(3) = resid/(eps*MAX(1,M)*dnorm)
      ELSE
         Result(3) = ZERO
      ENDIF
!
!     Copy D into DF again
!
      CALL DLACPY('Full',N,M,d,N,df,N)
!
!     Apply Q to D as QT*D
!
      CALL DGEMLQT('L','T',N,M,k,Nb,af,M,t,Nb,df,N,work,info)
!
!     Compute |QT*D - QT*D| / |D|
!
      CALL DGEMM('T','N',N,M,N,-ONE,q,N,d,N,ONE,df,N)
      resid = DLANGE('1',N,M,df,N,rwork)
      IF ( dnorm>ZERO ) THEN
         Result(4) = resid/(eps*MAX(1,M)*dnorm)
      ELSE
         Result(4) = ZERO
      ENDIF
!
!     Generate random n-by-m matrix D and a copy DF
!
      DO j = 1 , N
         CALL DLARNV(2,iseed,M,c(1,j))
      ENDDO
      cnorm = DLANGE('1',M,N,c,M,rwork)
      CALL DLACPY('Full',M,N,c,M,cf,M)
!
!     Apply Q to C as C*Q
!
      CALL DGEMLQT('R','N',M,N,k,Nb,af,M,t,Nb,cf,M,work,info)
!
!     Compute |C*Q - C*Q| / |C|
!
      CALL DGEMM('N','N',M,N,N,-ONE,c,M,q,N,ONE,cf,M)
      resid = DLANGE('1',N,M,df,N,rwork)
      IF ( cnorm>ZERO ) THEN
         Result(5) = resid/(eps*MAX(1,M)*dnorm)
      ELSE
         Result(5) = ZERO
      ENDIF
!
!     Copy C into CF again
!
      CALL DLACPY('Full',M,N,c,M,cf,M)
!
!     Apply Q to D as D*QT
!
      CALL DGEMLQT('R','T',M,N,k,Nb,af,M,t,Nb,cf,M,work,info)
!
!     Compute |C*QT - C*QT| / |C|
!
      CALL DGEMM('N','T',M,N,N,-ONE,c,M,q,N,ONE,cf,M)
      resid = DLANGE('1',M,N,cf,M,rwork)
      IF ( cnorm>ZERO ) THEN
         Result(6) = resid/(eps*MAX(1,M)*dnorm)
      ELSE
         Result(6) = ZERO
      ENDIF
!
!     Deallocate all arrays
!
      DEALLOCATE (a,af,q,l,rwork,work,t,c,d,cf,df)
!
      END SUBROUTINE dlqt04
