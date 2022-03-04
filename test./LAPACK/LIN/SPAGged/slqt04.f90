!*==slqt04.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b SLQT04
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE SLQT04(M,N,NB,RESULT)
!
!       .. Scalar Arguments ..
!       INTEGER M, N, NB, LDT
!       .. Return values ..
!       REAL RESULT(6)
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SLQT04 tests SGELQT and SGEMLQT.
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
!>          RESULT is REAL array, dimension (6)
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
      SUBROUTINE SLQT04(M,N,Nb,Result)
      IMPLICIT NONE
!*--SLQT0477
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     April 2012
!
!     .. Scalar Arguments ..
      INTEGER M , N , Nb , ldt
!     .. Return values ..
      REAL Result(6)
!
!  =====================================================================
!
!     ..
!     .. Local allocatable arrays
      REAL , ALLOCATABLE  ::  af(:,:) , q(:,:) , l(:,:) , rwork(:) ,    &
     &                        work(:) , t(:,:) , cf(:,:) , df(:,:) ,    &
     &                        a(:,:) , c(:,:) , d(:,:)
!
!     .. Parameters ..
      REAL ONE , ZERO
      PARAMETER (ZERO=0.0,ONE=1.0)
!     ..
!     .. Local Scalars ..
      INTEGER info , j , k , ll , lwork
      REAL anorm , eps , resid , cnorm , dnorm
!     ..
!     .. Local Arrays ..
      INTEGER iseed(4)
!     ..
!     .. External Functions ..
      REAL SLAMCH , SLANGE , SLANSY
      LOGICAL LSAME
      EXTERNAL SLAMCH , SLANGE , SLANSY , LSAME
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX , MIN
!     ..
!     .. Data statements ..
      DATA iseed/1988 , 1989 , 1990 , 1991/
!
      eps = SLAMCH('Epsilon')
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
         CALL SLARNV(2,iseed,M,a(1,j))
      ENDDO
      CALL SLACPY('Full',M,N,a,M,af,M)
!
!     Factor the matrix A in the array AF.
!
      CALL SGELQT(M,N,Nb,af,M,t,ldt,work,info)
!
!     Generate the n-by-n matrix Q
!
      CALL SLASET('Full',N,N,ZERO,ONE,q,N)
      CALL SGEMLQT('R','N',N,N,k,Nb,af,M,t,ldt,q,N,work,info)
!
!     Copy R
!
      CALL SLASET('Full',M,N,ZERO,ZERO,l,ll)
      CALL SLACPY('Lower',M,N,af,M,l,ll)
!
!     Compute |L - A*Q'| / |A| and store in RESULT(1)
!
      CALL SGEMM('N','T',M,N,N,-ONE,a,M,q,N,ONE,l,ll)
      anorm = SLANGE('1',M,N,a,M,rwork)
      resid = SLANGE('1',M,N,l,ll,rwork)
      IF ( anorm>ZERO ) THEN
         Result(1) = resid/(eps*MAX(1,M)*anorm)
      ELSE
         Result(1) = ZERO
      ENDIF
!
!     Compute |I - Q'*Q| and store in RESULT(2)
!
      CALL SLASET('Full',N,N,ZERO,ONE,l,ll)
      CALL SSYRK('U','C',N,N,-ONE,q,N,ONE,l,ll)
      resid = SLANSY('1','Upper',N,l,ll,rwork)
      Result(2) = resid/(eps*MAX(1,N))
!
!     Generate random m-by-n matrix C and a copy CF
!
      DO j = 1 , M
         CALL SLARNV(2,iseed,N,d(1,j))
      ENDDO
      dnorm = SLANGE('1',N,M,d,N,rwork)
      CALL SLACPY('Full',N,M,d,N,df,N)
!
!     Apply Q to C as Q*C
!
      CALL SGEMLQT('L','N',N,M,k,Nb,af,M,t,Nb,df,N,work,info)
!
!     Compute |Q*D - Q*D| / |D|
!
      CALL SGEMM('N','N',N,M,N,-ONE,q,N,d,N,ONE,df,N)
      resid = SLANGE('1',N,M,df,N,rwork)
      IF ( dnorm>ZERO ) THEN
         Result(3) = resid/(eps*MAX(1,M)*dnorm)
      ELSE
         Result(3) = ZERO
      ENDIF
!
!     Copy D into DF again
!
      CALL SLACPY('Full',N,M,d,N,df,N)
!
!     Apply Q to D as QT*D
!
      CALL SGEMLQT('L','T',N,M,k,Nb,af,M,t,Nb,df,N,work,info)
!
!     Compute |QT*D - QT*D| / |D|
!
      CALL SGEMM('T','N',N,M,N,-ONE,q,N,d,N,ONE,df,N)
      resid = SLANGE('1',N,M,df,N,rwork)
      IF ( dnorm>ZERO ) THEN
         Result(4) = resid/(eps*MAX(1,M)*dnorm)
      ELSE
         Result(4) = ZERO
      ENDIF
!
!     Generate random n-by-m matrix D and a copy DF
!
      DO j = 1 , N
         CALL SLARNV(2,iseed,M,c(1,j))
      ENDDO
      cnorm = SLANGE('1',M,N,c,M,rwork)
      CALL SLACPY('Full',M,N,c,M,cf,M)
!
!     Apply Q to C as C*Q
!
      CALL SGEMLQT('R','N',M,N,k,Nb,af,M,t,Nb,cf,M,work,info)
!
!     Compute |C*Q - C*Q| / |C|
!
      CALL SGEMM('N','N',M,N,N,-ONE,c,M,q,N,ONE,cf,M)
      resid = SLANGE('1',N,M,df,N,rwork)
      IF ( cnorm>ZERO ) THEN
         Result(5) = resid/(eps*MAX(1,M)*dnorm)
      ELSE
         Result(5) = ZERO
      ENDIF
!
!     Copy C into CF again
!
      CALL SLACPY('Full',M,N,c,M,cf,M)
!
!     Apply Q to D as D*QT
!
      CALL SGEMLQT('R','T',M,N,k,Nb,af,M,t,Nb,cf,M,work,info)
!
!     Compute |C*QT - C*QT| / |C|
!
      CALL SGEMM('N','T',M,N,N,-ONE,c,M,q,N,ONE,cf,M)
      resid = SLANGE('1',M,N,cf,M,rwork)
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
      END SUBROUTINE SLQT04