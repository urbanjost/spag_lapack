!*==sqrt05.f90  processed by SPAG 7.51RB at 17:34 on  4 Mar 2022
 
!> \brief \b sqrt05
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE sqrt05(M,N,L,NB,RESULT)
!
!       .. Scalar Arguments ..
!       INTEGER LWORK, M, N, L, NB, LDT
!       .. Return values ..
!       REAL RESULT(6)
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> sqrt05 tests STPQRT and STPMQRT.
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
!>          RESULT is REAL array, dimension (6)
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
!> \ingroup single_lin
!
!  =====================================================================
      SUBROUTINE sqrt05(M,N,L,Nb,Result)
      IMPLICIT NONE
!*--SQRT0585
!
!  -- LAPACK test routine (version 3.8.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     April 2012
!
!     .. Scalar Arguments ..
      INTEGER lwork , M , N , L , Nb , ldt
!     .. Return values ..
      REAL Result(6)
!
!  =====================================================================
!
!     ..
!     .. Local allocatable arrays
      REAL , ALLOCATABLE  ::  af(:,:) , q(:,:) , r(:,:) , rwork(:) ,    &
     &                        work(:) , t(:,:) , cf(:,:) , df(:,:) ,    &
     &                        a(:,:) , c(:,:) , d(:,:)
!
!     .. Parameters ..
      REAL ZERO , ONE
      PARAMETER (ZERO=0.0,ONE=1.0)
!     ..
!     .. Local Scalars ..
      INTEGER info , j , k , m2 , np1
      REAL anorm , eps , resid , cnorm , dnorm
!     ..
!     .. Local Arrays ..
      INTEGER iseed(4)
!     ..
!     .. External Subroutine ..
EXTERNAL    SGEMM,   SLARNV,   STPMQRT,   STPQRT,   SGEMQRT,   SSYRK
EXTERNAL    SLACPY,   SLASET
!     ..
!     .. External Functions ..
      REAL SLAMCH
      REAL SLANGE , SLANSY
      LOGICAL LSAME
EXTERNAL    SLAMCH,   SLANGE,   SLANSY,   LSAME
!     ..
!     .. Data statements ..
      DATA iseed/1988 , 1989 , 1990 , 1991/
!
      eps = SLAMCH('Epsilon')
      k = N
      m2 = M + N
      IF ( M>0 ) THEN
         np1 = N + 1
      ELSE
         np1 = 1
      ENDIF
      lwork = m2*m2*Nb
!
!     Dynamically allocate all arrays
!
      ALLOCATE (a(m2,N),af(m2,N),q(m2,m2),r(m2,m2),rwork(m2),work(lwork)&
     &          ,t(Nb,N),c(m2,N),cf(m2,N),d(N,m2),df(N,m2))
!
!     Put random stuff into A
!
      ldt = Nb
      CALL SLASET('Full',m2,N,ZERO,ZERO,a,m2)
      CALL SLASET('Full',Nb,N,ZERO,ZERO,t,Nb)
      DO j = 1 , N
         CALL SLARNV(2,iseed,j,a(1,j))
      ENDDO
      IF ( M>0 ) THEN
         DO j = 1 , N
            CALL SLARNV(2,iseed,M-L,a(N+1,j))
         ENDDO
      ENDIF
      IF ( L>0 ) THEN
         DO j = 1 , N
            CALL SLARNV(2,iseed,MIN(j,L),a(N+M-L+1,j))
         ENDDO
      ENDIF
!
!     Copy the matrix A to the array AF.
!
      CALL SLACPY('Full',m2,N,a,m2,af,m2)
!
!     Factor the matrix A in the array AF.
!
      CALL STPQRT(M,N,L,Nb,af,m2,af(np1,1),m2,t,ldt,work,info)
!
!     Generate the (M+N)-by-(M+N) matrix Q by applying H to I
!
      CALL SLASET('Full',m2,m2,ZERO,ONE,q,m2)
      CALL SGEMQRT('R','N',m2,m2,k,Nb,af,m2,t,ldt,q,m2,work,info)
!
!     Copy R
!
      CALL SLASET('Full',m2,N,ZERO,ZERO,r,m2)
      CALL SLACPY('Upper',m2,N,af,m2,r,m2)
!
!     Compute |R - Q'*A| / |A| and store in RESULT(1)
!
      CALL SGEMM('T','N',m2,N,m2,-ONE,q,m2,a,m2,ONE,r,m2)
      anorm = SLANGE('1',m2,N,a,m2,rwork)
      resid = SLANGE('1',m2,N,r,m2,rwork)
      IF ( anorm>ZERO ) THEN
         Result(1) = resid/(eps*anorm*MAX(1,m2))
      ELSE
         Result(1) = ZERO
      ENDIF
!
!     Compute |I - Q'*Q| and store in RESULT(2)
!
      CALL SLASET('Full',m2,m2,ZERO,ONE,r,m2)
      CALL SSYRK('U','C',m2,m2,-ONE,q,m2,ONE,r,m2)
      resid = SLANSY('1','Upper',m2,r,m2,rwork)
      Result(2) = resid/(eps*MAX(1,m2))
!
!     Generate random m-by-n matrix C and a copy CF
!
      DO j = 1 , N
         CALL SLARNV(2,iseed,m2,c(1,j))
      ENDDO
      cnorm = SLANGE('1',m2,N,c,m2,rwork)
      CALL SLACPY('Full',m2,N,c,m2,cf,m2)
!
!     Apply Q to C as Q*C
!
      CALL STPMQRT('L','N',M,N,k,L,Nb,af(np1,1),m2,t,ldt,cf,m2,cf(np1,1)&
     &             ,m2,work,info)
!
!     Compute |Q*C - Q*C| / |C|
!
      CALL SGEMM('N','N',m2,N,m2,-ONE,q,m2,c,m2,ONE,cf,m2)
      resid = SLANGE('1',m2,N,cf,m2,rwork)
      IF ( cnorm>ZERO ) THEN
         Result(3) = resid/(eps*MAX(1,m2)*cnorm)
      ELSE
         Result(3) = ZERO
      ENDIF
!
!     Copy C into CF again
!
      CALL SLACPY('Full',m2,N,c,m2,cf,m2)
!
!     Apply Q to C as QT*C
!
      CALL STPMQRT('L','T',M,N,k,L,Nb,af(np1,1),m2,t,ldt,cf,m2,cf(np1,1)&
     &             ,m2,work,info)
!
!     Compute |QT*C - QT*C| / |C|
!
      CALL SGEMM('T','N',m2,N,m2,-ONE,q,m2,c,m2,ONE,cf,m2)
      resid = SLANGE('1',m2,N,cf,m2,rwork)
      IF ( cnorm>ZERO ) THEN
         Result(4) = resid/(eps*MAX(1,m2)*cnorm)
      ELSE
         Result(4) = ZERO
      ENDIF
!
!     Generate random n-by-m matrix D and a copy DF
!
      DO j = 1 , m2
         CALL SLARNV(2,iseed,N,d(1,j))
      ENDDO
      dnorm = SLANGE('1',N,m2,d,N,rwork)
      CALL SLACPY('Full',N,m2,d,N,df,N)
!
!     Apply Q to D as D*Q
!
      CALL STPMQRT('R','N',N,M,N,L,Nb,af(np1,1),m2,t,ldt,df,N,df(1,np1),&
     &             N,work,info)
!
!     Compute |D*Q - D*Q| / |D|
!
      CALL SGEMM('N','N',N,m2,m2,-ONE,d,N,q,m2,ONE,df,N)
      resid = SLANGE('1',N,m2,df,N,rwork)
      IF ( cnorm>ZERO ) THEN
         Result(5) = resid/(eps*MAX(1,m2)*dnorm)
      ELSE
         Result(5) = ZERO
      ENDIF
!
!     Copy D into DF again
!
      CALL SLACPY('Full',N,m2,d,N,df,N)
!
!     Apply Q to D as D*QT
!
      CALL STPMQRT('R','T',N,M,N,L,Nb,af(np1,1),m2,t,ldt,df,N,df(1,np1),&
     &             N,work,info)
 
!
!     Compute |D*QT - D*QT| / |D|
!
      CALL SGEMM('N','T',N,m2,m2,-ONE,d,N,q,m2,ONE,df,N)
      resid = SLANGE('1',N,m2,df,N,rwork)
      IF ( cnorm>ZERO ) THEN
         Result(6) = resid/(eps*MAX(1,m2)*dnorm)
      ELSE
         Result(6) = ZERO
      ENDIF
!
!     Deallocate all arrays
!
      DEALLOCATE (a,af,q,r,rwork,work,t,c,d,cf,df)
      END SUBROUTINE sqrt05
