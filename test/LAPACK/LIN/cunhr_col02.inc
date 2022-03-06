!*==cunhr_col02.f90  processed by SPAG 7.51RB at 17:34 on  4 Mar 2022
!> \brief \b cunhr_col02
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE cunhr_col02( M, N, MB1, NB1, NB2, RESULT )
!
!       .. Scalar Arguments ..
!       INTEGER           M, N, MB1, NB1, NB2
!       .. Return values ..
!       REAL              RESULT(6)
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> cunhr_col02 tests CUNGTSQR_ROW and CUNHR_COL inside CGETSQRHRT
!> (which calls CLATSQR, CUNGTSQR_ROW and CUNHR_COL) using CGEMQRT.
!> Therefore, CLATSQR (part of CGEQR), CGEMQRT (part of CGEMQR)
!> have to be tested before this test.
!>
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
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          Number of columns in test matrix.
!> \endverbatim
!> \param[in] MB1
!> \verbatim
!>          MB1 is INTEGER
!>          Number of row in row block in an input test matrix.
!> \endverbatim
!>
!> \param[in] NB1
!> \verbatim
!>          NB1 is INTEGER
!>          Number of columns in column block an input test matrix.
!> \endverbatim
!>
!> \param[in] NB2
!> \verbatim
!>          NB2 is INTEGER
!>          Number of columns in column block in an output test matrix.
!> \endverbatim
!>
!> \param[out] RESULT
!> \verbatim
!>          RESULT is REAL array, dimension (6)
!>          Results of each of the six tests below.
!>
!>            A is a m-by-n test input matrix to be factored.
!>            so that A = Q_gr * ( R )
!>                               ( 0 ),
!>
!>            Q_qr is an implicit m-by-m unitary Q matrix, the result
!>            of factorization in blocked WY-representation,
!>            stored in CGEQRT output format.
!>
!>            R is a n-by-n upper-triangular matrix,
!>
!>            0 is a (m-n)-by-n zero matrix,
!>
!>            Q is an explicit m-by-m unitary matrix Q = Q_gr * I
!>
!>            C is an m-by-n random matrix,
!>
!>            D is an n-by-m random matrix.
!>
!>          The six tests are:
!>
!>          RESULT(1) = |R - (Q**H) * A| / ( eps * m * |A| )
!>            is equivalent to test for | A - Q * R | / (eps * m * |A|),
!>
!>          RESULT(2) = |I - (Q**H) * Q| / ( eps * m ),
!>
!>          RESULT(3) = | Q_qr * C - Q * C | / (eps * m * |C|),
!>
!>          RESULT(4) = | (Q_gr**H) * C - (Q**H) * C | / (eps * m * |C|)
!>
!>          RESULT(5) = | D * Q_qr - D * Q | / (eps * m * |D|)
!>
!>          RESULT(6) = | D * (Q_qr**H) - D * (Q**H) | / (eps * m * |D|),
!>
!>          where:
!>            Q_qr * C, (Q_gr**H) * C, D * Q_qr, D * (Q_qr**H) are
!>            computed using CGEMQRT,
!>
!>            Q * C, (Q**H) * C, D * Q, D * (Q**H)  are
!>            computed using CGEMM.
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
!> \date November 2020
!
!> \ingroup complex_lin
!
!  =====================================================================
      SUBROUTINE cunhr_col02(M,N,Mb1,Nb1,Nb2,Result)
      IMPLICIT NONE
!*--CUNHR_COL02124
!
!  -- LAPACK test routine (version 3.10.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2020
!
!     .. Scalar Arguments ..
      INTEGER M , N , Mb1 , Nb1 , Nb2
!     .. Return values ..
      REAL Result(6)
!
!  =====================================================================
!
!     ..
!     .. Local allocatable arrays
      COMPLEX , ALLOCATABLE  ::  a(:,:) , af(:,:) , q(:,:) , r(:,:) ,   &
     &                           work(:) , t1(:,:) , t2(:,:) , diag(:) ,&
     &                           c(:,:) , cf(:,:) , d(:,:) , df(:,:)
      REAL , ALLOCATABLE  ::  rwork(:)
!
!     .. Parameters ..
      REAL ZERO
      PARAMETER (ZERO=0.0E+0)
      COMPLEX CONE , CZERO
      real,parameter :: one=1.0
      PARAMETER (CONE=(1.0E+0,0.0E+0),CZERO=(0.0E+0,0.0E+0))
!     ..
!     .. Local Scalars ..
      LOGICAL testzeros
      INTEGER info , j , k , l , lwork , nb2_ub , nrb
      REAL anorm , eps , resid , cnorm , dnorm
!     ..
!     .. Local Arrays ..
      INTEGER iseed(4)
      COMPLEX workquery(1)
!     ..
!     .. External Functions ..
      REAL SLAMCH , CLANGE , CLANSY
EXTERNAL    SLAMCH,   CLANGE,   CLANSY
!     ..
!     .. External Subroutines ..
EXTERNAL    CLACPY,   CLARNV,   CLASET,   CGETSQRHRT,   CSCAL,   CGEMM,   CGEMQRT,   CHERK
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC CEILING , REAL , MAX , MIN
!     ..
!     .. Scalars in Common ..
      CHARACTER(LEN=32) SRNamt
!     ..
!     .. Common blocks ..
      COMMON /SRMNAMC/ SRNamt
!     ..
!     .. Data statements ..
      DATA iseed/1988 , 1989 , 1990 , 1991/
!
!     TEST MATRICES WITH HALF OF MATRIX BEING ZEROS
!
      testzeros = .FALSE.
!
      eps = SLAMCH('Epsilon')
      k = MIN(M,N)
      l = MAX(M,N,1)
!
!     Dynamically allocate local arrays
!
      ALLOCATE (a(M,N),af(M,N),q(l,l),r(M,l),rwork(l),c(M,N),cf(M,N),   &
     &          d(N,M),df(N,M))
!
!     Put random numbers into A and copy to AF
!
      DO j = 1 , N
         CALL CLARNV(2,iseed,M,a(1,j))
      ENDDO
      IF ( testzeros ) THEN
         IF ( M>=4 ) THEN
            DO j = 1 , N
               CALL CLARNV(2,iseed,M/2,a(M/4,j))
            ENDDO
         ENDIF
      ENDIF
      CALL CLACPY('Full',M,N,a,M,af,M)
!
!     Number of row blocks in CLATSQR
!
      nrb = MAX(1,CEILING(REAL(M-N)/REAL(Mb1-N)))
!
      ALLOCATE (t1(Nb1,N*nrb))
      ALLOCATE (t2(Nb2,N))
      ALLOCATE (diag(N))
!
!     Begin determine LWORK for the array WORK and allocate memory.
!
!     CGEMQRT requires NB2 to be bounded by N.
!
      nb2_ub = MIN(Nb2,N)
!
!
      CALL CGETSQRHRT(M,N,Mb1,Nb1,Nb2,af,M,t2,Nb2,workquery,-1,info)
!
      lwork = INT(workquery(1))
!
!     In CGEMQRT, WORK is N*NB2_UB if SIDE = 'L',
!                or  M*NB2_UB if SIDE = 'R'.
!
      lwork = MAX(lwork,nb2_ub*N,nb2_ub*M)
!
      ALLOCATE (work(lwork))
!
!     End allocate memory for WORK.
!
!
!     Begin Householder reconstruction routines
!
!     Factor the matrix A in the array AF.
!
      SRNamt = 'CGETSQRHRT'
      CALL CGETSQRHRT(M,N,Mb1,Nb1,Nb2,af,M,t2,Nb2,work,lwork,info)
!
!     End Householder reconstruction routines.
!
!
!     Generate the m-by-m matrix Q
!
      CALL CLASET('Full',M,M,CZERO,CONE,q,M)
!
      SRNamt = 'CGEMQRT'
      CALL CGEMQRT('L','N',M,M,k,nb2_ub,af,M,t2,Nb2,q,M,work,info)
!
!     Copy R
!
      CALL CLASET('Full',M,N,CZERO,CZERO,r,M)
!
      CALL CLACPY('Upper',M,N,af,M,r,M)
!
!     TEST 1
!     Compute |R - (Q**T)*A| / ( eps * m * |A| ) and store in RESULT(1)
!
      CALL CGEMM('C','N',M,N,M,-CONE,q,M,a,M,CONE,r,M)
!
      anorm = CLANGE('1',M,N,a,M,rwork)
      resid = CLANGE('1',M,N,r,M,rwork)
      IF ( anorm>ZERO ) THEN
         Result(1) = resid/(eps*MAX(1,M)*anorm)
      ELSE
         Result(1) = ZERO
      ENDIF
!
!     TEST 2
!     Compute |I - (Q**T)*Q| / ( eps * m ) and store in RESULT(2)
!
      CALL CLASET('Full',M,M,CZERO,CONE,r,M)
      CALL CHERK('U','C',M,M,-ONE,q,M,ONE,r,M)
      resid = CLANSY('1','Upper',M,r,M,rwork)
      Result(2) = resid/(eps*MAX(1,M))
!
!     Generate random m-by-n matrix C
!
      DO j = 1 , N
         CALL CLARNV(2,iseed,M,c(1,j))
      ENDDO
      cnorm = CLANGE('1',M,N,c,M,rwork)
      CALL CLACPY('Full',M,N,c,M,cf,M)
!
!     Apply Q to C as Q*C = CF
!
      SRNamt = 'CGEMQRT'
      CALL CGEMQRT('L','N',M,N,k,nb2_ub,af,M,t2,Nb2,cf,M,work,info)
!
!     TEST 3
!     Compute |CF - Q*C| / ( eps *  m * |C| )
!
      CALL CGEMM('N','N',M,N,M,-CONE,q,M,c,M,CONE,cf,M)
      resid = CLANGE('1',M,N,cf,M,rwork)
      IF ( cnorm>ZERO ) THEN
         Result(3) = resid/(eps*MAX(1,M)*cnorm)
      ELSE
         Result(3) = ZERO
      ENDIF
!
!     Copy C into CF again
!
      CALL CLACPY('Full',M,N,c,M,cf,M)
!
!     Apply Q to C as (Q**T)*C = CF
!
      SRNamt = 'CGEMQRT'
      CALL CGEMQRT('L','C',M,N,k,nb2_ub,af,M,t2,Nb2,cf,M,work,info)
!
!     TEST 4
!     Compute |CF - (Q**T)*C| / ( eps * m * |C|)
!
      CALL CGEMM('C','N',M,N,M,-CONE,q,M,c,M,CONE,cf,M)
      resid = CLANGE('1',M,N,cf,M,rwork)
      IF ( cnorm>ZERO ) THEN
         Result(4) = resid/(eps*MAX(1,M)*cnorm)
      ELSE
         Result(4) = ZERO
      ENDIF
!
!     Generate random n-by-m matrix D and a copy DF
!
      DO j = 1 , M
         CALL CLARNV(2,iseed,N,d(1,j))
      ENDDO
      dnorm = CLANGE('1',N,M,d,N,rwork)
      CALL CLACPY('Full',N,M,d,N,df,N)
!
!     Apply Q to D as D*Q = DF
!
      SRNamt = 'CGEMQRT'
      CALL CGEMQRT('R','N',N,M,k,nb2_ub,af,M,t2,Nb2,df,N,work,info)
!
!     TEST 5
!     Compute |DF - D*Q| / ( eps * m * |D| )
!
      CALL CGEMM('N','N',N,M,M,-CONE,d,N,q,M,CONE,df,N)
      resid = CLANGE('1',N,M,df,N,rwork)
      IF ( dnorm>ZERO ) THEN
         Result(5) = resid/(eps*MAX(1,M)*dnorm)
      ELSE
         Result(5) = ZERO
      ENDIF
!
!     Copy D into DF again
!
      CALL CLACPY('Full',N,M,d,N,df,N)
!
!     Apply Q to D as D*QT = DF
!
      SRNamt = 'CGEMQRT'
      CALL CGEMQRT('R','C',N,M,k,nb2_ub,af,M,t2,Nb2,df,N,work,info)
!
!     TEST 6
!     Compute |DF - D*(Q**T)| / ( eps * m * |D| )
!
      CALL CGEMM('N','C',N,M,M,-CONE,d,N,q,M,CONE,df,N)
      resid = CLANGE('1',N,M,df,N,rwork)
      IF ( dnorm>ZERO ) THEN
         Result(6) = resid/(eps*MAX(1,M)*dnorm)
      ELSE
         Result(6) = ZERO
      ENDIF
!
!     Deallocate all arrays
!
      DEALLOCATE (a,af,q,r,rwork,work,t1,t2,diag,c,d,cf,df)
!
!
!     End of cunhr_col02
!
      END SUBROUTINE cunhr_col02
