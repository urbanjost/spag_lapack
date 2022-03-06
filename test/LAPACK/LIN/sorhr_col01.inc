!*==sorhr_col01.f90  processed by SPAG 7.51RB at 17:34 on  4 Mar 2022
!> \brief \b sorhr_col01
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE sorhr_col01( M, N, MB1, NB1, NB2, RESULT )
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
!> sorhr_col01 tests SORGTSQR and SORHR_COL using SLATSQR, SGEMQRT.
!> Therefore, SLATSQR (part of SGEQR), SGEMQRT (part of SGEMQR)
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
!>            Q_qr is an implicit m-by-m orthogonal Q matrix, the result
!>            of factorization in blocked WY-representation,
!>            stored in SGEQRT output format.
!>
!>            R is a n-by-n upper-triangular matrix,
!>
!>            0 is a (m-n)-by-n zero matrix,
!>
!>            Q is an explicit m-by-m orthogonal matrix Q = Q_gr * I
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
!>            computed using SGEMQRT,
!>
!>            Q * C, (Q**H) * C, D * Q, D * (Q**H)  are
!>            computed using SGEMM.
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
!> \ingroup single_lin
!
!  =====================================================================
      SUBROUTINE sorhr_col01(M,N,Mb1,Nb1,Nb2,Result)
      IMPLICIT NONE
!*--SORHR_COL01123
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
      REAL , ALLOCATABLE  ::  a(:,:) , af(:,:) , q(:,:) , r(:,:) ,      &
     &                        rwork(:) , work(:) , t1(:,:) , t2(:,:) ,  &
     &                        diag(:) , c(:,:) , cf(:,:) , d(:,:) ,     &
     &                        df(:,:)
!
!     .. Parameters ..
      REAL ONE , ZERO
      PARAMETER (ZERO=0.0E+0,ONE=1.0E+0)
!     ..
!     .. Local Scalars ..
      LOGICAL testzeros
      INTEGER info , i , j , k , l , lwork , nb1_ub , nb2_ub , nrb
      REAL anorm , eps , resid , cnorm , dnorm
!     ..
!     .. Local Arrays ..
      INTEGER iseed(4)
      REAL workquery(1)
!     ..
!     .. External Functions ..
      REAL SLAMCH , SLANGE , SLANSY
EXTERNAL    SLAMCH,   SLANGE,   SLANSY
!     ..
!     .. External Subroutines ..
EXTERNAL    SLACPY,   SLARNV,   SLASET,   SLATSQR,   SORHR_COL
EXTERNAL    SORGTSQR,   SSCAL,   SGEMM,   SGEMQRT,   SSYRK
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
         CALL SLARNV(2,iseed,M,a(1,j))
      ENDDO
      IF ( testzeros ) THEN
         IF ( M>=4 ) THEN
            DO j = 1 , N
               CALL SLARNV(2,iseed,M/2,a(M/4,j))
            ENDDO
         ENDIF
      ENDIF
      CALL SLACPY('Full',M,N,a,M,af,M)
!
!     Number of row blocks in SLATSQR
!
      nrb = MAX(1,CEILING(REAL(M-N)/REAL(Mb1-N)))
!
      ALLOCATE (t1(Nb1,N*nrb))
      ALLOCATE (t2(Nb2,N))
      ALLOCATE (diag(N))
!
!     Begin determine LWORK for the array WORK and allocate memory.
!
!     SLATSQR requires NB1 to be bounded by N.
!
      nb1_ub = MIN(Nb1,N)
!
!     SGEMQRT requires NB2 to be bounded by N.
!
      nb2_ub = MIN(Nb2,N)
!
      CALL SLATSQR(M,N,Mb1,nb1_ub,af,M,t1,Nb1,workquery,-1,info)
      lwork = INT(workquery(1))
      CALL SORGTSQR(M,N,Mb1,Nb1,af,M,t1,Nb1,workquery,-1,info)
 
      lwork = MAX(lwork,INT(workquery(1)))
!
!     In SGEMQRT, WORK is N*NB2_UB if SIDE = 'L',
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
      SRNamt = 'SLATSQR'
      CALL SLATSQR(M,N,Mb1,nb1_ub,af,M,t1,Nb1,work,lwork,info)
!
!     Copy the factor R into the array R.
!
      SRNamt = 'SLACPY'
      CALL SLACPY('U',N,N,af,M,r,M)
!
!     Reconstruct the orthogonal matrix Q.
!
      SRNamt = 'SORGTSQR'
      CALL SORGTSQR(M,N,Mb1,Nb1,af,M,t1,Nb1,work,lwork,info)
!
!     Perform the Householder reconstruction, the result is stored
!     the arrays AF and T2.
!
      SRNamt = 'SORHR_COL'
      CALL SORHR_COL(M,N,Nb2,af,M,t2,Nb2,diag,info)
!
!     Compute the factor R_hr corresponding to the Householder
!     reconstructed Q_hr and place it in the upper triangle of AF to
!     match the Q storage format in SGEQRT. R_hr = R_tsqr * S,
!     this means changing the sign of I-th row of the matrix R_tsqr
!     according to sign of of I-th diagonal element DIAG(I) of the
!     matrix S.
!
      SRNamt = 'SLACPY'
      CALL SLACPY('U',N,N,r,M,af,M)
!
      DO i = 1 , N
         IF ( diag(i)==-ONE ) CALL SSCAL(N+1-i,-ONE,af(i,i),M)
      ENDDO
!
!     End Householder reconstruction routines.
!
!
!     Generate the m-by-m matrix Q
!
      CALL SLASET('Full',M,M,ZERO,ONE,q,M)
!
      SRNamt = 'SGEMQRT'
      CALL SGEMQRT('L','N',M,M,k,nb2_ub,af,M,t2,Nb2,q,M,work,info)
!
!     Copy R
!
      CALL SLASET('Full',M,N,ZERO,ZERO,r,M)
!
      CALL SLACPY('Upper',M,N,af,M,r,M)
!
!     TEST 1
!     Compute |R - (Q**T)*A| / ( eps * m * |A| ) and store in RESULT(1)
!
      CALL SGEMM('T','N',M,N,M,-ONE,q,M,a,M,ONE,r,M)
!
      anorm = SLANGE('1',M,N,a,M,rwork)
      resid = SLANGE('1',M,N,r,M,rwork)
      IF ( anorm>ZERO ) THEN
         Result(1) = resid/(eps*MAX(1,M)*anorm)
      ELSE
         Result(1) = ZERO
      ENDIF
!
!     TEST 2
!     Compute |I - (Q**T)*Q| / ( eps * m ) and store in RESULT(2)
!
      CALL SLASET('Full',M,M,ZERO,ONE,r,M)
      CALL SSYRK('U','T',M,M,-ONE,q,M,ONE,r,M)
      resid = SLANSY('1','Upper',M,r,M,rwork)
      Result(2) = resid/(eps*MAX(1,M))
!
!     Generate random m-by-n matrix C
!
      DO j = 1 , N
         CALL SLARNV(2,iseed,M,c(1,j))
      ENDDO
      cnorm = SLANGE('1',M,N,c,M,rwork)
      CALL SLACPY('Full',M,N,c,M,cf,M)
!
!     Apply Q to C as Q*C = CF
!
      SRNamt = 'SGEMQRT'
      CALL SGEMQRT('L','N',M,N,k,nb2_ub,af,M,t2,Nb2,cf,M,work,info)
!
!     TEST 3
!     Compute |CF - Q*C| / ( eps *  m * |C| )
!
      CALL SGEMM('N','N',M,N,M,-ONE,q,M,c,M,ONE,cf,M)
      resid = SLANGE('1',M,N,cf,M,rwork)
      IF ( cnorm>ZERO ) THEN
         Result(3) = resid/(eps*MAX(1,M)*cnorm)
      ELSE
         Result(3) = ZERO
      ENDIF
!
!     Copy C into CF again
!
      CALL SLACPY('Full',M,N,c,M,cf,M)
!
!     Apply Q to C as (Q**T)*C = CF
!
      SRNamt = 'SGEMQRT'
      CALL SGEMQRT('L','T',M,N,k,nb2_ub,af,M,t2,Nb2,cf,M,work,info)
!
!     TEST 4
!     Compute |CF - (Q**T)*C| / ( eps * m * |C|)
!
      CALL SGEMM('T','N',M,N,M,-ONE,q,M,c,M,ONE,cf,M)
      resid = SLANGE('1',M,N,cf,M,rwork)
      IF ( cnorm>ZERO ) THEN
         Result(4) = resid/(eps*MAX(1,M)*cnorm)
      ELSE
         Result(4) = ZERO
      ENDIF
!
!     Generate random n-by-m matrix D and a copy DF
!
      DO j = 1 , M
         CALL SLARNV(2,iseed,N,d(1,j))
      ENDDO
      dnorm = SLANGE('1',N,M,d,N,rwork)
      CALL SLACPY('Full',N,M,d,N,df,N)
!
!     Apply Q to D as D*Q = DF
!
      SRNamt = 'SGEMQRT'
      CALL SGEMQRT('R','N',N,M,k,nb2_ub,af,M,t2,Nb2,df,N,work,info)
!
!     TEST 5
!     Compute |DF - D*Q| / ( eps * m * |D| )
!
      CALL SGEMM('N','N',N,M,M,-ONE,d,N,q,M,ONE,df,N)
      resid = SLANGE('1',N,M,df,N,rwork)
      IF ( dnorm>ZERO ) THEN
         Result(5) = resid/(eps*MAX(1,M)*dnorm)
      ELSE
         Result(5) = ZERO
      ENDIF
!
!     Copy D into DF again
!
      CALL SLACPY('Full',N,M,d,N,df,N)
!
!     Apply Q to D as D*QT = DF
!
      SRNamt = 'SGEMQRT'
      CALL SGEMQRT('R','T',N,M,k,nb2_ub,af,M,t2,Nb2,df,N,work,info)
!
!     TEST 6
!     Compute |DF - D*(Q**T)| / ( eps * m * |D| )
!
      CALL SGEMM('N','T',N,M,M,-ONE,d,N,q,M,ONE,df,N)
      resid = SLANGE('1',N,M,df,N,rwork)
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
!     End of sorhr_col01
!
      END SUBROUTINE sorhr_col01
