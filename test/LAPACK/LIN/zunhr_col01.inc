!*==zunhr_col01.f90  processed by SPAG 7.51RB at 17:35 on  4 Mar 2022
!> \brief \b zunhr_col01
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE zunhr_col01( M, N, MB1, NB1, NB2, RESULT )
!
!       .. Scalar Arguments ..
!       INTEGER           M, N, MB1, NB1, NB2
!       .. Return values ..
!       DOUBLE PRECISION  RESULT(6)
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> zunhr_col01 tests ZUNGTSQR and ZUNHR_COL using ZLATSQR, ZGEMQRT.
!> Therefore, ZLATSQR (part of ZGEQR), ZGEMQRT (part of ZGEMQR)
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
!>          RESULT is DOUBLE PRECISION array, dimension (6)
!>          Results of each of the six tests below.
!>
!>            A is a m-by-n test input matrix to be factored.
!>            so that A = Q_gr * ( R )
!>                               ( 0 ),
!>
!>            Q_qr is an implicit m-by-m unitary Q matrix, the result
!>            of factorization in blocked WY-representation,
!>            stored in ZGEQRT output format.
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
!>            computed using ZGEMQRT,
!>
!>            Q * C, (Q**H) * C, D * Q, D * (Q**H)  are
!>            computed using ZGEMM.
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
!> \ingroup complex16_lin
!
!  =====================================================================
      SUBROUTINE zunhr_col01(M,N,Mb1,Nb1,Nb2,Result)
      IMPLICIT NONE
!*--ZUNHR_COL01123
!
!  -- LAPACK test routine (version 3.10.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2020
!
!     .. Scalar Arguments ..
      INTEGER M , N , Mb1 , Nb1 , Nb2
!     .. Return values ..
      DOUBLE PRECISION Result(6)
!
!  =====================================================================
!
!     ..
!     .. Local allocatable arrays
      COMPLEX*16 , ALLOCATABLE  ::  a(:,:) , af(:,:) , q(:,:) , r(:,:) ,&
     &                              work(:) , t1(:,:) , t2(:,:) ,       &
     &                              diag(:) , c(:,:) , cf(:,:) , d(:,:) &
     &                              , df(:,:)
      DOUBLE PRECISION , ALLOCATABLE  ::  rwork(:)
!
!     .. Parameters ..
      DOUBLE PRECISION,parameter :: ZERO=0.0D+0
      COMPLEX*16,parameter :: CONE=(1.0D+0,0.0D+0),CZERO=(0.0D+0,0.0D+0)
      double precision, parameter :: dpone=1.0d0
!     ..
!     .. Local Scalars ..
      LOGICAL testzeros
      INTEGER info , i , j , k , l , lwork , nb1_ub , nb2_ub , nrb
      DOUBLE PRECISION anorm , eps , resid , cnorm , dnorm
!     ..
!     .. Local Arrays ..
      INTEGER iseed(4)
      COMPLEX*16 workquery(1)
!     ..
!     .. External Functions ..
      DOUBLE PRECISION DLAMCH , ZLANGE , ZLANSY
EXTERNAL    DLAMCH,   ZLANGE,   ZLANSY
!     ..
!     .. External Subroutines ..
EXTERNAL    ZLACPY,   ZLARNV,   ZLASET,   ZLATSQR,   ZUNHR_COL
EXTERNAL    ZUNGTSQR,   ZSCAL,   ZGEMM,   ZGEMQRT,   ZHERK
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC CEILING , DBLE , MAX , MIN
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
      eps = DLAMCH('Epsilon')
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
!     Number of row blocks in ZLATSQR
!
      nrb = MAX(1,CEILING(DBLE(M-N)/DBLE(Mb1-N)))
!
      ALLOCATE (t1(Nb1,N*nrb))
      ALLOCATE (t2(Nb2,N))
      ALLOCATE (diag(N))
!
!     Begin determine LWORK for the array WORK and allocate memory.
!
!     ZLATSQR requires NB1 to be bounded by N.
!
      nb1_ub = MIN(Nb1,N)
!
!     ZGEMQRT requires NB2 to be bounded by N.
!
      nb2_ub = MIN(Nb2,N)
!
      CALL ZLATSQR(M,N,Mb1,nb1_ub,af,M,t1,Nb1,workquery,-1,info)
      lwork = INT(workquery(1))
      CALL ZUNGTSQR(M,N,Mb1,Nb1,af,M,t1,Nb1,workquery,-1,info)
 
      lwork = MAX(lwork,INT(workquery(1)))
!
!     In ZGEMQRT, WORK is N*NB2_UB if SIDE = 'L',
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
      SRNamt = 'ZLATSQR'
      CALL ZLATSQR(M,N,Mb1,nb1_ub,af,M,t1,Nb1,work,lwork,info)
!
!     Copy the factor R into the array R.
!
      SRNamt = 'ZLACPY'
      CALL ZLACPY('U',N,N,af,M,r,M)
!
!     Reconstruct the orthogonal matrix Q.
!
      SRNamt = 'ZUNGTSQR'
      CALL ZUNGTSQR(M,N,Mb1,Nb1,af,M,t1,Nb1,work,lwork,info)
!
!     Perform the Householder reconstruction, the result is stored
!     the arrays AF and T2.
!
      SRNamt = 'ZUNHR_COL'
      CALL ZUNHR_COL(M,N,Nb2,af,M,t2,Nb2,diag,info)
!
!     Compute the factor R_hr corresponding to the Householder
!     reconstructed Q_hr and place it in the upper triangle of AF to
!     match the Q storage format in ZGEQRT. R_hr = R_tsqr * S,
!     this means changing the sign of I-th row of the matrix R_tsqr
!     according to sign of of I-th diagonal element DIAG(I) of the
!     matrix S.
!
      SRNamt = 'ZLACPY'
      CALL ZLACPY('U',N,N,r,M,af,M)
!
      DO i = 1 , N
         IF ( diag(i)==-CONE ) CALL ZSCAL(N+1-i,-CONE,af(i,i),M)
      ENDDO
!
!     End Householder reconstruction routines.
!
!
!     Generate the m-by-m matrix Q
!
      CALL ZLASET('Full',M,M,CZERO,CONE,q,M)
!
      SRNamt = 'ZGEMQRT'
      CALL ZGEMQRT('L','N',M,M,k,nb2_ub,af,M,t2,Nb2,q,M,work,info)
!
!     Copy R
!
      CALL ZLASET('Full',M,N,CZERO,CZERO,r,M)
!
      CALL ZLACPY('Upper',M,N,af,M,r,M)
!
!     TEST 1
!     Compute |R - (Q**H)*A| / ( eps * m * |A| ) and store in RESULT(1)
!
      CALL ZGEMM('C','N',M,N,M,-CONE,q,M,a,M,CONE,r,M)
!
      anorm = ZLANGE('1',M,N,a,M,rwork)
      resid = ZLANGE('1',M,N,r,M,rwork)
      IF ( anorm>ZERO ) THEN
         Result(1) = resid/(eps*MAX(1,M)*anorm)
      ELSE
         Result(1) = ZERO
      ENDIF
!
!     TEST 2
!     Compute |I - (Q**H)*Q| / ( eps * m ) and store in RESULT(2)
!
      CALL ZLASET('Full',M,M,CZERO,CONE,r,M)
      CALL ZHERK('U','C',M,M,-dpONE,q,M,dpONE,r,M)

      resid = ZLANSY('1','Upper',M,r,M,rwork)
      Result(2) = resid/(eps*MAX(1,M))
!
!     Generate random m-by-n matrix C
!
      DO j = 1 , N
         CALL ZLARNV(2,iseed,M,c(1,j))
      ENDDO
      cnorm = ZLANGE('1',M,N,c,M,rwork)
      CALL ZLACPY('Full',M,N,c,M,cf,M)
!
!     Apply Q to C as Q*C = CF
!
      SRNamt = 'ZGEMQRT'
      CALL ZGEMQRT('L','N',M,N,k,nb2_ub,af,M,t2,Nb2,cf,M,work,info)
!
!     TEST 3
!     Compute |CF - Q*C| / ( eps *  m * |C| )
!
      CALL ZGEMM('N','N',M,N,M,-CONE,q,M,c,M,CONE,cf,M)
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
!     Apply Q to C as (Q**H)*C = CF
!
      SRNamt = 'ZGEMQRT'
      CALL ZGEMQRT('L','C',M,N,k,nb2_ub,af,M,t2,Nb2,cf,M,work,info)
!
!     TEST 4
!     Compute |CF - (Q**H)*C| / ( eps * m * |C|)
!
      CALL ZGEMM('C','N',M,N,M,-CONE,q,M,c,M,CONE,cf,M)
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
!     Apply Q to D as D*Q = DF
!
      SRNamt = 'ZGEMQRT'
      CALL ZGEMQRT('R','N',N,M,k,nb2_ub,af,M,t2,Nb2,df,N,work,info)
!
!     TEST 5
!     Compute |DF - D*Q| / ( eps * m * |D| )
!
      CALL ZGEMM('N','N',N,M,M,-CONE,d,N,q,M,CONE,df,N)
      resid = ZLANGE('1',N,M,df,N,rwork)
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
!     Apply Q to D as D*QT = DF
!
      SRNamt = 'ZGEMQRT'
      CALL ZGEMQRT('R','C',N,M,k,nb2_ub,af,M,t2,Nb2,df,N,work,info)
!
!     TEST 6
!     Compute |DF - D*(Q**H)| / ( eps * m * |D| )
!
      CALL ZGEMM('N','C',N,M,M,-CONE,d,N,q,M,CONE,df,N)
      resid = ZLANGE('1',N,M,df,N,rwork)
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
!     End of zunhr_col01
!
      END SUBROUTINE zunhr_col01
