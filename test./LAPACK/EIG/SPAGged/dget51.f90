!*==dget51.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b DGET51
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DGET51( ITYPE, N, A, LDA, B, LDB, U, LDU, V, LDV, WORK,
!                          RESULT )
!
!       .. Scalar Arguments ..
!       INTEGER            ITYPE, LDA, LDB, LDU, LDV, N
!       DOUBLE PRECISION   RESULT
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), U( LDU, * ),
!      $                   V( LDV, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!>      DGET51  generally checks a decomposition of the form
!>
!>              A = U B V'
!>
!>      where ' means transpose and U and V are orthogonal.
!>
!>      Specifically, if ITYPE=1
!>
!>              RESULT = | A - U B V' | / ( |A| n ulp )
!>
!>      If ITYPE=2, then:
!>
!>              RESULT = | A - B | / ( |A| n ulp )
!>
!>      If ITYPE=3, then:
!>
!>              RESULT = | I - UU' | / ( n ulp )
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] ITYPE
!> \verbatim
!>          ITYPE is INTEGER
!>          Specifies the type of tests to be performed.
!>          =1: RESULT = | A - U B V' | / ( |A| n ulp )
!>          =2: RESULT = | A - B | / ( |A| n ulp )
!>          =3: RESULT = | I - UU' | / ( n ulp )
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The size of the matrix.  If it is zero, DGET51 does nothing.
!>          It must be at least zero.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (LDA, N)
!>          The original (unfactored) matrix.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of A.  It must be at least 1
!>          and at least N.
!> \endverbatim
!>
!> \param[in] B
!> \verbatim
!>          B is DOUBLE PRECISION array, dimension (LDB, N)
!>          The factored matrix.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of B.  It must be at least 1
!>          and at least N.
!> \endverbatim
!>
!> \param[in] U
!> \verbatim
!>          U is DOUBLE PRECISION array, dimension (LDU, N)
!>          The orthogonal matrix on the left-hand side in the
!>          decomposition.
!>          Not referenced if ITYPE=2
!> \endverbatim
!>
!> \param[in] LDU
!> \verbatim
!>          LDU is INTEGER
!>          The leading dimension of U.  LDU must be at least N and
!>          at least 1.
!> \endverbatim
!>
!> \param[in] V
!> \verbatim
!>          V is DOUBLE PRECISION array, dimension (LDV, N)
!>          The orthogonal matrix on the left-hand side in the
!>          decomposition.
!>          Not referenced if ITYPE=2
!> \endverbatim
!>
!> \param[in] LDV
!> \verbatim
!>          LDV is INTEGER
!>          The leading dimension of V.  LDV must be at least N and
!>          at least 1.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (2*N**2)
!> \endverbatim
!>
!> \param[out] RESULT
!> \verbatim
!>          RESULT is DOUBLE PRECISION
!>          The values computed by the test specified by ITYPE.  The
!>          value is currently limited to 1/ulp, to avoid overflow.
!>          Errors are flagged by RESULT=10/ulp.
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
!> \date December 2016
!
!> \ingroup double_eig
!
!  =====================================================================
      SUBROUTINE DGET51(Itype,N,A,Lda,B,Ldb,U,Ldu,V,Ldv,Work,Result)
      IMPLICIT NONE
!*--DGET51152
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Itype , Lda , Ldb , Ldu , Ldv , N
      DOUBLE PRECISION Result
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION A(Lda,*) , B(Ldb,*) , U(Ldu,*) , V(Ldv,*) ,      &
     &                 Work(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO , ONE , TEN
      PARAMETER (ZERO=0.0D0,ONE=1.0D0,TEN=10.0D0)
!     ..
!     .. Local Scalars ..
      INTEGER jcol , jdiag , jrow
      DOUBLE PRECISION anorm , ulp , unfl , wnorm
!     ..
!     .. External Functions ..
      DOUBLE PRECISION DLAMCH , DLANGE
      EXTERNAL DLAMCH , DLANGE
!     ..
!     .. External Subroutines ..
      EXTERNAL DGEMM , DLACPY
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC DBLE , MAX , MIN
!     ..
!     .. Executable Statements ..
!
      Result = ZERO
      IF ( N<=0 ) RETURN
!
!     Constants
!
      unfl = DLAMCH('Safe minimum')
      ulp = DLAMCH('Epsilon')*DLAMCH('Base')
!
!     Some Error Checks
!
      IF ( Itype<1 .OR. Itype>3 ) THEN
         Result = TEN/ulp
         RETURN
      ENDIF
!
      IF ( Itype<=2 ) THEN
!
!        Tests scaled by the norm(A)
!
         anorm = MAX(DLANGE('1',N,N,A,Lda,Work),unfl)
!
         IF ( Itype==1 ) THEN
!
!           ITYPE=1: Compute W = A - UBV'
!
            CALL DLACPY(' ',N,N,A,Lda,Work,N)
            CALL DGEMM('N','N',N,N,N,ONE,U,Ldu,B,Ldb,ZERO,Work(N**2+1), &
     &                 N)
!
            CALL DGEMM('N','C',N,N,N,-ONE,Work(N**2+1),N,V,Ldv,ONE,Work,&
     &                 N)
!
         ELSE
!
!           ITYPE=2: Compute W = A - B
!
            CALL DLACPY(' ',N,N,B,Ldb,Work,N)
!
            DO jcol = 1 , N
               DO jrow = 1 , N
                  Work(jrow+N*(jcol-1)) = Work(jrow+N*(jcol-1))         &
     &               - A(jrow,jcol)
               ENDDO
            ENDDO
         ENDIF
!
!        Compute norm(W)/ ( ulp*norm(A) )
!
         wnorm = DLANGE('1',N,N,Work,N,Work(N**2+1))
!
         IF ( anorm>wnorm ) THEN
            Result = (wnorm/anorm)/(N*ulp)
         ELSEIF ( anorm<ONE ) THEN
            Result = (MIN(wnorm,N*anorm)/anorm)/(N*ulp)
         ELSE
            Result = MIN(wnorm/anorm,DBLE(N))/(N*ulp)
         ENDIF
!
      ELSE
!
!        Tests not scaled by norm(A)
!
!        ITYPE=3: Compute  UU' - I
!
         CALL DGEMM('N','C',N,N,N,ONE,U,Ldu,U,Ldu,ZERO,Work,N)
!
         DO jdiag = 1 , N
            Work((N+1)*(jdiag-1)+1) = Work((N+1)*(jdiag-1)+1) - ONE
         ENDDO
!
         Result = MIN(DLANGE('1',N,N,Work,N,Work(N**2+1)),DBLE(N))      &
     &            /(N*ulp)
      ENDIF
!
!
!     End of DGET51
!
      END SUBROUTINE DGET51
