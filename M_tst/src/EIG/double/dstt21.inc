!*==dstt21.f90  processed by SPAG 7.51RB at 17:37 on  4 Mar 2022
!> \brief \b dstt21
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE dstt21( N, KBAND, AD, AE, SD, SE, U, LDU, WORK,
!                          RESULT )
!
!       .. Scalar Arguments ..
!       INTEGER            KBAND, LDU, N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   AD( * ), AE( * ), RESULT( 2 ), SD( * ),
!      $                   SE( * ), U( LDU, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> dstt21 checks a decomposition of the form
!>
!>    A = U S U'
!>
!> where ' means transpose, A is symmetric tridiagonal, U is orthogonal,
!> and S is diagonal (if KBAND=0) or symmetric tridiagonal (if KBAND=1).
!> Two tests are performed:
!>
!>    RESULT(1) = | A - U S U' | / ( |A| n ulp )
!>
!>    RESULT(2) = | I - UU' | / ( n ulp )
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The size of the matrix.  If it is zero, dstt21 does nothing.
!>          It must be at least zero.
!> \endverbatim
!>
!> \param[in] KBAND
!> \verbatim
!>          KBAND is INTEGER
!>          The bandwidth of the matrix S.  It may only be zero or one.
!>          If zero, then S is diagonal, and SE is not referenced.  If
!>          one, then S is symmetric tri-diagonal.
!> \endverbatim
!>
!> \param[in] AD
!> \verbatim
!>          AD is DOUBLE PRECISION array, dimension (N)
!>          The diagonal of the original (unfactored) matrix A.  A is
!>          assumed to be symmetric tridiagonal.
!> \endverbatim
!>
!> \param[in] AE
!> \verbatim
!>          AE is DOUBLE PRECISION array, dimension (N-1)
!>          The off-diagonal of the original (unfactored) matrix A.  A
!>          is assumed to be symmetric tridiagonal.  AE(1) is the (1,2)
!>          and (2,1) element, AE(2) is the (2,3) and (3,2) element, etc.
!> \endverbatim
!>
!> \param[in] SD
!> \verbatim
!>          SD is DOUBLE PRECISION array, dimension (N)
!>          The diagonal of the (symmetric tri-) diagonal matrix S.
!> \endverbatim
!>
!> \param[in] SE
!> \verbatim
!>          SE is DOUBLE PRECISION array, dimension (N-1)
!>          The off-diagonal of the (symmetric tri-) diagonal matrix S.
!>          Not referenced if KBSND=0.  If KBAND=1, then AE(1) is the
!>          (1,2) and (2,1) element, SE(2) is the (2,3) and (3,2)
!>          element, etc.
!> \endverbatim
!>
!> \param[in] U
!> \verbatim
!>          U is DOUBLE PRECISION array, dimension (LDU, N)
!>          The orthogonal matrix in the decomposition.
!> \endverbatim
!>
!> \param[in] LDU
!> \verbatim
!>          LDU is INTEGER
!>          The leading dimension of U.  LDU must be at least N.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (N*(N+1))
!> \endverbatim
!>
!> \param[out] RESULT
!> \verbatim
!>          RESULT is DOUBLE PRECISION array, dimension (2)
!>          The values computed by the two tests described above.  The
!>          values are currently limited to 1/ulp, to avoid overflow.
!>          RESULT(1) is always modified.
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
      SUBROUTINE dstt21(N,Kband,Ad,Ae,Sd,Se,U,Ldu,Work,Result)
      IMPLICIT NONE
!*--DSTT21130
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Kband , Ldu , N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION Ad(*) , Ae(*) , Result(2) , Sd(*) , Se(*) ,      &
     &                 U(Ldu,*) , Work(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO , ONE
      PARAMETER (ZERO=0.0D0,ONE=1.0D0)
!     ..
!     .. Local Scalars ..
      INTEGER j
      DOUBLE PRECISION anorm , temp1 , temp2 , ulp , unfl , wnorm
!     ..
!     .. External Functions ..
      DOUBLE PRECISION DLAMCH , DLANGE , DLANSY
      EXTERNAL DLAMCH , DLANGE , DLANSY
!     ..
!     .. External Subroutines ..
      EXTERNAL DGEMM , DLASET , DSYR , DSYR2
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS , DBLE , MAX , MIN
!     ..
!     .. Executable Statements ..
!
!     1)      Constants
!
      Result(1) = ZERO
      Result(2) = ZERO
      IF ( N<=0 ) RETURN
!
      unfl = DLAMCH('Safe minimum')
      ulp = DLAMCH('Precision')
!
!     Do Test 1
!
!     Copy A & Compute its 1-Norm:
!
      CALL DLASET('Full',N,N,ZERO,ZERO,Work,N)
!
      anorm = ZERO
      temp1 = ZERO
!
      DO j = 1 , N - 1
         Work((N+1)*(j-1)+1) = Ad(j)
         Work((N+1)*(j-1)+2) = Ae(j)
         temp2 = ABS(Ae(j))
         anorm = MAX(anorm,ABS(Ad(j))+temp1+temp2)
         temp1 = temp2
      ENDDO
!
      Work(N**2) = Ad(N)
      anorm = MAX(anorm,ABS(Ad(N))+temp1,unfl)
!
!     Norm of A - USU'
!
      DO j = 1 , N
         CALL DSYR('L',N,-Sd(j),U(1,j),1,Work,N)
      ENDDO
!
      IF ( N>1 .AND. Kband==1 ) THEN
         DO j = 1 , N - 1
            CALL DSYR2('L',N,-Se(j),U(1,j),1,U(1,j+1),1,Work,N)
         ENDDO
      ENDIF
!
      wnorm = DLANSY('1','L',N,Work,N,Work(N**2+1))
!
      IF ( anorm>wnorm ) THEN
         Result(1) = (wnorm/anorm)/(N*ulp)
      ELSEIF ( anorm<ONE ) THEN
         Result(1) = (MIN(wnorm,N*anorm)/anorm)/(N*ulp)
      ELSE
         Result(1) = MIN(wnorm/anorm,DBLE(N))/(N*ulp)
      ENDIF
!
!     Do Test 2
!
!     Compute  UU' - I
!
      CALL DGEMM('N','C',N,N,N,ONE,U,Ldu,U,Ldu,ZERO,Work,N)
!
      DO j = 1 , N
         Work((N+1)*(j-1)+1) = Work((N+1)*(j-1)+1) - ONE
      ENDDO
!
      Result(2) = MIN(DBLE(N),DLANGE('1',N,N,Work,N,Work(N**2+1)))      &
     &            /(N*ulp)
!
!
!     End of dstt21
!
      END SUBROUTINE dstt21
