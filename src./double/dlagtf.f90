!*==dlagtf.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b DLAGTF computes an LU factorization of a matrix T-λI, where T is a general tridiagonal matrix, and λ a scalar, using partial pivoting with row interchanges.
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download DLAGTF + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlagtf.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlagtf.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlagtf.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLAGTF( N, A, LAMBDA, B, C, TOL, D, IN, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, N
!       DOUBLE PRECISION   LAMBDA, TOL
!       ..
!       .. Array Arguments ..
!       INTEGER            IN( * )
!       DOUBLE PRECISION   A( * ), B( * ), C( * ), D( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DLAGTF factorizes the matrix (T - lambda*I), where T is an n by n
!> tridiagonal matrix and lambda is a scalar, as
!>
!>    T - lambda*I = PLU,
!>
!> where P is a permutation matrix, L is a unit lower tridiagonal matrix
!> with at most one non-zero sub-diagonal elements per column and U is
!> an upper triangular matrix with at most two non-zero super-diagonal
!> elements per column.
!>
!> The factorization is obtained by Gaussian elimination with partial
!> pivoting and implicit row scaling.
!>
!> The parameter LAMBDA is included in the routine so that DLAGTF may
!> be used, in conjunction with DLAGTS, to obtain eigenvectors of T by
!> inverse iteration.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix T.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is DOUBLE PRECISION array, dimension (N)
!>          On entry, A must contain the diagonal elements of T.
!>
!>          On exit, A is overwritten by the n diagonal elements of the
!>          upper triangular matrix U of the factorization of T.
!> \endverbatim
!>
!> \param[in] LAMBDA
!> \verbatim
!>          LAMBDA is DOUBLE PRECISION
!>          On entry, the scalar lambda.
!> \endverbatim
!>
!> \param[in,out] B
!> \verbatim
!>          B is DOUBLE PRECISION array, dimension (N-1)
!>          On entry, B must contain the (n-1) super-diagonal elements of
!>          T.
!>
!>          On exit, B is overwritten by the (n-1) super-diagonal
!>          elements of the matrix U of the factorization of T.
!> \endverbatim
!>
!> \param[in,out] C
!> \verbatim
!>          C is DOUBLE PRECISION array, dimension (N-1)
!>          On entry, C must contain the (n-1) sub-diagonal elements of
!>          T.
!>
!>          On exit, C is overwritten by the (n-1) sub-diagonal elements
!>          of the matrix L of the factorization of T.
!> \endverbatim
!>
!> \param[in] TOL
!> \verbatim
!>          TOL is DOUBLE PRECISION
!>          On entry, a relative tolerance used to indicate whether or
!>          not the matrix (T - lambda*I) is nearly singular. TOL should
!>          normally be chose as approximately the largest relative error
!>          in the elements of T. For example, if the elements of T are
!>          correct to about 4 significant figures, then TOL should be
!>          set to about 5*10**(-4). If TOL is supplied as less than eps,
!>          where eps is the relative machine precision, then the value
!>          eps is used in place of TOL.
!> \endverbatim
!>
!> \param[out] D
!> \verbatim
!>          D is DOUBLE PRECISION array, dimension (N-2)
!>          On exit, D is overwritten by the (n-2) second super-diagonal
!>          elements of the matrix U of the factorization of T.
!> \endverbatim
!>
!> \param[out] IN
!> \verbatim
!>          IN is INTEGER array, dimension (N)
!>          On exit, IN contains details of the permutation matrix P. If
!>          an interchange occurred at the kth step of the elimination,
!>          then IN(k) = 1, otherwise IN(k) = 0. The element IN(n)
!>          returns the smallest positive integer j such that
!>
!>             abs( u(j,j) ) <= norm( (T - lambda*I)(j) )*TOL,
!>
!>          where norm( A(j) ) denotes the sum of the absolute values of
!>          the jth row of the matrix A. If no such j exists then IN(n)
!>          is returned as zero. If IN(n) is returned as positive, then a
!>          diagonal element of U is small, indicating that
!>          (T - lambda*I) is singular or nearly singular,
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -k, the kth argument had an illegal value
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
!> \ingroup auxOTHERcomputational
!
!  =====================================================================
      SUBROUTINE DLAGTF(N,A,Lambda,B,C,Tol,D,In,Info)
      USE F77KINDS                        
      USE S_DLAMCH
      USE S_XERBLA
      IMPLICIT NONE
!*--DLAGTF163
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ZERO = 0.0D+0
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER , INTENT(IN) :: N
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: A
      REAL(R8KIND) , INTENT(IN) :: Lambda
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: B
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: C
      REAL(R8KIND) , INTENT(IN) :: Tol
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: D
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: In
      INTEGER :: Info
!
! Local variable declarations rewritten by SPAG
!
      REAL(R8KIND) :: eps , mult , piv1 , piv2 , scale1 , scale2 ,      &
     &                temp , tl
      INTEGER :: k
!
! End of declarations rewritten by SPAG
!
!     ..
!     .. Array Arguments ..
!     ..
!
! =====================================================================
!
!     .. Parameters ..
!     ..
!     .. Local Scalars ..
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. External Functions ..
!     ..
!     .. External Subroutines ..
!     ..
!     .. Executable Statements ..
!
      Info = 0
      IF ( N<0 ) THEN
         Info = -1
         CALL XERBLA('DLAGTF',-Info)
         RETURN
      ENDIF
!
      IF ( N==0 ) RETURN
!
      A(1) = A(1) - Lambda
      In(N) = 0
      IF ( N==1 ) THEN
         IF ( A(1)==ZERO ) In(1) = 1
         RETURN
      ENDIF
!
      eps = DLAMCH('Epsilon')
!
      tl = MAX(Tol,eps)
      scale1 = ABS(A(1)) + ABS(B(1))
      DO k = 1 , N - 1
         A(k+1) = A(k+1) - Lambda
         scale2 = ABS(C(k)) + ABS(A(k+1))
         IF ( k<(N-1) ) scale2 = scale2 + ABS(B(k+1))
         IF ( A(k)==ZERO ) THEN
            piv1 = ZERO
         ELSE
            piv1 = ABS(A(k))/scale1
         ENDIF
         IF ( C(k)==ZERO ) THEN
            In(k) = 0
            piv2 = ZERO
            scale1 = scale2
            IF ( k<(N-1) ) D(k) = ZERO
         ELSE
            piv2 = ABS(C(k))/scale2
            IF ( piv2<=piv1 ) THEN
               In(k) = 0
               scale1 = scale2
               C(k) = C(k)/A(k)
               A(k+1) = A(k+1) - C(k)*B(k)
               IF ( k<(N-1) ) D(k) = ZERO
            ELSE
               In(k) = 1
               mult = A(k)/C(k)
               A(k) = C(k)
               temp = A(k+1)
               A(k+1) = B(k) - mult*temp
               IF ( k<(N-1) ) THEN
                  D(k) = B(k+1)
                  B(k+1) = -mult*D(k)
               ENDIF
               B(k) = temp
               C(k) = mult
            ENDIF
         ENDIF
         IF ( (MAX(piv1,piv2)<=tl) .AND. (In(N)==0) ) In(N) = k
      ENDDO
      IF ( (ABS(A(N))<=scale1*tl) .AND. (In(N)==0) ) In(N) = N
!
!
!     End of DLAGTF
!
      END SUBROUTINE DLAGTF
