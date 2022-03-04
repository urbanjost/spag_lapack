!*==zptcon.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b ZPTCON
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZPTCON + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zptcon.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zptcon.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zptcon.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZPTCON( N, D, E, ANORM, RCOND, RWORK, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, N
!       DOUBLE PRECISION   ANORM, RCOND
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   D( * ), RWORK( * )
!       COMPLEX*16         E( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZPTCON computes the reciprocal of the condition number (in the
!> 1-norm) of a complex Hermitian positive definite tridiagonal matrix
!> using the factorization A = L*D*L**H or A = U**H*D*U computed by
!> ZPTTRF.
!>
!> Norm(inv(A)) is computed by a direct method, and the reciprocal of
!> the condition number is computed as
!>                  RCOND = 1 / (ANORM * norm(inv(A))).
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in] D
!> \verbatim
!>          D is DOUBLE PRECISION array, dimension (N)
!>          The n diagonal elements of the diagonal matrix D from the
!>          factorization of A, as computed by ZPTTRF.
!> \endverbatim
!>
!> \param[in] E
!> \verbatim
!>          E is COMPLEX*16 array, dimension (N-1)
!>          The (n-1) off-diagonal elements of the unit bidiagonal factor
!>          U or L from the factorization of A, as computed by ZPTTRF.
!> \endverbatim
!>
!> \param[in] ANORM
!> \verbatim
!>          ANORM is DOUBLE PRECISION
!>          The 1-norm of the original matrix A.
!> \endverbatim
!>
!> \param[out] RCOND
!> \verbatim
!>          RCOND is DOUBLE PRECISION
!>          The reciprocal of the condition number of the matrix A,
!>          computed as RCOND = 1/(ANORM * AINVNM), where AINVNM is the
!>          1-norm of inv(A) computed in this routine.
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is DOUBLE PRECISION array, dimension (N)
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value
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
!> \ingroup complex16PTcomputational
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  The method used is described in Nicholas J. Higham, "Efficient
!>  Algorithms for Computing the Condition Number of a Tridiagonal
!>  Matrix", SIAM J. Sci. Stat. Comput., Vol. 7, No. 1, January 1986.
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE ZPTCON(N,D,E,Anorm,Rcond,Rwork,Info)
      USE F77KINDS                        
      USE S_IDAMAX
      USE S_XERBLA
      IMPLICIT NONE
!*--ZPTCON126
!
! PARAMETER definitions rewritten by SPAG
!
      REAL(R8KIND) , PARAMETER  ::  ONE = 1.0D+0 , ZERO = 0.0D+0
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER :: N
      REAL(R8KIND) , INTENT(IN) , DIMENSION(*) :: D
      COMPLEX(CX16KIND) , INTENT(IN) , DIMENSION(*) :: E
      REAL(R8KIND) , INTENT(IN) :: Anorm
      REAL(R8KIND) , INTENT(OUT) :: Rcond
      REAL(R8KIND) , INTENT(INOUT) , DIMENSION(*) :: Rwork
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      REAL(R8KIND) :: ainvnm
      INTEGER :: i , ix
!
! End of declarations rewritten by SPAG
!
!     ..
!     .. Array Arguments ..
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
!     ..
!     .. Local Scalars ..
!     ..
!     .. External Functions ..
!     ..
!     .. External Subroutines ..
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments.
!
      Info = 0
      IF ( N<0 ) THEN
         Info = -1
      ELSEIF ( Anorm<ZERO ) THEN
         Info = -4
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('ZPTCON',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      Rcond = ZERO
      IF ( N==0 ) THEN
         Rcond = ONE
         RETURN
      ELSEIF ( Anorm==ZERO ) THEN
         RETURN
      ENDIF
!
!     Check that D(1:N) is positive.
!
      DO i = 1 , N
         IF ( D(i)<=ZERO ) RETURN
      ENDDO
!
!     Solve M(A) * x = e, where M(A) = (m(i,j)) is given by
!
!        m(i,j) =  abs(A(i,j)), i = j,
!        m(i,j) = -abs(A(i,j)), i .ne. j,
!
!     and e = [ 1, 1, ..., 1 ]**T.  Note M(A) = M(L)*D*M(L)**H.
!
!     Solve M(L) * x = e.
!
      Rwork(1) = ONE
      DO i = 2 , N
         Rwork(i) = ONE + Rwork(i-1)*ABS(E(i-1))
      ENDDO
!
!     Solve D * M(L)**H * x = b.
!
      Rwork(N) = Rwork(N)/D(N)
      DO i = N - 1 , 1 , -1
         Rwork(i) = Rwork(i)/D(i) + Rwork(i+1)*ABS(E(i))
      ENDDO
!
!     Compute AINVNM = max(x(i)), 1<=i<=n.
!
      ix = IDAMAX(N,Rwork,1)
      ainvnm = ABS(Rwork(ix))
!
!     Compute the reciprocal condition number.
!
      IF ( ainvnm/=ZERO ) Rcond = (ONE/ainvnm)/Anorm
!
!
!     End of ZPTCON
!
      END SUBROUTINE ZPTCON
