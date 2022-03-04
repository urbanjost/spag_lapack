!*==sptcon.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b SPTCON
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SPTCON + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sptcon.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sptcon.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sptcon.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SPTCON( N, D, E, ANORM, RCOND, WORK, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, N
!       REAL               ANORM, RCOND
!       ..
!       .. Array Arguments ..
!       REAL               D( * ), E( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SPTCON computes the reciprocal of the condition number (in the
!> 1-norm) of a real symmetric positive definite tridiagonal matrix
!> using the factorization A = L*D*L**T or A = U**T*D*U computed by
!> SPTTRF.
!>
!> Norm(inv(A)) is computed by a direct method, and the reciprocal of
!> the condition number is computed as
!>              RCOND = 1 / (ANORM * norm(inv(A))).
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
!>          D is REAL array, dimension (N)
!>          The n diagonal elements of the diagonal matrix D from the
!>          factorization of A, as computed by SPTTRF.
!> \endverbatim
!>
!> \param[in] E
!> \verbatim
!>          E is REAL array, dimension (N-1)
!>          The (n-1) off-diagonal elements of the unit bidiagonal factor
!>          U or L from the factorization of A,  as computed by SPTTRF.
!> \endverbatim
!>
!> \param[in] ANORM
!> \verbatim
!>          ANORM is REAL
!>          The 1-norm of the original matrix A.
!> \endverbatim
!>
!> \param[out] RCOND
!> \verbatim
!>          RCOND is REAL
!>          The reciprocal of the condition number of the matrix A,
!>          computed as RCOND = 1/(ANORM * AINVNM), where AINVNM is the
!>          1-norm of inv(A) computed in this routine.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is REAL array, dimension (N)
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
!> \ingroup realPTcomputational
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
      SUBROUTINE SPTCON(N,D,E,Anorm,Rcond,Work,Info)
      IMPLICIT NONE
!*--SPTCON122
!
!  -- LAPACK computational routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Info , N
      REAL Anorm , Rcond
!     ..
!     .. Array Arguments ..
      REAL D(*) , E(*) , Work(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ONE , ZERO
      PARAMETER (ONE=1.0E+0,ZERO=0.0E+0)
!     ..
!     .. Local Scalars ..
      INTEGER i , ix
      REAL ainvnm
!     ..
!     .. External Functions ..
      INTEGER ISAMAX
      EXTERNAL ISAMAX
!     ..
!     .. External Subroutines ..
      EXTERNAL XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC ABS
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
         CALL XERBLA('SPTCON',-Info)
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
!     and e = [ 1, 1, ..., 1 ]**T.  Note M(A) = M(L)*D*M(L)**T.
!
!     Solve M(L) * x = e.
!
      Work(1) = ONE
      DO i = 2 , N
         Work(i) = ONE + Work(i-1)*ABS(E(i-1))
      ENDDO
!
!     Solve D * M(L)**T * x = b.
!
      Work(N) = Work(N)/D(N)
      DO i = N - 1 , 1 , -1
         Work(i) = Work(i)/D(i) + Work(i+1)*ABS(E(i))
      ENDDO
!
!     Compute AINVNM = max(x(i)), 1<=i<=n.
!
      ix = ISAMAX(N,Work,1)
      ainvnm = ABS(Work(ix))
!
!     Compute the reciprocal condition number.
!
      IF ( ainvnm/=ZERO ) Rcond = (ONE/ainvnm)/Anorm
!
!
!     End of SPTCON
!
      END SUBROUTINE SPTCON
