!*==sptsv.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief <b> SPTSV computes the solution to system of linear equations A * X = B for PT matrices</b>
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SPTSV + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/sptsv.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/sptsv.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/sptsv.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SPTSV( N, NRHS, D, E, B, LDB, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, LDB, N, NRHS
!       ..
!       .. Array Arguments ..
!       REAL               B( LDB, * ), D( * ), E( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SPTSV computes the solution to a real system of linear equations
!> A*X = B, where A is an N-by-N symmetric positive definite tridiagonal
!> matrix, and X and B are N-by-NRHS matrices.
!>
!> A is factored as A = L*D*L**T, and the factored form of A is then
!> used to solve the system of equations.
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
!> \param[in] NRHS
!> \verbatim
!>          NRHS is INTEGER
!>          The number of right hand sides, i.e., the number of columns
!>          of the matrix B.  NRHS >= 0.
!> \endverbatim
!>
!> \param[in,out] D
!> \verbatim
!>          D is REAL array, dimension (N)
!>          On entry, the n diagonal elements of the tridiagonal matrix
!>          A.  On exit, the n diagonal elements of the diagonal matrix
!>          D from the factorization A = L*D*L**T.
!> \endverbatim
!>
!> \param[in,out] E
!> \verbatim
!>          E is REAL array, dimension (N-1)
!>          On entry, the (n-1) subdiagonal elements of the tridiagonal
!>          matrix A.  On exit, the (n-1) subdiagonal elements of the
!>          unit bidiagonal factor L from the L*D*L**T factorization of
!>          A.  (E can also be regarded as the superdiagonal of the unit
!>          bidiagonal factor U from the U**T*D*U factorization of A.)
!> \endverbatim
!>
!> \param[in,out] B
!> \verbatim
!>          B is REAL array, dimension (LDB,NRHS)
!>          On entry, the N-by-NRHS right hand side matrix B.
!>          On exit, if INFO = 0, the N-by-NRHS solution matrix X.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of the array B.  LDB >= max(1,N).
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value
!>          > 0:  if INFO = i, the leading minor of order i is not
!>                positive definite, and the solution has not been
!>                computed.  The factorization has not been completed
!>                unless i = N.
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
!> \ingroup realPTsolve
!
!  =====================================================================
      SUBROUTINE SPTSV(N,Nrhs,D,E,B,Ldb,Info)
      IMPLICIT NONE
!*--SPTSV118
!
!  -- LAPACK driver routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Info , Ldb , N , Nrhs
!     ..
!     .. Array Arguments ..
      REAL B(Ldb,*) , D(*) , E(*)
!     ..
!
!  =====================================================================
!
!     .. External Subroutines ..
      EXTERNAL SPTTRF , SPTTRS , XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      Info = 0
      IF ( N<0 ) THEN
         Info = -1
      ELSEIF ( Nrhs<0 ) THEN
         Info = -2
      ELSEIF ( Ldb<MAX(1,N) ) THEN
         Info = -6
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('SPTSV ',-Info)
         RETURN
      ENDIF
!
!     Compute the L*D*L**T (or U**T*D*U) factorization of A.
!
      CALL SPTTRF(N,D,E,Info)
!
!        Solve the system A*X = B, overwriting B with X.
!
      IF ( Info==0 ) CALL SPTTRS(N,Nrhs,D,E,B,Ldb,Info)
!
!     End of SPTSV
!
      END SUBROUTINE SPTSV
