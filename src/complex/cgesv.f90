!*==cgesv.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief <b> CGESV computes the solution to system of linear equations A * X = B for GE matrices (simple driver) </b>
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CGESV + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/cgesv.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/cgesv.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/cgesv.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CGESV( N, NRHS, A, LDA, IPIV, B, LDB, INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, LDA, LDB, N, NRHS
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * )
!       COMPLEX            A( LDA, * ), B( LDB, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CGESV computes the solution to a complex system of linear equations
!>    A * X = B,
!> where A is an N-by-N matrix and X and B are N-by-NRHS matrices.
!>
!> The LU decomposition with partial pivoting and row interchanges is
!> used to factor A as
!>    A = P * L * U,
!> where P is a permutation matrix, L is unit lower triangular, and U is
!> upper triangular.  The factored form of A is then used to solve the
!> system of equations A * X = B.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of linear equations, i.e., the order of the
!>          matrix A.  N >= 0.
!> \endverbatim
!>
!> \param[in] NRHS
!> \verbatim
!>          NRHS is INTEGER
!>          The number of right hand sides, i.e., the number of columns
!>          of the matrix B.  NRHS >= 0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is COMPLEX array, dimension (LDA,N)
!>          On entry, the N-by-N coefficient matrix A.
!>          On exit, the factors L and U from the factorization
!>          A = P*L*U; the unit diagonal elements of L are not stored.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,N).
!> \endverbatim
!>
!> \param[out] IPIV
!> \verbatim
!>          IPIV is INTEGER array, dimension (N)
!>          The pivot indices that define the permutation matrix P;
!>          row i of the matrix was interchanged with row IPIV(i).
!> \endverbatim
!>
!> \param[in,out] B
!> \verbatim
!>          B is COMPLEX array, dimension (LDB,NRHS)
!>          On entry, the N-by-NRHS matrix of right hand side matrix B.
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
!>          > 0:  if INFO = i, U(i,i) is exactly zero.  The factorization
!>                has been completed, but the factor U is exactly
!>                singular, so the solution could not be computed.
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
!> \date June 2017
!
!> \ingroup complexGEsolve
!
!  =====================================================================
      SUBROUTINE CGESV(N,Nrhs,A,Lda,Ipiv,B,Ldb,Info)
      IMPLICIT NONE
!*--CGESV126
!
!  -- LAPACK driver routine (version 3.7.1) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     June 2017
!
!     .. Scalar Arguments ..
      INTEGER Info , Lda , Ldb , N , Nrhs
!     ..
!     .. Array Arguments ..
      INTEGER Ipiv(*)
      COMPLEX A(Lda,*) , B(Ldb,*)
!     ..
!
!  =====================================================================
!
!     .. External Subroutines ..
      EXTERNAL CGETRF , CGETRS , XERBLA
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
      ELSEIF ( Lda<MAX(1,N) ) THEN
         Info = -4
      ELSEIF ( Ldb<MAX(1,N) ) THEN
         Info = -7
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('CGESV ',-Info)
         RETURN
      ENDIF
!
!     Compute the LU factorization of A.
!
      CALL CGETRF(N,N,A,Lda,Ipiv,Info)
!
!        Solve the system A*X = B, overwriting B with X.
!
      IF ( Info==0 ) CALL CGETRS('No transpose',N,Nrhs,A,Lda,Ipiv,B,Ldb,&
     &                           Info)
!
!     End of CGESV
!
      END SUBROUTINE CGESV
