!*==zppsv.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief <b> ZPPSV computes the solution to system of linear equations A * X = B for OTHER matrices</b>
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZPPSV + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zppsv.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zppsv.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zppsv.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZPPSV( UPLO, N, NRHS, AP, B, LDB, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            INFO, LDB, N, NRHS
!       ..
!       .. Array Arguments ..
!       COMPLEX*16         AP( * ), B( LDB, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZPPSV computes the solution to a complex system of linear equations
!>    A * X = B,
!> where A is an N-by-N Hermitian positive definite matrix stored in
!> packed format and X and B are N-by-NRHS matrices.
!>
!> The Cholesky decomposition is used to factor A as
!>    A = U**H * U,  if UPLO = 'U', or
!>    A = L * L**H,  if UPLO = 'L',
!> where U is an upper triangular matrix and L is a lower triangular
!> matrix.  The factored form of A is then used to solve the system of
!> equations A * X = B.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          = 'U':  Upper triangle of A is stored;
!>          = 'L':  Lower triangle of A is stored.
!> \endverbatim
!>
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
!> \param[in,out] AP
!> \verbatim
!>          AP is COMPLEX*16 array, dimension (N*(N+1)/2)
!>          On entry, the upper or lower triangle of the Hermitian matrix
!>          A, packed columnwise in a linear array.  The j-th column of A
!>          is stored in the array AP as follows:
!>          if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
!>          if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.
!>          See below for further details.
!>
!>          On exit, if INFO = 0, the factor U or L from the Cholesky
!>          factorization A = U**H*U or A = L*L**H, in the same storage
!>          format as A.
!> \endverbatim
!>
!> \param[in,out] B
!> \verbatim
!>          B is COMPLEX*16 array, dimension (LDB,NRHS)
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
!>          > 0:  if INFO = i, the leading minor of order i of A is not
!>                positive definite, so the factorization could not be
!>                completed, and the solution has not been computed.
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
!> \ingroup complex16OTHERsolve
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  The packed storage scheme is illustrated by the following example
!>  when N = 4, UPLO = 'U':
!>
!>  Two-dimensional storage of the Hermitian matrix A:
!>
!>     a11 a12 a13 a14
!>         a22 a23 a24
!>             a33 a34     (aij = conjg(aji))
!>                 a44
!>
!>  Packed storage of the upper triangle of A:
!>
!>  AP = [ a11, a12, a22, a13, a23, a33, a14, a24, a34, a44 ]
!> \endverbatim
!>
!  =====================================================================
      SUBROUTINE ZPPSV(Uplo,N,Nrhs,Ap,B,Ldb,Info)
      USE F77KINDS                        
      USE S_LSAME
      USE S_XERBLA
      USE S_ZPPTRF
      USE S_ZPPTRS
      IMPLICIT NONE
!*--ZPPSV153
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Nrhs
      COMPLEX(CX16KIND) , DIMENSION(*) :: Ap
      COMPLEX(CX16KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      INTEGER , INTENT(INOUT) :: Info
!
! End of declarations rewritten by SPAG
!
!     ..
!     .. Array Arguments ..
!     ..
!
!  =====================================================================
!
!     .. External Functions ..
!     ..
!     .. External Subroutines ..
!     ..
!     .. Intrinsic Functions ..
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      Info = 0
      IF ( .NOT.LSAME(Uplo,'U') .AND. .NOT.LSAME(Uplo,'L') ) THEN
         Info = -1
      ELSEIF ( N<0 ) THEN
         Info = -2
      ELSEIF ( Nrhs<0 ) THEN
         Info = -3
      ELSEIF ( Ldb<MAX(1,N) ) THEN
         Info = -6
      ENDIF
      IF ( Info/=0 ) THEN
         CALL XERBLA('ZPPSV ',-Info)
         RETURN
      ENDIF
!
!     Compute the Cholesky factorization A = U**H *U or A = L*L**H.
!
      CALL ZPPTRF(Uplo,N,Ap,Info)
!
!        Solve the system A*X = B, overwriting B with X.
!
!
      IF ( Info==0 ) CALL ZPPTRS(Uplo,N,Nrhs,Ap,B,Ldb,Info)
!
!     End of ZPPSV
!
      END SUBROUTINE ZPPSV