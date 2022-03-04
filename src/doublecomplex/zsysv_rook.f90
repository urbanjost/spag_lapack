!*==zsysv_rook.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief <b> ZSYSV_ROOK computes the solution to system of linear equations A * X = B for SY matrices</b>
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZSYSV_ROOK + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zsysv_rook.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zsysv_rook.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zsysv_rook.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZSYSV_ROOK( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK,
!                         LWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            INFO, LDA, LDB, LWORK, N, NRHS
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * )
!       COMPLEX*16         A( LDA, * ), B( LDB, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZSYSV_ROOK computes the solution to a complex system of linear
!> equations
!>    A * X = B,
!> where A is an N-by-N symmetric matrix and X and B are N-by-NRHS
!> matrices.
!>
!> The diagonal pivoting method is used to factor A as
!>    A = U * D * U**T,  if UPLO = 'U', or
!>    A = L * D * L**T,  if UPLO = 'L',
!> where U (or L) is a product of permutation and unit upper (lower)
!> triangular matrices, and D is symmetric and block diagonal with
!> 1-by-1 and 2-by-2 diagonal blocks.
!>
!> ZSYTRF_ROOK is called to compute the factorization of a complex
!> symmetric matrix A using the bounded Bunch-Kaufman ("rook") diagonal
!> pivoting method.
!>
!> The factored form of A is then used to solve the system
!> of equations A * X = B by calling ZSYTRS_ROOK.
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
!> \param[in,out] A
!> \verbatim
!>          A is COMPLEX*16 array, dimension (LDA,N)
!>          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
!>          N-by-N upper triangular part of A contains the upper
!>          triangular part of the matrix A, and the strictly lower
!>          triangular part of A is not referenced.  If UPLO = 'L', the
!>          leading N-by-N lower triangular part of A contains the lower
!>          triangular part of the matrix A, and the strictly upper
!>          triangular part of A is not referenced.
!>
!>          On exit, if INFO = 0, the block diagonal matrix D and the
!>          multipliers used to obtain the factor U or L from the
!>          factorization A = U*D*U**T or A = L*D*L**T as computed by
!>          ZSYTRF_ROOK.
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
!>          Details of the interchanges and the block structure of D,
!>          as determined by ZSYTRF_ROOK.
!>
!>          If UPLO = 'U':
!>               If IPIV(k) > 0, then rows and columns k and IPIV(k)
!>               were interchanged and D(k,k) is a 1-by-1 diagonal block.
!>
!>               If IPIV(k) < 0 and IPIV(k-1) < 0, then rows and
!>               columns k and -IPIV(k) were interchanged and rows and
!>               columns k-1 and -IPIV(k-1) were inerchaged,
!>               D(k-1:k,k-1:k) is a 2-by-2 diagonal block.
!>
!>          If UPLO = 'L':
!>               If IPIV(k) > 0, then rows and columns k and IPIV(k)
!>               were interchanged and D(k,k) is a 1-by-1 diagonal block.
!>
!>               If IPIV(k) < 0 and IPIV(k+1) < 0, then rows and
!>               columns k and -IPIV(k) were interchanged and rows and
!>               columns k+1 and -IPIV(k+1) were inerchaged,
!>               D(k:k+1,k:k+1) is a 2-by-2 diagonal block.
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
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX*16 array, dimension (MAX(1,LWORK))
!>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The length of WORK.  LWORK >= 1, and for best performance
!>          LWORK >= max(1,N*NB), where NB is the optimal blocksize for
!>          ZSYTRF_ROOK.
!>
!>          TRS will be done with Level 2 BLAS
!>
!>          If LWORK = -1, then a workspace query is assumed; the routine
!>          only calculates the optimal size of the WORK array, returns
!>          this value as the first entry of the WORK array, and no error
!>          message related to LWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0: successful exit
!>          < 0: if INFO = -i, the i-th argument had an illegal value
!>          > 0: if INFO = i, D(i,i) is exactly zero.  The factorization
!>               has been completed, but the block diagonal matrix D is
!>               exactly singular, so the solution could not be computed.
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
!> \ingroup complex16SYsolve
!
!> \par Contributors:
!  ==================
!>
!> \verbatim
!>
!>   December 2016, Igor Kozachenko,
!>                  Computer Science Division,
!>                  University of California, Berkeley
!>
!>  September 2007, Sven Hammarling, Nicholas J. Higham, Craig Lucas,
!>                  School of Mathematics,
!>                  University of Manchester
!>
!> \endverbatim
!
!  =====================================================================
      SUBROUTINE ZSYSV_ROOK(Uplo,N,Nrhs,A,Lda,Ipiv,B,Ldb,Work,Lwork,    &
     &                      Info)
      IMPLICIT NONE
!*--ZSYSV_ROOK208
!
!  -- LAPACK driver routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Uplo
      INTEGER Info , Lda , Ldb , Lwork , N , Nrhs
!     ..
!     .. Array Arguments ..
      INTEGER Ipiv(*)
      COMPLEX*16 A(Lda,*) , B(Ldb,*) , Work(*)
!     ..
!
!  =====================================================================
!
!     .. Local Scalars ..
      LOGICAL lquery
      INTEGER lwkopt
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL XERBLA , ZSYTRF_ROOK , ZSYTRS_ROOK
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      Info = 0
      lquery = (Lwork==-1)
      IF ( .NOT.LSAME(Uplo,'U') .AND. .NOT.LSAME(Uplo,'L') ) THEN
         Info = -1
      ELSEIF ( N<0 ) THEN
         Info = -2
      ELSEIF ( Nrhs<0 ) THEN
         Info = -3
      ELSEIF ( Lda<MAX(1,N) ) THEN
         Info = -5
      ELSEIF ( Ldb<MAX(1,N) ) THEN
         Info = -8
      ELSEIF ( Lwork<1 .AND. .NOT.lquery ) THEN
         Info = -10
      ENDIF
!
      IF ( Info==0 ) THEN
         IF ( N==0 ) THEN
            lwkopt = 1
         ELSE
            CALL ZSYTRF_ROOK(Uplo,N,A,Lda,Ipiv,Work,-1,Info)
            lwkopt = Work(1)
         ENDIF
         Work(1) = lwkopt
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('ZSYSV_ROOK ',-Info)
         RETURN
      ELSEIF ( lquery ) THEN
         RETURN
      ENDIF
!
!     Compute the factorization A = U*D*U**T or A = L*D*L**T.
!
      CALL ZSYTRF_ROOK(Uplo,N,A,Lda,Ipiv,Work,Lwork,Info)
!
!        Solve the system A*X = B, overwriting B with X.
!
!        Solve with TRS_ROOK ( Use Level 2 BLAS)
!
!
      IF ( Info==0 ) CALL ZSYTRS_ROOK(Uplo,N,Nrhs,A,Lda,Ipiv,B,Ldb,Info)
!
      Work(1) = lwkopt
!
!
!     End of ZSYSV_ROOK
!
      END SUBROUTINE ZSYSV_ROOK
