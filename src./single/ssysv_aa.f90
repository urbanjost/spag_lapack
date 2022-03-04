!*==ssysv_aa.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief <b> SSYSV_AA computes the solution to system of linear equations A * X = B for SY matrices</b>
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SSYSV_AA + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssysv_aa.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssysv_aa.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssysv_aa.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SSYSV_AA( UPLO, N, NRHS, A, LDA, IPIV, B, LDB, WORK,
!                            LWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            N, NRHS, LDA, LDB, LWORK, INFO
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * )
!       REAL               A( LDA, * ), B( LDB, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SSYSV computes the solution to a real system of linear equations
!>    A * X = B,
!> where A is an N-by-N symmetric matrix and X and B are N-by-NRHS
!> matrices.
!>
!> Aasen's algorithm is used to factor A as
!>    A = U**T * T * U,  if UPLO = 'U', or
!>    A = L * T * L**T,  if UPLO = 'L',
!> where U (or L) is a product of permutation and unit upper (lower)
!> triangular matrices, and T is symmetric tridiagonal. The factored
!> form of A is then used to solve the system of equations A * X = B.
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
!>          A is REAL array, dimension (LDA,N)
!>          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
!>          N-by-N upper triangular part of A contains the upper
!>          triangular part of the matrix A, and the strictly lower
!>          triangular part of A is not referenced.  If UPLO = 'L', the
!>          leading N-by-N lower triangular part of A contains the lower
!>          triangular part of the matrix A, and the strictly upper
!>          triangular part of A is not referenced.
!>
!>          On exit, if INFO = 0, the tridiagonal matrix T and the
!>          multipliers used to obtain the factor U or L from the
!>          factorization A = U**T*T*U or A = L*T*L**T as computed by
!>          SSYTRF.
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
!>          On exit, it contains the details of the interchanges, i.e.,
!>          the row and column k of A were interchanged with the
!>          row and column IPIV(k).
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
!> \param[out] WORK
!> \verbatim
!>          WORK is REAL array, dimension (MAX(1,LWORK))
!>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The length of WORK.  LWORK >= MAX(1,2*N,3*N-2), and for
!>          the best performance, LWORK >= MAX(1,N*NB), where NB is
!>          the optimal blocksize for SSYTRF_AA.
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
!> \date November 2017
!
!> \ingroup realSYsolve
!
!  =====================================================================
      SUBROUTINE SSYSV_AA(Uplo,N,Nrhs,A,Lda,Ipiv,B,Ldb,Work,Lwork,Info)
      USE S_LSAME
      USE S_SSYTRF_AA
      USE S_SSYTRS_AA
      USE S_XERBLA
      IMPLICIT NONE
!*--SSYSV_AA169
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Nrhs
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      INTEGER , DIMENSION(*) :: Ipiv
      REAL , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      LOGICAL :: lquery
      INTEGER :: lwkopt , lwkopt_sytrf , lwkopt_sytrs
!
! End of declarations rewritten by SPAG
!
!     ..
!     .. Array Arguments ..
!     ..
!
!  =====================================================================
!
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
      ELSEIF ( Lwork<MAX(2*N,3*N-2) .AND. .NOT.lquery ) THEN
         Info = -10
      ENDIF
!
      IF ( Info==0 ) THEN
         CALL SSYTRF_AA(Uplo,N,A,Lda,Ipiv,Work,-1,Info)
         lwkopt_sytrf = INT(Work(1))
         CALL SSYTRS_AA(Uplo,N,Nrhs,A,Lda,Ipiv,B,Ldb,Work,-1,Info)
         lwkopt_sytrs = INT(Work(1))
         lwkopt = MAX(lwkopt_sytrf,lwkopt_sytrs)
         Work(1) = lwkopt
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('SSYSV_AA',-Info)
         RETURN
      ELSEIF ( lquery ) THEN
         RETURN
      ENDIF
!
!     Compute the factorization A = U**T*T*U or A = L*T*L**T.
!
      CALL SSYTRF_AA(Uplo,N,A,Lda,Ipiv,Work,Lwork,Info)
!
!        Solve the system A*X = B, overwriting B with X.
!
!
      IF ( Info==0 ) CALL SSYTRS_AA(Uplo,N,Nrhs,A,Lda,Ipiv,B,Ldb,Work,  &
     &                              Lwork,Info)
!
      Work(1) = lwkopt
!
!
!     End of SSYSV_AA
!
      END SUBROUTINE SSYSV_AA
