!*==zsysv_aa_2stage.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief <b> ZSYSV_AA_2STAGE computes the solution to system of linear equations A * X = B for SY matrices</b>
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZSYSV_AA_2STAGE + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zsysv_aasen_2stage.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zsysv_aasen_2stage.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zsysv_aasen_2stage.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!      SUBROUTINE ZSYSV_AA_2STAGE( UPLO, N, NRHS, A, LDA, TB, LTB,
!                                  IPIV, IPIV2, B, LDB, WORK, LWORK,
!                                  INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          UPLO
!       INTEGER            N, NRHS, LDA, LTB, LDB, LWORK, INFO
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * ), IPIV2( * )
!       COMPLEX*16         A( LDA, * ), TB( * ), B( LDB, *), WORK( * )
!       ..
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZSYSV_AA_2STAGE computes the solution to a complex system of
!> linear equations
!>    A * X = B,
!> where A is an N-by-N symmetric matrix and X and B are N-by-NRHS
!> matrices.
!>
!> Aasen's 2-stage algorithm is used to factor A as
!>    A = U**T * T * U,  if UPLO = 'U', or
!>    A = L * T * L**T,  if UPLO = 'L',
!> where U (or L) is a product of permutation and unit upper (lower)
!> triangular matrices, and T is symmetric and band. The matrix T is
!> then LU-factored with partial pivoting. The factored form of A
!> is then used to solve the system of equations A * X = B.
!>
!> This is the blocked version of the algorithm, calling Level 3 BLAS.
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
!>          On exit, L is stored below (or above) the subdiaonal blocks,
!>          when UPLO  is 'L' (or 'U').
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,N).
!> \endverbatim
!>
!> \param[out] TB
!> \verbatim
!>          TB is COMPLEX*16 array, dimension (LTB)
!>          On exit, details of the LU factorization of the band matrix.
!> \endverbatim
!>
!> \param[in] LTB
!> \verbatim
!>          LTB is INTEGER
!>          The size of the array TB. LTB >= 4*N, internally
!>          used to select NB such that LTB >= (3*NB+1)*N.
!>
!>          If LTB = -1, then a workspace query is assumed; the
!>          routine only calculates the optimal size of LTB,
!>          returns this value as the first entry of TB, and
!>          no error message related to LTB is issued by XERBLA.
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
!> \param[out] IPIV2
!> \verbatim
!>          IPIV2 is INTEGER array, dimension (N)
!>          On exit, it contains the details of the interchanges, i.e.,
!>          the row and column k of T were interchanged with the
!>          row and column IPIV(k).
!> \endverbatim
!>
!> \param[in,out] B
!> \verbatim
!>          B is COMPLEX*16 array, dimension (LDB,NRHS)
!>          On entry, the right hand side matrix B.
!>          On exit, the solution matrix X.
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
!>          WORK is COMPLEX*16 workspace of size LWORK
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The size of WORK. LWORK >= N, internally used to select NB
!>          such that LWORK >= N*NB.
!>
!>          If LWORK = -1, then a workspace query is assumed; the
!>          routine only calculates the optimal size of the WORK array,
!>          returns this value as the first entry of the WORK array, and
!>          no error message related to LWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value.
!>          > 0:  if INFO = i, band LU factorization failed on i-th column
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
!> \ingroup complex16SYcomputational
!
!  =====================================================================
      SUBROUTINE ZSYSV_AA_2STAGE(Uplo,N,Nrhs,A,Lda,Tb,Ltb,Ipiv,Ipiv2,B, &
     &                           Ldb,Work,Lwork,Info)
!
!  -- LAPACK computational routine (version 3.8.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     November 2017
!
      USE F77KINDS                        
      USE S_LSAME
      USE S_XERBLA
      USE S_ZSYTRF_AA_2STAGE
      USE S_ZSYTRS_AA_2STAGE
      IMPLICIT NONE
!*--ZSYSV_AA_2STAGE200
!
! Dummy argument declarations rewritten by SPAG
!
      CHARACTER :: Uplo
      INTEGER :: N
      INTEGER :: Nrhs
      COMPLEX(CX16KIND) , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      COMPLEX(CX16KIND) , DIMENSION(*) :: Tb
      INTEGER :: Ltb
      INTEGER , DIMENSION(*) :: Ipiv
      INTEGER , DIMENSION(*) :: Ipiv2
      COMPLEX(CX16KIND) , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      COMPLEX(CX16KIND) , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      INTEGER :: lwkopt
      LOGICAL :: tquery , upper , wquery
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
      upper = LSAME(Uplo,'U')
      wquery = (Lwork==-1)
      tquery = (Ltb==-1)
      IF ( .NOT.upper .AND. .NOT.LSAME(Uplo,'L') ) THEN
         Info = -1
      ELSEIF ( N<0 ) THEN
         Info = -2
      ELSEIF ( Nrhs<0 ) THEN
         Info = -3
      ELSEIF ( Lda<MAX(1,N) ) THEN
         Info = -5
      ELSEIF ( Ltb<(4*N) .AND. .NOT.tquery ) THEN
         Info = -7
      ELSEIF ( Ldb<MAX(1,N) ) THEN
         Info = -11
      ELSEIF ( Lwork<N .AND. .NOT.wquery ) THEN
         Info = -13
      ENDIF
!
      IF ( Info==0 ) THEN
         CALL ZSYTRF_AA_2STAGE(Uplo,N,A,Lda,Tb,-1,Ipiv,Ipiv2,Work,-1,   &
     &                         Info)
         lwkopt = INT(Work(1))
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('ZSYSV_AA_2STAGE',-Info)
         RETURN
      ELSEIF ( wquery .OR. tquery ) THEN
         RETURN
      ENDIF
!
!
!     Compute the factorization A = U**T*T*U or A = L*T*L**T.
!
      CALL ZSYTRF_AA_2STAGE(Uplo,N,A,Lda,Tb,Ltb,Ipiv,Ipiv2,Work,Lwork,  &
     &                      Info)
!
!        Solve the system A*X = B, overwriting B with X.
!
!
      IF ( Info==0 ) CALL ZSYTRS_AA_2STAGE(Uplo,N,Nrhs,A,Lda,Tb,Ltb,    &
     &     Ipiv,Ipiv2,B,Ldb,Info)
!
      Work(1) = lwkopt
!
!     End of ZSYSV_AA_2STAGE
!
      END SUBROUTINE ZSYSV_AA_2STAGE
