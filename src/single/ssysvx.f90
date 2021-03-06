!*==ssysvx.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief <b> SSYSVX computes the solution to system of linear equations A * X = B for SY matrices</b>
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SSYSVX + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssysvx.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssysvx.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssysvx.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SSYSVX( FACT, UPLO, N, NRHS, A, LDA, AF, LDAF, IPIV, B,
!                          LDB, X, LDX, RCOND, FERR, BERR, WORK, LWORK,
!                          IWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          FACT, UPLO
!       INTEGER            INFO, LDA, LDAF, LDB, LDX, LWORK, N, NRHS
!       REAL               RCOND
!       ..
!       .. Array Arguments ..
!       INTEGER            IPIV( * ), IWORK( * )
!       REAL               A( LDA, * ), AF( LDAF, * ), B( LDB, * ),
!      $                   BERR( * ), FERR( * ), WORK( * ), X( LDX, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SSYSVX uses the diagonal pivoting factorization to compute the
!> solution to a real system of linear equations A * X = B,
!> where A is an N-by-N symmetric matrix and X and B are N-by-NRHS
!> matrices.
!>
!> Error bounds on the solution and a condition estimate are also
!> provided.
!> \endverbatim
!
!> \par Description:
!  =================
!>
!> \verbatim
!>
!> The following steps are performed:
!>
!> 1. If FACT = 'N', the diagonal pivoting method is used to factor A.
!>    The form of the factorization is
!>       A = U * D * U**T,  if UPLO = 'U', or
!>       A = L * D * L**T,  if UPLO = 'L',
!>    where U (or L) is a product of permutation and unit upper (lower)
!>    triangular matrices, and D is symmetric and block diagonal with
!>    1-by-1 and 2-by-2 diagonal blocks.
!>
!> 2. If some D(i,i)=0, so that D is exactly singular, then the routine
!>    returns with INFO = i. Otherwise, the factored form of A is used
!>    to estimate the condition number of the matrix A.  If the
!>    reciprocal of the condition number is less than machine precision,
!>    INFO = N+1 is returned as a warning, but the routine still goes on
!>    to solve for X and compute error bounds as described below.
!>
!> 3. The system of equations is solved for X using the factored form
!>    of A.
!>
!> 4. Iterative refinement is applied to improve the computed solution
!>    matrix and calculate error bounds and backward error estimates
!>    for it.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] FACT
!> \verbatim
!>          FACT is CHARACTER*1
!>          Specifies whether or not the factored form of A has been
!>          supplied on entry.
!>          = 'F':  On entry, AF and IPIV contain the factored form of
!>                  A.  AF and IPIV will not be modified.
!>          = 'N':  The matrix A will be copied to AF and factored.
!> \endverbatim
!>
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
!>          of the matrices B and X.  NRHS >= 0.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is REAL array, dimension (LDA,N)
!>          The symmetric matrix A.  If UPLO = 'U', the leading N-by-N
!>          upper triangular part of A contains the upper triangular part
!>          of the matrix A, and the strictly lower triangular part of A
!>          is not referenced.  If UPLO = 'L', the leading N-by-N lower
!>          triangular part of A contains the lower triangular part of
!>          the matrix A, and the strictly upper triangular part of A is
!>          not referenced.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,N).
!> \endverbatim
!>
!> \param[in,out] AF
!> \verbatim
!>          AF is REAL array, dimension (LDAF,N)
!>          If FACT = 'F', then AF is an input argument and on entry
!>          contains the block diagonal matrix D and the multipliers used
!>          to obtain the factor U or L from the factorization
!>          A = U*D*U**T or A = L*D*L**T as computed by SSYTRF.
!>
!>          If FACT = 'N', then AF is an output argument and on exit
!>          returns the block diagonal matrix D and the multipliers used
!>          to obtain the factor U or L from the factorization
!>          A = U*D*U**T or A = L*D*L**T.
!> \endverbatim
!>
!> \param[in] LDAF
!> \verbatim
!>          LDAF is INTEGER
!>          The leading dimension of the array AF.  LDAF >= max(1,N).
!> \endverbatim
!>
!> \param[in,out] IPIV
!> \verbatim
!>          IPIV is INTEGER array, dimension (N)
!>          If FACT = 'F', then IPIV is an input argument and on entry
!>          contains details of the interchanges and the block structure
!>          of D, as determined by SSYTRF.
!>          If IPIV(k) > 0, then rows and columns k and IPIV(k) were
!>          interchanged and D(k,k) is a 1-by-1 diagonal block.
!>          If UPLO = 'U' and IPIV(k) = IPIV(k-1) < 0, then rows and
!>          columns k-1 and -IPIV(k) were interchanged and D(k-1:k,k-1:k)
!>          is a 2-by-2 diagonal block.  If UPLO = 'L' and IPIV(k) =
!>          IPIV(k+1) < 0, then rows and columns k+1 and -IPIV(k) were
!>          interchanged and D(k:k+1,k:k+1) is a 2-by-2 diagonal block.
!>
!>          If FACT = 'N', then IPIV is an output argument and on exit
!>          contains details of the interchanges and the block structure
!>          of D, as determined by SSYTRF.
!> \endverbatim
!>
!> \param[in] B
!> \verbatim
!>          B is REAL array, dimension (LDB,NRHS)
!>          The N-by-NRHS right hand side matrix B.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of the array B.  LDB >= max(1,N).
!> \endverbatim
!>
!> \param[out] X
!> \verbatim
!>          X is REAL array, dimension (LDX,NRHS)
!>          If INFO = 0 or INFO = N+1, the N-by-NRHS solution matrix X.
!> \endverbatim
!>
!> \param[in] LDX
!> \verbatim
!>          LDX is INTEGER
!>          The leading dimension of the array X.  LDX >= max(1,N).
!> \endverbatim
!>
!> \param[out] RCOND
!> \verbatim
!>          RCOND is REAL
!>          The estimate of the reciprocal condition number of the matrix
!>          A.  If RCOND is less than the machine precision (in
!>          particular, if RCOND = 0), the matrix is singular to working
!>          precision.  This condition is indicated by a return code of
!>          INFO > 0.
!> \endverbatim
!>
!> \param[out] FERR
!> \verbatim
!>          FERR is REAL array, dimension (NRHS)
!>          The estimated forward error bound for each solution vector
!>          X(j) (the j-th column of the solution matrix X).
!>          If XTRUE is the true solution corresponding to X(j), FERR(j)
!>          is an estimated upper bound for the magnitude of the largest
!>          element in (X(j) - XTRUE) divided by the magnitude of the
!>          largest element in X(j).  The estimate is as reliable as
!>          the estimate for RCOND, and is almost always a slight
!>          overestimate of the true error.
!> \endverbatim
!>
!> \param[out] BERR
!> \verbatim
!>          BERR is REAL array, dimension (NRHS)
!>          The componentwise relative backward error of each solution
!>          vector X(j) (i.e., the smallest relative change in
!>          any element of A or B that makes X(j) an exact solution).
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
!>          The length of WORK.  LWORK >= max(1,3*N), and for best
!>          performance, when FACT = 'N', LWORK >= max(1,3*N,N*NB), where
!>          NB is the optimal blocksize for SSYTRF.
!>
!>          If LWORK = -1, then a workspace query is assumed; the routine
!>          only calculates the optimal size of the WORK array, returns
!>          this value as the first entry of the WORK array, and no error
!>          message related to LWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] IWORK
!> \verbatim
!>          IWORK is INTEGER array, dimension (N)
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0: successful exit
!>          < 0: if INFO = -i, the i-th argument had an illegal value
!>          > 0: if INFO = i, and i is
!>                <= N:  D(i,i) is exactly zero.  The factorization
!>                       has been completed but the factor D is exactly
!>                       singular, so the solution and error bounds could
!>                       not be computed. RCOND = 0 is returned.
!>                = N+1: D is nonsingular, but RCOND is less than machine
!>                       precision, meaning that the matrix is singular
!>                       to working precision.  Nevertheless, the
!>                       solution and error bounds are computed because
!>                       there are a number of situations where the
!>                       computed solution can be more accurate than the
!>                       value of RCOND would suggest.
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
!> \date April 2012
!
!> \ingroup realSYsolve
!
!  =====================================================================
      SUBROUTINE SSYSVX(Fact,Uplo,N,Nrhs,A,Lda,Af,Ldaf,Ipiv,B,Ldb,X,Ldx,&
     &                  Rcond,Ferr,Berr,Work,Lwork,Iwork,Info)
      IMPLICIT NONE
!*--SSYSVX287
!
!  -- LAPACK driver routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     April 2012
!
!     .. Scalar Arguments ..
      CHARACTER Fact , Uplo
      INTEGER Info , Lda , Ldaf , Ldb , Ldx , Lwork , N , Nrhs
      REAL Rcond
!     ..
!     .. Array Arguments ..
      INTEGER Ipiv(*) , Iwork(*)
      REAL A(Lda,*) , Af(Ldaf,*) , B(Ldb,*) , Berr(*) , Ferr(*) ,       &
     &     Work(*) , X(Ldx,*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ZERO
      PARAMETER (ZERO=0.0E+0)
!     ..
!     .. Local Scalars ..
      LOGICAL lquery , nofact
      INTEGER lwkopt , nb
      REAL anorm
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      INTEGER ILAENV
      REAL SLAMCH , SLANSY
      EXTERNAL ILAENV , LSAME , SLAMCH , SLANSY
!     ..
!     .. External Subroutines ..
      EXTERNAL SLACPY , SSYCON , SSYRFS , SSYTRF , SSYTRS , XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      Info = 0
      nofact = LSAME(Fact,'N')
      lquery = (Lwork==-1)
      IF ( .NOT.nofact .AND. .NOT.LSAME(Fact,'F') ) THEN
         Info = -1
      ELSEIF ( .NOT.LSAME(Uplo,'U') .AND. .NOT.LSAME(Uplo,'L') ) THEN
         Info = -2
      ELSEIF ( N<0 ) THEN
         Info = -3
      ELSEIF ( Nrhs<0 ) THEN
         Info = -4
      ELSEIF ( Lda<MAX(1,N) ) THEN
         Info = -6
      ELSEIF ( Ldaf<MAX(1,N) ) THEN
         Info = -8
      ELSEIF ( Ldb<MAX(1,N) ) THEN
         Info = -11
      ELSEIF ( Ldx<MAX(1,N) ) THEN
         Info = -13
      ELSEIF ( Lwork<MAX(1,3*N) .AND. .NOT.lquery ) THEN
         Info = -18
      ENDIF
!
      IF ( Info==0 ) THEN
         lwkopt = MAX(1,3*N)
         IF ( nofact ) THEN
            nb = ILAENV(1,'SSYTRF',Uplo,N,-1,-1,-1)
            lwkopt = MAX(lwkopt,N*nb)
         ENDIF
         Work(1) = lwkopt
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('SSYSVX',-Info)
         RETURN
      ELSEIF ( lquery ) THEN
         RETURN
      ENDIF
!
      IF ( nofact ) THEN
!
!        Compute the factorization A = U*D*U**T or A = L*D*L**T.
!
         CALL SLACPY(Uplo,N,N,A,Lda,Af,Ldaf)
         CALL SSYTRF(Uplo,N,Af,Ldaf,Ipiv,Work,Lwork,Info)
!
!        Return if INFO is non-zero.
!
         IF ( Info>0 ) THEN
            Rcond = ZERO
            RETURN
         ENDIF
      ENDIF
!
!     Compute the norm of the matrix A.
!
      anorm = SLANSY('I',Uplo,N,A,Lda,Work)
!
!     Compute the reciprocal of the condition number of A.
!
      CALL SSYCON(Uplo,N,Af,Ldaf,Ipiv,anorm,Rcond,Work,Iwork,Info)
!
!     Compute the solution vectors X.
!
      CALL SLACPY('Full',N,Nrhs,B,Ldb,X,Ldx)
      CALL SSYTRS(Uplo,N,Nrhs,Af,Ldaf,Ipiv,X,Ldx,Info)
!
!     Use iterative refinement to improve the computed solutions and
!     compute error bounds and backward error estimates for them.
!
      CALL SSYRFS(Uplo,N,Nrhs,A,Lda,Af,Ldaf,Ipiv,B,Ldb,X,Ldx,Ferr,Berr, &
     &            Work,Iwork,Info)
!
!     Set INFO = N+1 if the matrix is singular to working precision.
!
      IF ( Rcond<SLAMCH('Epsilon') ) Info = N + 1
!
      Work(1) = lwkopt
!
!
!     End of SSYSVX
!
      END SUBROUTINE SSYSVX
