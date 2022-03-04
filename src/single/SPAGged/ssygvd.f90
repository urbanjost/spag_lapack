!*==ssygvd.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b SSYGVD
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SSYGVD + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssygvd.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssygvd.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssygvd.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SSYGVD( ITYPE, JOBZ, UPLO, N, A, LDA, B, LDB, W, WORK,
!                          LWORK, IWORK, LIWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          JOBZ, UPLO
!       INTEGER            INFO, ITYPE, LDA, LDB, LIWORK, LWORK, N
!       ..
!       .. Array Arguments ..
!       INTEGER            IWORK( * )
!       REAL               A( LDA, * ), B( LDB, * ), W( * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SSYGVD computes all the eigenvalues, and optionally, the eigenvectors
!> of a real generalized symmetric-definite eigenproblem, of the form
!> A*x=(lambda)*B*x,  A*Bx=(lambda)*x,  or B*A*x=(lambda)*x.  Here A and
!> B are assumed to be symmetric and B is also positive definite.
!> If eigenvectors are desired, it uses a divide and conquer algorithm.
!>
!> The divide and conquer algorithm makes very mild assumptions about
!> floating point arithmetic. It will work on machines with a guard
!> digit in add/subtract, or on those binary machines without guard
!> digits which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or
!> Cray-2. It could conceivably fail on hexadecimal or decimal machines
!> without guard digits, but we know of none.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] ITYPE
!> \verbatim
!>          ITYPE is INTEGER
!>          Specifies the problem type to be solved:
!>          = 1:  A*x = (lambda)*B*x
!>          = 2:  A*B*x = (lambda)*x
!>          = 3:  B*A*x = (lambda)*x
!> \endverbatim
!>
!> \param[in] JOBZ
!> \verbatim
!>          JOBZ is CHARACTER*1
!>          = 'N':  Compute eigenvalues only;
!>          = 'V':  Compute eigenvalues and eigenvectors.
!> \endverbatim
!>
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          = 'U':  Upper triangles of A and B are stored;
!>          = 'L':  Lower triangles of A and B are stored.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrices A and B.  N >= 0.
!> \endverbatim
!>
!> \param[in,out] A
!> \verbatim
!>          A is REAL array, dimension (LDA, N)
!>          On entry, the symmetric matrix A.  If UPLO = 'U', the
!>          leading N-by-N upper triangular part of A contains the
!>          upper triangular part of the matrix A.  If UPLO = 'L',
!>          the leading N-by-N lower triangular part of A contains
!>          the lower triangular part of the matrix A.
!>
!>          On exit, if JOBZ = 'V', then if INFO = 0, A contains the
!>          matrix Z of eigenvectors.  The eigenvectors are normalized
!>          as follows:
!>          if ITYPE = 1 or 2, Z**T*B*Z = I;
!>          if ITYPE = 3, Z**T*inv(B)*Z = I.
!>          If JOBZ = 'N', then on exit the upper triangle (if UPLO='U')
!>          or the lower triangle (if UPLO='L') of A, including the
!>          diagonal, is destroyed.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= max(1,N).
!> \endverbatim
!>
!> \param[in,out] B
!> \verbatim
!>          B is REAL array, dimension (LDB, N)
!>          On entry, the symmetric matrix B.  If UPLO = 'U', the
!>          leading N-by-N upper triangular part of B contains the
!>          upper triangular part of the matrix B.  If UPLO = 'L',
!>          the leading N-by-N lower triangular part of B contains
!>          the lower triangular part of the matrix B.
!>
!>          On exit, if INFO <= N, the part of B containing the matrix is
!>          overwritten by the triangular factor U or L from the Cholesky
!>          factorization B = U**T*U or B = L*L**T.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of the array B.  LDB >= max(1,N).
!> \endverbatim
!>
!> \param[out] W
!> \verbatim
!>          W is REAL array, dimension (N)
!>          If INFO = 0, the eigenvalues in ascending order.
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
!>          The dimension of the array WORK.
!>          If N <= 1,               LWORK >= 1.
!>          If JOBZ = 'N' and N > 1, LWORK >= 2*N+1.
!>          If JOBZ = 'V' and N > 1, LWORK >= 1 + 6*N + 2*N**2.
!>
!>          If LWORK = -1, then a workspace query is assumed; the routine
!>          only calculates the optimal sizes of the WORK and IWORK
!>          arrays, returns these values as the first entries of the WORK
!>          and IWORK arrays, and no error message related to LWORK or
!>          LIWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] IWORK
!> \verbatim
!>          IWORK is INTEGER array, dimension (MAX(1,LIWORK))
!>          On exit, if INFO = 0, IWORK(1) returns the optimal LIWORK.
!> \endverbatim
!>
!> \param[in] LIWORK
!> \verbatim
!>          LIWORK is INTEGER
!>          The dimension of the array IWORK.
!>          If N <= 1,                LIWORK >= 1.
!>          If JOBZ  = 'N' and N > 1, LIWORK >= 1.
!>          If JOBZ  = 'V' and N > 1, LIWORK >= 3 + 5*N.
!>
!>          If LIWORK = -1, then a workspace query is assumed; the
!>          routine only calculates the optimal sizes of the WORK and
!>          IWORK arrays, returns these values as the first entries of
!>          the WORK and IWORK arrays, and no error message related to
!>          LWORK or LIWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value
!>          > 0:  SPOTRF or SSYEVD returned an error code:
!>             <= N:  if INFO = i and JOBZ = 'N', then the algorithm
!>                    failed to converge; i off-diagonal elements of an
!>                    intermediate tridiagonal form did not converge to
!>                    zero;
!>                    if INFO = i and JOBZ = 'V', then the algorithm
!>                    failed to compute an eigenvalue while working on
!>                    the submatrix lying in rows and columns INFO/(N+1)
!>                    through mod(INFO,N+1);
!>             > N:   if INFO = N + i, for 1 <= i <= N, then the leading
!>                    minor of order i of B is not positive definite.
!>                    The factorization of B could not be completed and
!>                    no eigenvalues or eigenvectors were computed.
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
!> \ingroup realSYeigen
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  Modified so that no backsubstitution is performed if SSYEVD fails to
!>  converge (NEIG in old code could be greater than N causing out of
!>  bounds reference to A - reported by Ralf Meyer).  Also corrected the
!>  description of INFO and the test on ITYPE. Sven, 16 Feb 05.
!> \endverbatim
!
!> \par Contributors:
!  ==================
!>
!>     Mark Fahey, Department of Mathematics, Univ. of Kentucky, USA
!>
!  =====================================================================
      SUBROUTINE SSYGVD(Itype,Jobz,Uplo,N,A,Lda,B,Ldb,W,Work,Lwork,     &
     &                  Iwork,Liwork,Info)
      USE S_LSAME
      USE S_SPOTRF
      USE S_SSYEVD
      USE S_SSYGST
      USE S_STRMM
      USE S_STRSM
      USE S_XERBLA
      IMPLICIT NONE
!*--SSYGVD238
!
! PARAMETER definitions rewritten by SPAG
!
      REAL , PARAMETER  ::  ONE = 1.0E+0
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER :: Itype
      CHARACTER :: Jobz
      CHARACTER :: Uplo
      INTEGER :: N
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL , DIMENSION(*) :: W
      REAL , INTENT(INOUT) , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , INTENT(INOUT) , DIMENSION(*) :: Iwork
      INTEGER :: Liwork
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      INTEGER :: liopt , liwmin , lopt , lwmin
      LOGICAL :: lquery , upper , wantz
      CHARACTER :: trans
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
!     Test the input parameters.
!
      wantz = LSAME(Jobz,'V')
      upper = LSAME(Uplo,'U')
      lquery = (Lwork==-1 .OR. Liwork==-1)
!
      Info = 0
      IF ( N<=1 ) THEN
         liwmin = 1
         lwmin = 1
      ELSEIF ( wantz ) THEN
         liwmin = 3 + 5*N
         lwmin = 1 + 6*N + 2*N**2
      ELSE
         liwmin = 1
         lwmin = 2*N + 1
      ENDIF
      lopt = lwmin
      liopt = liwmin
      IF ( Itype<1 .OR. Itype>3 ) THEN
         Info = -1
      ELSEIF ( .NOT.(wantz .OR. LSAME(Jobz,'N')) ) THEN
         Info = -2
      ELSEIF ( .NOT.(upper .OR. LSAME(Uplo,'L')) ) THEN
         Info = -3
      ELSEIF ( N<0 ) THEN
         Info = -4
      ELSEIF ( Lda<MAX(1,N) ) THEN
         Info = -6
      ELSEIF ( Ldb<MAX(1,N) ) THEN
         Info = -8
      ENDIF
!
      IF ( Info==0 ) THEN
         Work(1) = lopt
         Iwork(1) = liopt
!
         IF ( Lwork<lwmin .AND. .NOT.lquery ) THEN
            Info = -11
         ELSEIF ( Liwork<liwmin .AND. .NOT.lquery ) THEN
            Info = -13
         ENDIF
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('SSYGVD',-Info)
         RETURN
      ELSEIF ( lquery ) THEN
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( N==0 ) RETURN
!
!     Form a Cholesky factorization of B.
!
      CALL SPOTRF(Uplo,N,B,Ldb,Info)
      IF ( Info/=0 ) THEN
         Info = N + Info
         RETURN
      ENDIF
!
!     Transform problem to standard eigenvalue problem and solve.
!
      CALL SSYGST(Itype,Uplo,N,A,Lda,B,Ldb,Info)
      CALL SSYEVD(Jobz,Uplo,N,A,Lda,W,Work,Lwork,Iwork,Liwork,Info)
      lopt = MAX(REAL(lopt),REAL(Work(1)))
      liopt = MAX(REAL(liopt),REAL(Iwork(1)))
!
      IF ( wantz .AND. Info==0 ) THEN
!
!        Backtransform eigenvectors to the original problem.
!
         IF ( Itype==1 .OR. Itype==2 ) THEN
!
!           For A*x=(lambda)*B*x and A*B*x=(lambda)*x;
!           backtransform eigenvectors: x = inv(L)**T*y or inv(U)*y
!
            IF ( upper ) THEN
               trans = 'N'
            ELSE
               trans = 'T'
            ENDIF
!
            CALL STRSM('Left',Uplo,trans,'Non-unit',N,N,ONE,B,Ldb,A,Lda)
!
         ELSEIF ( Itype==3 ) THEN
!
!           For B*A*x=(lambda)*x;
!           backtransform eigenvectors: x = L*y or U**T*y
!
            IF ( upper ) THEN
               trans = 'T'
            ELSE
               trans = 'N'
            ENDIF
!
            CALL STRMM('Left',Uplo,trans,'Non-unit',N,N,ONE,B,Ldb,A,Lda)
         ENDIF
      ENDIF
!
      Work(1) = lopt
      Iwork(1) = liopt
!
!
!     End of SSYGVD
!
      END SUBROUTINE SSYGVD
