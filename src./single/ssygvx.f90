!*==ssygvx.f90  processed by SPAG 7.51RB at 21:36 on  3 Mar 2022
!> \brief \b SSYGVX
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download SSYGVX + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ssygvx.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ssygvx.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ssygvx.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE SSYGVX( ITYPE, JOBZ, RANGE, UPLO, N, A, LDA, B, LDB,
!                          VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ, WORK,
!                          LWORK, IWORK, IFAIL, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          JOBZ, RANGE, UPLO
!       INTEGER            IL, INFO, ITYPE, IU, LDA, LDB, LDZ, LWORK, M, N
!       REAL               ABSTOL, VL, VU
!       ..
!       .. Array Arguments ..
!       INTEGER            IFAIL( * ), IWORK( * )
!       REAL               A( LDA, * ), B( LDB, * ), W( * ), WORK( * ),
!      $                   Z( LDZ, * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> SSYGVX computes selected eigenvalues, and optionally, eigenvectors
!> of a real generalized symmetric-definite eigenproblem, of the form
!> A*x=(lambda)*B*x,  A*Bx=(lambda)*x,  or B*A*x=(lambda)*x.  Here A
!> and B are assumed to be symmetric and B is also positive definite.
!> Eigenvalues and eigenvectors can be selected by specifying either a
!> range of values or a range of indices for the desired eigenvalues.
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
!> \param[in] RANGE
!> \verbatim
!>          RANGE is CHARACTER*1
!>          = 'A': all eigenvalues will be found.
!>          = 'V': all eigenvalues in the half-open interval (VL,VU]
!>                 will be found.
!>          = 'I': the IL-th through IU-th eigenvalues will be found.
!> \endverbatim
!>
!> \param[in] UPLO
!> \verbatim
!>          UPLO is CHARACTER*1
!>          = 'U':  Upper triangle of A and B are stored;
!>          = 'L':  Lower triangle of A and B are stored.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The order of the matrix pencil (A,B).  N >= 0.
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
!>          On exit, the lower triangle (if UPLO='L') or the upper
!>          triangle (if UPLO='U') of A, including the diagonal, is
!>          destroyed.
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
!> \param[in] VL
!> \verbatim
!>          VL is REAL
!>          If RANGE='V', the lower bound of the interval to
!>          be searched for eigenvalues. VL < VU.
!>          Not referenced if RANGE = 'A' or 'I'.
!> \endverbatim
!>
!> \param[in] VU
!> \verbatim
!>          VU is REAL
!>          If RANGE='V', the upper bound of the interval to
!>          be searched for eigenvalues. VL < VU.
!>          Not referenced if RANGE = 'A' or 'I'.
!> \endverbatim
!>
!> \param[in] IL
!> \verbatim
!>          IL is INTEGER
!>          If RANGE='I', the index of the
!>          smallest eigenvalue to be returned.
!>          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0.
!>          Not referenced if RANGE = 'A' or 'V'.
!> \endverbatim
!>
!> \param[in] IU
!> \verbatim
!>          IU is INTEGER
!>          If RANGE='I', the index of the
!>          largest eigenvalue to be returned.
!>          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0.
!>          Not referenced if RANGE = 'A' or 'V'.
!> \endverbatim
!>
!> \param[in] ABSTOL
!> \verbatim
!>          ABSTOL is REAL
!>          The absolute error tolerance for the eigenvalues.
!>          An approximate eigenvalue is accepted as converged
!>          when it is determined to lie in an interval [a,b]
!>          of width less than or equal to
!>
!>                  ABSTOL + EPS *   max( |a|,|b| ) ,
!>
!>          where EPS is the machine precision.  If ABSTOL is less than
!>          or equal to zero, then  EPS*|T|  will be used in its place,
!>          where |T| is the 1-norm of the tridiagonal matrix obtained
!>          by reducing C to tridiagonal form, where C is the symmetric
!>          matrix of the standard symmetric problem to which the
!>          generalized problem is transformed.
!>
!>          Eigenvalues will be computed most accurately when ABSTOL is
!>          set to twice the underflow threshold 2*DLAMCH('S'), not zero.
!>          If this routine returns with INFO>0, indicating that some
!>          eigenvectors did not converge, try setting ABSTOL to
!>          2*SLAMCH('S').
!> \endverbatim
!>
!> \param[out] M
!> \verbatim
!>          M is INTEGER
!>          The total number of eigenvalues found.  0 <= M <= N.
!>          If RANGE = 'A', M = N, and if RANGE = 'I', M = IU-IL+1.
!> \endverbatim
!>
!> \param[out] W
!> \verbatim
!>          W is REAL array, dimension (N)
!>          On normal exit, the first M elements contain the selected
!>          eigenvalues in ascending order.
!> \endverbatim
!>
!> \param[out] Z
!> \verbatim
!>          Z is REAL array, dimension (LDZ, max(1,M))
!>          If JOBZ = 'N', then Z is not referenced.
!>          If JOBZ = 'V', then if INFO = 0, the first M columns of Z
!>          contain the orthonormal eigenvectors of the matrix A
!>          corresponding to the selected eigenvalues, with the i-th
!>          column of Z holding the eigenvector associated with W(i).
!>          The eigenvectors are normalized as follows:
!>          if ITYPE = 1 or 2, Z**T*B*Z = I;
!>          if ITYPE = 3, Z**T*inv(B)*Z = I.
!>
!>          If an eigenvector fails to converge, then that column of Z
!>          contains the latest approximation to the eigenvector, and the
!>          index of the eigenvector is returned in IFAIL.
!>          Note: the user must ensure that at least max(1,M) columns are
!>          supplied in the array Z; if RANGE = 'V', the exact value of M
!>          is not known in advance and an upper bound must be used.
!> \endverbatim
!>
!> \param[in] LDZ
!> \verbatim
!>          LDZ is INTEGER
!>          The leading dimension of the array Z.  LDZ >= 1, and if
!>          JOBZ = 'V', LDZ >= max(1,N).
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
!>          The length of the array WORK.  LWORK >= max(1,8*N).
!>          For optimal efficiency, LWORK >= (NB+3)*N,
!>          where NB is the blocksize for SSYTRD returned by ILAENV.
!>
!>          If LWORK = -1, then a workspace query is assumed; the routine
!>          only calculates the optimal size of the WORK array, returns
!>          this value as the first entry of the WORK array, and no error
!>          message related to LWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] IWORK
!> \verbatim
!>          IWORK is INTEGER array, dimension (5*N)
!> \endverbatim
!>
!> \param[out] IFAIL
!> \verbatim
!>          IFAIL is INTEGER array, dimension (N)
!>          If JOBZ = 'V', then if INFO = 0, the first M elements of
!>          IFAIL are zero.  If INFO > 0, then IFAIL contains the
!>          indices of the eigenvectors that failed to converge.
!>          If JOBZ = 'N', then IFAIL is not referenced.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value
!>          > 0:  SPOTRF or SSYEVX returned an error code:
!>             <= N:  if INFO = i, SSYEVX failed to converge;
!>                    i eigenvectors failed to converge.  Their indices
!>                    are stored in array IFAIL.
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
!> \date June 2016
!
!> \ingroup realSYeigen
!
!> \par Contributors:
!  ==================
!>
!>     Mark Fahey, Department of Mathematics, Univ. of Kentucky, USA
!
!  =====================================================================
      SUBROUTINE SSYGVX(Itype,Jobz,Range,Uplo,N,A,Lda,B,Ldb,Vl,Vu,Il,Iu,&
     &                  Abstol,M,W,Z,Ldz,Work,Lwork,Iwork,Ifail,Info)
      USE S_ILAENV
      USE S_LSAME
      USE S_SPOTRF
      USE S_SSYEVX
      USE S_SSYGST
      USE S_STRMM
      USE S_STRSM
      USE S_XERBLA
      IMPLICIT NONE
!*--SSYGVX308
!
! PARAMETER definitions rewritten by SPAG
!
      REAL , PARAMETER  ::  ONE = 1.0E+0
!
! Dummy argument declarations rewritten by SPAG
!
      INTEGER :: Itype
      CHARACTER :: Jobz
      CHARACTER :: Range
      CHARACTER :: Uplo
      INTEGER :: N
      REAL , DIMENSION(Lda,*) :: A
      INTEGER :: Lda
      REAL , DIMENSION(Ldb,*) :: B
      INTEGER :: Ldb
      REAL :: Vl
      REAL :: Vu
      INTEGER :: Il
      INTEGER :: Iu
      REAL :: Abstol
      INTEGER :: M
      REAL , DIMENSION(*) :: W
      REAL , DIMENSION(Ldz,*) :: Z
      INTEGER :: Ldz
      REAL , DIMENSION(*) :: Work
      INTEGER :: Lwork
      INTEGER , DIMENSION(*) :: Iwork
      INTEGER , DIMENSION(*) :: Ifail
      INTEGER , INTENT(INOUT) :: Info
!
! Local variable declarations rewritten by SPAG
!
      LOGICAL :: alleig , indeig , lquery , upper , valeig , wantz
      INTEGER :: lwkmin , lwkopt , nb
      CHARACTER :: trans
!
! End of declarations rewritten by SPAG
!
!     ..
!     .. Array Arguments ..
!     ..
!
! =====================================================================
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
      upper = LSAME(Uplo,'U')
      wantz = LSAME(Jobz,'V')
      alleig = LSAME(Range,'A')
      valeig = LSAME(Range,'V')
      indeig = LSAME(Range,'I')
      lquery = (Lwork==-1)
!
      Info = 0
      IF ( Itype<1 .OR. Itype>3 ) THEN
         Info = -1
      ELSEIF ( .NOT.(wantz .OR. LSAME(Jobz,'N')) ) THEN
         Info = -2
      ELSEIF ( .NOT.(alleig .OR. valeig .OR. indeig) ) THEN
         Info = -3
      ELSEIF ( .NOT.(upper .OR. LSAME(Uplo,'L')) ) THEN
         Info = -4
      ELSEIF ( N<0 ) THEN
         Info = -5
      ELSEIF ( Lda<MAX(1,N) ) THEN
         Info = -7
      ELSEIF ( Ldb<MAX(1,N) ) THEN
         Info = -9
      ELSEIF ( valeig ) THEN
         IF ( N>0 .AND. Vu<=Vl ) Info = -11
      ELSEIF ( indeig ) THEN
         IF ( Il<1 .OR. Il>MAX(1,N) ) THEN
            Info = -12
         ELSEIF ( Iu<MIN(N,Il) .OR. Iu>N ) THEN
            Info = -13
         ENDIF
      ENDIF
      IF ( Info==0 ) THEN
         IF ( Ldz<1 .OR. (wantz .AND. Ldz<N) ) Info = -18
      ENDIF
!
      IF ( Info==0 ) THEN
         lwkmin = MAX(1,8*N)
         nb = ILAENV(1,'SSYTRD',Uplo,N,-1,-1,-1)
         lwkopt = MAX(lwkmin,(nb+3)*N)
         Work(1) = lwkopt
!
         IF ( Lwork<lwkmin .AND. .NOT.lquery ) Info = -20
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('SSYGVX',-Info)
         RETURN
      ELSEIF ( lquery ) THEN
         RETURN
      ENDIF
!
!     Quick return if possible
!
      M = 0
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
      CALL SSYEVX(Jobz,Range,Uplo,N,A,Lda,Vl,Vu,Il,Iu,Abstol,M,W,Z,Ldz, &
     &            Work,Lwork,Iwork,Ifail,Info)
!
      IF ( wantz ) THEN
!
!        Backtransform eigenvectors to the original problem.
!
         IF ( Info>0 ) M = Info - 1
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
            CALL STRSM('Left',Uplo,trans,'Non-unit',N,M,ONE,B,Ldb,Z,Ldz)
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
            CALL STRMM('Left',Uplo,trans,'Non-unit',N,M,ONE,B,Ldb,Z,Ldz)
         ENDIF
      ENDIF
!
!     Set WORK(1) to optimal workspace size.
!
      Work(1) = lwkopt
!
!
!     End of SSYGVX
!
      END SUBROUTINE SSYGVX
