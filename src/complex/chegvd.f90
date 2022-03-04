!*==chegvd.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b CHEGVD
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download CHEGVD + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/chegvd.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/chegvd.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/chegvd.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE CHEGVD( ITYPE, JOBZ, UPLO, N, A, LDA, B, LDB, W, WORK,
!                          LWORK, RWORK, LRWORK, IWORK, LIWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          JOBZ, UPLO
!       INTEGER            INFO, ITYPE, LDA, LDB, LIWORK, LRWORK, LWORK, N
!       ..
!       .. Array Arguments ..
!       INTEGER            IWORK( * )
!       REAL               RWORK( * ), W( * )
!       COMPLEX            A( LDA, * ), B( LDB, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> CHEGVD computes all the eigenvalues, and optionally, the eigenvectors
!> of a complex generalized Hermitian-definite eigenproblem, of the form
!> A*x=(lambda)*B*x,  A*Bx=(lambda)*x,  or B*A*x=(lambda)*x.  Here A and
!> B are assumed to be Hermitian and B is also positive definite.
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
!>          A is COMPLEX array, dimension (LDA, N)
!>          On entry, the Hermitian matrix A.  If UPLO = 'U', the
!>          leading N-by-N upper triangular part of A contains the
!>          upper triangular part of the matrix A.  If UPLO = 'L',
!>          the leading N-by-N lower triangular part of A contains
!>          the lower triangular part of the matrix A.
!>
!>          On exit, if JOBZ = 'V', then if INFO = 0, A contains the
!>          matrix Z of eigenvectors.  The eigenvectors are normalized
!>          as follows:
!>          if ITYPE = 1 or 2, Z**H*B*Z = I;
!>          if ITYPE = 3, Z**H*inv(B)*Z = I.
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
!>          B is COMPLEX array, dimension (LDB, N)
!>          On entry, the Hermitian matrix B.  If UPLO = 'U', the
!>          leading N-by-N upper triangular part of B contains the
!>          upper triangular part of the matrix B.  If UPLO = 'L',
!>          the leading N-by-N lower triangular part of B contains
!>          the lower triangular part of the matrix B.
!>
!>          On exit, if INFO <= N, the part of B containing the matrix is
!>          overwritten by the triangular factor U or L from the Cholesky
!>          factorization B = U**H*U or B = L*L**H.
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
!>          WORK is COMPLEX array, dimension (MAX(1,LWORK))
!>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The length of the array WORK.
!>          If N <= 1,                LWORK >= 1.
!>          If JOBZ  = 'N' and N > 1, LWORK >= N + 1.
!>          If JOBZ  = 'V' and N > 1, LWORK >= 2*N + N**2.
!>
!>          If LWORK = -1, then a workspace query is assumed; the routine
!>          only calculates the optimal sizes of the WORK, RWORK and
!>          IWORK arrays, returns these values as the first entries of
!>          the WORK, RWORK and IWORK arrays, and no error message
!>          related to LWORK or LRWORK or LIWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is REAL array, dimension (MAX(1,LRWORK))
!>          On exit, if INFO = 0, RWORK(1) returns the optimal LRWORK.
!> \endverbatim
!>
!> \param[in] LRWORK
!> \verbatim
!>          LRWORK is INTEGER
!>          The dimension of the array RWORK.
!>          If N <= 1,                LRWORK >= 1.
!>          If JOBZ  = 'N' and N > 1, LRWORK >= N.
!>          If JOBZ  = 'V' and N > 1, LRWORK >= 1 + 5*N + 2*N**2.
!>
!>          If LRWORK = -1, then a workspace query is assumed; the
!>          routine only calculates the optimal sizes of the WORK, RWORK
!>          and IWORK arrays, returns these values as the first entries
!>          of the WORK, RWORK and IWORK arrays, and no error message
!>          related to LWORK or LRWORK or LIWORK is issued by XERBLA.
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
!>          routine only calculates the optimal sizes of the WORK, RWORK
!>          and IWORK arrays, returns these values as the first entries
!>          of the WORK, RWORK and IWORK arrays, and no error message
!>          related to LWORK or LRWORK or LIWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value
!>          > 0:  CPOTRF or CHEEVD returned an error code:
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
!> \ingroup complexHEeigen
!
!> \par Further Details:
!  =====================
!>
!> \verbatim
!>
!>  Modified so that no backsubstitution is performed if CHEEVD fails to
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
      SUBROUTINE CHEGVD(Itype,Jobz,Uplo,N,A,Lda,B,Ldb,W,Work,Lwork,     &
     &                  Rwork,Lrwork,Iwork,Liwork,Info)
      IMPLICIT NONE
!*--CHEGVD253
!
!  -- LAPACK driver routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Jobz , Uplo
      INTEGER Info , Itype , Lda , Ldb , Liwork , Lrwork , Lwork , N
!     ..
!     .. Array Arguments ..
      INTEGER Iwork(*)
      REAL Rwork(*) , W(*)
      COMPLEX A(Lda,*) , B(Ldb,*) , Work(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      COMPLEX CONE
      PARAMETER (CONE=(1.0E+0,0.0E+0))
!     ..
!     .. Local Scalars ..
      LOGICAL lquery , upper , wantz
      CHARACTER trans
      INTEGER liopt , liwmin , lopt , lropt , lrwmin , lwmin
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
!     ..
!     .. External Subroutines ..
      EXTERNAL CHEEVD , CHEGST , CPOTRF , CTRMM , CTRSM , XERBLA
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX , REAL
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      wantz = LSAME(Jobz,'V')
      upper = LSAME(Uplo,'U')
      lquery = (Lwork==-1 .OR. Lrwork==-1 .OR. Liwork==-1)
!
      Info = 0
      IF ( N<=1 ) THEN
         lwmin = 1
         lrwmin = 1
         liwmin = 1
      ELSEIF ( wantz ) THEN
         lwmin = 2*N + N*N
         lrwmin = 1 + 5*N + 2*N*N
         liwmin = 3 + 5*N
      ELSE
         lwmin = N + 1
         lrwmin = N
         liwmin = 1
      ENDIF
      lopt = lwmin
      lropt = lrwmin
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
         Rwork(1) = lropt
         Iwork(1) = liopt
!
         IF ( Lwork<lwmin .AND. .NOT.lquery ) THEN
            Info = -11
         ELSEIF ( Lrwork<lrwmin .AND. .NOT.lquery ) THEN
            Info = -13
         ELSEIF ( Liwork<liwmin .AND. .NOT.lquery ) THEN
            Info = -15
         ENDIF
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('CHEGVD',-Info)
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
      CALL CPOTRF(Uplo,N,B,Ldb,Info)
      IF ( Info/=0 ) THEN
         Info = N + Info
         RETURN
      ENDIF
!
!     Transform problem to standard eigenvalue problem and solve.
!
      CALL CHEGST(Itype,Uplo,N,A,Lda,B,Ldb,Info)
      CALL CHEEVD(Jobz,Uplo,N,A,Lda,W,Work,Lwork,Rwork,Lrwork,Iwork,    &
     &            Liwork,Info)
      lopt = MAX(REAL(lopt),REAL(Work(1)))
      lropt = MAX(REAL(lropt),REAL(Rwork(1)))
      liopt = MAX(REAL(liopt),REAL(Iwork(1)))
!
      IF ( wantz .AND. Info==0 ) THEN
!
!        Backtransform eigenvectors to the original problem.
!
         IF ( Itype==1 .OR. Itype==2 ) THEN
!
!           For A*x=(lambda)*B*x and A*B*x=(lambda)*x;
!           backtransform eigenvectors: x = inv(L)**H *y or inv(U)*y
!
            IF ( upper ) THEN
               trans = 'N'
            ELSE
               trans = 'C'
            ENDIF
!
            CALL CTRSM('Left',Uplo,trans,'Non-unit',N,N,CONE,B,Ldb,A,   &
     &                 Lda)
!
         ELSEIF ( Itype==3 ) THEN
!
!           For B*A*x=(lambda)*x;
!           backtransform eigenvectors: x = L*y or U**H *y
!
            IF ( upper ) THEN
               trans = 'C'
            ELSE
               trans = 'N'
            ENDIF
!
            CALL CTRMM('Left',Uplo,trans,'Non-unit',N,N,CONE,B,Ldb,A,   &
     &                 Lda)
         ENDIF
      ENDIF
!
      Work(1) = lopt
      Rwork(1) = lropt
      Iwork(1) = liopt
!
!
!     End of CHEGVD
!
      END SUBROUTINE CHEGVD