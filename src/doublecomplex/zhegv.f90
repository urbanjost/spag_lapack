!*==zhegv.f90  processed by SPAG 7.51RB at 20:08 on  3 Mar 2022
!> \brief \b ZHEGV
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!> \htmlonly
!> Download ZHEGV + dependencies
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zhegv.f">
!> [TGZ]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zhegv.f">
!> [ZIP]</a>
!> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zhegv.f">
!> [TXT]</a>
!> \endhtmlonly
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZHEGV( ITYPE, JOBZ, UPLO, N, A, LDA, B, LDB, W, WORK,
!                         LWORK, RWORK, INFO )
!
!       .. Scalar Arguments ..
!       CHARACTER          JOBZ, UPLO
!       INTEGER            INFO, ITYPE, LDA, LDB, LWORK, N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   RWORK( * ), W( * )
!       COMPLEX*16         A( LDA, * ), B( LDB, * ), WORK( * )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZHEGV computes all the eigenvalues, and optionally, the eigenvectors
!> of a complex generalized Hermitian-definite eigenproblem, of the form
!> A*x=(lambda)*B*x,  A*Bx=(lambda)*x,  or B*A*x=(lambda)*x.
!> Here A and B are assumed to be Hermitian and B is also
!> positive definite.
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
!>          A is COMPLEX*16 array, dimension (LDA, N)
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
!>          B is COMPLEX*16 array, dimension (LDB, N)
!>          On entry, the Hermitian positive definite matrix B.
!>          If UPLO = 'U', the leading N-by-N upper triangular part of B
!>          contains the upper triangular part of the matrix B.
!>          If UPLO = 'L', the leading N-by-N lower triangular part of B
!>          contains the lower triangular part of the matrix B.
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
!>          W is DOUBLE PRECISION array, dimension (N)
!>          If INFO = 0, the eigenvalues in ascending order.
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
!>          The length of the array WORK.  LWORK >= max(1,2*N-1).
!>          For optimal efficiency, LWORK >= (NB+1)*N,
!>          where NB is the blocksize for ZHETRD returned by ILAENV.
!>
!>          If LWORK = -1, then a workspace query is assumed; the routine
!>          only calculates the optimal size of the WORK array, returns
!>          this value as the first entry of the WORK array, and no error
!>          message related to LWORK is issued by XERBLA.
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is DOUBLE PRECISION array, dimension (max(1, 3*N-2))
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0:  successful exit
!>          < 0:  if INFO = -i, the i-th argument had an illegal value
!>          > 0:  ZPOTRF or ZHEEV returned an error code:
!>             <= N:  if INFO = i, ZHEEV failed to converge;
!>                    i off-diagonal elements of an intermediate
!>                    tridiagonal form did not converge to zero;
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
!> \ingroup complex16HEeigen
!
!  =====================================================================
      SUBROUTINE ZHEGV(Itype,Jobz,Uplo,N,A,Lda,B,Ldb,W,Work,Lwork,Rwork,&
     &                 Info)
      IMPLICIT NONE
!*--ZHEGV185
!
!  -- LAPACK driver routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      CHARACTER Jobz , Uplo
      INTEGER Info , Itype , Lda , Ldb , Lwork , N
!     ..
!     .. Array Arguments ..
      DOUBLE PRECISION Rwork(*) , W(*)
      COMPLEX*16 A(Lda,*) , B(Ldb,*) , Work(*)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      COMPLEX*16 ONE
      PARAMETER (ONE=(1.0D+0,0.0D+0))
!     ..
!     .. Local Scalars ..
      LOGICAL lquery , upper , wantz
      CHARACTER trans
      INTEGER lwkopt , nb , neig
!     ..
!     .. External Functions ..
      LOGICAL LSAME
      INTEGER ILAENV
      EXTERNAL LSAME , ILAENV
!     ..
!     .. External Subroutines ..
      EXTERNAL XERBLA , ZHEEV , ZHEGST , ZPOTRF , ZTRMM , ZTRSM
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      wantz = LSAME(Jobz,'V')
      upper = LSAME(Uplo,'U')
      lquery = (Lwork==-1)
!
      Info = 0
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
         nb = ILAENV(1,'ZHETRD',Uplo,N,-1,-1,-1)
         lwkopt = MAX(1,(nb+1)*N)
         Work(1) = lwkopt
!
         IF ( Lwork<MAX(1,2*N-1) .AND. .NOT.lquery ) Info = -11
      ENDIF
!
      IF ( Info/=0 ) THEN
         CALL XERBLA('ZHEGV ',-Info)
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
      CALL ZPOTRF(Uplo,N,B,Ldb,Info)
      IF ( Info/=0 ) THEN
         Info = N + Info
         RETURN
      ENDIF
!
!     Transform problem to standard eigenvalue problem and solve.
!
      CALL ZHEGST(Itype,Uplo,N,A,Lda,B,Ldb,Info)
      CALL ZHEEV(Jobz,Uplo,N,A,Lda,W,Work,Lwork,Rwork,Info)
!
      IF ( wantz ) THEN
!
!        Backtransform eigenvectors to the original problem.
!
         neig = N
         IF ( Info>0 ) neig = Info - 1
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
            CALL ZTRSM('Left',Uplo,trans,'Non-unit',N,neig,ONE,B,Ldb,A, &
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
            CALL ZTRMM('Left',Uplo,trans,'Non-unit',N,neig,ONE,B,Ldb,A, &
     &                 Lda)
         ENDIF
      ENDIF
!
      Work(1) = lwkopt
!
!
!     End of ZHEGV
!
      END SUBROUTINE ZHEGV
