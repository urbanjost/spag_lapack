!*==cgeqls.f90  processed by SPAG 7.51RB at 17:34 on  4 Mar 2022
!> \brief \b cgeqls
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE cgeqls( M, N, NRHS, A, LDA, TAU, B, LDB, WORK, LWORK,
!                          INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, LDA, LDB, LWORK, M, N, NRHS
!       ..
!       .. Array Arguments ..
!       COMPLEX            A( LDA, * ), B( LDB, * ), TAU( * ),
!      $                   WORK( LWORK )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> Solve the least squares problem
!>     min || A*X - B ||
!> using the QL factorization
!>     A = Q*L
!> computed by CGEQLF.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix A.  M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix A.  M >= N >= 0.
!> \endverbatim
!>
!> \param[in] NRHS
!> \verbatim
!>          NRHS is INTEGER
!>          The number of columns of B.  NRHS >= 0.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is COMPLEX array, dimension (LDA,N)
!>          Details of the QL factorization of the original matrix A as
!>          returned by CGEQLF.
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the array A.  LDA >= M.
!> \endverbatim
!>
!> \param[in] TAU
!> \verbatim
!>          TAU is COMPLEX array, dimension (N)
!>          Details of the orthogonal matrix Q.
!> \endverbatim
!>
!> \param[in,out] B
!> \verbatim
!>          B is COMPLEX array, dimension (LDB,NRHS)
!>          On entry, the m-by-nrhs right hand side matrix B.
!>          On exit, the n-by-nrhs solution matrix X, stored in rows
!>          m-n+1:m.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of the array B. LDB >= M.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX array, dimension (LWORK)
!> \endverbatim
!>
!> \param[in] LWORK
!> \verbatim
!>          LWORK is INTEGER
!>          The length of the array WORK.  LWORK must be at least NRHS,
!>          and should be at least NRHS*NB, where NB is the block size
!>          for this environment.
!> \endverbatim
!>
!> \param[out] INFO
!> \verbatim
!>          INFO is INTEGER
!>          = 0: successful exit
!>          < 0: if INFO = -i, the i-th argument had an illegal value
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
!> \ingroup complex_lin
!
!  =====================================================================
      SUBROUTINE cgeqls(M,N,Nrhs,A,Lda,Tau,B,Ldb,Work,Lwork,Info)
      IMPLICIT NONE
!*--CGEQLS125
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER Info , Lda , Ldb , Lwork , M , N , Nrhs
!     ..
!     .. Array Arguments ..
      COMPLEX A(Lda,*) , B(Ldb,*) , Tau(*) , Work(Lwork)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      COMPLEX ONE
      PARAMETER (ONE=(1.0E+0,0.0E+0))
!     ..
!     .. External Subroutines ..
EXTERNAL    CTRSM,   CUNMQL
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input arguments.
!
      Info = 0
      IF ( M<0 ) THEN
         Info = -1
      ELSEIF ( N<0 .OR. N>M ) THEN
         Info = -2
      ELSEIF ( Nrhs<0 ) THEN
         Info = -3
      ELSEIF ( Lda<MAX(1,M) ) THEN
         Info = -5
      ELSEIF ( Ldb<MAX(1,M) ) THEN
         Info = -8
      ELSEIF ( Lwork<1 .OR. Lwork<Nrhs .AND. M>0 .AND. N>0 ) THEN
         Info = -10
      ENDIF
      IF ( Info/=0 ) THEN
         CALL xerbla('cgeqls',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( N==0 .OR. Nrhs==0 .OR. M==0 ) RETURN
!
!     B := Q' * B
!
      CALL CUNMQL('Left','Conjugate transpose',M,Nrhs,N,A,Lda,Tau,B,Ldb,&
     &            Work,Lwork,Info)
!
!     Solve L*X = B(m-n+1:m,:)
!
      CALL CTRSM('Left','Lower','No transpose','Non-unit',N,Nrhs,ONE,   &
     &           A(M-N+1,1),Lda,B(M-N+1,1),Ldb)
!
!
!     End of cgeqls
!
      END SUBROUTINE cgeqls
