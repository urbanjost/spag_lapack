!*==sgerqs.f90  processed by SPAG 7.51RB at 17:34 on  4 Mar 2022
!> \brief \b sgerqs
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE sgerqs( M, N, NRHS, A, LDA, TAU, B, LDB, WORK, LWORK,
!                          INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, LDA, LDB, LWORK, M, N, NRHS
!       ..
!       .. Array Arguments ..
!       REAL               A( LDA, * ), B( LDB, * ), TAU( * ),
!      $                   WORK( LWORK )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> Compute a minimum-norm solution
!>     min || A*X - B ||
!> using the RQ factorization
!>     A = R*Q
!> computed by SGERQF.
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
!>          The number of columns of the matrix A.  N >= M >= 0.
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
!>          A is REAL array, dimension (LDA,N)
!>          Details of the RQ factorization of the original matrix A as
!>          returned by SGERQF.
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
!>          TAU is REAL array, dimension (M)
!>          Details of the orthogonal matrix Q.
!> \endverbatim
!>
!> \param[in,out] B
!> \verbatim
!>          B is REAL array, dimension (LDB,NRHS)
!>          On entry, the right hand side vectors for the linear system.
!>          On exit, the solution vectors X.  Each solution vector
!>          is contained in rows 1:N of a column of B.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of the array B. LDB >= max(1,N).
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is REAL array, dimension (LWORK)
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
!> \ingroup single_lin
!
!  =====================================================================
      SUBROUTINE sgerqs(M,N,Nrhs,A,Lda,Tau,B,Ldb,Work,Lwork,Info)
      IMPLICIT NONE
!*--SGERQS125
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
      REAL A(Lda,*) , B(Ldb,*) , Tau(*) , Work(Lwork)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ZERO , ONE
      PARAMETER (ZERO=0.0E+0,ONE=1.0E+0)
!     ..
!     .. External Subroutines ..
EXTERNAL    SLASET,   SORMRQ,   STRSM
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX
!     ..
!     .. Executable Statements ..
!
!     Test the input parameters.
!
      Info = 0
      IF ( M<0 ) THEN
         Info = -1
      ELSEIF ( N<0 .OR. M>N ) THEN
         Info = -2
      ELSEIF ( Nrhs<0 ) THEN
         Info = -3
      ELSEIF ( Lda<MAX(1,M) ) THEN
         Info = -5
      ELSEIF ( Ldb<MAX(1,N) ) THEN
         Info = -8
      ELSEIF ( Lwork<1 .OR. Lwork<Nrhs .AND. M>0 .AND. N>0 ) THEN
         Info = -10
      ENDIF
      IF ( Info/=0 ) THEN
         CALL xerbla('sgerqs',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( N==0 .OR. Nrhs==0 .OR. M==0 ) RETURN
!
!     Solve R*X = B(n-m+1:n,:)
!
      CALL STRSM('Left','Upper','No transpose','Non-unit',M,Nrhs,ONE,   &
     &           A(1,N-M+1),Lda,B(N-M+1,1),Ldb)
!
!     Set B(1:n-m,:) to zero
!
      CALL SLASET('Full',N-M,Nrhs,ZERO,ZERO,B,Ldb)
!
!     B := Q' * B
!
      CALL SORMRQ('Left','Transpose',N,Nrhs,M,A,Lda,Tau,B,Ldb,Work,     &
     &            Lwork,Info)
!
!
!     End of sgerqs
!
      END SUBROUTINE sgerqs
