!*==dgelqs.f90  processed by SPAG 7.51RB at 17:34 on  4 Mar 2022
!> \brief \b dgelqs
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE dgelqs( M, N, NRHS, A, LDA, TAU, B, LDB, WORK, LWORK,
!                          INFO )
!
!       .. Scalar Arguments ..
!       INTEGER            INFO, LDA, LDB, LWORK, M, N, NRHS
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), TAU( * ),
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
!> using the LQ factorization
!>     A = L*Q
!> computed by DGELQF.
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
!>          A is DOUBLE PRECISION array, dimension (LDA,N)
!>          Details of the LQ factorization of the original matrix A as
!>          returned by DGELQF.
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
!>          TAU is DOUBLE PRECISION array, dimension (M)
!>          Details of the orthogonal matrix Q.
!> \endverbatim
!>
!> \param[in,out] B
!> \verbatim
!>          B is DOUBLE PRECISION array, dimension (LDB,NRHS)
!>          On entry, the m-by-nrhs right hand side matrix B.
!>          On exit, the n-by-nrhs solution matrix X.
!> \endverbatim
!>
!> \param[in] LDB
!> \verbatim
!>          LDB is INTEGER
!>          The leading dimension of the array B. LDB >= N.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is DOUBLE PRECISION array, dimension (LWORK)
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
!> \ingroup double_lin
!
!  =====================================================================
      SUBROUTINE dgelqs(M,N,Nrhs,A,Lda,Tau,B,Ldb,Work,Lwork,Info)
      IMPLICIT NONE
!*--DGELQS124
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
      DOUBLE PRECISION A(Lda,*) , B(Ldb,*) , Tau(*) , Work(Lwork)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO , ONE
      PARAMETER (ZERO=0.0D+0,ONE=1.0D+0)
!     ..
!     .. External Subroutines ..
EXTERNAL    DLASET,   DORMLQ,   DTRSM
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
         CALL xerbla('dgelqs',-Info)
         RETURN
      ENDIF
!
!     Quick return if possible
!
      IF ( N==0 .OR. Nrhs==0 .OR. M==0 ) RETURN
!
!     Solve L*X = B(1:m,:)
!
      CALL DTRSM('Left','Lower','No transpose','Non-unit',M,Nrhs,ONE,A, &
     &           Lda,B,Ldb)
!
!     Set B(m+1:n,:) to zero
!
      IF ( M<N ) CALL DLASET('Full',N-M,Nrhs,ZERO,ZERO,B(M+1,1),Ldb)
!
!     B := Q' * B
!
      CALL DORMLQ('Left','Transpose',N,Nrhs,M,A,Lda,Tau,B,Ldb,Work,     &
     &            Lwork,Info)
!
!
!     End of dgelqs
!
      END SUBROUTINE dgelqs
