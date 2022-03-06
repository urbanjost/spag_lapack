!*==sqlt02.f90  processed by SPAG 7.51RB at 17:34 on  4 Mar 2022
!> \brief \b sqlt02
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE sqlt02( M, N, K, A, AF, Q, L, LDA, TAU, WORK, LWORK,
!                          RWORK, RESULT )
!
!       .. Scalar Arguments ..
!       INTEGER            K, LDA, LWORK, M, N
!       ..
!       .. Array Arguments ..
!       REAL               A( LDA, * ), AF( LDA, * ), L( LDA, * ),
!      $                   Q( LDA, * ), RESULT( * ), RWORK( * ), TAU( * ),
!      $                   WORK( LWORK )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> sqlt02 tests SORGQL, which generates an m-by-n matrix Q with
!> orthonornmal columns that is defined as the product of k elementary
!> reflectors.
!>
!> Given the QL factorization of an m-by-n matrix A, sqlt02 generates
!> the orthogonal matrix Q defined by the factorization of the last k
!> columns of A; it compares L(m-n+1:m,n-k+1:n) with
!> Q(1:m,m-n+1:m)'*A(1:m,n-k+1:n), and checks that the columns of Q are
!> orthonormal.
!> \endverbatim
!
!  Arguments:
!  ==========
!
!> \param[in] M
!> \verbatim
!>          M is INTEGER
!>          The number of rows of the matrix Q to be generated.  M >= 0.
!> \endverbatim
!>
!> \param[in] N
!> \verbatim
!>          N is INTEGER
!>          The number of columns of the matrix Q to be generated.
!>          M >= N >= 0.
!> \endverbatim
!>
!> \param[in] K
!> \verbatim
!>          K is INTEGER
!>          The number of elementary reflectors whose product defines the
!>          matrix Q. N >= K >= 0.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is REAL array, dimension (LDA,N)
!>          The m-by-n matrix A which was factorized by sqlt01.
!> \endverbatim
!>
!> \param[in] AF
!> \verbatim
!>          AF is REAL array, dimension (LDA,N)
!>          Details of the QL factorization of A, as returned by SGEQLF.
!>          See SGEQLF for further details.
!> \endverbatim
!>
!> \param[out] Q
!> \verbatim
!>          Q is REAL array, dimension (LDA,N)
!> \endverbatim
!>
!> \param[out] L
!> \verbatim
!>          L is REAL array, dimension (LDA,N)
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the arrays A, AF, Q and L. LDA >= M.
!> \endverbatim
!>
!> \param[in] TAU
!> \verbatim
!>          TAU is REAL array, dimension (N)
!>          The scalar factors of the elementary reflectors corresponding
!>          to the QL factorization in AF.
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
!>          The dimension of the array WORK.
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is REAL array, dimension (M)
!> \endverbatim
!>
!> \param[out] RESULT
!> \verbatim
!>          RESULT is REAL array, dimension (2)
!>          The test ratios:
!>          RESULT(1) = norm( L - Q'*A ) / ( M * norm(A) * EPS )
!>          RESULT(2) = norm( I - Q'*Q ) / ( M * EPS )
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
      SUBROUTINE sqlt02(M,N,K,A,Af,Q,L,Lda,Tau,Work,Lwork,Rwork,Result)
      IMPLICIT NONE
!*--SQLT02139
!
!  -- LAPACK test routine (version 3.7.0) --
!  -- LAPACK is a software package provided by Univ. of Tennessee,    --
!  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
!     December 2016
!
!     .. Scalar Arguments ..
      INTEGER K , Lda , Lwork , M , N
!     ..
!     .. Array Arguments ..
      REAL A(Lda,*) , Af(Lda,*) , L(Lda,*) , Q(Lda,*) , Result(*) ,     &
     &     Rwork(*) , Tau(*) , Work(Lwork)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      REAL ZERO , ONE
      PARAMETER (ZERO=0.0E+0,ONE=1.0E+0)
      REAL ROGUE
      PARAMETER (ROGUE=-1.0E+10)
!     ..
!     .. Local Scalars ..
      INTEGER info
      REAL anorm , eps , resid
!     ..
!     .. External Functions ..
      REAL SLAMCH , SLANGE , SLANSY
EXTERNAL    SLAMCH,   SLANGE,   SLANSY
!     ..
!     .. External Subroutines ..
EXTERNAL    SGEMM,   SLACPY,   SLASET,   SORGQL,   SSYRK
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC MAX , REAL
!     ..
!     .. Scalars in Common ..
      CHARACTER*32 SRNamt
!     ..
!     .. Common blocks ..
      COMMON /SRNAMC/ SRNamt
!     ..
!     .. Executable Statements ..
!
!     Quick return if possible
!
      IF ( M==0 .OR. N==0 .OR. K==0 ) THEN
         Result(1) = ZERO
         Result(2) = ZERO
         RETURN
      ENDIF
!
      eps = SLAMCH('Epsilon')
!
!     Copy the last k columns of the factorization to the array Q
!
      CALL SLASET('Full',M,N,ROGUE,ROGUE,Q,Lda)
      IF ( K<M ) CALL SLACPY('Full',M-K,K,Af(1,N-K+1),Lda,Q(1,N-K+1),   &
     &                       Lda)
      IF ( K>1 ) CALL SLACPY('Upper',K-1,K-1,Af(M-K+1,N-K+2),Lda,       &
     &                       Q(M-K+1,N-K+2),Lda)
!
!     Generate the last n columns of the matrix Q
!
      SRNamt = 'SORGQL'
      CALL SORGQL(M,N,K,Q,Lda,Tau(N-K+1),Work,Lwork,info)
!
!     Copy L(m-n+1:m,n-k+1:n)
!
      CALL SLASET('Full',N,K,ZERO,ZERO,L(M-N+1,N-K+1),Lda)
      CALL SLACPY('Lower',K,K,Af(M-K+1,N-K+1),Lda,L(M-K+1,N-K+1),Lda)
!
!     Compute L(m-n+1:m,n-k+1:n) - Q(1:m,m-n+1:m)' * A(1:m,n-k+1:n)
!
      CALL SGEMM('Transpose','No transpose',N,K,M,-ONE,Q,Lda,A(1,N-K+1),&
     &           Lda,ONE,L(M-N+1,N-K+1),Lda)
!
!     Compute norm( L - Q'*A ) / ( M * norm(A) * EPS ) .
!
      anorm = SLANGE('1',M,K,A(1,N-K+1),Lda,Rwork)
      resid = SLANGE('1',N,K,L(M-N+1,N-K+1),Lda,Rwork)
      IF ( anorm>ZERO ) THEN
         Result(1) = ((resid/REAL(MAX(1,M)))/anorm)/eps
      ELSE
         Result(1) = ZERO
      ENDIF
!
!     Compute I - Q'*Q
!
      CALL SLASET('Full',N,N,ZERO,ONE,L,Lda)
      CALL SSYRK('Upper','Transpose',N,M,-ONE,Q,Lda,ONE,L,Lda)
!
!     Compute norm( I - Q'*Q ) / ( M * EPS ) .
!
      resid = SLANSY('1','Upper',N,L,Lda,Rwork)
!
      Result(2) = (resid/REAL(MAX(1,M)))/eps
!
!
!     End of sqlt02
!
      END SUBROUTINE sqlt02
