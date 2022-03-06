!*==slqt02.f90  processed by SPAG 7.51RB at 17:34 on  4 Mar 2022
!> \brief \b slqt02
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE slqt02( M, N, K, A, AF, Q, L, LDA, TAU, WORK, LWORK,
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
!> slqt02 tests SORGLQ, which generates an m-by-n matrix Q with
!> orthonornmal rows that is defined as the product of k elementary
!> reflectors.
!>
!> Given the LQ factorization of an m-by-n matrix A, slqt02 generates
!> the orthogonal matrix Q defined by the factorization of the first k
!> rows of A; it compares L(1:k,1:m) with A(1:k,1:n)*Q(1:m,1:n)', and
!> checks that the rows of Q are orthonormal.
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
!>          N >= M >= 0.
!> \endverbatim
!>
!> \param[in] K
!> \verbatim
!>          K is INTEGER
!>          The number of elementary reflectors whose product defines the
!>          matrix Q. M >= K >= 0.
!> \endverbatim
!>
!> \param[in] A
!> \verbatim
!>          A is REAL array, dimension (LDA,N)
!>          The m-by-n matrix A which was factorized by slqt01.
!> \endverbatim
!>
!> \param[in] AF
!> \verbatim
!>          AF is REAL array, dimension (LDA,N)
!>          Details of the LQ factorization of A, as returned by SGELQF.
!>          See SGELQF for further details.
!> \endverbatim
!>
!> \param[out] Q
!> \verbatim
!>          Q is REAL array, dimension (LDA,N)
!> \endverbatim
!>
!> \param[out] L
!> \verbatim
!>          L is REAL array, dimension (LDA,M)
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the arrays A, AF, Q and L. LDA >= N.
!> \endverbatim
!>
!> \param[in] TAU
!> \verbatim
!>          TAU is REAL array, dimension (M)
!>          The scalar factors of the elementary reflectors corresponding
!>          to the LQ factorization in AF.
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
!>          RESULT(1) = norm( L - A*Q' ) / ( N * norm(A) * EPS )
!>          RESULT(2) = norm( I - Q*Q' ) / ( N * EPS )
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
      SUBROUTINE slqt02(M,N,K,A,Af,Q,L,Lda,Tau,Work,Lwork,Rwork,Result)
      IMPLICIT NONE
!*--SLQT02138
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
EXTERNAL    SGEMM,   SLACPY,   SLASET,   SORGLQ,   SSYRK
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
      eps = SLAMCH('Epsilon')
!
!     Copy the first k rows of the factorization to the array Q
!
      CALL SLASET('Full',M,N,ROGUE,ROGUE,Q,Lda)
      CALL SLACPY('Upper',K,N-1,Af(1,2),Lda,Q(1,2),Lda)
!
!     Generate the first n columns of the matrix Q
!
      SRNamt = 'SORGLQ'
      CALL SORGLQ(M,N,K,Q,Lda,Tau,Work,Lwork,info)
!
!     Copy L(1:k,1:m)
!
      CALL SLASET('Full',K,M,ZERO,ZERO,L,Lda)
      CALL SLACPY('Lower',K,M,Af,Lda,L,Lda)
!
!     Compute L(1:k,1:m) - A(1:k,1:n) * Q(1:m,1:n)'
!
      CALL SGEMM('No transpose','Transpose',K,M,N,-ONE,A,Lda,Q,Lda,ONE, &
     &           L,Lda)
!
!     Compute norm( L - A*Q' ) / ( N * norm(A) * EPS ) .
!
      anorm = SLANGE('1',K,N,A,Lda,Rwork)
      resid = SLANGE('1',K,M,L,Lda,Rwork)
      IF ( anorm>ZERO ) THEN
         Result(1) = ((resid/REAL(MAX(1,N)))/anorm)/eps
      ELSE
         Result(1) = ZERO
      ENDIF
!
!     Compute I - Q*Q'
!
      CALL SLASET('Full',M,M,ZERO,ONE,L,Lda)
      CALL SSYRK('Upper','No transpose',M,N,-ONE,Q,Lda,ONE,L,Lda)
!
!     Compute norm( I - Q*Q' ) / ( N * EPS ) .
!
      resid = SLANSY('1','Upper',M,L,Lda,Rwork)
!
      Result(2) = (resid/REAL(MAX(1,N)))/eps
!
!
!     End of slqt02
!
      END SUBROUTINE slqt02
