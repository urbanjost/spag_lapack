!*==dlqt02.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b DLQT02
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DLQT02( M, N, K, A, AF, Q, L, LDA, TAU, WORK, LWORK,
!                          RWORK, RESULT )
!
!       .. Scalar Arguments ..
!       INTEGER            K, LDA, LWORK, M, N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   A( LDA, * ), AF( LDA, * ), L( LDA, * ),
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
!> DLQT02 tests DORGLQ, which generates an m-by-n matrix Q with
!> orthonornmal rows that is defined as the product of k elementary
!> reflectors.
!>
!> Given the LQ factorization of an m-by-n matrix A, DLQT02 generates
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
!>          A is DOUBLE PRECISION array, dimension (LDA,N)
!>          The m-by-n matrix A which was factorized by DLQT01.
!> \endverbatim
!>
!> \param[in] AF
!> \verbatim
!>          AF is DOUBLE PRECISION array, dimension (LDA,N)
!>          Details of the LQ factorization of A, as returned by DGELQF.
!>          See DGELQF for further details.
!> \endverbatim
!>
!> \param[out] Q
!> \verbatim
!>          Q is DOUBLE PRECISION array, dimension (LDA,N)
!> \endverbatim
!>
!> \param[out] L
!> \verbatim
!>          L is DOUBLE PRECISION array, dimension (LDA,M)
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
!>          TAU is DOUBLE PRECISION array, dimension (M)
!>          The scalar factors of the elementary reflectors corresponding
!>          to the LQ factorization in AF.
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
!>          The dimension of the array WORK.
!> \endverbatim
!>
!> \param[out] RWORK
!> \verbatim
!>          RWORK is DOUBLE PRECISION array, dimension (M)
!> \endverbatim
!>
!> \param[out] RESULT
!> \verbatim
!>          RESULT is DOUBLE PRECISION array, dimension (2)
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
!> \ingroup double_lin
!
!  =====================================================================
      SUBROUTINE DLQT02(M,N,K,A,Af,Q,L,Lda,Tau,Work,Lwork,Rwork,Result)
      IMPLICIT NONE
!*--DLQT02138
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
      DOUBLE PRECISION A(Lda,*) , Af(Lda,*) , L(Lda,*) , Q(Lda,*) ,     &
     &                 Result(*) , Rwork(*) , Tau(*) , Work(Lwork)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO , ONE
      PARAMETER (ZERO=0.0D+0,ONE=1.0D+0)
      DOUBLE PRECISION ROGUE
      PARAMETER (ROGUE=-1.0D+10)
!     ..
!     .. Local Scalars ..
      INTEGER info
      DOUBLE PRECISION anorm , eps , resid
!     ..
!     .. External Functions ..
      DOUBLE PRECISION DLAMCH , DLANGE , DLANSY
      EXTERNAL DLAMCH , DLANGE , DLANSY
!     ..
!     .. External Subroutines ..
      EXTERNAL DGEMM , DLACPY , DLASET , DORGLQ , DSYRK
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC DBLE , MAX
!     ..
!     .. Scalars in Common ..
      CHARACTER*32 SRNamt
!     ..
!     .. Common blocks ..
      COMMON /SRNAMC/ SRNamt
!     ..
!     .. Executable Statements ..
!
      eps = DLAMCH('Epsilon')
!
!     Copy the first k rows of the factorization to the array Q
!
      CALL DLASET('Full',M,N,ROGUE,ROGUE,Q,Lda)
      CALL DLACPY('Upper',K,N-1,Af(1,2),Lda,Q(1,2),Lda)
!
!     Generate the first n columns of the matrix Q
!
      SRNamt = 'DORGLQ'
      CALL DORGLQ(M,N,K,Q,Lda,Tau,Work,Lwork,info)
!
!     Copy L(1:k,1:m)
!
      CALL DLASET('Full',K,M,ZERO,ZERO,L,Lda)
      CALL DLACPY('Lower',K,M,Af,Lda,L,Lda)
!
!     Compute L(1:k,1:m) - A(1:k,1:n) * Q(1:m,1:n)'
!
      CALL DGEMM('No transpose','Transpose',K,M,N,-ONE,A,Lda,Q,Lda,ONE, &
     &           L,Lda)
!
!     Compute norm( L - A*Q' ) / ( N * norm(A) * EPS ) .
!
      anorm = DLANGE('1',K,N,A,Lda,Rwork)
      resid = DLANGE('1',K,M,L,Lda,Rwork)
      IF ( anorm>ZERO ) THEN
         Result(1) = ((resid/DBLE(MAX(1,N)))/anorm)/eps
      ELSE
         Result(1) = ZERO
      ENDIF
!
!     Compute I - Q*Q'
!
      CALL DLASET('Full',M,M,ZERO,ONE,L,Lda)
      CALL DSYRK('Upper','No transpose',M,N,-ONE,Q,Lda,ONE,L,Lda)
!
!     Compute norm( I - Q*Q' ) / ( N * EPS ) .
!
      resid = DLANSY('1','Upper',M,L,Lda,Rwork)
!
      Result(2) = (resid/DBLE(MAX(1,N)))/eps
!
!
!     End of DLQT02
!
      END SUBROUTINE DLQT02
