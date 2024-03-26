!*==zlqt02.f90  processed by SPAG 7.51RB at 17:35 on  4 Mar 2022
!> \brief \b zlqt02
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE zlqt02( M, N, K, A, AF, Q, L, LDA, TAU, WORK, LWORK,
!                          RWORK, RESULT )
!
!       .. Scalar Arguments ..
!       INTEGER            K, LDA, LWORK, M, N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   RESULT( * ), RWORK( * )
!       COMPLEX*16         A( LDA, * ), AF( LDA, * ), L( LDA, * ),
!      $                   Q( LDA, * ), TAU( * ), WORK( LWORK )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> zlqt02 tests ZUNGLQ, which generates an m-by-n matrix Q with
!> orthonornmal rows that is defined as the product of k elementary
!> reflectors.
!>
!> Given the LQ factorization of an m-by-n matrix A, zlqt02 generates
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
!>          A is COMPLEX*16 array, dimension (LDA,N)
!>          The m-by-n matrix A which was factorized by zlqt01.
!> \endverbatim
!>
!> \param[in] AF
!> \verbatim
!>          AF is COMPLEX*16 array, dimension (LDA,N)
!>          Details of the LQ factorization of A, as returned by ZGELQF.
!>          See ZGELQF for further details.
!> \endverbatim
!>
!> \param[out] Q
!> \verbatim
!>          Q is COMPLEX*16 array, dimension (LDA,N)
!> \endverbatim
!>
!> \param[out] L
!> \verbatim
!>          L is COMPLEX*16 array, dimension (LDA,M)
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
!>          TAU is COMPLEX*16 array, dimension (M)
!>          The scalar factors of the elementary reflectors corresponding
!>          to the LQ factorization in AF.
!> \endverbatim
!>
!> \param[out] WORK
!> \verbatim
!>          WORK is COMPLEX*16 array, dimension (LWORK)
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
!> \ingroup complex16_lin
!
!  =====================================================================
      SUBROUTINE zlqt02(M,N,K,A,Af,Q,L,Lda,Tau,Work,Lwork,Rwork,Result)
      IMPLICIT NONE
!*--ZLQT02138
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
      DOUBLE PRECISION Result(*) , Rwork(*)
      COMPLEX*16 A(Lda,*) , Af(Lda,*) , L(Lda,*) , Q(Lda,*) , Tau(*) ,  &
     &           Work(Lwork)
!     ..
!
!  =====================================================================
!
!     .. Parameters ..
      DOUBLE PRECISION ZERO , ONE
      PARAMETER (ZERO=0.0D+0,ONE=1.0D+0)
      COMPLEX*16 ROGUE
      PARAMETER (ROGUE=(-1.0D+10,-1.0D+10))
!     ..
!     .. Local Scalars ..
      INTEGER info
      DOUBLE PRECISION anorm , eps , resid
!     ..
!     .. External Functions ..
      DOUBLE PRECISION DLAMCH , ZLANGE , ZLANSY
EXTERNAL    DLAMCH,   ZLANGE,   ZLANSY
!     ..
!     .. External Subroutines ..
EXTERNAL    ZGEMM,   ZHERK,   ZLACPY,   ZLASET,   ZUNGLQ
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC DBLE , DCMPLX , MAX
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
      CALL ZLASET('Full',M,N,ROGUE,ROGUE,Q,Lda)
      CALL ZLACPY('Upper',K,N-1,Af(1,2),Lda,Q(1,2),Lda)
!
!     Generate the first n columns of the matrix Q
!
      SRNamt = 'ZUNGLQ'
      CALL ZUNGLQ(M,N,K,Q,Lda,Tau,Work,Lwork,info)
!
!     Copy L(1:k,1:m)
!
      CALL ZLASET('Full',K,M,DCMPLX(ZERO),DCMPLX(ZERO),L,Lda)
      CALL ZLACPY('Lower',K,M,Af,Lda,L,Lda)
!
!     Compute L(1:k,1:m) - A(1:k,1:n) * Q(1:m,1:n)'
!
      CALL ZGEMM('No transpose','Conjugate transpose',K,M,N,DCMPLX(-ONE)&
     &           ,A,Lda,Q,Lda,DCMPLX(ONE),L,Lda)
!
!     Compute norm( L - A*Q' ) / ( N * norm(A) * EPS ) .
!
      anorm = ZLANGE('1',K,N,A,Lda,Rwork)
      resid = ZLANGE('1',K,M,L,Lda,Rwork)
      IF ( anorm>ZERO ) THEN
         Result(1) = ((resid/DBLE(MAX(1,N)))/anorm)/eps
      ELSE
         Result(1) = ZERO
      ENDIF
!
!     Compute I - Q*Q'
!
      CALL ZLASET('Full',M,M,DCMPLX(ZERO),DCMPLX(ONE),L,Lda)
      CALL ZHERK('Upper','No transpose',M,N,-ONE,Q,Lda,ONE,L,Lda)
!
!     Compute norm( I - Q*Q' ) / ( N * EPS ) .
!
      resid = ZLANSY('1','Upper',M,L,Lda,Rwork)
!
      Result(2) = (resid/DBLE(MAX(1,N)))/eps
!
!
!     End of zlqt02
!
      END SUBROUTINE zlqt02
