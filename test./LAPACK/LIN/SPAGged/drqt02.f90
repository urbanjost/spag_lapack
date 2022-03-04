!*==drqt02.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b DRQT02
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DRQT02( M, N, K, A, AF, Q, R, LDA, TAU, WORK, LWORK,
!                          RWORK, RESULT )
!
!       .. Scalar Arguments ..
!       INTEGER            K, LDA, LWORK, M, N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   A( LDA, * ), AF( LDA, * ), Q( LDA, * ),
!      $                   R( LDA, * ), RESULT( * ), RWORK( * ), TAU( * ),
!      $                   WORK( LWORK )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> DRQT02 tests DORGRQ, which generates an m-by-n matrix Q with
!> orthonornmal rows that is defined as the product of k elementary
!> reflectors.
!>
!> Given the RQ factorization of an m-by-n matrix A, DRQT02 generates
!> the orthogonal matrix Q defined by the factorization of the last k
!> rows of A; it compares R(m-k+1:m,n-m+1:n) with
!> A(m-k+1:m,1:n)*Q(n-m+1:n,1:n)', and checks that the rows of Q are
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
!>          The m-by-n matrix A which was factorized by DRQT01.
!> \endverbatim
!>
!> \param[in] AF
!> \verbatim
!>          AF is DOUBLE PRECISION array, dimension (LDA,N)
!>          Details of the RQ factorization of A, as returned by DGERQF.
!>          See DGERQF for further details.
!> \endverbatim
!>
!> \param[out] Q
!> \verbatim
!>          Q is DOUBLE PRECISION array, dimension (LDA,N)
!> \endverbatim
!>
!> \param[out] R
!> \verbatim
!>          R is DOUBLE PRECISION array, dimension (LDA,M)
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
!>          to the RQ factorization in AF.
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
!>          RESULT(1) = norm( R - A*Q' ) / ( N * norm(A) * EPS )
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
      SUBROUTINE DRQT02(M,N,K,A,Af,Q,R,Lda,Tau,Work,Lwork,Rwork,Result)
      IMPLICIT NONE
!*--DRQT02139
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
      DOUBLE PRECISION A(Lda,*) , Af(Lda,*) , Q(Lda,*) , R(Lda,*) ,     &
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
      EXTERNAL DGEMM , DLACPY , DLASET , DORGRQ , DSYRK
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
!     Quick return if possible
!
      IF ( M==0 .OR. N==0 .OR. K==0 ) THEN
         Result(1) = ZERO
         Result(2) = ZERO
         RETURN
      ENDIF
!
      eps = DLAMCH('Epsilon')
!
!     Copy the last k rows of the factorization to the array Q
!
      CALL DLASET('Full',M,N,ROGUE,ROGUE,Q,Lda)
      IF ( K<N ) CALL DLACPY('Full',K,N-K,Af(M-K+1,1),Lda,Q(M-K+1,1),   &
     &                       Lda)
      IF ( K>1 ) CALL DLACPY('Lower',K-1,K-1,Af(M-K+2,N-K+1),Lda,       &
     &                       Q(M-K+2,N-K+1),Lda)
!
!     Generate the last n rows of the matrix Q
!
      SRNamt = 'DORGRQ'
      CALL DORGRQ(M,N,K,Q,Lda,Tau(M-K+1),Work,Lwork,info)
!
!     Copy R(m-k+1:m,n-m+1:n)
!
      CALL DLASET('Full',K,M,ZERO,ZERO,R(M-K+1,N-M+1),Lda)
      CALL DLACPY('Upper',K,K,Af(M-K+1,N-K+1),Lda,R(M-K+1,N-K+1),Lda)
!
!     Compute R(m-k+1:m,n-m+1:n) - A(m-k+1:m,1:n) * Q(n-m+1:n,1:n)'
!
      CALL DGEMM('No transpose','Transpose',K,M,N,-ONE,A(M-K+1,1),Lda,Q,&
     &           Lda,ONE,R(M-K+1,N-M+1),Lda)
!
!     Compute norm( R - A*Q' ) / ( N * norm(A) * EPS ) .
!
      anorm = DLANGE('1',K,N,A(M-K+1,1),Lda,Rwork)
      resid = DLANGE('1',K,M,R(M-K+1,N-M+1),Lda,Rwork)
      IF ( anorm>ZERO ) THEN
         Result(1) = ((resid/DBLE(MAX(1,N)))/anorm)/eps
      ELSE
         Result(1) = ZERO
      ENDIF
!
!     Compute I - Q*Q'
!
      CALL DLASET('Full',M,M,ZERO,ONE,R,Lda)
      CALL DSYRK('Upper','No transpose',M,N,-ONE,Q,Lda,ONE,R,Lda)
!
!     Compute norm( I - Q*Q' ) / ( N * EPS ) .
!
      resid = DLANSY('1','Upper',M,R,Lda,Rwork)
!
      Result(2) = (resid/DBLE(MAX(1,N)))/eps
!
!
!     End of DRQT02
!
      END SUBROUTINE DRQT02
