!*==dqrt02.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b DQRT02
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE DQRT02( M, N, K, A, AF, Q, R, LDA, TAU, WORK, LWORK,
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
!> DQRT02 tests DORGQR, which generates an m-by-n matrix Q with
!> orthonornmal columns that is defined as the product of k elementary
!> reflectors.
!>
!> Given the QR factorization of an m-by-n matrix A, DQRT02 generates
!> the orthogonal matrix Q defined by the factorization of the first k
!> columns of A; it compares R(1:n,1:k) with Q(1:m,1:n)'*A(1:m,1:k),
!> and checks that the columns of Q are orthonormal.
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
!>          A is DOUBLE PRECISION array, dimension (LDA,N)
!>          The m-by-n matrix A which was factorized by DQRT01.
!> \endverbatim
!>
!> \param[in] AF
!> \verbatim
!>          AF is DOUBLE PRECISION array, dimension (LDA,N)
!>          Details of the QR factorization of A, as returned by DGEQRF.
!>          See DGEQRF for further details.
!> \endverbatim
!>
!> \param[out] Q
!> \verbatim
!>          Q is DOUBLE PRECISION array, dimension (LDA,N)
!> \endverbatim
!>
!> \param[out] R
!> \verbatim
!>          R is DOUBLE PRECISION array, dimension (LDA,N)
!> \endverbatim
!>
!> \param[in] LDA
!> \verbatim
!>          LDA is INTEGER
!>          The leading dimension of the arrays A, AF, Q and R. LDA >= M.
!> \endverbatim
!>
!> \param[in] TAU
!> \verbatim
!>          TAU is DOUBLE PRECISION array, dimension (N)
!>          The scalar factors of the elementary reflectors corresponding
!>          to the QR factorization in AF.
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
!>          RESULT(1) = norm( R - Q'*A ) / ( M * norm(A) * EPS )
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
!> \ingroup double_lin
!
!  =====================================================================
      SUBROUTINE DQRT02(M,N,K,A,Af,Q,R,Lda,Tau,Work,Lwork,Rwork,Result)
      IMPLICIT NONE
!*--DQRT02138
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
      EXTERNAL DGEMM , DLACPY , DLASET , DORGQR , DSYRK
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
!     Copy the first k columns of the factorization to the array Q
!
      CALL DLASET('Full',M,N,ROGUE,ROGUE,Q,Lda)
      CALL DLACPY('Lower',M-1,K,Af(2,1),Lda,Q(2,1),Lda)
!
!     Generate the first n columns of the matrix Q
!
      SRNamt = 'DORGQR'
      CALL DORGQR(M,N,K,Q,Lda,Tau,Work,Lwork,info)
!
!     Copy R(1:n,1:k)
!
      CALL DLASET('Full',N,K,ZERO,ZERO,R,Lda)
      CALL DLACPY('Upper',N,K,Af,Lda,R,Lda)
!
!     Compute R(1:n,1:k) - Q(1:m,1:n)' * A(1:m,1:k)
!
      CALL DGEMM('Transpose','No transpose',N,K,M,-ONE,Q,Lda,A,Lda,ONE, &
     &           R,Lda)
!
!     Compute norm( R - Q'*A ) / ( M * norm(A) * EPS ) .
!
      anorm = DLANGE('1',M,K,A,Lda,Rwork)
      resid = DLANGE('1',N,K,R,Lda,Rwork)
      IF ( anorm>ZERO ) THEN
         Result(1) = ((resid/DBLE(MAX(1,M)))/anorm)/eps
      ELSE
         Result(1) = ZERO
      ENDIF
!
!     Compute I - Q'*Q
!
      CALL DLASET('Full',N,N,ZERO,ONE,R,Lda)
      CALL DSYRK('Upper','Transpose',N,M,-ONE,Q,Lda,ONE,R,Lda)
!
!     Compute norm( I - Q'*Q ) / ( M * EPS ) .
!
      resid = DLANSY('1','Upper',N,R,Lda,Rwork)
!
      Result(2) = (resid/DBLE(MAX(1,M)))/eps
!
!
!     End of DQRT02
!
      END SUBROUTINE DQRT02
