!*==zqrt02.f90  processed by SPAG 7.51RB at 20:37 on  3 Mar 2022
!> \brief \b ZQRT02
!
!  =========== DOCUMENTATION ===========
!
! Online html documentation available at
!            http://www.netlib.org/lapack/explore-html/
!
!  Definition:
!  ===========
!
!       SUBROUTINE ZQRT02( M, N, K, A, AF, Q, R, LDA, TAU, WORK, LWORK,
!                          RWORK, RESULT )
!
!       .. Scalar Arguments ..
!       INTEGER            K, LDA, LWORK, M, N
!       ..
!       .. Array Arguments ..
!       DOUBLE PRECISION   RESULT( * ), RWORK( * )
!       COMPLEX*16         A( LDA, * ), AF( LDA, * ), Q( LDA, * ),
!      $                   R( LDA, * ), TAU( * ), WORK( LWORK )
!       ..
!
!
!> \par Purpose:
!  =============
!>
!> \verbatim
!>
!> ZQRT02 tests ZUNGQR, which generates an m-by-n matrix Q with
!> orthonornmal columns that is defined as the product of k elementary
!> reflectors.
!>
!> Given the QR factorization of an m-by-n matrix A, ZQRT02 generates
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
!>          A is COMPLEX*16 array, dimension (LDA,N)
!>          The m-by-n matrix A which was factorized by ZQRT01.
!> \endverbatim
!>
!> \param[in] AF
!> \verbatim
!>          AF is COMPLEX*16 array, dimension (LDA,N)
!>          Details of the QR factorization of A, as returned by ZGEQRF.
!>          See ZGEQRF for further details.
!> \endverbatim
!>
!> \param[out] Q
!> \verbatim
!>          Q is COMPLEX*16 array, dimension (LDA,N)
!> \endverbatim
!>
!> \param[out] R
!> \verbatim
!>          R is COMPLEX*16 array, dimension (LDA,N)
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
!>          TAU is COMPLEX*16 array, dimension (N)
!>          The scalar factors of the elementary reflectors corresponding
!>          to the QR factorization in AF.
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
!> \ingroup complex16_lin
!
!  =====================================================================
      SUBROUTINE ZQRT02(M,N,K,A,Af,Q,R,Lda,Tau,Work,Lwork,Rwork,Result)
      IMPLICIT NONE
!*--ZQRT02138
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
      COMPLEX*16 A(Lda,*) , Af(Lda,*) , Q(Lda,*) , R(Lda,*) , Tau(*) ,  &
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
      EXTERNAL DLAMCH , ZLANGE , ZLANSY
!     ..
!     .. External Subroutines ..
      EXTERNAL ZGEMM , ZHERK , ZLACPY , ZLASET , ZUNGQR
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
!     Copy the first k columns of the factorization to the array Q
!
      CALL ZLASET('Full',M,N,ROGUE,ROGUE,Q,Lda)
      CALL ZLACPY('Lower',M-1,K,Af(2,1),Lda,Q(2,1),Lda)
!
!     Generate the first n columns of the matrix Q
!
      SRNamt = 'ZUNGQR'
      CALL ZUNGQR(M,N,K,Q,Lda,Tau,Work,Lwork,info)
!
!     Copy R(1:n,1:k)
!
      CALL ZLASET('Full',N,K,DCMPLX(ZERO),DCMPLX(ZERO),R,Lda)
      CALL ZLACPY('Upper',N,K,Af,Lda,R,Lda)
!
!     Compute R(1:n,1:k) - Q(1:m,1:n)' * A(1:m,1:k)
!
      CALL ZGEMM('Conjugate transpose','No transpose',N,K,M,DCMPLX(-ONE)&
     &           ,Q,Lda,A,Lda,DCMPLX(ONE),R,Lda)
!
!     Compute norm( R - Q'*A ) / ( M * norm(A) * EPS ) .
!
      anorm = ZLANGE('1',M,K,A,Lda,Rwork)
      resid = ZLANGE('1',N,K,R,Lda,Rwork)
      IF ( anorm>ZERO ) THEN
         Result(1) = ((resid/DBLE(MAX(1,M)))/anorm)/eps
      ELSE
         Result(1) = ZERO
      ENDIF
!
!     Compute I - Q'*Q
!
      CALL ZLASET('Full',N,N,DCMPLX(ZERO),DCMPLX(ONE),R,Lda)
      CALL ZHERK('Upper','Conjugate transpose',N,M,-ONE,Q,Lda,ONE,R,Lda)
!
!     Compute norm( I - Q'*Q ) / ( M * EPS ) .
!
      resid = ZLANSY('1','Upper',N,R,Lda,Rwork)
!
      Result(2) = (resid/DBLE(MAX(1,M)))/eps
!
!
!     End of ZQRT02
!
      END SUBROUTINE ZQRT02
